#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#define NOB_IMPLEMENTATION
#include <nob.h>

#include <utils.h>

/*
    --- NOTE: ---
    Color layout used by stbi_load (bytes in memory): [R][G][B][A]
    On little-endian machines a `uint32_t` view of the 4 bytes becomes 0xAABBGGRR.

    Extraction of each channel (works on this platform):
    RED:    (color >> 0) & 0xFF;
    GREEN:  (color >> 8) & 0xFF;
    BLUE:   (color >> 16) & 0xFF;
    ALPHA:  (color >> 24) & 0xFF;

    uint32_t packed = (a << 24) | (b << 16) | (g << 8) | r;
*/

typedef struct
{
    float *items;
    int width;
    int height;
    int stride;
} Mat;

typedef struct
{
    uint32_t *pixels;
    int width;
    int height;
    int stride;
} Img;

#define mat_at(mat, i, j) ((mat).items[(mat).stride * (i) + (j)])
#define img_at(img, i, j) ((img).pixels[(img).stride * (i) + (j)])

#ifdef _OPENMP
#define SEAM_CARVING_OMP_PARALLEL_FOR _Pragma("omp parallel for schedule(static)")
#else
#define SEAM_CARVING_OMP_PARALLEL_FOR
#endif

static inline Mat mat_alloc(int width, int height)
{
    Mat mat = {0};
    mat.items = (float *)malloc(sizeof(*mat.items) * width * height);
    assert(mat.items != NULL);
    mat.height = height;
    mat.width = width;
    mat.stride = width;
    return mat;
}

static inline void mat_free(Mat *mat)
{
    free(mat->items);
    mat->items = NULL;
}

static float rgb_to_lum(uint32_t rgb)
{
    float r = ((rgb >> (8 * 0)) & 0xFF) / 255.0f;
    float g = ((rgb >> (8 * 1)) & 0xFF) / 255.0f;
    float b = ((rgb >> (8 * 2)) & 0xFF) / 255.0f;
    return (0.2126f * r + 0.7152f * g + 0.0722f * b);
}

/* Plain `<` instead of fminf(): fminf() has NaN-propagation semantics that stop
   GCC/Clang from folding it to a bare minss/vminps unless -ffinite-math-only is
   on. This form always lowers to a single min instruction and vectorises. */
static inline float fmin2(float a, float b) { return a < b ? a : b; }
static inline float fmin3(float a, float b, float c) { return fmin2(fmin2(a, b), c); }

static inline float lum_at(const Img img, int y, int x) {
    // Compute luminance at a image pixel (y, x)
    if (x < 0 || x >= img.width || y < 0 || y >= img.height) return 0.0f;
    return rgb_to_lum(img_at(img, y, x));
}

// https://en.wikipedia.org/wiki/Sobel_operator#Formulation
static float sobel_at(const Img img, int cy, int cx) {
    const float p00 = lum_at(img, cy - 1, cx - 1);
    const float p01 = lum_at(img, cy - 1, cx    );
    const float p02 = lum_at(img, cy - 1, cx + 1);
    const float p10 = lum_at(img, cy,     cx - 1);
    const float p12 = lum_at(img, cy,     cx + 1);
    const float p20 = lum_at(img, cy + 1, cx - 1);
    const float p21 = lum_at(img, cy + 1, cx    );
    const float p22 = lum_at(img, cy + 1, cx + 1);

    const float sx = (p00 + 2.0f * p10 + p20) - (p02 + 2.0f * p12 + p22);
    const float sy = (p00 + 2.0f * p01 + p02) - (p20 + 2.0f * p21 + p22);
    return sqrtf(sx * sx + sy * sy);
}

static void sobel(const Img img, Mat grad)
{
    assert(img.pixels != NULL);
    assert(img.width == grad.width);
    assert(img.height == grad.height);

    SEAM_CARVING_OMP_PARALLEL_FOR
    for (int cy = 0; cy < img.height; ++cy)
    {
        for (int cx = 0; cx < img.width; ++cx)
        {
            mat_at(grad, cy, cx) = sobel_at(img, cy, cx);
        }
    }
}

static void energy(const Mat grad, Mat energy)
{
    assert(grad.items != NULL);
    assert(grad.width == energy.width);
    assert(grad.height == energy.height);

    const int width = grad.width;
    const int height = grad.height;

    if (width <= 0 || height <= 0) return;
    memcpy(&mat_at(energy, 0, 0), &mat_at(grad, 0, 0), sizeof(float) * (size_t)width);

    for (int y = 1; y < height; ++y)
    {
        const float *restrict g = &mat_at(grad, y, 0);
        const float *restrict p = &mat_at(energy, y - 1, 0);
        float *restrict c = &mat_at(energy, y, 0);

        if (width == 1)
        {
            c[0] = g[0] + p[0];
            continue;
        }

        c[0] = g[0] + fmin2(p[0], p[1]);

        for (int x = 1; x < width - 1; ++x)
        {
            c[x] = g[x] + fmin3(p[x - 1], p[x], p[x + 1]);
        }

        c[width - 1] = g[width - 1] + fmin2(p[width - 2], p[width - 1]);
    }
}

static inline void remove_pixel_img(Img img, int r, int c)
{
    uint32_t *row = &img_at(img, r, 0);
    memmove(row + c, row + c + 1, (size_t)(img.width - c - 1) * sizeof(uint32_t));
}

static inline void remove_pixel_mat(Mat mat, int r, int c)
{
    float *row = &mat_at(mat, r, 0);
    memmove(row + c, row + c + 1, (size_t)(mat.width - c - 1) * sizeof(float));
}

static void compute_seam(const Mat egy, int *seams)
{
    const int width = egy.width;
    const int height = egy.height;
    if (width <= 0 || height <= 0) return;

    float *dp = (float *)malloc(sizeof(*dp) * width);
    int *parents = (int *)malloc(sizeof(*parents) * width * height);
    assert(dp != NULL);
    assert(parents != NULL);

    for (int x = 0; x < width; ++x)
    {
        dp[x] = mat_at(egy, 0, x);
        parents[x] = -1;
    }

    for (int y = 1; y < height; ++y)
    {
        const int row_offset = y * width;
        for (int x = 0; x < width; ++x)
        {
            int best_col = x;
            float best_cost = dp[x];

            if (x > 0 && dp[x - 1] < best_cost)
            {
                best_col = x - 1;
                best_cost = dp[x - 1];
            }
            if (x + 1 < width && dp[x + 1] < best_cost)
            {
                best_col = x + 1;
                best_cost = dp[x + 1];
            }

            dp[x] = best_cost + mat_at(egy, y, x);
            parents[row_offset + x] = best_col;
        }
    }

    int best_col = 0;
    float best_cost = dp[0];
    for (int x = 1; x < width; ++x)
    {
        if (dp[x] < best_cost)
        {
            best_cost = dp[x];
            best_col = x;
        }
    }

    for (int y = height - 1; y >= 0; --y)
    {
        seams[y] = best_col;
        if (y > 0)
        {
            best_col = parents[y * width + best_col];
        }
    }

    free(dp);
    free(parents);
}

static void markout_sobel_patches(Mat grad, const int *seams)
{
    SEAM_CARVING_OMP_PARALLEL_FOR
    for (int cy = 0; cy < grad.height; ++cy)
    {
        int cx = seams[cy];
        for (int dy = -1; dy <= 1; ++dy)
        {
            for (int dx = -1; dx <= 1; ++dx)
            {
                int x = cx + dx;
                int y = cy + dy;
                if (0 <= x && x < grad.width && 0 <= y && y < grad.height)
                {
                    mat_at(grad, y, x) = NAN;
                }
            }
        }
    }
}

static bool dump_img(const char *fp, const Img img)
{
    assert(img.pixels != NULL);

    if (!stbi_write_png(fp, img.width, img.height, 4, img.pixels, img.stride * sizeof(uint32_t)))
    {
        LOG_ERROR("ERROR: Could not write the image file %s\n", fp);
        return false;
    }

    LOG_ERROR("OK: Saved the image file at %s\n", fp);
    return true;
}

int seam_carving_main(int argc, char *argv[])
{
    if (argc != 3 && argc != 4)
    {
        LOG_USAGE(argv[0], "<input_file_path> <n_seams_to_remove> [output_file_path]");
        return EXIT_FAILURE;
    }

    const char *file_path = argv[1];
    size_t seams_to_remove = atoi(argv[2]);
    const char *output_path = (argc == 4) ? argv[3] : "output.png";

    int width = 0;
    int height = 0;
    uint32_t *pixels = (uint32_t *)stbi_load(file_path, &width, &height, NULL, 4);
    if (!pixels)
    {
        LOG_ERROR("ERROR: Failed to load image %s\n", file_path);
        return EXIT_FAILURE;
    }

    Img img = {0};
    img.pixels = pixels;
    img.width = width;
    img.height = height;
    img.stride = width;

    CLAMP_ASSIGN(seams_to_remove, 0UL, width - 1UL);

    Mat grad = mat_alloc(width, height);
    Mat egy = mat_alloc(width, height);
    int *seams = (int *)malloc(sizeof(*seams) * height);
    assert(seams != NULL);

    sobel(img, grad);
    energy(grad, egy);

    for (size_t i = 0; i < seams_to_remove; ++i)
    {
        energy(grad, egy);
        compute_seam(egy, seams);
        markout_sobel_patches(grad, seams);

        for (int cy = 0; cy < grad.height; ++cy)
        {
            int cx = seams[cy];
            remove_pixel_mat(grad, cy, cx);
            remove_pixel_img(img, cy, cx);
        }

        img.width -= 1;
        grad.width -= 1;
        egy.width -= 1;

        SEAM_CARVING_OMP_PARALLEL_FOR
        for (int cy = 0; cy < grad.height; ++cy)
        {
            int s = seams[cy];
            for (int cx = s; cx < grad.width; ++cx)
            {
                if (!isnan(mat_at(grad, cy, cx))) break;
                mat_at(grad, cy, cx) = sobel_at(img, cy, cx);
            }
            for (int cx = s - 1; cx >= 0; --cx)
            {
                if (!isnan(mat_at(grad, cy, cx))) break;
                mat_at(grad, cy, cx) = sobel_at(img, cy, cx);
            }
        }
    }

    if (!dump_img(output_path, img))
    {
        stbi_image_free(pixels);
        mat_free(&egy);
        mat_free(&grad);
        free(seams);
        return EXIT_FAILURE;
    }

    stbi_image_free(pixels);
    mat_free(&egy);
    mat_free(&grad);
    free(seams);

    return EXIT_SUCCESS;
}

int main(int argc, char **argv)
{
    return seam_carving_main(argc, argv);
}
