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

static void luminance(const Img img, Mat lum)
{
    assert(lum.items != NULL);
    assert(lum.width == img.width);
    assert(lum.height == img.height);

    SEAM_CARVING_OMP_PARALLEL_FOR
    for (int y = 0; y < lum.height; ++y)
    {
        for (int x = 0; x < lum.width; ++x)
        {
            mat_at(lum, y, x) = rgb_to_lum(img_at(img, y, x));
        }
    }
}

static const float gx[3][3] = {
    {1.0f, 0.0f, -1.0f},
    {2.0f, 0.0f, -2.0f},
    {1.0f, 0.0f, -1.0f},
};

static const float gy[3][3] = {
    {1.0f, 2.0f, 1.0f},
    {0.0f, 0.0f, 0.0f},
    {-1.0f, -2.0f, -1.0f},
};

static float sobel_at(const Mat lum, int cy, int cx)
{
    float sx = 0.0f;
    float sy = 0.0f;
    for (int ky = -1; ky <= 1; ++ky)
    {
        for (int kx = -1; kx <= 1; ++kx)
        {
            int x = cx + kx;
            int y = cy + ky;
            float l = (0 <= x && x < lum.width && 0 <= y && y < lum.height) ? mat_at(lum, y, x) : 0.0f;
            sx += gx[ky + 1][kx + 1] * l;
            sy += gy[ky + 1][kx + 1] * l;
        }
    }
    return sqrtf(sx * sx + sy * sy);
}

static void sobel(const Mat lum, Mat grad)
{
    assert(lum.items != NULL);
    assert(lum.width == grad.width);
    assert(lum.height == grad.height);

    SEAM_CARVING_OMP_PARALLEL_FOR
    for (int cy = 0; cy < lum.height; ++cy)
    {
        for (int cx = 0; cx < lum.width; ++cx)
        {
            mat_at(grad, cy, cx) = sobel_at(lum, cy, cx);
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
    const int stride = grad.stride;

    for (int x = 0; x < width; ++x)
    {
        mat_at(energy, 0, x) = mat_at(grad, 0, x);
    }

    for (int y = 1; y < height; ++y)
    {
        float *prev_row = &energy.items[(y - 1) * stride];
        float *cur_row = &energy.items[y * stride];
        const float *grad_row = &grad.items[y * stride];

        if (width == 1)
        {
            cur_row[0] = grad_row[0] + prev_row[0];
            continue;
        }

        cur_row[0] = grad_row[0] + ((prev_row[0] < prev_row[1]) ? prev_row[0] : prev_row[1]);

        for (int x = 1; x < width - 1; ++x)
        {
            float best_prev = prev_row[x - 1];
            if (prev_row[x] < best_prev) best_prev = prev_row[x];
            if (prev_row[x + 1] < best_prev) best_prev = prev_row[x + 1];
            cur_row[x] = grad_row[x] + best_prev;
        }

        cur_row[width - 1] = grad_row[width - 1] + ((prev_row[width - 2] < prev_row[width - 1]) ? prev_row[width - 2] : prev_row[width - 1]);
    }
}

static void remove_pixel_img(Img img, int r, int c)
{
    uint32_t *pixels_row = &img_at(img, r, 0);
    memmove(pixels_row + c, pixels_row + c + 1, ((img.width - c - 1) * sizeof(uint32_t)));
}

static void remove_pixel_mat(Mat mat, int r, int c)
{
    float *pixels_row = &mat_at(mat, r, 0);
    memmove(pixels_row + c, pixels_row + c + 1, ((mat.width - c - 1) * sizeof(float)));
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

    Mat lum = mat_alloc(width, height);
    Mat grad = mat_alloc(width, height);
    Mat egy = mat_alloc(width, height);
    int *seams = (int *)malloc(sizeof(*seams) * height);
    assert(seams != NULL);

    luminance(img, lum);
    sobel(lum, grad);
    energy(grad, egy);

    for (size_t i = 0; i < seams_to_remove; ++i)
    {
        energy(grad, egy);
        compute_seam(egy, seams);
        markout_sobel_patches(grad, seams);

        for (int cy = 0; cy < grad.height; ++cy)
        {
            int cx = seams[cy];
            remove_pixel_mat(lum, cy, cx);
            remove_pixel_mat(grad, cy, cx);
            remove_pixel_img(img, cy, cx);
        }

        img.width -= 1;
        lum.width -= 1;
        grad.width -= 1;
        egy.width -= 1;

        SEAM_CARVING_OMP_PARALLEL_FOR
        for (int cy = 0; cy < grad.height; ++cy)
        {
            int s = seams[cy];
            for (int cx = s; cx < grad.width; ++cx)
            {
                if (!isnan(mat_at(grad, cy, cx))) break;
                mat_at(grad, cy, cx) = sobel_at(lum, cy, cx);
            }
            for (int cx = s - 1; cx >= 0; --cx)
            {
                if (!isnan(mat_at(grad, cy, cx))) break;
                mat_at(grad, cy, cx) = sobel_at(lum, cy, cx);
            }
        }
    }

    if (!dump_img(output_path, img))
    {
        stbi_image_free(pixels);
        mat_free(&lum);
        mat_free(&egy);
        mat_free(&grad);
        free(seams);
        return EXIT_FAILURE;
    }

    stbi_image_free(pixels);
    mat_free(&lum);
    mat_free(&egy);
    mat_free(&grad);
    free(seams);

    return EXIT_SUCCESS;
}

int main(int argc, char **argv)
{
    return seam_carving_main(argc, argv);
}
