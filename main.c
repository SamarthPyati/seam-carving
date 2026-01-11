#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <omp.h>

#include <stb_image.h>
#include <stb_image_write.h>
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


// float matrix (for luminance, edge gradient & energy function)
typedef struct  
{
    float *items;
    int width, height, stride;
} Mat;

// image 
typedef struct 
{
    uint32_t *pixels;
    int width, height, stride;
}Img;


#define mat_at(mat, i, j) ((mat).items[(mat).stride * (i) + (j)])
#define img_at(img, i, j) ((img).pixels[(img).stride * (i) + (j)])
#define mat_within(mat, i, j) (0 <= i && i < (mat).height && 0 <= j && j < (mat).width)

#define RED (uint32_t)0xFF0000FF

static inline Mat mat_alloc(int width, int height)
{
    Mat mat = {0};
    mat.items = (float *) malloc(sizeof(*mat.items) * width * height);
    assert(mat.items != NULL);    
    mat.height = height;
    mat.width = width;
    mat.stride = width;
    return mat;
}

static inline Img img_alloc(int width, int height)
{
    Img img = {0};
    img.pixels = (uint32_t *) malloc(sizeof(*img.pixels) * width * height);
    assert(img.pixels != NULL);
    img.height = height;
    img.width = width;
    img.stride = width;
    return img;
}

static inline void mat_free(Mat *mat)
{
    free(mat->items);
    mat->items = NULL;
}

static inline void img_free(Img *img)
{
    free(img->pixels);
    img->pixels = NULL;
}

typedef struct {
    uint8_t r;
    uint8_t g;
    uint8_t b;
    uint8_t a;
} color_t;


static color_t rgb_to_color_t(uint32_t rgb)
{
    uint8_t r = (rgb >> (8 * 0)) & 0xFF;
    uint8_t g = (rgb >> (8 * 1)) & 0xFF;
    uint8_t b = (rgb >> (8 * 2)) & 0xFF;
    uint8_t a = (rgb >> (8 * 3)) & 0xFF;
    return (color_t) {.r = r, .g = g, .b = b, .a = a};
}

static void print_color_t(color_t c)
{   
    if (c.a == (uint8_t)-1)    printf("rgb(%u, %u, %u)\n", c.r, c.g, c.b);
    else printf("rgb(%u, %u, %u, %u)\n", c.r, c.g, c.b, c.a);
}

static float rgb_to_lum(uint32_t rgb)
{   
    // REFR: https://stackoverflow.com/questions/596216/formula-to-determine-perceived-brightness-of-rgb-color
    // Luminance is normalized 
    float r = ((rgb >> (8 * 0)) & 0xFF) / 255.0;
    float g = ((rgb >> (8 * 1)) & 0xFF) / 255.0;
    float b = ((rgb >> (8 * 2)) & 0xFF) / 255.0;
    return (0.2126 * r + 0.7152 * g + 0.0722 * b);
}

static void min_and_max(Mat *values, float *min, float *max, int l, int r)
{   
    // Divide and conquer approach
    if (l == r) // only element present
    {
        *min = values->items[l];
        *max = values->items[l];
    }
    else if (l == r - 1)  // two elements present 
    {
        *max = (values->items[l] > values->items[r]) ? values->items[l] : values->items[r]; 
        *min = (values->items[l] < values->items[r]) ? values->items[l] : values->items[r]; 
    }
    else 
    {   
        // Divide the arrays into 2 
        float lmax, lmin, rmax, rmin;
        int mid = l + (r - l) / 2;
        min_and_max(values, &lmin, &lmax, l, mid);
        min_and_max(values, &rmin, &rmax, mid + 1, r);
        *max = (lmax > rmax) ? lmax : rmax;
        *min = (lmin < rmin) ? lmin : rmin;
    }
}

static void analyse_min_and_max(const char *prompt, Mat *values)
{   
    assert(values != NULL);
    float min, max;
    min_and_max(values, &min, &max, 0, values->width * values->height - 1);
    printf("%s: (MIN -> %f & MAX -> %f)\n", prompt, min, max);
}

static void normalize_pixels(Mat *values)
{
    int count = values->width * values->height;
    if (count <= 0) return;

    float min, max;
    min_and_max(values, &min, &max, 0, count - 1);

    if (max == min) {
        // Identical elements, avoid division by zero and instead set all to zero
        memset(values->items, 0.0f, sizeof(*values->items) * count);
        return;
    }

    for (int i = 0; i < count; ++i)
    {
        values->items[i] = (values->items[i] - min) / (max - min);
    }
} 

static bool dump_mat(const char *fp, Mat mat)
{
    uint32_t *pixels = NULL;
    bool result = true;

    pixels = malloc(sizeof(*pixels) * mat.width * mat.height);
    assert(pixels != NULL);

    normalize_pixels(&mat);
    for (int y = 0; y < mat.height; ++y) 
    {
        for (int x = 0; x < mat.width; ++x) 
        {
            int i = y * mat.width + x;
            uint32_t value = 255 * mat_at(mat, y, x);
            pixels[i] = 0xFF000000 | (value << (8*2)) | (value << (8*1)) | (value << (8*0));
        }
    }

    if (!stbi_write_png(fp, mat.width, mat.height, 4, pixels, mat.width * sizeof(*pixels))) 
    {
        LOG_ERROR("ERROR: could not save file %s", fp);
        nob_return_defer(false);
    }
    else LOG_ERROR("OK: Saved the image file at %s\n", fp);

defer:
    free(pixels);
    return result;
}

static void dump_img(const char *fp, Img img)
{       
    // Dumps the normalized images i.e with float values from 0.0 to 1.0 to a file path 
    assert(img.pixels != NULL);

    if (!stbi_write_png(fp, img.width, img.height, 4, img.pixels, img.stride * sizeof(uint32_t)))  
    {
        LOG_ERROR_AND_ABORT("ERROR: Could not write the image file %s\n", fp);
    }
    else LOG_ERROR("OK: Saved the image file at %s\n", fp);
}

static void luminance(Img img, Mat lum)
{   
    assert(lum.items != NULL);
    assert(lum.width == img.width);
    assert(lum.height == img.height);

    for (int y = 0; y < lum.height; ++y)
    {
        for (int x = 0; x < lum.width; ++x) {
            mat_at(lum, y, x) = rgb_to_lum(img_at(img, y, x));
        }
    }
}

const static float gx[3][3] = {
    {1.0, 0.0, -1.0},
    {2.0, 0.0, -2.0},
    {1.0, 0.0, -1.0}
};

const static float gy[3][3] = {
    { 1.0,  2.0,  1.0},
    { 0.0,  0.0,  0.0},
    {-1.0, -2.0, -1.0}
};

static float sobel_at(Mat lum, int cy, int cx)
{
    float sx = 0.0f;
    float sy = 0.0f;
    for (int ky = -1; ky <= 1; ++ky)
    {
        for (int kx = -1; kx <= 1; ++kx)
        {
            int x = cx + kx;
            int y = cy + ky;
            // ensure x & y stay in bounds of image dimensions
            float l = (0 <= x && x < lum.width && 0 <= y && y < lum.height) ? mat_at(lum, y, x) : 0.0f;
            // kx and ky start from -1 so offset by 1
            sx += gx[ky + 1][kx + 1] * l;
            sy += gy[ky + 1][kx + 1] * l;
        }
    }
    return sqrtf(sx * sx + sy * sy);
}

static void sobel(Mat lum, Mat grad) 
{       
    // Convolutional kernel for SOBEL filter 
    // REFR: https://en.wikipedia.org/wiki/Sobel_operator
    assert(lum.items != NULL);
    assert(lum.width == grad.width);
    assert(lum.height == grad.height);

    // Convolving the kernel and the image
    #pragma omp parallel for
    for (int cy = 0; cy < lum.height; ++cy)
    {
        #pragma omp parallel for
        for (int cx = 0; cx < lum.width; ++cx)
        {
            // magnitude of the gradient 
            mat_at(grad, cy, cx) = sobel_at(lum, cy, cx);
        }
    }
}   

static void energy(Mat grad, Mat energy)
{   
    assert(grad.items != NULL);
    assert(grad.width == energy.width);
    assert(grad.height == energy.height);

    // Calculate the energy 
    #pragma omp parallel for 
    for (int x = 0; x < grad.width; ++x) 
    {   
        // Copy the first row as it is, as there's no row above to calculate the energy
        mat_at(energy, 0, x) = mat_at(grad, 0, x);
    }
    
    #pragma omp parallel for 
    for (int cy = 1; cy < grad.height; ++cy)
    {
        #pragma omp parallel for 
        for (int cx = 0; cx < grad.width; ++cx)
        {       
            // WORKING: M(i, j) = e(i, j) + min(M(i − 1, j − 1), M(i − 1, j), M(i − 1, j + 1))
            // TEST: M(i, j) = e(i, j) + min(M(i − 1, j - 2), M(i − 1, j − 1), M(i − 1, j), M(i − 1, j + 1), M(i − 1, j + 2))
            
            float min = FLT_MAX;
            for (int dx = -1; dx <= 1; ++dx)
            {
                int x = cx + dx;
                // float value = (0 <= x && x < width) ? energy[(cy - 1) * width + x] : 0.0f;
                float value = (0 <= x && x < grad.width) ? mat_at(energy, (cy - 1), x) : FLT_MAX;
                // Why FLT_MAX?: to not make it seam as minimum energy seam
                if (value < min)        min = value;
            }
            // energy[cy * width + cx] = grad[cy * width + cx] + min;
            mat_at(energy, cy, cx) = mat_at(grad, cy, cx) + min;
        }
    }
}

#if 0
    const static float bk3[3][3] = {       // blur kernel (3 x 3) gaussian blur
        {1, 2, 1},
        {2, 4, 2},
        {1, 2, 1}
    };  // (/ 16) divide every element by 16
#else 
    static float bk5[5][5] = {      // blur kernel (5 x 5) gaussian blur (better approx)
        {1.0, 4.0 , 6.0 , 4.0 , 1.0},
        {4.0, 16.0, 24.0, 16.0, 4.0},
        {6.0, 24.0, 36.0, 24.0, 6.0},
        {4.0, 16.0, 24.0, 16.0, 4.0},
        {1.0, 4.0 , 6.0 , 4.0 , 1.0}
    };  // (/ 256) divide every element by 256
#endif 

static void gaussianBlur(Mat lum, Mat img)
{   
    
    int kdim = 5;                  // kernel dim in order of ranges of 2k + 1, -(2k + 1)

    // normalize_filter(&bk5[0][0], 5);

    int width = lum.width;
    int height = lum.height;

    for (int cy = 0; cy < height; ++cy)
    {
        for (int cx = 0; cx < width; ++cx)
        {
            float acc = 0.0f;
            for (int ky = -((kdim - 1) / 2); ky <= ((kdim - 1) / 2); ++ky)
            {
                for (int kx = -((kdim - 1) / 2); kx <= ((kdim - 1) / 2); ++kx)
                {
                    int x = cx + kx;
                    int y = cy + ky;
                    // ensure x & y stay in bounds of image dimensions
                    float l = (0 <= x && x < width && 0 <= y && y < height) ? mat_at(lum, y, x) : 0.0f;
                    // kx and ky run from -2..2 for kdim=5 so offset by 2
                    acc += bk5[ky + 2][kx + 2] * l;
                }   
            }
            mat_at(img, cy, cx) = acc;
        }
    }
}

static void remove_pixel_img(Img img, int r, int c)
{   
    // Remove a particular pixel from image ((r, c) pixel coordinate)
    uint32_t *pixels_row = &img_at(img, r, 0);   
    memmove(pixels_row + c, pixels_row + c + 1, ((img.width - c - 1) * sizeof(uint32_t)));      
}

static void remove_pixel_mat(Mat mat, int r, int c)
{   
    // Remove a particular pixel from image ((r, c) pixel coordinate)
    float *pixels_row = &mat_at(mat, r, 0);   
    memmove(pixels_row + c, pixels_row + c + 1, ((mat.width - c - 1) * sizeof(float)));      
}

static void compute_seam(Mat egy, int *seams)
{
    // --- What are we doing here? --- 
    // Iteration the rows from top to bottom and find the lowest energy pixel   
    // or just pick a random seam from the vicinity
    // i.e, the three neighbouring images on top of the low energy pixel 

    // 'abcdefghijklm' -> image row , lets say 'f' is seam 
    // - pixels_row + seam points to 'f'
    // - pixels_row + seam + 1 points to next pixel of 'f'
    // - copying entire array of &(pixels_row + seam + 1) to 1 step back and then reduce width

    // Finding the minimum energy seam or selecting a random one
    int y = egy.height - 1;  
    seams[y] = 0;
    for (int x = 0; x < egy.width; ++x) 
    {
        if (mat_at(egy, y, x) < mat_at(egy, y, seams[y])) seams[y] = x;
    }

    for (y = egy.height - 2; y >= 0; --y)
    {           
        // Since we are calculating this from bottom to top 
        // We add 1 to get the previous seam
        
    //  current_seam = previous_seam
        seams[y]     = seams[y + 1];
        for (int dx = -1; dx <= 1; ++dx)
        {
            int x = seams[y] + dx;                 
            if (0 <= x && x < egy.width && mat_at(egy, y, x) < mat_at(egy, y, seams[y]))
            {   
                seams[y] = x;
            }
        }
    }
}

// static void remove_seam(Img img, Mat lum, Mat egy, int *seams)
// {                
    // compute_seam(egy, seams); 
    // for (int y = 0; y < img.height; ++y)
    // {
    //     remove_pixel_img(img, y, seams[y]);
    //     remove_pixel_mat(lum, y, seams[y]);
    // }
// }

void markout_sobel_patches(Mat grad, int *seams)
{
    // -- Optimized the process of removing seams -- 
    // mark the pixels around the seam
    #pragma omp parallel for 
    for (int cy = 0; cy < grad.height; ++cy)
    {
        int cx = seams[cy];
        for (int dy = -1; dy <= 1; ++dy)
        {
            for (int dx = -1; dx <= 1; ++dx)
            {
                int x = cx + dx;
                int y = cy + dy;
                if (0 <= x && x < grad.width && 0 <= y && y < grad.height) {
                    mat_at(grad, y, x) = NAN;
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{    
    if (argc != 3)
    {
        LOG_USAGE(argv[0], "<image_file_path> <seams_to_remove>");
    }

    const char *file_path = argv[1];
    size_t seams_to_remove = atoi(argv[2]);

    int width, height;
    uint32_t *pixels = (uint32_t *) stbi_load(file_path, &width, &height, NULL, 4);
    if (!pixels)  
    {   
        LOG_ERROR("ERROR: Failed to load image %s\n", file_path);
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
    // stores the index of seam in every row, so any index of this array represent the row of image
    int *seams = malloc(sizeof(*seams) * height);

    luminance(img, lum);
    sobel(lum, grad);   
    energy(grad, egy);

    for (size_t i = 0; i < seams_to_remove; ++i) {
        // Compute energy at every iteration
        // sobel(lum, grad);
        energy(grad, egy);
        
        // Compute the low energy seams or the indices of the seam for each row & remove them 
        compute_seam(egy, seams);
        markout_sobel_patches(grad, seams);
        
        for (int cy = 0; cy < grad.height; ++cy) {
            // remove the seam 
            int cx = seams[cy];
            remove_pixel_mat(lum, cy, cx);
            remove_pixel_mat(grad, cy, cx);
            remove_pixel_img(img, cy, cx);
        }        


        // reduce the width for each removed seam 
        img.width  -= 1;
        lum.width  -= 1;
        grad.width -= 1;
        egy.width  -= 1;  

        // recompute the sobel only for the marked seam 
        #pragma omp parallel for 
        for (int cy = 0; cy < grad.height; ++cy) {
            int s = seams[cy];
            /* recompute to the right of seam until a non-NaN value is found */
            for (int cx = s; cx < grad.width; ++cx) {
                if (!isnan(mat_at(grad, cy, cx))) break;
                mat_at(grad, cy, cx) = sobel_at(lum, cy, cx);
            }
            /* recompute to the left of seam until a non-NaN value is found */
            for (int cx = s - 1; cx >= 0; --cx) {
                if (!isnan(mat_at(grad, cy, cx))) break;
                mat_at(grad, cy, cx) = sobel_at(lum, cy, cx);
            }
        }
    }   

    // dump_mat("gradient.png", grad);
    // dump_mat("energy.png", egy);
    dump_img("output.png", img);

    // analyse_min_and_max("Gradient", &grad);
    // analyse_min_and_max("Luminance", &lum);

    stbi_image_free(pixels);
    mat_free(&lum);
    mat_free(&egy);
    mat_free(&grad);

    return 0;
}
