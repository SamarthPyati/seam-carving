#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>

#include "external/stb_image.h"
#include "external/stb_image_write.h"
#include "external/nob.h"


/* 
    -- NOTE: -- 
    Color Information stored in uint32_t
    - > [BBBBBBBB][GGGGGGGG][RRRRRRRR][AAAAAAAA] or 0xAARRGGBB

    Extraction of each:
    RED:    (color >> 16) & 0xFF;
    GREEN:  (color >> 8)  & 0xFF;
    BLUE:   (color >> 0)  & 0xFF;
    
    uint32_t color = (a << 24) | (r << 16) | (g << 8) | b;
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


#define mat_at(mat, i, j) ((mat).items[(mat).width * (i) + (j)])
#define img_at(img, i, j) ((img).pixels[(img).width * (i) + (j)])

Mat mat_alloc(int width, int height)
{
    Mat m = {0};
    m.items = (float *) malloc(sizeof(*m.items) * width * height);
    assert(m.items != NULL);    
    m.height = height;
    m.width = width;
    m.stride = width;
    return m;
}

Img img_alloc(int width, int height)
{
    Img im = {0};
    im.pixels = (uint32_t *) malloc(sizeof(*im.pixels) * width * height);
    assert(im.pixels != NULL);
    im.height = height;
    im.width = width;
    im.stride = width;
    return im;
}

void mat_free(Mat *mat)
{
    free(mat->items);
}

void img_free(Img *img)
{
    free(img->pixels);
}


#define SPDA_ARR_LEN(xs) sizeof(xs) / sizeof(xs[0])
#define raise(err, type)                                        \
    do {                                                        \
        fprintf(stderr, "ERROR {%s} : %s\n", err, type);        \
    } while (0)

#define RED (uint32_t)0xFF0000FF;

typedef struct {
    uint8_t r;
    uint8_t g;
    uint8_t b;
    uint8_t a;
} color_t;


color_t rgb_to_color_t(uint32_t rgb)
{
    uint8_t r = (rgb >> (8 * 0)) & 0xFF;
    uint8_t g = (rgb >> (8 * 1)) & 0xFF;
    uint8_t b = (rgb >> (8 * 2)) & 0xFF;
    return (color_t) {.r = r, .g = g, .b = b};
}

void print_color_t(color_t c)
{   
    if (c.a == (uint8_t)-1)    printf("rgb(%u, %u, %u)\n", c.r, c.g, c.b);
    else printf("rgb(%u, %u, %u, %u)\n", c.r, c.g, c.b, c.a);
}

float rgb_to_lum(uint32_t rgb)
{   
    // REFR: https://stackoverflow.com/questions/596216/formula-to-determine-perceived-brightness-of-rgb-color
    // Luminance is normalized 
    float r = ((rgb >> (8 * 0)) & 0xFF) / 255.0;
    float g = ((rgb >> (8 * 1)) & 0xFF) / 255.0;
    float b = ((rgb >> (8 * 2)) & 0xFF) / 255.0;
    return (0.2126 * r + 0.7152 * g + 0.0722 * b);
}

void luminance(Img img, Mat lum)
{   
    assert(lum.items != NULL);
    assert(lum.width == img.width);
    assert(lum.height == img.height);

    for (int y = 0; y < img.height; ++y)
    {
        for (int x = 0; x < img.width; ++x) {
            mat_at(lum, y, x) = rgb_to_lum(img_at(img, y, x));
        }
    }
}

void min_and_max(Mat *values, float *min, float *max)
{
    *min = FLT_MAX;
    *max = FLT_MIN;
    for (int i = 0; i < (values->width * values->height); ++i)
    {
        if (values->items[i] < *min) *min = values->items[i];
        if (values->items[i] > *max) *max = values->items[i];
    }
}

void analyse_min_and_max(const char *prompt, Mat *values)
{   
    assert(values != NULL);
    float min, max;
    min_and_max(values, &min, &max);
    printf("%s: (MIN -> %f & MAX -> %f)\n", prompt, min, max);
}

void normalize_pixels(Mat *values)
{
    float min, max;
    min_and_max(values, &min, &max);
    for (int i = 0; i < values->width * values->height; ++i)
    {
        values->items[i] = (values->items[i] - min) / (max - min);
    }
}

bool dump_mat(const char *fp, Mat mat)
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
        fprintf(stderr, "ERROR: could not save file %s", fp);
        nob_return_defer(false);
    }
    else fprintf(stderr, "OK: Saved the image file at %s\n", fp);

defer:
    free(pixels);
    return result;
}

void dump_img(const char *fp, Img img)
{       
    // Dumps the normalized images i.e with float values from 0.0 to 1.0 to a file path 
    assert(img.pixels != NULL);

    if (!stbi_write_png(fp, img.width, img.height, 4, img.pixels, img.stride * sizeof(*img.pixels)))  
    {
        fprintf(stderr, "ERROR: Could not write the image file %s\n", fp);
        exit(EXIT_FAILURE);
    }
    else fprintf(stderr, "OK: Saved the image file at %s\n", fp);
}

void sobel(Mat lum, Mat grad) 
{       
    // Convolutional kernel for SOBEL filter 
    // REFR: https://en.wikipedia.org/wiki/Sobel_operator
    assert(lum.items != NULL);
    assert(lum.width == grad.width);
    assert(lum.height == grad.height);

    static float gx[3][3] = {
        {1, 0, -1},
        {2, 0, -2},
        {1, 0, -1}
    };

    static float gy[3][3] = {
        {1, 2, 1},
        {0, 0, 0},
        {-1, -2, -1}
    };

    // gradient or edge mapping of image
    int width = lum.width;
    int height = lum.height;
    
    // Convolving the kernel and the image
    for (int cy = 0; cy < height; ++cy)
    {
        for (int cx = 0; cx < width; ++cx)
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
                    float l = (0 <= x && x < width && 0 <= y && y < height) ? mat_at(lum, y, x) : 0.0f;
                    // kx and ky start from -1 so offset by 1
                    sx += gx[ky + 1][kx + 1] * l;
                    sy += gy[ky + 1][kx + 1] * l;
                }
            }
            // magnitude of the gradient 
            mat_at(grad, cy, cx) = sqrtf(sx * sx + sy * sy);
        }
    }
}   

void energy(Mat grad, Mat energy)
{   
    assert(grad.items != NULL);
    assert(grad.width == energy.width);
    assert(grad.height == energy.height);

    int width = energy.width;
    int height = energy.height;
    // Calculate the energy 
    for (int x = 0; x < width; ++x) 
    {   
        // Copy the first row as it is, as there's no row above to calculate the energy
        mat_at(energy, 0, x) = mat_at(grad, 0, x);
    }
    
    for (int cy = 1; cy < height; ++cy)
    {
        for (int cx = 0; cx < width; ++cx)
        {       
            // M(i, j) = e(i, j) + min(M(i − 1, j − 1), M(i − 1, j), M(i − 1, j + 1))
            float min = FLT_MAX;
            for (int dx = -1; dx <= 1; ++dx)
            {
                int x = cx + dx;
                // float value = (0 <= x && x < width) ? energy[(cy - 1) * width + x] : 0.0f;
                float value = (0 <= x && x < width) ? mat_at(energy, (cy - 1), x) : FLT_MAX;
                // Why FLT_MAX?: to not make it seam as minimum energy seam
                if (value < min)        min = value;
            }
            // energy[cy * width + cx] = grad[cy * width + cx] + min;
            mat_at(energy, cy, cx) = mat_at(grad, cy, cx) + min;
        }
    }
}

Mat gaussianBlur(Mat lum)
{   
    Mat blur = mat_alloc(lum.width, lum.height);

#if 0
    static float bk3[3][3] = {       // blur kernel (3 x 3) gaussian blur
        {1, 2, 1},
        {2, 4, 2},
        {1, 2, 1}
    };  // (/ 16) divide every element by 16
#else 
    static float bk5[5][5] = {      // blur kernel (5 x 5) gaussian blur (better approx)
        {1, 4 , 6 , 4 , 1},
        {4, 16, 24, 16, 4},
        {6, 24, 36, 24, 6},
        {4, 16, 24, 16, 4},
        {1, 4 , 6 , 4 , 1}
    };  // (/ 256) divide every element by 256
#endif 

    int kdim = 5;               // kernel dim in order of ranges of 2k + 1, -(2k + 1)

    for (int i = 0; i < kdim; ++i) {
        for (int j = 0; j < kdim; ++j) {
            bk5[i][j] = bk5[i][j] / 256;
        }
    }

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
                    float l = (0 <= x && x < width && 0 <= y && y < width) ? mat_at(lum, y, x) : 0.0f;
                    // kx and ky start from -1 so offset by 1
                    acc += bk5[ky + 1][kx + 1] * l;
                }   
            }
            mat_at(blur, cy, cx) = acc;
        }
    }
    return blur;
}


// Implementing bit masking 
typedef enum {
    DUMP_DUMP_IMg    = 0x00000400,  // 1 << 10
    DUMP_LUMINANCE     = 0x00000001,  // 1 << 0
    DUMP_GRADIENT      = 0x00000002,  // 1 << 1
    DUMP_ENERGY        = 0x00000004,  // 1 << 2
    DUMP_GAUSSIAN_BLUR = 0x00000008,  // 1 << 3
    DUMP_COUNTS        = 0x00000010   // 1 << 4
} dump_modes;

int main(int argc, char *argv[])
{       
   if (argc != 2)
    {
        fprintf(stderr, "USAGE: <c_exec_file> <image_path (.jpgs file only)>\n");
        exit(EXIT_FAILURE);
    }

    const char *inp_img_fp = argv[1];

    int _width, _height;
    uint32_t *pixels = (uint32_t *) stbi_load(inp_img_fp, &_width, &_height, NULL, 4);
    if (!pixels)  
    {
        fprintf(stderr, "ERROR: Failed to load image %s\n", inp_img_fp);
        exit(EXIT_FAILURE);
    }

    Img img = {
        .height = _height,
        .width = _width,
        .pixels = pixels,
        .stride = _width
    };

    Mat lum = mat_alloc(_width, _height);
    Mat grad = mat_alloc(_width, _height);
    Mat egy = mat_alloc(_width, _height);  

    // Iteration the rows from top to bottom and find the lowest energy pixel 
    // or just pick a random seam from the vicinity
    // i.e, the three neighbouring images on top of the low energy pixel 

    // What are we doing here?
    // 'abcdefghijklm' -> image row , lets say 'f' is seam 
    // - pixels_row + seam points to 'f'
    // - pixels_row + seam + 1 points to next pixel of 'f'
    // - copying entire array of &(pixels_row + seam + 1) to 1 step back and then reduce width

    size_t n_iters = 100;
    for (size_t i = 0; i < n_iters; ++i) {
        printf("Removing Seam %zu\n", i);
       
        // Update these at every iteration 
        luminance(img, lum);
        sobel(lum, grad);
        energy(grad, egy);

        // Finding the minimum energy seam or selecting a random one
        int seam = (rand() % img.width);
        int y = img.height - 1;  

        uint32_t *pixels_row = &img_at(img, y, 0);   // get the last row of the image
        memmove(pixels_row + seam, pixels_row + seam + 1, ((img.width - seam - 1) * sizeof(uint32_t)));       // removing the seam 

        for (y = img.height - 2; y >= 0; --y)
        {   
            for (int dx = -1; dx <= 1; ++dx)
            {
                int x = seam + dx;                 // find the appropriate seam col
                if (0 <= x && x < img.width && mat_at(egy, y, x) < mat_at(egy, y, seam))
                {   
                    seam = x;
                }
            }
            pixels_row = &img_at(img, y, 0);                                                         // get the row of the image
            memmove(pixels_row + seam, pixels_row + seam + 1, ((img.width - seam - 1) * sizeof(uint32_t)));    // remove the seam 
            // seam_image.pixels[y * seam_image.width + seam] = RED;     // Visualize the seams
        }

        // Everytime we remove column reduce the width 
        img.width  -= 1;
        lum.width  -= 1;
        grad.width -= 1;
        egy.width  -= 1;  
    }   

    dump_img("seam.png", img);

    stbi_image_free(pixels);
    mat_free(&lum);
    mat_free(&egy);
    mat_free(&grad);

    return 0;
}

