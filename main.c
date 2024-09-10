#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <math.h>

#define STB_IMAGE_IMPLEMENTATION
#include "external/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION 
#include "external/stb_image_write.h"


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


#define mat_at(mat, i, j) ((mat).items[(mat).width * i + j])
#define img_at(img, i, j) ((img).pixels[(img).width * i + j])

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

float *get_lum(uint32_t* pixels, int width, int height)
{   
    float *luminance = (float *)malloc (width * height * sizeof(*luminance));
    assert(luminance != NULL);

    for (int i = 0; i < height * width; ++i)
    {
        luminance[i] = rgb_to_lum(pixels[i]);
    }
    return luminance;
}

void min_and_max(float *values, int width, int height, float *min, float *max)
{
    *min = FLT_MAX;
    *max = FLT_MIN;
    for (size_t i = 0; i < (size_t)(width * height); ++i)
    {
        if (values[i] < *min) *min = values[i];
        if (values[i] > *max) *max = values[i];
    }
}

void analyse_min_and_max(const char *prompt, float *values, int width, int height)
{   
    assert(values != NULL);
    float min, max;
    min_and_max(values, width, height, &min, &max);
    printf("%s: (MIN -> %f & MAX -> %f)\n", prompt, min, max);
}

void normalize_pixels(float *values, int width, int height)
{
    float min, max;
    min_and_max(values, width, height, &min, &max);
    for (size_t i = 0; i < (size_t)(width * height); ++i)
    {
        values[i] = (values[i] - min) / (max - min);
    }
}

uint32_t *dump_image_dn(const char *fp, float *pixel_values, int width, int height)
{       
    // Dumps the normalized images i.e with float values from 0.0 to 1.0 to a file path 
    uint32_t *pixels = (uint32_t *)malloc(sizeof(*pixels) * width * height);
    assert(pixels != NULL);

    // normalize the image (if not normalized)
    normalize_pixels(pixel_values, width, height);

    for (size_t i = 0; i < (size_t)(width * height); ++i)
    {
        uint32_t val = (uint32_t)(pixel_values[i] * 255.0f);
        pixels[i] = (pixels[i] & 0xFF000000) | val << (8 * 2) | val << (8 * 1) | val << (8 * 0);
    }

    if (!stbi_write_jpg(fp, width, height, 4, pixels, 100))  
    {
        fprintf(stderr, "ERROR: Could not write the image file %s\n", fp);
        exit(EXIT_FAILURE);
    }
    else fprintf(stderr, "OK: Saved the image file at %s\n", fp);
    return pixels;
}

void dump_image(const char *fp, uint32_t *pixel_values, int width, int height)
{       
    // Dumps the normalized images i.e with float values from 0.0 to 1.0 to a file path 
    assert(pixel_values != NULL);

    if (!stbi_write_jpg(fp, width, height, 4, pixel_values, 100))  
    {
        fprintf(stderr, "ERROR: Could not write the image file %s\n", fp);
        exit(EXIT_FAILURE);
    }
    else fprintf(stderr, "OK: Saved the image file at %s\n", fp);
}



float *sobel(float *luminance, int width, int height) 
{       
    // Convolutional kernel for SOBEL filter 
    // REFR: https://en.wikipedia.org/wiki/Sobel_operator
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
    float *grad = (float *)malloc(sizeof(*grad) * width * height);
    assert(grad != NULL);

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
                    float l = (0 <= x && x < width && 0 <= y && y < width) ? luminance[y * width + x] : 0.0f;
                    // kx and ky start from -1 so offset by 1
                    sx += gx[ky + 1][kx + 1] * l;
                    sy += gy[ky + 1][kx + 1] * l;
                }
            }
            // magnitude of the gradient 
            grad[cy * width + cx] = sqrtf(sx * sx + sy * sy);
        }
    }
    return grad;
}

float *energy(float *grad, int width, int height)
{   
    float *energy = (float *) malloc(sizeof(*energy) * width * height);
    assert(energy != NULL);
    
    // Calculate the energy 
    for (int x = 0; x < width; ++x) 
    {   
        // Copy the first row as it is, as there's no row above to calculate the energy
        energy[0 * width + x] = grad[0 * width + x];
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
                float value = (0 <= x && x < width) ? energy[(cy - 1) * width + x] : 0.0f;
                if (value < min)        min = value;
            }
            energy[cy * width + cx] = grad[cy * width + cx] + min;
        }
    }
    return energy;
}

int main(int argc, char *argv[])
{   
    if (argc != 2)
    {
        fprintf(stderr, "USAGE: <c_exec_file> <image_path (.jpgs file only)>\n");
        exit(EXIT_FAILURE);
    }

    const char *inp_img_fp = argv[1];

    int width, height;
    uint32_t *pixels = (uint32_t *) stbi_load(inp_img_fp, &width, &height, NULL, 4);
    if (!pixels)  
    {
        fprintf(stderr, "ERROR: Failed to load image %s\n", inp_img_fp);
        exit(EXIT_FAILURE);
    }

    float *luminance = get_lum(pixels, width, height);
    dump_image_dn("luminance.jpg", luminance, width, height);

    float *grad = sobel(luminance, width, height);
    dump_image_dn("gradient.jpg", grad, width, height);
    
    float *energy_map = energy(luminance, width, height);
    dump_image_dn("energy.jpg", energy_map, width, height);

#if 0
    // Find the lowest energy column from the bottom most row 
    int minCol = 0;
    for (int x = 0; x < width; ++x)
    {
        if (energy_map[(height - 1) * width + x] < energy_map[(height - 1) * width + minCol]) minCol = x;
    }   
    pixels[(height - 1) * width + minCol] = red;
#else 
    int minCol = width / 2 ;        // column with minimun energy density
#endif

    // Iteration the rows from top to bottom and find the lowest energy pixel from the vicinity
    // i.e, the three neighbouring images on top of the low energy pixel 
    for (size_t i = 0; i < 100; ++i) {  // Calculate 100 seams at equal intervals
        minCol = (width / 100) * i;
        for (int y = height - 2; y >= 0; --y)
        {   
            for (int dx = -1; dx <= 1; ++dx)
            {
                int x = minCol + dx;    // find the appropriate col
                if (0 <= x && x < width && energy_map[y * width + x] < energy_map[y * width + minCol])
                {
                    minCol = x;
                }
            }
            pixels[y * width + minCol] = RED;
        }
    }

    dump_image("image_with_seam.jpg", pixels, width, height);
    
    stbi_image_free(pixels);
    free(luminance);
    free(energy_map);
    free(grad);

    return 0;
}