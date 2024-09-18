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


#define mat_at(mat, i, j) ((mat).items[(mat).stride * (i) + (j)])
#define img_at(img, i, j) ((img).pixels[(img).stride * (i) + (j)])

static Mat mat_alloc(int width, int height)
{
    Mat mat = {0};
    mat.items = (float *) malloc(sizeof(*mat.items) * width * height);
    assert(mat.items != NULL);    
    mat.height = height;
    mat.width = width;
    mat.stride = width;
    return mat;
}

static Img img_alloc(int width, int height)
{
    Img img = {0};
    img.pixels = (uint32_t *) malloc(sizeof(*img.pixels) * width * height);
    assert(img.pixels != NULL);
    img.height = height;
    img.width = width;
    img.stride = width;
    return img;
}

void mat_free(Mat *mat)
{
    free(mat->items);
}

void img_free(Img *img)
{
    free(img->pixels);
}

#define RED (uint32_t)0xFF0000FF;

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
    return (color_t) {.r = r, .g = g, .b = b};
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

static void min_and_max(Mat *values, float *min, float *max)
{
    *min = FLT_MAX;
    *max = FLT_MIN;
    for (int y = 0; y < values->height; ++y)
    {
        for (int x = 0; x < values->width; ++x)
        {
            if (mat_at(*values, y, x) < *min) *min = mat_at(*values, y, x);
            if (mat_at(*values, y, x) > *max) *max = mat_at(*values, y, x);
        }
    }
}

static void analyse_min_and_max(const char *prompt, Mat *values)
{   
    assert(values != NULL);
    float min, max;
    min_and_max(values, &min, &max);
    printf("%s: (MIN -> %f & MAX -> %f)\n", prompt, min, max);
}

static void normalize_pixels(Mat *values)
{
    float min, max;
    min_and_max(values, &min, &max);
    for (int i = 0; i < values->width * values->height; ++i)
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
        fprintf(stderr, "ERROR: could not save file %s", fp);
        nob_return_defer(false);
    }
    else fprintf(stderr, "OK: Saved the image file at %s\n", fp);

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
        fprintf(stderr, "ERROR: Could not write the image file %s\n", fp);
        exit(EXIT_FAILURE);
    }
    else fprintf(stderr, "OK: Saved the image file at %s\n", fp);
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

static float sobel_at(Mat lum, int cy, int cx)
{
    static float gx[3][3] = {
        {1.0, 0.0, -1.0},
        {2.0, 0.0, -2.0},
        {1.0, 0.0, -1.0}
    };

    static float gy[3][3] = {
        {1.0, 2.0, 1.0},
        {0.0, 0.0, 0.0},
        {-1.0, -2.0, -1.0}
    };

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
    for (int cy = 0; cy < lum.height; ++cy)
    {
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
    for (int x = 0; x < grad.width; ++x) 
    {   
        // Copy the first row as it is, as there's no row above to calculate the energy
        mat_at(energy, 0, x) = mat_at(grad, 0, x);
    }
    
    for (int cy = 1; cy < grad.height; ++cy)
    {
        for (int cx = 0; cx < grad.width; ++cx)
        {       
            // M(i, j) = e(i, j) + min(M(i − 1, j − 1), M(i − 1, j), M(i − 1, j + 1))
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

static Mat gaussianBlur(Mat lum)
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
        {1.0, 4.0 , 6.0 , 4.0 , 1.0},
        {4.0, 16.0, 24.0, 16.0, 4.0},
        {6.0, 24.0, 36.0, 24.0, 6.0},
        {4.0, 16.0, 24.0, 16.0, 4.0},
        {1.0, 4.0 , 6.0 , 4.0 , 1.0}
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

static void remove_seam(Img img, Mat lum, Mat egy)
{
    // Finding the minimum energy seam or selecting a random one
    int y = img.height - 1;  
    int seam = 0;
    for (int x = 0; x < img.width; ++x) 
    {
        if (mat_at(egy, y, x) < mat_at(egy, y, seam)) seam = x;
    }

    remove_pixel_img(img, y, seam);
    // Since lum is one to one mapping of image we can just remove the pixels instead of 
    // recomputing the lum everytime
    remove_pixel_mat(lum, y, seam);

    for (y = img.height - 2; y >= 0; --y)
    {       
        int original_seam = seam;
        for (int dx = -1; dx <= 1; ++dx)
        {
            int x = original_seam + dx;                 
            if (0 <= x && x < img.width && mat_at(egy, y, x) < mat_at(egy, y, seam))
            {   
                seam = x;
            }
        }
        remove_pixel_img(img, y, seam);
        remove_pixel_mat(lum, y, seam);
        // seam_image.pixels[y * seam_image.width + seam] = RED;     // Visualize the seams
    }
}

int main(int argc, char *argv[])
{       
   if (argc != 3)
    {
        fprintf(stderr, "USAGE: <c_exec_file> <image_path (.jpgs file only)> <seams_to_remove>\n");
        exit(EXIT_FAILURE);
    }

    const char *inp_img_fp = argv[1];
    size_t seams_to_remove = atoi(argv[2]);

    int _width, _height;
    uint32_t *pixels_ = (uint32_t *) stbi_load(inp_img_fp, &_width, &_height, NULL, 4);
    if (!pixels_)  
    {
        fprintf(stderr, "ERROR: Failed to load image %s\n", inp_img_fp);
        exit(EXIT_FAILURE);
    }

    Img img = img_alloc(_width, _height);
    img.pixels = pixels_;

    Mat lum = mat_alloc(_width, _height);
    Mat grad = mat_alloc(_width, _height);
    Mat egy = mat_alloc(_width, _height);  
    
    luminance(img, lum);
    // Iteration the rows from top to bottom and find the lowest energy pixel   
    // or just pick a random seam from the vicinity
    // i.e, the three neighbouring images on top of the low energy pixel 

    // What are we doing here?
    // 'abcdefghijklm' -> image row , lets say 'f' is seam 
    // - pixels_row + seam points to 'f'
    // - pixels_row + seam + 1 points to next pixel of 'f'
    // - copying entire array of &(pixels_row + seam + 1) to 1 step back and then reduce width

    
    for (size_t i = 0; i < seams_to_remove; ++i) {
        
        printf("Removing Seam %zu\n", i);
       
        // Update these at every iteration 
        sobel(lum, grad);   
        energy(grad, egy);

        remove_seam(img, lum, egy);

        // Everytime we remove seam, reduce the width 
        img.width  -= 1;
        lum.width  -= 1;
        grad.width -= 1;
        egy.width  -= 1;  
    }   

    dump_img("output.png", img);

    stbi_image_free(pixels_);
    mat_free(&lum);
    mat_free(&egy);
    mat_free(&grad);

    return 0;
}

