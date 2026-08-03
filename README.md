# Seam Carving

Seam carving is a content-aware image resizing technique that removes or inserts low-energy seams while preserving the most important visual regions of an image.

## Demo

<div align="center">
 <div style="display: flex; justify-content: center; gap: 20px;">
   <div>
     <img src="./demo-images/tower.jpg" height="200" alt="Original image">
     <p><b>Original image</b></p>
   </div>
   <div>
     <img src="./demo-images/tower_cropped.png" height="200" alt="Seam carved image">
     <p><b>Seam carved image</b></p>
   </div>
 </div>
</div>

## Features

- Content-aware image resizing
- Optimized seam-energy computation in C
- Clean source layout with dedicated source and include directories

## Project layout

- [src/main.c](src/main.c): seam carving implementation and CLI entry point
- [include/utils.h](include/utils.h): shared helper macros
- [nob.c](nob.c): lightweight build driver for the project
- [external/](external/): bundled third-party headers for image I/O

## Build and run

1. Clone the repository:
  ```sh
  git clone https://github.com/SamarthPyati/seam-carving.git
  ```
2. Change into the project directory:
  ```sh
  cd seam-carving
  ```
3. Build the executable:
  ```sh
  cc nob.c -o nob
  ```
4. Run the program:
  ```sh
  ./nob <image_file_path> <seams_to_remove> [output_file_path]
  ```

Example:

```sh
./nob demo-images/tower.jpg 5
./nob demo-images/tower.jpg 5 output.png
```

If no output path is supplied, the image is written to [output.png](output.png).

## Acknowledgements

- Seam carving paper: [Seam Carving for Content-Aware Image Resizing](https://dl.acm.org/doi/10.1145/1275808.1276390)

## Contributing

Feel free to fork the repository and submit pull requests.

