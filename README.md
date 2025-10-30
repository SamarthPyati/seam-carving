# Seam Carving

Seam Carving is an implementation of content-aware image resizing algorithm that selectively removes or inserts pixels in an image to change its dimensions while preserving important visual content without any application of extensive machine learning models. 

## Demo
<p align="center">
  <figure style="display: inline-block; margin: 0 10px; text-align: center; padding:20px">
    <img src="./demo-images/tower.jpg" height="200">
    <figcaption><b>Original Image</b></figcaption>
  </figure>
  <figure style="display: inline-block; margin: 0 10px; text-align: center; padding:20px">
    <img src="./demo-images/tower_cropped.png" height="200">
    <figcaption><b>Seam Carved Image</b></figcaption>
  </figure>
</p>

## Features

- Content-aware image resizing
- Efficient implementation in C
- Parellel Processing with OpenMP

## How to Run

1. Clone the repository:
   ```sh
   git clone https://github.com/SamarthPyati/seam-carving.git
   ```
2. Navigate to the project directory:
   ```sh
   cd seam-carving
   ```
3. Compile the program:
   ```sh
   ./nob
   ```
4. Run the program:
   ```sh
   ./build/main <image_file_path> <seam_to_remove>
   ```

## Acknowledgements
 Seam Carving Paper: [Seam carving for content-aware image resizing](https://dl.acm.org/doi/10.1145/1275808.1276390)

## Contributing

Feel free to fork this repository and submit pull requests. 
