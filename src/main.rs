use std::io::BufWriter;

use image::{Rgba, RgbaImage, imageops, ImageOutputFormat};

mod vector;
mod utils;

const WIDTH: u32 = 100;
const HEIGHT: u32 = 100;

fn main() {
    let mut img = RgbaImage::from_pixel(WIDTH, HEIGHT, Rgba([103, 60, 200, 255]));
    img.put_pixel(52, 41, Rgba([255, 0, 0, 255]));
    img = imageops::flip_vertical(&img);
    img.save("test.tga").unwrap();
}
