mod renderer;
mod utils;
mod vector;
mod model;

use std::path::Path;

use image::{imageops, Rgba, RgbaImage};
use renderer::draw_wireframe;

use crate::model::Model;

const WIDTH: u32 = 1000;
const HEIGHT: u32 = 1000;

fn main() {
    let mut img = RgbaImage::from_pixel(WIDTH, HEIGHT, Rgba([0, 0, 0, 255]));
    
    let obj = Model::from(Path::new("obj/body.obj")).expect("Error parsing .obj file!");
    
    draw_wireframe(obj, &mut img, Rgba([255, 255, 255, 255]), WIDTH, HEIGHT);

    imageops::flip_vertical_in_place(&mut img);
    img.save("body.tga").unwrap();
}
