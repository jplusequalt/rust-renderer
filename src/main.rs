mod renderer;
mod utils;
mod model;

use std::path::Path;
use image::{imageops, Rgba, RgbaImage};
use renderer::render_model;

use crate::model::Model;

const WIDTH: u32 = 800;
const HEIGHT: u32 = 800;

fn main() {
    let mut img = RgbaImage::from_pixel(WIDTH, HEIGHT, Rgba([0, 0, 0, 255]));
    
    let obj = Model::from(Path::new("obj/african_head.obj")).expect("Error parsing .obj file!");

    render_model(obj, &mut img, WIDTH, HEIGHT);

    imageops::flip_vertical_in_place(&mut img);
    img.save("head.tga").unwrap();
}
