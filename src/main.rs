mod renderer;
mod utils;
mod model;

use std::path::Path;
use nalgebra::Vector2;
use image::{imageops, Rgba, RgbaImage};
use renderer::{draw_wireframe, triangle};

use crate::model::Model;

const WIDTH: u32 = 200;
const HEIGHT: u32 = 200;

fn main() {
    let mut img = RgbaImage::from_pixel(WIDTH, HEIGHT, Rgba([0, 0, 0, 255]));
    
    // let obj = Model::from(Path::new("obj/body.obj")).expect("Error parsing .obj file!");
    
    // draw_wireframe(obj, &mut img, Rgba([255, 255, 255, 255]), WIDTH, HEIGHT);

    let mut v0 = Vector2::new(32f64, 90f64);
    let mut v1 = Vector2::new(60f64, 100f64);
    let mut v2 = Vector2::new(60f64, 140f64);
    
    let mut u0 = Vector2::new(40f64, 60f64); 
    let mut u1 = Vector2::new(91f64, 75f64);
    let mut u2 = Vector2::new(22f64, 0f64);

    let mut w0 = Vector2::new(180f64, 150f64);
    let mut w1 = Vector2::new(120f64, 160f64);
    let mut w2 = Vector2::new(130f64, 180f64);

    triangle(&mut v0, &mut v1, &mut v2, &mut img, Rgba([255, 0, 0, 255]));
    triangle(&mut u0, &mut u1, &mut u2, &mut img, Rgba([255, 255, 255, 255]));
    triangle(&mut w0, &mut w1, &mut w2, &mut img, Rgba([0, 255, 0, 255]));

    imageops::flip_vertical_in_place(&mut img);
    img.save("new.tga").unwrap();
}
