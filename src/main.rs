mod renderer;
mod utils;

use anyhow::Result;
use image::{imageops, Rgba, RgbaImage};
use nalgebra::Vector3;
use std::{path::Path, borrow::BorrowMut};

use renderer::render;
use utils::parse_model;

const WIDTH: u32 = 800;
const HEIGHT: u32 = 800;

fn main() {
    match parse_model(
        Path::new("assets/african_head"),
        "african_head.obj",
        "african_head_diffuse.tga",
    ) {
        Ok((model, texture)) => {
            let mut img = RgbaImage::from_pixel(WIDTH, HEIGHT, Rgba([0, 0, 0, 255]));

            let light_dir = Vector3::new(1., 1.0, 1.0).normalize();
            let camera = Vector3::new(0.0, 0.0, 1.0);
            let model_position = Vector3::new(0.0, 0.0, 0.0);
            let up = Vector3::new(0.0, 1.0, 0.0);

            render(&model, &mut img, texture, light_dir, camera, model_position, up);

            imageops::flip_vertical_in_place(&mut img);
            img.save("result.png").unwrap();
        }
        Err(e) => eprintln!("{e}")
    }
}
