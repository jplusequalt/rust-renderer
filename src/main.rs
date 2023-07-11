mod renderer;
mod utils;

use std::{fs::File, io::BufReader, path::Path};

use image::{imageops, Rgba, RgbaImage};
use nalgebra::Vector3;
use obj::{load_obj, Obj, TexturedVertex};
use renderer::{render};

const WIDTH: u32 = 800;
const HEIGHT: u32 = 800;

fn main() {
    let mut img = RgbaImage::from_pixel(WIDTH, HEIGHT, Rgba([0, 0, 0, 255]));

    let reader =
        BufReader::new(File::open("obj/african_head/african_head.obj").expect("Error reading obj file!"));
    let obj: Obj<TexturedVertex> = load_obj(reader).expect("Error parsing obj file!");
    let mut texture = image::open(Path::new("obj/african_head/african_head_diffuse.tga"))
        .expect("Error opening diffuse file!")
        .into_rgba8();

    imageops::flip_vertical_in_place(&mut texture);

    render(&obj, &mut img, texture);

    imageops::flip_vertical_in_place(&mut img);
    img.save("result.png").unwrap();
}
