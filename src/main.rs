use std::mem::swap;

use image::{imageops, ImageBuffer, Rgba, RgbaImage};

mod utils;
mod vector;

const WIDTH: u32 = 100;
const HEIGHT: u32 = 100;

// first attempt
// fn line(x0: i32, y0: i32, x1: i32, y1: i32, img: &mut ImageBuffer<Rgba<u8>, Vec<u8>>, c: Rgba<u8>) {
//     for t in 0..100 {
//         let t_f32 = t as f32 * 0.01;
//         let x = x0 as f32 + (x1 - x0) as f32 * t_f32;
//         let y = y0 as f32 + (y1 - y0) as f32 * t_f32;
//         img.put_pixel(x as u32, y as u32, c);
//     }
// }

// second attempt
// fn line(x0: i32, y0: i32, x1: i32, y1: i32, img: &mut ImageBuffer<Rgba<u8>, Vec<u8>>, c: Rgba<u8>) {
//     for x in x0..x1 {
//         let t = (x - x0) as f32 / (x1 - x0) as f32;
//         let y = y0 as f32 * (1.0 - t) + y1 as f32 * t;
//         img.put_pixel(x as u32, y as u32, c);
//     }
// }

// third attempt
// fn line(mut x0: i32, mut y0: i32, mut x1: i32, mut y1: i32, img: &mut ImageBuffer<Rgba<u8>, Vec<u8>>, c: Rgba<u8>) {
//     let mut steep = false;
//     if (x0 - x1).abs() < (y0 - y1).abs() { // if the line is steep, we transpose the image
//         (x0, y0) = (y0, x0);
//         (x1, y1) = (y1, x1);
//         steep = true;
//     }
//     if x0 > x1 { // make it left to right
//         (x0, x1) = (x1, x0);
//         (y0, y1) = (y1, y0);
//     }

//     for x in x0..x1 { // step along x is just the number of pixels to draw
//         let t = (x - x0) as f32 / (x1 - x0) as f32; // t goes from 0 to 1
//         let y = y0 as f32 * (1.0 - t) + y1 as f32 * t;
//         if steep {
//             img.put_pixel(y as u32, x as u32, c);
//         } else {
//             img.put_pixel(x as u32, y as u32, c);
//         }
//     }
// }

// fourth attempt
// fn line(
//     mut x0: i32,
//     mut y0: i32,
//     mut x1: i32,
//     mut y1: i32,
//     img: &mut ImageBuffer<Rgba<u8>, Vec<u8>>,
//     c: Rgba<u8>,
// ) {
//     let mut steep = false;
//     if (x0 - x1).abs() < (y0 - y1).abs() {
//         // if the line is steep, we transpose the image
//         (x0, y0) = (y0, x0);
//         (x1, y1) = (y1, x1);
//         steep = true;
//     }
//     if x0 > x1 {
//         // make it left to right
//         (x0, x1) = (x1, x0);
//         (y0, y1) = (y1, y0);
//     }

//     // dx is the same, so pull out of loop
//     let dx = x1 - x0;
//     let dy = y1 - y0;
//     let d_err = (dy as f32 / dx as f32).abs();
//     // gives the dist to the best line from current (x, y) pixel
//     let mut error = 0f32;
//     let mut y = y0;

//     // step along x is just the number of pixels to draw
//     for x in x0..x1 {
//         if steep {
//             img.put_pixel(y as u32, x as u32, c);
//         } else {
//             img.put_pixel(x as u32, y as u32, c);
//         }
//         error += d_err;
//         // error is > 1 pixel (0.5 is half the pixel)
//         if error > 0.5 {
//             y += if y1 > y0 { 1 } else { -1 };
//             error -= 1f32;
//         }
//     }
// }

// fifth and final attempt
fn line(
    mut x0: i32,
    mut y0: i32,
    mut x1: i32,
    mut y1: i32,
    img: &mut ImageBuffer<Rgba<u8>, Vec<u8>>,
    c: Rgba<u8>,
) {
    let mut steep = false;
    if (x0 - x1).abs() < (y0 - y1).abs() {
        // if the line is steep, we transpose the image
        (x0, y0) = (y0, x0);
        (x1, y1) = (y1, x1);
        steep = true;
    }
    if x0 > x1 {
        // make it left to right
        (x0, x1) = (x1, x0);
        (y0, y1) = (y1, y0);
    }
    
    // dx is the same, so pull out of loop
    let dx = x1 - x0;
    let dy = y1 - y0;
    let d_err2 = 2 * dy.abs();
    // gives the dist to the best line from current (x, y) pixel
    let mut error2 = 0;
    let mut y = y0;
    let y_incr = if y1 > y0 { 1 } else { -1 };
    
    // step along x is just the number of pixels to draw
    if steep {
        for x in x0..x1 {
            img.put_pixel(y as u32, x as u32, c);
            error2 += d_err2;
            // error is > 1 pixel (0.5 is half the pixel)
            if error2 > dx {
                y += y_incr;
                error2 -= dx * 2;
            }
        }
    } else {
        for x in x0..x1 {
            img.put_pixel(x as u32, y as u32, c);
            error2 += d_err2;
            // error is > 1 pixel (0.5 is half the pixel)
            if error2 > dx {
                y += y_incr;
                error2 -= dx * 2;
            }
        }
    }
}

fn main() {
    let mut img = RgbaImage::from_pixel(WIDTH, HEIGHT, Rgba([0, 0, 0, 255]));

    line(13, 20, 80, 40, &mut img, Rgba([255, 255, 255, 255]));
    line(20, 13, 40, 80, &mut img, Rgba([255, 0, 0, 255]));
    line(80, 40, 13, 20, &mut img, Rgba([255, 0, 0, 255]));

    imageops::flip_vertical_in_place(&mut img);
    img.save("test.tga").unwrap();
}
