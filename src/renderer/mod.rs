use image::{Rgba, RgbaImage};

use crate::model::Model;

// fifth and final attempt
pub fn line(
    mut x0: i32,
    mut y0: i32,
    mut x1: i32,
    mut y1: i32,
    img: &mut RgbaImage,
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
        for x in x0..=x1 {
            img.put_pixel(y as u32, x as u32, c);
            error2 += d_err2;
            // error is > dx
            if error2 > dx {
                y += y_incr;
                error2 -= dx * 2;
            }
        }
    } else {
        for x in x0..=x1 {
            img.put_pixel(x as u32, y as u32, c);
            error2 += d_err2;
            // error is > dx
            if error2 > dx {
                y += y_incr;
                error2 -= dx * 2;
            }
        }
    }
}

pub fn draw_wireframe(obj: Model, img: &mut RgbaImage, c: Rgba<u8>, width: u32, height: u32) {
    for i in 0..obj.num_faces() {
        let face = obj.face(i).unwrap();
        for j in 0..3 {
            let v0 = obj.vert(face[j] as usize).unwrap();
            let v1 = obj.vert(face[(j+1) % 3] as usize).unwrap();
            let x0 = ((v0.x + 1f64) * (width - 1) as f64) / 2.0;
            let y0 = ((v0.y + 1f64) * (height - 1) as f64) / 2.0;
            let x1 = ((v1.x + 1f64) * (width - 1) as f64) / 2.0;
            let y1 = ((v1.y + 1f64) * (height - 1) as f64) / 2.0;
            line(x0 as i32, y0 as i32, x1 as i32, y1 as i32, img, c);
        }
    }
}
