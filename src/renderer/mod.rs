use std::mem::swap;

use image::{Rgba, RgbaImage};
use nalgebra::{Vector2, Vector3};

use crate::model::Model;

// fifth and final attempt
pub fn line(v0: &mut Vector2<f64>, v1: &mut Vector2<f64>, img: &mut RgbaImage, c: Rgba<u8>) {
    let (mut x0, mut y0) = (v0.x as i32, v0.y as i32);
    let (mut x1, mut y1) = (v1.x as i32, v1.y as i32);
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
            let v1 = obj.vert(face[(j + 1) % 3] as usize).unwrap();

            let mut u0 = Vector2::new(
                ((v0.x + 1f64) * (width - 1) as f64) / 2.0,
                ((v0.y + 1f64) * (height - 1) as f64) / 2.0,
            );
            let mut u1 = Vector2::new(
                ((v1.x + 1f64) * (width - 1) as f64) / 2.0,
                ((v1.y + 1f64) * (height - 1) as f64) / 2.0,
            );

            line(&mut u0, &mut u1, img, c);
        }
    }
}

pub fn triangle(
    v0: &mut Vector2<f64>,
    v1: &mut Vector2<f64>,
    v2: &mut Vector2<f64>,
    img: &mut RgbaImage,
    c: Rgba<u8>,
) {
    if v0.y == v1.y && v0.y == v2.y {
        return ();
    }

    // // sort the verts in ascending order (v0 <= v1 <= v2)
    if v0.y > v1.y {
        swap(v0, v1);
    }
    if v0.y > v2.y {
        swap(v0, v2);
    }
    if v1.y > v2.y {
        swap(v1, v2);
    }

    // height of triangle
    let total_h = v2.y - v0.y;

    // boundary a will be from v0 -> v2
    // boundary b will be split in two parts: v0 -> v1, v1 -> v2
    for i in 0..(total_h as i32) {
        let upper_seg = (i > (v1.y - v0.y) as i32) || (v1.y == v0.y);

        let segment_h = if upper_seg {
            v2.y - v1.y
        } else {
            v1.y - v0.y
        };

        let alpha = i as f64 / total_h as f64;
        let beta = (i as f64 - (if upper_seg { v1.y - v0.y } else { 0.0 })) / segment_h;

        let mut a = *v0 + (*v2 - *v0) * alpha;
        let mut b = if upper_seg {
            *v1 + (*v2 - *v1) * beta
        } else {
            *v0 + (*v1 - *v0) * beta
        };

        if a.x > b.x  {
            swap(&mut a, &mut b);
        }

        for j in (a.x as u32)..=(b.x as u32) {
            img.put_pixel(j as u32, (v0.y + i as f64) as u32, c);
        }
    }
}
