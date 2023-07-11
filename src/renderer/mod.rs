#![allow(dead_code)]
use std::mem::swap;

use image::{Pixel, Rgba, RgbaImage};
use nalgebra::{Vector2, Vector3};
use obj::{Obj, TexturedVertex};

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

pub fn barycentric(points: &[Vector2<f32>; 3], p: &Vector3<f32>) -> Vector3<f32> {
    // create two vectors from the sides of the tri
    let u = Vector3::new(
        points[2].x - points[0].x,
        points[1].x - points[0].x,
        points[0].x - p.x,
    );
    let v = Vector3::new(
        points[2].y - points[0].y,
        points[1].y - points[0].y,
        points[0].y - p.y,
    );

    // compute the normal
    let n = u.cross(&v);

    // barycentric coords must all be > 0
    // if abs(n.z) is close to zero, then the triangle is degenerate
    // return a vec3 with negative coords
    if f32::abs(n.z) < 0.001 {
        return Vector3::new(-1.0, 1.0, 1.0);
    }

    Vector3::new(1.0 - (n.x + n.y) / n.z, n.y / n.z, n.x / n.z)
}

pub fn face_world_to_screen(
    obj: &Obj<TexturedVertex>,
    face: &[u16],
    width: f32,
    height: f32,
) -> [Vector2<i32>; 3] {
    // v0, v1, v2 are the vertices of the triangle defined by this face in world space
    let v0 = obj.vertices[face[0] as usize].position;
    let v1 = obj.vertices[face[1] as usize].position;
    let v2 = obj.vertices[face[2] as usize].position;

    let p0 = Vector2::new(
        ((v0[0] + 1.0) * width) as i32,
        ((v0[1] + 1.0) * height) as i32,
    );

    let p1 = Vector2::new(
        ((v1[0] + 1.0) * width) as i32,
        ((v1[1] + 1.0) * height) as i32,
    );

    let p2 = Vector2::new(
        ((v2[0] + 1.0) * width) as i32,
        ((v2[1] + 1.0) * height) as i32,
    );

    [p0, p1, p2]
}

pub fn face_world_coords(obj: &Obj<TexturedVertex>, face: &[u16]) -> [Vector3<f32>; 3] {
    // v0, v1, v2 are the vertices of the triangle defined by this face in world space
    let v0 = obj.vertices[face[0] as usize].position;
    let v1 = obj.vertices[face[1] as usize].position;
    let v2 = obj.vertices[face[2] as usize].position;

    let p0 = Vector3::new(v0[0], v0[1], v0[2]);
    let p1 = Vector3::new(v1[0], v1[1], v1[2]);
    let p2 = Vector3::new(v2[0], v2[1], v2[2]);

    [p0, p1, p2]
}

fn face_uv_coords(
    obj: &Obj<TexturedVertex>,
    face: &[u16],
    texture_width: f32,
    texture_height: f32,
) -> [Vector2<f32>; 3] {
    // v0, v1, v2 are the texture coordinates at the vertices of the triangle
    let v0 = obj.vertices[face[0] as usize].texture;
    let v1 = obj.vertices[face[1] as usize].texture;
    let v2 = obj.vertices[face[2] as usize].texture;

    let t0 = Vector2::new(
        (v0[0] * texture_width) - 1.0,
        (v0[1] * texture_height) - 1.0,
    );
    let t1 = Vector2::new(
        (v1[0] * texture_width) - 1.0,
        (v1[1] * texture_height) - 1.0,
    );
    let t2 = Vector2::new(
        (v2[0] * texture_width) - 1.0,
        (v2[1] * texture_height) - 1.0,
    );

    [t0, t1, t2]
}

fn calc_light_intesity(world_coords: &[Vector3<f32>; 3], light_dir: &Vector3<f32>) -> f32 {
    // calc normal from two edges of the triangle
    let v0 = Vector3::from(world_coords[2] - world_coords[0]);
    let v1 = Vector3::from(world_coords[1] - world_coords[0]);
    let norm = v0.cross(&v1).normalize();

    // light intensity is the dot product of the normal and light vector
    norm.dot(&light_dir)
}

#[inline]
fn edge_function(v0: &Vector2<i32>, v1: &Vector2<i32>, v2: &Vector2<i32>) -> i32 {
    (v2.x - v0.x) * (v1.y - v0.y) - (v2.y - v0.y) * (v1.x - v0.x)
}

pub fn render(obj: &Obj<TexturedVertex>, img: &mut RgbaImage, texture: RgbaImage) {
    // depth buffer
    // indicates for each pixel, how far away that pixel is from the camera
    let mut z_buffer = vec![f32::NEG_INFINITY; (img.width() * img.height()) as usize];

    let faces = &obj.indices[..obj.indices.len()];
    for face in faces.chunks(3) {
        // fetch data for this tri
        let screen_coords = face_world_to_screen(
            obj,
            face,
            (img.width() - 1) as f32 / 2.0,
            (img.height() - 1) as f32 / 2.0,
        );
        let world_coords = face_world_coords(obj, face);
        let uv_coords = face_uv_coords(obj, face, texture.width() as f32, texture.height() as f32);

        let light_dir = Vector3::new(0.0, 0.0, -1.0);
        let intensity = calc_light_intesity(&world_coords, &light_dir);

        if intensity > 0.0 {
            draw_triangle_barycentric(
                &screen_coords,
                &world_coords,
                &uv_coords,
                &texture,
                intensity,
                &mut z_buffer,
                img,
            );
        }
    }
}

fn draw_triangle_barycentric(
    screen_coords: &[Vector2<i32>; 3],
    world_coords: &[Vector3<f32>; 3],
    uv_coords: &[Vector2<f32>; 3],
    texture: &RgbaImage,
    light_intensity: f32,
    z_buffer: &mut Vec<f32>,
    img: &mut RgbaImage,
) {
    let mut bbox_min = Vector2::new(f32::MAX, f32::MAX);
    let mut bbox_max = Vector2::new(f32::MIN, f32::MIN);
    let clamped_bbox = Vector2::new(img.width() as f32 - 1.0, img.height() as f32 - 1.0);

    for i in 0..3 {
        bbox_min.x = f32::max(0.0, f32::min(bbox_min.x, screen_coords[i].x as f32));
        bbox_min.y = f32::max(0.0, f32::min(bbox_min.y, screen_coords[i].y as f32));

        bbox_max.x = f32::min(
            clamped_bbox.x,
            f32::max(bbox_max.x, screen_coords[i].x as f32),
        );
        bbox_max.y = f32::min(
            clamped_bbox.y,
            f32::max(bbox_max.y, screen_coords[i].y as f32),
        );
    }

    let area = edge_function(&screen_coords[0], &screen_coords[1], &screen_coords[2]) as f32;

    // Calculate if point2 of the bounding box is inside triangle
    for x in (bbox_min.x as i32)..=(bbox_max.x as i32) {
        for y in (bbox_min.y as i32)..=(bbox_max.y as i32) {
            let p = Vector2::new(x, y);

            // check if pixel is in triangle
            let mut w0 = edge_function(&screen_coords[1], &screen_coords[2], &p) as f32;
            let mut w1 = edge_function(&screen_coords[2], &screen_coords[0], &p) as f32;
            // Barycentric coordinates
            w0 /= area;
            w1 /= area;
            let w2 = 1.0 - w0 - w1;

            if w0 >= 0.0 && w1 >= 0.0 && w2 >= 0.0 {
                // interpolate z
                let z_value =
                    w0 * world_coords[0].z + w1 * world_coords[1].z + w2 * world_coords[2].z;

                // check if distance to intersection is closer than value
                // currently stored in the depth buffer
                if z_buffer[(x + y * (img.width() as i32 - 1)) as usize] <= z_value {
                    z_buffer[(x + y * (img.width() as i32 - 1)) as usize] = z_value;

                    // interpolate texture coords
                    let u = w0 * uv_coords[0].x + w1 * uv_coords[1].x + w2 * uv_coords[2].x;
                    let v = w0 * uv_coords[0].y + w1 * uv_coords[1].y + w2 * uv_coords[2].y;

                    // apply shading
                    let mut uv = texture.get_pixel(u as u32, v as u32).to_rgba();
                    uv.apply_without_alpha(|ch| ((ch as f32) * light_intensity) as u8);

                    img.put_pixel(x as u32, y as u32, uv);
                }
            }
        }
    }
}

pub fn triangle_from_verts(
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

        let segment_h = if upper_seg { v2.y - v1.y } else { v1.y - v0.y };

        let alpha = i as f64 / total_h as f64;
        let beta = (i as f64 - (if upper_seg { v1.y - v0.y } else { 0.0 })) / segment_h;

        let mut a = *v0 + (*v2 - *v0) * alpha;
        let mut b = if upper_seg {
            *v1 + (*v2 - *v1) * beta
        } else {
            *v0 + (*v1 - *v0) * beta
        };

        if a.x > b.x {
            swap(&mut a, &mut b);
        }

        for j in (a.x as u32)..=(b.x as u32) {
            img.put_pixel(j, (v0.y + i as f64) as u32, c);
        }
    }
}
