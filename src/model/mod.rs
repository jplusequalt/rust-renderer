use std::{
    fs::File,
    io,
    io::{BufRead, BufReader},
    path::Path,
};

use crate::vector::vec3::Vec3;

pub struct Model {
    verts: Vec<Vec3>,
    faces: Vec<Vec<i32>>,
}

impl Model {
    pub fn from(filename: &Path) -> io::Result<Self> {
        let file = File::open(filename)?;
        let content = BufReader::new(file);

        let mut verts: Vec<Vec3> = Vec::new();
        let mut faces: Vec<Vec<i32>> = Vec::new();
        content.lines().for_each(|line| {
            if let Ok(l) = line {
                if l.starts_with("v ") {
                    let split_verts = l
                        .split(" ")
                        .filter_map(|s| s.parse::<f64>().ok())
                        .collect::<Vec<f64>>();

                    verts.push(Vec3::from(split_verts[0], split_verts[1], split_verts[2]));
                } else if l.starts_with("f ") {
                    let split_faces = l
                        .split(" ")
                        .filter_map(|s| {
                            if let Ok(mut index) = s.split("/").nth(0).unwrap().parse::<i32>() {
                                index -= 1;
                                return Some(index);
                            }
                            None
                        })
                        .collect::<Vec<i32>>();
                    faces.push(split_faces);
                }
            }
        });

        Ok(Model { verts, faces })
    }

    pub fn num_verts(&self) -> usize {
        self.verts.len()
    }

    pub fn num_faces(&self) -> usize {
        self.faces.len()
    }

    pub fn face(&self, index: usize) -> Option<&Vec<i32>> {
        self.faces.get(index)
    }

    pub fn vert(&self, index: usize) -> Option<&Vec3> {
        self.verts.get(index)
    }
}
