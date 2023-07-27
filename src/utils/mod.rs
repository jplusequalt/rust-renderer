use std::{fs::File, io::BufReader, path::Path, env};

use anyhow::{Context, Result};
use image::{RgbaImage, imageops};
use obj::{load_obj, Obj, TexturedVertex};

pub fn parse_model(
    assets_root: &Path,
    obj_name: &str,
    diffuse_name: &str,
) -> Result<(Obj<TexturedVertex>, RgbaImage)> {
    let reader = BufReader::new(File::open(assets_root.join(obj_name))?);
    let model: Obj<TexturedVertex> = load_obj(reader)?;
    let mut texture = image::open(&Path::new(&assets_root.join(diffuse_name)))?.into_rgba8();
    
    imageops::flip_vertical_in_place(&mut texture);

    Ok((model, texture))
}
