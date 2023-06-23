pub struct Vec4 {
  pub x: f64,
  pub y: f64,
  pub z: f64,
  pub w: f64
}

impl Vec4 {
  pub fn from(x: f64, y: f64, z: f64, w: f64) -> Self {
    Vec4 { x, y, z, w }
  }
}