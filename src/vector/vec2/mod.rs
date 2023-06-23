pub struct Vec2 {
  pub x: f64,
  pub y: f64
}

impl Vec2 {
  pub fn new() -> Self {
    Vec2::from(0f64, 0f64)
  }

  pub fn from(x: f64, y: f64) -> Self {
    Vec2 { x, y }
  }
}