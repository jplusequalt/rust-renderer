use std::ops;

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Vec2 {
  pub x: f64,
  pub y: f64
}

impl Vec2 {
  pub fn new() -> Self {
      Vec2::from(0.0, 0.0)
  }

  pub fn from(x: f64, y: f64) -> Self {
      Vec2 { x, y }
  }

  #[inline]
  pub fn length(&self) -> f64 {
      f64::sqrt(Self::length_squared(&self))
  }

  #[inline]
  pub fn length_squared(&self) -> f64 {
      self.x * self.x + self.y * self.y
  }

  #[inline]
  pub fn dot(u: &Vec2, v: &Vec2) -> f64 {
      u.x * v.x + u.y * v.y 
  }

  pub fn reflect(v: &Vec2, n: &Vec2) -> Self {
      *v - *n * Self::dot(v, n) * 2f64
  }

  pub fn refract(uv: &Vec2, n: &Vec2, etai_over_etat: f64) -> Self {
      let cos_theta = f64::min(Self::dot(&-*uv, &n), 1.0);
      let r_out_perp = (*uv + (*n * cos_theta)) * etai_over_etat;
      let r_out_parallel = *n * -f64::sqrt(f64::abs(1.0 - r_out_perp.length_squared()));
      r_out_perp + r_out_parallel
  }

  pub fn near_zero(&self) -> bool {
      let s = 1e-8;
      (f64::abs(self.x) < s) && (f64::abs(self.y) < s)
  }

  #[inline]
  pub fn unit_vector(v: &Vec2) -> Self {
      *v / v.length()
  }
}

impl ops::Neg for Vec2 {
  type Output = Self;

  fn neg(self) -> Self::Output {
      Vec2 {
          x: -self.x,
          y: -self.y
      }
  }
}

impl ops::AddAssign for Vec2 {
  fn add_assign(&mut self, rhs: Self) {
      *self = Self {
          x: self.x + rhs.x,
          y: self.y + rhs.y
      }
  }
}

impl ops::MulAssign for Vec2 {
  fn mul_assign(&mut self, rhs: Self) {
      *self = Self {
          x: self.x * rhs.x,
          y: self.y * rhs.y
      }
  }
}

impl ops::DivAssign for Vec2 {
  fn div_assign(&mut self, rhs: Self) {
      *self = Self {
          x: self.x / rhs.x,
          y: self.y / rhs.y
      }
  }
}

impl ops::Add for Vec2 {
  type Output = Vec2;

  fn add(self, rhs: Self) -> Self::Output {
      Vec2 {
          x: self.x + rhs.x,
          y: self.y + rhs.y
      }
  }
}

impl ops::Sub for Vec2 {
  type Output = Vec2;

  fn sub(self, rhs: Self) -> Self::Output {
      Vec2 {
          x: self.x - rhs.x,
          y: self.y - rhs.y
      }
  }
}

impl ops::Mul<Vec2> for Vec2 {
  type Output = Vec2;
  fn mul(self, rhs: Vec2) -> Self::Output {
      Vec2 {
          x: self.x * rhs.x,
          y: self.y * rhs.y
      }
  }
}

impl ops::Mul<f64> for Vec2 {
  type Output = Vec2;
  fn mul(self, rhs: f64) -> Self::Output {
      Vec2 {
          x: self.x * rhs,
          y: self.y * rhs
      }
  }
}

impl ops::Div<f64> for Vec2 {
  type Output = Vec2;
  fn div(self, rhs: f64) -> Self::Output {
      self * (1.0 / rhs)
  }
}