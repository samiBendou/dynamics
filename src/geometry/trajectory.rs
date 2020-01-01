use std::fmt;
use std::fmt::Debug;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

use crate::geometry::vector::Vector2;

pub const TRAJECTORY_SIZE: usize = 256;

#[derive(Copy, Clone)]
pub struct Trajectory2 {
    positions: [Vector2; TRAJECTORY_SIZE],
    index: usize,
}

impl From<Vector2> for Trajectory2 {
    fn from(position: Vector2) -> Self {
        let position = position.clone();
        Trajectory2::new([position; TRAJECTORY_SIZE], 0)
    }
}

impl From<[Vector2; TRAJECTORY_SIZE]> for Trajectory2 {
    fn from(positions: [Vector2; TRAJECTORY_SIZE]) -> Self {
        Trajectory2::new(positions, 0)
    }
}

impl Trajectory2 {
    pub fn new(positions: [Vector2; TRAJECTORY_SIZE], index: usize) -> Trajectory2 {
        Trajectory2 { positions, index }
    }

    pub fn zeros() -> Trajectory2 {
        let position = Vector2::zeros();
        Trajectory2::new([position; TRAJECTORY_SIZE], 0)
    }

    pub fn reset0(&mut self) -> &mut Self {
        for position in self.positions.iter_mut() {
            position.reset0();
        }
        self
    }

    pub fn reset1(&mut self) -> &mut Self {
        for position in self.positions.iter_mut() {
            position.reset0();
        }
        self
    }

    pub fn reset(&mut self, position: &Vector2) -> &mut Self {
        for pos in self.positions.iter_mut() {
            *pos = *position;
        }
        self
    }


    pub fn push(&mut self, position: &Vector2) {
        self.positions[self.index] = *position;
        self.index = self.index_offset(0);
    }


    pub fn position(&self, i: usize) -> &Vector2 {
        &self.positions[self.index_offset(i)]
    }


    pub fn position_mut(&mut self, i: usize) -> &mut Vector2 {
        &mut self.positions[self.index_offset(i)]
    }


    fn index_offset(&self, i: usize) -> usize {
        (i + self.index + 1) % TRAJECTORY_SIZE
    }
}

impl Debug for Trajectory2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let mut buffer = String::new();
        for i in (TRAJECTORY_SIZE - 32)..TRAJECTORY_SIZE {
            buffer += format!("\n{:?}", self.positions[self.index_offset(i)]).as_str();
        }
        write!(f, "{}", buffer)
    }
}

impl Add<Trajectory2> for Trajectory2 {
    type Output = Trajectory2;

    fn add(self, rhs: Trajectory2) -> Self::Output {
        let mut ret = self;
        ret += rhs;
        ret
    }
}

impl Add<Vector2> for Trajectory2 {
    type Output = Trajectory2;

    fn add(self, rhs: Vector2) -> Self::Output {
        let mut ret = self;
        ret += rhs;
        ret
    }
}

impl AddAssign<Trajectory2> for Trajectory2 {
    fn add_assign(&mut self, rhs: Trajectory2) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[self.index_offset(i)] += rhs.positions[rhs.index_offset(i)];
        }
    }
}

impl AddAssign<Vector2> for Trajectory2 {
    fn add_assign(&mut self, rhs: Vector2) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[i] += rhs;
        }
    }
}

impl Sub<Trajectory2> for Trajectory2 {
    type Output = Trajectory2;

    fn sub(self, rhs: Trajectory2) -> Self::Output {
        let mut ret = self;
        ret -= rhs;
        ret
    }
}

impl Sub<Vector2> for Trajectory2 {
    type Output = Trajectory2;

    fn sub(self, rhs: Vector2) -> Self::Output {
        let mut ret = self;
        ret -= rhs;
        ret
    }
}

impl SubAssign<Trajectory2> for Trajectory2 {
    fn sub_assign(&mut self, rhs: Trajectory2) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[self.index_offset(i)] -= rhs.positions[rhs.index_offset(i)];
        }
    }
}

impl SubAssign<Vector2> for Trajectory2 {
    fn sub_assign(&mut self, rhs: Vector2) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[i] -= rhs;
        }
    }
}

impl Mul<f64> for Trajectory2 {
    type Output = Trajectory2;

    fn mul(self, rhs: f64) -> Self::Output {
        let mut ret = self;
        ret *= rhs;
        ret
    }
}

impl MulAssign<f64> for Trajectory2 {
    fn mul_assign(&mut self, rhs: f64) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[i] *= rhs;
        }
    }
}

impl Div<f64> for Trajectory2 {
    type Output = Trajectory2;

    fn div(self, rhs: f64) -> Self::Output {
        let mut ret = self;
        ret /= rhs;
        ret
    }
}

impl DivAssign<f64> for Trajectory2 {
    fn div_assign(&mut self, rhs: f64) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[i] /= rhs;
        }
    }
}