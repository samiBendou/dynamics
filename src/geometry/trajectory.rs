use std::fmt;
use std::fmt::Debug;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

use crate::geometry::common::{Initializer, Reset};
use crate::geometry::vector::{Vector2, Vector3};

pub type Trajectory2 = Trajectory<Vector2>;
pub type Trajectory3 = Trajectory<Vector3>;

pub const TRAJECTORY_SIZE: usize = 256;
pub const ZERO: Trajectory3 = Trajectory3 {
    positions: [Vector3 { x: 0., y: 0., z: 0. }; TRAJECTORY_SIZE],
    index: 0,
};

#[derive(Copy, Clone)]
pub struct Trajectory<T> {
    positions: [T; TRAJECTORY_SIZE],
    index: usize,
}

impl<T> From<T> for Trajectory<T> where
    T: Copy + Clone
{
    fn from(position: T) -> Self {
        Trajectory::new([position; TRAJECTORY_SIZE], 0)
    }
}

impl<T> From<[T; TRAJECTORY_SIZE]> for Trajectory<T> where
    T: Copy + Clone {
    fn from(positions: [T; TRAJECTORY_SIZE]) -> Self {
        Trajectory::new(positions, 0)
    }
}

impl<T> Trajectory<T> where
    T: Copy + Clone {
    #[inline]
    pub fn new(positions: [T; TRAJECTORY_SIZE], index: usize) -> Trajectory<T> {
        Trajectory { positions, index }
    }

    #[inline]
    pub fn push(&mut self, position: &T) {
        self.positions[self.index] = *position;
        self.index = self.index_offset(0);
    }

    #[inline]
    pub fn position(&self, i: usize) -> &T {
        &self.positions[self.index_offset(i)]
    }

    #[inline]
    pub fn position_mut(&mut self, i: usize) -> &mut T {
        &mut self.positions[self.index_offset(i)]
    }

    #[inline]
    fn index_offset(&self, i: usize) -> usize {
        (i + self.index + 1) % TRAJECTORY_SIZE
    }
}

impl<T> Initializer for Trajectory<T> where
    T: Initializer + Copy + Clone {
    #[inline]
    fn zeros() -> Self {
        Trajectory::new([T::zeros(); TRAJECTORY_SIZE], 0)
    }

    #[inline]
    fn ones() -> Self {
        Trajectory::new([T::ones(); TRAJECTORY_SIZE], 0)
    }
}

impl<T> Reset<T> for Trajectory<T> where
    T: Reset<T> + Copy + Clone {
    #[inline]
    fn reset0(&mut self) -> &mut Self {
        for position in self.positions.iter_mut() {
            position.reset0();
        }
        self
    }

    #[inline]
    fn reset1(&mut self) -> &mut Self {
        for position in self.positions.iter_mut() {
            position.reset1();
        }
        self
    }

    #[inline]
    fn reset(&mut self, val: &T) -> &mut Self {
        for pos in self.positions.iter_mut() {
            pos.reset(val);
        }
        self
    }
}

impl<T> Debug for Trajectory<T> where
    T: Debug + Copy + Clone {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let mut buffer = String::new();
        for i in (TRAJECTORY_SIZE - 32)..TRAJECTORY_SIZE {
            buffer += format!("\n{:?}", self.positions[self.index_offset(i)]).as_str();
        }
        write!(f, "{}", buffer)
    }
}

impl<T> Add<Trajectory<T>> for Trajectory<T> where
    T: AddAssign<T> + Copy + Clone {
    type Output = Trajectory<T>;

    #[inline]
    fn add(self, rhs: Trajectory<T>) -> Self::Output {
        let mut ret = self;
        ret += rhs;
        ret
    }
}

impl<T> Add<T> for Trajectory<T> where
    T: AddAssign<T> + Copy + Clone {
    type Output = Trajectory<T>;

    #[inline]
    fn add(self, rhs: T) -> Self::Output {
        let mut ret = self;
        ret += rhs;
        ret
    }
}

impl<T> AddAssign<Trajectory<T>> for Trajectory<T> where
    T: AddAssign<T> + Copy + Clone {
    #[inline]
    fn add_assign(&mut self, rhs: Trajectory<T>) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[self.index_offset(i)] += rhs.positions[rhs.index_offset(i)];
        }
    }
}

impl<T> AddAssign<T> for Trajectory<T> where
    T: AddAssign<T> + Copy + Clone {
    #[inline]
    fn add_assign(&mut self, rhs: T) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[i] += rhs;
        }
    }
}

impl<T> Sub<Trajectory<T>> for Trajectory<T> where
    T: SubAssign<T> + Copy + Clone {
    type Output = Trajectory<T>;

    #[inline]
    fn sub(self, rhs: Trajectory<T>) -> Self::Output {
        let mut ret = self;
        ret -= rhs;
        ret
    }
}

impl<T> Sub<T> for Trajectory<T> where
    T: SubAssign<T> + Copy + Clone {
    type Output = Trajectory<T>;

    #[inline]
    fn sub(self, rhs: T) -> Self::Output {
        let mut ret = self;
        ret -= rhs;
        ret
    }
}

impl<T> SubAssign<Trajectory<T>> for Trajectory<T> where
    T: SubAssign<T> + Copy + Clone {
    #[inline]
    fn sub_assign(&mut self, rhs: Trajectory<T>) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[self.index_offset(i)] -= rhs.positions[rhs.index_offset(i)];
        }
    }
}

impl<T> SubAssign<T> for Trajectory<T> where
    T: SubAssign<T> + Copy + Clone {
    #[inline]
    fn sub_assign(&mut self, rhs: T) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[i] -= rhs;
        }
    }
}

impl<T> Mul<f64> for Trajectory<T> where
    T: MulAssign<f64> + Copy + Clone {
    type Output = Trajectory<T>;

    #[inline]
    fn mul(self, rhs: f64) -> Self::Output {
        let mut ret = self;
        ret *= rhs;
        ret
    }
}

impl<T> MulAssign<f64> for Trajectory<T> where
    T: MulAssign<f64> + Copy + Clone {
    #[inline]
    fn mul_assign(&mut self, rhs: f64) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[i] *= rhs;
        }
    }
}

impl<T> Div<f64> for Trajectory<T> where
    T: DivAssign<f64> + Copy + Clone {
    type Output = Trajectory<T>;

    #[inline]
    fn div(self, rhs: f64) -> Self::Output {
        let mut ret = self;
        ret /= rhs;
        ret
    }
}

impl<T> DivAssign<f64> for Trajectory<T> where
    T: DivAssign<f64> + Copy + Clone {
    #[inline]
    fn div_assign(&mut self, rhs: f64) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[i] /= rhs;
        }
    }
}