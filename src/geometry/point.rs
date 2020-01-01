use std::fmt;
use std::fmt::Debug;
use std::ops::{AddAssign, DivAssign, Mul, MulAssign, SubAssign, Add, Sub, Div};

use crate::geometry::vector::{Angle, Split, Vector2, Vector4};
use crate::geometry::trajectory::Trajectory2;

pub trait State<T> {
    fn set_state(&mut self, state: &T) -> &mut Self;
}

#[derive(Copy, Clone)]
pub struct Point2 {
    pub position: Vector2,
    pub speed: Vector2,
    pub trajectory: Trajectory2,
}

impl From<Vector2> for Point2 {
    fn from(vector: Vector2) -> Self {
        Point2::new(vector, Vector2::zeros())
    }
}

impl From<Vector4> for Point2 {
    fn from(vector: Vector4) -> Self {
        Point2::new(vector.upper(), vector.lower())
    }
}


impl Point2 {
    pub fn new(position: Vector2, speed: Vector2) -> Point2 {
        Point2 {
            position,
            speed,
            trajectory: Trajectory2::from(position.clone()),
        }
    }

    #[inline]
    pub fn zeros() -> Point2 {
        Point2::new(Vector2::zeros(), Vector2::zeros())
    }

    pub fn reset0(&mut self) -> &mut Self {
        self.position.reset0();
        self.speed.reset0();
        self.trajectory.reset0();
        self
    }

    pub fn reset1(&mut self) -> &mut Self {
        self.position.reset1();
        self.speed.reset1();
        self.trajectory.reset1();
        self
    }

    pub fn reset(&mut self, position: &Vector2, speed: &Vector2) -> &mut Self {
        self.position = *position;
        self.speed = *speed;
        self.trajectory.reset(position);
        self
    }

    pub fn distance(&self, position: &Vector2) -> f64 {
        self.position.distance(*position)
    }

    pub fn set_origin(&mut self, origin: &Point2, old_origin: &Option<Point2>) -> &mut Self {
        if let Some(old_origin) = old_origin {
            *self += *old_origin;
            self.trajectory += old_origin.trajectory;
        }
        *self -= *origin;
        self.trajectory -= origin.trajectory;
        self
    }
}

impl Debug for Point2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        write!(
            f,
            "position: {:?}\nspeed: {:?}\n",
            self.position,
            self.speed,
        )
    }
}

impl PartialEq for Point2 {
    fn eq(&self, other: &Self) -> bool {
        self.position == other.position && self.speed == other.speed
    }

    fn ne(&self, other: &Self) -> bool {
        self.position != other.position || self.speed != other.speed
    }
}

impl Add<Vector2> for Point2 {
    type Output = Point2;

    fn add(self, rhs: Vector2) -> Self::Output {
        let mut ret = self;
        ret += rhs;
        ret
    }
}

impl Add<Point2> for Point2 {
    type Output = Point2;

    fn add(self, rhs: Point2) -> Self::Output {
        let mut ret = self;
        ret += rhs;
        ret
    }
}

impl Add<Vector4> for Point2 {
    type Output = Point2;

    fn add(self, rhs: Vector4) -> Self::Output {
        let mut ret = self;
        ret += rhs;
        ret
    }
}

impl AddAssign<Vector2> for Point2 {
    fn add_assign(&mut self, rhs: Vector2) {
        self.position += rhs;
    }
}

impl AddAssign<Point2> for Point2 {
    fn add_assign(&mut self, rhs: Point2) {
        self.position += rhs.position;
        self.speed += rhs.speed;
    }
}

impl AddAssign<Vector4> for Point2 {
    fn add_assign(&mut self, rhs: Vector4) {
        self.position.x += rhs.x;
        self.position.y += rhs.y;
        self.speed.x += rhs.z;
        self.speed.y += rhs.w;
    }
}

impl Sub<Vector2> for Point2 {
    type Output = Point2;

    fn sub(self, rhs: Vector2) -> Self::Output {
        let mut ret = self;
        ret -= rhs;
        ret
    }
}

impl Sub<Point2> for Point2 {
    type Output = Point2;

    fn sub(self, rhs: Point2) -> Self::Output {
        let mut ret = self;
        ret -= rhs;
        ret
    }
}

impl Sub<Vector4> for Point2 {
    type Output = Point2;

    fn sub(self, rhs: Vector4) -> Self::Output {
        let mut ret = self;
        ret -= rhs;
        ret
    }
}

impl SubAssign<Vector2> for Point2 {
    fn sub_assign(&mut self, rhs: Vector2) {
        self.position -= rhs;
    }
}

impl SubAssign<Point2> for Point2 {
    fn sub_assign(&mut self, rhs: Point2) {
        self.position -= rhs.position;
        self.speed -= rhs.speed;
    }
}

impl SubAssign<Vector4> for Point2 {
    fn sub_assign(&mut self, rhs: Vector4) {
        self.position.x -= rhs.x;
        self.position.y -= rhs.y;
        self.speed.x -= rhs.z;
        self.speed.y -= rhs.w;
    }
}

impl Mul<f64> for Point2 {
    type Output = Point2;

    fn mul(self, rhs: f64) -> Self::Output {
        let mut output = self;
        output *= rhs;
        output
    }
}

impl MulAssign<f64> for Point2 {
    fn mul_assign(&mut self, rhs: f64) {
        self.position *= rhs;
        self.speed *= rhs;
    }
}

impl Div<f64> for Point2 {
    type Output = Point2;

    fn div(self, rhs: f64) -> Self::Output {
        let mut output = self;
        output /= rhs;
        output
    }
}

impl DivAssign<f64> for Point2 {
    fn div_assign(&mut self, rhs: f64) {
        self.position /= rhs;
        self.speed /= rhs;
    }
}

