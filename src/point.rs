use geomath::prelude::*;
use geomath::point;
use geomath::vector::*;

#[derive(Copy, Clone)]
pub struct Point2 {
    pub state: point::Point2,
    pub gradient: Vector4,
    pub mass: f64,
}

#[derive(Copy, Clone)]
pub struct Point3 {
    pub state: point::Point3,
    pub gradient: Vector6,
    pub mass: f64,
}

macro_rules! impl_point {
    ($PointN:ident, $VectorN:ident, $Vector2N:ident) => {
        impl $PointN {
            #[inline]
            pub fn new(state: point::$PointN, mass: f64) -> $PointN {
                $PointN { state, mass, gradient: $Vector2N::zeros() }
            }

            #[inline]
            pub fn inertial(position: $VectorN, speed: $VectorN, mass: f64) -> $PointN {
                $PointN::new(point::$PointN::new(position, speed), mass)
            }

            #[inline]
            pub fn immobile(position: $VectorN, mass: f64) -> $PointN {
                $PointN::new(point::$PointN::from(position), mass)
            }

            #[inline]
            pub fn zeros(mass: f64) -> $PointN {
                $PointN::new(point::$PointN::zeros(), mass)
            }

            #[inline]
            pub fn kinetic_energy(&self) -> f64 {
                let speed = self.state.speed.magnitude();
                0.5 * self.mass * speed * speed
            }

            #[inline]
            pub fn angular_momentum(&self) -> f64 {
                self.state.position.area(&self.state.speed) * self.mass
            }

            #[inline]
            pub fn accelerate(&mut self, dt: f64) -> &mut Self {
                self.state += self.gradient * dt;
                self
            }
        }
    }
}

impl_point!(Point2, Vector2, Vector4);
impl_point!(Point3, Vector3, Vector6);