//!
//! A point encapsulates a kinematic point from [geomath](https://crates.io/crate/geomath) and
//! its physical properties such as the mass.
//!
//! * The kinematic state is stored in the field `state`, it contains speed, position and trajectory.
//! * The time derivative of this state is stored in the field `gradient`.
//!
//!
//! **Notes :**
//! * for now only the mass is implemented
//! * For a point having N degrees of freedom, The gradient will be of size 2N.
//!
//! ## Purpose
//! This abstraction aims to represent what's call a point in dynamics in order to simulate point's
//! dynamic.
//!
//! It also provides useful features such as momentum and energy computation.
//!
//! ## Conventions
//! * `gradient` is only a buffer for ODE solving, avoid modify it directly unless you need
//! another solver than the provided one
//! * Use SI units, notably, mass must be in kilograms

use geomath::prelude::*;
use geomath::point;
use geomath::vector::*;

#[derive(Copy, Clone)]
/// 2D point
pub struct Point2 {
    pub state: point::Point2,
    pub gradient: Vector4,
    pub mass: f64,
}

/// 3D point
#[derive(Copy, Clone)]
pub struct Point3 {
    pub state: point::Point3,
    pub gradient: Vector6,
    pub mass: f64,
}

macro_rules! impl_point {
    ($PointN:ident, $VectorN:ident, $Vector2N:ident) => {
        impl $PointN {
            /// Construct a point with given kinematic point and mass
            #[inline]
            pub fn new(state: point::$PointN, mass: f64) -> $PointN {
                $PointN { state, mass, gradient: $Vector2N::zeros() }
            }

            /// Construct a point with given position, speed and mass
            #[inline]
            pub fn inertial(position: $VectorN, speed: $VectorN, mass: f64) -> $PointN {
                $PointN::new(point::$PointN::new(position, speed), mass)
            }

            /// Construct a point with given position and mass and no speed
            #[inline]
            pub fn immobile(position: $VectorN, mass: f64) -> $PointN {
                $PointN::new(point::$PointN::from(position), mass)
            }

            /// Construct a point located at zero, with given mass and no speed
            #[inline]
            pub fn zeros(mass: f64) -> $PointN {
                $PointN::new(point::$PointN::zeros(), mass)
            }

            /// Get the kinetic energy of the point
            #[inline]
            pub fn kinetic_energy(&self) -> f64 {
                let speed = self.state.speed.magnitude();
                0.5 * self.mass * speed * speed
            }

            /// Get the norm of angular momentum
            #[inline]
            pub fn angular_momentum(&self) -> f64 {
                self.state.position.area(&self.state.speed) * self.mass
            }

            /// Accelerate the point according to its current `gradient` field
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