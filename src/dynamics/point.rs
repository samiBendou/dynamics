use crate::geometry::common::Angle;
use crate::geometry::point;
use crate::geometry::vector::{Vector2, Vector4};

#[derive(Copy, Clone)]
pub struct Point2 {
    pub state: point::Point2,
    pub gradient: Vector4,
    pub mass: f64,
}

impl Point2 {

    pub fn new(state: point::Point2, mass: f64) -> Point2 {
        Point2 { state, mass, gradient: Vector4::zeros() }
    }


    pub fn inertial(position: Vector2, speed: Vector2, mass: f64) -> Point2 {
        Point2::new(point::Point2::new(position, speed), mass)
    }


    pub fn immobile(position: Vector2, mass: f64) -> Point2 {
        Point2::new(point::Point2::from(position), mass)
    }


    pub fn zeros(mass: f64) -> Point2 {
        Point2::new(point::Point2::zeros(), mass)
    }

    pub fn kinetic_energy(&self) -> f64 {
        let speed = self.state.speed.magnitude();
        0.5 * self.mass * speed * speed
    }

    pub fn angular_momentum(&self) -> f64 {
        self.state.position.area(&self.state.speed) * self.mass
    }

    pub fn accelerate(&mut self, dt: f64) -> &mut Self {
        self.state += self.gradient * dt;
        self
    }
}