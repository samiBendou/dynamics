use crate::dynamics::{Body, point::Point3};
use crate::geometry::common::*;
use crate::geometry::common::Split;
use crate::geometry::vector::*;
use crate::units::consts::G_UNIV;

const RESISTANCE: f64 = 0.001;

pub fn nav_stokes(point: &Point3) -> Vector6 {
    let speed = point.state.speed.magnitude();
    let acceleration = point.state.speed * (-RESISTANCE / point.mass * speed);
    Vector6::concat(&point.state.speed, &acceleration)
}

pub fn gravity(point: &Point3, bodies: &Vec<Body>) -> Vector6 {
    let len = bodies.len();
    let mut acceleration = Vector3::zeros();
    let mut distance: Vector3;
    let mut magnitude: f64;
    for i in 0..len {
        distance = bodies[i].center.state.position - point.state.position;
        magnitude = distance.magnitude();
        if magnitude < std::f64::EPSILON {
            continue;
        }
        acceleration += distance * G_UNIV * bodies[i].center.mass / (magnitude * magnitude * magnitude);
    }
    Vector6::concat(&point.state.speed, &acceleration)
}
