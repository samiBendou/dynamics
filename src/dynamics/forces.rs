use geomath::common::*;
use geomath::vector::*;

use crate::dynamics::point::Point3;
use crate::units::consts::G_UNIV;

const RESISTANCE: f64 = 0.001;

pub fn nav_stokes(point: &Point3) -> Vector6 {
    let speed = point.state.speed.magnitude();
    let acceleration = point.state.speed * (-RESISTANCE / point.mass * speed);
    Vector6::concat(&point.state.speed, &acceleration)
}

pub fn gravity(point: &Point3, points: &Vec<Point3>) -> Vector6 {
    let len = points.len();
    let mut acceleration = Vector3::zeros();
    let mut distance: Vector3;
    let mut magnitude: f64;
    for i in 0..len {
        distance = points[i].state.position - point.state.position;
        magnitude = distance.magnitude();
        if magnitude < std::f64::EPSILON {
            continue;
        }
        acceleration += distance * G_UNIV * points[i].mass / (magnitude * magnitude * magnitude);
    }
    Vector6::concat(&point.state.speed, &acceleration)
}
