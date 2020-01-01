use crate::dynamics::{Cluster, point::Point2};
use crate::units::consts::G_UNIV;
use crate::geometry::vector::*;

const RESISTANCE: f64 = 0.001;

pub fn nav_stokes(point: &Point2) -> Vector4 {
    let speed = point.state.speed.magnitude();
    let acceleration = point.state.speed * (-RESISTANCE / point.mass * speed);
    Vector4::concat(&point.state.speed, &acceleration)
}

pub fn gravity(point: &Point2, cluster: &Cluster) -> Vector4 {
    let len = cluster.len();
    let mut acceleration = Vector2::zeros();
    let mut distance: Vector2;
    let mut magnitude: f64;
    for i in 0..cluster.len() {
        distance = cluster[i].center.state.position - point.state.position;
        magnitude = distance.magnitude();
        if magnitude < std::f64::EPSILON {
            continue;
        }
        acceleration += distance * G_UNIV * cluster[i].center.mass / (magnitude * magnitude * magnitude);
    }
    Vector4::concat(&point.state.speed, &acceleration)
}
