use crate::dynamics::{Cluster, point::Point2};
use crate::units::consts::G_UNIV;
use crate::vector::*;

const RESISTANCE: f64 = 0.001;

pub fn nav_stokes(point: &Point2) -> Vector2 {
    point.speed * (-RESISTANCE / point.mass * point.speed.magnitude())
}

pub fn gravity(point: &Point2, cluster: &Cluster) -> Vector4 {
    let mut result = Vector2::zeros();
    let mut distance: Vector2;
    let mut magnitude: f64;
    for i in 0..cluster.count() {
        distance = cluster[i].shape.center.position - point.position;
        magnitude = distance.magnitude();
        if magnitude < std::f64::EPSILON {
            continue;
        }
        result += distance * G_UNIV * cluster[i].mass / (magnitude * magnitude * magnitude);
    }
    Vector4::concat(&point.speed, &result)
}
