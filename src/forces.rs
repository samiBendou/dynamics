//!
//! This module provides functions to compute vectors for the commonly used force in point's dynamics
//!
//! ## Conventions
//! * The forces takes dynamics points as parameters as first parameters
//! * Forces returned are of size 2N for a point with N degrees of freedom

use geomath::prelude::*;
use geomath::vector::{*, self};

use crate::consts::G_UNIV;
use crate::point::Point3;

/// Simplified quadratic fluid drag.
///
/// This is a force that is quadratically opposed to the speed.
/// You change fluid parameters such as density or change cross sectional area by computing coherent values of
/// the resistance.
///
/// **Note :** The resistance is given by the following equation **R = 0.5 D S C** where
/// * **D** is the density of the fluid
/// * **S** is the cross sectional area
/// * **C** is the drag coefficient
///
pub fn quadratic_drag(point: &Point3, resistance: f64) -> Vector6 {
    let speed = point.state.speed.magnitude();
    let acceleration = point.state.speed * (-resistance / point.mass * speed);
    Vector6::concat(&point.state.speed, &acceleration)
}

/// Newton's n-body gravity attraction
///
/// This the sum of the 2-body newtonian forces felt on a specific point surrounded by masses at specific locations.
pub fn newton_gravity(point: &Point3, points: &Vec<Point3>) -> Vector6 {
    let len = points.len();
    let mut acceleration = vector::consts::ZEROS_3;
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
