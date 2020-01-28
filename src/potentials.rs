//!
//! This module provides functions to compute values for the commonly used potentials in point's dynamics
//!
//! ## Conventions
//! * The potentials takes dynamics points as parameters as first parameters

use crate::consts::G_UNIV;
use crate::point::Point3;

/// Newton's n-body gravitational potential
///
/// This the sum of the 2-body newtonian potential felt on a specific point surrounded by masses at specific locations.
pub fn newton_gravity(point: &Point3, points: &Vec<Point3>) -> f64 {
    let mut ret = 0.;
    let mut distance: f64;
    for i in 0..points.len() {
        distance = points[i].state.position % point.state.position;
        if distance < std::f64::EPSILON {
            continue;
        }
        ret += G_UNIV * points[i].mass / distance;
    }
    ret
}