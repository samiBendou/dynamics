use crate::consts::G_UNIV;
use crate::point::Point3;

pub fn gravity(point: &Point3, points: &Vec<Point3>) -> f64 {
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