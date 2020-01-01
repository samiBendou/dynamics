use crate::dynamics::{Cluster, point::Point2};
use crate::units::consts::G_UNIV;

pub fn gravity(point: &Point2, cluster: &Cluster) -> f64 {
    let len = cluster.len();
    let mut ret = 0.;
    let mut distance: f64;
    for i in 0..len {
        distance = cluster[i].center.state.position % point.state.position;
        if distance < std::f64::EPSILON {
            continue;
        }
        ret += G_UNIV * cluster[i].center.mass / distance;
    }
    ret
}