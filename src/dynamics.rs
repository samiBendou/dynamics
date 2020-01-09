use std::fmt;
use std::fmt::Debug;
use std::ops::{Index, IndexMut};

use crate::dynamics::orbital::Orbit;
use crate::dynamics::point::Point3;
use crate::dynamics::solver::Solver;
use crate::geometry;
use crate::geometry::common::*;
use crate::geometry::vector::{Vector3, Vector6};

pub mod point;
pub mod forces;
pub mod potentials;
pub mod orbital;
pub mod solver;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Frame {
    Zero,
    Current,
    Barycenter,
}

impl Frame {
    pub fn next(&mut self) {
        use Frame::*;
        *self = match self {
            Zero => Current,
            Current => Barycenter,
            Barycenter => Zero,
        }
    }
}

pub struct Body {
    pub name: String,
    pub center: Point3,
}

impl Body {
    pub fn new(name: &str, center: Point3) -> Body {
        Body { name: String::from(name), center }
    }

    pub fn orbital(name: &str, orbit: &Orbit, true_anomaly: f64, mass: f64) -> Body {
        let position = orbit.position_at(true_anomaly);
        let speed = orbit.speed_at(true_anomaly);
        Body::new(name, Point3::inertial(position, speed, mass))
    }
}

impl Debug for Body {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        write!(f, "name: {}\n{:?}", self.name, self.center.state)
    }
}

pub struct Cluster {
    pub bodies: Vec<Body>,
    barycenter: Point3,
}

impl Cluster {
    pub fn new(bodies: Vec<Body>) -> Self {
        Cluster {
            bodies,
            barycenter: Point3::zeros(0.),
        }
    }

    pub fn empty() -> Self {
        Cluster::new(vec![])
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.bodies.len() == 0
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.bodies.len()
    }

    #[inline]
    pub fn barycenter(&self) -> &Point3 {
        &self.barycenter
    }

    #[inline]
    pub fn kinetic_energy(&self) -> f64 {
        self.bodies.iter().map(|body| body.center.kinetic_energy()).sum()
    }

    #[inline]
    pub fn angular_momentum(&self) -> f64 {
        self.bodies.iter().map(|body| body.center.angular_momentum()).sum()
    }

    pub fn potential_energy<T>(&self, mut f: T) -> f64 where
        T: FnMut(&Cluster, usize) -> f64 {
        let len = self.bodies.len();
        let mut ret = 0.;
        for i in 0..len {
            ret += f(self, i);
        }
        ret * 0.5
    }

    pub fn max_distance(&self) -> (f64, usize) {
        let mut max_distance = 0.;
        let mut max_index: usize = 0;
        let mut distance: f64;
        let len = self.bodies.len();
        for i in 0..len {
            distance = self.bodies[i].center.state.distance(&self.barycenter.state);
            if distance > max_distance {
                max_distance = distance;
                max_index = i;
            }
        }
        (max_distance, max_index)
    }

    pub fn stats_distance_without(&self, index: Option<usize>) -> (f64, f64, Vec<f64>) {
        let len = self.bodies.len();
        let mut mean = 0.;
        let mut sum2 = 0.;
        let mut distances: Vec<f64> = Vec::with_capacity(len);
        let index = match index {
            None => len,
            Some(index) => index,
        };
        for i in 0..len {
            distances.push(self.bodies[i].center.state.distance(&self.barycenter.state));
            if i == index {
                continue;
            }
            mean += distances[i];
            sum2 += distances[i] * distances[i];
        }
        let len = len as f64;
        mean /= len;
        (mean, (sum2 / len - mean * mean).sqrt(), distances)
    }

    pub fn remove_aways(&mut self) -> &mut Self {
        let (max_distance, max_index) = self.max_distance();
        let (mean, deviation, _distances) = if self.bodies.len() < 3 {
            self.stats_distance_without(None)
        } else {
            self.stats_distance_without(Some(max_index))
        };
        if max_distance > mean + 10e2 * deviation {
            self.remove(max_index);
            self.update_barycenter();
        }
        self
    }

    pub fn push(&mut self, body: Body) -> &mut Self {
        self.bodies.push(body);
        self.update_barycenter()
    }

    pub fn pop(&mut self) -> Option<Body> {
        let ret = self.bodies.pop();
        self.update_barycenter();
        ret
    }

    pub fn remove(&mut self, i: usize) -> Body {
        let ret = self.bodies.remove(i);
        self.update_barycenter();
        ret
    }

    #[inline]
    pub fn reset0_at(&mut self, i: usize) -> &mut Self {
        self.bodies[i].center.state.reset0();
        self.bodies[i].center.state.trajectory.reset0();
        self.update_barycenter()
    }

    #[inline]
    pub fn reset_trajectory_at(&mut self, i: usize) -> &mut Self {
        let position = self.bodies[i].center.state.position;
        self.bodies[i].center.state.trajectory.reset(&position);
        self.update_barycenter()
    }

    #[inline]
    pub fn translate_at(&mut self, i: usize, direction: &Vector3) -> &mut Self {
        self.bodies[i].center.state.position += *direction;
        self.bodies[i].center.state.update_trajectory();
        self.update_barycenter()
    }

    #[inline]
    pub fn translate(&mut self, direction: &Vector3) -> &mut Self {
        self.barycenter.state.position += *direction;
        for body in self.bodies.iter_mut() {
            body.center.state.position += *direction;
        }
        self.update_trajectory()
    }

    #[inline]
    pub fn apply<T>(&mut self, solver: &mut Solver, f: T) -> &mut Self where
        T: FnMut(&Vec<Body>, usize) -> Vector6 {
        solver.step(&mut self.bodies, f);
        self.update_barycenter()
            .update_trajectory()
    }

    #[inline]
    pub fn set_absolute(&mut self, origin: &geometry::point::Point3) -> &mut Self {
        self.barycenter.state += *origin;
        for body in self.bodies.iter_mut() {
            body.center.state += *origin;
        }
        self
    }

    #[inline]
    pub fn set_relative(&mut self, origin: &geometry::point::Point3) -> &mut Self {
        self.barycenter.state -= *origin;
        *self.barycenter.state.trajectory.last_mut() -= *origin.trajectory.last();
        for body in self.bodies.iter_mut() {
            body.center.state -= *origin;
            *body.center.state.trajectory.last_mut() -= *origin.trajectory.last();
        }
        self
    }

    #[inline]
    pub fn reset_origin(&mut self, origin: &geometry::point::Point3, old_origin: &geometry::point::Point3) -> &mut Self {
        self.barycenter.state.reset_origin(origin, old_origin);
        for body in self.bodies.iter_mut() {
            body.center.state.reset_origin(origin, old_origin);
        }
        self
    }

    #[inline]
    pub fn update_barycenter(&mut self) -> &mut Self {
        self.barycenter.mass = 0.;
        self.barycenter.state.reset0();
        for body in self.bodies.iter() {
            self.barycenter.mass += body.center.mass;
            self.barycenter.state += body.center.state * body.center.mass;
        }
        self.barycenter.state /= self.barycenter.mass;
        self
    }

    #[inline]
    pub fn reset_trajectory(&mut self) -> &mut Self {
        self.barycenter.state.trajectory.reset(&self.barycenter.state.position);
        for body in self.bodies.iter_mut() {
            body.center.state.trajectory.reset(&body.center.state.position);
        }
        self
    }

    #[inline]
    fn update_trajectory(&mut self) -> &mut Self {
        self.barycenter.state.update_trajectory();
        for body in self.bodies.iter_mut() {
            body.center.state.update_trajectory();
        }
        self
    }
}

impl Debug for Cluster {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let mut buffer = String::from("");
        buffer.push_str(format!("{:?}\n\n", self.barycenter.state).as_str());
        for body in self.bodies.iter() {
            buffer.push_str(format!("{:?}\n\n", body.center.state).as_str());
        }
        write!(f, "{}", buffer)
    }
}

impl Index<usize> for Cluster {
    type Output = Body;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        &self.bodies[index]
    }
}

impl IndexMut<usize> for Cluster {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.bodies[index]
    }
}