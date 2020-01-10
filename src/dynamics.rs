use std::fmt;
use std::fmt::Debug;
use std::ops::{Index, IndexMut};

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

pub struct Cluster {
    pub points: Vec<Point3>,
    barycenter: Point3,
}

impl Cluster {
    pub fn new(bodies: Vec<Point3>) -> Self {
        Cluster {
            points: bodies,
            barycenter: Point3::zeros(0.),
        }
    }

    pub fn empty() -> Self {
        Cluster::new(vec![])
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.points.len() == 0
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.points.len()
    }

    #[inline]
    pub fn barycenter(&self) -> &Point3 {
        &self.barycenter
    }

    #[inline]
    pub fn kinetic_energy(&self) -> f64 {
        self.points.iter().map(|point| point.kinetic_energy()).sum()
    }

    #[inline]
    pub fn angular_momentum(&self) -> f64 {
        self.points.iter().map(|point| point.angular_momentum()).sum()
    }

    pub fn potential_energy<T>(&self, mut f: T) -> f64 where
        T: FnMut(&Vec<Point3>, usize) -> f64 {
        let len = self.points.len();
        let mut ret = 0.;
        for i in 0..len {
            ret += f(&self.points, i);
        }
        ret * 0.5
    }

    /*
    pub fn remove_aways(&mut self) -> &mut Self {
        let (max_distance, max_index) = self.max_distance();
        let (mean, deviation, _distances) = if self.points.len() < 3 {
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
    */

    pub fn push(&mut self, point: Point3) -> &mut Self {
        self.points.push(point);
        self.update_barycenter()
    }

    pub fn pop(&mut self) -> Option<Point3> {
        let ret = self.points.pop();
        self.update_barycenter();
        ret
    }

    pub fn remove(&mut self, i: usize) -> Point3 {
        let ret = self.points.remove(i);
        self.update_barycenter();
        ret
    }

    #[inline]
    pub fn reset0_at(&mut self, i: usize) -> &mut Self {
        self.points[i].state.reset0();
        self.points[i].state.trajectory.reset0();
        self.update_barycenter()
    }

    #[inline]
    pub fn reset_trajectory_at(&mut self, i: usize) -> &mut Self {
        let position = self.points[i].state.position;
        self.points[i].state.trajectory.reset(&position);
        self.update_barycenter()
    }

    #[inline]
    pub fn translate_at(&mut self, i: usize, direction: &Vector3) -> &mut Self {
        self.points[i].state.position += *direction;
        self.points[i].state.update_trajectory();
        self.update_barycenter()
    }

    #[inline]
    pub fn translate(&mut self, direction: &Vector3) -> &mut Self {
        self.barycenter.state.position += *direction;
        for body in self.points.iter_mut() {
            body.state.position += *direction;
        }
        self.update_trajectory()
    }

    #[inline]
    pub fn apply<T>(&mut self, solver: &mut Solver, f: T) -> &mut Self where
        T: FnMut(&Vec<Point3>, usize) -> Vector6 {
        solver.step(&mut self.points, f);
        self.update_barycenter()
            .update_trajectory()
    }

    #[inline]
    pub fn set_absolute(&mut self, origin: &geometry::point::Point3) -> &mut Self {
        self.barycenter.state += *origin;
        for body in self.points.iter_mut() {
            body.state += *origin;
        }
        self
    }

    #[inline]
    pub fn set_relative(&mut self, origin: &geometry::point::Point3) -> &mut Self {
        self.barycenter.state -= *origin;
        *self.barycenter.state.trajectory.last_mut() -= *origin.trajectory.last();
        for body in self.points.iter_mut() {
            body.state -= *origin;
            *body.state.trajectory.last_mut() -= *origin.trajectory.last();
        }
        self
    }

    #[inline]
    pub fn reset_origin(&mut self, origin: &geometry::point::Point3, old_origin: &geometry::point::Point3) -> &mut Self {
        self.barycenter.state.reset_origin(origin, old_origin);
        for body in self.points.iter_mut() {
            body.state.reset_origin(origin, old_origin);
        }
        self
    }

    #[inline]
    pub fn update_barycenter(&mut self) -> &mut Self {
        self.barycenter.mass = 0.;
        self.barycenter.state.reset0();
        for body in self.points.iter() {
            self.barycenter.mass += body.mass;
            self.barycenter.state += body.state * body.mass;
        }
        self.barycenter.state /= self.barycenter.mass;
        self
    }

    #[inline]
    pub fn reset_trajectory(&mut self) -> &mut Self {
        self.barycenter.state.trajectory.reset(&self.barycenter.state.position);
        for body in self.points.iter_mut() {
            body.state.trajectory.reset(&body.state.position);
        }
        self
    }

    #[inline]
    fn update_trajectory(&mut self) -> &mut Self {
        self.barycenter.state.update_trajectory();
        for body in self.points.iter_mut() {
            body.state.update_trajectory();
        }
        self
    }
}

impl Debug for Cluster {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let mut buffer = String::from("");
        buffer.push_str(format!("{:?}\n\n", self.barycenter.state).as_str());
        for body in self.points.iter() {
            buffer.push_str(format!("{:?}\n\n", body.state).as_str());
        }
        write!(f, "{}", buffer)
    }
}

impl Index<usize> for Cluster {
    type Output = Point3;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        &self.points[index]
    }
}

impl IndexMut<usize> for Cluster {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.points[index]
    }
}