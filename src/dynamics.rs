use std::fmt;
use std::fmt::Debug;
use std::ops::{Index, IndexMut};

use rand::Rng;

use crate::dynamics::orbital::Orbit;
use crate::dynamics::point::Point2;
use crate::dynamics::solver::{Method, Solver};
use crate::geometry;
use crate::geometry::point::ZERO;
use crate::geometry::vector::{Vector2, Vector4};

pub mod point;
pub mod forces;
pub mod potentials;
pub mod orbital;
pub mod solver;

pub const SPEED_SCALING_FACTOR: f64 = 5e-7;

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
    pub center: Point2,
}

impl Body {
    pub fn new(name: &str, center: Point2) -> Body {
        Body { name: String::from(name), center }
    }

    pub fn orbital(name: &str, orbit: &Orbit, true_anomaly: f64, mass: f64) -> Body {
        let position = orbit.position_at(true_anomaly);
        let speed = orbit.speed_at(true_anomaly);
        Body::new(name, Point2::inertial(position, speed, mass))
    }
}

impl Debug for Body {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        write!(f, "name: {}\n{:?}", self.name, self.center.state)
    }
}

pub struct Cluster {
    pub bodies: Vec<Body>,
    barycenter: Point2,
    origin: geometry::point::Point2,
    current: usize,
    frame: Frame,
    solver: Solver,
}

impl Cluster {
    pub fn new(bodies: Vec<Body>) -> Self {
        Cluster {
            bodies,
            barycenter: Point2::zeros(0.),
            origin: ZERO,
            current: 0,
            solver: Solver::new(1., Method::RungeKutta4),
            frame: Frame::Zero,
        }
    }

    pub fn empty() -> Self {
        Cluster::new(vec![])
    }

    pub fn orbital(cluster: &orbital::Cluster, true_anomalies: Vec<f64>) -> Self {
        let len = cluster.bodies.len();
        let mut bodies: Vec<Body> = Vec::with_capacity(len);
        let mut body;
        for i in 0..len {
            body = Body::orbital(
                &cluster.bodies[i].name,
                &cluster.bodies[i].orbit,
                true_anomalies[i],
                cluster.bodies[i].mass,
            );
            bodies.push(body);
        }
        Cluster::new(bodies)
    }

    pub fn orbital_at(cluster: &orbital::Cluster, true_anomaly: f64) -> Self {
        let mut true_anomalies = Vec::with_capacity(cluster.bodies.len());
        for _ in cluster.bodies.iter() {
            true_anomalies.push(true_anomaly)
        }
        Cluster::orbital(cluster, true_anomalies)
    }

    pub fn orbital_at_random(cluster: &orbital::Cluster) -> Self {
        let two_pi = 2. * std::f64::consts::PI;
        let mut rng = rand::thread_rng();
        let mut true_anomalies: Vec<f64> = Vec::with_capacity(cluster.bodies.len());
        for _ in cluster.bodies.iter() {
            true_anomalies.push(rng.gen_range(0., two_pi))
        }
        Cluster::orbital(cluster, true_anomalies)
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
    pub fn barycenter(&self) -> &Point2 {
        &self.barycenter
    }

    #[inline]
    pub fn current(&self) -> Option<&Body> {
        self.bodies.get(self.current)
    }

    #[inline]
    pub fn current_mut(&mut self) -> Option<&mut Body> {
        self.bodies.get_mut(self.current)
    }

    #[inline]
    pub fn current_index(&self) -> usize {
        self.current
    }

    #[inline]
    pub fn last(&self) -> Option<&Body> { self.bodies.last() }

    #[inline]
    pub fn last_mut(&mut self) -> Option<&mut Body> { self.bodies.last_mut() }

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
            distance = self.bodies[i].center.state.distance(&self.barycenter.state.position);
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
            distances.push(self.bodies[i].center.state.distance(&self.barycenter.state.position));
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

    #[inline]
    fn update_origin(&mut self) -> &mut Self {
        self.origin.position = self.current_origin().position;
        self.origin.speed = self.current_origin().speed;
        self
    }

    fn current_origin(&mut self) -> &geometry::point::Point2 {
        if self.is_empty() {
            return &ZERO;
        }
        match self.frame {
            Frame::Zero => &ZERO,
            Frame::Current => &self.bodies[self.current].center.state,
            Frame::Barycenter => &self.barycenter.state,
        }
    }

    pub fn reset_origin(&mut self, origin: geometry::point::Point2) -> &mut Self {
        if self.is_empty() {
            return self;
        }
        self.barycenter.state.update_origin(&origin, &ZERO);
        for body in self.bodies.iter_mut() {
            body.center.state.update_origin(&origin, &ZERO);
        }
        self.origin = origin;
        self
    }

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

    pub fn push(&mut self, body: Body) -> &mut Self {
        self.bodies.push(body);
        self.update_barycenter()
    }

    pub fn pop(&mut self) -> Option<Body> {
        let len = self.bodies.len();
        if self.current != 0 && self.current == len - 1 {
            self.current -= 1;
        }
        let body = self.bodies.pop();
        self.update_barycenter();
        body
    }

    pub fn remove(&mut self, index: usize) -> Body {
        let len = self.bodies.len();
        if index == len - 1 {
            self.pop().unwrap()
        } else {
            if self.current == len - 1 {
                self.current -= 1;
            }
            let ret = self.bodies.remove(index);
            self.update_barycenter();
            ret
        }
    }

    pub fn update_current_index(&mut self, increase: bool, bypass_last: bool) -> &mut Self {
        if increase {
            self.increase_current(bypass_last);
        } else {
            self.decrease_current();
        }
        if self.frame == Frame::Current {
            let origin = *self.current_origin();
            self.reset_origin(origin).update_barycenter();
        }
        self
    }

    #[inline]
    pub fn update_frame(&mut self) -> &mut Self {
        self.frame.next();
        let origin = *self.current_origin();
        self.reset_origin(origin).update_barycenter()
    }

    #[inline]
    pub fn reset0_current(&mut self) -> &mut Self {
        self.bodies[self.current].center.state.reset0();
        self.bodies[self.current].center.state.trajectory.reset0();
        self.update_barycenter()
    }

    #[inline]
    pub fn reset_current_trajectory(&mut self) -> &mut Self {
        let position = self.bodies[self.current].center.state.position;
        self.bodies[self.current].center.state.trajectory.reset(&position);
        self.update_barycenter()
    }

    #[inline]
    pub fn update_current_trajectory(&mut self) -> &mut Self {
        self.bodies[self.current].center.state.update_trajectory();
        self
    }

    #[inline]
    pub fn translate_current(&mut self, direction: &Vector2) -> &mut Self {
        self.bodies[self.current].center.state.position += *direction;
        self.update_barycenter()
    }

    pub fn translate(&mut self, direction: &Vector2) -> &mut Self {
        for body in self.bodies.iter_mut() {
            body.center.state.position += *direction;
        }
        self.update_barycenter()
    }

    #[inline]
    pub fn apply<T>(&mut self, dt: f64, iterations: u32, f: T) -> &mut Self where
        T: FnMut(&Vec<Body>, usize) -> Vector4 {
        self.set_absolute();
        self.solver.dt = dt;
        self.solver.step(&mut self.bodies, f, iterations);
        self.update_barycenter()
            .update_origin()
            .set_relative()
            .update_trajectory()
    }

    fn set_absolute(&mut self) -> &mut Self {
        self.barycenter.state += self.origin;
        for body in self.bodies.iter_mut() {
            body.center.state += self.origin;
        }
        self
    }

    fn set_relative(&mut self) -> &mut Self {
        self.barycenter.state -= self.origin;
        for body in self.bodies.iter_mut() {
            body.center.state -= self.origin;
        }
        self
    }

    fn update_trajectory(&mut self) -> &mut Self {
        self.barycenter.state.update_trajectory();
        self.origin.update_trajectory();
        for body in self.bodies.iter_mut() {
            body.center.state.update_trajectory();
        }
        self
    }

    pub fn reset_trajectory(&mut self) -> &mut Self {
        let origin_position = self.current_origin().position;
        self.barycenter.state.trajectory.reset(&self.barycenter.state.position);
        self.origin.trajectory.reset(&origin_position);
        for body in self.bodies.iter_mut() {
            body.center.state.trajectory.reset(&body.center.state.position);
        }
        self
    }

    fn decrease_current(&mut self) -> &mut Self {
        if self.current > 0 {
            self.current -= 1;
        }
        self
    }

    fn increase_current(&mut self, bypass_last: bool) -> &mut Self {
        let offset = if bypass_last { 2 } else { 1 };
        if self.current < self.bodies.len() - offset {
            self.current += 1;
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