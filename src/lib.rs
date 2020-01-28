//!
//! dynamics is a physics framework that provides tools for accurate dynamics simulation:
//!
//! * n-body solving
//! * forces and potentials
//! * energy and momentum
//! * 2-body orbitals
//! * common constants
//!
//! It relies on a straight forward API that mainly offers a simple solver for points dynamics with
//! an easy to use forces generator.
//! Furthermore many tools are provided to serialize, deserialize in JSON and interpolate orbits.
//!
//! It still uncomplete and other features and improvements are planned,
//! take a look at the [GitHub developement branch](https://github.com/samiBendou/dynamics)  for more details.
//!
//! **Note :** It uses [geomath](https://crates.io/crates/geomath) for vectors, matrices and points.
//!
//! # Get started
//! You can find a detailled documentation to get started with specifical features in each module of this crate.
//!
//! ## Conventions
//! Some simple conventions have been choosed either to gain clarity or performance.
//!
//! * Bodies are considered punctual and called *dynamic points*
//! * n-body problems are defined in the form **a = g** for each body
//! * **a** is the acceleration of the body
//! * **g** is the sum of the forces felt by the body divided by the mass
//!
//! ## Simulation flow
//! The simulation flow is very simple :
//! * Instanciate a cluster containing the dynamic points
//! * Instanciate a solver with choosen parameters
//! * Call the method [`apply`](struct.Cluster.html#method.apply) passing solver and closure as parameters
//! * The closure must return the value of the **g** vector to apply
//!
//! Each time you call the method [`apply`](struct.Cluster.html#method.apply),
//! the solver updates the kinematical state of the cluster according to its current state and the value of **g**.
//!
//!
//! ### Example
//! ```rust
//! use geomath::prelude::*;
//! use geomath::vector;
//! use dynamics::{solver::*, point::*, Cluster};
//! use vector::Vector6;
//!
//! // Set initial position and speed
//! let position = vector::consts::ZEROS_3;
//! let speed = vector::consts::ONES_3;
//!
//! // Create a point of mass 1 kg with given position and speed
//! let point = Point3::inertial(position, speed, 1.);
//!
//! // Create a cluster containing that unique point
//! let mut cluster = Cluster::new(vec![point]);
//!
//! // Initialize solver with dt = 0.1 and 1 iteration per step
//! let mut solver = Solver::new(0.1, 1, Method::RungeKutta4);
//!
//! // Apply an acceleration to the points of the cluster
//! cluster.apply(&mut solver, |points, index| {
//!     // Same constant downwards unit acceleration for all points
//!     Vector6::concat(&points[index].state.speed, &vector::consts::N_EZ_3)
//! });
//! ```
//!
//! # JSON format
//!
//! You can have an idea of what format is required to deserialize JSON stellar systems by looking at
//! [this example](https://github.com/samiBendou/nbodies/blob/master/data/solar_system.json) file.
//! It represents the eight main bodies of the sollar system.
//!

use std::fmt;
use std::fmt::Debug;
use std::ops::{Index, IndexMut};

use geomath;
use geomath::prelude::*;
use geomath::vector::{Vector3, Vector6};

use crate::point::Point3;
use crate::solver::Solver;

/// Helper structures
pub mod common;
/// Numerical constants
pub mod consts;
/// Dynamic points
pub mod point;
/// Commonly used forces
pub mod forces;
/// Commonly used potentials
pub mod potentials;
/// 2-body Kepler's orbit
pub mod orbital;
/// N-body solver
pub mod solver;

/// Cluster of dynamic points
///
/// Encapsulates an array of dynamic points in order to represent and simulate n-body systems.
pub struct Cluster {
    /// Array of dynamic points
    pub points: Vec<Point3>,
    barycenter: Point3,
}

impl Cluster {
    /// Construct a cluster from given array of points
    #[inline]
    pub fn new(points: Vec<Point3>) -> Self {
        Cluster { points, barycenter: Point3::zeros(0.) }
    }

    /// Construct an empty cluster
    #[inline]
    pub fn empty() -> Self {
        Cluster::new(vec![])
    }

    /// Check either the cluster is empty
    ///
    /// An empty cluster does not contain any point.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.points.len() == 0
    }

    /// Number of points in the cluster
    #[inline]
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Barycenter of the cluster
    ///
    /// The barycenter is the center of mass.
    #[inline]
    pub fn barycenter(&self) -> &Point3 {
        &self.barycenter
    }

    /// Total kinetic energy of the cluster
    #[inline]
    pub fn kinetic_energy(&self) -> f64 {
        self.points.iter().map(|point| point.kinetic_energy()).sum()
    }

    /// Total angular momentum of the cluster
    #[inline]
    pub fn angular_momentum(&self) -> f64 {
        self.points.iter().map(|point| point.angular_momentum()).sum()
    }

    /// Potential energy of the cluster from given potential energy
    ///
    /// The function `f` must compute the potential energy of a particular point of the cluster
    pub fn potential_energy<T>(&self, mut f: T) -> f64 where
        T: FnMut(&Vec<Point3>, usize) -> f64 {
        let len = self.points.len();
        let mut ret = 0.;
        for i in 0..len {
            ret += f(&self.points, i);
        }
        ret * 0.5
    }

    /// Push a point into the cluster
    #[inline]
    pub fn push(&mut self, point: Point3) -> &mut Self {
        self.points.push(point);
        self.update_barycenter()
    }

    /// Pop a point from the cluster
    #[inline]
    pub fn pop(&mut self) -> Option<Point3> {
        let ret = self.points.pop();
        self.update_barycenter();
        ret
    }

    /// Remove a point at given cluster's index
    #[inline]
    pub fn remove(&mut self, i: usize) -> Point3 {
        let ret = self.points.remove(i);
        self.update_barycenter();
        ret
    }

    /// Set zero position, speed and trajectory at given cluster's index
    #[inline]
    pub fn reset0_at(&mut self, i: usize) -> &mut Self {
        self.points[i].state.reset0();
        self.points[i].state.trajectory.reset0();
        self.update_barycenter()
    }

    /// Set state and trajectory of point at given cluster's index
    #[inline]
    pub fn reset_at(&mut self, i: usize, state: &geomath::point::Point3) -> &mut Self {
        self.points[i].state.reset(state);
        self.points[i].state.trajectory.reset(&state.position);
        self.update_barycenter()
    }

    /// Set position of point at given cluster's index
    #[inline]
    pub fn reset_position_at(&mut self, i: usize, position: &Vector3) -> &mut Self {
        self.points[i].state.position = *position;
        self.points[i].state.trajectory.reset(position);
        self.update_barycenter()
    }

    /// Set position of point at given cluster's index
    #[inline]
    pub fn reset_speed_at(&mut self, i: usize, speed: &Vector3) -> &mut Self {
        self.points[i].state.speed = *speed;
        self.update_barycenter()
    }

    /// Set trajectory at given cluster's index
    #[inline]
    pub fn reset_trajectory_at(&mut self, i: usize) -> &mut Self {
        let position = self.points[i].state.position;
        self.points[i].state.trajectory.reset(&position);
        self.update_barycenter()
    }

    /// Translate point a given cluster index along direction vector
    #[inline]
    pub fn translate_at(&mut self, i: usize, direction: &Vector3) -> &mut Self {
        self.points[i].state.position += *direction;
        self.points[i].state.update_trajectory();
        self.update_barycenter()
    }

    /// Translate the cluster along given direction vector
    #[inline]
    pub fn translate(&mut self, direction: &Vector3) -> &mut Self {
        self.barycenter.state.position += *direction;
        for body in self.points.iter_mut() {
            body.state.position += *direction;
        }
        self.update_trajectory()
    }

    /// Apply forces to the points of the cluster
    ///
    /// The state and trajectory of the points gets updated according to the equation function
    /// See [`solver`](solver/index.html) module for more details.
    #[inline]
    pub fn apply<T>(&mut self, solver: &mut Solver, f: T) -> &mut Self where
        T: FnMut(&Vec<Point3>, usize) -> Vector6 {
        solver.step(&mut self.points, f);
        self.update_barycenter()
            .update_trajectory()
    }

    /// Reset trajectories of all points of the cluster
    ///
    /// The trajectory of each point is set to its current position
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

    #[inline]
    fn update_barycenter(&mut self) -> &mut Self {
        self.barycenter.mass = 0.;
        self.barycenter.state.reset0();
        for body in self.points.iter() {
            self.barycenter.mass += body.mass;
            self.barycenter.state += body.state * body.mass;
        }
        self.barycenter.state /= self.barycenter.mass;
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