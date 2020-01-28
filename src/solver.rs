//!
//! A solver implements iterative methods for approximated ODE resolution between dynamic points.
//!
//! It encapsulates the solving step and the number of iterations in order to allow
//! the implementation of variable step methods.
//!
//! The model equation is **u' = f(u, u', t)** where **u = u(t)** is the state of a point
//! * The vector **u(t)** is represented as a dynamic point
//! * The function **f** is specific to each point but can represent interactions between theses
//!
//! ## Conventions

//!
//! The solver uses 2nd to 1rst order transformations for ODE systems.
//! This means, when the point has N degrees of freedom, to represent the equation **a = g**
//! the equation must be between vectors of size 2N:
//! * the first half of the components contains the speed of the point
//! * the last half of the components contains the N-sized **g** vector felt by the point
//!
//! Therefore **u** must also be of size 2N:
//! * The first half of the components is the position
//! * The second half of the components is the speed

use geomath::prelude::*;
use geomath::vector::Vector6;

use crate::point::Point3;

const FRAC_1_6: f64 = 1. / 6.;

/// Solving method provided
pub enum Method {
    /// Euler explicit method
    EulerExplicit,
    /// 4th order Runge Kutta method
    RungeKutta4,
}

impl Method {
    /// Method switching logic
    pub fn next(&mut self) {
        use crate::solver::Method::*;
        *self = match self {
            EulerExplicit => RungeKutta4,
            RungeKutta4 => EulerExplicit,
        }
    }
}

/// 3D Solver
pub struct Solver {
    pub dt: f64,
    pub iterations: u32,
    pub method: Method,
    state: Vec<Vector6>,
    tmp: Vec<[Vector6; 4]>,
}

impl From<Method> for Solver {
    fn from(method: Method) -> Self {
        Solver::new(1., 1, method)
    }
}

impl Solver {
    /// Construct a solver given step, number of iterations and solving method
    pub fn new(dt: f64, iterations: u32, method: Method) -> Solver {
        Solver { dt, iterations, method, state: Vec::new(), tmp: Vec::new() }
    }

    /// Updates the given points according to the given ODE
    ///
    /// * The points passed are supposed already initialized
    /// * The parameter `f` computes the value of the **f** function
    pub fn step<T>(&mut self, points: &mut Vec<Point3>, mut f: T) -> &mut Self where
        T: FnMut(&Vec<Point3>, usize) -> Vector6 {
        self.tmp = points.iter().map(|_point| [Vector6::zeros(); 4]).collect();
        self.state = points.iter().map(|_point| Vector6::zeros()).collect();
        for _ in 0..self.iterations {
            match self.method {
                Method::EulerExplicit => self.euler_explicit(points, &mut f),
                Method::RungeKutta4 => self.runge_kutta_4(points, &mut f),
            };
            for point in points.iter_mut() {
                point.accelerate(self.dt);
            }
        }
        self
    }

    fn euler_explicit<T>(&mut self, points: &mut Vec<Point3>, f: &mut T) -> &mut Self where
        T: FnMut(&Vec<Point3>, usize) -> Vector6 {
        for j in 0..points.len() {
            points[j].gradient = f(points, j);
        }
        self
    }

    fn runge_kutta_4<T>(&mut self, points: &mut Vec<Point3>, f: &mut T) -> &mut Self where
        T: FnMut(&Vec<Point3>, usize) -> Vector6 {
        let half_dt = 0.5 * self.dt;
        for j in 0..points.len() {
            self.tmp[j][0] = f(points, j);
            points[j].gradient = self.tmp[j][0];
        }
        for j in 0..points.len() {
            self.state[j] = points[j].state.to_vector();
            points[j].state.set_vector(&(self.state[j] + self.tmp[j][0] * half_dt));
        }
        for j in 0..points.len() {
            self.tmp[j][1] = f(points, j);
        }
        for j in 0..points.len() {
            points[j].gradient += self.tmp[j][1] * 2.;
            points[j].state.set_vector(&(self.state[j] + self.tmp[j][1] * half_dt));
        }
        for j in 0..points.len() {
            self.tmp[j][2] = f(points, j);
        }
        for j in 0..points.len() {
            points[j].gradient += self.tmp[j][2] * 2.;
            points[j].state.set_vector(&(self.state[j] + self.tmp[j][2] * self.dt));
        }
        for j in 0..points.len() {
            self.tmp[j][3] = f(points, j);
        }
        for j in 0..points.len() {
            points[j].gradient += self.tmp[j][3];
            points[j].gradient *= FRAC_1_6;
            points[j].state.set_vector(&self.state[j]);
        }
        self
    }
}