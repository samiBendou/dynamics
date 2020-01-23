use geomath::prelude::*;
use geomath::vector::Vector6;

use crate::point::Point3;

const FRAC_1_6: f64 = 1. / 6.;

pub enum Method {
    EulerExplicit,
    RungeKutta4,
}

impl Method {
    pub fn next(&mut self) {
        use crate::solver::Method::*;
        *self = match self {
            EulerExplicit => RungeKutta4,
            RungeKutta4 => EulerExplicit,
        }
    }
}

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
    pub fn new(dt: f64, iterations: u32, method: Method) -> Solver {
        Solver { dt, iterations, method, state: Vec::new(), tmp: Vec::new() }
    }

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