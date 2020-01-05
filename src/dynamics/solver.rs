use crate::dynamics::Body;
use crate::geometry::common::*;
use crate::geometry::vector::Vector4;

const FRAC_1_6: f64 = 1. / 6.;

pub enum Method {
    EulerExplicit,
    RungeKutta4,
}

impl Method {
    pub fn next(&mut self) {
        use crate::dynamics::solver::Method::*;
        *self = match self {
            EulerExplicit => RungeKutta4,
            RungeKutta4 => EulerExplicit,
        }
    }
}


pub struct Solver {
    pub dt: f64,
    pub method: Method,
    state: Vec<Vector4>,
    tmp: Vec<[Vector4; 4]>,
}

impl Solver {
    pub fn new(dt: f64, method: Method) -> Solver {
        Solver { dt, method, state: Vec::new(), tmp: Vec::new() }
    }

    pub fn step<T>(&mut self, bodies: &mut Vec<Body>, mut f: T, iterations: u32) -> &mut Self where
        T: FnMut(&Vec<Body>, usize) -> Vector4 {
        self.tmp = bodies.iter().map(|_body| [Vector4::zeros(); 4]).collect();
        self.state = bodies.iter().map(|_body| Vector4::zeros()).collect();

        for _ in 0..iterations {
            match self.method {
                Method::EulerExplicit => self.euler_explicit(bodies, &mut f),
                Method::RungeKutta4 => self.runge_kutta_4(bodies, &mut f),
            };
            for body in bodies.iter_mut() {
                body.center.accelerate(self.dt);
            }
        }
        self
    }

    fn euler_explicit<T>(&mut self, bodies: &mut Vec<Body>, f: &mut T) -> &mut Self where
        T: FnMut(&Vec<Body>, usize) -> Vector4 {
        for j in 0..bodies.len() {
            bodies[j].center.gradient = f(bodies, j);
        }
        self
    }

    fn runge_kutta_4<T>(&mut self, bodies: &mut Vec<Body>, f: &mut T) -> &mut Self where
        T: FnMut(&Vec<Body>, usize) -> Vector4 {
        let half_dt = 0.5 * self.dt;
        for j in 0..bodies.len() {
            self.tmp[j][0] = f(bodies, j);
            bodies[j].center.gradient = self.tmp[j][0];
        }
        for j in 0..bodies.len() {
            self.state[j] = bodies[j].center.state.vector();
            bodies[j].center.state.set_vector(&(self.state[j] + self.tmp[j][0] * half_dt));
        }
        for j in 0..bodies.len() {
            self.tmp[j][1] = f(bodies, j);
        }
        for j in 0..bodies.len() {
            bodies[j].center.gradient += self.tmp[j][1] * 2.;
            bodies[j].center.state.set_vector(&(self.state[j] + self.tmp[j][1] * half_dt));
        }
        for j in 0..bodies.len() {
            self.tmp[j][2] = f(bodies, j);
        }
        for j in 0..bodies.len() {
            bodies[j].center.gradient += self.tmp[j][2] * 2.;
            bodies[j].center.state.set_vector(&(self.state[j] + self.tmp[j][2] * self.dt));
        }
        for j in 0..bodies.len() {
            self.tmp[j][3] = f(bodies, j);
        }
        for j in 0..bodies.len() {
            bodies[j].center.gradient += self.tmp[j][3];
            bodies[j].center.gradient *= FRAC_1_6;
            bodies[j].center.state.set_vector(&self.state[j]);
        }
        self
    }
}