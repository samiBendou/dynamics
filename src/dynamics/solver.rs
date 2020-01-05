use crate::dynamics::Body;
use crate::geometry::common::Vector;
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
    state: Vector4,
    tmp: [Vector4; 4],
}

impl Solver {
    pub fn new(dt: f64, method: Method) -> Solver {
        Solver { dt, method, state: Vector4::zeros(), tmp: [Vector4::zeros(); 4] }
    }

    pub fn step<T>(&mut self, bodies: &mut Vec<Body>, mut f: T, iterations: u32) -> &mut Self where
        T: FnMut(&Vec<Body>, usize) -> Vector4 {
        let len = bodies.len();


        for _ in 0..iterations {
            for i in 0..len {
                match self.method {
                    Method::EulerExplicit => self.euler_explicit(bodies, i, &mut f),
                    Method::RungeKutta4 => self.runge_kutta_4(bodies, i, &mut f),
                };
            }
            for i in 0..len {
                bodies[i].center.accelerate(self.dt);
            }
        }
        self
    }

    pub fn euler_explicit<T>(&mut self, bodies: &mut Vec<Body>, i: usize, f: &mut T) -> &mut Self where
        T: FnMut(&Vec<Body>, usize) -> Vector4 {
        bodies[i].center.gradient = f(bodies, i);
        self
    }
    pub fn runge_kutta_4<T>(&mut self, bodies: &mut Vec<Body>, i: usize, f: &mut T) -> &mut Self where
        T: FnMut(&Vec<Body>, usize) -> Vector4 {
        let half_dt = 0.5 * self.dt;
        self.tmp[0] = f(bodies, i);
        self.state = bodies[i].center.state.vector();
        bodies[i].center.state.set_vector(&(self.state + self.tmp[0] * half_dt));
        self.tmp[1] = f(bodies, i);
        bodies[i].center.state.set_vector(&(self.state + self.tmp[1] * half_dt));
        self.tmp[2] = f(bodies, i);
        bodies[i].center.state.set_vector(&(self.state + self.tmp[2] * self.dt));
        self.tmp[3] = f(bodies, i);
        bodies[i].center.state.set_vector(&self.state);
        bodies[i].center.gradient = (self.tmp[0] + (self.tmp[1] + self.tmp[2]) * 2. + self.tmp[3]) * FRAC_1_6;
        self
    }
}