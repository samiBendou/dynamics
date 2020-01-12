use std::error::Error;
use std::fs::File;
use std::io::prelude::*;
use std::ops::{Index, IndexMut};
use std::path::Path;

use rand::prelude::*;
use serde::{Deserialize, Serialize};

use crate::common::{Average, random_color};
use crate::dynamics::point::Point3;
use crate::geometry::common::*;
use crate::geometry::common::coordinates::{Cartesian2, Cartesian3, Polar};
use crate::geometry::common::transforms::Rotation3;
use crate::geometry::matrix::Matrix3;
use crate::geometry::trajectory::{Trajectory3, TRAJECTORY_SIZE};
use crate::geometry::vector::Vector3;
use crate::units::consts::G_UNIV;

#[derive(Serialize, Deserialize, Debug, PartialEq, Copy, Clone)]
pub enum Kind {
    Artificial,
    Terrestrial,
    Giant,
    Star,
    Hole,
}

impl Kind {
    pub fn random() -> Kind {
        use Kind::*;
        let mut rng = rand::thread_rng();
        match rng.gen_range(0, 4) {
            1 => Terrestrial,
            2 => Giant,
            3 => Star,
            4 => Hole,
            _ => Artificial,
        }
    }

    pub fn random_mass(&self) -> f64 {
        let mut rng = rand::thread_rng();
        match self {
            Kind::Artificial => rng.gen_range(1., 1e6),
            Kind::Terrestrial => rng.gen_range(1e22, 1e25),
            Kind::Giant => rng.gen_range(1e25, 1e28),
            Kind::Star => rng.gen_range(1e28, 1e31),
            Kind::Hole => rng.gen_range(1e32, 1e31),
        }
    }


    pub fn random_radius(&self) -> f64 {
        let mut rng = rand::thread_rng();
        match self {
            Kind::Artificial => rng.gen_range(1., 100.),
            Kind::Terrestrial => rng.gen_range(1e6, 1e7),
            Kind::Giant => rng.gen_range(1e7, 1e8),
            Kind::Star => rng.gen_range(1e7, 1e9),
            Kind::Hole => rng.gen_range(1e6, 1e10),
        }
    }

    pub fn scaled_radius(&self, radius: f64) -> f64 {
        radius / 10f64.powf(radius.log10()) + match self {
            Kind::Artificial => 1.,
            Kind::Terrestrial => 20.,
            Kind::Giant => 22.,
            Kind::Star => 24.,
            Kind::Hole => 24.,
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Copy, Clone)]
pub struct Inclination {
    pub value: Average<f64>,
    pub argument: Average<f64>,
}

impl Inclination {
    pub fn zeros() -> Self {
        Inclination { value: Average::new(&0.0), argument: Average::new(&0.0) }
    }
}

#[derive(Serialize, Deserialize, Debug, Copy, Clone)]
pub struct Orbit {
    pub mu: f64,
    pub apoapsis: Average<f64>,
    pub periapsis: Average<f64>,
    pub argument: Average<f64>,
    pub inclination: Inclination,
}

impl Orbit {
    pub fn zeros() -> Self {
        Orbit {
            mu: 0.0,
            apoapsis: Average::new(&0.0),
            periapsis: Average::new(&0.0),
            argument: Average::new(&0.0),
            inclination: Inclination::zeros(),
        }
    }

    pub fn set_interpolation(&mut self, trajectory: &Trajectory3, barycenter: &Point3) -> &Self {
        let rad_to_deg = 180. / std::f64::consts::PI;
        let pi_frac_2 = std::f64::consts::FRAC_PI_2;
        let last_index = TRAJECTORY_SIZE - 1;
        let mut normal0 = (trajectory[last_index] - trajectory[last_index - 1])
            .cross(&(trajectory[last_index] - barycenter.state.trajectory[last_index]));
        normal0.normalize();
        let position0 = normal0.cross(&Vector3::unit_y());
        if normal0.magnitude() < std::f64::EPSILON {
            return self;
        }
        let mut distances: Vec<f64> = Vec::with_capacity(TRAJECTORY_SIZE);
        let mut anomalies: Vec<f64> = Vec::with_capacity(TRAJECTORY_SIZE);
        let mut inclinations: Vec<f64> = Vec::with_capacity(TRAJECTORY_SIZE);
        let mut position: Vector3;
        let mut normal: Vector3;
        for i in 0..TRAJECTORY_SIZE {
            position = trajectory[i] - barycenter.state.trajectory[i];
            normal = position.cross(&position0);
            normal.normalize();
            inclinations.push(position.angle(&Vector3::unit_z()));
            distances.push(position.magnitude());
            if normal | normal0 > 0. {
                anomalies.push(position.angle(&position0));
            } else {
                anomalies.push(2. * std::f64::consts::PI - position.angle(&position0));
            }
        }
        let (_, max_magnitude) = self.interpolate(&distances, &anomalies, |rhs, lhs| rhs < lhs);
        let (min_argument, min_magnitude) = self.interpolate(&distances, &anomalies, |rhs, lhs| rhs > lhs);
        let (asc_argument, asc_inclination) = self.interpolate(&inclinations, &anomalies, |rhs, lhs| rhs > lhs);

        self.argument.push(&((min_argument) * rad_to_deg));
        self.apoapsis.push(&max_magnitude);
        self.periapsis.push(&min_magnitude);
        self.inclination.argument.push(&((asc_argument - pi_frac_2) * rad_to_deg));
        self.inclination.value.push(&((pi_frac_2 - asc_inclination) * rad_to_deg));
        self.mu = G_UNIV * barycenter.mass;
        self
    }

    pub fn semi_minor(&self) -> f64 {
        (self.apoapsis.get() * self.periapsis.get()).sqrt()
    }

    pub fn semi_major(&self) -> f64 {
        0.5 * (self.apoapsis.get() + self.periapsis.get())
    }

    pub fn is_degenerated(&self) -> bool {
        let epsilon_f64 = std::f64::EPSILON;
        if self.semi_minor() < epsilon_f64 || self.semi_major() < epsilon_f64 {
            true
        } else {
            false
        }
    }

    pub fn eccentricity(&self) -> f64 {
        let ra = self.apoapsis.get();
        let rp = self.periapsis.get();
        let total_r = ra + rp;
        if total_r > 0. {
            (ra - rp) / total_r
        } else {
            0.
        }
    }

    pub fn radius_at(&self, true_anomaly: f64) -> f64 {
        let a = self.semi_major();
        let epsilon = self.eccentricity();
        a * (1. - epsilon * epsilon) / (1. + epsilon * true_anomaly.cos())
    }

    pub fn eccentric_anomaly_at(&self, true_anomaly: f64) -> f64 {
        let epsilon = self.eccentricity();
        (true_anomaly.sin() * (1. - epsilon * epsilon).sqrt() / (1. + epsilon * true_anomaly.cos())).atan()
    }

    pub fn flight_angle_at(&self, true_anomaly: f64) -> f64 {
        let epsilon = self.eccentricity();
        let ec = epsilon * true_anomaly.cos();
        ((1. + ec) / (1. + epsilon * epsilon + 2. * ec).sqrt()).min(1.0).acos()
    }

    pub fn position_at(&self, true_anomaly: f64) -> Vector3 {
        let rad_to_deg = std::f64::consts::PI / 180.;
        let argument = *self.inclination.argument.get() * rad_to_deg;
        let value = *self.inclination.value.get() * rad_to_deg;
        let mag = self.radius_at(true_anomaly);
        let rotation_node = Matrix3::from_rotation(value, Vector3::unit_x().set_rotation_z(argument));
        let rotation_z = Matrix3::from_rotation_z(*self.argument.get() * rad_to_deg);
        rotation_node * (rotation_z * Vector3::from_polar(mag, true_anomaly))
    }

    pub fn speed_at(&self, true_anomaly: f64) -> Vector3 {
        if self.is_degenerated() {
            return Vector3::zeros();
        }
        let pi_frac_2 = std::f64::consts::FRAC_PI_2;
        let rad_to_deg = std::f64::consts::PI / 180.;
        let argument = *self.inclination.argument.get() * rad_to_deg;
        let value = *self.inclination.value.get() * rad_to_deg;
        let mag = (self.mu * (2. / self.radius_at(true_anomaly) - 1. / self.semi_major())).sqrt();
        let phi = true_anomaly + pi_frac_2 - self.flight_angle_at(true_anomaly);
        let rotation_node = Matrix3::from_rotation(value, Vector3::unit_x().set_rotation_z(argument));
        let rotation_z = Matrix3::from_rotation_z(*self.argument.get() * rad_to_deg);
        rotation_node * (rotation_z * Vector3::from_polar(mag, phi))
    }

    fn interpolate<T>(&self, values: &Vec<f64>, angles: &Vec<f64>, mut f: T) -> (f64, f64) where
        T: FnMut(f64, f64) -> bool {
        let mut value = values[0];
        let mut index = 0;
        for i in 1..TRAJECTORY_SIZE - 1 {
            if angles[i] == std::f64::NAN {
                continue;
            }
            if f(value, values[i]) {
                value = values[i];
                index = i;
            }
        }
        (angles[index], values[index])
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Body {
    pub name: String,
    pub mass: f64,
    pub kind: Kind,
    pub color: [f32; 4],
    pub radius: f64,
    pub orbit: Orbit,
}

impl Body {
    pub fn new() -> Self {
        Body {
            name: "untitled".to_string(),
            mass: 1.0,
            kind: Kind::Artificial,
            color: [1., 0., 0., 1.],
            radius: 1e5,
            orbit: Orbit::zeros(),
        }
    }

    pub fn random() -> Self {
        let mut rng = rand::thread_rng();
        let kind = Kind::random();
        let mass = kind.random_mass();
        let radius = kind.random_radius();
        let id = rng.gen_range(0, std::i16::MAX);
        Body {
            name: format!("{:#?}-{}", kind, id),
            mass,
            kind,
            color: random_color(),
            radius,
            orbit: Orbit::zeros(),
        }
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Cluster {
    pub bodies: Vec<Body>
}

impl From<Vec<Body>> for Cluster {
    fn from(bodies: Vec<Body>) -> Self {
        Cluster { bodies }
    }
}

impl Cluster {
    pub fn from_file(path: &Path) -> Result<Self, Box<dyn Error>> {
        let mut file = File::open(path)?;
        let mut contents = String::new();
        file.read_to_string(&mut contents)?;
        let bodies: Vec<Body> = serde_json::from_str(&contents)?;
        Ok(Cluster::from(bodies))
    }

    pub fn update_orbits(&mut self, points: &Vec<Point3>, barycenter: &Point3) -> &mut Self {
        for i in 0..points.len() {
            self.bodies[i].orbit.set_interpolation(&points[i].state.trajectory, barycenter);
        }
        self
    }

    #[inline]
    pub fn push(&mut self, body: Body) -> &mut Self {
        self.bodies.push(body);
        self
    }

    #[inline]
    pub fn pop(&mut self) -> Option<Body> {
        self.bodies.pop()
    }

    #[inline]
    pub fn remove(&mut self, i: usize) -> Body {
        self.bodies.remove(i)
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

#[cfg(test)]
mod tests {
    mod cluster {
        use std::path::Path;

        use super::super::Cluster;

        #[test]
        fn simple_deserialize() {
            let path: &Path = Path::new("data/solar_system.json");
            let cluster = Cluster::from_file(path);
            println!("{:?}", cluster);
        }
    }
}

