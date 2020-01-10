use std::error::Error;
use std::fs::File;
use std::io::prelude::*;
use std::ops::{Index, IndexMut};
use std::path::Path;

use rand::prelude::*;
use serde::{Deserialize, Serialize};

use crate::common::random_color;
use crate::geometry::common::*;
use crate::geometry::common::coordinates::{Cartesian2, Spherical};
use crate::geometry::trajectory::{Trajectory3, TRAJECTORY_SIZE};
use crate::geometry::vector::Vector3;

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

#[derive(Serialize, Deserialize, Debug, PartialEq, Copy, Clone)]
pub struct Inclination {
    pub value: f64,
    pub argument: f64,
}

impl Inclination {
    pub fn zeros() -> Self {
        Inclination { value: 0.0, argument: 0.0 }
    }
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Copy, Clone)]
pub struct Orbit {
    pub mu: f64,
    pub apoapsis: f64,
    pub periapsis: f64,
    pub argument: f64,
    pub inclination: Inclination,
}

impl Orbit {
    pub fn zeros() -> Self {
        Orbit {
            mu: 0.0,
            apoapsis: 0.0,
            periapsis: 0.0,
            argument: 0.0,
            inclination: Inclination::zeros(),
        }
    }

    pub fn interpolate_trajectory(&self, trajectory: &Trajectory3, center: &Vector3) -> &Self {
        let position0 = trajectory[0] - *center;
        let magnitudes: Vec<f64> = trajectory.positions().iter()
            .map(|position| position.distance(center))
            .collect();
        let angles: Vec<f64> = trajectory.positions().iter()
            .map(|position| (*position - *center).angle(&Vector3::unit_x()))
            .collect();

        let mut max_indexes: Vec<usize> = Vec::new();
        let mut min_indexes: Vec<usize> = Vec::new();
        let mut max_magnitude = magnitudes[0];
        let mut min_magnitude = magnitudes[0];
        let mut sign = if angles[1] > angles[0] { true } else { false };
        let mut sign_changes = 0;
        let mut index = 0;
        max_indexes.push(0);
        min_indexes.push(0);
        for i in 1..TRAJECTORY_SIZE - 1 {
            println!("r({:.3}) = {:.3e}", angles[i] * 180. / std::f64::consts::PI, magnitudes[i]);
            if max_magnitude < magnitudes[i] {
                max_magnitude = magnitudes[i];
                max_indexes[index] = i;
            }
            if min_magnitude > magnitudes[i] {
                min_magnitude = magnitudes[i];
                min_indexes[index] = i;
            }
            if (angles[i] > angles[i - 1]) != sign {
                sign_changes += 1;
                sign = !sign;
            }

            if sign_changes == 2 {
                max_indexes.push(0);
                min_indexes.push(0);
                max_magnitude = 0.;
                min_magnitude = std::f64::INFINITY;
                index += 1;
                sign_changes = 0;
            }
        }
        let mut average_max_argument = 0.;
        let mut average_max_magnitude = 0.;
        for index in max_indexes.iter() {
            average_max_argument += angles[*index];
            average_max_magnitude += magnitudes[*index];
        }
        average_max_argument /= max_indexes.len() as f64;
        average_max_magnitude /= max_indexes.len() as f64;

        let mut average = Vector3::zeros();
        for i in 1..TRAJECTORY_SIZE - 1 {
            average += trajectory[i];
        }
        average /= TRAJECTORY_SIZE as f64;
        println!("avg arg:{:.2}, mag:{:.3e}", average_max_argument * 180. / std::f64::consts::PI, average_max_magnitude);
        self
    }

    pub fn semi_minor(&self) -> f64 {
        (self.apoapsis * self.periapsis).sqrt()
    }

    pub fn semi_major(&self) -> f64 {
        0.5 * (self.apoapsis + self.periapsis)
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
        let ra = self.apoapsis;
        let rp = self.periapsis;
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
        let pi_frac_2 = std::f64::consts::FRAC_PI_2;
        let argument = self.inclination.argument * std::f64::consts::PI / 180.;
        let value = self.inclination.value * std::f64::consts::PI / 180.;
        let mag = self.radius_at(true_anomaly);
        let theta = pi_frac_2 - value * (true_anomaly - argument) / pi_frac_2;
        Vector3::from_spherical(mag, true_anomaly + self.argument, theta)
    }

    pub fn speed_at(&self, true_anomaly: f64) -> Vector3 {
        if self.is_degenerated() {
            return Vector3::zeros();
        }
        let pi_frac_2 = std::f64::consts::FRAC_PI_2;
        let argument = self.inclination.argument * std::f64::consts::PI / 180.;
        let value = self.inclination.value * std::f64::consts::PI / 180.;
        let mag = (self.mu * (2. / self.radius_at(true_anomaly) - 1. / self.semi_major())).sqrt();
        let phi = true_anomaly + std::f64::consts::FRAC_PI_2 - self.flight_angle_at(true_anomaly);
        let theta = pi_frac_2 - value * (1. - (true_anomaly - argument) / pi_frac_2);
        Vector3::from_spherical(mag, phi + self.argument, theta)
    }
}

#[derive(Serialize, Deserialize, Debug, PartialEq)]
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
        let id = rng.gen_range(0, std::i32::MAX);
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

