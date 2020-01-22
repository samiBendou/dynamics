use std::fmt;
use std::fmt::Debug;
use std::ops::{AddAssign, DivAssign};

use geomath::trajectory::consts::TRAJECTORY_SIZE;
use rand::Rng;
use serde::{de, de::Visitor, Deserialize, Deserializer, Serialize, Serializer};

pub const AVERAGE_SIZE: usize = TRAJECTORY_SIZE - 1;

pub fn random_color() -> [f32; 4] {
    let mut rng = rand::thread_rng();
    [rng.gen(), rng.gen(), rng.gen(), 1.]
}

#[derive(Copy, Clone)]
pub struct Average<T> {
    index: usize,
    values: [T; AVERAGE_SIZE],
    average: T,
}

impl<T> Average<T> where
    T: Copy + AddAssign<T> + DivAssign<f64> + Serialize + Deserialize {
    pub fn new(val: &T) -> Average<T> {
        Average { index: 0, values: [*val; AVERAGE_SIZE], average: *val }
    }

    pub fn push(&mut self, val: &T) -> &mut Self {
        self.average = *val;
        for val in self.values.iter() {
            self.average += *val;
        }
        self.average /= (AVERAGE_SIZE + 1) as f64;
        self.values[self.index] = *val;
        self.index = (self.index + 1) % AVERAGE_SIZE;
        self
    }

    pub fn get(&self) -> &T {
        &self.average
    }
}

impl<T> Debug for Average<T> where
    T: Debug {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        write!(f, "{:?}", self.average)
    }
}

impl Serialize for Average<f64> {
    fn serialize<S>(&self, serializer: S) -> Result<<S as Serializer>::Ok, <S as Serializer>::Error> where S: Serializer {
        serializer.serialize_f64(self.average)
    }
}

impl Deserialize for Average<f64> {
    fn deserialize<D>(deserializer: D) -> Result<Average<f64>, D::Error> where
        D: Deserializer {
        Ok(deserializer.deserialize_f64(AverageF64Visitor)?)
    }
}


struct AverageF64Visitor;

impl Visitor for AverageF64Visitor {
    type Value = Average<f64>;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("a double ")
    }

    fn visit_f64<E>(self, value: f64) -> Result<Self::Value, E>
        where
            E: de::Error,
    {
        Ok(Average::new(&value))
    }
}