use crate::geometry::vector::Vector3;

pub trait Initializer {
    fn zeros() -> Self;
    fn ones() -> Self;
}

pub trait Reset<T> {
    fn reset0(&mut self) -> &mut Self;
    fn reset1(&mut self) -> &mut Self;
    fn reset(&mut self, val: &T) -> &mut Self;
}

pub trait Array<T> {
    fn array(&self) -> T;
    fn set_array(&mut self, arr: &T) -> &mut Self;
}

pub trait Vector<T> {
    fn vector(&self) -> T;
    fn set_vector(&mut self, arr: &T) -> &mut Self;
}

pub trait Split<T> {
    fn split(&self) -> [T; 2];
    fn concat(lhs: &T, rhs: &T) -> Self;
    fn upper(&self) -> T;
    fn lower(&self) -> T;
    fn set_upper(&mut self, val: &T) -> &mut Self;
    fn set_lower(&mut self, val: &T) -> &mut Self;
}

pub trait Metric {
    fn dot(&self, other: &Self) -> f64;
    fn distance2(&self, other: &Self) -> f64;
    fn distance(&self, other: &Self) -> f64;
    fn magnitude2(&self) -> f64;
    fn magnitude(&self) -> f64;
    fn normalize(&mut self) -> &mut Self;
}

pub trait Angle {
    fn cos(&self, _: &Self) -> f64;
    fn sin(&self, _: &Self) -> f64;
    fn angle(&self, _: &Self) -> f64;
    fn area(&self, _: &Self) -> f64;
    fn cross(&self, _: &Self) -> Vector3;
}