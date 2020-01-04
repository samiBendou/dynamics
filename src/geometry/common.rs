use crate::geometry::vector::Vector3;

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
    fn set_upper(&mut self, vector: &T) -> &mut Self;
    fn set_lower(&mut self, vector: &T) -> &mut Self;
}

pub trait Angle {
    fn cos(&self, _: &Self) -> f64;
    fn sin(&self, _: &Self) -> f64;
    fn angle(&self, _: &Self) -> f64;
    fn area(&self, _: &Self) -> f64;
    fn cross(&self, _: &Self) -> Vector3;
}