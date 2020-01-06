use std::fmt::{Debug, Error, Formatter};
use std::ops::{
    Add, AddAssign,
    BitOr, Div,
    DivAssign,
    Mul,
    MulAssign,
    Neg, Not,
    Rem, Sub,
    SubAssign,
};

use crate::geometry::common::*;
use crate::impl_vector;

#[derive(Copy, Clone)]
pub struct Matrix2 {
    pub xx: f64,
    pub xy: f64,
    pub yx: f64,
    pub yy: f64,
}

#[derive(Copy, Clone)]
pub struct Matrix3 {
    pub xx: f64,
    pub xy: f64,
    pub xz: f64,
    pub yx: f64,
    pub yy: f64,
    pub yz: f64,
    pub zx: f64,
    pub zy: f64,
    pub zz: f64,
}

impl_vector!(Matrix2 {xx, xy, yx, yy}, 4);
impl_vector!(Matrix3 {xx, xy, xz, yx, yy, yz, zx, zy, zz}, 9);

trait Algebra<T> {
    fn determinant(&self) -> f64;
    fn inverse(&self) -> Self;
    fn transposed(&self) -> Self;
    fn set_inverse(&mut self) -> &mut Self;
    fn set_transposed(&mut self) -> &mut Self;
}

impl MulAssign<Matrix3> for Matrix3 {
    fn mul_assign(&mut self, rhs: Matrix3) {
        self.xx = rhs.xx * self.xx + rhs.yx * self.xy + rhs.zx * self.xz;
        self.yx = rhs.xx * self.yx + rhs.yx * self.yy + rhs.zx * self.yz;
        self.zx = rhs.xx * self.zx + rhs.yx * self.zy + rhs.zx * self.zz;
        self.xy = rhs.xy * self.xx + rhs.yy * self.xy + rhs.zy * self.xz;
        self.yy = rhs.xy * self.yx + rhs.yy * self.yy + rhs.zy * self.yz;
        self.zy = rhs.xy * self.zx + rhs.yy * self.zy + rhs.zy * self.zz;
        self.zx = rhs.xz * self.xx + rhs.yz * self.xy + rhs.zz * self.xz;
        self.zy = rhs.xz * self.yx + rhs.yz * self.yy + rhs.zz * self.yz;
        self.zz = rhs.xz * self.zx + rhs.yz * self.zy + rhs.zz * self.zz;
    }
}

impl Algebra<Matrix3> for Matrix3 {
    fn determinant(&self) -> f64 {
        let dyx = self.zz * self.yy - self.zy * self.yz;
        let dyy = -self.zz * self.xy + self.zy * self.xz;
        let dyz = self.yz * self.xy - self.yy * self.xz;
        self.xx * dyx + self.yx * dyy + self.zx * dyz
    }

    fn inverse(&self) -> Self {
        let mut ret = *self;
        ret.set_inverse();
        ret
    }

    fn transposed(&self) -> Self {
        let mut ret = *self;
        ret.set_transposed();
        ret
    }

    fn set_inverse(&mut self) -> &mut Self {
        let dyx = self.zz * self.yy - self.zy * self.yz;
        let dyy = -self.zz * self.xy + self.zy * self.xz;
        let dyz = self.yz * self.xy - self.yy * self.xz;
        let mut det = self.xx * dyx + self.yx * dyy + self.zx * dyz;

        if det.abs() < std::f64::EPSILON {
            return self;
        }

        det = 1. / det;
        self.xx = dyx * det;
        self.yx = (-self.zz * self.yx + self.zx * self.yz) * det;
        self.zx = (self.zy * self.yx - self.zx * self.yy) * det;
        self.xy = dyy * det;
        self.yy = (self.zz * self.xx - self.zx * self.xz) * det;
        self.zy = (-self.zy * self.xx + self.zx * self.xy) * det;
        self.zx = dyz * det;
        self.zy = (-self.yz * self.xx + self.yx * self.xz) * det;
        self.zz = (self.yy * self.xx - self.yx * self.xy) * det;
        self
    }

    fn set_transposed(&mut self) -> &mut Self {
        let yx = self.yx;
        let zx = self.zx;
        let yz = self.zy;
        self.yx = self.xy;
        self.xy = yx;
        self.zx = self.xz;
        self.xz = zx;
        self.yz = self.zy;
        self.zy = yz;
        self
    }
}


#[cfg(test)]
mod tests {
    mod matrix3 {}
}