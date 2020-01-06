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

impl Mul<Matrix3> for Matrix3 {
    type Output = Matrix3;

    fn mul(self, rhs: Matrix3) -> Self::Output {
        let mut ret = self;
        ret *= rhs;
        ret
    }
}

impl MulAssign<Matrix3> for Matrix3 {
    fn mul_assign(&mut self, rhs: Matrix3) {
        let xx = self.xx;
        let yx = self.yx;
        let zx = self.zx;
        let xy = self.xy;
        let yy = self.yy;
        let zy = self.zy;
        let xz = self.xz;
        let yz = self.yz;
        let zz = self.zz;
        self.xx = rhs.xx * xx + rhs.yx * xy + rhs.zx * xz;
        self.yx = rhs.xx * yx + rhs.yx * yy + rhs.zx * yz;
        self.zx = rhs.xx * zx + rhs.yx * zy + rhs.zx * zz;
        self.xy = rhs.xy * xx + rhs.yy * xy + rhs.zy * xz;
        self.yy = rhs.xy * yx + rhs.yy * yy + rhs.zy * yz;
        self.zy = rhs.xy * zx + rhs.yy * zy + rhs.zy * zz;
        self.xz = rhs.xz * xx + rhs.yz * xy + rhs.zz * xz;
        self.yz = rhs.xz * yx + rhs.yz * yy + rhs.zz * yz;
        self.zz = rhs.xz * zx + rhs.yz * zy + rhs.zz * zz;
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
    mod matrix3 {
        use crate::geometry::common::*;

        use super::super::Matrix3;

        #[test]
        fn arithmetic() {
            let a = Matrix3::new(2., 0., 0., 0., 2., 0., 0., 0., 2.);
            let b = Matrix3::ones();
            let c = a * b;
            assert_eq!(c, b * 2.);
        }
    }
}