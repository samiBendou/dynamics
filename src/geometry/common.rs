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


pub mod coordinates {
    pub trait Cartesian2 {
        fn unit_x() -> Self;
        fn unit_neg_x() -> Self;
        fn unit_y() -> Self;
        fn unit_neg_y() -> Self;
    }

    pub trait Cartesian3 {
        fn unit_z() -> Self;
        fn unit_neg_z() -> Self;
    }

    pub trait Polar {
        fn from_polar(rho: f64, phi: f64) -> Self;
        fn set_polar(&mut self, rho: f64, phi: f64) -> &mut Self;
        fn unit_rho(phi: f64) -> Self;
        fn unit_phi(phi: f64) -> Self;
        fn rho(&self) -> f64;
        fn phi(&self) -> f64;
        fn set_rho(&mut self, rho: f64) -> &mut Self;
        fn set_phi(&mut self, phi: f64) -> &mut Self;
    }

    pub trait Cylindrical {
        fn from_cylindrical(rho: f64, phi: f64, z: f64) -> Self;
        fn set_cylindrical(&mut self, rho: f64, phi: f64, z: f64) -> &mut Self;
    }

    pub trait Spherical {
        fn from_spherical(radius: f64, phi: f64, theta: f64) -> Self;
        fn set_spherical(&mut self, radius: f64, phi: f64, theta: f64) -> &mut Self;
        fn unit_radius(phi: f64, theta: f64) -> Self;
        fn unit_theta(phi: f64, theta: f64) -> Self;
        fn theta(&self) -> f64;
        fn set_theta(&mut self, theta: f64) -> &mut Self;
    }

    pub trait Homogeneous<T> {
        fn from_homogeneous(vector: &T) -> Self;
        fn homogeneous(&self) -> T;
    }
}

pub mod transforms {
    use crate::geometry::common::Initializer;
    use crate::geometry::vector::{Vector2, Vector3};

    pub trait Cartesian2 where
        Self: std::marker::Sized + Copy + Clone + Initializer {
        fn from_left_up(middle: &Vector2, scale: f64) -> Self {
            let mut ret = Self::zeros();
            ret.set_left_up(middle, scale);
            ret
        }
        fn from_centered(middle: &Vector2, scale: f64) -> Self {
            let mut ret = Self::zeros();
            ret.set_centered(middle, scale);
            ret
        }

        fn left_up(&self, middle: &Vector2, scale: f64) -> Self {
            let mut ret = *self;
            ret.set_left_up(middle, scale);
            ret
        }
        fn centered(&self, middle: &Vector2, scale: f64) -> Self {
            let mut ret = *self;
            ret.set_centered(middle, scale);
            ret
        }
        fn set_left_up(&mut self, middle: &Vector2, scale: f64) -> &mut Self;

        fn set_centered(&mut self, middle: &Vector2, scale: f64) -> &mut Self;
    }

    pub trait Translation<T> where
        Self: std::marker::Sized + Copy + Clone + Initializer {
        fn from_translation(vector: &T) -> Self {
            let mut ret = Self::zeros();
            ret.set_translation(vector);
            ret
        }
        fn translation(&self, vector: &T) -> Self {
            let mut ret = *self;
            ret.set_translation(vector);
            ret
        }
        fn set_translation(&mut self, vector: &T) -> &mut Self;
    }

    pub trait Rigid<U, T> where
        Self: std::marker::Sized + Copy + Clone + Initializer {
        fn from_rigid(rotation: &U, vector: &T) -> Self {
            let mut ret = Self::zeros();
            ret.set_rigid(rotation, vector);
            ret
        }
        fn rigid(&self, rotation: &U, vector: &T) -> Self {
            let mut ret = *self;
            ret.set_rigid(rotation, vector);
            ret
        }
        fn set_rigid(&mut self, rotation: &U, vector: &T) -> &mut Self;
    }

    pub trait Similarity<U, T> where
        Self: std::marker::Sized + Copy + Clone + Initializer {
        fn from_similarity(scale: f64, rotation: &U, vector: &T) -> Self {
            let mut ret = Self::zeros();
            ret.set_similarity(scale, rotation, vector);
            ret
        }
        fn similarity(&self, scale: f64, rotation: &U, vector: &T) -> Self {
            let mut ret = *self;
            ret.set_similarity(scale, rotation, vector);
            ret
        }
        fn set_similarity(&mut self, scale: f64, rotation: &U, vector: &T) -> &mut Self;
    }

    pub trait Rotation2 where
        Self: std::marker::Sized + Copy + Clone + Initializer {
        fn from_rotation(angle: f64) -> Self {
            let mut ret = Self::zeros();
            ret.set_rotation(angle);
            ret
        }
        fn rotation(&self, angle: f64) -> Self {
            let mut ret = *self;
            ret.set_rotation(angle);
            ret
        }
        fn set_rotation(&mut self, angle: f64) -> &mut Self;
    }

    pub trait Rotation3 where
        Self: std::marker::Sized + Copy + Clone + Initializer {
        fn from_rotation(angle: f64, axis: &Vector3) -> Self {
            let mut ret = Self::zeros();
            ret.set_rotation(angle, axis);
            ret
        }
        fn from_rotation_x(angle: f64) -> Self {
            let mut ret = Self::zeros();
            ret.set_rotation_x(angle);
            ret
        }
        fn from_rotation_y(angle: f64) -> Self {
            let mut ret = Self::zeros();
            ret.set_rotation_y(angle);
            ret
        }
        fn from_rotation_z(angle: f64) -> Self {
            let mut ret = Self::zeros();
            ret.set_rotation_z(angle);
            ret
        }

        fn set_rotation(&mut self, angle: f64, axis: &Vector3) -> &mut Self;
        fn set_rotation_x(&mut self, angle: f64) -> &mut Self;
        fn set_rotation_y(&mut self, angle: f64) -> &mut Self;
        fn set_rotation_z(&mut self, angle: f64) -> &mut Self;
    }
}