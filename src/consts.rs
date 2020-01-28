//!
//! This modules defines all the commonly used numerical constants in physics. Two types of constants are
//! currently implemented :
//! * Fundamental constants
//! * Electromagnetic constants

pub const G_UNIV: f64 = 6.67430e-11;
pub const H_PLANCK: f64 = 6.62607015e-34;
pub const C_LIGHT: f64 = 299792458.;

pub const MU_0: f64 = 4. * std::f64::consts::PI * 1e-7;
pub const EPSILON_0: f64 = 1. / (MU_0 * C_LIGHT * C_LIGHT);