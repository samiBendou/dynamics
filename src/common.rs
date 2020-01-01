use rand::Rng;

#[macro_export]
macro_rules! assert_near {
    ($val: expr, $exp: expr, $tol: expr) => {
        assert!(
        ($val - $exp).abs() < $tol,
        "Approximation failed\nvalue: {}\nexpected: {}\n tolerance: {}",
        $val, $exp, $tol
        )
    }
}

pub fn random_color() -> [f32; 4] {
    let mut rng = rand::thread_rng();
    [rng.gen(), rng.gen(), rng.gen(), 1.]
}

