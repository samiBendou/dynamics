# dynamics
### A realistic physics engine

## Brief
dynamics is a framework that aims to provide common features needed when dealing with physics.
It has first been tough to be an efficient solver for n-body physics. 

It provides a physics engine that applies the 2nd law of Newton on a set of points with mass,
You specify the forces for each, and the points gets updated.

The framework also provides a helper module made to manipulate and simulate orbits easily using 2 body physics (Kepler's Laws).
 
In the future, other helper modules like this one will be implemented so that specific bodies (physics's particles, solids, ...)
can also be represented and simulated easily.

### Features
- n-body dynamics for points systems
- Common forces and potential (gravity, ...)
- Kepler's orbits representation, interpolation and simulation
- Solver for n-body problems (RK4, Euler, ...)

There is no documentation for the framework yet but I made an example app that represents the solar system
so you can easily get started for now. Find it [there](https://github.com/samiBendou/nbodies)

## Why another physics engine ?
I wanted to design a physics engine that allows to play with fundamental physics. Compared to other physics engine
this one is not yet providing strong joint or collision mechanisms but allows more advanced dynamics.
For example you can represent a particle that follows an electromagnetic or gravitational field easily.
You can represent complex springs systems or whatever idea of experience you have.

Furthermore, I needed an API that will ease the implementation of [nbodies](https://github.com/samiBendou/nbodies) 
3D simulator.

## Usage
API is made to be very straightforward, it's built using [geomath](https://github.com/samiBendou/geomath) framework.
```rust
let position = geomath::Vector3::zeros(); // set initial position
let speed = geomath::Vector3::ones(); // set initial speed
let solver = dynamics::Solver::new(0.1, 1, dynamics::Method::RungeKutta4) // initialize solver with dt = 0.1 and 1 iteration per step
let point = dynamics::Point3::inertial(position, speed, 1.); // create a point of mass 1 kg with position and speed 
let mut cluster = dynamics::Cluster::new(vec![point]); // create a cluster containing the point

// Apply an acceleration to the points of the cluster
cluster.apply(solver, |point, points| {
    Vector3::unit_neg_z() // Constant downwards acceleration for each point
});
```
The cluster can then easily be used to create an animation with points using piston for example.

## Install
For now since geomath is not yet on crates.io you need to bind the two frameworks manually.
Clone [geomath](https://github.com/samiBendou/geomath)'s repo and modify your Cargo.toml to match
the path of geomath with the path of the cloned repo.