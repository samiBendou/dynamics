# dynamics
### A physics toolbox with a realistic physics engine

## Brief
dynamics is a framework that aims to provide common features needed when dealing with physics.

It mainly consists on a little realistic physics engine that applies the 2nd law of Newton on a set of body,
You specify the forces for each body, and the hole gets updated.

The framework also provides an orbital module that is made to manipulate and simulate orbits easily. 
In the future, other helper modules like this one will be implemented so that specific bodies (physics's particles, solids, ...)
can also be represented and simulated easily.

### Features
- Point systems dynamics
- Forces and potentials computation
- Kepler's orbits representation, interpolation, serialization and simulation
- Generic solver for points systems (RK4, Euler, ...)

There is no documentation for the framework yet but I made an example app that represents the solar system
so you can easily get started for now. Find it [there](https://github.com/samiBendou/nbodies)

### Why another physics engine ?
I wanted to design a physics engine that allows to play with fundamental physics. Compared to other physics engine
this one is not yet providing strong joint or collision mechanisms but allows more advanced movement definition.
For example you can represent a particle that follows an electromagnetic or gravitational field easily.
You can represent complex springs systems or whatever idea you have

### Usage
API is made to be very straightforward, it's built using [geomath](https://github.com/samiBendou/nbodies) framework.
For now since geomath is not yet on crates.io you need to bind the two frameworks manually. You can do so just by cloning
geomath and putting it's root directory at the place dynamics's root directory is located.
```rust
let position = geomath::Vector3::zeros(); // set initial position
let speed = geomath::Vector3::ones(); // set initial speed
let solver = dynamics::Solver::new(0.1, 1, dynamics::Method::RungeKutta4) // initialize solver with dt = 0.1 and 1 iteration per step
let point = dynamics::Point3::inertial(position, speed, 1.); // create a point of mass 1 kg with position and speed 
let mut cluster = dynamics::Cluster::new(vec![point]); // create a cluster containing the point

// Apply an acceleration to the points of the cluster
cluster.apply(solver, |point, points| {
    Vector3::unit_neg_z() // Constant downwards acceleration
});
```
The cluster can then easily be used to create an animation with points using piston for example. 