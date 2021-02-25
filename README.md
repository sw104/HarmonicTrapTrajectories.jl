# HarmonicTrapTrajectories

This simple package implements some transport trajectories which may be used to
move atoms using harmonic traps in ultracold atom experiments.

It is a very simple package and was made for use with
[QuantumSimulation.jl](https://github.com/sw104/QuantumSimulation.jl), but is
quite independent.

The trajectories included are currently:

* Linear ramp.
* Minimum jerk trajectory.
* Lewis-Riesenfeld invariants reverse engineered trajectory with the minimum
  jerk trajectory as the classical trajectory.

This package is still in development and is liable to change significantly
without notice. Use at your own risk.
