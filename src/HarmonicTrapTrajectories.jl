module HarmonicTrapTrajectories
export Trajectory, LinearRamp, MinimumJerk, LewisRiesenfeld, trajectory;

"""
    abstract type Trajectory

Transport trajectories in harmonic traps.

Subtypes will contain at least the following parameters:
* `t` - total motional time for the trajectory.
* `d` - total distance moved during the motion.

See also: [`LinearRamp`](@ref), [`MinimumJerk`](@ref), [`LewisRiesenfeld`](@ref)
"""
abstract type Trajectory end;

"""
    LinearRamp(t::Number, d::Number) <: Trajectory

Simple linear trajectory where particle moves at constant speed between start
and end points.

See also: [`Trajectory`](@ref), [`trajectory`](@ref)
"""
struct LinearRamp <: Trajectory
  t::Number;
  d::Number;
end

"""
    MinimumJerk(t::Number, d::Number) <: Trajectory

Trajectory minimising jerk throughout the motion.

See also: [`Trajectory`](@ref), [`trajectory`](@ref)
"""
struct MinimumJerk <: Trajectory
  t::Number;
  d::Number;
end

"""
    LewisRiesenfeld{T}(t::Number, d::Number, classical::T) where T <: Trajectory

Use Lewis-Riesenfeld invariants reverse engineering method to get optimum
trajectory based on the `classical` trajectory.

Refer to E. Torrontegui et al. Phys. Rev. A 83, 013415, 2011 (doi:
10.1103/PhysRevA.83.013415) for details.

See also: [`Trajectory`](@ref), [`trajectory`](@ref)
"""
struct LewisRiesenfeld{T<:Trajectory} <: Trajectory
  t::Number;
  d::Number;

  "Classical trajectory."
  classical::T;

  function LewisRiesenfeld{T}(t::Number, d::Number, classical::T) where T <: Trajectory
    if t !== classical.t || d !== classical.d
      throw(ArgumentError("classical trajectory must have the same `t` and `d` parameters."));
    end

    new(t, d, classical);
  end
end

"""
    LewisRiesenfeld{T}(t::Number, d::Number) where T<:Trajectory

Implicitly construct classical trajectory object from the given type.
"""
LewisRiesenfeld{T}(t::Number, d::Number) where T<:Trajectory = LewisRiesenfeld{T}(t, d, T(t, d));

"""
    trajectory(trajectory::T, t::Number, start::Number, params::Vararg) where T<:Trajectory

Return the location at time `t` for a particle on `trajectory` starting at a
position given by `start`.

Additional parameters required by the `fractional_trajectory()` corresponding to
the relevant `trajectory` type should be passed as additional parameters.


See also: [`fractional_trajectory`](@ref)
"""
function trajectory(trajectory::T, t::Number, start::Number, params::Vararg) where T<:Trajectory
  return start + trajectory.d * fractional_trajectory(trajectory, t, params...);
end

"""
    fractional_trajectory(trajectory::T, t::Number) where T <: Trajectory
    fractional_trajectory(trajectory::LinearRamp, t::Number)
    fractional_trajectory(trajectory::MinimumJerk, t::Number)

Positions for particles on `trajectory` at the time `t` as a fraction of the
total distance which should be travelled on the trajectory.
"""
fractional_trajectory(trajectory::LinearRamp, t::Number) = t/trajectory.t;

function fractional_trajectory(trajectory::MinimumJerk, t::Number)
  s = t / trajectory.t;
  return 10s^3 - 15s^4 + 6s^5;
end

"""
    fractional_trajectory(trajectory::LewisRiesenfeld{MinimumJerk}, t::Number, ω::Number)

Implementation for use of method of reverse engineered Lewis-Riesenfeld invariants where the
classical trajectory is the `MinimumJerk` trajectory and the harmonic trap frequency is `ω`.
"""
function fractional_trajectory(trajectory::LewisRiesenfeld{MinimumJerk}, t::Number, ω::Number)
  s = t / trajectory.t;
  return (60s - 180s^2 + 120s^3) / (trajectory.t * ω)^2 + fractional_trajectory(trajectory.classical, t);
end
end # module
