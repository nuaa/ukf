/*
Copyright (C) 2013 Daniel Dyer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <stdint.h>
#include <cmath>


#include "types.h"
#include "state.h"

#define CONTROL_MAX_DIM 4
/* Disable dynamics model if velocity is less than 1m/s or greater than 100m/s */
#define DYNAMICS_MIN_V 1.0
#define DYNAMICS_MAX_V 100.0
/* Disable dynamics model if angular velocity exceeds ~120deg/s */
#define DYNAMICS_MAX_RATE (M_PI*0.63)
/* Airframe minimums */
#define AIRFRAME_MIN_MASS 0.1
#define AIRFRAME_MIN_MOMENT 1e-6

/*
Typedef for control vector. It is dynamic in order to allow different dynamics
models to accept different numbers of control inputs.
*/
typedef Eigen::Matrix<
    real_t,
    Eigen::Dynamic,
    1,
    0,
    CONTROL_MAX_DIM> ControlVector;

typedef Eigen::Matrix<real_t, 6, 1> AccelerationVector;

/*
Dynamics model base class. The public interface is via the evaluate() method,
which returns a 6-dimensional column vector containing linear acceleration
and angular acceleration. In the future, this will need to accept a control
input vector as a parameter.
*/
class DynamicsModel {
public:
    virtual ~DynamicsModel();
    virtual AccelerationVector evaluate(
    const State &in, const ControlVector &control) const = 0;
};

/*
Just a basic dynamics model which predicts only acceleration using a
centripetal acceleration model.
*/
class CentripetalModel: public DynamicsModel {
public:
    AccelerationVector evaluate(
    const State &in, const ControlVector &control) const;
};

/*
Fixed-wing flight dynamics model using coefficient build-up with quartic
*/
class FixedWingFlightDynamicsModel: public DynamicsModel {
    /*
    Prop area and prop Cve are used in calculating thrust based on RPM. The
    Cve relates control value to the exit velocity, and the prop area
    determines the ratio of thrust to velocity differential.
    */
    real_t prop_area, prop_cve;

    /* Use reciprocal of mass for performance */
    real_t mass_inv;

    /*
    Coefficient vector for side force:
        alpha^2 * qbar,
        alpha * qbar,
        beta^2 * qbar,
        beta * qbar,
        alpha^2 * beta * qbar,
        alpha * beta * qbar,
        yaw rate * qbar, roll rate * qbar,
        [control 0-3] * qbar
    */
    Eigen::Matrix<real_t, 8 + CONTROL_MAX_DIM, 1> c_side_force;

    /*
    Coefficient vectors for aerodynamic moments:
    CM (pitch) -
        alpha * qbar
        Q^2 * sign(Q) * qbar
        control[0-3] * qbar
    */
    Eigen::Matrix<real_t, 2 + CONTROL_MAX_DIM, 1> c_pitch_moment;

    /*
    CN (yaw) -
        beta * qbar
        R * qbar
        control[0-3] * qbar
    */
    Eigen::Matrix<real_t, 2 + CONTROL_MAX_DIM, 1> c_yaw_moment;

    /*
    CL (roll) -
        P * qbar
        control[0-3] * qbar
    */
    Eigen::Matrix<real_t, 1 + CONTROL_MAX_DIM, 1> c_roll_moment;

    /* Store inverse inertia tensor for performance */
    Eigen::Matrix<real_t, 3, 3> inertia_tensor_inv;

    /*
    Polynomial coefficients for varying alpha:
    C = v[0]x^4 + v[1]x^3 + v[2]x^2 + v[3]x^1 + v[4]
    coefficient of drag (alpha)
    coefficient of lift (alpha)
    */
    Eigen::Matrix<real_t, 5, 1> c_drag_alpha, c_lift_alpha;

    /* Motor force direction vector in body frame */
    Vector3r motor_thrust;

    /* Motor control channel index */
    ControlVector::Index motor_idx;

public:
    FixedWingFlightDynamicsModel(void) {
        c_drag_alpha.setZero();
        c_lift_alpha.setZero();
        c_side_force.setZero();
        c_yaw_moment.setZero();
        c_pitch_moment.setZero();
        c_roll_moment.setZero();

        prop_area = 0.0;
        prop_cve = 0.0;

        mass_inv = 0.0;
        inertia_tensor_inv.setZero();

        motor_idx = 0;
        motor_thrust << -1, 0, 0;
    }

    AccelerationVector evaluate(
    const State &in, const ControlVector &control) const;

    /* Airframe property and coefficient setters */
    void set_mass(real_t in_mass) {
        assert(std::abs(in_mass) > AIRFRAME_MIN_MASS);
        mass_inv = 1.0 / in_mass;
    }

    void set_inertia_tensor(
    const Eigen::Matrix<real_t, 3, 3>& in_inertia_tensor) {
        inertia_tensor_inv = in_inertia_tensor.inverse();
    }

    void set_prop_coeffs(real_t in_prop_area, real_t in_prop_cve) {
        prop_area = in_prop_area;
        prop_cve = in_prop_cve;
    }

    void set_drag_coeffs(const Eigen::Matrix<real_t, 5, 1>& in_coeffs) {
        c_drag_alpha << in_coeffs;
    }

    void set_lift_coeffs(const Eigen::Matrix<real_t, 5, 1>& in_coeffs) {
        c_lift_alpha << in_coeffs;
    }

    void set_side_coeffs(
    const Eigen::Matrix<real_t, 8 + CONTROL_MAX_DIM, 1>& in_coeffs) {
        c_side_force << in_coeffs;
    }

    void set_pitch_moment_coeffs(
    const Eigen::Matrix<real_t, 2 + CONTROL_MAX_DIM, 1>& in_coeffs) {
        c_pitch_moment << in_coeffs;
    }

    void set_roll_moment_coeffs(
    const Eigen::Matrix<real_t, 1 + CONTROL_MAX_DIM, 1>& in_coeffs) {
        c_roll_moment << in_coeffs;
    }

    void set_yaw_moment_coeffs(
    const Eigen::Matrix<real_t, 2 + CONTROL_MAX_DIM, 1>& in_coeffs) {
        c_yaw_moment << in_coeffs;
    }

    void set_motor_index(ControlVector::Index idx) {
        motor_idx = idx;
    }
};

#endif
