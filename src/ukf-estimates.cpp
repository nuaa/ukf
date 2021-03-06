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

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

#include "types.h"
#include "ukf.h"
#include "state.h"
#include "dynamics.h"
#include "debug.h"

/*
This implementation is split into two files, ukf.cpp and ukf-estimates.cpp,
in order to work around memory limitations in the TI CCS compiler.

See ukf.cpp for details.
*/

/*
Follows section 3.3.
*/
void UnscentedKalmanFilter::apriori_estimate(real_t dt, ControlVector c) {
    uint32_t i;

    /*
    Run the sigma points through the process model as per equation 37.
    No need to make a new sigma point distribution, we won't need the old one
    again.
    */
    for(i = 0; i < UKF_NUM_SIGMA; i++) {
        sigma_points.col(i) = integrator.integrate(
            State(sigma_points.col(i)), dt);

        /* Orientation often doesn't stay normalized after integration */
        Quaternionr q(sigma_points.col(i).segment<4>(9));
        q.normalize();
        sigma_points.col(i).segment<4>(9) << q.vec(), q.w();
    }

    /*
    This is where the dynamics model should be run, which is effectively the
    control input part of the process model.
    If there's no dynamics model selected, then the acceleration and angular
    acceleration should be zeroed as there's no way to predict them.
    */
    if(dynamics != NULL) {
        for(i = 0; i < UKF_NUM_SIGMA; i++) {
            AccelerationVector temp = dynamics->evaluate(
                sigma_points.col(i), c);
            sigma_points.col(i).segment<3>(6) = temp.segment<3>(0);
            sigma_points.col(i).segment<3>(16) = temp.segment<3>(3);
        }
    } else {
        /*
        FIXME: may not be necessary/helpful to clear out linear acceleration
        in the absence of a dynamics model.

        sigma_points.block<3, UKF_NUM_SIGMA>(6, 0) =
            Eigen::Matrix<real_t, 3, UKF_NUM_SIGMA>::Zero();
        */
        sigma_points.block<3, UKF_NUM_SIGMA>(16, 0) =
            Eigen::Matrix<real_t, 3, UKF_NUM_SIGMA>::Zero();
    }

    /*
    Now, calculate the mean as described in section 3.4, using the weights
    from the scaled unscented transform.
    For calculating a priori attitude estimate, we use the algorithm for
    computing the quaternion mean found in "Unscented Filtering in a Unit
    Quaternion Space for Spacecraft Attitude Estimation" by Yee-Jin Cheon.
    The following algorithm implements equation (41d) from that paper.
    */
    apriori_mean.segment<UKF_DIM+1>(0) =
        UKF_SIGMA_WMI*sigma_points.block<
            UKF_DIM+1, UKF_NUM_SIGMA-1>(0, 1).rowwise().sum() +
        UKF_SIGMA_WM0*sigma_points.col(0);

    /*
    Calculate the covariance as described in section 3.5.1. Note that there is
    another error with equation 63; we want to operate on the transformed
    sigma points, not the untransformed ones.
    For the attitude error, we use the same error vectors calculated using the
    gradient descent algorithm above.
    */
    w_prime.topRows<9>() = sigma_points.topRows<9>().colwise() -
        sigma_points.col(0).segment<9>(0);
    w_prime.bottomRows<12>() = sigma_points.bottomRows<12>().colwise() -
        sigma_points.col(0).segment<12>(13);

    /*
    The attitude part of this set of vectors is calculated using equation (45)
    from the paper mentioned above.
    */
    for(i = 0; i < UKF_NUM_SIGMA; i++) {
        Quaternionr err_q =
            (Quaternionr(sigma_points.col(i).segment<4>(9)) *
            Quaternionr(sigma_points.col(0).segment<4>(9)).conjugate());

        w_prime.block<3, 1>(9, i) = UKF_MRP_F *
            (err_q.vec() / (UKF_MRP_A + err_q.w()));
    }

    /* Calculate the covariance using equation 64. */
    apriori_covariance =
        UKF_SIGMA_WC0 * (w_prime.col(0) * w_prime.col(0).transpose());
    for(i = 1; i < UKF_NUM_SIGMA; i++) {
        apriori_covariance +=
            UKF_SIGMA_WCI * (w_prime.col(i) * w_prime.col(i).transpose());
    }
}

/*
This function implements equations 40, 41 and 42. Once again, there's an error
in equation 42; it's the a priori sigma points which should be transformed by
the measurement model, not the original sigma point distribution.
*/
void UnscentedKalmanFilter::measurement_estimate() {
    uint32_t i;
    Eigen::Matrix<
        real_t,
        Eigen::Dynamic,
        UKF_NUM_SIGMA,
        0,
        MEASUREMENT_MAX_DIM> measurement_sigma_points(
            sensor.size(), UKF_NUM_SIGMA);

    /*
    Run the a priori sigma points through the measurement model as per
    equation 40.
    */
    for(i = 0; i < UKF_NUM_SIGMA; i++) {
        measurement_sigma_points.col(i) = sensor.predict(sigma_points.col(i));
    }

    /* Now, calculate the mean of the measurement sigma point distribution. */
    Eigen::Matrix<real_t, UKF_NUM_SIGMA, 1> weights =
        Eigen::Matrix<real_t, UKF_NUM_SIGMA, 1>::Constant(UKF_SIGMA_WMI);
    weights[0] = UKF_SIGMA_WM0;
    measurement_estimate_mean = sensor.calculate_mean(
        measurement_sigma_points, weights);

    /*
    Calculating the measurement estimate covariance is also the same as for
    the a priori case, but without the quaternion stuff.
    */
    z_prime = sensor.calculate_deltas(
        measurement_sigma_points, measurement_estimate_mean);

    measurement_estimate_covariance =
        UKF_SIGMA_WC0 * (z_prime.col(0) * z_prime.col(0).transpose());
    for(i = 1; i < UKF_NUM_SIGMA; i++) {
        measurement_estimate_covariance +=
            UKF_SIGMA_WCI * (z_prime.col(i) * z_prime.col(i).transpose());
    }
}

/*
The final update step is calculated using section 3.6, specifically with
equations 74 and 75.
*/
void UnscentedKalmanFilter::aposteriori_estimate() {
    StateCovarianceVector update_temp;

    /*
    The update_temp vector will contain the attitude update as a vector, which
    we need to apply the to quaternion section of the state vector.
    */
    update_temp = kalman_gain * innovation;

    /*
    Update the attitude quaternion from the MRP vector, equation 45 from the
    Markley paper.
    */
    Vector3r d_p = update_temp.segment<3>(9);
    real_t x_2 = d_p.squaredNorm();
    real_t d_q_w = (-UKF_MRP_A * x_2 +
        UKF_MRP_F * std::sqrt(UKF_MRP_F_2 + (1.0 - UKF_MRP_A_2) * x_2)) /
        (UKF_MRP_F_2 + x_2);
    Vector3r d_q_xyz = ((1.0 / UKF_MRP_F) * (UKF_MRP_A + d_q_w)) * d_p;
    Quaternionr d_q;
    d_q.vec() = d_q_xyz;
    d_q.w() = d_q_w;

    Quaternionr update_q = d_q *
        Quaternionr(sigma_points.col(0).segment<4>(9));
    state.attitude() << update_q.vec(), update_q.w();

    AssertNormalized(Quaternionr(state.attitude()));

    /* The rest of the update step is very straightforward. */
    state.segment<9>(0) = apriori_mean.segment<9>(0) +
        update_temp.segment<9>(0);
    state.segment<12>(13) = apriori_mean.segment<12>(13) +
        update_temp.segment<12>(12);

    /*
    And the state covariance update from equation 75, no quaternion
    manipulation necessary.
    */
    state_covariance = apriori_covariance -
        (kalman_gain * innovation_covariance * kalman_gain.transpose());
}
