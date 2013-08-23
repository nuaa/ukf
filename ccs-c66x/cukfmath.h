#ifndef CUKFMATH_H_
#define CUKFMATH_H_

#define X 0
#define Y 1
#define Z 2
#define W 3

#ifdef UKF_USE_DSP_INTRINSICS
#undef cos
#define cos Cosdp
#undef sin
#define sin Sindp
#undef atan2
#define atan2 Atan2dp
#undef sqrt
#define sqrt Sqrtdp
#define sqrt_inv(x) Rsqrtdp(x)
#else
#define sqrt_inv(x) (1.0 / sqrt((x)))
#endif

static inline void _mul_quat_vec3(real_t res[3], real_t q[4], real_t v[3]) {
    /*
    Multiply a quaternion by a vector (i.e. transform a vectory by a
    quaternion)

    v' = q * v * conjugate(q), or:
    t = 2 * cross(q.xyz, v)
    v' = v + q.w * t + cross(q.xyz, t)

    http://molecularmusings.wordpress.com/2013/05/24/a-faster-quaternion-vector-multiplication/
    */

    assert(res);
    assert(q);
    assert(v);
    assert(res != v);

    real_t t[3];
    t[X] = 2.0 * (q[Y]*v[Z] - q[Z]*v[Y]);
    t[Y] = 2.0 * (q[Z]*v[X] - q[X]*v[Z]);
    t[Z] = 2.0 * (q[X]*v[Y] - q[Y]*v[X]);

    res[X] = v[X] + q[W]*t[X] + (q[Y]*t[Z] - q[Z]*t[Y]);
    res[Y] = v[Y] + q[W]*t[Y] + (q[Z]*t[X] - q[X]*t[Z]);
    res[Z] = v[Z] + q[W]*t[Z] + (q[X]*t[Y] - q[Y]*t[X]);
}

static inline void _mul_vec_scalar_add_vec(real_t v1[], real_t a, real_t v2[],
        size_t len) {
    assert(v1);
    assert(v2);
    for (size_t i = 0; i < len; i++) {
        v1[i] = v1[i] * a + v2[i];
    }
}

static inline void _add_vec_vec(real_t v1[], real_t v2[], size_t len) {
    assert(v1);
    assert(v2);
    for (size_t i = 0; i < len; i++) {
        v1[i] += v2[i];
    }
}

static inline void _inv_mat3x3(real_t out[9], real_t M[9]) {
    /*
    M = 0 1 2
        3 4 5
        6 7 8
    */
    assert(M);
    assert(out);
    real_t det = M[0] * (M[8]*M[4] - M[5]*M[7]) -
                 M[1] * (M[8]*M[3] - M[5]*M[6]) +
                 M[2] * (M[7]*M[3] - M[4]*M[6]);

    assert(fabs(det) > 1e-6);
    det = 1.0 / det;

    out[0] =  (M[8]*M[4] - M[5]*M[7]) * det;
    out[3] = -(M[8]*M[3] - M[5]*M[6]) * det;
    out[6] =  (M[7]*M[3] - M[4]*M[6]) * det;
    out[1] = -(M[8]*M[1] - M[2]*M[7]) * det;
    out[4] =  (M[8]*M[0] - M[2]*M[6]) * det;
    out[7] = -(M[7]*M[0] - M[1]*M[6]) * det;
    out[2] =  (M[5]*M[1] - M[2]*M[4]) * det;
    out[5] = -(M[5]*M[0] - M[2]*M[3]) * det;
    out[8] =  (M[4]*M[0] - M[1]*M[3]) * det;
}

#endif
