/* 
 * File:   extras.c
 * Author: stefano
 * 
 * Extras.(h|c) contains a few utilities that are not strictly related to the
 * core functionality.
 */

#include "extras.h"
#include "utils.h"
#include "stdlib.h"
#include "kernel.h"
#include <gsl/gsl_poly.h>
#include <gsl/gsl_statistics_double.h>
#ifndef JAVASCRIPT
#include "omp.h"
#else
#include "omp_shim.h"
#endif

double K_E_n(const gsl_matrix* alle, double n, double sigma) {
    double m_1 = MSUN_TO_MJUP(MGET(alle, 0, MASS));
    double m_2 = MGET(alle, 1, MASS);
    double m_3 = MGET(alle, 2, MASS);
    double m_12 = m_1 + m_2;

    double a_i = MGET(alle, 1, SMA);
    double a_o = MGET(alle, 2, SMA);

    double e_i = MGET(alle, 1, ECC);
    double e_o = MGET(alle, 2, ECC);

    double xi = acosh(1 / e_o) * sqrt(1 - e_o * e_o);
    double E_22 = 4. * SQRT_TWOPI / 3. * pow(1 - e_o*e_o, 0.75) / (e_o * e_o) * pow(sigma, 5. / 2.) * exp(-sigma * xi);
    double I_22 = 9. / 4. * m_3 / m_12 * POW_3(a_i / a_o) * E_22;

    double e_i_ind = sqrt(e_i * e_i + I_22 * I_22);
    double eps_o = sqrt(1 - e_i * e_i);

    double e_eq;
    if (e_i < 0.) {
        e_eq = (5. / 4.) * e_o * m_3 * (m_1 - m_2) * SQR(a_i / a_o) * sigma /
                (eps_o * fabs(m_1 * m_2 - m_12 * m_3 * (a_i / a_o) * eps_o * sigma));

    } else {

        double A = 0.75 * (m_3 / m_12) * POW_3(a_i / a_o) / POW_3(eps_o);
        double B = 15. / 64. * (m_3 / m_12) * (m_1 - m_2) / m_12 * POW_4(a_i / a_o) / POW_5(eps_o);
        double C = 3. / 4. * (m_1 * m_2 / SQR(m_12)) * SQR(a_i / a_o) / POW_4(eps_o);
        double D = 15. / 64. * (m_1 * m_2 / SQR(m_12)) * (m_1 - m_2) / m_12 * POW_3(a_i / a_o) * (1 + 4 * e_o * e_o) / (e_o * POW_6(eps_o));

        double a[9] = {-B*B, 2 * A*B, B * B + C * C - A*A,
            -2 * (A * B + 4 * C * D), A * A + 3 * C * C + 16 * D*D,
            -18 * C*D, 9. / 4. * C * C + 24 * D*D, -9 * C*D, 9 * D * D};

        gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(9);
        double z[16];
        gsl_poly_complex_solve(a, 9, w, z);

        for (int i = 0; i < 16; i += 2)
            if (z[i] > 0 && z[i] < 1) {
                e_eq = z[i];
                break;
            }

        gsl_poly_complex_workspace_free(w);
    }

    double alpha = fabs(1 - e_i / e_eq);
    double e_i_oct = (alpha < 1 ? (1 + alpha) * e_eq : e_i + 2 * e_eq);

    e_i = MAX(e_i_ind, e_i_oct);
    double s = -3 * e_i + 13. / 8. * POW_3(e_i) + 5. / 192. * POW_5(e_i) - 227. / 3072. * POW_7(e_i);
    double F = E_22 / (TWOPI * n);
    double M_i = m_3 / (m_12 + m_3);
    double M_o = (m_1 * m_2 / POW_2(m_12)) * pow(m_12 / (m_12 + m_3), 2. / 3.);
    double A_n = -9. / 2. * s * F * (M_i + M_o * pow(n, 2. / 3.));

    double E_n = 0.5 * POW_2(sigma - n) - 2. * A_n;
    /*
    PRINTNUM(A_n);
    PRINTNUM(E_n);
    PRINTNUM(I_22);
    PRINTNUM(n);
    PRINTNUM(sigma);
    PRINTNUM(s);
    PRINTNUM(F);
    PRINTNUM(e_i);
    PRINTNUM(e_i_oct);
    PRINTNUM(e_i_ind);
    PRINTNUM(e_eq);
     */
    return E_n;
}

/**
 * This function implements the stability criterion of Mardling, 2008a for coplanar,
 * three body systems. 
 * @param alle A matrix of _all_ elements in hierarchical coordinates, as 
 * returned by K_getAllElements_jacobi (or supplied by the user).
 * 
 * @return One of T_INAPPLICABLE, T_STABLE or T_UNSTABLE. The routine returns T_INAPPLICABLE
 * if one of the following is true: (1) the system contains more of three bodies (returns
 * T_STABLE if two-body system); (2) the mass conditions for which n:1 resonances dominate
 * (m_2/m_1 > 0.01 && m_2/m_1 > 0.01 or m_2/m1 > 0.05 || m_3/m_1 > 0.05). The routine returns 
 * T_UNSTABLE if the stability criterion (E_n < 0 and E_(n+1) < 0) is satisfied,
 * T_STABLE otherwise.
 */
int K_isMstable_coplanar(const gsl_matrix* alle) {

    if (MROWS(alle) > 3)
        return T_INAPPLICABLE;

    if (MROWS(alle) == 2)
        return T_STABLE;

    double m_1 = MSUN_TO_MJUP(MGET(alle, 0, MASS));
    double m_2 = MGET(alle, 1, MASS);
    double m_3 = MGET(alle, 2, MASS);

    // outside mass criterion
    if ((m_2 / m_1 < 0.01 || m_3 / m_1 < 0.01) && (m_2 / m_1 < 0.05 && m_3 / m_1 < 0.05))
        return T_INAPPLICABLE;

    // Elements of inner and outer binary
    double P_i = MGET(alle, 1, PER);
    double P_o = MGET(alle, 2, PER);
    double sigma = P_o / P_i;
    double n = floor(sigma);

    double E_n = K_E_n(alle, n, sigma);
    double E_n1 = K_E_n(alle, n + 1, sigma);

    if (E_n < 0 && E_n1 < 0)
        return T_UNSTABLE;
    else
        return T_STABLE;
}

double K_crossval_l1o(ok_kernel* k, int minalgo, int maxiter, double params[]) {
    int nd = K_getNdata(k);

    int np = omp_get_max_threads();
    ok_kernel * ks[np];
    double lh[np];
    double rms[nd];
    ok_progress prog = k->progress;

    for (int p = 0; p < np; p++) {
        ks[p] = K_clone(k);
        ks[p]->progress = NULL;
        lh[p] = 0.;
    }
    for (int i = 0; i < nd; i++) {
        rms[i] = 0;
    }

    bool invalid = false;

#pragma omp parallel for
    for (int i = 0; i < nd; i++) {
        if (invalid)
            continue;

        int p = omp_get_thread_num();
        gsl_matrix_memcpy(ks[p]->system->elements, k->system->elements);
        gsl_vector_memcpy(ks[p]->params, k->params);

        ks[p]->flags |= NEEDS_SETUP;
        K_calculate(ks[p]);

        double err = ks[p]->compiled[i][T_ERR];
        int set = (int) (ks[p]->compiled[i][T_SET]);


        ks[p]->compiled[i][T_ERR] = -1;

        K_minimize(ks[p], minalgo, maxiter, params);
        K_calculate(ks[p]);
        double n = K_getPar(ks[p], set + DATA_SETS_SIZE);
        double s = err * err + n * n;
        double diff = ks[p]->compiled[i][T_SVAL] - ks[p]->compiled[i][T_PRED];
        lh[p] += -0.5 * log(s) - 0.5 * diff * diff / s;
        ks[p]->compiled[i][T_ERR] = err;
        rms[i] = (diff * diff) / (err * err);

        if (prog != NULL && omp_get_thread_num() == 0) {
            char message[MAX_LINE];
            int ret = prog(i * np, nd, ks[p],
                    "K_crossVal_l1o");
            if (ret == PROGRESS_STOP) {
                invalid = true;
            }
        }
    }

    double lh2 = 0;
    for (int p = 0; p < np; p++) {
        lh2 += lh[p];
        K_free(ks[p]);
    }
    if (invalid)
        return INVALID_NUMBER;

    lh2 += -0.5 * (double) k->ndata * LOG_2PI;
    return -lh2;
}

double K_crossval_lno(ok_kernel* k, int no, int minalgo, int maxiter, double params[]) {
    int nd = K_getNdata(k);

    int np = omp_get_max_threads();
    ok_kernel * ks[np];
    double lh[np];
    double rms[nd];
    ok_progress prog = k->progress;

    for (int p = 0; p < np; p++) {
        ks[p] = K_clone(k);
        ks[p]->progress = NULL;
        lh[p] = 0.;
    }
    for (int i = 0; i < nd; i++) {
        rms[i] = 0;
    }

    bool invalid = false;



#pragma omp parallel for
    for (int i = 0; i < nd; i++) {
        if (invalid)
            continue;

        int p = omp_get_thread_num();
        gsl_matrix_memcpy(ks[p]->system->elements, k->system->elements);
        gsl_vector_memcpy(ks[p]->params, k->params);

        ks[p]->flags |= NEEDS_SETUP;
        K_calculate(ks[p]);

        double err = ks[p]->compiled[i][T_ERR];
        int set = (int) (ks[p]->compiled[i][T_SET]);


        ks[p]->compiled[i][T_ERR] = -1;

        K_minimize(ks[p], minalgo, maxiter, params);
        K_calculate(ks[p]);
        double n = K_getPar(ks[p], set + DATA_SETS_SIZE);
        double s = err * err + n * n;
        double diff = ks[p]->compiled[i][T_SVAL] - ks[p]->compiled[i][T_PRED];
        lh[p] += -0.5 * log(s) - 0.5 * diff * diff / s;
        ks[p]->compiled[i][T_ERR] = err;
        rms[i] = (diff * diff) / (err * err);

        if (prog != NULL && omp_get_thread_num() == 0) {
            char message[MAX_LINE];
            int ret = prog(i * np, nd, ks[p],
                    "K_crossVal_l1o");
            if (ret == PROGRESS_STOP) {
                invalid = true;
            }
        }
    }

    double lh2 = 0;
    for (int p = 0; p < np; p++) {
        lh2 += lh[p];
        K_free(ks[p]);
    }
    if (invalid)
        return INVALID_NUMBER;

    lh2 += -0.5 * (double) k->ndata * LOG_2PI;
    return -lh2;
}
