#include "periodogram.h"
#include "stdint.h"
#include "utils.h"
#include "math.h"
#ifndef JAVASCRIPT
#include "omp.h"
#else
#include "omp_shim.h"
#endif

#include "kernel.h"
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_sort.h"
#include "gsl/gsl_roots.h"

#define BUF_SDF 0
#define BUF_CDF 1
#define BUF_C 2
#define BUF_S 3
#define BUF_SIG 4

#define SDF(j) (MGET(buf, j, BUF_SDF))
#define CDF(j) (MGET(buf, j, BUF_CDF))
#define S(j) (MGET(buf, j, BUF_S))
#define C(j) (MGET(buf, j, BUF_C))
#define SIG(j) (MGET(buf, j, BUF_SIG))

double _baluev_tau(double z_1, void* params) {
    double* p = (double*) params;
    double ndata = p[0];
    double W = p[1];
    double target = p[2];
    double fap_single = pow(1. - 2. * z_1 / (double) ndata, 0.5 * (double) (ndata - 2.));
    return W * fap_single * sqrt(z_1) - target;
}

double _find_z(gsl_root_fsolver* s, gsl_function* F, double target, double x_lo, double x_hi) {
    ((double*) F->params)[2] = target;
    gsl_root_fsolver_set(s, F, x_lo, x_hi);


    int status = GSL_SUCCESS;
    do {
        if (gsl_root_fsolver_iterate(s) != GSL_SUCCESS)
            return INVALID_NUMBER;
        x_lo = gsl_root_fsolver_x_lower(s);
        x_hi = gsl_root_fsolver_x_upper(s);

        status = gsl_root_test_interval(x_lo, x_hi, 1e-5, 1e-3);
    } while (status == GSL_CONTINUE);
    return 0.5 * (x_lo + x_hi);
}

/**
 * Computes the Lomb-Scargle periodogram of the matrix "data". "data" should contain at least three
 * columns: time, measurement and measurement error. The periodogram is calculated in "samples" intervals
 * between "Pmin" and "Pmax", spaced logarithmically. 
 * 
 * The function returns a matrix of "samples" rows and several columns, including period, power (z) and 
 * an estimation of the upper bound for the false alarm probability. The estimation is calculated using 
 * the method of Baluev, 2008 (Baluev08). The column PS_Z_LS contains the unnormalized LS periodogram 
 * (z = 1/2 * (Chi^2_0 - Chi^2_SC)), while the column PS_Z contains z_1 = 1/2 * N_H * z / Chi^2_0 (z_1 in Baluev08). 
 * The FAP upper bound is estimated as ~ tau(z_1). (Another estimate of the FAP can be calculated by 
 * estimating the indep. frequencies through your own algorithm, or using the ok_periodogram_boot routine.)
 * 
 * @param data Input data containing the data; each row containing (t_i, x_i, sigma_i)
 * @param samples Number of frequencies sampled
 * @param Pmin Minimum period sampled
 * @param Pmax Maximum period sampled
 * @param method Method to compute periodogram (ignored)
 * @param timecol Time column (e.g. 0) in the matrix data
 * @param valcol Value column (e.g. 1) in the matrix data
 * @param sigmacol Sigma column (e.g. 2) in the matrix data
 * @param p If not NULL, it is used to return additional info for the periodogram and reuse matrices to save space/speed. If you pass
 * a value different than NULL, you are responsible for deallocating the workspace and its fields. p->buf is an array of
 * gsl_matrix*, sized the same as the value of omp_get_max_threads().
 * @return A matrix containing: {PS_TIME, PS_Z, PS_FAP, PS_Z_LS} (period, power, FAP upper limit, unnormalized
 * LS power). You are responsible for deallocating it.
 */
gsl_matrix* ok_periodogram_ls(const gsl_matrix* data, const unsigned int samples, const double Pmin, const double Pmax, const int method,
                              unsigned int timecol, unsigned int valcol, unsigned int sigcol, ok_periodogram_workspace* p) {

    gsl_matrix* ret = NULL;
    gsl_matrix* buf = NULL;
    gsl_vector* bufv = gsl_vector_alloc(data->size1);

    int ndata = data->size1;

    // If no pre-allocated buffers are passed through p, or p is null,
    // allocate those buffers.
    if (p != NULL) {
        if (p->per != NULL && MROWS(p->per) == samples && MCOLS(p->per) == PS_SIZE)
            ret = p->per;
        if (p->buf != NULL && MROWS(p->buf) == ndata && MCOLS(p->per) == 5)
            ret = p->buf;
    }

    ret = (ret != NULL ? ret : gsl_matrix_alloc(samples, PS_SIZE));
    buf = (buf != NULL ? buf : gsl_matrix_alloc(ndata, 5));

    double fmin = 1. / Pmax;
    double fmax = 1. / Pmin;
    double df = (fmax - fmin) / (double) samples;


    gsl_matrix_get_col(bufv, data, timecol);
    double W = 2. * M_PI * gsl_stats_sd(bufv->data, 1, ndata) / Pmin;
    gsl_matrix_get_col(bufv, data, valcol);
    double avg = gsl_stats_mean(bufv->data, 1, ndata);
    double z1_max = 0.;
    double xa[ndata];

    // pre-calculate cdf, sdf
    for (int i = 0; i < ndata; i++) {
        double t = MGET(data, i, timecol) - MGET(data, 0, timecol);
        MSET(buf, i, BUF_CDF, cos(2 * M_PI * df * t));
        MSET(buf, i, BUF_SDF, sin(2 * M_PI * df * t));
        MSET(buf, i, BUF_C, cos(2 * M_PI * fmin * t));
        MSET(buf, i, BUF_S, sin(2 * M_PI * fmin * t));
        MSET(buf, i, BUF_SIG, 1. / (MGET(data, i, sigcol) * MGET(data, i, sigcol)));
        xa[i] = MGET(data, i, valcol) - avg;
    }

    // Calculate periodogram by looping over all angular frequencies
    for (int i = 0; i < samples; i++) {
        // Current frequency
        double f = fmin + df * i;


        double w = 2 * M_PI*f;

        // Calculate tau(w)
        double s_2wt = 0.;
        double c_2wt = 0.;

        for (int j = 0; j < ndata; j++) {
            double cos_wt = C(j);
            double sin_wt = S(j);
            c_2wt += (1. - 2. * sin_wt * sin_wt) * SIG(j);
            s_2wt += (2. * sin_wt * cos_wt) * SIG(j);
        }

        double tau = atan2(s_2wt, c_2wt) / (2. * w);
        double numa = 0.;
        double numb = 0.;
        double dena = 0.;
        double denb = 0.;
        double numa_w = 0.;
        double numb_w = 0.;
        double dena_w = 0.;
        double denb_w = 0.;

        double coswtau = cos(w * tau);
        double sinwtau = sin(w * tau);
        double chi2_h = 0.;
        double chi2_h_w = 0;

        for (int j = 0; j < ndata; j++) {

            double sig = SIG(j);

            const double cos_wt = C(j);
            const double sin_wt = S(j);

            double cos_wdf = CDF(j);
            double sin_wdf = SDF(j);

            double c = cos_wt * coswtau + sin_wt * sinwtau;
            double s = sin_wt * coswtau - cos_wt * sinwtau;
            double x = xa[j];

            MSET(buf, j, BUF_C, cos_wt * cos_wdf - sin_wt * sin_wdf);
            MSET(buf, j, BUF_S, sin_wt * cos_wdf + cos_wt * sin_wdf);

            numa += x * c * sig;
            numb += x * s * sig;
            dena += c * c * sig;
            denb += s * s * sig;
            chi2_h += x * x * sig;

            numa_w += c;
            numb_w += s;
            dena_w += c*c;
            denb_w += s*s;

            chi2_h_w += 1;
        }


        double z = 0.5 * (numa * numa / dena + numb * numb / denb);
        double z_1 = z * ndata / chi2_h;

        double w_1 = 0.5 * (numa_w * numa_w / dena_w + numb_w * numb_w / denb_w) * ndata / chi2_h_w;

        double fap_single = pow(1. - 2. * z_1 / (double) ndata, 0.5 * (double) (ndata - 3.));
        double tau_z = W * fap_single * sqrt(z_1);

        MSET(ret, samples - i - 1, PS_TIME, 1. / f);
        MSET(ret, samples - i - 1, PS_Z, z_1);
        MSET(ret, samples - i - 1, PS_Z_LS, z);
        MSET(ret, samples - i - 1, PS_FAP, MIN(fap_single + tau_z, 1.));
        MSET(ret, samples - i - 1, PS_TAU, tau);
        MSET(ret, samples - i - 1, PS_WIN, w_1);

        z1_max = MAX(z1_max, z_1);
    }

    if (p != NULL && p->calc_z_fap) {
        gsl_root_fsolver * s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
        double pars[3];
        pars[0] = ndata;
        pars[1] = W;
        pars[2] = 0.;

        gsl_function F;
        F.function = _baluev_tau;
        F.params = pars;

        double zz = z1_max;
        while (_baluev_tau(zz, pars) > 1e-3)
            zz *= 2;

        p->z_fap_3 = _find_z(s, &F, 1e-3, 0.1, zz);
        p->z_fap_2 = _find_z(s, &F, 1e-2, 0.1, p->z_fap_3);
        p->z_fap_1 = _find_z(s, &F, 1e-1, 0.1, p->z_fap_2);


        gsl_root_fsolver_free(s);
        p->calc_z_fap = false;
    }

    if (p == NULL) {
        gsl_matrix_free(buf);
    } else {
        p->per = ret;
        p->buf = buf;
        p->zmax = z1_max;
    };

    gsl_vector_free(bufv);

    return ret;
}

double _kminimize(ok_kernel* k, int algo) {
    K_calculate(k);

    double chi2 = K_getChi2_nr(k);
    while (true) {
        K_minimize(k, algo, 10000, NULL);

        if (K_getChi2_nr(k) - chi2 < -0.01)
            chi2 = K_getChi2_nr(k);
        else
            break;
    }

    return K_getChi2_nr(k);
}

gsl_matrix* ok_periodogram_full(ok_kernel* k, int type, int algo, bool circular, unsigned int sample,
                                const unsigned int samples, const double Pmin, const double Pmax) {

    k = K_clone(k);
    K_calculate(k);

    // Input data for LS periodogram
    gsl_matrix* data = ok_buf_to_matrix(K_compileData(k), K_getNdata(k), DATA_SIZE);


    if (type == PS_TYPE_RESIDUALS) {
        // If residuals periodogram, subtract signal from data
        for (int i = 0; i < data->size1; i++)
            MSET(data, i, T_SVAL, MGET(data, i, T_SVAL) - MGET(data, i, T_PRED));
    } else if (type == PS_TYPE_DATA) {
        // If full periodogram, then start with no planets
        K_removePlanet(k, -1);
    }

    // Calculate LS periodogram
    gsl_matrix* ret = ok_periodogram_ls(data, samples, Pmin, Pmax, 0, T_TIME, T_SVAL, T_ERR, NULL);
    int np = K_getNplanets(k) + 1;

    // Number of minimizable offsets
    int no = 0;
    for (int i = 0; i < DATA_SETS_SIZE; i++)
        if (VIGET(k->parFlags, i) & MINIMIZE) {
            no++;
        }

    // Calculate baseline chi^2 (Chi^2_H)

    double Chi2_H = _kminimize(k, algo);

    // Normalizing factor for power
    double nd = 0.5 * (K_getNdata(k) - no);



    #pragma omp parallel for
    for (int r = 0; r < samples; r++) {
        double P = MGET(ret, r, PS_TIME);
        double K = sqrt(MGET(ret, r, PS_Z));

        ok_kernel* k2 = K_clone(k);
        K_calculate(k2);

        double args[] = {PER, P, DONE};
        K_addPlanet(k2, args);
        K_setElement(k2, np, SEMIAMP, K);

        K_setElementFlag(k2, np, PER, ACTIVE);

        if (circular) {
            K_setElementFlag(k2, np, ECC, ACTIVE);
            K_setElementFlag(k2, np, LOP, ACTIVE);
        }

        double Chi2_K = _kminimize(k2, algo);

        double z = nd * (Chi2_H - Chi2_K) / Chi2_H;
        MSET(ret, r, PS_Z, z);
        fflush(stdout);
    }

    return ret;

}

/**
 * Estimates the FAP by Monte Carlo bootstrapping of the original data. "trials" bootstrapped data sets are
 * generated using the random number generator "rng"; for each data set, the routine ok_periodogram_ls estimates
 * z_max, which is collected into an array and returned into "zmax". Bootstrapped datasets are built by selecting
 * with replacement from the input dataset, keeping times of observation fixed. 
 * @param data Input matrix containing the data; each row containing (t_i, x_i, sigma_i)
 * @param trials Number of bootstrap trials
 * @param samples Number of frequencies sampled
 * @param Pmin Minimum period sampled
 * @param Pmax Maximum period sampled
 * @param method Method used to compute periodogram (ignored)
 * @param timecol Time column (e.g. 0) in the matrix data
 * @param valcol Value column (e.g. 1) in the matrix data
 * @param sigmacol Sigma column (e.g. 2) in the matrix data
 * @param rng A pre-allocated random number generator
 * @param p If specified, returns additional info for the periodogram and reuses matrices to save space/speed. If you pass
 * a value different than NULL, you are responsible for deallocating the workspace and its fields. p->zm returns a sorted
 * vector of the maximum powers in each synthetic trial.
 * @param prog An ok_progress* callback; if different from NULL, can be used to stop or report progress.
 * @return A matrix containing: {PS_TIME, PS_Z, PS_FAP, PS_Z_LS} (period, power, bootstrapped FAP, unnormalized
 * LS power). You are responsible for deallocating it.

 */
gsl_matrix* ok_periodogram_boot(const gsl_matrix* data, const unsigned int trials, const unsigned int samples,
                                const double Pmin, const double Pmax, const int method,
                                const unsigned int timecol, const unsigned int valcol, const unsigned int sigcol,
                                const unsigned long int seed, ok_periodogram_workspace* p, ok_progress prog) {


    int nthreads = omp_get_max_threads();

    ok_periodogram_workspace * w[nthreads];
    gsl_matrix * mock[nthreads];
    gsl_rng * rng[nthreads];

    rng[0] = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng[0], seed);

    for (int i = 0; i < nthreads; i++) {
        w[i] = (ok_periodogram_workspace*) malloc(sizeof (ok_periodogram_workspace));
        w[i]->per = NULL;
        w[i]->buf = NULL;
        w[i]->calc_z_fap = false;
        mock[i] = ok_matrix_copy(data);
        if (i > 0) {
            rng[i] = gsl_rng_alloc(gsl_rng_default);
            gsl_rng_set(rng[i], seed + i);
        }
    }

    gsl_matrix* ret = ok_matrix_copy(ok_periodogram_ls(data, samples, Pmin, Pmax, method, timecol, valcol, sigcol, NULL));

    gsl_vector* zmax = (p != NULL && p->zm != NULL ? p->zm : gsl_vector_alloc(trials));

    bool abort = false;
    #pragma omp parallel for
    for (int i = 0; i < trials; i++) {
        if (!abort) {
            int nt = omp_get_thread_num();

            ok_bootstrap_matrix_mean(data, T_TIME, T_VAL, mock[nt], rng[nt]);
            ok_periodogram_ls(mock[nt], samples, Pmin, Pmax, method, timecol, valcol, sigcol, w[nt]);
            zmax->data[i] = w[nt]->zmax;

            if (nt == 0 && prog != NULL) {
                int ret = prog(i * nthreads, trials, NULL,
                               "ok_periodogram_boot");
                if (ret == PROGRESS_STOP) {
                    abort = true;
                    #pragma omp flush (abort)
                }
            }
        }
    }

    gsl_sort(zmax->data, 1, trials);

    for (int i = 0; i < ret->size1; i++) {
        if (MGET(ret, i, PS_Z) > zmax->data[trials - 1])
            MSET(ret, i, PS_FAP, 1. / (double) trials);
        else if (MGET(ret, i, PS_Z) < zmax->data[0])
            MSET(ret, i, PS_FAP, 1.);
        else {
            int idx = ok_bsearch(zmax->data, MGET(ret, i, PS_Z), trials);
            MSET(ret, i, PS_FAP, (double) (trials - idx) / (double) trials);
        }
    }



    for (int i = 0; i < nthreads; i++) {
        gsl_matrix_free(w[i]->buf);
        gsl_matrix_free(w[i]->per);
        free(w[i]);
        gsl_matrix_free(mock[i]);
        gsl_rng_free(rng[i]);
    }

    if (p != NULL) {
        if (p->zm != NULL)
            gsl_vector_free(zmax);
    }
    return ret;
}
