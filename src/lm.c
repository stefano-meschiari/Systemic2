#include <gsl/gsl_multifit_nlin.h>

#include "math.h"
#ifndef JAVASCRIPT
#include "omp.h"
#else
#include "omp_shim.h"
#endif

#include "utils.h"

#include "kernel.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_int.h>

static inline void K_lm_calc(ok_kernel* k, double* f) {
    k->flags |= NEEDS_SETUP;
    K_calculate(k);
    double** compiled = k->compiled;


    for (int i = 0; i < k->ndata; i++) {
        const double* comprow = compiled[i];
        int set = (int) comprow[T_SET];
        double n = K_getPar(k, set + DATA_SETS_SIZE);

        double diff = (comprow[T_PRED] - comprow[T_SVAL]);
        double s = comprow[T_ERR];

        f[i] = diff / sqrt(s * s + n * n);
        f[i + k->ndata] = sqrt(fabs(log(n * n + s * s))) * SIGN(log(n * n + s * s));
    }


}

int K_lm_f(const gsl_vector* x, void* params, gsl_vector* f) {
    ok_kernel_minimizer_pars* sp = (ok_kernel_minimizer_pars*) params;

    ok_kernel* k = sp->kernel;

    double** pars = sp->pars;

    for (int i = 0; i < x->size; i++)
        *(pars[i]) = x->data[i] * sp->steps[i];

    K_lm_calc(k, f->data);
    printf("%e\n", K_getLoglik(k));
    for (int i = 0; i < x->size; i++)
        printf("%e ", i, x->data[i]);
    printf("\n");
    if (IS_NOT_FINITE(K_getLoglik(k))) {


        exit(0);
    }

    return GSL_SUCCESS;
}

int K_minimize_lm(ok_kernel* k, int maxiter, double params[]) {
    K_calculate(k);

    double prev_chi2 = K_getLoglik(k);
    int verbose = 0;

    if (params != NULL) {
        int i = 0;
        while (true) {
            if (params[i] == DONE)
                break;
            if (params[i] == OPT_VERBOSE_DIAGS)
                verbose = (int) params[i + 1];
            i += 2;
        }
    }



    ok_kernel_minimizer_pars sp = K_getMinimizedVariables(k);
    int npars = sp.npars;
    int nx = 2 * k->ndata;

    double L = K_getLoglik(k);
    double epsilon = sqrt(GSL_DBL_EPSILON);
    double target_DL = 10;

    /*
    for (int i = 0; i < npars; i++) {
        for (int j = 0; j < 10; j++) {
            double v = *(sp.pars[i]);

     *(sp.pars[i]) = v + sp.steps[i];

            k->flags |= NEEDS_SETUP;
            K_calculate(k);
            double dL = fabs(K_getLoglik(k) - L);

     *(sp.pars[i]) = v;
            sp.steps[i] = 0.5 * sp.steps[i] * MAX(MIN(1 + target_DL / dL,
                    10), 0.1);
        }
    }
     */
    for (int i = 0; i < npars; i++)
        if (sp.type[i] == MA || sp.type == INC || sp.type == LOP || sp.type == NODE)
            sp.steps[i] = 180;
        else if (sp.type[i] == ECC)
            sp.steps[i] = 0.1;
        else if (sp.type[i] == -1)
            sp.steps[i] = 1.;
        else
            sp.steps[i] = *(sp.pars[i]);

    gsl_multifit_fdfsolver * s
            = gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmder, nx, npars);

    gsl_multifit_function_fdf fdf;
    fdf.f = &K_lm_f;
    fdf.df = NULL;
    fdf.fdf = NULL;
    fdf.n = nx;
    fdf.p = npars;
    fdf.params = &sp;

    gsl_vector* x = gsl_vector_alloc(npars);
    for (int i = 0; i < npars; i++) {
        VSET(x, i, *(sp.pars[i]) / sp.steps[i]);
        printf("%e %e\n", *(sp.pars[i]), x->data[i]);
    }

    gsl_multifit_fdfsolver_set(s, &fdf, x);

    double loglik = ok_vector_sum_2(s->f);


    int status = 0;
    int user_status = PROGRESS_CONTINUE;

    for (int i = 0; i < maxiter; i++) {
        int status = gsl_multifit_fdfsolver_iterate(s);
        double loglik = ok_vector_sum_2(s->f);
        printf("%d %e %s\n", i, loglik, gsl_strerror(status));
    }


    k->flags |= NEEDS_SETUP;
    K_calculate(k);

    return user_status;
}
