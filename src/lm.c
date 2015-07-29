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

#define PARTYPE_PAR -1
#define STEP_EPS (sqrt(1e-20))

typedef struct {
    ok_kernel* k;
    double** pars;
    double* best;
    double* stepscale;
    double** compiled;
    double* f0;
    double* f1;
    double* f2;
    double* f3;

    int* parstype;
    unsigned int ndata;
    unsigned int iterations;
    unsigned int every;
    unsigned int maxiterations;
    unsigned int status;
    unsigned int npars;
    double min_chi;
    double st;

    bool high_df;
} ok_lm_params;

static inline void K_lm_calc(ok_kernel* k, double* f, const int ndata, const double** compiled, ok_lm_params* params) {
    k->flags |= NEEDS_SETUP;
    K_calculate(k);

    for (int i = 0; i < k->ndata; i++) {
        const double* comprow = compiled[i];
        int set = (int) comprow[T_SET];
        double n = K_getPar(k, set + DATA_SETS_SIZE);

        double diff = (comprow[T_PRED] - comprow[T_SVAL]);
        double s = comprow[T_ERR];

        f[i] = diff / sqrt(s * s + n * n);
        f[i + k->ndata] = sqrt(log(n * n + s * s)) * SIGN(log(n * n + s * s));

    }
    double fa = 0;
    for (int i = 0; i < k->ndata; i++)
        fa += f[i] * f[i];

    if (K_getLoglik(k) < params->min_chi) {
        params->min_chi = K_getLoglik(k);
        for (int i = 0; i < params->npars; i++)
            params->best[i] = *(params->pars[i]);
    }
}

int K_lm_f(const gsl_vector* x, void* params, gsl_vector* f) {
    ok_lm_params* sp = (ok_lm_params*) params;

    if (sp->status != PROGRESS_CONTINUE) {
        for (int i = 0; i < f->size; i++)
            f->data[i] = INVALID_NUMBER;
        return GSL_EINVAL;
    }

    ok_kernel* k = sp->k;

    double** pars = sp->pars;
    double** compiled = sp->compiled;
    unsigned int ndata = sp->ndata;

    for (int i = 0; i < x->size; i++)
        *(pars[i]) = x->data[i];

    K_lm_calc(k, f->data, ndata, (const double**) compiled, sp);

    sp->iterations++;

    if (sp->iterations % sp->every == 0 && k->progress != NULL) {
        k->chi2 = sp->min_chi;
        if ((k->progress)(sp->iterations, sp->maxiterations, k, __func__) != PROGRESS_CONTINUE) {
            sp->status = PROGRESS_STOP;
            for (int i = 0; i < ndata; i++)
                f->data[i] = INVALID_NUMBER;
            return GSL_EINVAL;
        }
    }

    return GSL_SUCCESS;
}

int K_lm_jac(const gsl_vector * x, void * params, gsl_matrix * J) {
    ok_lm_params* sp = (ok_lm_params*) params;

    if (sp->status != PROGRESS_CONTINUE) {
        for (int i = 0; i < J->size1 * J->size2; i++)
            J->data[i] = INVALID_NUMBER;
        return GSL_SUCCESS;
    }

    double* inp = x->data;
    double* jac = J->data;
    int m = x->size;

    ok_kernel* k = sp->k;
    double** pars = sp->pars;
    double* scale = sp->stepscale;
    int* parstype = sp->parstype;
    int ndata = sp->ndata;
    double** compiled = sp->compiled;
    double* f0 = sp->f0;
    double* f1 = sp->f1;
    double* f2 = sp->f2;
    double* f3 = sp->f3;

    for (int i = 0; i < m; i++)
        *(pars[i]) = inp[i];


    for (int j = 0; j < m; j++) {
        // Forward differences
        if ((parstype[j] == ECC || parstype[j] == MASS) && (inp[j] - 2. * scale[j] < 0)) {
            *(pars[j]) = inp[j];
            K_lm_calc(k, f0, ndata, (const double**) compiled, sp);

            *(pars[j]) = inp[j] + scale[j];
            K_lm_calc(k, f1, ndata, (const double**) compiled, sp);

            *(pars[j]) = inp[j] + 2 * scale[j];
            K_lm_calc(k, f2, ndata, (const double**) compiled, sp);

            for (int i = 0; i < ndata; i++) {
                jac[i * m + j] = (-0.5 * f2[i] + 2 * f1[i] - 1.5 * f0[i]) / (scale[j]);
            }

            *(pars[j]) = inp[j];
        } else {
            // Central differences
            if (sp->high_df) {
                *(pars[j]) = inp[j] - 2. * scale[j];
                K_lm_calc(k, f0, ndata, (const double**) compiled, sp);

                *(pars[j]) = inp[j] - scale[j];
                K_lm_calc(k, f1, ndata, (const double**) compiled, sp);

                *(pars[j]) = inp[j] + scale[j];
                K_lm_calc(k, f2, ndata, (const double**) compiled, sp);

                *(pars[j]) = inp[j] + 2. * scale[j];
                K_lm_calc(k, f3, ndata, (const double**) compiled, sp);

                for (int i = 0; i < ndata; i++) {
                    jac[i * m + j] = (f0[i] + 8. * f2[i] - 8. * f1[i] - f3[i]) / (12. * scale[j]);
                }
            } else {
                *(pars[j]) = inp[j] - scale[j];
                K_lm_calc(k, f1, ndata, (const double**) compiled, sp);

                *(pars[j]) = inp[j] + scale[j];
                K_lm_calc(k, f2, ndata, (const double**) compiled, sp);

                for (int i = 0; i < ndata; i++) {
                    jac[i * m + j] = (f2[i] - f1[i]) / (double) (2 * scale[j]);
                }
            }

            *(pars[j]) = inp[j];
        }

    }

    sp->iterations++;

    if (sp->iterations % sp->every == 0 && k->progress != NULL) {
        k->chi2 = sp->min_chi;
        if ((k->progress)(sp->iterations, sp->maxiterations, k, __func__) != PROGRESS_CONTINUE) {
            sp->status = PROGRESS_STOP;
            for (int i = 0; i < J->size1 * J->size2; i++)
                jac[i] = INVALID_NUMBER;
            return GSL_EINVAL;
        }
    }

    return GSL_SUCCESS;
}

int K_lm_fdf(const gsl_vector * x, void * params, gsl_vector* f, gsl_matrix * J) {
    ok_lm_params* sp = (ok_lm_params*) params;

    if (sp->status != PROGRESS_CONTINUE) {
        for (int i = 0; i < J->size1 * J->size2; i++)
            J->data[i] = INVALID_NUMBER;
        return GSL_SUCCESS;
    }

    double* inp = x->data;
    double* jac = J->data;
    int m = x->size;

    ok_kernel* k = sp->k;
    double** pars = sp->pars;
    double* scale = sp->stepscale;
    int* parstype = sp->parstype;
    int ndata = sp->ndata;
    double** compiled = sp->compiled;
    double* f0 = sp->f0;
    double* f1 = sp->f1;
    double* f2 = sp->f2;
    double* f3 = sp->f3;

    for (int i = 0; i < m; i++)
        *(pars[i]) = inp[i];

    K_lm_calc(k, f->data, ndata, (const double**) compiled, sp);

    for (int j = 0; j < m; j++) {
        // Forward differences
        if ((parstype[j] == ECC || parstype[j] == MASS) && (inp[j] - 2. * scale[j] < 0)) {

            *(pars[j]) = inp[j] + scale[j];
            K_lm_calc(k, f1, ndata, (const double**) compiled, sp);

            *(pars[j]) = inp[j] + 2 * scale[j];
            K_lm_calc(k, f2, ndata, (const double**) compiled, sp);

            for (int i = 0; i < ndata; i++) {
                jac[i * m + j] = (-0.5 * f2[i] + 2 * f1[i] - 1.5 * f->data[i]) / (scale[j]);
            }

            *(pars[j]) = inp[j];
        } else {
            // Central differences
            if (sp->high_df) {
                *(pars[j]) = inp[j] - 2. * scale[j];
                K_lm_calc(k, f0, ndata, (const double**) compiled, sp);

                *(pars[j]) = inp[j] - scale[j];
                K_lm_calc(k, f1, ndata, (const double**) compiled, sp);

                *(pars[j]) = inp[j] + scale[j];
                K_lm_calc(k, f2, ndata, (const double**) compiled, sp);

                *(pars[j]) = inp[j] + 2. * scale[j];
                K_lm_calc(k, f3, ndata, (const double**) compiled, sp);

                for (int i = 0; i < ndata; i++) {
                    jac[i * m + j] = (f0[i] + 8. * f2[i] - 8. * f1[i] - f3[i]) / (12. * scale[j]);
                }
            } else {
                *(pars[j]) = inp[j] + scale[j];
                K_lm_calc(k, f0, ndata, (const double**) compiled, sp);

                *(pars[j]) = inp[j] - scale[j];
                K_lm_calc(k, f1, ndata, (const double**) compiled, sp);

                for (int i = 0; i < ndata; i++) {
                    jac[i * m + j] = (f0[i] - f1[i]) / (2. * scale[j]);
                }
            }

            *(pars[j]) = inp[j];
        }
    }

    sp->iterations++;

    if (sp->iterations % sp->every == 0 && k->progress != NULL) {
        k->chi2 = sp->min_chi;
        if ((k->progress)(sp->iterations, sp->maxiterations, k, __func__) != PROGRESS_CONTINUE) {
            sp->status = PROGRESS_STOP;
            for (int i = 0; i < J->size1 * J->size2; i++)
                jac[i] = INVALID_NUMBER;
            return GSL_EINVAL;
        }
    }

    return GSL_SUCCESS;
}

int K_minimize_lm(ok_kernel* k, int maxiter, double params[]) {
    double min_chi_par = 1e-4;
    K_calculate(k);
    double prev_chi2 = K_getLoglik(k);
    bool high_df = false;
    int max_iter_at_scale = 200;
    double initial_st = 1.;

    int max_kt = 1;
    for (int i = 0; i < k->ndata; i++)
        if (k->compiled[i][T_FLAG] == T_TIMING) {
            high_df = true;
            max_kt = 2;
            max_iter_at_scale = 500;
            break;
        }

    if (params != NULL) {
        int i = 0;
        while (true) {
            if (params[i] == DONE)
                break;
            if (params[i] == OPT_LM_MINCHI_PAR)
                min_chi_par = params[i + 1];
            else if (params[i] == OPT_LM_HIGH_DF)
                high_df = !((int) params[i + 1] == 0);
            else if (params[i] == OPT_LM_MAX_ITER_AT_SCALE)
                max_iter_at_scale = (int) params[i + 1];
            else if (params[i] == OPT_LM_INITIAL_SCALE)
                initial_st = params[i + 1];
            i += 2;
        }
    }
    unsigned int npars = 0;

    // Count all the parameters to minimize on
    for (int i = 1; i < k->system->nplanets + 1; i++)
        for (int j = 0; j < ELEMENTS_SIZE; j++)
            npars += (MIGET(k->plFlags, i, j) & MINIMIZE ? 1 : 0);
    for (int i = 0; i < k->parFlags->size; i++)
        npars += (VIGET(k->parFlags, i) & MINIMIZE ? 1 : 0);

    if (npars == 0)
        return 0;

    // Create a pointer table (flat array -> matrices)
    double** pars = (double**) malloc(sizeof (double*) * npars);
    double prevpars[npars];

    double* steps = (double*) malloc(npars * sizeof (double));
    double* stepscale = (double*) malloc(npars * sizeof (double));
    int* parstype = (int*) malloc(npars * sizeof (int));

    gsl_vector* x = gsl_vector_alloc(npars);


    int idx = 0;
    for (int i = 1; i < k->system->nplanets + 1; i++)
        for (int j = 0; j < ELEMENTS_SIZE; j++)
            if (MIGET(k->plFlags, i, j) & MINIMIZE) {
                pars[idx] = gsl_matrix_ptr(k->system->elements, i, j);
                x->data[idx] = MGET(k->system->elements, i, j);
                prevpars[idx] = x->data[idx];
                steps[idx] = stepscale[idx] = MGET(k->plSteps, i, j);
                parstype[idx] = j;

                if (steps[idx] < 1e-10) {
                    printf("Warning: step for element %d of planet %d is <= 0\n", j, i);
                }

                idx++;
            }

    for (int i = 0; i < k->parFlags->size; i++)
        if (VIGET(k->parFlags, i) & MINIMIZE) {
            pars[idx] = gsl_vector_ptr(k->params, i);
            x->data[idx] = VGET(k->params, i);
            prevpars[idx] = x->data[idx];
            steps[idx] = stepscale[idx] = VGET(k->parSteps, i);
            parstype[idx] = PARTYPE_PAR;

            if (steps[idx] < 1e-10)
                printf("Warning: step for parameter %d is <= 0\n", i);

            idx++;
        }



    gsl_multifit_fdfsolver * s
            = gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, 2 * k->ndata, npars);




    ok_lm_params sp;
    sp.k = k;
    sp.pars = pars;
    sp.best = (double*) malloc(sizeof (double) * npars);
    sp.stepscale = stepscale;
    sp.compiled = k->compiled;
    sp.f0 = (double*) malloc(sizeof (double)*2 * k->ndata);
    sp.f1 = (double*) malloc(sizeof (double)*2 * k->ndata);
    sp.f2 = (double*) malloc(sizeof (double)*2 * k->ndata);
    sp.f3 = (double*) malloc(sizeof (double)*2 * k->ndata);
    sp.parstype = parstype;
    sp.ndata = 2 * k->ndata;
    sp.iterations = 0;
    sp.maxiterations = maxiter;
    sp.npars = npars;
    sp.every = (k->intMethod == KEPLER ? 10 : 1);
    sp.status = PROGRESS_CONTINUE;
    sp.high_df = high_df;
    sp.min_chi = K_getLoglik(k);
    sp.st = initial_st;

    for (int i = 0; i < npars; i++)
        sp.best[i] = *(pars[i]);

    gsl_multifit_function_fdf fdf;
    fdf.f = &K_lm_f;
    fdf.df = &K_lm_jac;
    fdf.fdf = &K_lm_fdf;
    fdf.n = 2 * k->ndata;
    fdf.p = npars;
    fdf.params = &sp;

    gsl_multifit_fdfsolver_set(s, &fdf, x);

    bool improved = true;
    int status = 0;
    int kt = 0;
    int tot_iter = 0;
    int iter_at_scale = 0;


    bool last_ditch = false;
    while (improved || last_ditch) {
        k->flags |= NEEDS_SETUP;
        iter_at_scale = 0;

        while (true) {
            double chi2 = sp.min_chi;
            int status = gsl_multifit_fdfsolver_iterate(s);

            iter_at_scale++;
            tot_iter++;
            if (chi2 - sp.min_chi > min_chi_par)
                iter_at_scale = 0;

            if (status || iter_at_scale > max_iter_at_scale) {
                break;
            }
        }

        gsl_multifit_fdfsolver_set(s, &fdf, x);


        for (int i = 0; i < npars; i++) {
            *(pars[i]) = sp.best[i];
            x->data[i] = sp.best[i];
        }

        k->flags |= NEEDS_SETUP;
        K_calculate(k);

        if (fabs(prev_chi2 - sp.min_chi) / fabs(sp.min_chi) > 1e-2 && iter_at_scale > 1) {
            kt = 0;
            last_ditch = false;
        } else {
            sp.st *= 0.1;
        }

        improved = (kt < max_kt || fabs(sp.min_chi - prev_chi2) / fabs(sp.min_chi) > 1e-2);

        if (last_ditch)
            break;

        if (!improved && kt <= 3) {
            last_ditch = true;
            sp.st *= 0.1;
        }

        kt++;
        //printf("-> %d %d %d %e %e, last_ditch=%d\n", iter_at_scale, kt, improved, sp.st, sp.min_chi, last_ditch);
        prev_chi2 = K_getLoglik(k);

        for (int idx = 0; idx < npars; idx++)
            stepscale[idx] = steps[idx] * sp.st;

        if (sp.iterations > maxiter || sp.st < 1e-12)
            break;
    }

    for (int i = 0; i < npars; i++) {
        *(pars[i]) = sp.best[i];
        x->data[i] = sp.best[i];
    }

    k->flags |= NEEDS_SETUP;
    K_calculate(k);


    free(sp.stepscale);
    free(sp.f0);
    free(sp.f1);
    free(sp.f2);
    free(sp.f3);
    free(sp.pars);
    free(sp.parstype);
    free(sp.best);

    if (sp.status == PROGRESS_STOP)
        return PROGRESS_STOP;
    else
        return status;

}
