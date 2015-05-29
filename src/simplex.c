#include "math.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_vector.h"
#include "utils.h"
#include "kernel.h"

#define SIMPLEX_NON_FINITE_VALUE 1

typedef struct  {
    ok_kernel* k;
    double** pars;
    int flag;
} ok_simplex_params;

double K_simplex_f(const gsl_vector* cur, void* params) {
    
    ok_simplex_params* sp = (ok_simplex_params*) params;
    if (sp->flag != 0) {
        return 0.;
    }
    ok_kernel* k = sp->k;
    double** pars = sp->pars;
    
    for (int i = 0; i < cur->size; i++)
        *(pars[i]) = cur->data[i];
    
    k->flags |= NEEDS_SETUP;
    K_calculate(k);
    double res = k->minfunc(k);
    if (IS_NOT_FINITE(res)) {
        res = 0.;
        sp->flag = SIMPLEX_NON_FINITE_VALUE;
    }
    return res;
}

int K_minimize_simplex_iter(ok_kernel* k, int maxiter, double params[]) {
    double eps = 1e-8;
    double dminValue = 1e-4;
    //const int max_steps_wo_improvement = 10;
    
    int i = 0;
    if (params != NULL) {
        while (true) {
            if (params[i] == DONE)
                break;
            else if (round(params[i]) == OPT_EPS)
                eps = params[i+1];
            i+=2;
        }
    }
    int npars = 0;
    
    // Count all the parameters to minimize on
    for (int i = 1; i < k->system->nplanets + 1; i++)
        for (int j = 0; j < ELEMENTS_SIZE; j++)
            npars += (MIGET(k->plFlags, i, j) & MINIMIZE ? 1 : 0);
    for (int i = 0; i < k->parFlags->size; i++)
        npars += (VIGET(k->parFlags, i) & MINIMIZE ? 1 : 0);

    if (npars == 0)
        return 0;
    
    // Create a pointer table (flat array -> matrices)
    double** pars = (double**) malloc(sizeof(double*) * npars);
    
    int idx = 0;
    gsl_vector* x = gsl_vector_alloc(npars);
    gsl_vector* step = gsl_vector_alloc(npars);
    
    for (int i = 1; i < k->system->nplanets + 1; i++)
        for (int j = 0; j < ELEMENTS_SIZE; j++)
            if (MIGET(k->plFlags, i, j) & MINIMIZE) {
                
                pars[idx] = gsl_matrix_ptr(k->system->elements, i, j);
                
                x->data[idx] = MGET(k->system->elements, i, j);
                if ((j == PER) || (j == MASS))
                    step->data[idx] = MAX(MGET(k->plSteps, i, j), MGET(k->plSteps, i, j) * x->data[idx]);
                else
                    step->data[idx] = MGET(k->plSteps, i, j);
                
                
                if (step->data[idx] < 1e-10)
                    printf("Warning: step for element %d of planet %d is <= 0\n", j, i);
            
                idx++;
            }

    for (int i = 0; i < k->parFlags->size; i++)
        if (VIGET(k->parFlags, i) & MINIMIZE) {
            pars[idx] = gsl_vector_ptr(k->params, i);
            x->data[idx] = VGET(k->params, i);
            step->data[idx] = VGET(k->parSteps, i);
            
            if (step->data[idx] < 1e-10)
                printf("Warning: step for parameter %d is <= 0\n", i);
            
            idx++;
        }
        
    
    ok_simplex_params sp;
    sp.k = k;
    sp.pars = pars;
    sp.flag = 0;
    
    // Initialize minimizer
    gsl_multimin_fminimizer* s = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,
                npars);
    gsl_multimin_function f;
    f.n = npars;
    f.f = K_simplex_f;
    f.params = &sp;
    
    gsl_multimin_fminimizer_set(s, &f, x, step);
    
    
    int iter = 0;
    int steps_wo_improvement = 0;
    
    int status;
    
    // Loop until the desired eps (the size of the simplex) is reached.
    
    ok_progress pr = k->progress;
    const int every = (k->intMethod == 0 ? 30 : 1);
    double last_min_value = K_getMinValue(k);
    
    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if (status || sp.flag != 0)
            break;
        double size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, eps);
        if (status == GSL_SUCCESS)
            break;
        
        double min_value = gsl_multimin_fminimizer_minimum(s);
        /*if (last_min_value - min_value < dminValue) {
            steps_wo_improvement++;
            if (steps_wo_improvement > max_steps_wo_improvement)
                break;
        } else {
            last_min_value = min_value;
            steps_wo_improvement = 0;
        }*/
        if (pr != NULL && iter % every == 0) {
            k->chi2 = k->minfunc(k);
            if (pr(iter, maxiter, k, __func__) != PROGRESS_CONTINUE) {
                status = PROGRESS_STOP;
                break;
            }
        }
    } while (iter < maxiter);
    
    if (sp.flag != 0) {
        fprintf(stderr, "Non-finite value encountered by minimizer [%d].\n", sp.flag);
    }
    gsl_vector* best = gsl_multimin_fminimizer_x(s);
    K_simplex_f(best, &sp);
    K_calculate(k);
    gsl_multimin_fminimizer_free(s);
    free(pars);
    gsl_vector_free(step);
    gsl_vector_free(x);
    return status;
}

int K_minimize_simplex(ok_kernel* k, int maxiter, double params[]) {
    double dchi = 1e10;
    int status = PROGRESS_CONTINUE;
    int iter = 0;
    while (dchi > 1e-3 && status != PROGRESS_STOP) {
        double chi = k->minfunc(k);
        status = K_minimize_simplex_iter(k, maxiter, params);
        dchi = chi - k->minfunc(k);
        iter++;
    }
    return status;
}