
#include <gsl/gsl_randist.h>

#ifndef JAVASCRIPT
#include "omp.h"
#else
#include "omp_shim.h"
#endif

#include "math.h"
#include "utils.h"
#include "kernel.h"

#define DISTINCT(a, b, c, x) (a != b && a != c && b != c && a != x)
#define IS_ANGLE(b) ((b) == MA || (b) == LOP || (b) == INC || (b) == NODE)
#define IS_LOG(b) ((b) == MASS || (b) == PER)
typedef struct {
    double* pars;
    double* old;
    double chi;
} ok_de_cand;

const double ok_de_min[ELEMENTS_SIZE] = { 1e-3, 1e-4, 0., 0, 0., 0., -100000., 0, 0, 0 };
const double ok_de_max[ELEMENTS_SIZE] = { 1e4, 100, 360., 0.99, 360., 360., 360., 0, 0, 0};

void K_validate(ok_kernel*);

bool ok_de_isCrossing(ok_kernel* k) {
    K_validate(k);
    double a[k->system->nplanets+1];
    for (int i = 1; i < k->system->nplanets+1; i++)
        a[i] = K_getElement(k, i, SMA);
    
    for (int i = 1; i < k->system->nplanets+1; i++)
        for (int j = 1; j < k->system->nplanets+1; j++) 
            if (i != j) {
                int imin = (a[i] < a[j] ? i : j);
                int imax = (a[i] > a[j] ? i : j);
                if (a[imin] * (1.+MGET(k->system->elements, imin, ECC)) > 
                    a[imax] * (1.-MGET(k->system->elements, imax, ECC)))
                    return true;
            }
    return false;
}

int K_minimize_de(ok_kernel* k, int trials, double params[]) {
    int status = PROGRESS_CONTINUE;
    K_calculate(k);
    double F_min = 0.5;
    double F_max = 1.;
    double CR = 0.2;
    int NP_fac = 25;
    bool use_step = false;
    
    int threads = omp_get_max_threads();
    
    ok_kernel_minimizer_pars mpars = K_getMinimizedVariables(k);
    
    if (mpars.npars == 0) {
        FREE_MINIMIZER_PARS(mpars);
        return 0;
    }
    
    int idx = 0;
    while (params != NULL) {
        if (params[idx] == DONE)
            break;
        else if (params[idx] == OPT_DE_CR)
            CR = params[idx+1];
        else if (params[idx] == OPT_DE_F_MIN)
            F_min = params[idx+1];
        else if (params[idx] == OPT_DE_F_MAX)
            F_max = params[idx+1];
        else if (params[idx] == OPT_DE_NP_FAC)
            NP_fac = (int) params[idx+1];
        else if (params[idx] == OPT_DE_USE_STEPS)
            use_step = (((int) params[idx+1]) != 0);
        
        idx+=2;
    }
    
    const int npars = mpars.npars;
    const int ncand = NP_fac * npars;
    
    ok_de_cand cand[ncand];
    
    double chi = K_getChi2(k);
    
    for (int i = 0; i < mpars.npars; i++) {
        double min = mpars.min[i] + 1e-16;
        double max = mpars.max[i];
        
        if (IS_INVALID(min) || IS_INVALID(max)) {
            if (mpars.type[i] == -1) {
                min = (IS_INVALID(min) ? -1e4 : min);
                max = (IS_INVALID(max) ? 1e4 : max);
            } else {
                min = (IS_INVALID(min) ? ok_de_min[mpars.type[i]] : min);
                max = (IS_INVALID(max) ? ok_de_max[mpars.type[i]] : max);
            }
        }

        mpars.min[i] = min;
        mpars.max[i] = max;
    }
    
    
    if (use_step) {
        for (int i = 0; i < ncand; i++) {
            cand[i].pars = malloc(sizeof(double) * npars);
            cand[i].old = malloc(sizeof(double) * npars);
            bool out_of_range = true;
            
            do {
                out_of_range = false;
                for (int j = 0; j < npars; j++) {
                    double r =  gsl_ran_gaussian(k->rng, mpars.steps[j]);
                    cand[i].pars[j] = *(mpars.pars[j]) + r;
                    out_of_range |= (cand[i].pars[j] < mpars.min[j]) || (cand[i].pars[j] > mpars.max[j]);
                }
                
                memcpy(cand[i].old, cand[i].pars, sizeof(double) * npars);
            } while (out_of_range);
        }
        
    } else {
        for (int i = 0; i < ncand; i++) {
            cand[i].pars = malloc(sizeof(double) * npars);
            cand[i].old = malloc(sizeof(double) * npars);

            for (int j = 0; j < npars; j++) {
                double min = mpars.min[j];
                double max = mpars.max[j];
                double r = gsl_rng_uniform(k->rng);
                
                if (mpars.type[j] == PER || mpars.type[j] == MASS) {
                    cand[i].pars[j] = exp(log(min) + (log(max) - log(min)) * r);
                } else if (mpars.type[j] == -1) {
                    cand[i].pars[j] = gsl_ran_gaussian(k->rng, 1);
                } else {
                    cand[i].pars[j] = min + (max - min) * r;
                } 

            }
            
            memcpy(cand[i].old, cand[i].pars, sizeof(double) * npars);
        }
    }
    
    ok_kernel* k_t[threads];
    ok_kernel_minimizer_pars mpars_t[threads];
    
    
    for (int i = 0; i < threads; i++) {
        k_t[i] = K_clone(k);
        mpars_t[i] = K_getMinimizedVariables(k_t[i]);
    }
    
    bool stop = false;
    #pragma omp parallel for 
    for (int i = 0; i < ncand; i++) {
        int th = omp_get_thread_num();
        for (int j = 0; j < npars; j++)
            *(mpars_t[th].pars[j]) = cand[i].pars[j];
        
        k_t[th]->flags |= NEEDS_SETUP;
        K_calculate(k_t[th]);
        fflush(stdout);
        cand[i].chi = k->minfunc(k_t[th]);
        
        if (k->progress != NULL && th == 0) {
            status = k->progress(0, trials, k_t[0], __func__);
            if (status == PROGRESS_STOP)
                stop = true;
        };
        if (stop)
            continue;
    }
    
    
    for (int tr = 0; tr < trials; tr++) {
        if (stop)
            break;
        double chi2 = 0;
        #pragma omp parallel for
        for (int x = 0; x < ncand; x++) {
            int nt = omp_get_thread_num();
            ok_kernel* kd = k_t[nt];
            ok_kernel_minimizer_pars mp = mpars_t[nt];
            double** parsd = mp.pars;

            int R = gsl_rng_uniform_int(kd->rng, npars);
            int a = 0, b = 0, c = 0;
            double F = gsl_rng_uniform(kd->rng) * (F_max-F_min) + F_min;

            while (! DISTINCT(a, b, c, x)) {
                a = gsl_rng_uniform_int(kd->rng, npars);
                b = gsl_rng_uniform_int(kd->rng, npars);
                c = gsl_rng_uniform_int(kd->rng, npars);
            }

            bool out_of_range = false;
                       
            for (int j = 0; j < npars; j++) {
                if ((x != R) && (gsl_rng_uniform(kd->rng) > CR)) {
                    *(parsd[j]) = cand[x].old[j];
                    cand[x].pars[j] = cand[x].old[j];
                } else {
                    double y;
                    if (IS_LOG(mpars.type[j])) 
                        y = exp(log(cand[a].old[j]) + F * (log(cand[b].old[j])-log(cand[c].old[j])));
                    else
                        y = cand[a].old[j] + F * (cand[b].old[j]-cand[c].old[j]);
                    
                    if (IS_ANGLE(mpars.type[j]))
                        y = DEGRANGE(y);
                    
                    if (y < mpars.min[j] || y > mpars.max[j]) {
                        out_of_range = true;
                        break;
                    }
                    *(parsd[j]) = y;
                    cand[x].pars[j] = y;
                }
            }

            
            kd->flags |= NEEDS_SETUP;
            bool cross = ok_de_isCrossing(kd);
            if (! (cross || out_of_range))
                K_calculate(kd);
            
            double chi_x = k->minfunc(kd);
            
            if ((chi_x < cand[x].chi || IS_INVALID(cand[x].chi)) && ! ok_de_isCrossing(kd) && ! out_of_range && ! IS_INVALID(chi_x)) {
                cand[x].chi = k->minfunc(kd);
            } else {
                for (int j = 0; j < npars; j++)
                    cand[x].pars[j] = cand[x].old[j];
            }
        }
        
        double min_chi_new = cand[0].chi;
        for (int x = 0; x < ncand; x++) {
            memcpy(cand[x].old, cand[x].pars, sizeof(double) * npars);
            if (cand[x].chi <= min_chi_new) {
                min_chi_new = cand[x].chi;
                for (int j = 0; j < npars; j++)
                    *(mpars_t[0].pars[j]) = cand[x].pars[j];
            }
            
            chi2 += cand[x].chi;
        }
       
       
        
        if (k->progress != NULL) {
            k_t[0]->flags |= NEEDS_SETUP;
            K_calculate(k_t[0]);
            status = k->progress(tr, trials, k_t[0], __func__);
            if (status == PROGRESS_STOP)
                break;
        };
    }
    
    for (int i = 0; i < threads; i++) {
        K_free(k_t[i]);
        FREE_MINIMIZER_PARS(mpars_t[i]);
    }
    
    
    for (int i = 0; i < ncand; i++) {
        if (cand[i].chi <= chi) {
            chi = cand[i].chi;
            for (int j = 0; j < npars; j++) {
                *(mpars.pars[j]) = cand[i].pars[j];
            }
        }
        
        free(cand[i].pars);
        free(cand[i].old);
    }
    
    k->flags |= NEEDS_SETUP;
    K_calculate(k);
    
    FREE_MINIMIZER_PARS(mpars);
    return status;
}