#include <gsl/gsl_randist.h>
#ifndef JAVASCRIPT
#include "omp.h"
#else
#include "omp_shim.h"
#endif

#include "math.h"
#include "utils.h"
#include "kernel.h"


double* K_minimize_sa_iter(ok_kernel* k2, const int N, const double T_0, const double alpha,
    const double* steps, double* best_chi, int* stop) {
    
    ok_kernel* k = K_clone(k2);
    ok_kernel_minimizer_pars mpars = K_getMinimizedVariables(k);
    double** pars = mpars.pars;
    int npars = mpars.npars;
    
    double* best_pars = (double*) malloc(sizeof(double) * npars);
    
    double old_pars[npars];
    K_calculate(k);
    
    double E = k->minfunc(k);
    double E_best = E;
    
    int n = 0;
    while (n < N) {
        double T = T_0 * pow(1-(double)n/(double) N, alpha);
        
        for (int i = 0; i < npars; i++) {
            old_pars[i] = *(pars[i]);
            *(pars[i]) += gsl_ran_gaussian(k->rng, steps[i]);
        }
        
        k->flags |= NEEDS_SETUP;
        K_calculate(k);
        
        if (k->minfunc(k) < E_best) {
            E_best = E;
            for (int i = 0; i < npars; i++)
                best_pars[i] = *(pars[i]);
        }
        
        if (*stop == PROGRESS_STOP) 
            break;
        if (omp_get_thread_num() == 0 && k2->progress != NULL) {
            *stop = k2->progress(n, N, k,
                    "K_sa");
        }
        
        printf("%d %e %e\n", n, T, E_best);
        
        if (k->minfunc(k) < E) {
            E = k->minfunc(k);
            n++;
            continue;
        }
        else {
            double P = exp(-fabs(k->minfunc(k)-E)/T);
            
            if (gsl_rng_uniform(k->rng) < P) {
                n++;
                E = k->minfunc(k);
            } else {
                for (int i = 0; i < npars; i++) {
                    *(pars[i]) = old_pars[i];
                }
                continue;
            }
        }
        
        
        
        
    }
    
    FREE_MINIMIZER_PARS(mpars);
    
    K_free(k);
    
    *best_chi = E_best;
    
    return best_pars;
}

int K_minimize_sa(ok_kernel* k, int trials, double params[]) {
    int status = PROGRESS_CONTINUE;
    K_calculate(k);
    double alpha = 2.;
    double T_0 =  k->minfunc(k);
    bool auto_step = true;
    int chains = omp_get_max_threads();
    
    ok_kernel_minimizer_pars mpars = K_getMinimizedVariables(k);
    
    int idx = 0;
    while (params != NULL) {
        if (params[idx] == DONE)
            break;
        else if (params[idx] == OPT_SA_T0)
            T_0 = params[idx+1];
        else if (params[idx] == OPT_SA_ALPHA)
            alpha = params[idx+1];
        else if (params[idx] == OPT_SA_CHAINS)
            chains = (int) params[idx+1];
        else if (params[idx] == OPT_SA_AUTO_STEPS)
            auto_step = ((int) params[idx+1] != 0);
        
        idx+=2;
    }
    
    double* best_pars[chains];
    double best_chi[chains];
    
    
    double chi = K_getChi2(k);
    if (auto_step) {
        for (int i = 0; i < mpars.npars; i++) {
            double p = *(mpars.pars[i]);
            double dchi;
            
            do {
                *(mpars.pars[i]) += mpars.steps[i];
                k->flags |= NEEDS_SETUP;
                K_calculate(k);
                dchi = fabs(chi - k->minfunc(k));
                mpars.steps[i] = 0.5 * mpars.steps[i] * (1 + 0.1*T_0/dchi);
                *(mpars.pars[i]) = p;       
            } while (fabs(dchi-0.1*T_0)/(0.1*T_0) > 0.1);
            printf("%d %e\n", i, mpars.steps[i]);
        }
    }
    
    #pragma omp parallel for
    for (int ch = 0; ch < chains; ch++) {
        best_pars[ch] = K_minimize_sa_iter(k, trials, T_0, alpha,
                mpars.steps,
                &(best_chi[ch]), &status);        
    }
    
    for (int ch = 1; ch < chains; ch++) {
        if (best_chi[ch] < best_chi[0]) {
            best_chi[0] = best_chi[ch];
            memcpy(best_pars[0], best_pars[ch], mpars.npars * sizeof(double));
        }
        free(best_pars[ch]);
    }
    
    for (int i = 0; i < mpars.npars; i++) {
        *(mpars.pars[i]) = best_pars[0][i];
    }
    
    free(best_pars[0]);
    FREE_MINIMIZER_PARS(mpars);
    
    K_calculate(k);
    return status;
}