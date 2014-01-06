#include "bootstrap.h"
#include "kernel.h"
#include <gsl/gsl_randist.h>

#ifndef JAVASCRIPT
#include "omp.h"
#else
#include "omp_shim.h"
#endif

ok_list* K_bootstrap(ok_kernel* k, int trials, int warmup, int malgo, int miter, double mparams[]) {
    ok_progress prog = k->progress;
    gsl_matrix* dev = NULL;
    
    int nthreads = omp_get_max_threads();
    bool invalid = false;
    
    if (prog != NULL) {
        int ret = prog(0, (int)((double) trials / (double) nthreads), k,
                "K_bootstrap");
        if (ret == PROGRESS_STOP) {
            return NULL;
        }
    }
    
    K_minimize(k, malgo, trials, mparams);
    ok_kernel* k0 = K_clone(k);
    
    ok_list* kl = KL_alloc(trials, k0);
    
    ok_list* wu = NULL;
    
    if (warmup > 0) {
        wu = KL_alloc(warmup, K_clone(k));
        
        #pragma omp parallel for
        for (int i = 0; i < warmup; i++) {
            if (invalid)
                continue;
            
            ok_kernel* k2 = K_cloneFlags(k, SHARE_FLAGS | SHARE_STEPS | SHARE_RANGES);
            k2->flags |= NEEDS_COMPILE | BOOTSTRAP_DATA;
            k2->progress = NULL;
            K_calculate(k2);
            
            K_minimize(k2, malgo, miter, mparams);
            KL_set(wu, i, k2->system->elements, k2->params, k2->minfunc(k2), 0);
            
            if (prog != NULL && omp_get_thread_num() == 0) {
                int ret = prog(i * nthreads, warmup, k2,
                    "K_bootstrap_warmup");
                if (ret == PROGRESS_STOP) {
                    invalid = true;
                }
            }
            
            K_free(k2);
        }
        
        dev = KL_getElementsStats(wu, STAT_STDDEV);
    }
    
    
    
    #pragma omp parallel for
    for (int i = 0; i < trials; i++) {
        if (invalid)
            continue;
       
        ok_kernel * k2 = K_cloneFlags(k, SHARE_FLAGS | SHARE_STEPS | SHARE_RANGES);
        k2->flags |= NEEDS_COMPILE | BOOTSTRAP_DATA;
        k2->progress = NULL;
        
        if (warmup > 0) {
            for (int i = 1; i <= k2->system->nplanets; i++) {
                for (int j = 0; j < ELEMENTS_SIZE; j++) 
                    if (MIGET(k2->plFlags, i, j) & MINIMIZE) {
                        double pert = gsl_ran_gaussian(k->rng, MGET(dev, i, j));
                        MINC(k2->system->elements, i, j, pert * MGET(k2->system->elements, i, j));
                    }
            }
        }
        
        
        K_calculate(k2);
        
        K_minimize(k2, malgo, miter, mparams);
        
        KL_set(kl, i, K_getAllElements(k2), ok_vector_copy(k2->params), k2->minfunc(k2), 0);
        
        if (prog != NULL && omp_get_thread_num() == 0) {
            int ret = prog(i * nthreads, trials, k2,
                    "K_bootstrap");
            if (ret == PROGRESS_STOP) {
                invalid = true;
            }
        }
        
        K_free(k2);
    }
    
    if (invalid) {
        if (wu != NULL)
            KL_free(wu);
        KL_free(kl);
        return NULL;
    }
        
    
    KL_free(wu);
    return kl;
}
