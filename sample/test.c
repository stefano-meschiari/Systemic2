#include <stdio.h>
#include "systemic.h"
#include "kernel.h"
#include "mcmc.h"

int progress(int current, int max, void* state, const char* msg);

int main() {
    // Open the best fit for HD141399.fit (1)
    FILE* fid = fopen("private/hd141399.fit", "r");
    ok_kernel* k = K_load(fid, 0);
    
    // Set the callback function for k to function "progress" (2)
    K_setProgress(k, progress);
    
    // Activate the noise parameter (3)
    K_setParFlag(k, P_DATA_NOISE1, ACTIVE | MINIMIZE);
    K_setParFlag(k, P_DATA_NOISE2, ACTIVE | MINIMIZE);
    
    // Create two different starting states from perturbations
    // of k (4)
    K_perturb(k);
    
    ok_kernel* k2 = K_clone(k);
    K_perturb(k2);
    
    ok_kernel* kk[2] = {k, k2};
    
    // Call the MCMC routine, return the result in kl
    ok_list* kl = K_mcmc_mult(kk, 2, 1, 1000, 200, (double[3]) {OPT_MCMC_NMIN, 5000, DONE},
            1.1, NULL);
    
    // Calculate the median and mean absolute deviation
    // of the orbital elements... (5)
    printf("Final chain length: %d\n", kl->size);
    
    gsl_matrix* med = KL_getElementsStats(kl, STAT_MEDIAN);
    gsl_matrix* mad = KL_getElementsStats(kl, STAT_MAD);
    char* labels[5] = {"Period", "Mass", "MA", "Ecc", "Long. peri"};
    
    // ... and print it out.
    for (int i = 1; i < K_getNplanets(k)+1; i++) {
        for (int j = PER; j <= LOP; j++)
                printf("%s [%d]: %e +- %e\n", labels[j], i,
                        MGET(med, i, j), MGET(mad, i, j));                   
    }
    
    // Finally, save it to file.
    FILE* fid2 = fopen("private/hd141399.kl", "w");
    KL_save(kl, fid2);
    fclose(fid2);
    
}

int progress(int current, int max, void* state, const char* msg) {
    if (! strcmp(msg, "") == 0)
        printf("%s [%d/%d]\n", msg, current, max);
    return PROGRESS_CONTINUE;
}
