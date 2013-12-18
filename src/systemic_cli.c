/* 
 * File:   systemic_cli.c
 * Author: stefano
 *
 * Created on February 7, 2013, 10:01 AM
 */

#include <stdio.h>
#include <stdlib.h>
#include "systemic.h"
#include "kernel.h"
#include "bootstrap.h"
#include "assert.h"
/*
 * 
 */
int main(int argc, char** argv) {
    char* loadfn = NULL;
    char* savefn = NULL;
    ok_kernel* k = NULL;
    ok_list* kl = NULL;
    int trials = 5000;
    int iters = 3000;
    int algo = SIMPLEX;
    
    for (int i = 1; i < argc; i+=2) {
        char* opt = argv[i];
        if (i + 1 >= argc) {
            printf("Wrong number of arguments.\n");
            exit(0);
        }
        char* arg = argv[i+1];
        
        
        if (strcmp(opt, "-file") == 0)
            loadfn = arg;
        else if (strcmp(opt, "-trials") == 0)
            trials = atoi(arg);
        else if (strcmp(opt, "-iters") == 0)
            iters = atoi(arg);
        else if (strcmp(opt, "-algo") == 0)
            algo = atoi(arg);
        else if (strcmp(opt, "-save") == 0)
            savefn = arg;
        else if (strcmp(opt, "-action") == 0) {
            if (strcmp(arg, "load") == 0) {
                printf("Loading kernel from file %s\n", loadfn);
                assert(loadfn != NULL);
                FILE* fid = fopen(loadfn, "r");
                assert(fid != NULL);
                k = K_load(fid, 0);
                assert(k != NULL);
                loadfn = NULL;
                fclose(fid);
                K_calculate(k);
                printf("Loaded, %d planets, %d datasets, chi2 = %e, intmethod = %d\n", k->system->nplanets, k->nsets, k->chi2,
                    k->intMethod);
            } else if (strcmp(arg, "bootstrap") == 0) {
                assert(k != NULL);
                printf("Running bootstrap with %d trials, %d iters, %d algo, saving to %s\n",
                    trials, iters, algo, savefn);
                kl = K_bootstrap(k, trials, 0, algo, iters, NULL);
                
                if (savefn != NULL) {
                    FILE* fid = fopen(savefn, "w");
                    KL_save(kl, fid);
                    savefn = NULL;
                    fclose(fid);
                }
            }
        }
    }
    
    return (EXIT_SUCCESS);
}

