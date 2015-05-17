#include <stdio.h>
#include "systemic.h"
#include "kernel.h"
#include "mcmc.h"
#include "bootstrap.h"
#include "time.h"

int main() {
    FILE* fid = fopen("/Users/sm52286/Dropbox/Projects/HJST/psidraa.k", "r");
    
    ok_kernel* k = K_load(fid, 0);
    ok_kernel* kbuf[4];
    for (int i = 0; i < 4; i++)
        kbuf[i] = K_clone(k);
    K_mcmc_mult(kbuf, 4, 1, 100, 100, NULL, 1.1, NULL);
}
