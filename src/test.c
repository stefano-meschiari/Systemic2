#include <stdio.h>
#include "systemic.h"
#include "kernel.h"
#include "mcmc.h"
#include "bootstrap.h"
#include "time.h"

int main() {
    FILE* fid = fopen("/Users/sm52286/Dropbox/Projects/HJST/psidraa_b.k", "r");
    
    ok_kernel* k = K_load(fid, 0);
    K_calculate(k);
    printf("%e\n", k->chi2);
    K_setIntMethod(k, RK89);
    K_calculate(k);
    printf("%e\n", k->chi2);
    
}
