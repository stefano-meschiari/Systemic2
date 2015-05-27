#include <stdio.h>
#include "systemic.h"
#include "kernel.h"
#include "mcmc.h"
#include "bootstrap.h"
#include "time.h"

int main() {
    FILE* fid = fopen("/Users/sm52286/Dropbox/Projects/HJST/psidraa.k", "r");
    
    ok_kernel* k = K_load(fid, 0);
    
    
}
