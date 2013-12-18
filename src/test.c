#include <stdio.h>
#include "systemic.h"
#include "kernel.h"
#include "mcmc.h"
#include "bootstrap.h"
#include "time.h"

int main() {
    FILE* fid = fopen("test.fit", "r");
    ok_kernel* k = K_load(fid, 0);
    k = K_clone(k);
    printf("%s %s\n", K_getDataName(k, 0), K_getDataName(k, 1));
    
    
}
