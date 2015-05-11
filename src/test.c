#include <stdio.h>
#include "systemic.h"
#include "kernel.h"
#include "mcmc.h"
#include "bootstrap.h"
#include "time.h"

int main() {
    FILE* fid = fopen("private/test.fit", "r");
    ok_kernel* k = K_alloc();
    
    
    K_addDataFromSystem(k, "/Users/sm52286/Projects/Systemic2/datafiles/14Her.sys");
    FILE* fid2 = fopen("private/test2.fit", "w");
    
    K_save(k, fid2);
    K_free(k);
 
}
