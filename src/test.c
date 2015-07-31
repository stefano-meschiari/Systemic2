#include <stdio.h>
#include "systemic.h"
#include "kernel.h"
#include "mcmc.h"
#include "bootstrap.h"
#include "time.h"
#include "gd.h"

int main() {

    ok_kernel* k = K_load(fopen("/Users/sm52286/Desktop/test.k", "r"), 0);
    K_minimize_gd(k, 1000, NULL);

}
