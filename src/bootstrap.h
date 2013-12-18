/* 
 * File:   bootstrap.h
 * Author: stefano
 *
 * Created on March 1, 2012, 1:25 PM
 */

#ifndef BOOTSTRAP_H
#define	BOOTSTRAP_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "kl.h"

    ok_list* K_bootstrap(ok_kernel* k, int trials, int warmup, int malgo, int miter, double mparams[]);
    
#ifdef	__cplusplus
}
#endif

#endif	/* BOOTSTRAP_H */

