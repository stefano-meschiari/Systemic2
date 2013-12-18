/* 
 * File:   extras.h
 * Author: stefano
 *
 * Extras.(h|c) contains a few utilities that are not strictly related to the
 * core functionality.
 * 
 */


#ifndef EXTRAS_H
#define	EXTRAS_H


#include "systemic.h"

#ifdef	__cplusplus
extern "C" {
#endif

#define T_INAPPLICABLE -1
#define T_STABLE 0
#define T_UNSTABLE 1

    int K_isMstable_coplanar(const gsl_matrix* alle);

    double K_crossval_l1o(ok_kernel* k, int minalgo, int maxiter, double params[]);
    
#ifdef	__cplusplus
}
#endif

#endif	/* EXTRAS_H */

