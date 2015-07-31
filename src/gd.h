/* 
 * File:   gd.h
 * Author: sm52286
 *
 * Created on July 28, 2015, 2:46 PM
 */

#ifndef GD_H
#define	GD_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "kernel.h"
#include "systemic.h"

    int K_minimize_gd(ok_kernel* k, int maxiter, double params[]);


#ifdef	__cplusplus
}
#endif

#endif	/* GD_H */

