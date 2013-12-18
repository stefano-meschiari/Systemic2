/* 
 * File:   swift.h
 * Author: stefano
 *
 * Created on February 28, 2012, 1:38 PM
 */

#ifndef SWIFT_H
#define	SWIFT_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "systemic.h"
    

ok_system** ok_integrate_swift(ok_system* initial, const gsl_vector* times, const ok_integrator_options* options,
         ok_system** bag, int* error);

#ifdef	__cplusplus
}
#endif

#endif	/* SWIFT_H */

