/* 
 * File:   hermite.h
 * Author: sm52286
 *
 * Created on February 3, 2013, 3:08 PM
 */

// THIS IS WORK IN PROGRESS -- DO NOT USE

#ifndef HERMITE_H
#define	HERMITE_H

#include "systemic.h"

#ifdef	__cplusplus
extern "C" {
#endif

ok_system** ok_integrate_hermite(ok_system* initial, const gsl_vector* times, ok_integrator_options* options,
    ok_system** bag, int* error);



#ifdef	__cplusplus
}
#endif

#endif	/* HERMITE_H */

