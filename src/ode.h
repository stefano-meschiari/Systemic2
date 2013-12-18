/* 
 * File:   ode.h
 * Author: stefano
 *
 * Created on February 1, 2013, 12:38 PM
 */

#ifndef ODE_H
#define	ODE_H

#include "systemic.h"

#ifdef	__cplusplus
extern "C" {
#endif

ok_system** ok_integrate_ode(ok_system* initial, const gsl_vector* times, ok_integrator_options* options, 
        ok_system** bag, int* error);


#ifdef	__cplusplus
}
#endif

#endif	/* ODE_H */

