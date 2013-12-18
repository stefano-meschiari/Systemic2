//
//  Integration.h
//  Systemic Framework
//
//  Created by Stefano Meschiari on 12/3/10.

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include "math.h"
#include "utils.h"
#include "mercury.h"
#include "assert.h"
#include "string.h"
#include "stdbool.h"
#include "stdarg.h"
#include "systemic.h"

#define OK_SUCCESS 0
#define OK_NOCONV 1

extern double ok_min_distance;
extern ok_integrator_options defoptions;
//typedef struct ok_transit_wspace 

/// Allocates a new system variable.
ok_system* ok_alloc_system(int nplanets);

/// Frees the memory allocated for the system variable
void ok_free_system(ok_system* system);

void ok_free_systems(ok_system** system, unsigned int len);

/// Returns a deep copy of the given system
ok_system* ok_copy_system(const ok_system* orig);

/// Returns a resized copy of the given system
void ok_resize_system(ok_system* system, int npnew);

/// Calculates a matrix containing the cartesian coordinates
/// of the system, given its orbital elements in system->orbits. 
void ok_el2cart(ok_system* system, gsl_matrix* xyz);

/// Calculates a matrix containing the cartesian coordinates
/// of the system, given its cartesian coords in system->xyz. 
void ok_cart2el(ok_system* system, gsl_matrix* els, bool internal);

/// Sets the system up. Call this function once all 
/// system parameters have been filled.
void ok_setup(ok_system* system);

/// Default force calculation (Newtonian gravitation)
int ok_force(double t, const double y[], double f[], void *params);

/// Default force jacobian
int ok_jac(double t, const double y[], double *dfdy, 
    double dfdt[], void *params);

/// Integrate routine
ok_system** ok_integrate(ok_system* initial, const gsl_vector* times, ok_integrator_options* options, const int integrator,
        ok_system** bag, int* error);

/// Extract rv, in m/s
double ok_get_rv(ok_system* sys);
gsl_matrix* ok_get_rvs(ok_system** sys, int len);

/// Extracts the coords; each row is a time snapshot of the xyz matrix
gsl_matrix* ok_get_xyzs(ok_system** bag, int len);

/// Extracts the elements; each row is a time snapshot of the elements.
gsl_matrix* ok_get_els(ok_system** bag, int len, bool internal);

gsl_vector* ok_find_transits(ok_system** bag, const int len, const int pidx, const int intMethod, const double eps, const int flags[], int* error);
/// Moves the origin to COM
void ok_to_cm(ok_system* system, gsl_matrix* xyz);

/// Moves the origin to the central star
void ok_to_star(ok_system* system, gsl_matrix* xyz);

int ok_find_closest_transit(ok_system* sys, const int pidx, ok_integrator_options* options, const int intMethod, const double eps, const int type, double* timeout, int* error);

double ok_pcalc(double Mcenter, double Mp, double a);
double ok_acalc(double Mcenter, double Mp, double P);
