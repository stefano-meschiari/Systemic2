//
//  Integration.c
//  Systemic Framework
//
//

#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_odeiv2.h>
#include "integration.h"
#include "ode.h"
#include "odex.h"

#ifndef JAVASCRIPT
#include "swift.h"
#endif

int ok_last_error(ok_system* system) {
    if (system->flag & INTEGRATION_FAILURE_CLOSE_ENCOUNTER)
        return INTEGRATION_FAILURE_CLOSE_ENCOUNTER;
    if (system->flag & INTEGRATION_FAILURE_CLOSE_ENCOUNTER_STAR)
        return INTEGRATION_FAILURE_CLOSE_ENCOUNTER_STAR;
    
    return INTEGRATION_SUCCESS;
}


double ok_min_distance = RJUP / AU;
ok_system* ok_alloc_system(int nplanets) {
    ok_system* system = (ok_system*) malloc(sizeof(ok_system));
    system->nplanets = nplanets;
    system->epoch = INVALID_NUMBER;
    system->time = INVALID_NUMBER;
    system->elements = gsl_matrix_calloc(nplanets + 1, ELEMENTS_SIZE);
    MSET(system->elements, 0, MASS, 1.);
    system->xyz = NULL;
    system->orbits = NULL;
    system->flag = 0;
    return system;
}

void ok_free_system(ok_system* system) {
    if (system == NULL)
        return;
    
    gsl_matrix_free(system->elements);
    if (system->xyz != NULL)
        gsl_matrix_free(system->xyz);
    if (system->orbits != NULL)
        gsl_matrix_free(system->orbits);
    
    free(system);
}

void ok_free_systems(ok_system** system, const unsigned int len) {
    for (unsigned int i = 0; i < len; i++)
        ok_free_system(system[i]);
    free(system);
}

ok_system* ok_copy_system(const ok_system* orig) {
    ok_system* system = ok_alloc_system(orig->nplanets);
    MATRIX_MEMCPY(system->elements, orig->elements);
    if (orig->orbits != NULL)
        system->orbits = ok_matrix_copy(orig->orbits);
    if (orig->xyz != NULL)
        system->xyz = ok_matrix_copy(orig->xyz);
    system->epoch = orig->epoch;
    system->time = orig->time;
    system->flag = orig->flag;
    return system;
}


void ok_copy_system_to(const ok_system* orig, ok_system* dest) {
    MATRIX_MEMCPY(dest->elements, orig->elements);
    if (orig->orbits != NULL) {
        if (dest->orbits == NULL)
            dest->orbits = ok_matrix_copy(orig->orbits);
        else
            MATRIX_MEMCPY(dest->orbits, orig->orbits);
    }
    if (orig->xyz != NULL) {
        if (dest->xyz == NULL)
            dest->xyz = ok_matrix_copy(orig->xyz);
        else
            MATRIX_MEMCPY(dest->xyz, orig->xyz);
    }
    dest->epoch = orig->epoch;
    dest->time = orig->time;
}


void ok_resize_system(ok_system* system, const int npnew) {
    
    gsl_matrix* old = system->elements;
    assert(old != NULL);
    assert(npnew >= 1);
    
    gsl_matrix* new = gsl_matrix_calloc(npnew, ELEMENTS_SIZE); 
    for (int i = 0; i < min(npnew, old->size1); i++) {
        for (int j = 0; j < ELEMENTS_SIZE; j++)
            MSET(new, i, j, MGET(old, i, j));
    }
    
    gsl_matrix_free(old);
    system->elements = new;
}

double ok_acalc(const double P, const double Mcenter, const double Mp) {
    double a = cbrt(P * P * (Mcenter + Mp) / (4.*M_PI*M_PI));
    return a;
}

double ok_pcalc(const double a, const double Mcenter, const double Mp) {
    return sqrt(a * a * a / (Mcenter + Mp)) * 2. * M_PI;
}

void ok_setup(ok_system* system) {
    system->flag &= ~(INTEGRATION_FAILURE_CLOSE_ENCOUNTER | INTEGRATION_FAILURE_CLOSE_ENCOUNTER_STAR);
    
    if (system->orbits != NULL && MROWS(system->orbits) != system->nplanets + 1) {
        gsl_matrix_free(system->orbits);
        system->orbits = NULL;
    }
    if (system->orbits == NULL)
        system->orbits = gsl_matrix_calloc(system->nplanets + 1, ELEMENTS_SIZE + 1);
    
    if (system->xyz != NULL && MROWS(system->xyz) != system->nplanets + 1) {
        gsl_matrix_free(system->xyz);
        system->xyz = NULL;
    }
    if (system->xyz == NULL)
        system->xyz = gsl_matrix_calloc(system->nplanets + 1, 7);
    
    
    // Copy elements matrix into orbits matrix (for instance, to carry on flags
    // such as ORD)
    MATRIX_MEMCPY(system->orbits, system->elements);
    
    double Mcenter = MSUN_TO_INT(MGET(system->elements, 0, MASS));

    MSET(system->orbits, 0, MASS, Mcenter);
    
    for (int i = 1; i < system->nplanets + 1; i++) {
        double mass = MJUP_TO_INT(MGET(system->elements, i, MASS));
        double per = MGET(system->elements, i, PER);
        
        double a = ok_acalc(per, Mcenter, mass);
        
        double e = MGET(system->elements, i, ECC);
        double lop = MGET(system->elements, i, LOP);
        double ma = MGET(system->elements, i, MA);
        double node = MGET(system->elements, i, NODE);
        double inc = MGET(system->elements, i, INC);
        
        MSET(system->orbits, i, PER, per);
        MSET(system->orbits, i, MASS, mass);
        MSET(system->orbits, i, MA, TO_RAD(ma));
        MSET(system->orbits, i, ECC, e);
        MSET(system->orbits, i, LOP, TO_RAD(lop));
        MSET(system->orbits, i, INC, TO_RAD(inc));
        MSET(system->orbits, i, NODE, TO_RAD(node));
        MSET(system->orbits, i, SMA, a);
        
        if (system->flag & JACOBI) {
            Mcenter += mass;
        }
    }
    
    system->time = system->epoch;
    ok_el2cart(system, system->xyz);
}



void ok_el2cart(ok_system* system, gsl_matrix* xyz) {
    gsl_matrix_set_zero(xyz);
    
    double Mcent = MGET(system->orbits, 0, MASS);
    MSET(xyz, 0, 0, Mcent);
    
    double xc[6] = {0., 0., 0., 0., 0., 0.};
    
    for (int i = 1; i < system->nplanets + 1; i++) {
        
        double mass = MGET(system->orbits, i, MASS);
        
        double ma = MGET(system->orbits, i, MA);
        double lop = MGET(system->orbits, i, LOP);    
        
        double e = MGET(system->orbits, i, ECC);
        double inc = MGET(system->orbits, i, INC);
        double node = MGET(system->orbits, i, NODE);
        double q = MGET(system->orbits, i, SMA) * (1. - e);
        
        double mu = mass + Mcent;
        double x, y, z, u, v, w;
        
        mco_el2x__(mu, q, e, inc, lop, node, ma,
                   &x, &y, &z, &u, &v, &w);
        
        if (fabs(y) < 1e-23) 
            y = 0.;
        if (fabs(v) < 1e-23) 
            v = 0.;
        
        MSET(xyz, i, 0, mass);
        MSET(xyz, i, X, x + xc[0]);
        MSET(xyz, i, Y, y + xc[1]);
        MSET(xyz, i, Z, z + xc[2]);
        MSET(xyz, i, VX, u + xc[3]);
        MSET(xyz, i, VY, v + xc[4]);
        MSET(xyz, i, VZ, w + xc[5]);
        
        if (system->flag & JACOBI) {
            xc[0] += (x * mass + xc[0] * Mcent)/(Mcent + mass);
            xc[1] += (y * mass + xc[1] * Mcent)/(Mcent + mass);
            xc[2] += (z * mass + xc[2] * Mcent)/(Mcent + mass);
            xc[3] += (u * mass + xc[3] * Mcent)/(Mcent + mass);
            xc[4] += (v * mass + xc[4] * Mcent)/(Mcent + mass);
            xc[5] += (w * mass + xc[5] * Mcent)/(Mcent + mass);
            Mcent += mass;
        }
    }
    
}

void ok_to_cm(ok_system* system, gsl_matrix* xyz) {
    double com[7] = {0.,0.,0.,0.,0.,0.,0.};
    
    for (int i = 0; i < system->nplanets + 1; i++) {
        com[0] += MGET(xyz, i, 0);
        for (int j = 1; j < 7; j++)
            com[j] += MGET(xyz, i, 0) * MGET(xyz, i, j);
    }
    
    for (int i = 0; i < system->nplanets + 1; i++) {
        for (int j = 1; j < 7; j++)
            MINC(xyz, i, j, - com[j]/com[0]);
    }   
}

void ok_to_star(ok_system* system, gsl_matrix* xyz) {
    for (int i = 1; i < system->nplanets + 1; i++)
        for (int j = 1; j < 7; j++)
            MINC(xyz, i, j, -MGET(xyz, 0, j));
    for (int j = 1; j < 7; j++)
        MSET(xyz, 0, j, 0.);
}

void ok_cart2el(ok_system* system, gsl_matrix* els, bool internal) {
    
    gsl_matrix* xyz = system->xyz;
    double xc[6] = {MGET(xyz, 0, X), MGET(xyz, 0, Y), MGET(xyz, 0, Z),
                MGET(xyz, 0, VX), MGET(xyz, 0, VY), MGET(xyz, 0, VZ)};
    double Mcent = MGET(xyz, 0, 0);
    
    for (int i = 1; i < system->nplanets + 1; i++) {
        double mass = MGET(xyz, i, 0);
        double mu = mass + Mcent;
        double x, y, z, u, v, w;
        
        x = MGET(xyz, i, X) - xc[0];
        y = MGET(xyz, i, Y) - xc[1];
        z = MGET(xyz, i, Z) - xc[2];
        u = MGET(xyz, i, VX) - xc[3];
        v = MGET(xyz, i, VY) - xc[4];
        w = MGET(xyz, i, VZ) - xc[5];

        double q, e, inc, p, n, l;
        
        mco_x2el__(&mu, &x, &y, &z, &u, &v, &w, &q, &e, &inc, &p, &n, &l);
        
        double a = q / (1-e);
        double P = ok_pcalc(a, Mcent, mass);
        
        MSET(els, i, PER, P);
        MSET(els, i, MASS, mass);
        MSET(els, i, MA, l);
        MSET(els, i, ECC, (e > 1e-6 ? e : 0.));
        MSET(els, i, INC, inc);
        MSET(els, i, NODE, n);
        MSET(els, i, SMA, a);
        MSET(els, i, LOP, p);
        
        if (system->flag & JACOBI) {
            xc[0] = (xc[0] * Mcent + MGET(xyz, i, X) * mass) / (Mcent + mass);
            xc[1] = (xc[1] * Mcent + MGET(xyz, i, Y) * mass) / (Mcent + mass);
            xc[2] = (xc[2] * Mcent + MGET(xyz, i, Z) * mass) / (Mcent + mass);
            xc[3] = (xc[3] * Mcent + MGET(xyz, i, VX) * mass) / (Mcent + mass);
            xc[4] = (xc[4] * Mcent + MGET(xyz, i, VY) * mass) / (Mcent + mass);
            xc[5] = (xc[5] * Mcent + MGET(xyz, i, VZ) * mass) / (Mcent + mass);
            
            Mcent += mass;   
        }
    }
    
}

unsigned long int ok_force_counter = 0;

int ok_force(double t, const double y[], double f[], void* params) {
    ok_force_counter++;
    ok_system* system = (ok_system*) params;
    const int N = system->nplanets + 1;
    
    for (int i = 0; i < 7 * N; i++)
        f[i] = 0.;
    
    for (int i = 0; i < N; i++) {
        f[i * 7 + 1] = y[i * 7 + 4];
        f[i * 7 + 2] = y[i * 7 + 5];
        f[i * 7 + 3] = y[i * 7 + 6];
        
        for (int j = i + 1; j < N; j++) {
            double m1 = y[i * 7];
            double m2 = y[j * 7];
            
            double i_rsq = 1./(sqr(y[i * 7 + 1] - y[j * 7 + 1]) +
                sqr(y[i * 7 + 2] - y[j * 7 + 2]) +
                sqr(y[i * 7 + 3] - y[j * 7 + 3]));
            
            double i_r = sqrt(i_rsq);
            
            if (i_r < ok_min_distance) {
                if (i == 0 || j == 0) {
                    system->flag |= INTEGRATION_FAILURE_CLOSE_ENCOUNTER_STAR;
                    return INTEGRATION_FAILURE_CLOSE_ENCOUNTER_STAR;
                }
                else {
                    system->flag |= INTEGRATION_FAILURE_CLOSE_ENCOUNTER;
                    return INTEGRATION_FAILURE_CLOSE_ENCOUNTER;
                }
            }
            
            double a1 = -m2 * i_rsq * i_r;
            double a2 = -m1 * i_rsq * i_r;
            
            for (int d = 1; d <= 3; d++) {
                double w = (y[i * 7 + d] - y[j * 7 + d]);
                f[i * 7 + d + 3] += w * a1;
                f[j * 7 + d + 3] -= w * a2;
            }
        }        
    }

    /*
    static unsigned int calls = 0;
    calls++;
    if (calls % 1000000 == 0)
        printf("rk: %u\n", calls);
    */
    
    return INTEGRATION_SUCCESS;
}

int ok_force_on(double t, const double y[], double f[], void* params, int i) {
    
    ok_system* system = (ok_system*) params;
    
    const int N = system->nplanets + 1;
    f[0] = f[1] = f[2] = 0.;
    
    for (int j = 0; j < N; j++) {
        if (j == i)
            continue;
        
        double m2 = y[j * 7];

        double i_rsq = 1./(sqr(y[i * 7 + 1] - y[j * 7 + 1]) +
            sqr(y[i * 7 + 2] - y[j * 7 + 2]) +
            sqr(y[i * 7 + 3] - y[j * 7 + 3]));

        double i_r = sqrt(i_rsq);

        if (i_r < ok_min_distance) {
            if (i == 0 || j == 0)
                return INTEGRATION_FAILURE_CLOSE_ENCOUNTER_STAR;
            else
                return INTEGRATION_FAILURE_CLOSE_ENCOUNTER;
        }

        double a1 = -m2 * i_rsq * i_r;

        for (int d = 1; d <= 3; d++) {
            double w = (y[i * 7 + d] - y[j * 7 + d]);
            f[d-1] += w * a1;
        }
    }        
    
    return INTEGRATION_SUCCESS;
}


int ok_force_jerk(double t, const double y[], double f[], double jerk[], void* params) {
    ok_system* system = (ok_system*) params;
    
    const int N = system->nplanets + 1;
    for (int i = 0; i < 7*N; i++) {
        f[i] = 0.;
        if (i < 3*N+1)
            jerk[i] = 0.;
    }
    
    double t_min = 0; 
    
    for (int i = 0; i < N; i++) {
        const double* xi = y + i*7;
        double* fi = f + i*7;
        fi[1] = xi[4];
        fi[2] = xi[5];
        fi[3] = xi[6];
        
        for (int j = i + 1; j < N; j++) {
            const double* xj = y + j * 7;
            double* fj = f + j*7;
            const double i_rsq = 1./(sqr(xi[1] - xj[1]) +
                sqr(xi[2] - xj[2]) +
                sqr(xi[3] - xj[3]));
            
            const double vij2 = sqr(xi[4]-xj[4]) + sqr(xi[5]-xj[5]) + sqr(xi[6] - xj[6]);
            const double i_r = sqrt(i_rsq);
            const double i_r3 = i_rsq * i_r;
            
            if (i_r < ok_min_distance) {
                if (i == 0 || j == 0) {
                    system->flag |= INTEGRATION_FAILURE_CLOSE_ENCOUNTER_STAR;
                    return INTEGRATION_FAILURE_CLOSE_ENCOUNTER_STAR;
                }
                else {
                    system->flag |= INTEGRATION_FAILURE_CLOSE_ENCOUNTER;
                    return INTEGRATION_FAILURE_CLOSE_ENCOUNTER;
                }
            }
            
            const double a1 = -xj[0] * i_r3;
            const double a2 = -xi[0] * i_r3;
            
            const double jw =  i_rsq * ((xi[1] - xj[1]) * (xi[4] - xj[4]) +
                (xi[2]-xj[2]) * (xi[5]-xj[5]) + (xi[3]-xj[3]) * (xi[6]-xj[6]));
            
            for (int d = 1; d <= 3; d++) {
                const double w = (xi[d] - xj[d]);
                const double vd = (xi[d+3]-xj[d+3]);
                fi[d + 3] += w * a1;
                fj[d + 3] -= w * a2;
                
                jerk[i * 3 + d - 1] += a1 * (vd - 3.*jw * w);
                jerk[j * 3 + d - 1] -= a2 * (vd - 3.*jw * w);
            }
            
            t_min = MAX(MAX(MAX((i_rsq * vij2), t_min), -a1), -a2);
        }
    }
    
    for (int i = 0; i < N; i++)
        t_min = MAX(t_min, 
                (sqr(jerk[i * 3]) + sqr(jerk[i * 3 + 1]) + sqr(jerk[i * 3 + 2]))/
                (sqr(f[i*7 + 4]) + sqr(f[i * 7 + 5]) + sqr(f[i * 7 + 6])));
    
    jerk[3*N] = t_min;
    /*
    static unsigned int calls = 0;
    calls++;
    */
    return GSL_SUCCESS;
}


int ok_star_force(double t, const double y[], double f[], void* params) {
    ok_system* system = (ok_system*) params;
    
    const int N = system->nplanets + 1;
    
    for (int i = 0; i < 7 * N; i++)
        f[i] = 0.;
    
    for (int i = 0; i < N; i++) {
        f[i * 7 + 1] = y[i * 7 + 4];
        f[i * 7 + 2] = y[i * 7 + 5];
        f[i * 7 + 3] = y[i * 7 + 6];
        
        if (i != 0) {
            double m2 = y[i * 7];
            
            double i_rsq = 1./(sqr(y[i * 7 + 1] - y[1]) +
                sqr(y[i * 7 + 2] - y[2]) +
                sqr(y[i * 7 + 3] - y[3]));
            
            double i_r = sqrt(i_rsq);
            
            if (i_r < ok_min_distance) {
                return GSL_ERANGE;
            }
            
            double a1 = -m2 * i_rsq * i_r;
            
            for (int d = 1; d <= 3; d++) {
                double w = (y[i * 7 + d] - y[d]);
                f[i * 7 + d + 3] += w * a1;
            }
        }        
    }


    return GSL_SUCCESS;
}

int ok_star_force_on(double t, const double y[], double f[], void* params, int i) {
    if (i == 0) {
        f[0] = f[1] = f[2] = 0.;
        return INTEGRATION_SUCCESS;
    }
    
    double m2 = y[i * 7];

    double i_rsq = 1./(sqr(y[i * 7 + 1] - y[1]) +
        sqr(y[i * 7 + 2] - y[2]) +
        sqr(y[i * 7 + 3] - y[3]));

    double i_r = sqrt(i_rsq);

    if (i_r < ok_min_distance) {
        return GSL_ERANGE;
    }

    double a1 = -m2 * i_rsq * i_r;

    for (int d = 1; d <= 3; d++) {
        double w = (y[i * 7 + d] - y[d]);
        f[d-1] = w * a1;
    }

    return GSL_SUCCESS;
}

int ok_jac(double t, const double y[], double *dfdy, 
           double dfdt[], void *params) {
    ok_system* system = (ok_system*) params;
    const int N = system->nplanets + 1;
    
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 7 * N, 7 * N);
    gsl_matrix * jac = &dfdy_mat.matrix; 
    gsl_matrix_set_zero(jac);
    
    for (int i = 0; i < N; i++) {
        // the N * 7-th row [mass row] is always 0
        
        // df_x / dx = 0
        // df_x_i / dv_i = 1        
        MSET(jac, i * 7 + 1, i * 7 + 4, 1);
        MSET(jac, i * 7 + 2, i * 7 + 5, 1);
        MSET(jac, i * 7 + 3, i * 7 + 6, 1);        
        
        // df_v / dv = 0
        // df_v_i / dx =
        
        for (int j = i+1; j < N; j++) {
            const double inv_rij = 1./sqrt(sqr(y[i * 7 + 1] - y[j * 7 + 1]) + 
                                        sqr(y[i * 7 + 2] - y[j * 7 + 2]) + 
                                        sqr(y[i * 7 + 3] - y[j * 7 + 3]));
            const double inv_rij_3 = inv_rij * inv_rij * inv_rij;
            const double inv_rij_5 = inv_rij_3 * inv_rij * inv_rij;
            
            for (int m = 1; m <= 3; m++) 
                for (int n = 1; n <= 3; n++) {
                    double term = - (y[i * 7 + m] - y[j * 7 + m])
                         * (y[i * 7 + n] - y[j * 7 + n]) * inv_rij_5;
                    
                    if (m == n)
                        term += - inv_rij_3;
                    
                    MINC(jac, i * 7 + m + 3, i * 7 + n, y[j * 7] * term);
                    MINC(jac, j * 7 + m + 3, j * 7 + n, - y[i * 7] * term);   
                }

        }
    }
    return GSL_SUCCESS;
}

ok_integrator_options defoptions = { 1e-13, 1e-13, 0.15, 1., 1e-6, 2, true, &ok_force, &ok_jac, &ok_force_jerk, NULL, NULL, NULL };

/// This routine integrates the system in time. An array of snapshots taken at each time specified
/// by the times vector is returned. 
ok_system** ok_integrate_gsl(ok_system* initial, const gsl_vector* times, ok_integrator_options* options, const gsl_odeiv2_step_type * solver,
        ok_system** bag, int* error) {
    
    // Check that the system has been set-up
    assert(initial->xyz != NULL);
    // Check input arguments
    assert(times != NULL);
    assert(solver != NULL);
    
    const double startTime = initial->epoch;
    const int NDIMS = initial->nplanets + 1;
    
    // Allocate the return array of snapshots
    const int SAMPLES = times->size;
   
    if (bag == NULL) {
        bag = (ok_system**) calloc(SAMPLES, sizeof(ok_system*));
        for (int i = 0; i < SAMPLES; i++) {
            bag[i] = ok_copy_system(initial);
            bag[i]->epoch = initial->epoch;
            bag[i]->xyz = (bag[i]->xyz != NULL ? bag[i]->xyz : gsl_matrix_alloc(initial->nplanets+1, 7));
        }
    } else {
        for (int i = 0; i < SAMPLES; i++) {
            bag[i]->epoch = initial->epoch;
            bag[i]->flag = initial->flag;
            //MATRIX_MEMCPY(bag[i]->elements, initial->elements);
            //MATRIX_MEMCPY(bag[i]->orbits, initial->orbits);
            bag[i]->xyz = (bag[i]->xyz != NULL ? bag[i]->xyz : gsl_matrix_alloc(initial->nplanets+1, 7));
        }
    }

    // Initialize the GSL structures to solve the ODE.
    gsl_odeiv2_step *stepper = gsl_odeiv2_step_alloc(solver, NDIMS * 7);
    gsl_odeiv2_control * control = gsl_odeiv2_control_standard_new(options->abs_acc, options->rel_acc, 1., 1.);
    gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (NDIMS * 7);
    
    gsl_odeiv2_system eqns;
    eqns.dimension = NDIMS * 7;
    eqns.function = options->force;
    eqns.params = initial;
    
    // Allocates the temporary buffer to hold cartesian coordinates
    double prevTime = startTime;
    gsl_matrix* prevOrbits = initial->orbits;
    
    double h = 0.01;
    
    if (options->buffer == NULL || options->buffer->size < NDIMS*7) {
        if (options->buffer != NULL)
            gsl_vector_free(options->buffer);
        options->buffer = gsl_vector_alloc(NDIMS*7);
    }
    
    double* xyz = options->buffer->data;
    MATRIX_MEMCPY_TOARRAY(xyz, initial->xyz);
    
    ok_progress progress = options->progress;
    
    // Loop through the times vector
    for (int i = 0; i < SAMPLES; i++) {
        double time = times->data[i];

        // Integrate between prevTime and time
        if (fabs(time - prevTime) > 1e-10) {
            while (fabs(time - prevTime) > 1e-10) {
                h = copysign(h, time-prevTime);
                const int result = gsl_odeiv2_evolve_apply(e, control, stepper, &eqns, &prevTime, time, &h, xyz);
                
                if (result != GSL_SUCCESS) {
                    if (i == 0) {
                        if (error != NULL) {
                            *error = ok_last_error(bag[i]);
                        }
                        
                        for (int j = 0; j < SAMPLES; j++)
                            ok_free_system(bag[i]);
                        free(bag);
                        
                        
                        gsl_odeiv2_control_free(control);
                        gsl_odeiv2_evolve_free(e);
                        gsl_odeiv2_step_free(stepper);

                        
                        return NULL;
                    } else {
                        if (error != NULL) {
                            *error = ok_last_error(bag[i]);
                        }
                        
                        for (int j = i - 1; j < SAMPLES; j++) {
                            bag[j]->time = bag[j]->epoch = bag[i]->time;
                            gsl_matrix_set_all(bag[j]->xyz, INVALID_NUMBER);
                            gsl_matrix_set_all(bag[j]->orbits, INVALID_NUMBER);
                        }
                        
                        gsl_odeiv2_control_free(control);
                        gsl_odeiv2_evolve_free(e);
                        gsl_odeiv2_step_free(stepper);
                        
                        
                        return bag;
                    } 
                }
            }
        } else {
            if (i == 0)
                MATRIX_MEMCPY_TOARRAY(xyz, initial->xyz);
            else
                MATRIX_MEMCPY_TOARRAY(xyz, bag[i-1]->xyz);
        };
        
        
        // Set up the return vector
        bag[i]->time = bag[i]->epoch = time;
       
        MATRIX_MEMCPY_FROMARRAY(bag[i]->xyz, xyz);
        if (options == NULL || options->calc_elements) {
                MATRIX_MEMCPY(bag[i]->orbits, prevOrbits);
                ok_cart2el(bag[i], bag[i]->orbits, true);
        }
        prevOrbits = bag[i]->orbits;
        
        
        // Ensures that the force/jacobian routines are passed the most recent state of
        // the system
        eqns.params = bag[i];
        prevTime = time;
        
        if (progress != NULL) {
            
            int ret = progress(i, SAMPLES, NULL, "Integration");
            
            if (ret == PROGRESS_STOP) {
                for (int i = 0; i < SAMPLES; i++)
                        ok_free_system(bag[i]);
                
                free(bag);
                gsl_odeiv2_control_free(control);
                gsl_odeiv2_evolve_free(e);
                gsl_odeiv2_step_free(stepper);

                if (error != NULL) {
                    *error = INTEGRATION_FAILURE_STOPPED;
                }
                return NULL;
            }
        }
    }
    
    gsl_odeiv2_control_free(control);
    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_step_free(stepper);
    
    if (error != NULL)
        *error = INTEGRATION_SUCCESS;

    return bag;
}

ok_system** ok_integrate_kep(ok_system* initial, const gsl_vector* times, 
                             ok_integrator_options* options, ok_system** bag, int* error) {
    
    assert(initial->orbits != NULL);
    
    const double startTime = initial->epoch;
    const int NDIMS = initial->nplanets + 1;
    
    // Allocate the return array of snapshots
    const int SAMPLES = times->size;
    if (bag == NULL) {
        bag = (ok_system**) calloc(SAMPLES, sizeof(ok_system*));
        for (int i = 0; i < SAMPLES; i++) {
            bag[i] = ok_copy_system(initial);
            bag[i]->epoch = initial->epoch;
            bag[i]->xyz = (bag[i]->xyz != NULL ? bag[i]->xyz : gsl_matrix_alloc(initial->nplanets+1, 7));
        }
    } else {
        for (int i = 0; i < SAMPLES; i++) {
            bag[i]->epoch = initial->epoch;
            bag[i]->flag = initial->flag;
            MATRIX_MEMCPY(bag[i]->elements, initial->elements);
            MATRIX_MEMCPY(bag[i]->orbits, initial->orbits);
            bag[i]->xyz = (bag[i]->xyz != NULL ? bag[i]->xyz : gsl_matrix_alloc(initial->nplanets+1, 7));
        }
    }
    
    double prevTime = startTime;
    ok_system* prevSystem = initial;
    
    double n[NDIMS];
    for (int i = 1; i < NDIMS; i++)
        n[i] = 2*M_PI/MGET(prevSystem->orbits, i, PER);
    
    // Loop through the times vector
    for (int i = 0; i < SAMPLES; i++) {
        
        double dt = times->data[i] - prevTime;
        for (int j = 1; j < NDIMS; j++) {
            double ma = RADRANGE(MGET(prevSystem->orbits, j, MA) + n[j] * dt);
            MSET(bag[i]->orbits, j, MA, ma);
        }
        
        bag[i]->time = times->data[i];
        bag[i]->epoch = bag[i]->time;
        prevTime = times->data[i];
        prevSystem = bag[i];
        ok_el2cart(bag[i], bag[i]->xyz);
    }
    
    if (error != NULL)
        *error = INTEGRATION_SUCCESS;
    
    return bag;
}

ok_integrator* ok_integrators[4];

ok_system** ok_integrate(ok_system* initial, const gsl_vector* times, ok_integrator_options* options, const int integrator,
        ok_system** bag, int* error) {
    // Check that the system has been set-up
    assert(initial->xyz != NULL);
    // Check input arguments
    assert(times != NULL);
    
    assert(! IS_INVALID(initial->epoch));
    assert(options != NULL);
    
    if (IS_INVALID(initial->time))
        initial->time = initial->epoch;
    
    
    switch (integrator) {
        case KEPLER:
            return ok_integrate_kep(initial, times, options, bag, error);
            break;
        case RK45:
            return ok_integrate_gsl(initial, times, options, gsl_odeiv2_step_rkf45, bag, error);
            break;
        case RK89:
            return ok_integrate_gsl(initial, times, options, gsl_odeiv2_step_rk8pd, bag, error);
            break;
#ifndef JAVASCRIPT
        case ADAMS: 
            return ok_integrate_ode(initial, times, options, bag, error);
        case BULIRSCHSTOER: 
            return ok_integrate_odex(initial, times, options, bag, error);
        case SWIFTRMVS: 
            return ok_integrate_swift(initial, times, options, bag, error);
#endif         
    }
    
    return NULL;
}



double ok_get_rv(ok_system* sys) {
    assert(sys->xyz != NULL);    
    ok_to_cm(sys, sys->xyz);
    double rv = -MGET(sys->xyz, 0, VZ);
    rv = AUPDAY_TO_MPS(rv);
    return rv;
}

gsl_matrix* ok_get_rvs(ok_system** sys, int len) {
    gsl_matrix* ret = gsl_matrix_alloc(len, 2);
    
    for (int i = 0; i < len; i++) {
        MSET(ret, i, 0, sys[i]->time);
        MSET(ret, i, 1, ok_get_rv(sys[i]));
    }
    return ret;
}

gsl_matrix* ok_get_xyzs(ok_system** bag, int len) {
    int N = bag[0]->nplanets + 1;
    gsl_matrix* xyzs = gsl_matrix_calloc(len, 7 * N + 1);
    
    for (int i = 0; i < len; i++) {
        assert(bag[i]->xyz != NULL);
        
        MSET(xyzs, i, 0, bag[i]->time);
        
        for (int j = 0; j < N; j++)
            for (int p = 0; p < 7; p++)
                MSET(xyzs, i, j * 7 + p + 1, MGET(bag[i]->xyz, j, p));
    }
    
    return xyzs;
}

gsl_matrix* ok_get_els(ok_system** bag, int len, bool internal) {
    int N = bag[0]->nplanets + 1;
    gsl_matrix* els = gsl_matrix_calloc(len, N * ELEMENTS_SIZE + 1);
    
    for (int i = 0; i < len; i++) {
        assert(bag[i]->orbits != NULL);
        
        MSET(els, i, 0, bag[i]->time);
        
        for (int j = 0; j < N; j++)
            for (int p = 0; p < ELEMENTS_SIZE; p++) {
                double v = MGET(bag[i]->orbits, j, p);
                
                if (p == MASS && !internal)
                    v = INT_TO_MJUP(v);
                else if ((p == MA || p == LOP || p == INC || p == NODE) && !internal)
                    v = TO_DEG(v);
                
                MSET(els, i, j * ELEMENTS_SIZE + p + 1, v);
            }
    }
    
    return els;
}



/**
 * Returns the time of the transit closest to the time of the current state of the system for the specified planet.
 * @param state pointer to a ok_system structure containing cartesian coords + elements, time and epoch.
 * @param plidx index of the planet (1..n)
 * @param options NULL for the default integrator options, or pass a pointer to an ok_integrator_options structure
 * @param intMethod integration method (one of the builtin integrators: KEPLER, RK5, RK8, HERMITE, HERMITE_ACC, SWIFT_BS, SWIFT_RMVS or
 * a custom integrator)
 * @param eps tolerance on the time
 * @return Closest transit to the time specified in state
 */
int ok_find_closest_transit(ok_system* state, const int plidx, ok_integrator_options* options, const int intMethod, const double eps, const int type, double* timeout, int* error) {
    // Find planet tagged with the specified index (ORD column), since internally we might keep planets sorted
    // by period.
    
    if (IS_INVALID(state->time) || IS_INVALID(state->epoch)) {
        *(timeout) = INVALID_NUMBER;
        return 0;
    }
    
    int pidx = plidx;
    int retval = OK_SUCCESS;
    
    ok_to_star(state, state->xyz);
    options->iterations = 1;
    
    for (int i = 1; i < state->elements->size1; i++) {
        if (pidx == (int) MGET(state->elements, i, ORD)) {
            pidx = i;
            break;
        }
    }
    
    double t;
    

    // Calculate approximate time of closest transit. Angles in the orbits matrix are in radians,
    // with epoch stored in state->epoch.
    double aop = MGET(state->orbits, pidx, LOP) - MGET(state->orbits, pidx, NODE);
    double per = MGET(state->orbits, pidx, PER);
    double n = 2*M_PI/per;
    double tau = (n * state->time - MGET(state->orbits, pidx, MA))/n;
    double e = MGET(state->orbits, pidx, ECC);

    // Elliptic expansion; fail if e > 0.6. Should handle this case later.
    if (e > 0.6) {
        *timeout = INVALID_NUMBER;
        return OK_NOCONV;
    }

    double ftra = RADRANGE(0.5 * M_PI - aop);
    
    if (type == TDS_SECONDARY)
        ftra += M_PI;
    
    double M_1 = RADRANGE(ftra-2.*e*sin(ftra)+0.75*e*e*sin(2*ftra));
    t = M_1/n + tau;

    if (fabs(t - state->time - per) < fabs(t-state->time))
        t -= per;
    else if (fabs(t - state->time + per) < fabs(t-state->time))
        t += per;
    
    if (IS_INVALID(t) || IS_INVALID(state->epoch)) {
        *(timeout) = INVALID_NUMBER;
        return 0;
    }
        
    // Allocate working space
    double fout[3];
    
    gsl_vector* times = gsl_vector_alloc(1);
    
    ok_system** bag = NULL;
    VSET(times, 0, t);
    bag = ok_integrate(state, times, options, intMethod, bag, error);
    ok_system* buf = ok_copy_system(bag[0]);
    
    double diff = DBL_MAX;
    int steps = 0;
    
    while (diff > eps) {
        // Center on star
        ok_to_star(buf, buf->xyz);
        
        // Calculate cartesian accelerations
        if (intMethod != KEPLER)
            ok_force_on(buf->time, buf->xyz->data, fout, buf, pidx);
        else
            ok_star_force_on(buf->time, buf->xyz->data, fout, buf, pidx);
        
        double rs_1 = MGET(buf->xyz, pidx, X) * MGET(buf->xyz, pidx, VX) +
                MGET(buf->xyz, pidx, Y) * MGET(buf->xyz, pidx, VY);
        double rs_2 = MGET(buf->xyz, pidx, VX) * MGET(buf->xyz, pidx, VX) +
                MGET(buf->xyz, pidx, X) * fout[0] +
                MGET(buf->xyz, pidx, VY) * MGET(buf->xyz, pidx, VY) +
                MGET(buf->xyz, pidx, Y) * fout[1];
        
        double t_1 = t - rs_1/rs_2;
        diff = fabs(t-t_1);
        
        t = t_1;
        VSET(times, 0, t);
        bag = ok_integrate(buf, times, options, intMethod, bag, error);
        
        // Set the new system t
        ok_copy_system_to(bag[0], buf);
        
        
        if (steps > 30) {
            printf("breaking, %d %e %e %e %e\n", steps, t, sqrt(SQR(MGET(buf->xyz, pidx, X)) + SQR(MGET(buf->xyz, pidx, Y))), diff,
                    MGET(buf->orbits, pidx, SMA));
            retval = OK_NOCONV;
            t = INVALID_NUMBER;
            break;
        }
        
        steps++;
    };
    
    ok_free_system(buf);
    ok_free_system(bag[0]);
    ok_free_system(bag[1]);
    free(bag);
    gsl_vector_free(times);
    options->iterations = 2;
    *timeout = t;
    
    return retval;
}

gsl_vector* ok_find_transits(ok_system** bag, const int len, const int pidx, const int intMethod, const double eps, const int flag[], int* error) {
    gsl_vector* times = gsl_vector_alloc(len);
    ok_integrator_options o;
    memcpy(&o, &defoptions, sizeof(ok_integrator_options));
    o.calc_elements = false;
    
    for (int i = 0; i < len; i++) {
        ok_find_closest_transit(bag[i], pidx, &o, intMethod, eps, (flag == NULL ? TDS_PRIMARY : flag[i]), &(times->data[i]), error);
    }
    return times;
}
