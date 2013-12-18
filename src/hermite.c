// THIS IS WORK IN PROGRESS -- DO NOT USE

#include "hermite.h"
#include "systemic.h"
#include "utils.h"
#include "assert.h"
#include "integration.h"

typedef int (* force_jerk) (double, const double[], double[], double[], void*);

#define SWAP(a, b) double* a_##a = a; a = b; b = a_##a;

int ok_hermite_step(const int nbody, const double prevTime, double toTime,
        double* xyz, double* xyzp, double* xyz1, double* f, double* jerk, double* f1, double* jerk1, const double tpar, double* h,
        void* par, const force_jerk force, const int ITERATIONS) {
    
    double t = prevTime;
    double dt = *h;
    
    while (fabs(toTime - t) > 1e-15) {
        
        dt = copysign(MIN(dt, fabs(toTime - t)), toTime - prevTime);
        double dt_prop;
        jerk1[nbody*3] = jerk[nbody*3];
        
        for (int k = 0; k < ITERATIONS; k++) {
            
            dt_prop = 0.5 * tpar * (sqrt(1./jerk[nbody*3]) + sqrt(1./jerk1[nbody*3]));
            dt = copysign(MIN(dt_prop, fabs(toTime-t)), toTime-prevTime);
            
            double dt_2 = dt*dt;
            double dt_3 = dt*dt_2;
            
            // PREDICTOR STEP
            for (int i = 0; i < nbody; i++) {
                for (int d = 1; d <= 3; d++) {
                        xyzp[i * 7 + d] = xyz[i * 7 + d] + xyz[i * 7 + d + 3] * dt 
                                + 0.5 * f[i * 7 + d + 3] * dt_2 
                                + (1./6.) * jerk[i * 3 + d - 1] * dt_3;

                        xyzp[i * 7 + d + 3] = xyz[i * 7 + d + 3] 
                                + f[i * 7 + d + 3] * dt 
                                + 0.5 * jerk[i * 3 + d - 1] * dt_2;
                }
            }
            
            (*force)(t, xyzp, f1, jerk1, par);
            /*
            // CORRECTOR STEP
            for (int i = 0; i < nbody; i++) {
                for (int d = 1; d <= 3; d++) {
                        xyz1[i * 7 + d] = xyzp[i * 7 + d]
                                - (0.15) * dt_2 * (f[i * 7 + d + 3] - f1[i * 7 + d + 3])
                                - (1./60.) * dt_3 * (7. * jerk[i * 3 + d - 1] + 2. * jerk1[i * 3 + d - 1]);

                        xyz1[i * 7 + d + 3] = xyzp[i * 7 + d + 3]
                                + (-.5) * dt * (f[i * 7 + d + 3] - f1[i * 7 + d + 3])
                                - 1./12. * dt_2 * (5. * jerk[i * 3 + d - 1] + jerk1[i * 3 + d - 1]);
                }
            }

            (*force)(t, xyz1, f1, jerk1, par);

            // CORRECTOR STEP
            if (k == ITERATIONS - 1)
                for (int i = 0; i < nbody; i++) {
                    for (int d = 1; d <= 3; d++) {
                            xyz[i * 7 + d] = xyzp[i * 7 + d]
                                    - (.15) * dt_2 * (f[i * 7 + d + 3] - f1[i * 7 + d + 3])
                                    - (1./60.) * dt_3 * (7. * jerk[i * 3 + d - 1] + 2. * jerk1[i * 3 + d - 1]);

                            xyz[i * 7 + d + 3] = xyzp[i * 7 + d + 3]
                                    - 0.5 * dt * (f[i * 7 + d + 3] - f1[i * 7 + d + 3])
                                    - 1./12. * dt_2 * (5. * jerk[i * 3 + d - 1] + jerk1[i * 3 + d - 1]);
                    }
                }
            else
                for (int i = 0; i < nbody; i++) {
                    for (int d = 1; d <= 3; d++) {
                            xyz1[i * 7 + d] = xyzp[i * 7 + d]
                                    - (.15) * dt_2 * (f[i * 7 + d + 3] - f1[i * 7 + d + 3])
                                    - (1./60.) * dt_3 * (7. * jerk[i * 3 + d - 1] + 2. * jerk1[i * 3 + d - 1]);

                            xyz1[i * 7 + d + 3] = xyzp[i * 7 + d + 3]
                                    - 0.5 * dt * (f[i * 7 + d + 3] - f1[i * 7 + d + 3])
                                    - 1./12. * dt_2 * (5. * jerk[i * 3 + d - 1] + jerk1[i * 3 + d - 1]);
                    }
                }**/
            
            for (int i = 0; i < nbody; i++) {
                for (int d = 1; d<=3; d++) {
                    double a2 = -6.*(f[i * 7 + d + 3] - f1[i*7 + d + 3]) -
                        (4.*jerk[i * 3 + d - 1] + 2.*jerk1[i * 3 + d - 1]) * dt;
                    double a3 = 12.*(f[i*7 + d + 3] - f1[i * 7 + d + 3]) +
                        6.*(jerk[i*3+d-1] + jerk1[i*3+d-1]) * dt;
                    
                    xyz[i * 7 + d] = xyzp[i*7 + d] + (a2/24.+a3/120.) * dt_2;
                    xyz[i * 7 + d + 3] = xyzp[i*7 + d + 3] + (a2/6.+a3/24.) * dt;   
                }
            }
            
            break;
        }
        
        t += dt;
        (*force)(t, xyz, f, jerk, par); 
    }
    
    *h = dt;
    return INTEGRATION_SUCCESS;
}

ok_system** ok_integrate_hermite(ok_system* initial, const gsl_vector* times, ok_integrator_options* options,
        ok_system** bag, int* error) {
    
    // Check that the system has been set-up
    assert(initial->xyz != NULL);
    // Check input arguments
    assert(times != NULL);
    
    const double startTime = initial->epoch;
    
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
            bag[i]->xyz = (bag[i]->xyz != NULL ? bag[i]->xyz : gsl_matrix_alloc(initial->nplanets+1, 7));
        }
    }
    
    
    double prevTime = startTime;
    gsl_matrix* prevOrbits = initial->orbits;
    
    const int nbody = initial->nplanets + 1;
    const int BUFLENGTH = 3*(nbody * 7) +
                        2*(nbody * 7) +
                        2*(nbody * 3 + 1);
    
    if (options->buffer == NULL || options->buffer->size != BUFLENGTH) {
        if (options->buffer != NULL) {
            gsl_vector_free(options->buffer);
        }
        options->buffer = gsl_vector_alloc(BUFLENGTH);
    }
    
    double* xyz = options->buffer->data;
    MATRIX_MEMCPY_TOARRAY(xyz, initial->xyz);
    
    double* xyz_scratch = xyz + 7*nbody;
    MATRIX_MEMCPY_TOARRAY(xyz_scratch, initial->xyz);
    
    double* xyz_scratch2 = xyz_scratch + 7*nbody;
    MATRIX_MEMCPY_TOARRAY(xyz_scratch2, initial->xyz);
    
    double* jerk = xyz_scratch2 + 7 * nbody;
    double* f = jerk + 3 * nbody + 1;
    double* jerk1 = f + 7 * nbody;
    double* f1 = jerk1 + 3 * nbody;
    
    void* params = (void*) initial;
    double tpar = options->acc_par;
    
    options->force_jerk(startTime, xyz, f, jerk, params);
    
    double h = tpar * sqrt(1./jerk[nbody*3]);
    
    // Loop through the times vector
    for (int i = 0; i < SAMPLES; i++) {
        double time = times->data[i];
        

        // Integrate between prevTime and time
        if (fabs(time - prevTime) > 1e-10) {
            h = copysign(h, time-prevTime);

            int err = ok_hermite_step(nbody, prevTime, time, xyz, xyz_scratch, xyz_scratch2,
                    f, jerk, f1, jerk1, tpar, &h, params,
                    options->force_jerk, options->iterations);

            if (err != INTEGRATION_SUCCESS) {
                for (int i = 0; i < SAMPLES; i++)
                    ok_free_system(bag[i]);
                free(bag);  
                return NULL;
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
        
        prevTime = time;
        params = bag[i];
    }
      
    return bag;
}
