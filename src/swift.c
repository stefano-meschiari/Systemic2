/// Work in progress

#include "swift.h"
#include "assert.h"
#include "systemic.h"
#include "utils.h"
#include "integration.h"
#include "stdint.h"
#include "unistd.h"

#define NTPMAX 1001
#define NPLMAX 51
#define NSTATP 3
#define NSTAT NSTATP + NPLMAX - 1

#define NSTATR NSTAT

extern void rmvs_step_systemic_(int*, double*, int*, int*, double*, double*, double*,
        double*, double*, double*, double*, double*, double*,
        double*, double*, double*, double*, double*, double*,
        int*, double*, double*, int*);

ok_system** ok_integrate_swift(ok_system* initial, const gsl_vector* times, const ok_integrator_options* options,
    ok_system** bag, int* error) {
    
    const double startTime = initial->epoch;
    int NDIMS = initial->nplanets + 1;
    int* swiftError = malloc(sizeof(int));
    swiftError[0] = 0;
    
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
    
    double** r = (double**) malloc(sizeof(double*) * 6);
    double** r1 = (double**) malloc(sizeof(double*) * 6);
    double* mass;
    
    for (int i = 0; i < 6; i++) {
        r[i] = (double*) malloc(sizeof(double) * NDIMS * 2);
        r1[i] = r[i] + NDIMS;
    }
    mass = (double*) malloc(sizeof(double) * NDIMS);
    
    gsl_matrix* xyz = initial->xyz;
    
    for (int i = 0; i < NDIMS; i++) {
        mass[i] = MGET(xyz, i, 0);
        
        for (int j = 0; j < 6; j++)
            r[j][i] = MGET(xyz, i, j+1) - MGET(xyz, 0, j+1);
    }
    
    
    double dum[0];
    int dum0 = 0;
    double dum0d = 0.;
    
    const double dt = options->dt;
    double h;
    
    ok_progress progress = options->progress;
    
    int step = 0;
    // Loop through the times vector
    
    for (int i = 0; i < SAMPLES; i++) {
        double time = times->data[i];
        
        for (int j = 0; j < 6; j++)
            memcpy(r1[j], r[j], sizeof(double) * NDIMS);
        
    
        // Integrate between prevTime and time
        
        if (fabs(time - prevTime) > 1e-10 && swiftError[0] == 0) {
            double t = prevTime;
            
            while (fabs(t - time) > 1e-10) {
                h = copysign(MIN(dt, fabs(t-time)), time - prevTime);
                
                rmvs_step_systemic_(&step, &t, &NDIMS, &dum0, mass, &dum0d, &dum0d, 
                        r[0], r[1], r[2], r[3], r[4], r[5], 
                        dum, dum, dum, dum, dum, dum,
                        NULL, NULL, &h, swiftError);
                if (swiftError[0] != 0)
                    break;
                t += h;
            }
        } 
        if (swiftError[0] != 0) {
            
            for (int j = 0; j < NDIMS; j++)
                for (int i = 0; i < 6; i++)
                       r1[i][j] = INVALID_NUMBER;
        }
        
        // Set up the return vector
        bag[i]->time = time;
        
        
        for (int j = 0; j < NDIMS; j++) {
            MSET(bag[i]->xyz, j, 0, mass[j]);
            MSET(bag[i]->xyz, j, 1, r1[0][j]);
            MSET(bag[i]->xyz, j, 2, r1[1][j]);
            MSET(bag[i]->xyz, j, 3, r1[2][j]);
            MSET(bag[i]->xyz, j, 4, r1[3][j]);
            MSET(bag[i]->xyz, j, 5, r1[4][j]);
            MSET(bag[i]->xyz, j, 6, r1[5][j]);
        }
        
        
        // Converts the current cartesian coordinates back to orbital elements
        if (options == NULL || options->calc_elements) {
                ok_cart2el(bag[i], bag[i]->orbits, true);
        }
        
        if (progress != NULL) {
            int ret = progress(i, SAMPLES, NULL, "");
            
            if (ret == PROGRESS_STOP) {
                for (int i = 0; i < SAMPLES; i++)
                        ok_free_system(bag[i]);
                
                free(bag);
                bag = NULL;
                free(r[0]);
                free(r[1]);
                free(r[2]);
                free(r[3]);
                free(r[4]);
                free(r[5]);
                free(r);
                free(r1);
                free(mass);   
                
                if (error != NULL) {
                    *error = INTEGRATION_FAILURE_STOPPED;
                }
                return NULL;
            }
        }
        
        prevTime = time;
    }
    
    if (swiftError[0] != 0) {
        *error = INTEGRATION_FAILURE_SWIFT;
    }
    
    free(r[0]);
    free(r[1]);
    free(r[2]);
    free(r[3]);
    free(r[4]);
    free(r[5]);
    free(r);
    free(r1);
    free(mass);   

    return bag;
}

