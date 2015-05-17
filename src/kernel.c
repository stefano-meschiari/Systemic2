/*
 *  kernel.c
 *  Systemic2
 *
 *
 */

#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>
#ifndef JAVASCRIPT
#include "omp.h"
#else
#include "omp_shim.h"
#endif

#include <gsl/gsl_randist.h>
#include "mercury.h"
#include "kernel.h"
#include "integration.h"
#include "utils.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_rng.h"
#include "simplex.h"
#include "lm.h"
#include "sa.h"
#include "de.h"
#include "time.h"
#include <libgen.h>


ok_minimizer ok_minimizers[] = {K_minimize_simplex, K_minimize_lm, K_minimize_de, K_minimize_sa, NULL, NULL, NULL, NULL};
char * ok_orb_labels[ELEMENTS_SIZE] = {"P", "M", "MA", "E", "LOP", "I", "NODE", "RADIUS", "ORD",
    "UNUSED1_", "UNUSED2_", "UNUSED3_", "UNUSED4_"};
char * ok_all_orb_labels[ALL_ELEMENTS_SIZE] = {"P", "M", "MA", "E", "LOP", "I", "NODE", "RADIUS", "ORD",
    "UNUSED1_", "UNUSED2_", "UNUSED3_", "UNUSED4_", "SMA", "K", "TPERI", "TRUEANOM", "UNUSED5_", "UNUSED6_", "UNUSED7_"};


double ok_default_steps[ELEMENTS_SIZE] = {[PER] = 1e-3, [MASS] = 1e-3, [MA] = 1e-2,
    [ECC] = 1e-2, [LOP] = 1e-2, [INC] = 1e-2, [NODE] = 1e-2, [RADIUS] = 1e-2, 0, [PRECESSION_RATE] = 1e-4, 0, 0};

/**
 * An internal function used to make sure the internal state of the kernel is
 * consistent. This is used, for instance, to:
 * - ensure the various matrices have the right sizes (e.g., the number of rows of 
 * the minimization flags matrix and the number of rows of the elements have to be kept in
 * sync)
 * - sort planets, if needed (e.g. if the coordinate system is JACOBI, the bodies have
 * to be sorted by increasing periods)
 * - setting up the system member struct (which contains orbital elements, cartesian
 * coordinates, epoch)
 * - other bookkeeping.
 * 
 * It is usually called by the other functions (e.g. K_setElement) to ensure the
 * kernel is in a consistent state.
 * 
 * This function should not be used directly by the user or a binding, but it could
 * be useful when debugging.
 * 
 * @param k The kernel to check.
 */
void K_validate(ok_kernel* k);

/**
 * Allocates a new kernel with 0 planets and an empty dataset.
 * @return A new kernel object.
 */
ok_kernel* K_alloc() {
    ok_kernel* k = (ok_kernel*) malloc(sizeof (ok_kernel));
    memset(k, 0, sizeof (ok_kernel));
    k->chi2 = INVALID_NUMBER;
    k->chi2_rvs = INVALID_NUMBER;
    k->chi2_tts = INVALID_NUMBER;
    k->chi2_other = INVALID_NUMBER;
    k->rms_tts = INVALID_NUMBER;
    k->jitter = INVALID_NUMBER;
    k->rms = INVALID_NUMBER;
    k->nsets = 0;
    k->ndata = 0;
    k->nrvs = 0;
    k->ntts = 0;
    k->last_error = 0;
    k->intMethod = KEPLER;
    k->Mstar = 1.;
    k->minfunc = K_getChi2;
    k->system = ok_alloc_system(0);
    k->info = NULL;
    
    MSET(k->system->elements, 0, MASS, 1.);
    k->compiled = NULL;
    k->plRanges = (gsl_matrix**) calloc(2, sizeof (gsl_matrix*));
    k->parRanges = (gsl_vector**) calloc(2, sizeof (gsl_vector*));

    k->times = NULL;
    k->datasets = (gsl_matrix**) calloc(PARAMS_SIZE, sizeof (gsl_matrix*));
    k->params = gsl_vector_calloc(PARAMS_SIZE);

    k->integration = NULL;
    k->plFlags = NULL;
    k->parFlags = gsl_vector_int_calloc(PARAMS_SIZE);
    k->parSteps = gsl_vector_calloc(PARAMS_SIZE);
    for (int i = 0; i < PARAMS_SIZE; i++)
        VSET(k->parSteps, i, 1e-2);
    for (int i = 0; i <= P_DATA_NOISE10; i++)
        VSET(k->parSteps, i, 5e-2);
    VSET(k->parSteps, P_RV_TREND_QUADRATIC, 1e-4);
    k->intOptions = (ok_integrator_options*) malloc(sizeof (ok_integrator_options));
    memcpy(k->intOptions, &defoptions, sizeof (ok_integrator_options));
    k->intOptions->buffer = NULL;
    k->intOptions->ibuffer = NULL;
    k->intOptions->progress = NULL;

    k->flags = 0;

    k->progress = NULL;
    k->model_function = NULL;

    k->rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(k->rng, clock());
    return k;
}

/**
 * Frees a kernel object allocated with K_alloc.
 * @param k The kernel object to be freed.
 */
void K_free(ok_kernel* k) {
    if (k == NULL)
        return;

    if (!(k->flags & SHARE_FLAGS)) {
        gsl_matrix_int_free(k->plFlags);
        gsl_vector_int_free(k->parFlags);
    }
    if (!(k->flags & SHARE_STEPS)) {
        gsl_matrix_free(k->plSteps);
        gsl_vector_free(k->parSteps);
    }
    if (!(k->flags & SHARE_RANGES)) {
        gsl_matrix_free(k->plRanges[0]);
        gsl_matrix_free(k->plRanges[1]);
        gsl_vector_free(k->parRanges[0]);
        gsl_vector_free(k->parRanges[1]);
        free(k->parRanges);
        free(k->plRanges);
    }

    if (!(k->flags & SHARE_DATA)) {
        for (int i = 0; i < k->nsets; i++) {
            gsl_matrix_free(k->datasets[i]);
        }
        free(k->datasets);

        gsl_vector_free(k->times);

        if (k->integration != NULL) {
            for (int i = 0; i < k->integrationSamples; i++)
                ok_free_system(k->integration[i]);
            free(k->integration);
        }
        free(k->compiled);
        
        K_clearInfo(k);
            
    }
    if (k->intOptions->buffer != NULL)
        gsl_vector_free(k->intOptions->buffer);
    if (k->intOptions->ibuffer != NULL)
        gsl_vector_int_free(k->intOptions->ibuffer);

    assert(!(k->system->flag & FREED));
    k->system->flag = FREED;
    ok_free_system(k->system);
    gsl_vector_free(k->params);
    gsl_rng_free(k->rng);
    k->flags = FREED;
    free(k->intOptions);
    free(k);
}

int K_compile_sort(const void* elem1, const void* elem2) {
    double* e1 = *(double**) elem1;
    double* e2 = *(double**) elem2;

    return (e1[T_TIME] - e2[T_TIME] < 0 ? -1 : 1);
}

/**
 * Compiles the data by merging all datasets and sorting them by date; a pointer table
 * is returned, each entry pointing to a row of a dataset. Data from different datasets
 * is mixed.
 * @param k The kernel containing the datasets.
 * @return A (double**) pointer table, containing all the data points sorted by date
 */
double** K_compileData(ok_kernel* k) {
    if (k->nsets < 1) {
        k->ndata = 0;
        return NULL;
    }
    if (k->times != NULL)
        gsl_vector_free(k->times);

    int ndata = 0;
    int ncols = MCOLS(k->datasets[0]);

    k->nrvs = 0;
    k->ntts = 0;

    for (int i = 0; i < k->nsets; i++) {
        ndata += MROWS(k->datasets[i]);
        ncols = MIN(ncols, MCOLS(k->datasets[i]));

        if (MROWS(k->datasets[i]) > 0) {

            if ((int) MGET(k->datasets[i], 0, T_FLAG) == T_RV) {
                k->nrvs += MROWS(k->datasets[i]);
            } else if ((int) MGET(k->datasets[i], 0, T_FLAG) == T_TIMING)
                k->ntts += MROWS(k->datasets[i]);
        }
    }

    ncols++;

    double** compiled = k->compiled;
    if (k->compiled == NULL || k->ndata != ndata) {
        compiled = (double**) malloc(sizeof (double*) * ndata);
        free(k->compiled);
    }

    int nr = 0;

    for (int i = 0; i < k->nsets; i++) {
        int rows = MROWS(k->datasets[i]);
        for (int j = 0; j < rows; j++) {
            compiled[nr] = gsl_matrix_ptr(k->datasets[i], j, 0);
            compiled[nr][T_SET] = (double) i;
            nr++;
        }
    }

    if (k->flags & BOOTSTRAP_DATA) {
        k->compiled = (double**) malloc(sizeof (double*) * ndata);
        for (int i = 0; i < ndata; i++) {
            k->compiled[i] = compiled[gsl_rng_uniform_int(k->rng, ndata)];
        }
        free(compiled);
        compiled = k->compiled;
        k->flags &= ~BOOTSTRAP_DATA;
    } else
        k->compiled = compiled;

    qsort(compiled, ndata, sizeof (double*), K_compile_sort);

    k->times = gsl_vector_calloc(ndata);
    for (int i = 0; i < ndata; i++) {
        VSET(k->times, i, k->compiled[i][T_TIME]);
    }

    if (IS_INVALID(k->system->epoch)) {
        k->system->epoch = k->system->time = VGET(k->times, 0);
    }
    k->ndata = ndata;

    k->flags &= ~NEEDS_COMPILE;
    return k->compiled;
}

/**
 * Adds a new data table (in the form of a matrix) to the list of datasets used 
 * by the kernel.
 * @param k Kernel where the data table will be added
 * @param table A matrix containing several columns (depending on the type of data);
 * common columns are TIME (=0), VALUE (=1), SIGMA (=2). For RVs, for instance, this would
 * correspond to the time of observation, the measured radial velocity and the
 * error of the measurement. Note that the matrix is owned by the kernel after this call,
 * and might be freed or resized as needed.
 * @param name Name describing the dataset (e.g. "AlphaCen_KECK").
 * @param type Type of dataset. Can be one of T_RV, T_PHOTO (photometry) or T_TIMING
 * (a TTV dataset)
 * @return The matrix itself.
 */
gsl_matrix* K_addDataTable(ok_kernel* k, gsl_matrix* table, const char* name, int type) {
    if (MCOLS(table) < DATA_SIZE) {
        table = ok_matrix_resize(table, MROWS(table), DATA_SIZE);
    }

    k->datasets[k->nsets] = table;
    
    for (int i = 0; i < MROWS(table); i++)
        MSET(table, i, T_FLAG, type);

    k->nsets++;
    k->flags |= NEEDS_COMPILE;
    if (type == T_RV)
        VISET(k->parFlags, k->nsets - 1, ACTIVE | MINIMIZE);
    else
        VISET(k->parFlags, k->nsets - 1, ACTIVE);
    
    char dataTag[100];
    sprintf(dataTag, "DataFileName%d", k->nsets-1);
    K_setInfo(k, dataTag, name);
    
    return table;
}

/**
 * Adds a new data table, read from the file specified by "path". 
 * @param k The kernel to receive the dataset
 * @param path Path to an ASCII file. The ASCII file should contain the data
 * separated by newlines; each column is separated by whitespace. Lines starting
 * with # are treated as comments.
 * @param type Type of dataset. Can be one of T_RV, T_PHOTO (photometry) or T_TIMING
 * (a TTV dataset)
 * @return A matrix containing the data read from the file. This matrix is owned
 * by the kernel, so make a copy if you need to hold on to it.
 */
gsl_matrix* K_addDataFile(ok_kernel* k, const char* path, int type) {
    FILE* fid = fopen(path, "r");
    assert(fid != NULL);
    if (fid == NULL)
        return NULL;


    k->datasets[k->nsets] = ok_read_table(fid);

    for (int i = 0; i < MROWS(k->datasets[k->nsets]); i++)
        MSET(k->datasets[k->nsets], i, T_FLAG, type);

    if (type == T_TIMING) {
        for (int i = 0; i < MROWS(k->datasets[k->nsets]); i++) {
            MSET(k->datasets[k->nsets], i, T_TDS_FLAG, MGET(k->datasets[k->nsets], i,
                    2));
            MSET(k->datasets[k->nsets], i, T_ERR, MGET(k->datasets[k->nsets], i,
                    T_VAL));
            MSET(k->datasets[k->nsets], i, T_VAL, MGET(k->datasets[k->nsets], i,
                    T_TIME));
        }
    }

    k->nsets++;
    fclose(fid);

    k->flags |= NEEDS_COMPILE;
    if (type == T_RV)
        VISET(k->parFlags, k->nsets - 1, ACTIVE | MINIMIZE);
    else
        VISET(k->parFlags, k->nsets - 1, ACTIVE);

    VSET(k->parSteps, k->nsets - 1, 0.05);
    VSET(k->params, k->nsets - 1, 0.);
    
    char dataTag[100];
    sprintf(dataTag, "DataFileName%d", k->nsets-1);
    K_setInfo(k, dataTag, path);
    
    return k->datasets[k->nsets - 1];
}

/**
 * Removes the idx-th dataset (shuffling datasets and parameter values as needed)
 * @param k The kernel containing the dataset
 * @param idx Index of the dataset
 */
void K_removeData(ok_kernel* k, int idx) {

    if (idx == -1) {
        for (int i = 0; i < k->nsets; i++) {
            gsl_matrix_free(k->datasets[i]);
            VSET(k->params, i, 0.);
            k->parFlags->data[i] = 0.;
        }
        k->nsets = 0;
        k->ndata = 0;
        k->nrvs = 0;
        k->ntts = 0;
        k->flags |= NEEDS_COMPILE | NEEDS_SETUP;
        return;
    }

    gsl_matrix_free(k->datasets[idx]);
    k->datasets[idx] = NULL;
    
    char a[100];
    
    for (int i = idx + 1; i < k->nsets; i++) {
        k->datasets[i - 1] = k->datasets[i];
        k->params->data[i - 1] = k->params->data[i];
        k->parFlags->data[i - 1] = k->parFlags->data[i];
        k->parSteps->data[i - 1] = k->parSteps->data[i];
        sprintf(a, "DataFileName%d", i);
        
        char* fn = K_getInfo(k, a);
        sprintf(a, "DataFileName%d", i-1);
        K_setInfo(k, a, fn);
    }
    sprintf(a, "DataFileName%d", k->nsets-1);
    K_setInfo(k, a, NULL);
    
    k->parFlags->data[k->nsets - 1] = 0;
    k->parSteps->data[k->nsets - 1] = 1e-2;

    k->nsets--;
    k->flags |= NEEDS_COMPILE | NEEDS_SETUP;
    K_validate(k);
}

gsl_matrix* K_getData(ok_kernel* k, int idx) {
    return k->datasets[idx];
}

int K_getDataType(ok_kernel* k, int idx) {
    if (MROWS(k->datasets[idx]) == 0)
        return T_RV;
    else
        return MGET(k->datasets[idx], 0, T_FLAG);
}

int K_getDataSize(ok_kernel* k, int idx) {
    return MROWS(k->datasets[idx]);
}

gsl_matrix* K_getCompiledDataMatrix(ok_kernel* k) {
    if (k->compiled == NULL)
        K_compileData(k);

    gsl_matrix* res = gsl_matrix_alloc(k->ndata, DATA_SIZE);

    for (int i = 0; i < k->ndata; i++) {
        for (int j = 0; j < DATA_SIZE; j++)
            MSET(res, i, j, k->compiled[i][j]);
    }
    return res;
}

void K_setData(ok_kernel* k, int idx, gsl_matrix* data) {
    data = ok_matrix_resize(data, MROWS(data), DATA_SIZE);
    gsl_matrix_free(k->datasets[idx]);
    k->datasets[idx] = data;
    k->flags |= NEEDS_COMPILE | NEEDS_SETUP;
    K_validate(k);
}

void K_addPlanet(ok_kernel* k, const double elements[]) {

    k->system->nplanets++;
    int row = k->system->nplanets;

    ok_resize_system(k->system, k->system->nplanets + 1);

    k->plFlags = ok_matrix_int_resize(k->plFlags, k->system->nplanets + 1, MCOLS(k->system->elements));
    k->plSteps = ok_matrix_resize(k->plSteps, k->system->nplanets + 1, MCOLS(k->system->elements));

    MSET(k->system->elements, row, ORD, -1);
    MSET(k->system->elements, row, INC, 90.);

    if (elements != NULL) {
        int i = 0;
        while (true) {
            if (elements[i] == DONE)
                break;
            else
                MSET(k->system->elements, row, (int) elements[i], elements[i + 1]);
            i += 2;
        }
    }

    if (MGET(k->system->elements, row, ORD) < 0)
        MSET(k->system->elements, row, ORD, row);


    for (int i = 0; i < 5; i++)
        MISET(k->plFlags, row, i, ACTIVE | MINIMIZE);
    for (int i = 0; i < ELEMENTS_SIZE; i++) {
        MSET(k->plSteps, row, i, ok_default_steps[i]);
    }

    K_setElementRange(k, row, PER, 0.1, 1e4);
    K_setElementRange(k, row, MASS, 1e-4, 2e3);
    K_setElementRange(k, row, ECC, 0, 0.99);

    ok_sort_matrix(k->system->elements, ORD);
    k->flags |= NEEDS_SETUP;
    K_validate(k);
}

/**
 * Removes the idx-th body from the system (the 0-th body representing the central star).
 * Matrices are shuffled automatically to accomodate the deletion of the body.
 * @param k Kernel to be modified
 * @param idx Index of the body to be deleted
 */

void K_removePlanet(ok_kernel* k, int idx) {
    if (k->system->nplanets == 0)
        return;

    if (idx == -1) {
        while (k->system->nplanets > 0) {
            K_removePlanet(k, 0);
        }
        k->flags |= NEEDS_SETUP;
        return;
    }



    gsl_matrix_free(k->system->xyz);
    k->system->xyz = NULL;
    gsl_matrix_free(k->system->orbits);
    k->system->orbits = NULL;


    gsl_matrix* els = ok_matrix_remove_row(k->system->elements, idx);
    gsl_matrix_free(k->system->elements);
    k->system->elements = els;

    gsl_matrix_int* plFlags = ok_matrix_int_remove_row(k->plFlags, idx);
    gsl_matrix_int_free(k->plFlags);
    k->plFlags = plFlags;

    if (k->plRanges[0] != NULL) {
        gsl_matrix* plRanges0 = ok_matrix_remove_row(k->plRanges[0], idx);
        gsl_matrix_free(k->plRanges[0]);
        k->plRanges[0] = plRanges0;
        gsl_matrix* plRanges1 = ok_matrix_remove_row(k->plRanges[1], idx);
        gsl_matrix_free(k->plRanges[1]);
        k->plRanges[1] = plRanges1;
    }

    for (int i = 0; i < k->system->elements->size1; i++)
        MSET(k->system->elements, i, ORD, i);

    k->system->nplanets--;
    k->flags |= NEEDS_SETUP;
}

/**
 * Sets the element specified by col (one of PER, MASS, ECC, MA,
 * LOP, INC, NODE, RADIUS, ECCANOMALY, TRUEANOMALY) to the specified value for the row-th body.
 * @param k Kernel to be modified
 * @param row Index of the planet (0-th body is the central star)
 * @param col Index of the element
 * @param value New value
 */
void K_setElement(ok_kernel* k, int row, int col, double value) {

    if (row == -1) {
        for (int i = 1; i < MROWS(k->system->elements); i++)
            K_setElement(k, i, col, value);
    } else {
        if (col >= 0 && col < ELEMENTS_SIZE)
            MSET(k->system->elements, row, col, value);
        else if (col == TPERI) {
            double M = 360. / MGET(k->system->elements, row, PER) * (k->system->epoch - value);
            MSET(k->system->elements, row, MA, DEGRANGE(M));
        } else if (col == TRUEANOMALY) {
            double e = MGET(k->system->elements, row, ECC);
            double f = RADRANGE(TO_RAD(value));

            double E = 2 * atan2(sqrt(1 - e) * sin(0.5 * f),
                    sqrt(1 + e) * cos(0.5 * f));

            double M = E - e * sin(E);
            MSET(k->system->elements, row, MA, TO_DEG(M));

        } else if (col == SEMIAMP) {
            double Mcent = MSUN_TO_INT(k->Mstar);
            double per = MGET(k->system->elements, row, PER);
            double ecc = MGET(k->system->elements, row, ECC);

            if (k->system->flag & JACOBI) {
                for (int j = 1; j < k->system->nplanets + 1; j++)
                    if (MGET(k->system->elements, j, PER) < per)
                        Mcent += MJUP_TO_INT(MGET(k->system->elements, j, MASS));
            }
            double K = MPS_TO_AUPDAY(value);
            double xi = K * cbrt(per / (2. * M_PI)) * sqrt(1. - ecc * ecc);

            double m = xi * pow(Mcent, 2. / 3.);

            int steps = 0;
            while (abs(1 - xi * pow(Mcent + m, 2. / 3.) / m) > 1e-6) {
                m = xi * pow(Mcent + m, 2. / 3.);

                steps++;
                if (steps > 20)
                    break;
            }

            MSET(k->system->elements, row, MASS, INT_TO_MJUP(m));

        } else if (col == MEANLONGITUDE) {
            MSET(k->system->elements, row, MA, value - MGET(k->system->elements, row, LOP));
        } else if (col == SMA) {

            double mass = MGET(k->system->elements, row, MASS);
            double a = value;
            double per = MGET(k->system->elements, row, PER);
            double dper;
            double Mcent;
            do {
                Mcent = k->Mstar * MSUN;

                if (k->system->flag & JACOBI) {
                    for (int j = 1; j < k->system->nplanets + 1; j++)
                        if (MGET(k->system->elements, j, PER) < per)
                            Mcent += MGET(k->system->elements, j, MASS) * MJUP;
                } else
                    break;

                double per_new = YEAR_TO_DAY(sqrt(a * a * a * MSUN / (Mcent + mass * MJUP)));
                dper = per_new - per;
                per = per_new;

            } while (abs(dper) > 1e-10);

            per = YEAR_TO_DAY(sqrt(a * a * a * MSUN / (Mcent + mass * MJUP)));
            MSET(k->system->elements, row, PER, per);
        }
    }

    k->flags |= NEEDS_SETUP;
    K_validate(k);
}

/**
 * Returns the element specified by col (one of PER, MASS, ECC, MA,
 * LOP, INC, NODE, RADIUS, ECCANOMALY, TRUEANOMALY) for the row-th body.
 * @param k Kernel to be modified
 * @param row Index of the planet (0-th body is the central star)
 * @param col Index of the element
 * @return Value of the element
 */
double K_getElement(ok_kernel* k, int row, int col) {
    K_validate(k);
    if (col >= 0 && col < ELEMENTS_SIZE) {
        return MGET(k->system->elements, row, col);
    } else if (col == SMA) {
        double Mcent = k->Mstar * MSUN;
        double per = MGET(k->system->elements, row, PER);
        double mass = MGET(k->system->elements, row, MASS);

        if (k->system->flag & JACOBI) {
            for (int j = 1; j < k->system->nplanets + 1; j++)
                if (MGET(k->system->elements, j, PER) < per)
                    Mcent += MGET(k->system->elements, j, MASS) * MJUP;
        }

        double a = ok_acalc(per, Mcent / MSUN * K2, mass * MJUP / MSUN * K2);
        return a;
    } else if (col == SEMIAMP) {
        double Mcent = k->Mstar * MSUN;
        double per = MGET(k->system->elements, row, PER);
        double mass = MGET(k->system->elements, row, MASS);
        double ecc = MGET(k->system->elements, row, ECC);

        if (k->system->flag & JACOBI) {
            for (int j = 1; j < k->system->nplanets + 1; j++)
                if (MGET(k->system->elements, j, PER) < per)
                    Mcent += MGET(k->system->elements, j, MASS) * MJUP;
        }

        double a = ok_acalc(per, Mcent / MSUN * K2, mass * MJUP / MSUN * K2);

        double K = 2 * M_PI / per * mass * MJUP / (Mcent + mass * MJUP) * a / sqrt(1 - ecc * ecc);
        return K * AU / (100. * DAY);
    } else if (col == TPERI) {
        double n = 360./MGET(k->system->elements, row, PER);
        return k->system->epoch - MGET(k->system->elements, row, MA) / n;
    } else if (col == TRUEANOMALY) {
        double e = MGET(k->system->elements, row, ECC);
        double ma = TO_RAD(MGET(k->system->elements, row, MA));
        double E = mco_kep__(e, ma);

        double f = 2. * atan2(
                sqrt(1 + e) * sin(0.5 * E),
                sqrt(1 - e) * cos(0.5 * E));

        return TO_DEG(f);
    } else if (col == MEANLONGITUDE) {
        return DEGRANGE(MGET(k->system->elements, row, MA) + MGET(k->system->elements, row, LOP));
    }
    return INVALID_NUMBER;
}

/**
 * Returns the whole matrix of orbital elements. This is still owned by the kernel;
 * therefore, it is best to make a copy before the kernel is modified.
 * @param k Kernel
 * @return A matrix of the elements used internally (PER, MASS, MA, ECC, LOP, INC,
 * NODE, RADIUS, plus other reserved values)
 */
gsl_matrix* K_getElements(ok_kernel* k) {
    K_validate(k);
    return k->system->elements;
}

/**
 * Sets the orbital element matrix to the specified matrix. The matrix becomes
 * owned by the kernel.
 * @param k Kernel
 * @param elements A matrix of size nbody x ELEMENTS_SIZE (where nbody includes
 * the planets and the central star)
 */
void K_setElements(ok_kernel* k, gsl_matrix* elements) {
    if (elements != k->system->elements)
        gsl_matrix_free(k->system->elements);

    k->system->nplanets = MROWS(elements) - 1;
    k->system->elements = elements;
    k->flags |= NEEDS_SETUP;
    K_validate(k);
}

/**
 * Returns a matrix of the internal elements (like getElements), plus other derived
 * elements such as the semimajor axis (SMA), semiamplitude (K), true and eccentric
 * anomaly (TRUEANOMALY, ECCANOMALY) and others.
 * @param k Kernel
 * @return A matrix containing the orbital elements, augmented by additional
 * derived elements
 */
gsl_matrix* K_getAllElements(ok_kernel* k) {
    K_validate(k);

    gsl_matrix* ret = ok_matrix_copy(k->system->elements);
    ret = ok_matrix_resize(ret, MROWS(k->system->elements), ALL_ELEMENTS_SIZE);

    for (int i = 0; i < MROWS(ret); i++)
        for (int j = ELEMENTS_SIZE; j < ALL_ELEMENTS_SIZE; j++) {
            MSET(ret, i, j, K_getElement(k, i, j));
        }

    return ret;
}

/**
 * Sets the parameters of the kernel (the meaning of "parameter" depends on the index
 * of the parameter). Indices 0-6 are reserved for the datasets associated with the kernel;
 * for RV datasets, a parameter specifies a RV offset. Index TREND specifies a trend
 * subtracted from the compiled RV dataset (not yet implemented). 
 * Index NOISE specifies an amount of jitter (currently used internally by the MCMC routine).
 * @param k Kernel
 * @param pars A vector of size DATA_SIZE (which is subsequently owned by the kernel)
 */
void K_setPars(ok_kernel* k, gsl_vector* pars) {
    if (pars != k->params)
        gsl_vector_free(k->params);
    k->params = ok_vector_resize(pars, PARAMS_SIZE);
}

/**
 * Returns a vector of parameters for the kernel (see K_setPars).
 * @param k
 * @return A vector of parameters
 */
gsl_vector* K_getPars(ok_kernel* k) {
    return k->params;
}

/**
 * Sets the idx-th parameter to the specified value.
 * @param k Kernel
 * @param idx Index
 * @param val The specified value
 */

void K_setPar(ok_kernel* k, int idx, double val) {
    if (idx == -1) {
        for (int i = 0; i < k->params->size; i++)
            VSET(k->params, i, val);
    } else {
        VSET(k->params, idx, val);
    }
}

/**
 * Sets the idx-th parameter to the specified value.
 * @param k Kernel
 * @param idx Index
 * @return The current value
 */
double K_getPar(ok_kernel* k, int idx) {
    return VGET(k->params, idx);
}

void K_setParFlag(ok_kernel* k, int idx, int value) {
    VSET(k->parFlags, idx, value);
}

int K_getParFlag(ok_kernel* k, int idx) {
    return VGET(k->parFlags, idx);
}

double K_getParStep(ok_kernel* k, int idx) {
    return VGET(k->parSteps, idx);
}

void K_setParStep(ok_kernel* k, int idx, double value) {
    VSET(k->parSteps, idx, value);
}

void K_setParRange(ok_kernel* k, int par, double min, double max) {
    if (k->parRanges[0] == NULL) {
        k->parRanges[0] = gsl_vector_alloc(PARAMS_SIZE);
        gsl_vector_set_all(k->parRanges[0], INVALID_NUMBER);
    }
    if (k->parRanges[1] == NULL) {
        k->parRanges[1] = gsl_vector_alloc(PARAMS_SIZE);
        gsl_vector_set_all(k->parRanges[1], INVALID_NUMBER);
    }
    K_validate(k);

    VSET(k->parRanges[0], par, min);
    VSET(k->parRanges[1], par, max);
}

void K_getParRange(ok_kernel* k, int idx, double* min, double* max) {
    K_validate(k);

    if (k->parRanges[0] == NULL)
        *min = INVALID_NUMBER;
    else
        *min = VGET(k->parRanges[0], idx);
    if (k->parRanges[1] == NULL)
        *max = INVALID_NUMBER;
    else
        *max = VGET(k->parRanges[1], idx);
}

void K_validate(ok_kernel* k) {

    for (int i = k->nsets; i < DATA_SETS_SIZE; i++)
        VSET(k->parFlags, i, 0);


    if (k->flags & NEEDS_COMPILE) {
        K_compileData(k);
    }


    if (k->plFlags == NULL || (MROWS(k->plFlags) != MROWS(k->system->elements))) {
        k->plFlags = ok_matrix_int_resize(k->plFlags, MROWS(k->system->elements), MCOLS(k->system->elements));
    }

    if ((k->plSteps == NULL) || MROWS(k->plSteps) != MROWS(k->system->elements)) {
        k->plSteps = ok_matrix_resize(k->plSteps, MROWS(k->system->elements), MCOLS(k->system->elements));
    }


    if (k->plRanges[0] != NULL && MROWS(k->plRanges[0]) != MROWS(k->system->elements)) {
        k->plRanges[0] = ok_matrix_resize_pad(k->plRanges[0], MROWS(k->system->elements), ELEMENTS_SIZE, INVALID_NUMBER);
        k->plRanges[1] = ok_matrix_resize_pad(k->plRanges[1], MROWS(k->system->elements), ELEMENTS_SIZE, INVALID_NUMBER);
    }

    if (k->parRanges[0] != NULL)
        for (int i = 0; i < PARAMS_SIZE; i++) {
            if (!IS_INVALID(VGET(k->parRanges[0], i)))
                VSET(k->params, i, MAX(VGET(k->parRanges[0], i), VGET(k->params, i)));
            if (!IS_INVALID(VGET(k->parRanges[1], i)))
                VSET(k->params, i, MIN(VGET(k->parRanges[1], i), VGET(k->params, i)));
        }


    for (int i = 1; i < k->system->nplanets + 1; i++) {

        double ecc = MGET(k->system->elements, i, ECC);
        if (ecc < 0) {
            MSET(k->system->elements, i, ECC, -ecc);
            MINC(k->system->elements, i, MA, -180.);
            MINC(k->system->elements, i, LOP, +180.);
        }

        double mass = MGET(k->system->elements, i, MASS);

        if (mass < 0) {
            MSET(k->system->elements, i, MASS, -mass);
            MINC(k->system->elements, i, MA, -180.);
        }
        if (k->plRanges[0] != NULL) {
            for (int j = 0; j < ELEMENTS_SIZE; j++) {
                if (!IS_INVALID(MGET(k->plRanges[0], i, j)))
                    MSET(k->system->elements, i, j, MAX(MGET(k->system->elements, i, j), MGET(k->plRanges[0], i, j)));
                if (!IS_INVALID(MGET(k->plRanges[1], i, j)))
                    MSET(k->system->elements, i, j, MIN(MGET(k->system->elements, i, j), MGET(k->plRanges[1], i, j)));
            }
        }

    }

    for (int i = P_DATA_NOISE1; i <= P_DATA_NOISE10; i++)
        VSET(k->params, i, fabs(VGET(k->params, i)));

    MSET(k->system->elements, 0, MASS, k->Mstar);


    // If the coordinate system is jacobian, ensure the bodies are sorted by period.
    // If they were not sorted, then start a new setup (recalculate xyz coordinates + orbits).
    if (k->system->flag & JACOBI) {
        gsl_matrix* el = k->system->elements;
        ok_sort_matrix(el, PER);

        for (int i = 0; i < k->system->nplanets + 1; i++)
            if (MGET(k->system->elements, i, ORD) != i) {
                k->flags |= NEEDS_SETUP;
                MSET(k->system->elements, i, ORD, i);
            }
    }



    if (k->flags & NEEDS_SETUP)
        ok_setup(k->system);
}

int K_getActiveElements(ok_kernel* k) {
    int flags = 0;
    for (int i = 0; i < MROWS(k->plFlags); i++)
        for (int j = 0; j < MCOLS(k->plFlags); j++)
            if (MIGET(k->plFlags, i, j) & ACTIVE)
                flags++;

    return flags;
}

int K_getActivePars(ok_kernel* k) {
    int flags = 0;
    for (int i = 0; i < k->parFlags->size; i++)
        if (VIGET(k->parFlags, i) & ACTIVE)
            flags++;

    return flags;
}

int K_getNrPars(ok_kernel* k) {
    return K_getActivePars(k) + K_getActiveElements(k);
}

void K_calculate(ok_kernel* k) {
    K_validate(k);

    int ppars = K_getActiveElements(k);
    int vpars = K_getActivePars(k);

    double pars = k->ndata - ppars - vpars;
    k->nrpars = pars;

    bool integrate = (k->system->nplanets > 0);

    if (k->ndata <= 0)
        return;

    if ((!integrate && k->integration != NULL) || ((k->integration != NULL) && (k->times->size != k->integrationSamples ||
            MROWS(k->integration[0]->elements) != MROWS(k->system->elements)))) {
        for (int i = 0; i < k->integrationSamples; i++)
            ok_free_system(k->integration[i]);
        free(k->integration);
        k->integration = NULL;
        k->integrationSamples = 0;

    }

    for (int i = 0; i < k->ndata; i++) {
        k->compiled[i][T_PRED] = 0.;
        k->compiled[i][T_SCRATCH] = -1.;
    }

    if (integrate) {
        k->integration = ok_integrate(k->system, k->times, k->intOptions, k->intMethod, k->integration,
                &k->last_error);
    }

    k->chi2_rvs = 0.;
    k->chi2_tts = 0.;
    k->chi2_other = 0.;

    k->rms = 0.;
    k->rms_tts = 0.;
    k->jitter = 0.;
    k->nrvs = 0;
    k->ntts = 0;

    ok_integrator_options o;
    memcpy(&o, k->intOptions, sizeof (ok_integrator_options));
    o.calc_elements = false;


    double** compiled = k->compiled;
    int ndata = k->ndata;
    double epoch = k->system->epoch;

    if (k->model_function != NULL)
        k->model_function(k, compiled, ndata);

    for (int i = 0; i < ndata; i++) {
        int set = (int) compiled[i][T_SET];
        double n = K_getPar(k, set + DATA_SETS_SIZE);

        if ((int) compiled[i][T_FLAG] == T_RV) {
            if ((int) compiled[i][T_SCRATCH] < 0)
                compiled[i][T_PRED] += (integrate && k->integration != NULL ? ok_get_rv(k->integration[i]) : 0.0);

            compiled[i][T_SCRATCH] = 0;
            int j = (int) compiled[i][T_SET];
            compiled[i][T_SVAL] = compiled[i][T_VAL] - VGET(k->params, j) - VGET(k->params, P_RV_TREND) * (compiled[i][T_TIME] - epoch) - VGET(k->params, P_RV_TREND_QUADRATIC) * (compiled[i][T_TIME] - epoch) * (compiled[i][T_TIME] - epoch);
            double diff = compiled[i][T_SVAL] - compiled[i][T_PRED];
            double s = compiled[i][T_ERR];

            k->chi2_rvs += diff * diff / (s * s + n * n);
            k->rms += diff * diff;
            k->jitter += s*s;
            k->nrvs++;
        } else if ((int) compiled[i][T_FLAG] == T_TIMING) {
            int pidx = (int) compiled[i][T_TDS_PLANET];
            compiled[i][T_SVAL] = compiled[i][T_VAL];

            if (pidx <= 0)
                pidx = 1;
            if (pidx >= k->system->nplanets + 1)
                continue;

            if (integrate && k->integration != NULL && ((int) compiled[i][T_SCRATCH] < 0)) {
                double to = 0.;
                ok_find_closest_time_to_transit(k->integration[i],
                        pidx, &o, k->intMethod, o.eps_tr, compiled[i][T_TDS_FLAG], &to, &k->last_error);
                compiled[i][T_PRED] += to;
            }
            compiled[i][T_SCRATCH] = 0;

            double diff = compiled[i][T_SVAL] - compiled[i][T_PRED];
            double s = compiled[i][T_ERR];
            k->chi2_tts += diff * diff / (s * s + n * n);
            k->rms_tts += diff * diff;
            k->ntts++;
        } else if ((int) compiled[i][T_FLAG] == T_DUMMY) {
            // do nothing
        } else {
            double diff = compiled[i][T_SVAL] - compiled[i][T_PRED];
            double s = compiled[i][T_ERR];
            k->chi2_other += diff * diff / (s * s + n * n);
        }
    }

    k->chi2 = (k->chi2_rvs + k->chi2_tts + k->chi2_other) / pars;
    k->rms = sqrt(k->rms / (double) k->nrvs);
    k->rms_tts = sqrt(k->rms_tts / (double) k->ntts);
    k->jitter = sqrt(k->rms * k->rms - k->jitter / (double) k->nrvs);

    k->integrationSamples = k->times->size;
    k->flags &= ~NEEDS_COMPILE;
    k->flags &= ~NEEDS_SETUP;
}

ok_system** K_integrate(ok_kernel* k, gsl_vector* times, ok_system** bag, int* error) {
    ok_setup(k->system);
    return ok_integrate(k->system, times, k->intOptions, k->intMethod, bag, error);
}

ok_system** K_integrateProgress(ok_kernel* k, gsl_vector* times, ok_system** bag, int* error) {
    k->intOptions->progress = k->progress;
    ok_system** bag2 = K_integrate(k, times, bag, error);
    k->intOptions->progress = NULL;

    return bag2;
}

ok_system** K_integrateRange(ok_kernel* k, double from, double to, unsigned int samples, ok_system** bag, int* error) {
    gsl_vector* times = gsl_vector_alloc(samples);
    for (int i = 0; i < samples; i++)
        VSET(times, i, i * (to - from) / (samples - 1) + from);
    ok_system** sl = K_integrate(k, times, bag, error);
    gsl_vector_free(times);
    return sl;
}

gsl_matrix* K_integrateStellarVelocity(ok_kernel* k, double from, double to, unsigned int samples, ok_system** bag, int* error) {
    ok_system** sl = K_integrateRange(k, from, to, samples, NULL, error);
    gsl_matrix* m;

    if (sl == NULL && k->system->nplanets > 0)
        return NULL;
    else if (sl == NULL || k->system->nplanets == 0) {
        m = gsl_matrix_calloc(samples, 2);
        for (int i = 0; i < samples; i++) {
            MSET(m, i, 0, i * (to - from) / (samples - 1) + from);
        }
    } else
        m = ok_get_rvs(sl, samples);

    if (k->model_function != NULL) {
        double** dr = (double**) malloc(samples * sizeof (double*));
        for (int i = 0; i < samples; i++) {
            dr[i] = (double*) malloc(DATA_SIZE * sizeof (double));
            dr[i][0] = MGET(m, i, 0);
            dr[i][T_SET] = -1;
            dr[i][T_FLAG] = T_RV;
            dr[i][T_PRED] = 0.;
        }

        k->model_function(k, dr, samples);

        for (int i = 0; i < samples; i++) {
            MSET(m, i, 1, MGET(m, i, 1) + dr[i][T_PRED]);
            free(dr[i]);
        }
        free(dr);
    }

    ok_free_systems(sl, samples);
    return m;
}

K_GETSET_C(Mstar, Mstar, double)
K_GETSETW_C(epoch, Epoch, double, k->system->epoch)

K_GETSET_C(intOptions->abs_acc, IntAbsAcc, double)
K_GETSET_C(intOptions->rel_acc, IntRelAcc, double)
K_GETSET_C(intOptions->dt, IntDt, double)

K_GET_C(chi2, Chi2, double)

double K_getChi2_nr(ok_kernel* k) {
    return k->chi2_rvs + k->chi2_tts + k->chi2_other;
}

double K_getLoglik(ok_kernel* k) {
    double chi2 = k->chi2_rvs + k->chi2_tts + k->chi2_other;

    double A = 0;
    for (int i = 0; i < k->ndata; i++) {
        int set = (int) k->compiled[i][T_SET];
        double n = VGET(k->params, set + DATA_SETS_SIZE);

        A += log(SQR(k->compiled[i][T_ERR]) + n * n);
    }

    return 0.5 * A + 0.5 * chi2;
};

K_GET_C(rms, Rms, double)
K_GET_C(jitter, Jitter, double)

K_GET_C(nrvs, Nrvs, unsigned int)
K_GET_C(ntts, Ntts, unsigned int)
K_GET_C(nsets, Nsets, unsigned int)
K_GET_C(chi2_rvs, Chi2_rvs, double)
K_GET_C(chi2_tts, Chi2_tts, double)
K_GET_C(rms_tts, Rms_tts, double)

void K_setMinFunc(ok_kernel* k, ok_callback f) {
    if (f == NULL)
        k->minfunc = K_getChi2;
    else
        k->minfunc = f;
}

double K_getMinValue(ok_kernel* k) {
    return k->minfunc(k);
};

K_GET_C(minfunc, MinFunc, ok_callback)
K_GETSET_C(intMethod, IntMethod, int)
K_GETSET_C(intOptions, IntOptions, ok_integrator_options*)
K_GET_C(compiled, Compiled, double**)
K_GETSET_C(flags, Flags, unsigned int)
K_GETSET_C(progress, Progress, ok_progress)
K_GETSET_C(model_function, CustomModelFunction, ok_model_function)

K_GETSET_FREE_C(plSteps, ElementSteps, gsl_matrix*, gsl_matrix_free)
K_GETSET_FREE_C(parSteps, ParSteps, gsl_vector*, gsl_vector_free)
K_GETSET_FREE_C(plFlags, ElementFlags, gsl_matrix_int*, gsl_matrix_int_free)
K_GETSET_FREE_C(parFlags, ParFlags, gsl_vector_int*, gsl_vector_int_free)


void K_setElementFlag(ok_kernel* k, int row, int col, int value) {
    K_validate(k);
    MISET(k->plFlags, row, col, value);
}

int K_getElementFlag(ok_kernel* k, int row, int col) {
    K_validate(k);
    return MIGET(k->plFlags, row, col);
}

void K_setVoffFlag(ok_kernel* k, int idx, int value) {
    K_validate(k);
    VISET(k->parFlags, idx, value);
}

int K_getVoffFlag(ok_kernel* k, int idx) {
    K_validate(k);
    return VIGET(k->parFlags, idx);
}

void K_setElementStep(ok_kernel* k, int row, int col, double value) {
    K_validate(k);

    if (row == -1 && col == -1 && IS_INVALID(value)) {
        for (int i = 1; i < MROWS(k->system->elements); i++)
            for (int j = 0; j < ELEMENTS_SIZE; j++)
                MSET(k->plSteps, i, j, ok_default_steps[j]);
        return;
    }

    if (row == -1) {
        for (int i = 1; i < MROWS(k->system->elements); i++)
            MSET(k->plSteps, i, col, value);
    } else
        MSET(k->plSteps, row, col, value);
}

double K_getElementStep(ok_kernel* k, int row, int col) {
    K_validate(k);
    return MGET(k->plSteps, row, col);
}

void K_setElementRange(ok_kernel* k, int row, int col, double min, double max) {
    if (k->plRanges[0] == NULL) {
        k->plRanges[0] = gsl_matrix_alloc(MROWS(k->system->elements), ELEMENTS_SIZE);
        gsl_matrix_set_all(k->plRanges[0], INVALID_NUMBER);
    }
    if (k->plRanges[1] == NULL) {
        k->plRanges[1] = gsl_matrix_alloc(MROWS(k->system->elements), ELEMENTS_SIZE);
        gsl_matrix_set_all(k->plRanges[1], INVALID_NUMBER);
    }
    K_validate(k);

    if (row == -1) {
        for (int i = 1; i < MROWS(k->system->elements); i++) {
            MSET(k->plRanges[0], i, col, min);
            MSET(k->plRanges[1], i, col, max);
        }
    } else {
        MSET(k->plRanges[0], row, col, min);
        MSET(k->plRanges[1], row, col, max);
    }
}

void K_getElementRange(ok_kernel* k, int row, int col, double* min, double* max) {
    K_validate(k);

    if (k->plRanges[0] == NULL)
        *min = INVALID_NUMBER;
    else
        *min = MGET(k->plRanges[0], row, col);

    if (k->plRanges[1] == NULL)
        *max = INVALID_NUMBER;
    else
        *max = MGET(k->plRanges[1], row, col);
}

void K_setVoffStep(ok_kernel* k, int idx, double value) {
    K_validate(k);
    VSET(k->parSteps, idx, value);
}

double K_getVoffStep(ok_kernel* k, int idx) {
    K_validate(k);
    return VGET(k->parSteps, idx);
}

ok_kernel* K_clone(ok_kernel* k) {
    return K_cloneFlags(k, 0);
}

ok_kernel* K_cloneFlags(ok_kernel* k, unsigned int flags) {
    ok_kernel* k2 = (ok_kernel*) malloc(sizeof (ok_kernel));
    memcpy(k2, k, sizeof (ok_kernel));

    if (!(flags & SHARE_DATA)) {
        k2->datasets = (gsl_matrix**) malloc(sizeof (gsl_matrix*) * PARAMS_SIZE);

        for (int i = 0; i < k->nsets; i++) {
            k2->datasets[i] = ok_matrix_copy(k->datasets[i]);
        }

        k2->compiled = NULL;
        k2->flags |= NEEDS_COMPILE;
        k2->times = NULL;
        k2->integration = NULL;
        k2->integrationSamples = -1;
        k2->info = NULL;
        
        ok_info* ptr = k->info;
        ok_info* last = NULL;
        while (ptr != NULL) {
            ok_info* n = (ok_info*) malloc(sizeof(ok_info));
            n->next = n->tag = n->info = NULL;
            
            if (ptr->info != NULL && ptr->tag != NULL) {
                n->tag = strdup(ptr->tag);
                n->info = strdup(ptr->info);
            }
            if (last == NULL)
                k2->info = n;
            else
                last->next = n;
            last = n;
            ptr = ptr->next;
        }
    }
    if (!(flags & SHARE_RANGES)) {
        k2->plRanges = (gsl_matrix**) calloc(2, sizeof (gsl_matrix*));
        if (k->plRanges[0] != NULL) {
            k2->plRanges[0] = ok_matrix_copy(k->plRanges[0]);
            k2->plRanges[1] = ok_matrix_copy(k->plRanges[1]);
        }
        k2->parRanges = (gsl_vector**) calloc(2, sizeof (gsl_vector*));
        if (k->parRanges[0] != NULL) {
            k2->parRanges[0] = ok_vector_copy(k->parRanges[0]);
            k2->parRanges[1] = ok_vector_copy(k->parRanges[1]);
        }
    }

    if (!(flags & SHARE_STEPS)) {
        k2->plSteps = ok_matrix_copy(k->plSteps);
        k2->parSteps = ok_vector_copy(k->parSteps);
    }
    if (!(flags & SHARE_FLAGS)) {
        k2->plFlags = ok_matrix_int_copy(k->plFlags);
        k2->parFlags = ok_vector_int_copy(k->parFlags);
    }

    k2->system = ok_copy_system(k->system);

    k2->params = ok_vector_copy(k->params);

    k2->system->epoch = k->system->epoch;
    k2->system->time = k->system->time;
    k2->tag = k->tag;
    k2->rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(k2->rng, clock());
    k2->flags |= flags;
    k2->progress = k->progress;
    k2->model_function = k->model_function;
    k2->intOptions = (ok_integrator_options*) malloc(sizeof (ok_integrator_options));
    memcpy(k2->intOptions, k->intOptions, sizeof (ok_integrator_options));
    k2->intOptions->buffer = NULL;
    k2->intOptions->ibuffer = NULL;
    k2->intOptions->progress = NULL;
    return k2;
}

int K_minimize(ok_kernel* k, int algo, int maxiter, double params[]) {
    int ret = ok_minimizers[algo](k, maxiter, params);

    for (int i = 1; i < k->system->elements->size1; i++) {
        MSET(k->system->elements, i, MA, DEGRANGE(MGET(k->system->elements, i, MA)));
        MSET(k->system->elements, i, LOP, DEGRANGE(MGET(k->system->elements, i, LOP)));
        MSET(k->system->elements, i, INC, DEGRANGE(MGET(k->system->elements, i, INC)));
        MSET(k->system->elements, i, NODE, DEGRANGE(MGET(k->system->elements, i, NODE)));
    }

    return ret;
}

int K_1dminimize(ok_kernel* k, int algo, int maxiter, int row, int column, double params[]) {
    gsl_matrix_int* m = ok_matrix_int_copy(k->plFlags);
    gsl_vector_int* p = ok_vector_int_copy(k->parFlags);

    for (int i = 0; i < m->size1; i++)
        for (int j = 0; j < m->size2; j++)
            MISET(m, i, j, MIGET(m, i, j) & ~MINIMIZE);
    for (int i = 0; i < p->size; i++)
        VISET(p, i, VIGET(p, i) & ~MINIMIZE);

    gsl_matrix_int* old_plFlags = k->plFlags;
    gsl_vector_int* old_parFlags = k->parFlags;

    if (row > 0) {
        MISET(m, row, column, ACTIVE | MINIMIZE);
    } else if (row == -1) {
        VISET(p, column, ACTIVE | MINIMIZE);
    }


    k->plFlags = m;
    k->parFlags = p;
    int ret = K_minimize(k, algo, maxiter, params);
    k->plFlags = old_plFlags;
    k->parFlags = old_parFlags;
    gsl_matrix_int_free(m);
    gsl_vector_int_free(p);
    return ret;
}

void K_print(ok_kernel* k, FILE* f) {
    f = (f == NULL ? stdout : f);
    fprintf(f, "Elements:\n");
    ok_fprintf_matrix(k->system->elements, f, "%e ");
    fprintf(f, "Data:\n");
    for (int i = 0; i < k->nsets; i++) {
        fprintf(f, "data[%d] -> %e\n", i, k->params->data[i]);
    }
    
    fprintf(f, "Info:\n");
    ok_info* el = k->info;
    while (el != NULL) {
        fprintf(f, "%s = %s\n", el->tag, el->info);
        el = el->next;
    }
        
    
    fprintf(f, "Orbits:\n");
    ok_fprintf_matrix(k->system->orbits, f, "%e ");
    fprintf(f, "XYZ:\n");
    ok_fprintf_matrix(k->system->xyz, f, "%e ");
    fprintf(f, "Compiled data:\n");
    ok_fprintf_buf(k->compiled, f, "%e ", k->ndata, DATA_SIZE);
    fflush(f);
}

void K_setInfo(ok_kernel* k, const char* tag, const char* content) {
    ok_info* el = NULL;
    if (k->info != NULL) {
        ok_info* ptr = k->info;
        do {
            if (strcmp(tag, ptr->tag) == 0) {
                el = ptr;
                break;
            }
            ptr = ptr->next;
        } while (ptr != NULL);
    }
    
    char* info = (content != NULL ? strdup(content) : NULL);
    char* newtag = (tag != NULL ? strdup(tag) : NULL);
        
    
    if (el == NULL) {
        el = (ok_info*) malloc(sizeof(ok_info));
        
        el->tag = newtag;
        el->info = info;
        el->next = k->info;
        k->info = el;
    } else {
        
        if (el->tag != NULL) free(el->tag);
        if (el->info != NULL) free(el->info);
        el->tag = newtag;
        el->info = info;
    }
}

bool K_infoExists(ok_kernel* k, const char * tag) {
    return(K_getInfo(k, tag) != NULL);
    
}

char* K_getInfo(ok_kernel* k, const char* tag) {
    if (k->info == NULL)
        return NULL;
    ok_info* ptr = k->info;
    do {
        if (ok_str_iequals(ptr->tag, tag)) {
            return ptr->info;
        }
        ptr = ptr->next;
    } while (ptr != NULL);
    return NULL;
}

char* K_getInfoTag(ok_kernel* k, int n) {
    ok_info* el = k->info;
    int i = 0;
    while (el != NULL) {
        if (i == n)
            return el->tag;
        el = el->next;
        i++;
    }
    return "";
}

void K_clearInfo(ok_kernel* k) {
    ok_info* ptr = k->info;
    while (ptr != NULL) {
        ok_info* next = ptr->next;
        
        free(ptr->tag);
        free(ptr->info);
        free(ptr);
        ptr = next;
    }
    k->info = NULL;
}

void K_save_old(ok_kernel* k, const char* stem) {
    K_validate(k);


    char* fitName = ok_str_cat(stem, ".fit");
    char dataName[k->nsets][255];
    for (int i = 0; i < k->nsets; i++) {
        sprintf(dataName[i], "%s_%d.vels", stem, i);
    }
    char* sysName = ok_str_cat(stem, ".sys");

    FILE* fit = fopen(fitName, "w");
    fprintf(fit, "# InitialEpoch: %e\n\n", k->system->epoch);
    fprintf(fit, "Components %d\n", k->system->nplanets - 1);
    fprintf(fit, "PrimaryRVSet \"%s\"\n", dataName[0]);
    fprintf(fit, "OverallRVOffset %e\n", VGET(k->params, 0));
    fprintf(fit, "RelativeRVOffsets {\n");
    for (int i = 1; i < k->nsets; i++)
        fprintf(fit, "\t\"%s\" %e\n", dataName[i], VGET(k->params, 1));
    fprintf(fit, "\t\"Trend\" %e\n}\n", 0.);

    for (int i = 1; i < k->system->nplanets + 1; i++) {
        fprintf(fit, "\"\" {\n");
        fprintf(fit, "\tPeriod %e\n", K_getElement(k, i, PER));
        fprintf(fit, "\tMass %e\n", K_getElement(k, i, MASS));
        fprintf(fit, "\tMeanAnomaly %e\n", K_getElement(k, i, MA));
        fprintf(fit, "\tEccentricity %e\n", K_getElement(k, i, ECC));
        fprintf(fit, "\tLongOfPericenter %e\n", K_getElement(k, i, LOP));
        fprintf(fit, "\tInclination %e\n", K_getElement(k, i, INC));
        fprintf(fit, "\tNode %e\n", K_getElement(k, i, NODE));
        fprintf(fit, "\tRadius %e\n", K_getElement(k, i, RADIUS));
        fprintf(fit, "}\n");
    }
    fclose(fit);

    for (int i = 0; i < k->nsets; i++) {
        FILE* vels = fopen(dataName[i], "w");
        ok_fprintf_buf(k->compiled, vels, "%e ", k->ndata, DATA_SIZE);
        fclose(vels);
    }

    FILE* sys = fopen(sysName, "w");
    fprintf(sys, "Data {\n");
    for (int i = 0; i < k->nsets; i++)
        fprintf(sys, "\tRV[] \"%s\"\n", dataName[i]);
    fprintf(sys, "}\n\"%s\" {\n Mass %e\n}\n", stem, k->Mstar);
    fclose(sys);
    free(fitName);
    free(sysName);
}

bool K_save(ok_kernel* k, FILE* fid) {
    int digits = 15;
    int fract = 15;
    char fmt[] = "%15.15e ";

    fid = (fid != NULL ? fid : stdout);
    fprintf(fid, "@Kernel\n\n");
    fprintf(fid, "Epoch = %*.*e\n", digits, fract, k->system->epoch);
    fprintf(fid, "Mstar = %*.*e\n", digits, fract, k->Mstar);
    fprintf(fid, "Chi2 = %*.*e\n", digits, fract, k->chi2);
    fprintf(fid, "RMS = %*.*e\n", digits, fract, k->rms);
    fprintf(fid, "Jitter = %*.*e\n", digits, fract, k->jitter);
    fprintf(fid, "IntMethod = %d\n", k->intMethod);
    fprintf(fid, "Version = %.4f\n", SYSTEMIC_VERSION);

    fprintf(fid, "\nElements = %zu\n", k->system->elements->size1);
    ok_fprintf_matrix(k->system->elements, fid, fmt);

    fprintf(fid, "\nParams = \n");
    ok_fprintf_vector(k->params, fid, fmt);

    if (k->plFlags != NULL) {
        fprintf(fid, "\nPlFlags = %zu\n", k->plFlags->size1);
        ok_fprintf_matrix_int(k->plFlags, fid, "%d ");
    }

    fprintf(fid, "\nParFlags = \n");
    ok_fprintf_vector_int(k->parFlags, fid, "%d ");

    if (k->plSteps != NULL) {
        fprintf(fid, "\nPlSteps = %zu\n", k->plSteps->size1);
        ok_fprintf_matrix(k->plSteps, fid, fmt);
    }
    fprintf(fid, "\nParSteps = \n");
    ok_fprintf_vector(k->parSteps, fid, fmt);


    if (k->plRanges[0] != NULL) {
        fprintf(fid, "\nPlRanges_min = %zu\n", k->plRanges[0]->size1);
        ok_fprintf_matrix(k->plRanges[0], fid, fmt);
    }
    if (k->plRanges[1] != NULL) {
        fprintf(fid, "\nPlRanges_max = %zu\n", k->plRanges[1]->size1);
        ok_fprintf_matrix(k->plRanges[1], fid, fmt);
    }
    if (k->parRanges[0] != NULL) {
        fprintf(fid, "\nParRanges_min =\n");
        ok_fprintf_vector(k->parRanges[0], fid, fmt);
    }
    if (k->parRanges[1] != NULL) {
        fprintf(fid, "\nParRanges_max =\n");
        ok_fprintf_vector(k->parRanges[1], fid, fmt);
    }

    for (int i = 0; i < k->nsets; i++) {
        fprintf(fid, "\nData = %zu\n", k->datasets[i]->size1);
        ok_fprintf_matrix(k->datasets[i], fid, fmt);
    }

    fprintf(fid, "Flags = %d\n", k->flags);
    fprintf(fid, "System_Flags = %d\n", k->system->flag);


    ok_info* info = k->info;
    while (info != NULL) {
        if (info->tag != NULL && info->info != NULL)
            fprintf(fid, "$%s = %s\n", info->tag, info->info);
        info = info->next;
    }
    
    fprintf(fid, "\n@End\n");
    return true;
}

ok_kernel* K_load(FILE* fid, int skip) {
    fid = (fid != NULL ? fid : stdin);
    char line[MAX_LINE];
    for (int i = 0; i < skip; i++) {
        while (fgets(line, sizeof (line), fid))
            if (strcmp(line, "@End\n") == 0)
                break;
        if (feof(fid))
            return NULL;
    }

    bool found = false;
    char* ret;
    while (!found) {
        while ((ret = fgets(line, sizeof (line), fid)) != NULL) {
            if (strcmp(line, "@Kernel\n") == 0) {
                found = true;
                break;
            }
            if (feof(fid))
                return NULL;
        }

        if (ret == NULL) {
            return NULL;
        }

    }

    ok_kernel* k = K_alloc();


    double v;
    int r;
    k->nsets = 0;
    
    while (true) {
        if ((fgets(line, sizeof (line), fid) == NULL) || feof(fid) || strcmp(line, "@End\n") == 0)
            break;
        char tag[100] = {0};
        sscanf(line, "%s =", tag);

        if (strcmp(tag, "Epoch") == 0) {
            sscanf(line, "%*s = %le", &v);
            K_setEpoch(k, v);
        } else if (strcmp(tag, "Mstar") == 0) {
            sscanf(line, "%*s = %le", &v);
            K_setMstar(k, v);
        } else if (strcmp(tag, "Flags") == 0) {
            sscanf(line, "%*s = %d", &r);
            k->flags = r;
        } else if (strcmp(tag, "System_Flags") == 0) {
            sscanf(line, "%*s = %d", &r);
            k->system->flag = r;
        } else if (strcmp(tag, "IntMethod") == 0) {
            sscanf(line, "%*s = %d", &r);
            k->intMethod = r;
        } else if (strcmp(tag, "Elements") == 0) {
            sscanf(line, "%*s = %d", &r);
            gsl_matrix_free(k->system->elements);
            k->system->elements = gsl_matrix_alloc(r, ELEMENTS_SIZE);
            ok_fscanf_matrix(k->system->elements, fid);
            k->system->nplanets = r - 1;
        } else if (strcmp(tag, "Params") == 0) {
            ok_fscanf_vector(k->params, fid);
        } else if (strcmp(tag, "PlSteps") == 0) {
            sscanf(line, "%*s = %d", &r);
            gsl_matrix* ps = gsl_matrix_alloc(r, ELEMENTS_SIZE);
            ok_fscanf_matrix(ps, fid);
            gsl_matrix_free(k->plSteps);
            k->plSteps = ps;
        } else if (strcmp(tag, "ParSteps") == 0) {
            ok_fscanf_vector(k->parSteps, fid);
        } else if (strcmp(tag, "PlFlags") == 0) {
            sscanf(line, "%*s = %d", &r);
            gsl_matrix_int* plflags = gsl_matrix_int_alloc(r, ELEMENTS_SIZE);
            ok_fscanf_matrix_int(plflags, fid);
            gsl_matrix_int_free(k->plFlags);
            k->plFlags = plflags;
        } else if (strcmp(tag, "ParFlags") == 0) {
            ok_fscanf_vector_int(k->parFlags, fid);
        } else if (strcmp(tag, "PlRanges_min") == 0) {
            sscanf(line, "%*s = %d", &r);
            gsl_matrix* pr = gsl_matrix_alloc(r, ELEMENTS_SIZE);
            ok_fscanf_matrix(pr, fid);
            k->plRanges[0] = pr;
        } else if (strcmp(tag, "PlRanges_max") == 0) {
            sscanf(line, "%*s = %d", &r);
            gsl_matrix* pr = gsl_matrix_alloc(r, ELEMENTS_SIZE);
            ok_fscanf_matrix(pr, fid);
            k->plRanges[1] = pr;
        } else if (strcmp(tag, "ParRanges_min") == 0) {
            gsl_vector* pr = gsl_vector_alloc(PARAMS_SIZE);
            ok_fscanf_vector(pr, fid);
            k->parRanges[0] = pr;
        } else if (strcmp(tag, "ParRanges_max") == 0) {
            gsl_vector* pr = gsl_vector_alloc(PARAMS_SIZE);
            ok_fscanf_vector(pr, fid);
            k->parRanges[1] = pr;
        } else if (strcmp(tag, "Data") == 0) {
            sscanf(line, "%*s = %d", &r);
            gsl_matrix* d = gsl_matrix_alloc(r, DATA_SIZE);
            ok_fscanf_matrix(d, fid);
            k->datasets[k->nsets] = d;
            k->nsets++;
        } else if (strlen(tag) > 0 && tag[0] == '$') {
            ok_info* el = (ok_info*) malloc(sizeof(ok_info));
            el->next = NULL;
            el->tag = strdup(tag + 1);
            line[strlen(line) - 1] = '\0';
            el->info = (char*) malloc((strlen(line)+1) * sizeof(char));
            strcpy(el->info, line + strlen(tag) + 3);
            if (k->info == NULL) {
                k->info = el;
            } else {
                el->next = k->info;
                k->info = el;
            }
        }
    }

    k->flags = NEEDS_SETUP | NEEDS_COMPILE;
    K_validate(k);
    
    return k;

}

bool K_addDataFromSystem(ok_kernel* k, const char* filename) {
    K_removeData(k, -1);
    K_clearInfo(k);

    FILE* fid = fopen(filename, "r");
    if (fid == NULL)
        return false;
    
    char* buf = (char*) malloc(sizeof(char) * MAX_LINE);
    char token[MAX_LINE];
    
    char path[MAX_LINE];

    char* fn = strdup(filename);
    char* dn = dirname(fn);

    while (fgets(buf, MAX_LINE, fid) != 0) {
        char* line = ok_str_trim(buf);
        char* value = (char*) malloc(sizeof(char) * MAX_LINE);
        if (sscanf(line, "%s", token) == 1) {
            if ((strcmp(token, "RV[]") == 0) || (strcmp(token, "TD[]") == 0)) {
                if (sscanf(line, "%*s %s", value) == 1) {
                    
                    char* df = ok_str_trim(value);
                    
                    if (strlen(df) > 1 && df[0] == '"')
                        df++;
                    if (strlen(df) > 1 && df[strlen(df) - 1] == '"')
                        df[strlen(df) - 1] = '\0';
                    
                    char* ext = strrchr(df, '.');

                    if (!ok_file_readable(df))
                        sprintf(path, "%s/%s", dn, df);
                    else
                        strcpy(path, df);

                    if (ok_file_readable(path)) {
                        int type = T_RV;
                        if (strcmp(ext, ".vels") == 0)
                            type = T_RV;
                        else if (strcmp(ext, ".tds") == 0)
                            type = T_TIMING;

                        K_addDataFile(k, path, type);
                    }
                }
            } else if (strcmp(token, "Mass") == 0) {
                double mass = 1.;
                if (sscanf(line, "%*s %le", &mass) == 1) {
                    K_setMstar(k, mass);
                }
            } else if (strlen(token) > 1 && token[0] != '#' && token[0] != '{' && token[0] != '}' && token[0] != '"') {
                line = line + strlen(token);
                line = ok_str_trim(line);
                if (strlen(line) > 1) {
                    K_setInfo(k, token, line);
                }
            }
        }
        free(value);
    }

    free(buf);
    free(fn);
    fclose(fid);
    K_setInfo(k, "SystemFile", filename);
    return true;
}

void K_setSeed(ok_kernel* k, unsigned long int seed) {
    gsl_rng_set(k->rng, seed);
}

unsigned int K_getNplanets(ok_kernel* k) {
    return k->system->nplanets;
}

unsigned int K_getNdata(ok_kernel* k) {
    K_validate(k);
    return k->ndata;
}

void K_getRange(ok_kernel* k, double* from, double* to) {
    K_compileData(k);

    if (k->ndata < 1) {
        *from = INVALID_NUMBER;
        *to = INVALID_NUMBER;
    } else {
        *from = k->compiled[0][0];
        *to = k->compiled[k->ndata - 1][0];
    }
}

/*
gsl_matrix* K_combineData(ok_kernel* k, const int data_type) {
    gsl_matrix* filt;
    double** cd = K_compileData(k);
    
    if (data_type >= 0) 
        filt = ok_matrix_buf_filter(cd, k->ndata, DATA_SIZE, T_FLAG, data_type);
    else {
        filt = gsl_matrix_alloc(k->ndata, DATA_SIZE);
        
        for (int i = 0; i < k->ndata; i++)
            for (int j = 0; j < DATA_SIZE; j++)
                MSET(filt, i, j, cd[i][j]);
    }
    if (filt == NULL)
        return NULL;
    
    double* means = (double*) calloc(k->nsets, sizeof(double));
    double* n = (double*) calloc(k->nsets, sizeof(double));
    for (int i = 0; i < filt->size1; i++) {
        int set = (int) MGET(filt, i, T_SET);
        means[set] += MGET(filt, i, T_VAL);
        n[set] += 1.;
    }
    
    
    for (int i = 0; i < filt->size1; i++) {
        double val = MGET(filt, i, T_VAL);
        double pred = MGET(filt, i, T_PRED);
        int set = (int) MGET(filt, i, T_SET);
        MSET(filt, i, T_PRED, val - pred);
    }
    return filt;
}
 */

gsl_matrix* K_getXYZ(ok_kernel* k) {
    K_validate(k);
    return k->system->xyz;
}

/**
 * Sets the type of coordinates for this kernel. This will influence how orbital
 * elements get converted into cartesian coordinates, and vice versa. By default,
 * the coordinate system is ASTROCENTRIC.
 * @param k Kernel to modify
 * @param type One of ASTROCENTRIC or JACOBI.
 */
void K_setElementType(ok_kernel* k, int type) {
    int flag = k->system->flag;
    if (type == JACOBI) {
        k->system->flag |= JACOBI;
    } else {
        k->system->flag &= ~JACOBI;
    }
    if (flag != k->system->flag)
        K_validate(k);
}

/**
 * Returns the coordinate type used by the kernel (one of ASTROCENTRIC or JACOBI).
 * @param k A kernel 
 * @return The coordinate type used by the kernel (one of ASTROCENTRIC or JACOBI).
 */
int K_getElementType(ok_kernel* k) {
    if (k->system->flag & JACOBI)
        return JACOBI;
    else
        return ASTROCENTRIC;
}

/**
 * Perturbs the parameters (orbital elements + pars) of the specified kernel 
 * that are flagged MINIMIZE with a gaussian rng with variance given by the 
 * corresponding step.
 * @param k
 */
void K_perturb(ok_kernel* k) {
    for (int i = 0; i < ELEMENTS_SIZE; i++)
        for (int j = 0; j < K_getNplanets(k) + 1; j++) {
            if (MGET(k->plFlags, j, i) & MINIMIZE) {
                double step = gsl_ran_gaussian(k->rng, MGET(k->plSteps, j, i));
                MINC(k->system->elements, j, i, step);
            }
        }
    for (int i = 0; i < PARAMS_SIZE; i++)
        if (VGET(k->parFlags, i) & MINIMIZE) {
            double step = VGET(k->parSteps, i);
            VINC(k->params, i, gsl_ran_gaussian(k->rng, step));
        }
    K_validate(k);
}

void* ok_bridge_kernel_buf(void* buf, int n, ok_kernel* k) {
    if (buf == NULL) {
        return (malloc(sizeof (ok_kernel*) * n));
    } else if (n < 0) {
        free(buf);
        buf = NULL;
    } else {
        ((ok_kernel**) buf)[n] = k;
    }
    return buf;
}

void K_getMinimizedIndex(ok_kernel* k, int index, int* row, int* column) {
    int npars = 0;
    for (int i = 1; i < k->system->nplanets + 1; i++)
        for (int j = 0; j < ELEMENTS_SIZE; j++)
            if ((MIGET(k->plFlags, i, j) & MINIMIZE ? 1 : 0)) {
                if (index == npars) {
                    *row = i;
                    *column = j;
                    return;
                }
                npars++;
            }
    for (int i = 0; i < PARAMS_SIZE; i++)
        if ((VIGET(k->parFlags, i) & MINIMIZE ? 1 : 0)) {
            if (index == npars) {
                *row = -1;
                *column = i;
                return;
            }
            npars++;
        }

    *row = -1;
    *column = -1;
}

void K_setMinimizedValues(ok_kernel* k, double* values) {
    int npars = 0;
    for (int i = 1; i < k->system->nplanets + 1; i++)
        for (int j = 0; j < ELEMENTS_SIZE; j++)
            if ((MIGET(k->plFlags, i, j) & MINIMIZE ? 1 : 0)) {
                MSET(k->system->elements, i, j, values[npars]);
                npars++;
            }
    for (int i = 0; i < PARAMS_SIZE; i++)
        if ((VIGET(k->parFlags, i) & MINIMIZE ? 1 : 0)) {
            VSET(k->params, i, values[npars]);
            npars++;
        }
    k->flags |= NEEDS_COMPILE | NEEDS_SETUP;
}

void K_getMinimizedValues(ok_kernel* k, double* values) {
    int npars = 0;
    for (int i = 1; i < k->system->nplanets + 1; i++)
        for (int j = 0; j < ELEMENTS_SIZE; j++)
            if ((MIGET(k->plFlags, i, j) & MINIMIZE ? 1 : 0)) {
                values[npars] = MGET(k->system->elements, i, j);
                npars++;
            }
    for (int i = 0; i < PARAMS_SIZE; i++)
        if ((VIGET(k->parFlags, i) & MINIMIZE ? 1 : 0)) {
            values[npars] = VGET(k->params, i);
            npars++;
        }
    k->flags |= NEEDS_COMPILE | NEEDS_SETUP;
}

ok_kernel_minimizer_pars K_getMinimizedVariables(ok_kernel* k) {
    ok_kernel_minimizer_pars mpars;

    mpars.npars = 0;

    // Count all the parameters to minimize on
    for (int i = 1; i < k->system->nplanets + 1; i++)
        for (int j = 0; j < ELEMENTS_SIZE; j++)
            mpars.npars += (MIGET(k->plFlags, i, j) & MINIMIZE ? 1 : 0);
    for (int i = 0; i < k->parFlags->size; i++)
        mpars.npars += (VIGET(k->parFlags, i) & MINIMIZE ? 1 : 0);

    double** pars = (double**) malloc(sizeof (double*) * mpars.npars);
    double* steps = (double*) malloc(sizeof (double) * mpars.npars);
    int* type = (int*) malloc(sizeof (int) * mpars.npars);
    double* min = (double*) malloc(sizeof (double) * mpars.npars);
    double* max = (double*) malloc(sizeof (double) * mpars.npars);
    int idx = 0;
    for (int i = 1; i < k->system->nplanets + 1; i++)
        for (int j = 0; j < ELEMENTS_SIZE; j++)
            if (MIGET(k->plFlags, i, j) & MINIMIZE) {
                pars[idx] = gsl_matrix_ptr(k->system->elements, i, j);
                steps[idx] = MGET(k->plSteps, i, j);
                type[idx] = j;
                K_getElementRange(k, i, j, &(min[idx]), &(max[idx]));
                idx++;
            }
    for (int i = 0; i < k->parFlags->size; i++)
        if (VIGET(k->parFlags, i) & MINIMIZE) {
            pars[idx] = gsl_vector_ptr(k->params, i);
            steps[idx] = VGET(k->parSteps, i);
            type[idx] = -1;
            K_getParRange(k, i, &(min[idx]), &(max[idx]));
            idx++;
        }

    mpars.pars = pars;
    mpars.steps = steps;
    mpars.type = type;
    mpars.min = min;
    mpars.max = max;
    return mpars;
}
