#include <stdio.h>
#include "systemic.h"
#include "kernel.h"
#include "mcmc.h"
#include "bootstrap.h"
#include "time.h"
#include "periodogram.h"
#include "javascript.h"

/* Contains a simplified interface for calling from the web app */


double K_getDataAt(ok_kernel* k, int subset, int row, int column) {
    if (subset == ALL)
        return K_getCompiled(k)[row][column];
    else
        return MGET(K_getData(k, subset), row, column);
}

void K_setDataAt(ok_kernel* k, int subset, int row, int column, double val) {
    if (subset == ALL)
        K_getCompiled(k)[row][column] = val;
    else
        MSET(K_getData(k, subset), row, column, val);
}


double K_getRVLine(ok_kernel* k, int row, int col) {
    static gsl_matrix* rvline = NULL;
    if (row < 0) {
        int samples = -row;
        if (rvline != NULL)
            gsl_matrix_free(rvline);
        
        if (k->ndata == 0)
            return -1;
        
        double** comp = K_getCompiled(k);
        
        rvline = K_integrateStellarVelocity(k, comp[0][0],
                comp[k->ndata-1][0], 
                samples,
                NULL, NULL);
        return 1;
    } else {
        if (rvline == NULL)
            return INVALID_NUMBER;
        else
            return MGET(rvline, row, col);
    }
}



double K_getPhasedDataForPlanet(ok_kernel* k, int planet, int row, int column) {
    static gsl_matrix* phased_data = NULL;

    if (planet >= 1) {
        if (phased_data != NULL)
            gsl_matrix_free(phased_data);
        
        planet = MIN(planet, K_getNplanets(k));
        double mass = K_getElement(k, planet, MASS);
        double period = K_getElement(k, planet, PER);
        K_setElement(k, planet, MASS, 0);
        K_calculate(k);
        
        phased_data = K_getCompiledDataMatrix(k);
        double mint = MGET(phased_data, 0, T_TIME);
        for (int i = 0; i < MROWS(phased_data); i++) {
            double t = fmod((MGET(phased_data, i, T_TIME) - mint), period);
            double v = MGET(phased_data, i, T_SVAL)-MGET(phased_data, i, T_PRED);
            MSET(phased_data, i, T_TIME, t);
            MSET(phased_data, i, T_VAL, v);
        }
        
        ok_sort_matrix(phased_data, T_TIME);
        K_setElement(k, planet, MASS, mass);
        return 1;
    } else {
        return MGET(phased_data, row, column);
    }
}


double K_getPhasedRVLine(ok_kernel* k, int planet, int row, int column) {
    static gsl_matrix* phasedRVLine = NULL;
    if (planet >= 1) {
        if (k->ndata == 0)
            return -1;
        
        int np = K_getNplanets(k);
        double masses[np+1];
        double periods[np+1];
        for (int i = 1; i <= np; i++) {
            masses[i] = K_getElement(k, i, MASS);
            periods[i] = K_getElement(k, i, PER);
            if (i != planet) {
                K_setElement(k, i, MASS, 0.);
                K_setElement(k, i, PER, 10000.);
            }
        };
        
        double period = K_getElement(k, planet, PER);
        int samples = -row;
        if (phasedRVLine != NULL)
            gsl_matrix_free(phasedRVLine);
        
        double** comp = K_getCompiled(k);
        
        phasedRVLine = K_integrateStellarVelocity(k, comp[0][0],
                comp[k->ndata-1][0], 
                samples,
                NULL, NULL);
        
        double mint = MGET(phasedRVLine, 0, T_TIME);
        for (int i = 0; i < MROWS(phasedRVLine); i++) {
            double t = fmod((MGET(phasedRVLine, i, 0) - mint), period);
            MSET(phasedRVLine, i, 0, t);
        }
        ok_sort_matrix(phasedRVLine, 0);
        
        for (int i = 1; i <= np; i++) {
            K_setElement(k, i, MASS, masses[i]);
            K_setElement(k, i, PER, periods[i]);
        }
        
        return 1;
    } else {
        return MGET(phasedRVLine, row, column);
    }
}


double K_getPeriodogramAt(ok_kernel* k, int row, int col) {
    
    static int length;
    static int samples = 30000;
    
    static ok_periodogram_workspace* p = NULL;
    static gsl_matrix* ps = NULL;
    if (p == NULL) {
        p = (ok_periodogram_workspace*) malloc(sizeof(ok_periodogram_workspace));
        p->buf = NULL;
        p->per = NULL;
        p->calc_z_fap = true;
    }
    if (row == -1) {
        if (ps != NULL)
            gsl_matrix_free(ps);
        
        gsl_matrix* data = K_getCompiledDataMatrix(k);
        for (int i = 0; i < MROWS(data); i++)
            MSET(data, i, T_SVAL, MGET(data, i, T_SVAL)-MGET(data, i, T_PRED));

        gsl_matrix* ret = ok_periodogram_ls(data, samples, 1., 20000., 
                0, T_TIME, T_SVAL, T_ERR, p);
        
        ps = ok_resample_curve(ret, 0, 1, 30, 0.1);
        length = MROWS(ps);
        gsl_matrix_free(data);
        
        return (double) length;
    } else if (row == -2) {
        if (p == NULL || ps == NULL)
            return 0;
        if (col == 1)
            return p->z_fap_1;
        else if (col == 2)
            return p->z_fap_2;
        else if (col == 3)
            return p->z_fap_3;
        else
            return 0.;
    }  else {
        if (ps == NULL)
            return 0;
        return MGET(ps, row, col);
    }
}

time_t mtime;
int timeout;
int failed;

int progressWithTimeout(int current, int max, void* state, const char* function) {
    if (difftime(time(NULL), mtime) > timeout) {
        failed = timeout;
        return PROGRESS_STOP;
    }
    failed = 0;
    return PROGRESS_CONTINUE;
}

int K_minimizeWithTimeout(ok_kernel* k, int to) {
    mtime = time(NULL);
    timeout = to;
    failed = 0;
    k->progress = progressWithTimeout;
    K_minimize(k, SIMPLEX, 1000, NULL);
    k->progress = NULL;
    return failed;
}

int main() {
    
}