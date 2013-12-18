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

gsl_matrix* rvline = NULL;
double K_getRVLine(ok_kernel* k, int row, int col) {
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


double** ps_resample(gsl_matrix* ps, int* length) {
    
    double** list = (double**) malloc(MROWS(ps) * sizeof(double*));
    list[0] = gsl_matrix_ptr(ps, 0, 0);
    int len = 1;
    double last_t = list[0][0];
    const double min_dt = 0.25;
    
    for (int i = 1; i < MROWS(ps) - 1; i++) {
        if (MGET(ps, i, 1) > MGET(ps, i-1, 1) && MGET(ps, i, 1) > MGET(ps, i+1, 1) &&
                MGET(ps, i, 0) - last_t > min_dt) {
            list[len] = gsl_matrix_ptr(ps, i, 0);
            last_t = list[len][0];
            len++;
        } else if (MGET(ps, i, 1) < MGET(ps, i-1, 1) && MGET(ps, i, 1) < MGET(ps, i+1, 1) &&
                        MGET(ps, i, 0) - last_t > min_dt) {
            list[len] = gsl_matrix_ptr(ps, i, 0);
            last_t = list[len][0];
            len++;
        } else if (MGET(ps, i, 0) - last_t > 3 * min_dt) {
            list[len] = gsl_matrix_ptr(ps, i, 0);
            last_t = list[len][0];
            len++;
        }
    }
    list[len] = gsl_matrix_ptr(ps, MROWS(ps)-1, 0);
    len++;
    *length = len;
    return list;
}

double K_getPeriodogramAt(ok_kernel* k, int row, int col) {
    
    static int length;
    static double** buf = NULL;
    static gsl_matrix* ps = NULL;
    
    static ok_periodogram_workspace* p = NULL;
    if (p == NULL) {
        p = (ok_periodogram_workspace*) malloc(sizeof(ok_periodogram_workspace));
        p->buf = NULL;
        p->per = NULL;
    }
    if (row == -1) {
        if (buf != NULL) {
            free(buf);
        }
        if (ps != NULL)
            gsl_matrix_free(ps);
        
        if (p->buf != NULL) {
            gsl_matrix_free(p->buf);
        }
        p->buf = NULL;
        p->per = NULL;
                
        gsl_matrix* data = K_getCompiledDataMatrix(k);
        for (int i = 0; i < MROWS(data); i++)
            MSET(data, i, T_SVAL, MGET(data, i, T_SVAL)-MGET(data, i, T_PRED));

        ps = ok_periodogram_ls(data, 40000, 0.5, 20000., 
                0, T_TIME, T_SVAL, T_ERR, p);
        buf = ps_resample(ps, &length);
        
        gsl_matrix_free(data);
        return (double) length;
    } else if (row == -2) {
        if (col == 1)
            return p->z_fap_1;
        else if (col == 2)
            return p->z_fap_2;
        else if (col == 3)
            return p->z_fap_3;
        else
            return 0.;
    } else {
        if (buf == NULL)
            return INVALID_NUMBER;
        else
            return buf[row][col];
    }
}

time_t mtime;
int timeout;
int failed;

int progressWithTimeout(int current, int max, void* state, const char* function) {
    if (time(NULL) - mtime > timeout) {
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