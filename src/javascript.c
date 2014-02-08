#include <stdio.h>
#include "systemic.h"
#include "kernel.h"
#include "mcmc.h"
#include "bootstrap.h"
#include "time.h"
#include "periodogram.h"
#include "javascript.h"
#include "integration.h"

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
    static double tolerance[1] = { 1e-3 };
    static int target_points = 600;
    
    if (row < 0) {
        int samples = col;
        if (rvline != NULL) {
            gsl_matrix_free(rvline);
            rvline = NULL;
        }
        if (k->ndata == 0)
            return -1;
        
        double** comp = K_getCompiled(k);
        
        gsl_matrix* rvline_full = K_integrateStellarVelocity(k, comp[0][0],
                comp[k->ndata-1][0], 
                samples,
                NULL, NULL);
        int fac = 1;
        
        rvline = ok_resample_curve(rvline_full,
                0, 1, 1, target_points, 100, tolerance, 0, false);
        
        if (MROWS(rvline) > 1.5 * target_points) {
            gsl_matrix_free(rvline);
            rvline = ok_resample_curve(rvline_full,
                0, 1, 0.2, target_points, 100, tolerance, 0, false);
            fac = -1;
        }
        
        gsl_matrix_free(rvline_full);
        
        
        
        return fac*(int)MROWS(rvline);
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
        if (phased_data != NULL) {
            gsl_matrix_free(phased_data);
            phased_data = NULL;
        }
        double chi2 = k->chi2;
        double rms = k->rms;
        double jitter = k->jitter;
        double chi2_rvs = k->chi2_rvs;
        
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
        K_calculate(k);
        k->chi2 = chi2;
        k->rms = rms;
        k->jitter = jitter;
        k->chi2_rvs = chi2_rvs;
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
        if (phasedRVLine != NULL) {
            gsl_matrix_free(phasedRVLine);
            phasedRVLine = NULL;
        }
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
    static int samples = 15000;
    static double Pmin = 1.;
    static double Pmax = 20000.;
    static ok_periodogram_workspace* p = NULL;
    static gsl_matrix* ps = NULL;
    static const int top_freqs = 10;
    static double* top = NULL;
    static double tolerance[1] = {1e-3};
    
    if (p == NULL) {
        p = (ok_periodogram_workspace*) malloc(sizeof(ok_periodogram_workspace));
        p->buf = NULL;
        p->per = NULL;
        p->calc_z_fap = true;
    }
    if (row == JS_PS_GET_TOP_PERIODS) {
        return top[col];  
    } else if (row == JS_PS_GET_TOP_POWERS) {
        return top[col+top_freqs];
    } else if (row == JS_PS_GET_TOP_FAPS) {
        return top[col+2*top_freqs];
    } else if (row == JS_PS_SET_PMIN) {
        Pmin = (double) col;        
        return 0;
    } else if (row == JS_PS_SET_PMAX) {
        Pmax = (double) col;
        return 0;
    } else if (row == JS_PS_SETUP) {
        if (ps != NULL) {
            gsl_matrix_free(ps);
            ps = NULL;
        }
        if (top == NULL)
            top = (double*) malloc(top_freqs * 3 * sizeof(double));
        gsl_matrix* data = K_getCompiledDataMatrix(k);
        for (int i = 0; i < MROWS(data); i++)
            MSET(data, i, T_SVAL, MGET(data, i, T_SVAL)-MGET(data, i, T_PRED));

        gsl_matrix* ret = ok_periodogram_ls(data, samples, Pmin, Pmax, 
                0, T_TIME, T_SVAL, T_ERR, p);
        
        ps = ok_resample_curve(ret, 0, 1, 0.1, 800,
                100, tolerance, 0, true);
        length = MROWS(ps);
        
        ok_sort_matrix(ret, PS_Z);
        double dt = 0.5;
        int idx = MROWS(ret);
        int i = 0;
        while (idx > 0 && i < top_freqs) {
            idx--;
            bool skip = false;
            for (int n = 0; n < i; n++)
                if (fabs(top[n] - MGET(ret, idx, PS_TIME)) < dt)
                    skip = true;
            
            if (!skip) {
                top[i] = MGET(ret, idx, PS_TIME);
                top[i+top_freqs] = MGET(ret, idx, PS_Z);
                top[i+2*top_freqs] = MGET(ret, idx, PS_FAP);   
                i++;
            }
        }
        gsl_matrix_free(data);
        
        return (double) length;
    } else if (row == JS_PS_GET_FAPS_LEVELS) {
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
    K_minimize(k, SIMPLEX, 5000, NULL);
    k->progress = NULL;
    return failed;
}

double K_integrateForward(ok_kernel* k, const int mode, const double nyears,
        const int row, const int col) {
    
    static ok_system* sys = NULL;
    static double time;
    static double yearspast;
    static double dt;
    static gsl_vector* times = NULL;
    static gsl_matrix* els = NULL;
    
    if (K_getNplanets(k) == 0)
        return 0;
    
    if (mode == JS_I_START) {
        if (times == NULL)
            times = gsl_vector_alloc(2);
        dt = MIN(nyears, 10);
        time = K_getEpoch(k);
        yearspast = 0;
        sys = ok_copy_system(k->system);
        ok_setup(sys);
    } else if (mode == JS_I_STEP) {
        if (yearspast >= nyears)
            return JS_I_ENDREACHED;
        
        if (els != NULL) {
            gsl_matrix_free(els);
            els = NULL;
        }
        VSET(times, 0, time);
        VSET(times, 1, time + dt * 365.25);
        
        yearspast += dt;
        
        int error;
        ok_system** bag = ok_integrate(sys, times, k->intOptions, RK89, NULL, &error);
        els = ok_get_els(bag, 2, false);
        
        ok_free_system(sys);
        sys = ok_copy_system(bag[1]);
        ok_free_systems(bag, 2);
        time += dt * 365.25;
        
        if (error != INTEGRATION_SUCCESS)
            return -error;
        
        return yearspast;
    } else if (mode == JS_I_END) {
        if (els != NULL) {
            gsl_matrix_free(els);
            els = NULL;
        }
        if (sys != NULL) {
            ok_free_system(sys);
            sys = NULL;
        }
    } else if (mode == JS_I_GET) {
        if (col == SMA) {
                int idx = row * ELEMENTS_SIZE + PER + 1;
                double mstar = K_getMstar(k) * K2;
                double mp = K_getElement(k, row, MASS) * MJUP / MSUN * K2;
                double P = MGET(els, 1, idx);
                
                return ok_acalc(P, mstar, mp);
        } else {
                int idx = row * ELEMENTS_SIZE + col + 1;
                return MGET(els, 1, idx);
        }
    };
    
    
    return 0;
}

int main() {
    
}
