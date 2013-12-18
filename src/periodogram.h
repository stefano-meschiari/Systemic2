/* 
 * File:   periodogram.h
 * Author: stefano
 *
 * Created on March 26, 2012, 1:56 PM
 */

#ifndef PERIODOGRAM_H
#define	PERIODOGRAM_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "systemic.h"
#include "gsl/gsl_matrix.h"
    
#define PS_TIME 0
#define PS_Z 1
#define PS_FAP 2
#define PS_Z_LS 3
#define PS_TAU 4
#define PS_WIN 5
#define PS_SIZE 6
    
#define PS_TYPE_DATA 0
#define PS_TYPE_RESIDUALS 1
    
    typedef struct {
        double W;
        double zmax;
        double z2max;
        double fapmax;
        double z_fap_1;
        double z_fap_2;
        double z_fap_3;
        bool calc_z_fap;
        gsl_matrix* per;
        gsl_matrix* buf;
        gsl_vector* zm;
    } ok_periodogram_workspace;
    
    gsl_matrix* ok_periodogram_ls(const gsl_matrix* data, const unsigned int samples, const double Pmin, const double Pmax, const int method,
        unsigned int timecol, unsigned int valcol, unsigned int sigcol, ok_periodogram_workspace* p);

    gsl_matrix* ok_periodogram_boot(const gsl_matrix* data, const unsigned int trials, const unsigned int samples, 
        const double Pmin, const double Pmax, const int method,
        const unsigned int timecol, const unsigned int valcol, const unsigned int sigcol,
        const unsigned long int seed, ok_periodogram_workspace* p, ok_progress prog);

    
    gsl_matrix* ok_periodogram_full(ok_kernel* k, int type, int algo, bool circular, unsigned int sample,
        const unsigned int samples, const double Pmin, const double Pmax);

#ifdef	__cplusplus
}
#endif

#endif	/* PERIODOGRAM_H */

