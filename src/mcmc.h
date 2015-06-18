/* 
 * File:   mcmc.h
 * Author: stefano
 *
 * Created on March 6, 2012, 1:25 PM
 */

#ifndef MCMC_H
#define	MCMC_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "systemic.h"
#include "kernel.h"
#include "kl.h"

#define MCMC
    ok_list* K_mcmc_single(ok_kernel* k, unsigned int nsteps, unsigned int skip, unsigned int discard, const double dparams[], ok_list* cont, ok_callback2 merit_function, int tag, int* flag);
    ok_list* K_mcmc_mult(ok_kernel** k, unsigned int nchains, unsigned int ntemps, unsigned int skip, unsigned int discard, const double params[], double Rstop, ok_callback2 merit_function);
    double K_default_prior(const ok_kernel* k);
    void K_mcmc_likelihood_and_prior_default(ok_kernel* k, double* ret);

#ifdef	__cplusplus
}
#endif

#endif	/* MCMC_H */

