/* 
 * File:   simplex.h
 * Author: stefano
 *
 * Created on February 23, 2012, 4:39 PM
 */

#ifndef MINPACK_H
#define	MINPACK_H


/**
 * Attempts to converge to a local minimum of the kernel k, using the lev-mar
 * library. Many of the algorithm settings are determined by the kernel object:
 * minimization flags (whether a parameter is held fixed or varied; through K_setFlags
 * and K_setVFlags to MINIMIZE), the initial vector of parameters (the state of the system; set
 * by loading the data and the orbital parameters) the step of parameters and elements (giving
 * the initial step in the simplex algorithm) and finally the merit function (through 
 * K_setMinFunc; by default, the reduced Chi^2 of the system). The state of the kernel
 * object is modified; at the end of the routine, it contains the parameters 
 * corresponding to the local minimum. 
 * 
 * @param k The kernel object containing the state of the system.
 * @param step A N-elements vector containing the initial step size used by 
 * the simplex algorithm, where N is total number of parameters being minimized. For
 * instance, if there are S planets, T velocity offsets and 5 parameters to minimize per
 * planet, N = 5S + T. 
 * @return 
 */

int K_minimize_lm(ok_kernel* k, int maxiter, double params[]);

#endif	/* SIMPLEX_H */

