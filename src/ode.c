/* ode.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/
#include "utils.h"

#include "systemic.h"
#include "integration.h"
#include "ode.h"
#include "assert.h"
#include "f2c.h"
/* Table of constant values */

static int c__1 = 1;
static doublereal c_b11 = 1.;
#define oneu .222e-15


static inline double ok_d_sign(const double a, const double b) {
    double x = (a >= 0 ? a : -a);
    return (b >= 0 ? x : -x);
}

static inline int ok_i_sign(const int a, const int b) {
    int x = (a >= 0 ? a : -a);
    return (b >= 0 ? x : -x);
}

typedef int (*ok_force_function) (const double t, const double y[], double dydt[],  void * params);

/* Subroutine */ int ok_ode_ode_(ok_force_function f, const int neqn, doublereal *y, doublereal *t,
	 doublereal *tout, doublereal *relerr, doublereal *abserr, int *
	iflag, doublereal *work, int *iwork, void* params)
{
    /* Initialized data */

     int ialpha = 1;
     int ih = 89;
     int ihold = 90;
     int istart = 91;
     int itold = 92;
     int idelsn = 93;
     int ibeta = 13;
     int isig = 25;
     int iv = 38;
     int iw = 50;
     int ig = 62;
     int iphase = 75;
     int ipsi = 76;
     int ix = 88;
     
     doublereal g[13] = { 1. };
     doublereal rho[13] = { 1. };


    extern /* Subroutine */ int ok_ode_de_(ok_force_function, const int, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, int *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, logical *,
	     doublereal *, doublereal *, doublereal *, doublereal *, logical *
	    , doublereal *, doublereal *, int *, logical *, int *, 
	    int *, int *, void *, doublereal*, doublereal*);
     int ip, iyp, iwt, iyy, iphi;
     logical nornd, start, phase1;
     int iypout;


/*   double precision subroutine ode integrates a system of neqn */
/*   first order ordinary differential equations of the form: */
/*             dy(i)/dt = f(t,y(1),y(2),...,y(neqn)) */
/*             y(i) given at  t . */
/*   the subroutine integrates from  t  to  tout .  on return the */
/*   parameters in the call list are set for continuing the integration. */
/*   the user has only to define a new value  tout  and call  ode  again. */

/*   the differential equations are actually solved by a suite of codes */
/*   de ,  step , and  intrp .  ode  allocates virtual storage in the */
/*   arrays  work  and  iwork  and calls  de .  de  is a supervisor which */
/*   directs the solution.  it calls on the routines  step  and  intrp */
/*   to advance the integration and to interpolate at output points. */
/*   step  uses a modified divided difference form of the adams pece */
/*   formulas and local extrapolation.  it adjusts the order and step */
/*   size to control the local error per unit step in a generalized */
/*   sense.  normally each call to  step  advances the solution one step */
/*   in the direction of  tout .  for reasons of efficiency  de */
/*   integrates beyond  tout  internally, though never beyond */
/*   t+10*(tout-t), and calls  intrp  to interpolate the solution at */
/*   tout .  an option is provided to stop the integration at  tout  but */
/*   it should be used only if it is impossible to continue the */
/*   integration beyond  tout . */

/*   this code is completely explained and documented in the text, */
/*   computer solution of ordinary differential equations:  the initial */
/*   value problem  by l. f. shampine and m. k. gordon. */

/*   the parameters represent: */
/*      f -- double precision subroutine f(t,y,yp) to evaluate */
/*                derivatives yp(i)=dy(i)/dt */
/*      neqn -- number of equations to be integrated (int*4) */
/*      y(*) -- solution vector at t                 (real*8) */
/*      t -- independent variable                    (real*8) */
/*      tout -- point at which solution is desired   (real*8) */
/*      relerr,abserr -- relative and absolute error tolerances for local */
/*           error test (real*8).  at each step the code requires */
/*             dabs(local error) .le. dabs(y)*relerr + abserr */
/*           for each component of the local error and solution vectors */
/*      iflag -- indicates status of integration     (int*4) */
/*      work(*)  (real*8)  -- arrays to hold information internal to */
/*      iwork(*) (int*4)    which is necessary for subsequent calls */

/*   first call to ode -- */

/*   the user must provide storage in his calling program for the arrays */
/*   in the call list, */
/*      y(neqn), work(100+21*neqn), iwork(5), */
/*   declare  f  in an external statement, supply the double precision */
/*   subroutine f(t,y,yp)  to evaluate */
/*      dy(i)/dt = yp(i) = f(t,y(1),y(2),...,y(neqn)) */
/*   and initialize the parameters: */
/*      neqn -- number of equations to be integrated */
/*      y(*) -- vector of initial conditions */
/*      t -- starting point of integration */
/*      tout -- point at which solution is desired */
/*      relerr,abserr -- relative and absolute local error tolerances */
/*      iflag -- +1,-1.  indicator to initialize the code.  normal input */
/*           is +1.  the user should set iflag=-1 only if it is */
/*           impossible to continue the integration beyond  tout . */
/*   all parameters except  f ,  neqn  and  tout  may be altered by the */
/*   code on output so must be variables in the calling program. */

/*   output from  ode  -- */

/*      neqn -- unchanged */
/*      y(*) -- solution at  t */
/*      t -- last point reached in integration.  normal return has */
/*           t = tout . */
/*      tout -- unchanged */
/*      relerr,abserr -- normal return has tolerances unchanged.  iflag=3 */
/*           signals tolerances increased */
/*      iflag = 2 -- normal return.  integration reached  tout */
/*            = 3 -- integration did not reach  tout  because error */
/*                   tolerances too small.  relerr ,  abserr  increased */
/*                   appropriately for continuing */
/*            = 4 -- integration did not reach  tout  because more than */
/*                   500 steps needed */
/*            = 5 -- integration did not reach  tout  because equations */
/*                   appear to be stiff */
/*            = 6 -- invalid input parameters (fatal error) */
/*           the value of  iflag  is returned negative when the input */
/*           value is negative and the integration does not reach  tout , */
/*           i.e., -3, -4, -5. */
/*      work(*),iwork(*) -- information generally of no interest to the */
/*           user but necessary for subsequent calls. */

/*   subsequent calls to  ode -- */

/*   subroutine  ode  returns with all information needed to continue */
/*   the integration.  if the integration reached  tout , the user need */
/*   only define a new  tout  and call again.  if the integration did not */
/*   reach  tout  and the user wants to continue, he just calls again. */
/*   the output value of  iflag  is the appropriate input value for */
/*   subsequent calls.  the only situation in which it should be altered */
/*   is to stop the integration internally at the new  tout , i.e., */
/*   change output  iflag=2  to input  iflag=-2 .  error tolerances may */
/*   be changed by the user before continuing.  all other parameters must */
/*   remain unchanged. */

/* *********************************************************************** */
/* *  subroutines  de  and  step  contain machine dependent constants. * */
/* *  be sure they are set before using  ode .                          * */
/* *********************************************************************** */

    /* Parameter adjustments */
    //--y;
    //--work;
    //--iwork;

    /* Function Body */
    iyy = 100;
    iwt = iyy + neqn;
    ip = iwt + neqn;
    iyp = ip + neqn;
    iypout = iyp + neqn;
    iphi = iypout + neqn;
    if (abs(*iflag) == 1) {
	goto L1;
    }
    start = work[istart] > 0.;
    phase1 = work[iphase] > 0.;
    nornd = iwork[2] != -1;
L1:
    ok_ode_de_((ok_force_function)f, neqn, y, t, tout, relerr, abserr, iflag, &work[iyy], &
	    work[iwt], &work[ip], &work[iyp], &work[iypout], &work[iphi], &
	    work[ialpha], &work[ibeta], &work[isig], &work[iv], &work[iw], &
	    work[ig], &phase1, &work[ipsi], &work[ix], &work[ih], &work[ihold]
	    , &start, &work[itold], &work[idelsn], &iwork[1], &nornd, &iwork[
	    3], &iwork[4], &iwork[5], params, g, rho);
    work[istart] = -1.;
    if (start) {
	work[istart] = 1.;
    }
    work[iphase] = -1.;
    if (phase1) {
	work[iphase] = 1.;
    }
    iwork[2] = -1;
    if (nornd) {
	iwork[2] = 1;
    }
    return 0;
} /* ode_ */

/* Subroutine */ int ok_ode_de_(ok_force_function f, int neqn, doublereal *y, doublereal *t, 
	doublereal *tout, doublereal *relerr, doublereal *abserr, int *
	iflag, doublereal *yy, doublereal *wt, doublereal *p, doublereal *yp, 
	doublereal *ypout, doublereal *phi, doublereal *alpha, doublereal *
	beta, doublereal *sig, doublereal *v, doublereal *w, doublereal *g, 
	logical *phase1, doublereal *psi, doublereal *x, doublereal *h__, 
	doublereal *hold, logical *start, doublereal *told, doublereal *
	delsgn, int *ns, logical *nornd, int *k, int *kold, 
	int *isnold, void* params, doublereal* ginterp, doublereal* rhointerp)
{
    /* Initialized data */

    static int maxnum = 10000;

    /* System generated locals */
    int phi_dim1, phi_offset, i__1;
    doublereal d__1, d__2, d__3, d__4, d__5;

   
    
    /* Local variables */
     int l;
     doublereal del, eps;
     int isn, kle4;
     doublereal tend;
    extern /* Subroutine */ int ok_ode_step_(doublereal *, doublereal *, ok_force_function, 
	    const int, doublereal *, doublereal *, doublereal *, logical *, 
	    doublereal *, int *, int *, logical *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, logical *, int *, logical *, void*);
     logical crash, stiff;
    extern /* Subroutine */ int ok_ode_intrp_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, const int, int *, doublereal *,
	     doublereal *, doublereal*, doublereal*);
    
     doublereal fouru, absdel, abseps, releps;
     int nostep;


/*   ode  merely allocates storage for  de  to relieve the user of the */
/*   inconvenience of a long call list.  consequently  de  is used as */
/*   described in the comments for  ode . */

/*   this code is completely explained and documented in the text, */
/*   computer solution of ordinary differential equations:  the initial */
/*   value problem  by l. f. shampine and m. k. gordon. */


/* *********************************************************************** */
/* *  the only machine dependent constant is based on the machine unit   * */
/* *  roundoff error  u  which is the smallest positive number such that * */
/* *  1.0+u .gt. 1.0 .  u  must be calculated and  fouru=4.0*u  inserted * */
/* *  in the following data statement before using  de .  the routine    * */
/* *  machin  calculates  u .  fouru  and  twou=2.0*u  must also be      * */
/* *  inserted in subroutine  step  before calling  de .                 * */
/*     data fouru/.888d-15/                                              *** */
/* *********************************************************************** */

/*   the constant  maxnum  is the maximum number of steps allowed in one */
/*   call to  de .  the user may change this limit by altering the */
/*   following statement */
    /* Parameter adjustments */
    phi_dim1 = neqn;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    
    --ypout;
    --yp;
    --p;
    --wt;
    --yy;
    --y;
    --alpha;
    --beta;
    --sig;
    --v;
    --w;
    --g;
    --psi;

    /* Function Body */

/*            ***            ***            *** */
/*   test for improper parameters */

    fouru = 4.f * oneu;
    if (neqn < 1) {
	goto L10;
    }
    if (*t == *tout) {
	goto L10;
    }
    if (*relerr < 0. || *abserr < 0.) {
	goto L10;
    }
    eps = max(*relerr,*abserr);
    if (eps <= 0.) {
	goto L10;
    }
    if (*iflag == 0) {
	goto L10;
    }
    isn = ok_i_sign(c__1, *iflag);
    *iflag = abs(*iflag);
    if (*iflag == 1) {
	goto L20;
    }
    if (*t != *told) {
	goto L10;
    }
    if (*iflag >= 2 && *iflag <= 5) {
	goto L20;
    }
L10:
    *iflag = 6;
    return 0;

/*   on each call set interval of integration and counter for number of */
/*   steps.  adjust input error tolerances to define weight vector for */
/*   subroutine  step */

L20:
    del = *tout - *t;
    absdel = abs(del);
    tend = *t + del * 10.;
    if (isn < 0) {
	tend = *tout;
    }
    nostep = 0;
    kle4 = 0;
    stiff = FALSE_;
    releps = *relerr / eps;
    abseps = *abserr / eps;
    if (*iflag == 1) {
	goto L30;
    }
    if (*isnold < 0) {
	goto L30;
    }
    if (*delsgn * del > 0.) {
	goto L50;
    }

/*   on start and restart also set work variables x and yy(*), store the */
/*   direction of integration and initialize the step size */

L30:
    *start = TRUE_;
    *x = *t;
    i__1 = neqn;
    int dd = 0;
    
    
    for (l = 1; l <= i__1; ++l) {
/* L40: */
	yy[l] = y[l];
        dd++;
    }
    *delsgn = ok_d_sign(c_b11, del);
/* Computing MAX */
    d__3 = (d__1 = *tout - *x, abs(d__1)), d__4 = fouru * abs(*x);
    d__2 = max(d__3,d__4);
    d__5 = *tout - *x;
    *h__ = ok_d_sign(d__2, d__5);

/*   if already past output point, interpolate and return */

L50:
    if ((d__1 = *x - *t, abs(d__1)) < absdel) {
	goto L60;
    }
    ok_ode_intrp_(x, &yy[1], tout, &y[1], &ypout[1], neqn, kold, &phi[phi_offset], &
	    psi[1], ginterp, rhointerp);
    *iflag = 2;
    *t = *tout;
    *told = *t;
    *isnold = isn;
    return 0;

/*   if cannot go past output point and sufficiently close, */
/*   extrapolate and return */

L60:
    if (isn > 0 || (d__1 = *tout - *x, abs(d__1)) >= fouru * abs(*x)) {
	goto L80;
    }
    *h__ = *tout - *x;
    ok_force(*x, &yy[1], &yp[1], params);
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
/* L70: */
	y[l] = yy[l] + *h__ * yp[l];
    }
    *iflag = 2;
    *t = *tout;
    *told = *t;
    *isnold = isn;
    return 0;

/*   test for too many steps */

L80:
    if (nostep < maxnum) {
	goto L100;
    }
    *iflag = isn << 2;
    if (stiff) {
	*iflag = isn * 5;
    }
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
/* L90: */
	y[l] = yy[l];
    }
    *t = *x;
    *told = *t;
    *isnold = 1;
    return 0;

/*   limit step size, set weight vector and take a step */

L100:
/* Computing MIN */
    d__3 = abs(*h__), d__4 = (d__1 = tend - *x, abs(d__1));
    d__2 = min(d__3,d__4);
    *h__ = ok_d_sign(d__2, *h__);
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
/* L110: */
	wt[l] = releps * (d__1 = yy[l], abs(d__1)) + abseps;
    }
    ok_ode_step_(x, &yy[1], (ok_force_function)f, neqn, h__, &eps, &wt[1], start, hold, k, kold, &
	    crash, &phi[phi_offset], &p[1], &yp[1], &psi[1], &alpha[1], &beta[
	    1], &sig[1], &v[1], &w[1], &g[1], phase1, ns, nornd, params);

/*   test for tolerances too small */

    if (! crash) {
	goto L130;
    }
    *iflag = isn * 3;
    *relerr = eps * releps;
    *abserr = eps * abseps;
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
/* L120: */
	y[l] = yy[l];
    }
    *t = *x;
    *told = *t;
    *isnold = 1;
    return 0;

/*   augment counter on number of steps and test for stiffness */

L130:
    ++nostep;
    ++kle4;
    if (*kold > 4) {
	kle4 = 0;
    }
    if (kle4 >= 50) {
	stiff = TRUE_;
    }
    goto L50;
} /* de_ */

/* Subroutine */ int ok_ode_step_(doublereal *x, doublereal *y, ok_force_function f, const int
	neqn, doublereal *h__, doublereal *eps, doublereal *wt, logical *
	start, doublereal *hold, int *k, int *kold, logical *crash, 
	doublereal *phi, doublereal *p, doublereal *yp, doublereal *psi, 
	doublereal *alpha, doublereal *beta, doublereal *sig, doublereal *v, 
	doublereal *w, doublereal *g, logical *phase1, int *ns, logical *
	nornd, void* params)
{
    /* Initialized data */

    static doublereal two[13] = { 2.,4.,8.,16.,32.,64.,128.,256.,512.,1024.,
	    2048.,4096.,8192. };
    static doublereal gstr[13] = { .5,.0833,.0417,.0264,.0188,.0143,.0114,
	    .00936,.00789,.00679,.00592,.00524,.00468 };

    /* System generated locals */
    int phi_dim1, phi_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);
    
    /* Local variables */
     int i__, j, l;
     doublereal r__;
     int iq, im1, km1, km2, ip1, kp1, kp2;
     doublereal erk, err, tau, rho, sum;
     int nsm2, nsp1, nsp2;
     doublereal absh, hnew;
     int knew;
     doublereal xold, twou, erkm1, erkm2, erkp1, temp1, temp2, temp3, 
	    temp4, temp5, temp6, p5eps;
     int ifail;
     doublereal reali, round;
    
     doublereal fouru;
     int limit1, limit2;
     doublereal realns;

/*    1  hold,k,kold,crash,phi,p,yp,psi) */

/*   double precision subroutine  step */
/*   integrates a system of first order ordinary */
/*   differential equations one step, normally from x to x+h, using a */
/*   modified divided difference form of the adams pece formulas.  local */
/*   extrapolation is used to improve absolute stability and accuracy. */
/*   the code adjusts its order and step size to control the local error */
/*   per unit step in a generalized sense.  special devices are included */
/*   to control roundoff error and to detect when the user is requesting */
/*   too much accuracy. */

/*   this code is completely explained and documented in the text, */
/*   computer solution of ordinary differential equations:  the initial */
/*   value problem  by l. f. shampine and m. k. gordon. */


/*   the parameters represent: */
/*      x -- independent variable             (real*8) */
/*      y(*) -- solution vector at x          (real*8) */
/*      yp(*) -- derivative of solution vector at  x  after successful */
/*           step                             (real*8) */
/*      neqn -- number of equations to be integrated (int*4) */
/*      h -- appropriate step size for next step.  normally determined by */
/*           code                             (real*8) */
/*      eps -- local error tolerance.  must be variable  (real*8) */
/*      wt(*) -- vector of weights for error criterion   (real*8) */
/*      start -- logical variable set .true. for first step,  .false. */
/*           otherwise                        (logical*4) */
/*      hold -- step size used for last successful step  (real*8) */
/*      k -- appropriate order for next step (determined by code) */
/*      kold -- order used for last successful step */
/*      crash -- logical variable set .true. when no step can be taken, */
/*           .false. otherwise. */
/*   the arrays  phi, psi  are required for the interpolation subroutine */
/*   intrp.  the array p is internal to the code.  all are real*8 */

/*   input to  step */

/*      first call -- */

/*   the user must provide storage in his driver program for all arrays */
/*   in the call list, namely */

/*     dimension y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),psi(12) */

/*   the user must also declare  start  and  crash  logical variables */
/*   and  f  an external subroutine, supply the subroutine  f(x,y,yp) */
/*   to evaluate */
/*      dy(i)/dx = yp(i) = f(x,y(1),y(2),...,y(neqn)) */
/*   and initialize only the following parameters: */
/*      x -- initial value of the independent variable */
/*      y(*) -- vector of initial values of dependent variables */
/*      neqn -- number of equations to be integrated */
/*      h -- nominal step size indicating direction of integration */
/*           and maximum size of step.  must be variable */
/*      eps -- local error tolerance per step.  must be variable */
/*      wt(*) -- vector of non-zero weights for error criterion */
/*      start -- .true. */

/*   step  requires the l2 norm of the vector with components */
/*   local error(l)/wt(l)  be less than  eps  for a successful step.  the */
/*   array  wt  allows the user to specify an error test appropriate */
/*   for his problem.  for example, */
/*      wt(l) = 1.0  specifies absolute error, */
/*            = dabs(y(l))  error relative to the most recent value of */
/*                 the l-th component of the solution, */
/*            = dabs(yp(l))  error relative to the most recent value of */
/*                 the l-th component of the derivative, */
/*            = dmax1(wt(l),dabs(y(l)))  error relative to the largest */
/*                 magnitude of l-th component obtained so far, */
/*            = dabs(y(l))*relerr/eps + abserr/eps  specifies a mixed */
/*                 relative-absolute test where  relerr  is relative */
/*                 error,  abserr  is absolute error and  eps = */
/*                 dmax1(relerr,abserr) . */

/*      subsequent calls -- */

/*   subroutine  step  is designed so that all information needed to */
/*   continue the integration, including the step size  h  and the order */
/*   k , is returned with each step.  with the exception of the step */
/*   size, the error tolerance, and the weights, none of the parameters */
/*   should be altered.  the array  wt  must be updated after each step */
/*   to maintain relative error tests like those above.  normally the */
/*   integration is continued just beyond the desired endpoint and the */
/*   solution interpolated there with subroutine  intrp .  if it is */
/*   impossible to integrate beyond the endpoint, the step size may be */
/*   reduced to hit the endpoint since the code will not take a step */
/*   larger than the  h  input.  changing the direction of integration, */
/*   i.e., the sign of  h , requires the user set  start = .true. before */
/*   calling  step  again.  this is the only situation in which  start */
/*   should be altered. */

/*   output from  step */

/*      successful step -- */

/*   the subroutine returns after each successful step with  start  and */
/*   crash  set .false. .  x  represents the independent variable */
/*   advanced one step of length  hold  from its value on input and  y */
/*   the solution vector at the new value of  x .  all other parameters */
/*   represent information corresponding to the new  x  needed to */
/*   continue the integration. */

/*      unsuccessful step -- */

/*   when the error tolerance is too small for the machine precision, */
/*   the subroutine returns without taking a step and  crash = .true. . */
/*   an appropriate step size and error tolerance for continuing are */
/*   estimated and all other information is restored as upon input */
/*   before returning.  to continue with the larger tolerance, the user */
/*   just calls the code again.  a restart is neither required nor */
/*   desirable. */

/* *********************************************************************** */
/* *  the only machine dependent constants are based on the machine unit * */
/* *  roundoff error  u  which is the smallest positive number such that * */
/* *  1.0+u .gt. 1.0  .  the user must calculate  u  and insert          * */
/* *  twou=2.0*u  and  fouru=4.0*u  in the data statement before calling * */
/* *  the code.  the routine  machin  calculates  u .                    * */
/*     data twou,fouru/.444d-15,.888d-15/                                *** */
/* *********************************************************************** */
    /* Parameter adjustments */
    --yp;
    --p;
    phi_dim1 = neqn;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --wt;
    --y;
    --psi;
    --alpha;
    --beta;
    --sig;
    --v;
    --w;
    --g;

    /* Function Body */
/*     data g(1),g(2)/1.0,0.5/,sig(1)/1.0/ */


    twou = 2.f * oneu;
    fouru = twou * 2.f;
/*       ***     begin block 0     *** */
/*   check if step size or error tolerance is too small for machine */
/*   precision.  if first step, initialize phi array and estimate a */
/*   starting step size. */
/*                   *** */

/*   if step size is too small, determine an acceptable one */

    *crash = TRUE_;
    if (abs(*h__) >= fouru * abs(*x)) {
	goto L5;
    }
    d__1 = fouru * abs(*x);
    *h__ = ok_d_sign(d__1, *h__);
    return 0;
L5:
    p5eps = *eps * .5;

/*   if error tolerance is too small, increase it to an acceptable value */

    round = 0.;
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
/* L10: */
/* Computing 2nd power */
	d__1 = y[l] / wt[l];
	round += d__1 * d__1;
    }
    round = twou * sqrt(round);
    if (p5eps >= round) {
	goto L15;
    }
    *eps = round * 2.f * (fouru + 1.);
    return 0;
L15:
    *crash = FALSE_;
    g[1] = 1.;
    g[2] = .5;
    sig[1] = 1.;
    if (! (*start)) {
	goto L99;
    }

/*   initialize.  compute appropriate step size for first step */

    ok_force(*x, &y[1], &yp[1], params);
    sum = 0.;
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
	phi[l + phi_dim1] = yp[l];
	phi[l + (phi_dim1 << 1)] = 0.;
/* L20: */
/* Computing 2nd power */
	d__1 = yp[l] / wt[l];
	sum += d__1 * d__1;
    }
    sum = sqrt(sum);
    absh = abs(*h__);
    if (*eps < sum * 16. * *h__ * *h__) {
	absh = sqrt(*eps / sum) * .25;
    }
/* Computing MAX */
    d__2 = absh, d__3 = fouru * abs(*x);
    d__1 = max(d__2,d__3);
    *h__ = ok_d_sign(d__1, *h__);
    *hold = 0.;
    *k = 1;
    *kold = 0;
    *start = FALSE_;
    *phase1 = TRUE_;
    *nornd = TRUE_;
    if (p5eps > round * 100.) {
	goto L99;
    }
    *nornd = FALSE_;
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
/* L25: */
	phi[l + phi_dim1 * 15] = 0.;
    }
L99:
    ifail = 0;
/*       ***     end block 0     *** */

/*       ***     begin block 1     *** */
/*   compute coefficients of formulas for this step.  avoid computing */
/*   those quantities not changed when step size is not changed. */
/*                   *** */

L100:
    kp1 = *k + 1;
    kp2 = *k + 2;
    km1 = *k - 1;
    km2 = *k - 2;

/*   ns is the number of steps taken with size h, including the current */
/*   one.  when k.lt.ns, no coefficients change */

    if (*h__ != *hold) {
	*ns = 0;
    }
    if (*ns <= *kold) {
	++(*ns);
    }
    nsp1 = *ns + 1;
    if (*k < *ns) {
	goto L199;
    }

/*   compute those components of alpha(*),beta(*),psi(*),sig(*) which */
/*   are changed */

    beta[*ns] = 1.;
    realns = (doublereal) (*ns);
    alpha[*ns] = 1. / realns;
    temp1 = *h__ * realns;
    sig[nsp1] = 1.;
    if (*k < nsp1) {
	goto L110;
    }
    i__1 = *k;
    for (i__ = nsp1; i__ <= i__1; ++i__) {
	im1 = i__ - 1;
	temp2 = psi[im1];
	psi[im1] = temp1;
	beta[i__] = beta[im1] * psi[im1] / temp2;
	temp1 = temp2 + *h__;
	alpha[i__] = *h__ / temp1;
	reali = (doublereal) i__;
/* L105: */
	sig[i__ + 1] = reali * alpha[i__] * sig[i__];
    }
L110:
    psi[*k] = temp1;

/*   compute coefficients g(*) */

/*   initialize v(*) and set w(*).  g(2) is set in data statement */

    if (*ns > 1) {
	goto L120;
    }
    i__1 = *k;
    for (iq = 1; iq <= i__1; ++iq) {
	temp3 = (doublereal) (iq * (iq + 1));
	v[iq] = 1. / temp3;
/* L115: */
	w[iq] = v[iq];
    }
    goto L140;

/*   if order was raised, update diagonal part of v(*) */

L120:
    if (*k <= *kold) {
	goto L130;
    }
    temp4 = (doublereal) (*k * kp1);
    v[*k] = 1. / temp4;
    nsm2 = *ns - 2;
    if (nsm2 < 1) {
	goto L130;
    }
    i__1 = nsm2;
    for (j = 1; j <= i__1; ++j) {
	i__ = *k - j;
/* L125: */
	v[i__] -= alpha[j + 1] * v[i__ + 1];
    }

/*   update v(*) and set w(*) */

L130:
    limit1 = kp1 - *ns;
    temp5 = alpha[*ns];
    i__1 = limit1;
    for (iq = 1; iq <= i__1; ++iq) {
	v[iq] -= temp5 * v[iq + 1];
/* L135: */
	w[iq] = v[iq];
    }
    g[nsp1] = w[1];

/*   compute the g(*) in the work vector w(*) */

L140:
    nsp2 = *ns + 2;
    if (kp1 < nsp2) {
	goto L199;
    }
    i__1 = kp1;
    for (i__ = nsp2; i__ <= i__1; ++i__) {
	limit2 = kp2 - i__;
	temp6 = alpha[i__ - 1];
	i__2 = limit2;
	for (iq = 1; iq <= i__2; ++iq) {
/* L145: */
	    w[iq] -= temp6 * w[iq + 1];
	}
/* L150: */
	g[i__] = w[1];
    }
L199:
/*       ***     end block 1     *** */

/*       ***     begin block 2     *** */
/*   predict a solution p(*), evaluate derivatives using predicted */
/*   solution, estimate local error at order k and errors at orders k, */
/*   k-1, k-2 as if constant step size were used. */
/*                   *** */

/*   change phi to phi star */

    if (*k < nsp1) {
	goto L215;
    }
    i__1 = *k;
    for (i__ = nsp1; i__ <= i__1; ++i__) {
	temp1 = beta[i__];
	i__2 = neqn;
	for (l = 1; l <= i__2; ++l) {
/* L205: */
	    phi[l + i__ * phi_dim1] = temp1 * phi[l + i__ * phi_dim1];
	}
/* L210: */
    }

/*   predict solution and differences */

L215:
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
	phi[l + kp2 * phi_dim1] = phi[l + kp1 * phi_dim1];
	phi[l + kp1 * phi_dim1] = 0.;
/* L220: */
	p[l] = 0.;
    }
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	i__ = kp1 - j;
	ip1 = i__ + 1;
	temp2 = g[i__];
	i__2 = neqn;
	for (l = 1; l <= i__2; ++l) {
	    p[l] += temp2 * phi[l + i__ * phi_dim1];
/* L225: */
	    phi[l + i__ * phi_dim1] += phi[l + ip1 * phi_dim1];
	}
/* L230: */
    }
    if (*nornd) {
	goto L240;
    }
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
	tau = *h__ * p[l] - phi[l + phi_dim1 * 15];
	p[l] = y[l] + tau;
/* L235: */
	phi[l + (phi_dim1 << 4)] = p[l] - y[l] - tau;
    }
    goto L250;
L240:
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
/* L245: */
	p[l] = y[l] + *h__ * p[l];
    }
L250:
    xold = *x;
    *x += *h__;
    absh = abs(*h__);
    ok_force(*x, &p[1], &yp[1], params);

/*   estimate errors at orders k,k-1,k-2 */

    erkm2 = 0.;
    erkm1 = 0.;
    erk = 0.;
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
	temp3 = 1. / wt[l];
	temp4 = yp[l] - phi[l + phi_dim1];
	if (km2 < 0) {
	    goto L265;
	} else if (km2 == 0) {
	    goto L260;
	} else {
	    goto L255;
	}
L255:
/* Computing 2nd power */
	d__1 = (phi[l + km1 * phi_dim1] + temp4) * temp3;
	erkm2 += d__1 * d__1;
L260:
/* Computing 2nd power */
	d__1 = (phi[l + *k * phi_dim1] + temp4) * temp3;
	erkm1 += d__1 * d__1;
L265:
/* Computing 2nd power */
	d__1 = temp4 * temp3;
	erk += d__1 * d__1;
    }
    if (km2 < 0) {
	goto L280;
    } else if (km2 == 0) {
	goto L275;
    } else {
	goto L270;
    }
L270:
    erkm2 = absh * sig[km1] * gstr[km2 - 1] * sqrt(erkm2);
L275:
    erkm1 = absh * sig[*k] * gstr[km1 - 1] * sqrt(erkm1);
L280:
    temp5 = absh * sqrt(erk);
    err = temp5 * (g[*k] - g[kp1]);
    erk = temp5 * sig[kp1] * gstr[*k - 1];
    knew = *k;

/*   test if order should be lowered */

    if (km2 < 0) {
	goto L299;
    } else if (km2 == 0) {
	goto L290;
    } else {
	goto L285;
    }
L285:
    if (max(erkm1,erkm2) <= erk) {
	knew = km1;
    }
    goto L299;
L290:
    if (erkm1 <= erk * .5) {
	knew = km1;
    }

/*   test if step successful */

L299:
    if (err <= *eps) {
	goto L400;
    }
/*       ***     end block 2     *** */

/*       ***     begin block 3     *** */
/*   the step is unsuccessful.  restore  x, phi(*,*), psi(*) . */
/*   if third consecutive failure, set order to one.  if step fails more */
/*   than three times, consider an optimal step size.  double error */
/*   tolerance and return if estimated step size is too small for machine */
/*   precision. */
/*                   *** */

/*   restore x, phi(*,*) and psi(*) */

    *phase1 = FALSE_;
    *x = xold;
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp1 = 1. / beta[i__];
	ip1 = i__ + 1;
	i__2 = neqn;
	for (l = 1; l <= i__2; ++l) {
/* L305: */
	    phi[l + i__ * phi_dim1] = temp1 * (phi[l + i__ * phi_dim1] - phi[
		    l + ip1 * phi_dim1]);
	}
/* L310: */
    }
    if (*k < 2) {
	goto L320;
    }
    i__1 = *k;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L315: */
	psi[i__ - 1] = psi[i__] - *h__;
    }

/*   on third failure, set order to one.  thereafter, use optimal step */
/*   size */

L320:
    ++ifail;
    temp2 = .5;
    if ((i__1 = ifail - 3) < 0) {
	goto L335;
    } else if (i__1 == 0) {
	goto L330;
    } else {
	goto L325;
    }
L325:
    if (p5eps < erk * .25) {
	temp2 = sqrt(p5eps / erk);
    }
L330:
    knew = 1;
L335:
    *h__ = temp2 * *h__;
    *k = knew;
    if (abs(*h__) >= fouru * abs(*x)) {
	goto L340;
    }
    *crash = TRUE_;
    d__1 = fouru * abs(*x);
    *h__ = ok_d_sign(d__1, *h__);
    *eps += *eps;
    return 0;
L340:
    goto L100;
/*       ***     end block 3     *** */

/*       ***     begin block 4     *** */
/*   the step is successful.  correct the predicted solution, evaluate */
/*   the derivatives using the corrected solution and update the */
/*   differences.  determine best order and step size for next step. */
/*                   *** */
L400:
    *kold = *k;
    *hold = *h__;

/*   correct and evaluate */

    temp1 = *h__ * g[kp1];
    if (*nornd) {
	goto L410;
    }
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
	rho = temp1 * (yp[l] - phi[l + phi_dim1]) - phi[l + (phi_dim1 << 4)];
	y[l] = p[l] + rho;
/* L405: */
	phi[l + phi_dim1 * 15] = y[l] - p[l] - rho;
    }
    goto L420;
L410:
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
/* L415: */
	y[l] = p[l] + temp1 * (yp[l] - phi[l + phi_dim1]);
    }
L420:
    ok_force(*x, &y[1], &yp[1], params);

/*   update differences for next step */

    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
	phi[l + kp1 * phi_dim1] = yp[l] - phi[l + phi_dim1];
/* L425: */
	phi[l + kp2 * phi_dim1] = phi[l + kp1 * phi_dim1] - phi[l + kp2 * 
		phi_dim1];
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = neqn;
	for (l = 1; l <= i__2; ++l) {
/* L430: */
	    phi[l + i__ * phi_dim1] += phi[l + kp1 * phi_dim1];
	}
/* L435: */
    }

/*   estimate error at order k+1 unless: */
/*     in first phase when always raise order, */
/*     already decided to lower order, */
/*     step size not constant so estimate unreliable */

    erkp1 = 0.;
    if (knew == km1 || *k == 12) {
	*phase1 = FALSE_;
    }
    if (*phase1) {
	goto L450;
    }
    if (knew == km1) {
	goto L455;
    }
    if (kp1 > *ns) {
	goto L460;
    }
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
/* L440: */
/* Computing 2nd power */
	d__1 = phi[l + kp2 * phi_dim1] / wt[l];
	erkp1 += d__1 * d__1;
    }
    erkp1 = absh * gstr[kp1 - 1] * sqrt(erkp1);

/*   using estimated error at order k+1, determine appropriate order */
/*   for next step */

    if (*k > 1) {
	goto L445;
    }
    if (erkp1 >= erk * .5) {
	goto L460;
    }
    goto L450;
L445:
    if (erkm1 <= min(erk,erkp1)) {
	goto L455;
    }
    if (erkp1 >= erk || *k == 12) {
	goto L460;
    }

/*   here erkp1 .lt. erk .lt. dmax1(erkm1,erkm2) else order would have */
/*   been lowered in block 2.  thus order is to be raised */

/*   raise order */

L450:
    *k = kp1;
    erk = erkp1;
    goto L460;

/*   lower order */

L455:
    *k = km1;
    erk = erkm1;

/*   with new order determine appropriate step size for next step */

L460:
    hnew = *h__ + *h__;
    if (*phase1) {
	goto L465;
    }
    if (p5eps >= erk * two[*k]) {
	goto L465;
    }
    hnew = *h__;
    if (p5eps >= erk) {
	goto L465;
    }
    temp2 = (doublereal) (*k + 1);
    d__1 = p5eps / erk;
    d__2 = 1. / temp2;
    r__ = pow(d__1, d__2);
/* Computing MAX */
    d__1 = .5, d__2 = min(.9,r__);
    hnew = absh * max(d__1,d__2);
/* Computing MAX */
    d__2 = hnew, d__3 = fouru * abs(*x);
    d__1 = max(d__2,d__3);
    hnew = ok_d_sign(d__1, *h__);
L465:
    *h__ = hnew;
    return 0;
/*       ***     end block 4     *** */
} /* step_ */

/* Subroutine */ int ok_ode_intrp_(doublereal *x, doublereal *y, doublereal *xout, 
	doublereal *yout, doublereal *ypout, const int neqn, int *kold, 
	doublereal *phi, doublereal *psi, doublereal* rho, doublereal* g)
{
    /* Initialized data */

    
    /* System generated locals */
    int phi_dim1, phi_offset, i__1, i__2;

    /* Local variables */
     int i__, j, l;
     doublereal w[13], hi;
     int ki, jm1;
     doublereal eta;
     int kip1;
     doublereal term, temp1, temp2, temp3, gamma;
     int limit1;
     doublereal psijm1;


/*   the methods in subroutine  step  approximate the solution near  x */
/*   by a polynomial.  subroutine  intrp  approximates the solution at */
/*   xout  by evaluating the polynomial there.  information defining this */
/*   polynomial is passed from  step  so  intrp  cannot be used alone. */

/*   this code is completely explained and documented in the text, */
/*   computer solution of ordinary differential equations:  the initial */
/*   value problem  by l. f. shampine and m. k. gordon. */

/*   input to intrp -- */

/*   all floating point variables are double precision */
/*   the user provides storage in the calling program for the arrays in */
/*   the call list */
/*   and defines */
/*      xout -- point at which solution is desired. */
/*   the remaining parameters are defined in  step  and passed to  intrp */
/*   from that subroutine */

/*   output from  intrp -- */

/*      yout(*) -- solution at  xout */
/*      ypout(*) -- derivative of solution at  xout */
/*   the remaining parameters are returned unaltered from their input */
/*   values.  integration with  step  may be continued. */

    /* Parameter adjustments */
    phi_dim1 = neqn;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --ypout;
    --yout;
    --y;
    --psi;

    /* Function Body */

    hi = *xout - *x;
    ki = *kold + 1;
    kip1 = ki + 1;

/*   initialize w(*) for computing g(*) */

    i__1 = ki;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp1 = (doublereal) i__;
/* L5: */
	w[i__ - 1] = 1. / temp1;
    }
    term = 0.;

/*   compute g(*) */

    i__1 = ki;
    for (j = 2; j <= i__1; ++j) {
	jm1 = j - 1;
	psijm1 = psi[jm1];
	gamma = (hi + term) / psijm1;
	eta = hi / psijm1;
	limit1 = kip1 - j;
	i__2 = limit1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L10: */
	    w[i__ - 1] = gamma * w[i__ - 1] - eta * w[i__];
	}
	g[j - 1] = w[0];
	rho[j - 1] = gamma * rho[jm1 - 1];
/* L15: */
	term = psijm1;
    }

/*   interpolate */

    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
	ypout[l] = 0.;
/* L20: */
	yout[l] = 0.;
    }
    i__1 = ki;
    for (j = 1; j <= i__1; ++j) {
	i__ = kip1 - j;
	temp2 = g[i__ - 1];
	temp3 = rho[i__ - 1];
	i__2 = neqn;
	for (l = 1; l <= i__2; ++l) {
	    yout[l] += temp2 * phi[l + i__ * phi_dim1];
/* L25: */
	    ypout[l] += temp3 * phi[l + i__ * phi_dim1];
	}
/* L30: */
    }
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
/* L35: */
	yout[l] = y[l] + hi * yout[l];
    }
    return 0;
} /* intrp_ */


ok_system** ok_integrate_ode(ok_system* initial, const gsl_vector* times, ok_integrator_options* options,
        ok_system** bag, int* error) {
    
    // Check that the system has been set-up
    assert(initial->xyz != NULL);
    // Check input arguments
    assert(times != NULL);
    
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
    
    int DIMENSIONS = NDIMS * 7;
    
    // Allocates the temporary buffer to hold cartesian coordinates
    double prevTime = startTime;
    gsl_matrix* prevOrbits = initial->orbits;
    
    double h = 0.01;
    gsl_matrix* xyz = gsl_matrix_alloc(initial->nplanets+1, 7);
    MATRIX_MEMCPY(xyz, initial->xyz);
    
    doublereal relerr = options->rel_acc * 100;
    doublereal abserr = options->abs_acc * 100;
    doublereal work[100 + 21 * DIMENSIONS + 1];
    
    int iwork[6];
    void* params = (void*) initial;
    
    // Loop through the times vector
    for (int i = 0; i < SAMPLES; i++) {
        double time = times->data[i];
        int iflag = 1;

        // Integrate between prevTime and time
        if (fabs(time - prevTime) > 1e-10) {
            while (fabs(time - prevTime) > 1e-10) {
                h = copysign(h, time-prevTime);
                ok_ode_ode_(&ok_force, DIMENSIONS, xyz->data, &prevTime, &time,
                    &relerr, &abserr, &iflag, work, iwork, params);
                
                if (iflag != 2) {
                    for (int i = 0; i < SAMPLES; i++)
                        ok_free_system(bag[i]);
                    free(bag);
                                        
                    return NULL;
                }
            }
        } else {
            if (i == 0)
                MATRIX_MEMCPY(xyz, initial->xyz);
            else
                MATRIX_MEMCPY(xyz, bag[i-1]->xyz);
        };
        
        
        // Set up the return vector
        bag[i]->time = bag[i]->epoch = time;
       
        MATRIX_MEMCPY(bag[i]->xyz, xyz);
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
    
    
    gsl_matrix_free(xyz);
    return bag;
}