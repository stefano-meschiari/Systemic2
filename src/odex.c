/* odex.f -- translated by f2c (version 20100827).
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
#include "gsl/gsl_vector.h"
#include "gsl/gsl_vector_int.h"
/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;
static doublereal c_b67 = 1.;


int ok_odex_odxcor_(integer *n, S_fp fcn, doublereal *x, doublereal *
	y, doublereal *xend, doublereal *hmax, doublereal *h__, doublereal *
	rtol, doublereal *atol, integer *itol, integer *km, 
	integer *iout, integer *idid, integer *nmax, doublereal *uround, 
	doublereal *dy, doublereal *yh1, doublereal *yh2, doublereal *dz, 
	doublereal *scal, doublereal *fsafe, doublereal *ysafe, doublereal *t,
	 doublereal *hh, doublereal *w, doublereal *a, doublereal *dens, 
	integer *ncom, integer *icomp, integer *nj, integer *ipoint, integer *
	nsequ, integer *mstab, integer *jstab, integer *lfsafe, doublereal *
	safe1, doublereal *safe2, doublereal *safe3, doublereal *fac1, 
	doublereal *fac2, doublereal *fac3, doublereal *fac4, integer *iderr, 
	doublereal *errfac, integer *mudif, integer *nrd, integer *nfcn, integer *nstep, integer *naccpt, 
	integer *nrejct, void* params);

int ok_odex_midex_(integer *j, doublereal *x, doublereal *y, 
	doublereal *h__, doublereal *hmax, integer *n, S_fp fcn, doublereal *
	dy, doublereal *yh1, doublereal *yh2, doublereal *dz, doublereal *t, 
	integer *nj, doublereal *hh, doublereal *w, doublereal *err, 
	doublereal *fac, doublereal *a, doublereal *safe1, doublereal *uround,
	 doublereal *fac1, doublereal *fac2, doublereal *safe2, doublereal *
	scal, logical *atov, doublereal *safe3, logical *reject, integer *km, 
	doublereal *rtol, doublereal *atol, integer *itol, integer *mstab, 
	integer *jstab, doublereal *errold, doublereal *fsafe, integer *
	lfsafe, integer *iout, integer *ipt, doublereal *ysafe, integer *
	icomp, integer *nrd, integer *nfcn, void* params);

int ok_odex_interp_(integer n, doublereal *y, integer imit);

/* Subroutine */ int ok_odex_(integer *n, U_fp fcn, doublereal *x, doublereal *y,
	 doublereal *xend, doublereal *h__, doublereal *rtol, doublereal *
	atol, integer *itol, integer *iout, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *idid, void* params)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    integer i__, km, iea, iet, nrd, iew;
    doublereal fac1, fac2, fac3, fac4;
    integer iehh, ieco, nfcn, ienj, iefs, icom, ieip, iedy, iedz;
    doublereal hmax;
    integer ncom, nmax, ieys;
    doublereal safe1, safe2, safe3;
    integer ieyh1, ieyh2, iefac, jstab, mudif, iderr, mstab;
    logical arret;
    integer nstep, nsequ, lfsafe, iescal, naccpt, nrejct, nrdens;
    
    integer istore;
    doublereal uround;

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 6, 0, 0, 0 };
    static cilist io___9 = { 0, 6, 0, 0, 0 };
    static cilist io___11 = { 0, 6, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, 0, 0 };
    static cilist io___16 = { 0, 6, 0, 0, 0 };
    static cilist io___18 = { 0, 6, 0, 0, 0 };
    static cilist io___20 = { 0, 6, 0, 0, 0 };
    static cilist io___23 = { 0, 6, 0, 0, 0 };
    static cilist io___26 = { 0, 6, 0, 0, 0 };
    static cilist io___48 = { 0, 6, 0, 0, 0 };
    static cilist io___52 = { 0, 6, 0, 0, 0 };


/* ---------------------------------------------------------- */
/*     NUMERICAL SOLUTION OF A SYSTEM OF FIRST 0RDER */
/*     ORDINARY DIFFERENTIAL EQUATIONS  Y'=F(X,Y). */
/*     THIS IS AN EXTRAPOLATION-ALGORITHM (GBS), BASED ON THE */
/*     EXPLICIT MIDPOINT RULE (WITH STEPSIZE CONTROL, */
/*     ORDER SELECTION AND DENSE OUTPUT). */

/*     AUTHORS: E. HAIRER AND G. WANNER */
/*              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES */
/*              CH-1211 GENEVE 24, SWITZERLAND */
/*              E-MAIL:  Ernst.Hairer@math.unige.ch */
/*                       Gerhard.Wanner@math.unige.ch */
/*              DENSE OUTPUT WRITTEN BY E. HAIRER AND A. OSTERMANN */

/*     THIS CODE IS DESCRIBED IN SECTION II.9 OF THE BOOK: */
/*         E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY */
/*         DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION. */
/*         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS, */
/*         SPRINGER-VERLAG (1993) */

/*     VERSION SEPTEMBER 30, 1995 */
/*         SMALL CORRECTIONS ON OCTOBER 11, 2009 */

/*     INPUT PARAMETERS */
/*     ---------------- */
/*     N           DIMENSION OF THE SYSTEM */

/*     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE */
/*                 VALUE OF F(X,Y): */
/*                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR) */
/*                    DOUBLE PRECISION X,Y(N),F(N) */
/*                    F(1)=...   ETC. */

/*     X           INITIAL X-VALUE */

/*     Y(N)        INITIAL VALUES FOR Y */

/*     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE) */

/*     H           INITIAL STEP SIZE GUESS; */
/*                 H=1.D0/(NORM OF F'), USUALLY 1.D-1 OR 1.D-3, IS GOOD. */
/*                 THIS CHOICE IS NOT VERY IMPORTANT, THE CODE QUICKLY */
/*                 ADAPTS ITS STEP SIZE. WHEN YOU ARE NOT SURE, THEN */
/*                 STUDY THE CHOSEN VALUES FOR A FEW */
/*                 STEPS IN SUBROUTINE "SOLOUT". */
/*                 (IF H=0.D0, THE CODE PUTS H=1.D-4). */

/*     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY */
/*                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N. */

/*     ITOL        SWITCH FOR RTOL AND ATOL: */
/*                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS. */
/*                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF */
/*                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL */
/*                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS. */
/*                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW */
/*                     RTOL(I)*ABS(Y(I))+ATOL(I). */

/*     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE */
/*                 NUMERICAL SOLUTION DURING INTEGRATION. */
/*                 IF IOUT.GE.1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP. */
/*                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0. */
/*                 IT MUST HAVE THE FORM */
/*                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,CON,NCON,ICOMP,ND, */
/*                                       RPAR,IPAR,IRTRN) */
/*                    DIMENSION X,Y(N),CON(NCON),ICOMP(ND) */
/*                    .... */
/*                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH */
/*                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS */
/*                    THE FIRST GRID-POINT). */
/*                 "XOLD" IS THE PRECEEDING GRID-POINT. */
/*                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN */
/*                    IS SET <0, ODEX WILL RETURN TO THE CALLING PROGRAM. */

/*          -----  CONTINUOUS OUTPUT (IF IOUT=2): ----- */
/*                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION */
/*                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH */
/*                 THE DOUBLE PRECISION FUNCTION */
/*                    >>>   CONTEX(I,S,CON,NCON,ICOMP,ND)   <<< */
/*                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH */
/*                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE */
/*                 S SHOULD LIE IN THE INTERVAL [XOLD,X]. */

/*     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT: */
/*                    IOUT=0: SUBROUTINE IS NEVER CALLED */
/*                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT */
/*                    IOUT=2: DENSE OUTPUT IS PERFORMED IN SOLOUT */

/*     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK". */
/*                 SERVES AS WORKING SPACE FOR ALL VECTORS. */
/*                 "LWORK" MUST BE AT LEAST */
/*                    N*(KM+5)+5*KM+20+(2*KM*(KM+2)+5)*NRDENS */
/*                 WHERE NRDENS=IWORK(8) (SEE BELOW) AND */
/*                        KM=9                IF IWORK(2)=0 */
/*                        KM=IWORK(2)         IF IWORK(2).GT.0 */
/*                 WORK(1),...,WORK(20) SERVE AS PARAMETERS */
/*                 FOR THE CODE. FOR STANDARD USE, SET THESE */
/*                 PARAMETERS TO ZERO BEFORE CALLING. */

/*     LWORK       DECLARED LENGTH OF ARRAY "WORK". */

/*     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK". */
/*                 "LIWORK" MUST BE AT LEAST */
/*                               2*KM+21+NRDENS */
/*                 IWORK(1),...,IWORK(20) SERVE AS PARAMETERS */
/*                 FOR THE CODE. FOR STANDARD USE, SET THESE */
/*                 PARAMETERS TO ZERO BEFORE CALLING. */

/*     LIWORK      DECLARED LENGTH OF ARRAY "IWORK". */

/*     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH */
/*                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING */
/*                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES. */

/* ----------------------------------------------------------------------- */

/*     SOPHISTICATED SETTING OF PARAMETERS */
/*     ----------------------------------- */
/*              SEVERAL PARAMETERS (WORK(1),...,IWORK(1),...) ALLOW */
/*              TO ADAPT THE CODE TO THE PROBLEM AND TO THE NEEDS OF */
/*              THE USER. FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES. */

/*    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 2.3D-16. */

/*    WORK(2)   MAXIMAL STEP SIZE, DEFAULT XEND-X. */

/*    WORK(3)   STEP SIZE IS REDUCED BY FACTOR WORK(3), IF THE */
/*              STABILITY CHECK IS NEGATIVE, DEFAULT 0.5. */

/*    WORK(4), WORK(5)   PARAMETERS FOR STEP SIZE SELECTION */
/*              THE NEW STEP SIZE FOR THE J-TH DIAGONAL ENTRY IS */
/*              CHOSEN SUBJECT TO THE RESTRICTION */
/*                 FACMIN/WORK(5) <= HNEW(J)/HOLD <= 1/FACMIN */
/*              WHERE FACMIN=WORK(4)**(1/(2*J-1)) */
/*              DEFAULT VALUES: WORK(4)=0.02D0, WORK(5)=4.D0 */

/*    WORK(6), WORK(7)   PARAMETERS FOR THE ORDER SELECTION */
/*              STEP SIZE IS DECREASED IF    W(K-1) <= W(K)*WORK(6) */
/*              STEP SIZE IS INCREASED IF    W(K) <= W(K-1)*WORK(7) */
/*              DEFAULT VALUES: WORK(6)=0.8D0, WORK(7)=0.9D0 */

/*    WORK(8), WORK(9)   SAFETY FACTORS FOR STEP CONTROL ALGORITHM */
/*             HNEW=H*WORK(9)*(WORK(8)*TOL/ERR)**(1/(J-1)) */
/*             DEFAULT VALUES: WORK(8)=0.65D0, */
/*                        WORK(9)=0.94D0  IF "HOPE FOR CONVERGENCE" */
/*                        WORK(9)=0.90D0  IF "NO HOPE FOR CONVERGENCE" */

/*    IWORK(1)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS. */
/*              THE DEFAULT VALUE (FOR IWORK(1)=0) IS 10000. */

/*    IWORK(2)  THE MAXIMUM NUMBER OF COLUMNS IN THE EXTRAPOLATION */
/*              TABLE. THE DEFAULT VALUE (FOR IWORK(2)=0) IS 9. */
/*              IF IWORK(2).NE.0 THEN IWORK(2) SHOULD BE .GE.3. */

/*    IWORK(3)  SWITCH FOR THE STEP SIZE SEQUENCE (EVEN NUMBERS ONLY) */
/*              IF IWORK(3).EQ.1 THEN 2,4,6,8,10,12,14,16,... */
/*              IF IWORK(3).EQ.2 THEN 2,4,8,12,16,20,24,28,... */
/*              IF IWORK(3).EQ.3 THEN 2,4,6,8,12,16,24,32,... */
/*              IF IWORK(3).EQ.4 THEN 2,6,10,14,18,22,26,30,... */
/*              IF IWORK(3).EQ.5 THEN 4,8,12,16,20,24,28,32,... */
/*              THE DEFAULT VALUE IS IWORK(3)=1 IF IOUT.LE.1; */
/*              THE DEFAULT VALUE IS IWORK(3)=4 IF IOUT.GE.2. */

/*    IWORK(4)  STABILITY CHECK IS ACTIVATED AT MOST IWORK(4) TIMES IN */
/*              ONE LINE OF THE EXTRAP. TABLE, DEFAULT IWORK(4)=1. */

/*    IWORK(5)  STABILITY CHECK IS ACTIVATED ONLY IN THE LINES */
/*              1 TO IWORK(5) OF THE EXTRAP. TABLE, DEFAULT IWORK(5)=1. */

/*    IWORK(6)  IF  IWORK(6)=0  ERROR ESTIMATOR IN THE DENSE */
/*              OUTPUT FORMULA IS ACTIVATED. IT CAN BE SUPPRESSED */
/*              BY PUTTING IWORK(6)=1. */
/*              DEFAULT IWORK(6)=0  (IF IOUT.GE.2). */

/*    IWORK(7)  DETERMINES THE DEGREE OF INTERPOLATION FORMULA */
/*              MU = 2 * KAPPA - IWORK(7) + 1 */
/*              IWORK(7) SHOULD LIE BETWEEN 1 AND 6 */
/*              DEFAULT IWORK(7)=4  (IF IWORK(7)=0). */

/*    IWORK(8)  = NRDENS = NUMBER OF COMPONENTS, FOR WHICH DENSE OUTPUT */
/*              IS REQUIRED */

/*    IWORK(21),...,IWORK(NRDENS+20) INDICATE THE COMPONENTS, FOR WHICH */
/*              DENSE OUTPUT IS REQUIRED */

/* ----------------------------------------------------------------------C */
/*     OUTPUT PARAMETERS */
/*     ----------------- */
/*     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED */
/*                 (AFTER SUCCESSFUL RETURN X=XEND). */

/*     Y(N)        NUMERICAL SOLUTION AT X */

/*     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP */

/*     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN: */
/*                   IDID=1  COMPUTATION SUCCESSFUL, */
/*                   IDID=-1 COMPUTATION UNSUCCESSFUL. */

/*   IWORK(17)  NFCN    NUMBER OF FUNCTION EVALUATIONS */
/*   IWORK(18)  NSTEP   NUMBER OF COMPUTED STEPS */
/*   IWORK(19)  NACCPT  NUMBER OF ACCEPTED STEPS */
/*   IWORK(20)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST), */
/*                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED) */
/* ----------------------------------------------------------------------- */
/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
/*          DECLARATIONS */
/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
/* *** *** *** *** *** *** *** */
/*        SETTING THE PARAMETERS */
/* *** *** *** *** *** *** *** */
    /* Parameter adjustments */
    --y;
    --rtol;
    --atol;
    --work;
    --iwork;

    /* Function Body */
    nfcn = 0;
    nstep = 0;
    naccpt = 0;
    nrejct = 0;
    arret = FALSE_;
/* -------- NMAX , THE MAXIMAL NUMBER OF STEPS ----- */
    if (iwork[1] == 0) {
	nmax = 10000;
    } else {
	nmax = iwork[1];
	if (nmax <= 0) {
	    s_wsle(&io___7);
	    do_lio(&c__9, &c__1, " WRONG INPUT IWORK(1)=", (ftnlen)22);
	    do_lio(&c__3, &c__1, (char *)&iwork[1], (ftnlen)sizeof(integer));
	    e_wsle();
	    arret = TRUE_;
	}
    }
/* -------- KM     MAXIMUM NUMBER OF COLUMNS IN THE EXTRAPOLATION */
    if (iwork[2] == 0) {
	km = 9;
    } else {
	km = iwork[2];
	if (km <= 2) {
	    s_wsle(&io___9);
	    do_lio(&c__9, &c__1, " CURIOUS INPUT IWORK(2)=", (ftnlen)24);
	    do_lio(&c__3, &c__1, (char *)&iwork[2], (ftnlen)sizeof(integer));
	    e_wsle();
	    arret = TRUE_;
	}
    }
/* -------- NSEQU     CHOICE OF STEP SIZE SEQUENCE */
    nsequ = iwork[3];
    if (iwork[3] == 0 && *iout <= 1) {
	nsequ = 1;
    }
    if (iwork[3] == 0 && *iout >= 2) {
	nsequ = 4;
    }
    if (nsequ <= 0 || nsequ >= 6) {
	s_wsle(&io___11);
	do_lio(&c__9, &c__1, " CURIOUS INPUT IWORK(3)=", (ftnlen)24);
	do_lio(&c__3, &c__1, (char *)&iwork[3], (ftnlen)sizeof(integer));
	e_wsle();
	arret = TRUE_;
    }
    if (nsequ <= 3 && *iout >= 2) {
	s_wsle(&io___12);
	do_lio(&c__9, &c__1, " IWORK(3) NOT COMPATIBLE WITH IOUT", (ftnlen)34)
		;
	e_wsle();
	arret = TRUE_;
    }
/* -------- MSTAB     PARAMETER FOR STABILITY CHECK */
    if (iwork[4] == 0) {
	mstab = 1;
    } else {
	mstab = iwork[4];
    }
/* -------- JSTAB     PARAMETER FOR STABILITY CHECK */
    if (iwork[5] == 0) {
	jstab = 2;
    } else {
	jstab = iwork[5];
    }
/* -------- IDERR  PARAMETER FOR ERROR ESTIMATION IN DENSE OUTPUT */
    if (iwork[6] == 0) {
	if (*iout <= 1) {
	    iderr = 1;
	}
	if (*iout >= 2) {
	    iderr = 0;
	}
    } else {
	iderr = iwork[6];
	if (*iout <= 1) {
	    s_wsle(&io___16);
	    do_lio(&c__9, &c__1, " ERROR ESTIMATION IN DENSE OUTPUT", (ftnlen)
		    33);
	    do_lio(&c__9, &c__1, " NOT POSSIBLE, WRONG IWORK(6)=", (ftnlen)30)
		    ;
	    do_lio(&c__3, &c__1, (char *)&iwork[6], (ftnlen)sizeof(integer));
	    e_wsle();
	    arret = TRUE_;
	}
    }
/* -------- MUDIF */
    if (iwork[7] == 0) {
	mudif = 4;
    } else {
	mudif = iwork[7];
	if (mudif <= 0 || mudif >= 7) {
	    s_wsle(&io___18);
	    do_lio(&c__9, &c__1, " WRONG INPUT IWORK(7)=", (ftnlen)22);
	    do_lio(&c__3, &c__1, (char *)&iwork[7], (ftnlen)sizeof(integer));
	    e_wsle();
	    arret = TRUE_;
	}
    }
/* -------- NRDENS   NUMBER OF DENSE OUTPUT COMPONENTS */
    nrdens = iwork[8];
    if (nrdens < 0 || nrdens > *n) {
	s_wsle(&io___20);
	do_lio(&c__9, &c__1, " CURIOUS INPUT IWORK(8)=", (ftnlen)24);
	do_lio(&c__3, &c__1, (char *)&iwork[8], (ftnlen)sizeof(integer));
	e_wsle();
	arret = TRUE_;
    }
    if (nrdens == *n) {
	i__1 = nrdens;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L17: */
	    iwork[i__ + 20] = i__;
	}
    }
/* -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0 */
    if (work[1] == 0.) {
	uround = 2.3e-16;
    } else {
	uround = work[1];
	if (uround <= 1e-35 || uround >= 1.) {
	    s_wsle(&io___23);
	    do_lio(&c__9, &c__1, " WHICH MACHINE DO YOU HAVE? YOUR UROUND WA"
		    "S:", (ftnlen)44);
	    do_lio(&c__5, &c__1, (char *)&work[1], (ftnlen)sizeof(doublereal))
		    ;
	    e_wsle();
	    arret = TRUE_;
	}
    }
/* -------- MAXIMAL STEP SIZE */
    if (work[2] == 0.) {
	hmax = *xend - *x;
    } else {
	hmax = abs(work[2]);
    }
/* -------- STEP SIZE REDUCTION FACTOR */
    if (work[3] == 0.) {
	safe3 = .5;
    } else {
	safe3 = work[3];
	if (safe3 <= uround || safe3 >= 1.) {
	    s_wsle(&io___26);
	    do_lio(&c__9, &c__1, " CURIOUS INPUT WORK(3)=", (ftnlen)23);
	    do_lio(&c__5, &c__1, (char *)&work[3], (ftnlen)sizeof(doublereal))
		    ;
	    e_wsle();
	    arret = TRUE_;
	}
    }
/* -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION */
    if (work[4] == 0.) {
	fac1 = .02;
    } else {
	fac1 = work[4];
    }
    if (work[5] == 0.) {
	fac2 = 4.;
    } else {
	fac2 = work[5];
    }
/* -------  FAC3, FAC4   PARAMETERS FOR THE ORDER SELECTION */
    if (work[6] == 0.) {
	fac3 = .8;
    } else {
	fac3 = work[6];
    }
    if (work[7] == 0.) {
	fac4 = .9;
    } else {
	fac4 = work[7];
    }
/* ------- SAFE1, SAFE2 SAFETY FACTORS FOR STEP SIZE PREDICTION */
    if (work[8] == 0.) {
	safe1 = .65;
    } else {
	safe1 = work[8];
    }
    if (work[9] == 0.) {
	safe2 = .94;
    } else {
	safe2 = work[9];
    }
/* ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK ----- */
    lfsafe = (km << 1) * km + km;
    iedy = 21;
    ieyh1 = iedy + *n;
    ieyh2 = ieyh1 + *n;
    iedz = ieyh2 + *n;
    iescal = iedz + *n;
    iet = iescal + *n;
    iefs = iet + km * *n;
    ieys = iefs + lfsafe * nrdens;
    iehh = ieys + km * nrdens;
    iew = iehh + km;
    iea = iew + km;
    iefac = iea + km;
/* ------ TOTAL STORAGE REQUIREMENT ----------- */
    ieco = iefac + (km << 1);
    istore = ieco + ((km << 1) + 5) * nrdens - 1;
    if (istore > *lwork) {
	s_wsle(&io___48);
	do_lio(&c__9, &c__1, " INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=", (
		ftnlen)43);
	do_lio(&c__3, &c__1, (char *)&istore, (ftnlen)sizeof(integer));
	e_wsle();
	arret = TRUE_;
    }
/* ------- ENTRY POINTS FOR INTEGER WORKSPACE ----- */
    icom = 21;
    ienj = icom + nrdens;
/* --------- TOTAL REQUIREMENT --------------- */
    ieip = ienj + km;
    istore = ieip + km;
    if (istore > *liwork) {
	s_wsle(&io___52);
	do_lio(&c__9, &c__1, " INSUFF. STORAGE FOR IWORK, MIN. LIWORK=", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&istore, (ftnlen)sizeof(integer));
	e_wsle();
	arret = TRUE_;
    }
/* ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1 */
    if (arret) {
	*idid = -1;
	return 0;
    }
/* -------- CALL TO CORE INTEGRATOR ------------ */
    nrd = max(1,nrdens);
/* Computing MAX */
    i__1 = 1, i__2 = ((km << 1) + 5) * nrdens;
    ncom = max(i__1,i__2);
    ok_odex_odxcor_(n, (U_fp)fcn, x, &y[1], xend, &hmax, h__, &rtol[1], &atol[1], 
	    itol, &km,  iout, idid, &nmax, &uround, &work[iedy], 
	    &work[ieyh1], &work[ieyh2], &work[iedz], &work[iescal], &work[
	    iefs], &work[ieys], &work[iet], &work[iehh], &work[iew], &work[
	    iea], &work[ieco], &ncom, &iwork[icom], &iwork[ienj], &iwork[ieip]
	    , &nsequ, &mstab, &jstab, &lfsafe, &safe1, &safe2, &safe3, &fac1, 
	    &fac2, &fac3, &fac4, &iderr, &work[iefac], &mudif, &nrd, &nfcn, 
            &nstep, &naccpt, &nrejct, params);
    iwork[17] = nfcn;
    iwork[18] = nstep;
    iwork[19] = naccpt;
    iwork[20] = nrejct;
/* ----------- RETURN ----------- */
    return 0;
} /* odex_ */




/*  ----- ... AND HERE IS THE CORE INTEGRATOR  ---------- */

/* Subroutine */ int ok_odex_odxcor_(integer *n, S_fp fcn, doublereal *x, doublereal *
	y, doublereal *xend, doublereal *hmax, doublereal *h__, doublereal *
	rtol, doublereal *atol, integer *itol, integer *km, 
	integer *iout, integer *idid, integer *nmax, doublereal *uround, 
	doublereal *dy, doublereal *yh1, doublereal *yh2, doublereal *dz, 
	doublereal *scal, doublereal *fsafe, doublereal *ysafe, doublereal *t,
	 doublereal *hh, doublereal *w, doublereal *a, doublereal *dens, 
	integer *ncom, integer *icomp, integer *nj, integer *ipoint, integer *
	nsequ, integer *mstab, integer *jstab, integer *lfsafe, doublereal *
	safe1, doublereal *safe2, doublereal *safe3, doublereal *fac1, 
	doublereal *fac2, doublereal *fac3, doublereal *fac4, integer *iderr, 
	doublereal *errfac, integer *mudif, integer *nrd, integer *nfcn, integer *nstep, integer *naccpt, 
	integer *nrejct, void* params)
{
    /* Format strings */
    static char fmt_979[] = "(\002 EXIT OF ODEX AT X=\002,d14.7,\002   H="
	    "\002,d14.7)";

    /* System generated locals */
    integer t_dim1, t_offset, fsafe_dim1, fsafe_offset, ysafe_dim1, 
	    ysafe_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), d_lg10(doublereal *), sqrt(
	    doublereal), pow_di(doublereal *, integer *), pow_dd(doublereal *,
	     doublereal *);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    integer i__, j, k, l, kc, kk, mu;
    doublereal fac;
    real hhh;
    integer kmi, kln;
    doublereal err;
    integer krn, ipt, kbeg, lbeg, lend;
    logical last;
    integer kmit;
    doublereal prod;
    logical atov;
    doublereal xold;
    integer kopt;
    doublereal errx;
    integer njadd;
    doublereal facnj;
    
    real xoldd;
    integer irtrn;
    doublereal dblenj;
    logical reject;
    doublereal factor, hoptde, errold, posneg;
    
    doublereal errint;

    /* Fortran I/O blocks */
    static cilist io___91 = { 0, 6, 0, fmt_979, 0 };


/* ---------------------------------------------------------- */
/*     CORE INTEGRATOR FOR ODEX */
/*     PARAMETERS SAME AS IN ODEX WITH WORKSPACE ADDED */
/* ---------------------------------------------------------- */
/*         DECLARATIONS */
/* ---------------------------------------------------------- */
/* --- DEFINE THE STEP SIZE SEQUENCE */
    /* Parameter adjustments */
    --scal;
    --dz;
    --yh2;
    --yh1;
    --dy;
    --y;
    --rtol;
    --atol;
    --errfac;
    --ipoint;
    --nj;
    --a;
    --w;
    --hh;
    t_dim1 = *km;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --dens;
    --icomp;
    ysafe_dim1 = *km;
    ysafe_offset = 1 + ysafe_dim1;
    ysafe -= ysafe_offset;
    fsafe_dim1 = *lfsafe;
    fsafe_offset = 1 + fsafe_dim1;
    fsafe -= fsafe_offset;

    /* Function Body */
    if (*nsequ == 1) {
	i__1 = *km;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L1: */
	    nj[i__] = i__ << 1;
	}
    }
    if (*nsequ == 2) {
	nj[1] = 2;
	i__1 = *km;
	for (i__ = 2; i__ <= i__1; ++i__) {
/* L2: */
	    nj[i__] = (i__ << 2) - 4;
	}
    }
    if (*nsequ == 3) {
	nj[1] = 2;
	nj[2] = 4;
	nj[3] = 6;
	i__1 = *km;
	for (i__ = 4; i__ <= i__1; ++i__) {
/* L11: */
	    nj[i__] = nj[i__ - 2] << 1;
	}
    }
    if (*nsequ == 4) {
	i__1 = *km;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L3: */
	    nj[i__] = (i__ << 2) - 2;
	}
    }
    if (*nsequ == 5) {
	i__1 = *km;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L6: */
	    nj[i__] = i__ << 2;
	}
    }
/* --- DEFINE THE A(I) FOR ORDER SELECTION */
    a[1] = nj[1] + 1.;
    i__1 = *km;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L4: */
	a[i__] = a[i__ - 1] + nj[i__];
    }
/* --- INITIAL SCALING */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*itol == 0) {
	    scal[i__] = atol[1] + rtol[1] * (d__1 = y[i__], abs(d__1));
	} else {
	    scal[i__] = atol[i__] + rtol[i__] * (d__1 = y[i__], abs(d__1));
	}
/* L8: */
    }
/* --- INITIAL PREPARATIONS */
    d__1 = *xend - *x;
    posneg = d_sign(&c_b67, &d__1);
/* Computing MAX */
/* Computing MIN */
    d__1 = rtol[1] + 1e-40;
    i__3 = *km - 1, i__4 = (integer) (-d_lg10(&d__1) * .6 + 1.5);
    i__1 = 2, i__2 = min(i__3,i__4);
    k = max(i__1,i__2);
    *hmax = abs(*hmax);
/* Computing MAX */
    d__1 = abs(*h__);
    *h__ = max(d__1,1e-4);
/* Computing MIN */
    d__2 = min(*h__,*hmax), d__3 = (d__1 = *xend - *x, abs(d__1)) / 2.;
    *h__ = posneg * min(d__2,d__3);
    if (*iout >= 1) {
	if (*iout >= 2) {
	    ipoint[1] = 0;
	    i__1 = *km;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		njadd = (i__ << 2) - 2;
		if (nj[i__] > njadd) {
		    ++njadd;
		}
/* L5: */
		ipoint[i__ + 1] = ipoint[i__] + njadd;
	    }
	    i__1 = *km << 1;
	    for (mu = 1; mu <= i__1; ++mu) {
		errx = sqrt(mu / (mu + 4.)) * .5;
/* Computing 2nd power */
		d__1 = mu + 4.;
		prod = 1. / (d__1 * d__1);
		i__2 = mu;
		for (j = 1; j <= i__2; ++j) {
/* L7: */
		    prod = prod * errx / j;
		}
/* L9: */
		errfac[mu] = prod;
	    }
	    ipt = 0;
	}
	irtrn = 0;
	xold = *x;
	i__1 = *naccpt + 1;
	if (irtrn < 0) {
	    goto L120;
	}
    }
    err = 0.;
    errold = 1e10;
    hoptde = posneg * *hmax;
    w[1] = 0.;
    reject = FALSE_;
    last = FALSE_;
L10:
    atov = FALSE_;
/* --- IS XEND REACHED IN THE NEXT STEP? */
    if ((d__1 = *xend - *x, abs(d__1)) * .1 <= abs(*x) * *uround) {
	goto L110;
    }
/* Computing MIN */
    d__2 = abs(*h__), d__3 = (d__1 = *xend - *x, abs(d__1)), d__2 = min(d__2,
	    d__3), d__2 = min(d__2,*hmax), d__3 = abs(hoptde);
    *h__ = posneg * min(d__2,d__3);
    if ((*x + *h__ * 1.01 - *xend) * posneg > 0.) {
	*h__ = *xend - *x;
	last = TRUE_;
    }
    if (*nstep == 0 || *iout != 2) {
	(*fcn)(*x, &y[1], &dz[1], params);
    }
    ++(*nfcn);
/* --- THE FIRST AND LAST STEP */
    if (*nstep == 0 || last) {
	ipt = 0;
	++(*nstep);
	i__1 = k;
	for (j = 1; j <= i__1; ++j) {
	    kc = j;
	    ok_odex_midex_(&j, x, &y[1], h__, hmax, n, (S_fp)fcn, &dy[1], &yh1[1], &
		    yh2[1], &dz[1], &t[t_offset], &nj[1], &hh[1], &w[1], &err,
		     &fac, &a[1], safe1, uround, fac1, fac2, safe2, &scal[1], 
		    &atov, safe3, &reject, km, &rtol[1], &atol[1], itol, 
		    mstab, jstab, &errold, &fsafe[fsafe_offset], lfsafe, iout,
		     &ipt, &ysafe[ysafe_offset], &icomp[1], nrd, nfcn, params);
	    if (atov) {
		goto L10;
	    }
/* L20: */
	    if (j > 1 && err <= 1.) {
		goto L60;
	    }
	}
	goto L55;
    }
/* --- BASIC INTEGRATION STEP */
L30:
    ipt = 0;
    ++(*nstep);
    if (*nstep >= *nmax) {
	goto L120;
    }
    kc = k - 1;
    i__1 = kc;
    for (j = 1; j <= i__1; ++j) {
	ok_odex_midex_(&j, x, &y[1], h__, hmax, n, (S_fp)fcn, &dy[1], &yh1[1], &yh2[1]
		, &dz[1], &t[t_offset], &nj[1], &hh[1], &w[1], &err, &fac, &a[
		1], safe1, uround, fac1, fac2, safe2, &scal[1], &atov, safe3, 
		&reject, km, &rtol[1], &atol[1], itol, mstab, jstab, &errold, 
		&fsafe[fsafe_offset], lfsafe, iout, &ipt, &ysafe[ysafe_offset]
		, &icomp[1], nrd, nfcn, params);
	if (atov) {
	    goto L10;
	}
/* L40: */
    }
/* --- CONVERGENCE MONITOR */
    if (k == 2 || reject) {
	goto L50;
    }
    if (err <= 1.) {
	goto L60;
    }
/* Computing 2nd power */
    d__1 = nj[k + 1] * nj[k] / 4.;
    if (err > d__1 * d__1) {
	goto L100;
    }
L50:
    ok_odex_midex_(&k, x, &y[1], h__, hmax, n, (S_fp)fcn, &dy[1], &yh1[1], &yh2[1], &
	    dz[1], &t[t_offset], &nj[1], &hh[1], &w[1], &err, &fac, &a[1], 
	    safe1, uround, fac1, fac2, safe2, &scal[1], &atov, safe3, &reject,
	     km, &rtol[1], &atol[1], itol, mstab, jstab, &errold, &fsafe[
	    fsafe_offset], lfsafe, iout, &ipt, &ysafe[ysafe_offset], &icomp[1]
	    , nrd, nfcn, params);
    if (atov) {
	goto L10;
    }
    kc = k;
    if (err <= 1.) {
	goto L60;
    }
/* --- HOPE FOR CONVERGENCE IN LINE K+1 */
L55:
/* Computing 2nd power */
    d__1 = nj[k + 1] / 2.;
    if (err > d__1 * d__1) {
	goto L100;
    }
    kc = k + 1;
    ok_odex_midex_(&kc, x, &y[1], h__, hmax, n, (S_fp)fcn, &dy[1], &yh1[1], &yh2[1], &
	    dz[1], &t[t_offset], &nj[1], &hh[1], &w[1], &err, &fac, &a[1], 
	    safe1, uround, fac1, fac2, safe2, &scal[1], &atov, safe3, &reject,
	     km, &rtol[1], &atol[1], itol, mstab, jstab, &errold, &fsafe[
	    fsafe_offset], lfsafe, iout, &ipt, &ysafe[ysafe_offset], &icomp[1]
	    , nrd, nfcn, params);
    if (atov) {
	goto L10;
    }
    if (err > 1.) {
	goto L100;
    }
/* --- STEP IS ACCEPTED */
L60:
    xold = *x;
    *x += *h__;
    if (*iout >= 2) {
/* ---  KMIT = MU OF THE PAPER */
	kmit = (kc << 1) - *mudif + 1;
	i__1 = *nrd;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L69: */
	    dens[i__] = y[icomp[i__]];
	}
	xoldd = xold;
	hhh = *h__;
	i__1 = *nrd;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L76: */
	    dens[*nrd + i__] = *h__ * dz[icomp[i__]];
	}
	kln = *nrd << 1;
	i__1 = *nrd;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L176: */
	    dens[kln + i__] = t[icomp[i__] * t_dim1 + 1];
	}
/* --- COMPUTE SOLUTION AT MID-POINT ---- */
	i__1 = kc;
	for (j = 2; j <= i__1; ++j) {
	    dblenj = (doublereal) nj[j];
	    for (l = j; l >= 2; --l) {
/* Computing 2nd power */
		d__1 = dblenj / nj[l - 1];
		factor = d__1 * d__1 - 1.;
		i__2 = *nrd;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ysafe[l - 1 + i__ * ysafe_dim1] = ysafe[l + i__ * 
			    ysafe_dim1] + (ysafe[l + i__ * ysafe_dim1] - 
			    ysafe[l - 1 + i__ * ysafe_dim1]) / factor;
/* L473: */
		}
	    }
	}
	krn = *nrd << 2;
	i__2 = *nrd;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L474: */
	    dens[krn + i__] = ysafe[i__ * ysafe_dim1 + 1];
	}
/* --- COMPUTE FIRST DERIVATIVE AT RIGHT END ---- */
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L478: */
	    yh1[i__] = t[i__ * t_dim1 + 1];
	}
	(*fcn)(*x, &yh1[1], &yh2[1], params);
	krn = *nrd * 3;
	i__2 = *nrd;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L274: */
	    dens[krn + i__] = yh2[icomp[i__]] * *h__;
	}
/* --- THE LOOP --- */
	i__2 = kmit;
	for (kmi = 1; kmi <= i__2; ++kmi) {
/* --- COMPUTE KMI-TH DERIVATIVE AT MID-POINT ---- */
	    kbeg = (kmi + 1) / 2;
	    i__1 = kc;
	    for (kk = kbeg; kk <= i__1; ++kk) {
		d__1 = nj[kk] / 2.;
		i__3 = kmi - 1;
		facnj = pow_di(&d__1, &i__3);
		ipt = ipoint[kk + 1] - (kk << 1) + kmi;
		i__3 = *nrd;
		for (i__ = 1; i__ <= i__3; ++i__) {
/* L371: */
		    ysafe[kk + i__ * ysafe_dim1] = fsafe[ipt + i__ * 
			    fsafe_dim1] * facnj;
		}
/* L375: */
	    }
	    i__1 = kc;
	    for (j = kbeg + 1; j <= i__1; ++j) {
		dblenj = (doublereal) nj[j];
		i__3 = kbeg + 1;
		for (l = j; l >= i__3; --l) {
/* Computing 2nd power */
		    d__1 = dblenj / nj[l - 1];
		    factor = d__1 * d__1 - 1.;
		    i__4 = *nrd;
		    for (i__ = 1; i__ <= i__4; ++i__) {
			ysafe[l - 1 + i__ * ysafe_dim1] = ysafe[l + i__ * 
				ysafe_dim1] + (ysafe[l + i__ * ysafe_dim1] - 
				ysafe[l - 1 + i__ * ysafe_dim1]) / factor;
/* L373: */
		    }
		}
	    }
	    krn = (kmi + 4) * *nrd;
	    i__4 = *nrd;
	    for (i__ = 1; i__ <= i__4; ++i__) {
/* L374: */
		dens[krn + i__] = ysafe[kbeg + i__ * ysafe_dim1] * *h__;
	    }
	    if (kmi == kmit) {
		goto L180;
	    }
/* --- COMPUTE DIFFERENCES */
	    i__4 = kc;
	    for (kk = (kmi + 2) / 2; kk <= i__4; ++kk) {
		lbeg = ipoint[kk + 1];
		lend = ipoint[kk] + kmi + 1;
		if (kmi == 1 && *nsequ == 4) {
		    lend += 2;
		}
		i__3 = lend;
		for (l = lbeg; l >= i__3; l += -2) {
		    i__1 = *nrd;
		    for (i__ = 1; i__ <= i__1; ++i__) {
/* L64: */
			fsafe[l + i__ * fsafe_dim1] -= fsafe[l - 2 + i__ * 
				fsafe_dim1];
		    }
		}
		if (kmi == 1 && *nsequ == 4) {
		    l = lend - 2;
		    i__1 = *nrd;
		    for (i__ = 1; i__ <= i__1; ++i__) {
/* L65: */
			fsafe[l + i__ * fsafe_dim1] -= dz[icomp[i__]];
		    }
		}
/* L66: */
	    }
/* --- COMPUTE DIFFERENCES */
	    i__4 = kc;
	    for (kk = (kmi + 2) / 2; kk <= i__4; ++kk) {
		lbeg = ipoint[kk + 1] - 1;
		lend = ipoint[kk] + kmi + 2;
		i__1 = lend;
		for (l = lbeg; l >= i__1; l += -2) {
		    i__3 = *nrd;
		    for (i__ = 1; i__ <= i__3; ++i__) {
/* L164: */
			fsafe[l + i__ * fsafe_dim1] -= fsafe[l - 2 + i__ * 
				fsafe_dim1];
		    }
		}
/* L166: */
	    }
L180:
	    ;
	}
	ok_odex_interp_(*nrd, &dens[1], kmit);
/* --- ESTIMATION OF INTERPOLATION ERROR */
	if (*iderr == 0 && kmit >= 1) {
	    errint = 0.;
	    i__2 = *nrd;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L187: */
/* Computing 2nd power */
		d__1 = dens[(kmit + 4) * *nrd + i__] / scal[icomp[i__]];
		errint += d__1 * d__1;
	    }
	    errint = sqrt(errint / *nrd) * errfac[kmit];
/* Computing MAX */
	    d__2 = 1. / (kmit + 4);
	    d__1 = pow_dd(&errint, &d__2);
	    hoptde = *h__ / max(d__1,.01);
	    if (errint > 10.) {
		*h__ = hoptde;
		*x = xold;
		++(*nrejct);
		reject = TRUE_;
		goto L10;
	    }
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L189: */
	    dz[i__] = yh2[i__];
	}
    }
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L70: */
	y[i__] = t[i__ * t_dim1 + 1];
    }
    ++(*naccpt);
    if (*iout >= 1) {
	i__2 = *naccpt + 1;
	if (irtrn < 0) {
	    goto L120;
	}
    }
/* --- COMPUTE OPTIMAL ORDER */
    if (kc == 2) {
/* Computing MIN */
	i__2 = 3, i__4 = *km - 1;
	kopt = min(i__2,i__4);
	if (reject) {
	    kopt = 2;
	}
	goto L80;
    }
    if (kc <= k) {
	kopt = kc;
	if (w[kc - 1] < w[kc] * *fac3) {
	    kopt = kc - 1;
	}
	if (w[kc] < w[kc - 1] * *fac4) {
/* Computing MIN */
	    i__2 = kc + 1, i__4 = *km - 1;
	    kopt = min(i__2,i__4);
	}
    } else {
	kopt = kc - 1;
	if (kc > 3 && w[kc - 2] < w[kc - 1] * *fac3) {
	    kopt = kc - 2;
	}
	if (w[kc] < w[kopt] * *fac4) {
/* Computing MIN */
	    i__2 = kc, i__4 = *km - 1;
	    kopt = min(i__2,i__4);
	}
    }
/* --- AFTER A REJECTED STEP */
L80:
    if (reject) {
	k = min(kopt,kc);
/* Computing MIN */
	d__2 = abs(*h__), d__3 = (d__1 = hh[k], abs(d__1));
	*h__ = posneg * min(d__2,d__3);
	reject = FALSE_;
	goto L10;
    }
/* --- COMPUTE STEPSIZE FOR NEXT STEP */
    if (kopt <= kc) {
	*h__ = hh[kopt];
    } else {
	if (kc < k && w[kc] < w[kc - 1] * *fac4) {
	    *h__ = hh[kc] * a[kopt + 1] / a[kc];
	} else {
	    *h__ = hh[kc] * a[kopt] / a[kc];
	}
    }
    k = kopt;
    *h__ = posneg * abs(*h__);
    goto L10;
/* --- STEP IS REJECTED */
L100:
/* Computing MIN */
    i__2 = min(k,kc), i__4 = *km - 1;
    k = min(i__2,i__4);
    if (k > 2 && w[k - 1] < w[k] * *fac3) {
	--k;
    }
    ++(*nrejct);
    *h__ = posneg * hh[k];
    reject = TRUE_;
    goto L30;
/* --- SOLUTION EXIT */
L110:
    *idid = 1;
    return 0;
/* --- FAIL EXIT */
L120:
    s_wsfe(&io___91);
    do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*h__), (ftnlen)sizeof(doublereal));
    e_wsfe();
    *idid = -1;
    return 0;
} /* odxcor_ */


/* Subroutine */ int ok_odex_midex_(integer *jj, doublereal *x, doublereal *y, 
	doublereal *h__, doublereal *hmax, integer *nn, S_fp fcn, doublereal *
	dy, doublereal *yh1, doublereal *yh2, doublereal *dz, doublereal *t, 
	integer *nj, doublereal *hh, doublereal *w, doublereal *err, 
	doublereal *fac, doublereal *a, doublereal *safe1, doublereal *uround,
	 doublereal *fac1, doublereal *fac2, doublereal *safe2, doublereal *
	scal, logical *atov, doublereal *safe3, logical *reject, integer *km, 
	doublereal *rtol, doublereal *atol, integer *itol, integer *mstab, 
	integer *jstab, doublereal *errold, doublereal *fsafe, integer *
	lfsafe, integer *iout, integer *ipt, doublereal *ysafe, integer *
	icomp, integer *nrd, integer *nfcn, void* params)
{
    /* System generated locals */
    integer t_dim1, t_offset, fsafe_dim1, fsafe_offset, ysafe_dim1, 
	    ysafe_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, l, m;
    doublereal hj;
    integer mm;
    doublereal ys, t1i, del1, del2, expo, quot;
    integer njmid;
    doublereal facmin, dblenj;

    integer j = *jj;
    integer n = *nn;
    
/* --- THIS SUBROUTINE COMPUTES THE J-TH LINE OF THE */
/* --- EXTRAPOLATION TABLE AND PROVIDES AN ESTIMATION */
/* --- OF THE OPTIMAL STEPSIZE */
    /* Parameter adjustments */
    --scal;
    --dz;
    --yh2;
    --yh1;
    --dy;
    --y;
    --a;
    --w;
    --hh;
    --nj;
    t_dim1 = *km;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --rtol;
    --atol;
    --icomp;
    ysafe_dim1 = *km;
    ysafe_offset = 1 + ysafe_dim1;
    ysafe -= ysafe_offset;
    fsafe_dim1 = *lfsafe;
    fsafe_offset = 1 + fsafe_dim1;
    fsafe -= fsafe_offset;
    

    /* Function Body */
    hj = *h__ / nj[j];
/* --- EULER STARTING STEP */
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yh1[i__] = y[i__];
/* L30: */
	yh2[i__] = y[i__] + hj * dz[i__];
    }
/* --- EXPLICIT MIDPOINT RULE */
    m = nj[j] - 1;
    njmid = nj[j] / 2;
    i__1 = m;
    for (mm = 1; mm <= i__1; ++mm) {
	if (*iout >= 2 && mm == njmid) {
	    i__2 = *nrd;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L31: */
		ysafe[j + i__ * ysafe_dim1] = yh2[icomp[i__]];
	    }
	}
	d__1 = *x + hj * mm;
	(*fcn)(d__1, &yh2[1], &dy[1], params);
	if (*iout >= 2 && (i__2 = mm - njmid, abs(i__2)) <= (j << 1) - 1) {
	    ++(*ipt);
	    i__2 = *nrd;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L32: */
		fsafe[*ipt + i__ * fsafe_dim1] = dy[icomp[i__]];
	    }
	}
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ys = yh1[i__];
	    yh1[i__] = yh2[i__];
/* L34: */
	    yh2[i__] = ys + hj * 2. * dy[i__];
	}
	if (mm <= *mstab && j <= *jstab) {
/* --- STABILITY CHECK */
	    del1 = 0.;
	    i__2 = n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L21: */
/* Computing 2nd power */
		d__1 = dz[i__] / scal[i__];
		del1 += d__1 * d__1;
	    }
	    del2 = 0.;
	    i__2 = n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L26: */
/* Computing 2nd power */
		d__1 = (dy[i__] - dz[i__]) / scal[i__];
		del2 += d__1 * d__1;
	    }
	    quot = del2 / max(*uround,del1);
	    if (quot > 4.) {
		++(*nfcn);
		goto L79;
	    }
	}
/* L35: */
    }
/* --- FINAL SMOOTHING STEP */
    d__1 = *x + *h__;
    (*fcn)(d__1, &yh2[1], &dy[1], params);
    if (*iout >= 2 && njmid <= (j << 1) - 1) {
	++(*ipt);
	i__1 = *nrd;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L39: */
	    fsafe[*ipt + i__ * fsafe_dim1] = dy[icomp[i__]];
	}
    }
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
	t[j + i__ * t_dim1] = (yh1[i__] + yh2[i__] + hj * dy[i__]) / 2.;
    }
    *nfcn += nj[j];
/* --- POLYNOMIAL EXTRAPOLATION */
    if (j == 1) {
	return 0;
    }
    dblenj = (doublereal) nj[j];
    for (l = j; l >= 2; --l) {
/* Computing 2nd power */
	d__1 = dblenj / nj[l - 1];
	*fac = d__1 * d__1 - 1.;
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    t[l - 1 + i__ * t_dim1] = t[l + i__ * t_dim1] + (t[l + i__ * 
		    t_dim1] - t[l - 1 + i__ * t_dim1]) / *fac;
/* L60: */
	}
    }
    *err = 0.;
/* --- SCALING */
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__3 = (d__1 = y[i__], abs(d__1)), d__4 = (d__2 = t[i__ * t_dim1 + 1],
		 abs(d__2));
	t1i = max(d__3,d__4);
	if (*itol == 0) {
	    scal[i__] = atol[1] + rtol[1] * t1i;
	} else {
	    scal[i__] = atol[i__] + rtol[i__] * t1i;
	}
/* L65: */
/* Computing 2nd power */
	d__1 = (t[i__ * t_dim1 + 1] - t[i__ * t_dim1 + 2]) / scal[i__];
	*err += d__1 * d__1;
    }
    *err = sqrt(*err / n);
    if (*err * *uround >= 1.) {
	goto L79;
    }
    if (j > 2 && *err >= *errold) {
	goto L79;
    }
/* Computing MAX */
    d__1 = *err * 4;
    *errold = max(d__1,1.);
/* --- COMPUTE OPTIMAL STEPSIZES */
    expo = 1. / ((j << 1) - 1);
    facmin = pow(*fac1, expo);
/* Computing MIN */
/* Computing MAX */
    d__5 = *err / *safe1;
    d__3 = facmin, d__4 = pow(d__5, expo) / *safe2;
    d__1 = *fac2 / facmin, d__2 = max(d__3,d__4);
    *fac = min(d__1,d__2);
    *fac = 1. / *fac;
/* Computing MIN */
    d__1 = abs(*h__) * *fac;
    hh[j] = min(d__1,*hmax);
    w[j] = a[j] / hh[j];
    return 0;
L79:
    *atov = TRUE_;
    *h__ *= *safe3;
    *reject = TRUE_;
    return 0;
} /* midex_ */


/* Subroutine */ int ok_odex_interp_(integer n, doublereal *y, integer imit)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    doublereal a[31];
    integer i__;
    doublereal y0, y1;
    integer im;
    doublereal ph0, ph1, ph2, ph3, yp0, yp1, fac1, fac2, aspl, bspl, ydiff;

/* --- COMPUTES THE COEFFICIENTS OF THE INTERPOLATION FORMULA */
/* --- BEGIN WITH HERMITE INTERPOLATION */
    /* Parameter adjustments */
    --y;

    /* Function Body */
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y0 = y[i__];
	y1 = y[(n << 1) + i__];
	yp0 = y[n + i__];
	yp1 = y[n * 3 + i__];
	ydiff = y1 - y0;
	aspl = -yp1 + ydiff;
	bspl = yp0 - ydiff;
	y[n + i__] = ydiff;
	y[(n << 1) + i__] = aspl;
	y[n * 3 + i__] = bspl;
	if (imit < 0) {
	    goto L100;
	}
/* --- COMPUTE THE DERIVATIVES OF HERMITE AT MIDPOINT */
	ph0 = (y0 + y1) * .5 + (aspl + bspl) * .125;
	ph1 = ydiff + (aspl - bspl) * .25;
	ph2 = -(yp0 - yp1);
	ph3 = (bspl - aspl) * 6.;
/* --- COMPUTE THE FURTHER COEFFICIENTS */
	if (imit < 1) {
	    goto L20;
	}
	a[1] = (y[n * 5 + i__] - ph1) * 16.;
	if (imit < 3) {
	    goto L20;
	}
	a[3] = (y[n * 7 + i__] - ph3 + a[1] * 3) * 16.;
	if (imit < 5) {
	    goto L20;
	}
	i__2 = imit;
	for (im = 5; im <= i__2; im += 2) {
	    fac1 = im * (im - 1) / 2.;
	    fac2 = fac1 * (im - 2) * (im - 3) * 2.;
/* L10: */
	    a[im] = (y[(im + 4) * n + i__] + fac1 * a[im - 2] - fac2 * a[im 
		    - 4]) * 16.;
	}
L20:
	a[0] = (y[(n << 2) + i__] - ph0) * 16.;
	if (imit < 2) {
	    goto L60;
	}
	a[2] = (y[n * 6 + i__] - ph2 + a[0]) * 16.;
	if (imit < 4) {
	    goto L60;
	}
	i__2 = imit;
	for (im = 4; im <= i__2; im += 2) {
	    fac1 = im * (im - 1) / 2.;
	    fac2 = (doublereal) (im * (im - 1) * (im - 2) * (im - 3));
/* L30: */
	    a[im] = (y[n * (im + 4) + i__] + a[im - 2] * fac1 - a[im - 4] * 
		    fac2) * 16.;
	}
L60:
	i__2 = imit;
	for (im = 0; im <= i__2; ++im) {
/* L70: */
	    y[n * (im + 4) + i__] = a[im];
	}
L100:
	;
    }
    return 0;
} /* interp_ */


ok_system** ok_integrate_odex(ok_system* initial, const gsl_vector* times, ok_integrator_options* options,
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
    
    integer DIMENSIONS = NDIMS * 7;
    
    int km = 9;
    integer liwork = 2 * km + 21; 
    integer* iwork = (integer*) calloc(liwork, sizeof(integer));
    iwork[1] = km;
    //iwork[6] = 6;
    //iwork[2] = 5;
    
    integer lwork = DIMENSIONS*(km+5) + 5*km + 30;
    double BUFSIZE = lwork + DIMENSIONS + 1;
    if (options->buffer == NULL || options->buffer->size < BUFSIZE) {
        if (options->buffer != NULL)
            gsl_vector_free(options->buffer);
        options->buffer = gsl_vector_alloc(BUFSIZE);
    }
    
    memset(options->buffer->data, 0, sizeof(double)*BUFSIZE);
    double* work = options->buffer->data;
    
    double prevTime = startTime;
    gsl_matrix* prevOrbits = initial->orbits;
    
    double h = 0.01;
    double* xyz = work + lwork;
    
    MATRIX_MEMCPY_TOARRAY(xyz, initial->xyz);
    
    doublereal relerr = options->rel_acc;
    doublereal abserr = options->abs_acc;
    integer itol = 0;
    integer iout = 0;
    
    void* params = (void*) initial;
    ok_progress progress = options->progress;
    
    // Loop through the times vector
    for (int i = 0; i < SAMPLES; i++) {
        double time = times->data[i];
        
        double fromTime = prevTime;
        double toTime = time;
        
        // Integrate between prevTime and time
        if (fabs(time - prevTime) > 1e-12) {
            bool invert = (time - prevTime < 0);
            
            if (invert)
                for (int i = 0; i < NDIMS; i++) {
                    xyz[i*7 + 4] *= -1;
                    xyz[i*7 + 5] *= -1;
                    xyz[i*7 + 6] *= -1;
                    
                    fromTime = time;
                    toTime = prevTime;
                }
                
                
                integer idid = 1;
                 
                ok_odex_(&DIMENSIONS, &ok_force,
                        &fromTime, xyz, 
                        &toTime, &h,
                        &relerr, &abserr,
                        &itol, &iout,
                        work, &lwork,
                        iwork, &liwork, &idid,
                        params);
                
                if (invert)
                for (int i = 0; i < NDIMS; i++) {
                    xyz[i*7 + 4] *= -1;
                    xyz[i*7 + 5] *= -1;
                    xyz[i*7 + 6] *= -1;
                }
                
                if (idid != 1 || (((ok_system*) params)->flag & INTEGRATION_FAILURE_CLOSE_ENCOUNTER) 
                || (((ok_system*) params)->flag & INTEGRATION_FAILURE_CLOSE_ENCOUNTER_STAR)) {
                    if (i == 0) {
                        for (int j = 0; j < SAMPLES; j++)
                            ok_free_system(bag[i]);
                        free(bag);
                        free(iwork);
                        
                        if (error != NULL) {
                            *error = INTEGRATION_FAILURE_SMALL_TIMESTEP;
                        }
                        
                        return NULL;
                    } else {
                        for (int j = i - 1; j < SAMPLES; j++) {
                            bag[j]->time = bag[j]->epoch = bag[i]->time;
                            gsl_matrix_set_all(bag[j]->xyz, INVALID_NUMBER);
                            gsl_matrix_set_all(bag[j]->orbits, INVALID_NUMBER);
                        }
                        
                        if (error != NULL) {
                            *error = INTEGRATION_FAILURE_SMALL_TIMESTEP;
                        }
                        free(iwork);
                        return bag;
                    }
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
        
        if (progress != NULL) {
            int ret = progress(i, SAMPLES, NULL, "");
            
            if (ret == PROGRESS_STOP) {
                for (int i = 0; i < SAMPLES; i++)
                        ok_free_system(bag[i]);
                
                free(bag);
                free(iwork);
                
                if (error != NULL) {
                    *error = INTEGRATION_FAILURE_STOPPED;
                }
                return NULL;
            }
        }
    }
    
    free(iwork);
    return bag;
}