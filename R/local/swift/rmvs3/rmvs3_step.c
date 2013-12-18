/* rmvs3_step.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__10 = 10;

/* ************************************************************************* */
/*                            RMVS3_STEP.F */
/* ************************************************************************* */
/* This subroutine takes a step in helio coord. */
/* both massive and test particles */
/* INCLUDES close approuches between test particles and planets */

/* VERSION 3 of RMVS */

/*             Input: */
/*                 i1st          ==>  = 0 if first step; = 1 not (int scalar) */
/*                 time          ==>  current time (real scalar) */
/*                 nbod          ==>  number of massive bodies (int scalar) */
/*                 ntp            ==>  number of massive bodies (int scalar) */
/*                 mass          ==>  mass of bodies (real array) */
/*                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4 */
/*                                     (real scalars) */
/*                 xh,yh,zh      ==>  initial planet position in helio coord */
/*                                    (real arrays) */
/*                 vxh,vyh,vzh   ==>  initial planet velocity in helio coord */
/*                                    (real arrays) */
/*                 xht,yht,zht    ==>  initial tp position in helio coord */
/*                                      (real arrays) */
/*                 vxht,vyht,vzht ==>  initial tp velocity in helio coord */
/*                                        (real arrays) */
/*                 istat           ==>  status of the test paricles */
/*                                      (2d integer array) */
/*                                      istat(i,1) = 0 ==> active:  = 1 not */
/*                                      istat(i,2) = -1 ==> Danby did not work */
/*                 rstat           ==>  status of the test paricles */
/*                                      (2d real array) */
/*                 dt            ==>  time step */
/*             Output: */
/*                 xh,yh,zh      ==>  final planet position in helio coord */
/*                                       (real arrays) */
/*                 vxh,vyh,vzh   ==>  final planet velocity in helio coord */
/*                                       (real arrays) */
/*                 xht,yht,zht    ==>  final tp position in helio coord */
/*                                       (real arrays) */
/*                 vxht,vyht,vzht ==>  final tp position in helio coord */
/*                                       (real arrays) */


/* Remarks: Based on rmvs2_step.f */
/* Authors:  Hal Levison */
/* Date:    7/10/96 */
/* Last revision: */
/* Subroutine */ int rmvs3_step__(integer *i1st, doublereal *time, integer *
	nbod, integer *ntp, doublereal *mass, doublereal *j2rp2, doublereal *
	j4rp4, doublereal *xh, doublereal *yh, doublereal *zh, doublereal *
	vxh, doublereal *vyh, doublereal *vzh, doublereal *xht, doublereal *
	yht, doublereal *zht, doublereal *vxht, doublereal *vyht, doublereal *
	vzht, integer *istat, doublereal *rstat, doublereal *dt)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    extern /* Subroutine */ int step_kdk__(integer *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *);
    static integer istattmp[53053]	/* was [1001][53] */, i__, j;
    extern /* Subroutine */ int rmvs3_chk__(integer *, integer *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *);
    static integer np;
    static doublereal rts;
    extern /* Subroutine */ int step_kdk_pl__(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), step_kdk_tp__(integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *)
	    ;
    static integer ienc[1001], nenc[51];
    static doublereal xbeg[51], ybeg[51], zbeg[51], peri[1001], xend[51], 
	    yend[51], zend[51], xtmp[510]	/* was [51][10] */, ytmp[510]	
	    /* was [51][10] */, ztmp[510]	/* was [51][10] */;
    extern /* Subroutine */ int rmvs3_interp__(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *);
    static integer i1sto, icflg;
    static doublereal vxbeg[51], vybeg[51], vzbeg[51], vxend[51], vyend[51], 
	    vzend[51], vxtmp[510]	/* was [51][10] */, vytmp[510]	/* 
	    was [51][10] */, vztmp[510]	/* was [51][10] */;
    static integer i1stpl, i1sttp, itpenc[51051]	/* was [1001][51] */, 
	    isperi[1001];
    extern /* Subroutine */ int rmvs3_step_out__(integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *);

/* ************************************************************************* */
/*                        SWIFT.INC */
/* ************************************************************************* */
/* Include file for SWIFT */

/* Author:  Hal Levison */
/* Date:    2/2/93 */
/* Last revision: 3/7/93 */
/* ...   Version of Swift */
/* you got it baby */
/* ...   Maximum array size */
/*       parameter  (NPLMAX = 21)   ! max number of planets, including the Sun */
/* max number of planets, including the */
/* ...   Size of the test particle integer status flag */
/* max number of test particles */
/* Number of status parameters */
/* Number of status parameters */
/* ...   Size of the test particle integer status flag */
/* include one for @ pl */
/* ...   convergence criteria for danby */
/* io_init_tp assumes NSTAT==NSTATR */
/* ...    loop limits in the Laguerre attempts */
/* ...    A small number */
/* ...    trig stuff */
/* ------------------------------------------------------------------------- */
/* ...  Inputs Only: */
/* ************************************************************************* */
/*                        RMVS.INC */
/* ************************************************************************* */
/* Include file for the encounter subroutines */

/* Author:  Hal Levison */
/* Date:    3/7/93 */
/* Last revision: */
/* ...	scale factor for hill's sphere to take shorter time step */
/*       parameter (RHSCALE=5.0) */
/* ...	scale factor for hill's sphere to go to planet centric coord */
/*       parameter (RHPSCALE=1.5) */
/* ..    ratio of the number of time steps in the */
/* ...        encounter (helocentric) vs normal */
/*       parameter (NTENC=100) */
/* ..    ratio of the number of time steps in the */
/* ...        planetcentric encounter  vs heliocentric */
/*       parameter (NTPHENC=5) */
/* ..    ratio of the number of time steps in the */
/* ...        encounter (planetcentric) vs normal */
/* ------------------------------------------------------------------------- */
/* ...  Inputs and Outputs: */
/* ...  Internals */
/* ---- */
/* ...  Executable code */
    /* Parameter adjustments */
    --vzh;
    --vyh;
    --vxh;
    --zh;
    --yh;
    --xh;
    --mass;
    --vzht;
    --vyht;
    --vxht;
    --zht;
    --yht;
    --xht;
    istat -= 1002;
    rstat -= 1002;

    /* Function Body */
    i1sttp = *i1st;
    i1sto = *i1st;
/* ...  are there any encounters? */
    rts = 12.25;
    rmvs3_chk__(nbod, ntp, &mass[1], &xh[1], &yh[1], &zh[1], &vxh[1], &vyh[1],
	     &vzh[1], &xht[1], &yht[1], &zht[1], &vxht[1], &vyht[1], &vzht[1],
	     &istat[1002], dt, &rts, &icflg, nenc, itpenc, ienc);
/* ...     nenc and itpenc not used here! */
/* .... if not just do a normal step and leave */
    if (icflg == 0) {
	step_kdk__(i1st, time, nbod, ntp, &mass[1], j2rp2, j4rp4, &xh[1], &yh[
		1], &zh[1], &vxh[1], &vyh[1], &vzh[1], &xht[1], &yht[1], &zht[
		1], &vxht[1], &vyht[1], &vzht[1], &istat[1002], &rstat[1002], 
		dt);
	i__1 = *ntp;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (istat[i__ + 1001] == 0) {
		istat[i__ + 2002] = 0;
	    }
	}
	return 0;
/*  NOTE AN EXIT */
    }
/* ...  ENCOUNTER STUFF FROM HERE ON!!!!! */
/* ...  save initial x and v of planets if there are planocentric enc */
    i__1 = *nbod;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xbeg[i__ - 1] = xh[i__];
	ybeg[i__ - 1] = yh[i__];
	zbeg[i__ - 1] = zh[i__];
	vxbeg[i__ - 1] = vxh[i__];
	vybeg[i__ - 1] = vyh[i__];
	vzbeg[i__ - 1] = vzh[i__];
    }
/* ... do a full step for the planets */
    i1stpl = *i1st;
    step_kdk_pl__(&i1stpl, nbod, &mass[1], j2rp2, j4rp4, &xh[1], &yh[1], &zh[
	    1], &vxh[1], &vyh[1], &vzh[1], dt);
/* ...  save the final position and velocity of planets */
    i__1 = *nbod;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xend[i__ - 1] = xh[i__];
	yend[i__ - 1] = yh[i__];
	zend[i__ - 1] = zh[i__];
	vxend[i__ - 1] = vxh[i__];
	vyend[i__ - 1] = vyh[i__];
	vzend[i__ - 1] = vzh[i__];
    }
/* ...  Now do the interpolation for intermediate steps */
    rmvs3_interp__(nbod, xbeg, ybeg, zbeg, vxbeg, vybeg, vzbeg, xend, yend, 
	    zend, vxend, vyend, vzend, dt, &mass[1], &c__10, xtmp, ytmp, ztmp,
	     vxtmp, vytmp, vztmp);
/* ...  do the integration */
    rmvs3_step_out__(i1st, nbod, ntp, &mass[1], j2rp2, j4rp4, xbeg, ybeg, 
	    zbeg, vxbeg, vybeg, vzbeg, xtmp, ytmp, ztmp, vxtmp, vytmp, vztmp, 
	    &xht[1], &yht[1], &zht[1], &vxht[1], &vyht[1], &vzht[1], &istat[
	    1002], ienc, dt, isperi, peri);
/* ...  As of this point all the test particles that are involved in an */
/* ...  encounter have been moved.  But not the ones that have not. */
/* ...  so move those,  BUT NOT the onces in the encounter */
/* ...  make a temporary istat array so only they are active */
    i__1 = *ntp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (istat[i__ + 1001] == 0) {
	    if (ienc[i__ - 1] != 0) {
		istattmp[i__ - 1] = 1;
/* don't integrate it */
	    } else {
		istattmp[i__ - 1] = 0;
/* integrate it */
	    }
	    for (j = 2; j <= 53; ++j) {
		istattmp[i__ + j * 1001 - 1002] = 0;
	    }
	} else {
	    istattmp[i__ - 1] = 1;
/* don't integrate it */
	}
    }
/* ...  do a full step */
    i1sto = 0;
/* we need to recalculate accel arrays */
    step_kdk_tp__(&i1sto, nbod, ntp, &mass[1], j2rp2, j4rp4, xbeg, ybeg, zbeg,
	     xend, yend, zend, &xht[1], &yht[1], &zht[1], &vxht[1], &vyht[1], 
	    &vzht[1], istattmp, dt);
/* ...  fix up the istat array */
    i__1 = *ntp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (istattmp[i__ + 1000] != 0) {
/* danby screwed up */
	    istat[i__ + 1001] = 1;
	    for (j = 2; j <= 53; ++j) {
		istat[i__ + j * 1001] = istattmp[i__ + j * 1001 - 1002];
	    }
	}
    }
/* ...  put the enc info into istat */
    i__1 = *ntp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (istat[i__ + 1001] == 0) {
	    if (ienc[i__ - 1] < 0) {
		istat[i__ + 2002] = ienc[i__ - 1];
		istat[i__ + 3003] = (i__2 = ienc[i__ - 1], abs(i__2));
		if (isperi[i__ - 1] == 0) {
		    rstat[i__ + 3003] = peri[i__ - 1];
		    np = (i__2 = ienc[i__ - 1], abs(i__2)) + 2;
		    if (rstat[i__ + np * 1001] == 0.) {
			rstat[i__ + np * 1001] = peri[i__ - 1];
		    } else {
/* Computing MIN */
			d__1 = rstat[i__ + np * 1001], d__2 = peri[i__ - 1];
			rstat[i__ + np * 1001] = min(d__1,d__2);
		    }
		} else {
		    rstat[i__ + 3003] = 0.;
		}
	    } else if (ienc[i__ - 1] > 0) {
		istat[i__ + 2002] = ienc[i__ - 1];
	    } else {
		istat[i__ + 2002] = 0;
	    }
	}
    }
/* ...  we MUST make sure that the saved accel arrays are ok */
/* ...  calculate them again */
    *i1st = 0;
    return 0;
} /* rmvs3_step__ */

