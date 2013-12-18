/* swift_bs.f -- translated by f2c (version 20100827).
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

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__0 = 0;

/* ********************************************************************** */
/* 		      SWIFT_BS.F */
/* ********************************************************************** */

/*                 INCLUDES CLOSE ENCOUNTERS */
/*                 To run, need 3 input files. The code prompts for */
/*                 the file names, but examples are : */

/*                   parameter file like       param.in */
/* 		    planet file like          pl.in */
/*                   test particle file like   tp.in */

/* Authors:  Hal Levison \& Martin Duncan */
/* Date:    8/5/93 */
/* Last revision: 12/27/96 */
/* ************************************************************************* */
/*                        SWIFT.INC */
/* ************************************************************************* */
/* Include file for SWIFT */

/* Author:  Hal Levison */
/* Date:    2/2/93 */
/* Last revision: 3/7/93 */
/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_999[] = "(a)";
    static char fmt_998[] = "(\002 Time = \002,1p1e12.5,\002: fraction done "
	    "= \002,0pf5.3,\002: Number of active tp =\002,i4)";

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_rsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_rsfe(void), s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static char inplfile[80], intpfile[80];
    extern /* Subroutine */ int io_discard_write__(integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, char *, char *, 
	    integer *, ftnlen, ftnlen), io_write_frame_r__(doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, char *, integer *, char *, 
	    ftnlen, ftnlen);
    static doublereal t;
    static char inparfile[80];
    static doublereal t0;
    extern /* Subroutine */ int anal_jacobi_write__(doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, char *, ftnlen);
    static char fopenstat[80];
    extern /* Subroutine */ int util_exit__(integer *), anal_energy_write__(
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, char *, doublereal *, 
	    ftnlen);
    static doublereal dt, xh[51], yh[51], zh[51];
    extern /* Subroutine */ int io_init_pl__(char *, shortlogical *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, ftnlen), io_dump_pl__(char *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, shortlogical *, integer 
	    *, doublereal *, doublereal *, doublereal *, ftnlen), 
	    io_init_tp__(char *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, ftnlen), io_dump_tp__(char *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, ftnlen);
    static integer iub, iud, iue, iuj, ntp;
    static doublereal xht[1001], yht[1001], zht[1001], vxh[51], vyh[51], vzh[
	    51];
    static integer i1st;
    static doublereal eoff;
    static integer nbod;
    static doublereal j2rp2, mass[51], j4rp4, rmin, qmin, rmax, vxht[1001], 
	    vyht[1001], vzht[1001], tout;
    extern /* Subroutine */ int util_version__(void);
    static doublereal tfrac;
    static integer nleft, istat[53053]	/* was [1001][53] */;
    static doublereal tdump, rmaxu, rstat[53053]	/* was [1001][53] */, 
	    dtout, rplsq[51];
    extern /* Subroutine */ int io_init_param__(char *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    shortlogical *, char *, char *, ftnlen, ftnlen, ftnlen);
    static doublereal tstop;
    extern /* Subroutine */ int io_dump_param__(char *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    shortlogical *, char *, ftnlen, ftnlen);
    static shortlogical lclose;
    static doublereal dtdump;
    extern /* Subroutine */ int io_write_frame__(doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, char *, integer *, char *, ftnlen, 
	    ftnlen);
    static integer iflgchk;
    extern /* Subroutine */ int discard_(doublereal *, doublereal *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, shortlogical *, doublereal *, integer *, doublereal 
	    *), bs_step__(integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *);
    static char outfile[80];

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };
    static cilist io___2 = { 0, 5, 0, fmt_999, 0 };
    static cilist io___17 = { 0, 6, 0, 0, 0 };
    static cilist io___18 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 5, 0, fmt_999, 0 };
    static cilist io___32 = { 0, 6, 0, 0, 0 };
    static cilist io___33 = { 0, 5, 0, fmt_999, 0 };
    static cilist io___54 = { 0, 6, 0, 0, 0 };
    static cilist io___56 = { 0, 6, 0, fmt_998, 0 };


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
/* ----- */
/* ...    Executable code */
/* ...    print version number */
    util_version__();
/* Get data for the run and the test particles */
    s_wsle(&io___1);
    do_lio(&c__9, &c__1, "Enter name of parameter data file : ", (ftnlen)36);
    e_wsle();
    s_rsfe(&io___2);
    do_fio(&c__1, inparfile, (ftnlen)80);
    e_rsfe();
    io_init_param__(inparfile, &t0, &tstop, &dt, &dtout, &dtdump, &iflgchk, &
	    rmin, &rmax, &rmaxu, &qmin, &lclose, outfile, fopenstat, (ftnlen)
	    80, (ftnlen)80, (ftnlen)80);
/* Prompt and read name of planet data file */
    s_wsle(&io___17);
    do_lio(&c__9, &c__1, " ", (ftnlen)1);
    e_wsle();
    s_wsle(&io___18);
    do_lio(&c__9, &c__1, "Enter name of planet data file : ", (ftnlen)33);
    e_wsle();
    s_rsfe(&io___19);
    do_fio(&c__1, inplfile, (ftnlen)80);
    e_rsfe();
    io_init_pl__(inplfile, &lclose, &iflgchk, &nbod, mass, xh, yh, zh, vxh, 
	    vyh, vzh, rplsq, &j2rp2, &j4rp4, (ftnlen)80);
/* Get data for the run and the test particles */
    s_wsle(&io___32);
    do_lio(&c__9, &c__1, "Enter name of test particle data file : ", (ftnlen)
	    40);
    e_wsle();
    s_rsfe(&io___33);
    do_fio(&c__1, intpfile, (ftnlen)80);
    e_rsfe();
    io_init_tp__(intpfile, &ntp, xht, yht, zht, vxht, vyht, vzht, istat, 
	    rstat, (ftnlen)80);
/* Initialize initial time and times for first output and first dump */
    t = t0;
    tout = t0 + dtout;
    tdump = t0 + dtdump;
    iub = 20;
    iuj = 30;
    iud = 40;
    iue = 60;
/* ...    Do the initial io write */
    if (bit_test(iflgchk,0)) {
/* bit 0 is set */
	io_write_frame__(&t0, &nbod, &ntp, mass, xh, yh, zh, vxh, vyh, vzh, 
		xht, yht, zht, vxht, vyht, vzht, istat, outfile, &iub, 
		fopenstat, (ftnlen)80, (ftnlen)80);
    }
    if (bit_test(iflgchk,1)) {
/* bit 1 is set */
	io_write_frame_r__(&t0, &nbod, &ntp, mass, xh, yh, zh, vxh, vyh, vzh, 
		xht, yht, zht, vxht, vyht, vzht, istat, outfile, &iub, 
		fopenstat, (ftnlen)80, (ftnlen)80);
    }
    if (bit_test(iflgchk,2)) {
/* bit 2 is set */
	eoff = 0.;
	anal_energy_write__(&t0, &nbod, mass, &j2rp2, &j4rp4, xh, yh, zh, vxh,
		 vyh, vzh, &iue, fopenstat, &eoff, (ftnlen)80);
    }
    if (bit_test(iflgchk,3)) {
/* bit 3 is set */
	anal_jacobi_write__(&t0, &nbod, &ntp, mass, xh, yh, zh, vxh, vyh, vzh,
		 xht, yht, zht, vxht, vyht, vzht, istat, &c__2, &iuj, 
		fopenstat, (ftnlen)80);
    }
/* ...    must initize discard io routine */
    if (bit_test(iflgchk,4)) {
/* bit 4 is set */
	io_discard_write__(&c__0, &t, &nbod, &ntp, xh, yh, zh, vxh, vyh, vzh, 
		xht, yht, zht, vxht, vyht, vzht, istat, rstat, &iud, "discar"
		"d.out", fopenstat, &nleft, (ftnlen)11, (ftnlen)80);
    }
    nleft = ntp;
    i1st = 0;
/* ***************here's the big loop ************************************* */
    s_wsle(&io___54);
    do_lio(&c__9, &c__1, " ************** MAIN LOOP ****************** ", (
	    ftnlen)45);
    e_wsle();
    while(t <= tstop && (ntp == 0 || nleft > 0)) {
	bs_step__(&i1st, &t, &nbod, &ntp, mass, &j2rp2, &j4rp4, xh, yh, zh, 
		vxh, vyh, vzh, xht, yht, zht, vxht, vyht, vzht, istat, rstat, 
		&dt);
	t += dt;
	if (bit_test(iflgchk,4)) {
/* bit 4 is set */
	    discard_(&t, &dt, &nbod, &ntp, mass, xh, yh, zh, vxh, vyh, vzh, 
		    xht, yht, zht, vxht, vyht, vzht, &rmin, &rmax, &rmaxu, &
		    qmin, &lclose, rplsq, istat, rstat);
	    io_discard_write__(&c__1, &t, &nbod, &ntp, xh, yh, zh, vxh, vyh, 
		    vzh, xht, yht, zht, vxht, vyht, vzht, istat, rstat, &iud, 
		    "discard.out", fopenstat, &nleft, (ftnlen)11, (ftnlen)80);
	} else {
	    nleft = ntp;
	}
/* if it is time, output orb. elements, */
	if (t >= tout) {
	    if (bit_test(iflgchk,0)) {
/* bit 0 is set */
		io_write_frame__(&t, &nbod, &ntp, mass, xh, yh, zh, vxh, vyh, 
			vzh, xht, yht, zht, vxht, vyht, vzht, istat, outfile, 
			&iub, fopenstat, (ftnlen)80, (ftnlen)80);
	    }
	    if (bit_test(iflgchk,1)) {
/* bit 1 is set */
		io_write_frame_r__(&t, &nbod, &ntp, mass, xh, yh, zh, vxh, 
			vyh, vzh, xht, yht, zht, vxht, vyht, vzht, istat, 
			outfile, &iub, fopenstat, (ftnlen)80, (ftnlen)80);
	    }
	    tout += dtout;
	}
/* If it is time, do a dump */
	if (t >= tdump) {
	    tfrac = (t - t0) / (tstop - t0);
	    s_wsfe(&io___56);
	    do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&tfrac, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&nleft, (ftnlen)sizeof(integer));
	    e_wsfe();
	    io_dump_pl__("dump_pl.dat", &nbod, mass, xh, yh, zh, vxh, vyh, 
		    vzh, &lclose, &iflgchk, rplsq, &j2rp2, &j4rp4, (ftnlen)11)
		    ;
	    io_dump_tp__("dump_tp.dat", &ntp, xht, yht, zht, vxht, vyht, vzht,
		     istat, rstat, (ftnlen)11);
	    io_dump_param__("dump_param.dat", &t, &tstop, &dt, &dtout, &
		    dtdump, &iflgchk, &rmin, &rmax, &rmaxu, &qmin, &lclose, 
		    outfile, (ftnlen)14, (ftnlen)80);
	    tdump += dtdump;
	    if (bit_test(iflgchk,2)) {
/* bit 2 is set */
		anal_energy_write__(&t, &nbod, mass, &j2rp2, &j4rp4, xh, yh, 
			zh, vxh, vyh, vzh, &iue, fopenstat, &eoff, (ftnlen)80)
			;
	    }
	    if (bit_test(iflgchk,3)) {
/* bit 3 is set */
		anal_jacobi_write__(&t, &nbod, &ntp, mass, xh, yh, zh, vxh, 
			vyh, vzh, xht, yht, zht, vxht, vyht, vzht, istat, &
			c__2, &iuj, fopenstat, (ftnlen)80);
	    }
	}
    }
/* ********** end of the big loop from time 't0' to time 'tstop' */
/* Do a final dump for possible resumption later */
    io_dump_pl__("dump_pl.dat", &nbod, mass, xh, yh, zh, vxh, vyh, vzh, &
	    lclose, &iflgchk, rplsq, &j2rp2, &j4rp4, (ftnlen)11);
    io_dump_tp__("dump_tp.dat", &ntp, xht, yht, zht, vxht, vyht, vzht, istat, 
	    rstat, (ftnlen)11);
    io_dump_param__("dump_param.dat", &t, &tstop, &dt, &dtout, &dtdump, &
	    iflgchk, &rmin, &rmax, &rmaxu, &qmin, &lclose, outfile, (ftnlen)
	    14, (ftnlen)80);
    util_exit__(&c__0);
    return 0;
} /* MAIN__ */

