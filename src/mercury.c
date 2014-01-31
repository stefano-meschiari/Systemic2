
/* This is a hand-optimized version of mercury.f's orbit routines */

#include "f2c.h"
#include <math.h>
#include "stdio.h"
/* Table of constant values */

#define pi (3.14159265358979311600)
#define twopi (6.28318530717958623200e+00)
#define piby2 (1.57079632679489655800e+00)
#define pi3by2 (4.71238898038468967400e+0)
static doublereal c_b2 = 1.;
static doublereal c_b4 = twopi;
static integer c__9 = 9;
static integer c__1 = 1;
static doublereal c_b14 = (1./3.);
static integer c__5 = 5;
static doublereal c_b34 = (1./3.);


#define mco_sine__(x, sx, cx) { sx = sin((*x)); cx = cos((*x)); }
#define mco_sine_noptr__(x, sx, cx) { sx = sin((x)); cx = cos((x)); }


/*
doublereal mco_kep__(doublereal e, doublereal M) {
        //represent function and its derivatives
        double f, f1, f2, f3;
        double d1, d2, d3;
        double accuracy_rel = 1e-7;
        double accuracy_abs = 1e-23;
        if (e == 0.) {
            return M;
        }
        double E = M + .85 * e;
        if (M > pi) {
            E = M - .85 * e;
        }
        f = E - e * sin(E) - M;
        while (fabs(f / M) > accuracy_rel && fabs(f) > accuracy_abs) {
            f1 = 1 - e * cos(E);
            f2 = -1 * f - M + E;
            f3 = -1 * f1 + 1;
            d1 = -1 * f / f1;
            d2 = -1 * f / (f1 + .5 * d1 * f2);
            d3 = -1 * f / (f1 + .5 * d2 * f2 + d2 * d2 * f3 / 6);
            E = E + d3;
            f = E - e * sin(E) - M;
        }
        return E;
    }

*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_KEP.FOR    (ErikSoft  7 July 1999) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Solves Kepler's equation for eccentricities less than one. */
/* Algorithm from A. Nijenhuis (1991) Cel. Mech. Dyn. Astron. 51, 319-330. */

/*  e = eccentricity */
/*  l = mean anomaly      (radians) */
/*  u = eccentric anomaly (   "   ) */

/* ------------------------------------------------------------------------------ */

/*<       function mco_kep (e,oldl) >*/
doublereal mco_kep__(doublereal e, doublereal oldl)
{
    /* System generated locals */
    doublereal ret_val;
    
    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), exp(doublereal);
    
    /* Local variables */
    doublereal l, p, q, x, f0, f1, f2, f3, p2, u1, u2, x2, z1, z2, z3, 
    cc, sn, ss;
    logical big;
    doublereal ome, dsn;
    logical bigg, flag__;
    doublereal sign;
    
    /*<       implicit none >*/
    
    /* Input/Outout */
    /*<       real*8 oldl,e,mco_kep >*/
    
    /* Local */
    /*<       real*8 l,pi,twopi,piby2,u1,u2,ome,sign >*/
    /*<       real*8 x,x2,sn,dsn,z1,z2,z3,f0,f1,f2,f3 >*/
    /*<       real*8 p,q,p2,ss,cc >*/
    /*<       logical flag,big,bigg >*/
    
    /* ------------------------------------------------------------------------------ */
    
    
    /* Reduce mean anomaly to lie in the range 0 < l < pi */
    /*<       if (oldl.ge.0) then >*/
    if (oldl >= 0.) {
        /*<         l = mod(oldl, twopi) >*/
        l = fmod(oldl, twopi);
        /*<       else >*/
    } else {
        /*<         l = mod(oldl, twopi) + twopi >*/
        l = fmod(oldl, twopi) + twopi;
        /*<       end if >*/
    }
    /*<       sign = 1.d0 >*/
    sign = 1.;
    /*<       if (l.gt.pi) then >*/
    if (l > pi) {
        /*<         l = twopi - l >*/
        l = twopi - l;
        /*<         sign = -1.d0 >*/
        sign = -1.;
        /*<       end if >*/
    }
    
    /*<       ome = 1.d0 - e >*/
    ome = 1. - e;
    
    /*<       if (l.ge..45d0.or.e.lt..55d0) then >*/
    if (l >= .45 || e < .55) {
        
        /* Regions A,B or C in Nijenhuis */
        /* ----------------------------- */
        
        /* Rough starting value for eccentric anomaly */
        /*<         if (l.lt.ome) then >*/
        if (l < ome) {
            /*<           u1 = ome >*/
            u1 = ome;
            /*<         else >*/
        } else {
            /*<           if (l.gt.(pi-1.d0-e)) then >*/
            if (l > pi - 1. - e) {
                /*<             u1 = (l+e*pi)/(1.d0+e) >*/
                u1 = (l + e * pi) / (e + 1.);
                /*<           else >*/
            } else {
                /*<             u1 = l + e >*/
                u1 = l + e;
                /*<           end if >*/
            }
            /*<         end if >*/
        }
        
        /* Improved value using Halley's method */
        /*<         flag = u1.gt.piby2 >*/
        flag__ = u1 > piby2;
        /*<         if (flag) then >*/
        if (flag__) {
            /*<           x = pi - u1 >*/
            x = pi - u1;
            /*<         else >*/
        } else {
            /*<           x = u1 >*/
            x = u1;
            /*<         end if >*/
        }
        /*<         x2 = x*x >*/
        x2 = x * x;
        /*<         sn = x*(1.d0 + x2*(-.16605 + x2*.00761) ) >*/
        sn = x * (x2 * (x2 * .00761 - .16605) + 1.);
        /*<         dsn = 1.d0 + x2*(-.49815 + x2*.03805) >*/
        dsn = x2 * (x2 * .03805 - .49815) + 1.;
        /*<         if (flag) dsn = -dsn >*/
        if (flag__) {
            dsn = -dsn;
        }
        /*<         f2 = e*sn >*/
        f2 = e * sn;
        /*<         f0 = u1 - f2 - l >*/
        f0 = u1 - f2 - l;
        /*<         f1 = 1.d0 - e*dsn >*/
        f1 = 1. - e * dsn;
        /*<         u2 = u1 - f0/(f1 - .5d0*f0*f2/f1) >*/
        u2 = u1 - f0 / (f1 - f0 * .5 * f2 / f1);
        /*<       else >*/
    } else {
        
        /* Region D in Nijenhuis */
        /* --------------------- */
        
        /* Rough starting value for eccentric anomaly */
        /*<         z1 = 4.d0*e + .5d0 >*/
        z1 = e * 4. + .5;
        /*<         p = ome / z1 >*/
        p = ome / z1;
        /*<         q = .5d0 * l / z1 >*/
        q = l * .5 / z1;
        /*<         p2 = p*p >*/
        p2 = p * p;
        /*<         z2 = exp( log( dsqrt( p2*p + q*q ) + q )/1.5 ) >*/
        z2 = exp(log(sqrt(p2 * p + q * q) + q) / 1.5);
        /*<         u1 = 2.d0*q / ( z2 + p + p2/z2 ) >*/
        u1 = q * 2. / (z2 + p + p2 / z2);

        /* Improved value using Newton's method */
        /*<         z2 = u1*u1 >*/
        z2 = u1 * u1;
        /*<         z3 = z2*z2 >*/
        z3 = z2 * z2;
        /*<         u2 = u1 - .075d0*u1*z3 / (ome + z1*z2 + .375d0*z3) >*/
        u2 = u1 - u1 * .075 * z3 / (ome + z1 * z2 + z3 * .375);
        /*<         u2 = l + e*u2*( 3.d0 - 4.d0*u2*u2 ) >*/
        u2 = l + e * u2 * (3. - u2 * 4. * u2);
        /*<       end if >*/
    }
    
    /* Accurate value using 3rd-order version of Newton's method */
    /* N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy! */
    
    /* First get accurate values for u2 - sin(u2) and 1 - cos(u2) */
    /*<       bigg = (u2.gt.piby2) >*/
    bigg = u2 > piby2;
    /*<       if (bigg) then >*/
    if (bigg) {
        /*<         z3 = pi - u2 >*/
        z3 = pi - u2;
        /*<       else >*/
    } else {
        /*<         z3 = u2 >*/
        z3 = u2;
        /*<       end if >*/
    }
    
    /*<       big = (z3.gt.(.5d0*piby2)) >*/
    big = z3 > (piby2 * .5);
    /*<       if (big) then >*/
    if (big) {
        /*<         x = piby2 - z3 >*/
        x = piby2 - z3;
        /*<       else >*/
    } else {
        /*<         x = z3 >*/
        x = z3;
        /*<       end if >*/
    }
    
    /*<       x2 = x*x >*/
    x2 = x * x;
    /*<       ss = 1.d0 >*/
    ss = 1.;
    /*<       cc = 1.d0 >*/
    cc = 1.;
    
    /*<        >*/
    ss = x * x2 / 6. * (1. - x2 / 20. * (1. - x2 / 42. * (1. - x2 / 
                                                               72. * (1. - x2 / 110. * (1. - x2 / 156. * (1. - x2 / 210. *
                                                                                                               (1. - x2 / 272.)))))));
    /*<        >*/
    cc = x2 / 2. * (1. - x2 / 12. * (1. - x2 / 30. * (1. - x2 / 56. * (
                                                                              1. - x2 / 90. * (1. - x2 / 132. * (1. - x2 / 182. * (1. - 
                                                                                                                                         x2 / 240. * (1. - x2 / 306.))))))));
    
    /*<       if (big) then >*/
    if (big) {
        /*<         z1 = cc + z3 - 1.d0 >*/
        z1 = cc + z3 - 1.;
        /*<         z2 = ss + z3 + 1.d0 - piby2 >*/
        z2 = ss + z3 + 1. - piby2;
        /*<       else >*/
    } else {
        /*<         z1 = ss >*/
        z1 = ss;
        /*<         z2 = cc >*/
        z2 = cc;
        /*<       end if >*/
    }
    
    /*<       if (bigg) then >*/
    if (bigg) {
        /*<         z1 = 2.d0*u2 + z1 - pi >*/
        z1 = u2 * 2. + z1 - pi;
        /*<         z2 = 2.d0 - z2 >*/
        z2 = 2. - z2;
        /*<       end if >*/
    }
    
    /*<       f0 = l - u2*ome - e*z1 >*/
    f0 = l - u2 * ome - e * z1;
    /*<       f1 = ome + e*z2 >*/
    f1 = ome + e * z2;
    /*<       f2 = .5d0*e*(u2-z1) >*/
    f2 = e * .5 * (u2 - z1);
    /*<       f3 = e/6.d0*(1.d0-z2) >*/
    f3 = e / 6. * (1. - z2);
    /*<       z1 = f0/f1 >*/
    z1 = f0 / f1;
    /*<       z2 = f0/(f2*z1+f1) >*/
    z2 = f0 / (f2 * z1 + f1);
    /*<       mco_kep = sign*( u2 + f0/((f3*z1+f2)*z2+f1) ) >*/
    ret_val = sign * (u2 + f0 / ((f3 * z1 + f2) * z2 + f1));
    
    /*<       return >*/
    return ret_val;
    /*<       end >*/
} /* mco_kep__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_X2EL.FOR    (ErikSoft  6 May 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates Keplerian orbital elements given relative coordinates and */
/* velocities, and MU = G times the sum of the masses. */

/* The elements are: q = perihelion distance */
/*                   e = eccentricity */
/*                   i = inclination */
/*                   p = longitude of perihelion (NOT argument of perihelion!!) */
/*                   n = longitude of ascending node */
/*                   l = mean anomaly (or mean longitude if e < 1.e-8) */

/* ------------------------------------------------------------------------------ */

/*<       subroutine mco_x2el (mu,x,y,z,u,v,w,q,e,i,p,n,l) >*/
/* Subroutine */ int mco_x2el__(doublereal *mu, doublereal *x, doublereal *y, 
                                doublereal *z__, doublereal *u, doublereal *v, doublereal *w, 
                                doublereal *q, doublereal *e, doublereal *i__, doublereal *p, 
                                doublereal *n, doublereal *l)
{
    /* System generated locals */
    doublereal d__1;
    
    /* Builtin functions */
    double sqrt(doublereal), acos(doublereal), atan2(doublereal, doublereal), 
    d_sign(doublereal *, doublereal *), sin(doublereal), log(
                                                             doublereal), sinh(doublereal), d_mod(doublereal *, doublereal *);
    
    /* Local variables */
    doublereal f, h__, r__, s, h2, v2, ce, cf, ci, hx, hy, hz, to, rv, 
    tmp2, bige, temp, true__;
    
    
    /*<       implicit none >*/
    /*<       integer NMAX, CMAX, NMESS >*/
    /*<       real*8 HUGE >*/
    /*<       parameter (NMAX = 2000) >*/
    /*<       parameter (CMAX = 50) >*/
    /*<       parameter (NMESS = 200) >*/
    /*<       parameter (HUGE = 9.9d29) >*/
    /* Constants: */
    
    /* DR = conversion factor from degrees to radians */
    /* K2 = Gaussian gravitational constant squared */
    /* AU = astronomical unit in cm */
    /* MSUN = mass of the Sun in g */
    
    /*<       real*8 PI,TWOPI,PIBY2,DR,K2,AU,MSUN >*/
    
    /*<       parameter (PI = 3.141592653589793d0) >*/
    /*<       parameter (TWOPI = PI * 2.d0) >*/
    /*<       parameter (PIBY2 = PI * .5d0) >*/
    /*<       parameter (DR = PI / 180.d0) >*/
    /*<       parameter (K2 = 2.959122082855911d-4) >*/
    /*<       parameter (AU = 1.4959787e13) >*/
    /*<       parameter (MSUN = 1.9891e33) >*/
    
    /* Input/Output */
    /*<       real*8 mu,q,e,i,p,n,l,x,y,z,u,v,w >*/
    
    /* Local */
    /*<       real*8 hx,hy,hz,h2,h,v2,r,rv,s,true >*/
    /*<       real*8 ci,to,temp,tmp2,bige,f,cf,ce >*/
    
    /* ------------------------------------------------------------------------------ */
    
    /*<       hx = y * w  -  z * v >*/
    hx = *y * *w - *z__ * *v;
    /*<       hy = z * u  -  x * w >*/
    hy = *z__ * *u - *x * *w;
    /*<       hz = x * v  -  y * u >*/
    hz = *x * *v - *y * *u;
    /*<       h2 = hx*hx + hy*hy + hz*hz >*/
    h2 = hx * hx + hy * hy + hz * hz;
    /*<       v2 = u * u  +  v * v  +  w * w >*/
    v2 = *u * *u + *v * *v + *w * *w;
    /*<       rv = x * u  +  y * v  +  z * w >*/
    rv = *x * *u + *y * *v + *z__ * *w;
    /*<       r = sqrt(x*x + y*y + z*z) >*/
    r__ = sqrt(*x * *x + *y * *y + *z__ * *z__);
    /*<       h = sqrt(h2) >*/
    h__ = sqrt(h2);
    /*<       s = h2 / mu >*/
    s = h2 / *mu;
    
    /* Inclination and node */
    /*<       ci = hz / h >*/
    ci = hz / h__;
    /*<       if (abs(ci).lt.1) then >*/
    if (abs(ci) < 1.) {
        /*<         i = acos (ci) >*/
        *i__ = acos(ci);
        /*<         n = atan2 (hx,-hy) >*/
        *n = atan2(hx, -hy);
        /*<         if (n.lt.0) n = n + TWOPI >*/
        if (*n < 0.) {
            *n += twopi;
        }
        /*<       else >*/
    } else {
        /*<         if (ci.gt.0) i = 0.d0 >*/
        if (ci > 0.) {
            *i__ = 0.;
        }
        /*<         if (ci.lt.0) i = PI >*/
        if (ci < 0.) {
            *i__ = pi;
        }
        /*<         n = 0.d0 >*/
        *n = 0.;
        /*<       end if >*/
    }
    
    /* Eccentricity and perihelion distance */
    /*<       temp = 1.d0 + s*(v2/mu - 2.d0/r) >*/
    temp = s * (v2 / *mu - 2. / r__) + 1.;
    /*<       if (temp.le.0) then >*/
    if (temp <= 0.) {
        /*<         e = 0.d0 >*/
        *e = 0.;
        /*<       else >*/
    } else {
        /*<         e = sqrt (temp) >*/
        *e = sqrt(temp);
        /*<       end if >*/
    }
    /*<       q = s / (1.d0 + e) >*/
    *q = s / (*e + 1.);
    
    /* True longitude */
    /*<       if (hy.ne.0) then >*/
    if (hy != 0.) {
        /*<         to = -hx/hy >*/
        to = -hx / hy;
        /*<         temp = (1.d0 - ci) * to >*/
        temp = (1. - ci) * to;
        /*<         tmp2 = to * to >*/
        tmp2 = to * to;
        /*<         true = atan2((y*(1.d0+tmp2*ci)-x*temp),(x*(tmp2+ci)-y*temp)) >*/
        true__ = atan2(*y * (tmp2 * ci + 1.) - *x * temp, *x * (tmp2 + ci) - *
                       y * temp);
        /*<       else >*/
    } else {
        /*<         true = atan2(y * ci, x) >*/
        true__ = atan2(*y * ci, *x);
        /*<       end if >*/
    }
    /*<       if (ci.lt.0) true = true + PI >*/
    if (ci < 0.) {
        true__ += pi;
    }
    
    /*<       if (e.lt.1.d-8) then >*/
    if (*e < 1e-6) {
        /*<         p = 0.d0 >*/
        *p = 0.;
        /*<         l = true >*/
        *l = true__;
        /*<       else >*/
    } else {
        /*<         ce = (v2*r - mu) / (e*mu) >*/
        ce = (v2 * r__ - *mu) / (*e * *mu);
        
        /* Mean anomaly for ellipse */
        /*<         if (e.lt.1) then >*/
        if (*e < 1.) {
            /*<           if (abs(ce).gt.1) ce = sign(1.d0,ce) >*/
            if (abs(ce) > 1.) {
                ce = d_sign(&c_b2, &ce);
            }
            /*<           bige = acos(ce) >*/
            bige = acos(ce);
            /*<           if (rv.lt.0) bige = TWOPI - bige >*/
            if (rv < 0.) {
                bige = twopi - bige;
            }
            /*<           l = bige - e*sin(bige) >*/
            *l = bige - *e * sin(bige);
            /*<         else >*/
        } else {
            
            /* Mean anomaly for hyperbola */
            /*<           if (ce.lt.1) ce = 1.d0 >*/
            if (ce < 1.) {
                ce = 1.;
            }
            /*<           bige = log( ce + sqrt(ce*ce-1.d0) ) >*/
            bige = log(ce + sqrt(ce * ce - 1.));
            /*<           if (rv.lt.0) bige = TWOPI - bige >*/
            if (rv < 0.) {
                bige = twopi - bige;
            }
            /*<           l = e*sinh(bige) - bige >*/
            *l = *e * sinh(bige) - bige;
            /*<         end if >*/
        }
        
        /* Longitude of perihelion */
        /*<         cf = (s - r) / (e*r) >*/
        cf = (s - r__) / (*e * r__);
        /*<         if (abs(cf).gt.1) cf = sign(1.d0,cf) >*/
        if (abs(cf) > 1.) {
            cf = d_sign(&c_b2, &cf);
        }
        /*<         f = acos(cf) >*/
        f = acos(cf);
        /*<         if (rv.lt.0) f = TWOPI - f >*/
        if (rv < 0.) {
            f = twopi - f;
        }
        /*<         p = true - f >*/
        *p = true__ - f;
        /*<         p = mod (p + TWOPI + TWOPI, TWOPI) >*/
        d__1 = *p + twopi + twopi;
        *p = d_mod(&d__1, &c_b4);
        /*<       end if >*/
    }
    
    /*<       if (l.lt.0) l = l + TWOPI >*/
    if (*l < 0.) {
        *l += twopi;
    }
    /*<       if (l.gt.TWOPI) l = mod (l, TWOPI) >*/
    if (*l > twopi) {
        *l = d_mod(l, &c_b4);
    }
    
    /* ------------------------------------------------------------------------------ */
    
    /*<       return >*/
    return 0;
    /*<       end >*/
} /* mco_x2el__ */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_SINE.FOR    (ErikSoft  17 April 1997) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates sin and cos of an angle X (in radians). */

/* ------------------------------------------------------------------------------ */

/*<       subroutine mco_sine (x,sx,cx) >*/

                   
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_SINH.FOR    (ErikSoft  12 June 1998) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Calculates sinh and cosh of an angle X (in radians) */

/* ------------------------------------------------------------------------------ */

/*<       subroutine mco_sinh (x,sx,cx) >*/
/* Subroutine */ int mco_sinh__(doublereal *x, doublereal *sx, doublereal *cx)
{
    /* Builtin functions */
    double sinh(doublereal), sqrt(doublereal);
    
    
    /*<       implicit none >*/
    
    /* Input/Output */
    /*<       real*8 x,sx,cx >*/
    
    /* ------------------------------------------------------------------------------ */
    
    /*<       sx = sinh(x) >*/
    *sx = sinh(*x);
    /*<       cx = sqrt (1.d0 + sx*sx) >*/
    *cx = sqrt(*sx * *sx + 1.);
    
    /* ------------------------------------------------------------------------------ */
    
    /*<       return >*/
    return 0;
    /*<       end >*/
} /* mco_sinh__ */

/* ********************************************************************** */
/*                    ORBEL_FGET.F */
/* ********************************************************************** */
/*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach. */

/*             Input: */
/*                           e ==> eccentricity anomaly. (real scalar) */
/*                        capn ==> hyperbola mean anomaly. (real scalar) */
/*             Returns: */
/*                  orbel_fget ==>  eccentric anomaly. (real scalar) */

/*     ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of */
/*           Cel. Mech. ".  Quartic convergence from Danby's book. */
/*     REMARKS: */
/*     AUTHOR: M. Duncan */
/*     DATE WRITTEN: May 11, 1992. */
/*     REVISIONS: 2/26/93 hfl */
/* ********************************************************************** */
/*< 	real*8 function orbel_fget(e,capn) >*/
doublereal orbel_fget__(doublereal *e, doublereal *capn)
{
    /* System generated locals */
    doublereal ret_val;
    
    /* Builtin functions */
    double log(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
    e_wsle(void);
    
    /* Local variables */
    doublereal f;
    integer i__;
    doublereal x, fp, dx, ech, esh, chx, fpp, tmp, shx, fppp;
    extern /* Subroutine */ int orbel_schget__(doublereal *, doublereal *, 
                                               doublereal *);
    
    /* Fortran I/O blocks */
    cilist io___60 = { 0, 6, 0, 0, 0 };
    
    
    /*<         implicit NONE     >*/
    /* ...   Version of Swift */
    /*<        real*8 VER_NUM >*/
    /*<        parameter(VER_NUM=2.0d0) >*/
    /* ...   Maximum array size */
    /*<        integer  NPLMAX, NTPMAX >*/
    /*       parameter  (NPLMAX = 21)   ! max number of planets, including the Sun */
    /*<        parameter  (NPLMAX = 51)   ! max number of planets, including the >*/
    /*<        parameter  (NTPMAX = 1001) ! max number of test particles >*/
    /* ...   Size of the test particle integer status flag */
    /*<         integer NSTATP            ! Number of status parameters >*/
    /*<         parameter  (NSTATP = 3) >*/
    /*<         integer NSTAT            ! Number of status parameters >*/
    /*<         parameter  (NSTAT = NSTATP + NPLMAX - 1)  ! include one for @ pl >*/
    /* ...   Size of the test particle integer status flag */
    /*<         integer NSTATR >*/
    /*<         parameter  (NSTATR = NSTAT)  ! io_init_tp assumes NSTAT==NSTATR >*/
    /* ...   convergence criteria for danby */
    /*<         real*8 DANBYAC , DANBYB >*/
    /*<         parameter (DANBYAC= 1.0d-14, DANBYB = 1.0d-13) >*/
    /* ...    loop limits in the Laguerre attempts */
    /*<         integer NLAG1, NLAG2 >*/
    /*<         parameter(NLAG1 = 50, NLAG2 = 400) >*/
    /* ...    A small number */
    /*<         real*8 TINY >*/
    /*<         PARAMETER(TINY=4.D-15) >*/
    /* ...    trig stuff */
    /*<         real*8 PI,TWOPI,PIBY2,DEGRAD >*/
    /*<         parameter (PI = 3.14159265358979D0) >*/
    /*<         parameter (TWOPI = 2.0D0 * PI) >*/
    /*<         parameter (PIBY2 = PI/2.0D0) >*/
    /*<         parameter (DEGRAD = 180.0D0 / PI) >*/
    /* ...  Inputs Only: */
    /*< 	real*8 e,capn >*/
    /* ...  Internals: */
    /*< 	integer i,IMAX >*/
    /*< 	real*8 tmp,x,shx,chx >*/
    /*< 	real*8 esh,ech,f,fp,fpp,fppp,dx >*/
    /*< 	PARAMETER (IMAX = 10) >*/
    /* ---- */
    /* ...  Executable code */
    /* Function to solve "Kepler's eqn" for F (here called */
    /* x) for given e and CAPN. */
    /*  begin with a guess proposed by Danby */
    /*< 	if( capn .lt. 0.d0) then >*/
    if (*capn < 0.) {
        /*< 	   tmp = -2.d0*capn/e + 1.8d0 >*/
        tmp = *capn * -2. / *e + 1.8;
        /*< 	   x = -log(tmp) >*/
        x = -log(tmp);
        /*< 	else >*/
    } else {
        /*< 	   tmp = +2.d0*capn/e + 1.8d0 >*/
        tmp = *capn * 2. / *e + 1.8;
        /*< 	   x = log( tmp) >*/
        x = log(tmp);
        /*< 	endif >*/
    }
    /*< 	orbel_fget = x >*/
    ret_val = x;
    /*< 	do i = 1,IMAX >*/
    for (i__ = 1; i__ <= 10; ++i__) {
        /*< 	  call orbel_schget(x,shx,chx) >*/
        orbel_schget__(&x, &shx, &chx);
        /*< 	  esh = e*shx >*/
        esh = *e * shx;
        /*< 	  ech = e*chx >*/
        ech = *e * chx;
        /*< 	  f = esh - x - capn >*/
        f = esh - x - *capn;
        /* 	  write(6,*) 'i,x,f : ',i,x,f */
        /*< 	  fp = ech - 1.d0   >*/
        fp = ech - 1.;
        /*< 	  fpp = esh  >*/
        fpp = esh;
        /*< 	  fppp = ech  >*/
        fppp = ech;
        /*< 	  dx = -f/fp >*/
        dx = -f / fp;
        /*< 	  dx = -f/(fp + dx*fpp/2.d0) >*/
        dx = -f / (fp + dx * fpp / 2.);
        /*< 	  dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0) >*/
        dx = -f / (fp + dx * fpp / 2. + dx * dx * fppp / 6.);
        /*< 	  orbel_fget = x + dx >*/
        ret_val = x + dx;
        /*   If we have converged here there's no point in going on */
        /*< 	  if(abs(dx) .le. TINY) RETURN >*/
        if (abs(dx) <= 4e-15) {
            return ret_val;
        }
        /*< 	  x = orbel_fget >*/
        x = ret_val;
        /*< 	enddo	 >*/
    }
    /*< 	write(6,*) 'FGET : RETURNING WITHOUT COMPLETE CONVERGENCE'  >*/
    s_wsle(&io___60);
    do_lio(&c__9, &c__1, "FGET : RETURNING WITHOUT COMPLETE CONVERGENCE", (
                                                                           ftnlen)45);
    e_wsle();
    /*< 	return >*/
    return ret_val;
    /*< 	end   ! orbel_fget >*/
} /* orbel_fget__ */

/* ------------------------------------------------------------------ */
/* ********************************************************************** */
/*                    ORBEL_FLON.F */
/* ********************************************************************** */
/*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach. */

/*             Input: */
/*                           e ==> eccentricity anomaly. (real scalar) */
/*                        capn ==> hyperbola mean anomaly. (real scalar) */
/*             Returns: */
/*                  orbel_flon ==>  eccentric anomaly. (real scalar) */

/*     ALGORITHM: Uses power series for N in terms of F and Newton,s method */
/*     REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6) */
/*     AUTHOR: M. Duncan */
/*     DATE WRITTEN: May 26, 1992. */
/*     REVISIONS: */
/* ********************************************************************** */
/*< 	real*8 function orbel_flon(e,capn) >*/
doublereal orbel_flon__(doublereal *e, doublereal *capn)
{
    /* System generated locals */
    doublereal ret_val, d__1;
    
    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
    e_wsle(void);
    double sinh(doublereal);
    
    /* Local variables */
    doublereal a, b, f;
    integer i__;
    doublereal x, a0, a1, b1, x2, fp, dx, sq, biga, bigb, diff;
    integer iflag;
    
    /* Fortran I/O blocks */
    cilist io___76 = { 0, 6, 0, 0, 0 };
    cilist io___78 = { 0, 6, 0, 0, 0 };
    cilist io___79 = { 0, 6, 0, 0, 0 };
    
    
    /*<         implicit NONE     >*/
    /* ...   Version of Swift */
    /*<        real*8 VER_NUM >*/
    /*<        parameter(VER_NUM=2.0d0) >*/
    /* ...   Maximum array size */
    /*<        integer  NPLMAX, NTPMAX >*/
    /*       parameter  (NPLMAX = 21)   ! max number of planets, including the Sun */
    /*<        parameter  (NPLMAX = 51)   ! max number of planets, including the >*/
    /*<        parameter  (NTPMAX = 1001) ! max number of test particles >*/
    /* ...   Size of the test particle integer status flag */
    /*<         integer NSTATP            ! Number of status parameters >*/
    /*<         parameter  (NSTATP = 3) >*/
    /*<         integer NSTAT            ! Number of status parameters >*/
    /*<         parameter  (NSTAT = NSTATP + NPLMAX - 1)  ! include one for @ pl >*/
    /* ...   Size of the test particle integer status flag */
    /*<         integer NSTATR >*/
    /*<         parameter  (NSTATR = NSTAT)  ! io_init_tp assumes NSTAT==NSTATR >*/
    /* ...   convergence criteria for danby */
    /*<         real*8 DANBYAC , DANBYB >*/
    /*<         parameter (DANBYAC= 1.0d-14, DANBYB = 1.0d-13) >*/
    /* ...    loop limits in the Laguerre attempts */
    /*<         integer NLAG1, NLAG2 >*/
    /*<         parameter(NLAG1 = 50, NLAG2 = 400) >*/
    /* ...    A small number */
    /*<         real*8 TINY >*/
    /*<         PARAMETER(TINY=4.D-15) >*/
    /* ...    trig stuff */
    /*<         real*8 PI,TWOPI,PIBY2,DEGRAD >*/
    /*<         parameter (PI = 3.14159265358979D0) >*/
    /*<         parameter (TWOPI = 2.0D0 * PI) >*/
    /*<         parameter (PIBY2 = PI/2.0D0) >*/
    /*<         parameter (DEGRAD = 180.0D0 / PI) >*/
    /* ...  Inputs Only: */
    /*< 	real*8 e,capn >*/
    /* ...  Internals: */
    /*< 	integer iflag,i,IMAX >*/
    /*< 	real*8 a,b,sq,biga,bigb >*/
    /*< 	real*8 x,x2 >*/
    /*< 	real*8 f,fp,dx >*/
    /*< 	real*8 diff >*/
    /*< 	real*8 a0,a1,a3,a5,a7,a9,a11 >*/
    /*< 	real*8 b1,b3,b5,b7,b9,b11 >*/
    /*< 	PARAMETER (IMAX = 10) >*/
    /*< 	PARAMETER (a11 = 156.d0,a9 = 17160.d0,a7 = 1235520.d0) >*/
    /*< 	PARAMETER (a5 = 51891840.d0,a3 = 1037836800.d0) >*/
    /*< 	PARAMETER (b11 = 11.d0*a11,b9 = 9.d0*a9,b7 = 7.d0*a7) >*/
    /*< 	PARAMETER (b5 = 5.d0*a5, b3 = 3.d0*a3) >*/
    /* ---- */
    /* ...  Executable code */
    /* Function to solve "Kepler's eqn" for F (here called */
    /* x) for given e and CAPN. Only good for smallish CAPN */
    /*< 	iflag = 0 >*/
    iflag = 0;
    /*< 	if( capn .lt. 0.d0) then >*/
    if (*capn < 0.) {
        /*< 	   iflag = 1 >*/
        iflag = 1;
        /*< 	   capn = -capn >*/
        *capn = -(*capn);
        /*< 	endif >*/
    }
    /*< 	a1 = 6227020800.d0 * (1.d0 - 1.d0/e) >*/
    a1 = (1. - 1. / *e) * 6227020800.;
    /*< 	a0 = -6227020800.d0*capn/e >*/
    a0 = *capn * -6227020800. / *e;
    /*< 	b1 = a1 >*/
    b1 = a1;
    /*  Set iflag nonzero if capn < 0., in which case solve for -capn */
    /*  and change the sign of the final answer for F. */
    /*  Begin with a reasonable guess based on solving the cubic for small F */
    /*< 	a = 6.d0*(e-1.d0)/e >*/
    a = (*e - 1.) * 6. / *e;
    /*< 	b = -6.d0*capn/e >*/
    b = *capn * -6. / *e;
    /*< 	sq = sqrt(0.25*b*b +a*a*a/27.d0) >*/
    sq = sqrt(b * .25 * b + a * a * a / 27.);
    /*< 	biga = (-0.5*b + sq)**0.3333333333333333d0 >*/
    d__1 = b * -.5 + sq;
    biga = pow_dd(&d__1, &c_b14);
    /*< 	bigb = -(+0.5*b + sq)**0.3333333333333333d0 >*/
    d__1 = b * .5 + sq;
    bigb = -pow_dd(&d__1, &c_b14);
    /*< 	x = biga + bigb >*/
    x = biga + bigb;
    /* 	write(6,*) 'cubic = ',x**3 +a*x +b */
    /*< 	orbel_flon = x >*/
    ret_val = x;
    /* If capn is tiny (or zero) no need to go further than cubic even for */
    /* e =1. */
    /*< 	if( capn .lt. TINY) go to 100 >*/
    if (*capn < 4e-15) {
        goto L100;
    }
    /*< 	do i = 1,IMAX >*/
    for (i__ = 1; i__ <= 10; ++i__) {
        /*< 	  x2 = x*x >*/
        x2 = x * x;
        /*< 	  f = a0 +x*(a1+x2*(a3+x2*(a5+x2*(a7+x2*(a9+x2*(a11+x2)))))) >*/
        f = a0 + x * (a1 + x2 * (x2 * (x2 * (x2 * (x2 * (x2 + 156.) + 17160.) 
                                             + 1235520.) + 51891840.) + 1037836800.));
        /*< 	  fp = b1 +x2*(b3+x2*(b5+x2*(b7+x2*(b9+x2*(b11 + 13.d0*x2)))))    >*/
        fp = b1 + x2 * (x2 * (x2 * (x2 * (x2 * (x2 * 13. + 1716.) + 154440.) 
                                    + 8648640.) + 259459200.) + 3113510400.);
        /*< 	  dx = -f/fp >*/
        dx = -f / fp;
        /* 	  write(6,*) 'i,dx,x,f : ' */
        /* 	  write(6,432) i,dx,x,f */
        /*< 432	  format(1x,i3,3(2x,1p1e22.15)) >*/
        /* L432: */
        /*< 	  orbel_flon = x + dx >*/
        ret_val = x + dx;
        /*   If we have converged here there's no point in going on */
        /*< 	  if(abs(dx) .le. TINY) go to 100 >*/
        if (abs(dx) <= 4e-15) {
            goto L100;
        }
        /*< 	  x = orbel_flon >*/
        x = ret_val;
        /*< 	enddo	 >*/
    }
    /* Abnormal return here - we've gone thru the loop */
    /* IMAX times without convergence */
    /*< 	if(iflag .eq. 1) then >*/
    if (iflag == 1) {
        /*< 	   orbel_flon = -orbel_flon >*/
        ret_val = -ret_val;
        /*< 	   capn = -capn >*/
        *capn = -(*capn);
        /*< 	endif >*/
    }
    /*< 	write(6,*) 'FLON : RETURNING WITHOUT COMPLETE CONVERGENCE'  >*/
    s_wsle(&io___76);
    do_lio(&c__9, &c__1, "FLON : RETURNING WITHOUT COMPLETE CONVERGENCE", (
                                                                           ftnlen)45);
    e_wsle();
    /*< 	  diff = e*sinh(orbel_flon) - orbel_flon - capn >*/
    diff = *e * sinh(ret_val) - ret_val - *capn;
    /*< 	  write(6,*) 'N, F, ecc*sinh(F) - F - N : ' >*/
    s_wsle(&io___78);
    do_lio(&c__9, &c__1, "N, F, ecc*sinh(F) - F - N : ", (ftnlen)28);
    e_wsle();
    /*< 	  write(6,*) capn,orbel_flon,diff >*/
    s_wsle(&io___79);
    do_lio(&c__5, &c__1, (char *)&(*capn), (ftnlen)sizeof(doublereal));
    do_lio(&c__5, &c__1, (char *)&ret_val, (ftnlen)sizeof(doublereal));
    do_lio(&c__5, &c__1, (char *)&diff, (ftnlen)sizeof(doublereal));
    e_wsle();
    /*< 	return >*/
    return ret_val;
    /*  Normal return here, but check if capn was originally negative */
    /*< 100	if(iflag .eq. 1) then >*/
L100:
    if (iflag == 1) {
        /*< 	   orbel_flon = -orbel_flon >*/
        ret_val = -ret_val;
        /*< 	   capn = -capn >*/
        *capn = -(*capn);
        /*< 	endif >*/
    }
    /*< 	return >*/
    return ret_val;
    /*< 	end     ! orbel_flon >*/
} /* orbel_flon__ */

/* ------------------------------------------------------------------ */
/* ********************************************************************** */
/* 	                  ORBEL_SCGET.F */
/* ********************************************************************** */
/*     PURPOSE:  Given an angle, efficiently compute sin and cos. */

/*        Input: */
/*             angle ==> angle in radians (real scalar) */

/*        Output: */
/*             sx    ==>  sin(angle)  (real scalar) */
/*             cx    ==>  cos(angle)  (real scalar) */

/*     ALGORITHM: Obvious from the code */
/*     REMARKS: The HP 700 series won't return correct answers for sin */
/*       and cos if the angle is bigger than 3e7. We first reduce it */
/*       to the range [0,2pi) and use the sqrt rather than cos (it's faster) */
/*       BE SURE THE ANGLE IS IN RADIANS - NOT DEGREES! */
/*     AUTHOR:  M. Duncan. */
/*     DATE WRITTEN:  May 6, 1992. */
/*     REVISIONS: */
/* ********************************************************************** */
/*< 	subroutine orbel_scget(angle,sx,cx) >*/
/* Subroutine */ int orbel_scget__(doublereal *angle, doublereal *sx, 
                                   doublereal *cx)
{
    /* Builtin functions */
    double sin(doublereal), sqrt(doublereal);
    
    /* Local variables */
    doublereal x;
    integer nper;
    
    /*<         implicit NONE     >*/
    /* ...   Version of Swift */
    /*<        real*8 VER_NUM >*/
    /*<        parameter(VER_NUM=2.0d0) >*/
    /* ...   Maximum array size */
    /*<        integer  NPLMAX, NTPMAX >*/
    /*       parameter  (NPLMAX = 21)   ! max number of planets, including the Sun */
    /*<        parameter  (NPLMAX = 51)   ! max number of planets, including the >*/
    /*<        parameter  (NTPMAX = 1001) ! max number of test particles >*/
    /* ...   Size of the test particle integer status flag */
    /*<         integer NSTATP            ! Number of status parameters >*/
    /*<         parameter  (NSTATP = 3) >*/
    /*<         integer NSTAT            ! Number of status parameters >*/
    /*<         parameter  (NSTAT = NSTATP + NPLMAX - 1)  ! include one for @ pl >*/
    /* ...   Size of the test particle integer status flag */
    /*<         integer NSTATR >*/
    /*<         parameter  (NSTATR = NSTAT)  ! io_init_tp assumes NSTAT==NSTATR >*/
    /* ...   convergence criteria for danby */
    /*<         real*8 DANBYAC , DANBYB >*/
    /*<         parameter (DANBYAC= 1.0d-14, DANBYB = 1.0d-13) >*/
    /* ...    loop limits in the Laguerre attempts */
    /*<         integer NLAG1, NLAG2 >*/
    /*<         parameter(NLAG1 = 50, NLAG2 = 400) >*/
    /* ...    A small number */
    /*<         real*8 TINY >*/
    /*<         PARAMETER(TINY=4.D-15) >*/
    /* ...    trig stuff */
    /*<         real*8 PI,TWOPI,PIBY2,DEGRAD >*/
    /*<         parameter (PI = 3.14159265358979D0) >*/
    /*<         parameter (TWOPI = 2.0D0 * PI) >*/
    /*<         parameter (PIBY2 = PI/2.0D0) >*/
    /*<         parameter (DEGRAD = 180.0D0 / PI) >*/
    /* ...  Inputs Only: */
    /*<         real*8 angle >*/
    /* ...  Output: */
    /*< 	real*8 sx,cx >*/
    /* ... Internals: */
    /*< 	integer nper >*/
    /*< 	real*8 x >*/
    /*< 	real*8 PI3BY2 >*/
    /*< 	parameter(PI3BY2 = 1.5d0*PI) >*/
    /* ---- */
    /* ...  Executable code */
    /*<         nper = angle/TWOPI >*/
    nper = (integer) (*angle / twopi);
    /*< 	x = angle - nper*TWOPI >*/
    x = *angle - nper * twopi;
    /*< 	if(x.lt.0.d0) then >*/
    if (x < 0.) {
        /*<            x = x + TWOPI >*/
        x += twopi;
        /*<         endif >*/
    }
    /*< 	sx = sin(x) >*/
    *sx = sin(x);
    /*< 	cx= sqrt(1.d0 - sx*sx) >*/
    *cx = sqrt(1. - *sx * *sx);
    /*< 	if( (x .gt. PIBY2) .and. (x .lt.PI3BY2)) then >*/
    if (x > piby2 && x < pi3by2) {
        /*<            cx = -cx >*/
        *cx = -(*cx);
        /*<         endif >*/
    }
    /*< 	return >*/
    return 0;
    /*< 	end   ! orbel_scget >*/
} /* orbel_scget__ */

/* ------------------------------------------------------------------- */
/* ********************************************************************** */
/* 	                  ORBEL_SCHGET.F */
/* ********************************************************************** */
/*     PURPOSE:  Given an angle, efficiently compute sinh and cosh. */

/*        Input: */
/*             angle ==> angle in radians (real scalar) */

/*        Output: */
/*             shx    ==>  sinh(angle)  (real scalar) */
/*             chx    ==>  cosh(angle)  (real scalar) */

/*     ALGORITHM: Obvious from the code */
/*     REMARKS: Based on the routine SCGET for sine's and cosine's. */
/*       We use the sqrt rather than cosh (it's faster) */
/*       BE SURE THE ANGLE IS IN RADIANS AND IT CAN'T BE LARGER THAN 300 */
/*       OR OVERFLOWS WILL OCCUR! */
/*     AUTHOR:  M. Duncan. */
/*     DATE WRITTEN:  May 6, 1992. */
/*     REVISIONS: */
/* ********************************************************************** */
/*< 	subroutine orbel_schget(angle,shx,chx) >*/
/* Subroutine */ int orbel_schget__(doublereal *angle, doublereal *shx, 
                                    doublereal *chx)
{
    /* Builtin functions */
    double sinh(doublereal), sqrt(doublereal);
    
    /*<         implicit NONE     >*/
    /* ...   Version of Swift */
    /*<        real*8 VER_NUM >*/
    /*<        parameter(VER_NUM=2.0d0) >*/
    /* ...   Maximum array size */
    /*<        integer  NPLMAX, NTPMAX >*/
    /*       parameter  (NPLMAX = 21)   ! max number of planets, including the Sun */
    /*<        parameter  (NPLMAX = 51)   ! max number of planets, including the >*/
    /*<        parameter  (NTPMAX = 1001) ! max number of test particles >*/
    /* ...   Size of the test particle integer status flag */
    /*<         integer NSTATP            ! Number of status parameters >*/
    /*<         parameter  (NSTATP = 3) >*/
    /*<         integer NSTAT            ! Number of status parameters >*/
    /*<         parameter  (NSTAT = NSTATP + NPLMAX - 1)  ! include one for @ pl >*/
    /* ...   Size of the test particle integer status flag */
    /*<         integer NSTATR >*/
    /*<         parameter  (NSTATR = NSTAT)  ! io_init_tp assumes NSTAT==NSTATR >*/
    /* ...   convergence criteria for danby */
    /*<         real*8 DANBYAC , DANBYB >*/
    /*<         parameter (DANBYAC= 1.0d-14, DANBYB = 1.0d-13) >*/
    /* ...    loop limits in the Laguerre attempts */
    /*<         integer NLAG1, NLAG2 >*/
    /*<         parameter(NLAG1 = 50, NLAG2 = 400) >*/
    /* ...    A small number */
    /*<         real*8 TINY >*/
    /*<         PARAMETER(TINY=4.D-15) >*/
    /* ...    trig stuff */
    /*<         real*8 PI,TWOPI,PIBY2,DEGRAD >*/
    /*<         parameter (PI = 3.14159265358979D0) >*/
    /*<         parameter (TWOPI = 2.0D0 * PI) >*/
    /*<         parameter (PIBY2 = PI/2.0D0) >*/
    /*<         parameter (DEGRAD = 180.0D0 / PI) >*/
    /* ...  Inputs Only: */
    /*<         real*8 angle >*/
    /* ...  Output: */
    /*< 	real*8 shx,chx >*/
    /* ---- */
    /* ...  Executable code */
    /*< 	shx = sinh(angle) >*/
    *shx = sinh(*angle);
    /*< 	chx= sqrt(1.d0 + shx*shx) >*/
    *chx = sqrt(*shx * *shx + 1.);
    /*< 	return >*/
    return 0;
    /*< 	end   ! orbel_schget >*/
} /* orbel_schget__ */

/* --------------------------------------------------------------------- */
/* ********************************************************************** */
/*                    ORBEL_ZGET.F */
/* ********************************************************************** */
/*     PURPOSE:  Solves the equivalent of Kepler's eqn. for a parabola */
/*          given Q (Fitz. notation.) */

/*             Input: */
/*                           q ==>  parabola mean anomaly. (real scalar) */
/*             Returns: */
/*                  orbel_zget ==>  eccentric anomaly. (real scalar) */

/*     ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech." */
/*     REMARKS: For a parabola we can solve analytically. */
/*     AUTHOR: M. Duncan */
/*     DATE WRITTEN: May 11, 1992. */
/*     REVISIONS: May 27 - corrected it for negative Q and use power */
/* 	      series for small Q. */
/* ********************************************************************** */
/*< 	real*8 function orbel_zget(q) >*/
doublereal orbel_zget__(doublereal *q)
{
    /* System generated locals */
    doublereal ret_val, d__1;
    
    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);
    
    /* Local variables */
    doublereal x, tmp;
    integer iflag;
    
    /*<             implicit NONE     >*/
    /* ...   Version of Swift */
    /*<        real*8 VER_NUM >*/
    /*<        parameter(VER_NUM=2.0d0) >*/
    /* ...   Maximum array size */
    /*<        integer  NPLMAX, NTPMAX >*/
    /*       parameter  (NPLMAX = 21)   ! max number of planets, including the Sun */
    /*<        parameter  (NPLMAX = 51)   ! max number of planets, including the >*/
    /*<        parameter  (NTPMAX = 1001) ! max number of test particles >*/
    /* ...   Size of the test particle integer status flag */
    /*<         integer NSTATP            ! Number of status parameters >*/
    /*<         parameter  (NSTATP = 3) >*/
    /*<         integer NSTAT            ! Number of status parameters >*/
    /*<         parameter  (NSTAT = NSTATP + NPLMAX - 1)  ! include one for @ pl >*/
    /* ...   Size of the test particle integer status flag */
    /*<         integer NSTATR >*/
    /*<         parameter  (NSTATR = NSTAT)  ! io_init_tp assumes NSTAT==NSTATR >*/
    /* ...   convergence criteria for danby */
    /*<         real*8 DANBYAC , DANBYB >*/
    /*<         parameter (DANBYAC= 1.0d-14, DANBYB = 1.0d-13) >*/
    /* ...    loop limits in the Laguerre attempts */
    /*<         integer NLAG1, NLAG2 >*/
    /*<         parameter(NLAG1 = 50, NLAG2 = 400) >*/
    /* ...    A small number */
    /*<         real*8 TINY >*/
    /*<         PARAMETER(TINY=4.D-15) >*/
    /* ...    trig stuff */
    /*<         real*8 PI,TWOPI,PIBY2,DEGRAD >*/
    /*<         parameter (PI = 3.14159265358979D0) >*/
    /*<         parameter (TWOPI = 2.0D0 * PI) >*/
    /*<         parameter (PIBY2 = PI/2.0D0) >*/
    /*<         parameter (DEGRAD = 180.0D0 / PI) >*/
    /* ...  Inputs Only: */
    /*< 	real*8 q >*/
    /* ...  Internals: */
    /*< 	integer iflag >*/
    /*< 	real*8 x,tmp >*/
    /* ---- */
    /* ...  Executable code */
    /*< 	iflag = 0 >*/
    iflag = 0;
    /*< 	if(q.lt.0.d0) then >*/
    if (*q < 0.) {
        /*< 	  iflag = 1 >*/
        iflag = 1;
        /*< 	  q = -q >*/
        *q = -(*q);
        /*< 	endif >*/
    }
    /*< 	if (q.lt.1.d-3) then >*/
    if (*q < .001) {
        /*< 	   orbel_zget = q*(1.d0 - (q*q/3.d0)*(1.d0 -q*q)) >*/
        ret_val = *q * (1. - *q * *q / 3. * (1. - *q * *q));
        /*< 	else >*/
    } else {
        /*< 	   x = 0.5d0*(3.d0*q + sqrt(9.d0*(q**2) +4.d0)) >*/
        /* Computing 2nd power */
        d__1 = *q;
        x = (*q * 3. + sqrt(d__1 * d__1 * 9. + 4.)) * .5;
        /*< 	   tmp = x**(1.d0/3.d0) >*/
        tmp = pow_dd(&x, &c_b34);
        /*< 	   orbel_zget = tmp - 1.d0/tmp >*/
        ret_val = tmp - 1. / tmp;
        /*< 	endif >*/
    }
    /*< 	if(iflag .eq.1) then >*/
    if (iflag == 1) {
        /*<            orbel_zget = -orbel_zget >*/
        ret_val = -ret_val;
        /*< 	   q = -q >*/
        *q = -(*q);
        /*< 	endif >*/
    }
    /*< 	return >*/
    return ret_val;
    /*< 	end    ! orbel_zget >*/
} /* orbel_zget__ */

/* ---------------------------------------------------------------------- */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_EL2X.FOR    (ErikSoft  7 July 1999) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates Cartesian coordinates and velocities given Keplerian orbital */
/* elements (for elliptical, parabolic or hyperbolic orbits). */

/* Based on a routine from Levison and Duncan's SWIFT integrator. */

/*  mu = grav const * (central + secondary mass) */
/*  q = perihelion distance */
/*  e = eccentricity */
/*  i = inclination                 ) */
/*  p = longitude of perihelion !!! )   in */
/*  n = longitude of ascending node ) radians */
/*  l = mean anomaly                ) */

/*  x,y,z = Cartesian positions  ( units the same as a ) */
/*  u,v,w =     "     velocities ( units the same as sqrt(mu/a) ) */

/* ------------------------------------------------------------------------------ */

/*<       subroutine mco_el2x (mu,q,e,i,p,n,l,x,y,z,u,v,w) >*/
/* Subroutine */ int mco_el2x__(doublereal mu, doublereal q, doublereal e, 
                                doublereal i__, doublereal p, doublereal n, doublereal l, 
                                doublereal *x, doublereal *y, doublereal *z__, doublereal *u, 
                                doublereal *v, doublereal *w) 
{
    /* Builtin functions */
    double sqrt(doublereal);
    
    /* Local variables */
    extern /* Subroutine */ int mco_sinh__(doublereal *, doublereal *, doublereal *);
                                                                     
    doublereal a, g, z1, z2, z3, z4, d11, d12, ce, d13, cg, d21, ci, 
    d22, d23, cn, se, sg, si, sn;
    extern doublereal orbel_zget__(doublereal *);
    doublereal temp, romes;
    extern doublereal orbel_fhybrid__(doublereal *, doublereal *);
    
    
    /*<       implicit none >*/
    /*<       integer NMAX, CMAX, NMESS >*/
    /*<       real*8 HUGE >*/
    /*<       parameter (NMAX = 2000) >*/
    /*<       parameter (CMAX = 50) >*/
    /*<       parameter (NMESS = 200) >*/
    /*<       parameter (HUGE = 9.9d29) >*/
    /* Constants: */
    
    /* DR = conversion factor from degrees to radians */
    /* K2 = Gaussian gravitational constant squared */
    /* AU = astronomical unit in cm */
    /* MSUN = mass of the Sun in g */
    
    /*<       real*8 PI,TWOPI,PIBY2,DR,K2,AU,MSUN >*/
    
    /*<       parameter (PI = 3.141592653589793d0) >*/
    /*<       parameter (TWOPI = PI * 2.d0) >*/
    /*<       parameter (PIBY2 = PI * .5d0) >*/
    /*<       parameter (DR = PI / 180.d0) >*/
    /*<       parameter (K2 = 2.959122082855911d-4) >*/
    /*<       parameter (AU = 1.4959787e13) >*/
    /*<       parameter (MSUN = 1.9891e33) >*/
    
    /* Input/Output */
    /*<       real*8 mu,q,e,i,p,n,l,x,y,z,u,v,w >*/
    
    /* Local */
    /*<       real*8 g,a,ci,si,cn,sn,cg,sg,ce,se,romes,temp >*/
    /*<       real*8 z1,z2,z3,z4,d11,d12,d13,d21,d22,d23 >*/
    /*<       real*8 mco_kep, orbel_fhybrid, orbel_zget >*/
    
    /* ------------------------------------------------------------------------------ */
    
    /* Change from longitude of perihelion to argument of perihelion */
    /*<       g = p - n >*/
    g = p - n;
    
    /* Rotation factors */
    /*<       call mco_sine (i,si,ci) >*/
    mco_sine_noptr__(i__, si, ci);
    /*<       call mco_sine (g,sg,cg) >*/
    mco_sine_noptr__(g, sg, cg);
    /*<       call mco_sine (n,sn,cn) >*/
    mco_sine_noptr__(n, sn, cn);
    /*<       z1 = cg * cn >*/
    z1 = cg * cn;
    /*<       z2 = cg * sn >*/
    z2 = cg * sn;
    /*<       z3 = sg * cn >*/
    z3 = sg * cn;
    /*<       z4 = sg * sn >*/
    z4 = sg * sn;
    /*<       d11 =  z1 - z4*ci >*/
    d11 = z1 - z4 * ci;
    /*<       d12 =  z2 + z3*ci >*/
    d12 = z2 + z3 * ci;
    /*<       d13 = sg * si >*/
    d13 = sg * si;
    /*<       d21 = -z3 - z2*ci >*/
    d21 = -z3 - z2 * ci;
    /*<       d22 = -z4 + z1*ci >*/
    d22 = -z4 + z1 * ci;
    /*<       d23 = cg * si >*/
    d23 = cg * si;
    
    /* Semi-major axis */
    /*<       a = q / (1.d0 - e) >*/
    a = q / (1. - e);
    
    /* Ellipse */
    /*<       if (e.lt.1.d0) then >*/
    if (e < 1.) {
        /*<         romes = sqrt(1.d0 - e*e) >*/
        romes = sqrt(1. - e * e);
        /*<         temp = mco_kep (e,l) >*/
        temp = mco_kep__(e, l);
        /*<         call mco_sine (temp,se,ce) >*/
        mco_sine_noptr__(temp, se, ce);
        /*<         z1 = a * (ce - e) >*/
        z1 = a * (ce - e);
        /*<         z2 = a * romes * se >*/
        z2 = a * romes * se;
        /*<         temp = sqrt(mu/a) / (1.d0 - e*ce) >*/
        temp = sqrt(mu / a) / (1. - e * ce);
        /*<         z3 = -se * temp >*/
        z3 = -se * temp;
        /*<         z4 = romes * ce * temp >*/
        z4 = romes * ce * temp;
        /*<       else >*/
    } else {
        /* Parabola */
        /*<         if (e.eq.1.d0) then >*/
        if (e == 1.) {
            /*<           ce = orbel_zget(l) >*/
            ce = orbel_zget__(&l);
            /*<           z1 = q * (1.d0 - ce*ce) >*/
            z1 = q * (1. - ce * ce);
            /*<           z2 = 2.d0 * q * ce >*/
            z2 = q * 2. * ce;
            /*<           z4 = sqrt(2.d0*mu/q) / (1.d0 + ce*ce) >*/
            z4 = sqrt(mu * 2. / q) / (ce * ce + 1.);
            /*<           z3 = -ce * z4 >*/
            z3 = -ce * z4;
            /*<         else >*/
        } else {
            /* Hyperbola */
            /*<           romes = sqrt(e*e - 1.d0) >*/
            romes = sqrt(e * e - 1.);
            /*<           temp = orbel_fhybrid(e,l) >*/
            temp = orbel_fhybrid__(&e, &l);
            /*<           call mco_sinh (temp,se,ce) >*/
            mco_sinh__(&temp, &se, &ce);
            /*<           z1 = a * (ce - e) >*/
            z1 = a * (ce - e);
            /*<           z2 = -a * romes * se >*/
            z2 = -a * romes * se;
            /*<           temp = sqrt(mu/abs(a)) / (e*ce - 1.d0) >*/
            temp = sqrt(mu / abs(a)) / (e * ce - 1.);
            /*<           z3 = -se * temp >*/
            z3 = -se * temp;
            /*<           z4 = romes * ce * temp >*/
            z4 = romes * ce * temp;
            /*<         end if >*/
        }
        /*<       endif >*/
    }
    
    /*<       x = d11*z1 + d21*z2 >*/
    *x = d11 * z1 + d21 * z2;
    /*<       y = d12*z1 + d22*z2 >*/
    *y = d12 * z1 + d22 * z2;
    /*<       z = d13*z1 + d23*z2 >*/
    *z__ = d13 * z1 + d23 * z2;
    /*<       u = d11*z3 + d21*z4 >*/
    *u = d11 * z3 + d21 * z4;
    /*<       v = d12*z3 + d22*z4 >*/
    *v = d12 * z3 + d22 * z4;
    /*<       w = d13*z3 + d23*z4 >*/
    *w = d13 * z3 + d23 * z4;
    
    /* ------------------------------------------------------------------------------ */
    
    /*<       return >*/
    return 0;
    /*<       end >*/
} 
                 
                 /* mco_el2x__ */

/* ********************************************************************** */
/*                    ORBEL_FHYBRID.F */
/* ********************************************************************** */
/*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach. */

/*             Input: */
/*                           e ==> eccentricity anomaly. (real scalar) */
/*                           n ==> hyperbola mean anomaly. (real scalar) */
/*             Returns: */
/*               orbel_fhybrid ==>  eccentric anomaly. (real scalar) */

/*     ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON */
/* 	         For larger N, uses FGET */
/*     REMARKS: */
/*     AUTHOR: M. Duncan */
/*     DATE WRITTEN: May 26,1992. */
/*     REVISIONS: */
/*     REVISIONS: 2/26/93 hfl */
/* ********************************************************************** */
/*< 	real*8 function orbel_fhybrid(e,n) >*/
doublereal orbel_fhybrid__(doublereal *e, doublereal *n)
{
    /* System generated locals */
    doublereal ret_val;
    
    /* Local variables */
    extern doublereal orbel_fget__(doublereal *, doublereal *), orbel_flon__(
                                                                             doublereal *, doublereal *);
    doublereal abn;
    
    /*<         implicit NONE     >*/
    /* ...   Version of Swift */
    /*<        real*8 VER_NUM >*/
    /*<        parameter(VER_NUM=2.0d0) >*/
    /* ...   Maximum array size */
    /*<        integer  NPLMAX, NTPMAX >*/
    /*       parameter  (NPLMAX = 21)   ! max number of planets, including the Sun */
    /*<        parameter  (NPLMAX = 51)   ! max number of planets, including the >*/
    /*<        parameter  (NTPMAX = 1001) ! max number of test particles >*/
    /* ...   Size of the test particle integer status flag */
    /*<         integer NSTATP            ! Number of status parameters >*/
    /*<         parameter  (NSTATP = 3) >*/
    /*<         integer NSTAT            ! Number of status parameters >*/
    /*<         parameter  (NSTAT = NSTATP + NPLMAX - 1)  ! include one for @ pl >*/
    /* ...   Size of the test particle integer status flag */
    /*<         integer NSTATR >*/
    /*<         parameter  (NSTATR = NSTAT)  ! io_init_tp assumes NSTAT==NSTATR >*/
    /* ...   convergence criteria for danby */
    /*<         real*8 DANBYAC , DANBYB >*/
    /*<         parameter (DANBYAC= 1.0d-14, DANBYB = 1.0d-13) >*/
    /* ...    loop limits in the Laguerre attempts */
    /*<         integer NLAG1, NLAG2 >*/
    /*<         parameter(NLAG1 = 50, NLAG2 = 400) >*/
    /* ...    A small number */
    /*<         real*8 TINY >*/
    /*<         PARAMETER(TINY=4.D-15) >*/
    /* ...    trig stuff */
    /*<         real*8 PI,TWOPI,PIBY2,DEGRAD >*/
    /*<         parameter (PI = 3.14159265358979D0) >*/
    /*<         parameter (TWOPI = 2.0D0 * PI) >*/
    /*<         parameter (PIBY2 = PI/2.0D0) >*/
    /*<         parameter (DEGRAD = 180.0D0 / PI) >*/
    /* ...  Inputs Only: */
    /*< 	real*8 e,n >*/
    /* ...  Internals: */
    /*< 	real*8 abn >*/
    /*<         real*8 orbel_flon,orbel_fget >*/
    /* ---- */
    /* ...  Executable code */
    /*< 	abn = n >*/
    abn = *n;
    /*< 	if(n.lt.0.d0) abn = -abn >*/
    if (*n < 0.) {
        abn = -abn;
    }
    /*< 	if(abn .lt. 0.636d0*e -0.6d0) then >*/
    if (abn < *e * .636 - .6) {
        /*< 	  orbel_fhybrid = orbel_flon(e,n) >*/
        ret_val = orbel_flon__(e, n);
        /*< 	else  >*/
    } else {
        /*< 	  orbel_fhybrid = orbel_fget(e,n) >*/
        ret_val = orbel_fget__(e, n);
        /*< 	endif    >*/
    }
    /*< 	return >*/
    return ret_val;
    /*< 	end  ! orbel_fhybrid >*/
} /* orbel_fhybrid__ */

