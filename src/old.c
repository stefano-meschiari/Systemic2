
/* Subroutine */ int mco_sine__(doublereal *x, doublereal *sx, doublereal *cx)
{
    /* Builtin functions */
    double d_mod(doublereal *, doublereal *), cos(doublereal), sqrt(
                                                                    doublereal);
    
    /* Local variables */
    doublereal pi, twopi;
    
    
    /*<       implicit none >*/
    
    /* Input/Output */
    /*<       real*8 x,sx,cx >*/
    
    /* Local */
    /*<       real*8 pi,twopi >*/
    
    /* ------------------------------------------------------------------------------ */
    
    /*<       pi = 3.141592653589793d0 >*/
    pi = 3.141592653589793;
    /*<       twopi = 2.d0 * pi >*/
    twopi = pi * 2.;
    
    /*<       if (x.gt.0) then >*/
    if (*x > 0.) {
        /*<         x = mod(x,twopi) >*/
        *x = d_mod(x, &twopi);
        /*<       else >*/
    } else {
        /*<         x = mod(x,twopi) + twopi >*/
        *x = d_mod(x, &twopi) + twopi;
        /*<       end if >*/
    }
    
    /*<       cx = cos(x) >*/
    *cx = cos(*x);
    
    /*<       if (x.gt.pi) then >*/
    if (*x > pi) {
        /*<         sx = -sqrt(1.d0 - cx*cx) >*/
        *sx = -sqrt(1. - *cx * *cx);
        /*<       else >*/
    } else {
        /*<         sx =  sqrt(1.d0 - cx*cx) >*/
        *sx = sqrt(1. - *cx * *cx);
        /*<       end if >*/
    }
    
    /* ------------------------------------------------------------------------------ */
    
    /*<       return >*/
    return 0;
    /*<       end >*/
} /* mco_sine__ */

