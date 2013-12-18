//
//  mercury.h
//  Systemic Console
//
//  Created by Stefano Meschiari on 12/3/10.
//  Copyright (c) 2010 __MyCompanyName__. All rights reserved.
//
#include "f2c.h"

int mco_el2x__(doublereal mu, doublereal q, doublereal e, 
               doublereal i__, doublereal p, doublereal n, doublereal l, 
               doublereal* x, doublereal* y, doublereal* z__, doublereal* u, 
               doublereal* v, doublereal* w);

int mco_x2el__(doublereal* mu, doublereal* x, doublereal* y, 
               doublereal* z__, doublereal* u, doublereal* v, doublereal* w, 
               doublereal* q, doublereal* e, doublereal* i__, doublereal* p, 
               doublereal* n, doublereal* l);
doublereal mco_kep__(doublereal e, doublereal oldl);