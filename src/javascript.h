#define ELEMENT 0
#define PARAM 1
#define DATA 2

#define ALL -1
#define LENGTH -2

#define JS_I_START 0
#define JS_I_STEP 1
#define JS_I_GET 2
#define JS_I_END 3
#define JS_I_ENDREACHED -1

double K_getDataAt(ok_kernel* k, int subset, int row, int column);
void K_setDataAt(ok_kernel* k, int subset, int row, int column, double val);
double K_getRVLine(ok_kernel* k, int row, int col);
double K_getPeriodogramAt(ok_kernel* k, int row, int col);
int K_minimizeWithTimeout(ok_kernel* k, int timeout);
double K_getPhasedRVLine(ok_kernel* k, int planet, int row, int column);
double K_getPhasedDataForPlanet(ok_kernel* k, int planet, int row, int column);
double K_integrateForward(ok_kernel* k, const int mode, const double nyears,
        const int row, const int col);