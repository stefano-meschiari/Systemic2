#define ELEMENT 0
#define PARAM 1
#define DATA 2

#define ALL -1
#define LENGTH -2

#define START 0
#define CONTINUE 1
#define GET 2

double K_getDataAt(ok_kernel* k, int subset, int row, int column);
void K_setDataAt(ok_kernel* k, int subset, int row, int column, double val);
double K_getRVLine(ok_kernel* k, int row, int col);
double K_getPeriodogramAt(ok_kernel* k, int row, int col);
int K_minimizeWithTimeout(ok_kernel* k, int timeout);
double K_getPhasedRVLine(ok_kernel* k, int planet, int row, int column);
double K_getPhasedDataForPlanet(ok_kernel* k, int planet, int row, int column);