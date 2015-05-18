/*
 *  kernel.h
 *  Systemic2
 *
 *
 */

#include "systemic.h"
#include "utils.h"

#define MV_VALUE 0
#define MV_MIN 1
#define MV_MAX 2
#define MV_STEP 3
#define MV_PARINDEX 4
#define MV_PARTYPE 5


#define MV_TYPE_ELEMENT 0
#define MV_TYPE_PAR 1


#define K_GETSET_H(name, Name, type) \
        void K_set##Name(ok_kernel* k, type value);\
        type K_get##Name(ok_kernel* k);

#define K_GETSET_C(name, Name, type) \
        void K_set##Name(ok_kernel* k, type value) { k-> name = value; K_validate(k); k->flags |= NEEDS_SETUP; } \
        type K_get##Name(ok_kernel* k) { return k-> name; }

#define K_GETSET_FREE_C(name, Name, type, freefcn) \
        void K_set##Name(ok_kernel* k, type name) { if (k-> name == name) return; freefcn(k->name); k-> name = name;  } \
        type K_get##Name(ok_kernel* k) { return k-> name; }


#define K_GETSETW_C(name, Name, type, where) \
        void K_set##Name(ok_kernel* k, type name) { where = name; } \
        type K_get##Name(ok_kernel* k) { return where; }

#define K_GET_H(name, Name, type) \
        type K_get##Name(ok_kernel* k);

#define K_GET_C(name, Name, type) \
        type K_get##Name(ok_kernel* k) { return k-> name; }


extern double ok_default_steps[];

// MEMORY MANAGEMENT
// Allocates a new kernel
ok_kernel* K_alloc();
// Frees the memory associated with the kernel
void K_free(ok_kernel* k);
// Returns a new copy of the current kernel (that you manage)
ok_kernel* K_clone(ok_kernel* k);
// Returns a new copy of the current kernel (that you manage), with data sharing options
ok_kernel* K_cloneFlags(ok_kernel* k, unsigned int shareFlags);

// RV MANAGEMENT
// add a new rv from an ASCII file; returns NULL on error
gsl_matrix* K_addDataFile(ok_kernel* k, const char* path, int type);
// add a new rv from a matrix
gsl_matrix* K_addDataTable(ok_kernel* k, gsl_matrix* rvtable, const char* name, int type);
void K_removeData(ok_kernel* k, int idx);
gsl_matrix* K_getData(ok_kernel* k, int idx);
void K_setData(ok_kernel* k, int idx, gsl_matrix* data);
const char* K_getDataName(ok_kernel* k, int idx);
bool K_addDataFromSystem(ok_kernel* k, const char* filename);
int K_getDataType(ok_kernel* k, int idx);

int K_getDataSize(ok_kernel* k, int idx);

// returns a compiled view of the data
double** K_compileData(ok_kernel* k);
// returns a matrix of the residuals
gsl_matrix* K_getCompiledDataMatrix(ok_kernel* k);
// returns a matrix of the data suitable for the periodogram
//gsl_matrix* K_combineData(ok_kernel* k, const int data_type);

// ELEMENTS MANAGEMENT
void K_addPlanet(ok_kernel* k, const double elements[]);
void K_removePlanet(ok_kernel* k, int idx);
void K_setElement(ok_kernel* k, int row, int col, double value);
double K_getElement(ok_kernel* k, int row, int col);
void K_setElements(ok_kernel* k, gsl_matrix* elements);
gsl_matrix* K_getElements(ok_kernel* k);
void K_setElementFlag(ok_kernel* k, int row, int col, int value);
int K_getElementFlag(ok_kernel* k, int row, int col);
void K_setElementStep(ok_kernel* k, int row, int col, double value);
double K_getElementStep(ok_kernel* k, int row, int col);
void K_setElementRange(ok_kernel* k, int row, int col, double min, double max);
void K_getElementRange(ok_kernel* k, int row, int col, double* min, double* max);
gsl_matrix* K_getAllElements(ok_kernel* k);
int K_getActivePars(ok_kernel* k);
int K_getActiveElements(ok_kernel* k);
int K_getNrPars(ok_kernel* k);

void K_setElementType(ok_kernel* k, int type);
int K_getElementType(ok_kernel* k);
double K_getMinValue(ok_kernel* k);

void K_setMinimizedValues(ok_kernel* k, double* values);
void K_getMinimizedValues(ok_kernel* k, double* values);

gsl_matrix* K_getXYZ(ok_kernel* k);


// STELLAR PARAMETERS
K_GETSET_H(Mstar, Mstar, double)
        
// SYSTEM PARAMETERS
K_GETSET_H(epoch, Epoch, double)
K_GET_H(chi2, Chi2, double)

double K_getChi2_nr(ok_kernel* k);
double K_getLoglik(ok_kernel* k);

K_GET_H(rms, Rms, double)
K_GET_H(jitter, Jitter, double)
K_GET_H(compiled, Compiled, double**)
        
K_GET_H(nrvs, Nrvs, unsigned int)
K_GET_H(ntts, Ntts, unsigned int)
K_GET_H(nsets, Nsets, unsigned int)
K_GET_H(chi2_rvs, Chi2_rvs, double)
K_GET_H(chi2_tts, Chi2_tts, double)        
K_GET_H(rms_tts, Rms_tts, double)        
        
K_GETSET_H(flags, Flags, unsigned int)
K_GETSET_H(minfunc, MinFunc, ok_callback)
K_GETSET_H(plSteps, ElementSteps, gsl_matrix*)
K_GETSET_H(parSteps, ParSteps, gsl_vector*)
K_GETSET_H(plFlags, ElementFlags, gsl_matrix_int*)
K_GETSET_H(parFlags, ParFlags, gsl_vector_int*)

K_GETSET_H(intMethod, IntMethod, int)
K_GETSET_H(progress, Progress, ok_progress)
K_GETSET_H(intOptions, IntOptions, ok_integrator_options*)
K_GETSET_H(model_function, CustomModelFunction, ok_model_function)

K_GETSET_H(intOptions->absacc, IntAbsAcc, double)
K_GETSET_H(intOptions->relacc, IntRelAcc, double)
K_GETSET_H(intOptions->dt, IntDt, double)
        
unsigned int K_getNplanets(ok_kernel* k);
unsigned int K_getNdata(ok_kernel* k);
void K_getRange(ok_kernel* k, double* from, double* to);

void K_perturb(ok_kernel* k);

// PARAMS MANAGEMENT
void K_setPars(ok_kernel* k, gsl_vector* pars);
gsl_vector* K_getPars(ok_kernel* k);
void K_setPar(ok_kernel* k, int idx, double val);
double K_getPar(ok_kernel* k, int idx);
void K_setParFlag(ok_kernel* k, int idx, int value);
int K_getParFlag(ok_kernel* k, int idx);
void K_setParStep(ok_kernel* k, int idx, double value);
double K_getParStep(ok_kernel* k, int idx);
void K_setParRange(ok_kernel* k, int idx, double min, double max);
void K_getParRange(ok_kernel* k, int idx, double* min, double* max);

void K_getMinimizedIndex(ok_kernel* k, int index, int* row, int* column);

// Load/Save
bool K_save(ok_kernel* k, FILE* fid);
ok_kernel* K_load(FILE* fid, int skip);
bool K_addDataFromSystem(ok_kernel* k, const char* filename);

// OPERATIONS
// Calculate chi^2
void K_calculate(ok_kernel* k);
int K_minimize(ok_kernel* k, int algo, int maxiter, double params[]);
int K_1dminimize(ok_kernel* k, int algo, int maxiter, int row, int column, double params[]);

ok_system** K_integrate(ok_kernel* k, gsl_vector* times, ok_system** bag, int* error);
ok_system** K_integrateRange(ok_kernel* k, double from, double to, unsigned int samples, ok_system** bag, int* error);
gsl_matrix* K_integrateStellarVelocity(ok_kernel* k, double from, double to, unsigned int samples, ok_system** bag, int* error);
ok_system** K_integrateProgress(ok_kernel* k, gsl_vector* times, ok_system** bag, int* error);

#ifndef PARSE
ok_kernel_minimizer_pars K_getMinimizedVariables(ok_kernel* k);
#endif

void K_setInfo(ok_kernel* k, const char* tag, const char* info);
char* K_getInfoTag(ok_kernel* k, int n);
char* K_getInfo(ok_kernel* k, const char* tag);
void K_clearInfo(ok_kernel* k);
bool K_infoExists(ok_kernel* k, const char * tag);

// DEBUG
void K_print(ok_kernel* k, FILE* f);
void K_save_old(ok_kernel* k, const char* stem);
void K_setSeed(ok_kernel* k, unsigned long int seed);

// BRIDGE UTILITIES
void* ok_bridge_kernel_buf(void* buf, int n, ok_kernel* k);

