/*
 *  types.h
 *  Systemic2
 *
 *
 */

#include "stdlib.h"
#include "stdbool.h"
#include "stdio.h"
#include "string.h"
#include "float.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

#ifndef types_h
#define types_h

#define INVALID_NUMBER (NAN)
#define IS_INVALID(x) (isnan(x))

#define SYSTEMIC_VERSION 2.1600

#define MAX_LINE 8192

#define T_RV 0
#define T_PHOTO 1
#define T_TIMING 2
#define T_DUMMY 99

// Number of columns for data matrices. Columns 0-5 can be used to store data,
// while columns 6-10 are used internally.
#define DATA_SIZE 11

#define DATA_SETS_SIZE 10
// Columns of data matrices
// Time (usually JD)
#define T_TIME 0
// Value (e.g. RV, flux, transit timing, etc.)
#define T_VAL 1
// Error
#define T_ERR 2
// Timing datasets: flag to associate a transit to a specific planet number. Set to 1 (first planet) by default.
#define T_TDS_PLANET 3
// Timing type: primary or secondary
#define T_TDS_FLAG 4
// Type of data
#define T_FLAG 6
// Shifted value (T_VAL shifted by the RV offset, otherwise same as T_VAL)
#define T_SVAL 7
// Computed value
#define T_PRED 8
// Dataset number
#define T_SET 9
// Reserved for GUI to do its dirty stuff
#define T_GUI 10
#define T_SCRATCH 10

#define X 1
#define Y 2
#define Z 3
#define VX 4
#define VY 5
#define VZ 6

#define STAT_MEAN 0
#define STAT_MEDIAN -1
#define STAT_STDDEV -2
#define STAT_MAD -3
#define STAT_QUART_25 25
#define STAT_QUART_50 50
#define STAT_QUART_75 75

#define ENSURE_SIZE(matrix, rows, columns) assert((matrix->size1 == rows) && (matrix->size2 == columns))

#define KEPLER 0
#define RK45 1
#define RK89 2
#define ADAMS 3
#define BULIRSCHSTOER 4
#define SWIFTRMVS 5

//#define K2 2.959122082855911e-4

#define AU 1.4959787e13
#define MSUN 1.98892e33
#define MJUP 1.8986e30
#define RJUP 7.1e9
#define GGRAV 6.67384e-8

#define DAY 8.64e4
#define TWOPI 6.2831853072e+00
#define SQRT_TWOPI 2.5066282746e+00
#define K2  ((GGRAV * MSUN * DAY * DAY) / (AU*AU*AU))
#define YEAR 31556926.

#define PER 0
#define MASS 1
#define MA 2
#define ECC 3
#define LOP 4
#define INC 5
#define NODE 6
#define RADIUS 7
#define ORD 8
// Extra elements that are calculated on-the-fly
#define SMA 13
#define SEMIAMP 14
#define TPERI 15
#define TRUEANOMALY 16
#define MEANLONGITUDE 17

#define ELEMENTS_SIZE 13
#define ALL_ELEMENTS_SIZE 20

#define P_DATA1 0
#define P_DATA2 1
#define P_DATA3 2
#define P_DATA4 3
#define P_DATA5 4
#define P_DATA6 5
#define P_DATA7 6
#define P_DATA8 7
#define P_DATA9 8
#define P_DATA10 9

#define P_DATA_NOISE1 10
#define P_DATA_NOISE2 11
#define P_DATA_NOISE3 12
#define P_DATA_NOISE4 13
#define P_DATA_NOISE5 14
#define P_DATA_NOISE6 15
#define P_DATA_NOISE7 16
#define P_DATA_NOISE8 17
#define P_DATA_NOISE9 18
#define P_DATA_NOISE10 19

#define P_RV_TREND 20
#define P_RV_TREND_QUADRATIC 21

#define PARAMS_SIZE 100

extern char * ok_orb_labels[ELEMENTS_SIZE];
extern char * ok_all_orb_labels[ALL_ELEMENTS_SIZE];

#define OPT_EPS 0
#define OPT_ECC_LAST 1

#define OPT_MCMC_SKIP_STEPS 0
#define OPT_MCMC_SKIP 1
#define OPT_MCMC_RSTOP_SINGLE 2
#define OPT_MCMC_RETURN_ALL 3
#define OPT_MCMC_BETA 4
#define OPT_MCMC_NSTOP 5
#define OPT_MCMC_TEMPFAC 6
#define OPT_MCMC_VERBOSE_DIAGS 7
#define OPT_MCMC_ACCRATIO 8
#define OPT_MCMC_NMIN 9

#define OPT_LM_MINCHI_PAR 10
#define OPT_LM_HIGH_DF 11
#define OPT_LM_MAX_ITER_AT_SCALE 12
#define OPT_LM_INITIAL_SCALE 13

#define OPT_SA_T0 20
#define OPT_SA_ALPHA 21
#define OPT_SA_CHAINS 22
#define OPT_SA_AUTO_STEPS 23

#define OPT_DE_CR 30
#define OPT_DE_NP_FAC 31
#define OPT_DE_F_MIN 32
#define OPT_DE_F_MAX 33
#define OPT_DE_USE_STEPS 34



#define PROGRESS_CONTINUE 0
#define PROGRESS_STOP 1

#define STATUS_SUCCESS GSL_SUCCESS

#define TDS_PRIMARY 1
#define TDS_SECONDARY 2

#define DONE -1

#define ACTIVE (1 << 1)
#define MINIMIZE (1 << 2)

#define NEEDS_COMPILE (1 << 1)
#define NEEDS_SETUP (1 << 2)
#define BOOTSTRAP_DATA (1 << 3)
#define RETAIN (1 << 4)
#define SHARE_FLAGS (1 << 5)
#define SHARE_STEPS (1 << 6)
#define SHARE_DATA (1 << 7)
#define SHARE_RANGES (1 << 11)
#define FREED (1 << 8)
#define DUPLICATE (1 << 9)
#define GUI_RESERVED_1 (1 << 12)
#define GUI_RESERVED_2 (1 << 13)
#define GUI_RESERVED_3 (1 << 14)
#define GUI_RESERVED_4 (1 << 15)

#define ASTROCENTRIC 0
#define JACOBI (1 << 10)

#define SIMPLEX 0
#define LM 1
#define DIFFEVOL 2
#define SA 3

#define INTEGRATION_SUCCESS 0
#define INTEGRATION_FAILURE_SMALL_TIMESTEP (1 << 11)
#define INTEGRATION_FAILURE_INCREASE_TOLERANCE (1 << 12)
#define INTEGRATION_FAILURE_CLOSE_ENCOUNTER (1 << 13)
#define INTEGRATION_FAILURE_CLOSE_ENCOUNTER_STAR (1 << 14)
#define INTEGRATION_FAILURE_STOPPED (1 << 15)
#define INTEGRATION_FAILURE_SWIFT (1 << 16)
typedef struct ok_system {
    /// The number of planets
    int nplanets;
    /// Initial epoch, to which the elements
    /// matrix below refers
    double epoch;
    /// Initial elements at epoch, in astrocentric 
    /// coordinates (by default)
    /// Units: period in days, mass in Mjup,
    /// MA, inc, node, lop in degrees.
    gsl_matrix* elements;

    /// The time to which the coordinates below 
    /// refer.
    double time;
    /// Cartesian coordinates. Internal units 
    /// are used.
    gsl_matrix* xyz;
    /// Current orbital elements. Internal units
    /// are used.
    gsl_matrix* orbits;
    /// flags
    int flag;
} ok_system;

typedef struct ok_kernel ok_kernel;
typedef double(*ok_callback)(ok_kernel*);
typedef void(*ok_callback2)(ok_kernel*, double*);
typedef void(*ok_icallback)(int);

#define INTEGRATORS_SIZE 4

typedef int(*ok_progress)(int current, int max, void* state, const char* function);


typedef struct ok_integrator_options {
    // Absolute and relative accuracy (used by the RK integrators and SWIFT_BS)
    double abs_acc;
    double rel_acc;
    // An accuracy factor used to determine the time step in HERMITE and HERMITE_ACC
    double acc_par;
    // Fixed time step for the other SWIFT integrators
    double dt;
    // Accuracy parameter for the transits
    double eps_tr;
    
    int iterations;
    
    bool calc_elements;
    /// A pointer to a force function that takes in time,
    /// the current state as an array, the output force
    /// array and a pointer to the current system state.
    int (* force) (double, const double[], double[], void*);
    int (* jac) (double, const double[], double*, double[], void*);
    int (* force_jerk) (double, const double[], double[], double[], void*);
    
    // A buffer, used to hold memory to be reused between
    // integrator calls.
    gsl_vector* buffer;
    gsl_vector_int* ibuffer;
    
    ok_progress progress;
} ok_integrator_options;


typedef void(*ok_model_function)(ok_kernel*, double** data, int ndata);
typedef ok_system**(*ok_integrator)(ok_system*, const gsl_vector*, ok_integrator_options*, ok_system** bag);
typedef int(*ok_minimizer)(ok_kernel*, int, double[]);


struct ok_kernel {
    // initial conditions
    ok_system* system;
    // reduced chi^2
    double chi2;
    // chi^2
    double chi2_rvs;
    // chi^2 tts
    double chi2_tts;
    // chi^2 for data that is not flagged T_RV or T_TIMING
    double chi2_other;
    // rms of residuals
    double rms;
    // rms of tts
    double rms_tts;
    // jitter
    double jitter;
    // number of radial velocities
    int nsets;
    // list of radial velocities
    gsl_matrix** datasets;
    
    // compiled data (merged and sorted); each row is a pointer to a row
    // of data (RV, photometry, transit timing, etc.)
    double** compiled;
    // velocity offsets
    gsl_vector* params;
    // times
    gsl_vector* times;
    // result of last integration
    ok_system** integration;
    // number of integration samples
    int integrationSamples;

    // minimization flags
    gsl_matrix_int* plFlags;
    gsl_vector_int* parFlags;
    gsl_matrix* plSteps;
    gsl_vector* parSteps;
    gsl_matrix** plRanges;
    gsl_vector** parRanges;

    // stellar parameters
    double Mstar;

    // total number of data points
    int ndata;
    int nrvs;
    int ntts;
    int nrpars;
    
    // integration
    int intMethod;
    ok_integrator_options* intOptions;

    // function to minimize; usually, chi^2
    ok_callback minfunc;

    // flags
    unsigned int flags;

    // private random number generator
    gsl_rng* rng;

    // tagging data
    void* tag;
    double tagValue;
    
    // progress callback
    ok_progress progress;
    
    // dataset names
    char datanames[DATA_SETS_SIZE][MAX_LINE];
    
    // custom model function
    ok_model_function model_function;
    int last_error;
    
    
};

typedef struct ok_list_item {
    gsl_matrix* elements;
    gsl_vector* params;
    double merit;
    double merit_pr;
    double merit_li;
    int tag;
} ok_list_item;
    

typedef struct ok_list {
    ok_kernel* prototype;
    ok_list_item** kernels;
    int size;
    void* diags;
    int type;
} ok_list;

typedef struct ok_kernel_minimizer_pars {
    double** pars;
    double* steps;
    double* min;
    double* max;
    int* type;
    int npars;
} ok_kernel_minimizer_pars;

#define FREE_MINIMIZER_PARS(p) do { free(p.pars); free(p.steps); free(p.min); free(p.max); free(p.type); } while (0);

#endif
