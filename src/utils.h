//
//  utils.h
//  Systemic Console

//

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_rng.h"
#include "math.h"
#include "stdbool.h"
#include "assert.h"
#include "gsl/gsl_sort_int.h"

#ifndef UTILS_H
#define UTILS_H

#ifndef M_PI
        #define M_PI 3.14159265e+00
#endif

#define DEBUG 

#ifdef __APPLE__
#define ok_qsort_r qsort_r
#else
typedef int		 cmp_t(void *, const void *, const void *);
extern void ok_qsort_r(void *a, size_t n, size_t es, void *thunk, cmp_t *cmp);
#endif

#define MIN(a,b) (a < b ? a : b)
#define MAX(a,b) (a > b ? a : b)

#define RANGE(x, a, b) (MIN(MAX(x, a), b))

#define MGET(m, i, j) (m->data[(i)*m->tda+(j)])
#define MSET(m, i, j, x) m->data[(i)*m->tda+(j)] = x 
#define MSET_DEBUG(m, i, j, x) assert((i)*m->tda+(j)<m->size1*m->size2); m->data[(i)*m->tda+(j)] = x 
#define MIGET gsl_matrix_int_get
#define MISET gsl_matrix_int_set

#define MINC(m, i, j, c) MSET(m, i, j, MGET(m, i, j) + c)
#define MROWS(m) (m->size1)
#define MCOLS(m) (m->size2)

#define VGET(v, i) (v->data[(i)])
#define VSET(v, i, x) v->data[(i)] = (x)
#define VINC(v, i, x) VSET(v, i, VGET(v, i) + x)

#define VIGET gsl_vector_int_get
#define VISET gsl_vector_int_set

#define ROW gsl_matrix_row
#define COL gsl_matrix_col

#define SQR(x) ((x)*(x))
#define POW_2(x) ((x)*(x))
#define POW_3(x) ((x)*POW_2(x))
#define POW_4(x) ((x)*POW_3(x))
#define POW_5(x) ((x)*POW_4(x))
#define POW_6(x) ((x)*POW_5(x))
#define POW_7(x) ((x)*POW_6(x))


#define PRINTNUM(v) printf("" #v " = %e\n", v)

#define VLEN(X) (sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]))
#define DELTA(I, J) (I == J ? 1 : 0)

#define MSUN_TO_INT(M) ((M) * K2)
#define INT_TO_MSUN(M) ((M) / K2)

#define MJUP_TO_INT(M) ((M) * MJUP / MSUN * K2)
#define INT_TO_MJUP(M) ((M) / K2 * MSUN / MJUP)

#define MSUN_TO_MJUP(M) ((M) * MSUN / MJUP)
#define MJUP_TO_MSUN(M) ((M) * MJUP / MSUN)

#define GM_TO_INT(M) ((M) / MSUN * K2)
#define INT_TO_GM(M) ((M) / K2 * MSUN)

#define DAY_TO_YEAR(T) ((T) * DAY / YEAR)
#define YEAR_TO_DAY(T) ((T) / DAY * YEAR)

#define CM_TO_INT(L) ((L) / AU)
#define INT_TO_CM(L) ((L) * AU)

#define AUPDAY_TO_MPS(V) ((V) * AU / 100 / DAY)
#define MPS_TO_AUPDAY(V) ((V) * DAY / (AU / 100))

#define sqr(x) ((x)*(x))
#define TO_DEG(ANGLE) (DEGRANGE((ANGLE) / TWOPI * 360.))
#define TO_RAD(ANGLE) (RADRANGE((ANGLE) / 360. * TWOPI))

#define MATRIX_MEMCPY(dest, src) memcpy(dest->data, src->data, sizeof(double)*src->size1*src->size2)
#define VECTOR_MEMCPY(dest, src) memcpy(dest->data, src->data, sizeof(double)*src->size)
#define MATRIX_MEMCPY_TOARRAY(dest, src) memcpy(dest, src->data, sizeof(double)*src->size1*src->size2)
#define MATRIX_MEMCPY_FROMARRAY(dest, src) memcpy(dest->data, src, sizeof(double)*dest->size1*dest->size2)



typedef struct {
    void* data;
    int tag[5];
} ok_tagdata;


double DEGRANGE(double angle);
double RADRANGE(double angle);


void ok_sort_small_matrix(gsl_matrix* matrix, const int column);
void ok_sort_matrix(gsl_matrix* matrix, const int column);
void ok_rsort_matrix(gsl_matrix* matrix, const int column);

void ok_fprintf_matrix(gsl_matrix* matrix, FILE* file, const char* format);
void ok_fprintf_vector(gsl_vector* v, FILE* file, const char* format);
void ok_fprintf_matrix_int(gsl_matrix_int* matrix, FILE* file, const char* format);
void ok_fprintf_vector_int(gsl_vector_int* v, FILE* file, const char* format);
void ok_fprintf_buf(double** buf, FILE* file, const char* format, int rows, int columns);

void ok_fscanf_matrix(gsl_matrix* matrix, FILE* file);
void ok_fscanf_matrix_int(gsl_matrix_int* matrix, FILE* file);
void ok_fscanf_vector(gsl_vector* matrix, FILE* file);
void ok_fscanf_vector_int(gsl_vector_int* matrix, FILE* file);


bool ok_save_matrix(gsl_matrix* matrix, FILE* fid, const char* format);
bool ok_save_matrix_bin(gsl_matrix* matrix, FILE* fid);
bool ok_save_buf(double** matrix, FILE* fid, const char* format, int rows, int cols);
bool ok_save_buf_bin(double** matrix, FILE* fid, int rows, int cols);

gsl_matrix* ok_read_table(FILE* fid);

gsl_vector* ok_vector_resize(gsl_vector* v, int len);
gsl_matrix* ok_matrix_resize(gsl_matrix* v, int rows, int cols);
gsl_vector* ok_vector_resize_pad(gsl_vector* v, int len, double pad);
gsl_matrix* ok_matrix_resize_pad(gsl_matrix* v, int rows, int cols, double pad);

gsl_vector_int* ok_vector_int_resize(gsl_vector_int* v, int len);
gsl_matrix_int* ok_matrix_int_resize(gsl_matrix_int* v, int rows, int cols);

gsl_matrix* ok_matrix_copy(const gsl_matrix* src);
gsl_matrix* ok_matrix_copy_sub(const gsl_matrix* src, int row1, int nrows, int col1, int ncols);
gsl_matrix_int* ok_matrix_int_copy(const gsl_matrix_int* src);
gsl_vector* ok_vector_copy(const gsl_vector* src);
gsl_vector_int* ok_vector_int_copy(const gsl_vector_int* src);
gsl_matrix* ok_buf_to_matrix(double** buf, int rows, int cols);
void ok_buf_col(double** buf, double* vector, int col, int nrows);

void ok_matrix_column_range(gsl_matrix* m, int col, double* min, double* max);

gsl_matrix* ok_matrix_remove_row(gsl_matrix* m, int row);
gsl_matrix* ok_matrix_remove_column(gsl_matrix* m, int col);

gsl_matrix_int* ok_matrix_int_remove_row(gsl_matrix_int* m, int row);
gsl_matrix_int* ok_matrix_int_remove_column(gsl_matrix_int* m, int col);
gsl_vector* ok_vector_remove(gsl_vector* m, int idx);
gsl_vector_int* ok_vector_int_remove(gsl_vector_int* m, int idx);

void ok_matrix_fill(gsl_matrix* src, gsl_matrix* dest);

void ok_bootstrap_matrix(const gsl_matrix* matrix, gsl_matrix* out, const int sortcol, gsl_rng* r);
void ok_bootstrap_matrix_mean(const gsl_matrix* matrix, int timecol, int valcol, gsl_matrix* out, gsl_rng* r);

gsl_matrix* ok_matrix_filter(gsl_matrix* matrix, const int column, const double filter);
gsl_matrix* ok_matrix_buf_filter(double** matrix, const int rows, const int columns, const int column, const double filter);

int ok_bsearch(double* v, double val, int len);

double ok_average_angle(const double* v, const int length, const bool isRadians);
double ok_median_angle(const double* v, const int length, const bool isRadians);
double ok_stddev_angle(const double* v, const int length, const bool isRadians);
double ok_mad_angle(double* v, const int length, const double med, const bool isRadians);
double ok_mad(double* v, const int length, const double med);

char* ok_str_copy(const char* src);
char* ok_str_cat(const char* a1, const char* a2);

void ok_avevar(const double* v, int len, double* ave, double* var);

gsl_matrix* ok_ptr_to_matrix(double* v, unsigned int rows, unsigned int cols);
gsl_vector* ok_ptr_to_vector(double* v, unsigned int len);
gsl_matrix_int* ok_iptr_to_imatrix(int* v, unsigned int rows, unsigned int cols);
gsl_vector_int* ok_iptr_to_ivector(int* v, unsigned int len);

void ok_block_to_ptr(void* vv, double* out);

void ok_buf_to_ptr(double** v,  unsigned int rows, unsigned int cols, double* out);

void ok_buf_add_to_col(double** buf, double* col_vector, int col, int nrows);

unsigned int ok_vector_len(void* v);
unsigned int ok_matrix_rows(void* v);
unsigned int ok_matrix_cols(void* v);

gsl_block* ok_vector_block(void* v);
gsl_block* ok_matrix_block(void* v);

gsl_matrix* ok_resample_curve(gsl_matrix* curve, const int xcol, const int ycol, const double peaks_frac, const int target_points,
    const int target_tolerance, double* start_tolerance, const int max_steps, const bool log_x);

bool ok_file_readable(char* fn);

typedef struct {
    int length;
    int max_length;
    int* v;
} ok_rivector;

ok_rivector* ok_rivector_alloc(const int maxlength);
#define ok_rivector_data(va) (va->v)
#define ok_rivector_sizeof(va) (sizeof(int))
#define ok_rivector_push(va, i) do { assert(va->length < va->max_length); va->v[va->length++] = i; } while (0)
#define ok_rivector_append(vdest, vsource) do { for (int i = 0; i < vsource->length; i++) ok_rivector_push(vdest, vsource->v[i]); } while (0)
#define ok_rivector_pop(va) (va->v[--va->length]);
#define ok_rivector_first(va) (va->v[0])
#define ok_rivector_last(va) (va->v[v->length-1])
#define ok_rivector_length(va) (va->length)
#define ok_rivector_reset(va) do { va->length = 0; } while (0)
#define ok_rivector_reset_to(va, len) do { va->length = len; } while (0)
#define ok_rivector_free(va) do { free(va->v); free(va); } while (0)
#define ok_rivector_sort(va) do { gsl_sort_int(va->v, 1, va->length); } while (0)
#define ok_rivector_foreach(va, val) for (int __ri_idx ## val = 0, val = va->v[0]; __ri_idx ## val < va->length; __ri_idx ## val++, val=va->v[__ri_idx ## val])
#define ok_rivector_foreach_i(va, val, idx) for (int idx = 0, val = va->v[0]; idx < va->length; idx++, val=va->v[idx])
#endif