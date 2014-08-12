//
//  utils.c
//  Systemic Console
//
//

#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_randist.h>
#include "systemic.h"
#include "utils.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "assert.h"
#include "kernel.h"
#include "time.h"

double DEGRANGE(double angle) {
    return (fmod((angle < 0. ? angle + (floor(-angle/360.) + 1.) * 360. : angle), 360.));
}

double RADRANGE(double angle) {
    return (fmod((angle < 0. ? angle + (floor(-angle/TWOPI) + 1.) * TWOPI : angle), TWOPI));
}

double ok_average_angle(const double* v, const int length, const bool isRadians) {
    
    double cos_avg = cos((isRadians ? v[0] : TO_RAD(v[0])));
    double sin_avg = sin((isRadians ? v[0] : TO_RAD(v[0])));

    for (int i = 1; i < length; i++) {
        cos_avg += cos((isRadians ? v[i] : TO_RAD(v[i])));
        sin_avg += sin((isRadians ? v[i] : TO_RAD(v[i])));
    }

    double avg = atan2(sin_avg / (double) length, cos_avg / (double) length);

    if (avg < 0) {
        avg += 2 * M_PI;
    }

    if (!isRadians) {
        avg = TO_DEG(avg);
    }

    return avg;
}


double ok_median_angle(const double* v, int length, bool isRadians) {
    
    gsl_vector* cos_avg = gsl_vector_alloc(length);
    gsl_vector* sin_avg = gsl_vector_alloc(length);
    

    for (int i = 0; i < length; i++) {
        cos_avg->data[i] = cos((isRadians ? v[i] : TO_RAD(v[i])));
        sin_avg->data[i] = sin((isRadians ? v[i] : TO_RAD(v[i])));
    }

    gsl_sort_vector(cos_avg);
    gsl_sort_vector(sin_avg);
    
    double avg = atan2(gsl_stats_median_from_sorted_data(sin_avg->data, 1, length), 
            gsl_stats_median_from_sorted_data(cos_avg->data, 1, length));

    if (avg < 0) {
        avg += 2 * M_PI;
    }

    if (!isRadians) {
        avg = TO_DEG(avg);
    }

    gsl_vector_free(cos_avg);
    gsl_vector_free(sin_avg);
    return avg;
}

double ok_stddev_angle(const double* v, const int length, const bool isRadians)  {
    if (length <= 1) {
        return 0;
    }

    double avg = ok_average_angle(v, length, isRadians);
    double avg_r = (isRadians ? avg : TO_RAD(avg));
        
    double dev = 0;

    for (int i = 0; i < length; i++) {
        double val = (isRadians ? v[i] : TO_RAD(v[i]));
        val = fmod(val, 2 * M_PI);
        if (val < 0)
            val += 2 * M_PI;

        double diff = MIN(RADRANGE(fabs(val - avg_r)), RADRANGE(fabs(val - avg_r + TWOPI)));
        diff = MIN(diff, RADRANGE(fabs(val - avg_r - TWOPI)));
        dev += diff * diff;
     
    }

    dev = sqrt(1. / (double) (length) * dev);


    if (! isRadians) {
        dev = TO_DEG(dev);
    }


    return dev;
}


double ok_mad_angle(double* v, const int length, const double med, const bool isRadians)  {
    if (length <= 1) {
        return 0;
    }

    
    double med_r = (isRadians ? med : TO_RAD(med));

    for (int i = 0; i < length; i++) {
        double val = (isRadians ? v[i] : TO_RAD(v[i]));
        val = fmod(val, 2 * M_PI);
        if (val < 0)
            val += 2 * M_PI;

        double diff = MIN(RADRANGE(fabs(val - med_r)), RADRANGE(fabs(val - med_r + TWOPI)));
        diff = MIN(diff, RADRANGE(fabs(val - med_r - TWOPI)));
        v[i] = fabs(diff);
    }
    
    gsl_sort(v, 1, length);

    double mad = gsl_stats_median_from_sorted_data(v, 1, length);
    return (isRadians ? mad : TO_DEG(mad));
}

double ok_mad(double* v, const int length, const double med) {
    
    for (int i = 0; i < length; i++)
        v[i] = fabs(v[i] - med);
    
    gsl_sort(v, 1, length);
    return gsl_stats_median_from_sorted_data(v, 1, length);
    
}


/*  */
void ok_sort_small_matrix(gsl_matrix* matrix, const int column) {
    const int nrows = matrix->size1;
    int swaps = 1;
    while (swaps > 0) {
        swaps = 0;
        for (int i = 1; i < nrows; i++) {
            if (MGET(matrix, i, column) < MGET(matrix, i-1, column)) {
                gsl_matrix_swap_rows(matrix, i, i-1);
                swaps++;
            }
        }
    }
}

void ok_bootstrap_matrix(const gsl_matrix* matrix, gsl_matrix* out, const int sortcol, gsl_rng* r) {
    if (!out)
        out = gsl_matrix_alloc(matrix->size1, matrix->size2);
    
    for (int i = 0; i < matrix->size1; i++) {
        int k = gsl_rng_uniform_int(r, matrix->size1);
        
        for (int j = 0; j < matrix->size2; j++) {
            MSET(out, i, j, MGET(matrix, k, j));
        }
    }
    
    if (sortcol >= 0)
        ok_sort_matrix(out, sortcol);
}

void ok_bootstrap_matrix_mean(const gsl_matrix* matrix, int timecol, int valcol, gsl_matrix* out, gsl_rng* r) {
    if (!out)
        out = gsl_matrix_alloc(matrix->size1, matrix->size2);
    
    double mean = 0.;
    //for (int i = 0; i < matrix->size1; i++)
      //  mean += MGET(matrix, i, valcol);
    //mean /= (double) matrix->size1;
    int a[matrix->size1];
    for (int i = 0; i < matrix->size1; i++)
        a[i] = i;
    gsl_ran_shuffle(r, a, matrix->size1, sizeof(int));
    
    for (int i = 0; i < matrix->size1; i++) {
        int k = a[i];
        
        for (int j = 0; j < matrix->size2; j++) {
            if (j == timecol) 
                MSET(out, i, j, MGET(matrix, i, j));
            else if (j == valcol) {
                MSET(out, i, j, MGET(matrix, k, j) - mean);
            } else 
                MSET(out, i, j, MGET(matrix, k, j));
        }
    }
    
}



void ok_fprintf_matrix(gsl_matrix* matrix, FILE* file, const char* format) {
    if (matrix == NULL) return;
    file = (file == NULL ? stdout : file);
    for (int i = 0; i < matrix->size1; i++) {
        for (int j = 0; j < matrix->size2; j++)
            fprintf(file, format, MGET(matrix, i, j));
        fprintf(file, "\n");
    }
}

void ok_fprintf_matrix_int(gsl_matrix_int* matrix, FILE* file, const char* format) {
    if (matrix == NULL) return;
    file = (file == NULL ? stdout : file);
    for (int i = 0; i < matrix->size1; i++) {
        for (int j = 0; j < matrix->size2; j++)
            fprintf(file, format, MIGET(matrix, i, j));
        fprintf(file, "\n");
    }
}

bool ok_save_matrix(gsl_matrix* matrix, FILE* fid, const char* format) {
    fid = (fid == NULL ? stdout : fid);
    ok_fprintf_matrix(matrix, fid, format);
    return true;
}

bool ok_save_matrix_bin(gsl_matrix* matrix, FILE* fid) {
    fid = (fid == NULL ? stdout : fid);
    gsl_matrix_fwrite(fid, matrix);
    return true;
}


bool ok_save_buf(double** matrix, FILE* fid, const char* format, int rows, int cols) {
    fid = (fid == NULL ? stdout : fid);
    ok_fprintf_buf(matrix, fid, format, rows, cols);
    return true;
}

bool ok_save_buf_bin(double** matrix, FILE* fid, int rows, int cols) {
    fid = (fid == NULL ? stdout : fid);
    for (int i = 0; i < rows; i++)
        fwrite(matrix[i], sizeof(double), cols, fid);
    return true;
}

void ok_fprintf_vector(gsl_vector* v, FILE* file, const char* format) {
    if (v == NULL) return;
    file = (file == NULL ? stdout : file);
    for (int i = 0; i < v->size; i++) 
       fprintf(file, format, VGET(v, i));
    fprintf(file, "\n");
}

void ok_fprintf_vector_int(gsl_vector_int* v, FILE* file, const char* format) {
    if (v == NULL) return;
    file = (file == NULL ? stdout : file);
    for (int i = 0; i < v->size; i++) 
       fprintf(file, format, VIGET(v, i));
    fprintf(file, "\n");
}

void ok_fprintf_buf(double** buf, FILE* file, const char* format, int rows, int columns) {
    if (buf == NULL) return;
    file = (file == NULL ? stdout : file);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++)
            fprintf(file, format, buf[i][j]);
        fprintf(file, "\n");
    }
}

gsl_vector* ok_vector_resize(gsl_vector* v, int len) {
    gsl_vector* nv = gsl_vector_calloc(len);
    if (v != NULL) {
        for (int i = 0; i < MIN(len, v->size); i++)
                VSET(nv, i, VGET(v, i));
        gsl_vector_free(v);
    }
    return nv;
}

gsl_vector* ok_vector_resize_pad(gsl_vector* v, int len, double pad) {
    gsl_vector* nv = gsl_vector_calloc(len);
    gsl_vector_set_all(v, pad);
    if (v != NULL) {
        for (int i = 0; i < MIN(len, v->size); i++)
                VSET(nv, i, VGET(v, i));
        gsl_vector_free(v);
    }
    return nv;
}

gsl_vector_int* ok_vector_int_resize(gsl_vector_int* v, int len) {
    gsl_vector_int* nv = gsl_vector_int_calloc(len);
    if (v != NULL) {
        for (int i = 0; i < MIN(len, v->size); i++)
                VISET(nv, i, VIGET(v, i));
        gsl_vector_int_free(v);
    }
    return nv;
}

gsl_matrix* ok_matrix_resize(gsl_matrix* v, int rows, int cols) {
    gsl_matrix* m = gsl_matrix_calloc(rows, cols);
    if (v != NULL) {
        for (int i = 0; i < MIN(rows, MROWS(v)); i++)
            for (int j = 0; j < MIN(cols, MCOLS(v)); j++)
                MSET(m, i, j, MGET(v, i, j));
        gsl_matrix_free(v);
    }
    return m;
}

gsl_matrix* ok_matrix_resize_pad(gsl_matrix* v, int rows, int cols, double pad) {
    gsl_matrix* m = gsl_matrix_calloc(rows, cols);
    gsl_matrix_set_all(m, pad);
    if (v != NULL) {
        for (int i = 0; i < MIN(rows, MROWS(v)); i++)
            for (int j = 0; j < MIN(cols, MCOLS(v)); j++)
                MSET(m, i, j, MGET(v, i, j));
        gsl_matrix_free(v);
    }
    return m;
}

gsl_matrix_int* ok_matrix_int_resize(gsl_matrix_int* v, int rows, int cols) {
    gsl_matrix_int* m = gsl_matrix_int_calloc(rows, cols);
    if (v != NULL) {
        for (int i = 0; i < MIN(rows, MROWS(v)); i++)
            for (int j = 0; j < MIN(cols, MCOLS(v)); j++)
                MISET(m, i, j, MIGET(v, i, j));
        gsl_matrix_int_free(v);
    }
    return m;
}

gsl_matrix* ok_read_table(FILE* fid) {
	const int MAXLEN = 4000;
	char line[MAXLEN];
    
    double JD = 0.0;
    int tds_plFlag = 1;
    
	fgets(line, MAXLEN, fid);
	while (line[0] == '#') {
        char tag[100] = {0};
        sscanf(&line[1], "%s =", tag);
        
        if (strcmp(tag, "JD") == 0) {
            sscanf(&line[1], "%*s = %le", &JD);
        } else if (strcmp(tag, "Planet") == 0) {
            sscanf(&line[1], "%*s = %d", &tds_plFlag);
        }
        
		fgets(line, MAXLEN, fid);
    }
		
	int nrows = 1;
        
	while (fgets(line, MAXLEN, fid) != NULL) {
		if (line[0] != '#')
			nrows++;
	}
	
	fseek(fid, 0, SEEK_SET); 
	
	gsl_matrix* ret = gsl_matrix_calloc(nrows, DATA_SIZE);
        double v[6];
	int nr = 0;
	while (fgets(line, MAXLEN, fid)!=NULL) {
		if (line[0] == '#')
			continue;
        for (int j = 0; j < 6; j++)
            v[j] = 0.;

        sscanf(line, "%le %le %le %le %le %le",
                &(v[0]), &(v[1]), &(v[2]),
                &(v[3]), &(v[4]), &(v[5]));
        for (int j = 0; j < 6; j++)
            MSET(ret, nr, j, v[j]);

		nr++;
	}
    
    for (int i = 0; i < nrows; i++) {
        MINC(ret, i, T_TIME, JD);
        if (tds_plFlag > 0)
            MSET(ret, i, T_TDS_PLANET, tds_plFlag);
    }
	
	return ret;
}

gsl_matrix* ok_matrix_copy(const gsl_matrix* src) {
    if (src == NULL) return NULL;
    gsl_matrix* m = gsl_matrix_alloc(MROWS(src), MCOLS(src));
    MATRIX_MEMCPY(m, src);
    return m;
}

gsl_matrix* ok_matrix_copy_sub(const gsl_matrix* src, int row1, int nrows, int col1, int ncols) {
    if (src == NULL) return NULL;
    gsl_matrix* m = gsl_matrix_alloc(nrows-row1, ncols-col1);
    
    for (int i = row1; i < nrows; i++)
        for (int j = col1; j < ncols; j++)
            MSET(m, i-row1, j-col1, MGET(src, i, j));
    
    return m;
}



gsl_matrix_int* ok_matrix_int_copy(const gsl_matrix_int* src) {
    if (src == NULL) return NULL;
    gsl_matrix_int* m = gsl_matrix_int_alloc(MROWS(src), MCOLS(src));
    gsl_matrix_int_memcpy(m, src);
    return m;
}
gsl_vector* ok_vector_copy(const gsl_vector* src) {
    if (src == NULL) return NULL;
    gsl_vector* v = gsl_vector_alloc(src->size);
    gsl_vector_memcpy(v, src);
    return v;
}
gsl_vector_int* ok_vector_int_copy(const gsl_vector_int* src) {
    if (src == NULL) return NULL;
    gsl_vector_int* v = gsl_vector_int_alloc(src->size);
    gsl_vector_int_memcpy(v, src);
    return v;
}

char* ok_str_copy(const char* src) {
    if (src == NULL) return NULL;
    char* d = (char*) malloc(strlen(src) * sizeof(char*));
    strcpy(d, src);
    return d;
}

char* ok_str_cat(const char* a1, const char* a2) {
    if (a1 == NULL) return ok_str_copy(a2);
    else if (a2 == NULL) return ok_str_copy(a1);
    
    char* a3 = (char*) malloc((strlen(a1) + strlen(a2) + 1) * sizeof(char));
    strcpy(a3, a1);
    strcat(a3, a2);
    return a3;
}

int ok_bsearch(double* v, double val, int len) {
    
    int min = 0;
    int max = len;
    
    while (max >= min) {
        int mid = (int)(0.5 * (min + max));
        if (v[mid] < val)
            min = mid + 1;
        else if (v[mid] > val)
            max = mid - 1;
        else if (v[mid] == val)
            return mid;
    }
    
    return min;
}

void ok_avevar(const double* v, int len, double* ave, double* var) {
    *ave = 0.;
    for (int i = 0; i < len; i++)
        *ave += v[i];
    *ave /= (double)len;
    *var = 0.;
    double sum2 = 0.;
    
    for (int i = 0; i < len; i++) {
        *var += SQR(v[i] - *ave);
        sum2 += v[i]-*ave;
    }
    
    *var = (*var - sum2/(double) len) / (double)(len-1);
}

int _ok_compare(void* col, const void* row1, const void* row2) {
    int c = *((int*)col);
    double* d1 = (double*) row1;
    double* d2 = (double*) row2;
    
    if (d1[c] < d2[c])
        return -1;
    else if (d1[c] > d2[c]) 
        return 1;
    else
        return 0;
}

void ok_sort_matrix(gsl_matrix* matrix, const int column) {
    ok_qsort_r(matrix->data, matrix->size1, matrix->size2 * sizeof(double), (void*)&column, _ok_compare);
}

static int _ok_rcompare(void* col, const void* row1, const void* row2) {
    int c = *((int*)col);
    double* d1 = (double*) row1;
    double* d2 = (double*) row2;
    
    if (d1[c] < d2[c])
        return 1;
    else if (d1[c] > d2[c]) 
        return -1;
    else
        return 0;
}

void ok_rsort_matrix(gsl_matrix* matrix, const int column) {
    ok_qsort_r(matrix->data, matrix->size1, matrix->size2 * sizeof(double*), (void*)&column, _ok_compare);
}



gsl_matrix* ok_matrix_filter(gsl_matrix* matrix, const int column, const double filter) {
    int rows = 0;
    for (int i = 0; i < matrix->size1; i++)
        if (fabs(MGET(matrix, i, column) - filter) < 1e-10)
            rows++;
    gsl_matrix* ret = gsl_matrix_alloc(rows, matrix->size2);
    
    rows = 0;
    for (int i = 0; i < matrix->size1; i++)
        if (fabs(MGET(matrix, i, column) - filter) < 1e-10) {
            rows++;
            for (int j = 0; j < matrix->size2; j++)
                MSET(ret, rows, j, MGET(matrix, i, j));
        }
    
    return ret;
    
}

gsl_matrix* ok_matrix_buf_filter(double** matrix, const int rows, const int columns, const int column, const double filter) {
    int nrows = 0;
    for (int i = 0; i < rows; i++) 
        if (fabs(matrix[i][column] - filter) < 1e-10) 
            nrows++;
 
    gsl_matrix* ret = gsl_matrix_alloc(nrows, columns);
    
    nrows = 0;
    for (int i = 0; i < rows; i++) {
        if (fabs(matrix[i][column] - filter) < 1e-10) {
            for (int j = 0; j < columns; j++)
                MSET(ret, nrows, j, matrix[i][j]);
            nrows++;
        }
    }
    return ret;
    
}

// These functions are defined to avoid 
// emscripten's scanf bug (does not recognize "nan" as 
// a valid number)
void ok_fscanf_double(FILE* file, double* out) {
#ifndef JAVASCRIPT
    fscanf(file, "%le", out);
#else
    char token[2000];
    fscanf(file, "%s", token);
    if (strcmp(token, "nan") == 0)
        *out = INVALID_NUMBER;
    else
        *out = atof(token);
#endif
}

void ok_fscanf_int(FILE* file, int* out) {
#ifndef JAVASCRIPT
    fscanf(file, "%d", out);
#else
    char token[2000];
    fscanf(file, "%s", token);
    if (strcmp(token, "nan") == 0)
        *out = INVALID_NUMBER;
    else
        *out = atoi(token);
#endif
}

void ok_fscanf_matrix(gsl_matrix* matrix, FILE* file) {
    
    double* data = matrix->data;
    
    for (int i = 0; i < matrix->size1; i++)
        for (int j = 0; j < matrix->size2; j++) {
            ok_fscanf_double(file, &(data[i*matrix->size2 + j]));
        }
}

void ok_fscanf_matrix_int(gsl_matrix_int* matrix, FILE* file) {
    int* data = matrix->data;
    for (int i = 0; i < matrix->size1; i++)
        for (int j = 0; j < matrix->size2; j++)
            ok_fscanf_int(file, &(data[i*matrix->size2 + j]));
}
void ok_fscanf_vector(gsl_vector* vector, FILE* file) {
    double* data = vector->data;
    for (int i = 0; i < vector->size; i++)
        ok_fscanf_double(file, &(data[i]));
}
void ok_fscanf_vector_int(gsl_vector_int* vector, FILE* file) {
    int* data = vector->data;
    for (int i = 0; i < vector->size; i++)
        ok_fscanf_int(file, &(data[i]));
}


gsl_matrix* ok_buf_to_matrix(double** buf, int rows, int cols) {
    gsl_matrix* m = gsl_matrix_alloc(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            MSET(m, i, j, buf[i][j]);
    return m;
}

gsl_matrix* ok_matrix_remove_row(gsl_matrix* m, int row) {
    gsl_matrix* n = gsl_matrix_alloc(m->size1-1, m->size2);
    
    for (int i = 0; i < m->size1; i++)
        for (int j = 0; j < m->size2; j++) {
            if (i < row)
                MSET(n, i, j, MGET(m, i, j));
            else if (i > row)
                MSET(n, i - 1, j, MGET(m, i, j));
        }
    
    return n;
}
gsl_matrix* ok_matrix_remove_column(gsl_matrix* m, int col) {
    gsl_matrix* n = gsl_matrix_alloc(m->size1-1, m->size2);
    
    for (int i = 0; i < m->size1; i++)
        for (int j = 0; j < m->size2; j++) {
            if (j < col)
                MSET(n, i, j, MGET(m, i, j));
            else if (i > col)
                MSET(n, i, j-1, MGET(m, i, j));
        }
    
    return n;
}

gsl_matrix_int* ok_matrix_int_remove_row(gsl_matrix_int* m, int row) {
    gsl_matrix_int* n = gsl_matrix_int_alloc(m->size1-1, m->size2);
    
    for (int i = 0; i < m->size1; i++)
        for (int j = 0; j < m->size2; j++) {
            if (i < row)
                MISET(n, i, j, MIGET(m, i, j));
            else if (i > row)
                MISET(n, i - 1, j, MIGET(m, i, j));
        }
    
    return n;
}

gsl_matrix_int* ok_matrix_int_remove_column(gsl_matrix_int* m, int col) {
    gsl_matrix_int* n = gsl_matrix_int_alloc(m->size1-1, m->size2);
    
    for (int i = 0; i < m->size1; i++)
        for (int j = 0; j < m->size2; j++) {
            if (j < col)
                MISET(n, i, j, MIGET(m, i, j));
            else if (i > col)
                MISET(n, i, j-1, MIGET(m, i, j));
        }
    
    return n;
}


gsl_vector* ok_vector_remove(gsl_vector* m, int idx) {
    gsl_vector* n = gsl_vector_alloc(m->size);
    
    for (int i = 0; i < m->size; i++) {
        if (i < idx)
            VSET(n, i, VGET(m, i));
        else if (i > idx)
            VSET(n, i - 1, VGET(m, i));
    }
    
    return n;
}

gsl_vector_int* ok_vector_int_remove(gsl_vector_int* m, int idx) {
    gsl_vector_int* n = gsl_vector_int_alloc(m->size);
    
    for (int i = 0; i < m->size; i++) {
        if (i < idx)
            VISET(n, i, VIGET(m, i));
        else if (i > idx)
            VISET(n, i - 1, VIGET(m, i));
    }
    
    return n;
}

void ok_matrix_fill(gsl_matrix* src, gsl_matrix* dest) {
    
    for (int i = 0; i < MIN(MROWS(src), MROWS(dest)); i++)
        for (int j = 0; j < MIN(MCOLS(src), MCOLS(dest)); j++)
            MSET(dest, i, j, MGET(src, i, j));
    
}

void ok_parcall(const ok_icallback fnc, const int n) {
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        fnc(i);
    }
}

gsl_matrix* ok_ptr_to_matrix(double* v, unsigned int rows, unsigned int cols) {
    gsl_matrix* m = gsl_matrix_alloc(rows, cols);
    memcpy(m->data, v, sizeof(double)*rows*cols);
    return m;
}

gsl_vector* ok_ptr_to_vector(double* v, unsigned int len) {
    gsl_vector* m = gsl_vector_alloc(len);
    memcpy(m->data, v, sizeof(double)*len);
    return m;
}

gsl_matrix_int* ok_iptr_to_imatrix(int* v, unsigned int rows, unsigned int cols) {
    gsl_matrix_int* m = gsl_matrix_int_alloc(rows, cols);
    memcpy(m->data, v, sizeof(int)*rows*cols);
    return m;
}
gsl_vector_int* ok_iptr_to_ivector(int* v, unsigned int len) {
    gsl_vector_int* m = gsl_vector_int_alloc(len);
    memcpy(m->data, v, sizeof(int)*len);
    return m;
}

void ok_block_to_ptr(void* vv, double* out) {
    gsl_block* v = (gsl_block*) vv;
    for (int i = 0; i < v->size; i++)
        out[i] = v->data[i];
}

void ok_buf_to_ptr(double** v,  unsigned int rows, unsigned int cols, double* out) {
    int k = 0;
    for (int i = 0; i < rows; i++) {
        double* p = v[i];
        for (int j = 0; j < cols; j++) {
            out[k] = p[j];
            k++;
        }
    }
}


void ok_buf_add_to_col(double** buf, double* col_vector, int col, int nrows) {
    for (int i = 0; i < nrows; i++)
        buf[i][col] += col_vector[i];
}

void ok_buf_col(double** buf, double* vector, int col, int nrows) {
    for (int i = 0; i < nrows; i++)
        vector[i] = buf[i][col];
}

unsigned int ok_vector_len(void* v) {
    return ((gsl_vector*) v)->size;
}
unsigned int ok_matrix_rows(void* v) {
    return ((gsl_matrix*) v)->size1;
}
unsigned int ok_matrix_cols(void* v) {
    return ((gsl_matrix*) v)->size2;
}

unsigned int ok_vector_len(void* v);
unsigned int ok_matrix_rows(void* v);
unsigned int ok_matrix_cols(void* v);

gsl_block* ok_vector_block(void* v) {
    return ((gsl_vector*) v)->block;
}
gsl_block* ok_matrix_block(void* v) {
    return ((gsl_matrix*) v)->block;
}

bool ok_file_readable(char* fn) {
    FILE* fid = fopen(fn, "r");
    if (fid == NULL)
        return false;
    fclose(fid);
    return true;
}

ok_rivector* ok_rivector_alloc(const int max_length) {
    ok_rivector* v = (ok_rivector*) malloc(sizeof(ok_rivector));
    v->v = (int*) malloc(max_length*sizeof(int));
    v->max_length = max_length;
    v->length = 0;
    return v;
}

void ok_matrix_column_range(gsl_matrix* m, int col, double* min, double* max) {
    *min = MGET(m, 0, col);
    *max = MGET(m, 0, col);
    
    for (int i = 1; i < MROWS(m); i++) {
        *min = MIN(*min, MGET(m, i, col));
        *max = MAX(*max, MGET(m, i, col));
    }
};

// Re-sample curve based on a triangle area criterion.
// Adapted from http://ariel.chronotext.org/dd/defigueiredo93adaptive.pdf
static inline double _ok_fastrand(unsigned int* seed) {
  *seed = (214013* *seed+2531011);
  *seed = (*seed>>16)&0x7FFF;
  return ((double) *seed / (1.+(double)RAND_MAX));
}

static void _ok_reduce_curve(gsl_matrix* curve, const int timecol,
        const int valcol, const double area_tol, ok_rivector* list, int a,
        int b, unsigned int* seed, const bool log) {
    
    if (b - a <= 1)
        return;
    
    int n = floor((_ok_fastrand(seed) * (0.55 - 0.45) + 0.45) * (b-a) + a);
    
    assert(n >= 0);
    assert(n < MROWS(curve));
    double x1, x2;
    
    if (!log) {
        x1 = MGET(curve, a, timecol) - MGET(curve, n, timecol);
        x2 = MGET(curve, b, timecol) - MGET(curve, n, timecol);
    } else {
        x1 = log10(MGET(curve, a, timecol)) - log10(MGET(curve, n, timecol));
        x2 = log10(MGET(curve, b, timecol)) - log10(MGET(curve, n, timecol));
    };
    
    double y1 = MGET(curve, a, valcol) - MGET(curve, n, valcol);
    double y2 = MGET(curve, b, valcol) - MGET(curve, n, valcol);
    double area = fabs(x1*y2 - x2*y1);
        
    if (area > area_tol) {
        if (a != n && b != n)
            ok_rivector_push(list, n);
        
        _ok_reduce_curve(curve, timecol, valcol, area_tol, list,
                a, n, seed, log);
        _ok_reduce_curve(curve, timecol, valcol, area_tol, list,
                n, b, seed, log);
    } else {
        if (a != n && b != n)
                ok_rivector_push(list, n);
    }
}

static int _ok_rsort_peaks(void* curve_v, const void* a_v, const void* b_v) {
    gsl_matrix* curve = (gsl_matrix*) ((ok_tagdata*) curve_v)->data;
    const int col = (int) ((ok_tagdata*) curve_v)->tag[0];
    
    const int a = *((const int*) a_v);
    const int b = *((const int*) b_v);
    return MGET(curve, b, col) - MGET(curve, a, col);
}
/**
 * Resamples a curve (represented by two columns of the matrix "curve",
 * "xcol" and "ycol" -- the x and y coordinates, respectively) 
 * in order to obtain a curve that retains its "points of interest"
 * (e.g., peaks, points where curvature is high, etc.) while being
 * represented with far fewer samplings. The sampling is then denser
 * where the curve changes rapidly, and sparser where the curve is
 * changing more smoothly.
 * 
 * @param curve Matrix with at least two columns
 * @param xcol "X" column of the curve
 * @param ycol "Y" column of the curve
 * @param peaks_frac Include at least a fraction peaks_frac of all the peaks 
 * in the curve. A peak is a point where Y[i] > Y[i-1] && Y[i] > Y[i+1]. The
 * regular algorithm can miss very steep or narrow peaks, so my algorithm
 * attempts to prevent this. 
 * @param target_points A target amount of points to sample the curve with. The 
 * algorithm will only *try* to meet this target.
 * @param target_tolerance A target tolerance (i.e. the algorithm will try
 * to return target_points +- target_tolerance amount of points)
 * @param start_tolerance An estimated area tolerance; start with a small number
 * (e.g. 1e-3). The routine will return a new estimate for this number 
 * that produces the target number of points. A smaller tolerance
 * yields more points (finer sampling); a tolerance of 0 would
 * return the original curve.
 * @param max_steps The maximum number of iterations to converge
 * @param log_x whether the x-axis is log (e.g. for periodograms)
 * @return A new matrix with two columns containing the resampled curve.
 * (The caller is responsible with freeing both the old and the new matrix.)
 */

gsl_matrix* ok_resample_curve(gsl_matrix* curve, const int xcol, const int ycol, const double peaks_frac, const int target_points,
    const int target_tolerance, double* start_tolerance, const int max_steps, const bool log_x) {
    const int max_points = MROWS(curve);
    
    ok_rivector* peaks = ok_rivector_alloc(max_points);
    ok_rivector* list = ok_rivector_alloc(max_points);
    
    for (int i = 1; i < MROWS(curve)-1; i++)
        if (MGET(curve, i, ycol) > MGET(curve, i-1, ycol) &&
                MGET(curve, i, ycol) > MGET(curve, i+1, ycol)) {
            ok_rivector_push(peaks, i);
        }
    
    ok_tagdata t_curve;
    t_curve.data = curve;
    t_curve.tag[0] = ycol;
    
    ok_qsort_r(ok_rivector_data(peaks), ok_rivector_length(peaks),
            ok_rivector_sizeof(peaks), &t_curve, _ok_rsort_peaks);
    
    int top_peaks_count = ceil(peaks_frac * ok_rivector_length(peaks));
    
    ok_rivector_reset_to(peaks, top_peaks_count);
    ok_rivector_push(peaks, 0);
    ok_rivector_push(peaks, MROWS(curve)-1);
    ok_rivector_sort(peaks);
    
    double xmin, xmax, ymin, ymax;
    xmin = MGET(curve, 0, xcol);
    xmax = MGET(curve, MROWS(curve)-1, xcol);
    ok_matrix_column_range(curve, ycol, &ymin, &ymax);
    if (log_x) {
        xmin = log10(xmin);
        xmax = log10(xmax);
    }
        
    const double area = fabs(xmax-xmin) * fabs(ymax - ymin);
    double tol = (start_tolerance == NULL ? 1e-3 : *start_tolerance);
    double area_tol = area * tol;
    
    unsigned int seed = time(NULL);
    
    
    ok_rivector_append(list, peaks);
    
    const int len = ok_rivector_length(peaks);
    
    for (int i = 1; i < len; i++) {
        _ok_reduce_curve(curve, xcol, ycol, area_tol,
            list, list->v[i-1], list->v[i], &seed,
                log_x);
    }
    
    int steps = max_steps;
    if (steps > 0 && len < 0.9 * max_steps) {
        double x0 = tol;
        double y0 = list->length;
        
        double x1 = *start_tolerance * 0.1;
        double y1 = -1;
        
        while (steps > 0) {
            area_tol = x1 * area;
            ok_rivector_reset_to(list, len);
            for (int i = 1; i < len; i++) {
                _ok_reduce_curve(curve, xcol, ycol, area_tol,
                    list, list->v[i-1], list->v[i], &seed,
                        log_x);
            }
            
            y1 = list->length;
            
            double xn = x1 - (y1 - target_points) * (x1-x0)/((y1 - target_points) - (y0 - target_points));
            if (xn < 0)
                xn = 0.25 * MIN(x0, x1);
            
            x0 = x1;
            y0 = y1;
            x1 = xn;
            
            if (fabs(y1 - target_points) < target_tolerance)
                break;
            
            steps--;
        }   
        
        *start_tolerance = x1;
    } 
    
    ok_rivector_sort(list);
    
    gsl_matrix* new_curve = gsl_matrix_alloc(ok_rivector_length(list),
            MCOLS(curve));
    
    ok_rivector_foreach_i(list, listidx, row) {
        for (int j = 0; j < MCOLS(curve); j++)
            MSET(new_curve, row, j, MGET(curve, listidx, j));
    }
    
    ok_rivector_free(peaks);
    ok_rivector_free(list);
    return new_curve;
}

