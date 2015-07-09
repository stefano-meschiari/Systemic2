
#include "kl.h"
#include "kernel.h"
#include "integration.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_sort_vector.h"
/*
 * ok_list is an object containing a list of orbital elements and 
 * parameters. It is conceptually similar to an array of ok_kernel* 
 * objects, except that it does not duplicate all the data associated
 * with each ok_kernel* object to avoid wasting memory (this is important,
 * for instance, when returning a chain of states from MCMC, when you
 * can have millions of objects!). They are typically returned by
 * the MCMC or bootstrap algorithms.
 * 
 */

/**
 * Returns a new ok_list of size size, using the specified kernel as a
 * prototype.
 * @param size size of the list
 * @param prototype prototype kernel.
 * @return a freshly allocated list.
 */

ok_list* KL_alloc(const int size, ok_kernel* prototype) {
    ok_list* kl = (ok_list*) malloc(sizeof (ok_list));

    kl->prototype = prototype;
    kl->kernels = (ok_list_item**) calloc(size, sizeof (ok_list_item*));
    kl->size = size;
    kl->diags = NULL;

    return kl;
}

/**
 * Appends two lists together (the src list to the end of the dest list).
 * The src list is subsequently freed.
 * 
 * @param dest List to append to
 * @param src List to be appended
 */

void KL_append(ok_list* dest, ok_list* src) {
    if (dest->prototype != src->prototype)
        K_free(src->prototype);

    int size = dest->size + src->size;

    ok_list_item** destItems = dest->kernels;
    dest->kernels = (ok_list_item**) calloc(size, sizeof (ok_list_item*));

    for (int i = 0; i < dest->size; i++)
        dest->kernels[i] = destItems[i];
    for (int i = 0; i < src->size; i++)
        dest->kernels[i + dest->size] = src->kernels[i];

    dest->size = size;

    free(src->kernels);
    free(src);
}

/**
 * Frees a list.
 * 
 * @param kl list to be freed
 */
void KL_free(ok_list* kl) {
    if (kl == NULL)
        return;

    for (int i = 0; i < kl->size; i++) {
        free(kl->kernels[i]);
    }

    if (kl->prototype != NULL)
        K_free(kl->prototype);
    if (kl->kernels != NULL)
        free(kl->kernels);

    free(kl);
}

/**
 * Loads a list previously saved to a file using KL_save.
 * 
 * @param fid opened file handle 
 * @param skip set to 0
 * @return a list containing the data read from the file
 */
ok_list* KL_load(FILE* fid, int skip) {
    char line[18192];
    int np = -1;

    int trials = -1;
    double Mstar = 0.;
    double epoch = 0.;

    for (int i = 0; i < skip; i++) {
        while (fgets(line, sizeof (line), fid))
            if (strcmp(line, "#End\n") == 0)
                break;
        if (feof(fid))
            return NULL;
    }

    bool found = false;
    char* ret = NULL;

    while (!found) {
        while ((ret = fgets(line, sizeof (line), fid)) != NULL) {
            if (strcmp(line, "#KernelList\n") == 0) {
                found = true;
                break;
            }
            if (feof(fid))
                return NULL;
        }

        if (ret == NULL) {
            return NULL;
        }

    }

    ok_list* kl = NULL;

    while (true) {
        //long int ft = ftell(fid);
        if ((fgets(line, sizeof (line), fid) == NULL) || feof(fid))
            break;

        if (line[0] == '#') {
            char tag[100] = {0};
            sscanf(line + 2, "%s = ", tag);

            if (strcmp(tag, "Planets") == 0) {
                sscanf(line + 2, "%*s = %d", &np);
            } else if (strcmp(tag, "Trials") == 0) {
                sscanf(line + 2, "%*s = %d", &trials);
                kl = KL_alloc(trials, NULL);
            } else if (strcmp(tag, "Epoch") == 0) {
                sscanf(line + 2, "%*s = %le", &epoch);
            } else if (strcmp(tag, "Mstar") == 0) {
                sscanf(line + 2, "%*s = %le", &Mstar);
            } else if (strcmp(tag, "Type") == 0) {
                sscanf(line + 2, "%*s = %d", &(kl->type));
            }

        } else {

            kl->prototype = K_alloc(np);
            K_setEpoch(kl->prototype, epoch);
            K_setMstar(kl->prototype, Mstar);


            for (int tr = 0; tr < trials; tr++) {
                gsl_matrix* elements = gsl_matrix_calloc(np + 1, ALL_ELEMENTS_SIZE);
                gsl_vector* pars = gsl_vector_calloc(PARAMS_SIZE);


                for (int i = 0; i < ALL_ELEMENTS_SIZE; i++) {
                    for (int j = 1; j <= np; j++) {
                        double v = 1e-10;
                        fscanf(fid, "%le", &v);
                        MSET(elements, j, i, v);
                    }
                }


                for (int i = 0; i < PARAMS_SIZE; i++) {
                    double v;
                    fscanf(fid, "%le", &v);
                    VSET(pars, i, v);
                }

                double merit;
                fscanf(fid, "%le", &merit);

                ok_list_item* it = KL_set(kl, tr, elements, pars, merit, 0);
                fscanf(fid, "%le", &it->merit_li);
                fscanf(fid, "%le", &it->merit_pr);
                double tag;
                fscanf(fid, "%le", &tag);
                it->tag = (int) tag;

            }
            break;
        }
    }

    return kl;
}

/**
 * Sets the idx-th element of the list kl to the following values.
 * 
 * @param kl List 
 * @param idx Index of the list element to be modified
 * @param elements Orbital elements (matrix) for the idx-th element
 * @param pars Parameters (vector) for the idx-th element
 * @param merit Merit value for the idx-th element
 * @param tag A tag value (can be used by the user to distinguish elements)
 * @return The item just modified
 */
ok_list_item* KL_set(ok_list* kl, const int idx, gsl_matrix* elements, gsl_vector* pars, double merit, int tag) {
    assert(idx < kl->size);
    if (kl->kernels[idx] != NULL) {
        gsl_matrix_free(kl->kernels[idx]->elements);
        gsl_vector_free(kl->kernels[idx]->params);
        free(kl->kernels[idx]);
    }

    ok_list_item* it = (ok_list_item*) malloc(sizeof (ok_list_item));
    it->elements = elements;
    it->params = pars;
    it->merit = merit;
    it->tag = tag;
    kl->kernels[idx] = it;
    return it;
}

/**
 * Returns a vector of the orbital element el for planet pl, built from
 * the list
 * @param kl
 * @param pl
 * @param el
 * @return 
 */
gsl_vector* KL_getElements(const ok_list* kl, const int pl, const int el) {
    gsl_vector* v = gsl_vector_alloc(kl->size);

    for (int i = 0; i < kl->size; i++)
        VSET(v, i, MGET(kl->kernels[i]->elements, pl, el));
    return v;
}

double KL_getElement(const ok_list* kl, const int index, const int pl, const int el) {

    return MGET(kl->kernels[index]->elements, pl, el);
}

double KL_getPar(const ok_list* kl, const int index, const int what) {

    return VGET(kl->kernels[index]->params, what);
}

/**
 * Returns a vector of the vo-th parameter, built from
 * the list
 * @param kl
 * @param pl
 * @param el
 * @return 
 */
gsl_vector* KL_getPars(const ok_list* kl, const int vo) {
    gsl_vector* v = gsl_vector_alloc(kl->size);

    for (int i = 0; i < kl->size; i++)
        VSET(v, i, VGET(kl->kernels[i]->params, vo));

    return v;
}

/**
 * Get a summary statistic for the orbital elements; for instance,
 * the median value calculated over all the elements of the list.
 * @param kl List
 * @param what Can be one of: STAT_MEAN, STAT_MEDIAN, STAT_STDDEV, STAT_MAD. 
 *      Summary statistic is calculated correctly for angle parameters.
 * @return A matrix whose entries are the summary statistic for the 
 * corresponding orbital element.
 */
gsl_matrix* KL_getElementsStats(const ok_list* kl, const int what) {

    int npl = MROWS(kl->kernels[0]->elements);
    if (npl == 0)
        return NULL;

    gsl_vector* v = gsl_vector_alloc(kl->size);

    gsl_matrix* m = gsl_matrix_alloc(npl, ALL_ELEMENTS_SIZE);
    gsl_matrix_set_all(m, 0.);


    for (int i = 0; i < npl; i++)
        for (int j = 0; j < ALL_ELEMENTS_SIZE; j++) {
            for (int n = 0; n < kl->size; n++) {
                VSET(v, n, MGET(kl->kernels[n]->elements, i, j));
            }

            switch (what) {
                case STAT_MEAN:
                    if (j == MA || j == LOP || j == INC || j == NODE || j == TRUEANOMALY)
                        MSET(m, i, j, ok_average_angle(v->data, v->size, false));
                    else
                        MSET(m, i, j, gsl_stats_mean(v->data, 1, v->size));
                    break;
                case STAT_STDDEV:
                    if (j == MA || j == LOP || j == INC || j == NODE || j == TRUEANOMALY) {
                        MSET(m, i, j, ok_stddev_angle(v->data, v->size, false));
                    } else
                        MSET(m, i, j, gsl_stats_sd(v->data, 1, v->size));
                    break;
                case STAT_MEDIAN:
                    if (j == MA || j == LOP || j == INC || j == NODE || j == TRUEANOMALY)
                        MSET(m, i, j, ok_median_angle(v->data, v->size, false));
                    else {
                        gsl_sort_vector(v);
                        MSET(m, i, j, gsl_stats_median_from_sorted_data(v->data, 1, v->size));
                    }
                    break;
                case STAT_MAD:
                    if (j == MA || j == LOP || j == INC || j == NODE || j == TRUEANOMALY) {
                        double med = ok_median_angle(v->data, v->size, false);
                        MSET(m, i, j, 1.4826 * ok_mad_angle(v->data, v->size, med, false));
                    } else {
                        gsl_sort_vector(v);
                        double med = gsl_stats_median_from_sorted_data(v->data, 1, v->size);

                        MSET(m, i, j, 1.4826 * ok_mad(v->data, v->size, med));
                    }
                    break;
                default:
                    // percentiles
                    gsl_sort_vector(v);
                    MSET(m, i, j, gsl_stats_quantile_from_sorted_data(v->data, 1, v->size, (double) (what) / 100.));
            };
        }
    gsl_vector_free(v);
    return m;
}

/**
 * Get a summary statistic for the parameters; for instance,
 * the median value calculated over all the elements of the list.
 * @param kl List
 * @param what Can be one of: STAT_MEAN, STAT_MEDIAN, STAT_STDDEV, STAT_MAD. 
 * @return A vector whose entries are the summary statistic for the 
 * corresponding orbital parameter.
 */
gsl_vector* KL_getParsStats(const ok_list* kl, const int what) {

    gsl_vector* v = gsl_vector_alloc(kl->size);
    gsl_vector* ret = gsl_vector_calloc(PARAMS_SIZE + 1);


    for (int j = 0; j < PARAMS_SIZE + 1; j++) {
        if (j == PARAMS_SIZE)
            for (int n = 0; n < kl->size; n++) {
                VSET(v, n, kl->kernels[n]->merit);
            } else
            for (int n = 0; n < kl->size; n++) {
                VSET(v, n, VGET(kl->kernels[n]->params, j));
            }

        switch (what) {
            case STAT_MEAN:
                VSET(ret, j, gsl_stats_mean(v->data, 1, v->size));
                break;
            case STAT_STDDEV:
                VSET(ret, j, gsl_stats_sd(v->data, 1, v->size));
                break;
            case STAT_MEDIAN:
                gsl_sort_vector(v);
                VSET(ret, j, gsl_stats_median_from_sorted_data(v->data, 1, v->size));
                break;
            case STAT_MAD:
                gsl_sort_vector(v);
                double med = gsl_stats_median_from_sorted_data(v->data, 1, v->size);
                VSET(ret, j, 1.4826 * ok_mad(v->data, v->size, med));
                break;
            default:
                // percentiles
                gsl_sort_vector(v);
                VSET(v, j, gsl_stats_quantile_from_sorted_data(v->data, 1, v->size, (double) (what) / 100.));
        };
    };

    gsl_vector_free(v);
    return ret;
}

int KL_getSize(const ok_list* kl) {
    return kl->size;
}

int KL_getNplanets(const ok_list* kl) {
    return MROWS(kl->kernels[0]->elements) - 1;
}

void KL_removeAtIndex(ok_list* kl, const int idx) {
    gsl_matrix_free(kl->kernels[idx]->elements);
    gsl_vector_free(kl->kernels[idx]->params);
    free(kl->kernels[idx]);

    for (int i = idx + 1; i < kl->size; i++)
        kl->kernels[i - 1] = kl->kernels[i];
    kl->size--;
}

void KL_fprintf(const ok_list* kl, FILE* out, const char* fmt, const char* lfmt) {
    lfmt = (lfmt != NULL ? lfmt : "%10s%d");

    int np = MROWS(kl->kernels[0]->elements) - 1;
    int vo = PARAMS_SIZE;

    fprintf(out, "# Planets = %d\n", np);
    fprintf(out, "# Trials = %d\n", kl->size);
    fprintf(out, "# Mstar = %e\n", K_getMstar(kl->prototype));
    fprintf(out, "# Epoch = %e\n", K_getEpoch(kl->prototype));



    for (int i = 0; i < ALL_ELEMENTS_SIZE; i++)
        for (int j = 1; j <= np; j++)
            fprintf(out, lfmt, ok_all_orb_labels[i], j);
    for (int i = 0; i < vo; i++)
        fprintf(out, lfmt, "PARAM", i);



    fprintf(out, "\n");

    for (int m = 0; m < kl->size; m++) {

        gsl_matrix* ae = kl->kernels[m]->elements;
        for (int i = 0; i < ALL_ELEMENTS_SIZE; i++)
            for (int j = 1; j <= np; j++)
                fprintf(out, fmt, MGET(ae, j, i));

        for (int i = 0; i < vo; i++)
            fprintf(out, fmt, VGET(kl->kernels[m]->params, i));

        fprintf(out, fmt, kl->kernels[m]->merit);
        fprintf(out, fmt, kl->kernels[m]->merit_pr);
        fprintf(out, fmt, kl->kernels[m]->merit_li);
        fprintf(out, fmt, kl->kernels[m]->tag);
        fprintf(out, " \n");

    }

}

/**
 * Save list to a file
 * @param kl List
 * @param out File handle (already opened for writing)
 */
void KL_save(const ok_list* kl, FILE* out) {
    fprintf(out, "#KernelList\n");
    KL_fprintf(kl, out, "%10.8e ", "%18s%d");
    fprintf(out, "#End\n");
}

/**
 * Flattens a list into a double array. Used by the R wrapper.
 * @param kl
 * @param out
 */
void KL_to_ptr(const ok_list* kl, double* out) {
    int idx = 0;
    int np = KL_getNplanets(kl);

    for (int row = 0; row < kl->size; row++) {
        for (int i = 1; i <= np; i++) {
            for (int j = 0; j < ALL_ELEMENTS_SIZE; j++)
                out[idx++] = MGET(kl->kernels[row]->elements, i, j);
        }


        for (int i = 0; i < PARAMS_SIZE; i++)
            out[idx++] = VGET(kl->kernels[row]->params, i);

        out[idx++] = kl->kernels[row]->merit;
        out[idx++] = kl->kernels[row]->merit_pr;
        out[idx++] = kl->kernels[row]->merit_li;
    }
}

