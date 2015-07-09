/* 
 * File:   kl.h
 * Author: stefano
 *
 * Created on March 1, 2012, 1:13 PM
 */

#ifndef KL_H
#define	KL_H

#ifdef	__cplusplus
extern "C" {
#endif
#include "systemic.h"

    ok_list* KL_alloc(const int size, ok_kernel* prototype);
    void KL_free(ok_list* list);
    ok_list* KL_load(FILE* fid, int skip);
    void KL_save(const ok_list* kl, FILE* out);
    void KL_append(ok_list* dest, ok_list* src);
    gsl_vector* KL_getParsStats(const ok_list* kl, const int what);
    gsl_vector* KL_getElements(const ok_list* kl, const int pl, const int el);
    gsl_vector* KL_getPars(const ok_list* kl, const int vo);
    gsl_matrix* KL_getElementsStats(const ok_list* kl, const int what);
    ok_list_item* KL_set(ok_list* kl, const int idx, gsl_matrix* elements, gsl_vector* pars, double merit, int tag);
    int KL_getSize(const ok_list* kl);
    void KL_removeAtIndex(ok_list* kl, const int idx);
    void KL_fprintf(const ok_list* kl, FILE* out, const char* fmt, const char* lfmt);
    void KL_to_ptr(const ok_list* kl, double* out);
    int KL_getNplanets(const ok_list* kl);
    double KL_getElement(const ok_list* kl, const int index, const int pl, const int el);
    double KL_getPar(const ok_list* kl, const int index, const int what);

#ifdef	__cplusplus
}
#endif

#endif	/* KL_H */

