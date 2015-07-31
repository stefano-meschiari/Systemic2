
#include <omp.h>

#include "systemic.h"
#include "kernel.h"
#include "utils.h"

#define COPY(from, to, n) for (int __i=0; __i < n; __i++) to[i] = from[i];

int K_minimize_gd(ok_kernel* k, int maxiter, double params[]) {
    int user_status = PROGRESS_CONTINUE;
    ok_kernel_minimizer_pars sp = K_getMinimizedVariables(k);

    double eta = 0.1;
    double alpha = 0.1;
    double eps = 1e-6;

    K_calculate(k);
    double L = K_getLoglik(k);

    double target_DL = 1;

    double dF[sp.npars];
    double dW[sp.npars];
    double dW1[sp.npars];
    double dW2[sp.npars];
    double m[sp.npars];

    int recompute_every = 3;
    int computed = 4;

    for (int i = 0; i < sp.npars; i++) {
        dF[i] = dW[i] = dW1[i] = dW2[i] = 0.;


        for (int j = 0; j < 10; j++) {
            double v = *(sp.pars[i]);

            *(sp.pars[i]) = v + sp.steps[i];

            k->flags |= NEEDS_SETUP;
            K_calculate(k);
            double dL = fabs(K_getLoglik(k) - L);

            *(sp.pars[i]) = v;
            sp.steps[i] = 0.5 * sp.steps[i] * MAX(MIN(1 + target_DL / dL,
                    10), 0.1);
        }
        m[i] = *(sp.pars[i]) / sp.steps[i];
    }



    double L1, L2;
    double grad_eps_abs = 1e-4;
    double grad_eps_rel = 1e-5;
    double grad0 = -1;

    for (int i = 0; i < maxiter; i++) {
        if (k->progress != NULL) {
            k->chi2 = L;
            if (k->progress(i, maxiter, k, __func__) != PROGRESS_CONTINUE) {
                k->flags |= NEEDS_SETUP;
                K_calculate(k);
                return (PROGRESS_STOP);
            }
        }
        if (computed >= recompute_every) {
            computed = 0;
            for (int j = 0; j < sp.npars; j++) {
                *(sp.pars[j]) = (m[j] + eps) * sp.steps[j];
                k->flags |= NEEDS_SETUP;
                K_calculate(k);

                dF[j] = (K_getLoglik(k) - L) / eps;
                dW1[j] = eta * dF[j];
                dW2[j] = eta * dF[j] + alpha * dW[j];
                *(sp.pars[j]) = m[j] * sp.steps[j];

                if (IS_NOT_FINITE(dF[j])) {
                    printf("The gradient on parameter %d is not finite, stopping...\n",
                            j);
                    return 0;
                }
            }
        } else
            computed++;

        if (grad0 < 0)
            grad0 = ok_ptr_sum_2(dF, sp.npars);

        double grad = ok_ptr_sum_2(dF, sp.npars);

        if (grad < grad_eps_abs + grad_eps_rel * grad0)
            break;

        for (int j = 0; j < sp.npars; j++)
            *(sp.pars[j]) = (m[j] - dW1[j]) * sp.steps[j];
        k->flags |= NEEDS_SETUP;
        K_calculate(k);
        L1 = K_getLoglik(k);

        for (int j = 0; j < sp.npars; j++)
            *(sp.pars[j]) = (m[j] - dW2[j]) * sp.steps[j];
        k->flags |= NEEDS_SETUP;
        K_calculate(k);
        L2 = K_getLoglik(k);

        double Lp;

        if (L1 < L2) {
            alpha *= 0.9;
            for (int j = 0; j < sp.npars; j++) {
                dW[j] = dW1[j];
                *(sp.pars[j]) = (m[j] - dW1[j]) * sp.steps[j];
            }
            Lp = L1;
        } else {
            alpha *= 1.1;
            for (int j = 0; j < sp.npars; j++)
                dW[j] = dW2[j];
            Lp = L2;
        }

        if (Lp >= L) {
            eta *= 0.9;
            computed = recompute_every + 1;
            for (int j = 0; j < sp.npars; j++)
                *(sp.pars[j]) = m[j] * sp.steps[j];
        } else {
            eta *= 1.01;
            L = Lp;
            for (int j = 0; j < sp.npars; j++)
                m[j] = m[j] - dW[j];
        }

    }


    return user_status;
}
