
#include <omp.h>

#include "systemic.h"
#include "kernel.h"
#include "utils.h"

#define LARGE_GRAD_JUMP 1000
#define COPY(from, to, n) for (int __i=0; __i < n; __i++) to[i] = from[i];

double K_recalculate_steps(ok_kernel* k, ok_kernel_minimizer_pars* sp, double m[],
        double target_DL, int verbose) {
    K_calculate(k);
    double L = K_getLoglik(k);
    for (int i = 0; i < sp->npars; i++) {
        for (int j = 0; j < 10; j++) {
            double v = *(sp->pars[i]);
            if (verbose > 0)
                printf("i = %d v = %e step = %e\n",
                    i, v, sp->steps[i]);
            *(sp->pars[i]) = v + sp->steps[i];

            k->flags |= NEEDS_SETUP;
            K_calculate(k);
            double dL = fabs(K_getLoglik(k) - L);

            *(sp->pars[i]) = v;
            sp->steps[i] = 0.5 * sp->steps[i] * MAX(MIN(1 + target_DL / dL,
                    10), 0.1);

        }
        m[i] = *(sp->pars[i]) / sp->steps[i];
    }
    return L;
}

int K_minimize_gd(ok_kernel* k, int maxiter, double params[]) {

    int user_status = PROGRESS_CONTINUE;
    ok_kernel_minimizer_pars sp = K_getMinimizedVariables(k);
    int verbose = 0;
    if (params != 0) {
        int i = 0;
        while (params[i] != DONE) {
            if (params[i] == OPT_VERBOSE_DIAGS) {
                verbose = (int) params[i + 1];
            }
            i++;
        }

    }

    if (verbose > 0)
        printf("Starting up...\n");

    double eta = 0.1;
    double alpha = 0.1;
    double eps = 1e-7;

    if (verbose > 0)
        printf("Calculating...\n");
    if (verbose > 0)
        printf("Done.\n");
    double target_DL = 1;

    double dF[sp.npars];
    double dW[sp.npars];
    double dW1[sp.npars];
    double dW2[sp.npars];
    double m[sp.npars];

    int last_improved = 0;
    int last_improved_max = 100;

    int recompute_every = 3;
    int computed = 4;

    for (int i = 0; i < sp.npars; i++) {
        dF[i] = dW[i] = dW1[i] = dW2[i] = 0.;
    }

    double L = K_recalculate_steps(k, &sp, m, target_DL, verbose);

    double L1, L2;
    double grad_eps_abs = 1e-4;
    double grad_eps_rel = 1e-6;
    double grad0 = -1;
    int succ_steps = 0;
    double grad;

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
                *(sp.pars[j]) = m[j] * sp.steps[j];

                if (IS_NOT_FINITE(dF[j])) {
                    printf("The gradient on parameter %d is not finite, stopping...\n",
                            j);
                    return 0;
                }
            }
        } else
            computed++;

        for (int j = 0; j < sp.npars; j++) {
            dW1[j] = eta * dF[j];
            dW2[j] = eta * dF[j] + alpha * dW[j];
        }

        if (grad0 < 0) {
            grad0 = sqrt(ok_ptr_sum_2(dF, sp.npars));
            grad = grad0;
        }
        double prev_grad = grad;
        grad = sqrt(ok_ptr_sum_2(dF, sp.npars));

        if (grad < grad_eps_abs + grad_eps_rel * grad0 || last_improved > last_improved_max) {
            break;
        }

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

        for (int j = 0; j < sp.npars; j++) {
            if (sp.)
            }

        if (Lp >= L || IS_NOT_FINITE(Lp)) {
            eta *= 0.9;
            if (computed == 0)
                computed--;
            else
                computed = recompute_every + 1;
            last_improved++;
            for (int j = 0; j < sp.npars; j++)
                *(sp.pars[j]) = m[j] * sp.steps[j];
        } else {
            eta *= 1.01;
            last_improved = 0;
            L = Lp;
            succ_steps++;
            for (int j = 0; j < sp.npars; j++)
                m[j] = m[j] - dW[j];
        }
        if (verbose)
            printf("step = %d [%d], L = %e [%e], eta = %e, alpha = %e, grad = %e, eps = %e\n",
                i, succ_steps, L, Lp, eta, alpha, grad, eps);
    }

    for (int j = 0; j < sp.npars; j++)
        *(sp.pars[j]) = m[j] * sp.steps[j];
    k->flags |= NEEDS_SETUP;
    K_calculate(k);
    return user_status;
}
