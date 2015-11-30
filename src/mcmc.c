#include <gsl/gsl_statistics_double.h>
#ifndef JAVASCRIPT
#include "omp.h"
#else
#include "omp_shim.h"
#endif


#include "mcmc.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#include "math.h"
#include "assert.h"
#include "kl.h"

#define ASSERTDO(x, action) if (!(x)) { action; assert((x)); } 

double K_getElMin(ok_kernel* k, int row, int column, double defMin) {
    double min, max;
    K_getElementRange(k, row, column, &min, &max);
    return (IS_INVALID(min) ? defMin : min);
}

double K_getElMax(ok_kernel* k, int row, int column, double defMax) {
    double min, max;
    K_getElementRange(k, row, column, &min, &max);
    return (IS_INVALID(max) ? defMax : max);
}

double K_getParMin(ok_kernel* k, int column, double defMin) {
    double min, max;
    K_getParRange(k, column, &min, &max);
    return (IS_INVALID(min) ? defMin : min);
}

double K_getParMax(ok_kernel* k, int column, double defMax) {
    double min, max;
    K_getParRange(k, column, &min, &max);
    return (IS_INVALID(max) ? defMax : max);
}

#define UNIFORM_PAR_PRIOR(PAR, DMIN, DMAX) \
        if (K_getElementFlag(k, i, PAR) & MINIMIZE) { \
                double Pmin = K_getElMin(k, i, PAR, DMIN); \
                double Pmax = K_getElMax(k, i, PAR, DMAX); \
                prior /= (Pmax - Pmin); \
        } \


double K_default_prior(ok_kernel* k) {
    double prior = 1;

    for (int i = 1; i < MROWS(k->system->elements); i++) {
        if (K_getElementFlag(k, i, PER) & MINIMIZE) {
            double Pmin = K_getElMin(k, i, PER, 0.2);
            double Pmax = K_getElMax(k, i, PER, 20000.);

            prior /= (MGET(k->system->elements, i, PER) + Pmin) * log(Pmax / Pmin);
        }
        if (K_getElementFlag(k, i, MASS) & MINIMIZE) {
            double Mmin = K_getElMin(k, i, MASS, 1e-5);
            double Mmax = K_getElMax(k, i, MASS, 50.);

            prior /= (MGET(k->system->elements, i, MASS) + Mmin) * log(Mmax / Mmin);
        }

        UNIFORM_PAR_PRIOR(MA, 0, 2 * M_PI);
        UNIFORM_PAR_PRIOR(INC, 0, 2 * M_PI);
        UNIFORM_PAR_PRIOR(NODE, 0, 2 * M_PI);
        UNIFORM_PAR_PRIOR(LOP, 0, 2 * M_PI);
        UNIFORM_PAR_PRIOR(ECC, 0, 0.99);
    }

    for (int i = P_DATA1; i <= P_DATA10; i++)
        prior /= K_getParMax(k, i, 1000.) - K_getParMin(k, i, -1000.);

    for (int i = P_DATA_NOISE1; i <= P_DATA_NOISE10; i++) {
        if (K_getParFlag(k, i) & MINIMIZE) {
            double smax = K_getParMax(k, i, 100.);
            prior /= (fabs(k->params->data[i]) + 0.3) * log((0.3 + smax) / 0.3);
        }
    }
    return prior;
}

void K_mcmc_likelihood_and_prior_default(ok_kernel* k, double* ret) {
    k->flags |= NEEDS_SETUP;
    K_calculate(k);

    double chi2 = K_getChi2_nr(k);

    double prior = K_default_prior(k);
    double A = 0;
    int nd = 0;
    for (int i = 0; i < k->ndata; i++) {
        if (k->compiled[i][T_ERR] >= 0.) {
            int set = (int) k->compiled[i][T_SET];
            double n = K_getPar(k, set + DATA_SETS_SIZE);

            A += log(SQR(k->compiled[i][T_ERR]) + n * n);
            nd++;
        }
    }

    ret[0] = -0.5 * A - 0.5 * chi2 - 0.5 * nd * LOG_2PI;
    ret[1] = log(prior);
}

#define STATE_STEPS 0
#define STATE_MAIN 1
#define STATE_SKIP 2

/**
 * Launches multiple parallel MCMC chains until convergence is achieved; returns a kernel list. The steps are automatically
 * derived by the routine to have a 44% acceptance rate on each minimized parameter. 
 * 
 * @param k Kernel to be used as the starting point. Set minimization flag to MINIMIZE to decide what parameters to vary.
 * @param nchains Number of chains to run in parallel. The ensemble of chains is used to determine convergence
 * @param ntemps Number of temperatures to run in parallel. 
 * @param skip Skip the first 'skip' elements of the chain
 * @param discard Only retain every 'discard'-th element of the chain; the others will be discarded
 * @param params Additional parameters 
 * @param Rstop Chains are considered converged when R < Rstop (usually < 1.2)
 * @param merit_function A function that returns the log of the merit of a given step; set to NULL for default. The default merit function
 * returns log(1/sqrt(A)) - 0.5*chi^2 + log(prior). 
 * @return A chain of systems 
 */
ok_list* K_mcmc_mult(ok_kernel** k, unsigned int nchains, unsigned int ntemps, unsigned int skip, unsigned int discard, const double params[], double Rstop, ok_callback2 merit_function) {
    int Nsteps = 40000;
    int Nmin = 5000;

    int Nstop = -1;

    double tempfac = (1. - 0.1) / ntemps;
    int verbose = 0;

    assert(nchains >= 1);
    assert(ntemps >= 1);


    ok_list * kls[nchains][ntemps];

    ok_progress progress = k[0]->progress;

    double R_single = -1;
    bool return_all = true;
    double acc_ratio = 0.44;
    int save_every = -1;

    bool skip_steps = false;

    int optIdx = 0;
    while (params != NULL) {

        if (params[optIdx] == DONE)
            break;
        else if (params[optIdx] == OPT_MCMC_RSTOP_SINGLE) {
            R_single = params[optIdx + 1];
        } else if (params[optIdx] == OPT_MCMC_RETURN_ALL) {
            return_all = !(params[optIdx + 1] == 0.);
        } else if (params[optIdx] == OPT_MCMC_NSTOP) {
            Nstop = params[optIdx + 1];
        } else if (params[optIdx] == OPT_MCMC_NMIN) {
            Nmin = params[optIdx + 1];
        } else if (params[optIdx] == OPT_MCMC_TEMPFAC) {
            tempfac = params[optIdx + 1];
        } else if (params[optIdx] == OPT_MCMC_VERBOSE_DIAGS) {
            verbose = (int) params[optIdx + 1];
        } else if (params[optIdx] == OPT_MCMC_ACCRATIO) {
            acc_ratio = params[optIdx + 1];
        } else if (params[optIdx] == OPT_MCMC_SKIP_STEPS) {
            skip_steps = ((int) params[optIdx + 1]) != 0;
        } else if (params[optIdx] == OPT_MCMC_SAVE_EVERY) {
            save_every = (int) params[optIdx + 1];
        }
        optIdx += 2;
    }


    double glOpts[ntemps][9];
    glOpts[0][0] = OPT_MCMC_BETA;
    glOpts[0][1] = 1.;
    glOpts[0][2] = OPT_MCMC_VERBOSE_DIAGS;
    glOpts[0][3] = verbose;
    glOpts[0][4] = OPT_MCMC_ACCRATIO;
    glOpts[0][5] = acc_ratio;
    glOpts[0][6] = OPT_MCMC_SKIP_STEPS;
    glOpts[0][7] = (skip_steps ? 1 : 0);
    glOpts[0][8] = DONE;

    for (int i = 1; i < ntemps; i++) {
        glOpts[i][0] = OPT_MCMC_BETA;
        glOpts[i][1] = glOpts[i - 1][1] - tempfac;
        glOpts[i][2] = OPT_MCMC_VERBOSE_DIAGS;
        glOpts[i][3] = verbose;
        glOpts[i][4] = OPT_MCMC_ACCRATIO;
        glOpts[i][5] = acc_ratio;
        glOpts[0][6] = OPT_MCMC_SKIP_STEPS;
        glOpts[0][7] = (skip_steps ? 1 : 0);
        glOpts[i][8] = DONE;
    };
    bool stopped = false;

    #pragma omp parallel for
    for (int n = 0; n < nchains; n++) {
        ok_kernel* k2 = K_clone(k[n]);

        K_calculate(k2);
        int flag;
        for (int j = 0; j < ntemps; j++) {
            kls[n][j] = K_mcmc_single(k2, 1, 0, discard, glOpts[j], NULL, merit_function, n,
                                      &flag);
            kls[n][j]->prototype = k2;
            kls[n][j]->kernels[0]->elements = K_getAllElements(k[n]);
            kls[n][j]->kernels[0]->params = ok_vector_copy(k[n]->params);

            if (flag == PROGRESS_STOP)
                stopped = true;
        }
    }

    if (stopped)
        return NULL;

    bool conv = false;
    bool conv_single = false;

    int iter = 0;
    int save = 0;

    int npars = 0;

    for (int i = 1; i < k[0]->system->nplanets + 1; i++)
        for (int j = 0; j < ELEMENTS_SIZE; j++)
            if (K_getElementFlag(k[0], i, j) & MINIMIZE)
                npars++;

    for (int i = 0; i < PARAMS_SIZE; i++)
        if (K_getParFlag(k[0], i) & MINIMIZE)
            npars++;

    double devs[npars][nchains];
    double avgs[npars][nchains];
    double devs_90[npars][nchains];
    double avgs_90[npars][nchains];
    double devs_2[npars][nchains];
    double avgs_2[npars][nchains];
    double vals[npars][nchains];


    double W[npars];
    double B[npars];
    double R[npars];

    double W_90[npars];
    double B_90[npars];
    double R_90[npars];

    bool isAngle[npars];
    int parType[npars];
    int parLabel[npars];

    for (int i = 0; i < npars; i++)
        isAngle[i] = false;

    int np = 0;
    for (int i = 1; i < k[0]->system->nplanets + 1; i++)
        for (int j = 0; j < ELEMENTS_SIZE; j++) {
            if (K_getElementFlag(k[0], i, j) & MINIMIZE) {
                isAngle[np] = (j == MA || j == LOP || j == INC || j == NODE || j == TRUEANOMALY);
                parType[np] = ELEMENT;
                parLabel[np] = j;
                np++;
            }

        }

    for (int i = 0; i < PARAMS_SIZE; i++)
        if (K_getParFlag(k[0], i) & MINIMIZE) {
            parType[np] = PARAMETER;
            parLabel[np] = i;
            np++;
        }


    double Rmax = 0;
    double Rmax_90 = 0;
    double Rsingle_max = 0;
    while ((!(conv || conv_single)) || (Nmin > kls[0][0]->size)) {

        bool stopped = false;



        #pragma omp parallel for 
        for (int n = 0; n < nchains * ntemps; n++) {

            int ntem = n % ntemps;
            int ncha = n - (n % ntemps);
            int flag;
            ok_list* kl = K_mcmc_single(kls[ncha][0]->prototype, Nsteps, (iter == 0 ? skip : 0),
                                        discard, glOpts[ntem], kls[ncha][ntem], merit_function, ncha, &flag);

            if (flag == PROGRESS_STOP)
                stopped = true;
            KL_append(kls[ncha][ntem], kl);
        }



        #pragma omp parallel for
        for (int n = 0; n < nchains; n++) {
            int size = kls[n][0]->size;

            gsl_matrix* dev = KL_getElementsStats(kls[n][0], STAT_STDDEV);
            gsl_vector* dev_p = KL_getParsStats(kls[n][0], STAT_STDDEV);
            gsl_matrix* avg = KL_getElementsStats(kls[n][0], STAT_MEAN);
            gsl_vector* avg_p = KL_getParsStats(kls[n][0], STAT_MEAN);

            int np = 0;
            for (int i = 1; i < k[0]->system->nplanets + 1; i++)
                for (int j = 0; j < ELEMENTS_SIZE; j++)
                    if (K_getElementFlag(k[0], i, j) & MINIMIZE) {
                        devs[np][n] = MGET(dev, i, j);
                        avgs[np][n] = MGET(avg, i, j);

                        np++;
                    }
            for (int i = 0; i < PARAMS_SIZE; i++)
                if (K_getParFlag(k[0], i) & MINIMIZE) {
                    devs[np][n] = VGET(dev_p, i);
                    avgs[np][n] = VGET(avg_p, i);

                    np++;
                }

            kls[n][0]->size = (int) (0.9 * size);

            gsl_matrix* dev_90 = KL_getElementsStats(kls[n][0], STAT_STDDEV);
            gsl_vector* dev_90_p = KL_getParsStats(kls[n][0], STAT_STDDEV);
            gsl_matrix* avg_90 = KL_getElementsStats(kls[n][0], STAT_MEAN);
            gsl_vector* avg_90_p = KL_getParsStats(kls[n][0], STAT_MEAN);

            np = 0;
            for (int i = 1; i < k[0]->system->nplanets + 1; i++)
                for (int j = 0; j < ELEMENTS_SIZE; j++)
                    if (K_getElementFlag(k[0], i, j) & MINIMIZE) {
                        devs_90[np][n] = MGET(dev_90, i, j);
                        avgs_90[np][n] = MGET(avg_90, i, j);
                        vals[np][n] = KL_getElement(kls[n][0], size - 1, i, j);

                        np++;

                    }
            for (int i = 0; i < PARAMS_SIZE; i++)
                if (K_getParFlag(k[0], i) & MINIMIZE) {
                    devs_90[np][n] = VGET(dev_90_p, i);
                    avgs_90[np][n] = VGET(avg_90_p, i);
                    vals[np][n] = KL_getPar(kls[n][0], size - 1, i);
                    np++;
                }

            kls[n][0]->size = (int) (0.5 * size);

            gsl_matrix* dev_2 = KL_getElementsStats(kls[n][0], STAT_STDDEV);
            gsl_vector* dev_p_2 = KL_getParsStats(kls[n][0], STAT_STDDEV);
            gsl_matrix* avg_2 = KL_getElementsStats(kls[n][0], STAT_MEAN);
            gsl_vector* avg_p_2 = KL_getParsStats(kls[n][0], STAT_MEAN);


            np = 0;
            for (int i = 1; i < k[0]->system->nplanets + 1; i++)
                for (int j = 0; j < ELEMENTS_SIZE; j++)
                    if (K_getElementFlag(k[0], i, j) & MINIMIZE) {
                        devs_2[np][n] = MGET(dev_2, i, j);
                        avgs_2[np][n] = MGET(avg_2, i, j);
                        np++;
                    }
            for (int i = 0; i < PARAMS_SIZE; i++)
                if (K_getParFlag(k[0], i) & MINIMIZE) {
                    devs_2[np][n] = VGET(dev_p_2, i);
                    avgs_2[np][n] = VGET(avg_p_2, i);
                    np++;
                }

            kls[n][0]->size = size;

            gsl_matrix_free(dev);
            gsl_vector_free(dev_p);
            gsl_matrix_free(avg);
            gsl_vector_free(avg_p);

            gsl_matrix_free(dev_2);
            gsl_vector_free(dev_p_2);
            gsl_matrix_free(avg_2);
            gsl_vector_free(avg_p_2);

            gsl_matrix_free(dev_90);
            gsl_vector_free(dev_90_p);
            gsl_matrix_free(avg_90);
            gsl_vector_free(avg_90_p);
        }

        Nsteps = discard * 500;
        conv = true;
        conv_single = (R_single < 0 ? false : true);

        int conv_params = 0;
        int conv_single_param = -1;
        int conv_single_chain = -1;
        int Rmax_param = -1;
        Rmax = 0.;
        Rmax_90 = 0.;
        Rsingle_max = 0.;
        for (int np = 0; np < npars; np++) {
            if (isAngle[np]) {
                W[np] = ok_average_angle(devs[np], nchains, false);
                B[np] = ok_stddev_angle(avgs[np], nchains, false);
                W_90[np] = ok_average_angle(devs_90[np], nchains, false);
                B_90[np] = ok_stddev_angle(avgs_90[np], nchains, false);
            } else {
                W[np] = gsl_stats_mean(devs[np], 1, nchains);
                B[np] = gsl_stats_sd(avgs[np], 1, nchains);
                W_90[np] = gsl_stats_mean(devs_90[np], 1, nchains);
                B_90[np] = gsl_stats_sd(avgs_90[np], 1, nchains);
            }

            R[np] = sqrt((W[np] + B[np]) / W[np]);
            R_90[np] = sqrt((W_90[np] + B_90[np]) / W_90[np]);

            if (Rmax < R[np])
                Rmax_param = np;

            Rmax = MAX(Rmax, R[np]);
            Rmax_90 = MAX(Rmax_90, R_90[np]);

            if (R[np] < Rstop && R_90[np] < Rstop)
                conv_params++;


            for (int j = 0; j < nchains; j++) {
                double R_single_dev = fabs(1 - fabs(devs[np][j] / devs_2[np][j]));


                if (R_single_dev > Rsingle_max) {
                    conv_single_param = np;
                    conv_single_chain = j;
                }
                Rsingle_max = MAX(Rsingle_max, R_single_dev);
            }

            conv = conv && (R[np] < Rstop) && (R_90[np] < Rstop);
        }

        if (R_single > 0)
            conv_single = Rsingle_max < R_single;

        if (progress != NULL) {
            char prog[400];
            sprintf(prog, "[%d] R[%d] = %.2e [1/2 = %.2e], Rsing_max = %.2e [par = %d, chain = %d, v = %e], Rstop = %.2e, size = %d [%d]",
                    kls[0][0]->size, Rmax_param, Rmax, Rmax_90, Rsingle_max, conv_single_param, conv_single_chain, vals[conv_single_param][conv_single_chain], Rstop, kls[0][0]->size, Nstop);

            double p = 100. * (1 - fabs(Rmax - Rstop) / Rstop);

            int progret = progress((int) p, 100, NULL, prog);

            if (progret == PROGRESS_STOP || progret == PROGRESS_BREAK)
                conv = true;
        }
        if (verbose > 1) {
            for (int i = 1; i < k[0]->system->nplanets + 1; i++) {
                for (int j = 0; j < ELEMENTS_SIZE; j++)
                    if (K_getElementFlag(k[0], i, j) & MINIMIZE) {
                        printf("%e ", KL_getElement(kls[conv_single_chain][0],
                                                    kls[conv_single_chain][0]->size - 2,
                                                    i, j));
                    }
                printf("\n");
            }
            for (int i = 0; i < PARAMS_SIZE; i++)
                if (K_getParFlag(k[0], i) & MINIMIZE) {
                    printf("%e ", KL_getPar(kls[conv_single_chain][0],
                                            kls[conv_single_chain][0]->size - 2,
                                            i));
                }
            printf("\n");
        }

        save++;

        if (save_every > 0 && save == save_every) {
            save = 0;
            FILE* fid = fopen("mcmc_stats.txt", "w");
            fprintf(fid, "size = %d\n", kls[0][0]->size);
            fprintf(fid, "R = %e\n", Rmax);
            fprintf(fid, "R_90 = %e\n", Rmax_90);
            fclose(fid);

            for (int i = 0; i < nchains; i++) {
                char fn[80];
                sprintf(fn, "mcmc_%d.txt", i);
                fid = fopen(fn, "w");
                KL_save(kls[i][0], fid);
                fclose(fid);
            }
        }

        if (verbose > 1) {
            printf("Rmax = %e [Rmax_90 = %e], Rsingle_max = %e, Chain length = %d\n", Rmax, Rmax_90, Rsingle_max, kls[0][0]->size);
        }

        if (kls[0][0]->size > Nstop && Nstop > 0) {
            conv = true;
        }

        if (ntemps > 1) {
            for (int n = 0; n < nchains; n++) {
                int j = gsl_rng_uniform_int(k[0]->rng, ntemps - 1);

                double beta_j = glOpts[j][1];
                double beta_j1 = glOpts[j + 1][1];

                double li_j = kls[n][j]->kernels[kls[n][j]->size - 1]->merit_li;
                double li_j1 = kls[n][j + 1]->kernels[kls[n][j + 1]->size - 1]->merit_li;

                double r = (beta_j * li_j1 + beta_j1 * li_j - beta_j * li_j - beta_j1 * li_j1);

                assert(!isnan(r));

                static int iter = 0;
                if (j == 0) iter++;


                if ((r > 0) || (gsl_rng_uniform(k[0]->rng) < MIN(exp(r), 1))) {
                    ok_list_item* it = kls[n][j]->kernels[kls[n][j]->size - 1];

                    kls[n][j]->kernels[kls[n][j]->size - 1] = kls[n][j + 1]->kernels[kls[n][j + 1]->size - 1];
                    kls[n][j + 1]->kernels[kls[n][j + 1]->size - 1] = it;
                }
            }
        }


        iter++;
        if (stopped)
            break;
    }

    if (verbose > 0) {
        printf("Final length: %d, final R_max = %e, final Rsingle_max = %e\n",
               kls[0][0]->size, Rmax, Rsingle_max);
    }

    if (return_all) {
        for (int i = 1; i < nchains; i++) {
            KL_append(kls[0][0], kls[i][0]);
            for (int j = 1; j < ntemps; j++)
                KL_free(kls[i][j]);
        }
    } else {
        for (int i = 1; i < nchains; i++)
            for (int j = 0; j < ntemps; j++)
                KL_free(kls[i][j]);
    }

    return kls[0][0];
}

ok_list* K_mcmc_single(ok_kernel* k2, unsigned int nsteps, unsigned int skip, unsigned int discard, const double dparams[], ok_list* cont, ok_callback2 merit_function, int tag,
                       int* flag) {


    gsl_matrix* plSteps = k2->plSteps;
    gsl_vector* parSteps = k2->parSteps;
    gsl_matrix_int* plFlags = k2->plFlags;
    gsl_vector_int* parFlags = k2->parFlags;

    gsl_matrix* oldEls;
    gsl_vector* oldPars;

    if (cont != NULL && cont->size > 0) {
        oldEls = ok_matrix_copy_sub(cont->kernels[cont->size - 1]->elements,
                                    0, MROWS(cont->kernels[cont->size - 1]->elements),
                                    0, ELEMENTS_SIZE);

        oldPars = ok_vector_copy(cont->kernels[cont->size - 1]->params);
        MATRIX_MEMCPY(k2->system->elements, oldEls);
        VECTOR_MEMCPY(k2->params, oldPars);
        k2->flags |= NEEDS_SETUP;
    } else {
        oldEls = ok_matrix_copy(K_getElements(k2));
        oldPars = ok_vector_copy(k2->params);
    }

    bool skipStepsConvergence = false;
    int optIdx = 0;
    double beta = 1.;
    int verbose = 2;
    double acc_ratio = 0.25;
    int progress_every = (k2->intMethod == KEPLER ? 2000 : 2);

    while (dparams != NULL) {
        if (dparams[optIdx] == DONE)
            break;
        else if (dparams[optIdx] == OPT_MCMC_SKIP_STEPS) {
            skipStepsConvergence = ((int) dparams[optIdx + 1]) != 0;
        } else if (dparams[optIdx] == OPT_MCMC_BETA)
            beta = dparams[optIdx + 1];
        else if (dparams[optIdx] == OPT_MCMC_VERBOSE_DIAGS)
            verbose = dparams[optIdx + 1];
        else if (dparams[optIdx] == OPT_MCMC_ACCRATIO)
            acc_ratio = dparams[optIdx + 1];

        optIdx += 2;
    }


    int nbodies = plSteps->size1;
    int npar = 0;

    for (int j = ELEMENTS_SIZE; j < nbodies * ELEMENTS_SIZE; j++)
        if (plFlags->data[j] & MINIMIZE)
            npar++;

    for (int j = 0; j < PARAMS_SIZE; j++)
        if (parFlags->data[j] & MINIMIZE)
            npar++;

    double acc = 0;
    double acc_par[npar];
    memset(acc_par, 0, sizeof (double) * npar);
    double n_par[npar];
    memset(n_par, 0, sizeof (double) * npar);
    bool conv_par[npar];

    double* steps[npar];
    if (discard < 0)
        discard = 10 * npar;

    ok_list* kl = KL_alloc((int) (MAX(nsteps / discard, 1)), NULL);

    int kpar = 0;
    for (int j = ELEMENTS_SIZE; j < nbodies * ELEMENTS_SIZE; j++)
        if (plFlags->data[j] & MINIMIZE) {
            steps[kpar] = &(plSteps->data[j]);
            conv_par[kpar] = false;
            kpar++;
        }
    for (int j = 0; j < PARAMS_SIZE; j++)
        if (parFlags->data[j] & MINIMIZE) {
            steps[kpar] = &(parSteps->data[j]);
            conv_par[kpar] = false;
            kpar++;
        }

    double prevMerit[2];
    double merit[2];


    (merit_function == NULL ? K_mcmc_likelihood_and_prior_default(k2, prevMerit) : merit_function(k2, prevMerit));

    int state = ((cont == NULL && !skipStepsConvergence) ? STATE_STEPS : STATE_SKIP);

    ok_list_item* it = KL_set(kl, 0, K_getAllElements(k2), ok_vector_copy(oldPars), prevMerit[0] + prevMerit[1], tag);
    it->merit_pr = prevMerit[1];
    it->merit_li = prevMerit[0];

    ok_progress progress = k2->progress;
    *flag = PROGRESS_CONTINUE;

    int it_idx = 0;

    for (unsigned int i = 0; (i <= nsteps) || (state != STATE_MAIN); i++) {
        int sub = 0;
        if (state == STATE_STEPS) {
            sub = gsl_rng_uniform_int(k2->rng, npar);
            n_par[sub] += 1.;
        }

        int par = 0;
        bool invalid = false;

        for (int pl = 1; pl < nbodies; pl++) {
            for (int j = 0; j < ELEMENTS_SIZE; j++) {
                if (MIGET(plFlags, pl, j) & MINIMIZE) {

                    if (state == STATE_MAIN || state == STATE_SKIP || (state == STATE_STEPS && par == sub))
                        MSET(k2->system->elements, pl, j, MGET(oldEls, pl, j) + gsl_ran_gaussian(k2->rng, MGET(plSteps, pl, j)));

                    if ((j == MA) || (j == LOP) || (j == INC) || (j == NODE)) {
                        MSET(k2->system->elements, pl, j, DEGRANGE(MGET(k2->system->elements, pl, j)));
                    }

                    if (k2->plRanges[0] != NULL) {
                        if ((!IS_INVALID(MGET(k2->plRanges[0], pl, j))) && (MGET(k2->system->elements, pl, j) < MGET(k2->plRanges[0], pl, j)))
                            invalid = true;
                        if ((!IS_INVALID(MGET(k2->plRanges[1], pl, j))) && (MGET(k2->system->elements, pl, j) > MGET(k2->plRanges[1], pl, j)))
                            invalid = true;
                    }

                    par++;
                }
            }
        }

        for (int j = 0; j < PARAMS_SIZE; j++) {
            if (parFlags->data[j] & MINIMIZE) {
                if (state == STATE_MAIN || state == STATE_SKIP || (state == STATE_STEPS && par == sub))
                    k2->params->data[j] = oldPars->data[j] + gsl_ran_gaussian(k2->rng, parSteps->data[j]);

                if (j >= P_DATA_NOISE1 && j <= P_DATA_NOISE10)
                    k2->params->data[j] = fabs(k2->params->data[j]);

                if (k2->parRanges[0] != NULL) {
                    if (!IS_INVALID(VGET(k2->parRanges[0], j)) && k2->params->data[j] < k2->parRanges[0]->data[j])
                        invalid = true;
                    if (!IS_INVALID(VGET(k2->parRanges[1], j)) && k2->params->data[j] > k2->parRanges[1]->data[j])
                        invalid = true;
                }
                par++;
            }
        }


        if (invalid) {
            if (i % discard == 0 && state == STATE_MAIN) {
                MATRIX_MEMCPY(k2->system->elements, oldEls);
                VECTOR_MEMCPY(k2->params, oldPars);

                k2->flags |= NEEDS_SETUP;
                it = KL_set(kl, it_idx, K_getAllElements(k2), ok_vector_copy(oldPars), prevMerit[0] + prevMerit[1], tag);
                it->merit_pr = prevMerit[1];
                it->merit_li = prevMerit[0];
                it_idx++;
            }
            continue;
        }

        k2->flags |= NEEDS_SETUP;

        (merit_function == NULL ? K_mcmc_likelihood_and_prior_default(k2, merit) : merit_function(k2, merit));
        ASSERTDO(!isnan(merit[0]), ok_fprintf_matrix(k2->system->elements, stdout, "%e "));
        ASSERTDO(!isnan(merit[1]), ok_fprintf_matrix(k2->system->elements, stdout, "%e "));

        assert(!isinf(merit[0]));
        assert(!isinf(merit[1]));

        double al = MIN(exp(merit[1] - prevMerit[1] + beta * (merit[0] - prevMerit[0])), 1.);
        double u = gsl_rng_uniform(k2->rng);



        if (u < al) {
            prevMerit[0] = merit[0];
            prevMerit[1] = merit[1];

            MATRIX_MEMCPY(oldEls, k2->system->elements);
            VECTOR_MEMCPY(oldPars, k2->params);

            acc += 1.;

            if (state == STATE_STEPS) {
                acc_par[sub] += 1.;
            }

            k2->flags |= NEEDS_SETUP;
        } else {
            MATRIX_MEMCPY(k2->system->elements, oldEls);
            VECTOR_MEMCPY(k2->params, oldPars);

            k2->flags |= NEEDS_SETUP;
        }

        if (i % discard == 0 && state == STATE_MAIN) {
            if (it_idx >= kl->size) {
                break;
            }
            it = KL_set(kl, it_idx, K_getAllElements(k2), ok_vector_copy(oldPars), prevMerit[0] + prevMerit[1], tag);
            it->merit_pr = prevMerit[1];
            it->merit_li = prevMerit[0];
            it_idx++;
        }

        if (state == STATE_SKIP && i > skip) {
            i = 0;
            state = STATE_MAIN;
        }

        int n_conv = 0;
        if (*flag == PROGRESS_STOP)
            break;

        if (n_par[sub] > 2000 && state == STATE_STEPS) {
            if (*steps[sub] < 1e-9)
                *steps[sub] = 1e-9;

            double oldstep = *steps[sub];


            double v = (acc_par[sub] / n_par[sub]) / acc_ratio;
            v = MIN(v, 5);

            *(steps[sub]) *= 0.5 * (1 + v);


            if (verbose > 2 && omp_get_thread_num() == 0) {
                printf("%d %e -> %e [%e]\n", sub, oldstep, *steps[sub], (acc_par[sub] / n_par[sub]));
                fflush(stdout);
            }
            if (fabs(acc_par[sub] / n_par[sub] - acc_ratio) < 0.1 * acc_ratio) {
                conv_par[sub] = true;
            }

            bool conv = true;
            for (int kpar = 0; kpar < npar; kpar++) {
                conv = conv && conv_par[kpar];
                n_conv++;
            }

            if (conv) {
                i = 0;
                state = STATE_SKIP;

                //printf("Steps computed:\n");
                //for (int kpar = 0; kpar < npar; kpar++) {
                //    printf("[%d] %.2e\n", kpar, *(steps[kpar]));
                //}
                //printf("\n");
            }



            acc_par[sub] = n_par[sub] = 0.;
        }

        if (progress != NULL && (i % progress_every == 0) && omp_get_thread_num() == 0) {
            int ret;
            if (state == STATE_STEPS) {
                char msg[500];
                double acc_max = 0;
                int par_max = 0;
                for (int i = 0; i < npar; i++) {
                    double acc = fabs(acc_par[i] / n_par[i]);
                    if (acc > acc_max) {
                        acc_max = acc;
                        par_max = i;
                    }
                }
                sprintf(msg, "Computing step sizes [cur = %.2f, targ = %.2f, par = %d, step = %.2e, n = %.0f]",
                        acc_max, acc_ratio, par_max, *(steps[par_max]), n_par[par_max]);
                ret = progress(n_conv, npar, NULL, msg);
            } else
                ret = progress(i, nsteps, NULL, "");
            *flag = ret;
            if (ret == PROGRESS_STOP)
                break;
        }
    }

    gsl_matrix_free(oldEls);
    gsl_vector_free(oldPars);

    return kl;
}
