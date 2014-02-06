# Constants
K_SYSTEMIC_VERSION <- 2.1400
K_MAX_LINE <- 8192
K_T_RV <- 0
K_T_PHOTO <- 1
K_T_TIMING <- 2
K_T_DUMMY <- 99
K_DATA_SIZE <- 11
K_DATA_SETS_SIZE <- 10
K_T_TIME <- 0
K_T_VAL <- 1
K_T_ERR <- 2
K_T_TDS_PLANET <- 3
K_T_TDS_FLAG <- 4
K_T_FLAG <- 6
K_T_SVAL <- 7
K_T_PRED <- 8
K_T_SET <- 9
K_T_GUI <- 10
K_T_SCRATCH <- 10
K_X <- 1
K_Y <- 2
K_Z <- 3
K_VX <- 4
K_VY <- 5
K_VZ <- 6
K_STAT_MEAN <- 0
K_STAT_MEDIAN <- -1
K_STAT_STDDEV <- -2
K_STAT_MAD <- -3
K_STAT_QUART_25 <- 25
K_STAT_QUART_50 <- 50
K_STAT_QUART_75 <- 75
K_KEPLER <- 0
K_RK45 <- 1
K_RK89 <- 2
K_ADAMS <- 3
K_BULIRSCHSTOER <- 4
K_SWIFTRMVS <- 5
K_AU <- 1.4959787e13
K_MSUN <- 1.98892e33
K_MJUP <- 1.8986e30
K_RJUP <- 7.1e9
K_GGRAV <- 6.67384e-8
K_DAY <- 8.64e4
K_TWOPI <- 6.2831853072e+00
K_SQRT_TWOPI <- 2.5066282746e+00
K_YEAR <- 31556926.
K_PER <- 0
K_MASS <- 1
K_MA <- 2
K_ECC <- 3
K_LOP <- 4
K_INC <- 5
K_NODE <- 6
K_RADIUS <- 7
K_ORD <- 8
K_SMA <- 13
K_SEMIAMP <- 14
K_TPERI <- 15
K_TRUEANOMALY <- 16
K_MEANLONGITUDE <- 17
K_ELEMENTS_SIZE <- 13
K_ALL_ELEMENTS_SIZE <- 20
K_P_DATA1 <- 0
K_P_DATA2 <- 1
K_P_DATA3 <- 2
K_P_DATA4 <- 3
K_P_DATA5 <- 4
K_P_DATA6 <- 5
K_P_DATA7 <- 6
K_P_DATA8 <- 7
K_P_DATA9 <- 8
K_P_DATA10 <- 9
K_P_DATA_NOISE1 <- 10
K_P_DATA_NOISE2 <- 11
K_P_DATA_NOISE3 <- 12
K_P_DATA_NOISE4 <- 13
K_P_DATA_NOISE5 <- 14
K_P_DATA_NOISE6 <- 15
K_P_DATA_NOISE7 <- 16
K_P_DATA_NOISE8 <- 17
K_P_DATA_NOISE9 <- 18
K_P_DATA_NOISE10 <- 19
K_P_RV_TREND <- 20
K_PARAMS_SIZE <- 100
K_OPT_EPS <- 0
K_OPT_ECC_LAST <- 1
K_OPT_MCMC_SKIP_STEPS <- 0
K_OPT_MCMC_SKIP <- 1
K_OPT_MCMC_RSTOP_SINGLE <- 2
K_OPT_MCMC_RETURN_ALL <- 3
K_OPT_MCMC_BETA <- 4
K_OPT_MCMC_NSTOP <- 5
K_OPT_MCMC_TEMPFAC <- 6
K_OPT_MCMC_VERBOSE_DIAGS <- 7
K_OPT_MCMC_ACCRATIO <- 8
K_OPT_MCMC_NMIN <- 9
K_OPT_LM_MINCHI_PAR <- 10
K_OPT_LM_HIGH_DF <- 11
K_OPT_LM_MAX_ITER_AT_SCALE <- 12
K_OPT_LM_INITIAL_SCALE <- 13
K_OPT_SA_T0 <- 20
K_OPT_SA_ALPHA <- 21
K_OPT_SA_CHAINS <- 22
K_OPT_SA_AUTO_STEPS <- 23
K_OPT_DE_CR <- 30
K_OPT_DE_NP_FAC <- 31
K_OPT_DE_F_MIN <- 32
K_OPT_DE_F_MAX <- 33
K_OPT_DE_USE_STEPS <- 34
K_PROGRESS_CONTINUE <- 0
K_PROGRESS_STOP <- 1
K_TDS_PRIMARY <- 1
K_TDS_SECONDARY <- 2
K_DONE <- -1
K_ASTROCENTRIC <- 0
K_SIMPLEX <- 0
K_LM <- 1
K_DIFFEVOL <- 2
K_SA <- 3
K_INTEGRATION_SUCCESS <- 0
K_INTEGRATORS_SIZE <- 4
K_M_PI <- 3.14159265e+00
K_OK_SUCCESS <- 0
K_OK_NOCONV <- 1
K_PS_TIME <- 0
K_PS_Z <- 1
K_PS_FAP <- 2
K_PS_Z_LS <- 3
K_PS_TAU <- 4
K_PS_WIN <- 5
K_PS_SIZE <- 6
K_PS_TYPE_DATA <- 0
K_PS_TYPE_RESIDUALS <- 1
K_T_INAPPLICABLE <- -1
K_T_STABLE <- 0
K_T_UNSTABLE <- 1
K_ACTIVE <- 2
K_MINIMIZE <- 4
K_NEEDS_COMPILE <- 2
K_NEEDS_SETUP <- 4
K_BOOTSTRAP_DATA <- 8
K_RETAIN <- 16
K_SHARE_FLAGS <- 32
K_SHARE_STEPS <- 64
K_SHARE_DATA <- 128
K_SHARE_RANGES <- 2048
K_FREED <- 256
K_DUPLICATE <- 512
K_GUI_RESERVED_1 <- 4096
K_GUI_RESERVED_2 <- 8192
K_GUI_RESERVED_3 <- 16384
K_GUI_RESERVED_4 <- 32768
K_JACOBI <- 1024
K_INTEGRATION_FAILURE_SMALL_TIMESTEP <- 2048
K_INTEGRATION_FAILURE_INCREASE_TOLERANCE <- 4096
K_INTEGRATION_FAILURE_CLOSE_ENCOUNTER <- 8192
K_INTEGRATION_FAILURE_CLOSE_ENCOUNTER_STAR <- 16384
K_INTEGRATION_FAILURE_STOPPED <- 32768
K_INTEGRATION_FAILURE_SWIFT <- 65536
# Structs 
parseStructInfos("gsl_block{L*d}size data;")
parseStructInfos("gsl_matrix{LLL*d*<gsl_block>i}size1 size2 tda data block owner;")
parseStructInfos("gsl_block_int{L*i}size data;")
parseStructInfos("gsl_matrix_int{LLL*i*<gsl_block_int>i}size1 size2 tda data block owner;")
parseStructInfos("gsl_vector_int{LL*i*<gsl_block_int>i}size stride data block owner;")
parseStructInfos("gsl_vector{LL*d*<gsl_block>i}size stride data block owner;")
parseStructInfos("ok_kernel{pddddddipppppippppppdiiiiippIppdpp}system chi2 chi2_rvs chi2_tts rms rms_tts jitter nsets datasets compiled params times integration integrationSamples plFlags parFlags plSteps parSteps plRanges parRanges Mstar ndata nrvs ntts npars intMethod intOptions minfunc flags rng tag tagValue progress datanames;")
# Function signatures
.lib <- dynbind(c("libsystemic.so", "libsystemic.dylib"), paste(sep=";",
# double DEGRANGE(double angle)
"DEGRANGE(d)d",
# double RADRANGE(double angle)
"RADRANGE(d)d",
# void ok_sort_small_matrix(gsl_matrix* matrix, const int column)
"ok_sort_small_matrix(*<gsl_matrix>i)v",
# void ok_sort_matrix(gsl_matrix* matrix, const int column)
"ok_sort_matrix(*<gsl_matrix>i)v",
# void ok_rsort_matrix(gsl_matrix* matrix, const int column)
"ok_rsort_matrix(*<gsl_matrix>i)v",
# void ok_fprintf_matrix(gsl_matrix* matrix, FILE* file, const char* format)
"ok_fprintf_matrix(*<gsl_matrix>*<FILE>Z)v",
# void ok_fprintf_vector(gsl_vector* v, FILE* file, const char* format)
"ok_fprintf_vector(*<gsl_vector>*<FILE>Z)v",
# void ok_fprintf_matrix_int(gsl_matrix_int* matrix, FILE* file, const char* format)
"ok_fprintf_matrix_int(*<gsl_matrix_int>*<FILE>Z)v",
# void ok_fprintf_vector_int(gsl_vector_int* v, FILE* file, const char* format)
"ok_fprintf_vector_int(*<gsl_vector_int>*<FILE>Z)v",
# void ok_fprintf_buf(double** buf, FILE* file, const char* format, int rows, int columns)
"ok_fprintf_buf(p*<FILE>Zii)v",
# void ok_fscanf_matrix(gsl_matrix* matrix, FILE* file)
"ok_fscanf_matrix(*<gsl_matrix>*<FILE>)v",
# void ok_fscanf_matrix_int(gsl_matrix_int* matrix, FILE* file)
"ok_fscanf_matrix_int(*<gsl_matrix_int>*<FILE>)v",
# void ok_fscanf_vector(gsl_vector* matrix, FILE* file)
"ok_fscanf_vector(*<gsl_vector>*<FILE>)v",
# void ok_fscanf_vector_int(gsl_vector_int* matrix, FILE* file)
"ok_fscanf_vector_int(*<gsl_vector_int>*<FILE>)v",
# bool ok_save_matrix(gsl_matrix* matrix, FILE* fid, const char* format)
"ok_save_matrix(*<gsl_matrix>*<FILE>Z)B",
# bool ok_save_matrix_bin(gsl_matrix* matrix, FILE* fid)
"ok_save_matrix_bin(*<gsl_matrix>*<FILE>)B",
# bool ok_save_buf(double** matrix, FILE* fid, const char* format, int rows, int cols)
"ok_save_buf(p*<FILE>Zii)B",
# bool ok_save_buf_bin(double** matrix, FILE* fid, int rows, int cols)
"ok_save_buf_bin(p*<FILE>ii)B",
# gsl_matrix* ok_read_table(FILE* fid)
"ok_read_table(*<FILE>)*<gsl_matrix>",
# gsl_vector* ok_vector_resize(gsl_vector* v, int len)
"ok_vector_resize(*<gsl_vector>i)*<gsl_vector>",
# gsl_matrix* ok_matrix_resize(gsl_matrix* v, int rows, int cols)
"ok_matrix_resize(*<gsl_matrix>ii)*<gsl_matrix>",
# gsl_vector* ok_vector_resize_pad(gsl_vector* v, int len, double pad)
"ok_vector_resize_pad(*<gsl_vector>id)*<gsl_vector>",
# gsl_matrix* ok_matrix_resize_pad(gsl_matrix* v, int rows, int cols, double pad)
"ok_matrix_resize_pad(*<gsl_matrix>iid)*<gsl_matrix>",
# gsl_vector_int* ok_vector_int_resize(gsl_vector_int* v, int len)
"ok_vector_int_resize(*<gsl_vector_int>i)*<gsl_vector_int>",
# gsl_matrix_int* ok_matrix_int_resize(gsl_matrix_int* v, int rows, int cols)
"ok_matrix_int_resize(*<gsl_matrix_int>ii)*<gsl_matrix_int>",
# gsl_matrix* ok_matrix_copy(const gsl_matrix* src)
"ok_matrix_copy(*<gsl_matrix>)*<gsl_matrix>",
# gsl_matrix* ok_matrix_copy_sub(const gsl_matrix* src, int row1, int nrows, int col1, int ncols)
"ok_matrix_copy_sub(*<gsl_matrix>iiii)*<gsl_matrix>",
# gsl_matrix_int* ok_matrix_int_copy(const gsl_matrix_int* src)
"ok_matrix_int_copy(*<gsl_matrix_int>)*<gsl_matrix_int>",
# gsl_vector* ok_vector_copy(const gsl_vector* src)
"ok_vector_copy(*<gsl_vector>)*<gsl_vector>",
# gsl_vector_int* ok_vector_int_copy(const gsl_vector_int* src)
"ok_vector_int_copy(*<gsl_vector_int>)*<gsl_vector_int>",
# gsl_matrix* ok_buf_to_matrix(double** buf, int rows, int cols)
"ok_buf_to_matrix(pii)*<gsl_matrix>",
# void ok_buf_col(double** buf, double* vector, int col, int nrows)
"ok_buf_col(p*dii)v",
# void ok_matrix_column_range(gsl_matrix* m, int col, double* min, double* max)
"ok_matrix_column_range(*<gsl_matrix>i*d*d)v",
# gsl_matrix* ok_matrix_remove_row(gsl_matrix* m, int row)
"ok_matrix_remove_row(*<gsl_matrix>i)*<gsl_matrix>",
# gsl_matrix* ok_matrix_remove_column(gsl_matrix* m, int col)
"ok_matrix_remove_column(*<gsl_matrix>i)*<gsl_matrix>",
# gsl_matrix_int* ok_matrix_int_remove_row(gsl_matrix_int* m, int row)
"ok_matrix_int_remove_row(*<gsl_matrix_int>i)*<gsl_matrix_int>",
# gsl_matrix_int* ok_matrix_int_remove_column(gsl_matrix_int* m, int col)
"ok_matrix_int_remove_column(*<gsl_matrix_int>i)*<gsl_matrix_int>",
# gsl_vector* ok_vector_remove(gsl_vector* m, int idx)
"ok_vector_remove(*<gsl_vector>i)*<gsl_vector>",
# gsl_vector_int* ok_vector_int_remove(gsl_vector_int* m, int idx)
"ok_vector_int_remove(*<gsl_vector_int>i)*<gsl_vector_int>",
# void ok_matrix_fill(gsl_matrix* src, gsl_matrix* dest)
"ok_matrix_fill(*<gsl_matrix>*<gsl_matrix>)v",
# void ok_bootstrap_matrix(const gsl_matrix* matrix, gsl_matrix* out, const int sortcol, gsl_rng* r)
"ok_bootstrap_matrix(*<gsl_matrix>*<gsl_matrix>i*<gsl_rng>)v",
# void ok_bootstrap_matrix_mean(const gsl_matrix* matrix, int timecol, int valcol, gsl_matrix* out, gsl_rng* r)
"ok_bootstrap_matrix_mean(*<gsl_matrix>ii*<gsl_matrix>*<gsl_rng>)v",
# gsl_matrix* ok_matrix_filter(gsl_matrix* matrix, const int column, const double filter)
"ok_matrix_filter(*<gsl_matrix>id)*<gsl_matrix>",
# gsl_matrix* ok_matrix_buf_filter(double** matrix, const int rows, const int columns, const int column, const double filter)
"ok_matrix_buf_filter(piiid)*<gsl_matrix>",
# int ok_bsearch(double* v, double val, int len)
"ok_bsearch(*ddi)i",
# double ok_average_angle(const double* v, const int length, const bool isRadians)
"ok_average_angle(*diB)d",
# double ok_median_angle(const double* v, const int length, const bool isRadians)
"ok_median_angle(*diB)d",
# double ok_stddev_angle(const double* v, const int length, const bool isRadians)
"ok_stddev_angle(*diB)d",
# double ok_mad_angle(double* v, const int length, const double med, const bool isRadians)
"ok_mad_angle(*didB)d",
# double ok_mad(double* v, const int length, const double med)
"ok_mad(*did)d",
# char* ok_str_copy(const char* src)
"ok_str_copy(Z)*<char>",
# char* ok_str_cat(const char* a1, const char* a2)
"ok_str_cat(ZZ)*<char>",
# void ok_avevar(const double* v, int len, double* ave, double* var)
"ok_avevar(*di*d*d)v",
# gsl_matrix* ok_ptr_to_matrix(double* v, unsigned int rows, unsigned int cols)
"ok_ptr_to_matrix(*dII)*<gsl_matrix>",
# gsl_vector* ok_ptr_to_vector(double* v, unsigned int len)
"ok_ptr_to_vector(*dI)*<gsl_vector>",
# gsl_matrix_int* ok_iptr_to_imatrix(int* v, unsigned int rows, unsigned int cols)
"ok_iptr_to_imatrix(*iII)*<gsl_matrix_int>",
# gsl_vector_int* ok_iptr_to_ivector(int* v, unsigned int len)
"ok_iptr_to_ivector(*iI)*<gsl_vector_int>",
# void ok_block_to_ptr(void* vv, double* out)
"ok_block_to_ptr(p*d)v",
# void ok_buf_to_ptr(double** v, unsigned int rows, unsigned int cols, double* out)
"ok_buf_to_ptr(pII*d)v",
# void ok_buf_add_to_col(double** buf, double* col_vector, int col, int nrows)
"ok_buf_add_to_col(p*dii)v",
# unsigned int ok_vector_len(void* v)
"ok_vector_len(p)I",
# unsigned int ok_matrix_rows(void* v)
"ok_matrix_rows(p)I",
# unsigned int ok_matrix_cols(void* v)
"ok_matrix_cols(p)I",
# gsl_block* ok_vector_block(void* v)
"ok_vector_block(p)*<gsl_block>",
# gsl_block* ok_matrix_block(void* v)
"ok_matrix_block(p)*<gsl_block>",
# gsl_matrix* ok_resample_curve(gsl_matrix* curve, const int xcol, const int ycol, const double peaks_frac, const int target_points,     const int target_tolerance, double* start_tolerance, const int max_steps, const bool log_x)
"ok_resample_curve(*<gsl_matrix>iidii*diB)*<gsl_matrix>",
# bool ok_file_readable(char* fn)
"ok_file_readable(*c)B",
# ok_rivector* ok_rivector_alloc(const int maxlength)
"ok_rivector_alloc(i)*<ok_rivector>",
# ok_kernel* K_alloc()
"K_alloc()*<ok_kernel>",
# void K_free(ok_kernel* k)
"K_free(p)v",
# ok_kernel* K_clone(ok_kernel* k)
"K_clone(p)*<ok_kernel>",
# ok_kernel* K_cloneFlags(ok_kernel* k, unsigned int shareFlags)
"K_cloneFlags(pI)*<ok_kernel>",
# gsl_matrix* K_addDataFile(ok_kernel* k, const char* path, int type)
"K_addDataFile(pZi)*<gsl_matrix>",
# gsl_matrix* K_addDataTable(ok_kernel* k, gsl_matrix* rvtable, const char* name, int type)
"K_addDataTable(p*<gsl_matrix>Zi)*<gsl_matrix>",
# void K_removeData(ok_kernel* k, int idx)
"K_removeData(pi)v",
# gsl_matrix* K_getData(ok_kernel* k, int idx)
"K_getData(pi)*<gsl_matrix>",
# void K_setData(ok_kernel* k, int idx, gsl_matrix* data)
"K_setData(pi*<gsl_matrix>)v",
# const char* K_getDataName(ok_kernel* k, int idx)
"K_getDataName(pi)Z",
# bool K_addDataFromSystem(ok_kernel* k, const char* filename)
"K_addDataFromSystem(pZ)B",
# int K_getDataType(ok_kernel* k, int idx)
"K_getDataType(pi)i",
# int K_getDataSize(ok_kernel* k, int idx)
"K_getDataSize(pi)i",
# double** K_compileData(ok_kernel* k)
"K_compileData(p)p",
# gsl_matrix* K_getCompiledDataMatrix(ok_kernel* k)
"K_getCompiledDataMatrix(p)*<gsl_matrix>",
# void K_addPlanet(ok_kernel* k, const double elements[])
"K_addPlanet(p*d)v",
# void K_removePlanet(ok_kernel* k, int idx)
"K_removePlanet(pi)v",
# void K_setElement(ok_kernel* k, int row, int col, double value)
"K_setElement(piid)v",
# double K_getElement(ok_kernel* k, int row, int col)
"K_getElement(pii)d",
# void K_setElements(ok_kernel* k, gsl_matrix* elements)
"K_setElements(p*<gsl_matrix>)v",
# gsl_matrix* K_getElements(ok_kernel* k)
"K_getElements(p)*<gsl_matrix>",
# void K_setElementFlag(ok_kernel* k, int row, int col, int value)
"K_setElementFlag(piii)v",
# int K_getElementFlag(ok_kernel* k, int row, int col)
"K_getElementFlag(pii)i",
# void K_setElementStep(ok_kernel* k, int row, int col, double value)
"K_setElementStep(piid)v",
# double K_getElementStep(ok_kernel* k, int row, int col)
"K_getElementStep(pii)d",
# void K_setElementRange(ok_kernel* k, int row, int col, double min, double max)
"K_setElementRange(piidd)v",
# void K_getElementRange(ok_kernel* k, int row, int col, double* min, double* max)
"K_getElementRange(pii*d*d)v",
# gsl_matrix* K_getAllElements(ok_kernel* k)
"K_getAllElements(p)*<gsl_matrix>",
# int K_getActivePars(ok_kernel* k)
"K_getActivePars(p)i",
# int K_getActiveElements(ok_kernel* k)
"K_getActiveElements(p)i",
# int K_getNrPars(ok_kernel* k)
"K_getNrPars(p)i",
# void K_setElementType(ok_kernel* k, int type)
"K_setElementType(pi)v",
# int K_getElementType(ok_kernel* k)
"K_getElementType(p)i",
# gsl_matrix* K_getXYZ(ok_kernel* k)
"K_getXYZ(p)*<gsl_matrix>",
# void K_setMstar(ok_kernel* k, double value)
"K_setMstar(pd)v",
# double K_getMstar(ok_kernel* k)
"K_getMstar(p)d",
# void K_setEpoch(ok_kernel* k, double value)
"K_setEpoch(pd)v",
# double K_getEpoch(ok_kernel* k)
"K_getEpoch(p)d",
# double K_getChi2(ok_kernel* k)
"K_getChi2(p)d",
# double K_getChi2_nr(ok_kernel* k)
"K_getChi2_nr(p)d",
# double K_getLoglik(ok_kernel* k)
"K_getLoglik(p)d",
# double K_getRms(ok_kernel* k)
"K_getRms(p)d",
# double K_getJitter(ok_kernel* k)
"K_getJitter(p)d",
# double** K_getCompiled(ok_kernel* k)
"K_getCompiled(p)p",
# unsigned int K_getNrvs(ok_kernel* k)
"K_getNrvs(p)I",
# unsigned int K_getNtts(ok_kernel* k)
"K_getNtts(p)I",
# unsigned int K_getNsets(ok_kernel* k)
"K_getNsets(p)I",
# double K_getChi2_rvs(ok_kernel* k)
"K_getChi2_rvs(p)d",
# double K_getChi2_tts(ok_kernel* k)
"K_getChi2_tts(p)d",
# double K_getRms_tts(ok_kernel* k)
"K_getRms_tts(p)d",
# void K_setFlags(ok_kernel* k, unsigned int value)
"K_setFlags(pI)v",
# unsigned int K_getFlags(ok_kernel* k)
"K_getFlags(p)I",
# void K_setMinFunc(ok_kernel* k, ok_callback value)
"K_setMinFunc(pp)v",
# ok_callback K_getMinFunc(ok_kernel* k)
"K_getMinFunc(p)p",
# void K_setElementSteps(ok_kernel* k, gsl_matrix* value)
"K_setElementSteps(p*<gsl_matrix>)v",
# gsl_matrix* K_getElementSteps(ok_kernel* k)
"K_getElementSteps(p)*<gsl_matrix>",
# void K_setParSteps(ok_kernel* k, gsl_vector* value)
"K_setParSteps(p*<gsl_vector>)v",
# gsl_vector* K_getParSteps(ok_kernel* k)
"K_getParSteps(p)*<gsl_vector>",
# void K_setElementFlags(ok_kernel* k, gsl_matrix_int* value)
"K_setElementFlags(p*<gsl_matrix_int>)v",
# gsl_matrix_int* K_getElementFlags(ok_kernel* k)
"K_getElementFlags(p)*<gsl_matrix_int>",
# void K_setParFlags(ok_kernel* k, gsl_vector_int* value)
"K_setParFlags(p*<gsl_vector_int>)v",
# gsl_vector_int* K_getParFlags(ok_kernel* k)
"K_getParFlags(p)*<gsl_vector_int>",
# void K_setIntMethod(ok_kernel* k, int value)
"K_setIntMethod(pi)v",
# int K_getIntMethod(ok_kernel* k)
"K_getIntMethod(p)i",
# void K_setProgress(ok_kernel* k, ok_progress value)
"K_setProgress(pp)v",
# ok_progress K_getProgress(ok_kernel* k)
"K_getProgress(p)p",
# void K_setIntOptions(ok_kernel* k, ok_integrator_options* value)
"K_setIntOptions(p*<ok_integrator_options>)v",
# ok_integrator_options* K_getIntOptions(ok_kernel* k)
"K_getIntOptions(p)*<ok_integrator_options>",
# void K_setCustomModelFunction(ok_kernel* k, ok_model_function value)
"K_setCustomModelFunction(pp)v",
# ok_model_function K_getCustomModelFunction(ok_kernel* k)
"K_getCustomModelFunction(p)p",
# void K_setIntAbsAcc(ok_kernel* k, double value)
"K_setIntAbsAcc(pd)v",
# double K_getIntAbsAcc(ok_kernel* k)
"K_getIntAbsAcc(p)d",
# void K_setIntRelAcc(ok_kernel* k, double value)
"K_setIntRelAcc(pd)v",
# double K_getIntRelAcc(ok_kernel* k)
"K_getIntRelAcc(p)d",
# void K_setIntDt(ok_kernel* k, double value)
"K_setIntDt(pd)v",
# double K_getIntDt(ok_kernel* k)
"K_getIntDt(p)d",
# unsigned int K_getNplanets(ok_kernel* k)
"K_getNplanets(p)I",
# unsigned int K_getNdata(ok_kernel* k)
"K_getNdata(p)I",
# void K_getRange(ok_kernel* k, double* from, double* to)
"K_getRange(p*d*d)v",
# void K_perturb(ok_kernel* k)
"K_perturb(p)v",
# void K_setPars(ok_kernel* k, gsl_vector* pars)
"K_setPars(p*<gsl_vector>)v",
# gsl_vector* K_getPars(ok_kernel* k)
"K_getPars(p)*<gsl_vector>",
# void K_setPar(ok_kernel* k, int idx, double val)
"K_setPar(pid)v",
# double K_getPar(ok_kernel* k, int idx)
"K_getPar(pi)d",
# void K_setParFlag(ok_kernel* k, int idx, int value)
"K_setParFlag(pii)v",
# int K_getParFlag(ok_kernel* k, int idx)
"K_getParFlag(pi)i",
# void K_setParStep(ok_kernel* k, int idx, double value)
"K_setParStep(pid)v",
# double K_getParStep(ok_kernel* k, int idx)
"K_getParStep(pi)d",
# void K_setParRange(ok_kernel* k, int idx, double min, double max)
"K_setParRange(pidd)v",
# void K_getParRange(ok_kernel* k, int idx, double* min, double* max)
"K_getParRange(pi*d*d)v",
# bool K_save(ok_kernel* k, FILE* fid)
"K_save(p*<FILE>)B",
# ok_kernel* K_load(FILE* fid, int skip)
"K_load(*<FILE>i)*<ok_kernel>",
# bool K_addDataFromSystem(ok_kernel* k, const char* filename)
"K_addDataFromSystem(pZ)B",
# void K_calculate(ok_kernel* k)
"K_calculate(p)v",
# int K_minimize(ok_kernel* k, int algo, int maxiter, double params[])
"K_minimize(pii*d)i",
# int K_1dminimize(ok_kernel* k, int algo, int maxiter, int row, int column, double params[])
"K_1dminimize(piiii*d)i",
# ok_system** K_integrate(ok_kernel* k, gsl_vector* times, ok_system** bag, int* error)
"K_integrate(p*<gsl_vector>p*i)p",
# ok_system** K_integrateRange(ok_kernel* k, double from, double to, unsigned int samples, ok_system** bag, int* error)
"K_integrateRange(pddIp*i)p",
# gsl_matrix* K_integrateStellarVelocity(ok_kernel* k, double from, double to, unsigned int samples, ok_system** bag, int* error)
"K_integrateStellarVelocity(pddIp*i)*<gsl_matrix>",
# ok_system** K_integrateProgress(ok_kernel* k, gsl_vector* times, ok_system** bag, int* error)
"K_integrateProgress(p*<gsl_vector>p*i)p",
# void K_print(ok_kernel* k, FILE* f)
"K_print(p*<FILE>)v",
# void K_save_old(ok_kernel* k, const char* stem)
"K_save_old(pZ)v",
# void K_setSeed(ok_kernel* k, unsigned long int seed)
"K_setSeed(pL)v",
# void* ok_bridge_kernel_buf(void* buf, int n, ok_kernel* k)
"ok_bridge_kernel_buf(pip)p",
# ok_system* ok_alloc_system(int nplanets)
"ok_alloc_system(i)*<ok_system>",
# void ok_free_system(ok_system* system)
"ok_free_system(*<ok_system>)v",
# void ok_free_systems(ok_system** system, unsigned int len)
"ok_free_systems(pI)v",
# ok_system* ok_copy_system(const ok_system* orig)
"ok_copy_system(*<ok_system>)*<ok_system>",
# void ok_resize_system(ok_system* system, int npnew)
"ok_resize_system(*<ok_system>i)v",
# void ok_el2cart(ok_system* system, gsl_matrix* xyz)
"ok_el2cart(*<ok_system>*<gsl_matrix>)v",
# void ok_cart2el(ok_system* system, gsl_matrix* els, bool internal)
"ok_cart2el(*<ok_system>*<gsl_matrix>B)v",
# void ok_setup(ok_system* system)
"ok_setup(*<ok_system>)v",
# int ok_force(double t, const double y[], double f[], void *params)
"ok_force(d*d*dv)i",
# int ok_jac(double t, const double y[], double *dfdy,     double dfdt[], void *params)
"ok_jac(d*dd*dv)i",
# ok_system** ok_integrate(ok_system* initial, const gsl_vector* times, ok_integrator_options* options, const int integrator,         ok_system** bag, int* error)
"ok_integrate(*<ok_system>*<gsl_vector>*<ok_integrator_options>ip*i)p",
# double ok_get_rv(ok_system* sys)
"ok_get_rv(*<ok_system>)d",
# gsl_matrix* ok_get_rvs(ok_system** sys, int len)
"ok_get_rvs(pi)*<gsl_matrix>",
# gsl_matrix* ok_get_xyzs(ok_system** bag, int len)
"ok_get_xyzs(pi)*<gsl_matrix>",
# gsl_matrix* ok_get_els(ok_system** bag, int len, bool internal)
"ok_get_els(piB)*<gsl_matrix>",
# gsl_vector* ok_find_transits(ok_system** bag, const int len, const int pidx, const int intMethod, const double eps, const int flags[], int* error)
"ok_find_transits(piiid*i*i)*<gsl_vector>",
# void ok_to_cm(ok_system* system, gsl_matrix* xyz)
"ok_to_cm(*<ok_system>*<gsl_matrix>)v",
# void ok_to_star(ok_system* system, gsl_matrix* xyz)
"ok_to_star(*<ok_system>*<gsl_matrix>)v",
# int ok_find_closest_transit(ok_system* sys, const int pidx, ok_integrator_options* options, const int intMethod, const double eps, const int type, double* timeout, int* error)
"ok_find_closest_transit(*<ok_system>i*<ok_integrator_options>idi*d*i)i",
# double ok_pcalc(const double a, const double Mcenter, const double Mp)
"ok_pcalc(ddd)d",
# double ok_acalc(const double P, const double Mcenter, const double Mp)
"ok_acalc(ddd)d",
# ok_list* KL_alloc(const int size, ok_kernel* prototype)
"KL_alloc(ip)*<ok_list>",
# void KL_free(ok_list* list)
"KL_free(*<ok_list>)v",
# ok_list* KL_load(FILE* fid, int skip)
"KL_load(*<FILE>i)*<ok_list>",
# void KL_save(const ok_list* kl, FILE* out)
"KL_save(*<ok_list>*<FILE>)v",
# void KL_append(ok_list* dest, ok_list* src)
"KL_append(*<ok_list>*<ok_list>)v",
# gsl_vector* KL_getParsStats(const ok_list* kl, const int what)
"KL_getParsStats(*<ok_list>i)*<gsl_vector>",
# gsl_vector* KL_getElements(const ok_list* kl, const int pl, const int el)
"KL_getElements(*<ok_list>ii)*<gsl_vector>",
# gsl_vector* KL_getPars(const ok_list* kl, const int vo)
"KL_getPars(*<ok_list>i)*<gsl_vector>",
# gsl_matrix* KL_getElementsStats(const ok_list* kl, const int what)
"KL_getElementsStats(*<ok_list>i)*<gsl_matrix>",
# ok_list_item* KL_set(ok_list* kl, const int idx, gsl_matrix* elements, gsl_vector* pars, double merit, int tag)
"KL_set(*<ok_list>i*<gsl_matrix>*<gsl_vector>di)*<ok_list_item>",
# int KL_getSize(const ok_list* kl)
"KL_getSize(*<ok_list>)i",
# void KL_removeAtIndex(ok_list* kl, const int idx)
"KL_removeAtIndex(*<ok_list>i)v",
# void KL_fprintf(const ok_list* kl, FILE* out, const char* fmt, const char* lfmt)
"KL_fprintf(*<ok_list>*<FILE>ZZ)v",
# void KL_to_ptr(const ok_list* kl, double* out)
"KL_to_ptr(*<ok_list>*d)v",
# int KL_getNplanets(const ok_list* kl)
"KL_getNplanets(*<ok_list>)i",
# ok_list* K_mcmc_single(ok_kernel* k, unsigned int nsteps, unsigned int skip, unsigned int discard, const double dparams[], ok_list* cont, ok_callback2 merit_function, int tag, int* flag)
"K_mcmc_single(pIII*d*<ok_list>pi*i)*<ok_list>",
# ok_list* K_mcmc_mult(ok_kernel** k, unsigned int nchains, unsigned int ntemps, unsigned int skip, unsigned int discard, const double params[], double Rstop, ok_callback2 merit_function)
"K_mcmc_mult(pIIII*ddp)*<ok_list>",
# ok_list* K_bootstrap(ok_kernel* k, int trials, int warmup, int malgo, int miter, double mparams[])
"K_bootstrap(piiii*d)*<ok_list>",
# gsl_matrix* ok_periodogram_ls(const gsl_matrix* data, const unsigned int samples, const double Pmin, const double Pmax, const int method,         unsigned int timecol, unsigned int valcol, unsigned int sigcol, ok_periodogram_workspace* p)
"ok_periodogram_ls(*<gsl_matrix>IddiIII*<ok_periodogram_workspace>)*<gsl_matrix>",
# gsl_matrix* ok_periodogram_boot(const gsl_matrix* data, const unsigned int trials, const unsigned int samples,         const double Pmin, const double Pmax, const int method,         const unsigned int timecol, const unsigned int valcol, const unsigned int sigcol,         const unsigned long int seed, ok_periodogram_workspace* p, ok_progress prog)
"ok_periodogram_boot(*<gsl_matrix>IIddiIIIL*<ok_periodogram_workspace>p)*<gsl_matrix>",
# gsl_matrix* ok_periodogram_full(ok_kernel* k, int type, int algo, bool circular, unsigned int sample,         const unsigned int samples, const double Pmin, const double Pmax)
"ok_periodogram_full(piiBIIdd)*<gsl_matrix>",
# int K_isMstable_coplanar(const gsl_matrix* alle)
"K_isMstable_coplanar(*<gsl_matrix>)i",
# double K_crossval_l1o(ok_kernel* k, int minalgo, int maxiter, double params[])
"K_crossval_l1o(pii*d)d",
# int mco_el2x__(doublereal mu, doublereal q, doublereal e,                doublereal i__, doublereal p, doublereal n, doublereal l,                doublereal* x, doublereal* y, doublereal* z__, doublereal* u,                doublereal* v, doublereal* w)
"mco_el2x__(ddddddd*d*d*d*d*d*d)i",
# int mco_x2el__(doublereal* mu, doublereal* x, doublereal* y,                doublereal* z__, doublereal* u, doublereal* v, doublereal* w,                doublereal* q, doublereal* e, doublereal* i__, doublereal* p,                doublereal* n, doublereal* l)
"mco_x2el__(*d*d*d*d*d*d*d*d*d*d*d*d*d)i",
# doublereal mco_kep__(doublereal e, doublereal oldl)
"mco_kep__(dd)d",
# gsl_matrix* gsl_matrix_alloc(size_t n1, size_t n2)
"gsl_matrix_alloc(LL)*<gsl_matrix>",
# gsl_matrix* gsl_matrix_calloc(size_t n1, size_t n2)
"gsl_matrix_calloc(LL)*<gsl_matrix>",
# void gsl_matrix_free(gsl_matrix* m)
"gsl_matrix_free(*<gsl_matrix>)v",
# void gsl_vector_free(gsl_vector* v)
"gsl_vector_free(*<gsl_vector>)v",
# int gsl_matrix_fwrite(FILE* stream, const gsl_matrix* m)
"gsl_matrix_fwrite(*<FILE>*<gsl_matrix>)i",
# int gsl_matrix_fread(FILE* stream, gsl_matrix* m)
"gsl_matrix_fread(*<FILE>*<gsl_matrix>)i",
# int gsl_matrix_fprintf(FILE* stream, const gsl_matrix* m, const char* format)
"gsl_matrix_fprintf(*<FILE>*<gsl_matrix>Z)i",
# int gsl_matrix_fscanf(FILE* stream, gsl_matrix* m)
"gsl_matrix_fscanf(*<FILE>*<gsl_matrix>)i",
# double gsl_matrix_get(const gsl_matrix* m, size_t i, size_t j)
"gsl_matrix_get(*<gsl_matrix>LL)d",
# void gsl_matrix_set(gsl_matrix* m, size_t i, size_t j, double x)
"gsl_matrix_set(*<gsl_matrix>LLd)v",
# void free(void* v)
"free(p)v",
# FILE* fopen(const char* filename, const char* mode)
"fopen(ZZ)*<FILE>",
# int fclose(FILE* stream)
"fclose(*<FILE>)i",
# size_t fwrite(const void* ptr, size_t size, size_t count, FILE* stream)
"fwrite(pLL*<FILE>)L",
# int fputs(const char* str, FILE* stream)
"fputs(Z*<FILE>)i",
""))

TIME <- K_T_TIME + 1
VAL <- K_T_VAL + 1
ERR <- K_T_ERR + 1
PRED <- K_T_PRED + 1
SVAL <- K_T_SVAL + 1
FLAG <- K_T_FLAG + 1
SET <- K_T_SET + 1
TDS_PLANET <- K_T_TDS_PLANET + 1
TDS_FLAG <- K_T_TDS_FLAG + 1

DATA1 <- K_P_DATA1 + 1
DATA2 <- K_P_DATA2 + 1
DATA3 <- K_P_DATA3 + 1
DATA4 <- K_P_DATA4 + 1
DATA5 <- K_P_DATA5 + 1
DATA6 <- K_P_DATA6 + 1
DATA7 <- K_P_DATA7 + 1
DATA8 <- K_P_DATA8 + 1
DATA9 <- K_P_DATA9 + 1
DATA10 <- K_P_DATA10 + 1

DATA.NOISE1 <- K_P_DATA_NOISE1 + 1
DATA.NOISE2 <- K_P_DATA_NOISE2 + 1
DATA.NOISE3 <- K_P_DATA_NOISE3 + 1
DATA.NOISE4 <- K_P_DATA_NOISE4 + 1
DATA.NOISE5 <- K_P_DATA_NOISE5 + 1
DATA.NOISE6 <- K_P_DATA_NOISE6 + 1
DATA.NOISE6 <- K_P_DATA_NOISE6 + 1
DATA.NOISE7 <- K_P_DATA_NOISE7 + 1
DATA.NOISE8 <- K_P_DATA_NOISE8 + 1
DATA.NOISE9 <- K_P_DATA_NOISE9 + 1
DATA.NOISE10 <- K_P_DATA_NOISE10 + 1

RV.TREND <- K_P_RV_TREND + 1

RV <- K_T_RV 
TIMING <- K_T_TIMING

TDS_PRIMARY <- K_TDS_PRIMARY
TDS_SECONDARY <- K_TDS_SECONDARY

PER <- K_PER + 1
MASS <- K_MASS + 1
MA <- K_MA + 1
ECC <- K_ECC + 1
LOP <- K_LOP + 1
INC <- K_INC + 1
NODE <- K_NODE + 1
RADIUS <- K_RADIUS + 1
ORD <- K_ORD + 1
SMA <- K_SMA + 1
SEMIAMP <- K_SEMIAMP + 1
TPERI <- K_TPERI + 1
TRUEANOMALY <- K_TRUEANOMALY + 1
DONE <- K_DONE

TDS_PLANET <- K_T_TDS_PLANET + 1

ELEMENTS_SIZE <- K_ELEMENTS_SIZE
ALL_ELEMENTS_SIZE <- K_ALL_ELEMENTS_SIZE
DATA_SIZE <- K_DATA_SIZE
PARAMS_SIZE <- K_PARAMS_SIZE
DATA_SETS_SIZE <- K_DATA_SETS_SIZE

MINIMIZE <- K_MINIMIZE
ACTIVE <- K_ACTIVE
INACTIVE <- 0

SIMPLEX <- K_SIMPLEX
LM <- K_LM
SA <- K_SA
DIFFEVOL <- K_DIFFEVOL

ASTROCENTRIC <- K_ASTROCENTRIC
JACOBI <- K_JACOBI

T_RV <- K_T_RV
T_TIMING <- K_T_TIMING
T_PHOTO <- K_T_PHOTO
T_DUMMY <- K_T_DUMMY

PS_TIME <- K_PS_TIME + 1
PS_Z <- K_PS_Z + 1
PS_FAP <- K_PS_FAP + 1
PS_Z_LS <- K_PS_Z_LS + 1
PS_TAU <- K_PS_TAU + 1
PS_WIN <- K_PS_WIN + 1

KEPLER <- K_KEPLER
RK45 <- K_RK45
RK89 <- K_RK89
ADAMS <- K_ADAMS
SWIFTRMVS <- K_SWIFTRMVS
BULIRSCHSTOER <- K_BULIRSCHSTOER

AU <- K_AU
MSUN <- K_MSUN
MJUP <- K_MJUP
DAY <- K_DAY
YEAR <- K_YEAR
GGRAV <- K_GGRAV
K2 <- ((GGRAV * MSUN * DAY * DAY) / (AU*AU*AU))

SYSTEMIC.VERSION <- K_SYSTEMIC_VERSION

# Labels

.elements <- c("period", "mass", "ma", "ecc", "lop", "inc", "node", "radius", "ord", "u1", "u2", "u3", "u4")
.allelements <- c(.elements, "a", "k", "tperi", "trueanomaly", "meanlongitude", "j1", "j2")
.params <- c(sprintf("par%d", 1:PARAMS_SIZE))
.params[1:DATA_SETS_SIZE] <- sprintf("data%d", 1:DATA_SETS_SIZE)
.params[(DATA_SETS_SIZE+1):(2*DATA_SETS_SIZE)] <- sprintf("data.noise%d", 1:DATA_SETS_SIZE)
.params[RV.TREND] <- "rv.trend"

.data <- sprintf("V%d", 1:DATA_SIZE)
.data[TIME] <- 'time'
.data[VAL] <- 'val'
.data[ERR] <- 'err'
.data[PRED] <- 'pred'
.data[SVAL] <- 'sval'
.data[FLAG] <- 'flag'
.data[SET] <- 'set'

.periodogram <- c("period", "power", "fap", "ls_power", "tau", "window")
.elements.labels <- c(period="Period", mass="Mass", ma="Mean anomaly", ecc="Eccentricity", lop="Longitude of periastron", 
	inc="Inclination", node="Node", radius="Radius", ord="Label", u1="u1", u2="u2", u3="u3", u4="u4",
	a="Semi-major axis", k="Semiamplitude", tperi="Time of periastron passage", trueanomaly="True anomaly", 
	meanlongitude="Mean longitude", j1="j1", j2="j2")

.integration.errors <- c() 
.integration.errors[K_INTEGRATION_SUCCESS+1] <- "Integration was successful"
.integration.errors[K_INTEGRATION_FAILURE_SMALL_TIMESTEP+1] <- "Timestep decreased too much or close encounter"
.integration.errors[K_INTEGRATION_FAILURE_INCREASE_TOLERANCE+1] <- "Increase tolerance to achieve convergence"
.integration.errors[K_INTEGRATION_FAILURE_CLOSE_ENCOUNTER+1] <- "Close encounter between two planets"
.integration.errors[K_INTEGRATION_FAILURE_CLOSE_ENCOUNTER_STAR+1] <- "Close encounter with the star"
.integration.errors[K_INTEGRATION_FAILURE_STOPPED+1] <- "Integration stopped"
.integration.errors[K_INTEGRATION_FAILURE_SWIFT+1] <- "Error during SWIFT integration (maybe planet was lost?)"

