
options(systemic.auto=FALSE)
options(systemic.psamples=5e4)
options(systemic.pmin=0.5)
options(systemic.pmax=1e4)

.properties <- list(
    nplanets = K_getNplanets,
    chi2 = K_getChi2,
    rms = K_getRms,
    chi2nr = K_getChi2_nr,
    loglik = K_getLoglik,
    jitter = K_getJitter,
    mstar = K_getMstar,
    ndata = K_getNdata,
    epoch = K_getEpoch,
    nplanets = K_getNplanets,
    int.method = K_getIntMethod,
    ntts = K_getNtts,
    nrvs = K_getNrvs,
    nsets = K_getNsets,
    chi2rvs = K_getChi2_rvs,
    chi2tts = K_getChi2_tts,
    rmstts = K_getRms_tts,
    nrpars = K_getNrPars,
    element.type = K_getElementType,
    abs.acc = K_getIntAbsAcc,
    rel.acc = K_getIntRelAcc,	
    dt = K_getIntDt,
    trange = function(h) {
        a <- NaN
        b <- NaN

        K_getRange(h, a, b)
        return(c(a, b))
    },
    ks.pvalue = function(h) {
        nd <- K_getNdata(h)
        if (nd <= 0)
            return(NaN)
        m <- .buf_to_R(K_compileData(h), K_getNdata(h), DATA_SIZE)
        p <- .gsl_vector_to_R(K_getPars(h))
        n <- p[m[, SET]+1+10]^2
        err <- sqrt(m[, ERR]^2 + n)
        return (suppressWarnings(ks.test((m[, PRED]-m[, SVAL])/err, 'pnorm')$p.value))
    }
)



systemic.names <- c(period='Period', mass='Mass', ma='Mean anomaly', ecc='Eccentricity',
                   lop='Longitude of pericenter', inc='Inclination', node='Node',
                   a='Semi-major axis', k='Semiamplitude', tperi='Periastron passage time', rv.trend='Linear trend',
                   rv.trend.quadratic='Quadratic trend', mstar='Stellar mass',
                   chi2='Reduced Chi-square', jitter='Stellar jitter', rms='RMS', epoch='Epoch', ndata='Data points', trange='Span of observations', data.noise1='Noise for dataset 1', data.noise2='Noise for dataset 2')
systemic.units <- c(period='[days]', mass='[M_{jup}]', ma='[deg]', ecc='',
                   lop='[deg]', inc='[deg]', node='[deg]',
                   a='[AU]', k='[m/s]', tperi='[JD]',
                   rv.trend='[m/s]', rv.trend.quadratic='[m/s^2]',
                   mstar = '[M_{sun}]', chi2='', jitter='[m/s]', rms='[m/s]',
                   epoch = '[JD]', ndata='', trange='[JD]', data.noise1='[m/s]', data.noise2='[m/s]')


ELEMENT <- 0
PARAMETER <- 1
PROPERTY <- 2

systemic.type <- rep(ELEMENT, times=length(systemic.names))
names(systemic.type) <- names(systemic.units)
systemic.type[names(systemic.type) %in% .params] <- PARAMETER
systemic.type[names(systemic.type) %in% names(.properties)] <- PROPERTY

.free.gsl_matrix <- function(m) {
    gsl_matrix_free(m)
}

.free.gsl_vector <- function(v) {
    gsl_vector_free(v)
}

.free.kernel <- function(h) {
    K_free(h)
}

.stop.ifnot <- function(b, msg) {
    if (!b)
        stop(msg)
}

# Used internally to convert from a gsl_matrix (double) to an R matrix
.gsl_matrix_to_R <- function(m, free = FALSE, .keep.h = FALSE) {
    dm <- c(ok_matrix_rows(m), ok_matrix_cols(m))
    len <- dm[1] * dm[2]

    v <- numeric(len)
    ok_block_to_ptr(ok_matrix_block(m), v)
    
    if (free) { 
        reg.finalizer(m, .free.gsl_matrix)
    }
    
    mat <- matrix(v, nrow=dm[1], ncol=dm[2], byrow=T)
    if (.keep.h)
        attr(mat, 'h') <- m
    
    return(mat)
}

.check_kernel <- function(k) {
    stopifnot(class(k) == "kernel")
    if (is.nullptr(k$h)) {
        stop("This kernel object is invalid (null pointer). Maybe you saved it using save() instead of save.systemic() ?")
    }
}

.gsl_vector_to_R <- function(v, free = F, .keep.h = F) {
    len <- ok_vector_len(v)
    v2 <- numeric(len)
    ok_block_to_ptr(ok_vector_block(v), v2)
    
    if (.keep.h)
        attr(v2, 'h') <- v
    
    if (free) { 
        reg.finalizer(v, .free.gsl_vector)
    }
    return(v2)
}


# Used internally to convert a double** matrix (with layout d[ROWS][COLUMNS]) to an R matrix
.buf_to_R <- function(b, r, c) {
    v <- numeric(r*c)
    ok_buf_to_ptr(b, r, c, v)
    return(matrix(v, nrow=r, ncol=c, byrow=T))	
}

# Used internally to convert from an R matrix to a gsl_matrix (double)
.R_to_gsl_matrix <- function(m, gc = F) {
    h <- ok_ptr_to_matrix(as.numeric(t(m)), nrow(m), ncol(m))
    if (gc) reg.finalizer(h, .free.gsl_matrix)
    return(h)
}

.label_to_index <- function(v) {
    if (is.null(v) || is.na(v))
        stop("Column is NULL or NA")
    if (is.numeric(v) || v == "all") return(v)

    a <- match(v, .allelements)
    if (is.na(a))
        a <- match(v, .params)
    if (is.na(a))
        stop(paste("Unrecognized label ", v, ", allowed: ", paste(.elements, collapse=' '), paste(.params, collapse=' ')))
    else
        return(a)
}

.mute <- function(k) {
    if (! is.null(k[['auto_old']]))
        return()
    k[['auto_old']] <- k[['auto']]
    k[['auto']] <- FALSE
}

.unmute <- function(k) {
    if (is.null(k[['auto_old']]))
        return()
    k[['auto']] <- k[['auto_old']]
    k[['auto_old']] <- NULL
}



knew <- function() {
    ## Creates a new kernel object. [1]
    #
    # Returns:
    # A new kernel object. A kernel object contains both the fit parameters and the data.
    #
    # Example:
    # myfit <- knew()
    # myfit[] <- c(period = 300, mass = 1, ma = 0, ecc = 0.1, lop = 30)
    
    k <- new.env()
    k$h <- K_alloc()
    k$auto <- getOption("systemic.auto", F)
    k$min.method <- SIMPLEX
    k$last.error.code <- integer(1)
    k$.new <- TRUE
    k$min.func <- "chi2"
    class(k) <- "kernel"
    k$epoch <- 2450000
    k$.epoch.set <- FALSE
    
    reg.finalizer(k$h, .free.kernel)
    kupdate(k)
    return(k)
}

kadd.planet <- function(k, period = 300, mass = 1, ma = 0, ecc = 0, lop = 0,
                       inc = 90, node = 0, tperi = NA, K = NA, a = NA) {
		
    ## Adds a new body with the given elements (you can specify period, mass, ecc, lop, inc, node) [2]
    # 
    # Alternative syntax:
    # k[] <- c(period = 300, mass = 1, ...)
    #
    # Args:
    # - k: the kernel to add the planet to
    # - period: period of planet (days)
    # - mass: mass of planet (Jupiter masses)
    # - ma: mean anomaly (deg)
    # - ecc: eccentricity
    # - lop: longitude of pericenter (deg)
    # - inc: inclination (deg)
    # - node: node (deg)
    # - tperi: time of passage through pericenter (days)
    # - K: semi-amplitude (m/s)
    # - a: semi-major axis (AU)
		#
    # Example:
    #
    # kadd.planet(k, period=300, mass=1, ecc=0.5, lop=100)
    #
    # equivalently:
    #
    # k[] <- c(period = 300, mass = 1, ecc = 0.5, lop = 100)
    
    .check_kernel(k)
    els <- c(PER-1, period, MASS-1, mass, MA-1, ma, ECC-1, ecc, LOP-1, lop,
             INC-1, inc, NODE-1, node, DONE)
    K_addPlanet(k$h, els)
    np <- K_getNplanets(k$h)
    K_setElementRange(k$h, np, K_ECC, 0, 0.95)
    K_setElementRange(k$h, np, K_MASS, 1e-4, 1e3)

    if (!is.na(tperi))
        k[np, 'tperi'] <- tperi
    if (!is.na(K))
        k[np, 'k'] <- K
    if (!is.na(a))
        k[np, 'a'] <- a
    
    on.exit(if (k$auto) kupdate(k))
    return()
}

kremove.planet <- function(k, idx) {
    ## Removes the idx-th planet [2]
    #
    # Alternative syntax: k[idx] <- NULL
    #
    # Args:
    # - k: kernel to remove the planet from
    # - idx: index of the planet to remove (starting at 1)
    on.exit(if (k$auto) kupdate(k))
    .check_kernel(k)
    if (idx == 'all') {
        if (k$nplanets > 0) {
            for (i in 1:k$nplanets) {
                kremove.planet(k, 1)
            }
        }
        return()
    }
    stopifnot(idx >= 1 && idx <= k$nplanets)

    K_removePlanet(k$h, idx)
    K_compileData(k$h)
   
}


kels <- function(k, keep.first = FALSE) {
    ## Returns a matrix of orbital elements [2]
    #
    # Alternative syntax: k[]
    # 
    # Args: 
    #  k: kernel to read the orbital elements from
    #
    # Returns:
    #  A matrix of orbital elements; each row is a planet's orbital elements.
    #
    # Example:
    #
    # elements <- kels(k)
    # print(elements)
    .check_kernel(k)
    m <- .gsl_matrix_to_R(K_getElements(k$h))
    if (! keep.first) {
        m = m[-1,, drop=FALSE]
    }
    colnames(m) <- .elements
    return(m)
}

`kels<-` <- function(k, value) {
    ## Sets the orbital elements of the kernel [2]
    #
    # Args:
    # - k: kernel where the orbital elements should be set
    # - value: a matrix of orbital elements
    #
    # Example:
    # elements <- kels(k)
    # elements[, PER] <- 2 * elements[, PER]  # Multiply all periods by two
    # kels(k) <- elements
    elements <- value
    .check_kernel(k)	
    stopifnot(class(elements) == "matrix",  nrow(value) >= 1)
    
    if (ncol(elements) > ELEMENTS_SIZE)
        elements <- value[,1:ELEMENTS_SIZE,drop=F]
    
    els <- kels(k, keep.first = T)
    elements <- rbind(els[1,], value)
    K_setElements(k$h, .R_to_gsl_matrix(elements, gc=F))
    if (k$auto) kupdate(k)
    return(k)
}


`[.kernel` <- function(k, idx1=NA, idx2=NA) {    
    ## Subsetting a kernel object with brackets returns the values of the orbital elements. [2] 
    # * k[] returns a matrix of elements (like kels)
    # * k[n, ] returns the elements for the n-th planet
    # * k[, m] returns the m-th element for all planets
    # * k[n, m] returns the m-th element for the n-th planet
    # * k['par'] returns a vector of all the parameters (e.g. offsets, linear trend, etc.)
    # * k['par', j] returns the j-th parameter
    # FUN: kernel[idx1, idx2]
    .check_kernel(k)	
    if (is.na(idx1) && is.na(idx2)) {
        return(kallels(k))
    } else if (is.na(idx2)) {
        
        if (idx1 == "par" || idx1 == -1) {
            return(kpars(k))
        }

        if (idx1 == "minimized") {
            return(apply(kminimized.indices(k), 2, function(l) k[l[1], l[2]]))
        }
        
        stopifnot(idx1 <= k$nplanets);
        return(kallels(k)[idx1,,drop=FALSE])
    } else if (is.na(idx1)) {
        idx2 = .label_to_index(idx2)
        return(kallels(k)[,idx2,drop=FALSE])
    } else {
        idx2 = .label_to_index(idx2)				
        if (idx1 == "par" || idx1 == -1) {
            stopifnot(idx2 <= PARAMS_SIZE)
            return(K_getPar(k$h, idx2-1))
        } else {
            stopifnot(idx1 <= k$nplanets);		
            return(K_getElement(k$h, idx1, idx2-1))
        }
    }
}

`[<-.kernel` <- function(k, idx1=NA, idx2=NA, value) {
    ## The parameters of a kernel object can be set by subsetting it with brackets and assigning values.  [2]
    #
    # - k[] <- matrix sets the orbital elements to the specified matrix
    # - k[] <- c(period=..., mass=..., ...) adds a new planet with the specified orbital elements 
    # - k[n] <- NULL removes the n-th planet
    # - k[n, m] <- v sets the m-th element for the n-th planet to v
    # - k[, m] <- v sets the m-th element to v for all planets
    # - k['par'] <- v sets parameter values to the vector v
    # - k['par', n] <- v sets the n-th parameter to v
    #
    # FUN: k[idx1, idx2] <- value
    
    .check_kernel(k)	
    if (is.na(idx1) && is.na(idx2)) {
        if (class(value) == "matrix") {
            kels(k) <- value
            return(k)
        } else {
            if (is.null(value)) {
                kremove.planet(k, k$nplanets)
                return(k)	
            }
            do.call(kadd.planet, append(list(k), as.list(value)))
            return(k)
        }
    } else if (! is.na(idx1) && idx1 < 0) {
        .mute(k)
        for (i in 1:k$nplanets) {
            if (i != -idx1) k[i, idx2] <- value
        }
        .unmute(k)		
        if (k$auto) kupdate(k)
        return(k)
    } else if (is.na(idx2)) {
        if (is.null(value) || is.na(value)) {
            kremove.planet(k, idx1)		
        } else if (idx1 == "par" || idx1 == -1) {
            if (length(value) == 1) {
                K_setPar(k$h, -1, value)	
            } else {
                for (i in 1:length(value))
                    K_setPar(k$h, i-1, value[[i]])
            }
        } else if (idx1 == "minimized") {
            nminimized.pars <- ncol(kminimized.indices(k))
            stopifnot(nminimized.pars == length(value))
            K_setMinimizedValues(k$h, value)
        } else {
            for (i in PER:RADIUS)
                K_setElement(k$h, idx1, i-1, value[[i]])
            
        }
        if (k$auto) kupdate(k)
        return(k)

    } else if (is.na(idx1)) {
        idx2 = .label_to_index(idx2)

        if (length(value) == 1) {
            K_setElement(k$h, -1, idx2-1, value)
        } else {
            for (i in 1:min(length(value), k$nplanets))
                K_setElement(k$h, i, idx2-1, value[[i]])
        }

        if (k$auto) kupdate(k)
        return(k)
    } else {
        if (idx1 == "par") {
            idx2 = .label_to_index(idx2)
            stopifnot(idx2 <= PARAMS_SIZE)
            K_setPar(k$h, idx2 - 1, value)
        } else {
            stopifnot(idx1 <= k$nplanets);		
            idx2 = .label_to_index(idx2)		
            K_setElement(k$h, idx1, idx2-1, value)
        }
        if (k$auto) kupdate(k)
        return (k)
    }
}



`$.kernel` <- function(k, idx) {
    ## Use the $ operator to access the following properties of the kernel. [6]
    #
    # Read-only properties:
    # * k$nplanets	Number of planets
    # * k$chi2 		Current reduced chi^2 value (normalized by (k$ndata - k$nrpars))
    # * k$chi2nr		Non-reduced chi^2
    # * k$rms		Current RMS value
    # * k$jitter	Current jitter value
    # * k$loglik	Current log likelihood (multiplied by -1)
    # * k$bic Current value of BIC (multiplied by -1)
    # * k$ks.pvalue	Current p-value of the KS test comparing normalized residuals to a unit gaussian
    # * k$ndata		Number of data points
    # * k$nrvs		Number of RV data points
    # * k$ntts		Number of central transits
    # * k$nsets		Number of data sets loaded
    # * k$chi2rvs	Chi^2 (RVs only)
    # * k$chi2tts	Chi^2 (central transits only)
    # * k$trange	Time range of the compiled dataset
    # * k$epoch Epoch of the fit (JD)
    # * k$mstar Mass of the star (Msun)
    # * k$int.method	Integration method (possible values: KEPLER, RK45, RK89)
    # * k$element.type Coordinate type (possible values: ASTROCENTRIC, JACOBI)
    # * k$min.func Function minimized by @kminimize. Possible values are "chi2" (default), "rms", or a function that takes a kernel and returns a number.
    # * k$nrpars	"Degrees of freedom" parameter used to calculate reduced chi^2. It is equal to the number of all the parameters that are marked as ACTIVE or MINIMIZE
    # FUN: k$property
    if (! is.null(.properties[[idx]])) {
        .check_kernel(k)			
        return(.properties[[idx]](k$h))
    } else if (idx == "custom.model") {
        return(k[["custom.model.function"]])
    } else if (idx == "last.error") {
        return(.integration.errors[[k$last.error.code+1]])	
    } else if (idx == "min.func") {
        return(k[['min.func']])
    } else if (idx == 'bic') {
        return(2*k$loglik - k$nrpars*(log(k$ndata)))
    } else if (idx == 'datanames') {
        if (k$nsets == 0) return(c())
        c <- sprintf("data-%d", 1:k$nsets)
        for (i in 1:k$nsets)
            if (K_infoExists(k$h, sprintf("DataFileName%d", i-1)))
                c[i] <- K_getInfo(k$h, sprintf("DataFileName%d", i-1))
        return(c)
    } else if (idx == 'filename') {
        if (K_infoExists(k$h, "FileName"))
            return(K_getInfo(k$h, "FileName"))
        else
            return("")
    } else
        return(k[[idx]])	
}


kget <- function(k, idx) {
    return(`$.kernel`(k, idx))
}

kprop <- function(k, idx = 'all') {    
    idx <- as.character(idx)

    if (idx != 'all') {
        if (K_infoExists(k$h, idx))
            return(K_getInfo(k$h, idx))
        else
            return(NULL)
    } else {
        i <- 0
        s <- K_getInfoTag(k$h, i)
        out <- c()
        while (s != "") {
            out[s] <- K_getInfo(k$h, s)
            i <- i +1
            s <- K_getInfoTag(k$h, i)
        }
        return(out)
    }
}

`kprop<-` <- function(k, idx, value) {
    idx <- as.character(idx)
    value <- as.character(value)
    K_setInfo(k$h, idx, value)
    return(k)
}

`$<-.kernel` <- function(k, idx, value) {
    ## Use the $ operator to set the following properties of the kernel. [6]
    #
    # Settable properties:
    # * k$int.method	Integration method (possible values: KEPLER, RK45, RK89)
    # * k$element.type Coordinate type (possible values: ASTROCENTRIC, JACOBI)
    # * k$epoch		Epoch in JD
    # * k$mstar		Mass of the star in solar masses
    # * k$minimize.func	Function or property string (one of chi2, rms, jitter, chi2nr, chi2rvs, chi2tts, loglik) which specifies the value to minimize
    # * k$min.func Function minimized by @kminimize. Possible values are "chi2" (default), "rms", or a function that takes a kernel and returns a number.
    # FUN: k$property <- value

    .check_kernel(k)	
    if (idx == "mstar") {
        K_setMstar(k$h, value)
        if (k$auto) kupdate(k)
    } else if (idx == "epoch") {
        k$.epoch.set <- TRUE
        K_setEpoch(k$h, value)
        if (k$auto) kupdate(k)	
    } else if (idx == "abs.acc") {
        K_setIntAbsAcc(k$h, value)	
        if (k$auto) kupdate(k)			
    } else if (idx == "rel.acc") {
        K_setIntRelAcc(k$h, value)	
        if (k$auto) kupdate(k)			
    } else if (idx == "dt") {
        K_setIntDt(k$h, value)	
        if (k$auto) kupdate(k)			
    } else if (idx == "int.method") {
        K_setIntMethod(k$h, value)
        if (k$auto) kupdate(k)
    } else if (idx == "min.method") {
        k[['min.method']] <- value
        if (k$auto) kupdate(k, calculate=FALSE)
    } else if (idx == "element.type") {
        if (value == "astrocentric")
            value = ASTROCENTRIC
        else if (value == "jacobi")
            value = JACOBI
        
        K_setElementType(k$h, value)	
        if (k$auto) kupdate(k)
    } else if (idx == "custom.model") {
        if (is.null(value)) {
            K_setCustomModelFunction(k[['h']], NULL)
            k[['custom.model.function']] <- NULL
            k[['custom.model.function.cb']] <- NULL
            print("Null value")
        } else {
            stopifnot(class(value) == "function")
            k[['custom.model.function']] <- value
            k[['custom.model.function.cb']] <- new.callback("ppi)v", 
                                                            function(khandle, bufhandle, ndata) {
                                                                d <- .buf_to_R(bufhandle, ndata, DATA_SIZE)
                                                                k <- new.env()
                                                                k[['h']] <- khandle
                                                                class(k) <- "kernel"
                                                                
                                                                times <- numeric(ndata)
                                                                flags <- numeric(ndata)					
                                                                ok_buf_col(bufhandle, times, K_T_TIME, ndata)
                                                                ok_buf_col(bufhandle, flags, K_T_FLAG, ndata)					
                                                                ret <- value(k, times, flags)

                                                                stopifnot(length(ret) == ndata)
                                                                ok_buf_add_to_col(bufhandle, ret, K_T_PRED, ndata)
                                                                d <- .buf_to_R(bufhandle, ndata, DATA_SIZE)
                                                            })
            K_setCustomModelFunction(k[['h']], k[['custom.model.function.cb']])
        }
        if (k$auto) kupdate(k)
    } else if (idx == "min.func") {
        if (is.character(value)) {
            if (is.null(.properties[[value]])) stop(sprintf("Could not find property %s", value))
            if (value == "chi2")
                K_setMinFunc(k[['h']], NULL)
            else
                K_setMinFunc(k[['h']], .dynsym(.lib$libhandle, paste("K_get", toupper(substring(value, 1, 1)),
                                                                     substring(value, 2), sep="")))
            k[['min.func']] <- value
        } else if (is.function(value)){ 
            k[['min.func']] <- function(h) {
                return(value(k))
            }

            k[['.min.func']] <- new.callback("p)d", k[['min.func']])
            
            K_setMinFunc(k[['h']], k[['.min.func']])
        } else 
            stop("The second argument should either be a string (chi2, chi2_nr, rms, jitter, loglik) or a function to minimize")
        
    } else if (idx == "progress") {
        stopifnot(is.function(value))
        k[['.progress']] <- function(cur, max, state, str) {
            ret <- value(k, list(cur=cur, max=max, job=.job, str=str))
            
            if (!is.numeric(ret))
                return(K_PROGRESS_CONTINUE)
            else
                return(ret)
        }
        
        k[['.progress.cb']] <- new.callback("iipZ)i", k[['.progress']])
        k[['progress']] <- value
        K_setProgress(k[['h']], k[['.progress.cb']])
    } else if (idx == 'filename') {
        K_setInfo(k[['h']], 'FileName', value)
    } else 
        k[[idx]] <- value
    return(k)
}

kset <- function(k, idx, val) {
    return(`$<-.kernel`(k, idx, val))
}


kallels <- function(k, keep.first = F) {
    ## Returns a matrix of orbital elements, including derived orbital elements [2]
    #
    # Args: 
    #	- k: kernel
    #
    # Returns:
    #	A matrix with all orbital elements (including semi-major axis, K, etc.)
    .check_kernel(k)
    m <- .gsl_matrix_to_R(K_getAllElements(k$h), free = TRUE)
    if (! keep.first) {
        m = m[-1,, drop=FALSE]
    }
    colnames(m) <- .allelements
    return(m)
}

kpars <- function(k) {
    ## Returns a vector of parameters [2]
    #
    # Args:
    #	- k: kernel
    #
    # Returns:
    #	A vector of parameters (e.g. offsets, trends)
    .check_kernel(k)
    v <- .gsl_vector_to_R(K_getPars(k$h))
    names(v) <- .params
    return(v)
}

kadd.data <- function(k, data, type = NA) {
    ## Adds a new dataset to the kernel. [5]
    #
    # Args:
    #	- k: kernel to add data to
    #	- data: either the path to a data file in textual format, or a matrix containing data
    #	- type: type of data contained in the data argument. One of the RV or TIMING constants
    #
    # Example:
    # fit <- knew()
    # kadd.data(fit, '~/data/mydata.vels')
    # print(fit$ndata)
    .check_kernel(k)
    if (! k$.epoch.set)
        k$epoch <- NaN
    
    if (k$nsets >= DATA_SETS_SIZE)
        stop(sprintf("You cannot add more than %d data sets.", DATA_SETS_SIZE))
    
    if (length(data) > 1 && ! class(data) == "matrix") {
        .mute(k)
        
        for (i in 1:length(data)) 
            kadd.data(k, data[i], type)
        
        .unmute(k)
        if (k$auto) kupdate(k)
        return(invisible())
    }
    
    if (class(data) == "character") {
        data <- normalizePath(data, mustWork=FALSE)
        if (! file.exists(data))
            stop(paste("Could not open", data))
        if (grepl(".vels", data)) {
            type <- T_RV
        } else if (grepl(".tds", data)) {
            type <- T_TIMING			
        } else if (grepl(".tt", data)) {
            type <- T_TIMING			
        }

        if (is.na(type))
            stop("Specify the type of data (one of RV or TIMING)")
        
        if (! is.nullptr(K_addDataFile(k$h, data, type))) {
        } else {
            stop("Could not open datafile")
        }
    } else if (class(data) == "matrix") {
        if (is.na(type))
            stop("Specify the type of data (one of RV or TIMING)")
        K_addDataTable(k$h, .R_to_gsl_matrix(data, gc=F), "", type)
        name <- deparse(substitute(data))
        
        if (type == RV)
            name <- paste(name, ".vels")
        else if (type == TIMING)
            name <- paste(name, ".tds")
        
    } else {
        stop("Data should be a string or a matrix")
    }
    
    if (k$auto) kupdate(k) else kcalculate(k)
    
    invisible()
}

kremove.data <- function(k, idx = -1) {
    # Removes the idx-th dataset [5]
    #
    # Args:
    #	k: kernel
    #	idx: the index of the dataset, or "all" if you want to remove all data.
    .check_kernel(k)
    
    if (idx == -1 || idx == "all") {
        K_removeData(k$h, -1)
    } else {
        stopifnot(idx <= k$nsets && idx >= 1)
        K_removeData(k$h, idx-1)
    }

    K_compileData(k$h)
    
    if (k$auto) kupdate(k)
    invisible()
}


.kdata.names <- c("TIME", "VAL", "ERR", "TDS_PLANET", "TDS_FLAG", "", "FLAG", "SVAL", "PRED", "SET", "")
kdata <- function(k, idx = 'all') {
    ## Returns a matrix containing the idx-th dataset (or all the data) [5].
    # The columns are:
    #
    # - TIME (time of observation)
    # - VAL (measurement)
    # - ERR (uncertainty on the measurement)
    # - SVAL (measurement, shifted vertically by the offset parameters)
    # - PRED (value of the fit at the time of measurement)
    # - SET (index of the dataset that contains the data point)
    #
    # Args:
    #	- k: kernel
    #	- idx: either the index of the dataset you are interested in, or 'all' to get all the data compiled from the loaded datasets.
    #
    # Examples:
    # k <- knew()
    # kadd.data(k, '~/Data/data.vels')
    # kadd.data(k, '~/Data/data2.vels')
    # kcalculate(k)
    # data2 <- kdata(k, 2)  # Returns the data originally read
    # data2[, ERR] <- 2 * data2[, ERR] # Inflates the uncertainty by a factor of 2
    # kdata(k, 2) <- data2 # Put data back into kernel
    .check_kernel(k)

    if (idx == 'all') {
        m <- .buf_to_R(K_compileData(k$h), k$ndata, DATA_SIZE)
        colnames(m) <- .kdata.names
        return(m)
    } else if (idx == 'res') {
        m <- .buf_to_R(K_compileData(k$h), k$ndata, DATA_SIZE)
        m[, VAL] <- m[, K_T_SVAL+1] - m[, K_T_PRED+1]
        m <- m[,c(TIME, VAL, ERR)]
        colnames(m) <- c("TIME", "VAL", "ERR")
        return(m)
    } else {
        stopifnot(class(idx) == "numeric")
        stopifnot(idx <= k$nsets && idx >= 1)
        m <- .gsl_matrix_to_R(K_getData(k$h, idx-1))
        colnames(m) <- .kdata.names
        return(m)
    }
}

`kdata<-` <- function(k, idx="all", value) {
    ## Sets a matrix containing the idx-th dataset (or all the data) [5]. See also @kdata.
    # Args:
    #	- k: kernel
    #	- idx: either the index of the dataset you are interested in, or 'all' to get all the data compiled from the loaded datasets.
    # - value: a matrix of observations.
    #
    # Examples:
    # k <- knew()
    # kadd.data(k, '~/Data/data.vels')
    # kadd.data(k, '~/Data/data2.vels')
    # kcalculate(k)
    # data2 <- kdata(k, 2)  # Returns the data originally read
    # data2[, ERR] <- 2 * data2[, ERR] # Inflates the uncertainty by a factor of 2
    # kdata(k, 2) <- data2 # Put data back into kernel
    #
    
    .check_kernel(k)
    stopifnot(class(value) == "matrix")
    if (idx != "all") {
        stopifnot(class(idx) == "numeric")
        stopifnot(idx > 0 && idx <= k$nsets)
        K_setData(k$h, idx-1, .R_to_gsl_matrix(value))
    } else {
        u <- unique(value[, SET])
        if (length(u) != k$nsets) 
            warning(sprintf("Dataset(s) index %s have no points left within the supplied matrix; they will be removed and the indices of the other datasets will be shifted accordingly.",
                            Reduce(paste, setdiff(1:k$nsets, u))-1))
        
        id <- 1
        for (i in (1:k$nsets)) {
            if (! (i-1) %in% u) {
                kremove.data(k, id)
            } else {
                m <- value[value[, SET] == (i - 1), ]
                m[, SET] <- id-1
                K_setData(k$h, id-1, .R_to_gsl_matrix(m))
                id <- id+1
            }
        }
    }
    if (k$auto) kupdate(k)
    return(k)
}


kcalculate <- function(k) {
    ## Recalculates or updates the statistics for the kernel (e.g. chi^2, rms, jitter, radial velocity response, etc.). [6]
    #
    # Args:
    #	- k: kernel
    #
    # Example:
    # k <- knew()
    # kadd.data(k, "~/Data/mydata.vels")
    # print(k$chi2)  # Prints the chi^2 with the data alone
    # k[] <- c(period=300, mass=1) # Adds a Jupiter-mass planet at P = 300 d
    # kcalculate(k) # Updates stats
    # print(k$chi2) # Prints new chi^2
    .check_kernel(k)
    K_calculate(k$h)
    return(K_getChi2(k$h))
}

kupdate <- function(k, calculate=TRUE) {
    if (calculate)
        kcalculate(k)
    
    if (length(k$hooks) > 0)
        for (i in 1:length(k$hooks))
            k$hooks[[i]](k)
}

kclone <- function(k) {
    ## Returns a new, independent copy of the kernel, with the same data and planets loaded. [1]
    # 
    # Args:
    #	k: a kernel
    #
    # Returns:
    #	A new kernel object that clones the input kernel object.
    .check_kernel(k)
    k2 <- new.env()
    
    for (n in ls(envir=k))
        k2[[n]] <- k[[n]]
    
    k2$h <- K_clone(k$h)
    k2$auto <- k$auto
    k2$.new <- TRUE
    k2$errors <- k$errors
    
    reg.finalizer(k2$h, .free.kernel)
    class(k2) <- "kernel"
    
    if (k2$auto) kupdate(k2)
    return(k2)
}

.extract <- function(line, header) {
    return(as.numeric(gsub(header, "", line)))
}

.trim <- function(line) {
    return(gsub("[ \t\n]+$", "", gsub("^[ \t\n]+", "", line)))
}

kload.old <- function(file, datafiles.dir=paste(dirname(file), "/datafiles")) {
    ## Loads an old-style fit (from the previous version of Systemic). [1]
    #
    # Args:
    #	file: path to the fit to be loaded
    #	datafiles.dir: path to the datafiles (.sys and .vels files)
    #
    # Returns:
    #	A new kernel object with data and parameters loaded according to the contents of file
    dir <- getwd()
    on.exit(setwd(dir))

    cat("Opening old-style systemic fit...\n")
    lines <- readLines(file)
    k <- knew()
    par <- 1
    rrvoffs <- FALSE
    for (line in lines) {
        line <- .trim(line)
        if (startsWith(line, "# InitialEpoch:")) {
            k$epoch <- .extract(line, "# InitialEpoch:")
        } else if (startsWith(line, "Parent")) {
            datafile <- .trim(gsub("\"", "", gsub("Parent", "", line)))
            setwd(datafiles.dir)
            if (file.exists(basename(datafile)))
                kload.datafile(k, basename(datafile))
            else {
                stop(sprintf("Could not load %s. You can specify where to locate datafiles using the datafiles.dir parameter.", datafile))
            }
            setwd(dir)
        } else if (startsWith(line, "OverallRVOffset")) {
            k['par', par] <- .extract(line, "OverallRVOffset")
            par <- par + 1
        } else if (startsWith(line, "RelativeRVOffsets")) {
            rrvoffs <- TRUE
        } else if (startsWith(line, "}"))
              rrvoffs <- FALSE
          else if (startsWith(line, "Period")) 
              kadd.planet(k, period=.extract(line, "Period"))
          else if (startsWith(line, "Mass")) 
              k[k$nplanets, "mass"] <- .extract(line, "Mass")
          else if (startsWith(line, "MeanAnomaly")) 
              k[k$nplanets, "ma"] <- .extract(line, "MeanAnomaly")
          else if (startsWith(line, "Eccentricity")) 
              k[k$nplanets, "ecc"] <- .extract(line, "Eccentricity")
          else if (startsWith(line, "LongOfPericenter")) 
              k[k$nplanets, "lop"] <- .extract(line, "LongOfPericenter")
          else if (startsWith(line, "Inclination")) 
              k[k$nplanets, "inc"] <- .extract(line, "Inclination")
          else if (startsWith(line, "Node")) 
              k[k$nplanets, "node"] <- .extract(line, "Node")
          else if (startsWith(line, "\"Trend\"")) 
              k['par', RV.TREND] <- .extract(line, "\"Trend\"")
          else if (rrvoffs) {
              k['par', par] <- .extract(line, "\".+\"")
              par <- par + 1
          }
    }
    
    if (k$nsets > 1) {
        for (i in 2:k$nsets)
            k['par', i] <- k['par', i] + k['par', 1]
    }
    for (i in 1:k$nplanets)
        k[i, 'inc'] <- k[i, 'inc']+90
    
    k$element.type <- JACOBI
    
    return(k)
}


kload.system <- function(k, datafile) {
    ## Loads the datafiles contained in a .sys file (see, e.g., the .sys files contained in the datafiles folder). [1]
    #
    # Args:
    #	- k: kernel to load data into
    #	- datafile: path to the .sys file

    dir <- getwd()
    
    if (! file.exists(datafile))
        stop(paste("Cannot open", datafile, ", current directory:", getwd()))

    K_addDataFromSystem(k$h, datafile)
    kupdate(k)
}

kload.datafile <- kload.system
kload <- function(file, skip = 0, chdir=TRUE) {
    ## Loads a kernel (previously saved with @ksave) from disk. [1]
    #
    # Args:
    #	- file: path to the file
    #	- skip: if multiple kernels are saved in a single file, the index of the kernel to read (usually only one kernel is saved per file, so specify 0)
    stopifnot(class(file) == "character")
    file <- normalizePath(file, mustWork=FALSE)
    header <- readLines(file, n=1)
    k <- NULL
    
    
    fid <- fopen(file, "r")
    
    if (is.nullptr(fid)) {
        stop(paste("Could not open file ", file))
    }
    k <- new.env()
    k$auto <- getOption("systemic.auto", F)
    k$min.method <- LM
    k$.new <- TRUE
    k$last.error.code <- integer(1)
    k$min.func <- 'chi2'
    
    if (substr(header, 1, 1) != "@") {
        k <- kload.old(file, chdir)
    } else {
        h <- K_load(fid, skip)
        if (is.nullptr(h)) {
            stop(paste("Could not parse file ", file))
            return(NULL)
        }
        
        k$h <- h
        reg.finalizer(k$h, .free.kernel)			
    }
    
    

    class(k) <- "kernel"	
    k$.epoch.set <- TRUE
    fclose(fid)
    kupdate(k)
    k$filename <- file
    return(k)
}

ksave <- function(k, file) {
    ## Saves a kernel (or a list of kernels) to a file. [1]
    #
    # Args:
    #	- k: a single kernel, or a list of kernels to save
    #	- file: file to save to. Multiple kernels might be saved to the same file.
    
    stopifnot(class(file) == "character")
    file <- normalizePath(file, mustWork=FALSE)
    
    if (class(k) == "kernel") {
        .check_kernel(k)
        
        fid <- fopen(file, "w")
        if (is.nullptr(fid)) {
            stop(paste("Could not open file ", file))
        }
        k$filename <- file
        K_save(k$h, fid)
        fclose(fid)
        
        return(invisible(k))
    } else if (class(k) == "list") {
        fid <- fopen(file, "w")
        if (is.nullptr(fid)) {
            stop(paste("Could not open file ", file))
        }

        for (ki in k) {
            .check_kernel(ki)
            ki$filename <- file
            K_save(ki$h, fid)
        }
        
        fclose(fid)
    }
}

kfind.peaks <- function(mat, column='power') {
    subm <- mat[which(diff(sign(diff(mat[, column]))) == -2)+1, ]
    if (mat[nrow(mat), column] > mat[nrow(mat)-1, column])
        subm <- rbind(subm, mat[nrow(mat),])
return(subm[order(subm[,column], decreasing=TRUE),])
}

print.periodogram <- function(x, what='peaks') {
    cat(sprintf("# Periodogram: pmin = %.2e, pmax = %.2e, samples = %d\n",
                attr(x, 'pmin'), attr(x, 'pmax'), attr(x, 'samples')))
    if (what == 'peaks') {
        cat(sprintf("# Peaks (sorted by power; window.cutoff = %e):\n", attr(x, 'window.cutoff')))
        print(attr(x, 'peaks')[, c('period', 'power', 'fap', 'window')])
        
        cat('# To print the full periodogram instead of power peaks, use\n',
            '# print(var, what="periodogram")\n',
            '# To save the periodogram, use write.f(var, file="p.txt")\n', sep="")
    } else if (what == 'periodogram') {
        print.matrix(x)
        cat('# To print the power peaks instead, use\n',
            '# print(var)\n',
            '# To save the periodogram, use write.f(var, file="myfile.txt")\n', sep="")
    }
    if (is.null(attr(x, 'is.boot')))
        cat("\n?(Note: The FAPs in this periodogram are roughly estimated analytically. Use kperiodogram.boot for a bootstrapped estimate.)\n")
    else
        cat(sprintf("\n# Trials: %d\n", attr(x, 'trials')))
}

kperiodogram <- function(k, per_type = "all", samples = getOption("systemic.psamples", 50000), pmin = getOption("systemic.pmin", 0.5), pmax = getOption("systemic.pmax", 1e4), data.flag = T_RV, timing.planet = NULL, val.col = SVAL, time.col = TIME, err.col = ERR, pred.col = PRED, plot = FALSE, print = FALSE, peaks = 25,
                         overplot.window=TRUE, .keep.h = FALSE) {
    ## Returns a periodogram of the supplied time series. [7]
    #
    # If the first parameter is a kernel, then this function will return 
    # periodogram of the loaded datasets (if per_type = "all") or the 
    # residuals (if per_type = "res"). If the first parameter is a matrix,
    # then this function will return a periodogram of the columns of the
    # matrix.
    #
    # The periodogram function will need three columns of data: a timestamp
    # column, a value column (e.g. the RV amplitude at that timestamp), and
    # an uncertainty column. The default indexes for those columns are TIME
    # [1], SVAL [8] and ERR [3]. 
    #
    # The false alarm probability returned is only an analytical
    # estimate. For more accurate false alarm probabilities, use 
    # kperiodogram.boot (which uses a bootstrap method to estimate the 
    # false alarm probabilities).
    #
    # Args:
    #	- k: either a kernel object, or a matrix
    #	- per_type: if k is a kernel object, one of "all" (periodogram of
    # 	the full compiled dataset) or "res" (periodogram of the 
    # 	residuals)
    #	- samples: number of periods (frequencies) at which to sample the
    #	periodogram
    # - peaks: identifies the N tallest peaks in the periodogram
    #	- pmin: minimum period at which to sample the periodogram
    #	- pmax: maximum period at which to sample the periodogram
    #	- data.flag: type of data (T_RV or T_TIMING)
    #	- timing.planet: if the data is of type T_TIMING, specifies which 
    #	transits to use to calculate the periodogram
    #	- val.col: the column to consider as the "value" column (by default, the SVAL column)
    #	- time.col: the column to consider as the "time" column (by default, the TIME column)
    #	- err.col: the column to consider as the "uncertainty" column (by default, the ERR column)
    #	- pred.col: the column to consider as the "model value" column 
    #	(used if per_type = "res") (by default, the PRED column)
    #	- plot: if TRUE, plot the periodogram after the calculation
    #	- overplot.window: if TRUE, overplot the periodogram of the sampling 
    #	window
    # - print: if TRUE, pretty-prints the periodogram sorted by power.
    # 
    # Returns:
    #	 A matrix with columns containing, respectively: period, power 
    # 	at that period, (analytical) false alarm probability, 
    #	unnormalized power, tau, power of the sampling window at that 
    #	period
    d <- NULL

    
    if (! per_type %in% c('all', 'res'))
        stop("per_type should be one of ['all', 'res']") 
    
    if (class(k) == "kernel") {
        .check_kernel(k)

        if (k$ndata == 0)
            stop("No data contained in the kernel.")

        d <- kdata(k)
        n <- kpars(k)[DATA.NOISE1:DATA.NOISE10]
        d[, err.col] <- sqrt(d[, err.col]^2 + n[d[, SET]+1]^2)
    } else {
        d <- k
        if (val.col > ncol(d))
            val.col <- 2
        stopifnot(nrow(d) > 1 && ncol(d) >= 3)
    }
    
    if (per_type == "res") {
        d[, val.col] <- d[, val.col] - d[, pred.col]
    } 
    
    if (! is.null(data.flag) && FLAG <= ncol(d)) {
        d <- d[d[, FLAG] == data.flag, ]
        
        if (data.flag == T_TIMING && per_type == "all") {
            if (is.null(timing.planet))
                stop("Please specify the option timing.planet = n to select the timing data associated with the n-th planet")
            val.col <- 2
            d <- d[d[, TDS_PLANET] == timing.planet && d[, TDS_FLAG] != TDS_SECONDARY, ]
            
            if (k$nplanets < timing.planet) {
                p <- diff(d[, time.col])
                minp <- median(p[p < min(p) * (1+0.4)])
            }
            else
                minp <- k[timing.planet, 'period']
            
            idx <- floor(d[, val.col] / minp)
            fit <- lm(d[, val.col] ~ idx, weights=1./d[, err.col])
            
            d[, val.col] <- fit$residuals
        }
    }


    stopifnot(nrow(d) > 0)

    m <- .R_to_gsl_matrix(d)
    per <- ok_periodogram_ls(m, samples, pmin, pmax, 0, time.col-1, val.col-1, err.col-1, NULL)

    if (samples > 1e4) {
        .periodogram.tol <- double(1)
        .periodogram.tol[1] <- 1e-4

        resampled <- .gsl_matrix_to_R(ok_resample_curve(per, 0, 1, 1, 10000,
                                                        2000, .periodogram.tol, 5, TRUE), free=TRUE)
        colnames(resampled) <- .periodogram
    } else 
        resampled <- NULL
    
    
    mat <- .gsl_matrix_to_R(per, free = TRUE, .keep.h = .keep.h)
    
    colnames(mat) <- .periodogram
    class(mat) <- "periodogram"
    
    if (plot)
        plot(mat, overplot.window=overplot.window)
    peaks.m <- kfind.peaks(mat)
    if (nrow(peaks.m) > 0) {
        mfap <- mat[mat[,'fap'] < 1, , drop=FALSE]
        if (nrow(mfap) > 5) {
            window.cutoff <- 10^approxfun(log10(mfap[,'fap']), log10(mfap[,'power']))(log10(1e-6))
            if (!is.na(window.cutoff)) {
                peaks.m <- peaks.m[peaks.m[,'window'] < window.cutoff, , drop=FALSE]
                if (print)
                    cat(sprintf("# Peaks with power in window > %e ignored\n", window.cutoff))
                attr(mat, 'window.cutoff') <- window.cutoff
            }
        }
        peaks <- min(peaks, nrow(peaks.m))
        attr(mat, "peaks") <- peaks.m[1:peaks, , drop=FALSE]
    }

    

    attr(mat, 'pmin') <- pmin
    attr(mat, 'pmax') <- pmax
    attr(mat, 'samples') <- samples
    attr(mat, 'resampled') <- resampled
    
    if (print) {
        print(mat, 'peaks')
        return(invisible(mat))
    } else {
        return(mat)
    }
}

kperiodogram.boot <- function(k, per_type = "all", trials = 1e5, samples = getOption("systemic.psamples", 50000), pmin = getOption("systemic.pmin", 0.5), pmax = getOption("systemic.pmax", 1e4), data.flag = T_RV, timing.planet = NULL, val.col = SVAL, time.col = TIME, err.col = ERR, seed = sample(1:1e4, 1), plot = FALSE, print = FALSE,
                              overplot.window=TRUE, peaks=25) {
    ## Returns a periodogram of the supplied time series, where the false alarm probabilities are estimated using a bootstrap method. [7]
    #
    # If the first parameter is a kernel, then this function will return 
    # periodogram of the loaded datasets (if per_type = "all") or the 
    # residuals (if per_type = "res"). If the first parameter is a matrix,
    # then this function will return a periodogram of the columns of the
    # matrix.
    #
    # The periodogram function will need three columns of data: a timestamp
    # column, a value column (e.g. the RV amplitude at that timestamp), and
    # an uncertainty column. The default indexes for those columns are TIME
    # [1], SVAL [8] and ERR [3]. 
    #
    # The false alarm probabilities are estimated by computing "trials" 
    # periodograms of gaussian noise. The routine is automatically 
    # parallelized.
    #
    # Args:
    # The arguments are the same as the @kperiodogram function, plus the following:
    #
    #	- trials: number of periodograms to use in the boostrap estimation
    #	- samples: number of periods (frequencies) at which to sample the
    #	periodogram
    #
    # Returns:
    #	A matrix with columns containing, respectively: period, power 
    # at that period, false alarm probability (calculated in the bootstrap procedure), 
    #	unnormalized power, tau, power of the sampling window at that 
    #	period
    
    .job <<- "Bootstrap periodogram"
    
    d <- NULL
    rng <- NULL
    
    if (class(k) == "kernel") {
        .check_kernel(k)
        stopifnot(k$ndata > 0)			
        d <- kdata(k)
        n <- kpars(k)[DATA.NOISE1:DATA.NOISE10]
        d[, err.col] <- sqrt(d[, err.col]^2 + n[d[, SET]+1]^2)
    } else {
        stopifnot(nrow(d) > 1 && ncol(d) >= 3)		
        d <- k
    }
    
    if (per_type == "res") {
        val.col <- PRED
        
        d[, PRED] <- d[, SVAL] - d[, PRED]
    }
    
    if (! is.null(data.flag)) {
        d <- d[d[, FLAG] == data.flag, ]
        
        if (data.flag == T_TIMING && per_type == "all") {
            if (is.null(timing.planet))
                stop("Please specify the option timing.planet = n to select the timing data associated with the n-th planet")
            val.col <- 2
            d <- d[d[, TDS_PLANET] == timing.planet, ]
						
            if (k$nplanets < timing.planet) {
                p <- diff(d[, time.col])
                minp <- median(p[p < min(p) * (1+0.4)])
            }
            else
                minp <- k[timing.planet, 'period']
            
            idx <- floor(d[, val.col] / minp)
            fit <- lm(d[, val.col] ~ idx, weights=1./d[, err.col])
            
            d[, val.col] <- fit$residuals
        }
    }

    
    m <- .R_to_gsl_matrix(d)
    per <- ok_periodogram_boot(m, trials, samples, pmin, pmax, 0, time.col-1, val.col-1, err.col-1, seed, NULL, K_getProgress(k$h))
    .job <<- "Bootstrap periodogram"


    if (samples > 1e4) {
        .periodogram.tol <- double(1)
        .periodogram.tol[1] <- 1e-5

        resampled <- .gsl_matrix_to_R(ok_resample_curve(per, 0, 1, 1, 10000,
                                                        2000, .periodogram.tol, 5, TRUE), free=TRUE)
        colnames(resampled) <- .periodogram
    } else
        resampled <- NULL
    
    m <- .gsl_matrix_to_R(per, free = T)
    colnames(m) <- .periodogram

    class(m) <- "periodogram"
    
    if (plot)
        plot(m, overplot.window=overplot.window)

    peaks.m <- kfind.peaks(m)
    if (nrow(peaks.m) > 0) {
        peaks <- min(peaks, nrow(peaks.m))
        attr(m, "peaks") <- peaks.m[1:peaks, ]
    }
    attr(m, "is.boot") <- TRUE

    attr(m, 'pmin') <- pmin
    attr(m, 'pmax') <- pmax
    attr(m, 'samples') <- samples
    attr(m, 'trials') <- trials
    attr(m, 'resampled') <- resampled
    
    if (print) {
        print(m, "peaks")
        return(invisible(m))
    } else {
        return(m)
    }

    return(m)
}

kflag <- function(k, row, column) {
    ## Returns/sets the flag (active, minimized or inactive) for an orbital element or parameter. [2]
    # This corresponds to clicking one of the "semaphore" buttons next to the fit parameters in the user interface.
    #
    # The value set can be one of INACTIVE (not minimized over, not counted 
    # as a parameter in reduced chi^2), ACTIVE (not minimized over, counted 
    # as a parameter in reduced chi^2), MINIMIZE (minimized over). By 
    # default, many orbital elements are ACTIVE + MINIMIZE (these appear 
    # as green in the GUI).
    #
    # Args:
    #	- k: kernel object
    #	- row: either the planet index, or "par" to specify a parameter 
    #	(e.g. an RV offset)
    #	- col: either the orbital element index (one of `constants`)
    #	- value: one of INACTIVE, ACTIVE, MINIMIZE, or ACTIVE + MINIMIZE
    #
    # Examples:
    # print(kflag(k, 1, 'period')) # prints the flag of the period parameter for the first planet
    # kflag(k, 1, 'period')) <- ACTIVE + MINIMIZE # marks the parameter as ACTIVE and MINIMIZE
    # kminimize(k) 
    .check_kernel(k)
		
    column <- .label_to_index(column)

    if (row == "par") {
        return(K_getParFlag(k$h, column - 1))
    } else {
        stopifnot(row <= k$nplanets)
        stopifnot(column <= ELEMENTS_SIZE && column >= 1)
        return(K_getElementFlag(k$h, row, column - 1))
    }
}

`kflag<-` <- function(k, row, column, value) {
    .check_kernel(k)
    
    column <- .label_to_index(column)
    
    if (length(row) == 1 && row == "par") {
        K_setParFlag(k$h, column-1, value)
    } else {
        stopifnot(row <= k$nplanets)		
        .stop.ifnot(column <= ELEMENTS_SIZE && column >= 1, sprintf("Index %d out of range [%d, %d]", column, 1, ELEMENTS_SIZE))						
        K_setElementFlag(k$h, row, column-1, value)
    }
    
    if (k$auto) {
        kupdate(k, calculate=FALSE)	
    }
    return(k)
}

kflags <- function(k, what = "all", type='standard') {
    .check_kernel(k)
    
    m <- NULL
    if (what == "els") {
        stopifnot(k$nplanets > 0)
        g <- expand.grid(1:k$nplanets, 1:ELEMENTS_SIZE)
        m <- matrix(mapply(function(i, j) { return(K_getElementFlag(k$h, i, j-1)) }, g[[1]], g[[2]]), nrow=k$nplanets)
        colnames(m) <- .elements
    } else if (what == "par") {
        m <- sapply(1:PARAMS_SIZE, function(idx) return(K_getParFlag(k$h, idx-1)))
        names(m) <- .params
    } else if (what == "all") {
        m <- list()
        m$els <- kflags(k, what="els", type=type)
        m$par <- kflags(k, what="par", type=type)		
        return(m)
    }
    
    if (type == "human") {
        m[m == 6] <- "M+A"
        m[m == "2"] <- "A"
        m[m == "4"] <- "M"
        m[m == "0"] <- ""
    }
    
    return(m)
}

`kflags<-` <- function(k, value) {
    flags <- value
    if (class(flags) == "list") {
        .check_kernel(k)
        stopifnot(k$nplanets == nrow(flags$els))
        g <- expand.grid(1:k$nplanets, 1:ELEMENTS_SIZE)
        mapply(function(i, j) { return(K_setElementFlag(k$h, i, j-1, flags$els[i, j])) }, g[[1]], g[[2]])
        sapply(1:PARAMS_SIZE, function(idx) K_setParFlag(k$h, idx-1, flags$par[idx]))
        if (k$auto)
            kupdate(k, calculate=FALSE)
    } 
    return(k)
}

kminimize <- function(k, iters = 5000, algo = NA, de.CR = 0.2,
                      de.NPfac = 10, de.Fmin = 0.5, de.Fmax = 1.0, de.use.steps = FALSE,
                      sa.T0 = k$chi2, sa.alpha=2, sa.auto=TRUE, sa.chains=4, repeat.steps = 10) {
    ## Minimizes the chi^2 of the fit. [3]
    #
    # kminimize uses one of the built-in algorithms to minimize the
    # reduced chi^2 of the fit. The algorithm is one of the following:
    #
    # - SIMPLEX uses the Nelder-Meade algorithm (as implemented in GSL)
    # to search for a local minimum.
    # - LM uses the Levenberg-Marquardt algorithm (as implemented in GSL)
    # to search for a local minimum.
    # - SA uses a simple implementation of the simulated annealing
    # algorithm.
    # - DE uses a simple implementation of the differential evolution
    # algorithm.
    #
    # The minimization algorithms may use the parameter steps set by
    # @kstep as initial scale parameters to explore the chi^2 landscape.
    #
    # The target function to be minimized is defined by k$min.func, which
    # has a value of "chi2" by default.
    #
    # Args:
    # - k: kernel to minimize
    # - iters: maximum number of iterations
    # - algo: one of SIMPLEX, LM, SA or DE. If none is specified, uses
    # the value in k$min.method
    # - sa.T0: for SA, the initial temperature of the annealer
    # - sa.alpha: the index of the annealer (T = T0 (1 - (n/N)^alpha))
    # - sa.auto: automatically derive steps that produce a variation of chi^2 = 10% T0
    # - de.CR: crossover probability for DE
    # - de.Fmin, de.Fmax: differential weight for DE
    # - repeat.steps: repeats the minimization algorithm if there is a change in chi^2 for max number of steps
    
    .check_kernel(k)
    
    if (is.na(algo))
        algo <- k$min.method
    
    p <- kflags(k, 'par')[DATA.NOISE1:DATA.NOISE10]
    if (any(bitAnd(p, MINIMIZE)==MINIMIZE) && k[['min.func']] == 'chi2')
        warning("Minimizing with one of the DATA.NOISE parameters selected")
    
    opts = c(K_OPT_SA_T0, sa.T0, K_OPT_SA_ALPHA, sa.alpha,
        K_OPT_SA_AUTO_STEPS, sa.auto, K_OPT_SA_CHAINS, sa.chains,
        K_OPT_DE_CR, de.CR, K_OPT_DE_NP_FAC, de.NPfac,
        K_OPT_DE_F_MIN, de.Fmin, K_OPT_DE_F_MAX, de.Fmax,
        K_OPT_DE_USE_STEPS, if (de.use.steps) 1 else 0, K_DONE)
    
    .job <<- "Minimization"
    stopifnot(k$ndata > 0)
    
    on.exit(if (k$auto) kupdate(k))
    old <- K_getMinValue(k$h)
    for (i in 1:repeat.steps) {
        K_minimize(k$h, algo, iters, as.numeric(opts))
        if (abs(K_getMinValue(k$h) - old) < 1e-6)
            break
        old <- K_getMinValue(k$h)
    }
    return(k$chi2)
}


kcrossval.l1o <- function(k, iters = 5000, algo = NA, type=NA) {
    ## Runs the "leave-one-out" cross validation algorithm. [4]
    #
    # Args:
    # - k: the kernel to run the routine on.
    # - iters: number of iterations (see @kminimize)
    # - algo: minimization algorithm (see @kminimize)
    # - type: if NA, returns the result of the cross-validation algorithm; otherwise, runs the cross-validation algorithm by removing planets of progressively smaller k.
    .check_kernel(k)
    .job <<- "Cross validation"
    if (is.na(algo))
        algo <- k$min.method
    stopifnot(k$ndata > 0)
    
    if (is.na(type))
        return(K_crossval_l1o(k$h, algo, iters, NULL))
    else {
        cat(paste("?This routine will attempt to automatically remove and refit planets one by one (starting from the smallest K),",
                  "using the results of the leave-1-out cross-validation routine to compare between models.",
                  "If the output of a fit with (n+1) planets is larger than that of a fit with n planets,",
                  "it is likely the (n+1)-planets model is worse than the n-planets model \n"))
        .k <- kclone(k)
        kels(.k) <- kels(.k)[sort.list(kallels(.k)[, 'k'], decreasing=TRUE),,drop=FALSE]
        .k$auto <- F

        cat("# ", .k$nplanets, " planets\n")
        ret <- c(kcrossval.l1o(.k))
        print(ret[1])
        
        for (i in .k$nplanets:1) {
            kremove.planet(.k, i)
            cat("# ", .k$nplanets, " planets\n")
            invisible(kminimize(.k))
            ret[length(ret)+1] <- kcrossval.l1o(.k)
            print(ret[length(ret)])
        }
        return(ret)
    }
}

kminimize1d <- function(k, row, column, algo = NA, iters = 5000) {
    ## Minimizes the specified parameter [3]
    #
    # Args:
    # - k: kernel
    # - row: either the index of the planet, or "par" to minimize a data parameter
    # - column: the planet parameter (e.g. 'period', 'ma', 'ecc') or the data parameter
    # - algo, iters: see @kminimize.
    .check_kernel(k)
    stopifnot(k$ndata > 0)
    
    if (is.na(algo))
        algo <- k$min.method
    
    column <- .label_to_index(column)
    .job <<- "1-D Minimization"
    if (row == "par") row = 0
    K_1dminimize(k$h, algo, iters, row, column - 1, NULL)
    
    on.exit(if (k$auto) kupdate(k))
    
    return(k$chi2)
}

kminimized.indices <- function(k) {
    idx1 <- integer(1)
    idx2 <- integer(1)
    v <- sapply(1:k$nrpars, function(i) {
        K_getMinimizedIndex(k$h, i-1, idx1, idx2)
        return(c(idx1, idx2+1))
    })
    return(v)
}

krange <- function(k, row, column) {
    ## Returns the allowed range of the given parameter. [2]
    #
    # See @krange<-.
    .check_kernel(k)

    a <- NaN
    b <- NaN	

    column <- .label_to_index(column)
    
    if (row == "par" || row == -1) {
        K_getParRange(k$h, column - 1, a, b)
        return(c(min=a, max=b))
    } else {
        stopifnot(row > 0 && row <= k$nplanets)
        K_getElementRange(k$h, row, column - 1, a, b)
        return(c(min=a, max=b))
    }
}

`krange<-` <- function(k, row, column, value) {
    ## Sets the allowed range of the given parameter. [2]
    #
    # A parameter range forces a parameter to lie within a specified
    # interval. The interval is specified as a two-element vector,
    # c(min, max). If either (or both) min or max is NaN, then
    # the minimum/maximum is not enforced.
    #
    # Args:
    # - k: kernel
    # - row: Either the index of the planet, or 'par' to specify the range of a data parameter.
    # - column: Either the name of the planet parameter (e.g. 'period', 'ma', etc.) or the name of the data parameter.
    # - value: A two-element vector c(min, max). min and max can be NaN.
    #
    # Example:
    # krange(k, 1, 'period') <- c(1, 100) # Constrain between 1 and 100 days
    # krange(k, 1, 'mass') <- c(0.1, NaN) # The mass must be larger than 0.1 Jupiter masses
    .check_kernel(k)
    
    column <- .label_to_index(column)
		
    value[1] <- if (is.na(value[1])) NaN else value[1]
    value[2] <- if (is.na(value[2])) NaN else value[2]
    
    if (row == "par" || row == -1) {
        K_setParRange(k$h, column - 1, value[1], value[2])		
    } else {
        stopifnot(row > 0 && row <= k$nplanets)
        if (row == "all")
            for (i in 1:k$nplanets) 
                K_setElementRange(k$h, row, i-1, value[1], value[2])
        else
            K_setElementRange(k$h, row, column - 1, value[1], value[2])
    }
    
    if (k$auto) kupdate(k)
    
    return(k)
}

kranges <- function(k) {
    .check_kernel(k)

    min <- kels(k)
    max <- kels(k)
    for (i in 1:k$nplanets)
        for (j in 1:ELEMENTS_SIZE) {
            r <- krange(k, i, j)
            min[i, j] <- r[1]
            max[i, j] <- r[2]
        }

    return(list(min=min, max=max))
}

ksteps <- function(k, row, column) {
    .check_kernel(k)
    
    ks <- kels(k)
    
    for (i in 1:k$nplanets)
        for (j in 1:ELEMENTS_SIZE)
            ks[i, j] = K_getElementStep(k$h, i, j-1)
    return(ks)
}

kstep <- function(k, row, column) {
    ## Returns the "step" for a given parameter.
    #
    # See @`kstep<-`.
    .check_kernel(k)
    
    column <- .label_to_index(column)
    
    if (row == "par" || row == -1) {
        return(K_getParStep(k$h, column - 1))
    } else {
        stopifnot(row > 0 && row <= k$nplanets)
        return(K_getElementStep(k$h, row, column - 1))
    }
}

`kstep<-` <- function(k, row, column, value) {
    ## Sets the "step" for a given parameter.
    #
    # The "step" for a parameter specifies the typical step taken by
    # minimization routines for the parameter. For the SIMPLEX algorithm,
    # it represents the initial size of the simplex; for the LM algorithm,
    # it specifies the step over which gradients are calculated.
    #
    # Args:
    # - k: kernel
    # - row: either the planet index, or 'par' to specify the step of a data parameter.
    # - column: Either the name of the planet parameter (e.g. 'period', 'ma', etc.) or the name of the data parameter.
    # - value: The value of the step.
    .check_kernel(k)
    
    column <- .label_to_index(column)
    
    if (row == "par") {
        K_setParStep(k$h, column - 1, value)
        return(k)
    } else {

        if (row == "all") {
            for (i in 1:k$nplanets) {
                K_setElementStep(k$h, i - 1, column - 1, value)
            }
        } else {
            stopifnot(row > 0 && row <= k$nplanets)			
            K_setElementStep(k$h, row, column - 1, value)
        }
    }

    if (k$auto)
        kupdate(k, calculate=FALSE)				
    return(k)

}

kselect <- function(k, row = "all", column = "all") {
    ## Selects (make available for minimization) the given parameters.
    #
    # This corresponds to the action of clicking on the semaphore buttons
    # in the user interface to select a parameter for minimization.
    #
    # Args:
    # - k: the kernel
    # - row: either the planet index, or 'par' to specify a data parameter.
    # - column: either the name of the planet parameter (e.g. 'period', 'ma', etc.) or the name of the data parameter.
    #
    # Example:
    # kselect(k, 1, 'period') # Selects the period parameter of the 1st planet
    # kselect(k, 'par', 1) # Selects the first data parameter (the offset of the first parameter)
    .check_kernel(k)
    column <- .label_to_index(column)
    
    if (length(row) == 1 && row == "par") {
        if (length(column) == 1 && column == "all")
            for (i in 1:k$nsets)
                kselect(k, "par", i)
        else
            kflag(k, "par", column) <- bitOr(kflag(k, "par", column), bitOr(MINIMIZE, ACTIVE))
        
    } else {
        
        if (length(row) == 1 &&  row == "all") {
            if (k$nplanets > 0) {
                for (i in 1:k$nplanets) {
                    kselect(k, i, column)
                }
            }
            
            
        } else {
            stopifnot(row > 0 && row <= k$nplanets)
            if (length(column) == 1 && column == "all") {
                activateable <- which(bitAnd(kflags(k)$els[row,], ACTIVE) == ACTIVE)
                for (j in activateable) {
                    kselect(k, row, j)
                }
            } else {
                kflag(k, row, column) <- bitOr(kflag(k, row, column), bitOr(MINIMIZE, ACTIVE))
            }
        }
    }
    
    if (k$auto)
        kupdate(k, calculate=FALSE)
}

kdeselect <- function(k, row = "all", column = "all") {
    ## Deselects (exclude from minimization) the given parameters.
    #
    # See @kselect.
    .check_kernel(k)
    column <- .label_to_index(column)	

    
    if (length(row) == 1 && row == "par") {
        if (length(column) == 1 && column == "all")
            for (i in 1:PARAMS_SIZE)
                kdeselect(k, "par", i)
        else
            kflag(k, "par", column) <- bitAnd(kflag(k, "par", column), bitFlip(MINIMIZE))
    } else {
        if (length(row) == 1 && row == "all") {
            if (k$nplanets > 0)
                for (i in 1:k$nplanets) {
                    kdeselect(k, i, column)
                }
        } else {
            stopifnot(row > 0 && row <= k$nplanets)
            if (length(column) == 1 && column == "all") {
                for (j in PER:LOP) {
                    kdeselect(k, row, j)
                }
            } else {
                kflag(k, row, column) <- bitAnd(kflag(k, row, column), bitFlip(MINIMIZE))
            }
        }
    }
    
    if (k$auto)
        kupdate(k, calculate=FALSE)	
}

kactivate <- function(k, row = "all", column = "all") {
    .check_kernel(k)
    column <- .label_to_index(column)

		
    if (length(row) == 1 && row == "par") {
        if (length(column) == 1 && column == "all")
            for (i in 1:k$nsets)
                kactivate(k, "par", i)
        else
            kflag(k, "par", column) <- bitOr(kflag(k, "par", column), ACTIVE)
        
    } else {
        if (length(row) == 1 &&  row == "all") {
            if (k$nplanets > 0)
                for (i in 1:k$nplanets) {
                    kactivate(k, i, column)
                }
            
        } else {
            stopifnot(row > 0 && row <= k$nplanets)
            if (length(column) == 1 && column == "all") {
                for (j in PER:LOP) {
                    kactivate(k, row, j)
                }
            } else {
                kflag(k, row, column) <- bitOr(kflag(k, row, column), ACTIVE)
            }
        }
    }
    
    if (k$auto)
        kupdate(k, calculate=FALSE)
}

kdeactivate <- function(k, row = "all", column = "all") {
    .check_kernel(k)
    column <- .label_to_index(column)	

    
    if (length(row) == 1 && row == "par") {
        if (length(column) == 1 && column == "all")
            for (i in 1:PARAMS_SIZE)
                kdeactivate(k, "par", i)
        else
            kflag(k, "par", column) <- bitAnd(kflag(k, "par", column), bitFlip(bitOr(ACTIVE, MINIMIZE)))
    } else {
        if (length(row) == 1 && row == "all") {
            for (i in 1:k$nplanets) {
                kdeactivate(k, i, column)
            }
        } else {
            stopifnot(row > 0 && row <= k$nplanets)
            if (length(column) == 1 && column == "all") {
                for (j in PER:LOP) {
                    kdeactivate(k, row, j)
                }
            } else {
                kflag(k, row, column) <- bitAnd(kflag(k, row, column), bitFlip(bitOr(ACTIVE, MINIMIZE)))
            }
        }
    }
    
    if (k$auto)
        kupdate(k, calculate=TRUE)	
}

kbootstrap <- function(k, algo = NA, trials = 5000, warmup = 0, min_iter = 2000, plot = FALSE, print = FALSE, save=NA) {
    ## Runs the bootstrap routine on the given kernel. [4]
    #
    # This function runs the bootstrap algorithm to estimate the
    # uncertainty on the parameters.
    #
    # Args:
    # - k: the kernel to run bootstrap on
    # - algo: the algorithm to use in minimization passes (see @kminimize)
    # - trials: the number of resampling trials
    # - plot: plots the resulting uncertainty object
    # - print: prints the resulting uncertainty object
    .check_kernel(k)
    stopifnot(k$ndata > 0)
    
    .job <<- "Bootstrap"

    if (is.na(algo))
        algo <- k$min.method
    
    kl <- K_bootstrap(k$h, trials, warmup, algo, min_iter, NULL)
    if (is.nullptr(kl)) {
        return(NULL)
    } else {
        a <- .klnew(kl, k, type="bootstrap", desc=sprintf("algo = %d, trials = %d", algo, trials))
        if (plot)
            plot(a)
        if (print)
            print(a)
        if (!is.na(save))
            save(a, file=save)
        return(a)
    }
}

K_LIMITS <- list()
K_LIMITS[[PER]] <- c(1e-2, 4e4)
K_LIMITS[[MASS]] <- c(1e-2, 100)
K_LIMITS[[MA]] <- c(0, 360)
K_LIMITS[[ECC]] <- c(0, 0.99)
K_LIMITS[[LOP]] <- c(0, 360)
K_LIMITS[[INC]] <- c(0, 360)
K_LIMITS[[NODE]] <- c(0, 360)
K_LIMITS[['par']] <- c(-100, 100)

kmcmc <- function(k, chains= 2, temps = 1, start = "perturb", noise=TRUE, skip.first = 1000, discard = k$nrpars * 10, R.stop = 1.1, 
                  min.length = 5000, max.iters = -1, auto.steps = TRUE, acc.ratio = 0.44, plot = FALSE, print = FALSE, save=NA,
                  debug.verbose.level = 1, random.log=TRUE) {
    ## Runs the MCMC routine on the given kernel. [4]
    #
    # This function runs a simple implementation of MCMC on the kernel
    # to estimate uncertainty on the parameters. 
    #
    # Args:
    # - k: the kernel to run bootstrap on, or a list of kernels with different parameters which represent the starting initial conditions.
    # - chains: number of chains to run in parallel
    # - skip.first: discard the first iterations
    # - R.stop: the Gelman-Rubin statistic used to estimate when to stop the routine.
    # - discard: only retain every n-th element of the chain
    # - min.length: minimum number of iterations
    # - acc.ratio: the acceptance ratio
    # - print: prints the resulting uncertainty object
    # - plot: plots the resulting uncertainty object
    .job <<- "MCMC"

    stopifnot(discard > 1)
    stopifnot(temps == 1)
    
    ka <- list()
    if (class(k) == "kernel") {
        .check_kernel(k)
        stopifnot(k$ndata > 0)	
        
        if (start == "random") {
            for (ch in 1:chains) {
                kch <- kclone(k)
                for (i in 1:k$nplanets)
                    for (j in 1:ELEMENTS_SIZE) 
                        if (bitAnd(kflag(k, i, j), MINIMIZE) == MINIMIZE) {
                            r <- krange(k, i, j)
                            min <- if (is.nan(r[1])) K_LIMITS[[j]][1] else r[1]
                            max <- if (is.nan(r[2])) K_LIMITS[[j]][2] else r[2]						
                            if (j == PER || j == MASS && random.log)
                                kch[i, j] <- exp(log(min) + runif(1) * (log(max)-log(min)))
                            else
                                kch[i, j] <- min + runif(1) * (max-min)
                        }
                for (i in 1:PARAMS_SIZE)
                    if (bitAnd(kflag(k, 'par', i), MINIMIZE) == MINIMIZE) {
                        r <- krange(k, 'par', i)
                        min <- if (is.nan(r[1])) K_LIMITS[['par']][1] else r[1]
                        max <- if (is.nan(r[2])) K_LIMITS[['par']][2] else r[2]						
                        kch['par', i] <- min + runif(1) * (max-min)
                    }
                ka[[ch]] <- kch
            }
            
        } else if (start == "perturb") {
            for (ch in 1:chains) {
                ka[[ch]] <- kperturb(k)
            }
        }
    } else if (class(k) == "list") {
        for (ch in 1:length(k)) {
            .check_kernel(k[[i]])
            stopifnot(k$ndata > 0)
            ka[[ch]] <- kclone(k[[ch]])
        }
        chains <- length(k)
        start <- "user-defined"
    }
    
    kbuf <- ok_bridge_kernel_buf(NULL, chains, NULL)
    
    for (i in 1:chains) {
        if (noise)
            for (j in 1:k$nsets) kselect(ka[[i]], 'par', j + DATA_SETS_SIZE)
        K_setProgress(ka[[i]]$h, K_getProgress(k$h))
        ok_bridge_kernel_buf(kbuf, i-1, ka[[i]]$h)
    }
    
    opts <- c(K_OPT_MCMC_NSTOP, max.iters,
              K_OPT_MCMC_NMIN, min.length,
              K_OPT_MCMC_VERBOSE_DIAGS, debug.verbose.level,
              K_OPT_MCMC_ACCRATIO, acc.ratio,
              K_OPT_MCMC_SKIP_STEPS, if (auto.steps) 0 else 1,
              DONE)
		
    kl <- K_mcmc_mult(kbuf, chains, temps, skip.first, discard, opts, R.stop, NULL)
    print(is.nullptr(kl))
    if (is.nullptr(kl)) {
        ok_bridge_kernel_buf(kbuf, -chains, NULL);
        return(NULL)
    } else {
        ok_bridge_kernel_buf(kbuf, -chains, NULL);		
        a <- .klnew(kl, k, type="mcmc", desc=sprintf("chains = %d, R.stop = %e, start = %s, noise=%s, skip = %d, discard = %d, tot. length = %d", chains, R.stop, start, 
                                            noise, skip.first, discard, KL_getSize(kl)), flags=kflags(ka[[1]], 'par'))
        
        if (plot)
            plot(a)
        if (print)
            print(a)
        if (!is.na(save))
            save(a, file=save)
        return(a)
    }
}



kxyz <- function(k, internal=TRUE) {
    ## Returns the cartesian coordinates of the bodies in the system. [9]
    #
    # Coordinates are returned in internal units (K2, AU, and day) if internal = TRUE,
    # else returned in cgs.
    .check_kernel(k)
    int <- .gsl_matrix_to_R(K_getXYZ(k$h))
    
    colnames(int) <- c('mass', 'x', 'y', 'z', 'u', 'v', 'w')
    if (internal)
        return(int)
    else {
        int[,1] <- int[,1] / K2 * MSUN
        int[,2:4] <- int[,2:4] * AU
        int[,5:7] <- int[,5:7] * AU / DAY
        return(int)
    }
        
}

keltype <- function(k) {
    .check_kernel(k)
    return (K_getElementType(k$h))
}

`keltype` <- function(k, value) {
    ## Sets the orbital elements format. [2]
    #
    # Args:
    # - k: kernel
    # - value: one of ASTROCENTRIC for astrocentric elements, or JACOBI for Jacobi elements.
    .check_kernel(k)
    
    K_setElementType(k$h, value)
    if (k$auto) kupdate(k)
    return(k)
}

kperturb <- function(k) {
    .check_kernel(k)
    k2 <- kclone(k)
    K_perturb(k2$h)
    if (k2$auto) kupdate(k2)
    return(k2)
}


ktrange <- function(k) {
    .check_kernel(k)
    d <- kdata(k)
    return(c(min(d[,TIME]), max(d[,TIME])))
}

krvcurve <- function(k, times=seq(from=min(ktrange(k)), to=max(ktrange(k)), length.out=getOption("systemic.rvsamples", 5000))) {
    ## Calculates the radial velocity curve over the specified time vector.
    #
    # Args:
    # - k: the kernel
    # - times: a vector of times where to sample the radial velocity curve.
    .check_kernel(k)
    k2 <- kclone(k)
    m <- matrix(0, nrow=length(times), ncol=DATA_SIZE)
    m[, TIME] <- times
    kremove.data(k2)
    kadd.data(k2, m, RV)
    kcalculate(k2)
    d <- kdata(k2)
    return(cbind(d[,TIME], d[,PRED]))
}

kintegrate <- function(k, times, int.method=k$int.method, transits = FALSE, plot=FALSE, print=FALSE, dt=k$dt) {
    .check_kernel(k)
    if (is.nan(k$epoch))
        stop("Set the epoch of the kernel.")
    if (k$nplanets < 1)
        stop("No planets.")
    
    if (length(times) == 1) {
        times <- seq(k$epoch, k$epoch + times, length.out=1000)
    }
    .job <<- sprintf("Integrating for %.2e years", (max(times)-min(times))/365.25)
    .auto <- k$auto
    .int.method <- k$int.method
    
    k$auto <- FALSE
    on.exit(k$auto <- .auto)
    k$int.method <- int.method
    k$dt <- dt
    stopifnot(length(times) >= 2, ! any(is.na(times)))
    nt <- length(times)
    v <- ok_ptr_to_vector(times, nt)
    
    sl <- K_integrateProgress(k$h, v, NULL, k$last.error.code)
    
    if (is.nullptr(sl) || k$last.error.code != K_INTEGRATION_SUCCESS) {
        warning("Error during integration: ", k$last.error, " [", k$last.error.code, "]")
        if (is.nullptr(sl))
            return (NULL)
    }
		
    ret <- list()
    ret$times <- times
    ret$int.method <- k$int.method
    ret$nplanets <- k$nplanets
    ret$rvs <- .gsl_matrix_to_R(ok_get_rvs(sl, nt), free=TRUE)
    els <- .gsl_matrix_to_R(ok_get_els(sl, nt, FALSE), free=TRUE)
    xyz <- .gsl_matrix_to_R(ok_get_xyzs(sl, nt), free=TRUE)

    if (transits) {
        ret$transits <- list()
        for (i in 1:k$nplanets) {
            ret$transits[[i]] <- .gsl_vector_to_R(ok_find_transits(sl, nt, i,
                                                                   k$int.method, 1e-6, NULL, k$last.error.code), free=TRUE)
        }
    }
    
    ret$els <- list()
    ret$xyz <- list()
    mstar <- k$mstar * K2
    
    for (i in 1:k$nplanets) {
        ret$els[[i]] <- els[, (i * ELEMENTS_SIZE + 2) : ((i+1) * ELEMENTS_SIZE + 1)]
        
        mass <- k[i, 'mass'] * MJUP/MSUN * K2
        a <- sapply(ret$els[[i]][, PER], ok_acalc, mstar, mass)
        ret$els[[i]] <- cbind(ret$els[[i]], a, a*(1-ret$els[[i]][, ECC]), ret$times)

        colnames(ret$els[[i]]) <- c(.elements, 'a', 'q', 'time')
        ret$xyz[[i]] <- cbind(xyz[, (i * 7 + 2) : ((i+1) * 7 + 1)], ret$times)
        colnames(ret$xyz[[i]]) <- c("m", "x", "y", "z", "u", "v", "w", "time")
    }

    tmax <- max(ret$times)
    ret$survival.time <- min(sapply(1:k$nplanets, function(n) {
        return(min(ret$times[ret$els[[n]][,'ecc'] > 1 | is.nan(ret$els[[n]][,'period']) | is.nan(ret$els[[n]][,'ecc'])], tmax))
    })) - min(ret$times)
    
    
    ret$xyz[['star']] <- xyz[, 1:7]
    
    colnames(ret$rvs) <- c("time", "rv")
    
    ok_free_systems(sl, nt)
    class(ret) <- "integration"
    k$integration <- ret
    
    if (plot)
        plot(ret)
    
    k$int.method <- .int.method
    k$auto <- .auto

    if (print) {
        print(ret)
        return(invisible(ret))
    } else 
        return(ret)
}

kel2x <- function(mu=1, q=1, e=0, i=0, p=0, n=0, l=0) {
    x <- numeric(1)
    y <- numeric(1)
    z <- numeric(1)
    u <- numeric(1)
    v <- numeric(1)
    w <- numeric(1)
    e <- c(q, e, i, p, n, l)
    
    mco_el2x__(mu, e[1], e[2], e[3], e[4], e[5], e[6],
               x, y, z, u, v, w)
    x <- c(x[1], y[1], z[1], u[1], v[1], w[1])
    names(x) <- c('x', 'y', 'z', 'u', 'v', 'w')
    return (x)
}

kx2el <- function(mu = 1, x=1, y=0, z=0, u=0, v=1, w=0) {
    q <- numeric(1)
    e <- numeric(1)
    i <- numeric(1)
    p <- numeric(1)
    n <- numeric(1)
    l <- numeric(1)
    x <- c(x, y, z, u, v, w)
    mco_x2el__(mu, x[1], x[2], x[3], x[4], x[5], x[6],
               q, e, i, p, n, l)
    q <- c(q, e, i, p, n, l)
    names(q) <- c('q', 'e', 'i', 'p', 'n', 'l')
    return (q)
}



.kl.stats.names <- c("bestfit", "median", "mad")
.klnew <- function(klptr, k, type="", desc="", flags=kflags(k, 'par')) {
    stopifnot(! is.nullptr(klptr))
    np <- KL_getNplanets(klptr)
    size <- KL_getSize(klptr)
    cols <- np * ALL_ELEMENTS_SIZE + PARAMS_SIZE + 3
    
    v <- numeric(size * cols)
    KL_to_ptr(klptr, v)
    
    m <- matrix(v, nrow = size, ncol = cols, byrow = TRUE)

    b <- list()
    
    .els <- kallels(k)
    .els.med <- .gsl_matrix_to_R(KL_getElementsStats(klptr, K_STAT_MEDIAN), free = TRUE)
    .els.mad <- .gsl_matrix_to_R(KL_getElementsStats(klptr, K_STAT_MAD), free = TRUE)
    
    stats <- list()
    if (np > 0) {
        for (i in 1:np) {
            range <- ((i-1) * ALL_ELEMENTS_SIZE + 1) : (i* ALL_ELEMENTS_SIZE)

            b[[i]] <- m[, range, drop = FALSE]

            colnames(b[[i]]) <- .allelements
            
            stats[[i]] <- matrix(ncol = 3, nrow = ALL_ELEMENTS_SIZE)
            stats[[i]][, 1] <- .els[i, , drop = FALSE]
            
            stats[[i]][, 2] <- .els.med[i+1, , drop = FALSE]
            stats[[i]][, 3] <- .els.mad[i+1, , drop = FALSE]
            rownames(stats[[i]]) <- .allelements	
            colnames(stats[[i]]) <- .kl.stats.names			
        }
    }

    .pars <- kpars(k)
    
    b$nplanets <- np
    b$size <- size
    b$stats <- stats
    b$nsets <- k$nsets

    b$params <- m[, (np*ALL_ELEMENTS_SIZE + 1): (ncol(m)-3)]
    b$params.stats <- matrix(ncol = 3, nrow = PARAMS_SIZE)
    
    
    for (i in 1:PARAMS_SIZE) {
        b$params.stats[i, 1] <- .pars[i]
        b$params.stats[i, 2] <- median(b$params[, i])
        b$params.stats[i, 3] <- mad(b$params[, i])
    }
    
    colnames(b$params) <- .params
    rownames(b$params.stats) <- .params	
    colnames(b$params.stats) <- .kl.stats.names
    
    b$chi2 <- m[, ncol(m)-2]
    b$merit <- m[, ncol(m)-2]
    b$prior <- m[, ncol(m)-1]
    b$likelihood <- m[, ncol(m)]
    
    b$fit.els <- kallels(k)
    b$fit.params <- kpars(k)
    b$element.type <- k$element.type
    b$length <- b$size
    b$type <- type
    b$desc <- desc
    b$par.flags <- flags
    class(b) <- "error.est"
    
    b$.matrix <- m
    
    KL_free(klptr)
    k$errors <- b
    return(b)
}

# Forces printing of string
print.ktor <- function(s) {
    cat(s)
}

ktor <- function(k) {
    name <- deparse(substitute(k))
    .check_kernel(k)

    tr <- sprintf("%s <- knew()", name)
    for (prop in c("mstar", "int.method", "epoch", "element.type"))  {
        tr <- c(tr, sprintf("%s$%s <- %s", name, prop, kget(k, prop)))
    }

    if (k$nplanets > 0)
        for (i in 1:k$nplanets)
            tr <- c(tr, sprintf("%s[] <- c(period=%s, mass=%s, ma=%s, ecc=%s, lop=%s, inc=%s, node=%s)", name, k[i, 'period'], k[i, 'mass'], k[i, 'ma'], k[i, 'ecc'], k[i, 'lop'], k[i, 'inc'], k[i, 'node']))

    tr <- paste(tr, collapse='\n')
    class(tr) <- c("ktor", "character")
    return(tr)
}


print.error.est <- function(e, all.pars = FALSE) {
    cat(sprintf("# %s, %s\n", e$type, e$desc))
    if (e$nplanets > 0)
        for (i in 1:e$nplanets) {
            cat(sprintf("# Planet %d ($stats[%d])\n", i, i))
            print(e$stats[[i]])
            cat("\n")
        }
    cat("# Parameters ($params.stats)\n")
    if (! all.pars) {
        if (is.null(e$par.flags))
            print(e$params.stats[1:e$nsets, , drop=F])
        else 
            print(e$params.stats[bitAnd(e$par.flags, MINIMIZE) == MINIMIZE, , drop=F ])
        
        cat("\n(Use print(name, all.pars = TRUE) to see all parameters)\n")
    } else 
        print(e$params.stats[, ])
}


print.kernel <- function(k, all.pars = FALSE) {
    .check_kernel(k)
    cat(sprintf("# %d planets, %d datasets, %d data points\n# intMethod: %d, epoch: %f\n",
                k$nplanets, k$nsets, k$ndata, k$int.method, k$epoch))
    cat(sprintf("# Chi^2 = %.2f, RMS = %.2f, jitter = %.2f\n",
                k$chi2, k$rms, k$jitter))
		
    if (k$nplanets > 0) {
        cat("\n# Orbital elements\n")
        print(kallels(k))
        cat("\n# Flags\n")
        print(kflags(k, what="els", type='human'))
    }
    cat("\n# Parameters\n")
    if (! all.pars)
        print(kpars(k)[1:10])
    else
        print(kpars(k))
		
    cat("\n# Flags\n")
    if (! all.pars)
        print(kflags(k, what='par', type='human')[1:10])	
    else
        print(kflags(k, what='par', type='human'))	
    if (! is.null(k$errors)) {
        cat("\n# Last error estimation run results:\n")
        print(k$errors)
        cat("\n")
    }
}

print.integration <- function(int) {

    cat(sprintf("# Integration: %.2f days [%.2f years], survival time: %.2f years, int.method: %d\n",
                max(int$times)-min(int$times),
                (max(int$times)-min(int$times))/365.25,
                int$survival.time/365.25,
                int$int.method))
    
    for (i in 1:int$nplanets) {
        cat("# Planet ", i, "\n")
        cat(sprintf("Pmin: %.2f, Pmax: %.2f, dev(P)/med(P): %.2e\n",
                    min(int$els[[i]][, 'period']), 
                    max(int$els[[i]][, 'period']), 
                    sd(int$els[[i]][, 'period'])/
                        median(int$els[[i]][, 'period'])))
        cat(sprintf("emin: %.2f, emax: %.2f, dev(e)/med(e): %.2e\n",
                    min(int$els[[i]][, 'ecc']), 
                    max(int$els[[i]][, 'ecc']), 
                    sd(int$els[[i]][, 'ecc'])/
                        median(int$els[[i]][, 'ecc'])))
        
        if (max(abs(diff(int$els[[i]][, 'inc']))) > 1e-5) {
            cat(sprintf("imin: %.2f, imax: %.2f, dev(i)/med(i): %.2e\n",
                        min(int$els[[i]][, 'inc']), 
                        max(int$els[[i]][, 'inc']), 
                        sd(int$els[[i]][, 'inc'])/
                            median(int$els[[i]][, 'inc'])))
        }
    }
    
    cat("\n?To access the elements evolution for planet i, use var$els[[i]]\n")
    cat("?To access the cartesian evolution, use var$xyz[[i]]\n")
    if (!is.null(int))
        cat("?To access the transit times, use var$transits[[i]]\n")
    cat("?To plot the results, use plot(var, what=c('a', 'ecc', 'lop', ...)\n")
}

save.systemic <- function(file) {
    for (kn in ls(envir=globalenv())) {
        ka <- get(kn, envir=globalenv())
        
        if ("kernel" %in% class(ka)) {
            fn <- tempfile()
            p <- fopen(fn, "w")
            K_save(ka$h, p)
            fclose(p)
            con <- file(fn)
            ka$representation <- readLines(con)
            close(con)
            assign(kn, ka, envir=globalenv())
        }
    }
    
    save(list = ls(envir=globalenv()), file=file)
    return(invisible())
}

load.systemic <- function(file) {
    load(file, envir=globalenv())
    
    for (kn in ls(envir=globalenv())) {
        k <- get(kn, envir=globalenv())
        if (class(k) == "kernel") {
            fn <- tempfile()
            con <- file(fn, "w")
            writeLines(k$representation, con)
            close(con)

            k2 <- kload(fn)
            k[['h']] <- k2[['h']]
            k$min.func <- k$min.func
            assign(kn, k, envir=globalenv())
        }
    }
    invisible()
}



## [1] Creating/opening/saving kernel objects
## [2] Setting fit parameters (planets, offsets, etc.)
## [5] Loading and manipulating data
## [6] Fit parameters (chi^2, jitter, etc.)
## [3] Minimization
## [4] Error estimation
## [7] Periodogram
## [8] Constants
## [9] Other
