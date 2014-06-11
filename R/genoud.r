default.ranges <- list()
default.ranges[[PER]] = c(0.1, 1e4)
default.ranges[[MASS]] = c(1e-4, 1e3)
default.ranges[[MA]] = c(0, 360) * 100
default.ranges[[ECC]] = c(0, 0.95)
default.ranges[[LOP]] = c(0, 360) * 100
default.ranges[[INC]] = c(0, 360) * 100
default.ranges[[NODE]] = c(0, 360) * 100
default.ranges[['par']] = c(-1e3, 1e3)

kminimize.domain <- function(k, log.period=FALSE, log.mass=FALSE) {
    indices <- kminimized.indices(k)
    
    return(t(sapply(1:k$nrpars, function(i) {
        r <- krange(k, indices[1, i], indices[2, i])

        if (log.period && indices[1, i] != -1 && indices[2, i] == PER)
            r <- log10(r)
        if (log.mass && indices[1, i] != -1 && indices[2, i] == MASS)
            r <- log10(r)
        
        
        if (is.nan(r[1]))
            if (indices[1, i] == -1) r[1] <- default.ranges[['par']][1] else r[1] <- default.ranges[[indices[2,i]]][1]
        if (is.nan(r[2]))
            if (indices[1, i] == -1) r[2] <- default.ranges[['par']][2] else r[2] <- default.ranges[[indices[2,i]]][2]
        
        return(r)
    })))
}

in.domain <- function(x, Domain, indices) all(x >= Domain[,1] & x <= Domain[,2])


kminimize.genoud <- function(k, minimize.function='default', log.period=TRUE, log.mass=TRUE, max.generations=100, ...) {
    .require.library('rgenoud')
    .check_kernel(k)
    stopifnot(k$nplanets > 0)
    indices <- kminimized.indices(k)
    per.indices <- indices[1,] != -1 & indices[2,] == PER
    mass.indices <- indices[1,] != -1 & indices[2,] == MASS
    on.exit(kupdate(k, calculate=TRUE))
    
    kminimize.default <- function(v, ...) {
        if (log.period)
            v[per.indices] = 10^v[per.indices]
        if (log.mass)
            v[mass.indices] = 10^v[mass.indices]

        K_setMinimizedValues(k$h, v)
        K_calculate(k$h)
        if (minimize.function == 'default')
            return(K_getMinValue(k$h))
        else
            minimize.function(k, v)
    }

    Domain <- kminimize.domain(k, log.period=log.period, log.mass=log.mass)
        
    v <- genoud(minimize.function, k$nrpars, Domains=Domain,  max.generations=max.generations, ...)
    
    return(invisible(v))
}

kminimize.de <- function(k, minimize.function='default', log.period=TRUE, log.mass=TRUE, population=min(40, 10*k$nrpars),
                         max.iterations=1000, F = 'dither', CR = 0.9, plot=NULL,
                         wait=10, check.function=NULL, mc.cores=getOption("mc.cores", 1), save=NULL, save.trials=NULL, min.f.spread=1e-3, ...) {
    .check_kernel(k)
    stopifnot(k$nplanets > 0)
    if (!is.null(plot)) {
        if (length(plot) %% 2 != 0)
            stop("Use pairs of parameter indices for the plot parameter, e.g. kminimize.de(..., plot=c(1, 2, 5, 6))")
        par(mfrow=c(length(plot)/2, 1))
    } 
    indices <- kminimized.indices(k)
    per.indices <- indices[1,] != -1 & indices[2,] == PER
    mass.indices <- indices[1,] != -1 & indices[2,] == MASS
    
    Domain <- kminimize.domain(k, log.period=log.period, log.mass=log.mass)
    if (is.character(minimize.function) && minimize.function == 'default')
        minimize.function <- function(k, ...) {
            return(K_getMinValue(k$h))
        }
    
    set.values <- function(v, ..., clone=mc.cores != 1) {
        if (clone)
            k2 <- kclone(k)
        else
            k2 <- k
        
        if (log.period)
            v[per.indices] = 10^v[per.indices]
        if (log.mass)
            v[mass.indices] = 10^v[mass.indices]

        K_setMinimizedValues(k2$h, v)
        K_calculate(k2$h)
        if (!is.null(check.function)) 
            if (!check.function(k2))
                return(NaN);

        return(minimize.function(k2, v))
    }

    
    x <- lapply(1:population, function(...) {
        return(runif(k$nrpars, Domain[,1], Domain[,2]))
    })
    
    f <- unlist(mclapply(x, set.values, mc.cores=mc.cores))
    print(f)
    x <- lapply(1:length(x), function(i) c(x[[i]], f[i]))

    on.exit({ set.values(x[[which.min(f)]], clone=FALSE); kupdate(k, calculate=TRUE) })

    dither <- F=='dither'
    print(summary(f))

    iters <- c(0, min(f))
    for (reps in 1:max.iterations) {
        if (! any(is.nan(f))) 
            min.f <- which.min(f)
        else
            min.f <- 1
        
        x <- mclapply(1:length(x), function(me, ...) {
            
            repeat {
                idx <- sample(1:length(x), 3)
                if (! any(idx == me)) break
            }

            R <- sample(1:k$nrpars, 1)
            if (dither)
                    F <- runif(1, 0.5, 1)
                
            v <- sapply(1:k$nrpars, function(i) {
                
                if (runif(1) < CR || i == R)
                    return(x[[idx[1]]][i] + F * (x[[idx[2]]][i] - x[[idx[3]]][i]))
                else
                    return(x[[me]][i])
            })

            if (!in.domain(v, Domain))
                return(x[[me]])
            
            fnew <- set.values(v, check.function=check.function)
            if (is.nan(fnew))
                return(x[[me]])
            
            ## if (me != min.f && !is.nan(f[me]) && fnew > f[me]) {
            ##     spread <- 1/sd(f, na.rm=TRUE)
            ##     r <- (exp(-(fnew-f[me])/spread) > runif(1))

            ##     if (r) {
            ##         print(c(spread, fnew-f[me]))
            ##         stop()
            ##     }
            ## } else
            ##     r <- FALSE
            
            r <- FALSE
            if (is.nan(f[me]) || (fnew < f[me]) || r) {
                return(c(v, fnew))
            }
            else
                return(x[[me]])
        }, mc.cores=mc.cores)

        f <- sapply(x, function(v) v[length(v)])

        iters <- rbind(reps, min.f)
        if (!is.null(save.trials)) {
            de.pars <- list(pop=x, iters=iters)
            save(de.pars, file=save.trials)
        }
        print(summary(f))
        if (!is.null(plot)) {
            pm <- sapply(x, function(p) return(p[plot]))
            col <- rep('black', length(x))
            if (!all(is.nan(f)))
                col[which.min(f)] <- 'red'
            plot(pm[1,], pm[2,], col=col, pch=19)
            if (length(plot) > 2)
                for (i in seq(3, length(plot), 2))
                    plot(pm[i,], pm[i+1,], col=col, pch=19)
        }
        if (!all(is.nan(f))) {
            a <- x[[which.min(f)]]
            if (log.period)
                a[per.indices] <- 10^a[per.indices]
            if (log.mass)
                a[mass.indices] <- 10^a[mass.indices]
            print(x[[which.min(f)]])
            if (!is.null(save)) {
                set.values(x[[which.min(f)]], clone=FALSE)
                ksave(k, save)
            }
        }
       
        else {
            print("All solutions are NaN for now (either the integrator returned NaN for all trial solutions, or no trial solutions satisfy check.function)")
            
        }
    }

    return(x[[which.min(f)]])
}
