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

        return(K_getMinValue(k$h))
    }

    if (minimize.function == 'default')
        minimize.function = kminimize.default

    Domain <- kminimize.domain(k, log.period=log.period, log.mass=log.mass)
        
    v <- genoud(minimize.function, k$nrpars, Domains=Domain,  max.generations=max.generations, ...)
    
    return(invisible(v))
}

kminimize.de <- function(k, minimize.function='default', log.period=TRUE, log.mass=TRUE, population=50*k$nrpars,
                         max.iterations=1000, F = 0.5, CR = 1, plot.status=NULL,
                         wait=10, check.function=NULL, mc.cores=getOption("mc.cores", 1), ...) {
    .check_kernel(k)
    stopifnot(k$nplanets > 0)
    indices <- kminimized.indices(k)
    per.indices <- indices[1,] != -1 & indices[2,] == PER
    mass.indices <- indices[1,] != -1 & indices[2,] == MASS
    
    Domain <- kminimize.domain(k, log.period=log.period, log.mass=log.mass)
    print(Domain)
    print(indices)

    kminimize.default <- function(v, ..., clone=mc.cores != 1) {
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

        return(K_getMinValue(k2$h))
    }

   

    if (minimize.function=='default')
        minimize.function <- kminimize.default
    
    x <- lapply(1:population, function(...) {
        return(runif(k$nrpars, Domain[,1], Domain[,2]))
    })
    
    f <- unlist(mclapply(x, minimize.function, mc.cores=mc.cores))
    x <- lapply(1:length(x), function(i) c(x[[i]], f[i]))
    
    on.exit({ minimize.function(x[[which.min(f)]], clone=FALSE); kupdate(k, calculate=TRUE) })
    
    
    print(summary(f))
    if (!is.null(plot.status)) {
        pm <- sapply(x, function(p) return(p[plot.status]))
        xrange <- range(pm[1,])
        yrange <- range(pm[2,])
    }

    
    for (reps in 1:max.iterations) {
        x <- mclapply(1:length(x), function(me) {
            
            repeat {
                idx <- sample(1:length(x), 3)
                if (! any(idx == me)) break
            }

            R <- sample(1:k$nrpars, 1)
            v <- sapply(1:k$nrpars, function(i) {
                if (runif(1) < CR || i == R)
                    return(x[[idx[1]]][i] + F * (x[[idx[2]]][i] - x[[idx[3]]][i]))
                else
                    return(x[[me]][i])
            })

            if (!in.domain(v, Domain))
                return(x[[me]])
            
            fnew <- minimize.function(v, check.function=check.function)
            if (is.nan(fnew))
                return(x[[me]])

            if (is.nan(f[me]) || fnew < f[me]) {
                return(c(v, fnew))
            }
            else
                return(x[[me]])
        }, mc.cores=mc.cores)

        f <- sapply(x, function(v) v[length(v)])
        
        print(summary(f))
        if (!is.null(plot.status)) {
            pm <- sapply(x, function(p) return(p[plot.status]))
            col <- rep('black', length(x))
            col[which.min(f)] <- 'red'
            plot(pm[1,], pm[2,], col=col, pch=19)
        }
        print(x[[which.min(f)]])
    }

    return(x[[which.min(f)]])
}
