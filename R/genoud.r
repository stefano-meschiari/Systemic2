
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

print.de <- function(x) {
    cat(sprintf("%-30s: %e\n", "Final value", x$min.f))
    cat(sprintf("%-30s: %.0f\n", "Total # of iterations", x$niters))
    cat(sprintf("%-30s: %.0f\n", "Total # of evaluations", x$niters * x$population))
    cat(sprintf("%-30s: %s\n", "Reason for stopping:", x$reason))                
}

kminimize.de <- function(k, minimize.function='default', domain=NULL, log.period=TRUE, log.mass=TRUE, population=10*k$nrpars,
                         max.iterations=1000, nochange = 500, target.val=NA, target.val.tol = 0.05, F = 'dither', CR = 0.9, plot=NULL, 
                         check.function=NULL, accept.function=NULL, mc.cores=getOption("mc.cores", 1), save=NULL, save.trials=NULL, save.perf=NULL, plot.perf=FALSE,
                          ..., initial.pop=NULL, use.k=FALSE, type='rand', break.if=NULL, silent=TRUE, log='y', function.best=NULL) {
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
    vinitial <- k['minimized']
    reason <- sprintf("Reached max.iterations [%d]", max.iterations)

    if (is.null(domain))
        Domain <- kminimize.domain(k, log.period=log.period, log.mass=log.mass)
    else
        Domain <- domain
    
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
        if (!is.null(accept.function))
            if (!accept.function(v))
                return(NaN)
        if (!is.null(check.function)) 
            if (!check.function(k2))
                return(NaN);
        
        return(minimize.function(k2, v))
    }

    if (is.null(initial.pop)) {
        if (!silent)
            cat("Initial population is NULL, trying to find suitable population...\n")
        
        x <- lapply(1:population, function(i, ...) {
            while (TRUE) {
                vals <- runif(k$nrpars, Domain[,1], Domain[,2])
                vals2 <- vals
                if (log.period)
                    vals2[per.indices] <- 10^vals2[per.indices]
                if (log.mass)
                    vals2[mass.indices] <- 10^vals2[mass.indices]
                if (!is.null(accept.function)) {
                    if (!accept.function(vals2)) {
                        next
                    } else {
                        break
                    }
                } else
                    break
            }
            return(vals)
        })
        if (!silent)
            cat("Done.\n")
        initial.pop <- x
    } else {
        x <- initial.pop
    }

    if (use.k) {
        x[[1]] <- vinitial
        if (log.period)
            x[[1]][per.indices] <- log10(x[[1]][per.indices])
        if (log.mass)
            x[[1]][mass.indices] <- log10(x[[1]][mass.indices])
        
    }

    if (mc.cores == 1)
        apply.function <- lapply
    else
        apply.function <- mclapply

    f <- unlist(apply.function(x, set.values, mc.cores=mc.cores))
    x <- apply.function(1:length(x), function(i) c(x[[i]], f[i]))
    if (!silent)
        cat("Starting...\n")
    initial.pop <- x

    on.exit({ set.values(x[[which.min(f)]], clone=FALSE); kupdate(k, calculate=TRUE) })

    dither <- F=='dither'

    iters <- c(0, min(f))
    no.imp <- 0
    min.f <- 1
    perf <- c()
    
    for (reps in 1:max.iterations) {
        last.min.f <- min.f
        if (! any(is.nan(f))) 
            min.f <- which.min(f)
        else
            min.f <- 1
        
        if (f[min.f] != last.min.f)
            no.imp <- 0
        else
            no.imp <- no.imp + 1
        
        last.min.f <- min(f)

        ## if ((max(f) - min(f))/min(f) < 0.5) {
        ##     cat(sprintf("Introducing new element [%e, %e]\n", min(f), max(f)))
        ##     idx <- sample(1:population, 1)
        ##     x[[idx]] <- initial.pop[[idx]]
        ##     f[idx] <- x[[idx]][length(x)]
        ##     print(x[[idx]])
        ##     print(f[idx])
        ## }

        
        ## if (reps %% 1000 == 0) {
        ##     set.values(x[[min.f]], clone=FALSE)
        ##     kminimize(k)
        ##     i <- sample(1:population, 1)
        ##     f[i] <- minimize.function(k, v)
        ##     x[[i]] <- c(k['minimized'], f[i])
        ## }

        if (no.imp > nochange) {
            reason <- sprintf("No improvement in the last %d of %d total iterations, assuming convergence.", nochange, nrow(iters))
            break
        }


        if (!is.na(target.val) && ((abs((target.val - min(f))/target.val) < target.val.tol) || (target.val > min(f)))) {
            reason <- sprintf("Min.f [%e] is within target.val.tol [%e] of target.val [%e].", min.f, target.val.tol, target.val)
            break
        }
        
        if (dither)
            F <- runif(1, 0.5, 1)
        x <- apply.function(1:length(x), function(me, ...) {
            
            repeat {
                idx <- sample(1:length(x), 3)
                if (! any(idx == me)) break
            }

            R <- sample(1:k$nrpars, 1)

                
            v <- sapply(1:k$nrpars, function(i) {
                
                if (runif(1) < CR || i == R) {
                    if (type == 'rand') {
                        return(x[[idx[1]]][i] + F * (x[[idx[2]]][i] - x[[idx[3]]][i]))
                    } else if (type == 'best') {
                        return(x[[min.f]][i] + F * (x[[idx[2]]][i] - x[[idx[3]]][i]))
                    } else if (type == "current-to-best") {
                        return(x[[me]][i] + F * (x[[min.f]][i] - x[[me]][i]) + F * (x[[idx[2]]][i] - x[[idx[3]]][i]))
                    }
                }
                else
                    return(x[[me]][i])
            })

            #if (!in.domain(v, Domain))
            #    return(x[[me]])
            
            fnew <- set.values(v, check.function=check.function)
            if (is.nan(fnew))
                return(x[[me]])
            
            if (is.nan(f[me]) || (fnew < f[me])) {
                return(c(v, fnew))
            }
            else
                return(x[[me]])
        }, mc.cores=mc.cores)

        f <- sapply(x, function(v) v[length(v)])

        iters <- rbind(iters, c(reps, min(f)))
        if (!is.null(save.trials)) {
            de.pars <- list(pop=x, iters=iters)
            save(de.pars, file=save.trials)
        }

        if (!silent) {
            print(summary(f))
        }
        
        if (!is.null(save.perf)) {
            p <- cbind(sapply(x, function(v) v[save.perf]), f)
            perf <- rbind(perf, p)
            if (plot.perf) {
                if (nrow(perf) == nrow(p))
                    plot(perf[,1], perf[,2], log='y', ylim=c(1, max(perf[,2])))
                else
                    points(p[,1], p[,2])
            }
        }
        if (!is.null(plot)) {
            pm <- sapply(x, function(p) return(p[plot]))
            col <- rep('black', length(x))
            if (!all(is.nan(f)))
                col[which.min(f)] <- 'red'
            plot(pm[1,], pm[2,], col=col, pch=19, log=log)
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
            if (!silent)
                print(c(reps, x[[which.min(f)]]))
            if (!is.null(save)) {
                set.values(x[[which.min(f)]], clone=FALSE)
                ksave(k, save)
            }
        }
       
        else {
            if (!silent)
                print("All solutions are NaN for now (either the integrator returned NaN for all trial solutions, or no trial solutions satisfy check.function)")
            
        }

        if (!is.null(break.if)) {
            set.values(x[[which.min(f)]], clone=FALSE)            
            kupdate(k)
            if (break.if(k)) {
                break
            }
        }
        if (!is.null(function.best)) {
            set.values(x[[which.min(f)]], clone=FALSE)            
            kupdate(k)
            function.best(k)
        }
    }

    set.values(x[[which.min(f)]], clone=FALSE)            
    kupdate(k)

    a <- list(best=x[[which.min(f)]], min.f=min(f), perf=perf, niters=nrow(iters), iters=iters, reason=reason, population=population)
    class(a) <- 'de'
    
    return(a)
}
