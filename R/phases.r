phases <- function(k, from, to, phases=c(0, 0.25, 0.5, 0.75), save=NULL, print=TRUE,  plIndex = 1) {
    .check_kernel(k)
    if (to < from)
        stop("Incorrect time interval passed (*from* should be smaller than *to*)") 
    k1 <- kclone(k)

    P <- k[plIndex, 'period']
    M <- seq(0, 360, length.out=2e4)

    phases.deg <- phases * 360 + 90 - k[1, 'lop']
    phases.deg[phases.deg < 0] <- phases.deg[phases.deg < 0] + 360
    times <- c()
    for (phase in phases.deg) {
        k1[plIndex, 'trueanomaly'] <- phase
        
        t <- k1[plIndex, 'ma'] * P/360 + k[1, 'tperi']
        times <- c(times, t)
    }
            
    out <- list()

    out$from <- from
    out$to <- to
    out$phases <- list()
    
    
    for (i in 1:length(phases)) {
        phase <- phases[i]
        t.phase <- c()
        d <- floor((from-times[i])/P) * P - P
        t <- times[i] + d
        while (t <= to) {
            if (t >= from) {
                t.phase <- c(t.phase, t)
            }
            t <- t + P
        }
        out$phases[[paste(phase)]] <- t.phase
    }

    if (!is.null(save)) {
        cat(file=save, toString(out))
    }

    class(out) <- 'phases'
    return(out)
}

toString.phases <- function(p) {

    a <- c(sprintf("# Between JD = %f and JD = %f:\n", p$from, p$to))

    for (phase in names(p$phases)) {
        a <- c(a, sprintf("# %s\n", phase))
        if (is.null(p$errors)) {
            a <- c(a, sprintf("%10.6f\n", p$phases[[phase]]))
        } else {
            a <- c(a, sprintf("%10.6f\t%e\n", p$phases[[phase]], p$errors[[phase]]))
        }
    }
    return(paste(a, collapse='', sep=''))
}

print.phases <- function(p) {
    cat(toString(p))
}

phases.avg <- function(plist) {
    pout <- plist[[1]]
    
    for (phase in names(plist[[1]]$phases)) {        
        ts <- plist[[1]]$phases[[phase]]
        pout$phases[[phase]] <- sapply(1:length(ts), function(i) {
            return(mean(sapply(plist, function(p) {
                return(p$phases[[phase]][i])
            })))
        })
        pout$errors[[phase]] <- sapply(1:length(ts), function(i) {
            return(sd(sapply(plist, function(p) {
                return(p$phases[[phase]][i])
            })))
        })
    }

    return(pout)
}
