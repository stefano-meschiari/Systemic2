.peaks.reduce <- function(p, n = 4, peak.tolerance) {
    a <- p[1,,drop=FALSE]
    i <- 2
    while (nrow(a) <= n && i < nrow(p)) {
        if (min(abs(a[,1] - p[i, 1]) / p[i,1]) > peak.tolerance)
            a <- rbind(a, p[i,])
        i <- i + 1
    }
    return(a)
}

kauto <- function(k, trials=4, samples = getOption("systemic.psamples", 50000), pmin = getOption("systemic.pmin", 0.5), pmax = getOption("systemic.pmax", 1e4), fap.limit=1e-3, fit.trends=TRUE, peak.tolerance=0.05) {

    out <- list()
    k <- kclone(k)
    if (k$nplanets > 0) {
        k[] <- NULL
    }
    kselect(k, 'par', 'all')
    kminimize(k)

    
    p <- kperiodogram(k, 'res', samples, pmin, pmax)
    if (min(p[, 'fap']) > fap.limit) {
        break
    }
    print('hi')
    
    p <- kperiodogram.boot(k, 'res', samples, pmin, pmax)
    if (min(p[, 'fap']) > fap.limit) {
        break
    }

    peaks <- attr(p, 'peaks')
    peaks <- peaks[, peaks[,3] < fap.limit]
    peaks <- .peaks.reduce(peaks, trials, peak.tolerance)
    
    
    return(peaks)
}
