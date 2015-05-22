write.f <- function(m, file="", col.names=TRUE, format="%18.10e", sformat="%18s", comments=NULL) {
    f <- file(file, "w")

    if (!isOpen(f))
        stop(sprintf("Could not open %s", file))

    on.exit(close(f))
    if (!is.null(comments)) 
        writeLines(sprintf("# %s", comments))
    
    if (is.matrix(m) || is.data.frame(m)) {
        if (col.names && !is.null(colnames(m))) {
            cat('# ', file=f)
            writeLines(Reduce(paste, sprintf(sformat, colnames(m))), con=f)
        }
        for (r in 1:nrow(m)) 
            writeLines(Reduce(paste, sprintf(format, m[r, ])), con=f)
    } else if (is.vector(m)) {
        writeLines(sprintf(format, m), con=f)        
    } else {
        stop(sprintf("I don't know how to write out object of class %s; try to use the write, dump or save functions instead"))
    }

}


write.fmatrix <- function(m, file="", col.names=TRUE, format="%18.10e", sformat="%18s") {
    warning("write.fmatrix is deprecated, using write.f instead")
    write.f(m, file, col.names, format, sformat)
}

longest.prefix <- function(c) {
    stopifnot(class(c) == "character")
    if (length(c) == 0)
        return(NULL)
    else if (length(c) == 1)
        return(c[1])
    else
        return(paste(Reduce(function(a, b) {
            if (is.null(a) || is.null(b))
                return(NULL)
            if (class(a) == "character" && length(a) == 1)
                a <- strsplit(a, '')[[1]]
            b <- strsplit(b, '')[[1]]
            
            ret <- c()

            for (i in 1:(min(length(a), length(b)))) {
                if (a[i] == b[i])
                    ret <- c(ret, a[i])
                else
                    break
            }

            return(ret)
        }, c), collapse=''))
}

kscan.error.est <- function(k, obj, sample.length=obj$length, time=1e3*365.25, samples=NULL, int.method=SWIFTRMVS, dt=min(k[,'period'])*0.01, criterion=NULL, print=TRUE, progress=NULL, save=NULL) {

    stable <- c()
    
    times <- c()
    if (is.null(samples))
        samples <- sample(1:obj$length, sample.length)

    k <- kclone(k)
    k$silent <- TRUE
    k$errors <- NULL
    for (i in samples) {
        for (j in 1:k$nplanets)
            k[j,] <- obj[[j]][i, ]
        io <- kintegrate(k, time, int.method=int.method, dt=dt)
        if (!is.null(criterion))
            st <- criterion(k, io)
        else 
            st <- !(io$survival.time < time)

        stable <- c(stable, st)
        if (exists('gui.progress'))
            gui.progress(length(stable), sample.length, io$survival.time/365.25, sprintf("Integrating..."))
                                              
        times <- c(times, io$survival.time)
        if (!is.null(progress))
            progress(length(stable), sample.length, stable, times)

    }

    if (print)
        cat(sprintf("%d/%d samples are unstable [fraction = %.2e]", sum(!stable), length(stable),
                    sum(!stable)/length(stable)))

    kscan.ret <- list(samples=samples, stable=stable, survival.times=times)
    if (!is.null(save)) {
        save(kscan.ret, file=save)
    }
    
    return(invisible(kscan.ret))
    
}

kclone.from <- function(k, obj, index) {
    .check_kernel(k)
    stopifnot(class(obj) == "error.est")

    k <- kclone(k)

    for (j in 1:k$nplanets) 
        k[j,] <- obj[[j]][index, ]
    
    k['par', ] <- obj$params[index, ]
 
    return(k)
}

kauto.steps <- function(k1, delta.chi=0.1, max.iters=10, verbose=TRUE) {
    k <- kclone(k1)
    a <- k['minimized']
    chi2.0 <- k$chi2
    for (i in 1:length(a)) {
        for (j in 1:max.iters) {
            k['minimized'][i] <- a[i] + k['minimized.steps']
            kcalculate(k)
            dchi <- abs(k$chi2 - chi2.0)
            iters <- iters + 1
            k['minimized.steps'][i] <- 0.5 * k['minimized.steps'][i] * (1 + delta.chi/dchi)
            if (verbose)
                cat(sprintf("[%d] dchi = %e, step = %e\n", k['minimized.steps'][i]))
        }
        k['minimized'] <- a
    }
}
