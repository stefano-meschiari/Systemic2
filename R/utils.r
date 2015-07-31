str_rep <- function (ch, times)  
    str_c(rep(ch, times), collapse = "")

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

kauto.steps <- function(k1, delta.func=0.1, max.iters=10, verbose=TRUE) {
    k <- kclone(k1)
    a <- k['minimized']
    p <- k['minimized.steps']
    
    chi2.0 <- K_getMinValue(k1$h)
    for (i in 1:length(a)) {
        iters <- 0
        for (j in 1:max.iters) {
            k['minimized'][i] <- a[i] + p[i]
            kcalculate(k)
            chi2 <- K_getMinValue(k$h)
            dchi <- abs(chi2 - chi2.0)
            iters <- iters + 1
            p[i] <- 0.5 * p[i] * max(min(1 + delta.func/dchi, 10), 0.1)
            if (verbose)
                cat(sprintf("[%d] v = %e, dfunc = %e, step = %e\n", i, a[i], dchi, p[i]))
        }
        k['minimized'] <- a
    }
    idx <- kminimized.indices(k)
    for (i in 1:ncol(idx)) {
        kstep(k1, idx[1, i], idx[2, i]) <- p[i]
    }
    
}

exoplanet.archive <- function(name) {
    if (is.na(name)) {
        stop("Specify the star name (either assign it to a kernel using k$starname <- NAME, or pass it to this function directly)")
    }
    message(str_c("Using star name = ", name))
   browseURL(sprintf("http://exoplanetarchive.ipac.caltech.edu/cgi-bin/ExoOverview/nph-ExoOverview?objname=%s&type=PLANET%%20HOST&label&aliases&exo&iden&orb&ppar&tran&disc&ospar", name))
}

simbad <- function(name=NA) {
    if (is.na(name)) {
        stop("Specify the star name (either assign it to a kernel using k$starname <- 'NAME', or pass the star name to this function directly)")
    }
    message(str_c("Using star name = ", name))
    browseURL(sprintf("http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%s", name))
}

exoplanets.org <- function(name=NA) {
    if (is.na(name)) {
        stop("Specify the star name (either assign it to a kernel using k$starname <- 'NAME', or pass the star name to this function directly")
    }

    cat("Locating star in database...\n")
    exoarchive <- sprintf("http://exoplanetarchive.ipac.caltech.edu/cgi-bin/ExoOverview/nph-ExoOverview?objname=%s&type=PLANET%%20HOST&label&aliases&exo&iden&orb&ppar&tran&disc&ospar", name)


    a <- paste0(collapse='\n', readLines(url(exoarchive)))
    url <- str_match(a, '(http://exoplanets.org/.*?) ')[1,2]
    if (is.na(url)) {
        stop("Could not locate correct link to exoplanets.org")
    } else {
        name <- str_match(url, 'detail/(.*?)_b')[1,2]
        cat('(called ', name, ' on exoplanets.org.)\n')
        for (l in letters[3:20]) {
            url <- paste0(url, ',', name, '_', l)
        }
        browseURL(url)
    }
}

exoplanet.eu <- function(name) {
    if (is.na(name)) {
        stop("Specify the star name (either assign it to a kernel using k$starname <- 'NAME', or pass the star name to this function directly")
    }

    cat("Locating star in database...\n")
    exoarchive <- sprintf("http://exoplanetarchive.ipac.caltech.edu/cgi-bin/ExoOverview/nph-ExoOverview?objname=%s&type=PLANET%%20HOST&label&aliases&exo&iden&orb&ppar&tran&disc&ospar", name)


    a <- paste0(collapse='\n', readLines(url(exoarchive)))
    url <- str_match(a, '(http://exoplanet.eu/.*?) ')[1,2]
    if (is.na(url)) {
        stop("Could not locate correct link to exoplanet.eu")
    } else {
        name <- str_match(url, 'detail/(.*?)_b')[1,2]
        cat('(called ', name, ' on exoplanets.eu.)\n')
        browseURL(url)
    }
}
    

kclone.boot <- function(k) {
  k2 <- kclone(k)
  d <- kdata(k)
  
  while (TRUE) {
    m <- d[sample(1:nrow(d), size=nrow(d), replace=TRUE), ]
    if (length(unique(m[, SET])) == k$nsets)
      break
  }
  kdata(k2) <- m
  kcalculate(k2)
  return(k2)
}

kboot.cv <- function(k, training.len=0.9) {
  k2 <- kclone(k)
  d <- kdata(k)
  
  while (TRUE) {
    training.indices <- sample(1:nrow(d), size=ceil(training.len*nrow(d)))
    if (length(unique(d[training.indices, SET])) == k$nsets)
      break
  }

  d2 <- d
  d2[-training.indices, ERR] <- -1
  kdata(k2) <- d2
  kminimize(k2)

  d[training.indices, ERR] <- -1
  kdata(k2) <- d
  kcalculate(k2)
  k2$loglik
}
