
systemic.par <- list()
systemic.par$font.lab <- 2
systemic.par$tck <- 0.02
systemic.par$lwd <- 2

systemic.plot.style <- function() {
    par(systemic.par)
}


plot.kernel <- function(k, type = "rv", wrap=NA, plot.residuals=TRUE, transiting.planet = NA, transiting.per = NA, xlim = NA, breaks=NA, plot.gaussian=TRUE, density=FALSE, ...) {
    .check_kernel(k)
    oldpar <- par(no.readonly=TRUE)
    on.exit(suppressWarnings(par(oldpar)))
    
    par(systemic.par)
    if (is.nan(k$epoch)) {
        stop("No epoch set")
    }
    
    if (type == "rv") {
        rows <- if (plot.residuals) 2 else 1
        
        par(mfrow=c(rows, 1), mar=c(4, 4, 2, 2))
        
        data <- kdata(k)
        data <- data[data[, FLAG] == RV, ]
        rvsamples <- getOption("systemic.rvsamples", 5000)
        trange <- k$trange
        if (is.nan(trange[1])) {
            trange <- c(1, 1000)
        }
        
        
        sl <- K_integrateRange(k$h, trange[1], trange[2], rvsamples, NULL, k$last.error.code)
        m <- .gsl_matrix_to_R(ok_get_rvs(sl, rvsamples))

        ylim <- c(min(data[,SVAL], m[,VAL]), max(data[, SVAL], m[,VAL]))

        if (! is.na(wrap)) {
            if (wrap == T) wrap <- k[1, 'period']
            data[, TIME] <- data[, TIME] %% wrap
            m[, TIME] <- m[, TIME] %% wrap
            m <- m[order(m[, TIME]), ]
        }

        suppressWarnings(plotCI(data[,TIME], data[, SVAL], data[, ERR], xlab="Time [JD]", ylab="Radial velocity [m/s]", ylim=ylim, col=data[,SET]+2, pch=20, gap=0, ...))
        lines(m[,TIME], m[,VAL])
        axis(3, labels=FALSE)
        axis(4, labels=FALSE)
        
        if (plot.residuals) {
            suppressWarnings(plotCI(data[,TIME], data[,SVAL] - data[, PRED], data[,ERR], xlab="Time [JD]", ylab="Residuals [m/s]",  col=data[,SET]+2, pch=20, gap = 0, ...))
            axis(3, labels=FALSE)
            axis(4, labels=FALSE)
        }
        
    } else if (type == "allrv") {
        stopifnot(k$nplanets > 0)
        np <- k$nplanets
        rows <- if (plot.residuals) np+1 else np
        par(mfrow=c(rows, 1), mar=c(4, 4, 2, 2))
        
        k <- kclone(k)
        rvsamples <- getOption("systemic.rvsamples", 5000)
        
        kcalculate(k)
        data <- kdata(k)
        data <- data[data[, FLAG] == RV, ]		
        trange <- c(min(data[,TIME]), max(data[,TIME]))
        
        ret <- list()
        
        for (i in 1:np) {
            masses <- k[, 'mass']
            k[i, 'mass'] <- 0

            kcalculate(k)
            data_i <- kdata(k)
            data_i <- data_i[data_i[, FLAG] == RV, ]			
            k[,'mass'] <- masses
            k[-i, 'mass'] <- 0
						
            sl <- K_integrateRange(k$h, trange[1], trange[2], rvsamples, NULL, k$last.error.code)
            m <- .gsl_matrix_to_R(ok_get_rvs(sl, rvsamples))
            
            ylim <- c(min(data_i[, SVAL] - data_i[, PRED], m[,VAL]), max(data_i[, SVAL] - data_i[, PRED], m[,VAL]))
            
            if (! is.na(wrap)) {
                data_i[, TIME] <- data_i[, TIME] %% k[i, 'period']
                m[, TIME] <- m[, TIME] %% k[i, 'period']
                m <- m[order(m[, TIME]), ]
            }
            xlim <- c(min(data_i[, TIME]), max(data_i[, TIME]))
            plotCI(data_i[,TIME], data_i[, SVAL] - data_i[,PRED], data_i[, ERR], xlab="Time [JD]", ylab=sprintf("RV, Planet %d [m/s]", i), ylim=ylim, col=data_i[,SET]+2, xlim=xlim, pch=20, gap = 0)
            lines(m[,TIME], m[,VAL], xlim=xlim)
            
            axis(3, labels=FALSE)
            axis(4, labels=FALSE)
            k[,'mass'] <- masses
        }
        
        if (plot.residuals) {
            plotCI(data[,TIME], data[,SVAL] - data[, PRED], data[,ERR], xlab="Time [JD]", ylab="Residuals [m/s]",  col=data[,SET]+2, gap = 0, pch=20)
            
            axis(3, labels=FALSE)
            axis(4, labels=FALSE)
        }
        

    } else if (type == "periodogram") {
        par(mfrow=c(2, 1), mar=c(4, 4, 1, 1))
        
        p <- kperiodogram(k, samples=getOption("systemic.psamples", 3e4), pmin=getOption("systemic.pmin", 0.5), pmax=getOption("systemic.pmax", 2e4))
        
        pr <- kperiodogram(k, per_type="res", samples=getOption("systemic.psamples", 3e4), pmin=getOption("systemic.pmin", 0.5), pmax=getOption("systemic.pmax", 2e4))
        
        plot(p, overplot.window = TRUE)
        plot(pr, overplot.window = TRUE)

    } else if (type == "residuals") {
        oldpar <- par()
        par(mfrow=c(ceiling(k$nsets/2), 2), mar=c(4, 4, 1, 1))
        data <- kdata(k)
        kcalculate(k)
        for (i in 0:k$nsets-1) {
            data_i <- data[data[,SET] == i, ]
            
            if (! is.na(wrap)) {
                data[, TIME] <- data[, TIME] %% wrap
            }
            
            if (nrow(data_i) > 0) {
                plot(data_i[,TIME], data_i[,SVAL]-data_i[,PRED], xlab="Time [JD]", ylab=sprintf("Radial velocity, dataset %d [m/s]", i), col=i+2)
                
                axis(3, labels=FALSE)
                axis(4, labels=FALSE)
            }
        }

    } else if (type == "ttv") {
        stopifnot(k$nplanets > 0)
        k <- kclone(k)
        data <- kdata(k)
        data <- data[data[, FLAG] == TIMING, ]
        data <- data[data[, TDS_FLAG] != TDS_SECONDARY, ]
        
        if (nrow(data) < 1)
            stop("Not enough timing data to plot")
        
        if (is.na(transiting.planet)) 
            transiting.planet = unique(data[, TDS_PLANET])

        par(mfrow=c(length(transiting.planet), 1), mar=c(4, 4, 1, 1))
        for (pl in transiting.planet) {
            
            kd <- data[data[, TDS_PLANET] == pl, ]			
            if (! is.na(transiting.per))
                per <- transiting.per[pl]
            else {
                per <- k[pl, 'period']
                perm <- min(diff(kd[, 1]))
                idx <- diff(kd[, 1]) < 1.4 * perm
                per <- median(diff(kd[, 1])[idx])
            }

            trange <- c(min(kd[, TIME]), max(kd[, TIME]))
            
            tsamp <- seq(from=trange[1], to=trange[2], by=per)
            ret <- kintegrate(k, tsamp, transits = TRUE)
            idx <- floor(kd[, TIME] / per)
            idx2 <- floor(tsamp / per)
            
            tim <- kd[, TIME] 
            fit <- lm(tim ~ idx, weights=1./kd[, ERR])
            omc <- fit$residuals
            
            if (max(abs(omc)) > 0.2 * per) {
                warning(sprintf("Potentially misdetected period for linear ephemeris fit [P = %e d] (or huge TTVs?)\nSpecify a better initial starting guess by supplying the transiting.per argument to plot", fit$coefficients[2]))
            }
            
            suppressWarnings(plotCI(kd[, TIME], omc, kd[, ERR], xlab="Time [d]", ylab="O-C [d]", pch=20, gap=0))
            
            lines(tsamp, (-fit$coefficients[1] - fit$coefficients[2] * idx2) + ret$transits[[pl]], col="red")
            
            axis(3, labels=FALSE)
            axis(4, labels=FALSE)
        }
        

    } else if (type == "normres") {
        kd <- kdata(k)
        
        err <- sqrt(kd[, ERR]^2 + (kpars(k)[kd[, SET]+1+10])^2)
        nr <- (kd[, PRED]-kd[, SVAL])/err
        
        if (is.na(breaks)) {
            h2 <- hist(nr, plot=FALSE)
            breaks <- h2$breaks
        }
        xlim <- if (is.na(xlim)) c(min(nr), max(nr)) else xlim
        ylim <- c(0, -1e10)
        
        h <- list()
        cat("# Median, mad and results of K-S test compared to a unit gaussian\n")
        m <- matrix(nrow=k$nsets, ncol=4)
        for (set in 1:k$nsets) {
            ns <- nr[kd[, SET]==(set-1)]
            m[set, 1] <- median(ns)
            m[set, 2] <- mad(ns)
            ks <- ks.test(ns, "pnorm")
            m[set, 3] <- ks$p.value
            m[set, 4] <- ks$statistic
            hs <- if (density) density(ns) else hist(ns, breaks = breaks, plot=FALSE)
            ylim <- if (density) c(0, max(ylim[2], max(hs$y))) else c(0, max(ylim[2], max(hs$density)))
            h[[set]] <- hs
        }
        ylim[2] <- 1.1*ylim[2]
        add <- FALSE
        b <- 1
        for (hs in h) {
            if (! density) {
                if (! add)
                    plot(c(hs$breaks, hs$breaks[length(hs$breaks)]), c(0, hs$density, 0), col=b, xlim=xlim, ylim=ylim, type="S", lwd=2, xlab="Normalized residuals", ylab="Density")
                else
                    lines(c(hs$breaks, hs$breaks[length(hs$breaks)]), c(0, hs$density, 0), col=b, xlim=xlim, ylim=ylim, type="S", lwd=2)
            } else {
                if (! add)
                    plot(hs$x, hs$y, xlim=xlim, ylim=ylim, lwd=2, col=b, xlab="Normalized residuals", ylab="Density", type="l", ...)
                else
                    lines(hs$x, hs$y, xlim=xlim, ylim=ylim, lwd=2, col=b, xlab="Normalized residuals", ylab="Density", ...)
            }
            add <- TRUE
            b <- b+1
        }
        x <- seq(min(xlim), max(xlim), length.out=40)
        lines(x, dnorm(x), lty="dotted")
        
        legend("topright", col=1:k$nsets, legend=sapply(k$datanames, basename), lty=rep(1, k$nsets), box.col="white")
        colnames(m) <- c("median", "mad", "p.value", "statistic")
        rownames(m) <- sapply(k$datanames, basename)
        print(m)
        cat("# Specify the breaks parameter to use a different number of bins\n")
    } else if (type == "orbits") {
        t <- seq(0, 2*pi, length.out=1000)

        xlim <- c(-max(k[, 'a'] * (1+k[, 'ecc'])), max(k[, 'a'] * (1+k[, 'ecc'])))

        plot(0, 0, pch=19, xlim=xlim, ylim=xlim, xlab="", ylab="")

        for (i in 1:k$nplanets) {
            f <- k[i, 'trueanomaly'] * pi/180
            pom <- k[i, 'lop'] * pi/180
            e <- k[i, 'ecc']
            a <- k[i, 'a']
            r <- a*(1-e^2)/(1+e*cos(t-pom))
            q <- a*(1-e)
            rf <- a*(1-e^2)/(1+e*cos(f))
            lines(r*cos(t), r*sin(t), lwd=2, type='l')
            lines(c(0, q*cos(pom)), c(0, q*sin(pom)), lty="dotted")
            points(rf*cos(f+pom), rf*sin(f+pom), pch=19)
        }
    }
    
    axis(3, labels=FALSE)
   
    par(mfrow=c(1, 1))
}

plot.periodogram <- function(p, overplot.window = F, what = 'power') {
    oldpar <- par(no.readonly=TRUE)
    on.exit(suppressWarnings(par(oldpar)))
    
    par(systemic.par)
    ymax = max(p[, what])
    if (overplot.window)
        ymax = max(ymax, p[, 'window'])
    
    plot(p[, 'period'], p[, what], type="l", log="x", xlab="Period [d]", ylab="Power", ylim=c(0, ymax))
    
    if (overplot.window) {
        lines(p[, 'period'], p[, 'window'], col="red")
    }
    
    if (sum(p[, 'fap'] < 0.1) > 3) {
        
        faps <- approx(log(p[,'fap']), p[,what], log(c(1e-1, 1e-2, 1e-3)))
        
        xmin = min(p[,'period'])
        xmax = max(p[,'period'])
        for (i in 1:3)	{
            if (! is.na(faps$y[i])) {
                lines(c(xmin, xmax), c(faps$y[i], faps$y[i]), lty=i+1)
            }
        }
    }	

    axis(3, labels=FALSE)
    axis(4, labels=FALSE)
    
    invisible()
}



.cap <- function(str) {
    return(paste(toupper(substring(str, 1, 1)), substring(str, 2), sep=""))
}

.grayscale.pal <- c(rgb(0, 0, 0), rgb(0.2, 0.2, 0.2), rgb(0.4, 0.4, 0.4), rgb(0.6, 0.6, 0.6), rgb(0.8, 0.8, 0.8), rgb(1, 1, 1))

                                        # Replacement for ci2d
.ci2d <- function (x, y = NULL, nbins = 400, method = c("bkde2D", "hist2d"), bandwidth, factor = 1, ci.levels = c(0.5, 0.75, 0.9, 0.95, 0.975), show = c("filled.contour", "contour", "image", "none"), col = topo.colors(length(breaks) - 1), show.points = FALSE, pch = par("pch"), points.col = "red", xlab, ylab, xlim, ylim, extra.points=NULL, ...) {

    oldpar <- par(no.readonly=TRUE)
    on.exit(suppressWarnings(par(oldpar)))

    par(systemic.par)
    show <- match.arg(show)
    method <- match.arg(method)
    breaks <- unique(c(0, ci.levels, 1))
    if (missing(xlab))  {
        xlab <- if (missing(x)) 
            ""
        else deparse(substitute(x))
    }
    if (missing(ylab)) {
        ylab <- if (missing(y)) 
            ""
        else deparse(substitute(y))
    }
    if (!is.null(y)) 
        x <- cbind(x, y)
    if (method == "hist2d") {
        h2d <- hist2d(x, show = FALSE, nbins = nbins, ...)
        h2d$density <- h2d$counts/sum(h2d$counts, na.rm = TRUE)
    }
    else if (method == "bkde2D") {
        if (length(nbins) == 1) 
            nbins <- c(nbins, nbins)
        if (missing(bandwidth)) {
            h.x = dpik(x[, 1])
            h.y = dpik(x[, 2])
            bandwidth <- c(h.x, h.y)
        }
        est <- bkde2D(x, bandwidth = bandwidth * factor, gridsize = nbins, 
                      ...)
        h2d <- list()
        h2d$x <- est$x1
        h2d$y <- est$x2
        h2d$counts <- est$fhat
        h2d$nobs <- nrow(x)
        h2d$density <- est$fhat/sum(est$fhat)
    }
    else stop("Unknown method: '", method, "'")
    uniqueVals <- rev(unique(sort(h2d$density)))
    cumProbs <- sapply(uniqueVals, function(val) sum(h2d$density[h2d$density >= 
                                                                 val]))
    names(cumProbs) <- uniqueVals
    h2d$cumDensity <- matrix(nrow = nrow(h2d$density), ncol = ncol(h2d$density))
    h2d$cumDensity[] <- cumProbs[as.character(h2d$density)]
    if (show == "image") {
        image(h2d$x, h2d$y, h2d$cumDensity, xlab = xlab, ylab = ylab, 
              breaks = breaks, col = col, xlim=xlim, ylim=ylim)
        if (show.points) 
            points(x[, 1], x[, 2], pch = '.', col = points.col)
    }
    else if (show == "filled.contour") {
        if (show.points) 
            plot.title <- function() {
                points(x[, 1], x[, 2], pch = pch, col = points.col)
                if (! is.null(extra.points))
                    points(extra.points[, 1], extra.points[, 2], pch=19, col='red')
            }
        else plot.title <- function() {
            if (! is.null(extra.points))
                points(extra.points[, 1], extra.points[, 2], pch=19, col='red')
        }
        filled.contour(h2d$x, h2d$y, h2d$cumDensity, levels = breaks, 
                       col = col, xlab = xlab, ylab = ylab, plot.title = plot.title(), 
                       key.title = title(""), key.axes = axis(4, 
                                                  at = breaks), xlim=xlim, ylim=ylim, cex.axis=0.5)
    }
    else if (show == "contour") {
        tmpBreaks <- breaks[breaks < 1]
        contour(h2d$x, h2d$y, h2d$cumDensity, levels = tmpBreaks, 
                labels = tmpBreaks, xlab = xlab, ylab = ylab, nlevels = length(tmpBreaks), 
                col = col, xlim=xlim, ylim=ylim)
        
        if (show.points) 
            points(x[, 1], x[, 2], pch = pch, col = points.col)
        if (! is.null(extra.points))
            points(extra.points[, 1], extra.points[, 2], pch=19, col='red')
    }
    h2d$contours <- contourLines(h2d$x, h2d$y, h2d$cumDensity, 
                                 levels = breaks, nlevels = length(breaks))
    names(h2d$contours) <- sapply(h2d$contours, function(x) x$level)
    h2d$contours <- lapply(h2d$contours, function(J) data.frame(x = J$x, 
                                                                y = J$y))
    h2d$call <- match.call()
    class(h2d) <- "ci2d"
    invisible(h2d)
}


plot.error.est <- function(e, type="histogram", px=list(1, "period"), py=list(1, "mass"), dev.factor = 5, planet, xlab, ylab, main, xlim, ylim, pch=16, col=ifelse(type=="histogram", 'white', 'black'), show.points = FALSE, points.col=20, breaks=c(0.5, 0.75, 0.9, 0.95, 0.99), cut.outliers = 12, scatter.bins = 8, ...) {
    oldpar <- par(no.readonly=TRUE)
    on.exit(suppressWarnings(par(oldpar)))
    x <- px
    y <- py
    bfx <- 0
    bfy <- 0
    datax <- NULL
    datay <- NULL
    medx <- 0
    medy <- 0
    devx <- 0
    devy <- 0
    labx <- ""
    laby <- ""
    par(systemic.par)
    if (x[[1]] != "par") {
        bfx <- e$fit.els[x[[1]], x[[2]]]
        datax <- e[[x[[1]]]][, x[[2]]]
        medx <- e$stats[[x[[1]]]][x[[2]], 'median']
        devx <- e$stats[[x[[1]]]][x[[2]], 'mad']
        labx <- sprintf("%s of planet %d", .elements.labels[[x[[2]]]], x[[1]])
    }
    else {
        bfx <- e$fit.params[x[[2]]]
        datax <- e$params[, x[[2]]]
        medx <- e$params.stats[x[[2]], 'median']
        devx <- e$params.stats[x[[2]], 'mad']
        labx <- sprintf("Parameter %s", as.character(x[[2]]))
    }

    if (y[[1]] != "par") {
        bfy <- e$fit.els[y[[1]], y[[2]]]
        datay <- e[[y[[1]]]][, y[[2]]]
        medy <- e$stats[[y[[1]]]][y[[2]], 'median']
        devy <- e$stats[[y[[1]]]][y[[2]], 'mad']		
        laby <- sprintf("%s of planet %d", .elements.labels[[y[[2]]]], y[[1]])
    }
    else {
        bfy <- e$fit.params[y[[2]]]
        datay <- e$params[, y[[2]]]
        medy <- e$params.stats[y[[2]], 'median']
        devy <- e$params.stats[y[[2]], 'mad']
        laby <- sprintf("Parameter %s", as.character(y[[2]]))		
    }


    if (!missing(xlab)) {
        labx <- xlab
    }
    if (!missing(ylab)) {
        laby <- ylab
    }

    if (dev.factor > 0 && type=="histogram") {
        datax <- datax[datax > medx - dev.factor * devx & datax < medx + dev.factor * devx]
    } else {
        limx <- c(min(datax),max(datax))
        limy <- c(min(datay),max(datay))
        
        if (dev.factor > 0) {
            limx <- c(max(limx[1], medx - dev.factor*devx), min(limx[2], medx + dev.factor * devx))
            limy <- c(max(limy[1], medy - dev.factor*devy), min(limy[2], medy + dev.factor * devy))
        }
    }

    if (!missing(xlim))
        limx <- xlim
    if (!missing(ylim))
        limy <- ylim

    if (type != "histogram" && type != "scatter" && type != "smoothScatter") {

        if (is.numeric(cut.outliers)) {
            idx <- datax > medx - cut.outliers * devx & datax < medx + cut.outliers * devx & datay > medy - cut.outliers*devy & datay < medy + cut.outliers*devy
            datax1 <- datax[idx]
            datay1 <- datay[idx]
            
            datax <- datax1
            datay <- datay1
        }
        bins <- c(as.integer(diff(range(datax))/devx * scatter.bins), as.integer(diff(range(datay))/devy * scatter.bins))
    }


    if (type == "histogram") {
        a <- hist(datax, freq=FALSE, xlab=labx, main=ifelse(missing(main), labx, main), col=col, ...)
        points(c(bfx), c(0.), pch=19,  col='red')
        return(invisible())
    } else if (type == "scatter") {
        plot(datax, datay, xlab=labx, ylab=laby, xlim=limx, ylim=limy, pch=pch, col=col,  ...)
        points(c(bfx), c(bfy), pch=19, col='red')
        
        axis(3, labels=FALSE)
        axis(4, labels=FALSE)
    } else if (type == "smoothScatter") {
        smoothScatter(datax, datay, nrpoints=0, xlab=labx, ylab=laby, col=col, xlim=limx, ylim=limy, pch=pch, ...)
        points(c(bfx), c(bfy), pch=19, col='red')
        
        axis(3, labels=FALSE)
        axis(4, labels=FALSE)
    } else if (type == "filled.contour") {
        
        fcol <- if (missing(col)) .grayscale.pal else col

        .ci2d(datax, datay, show="filled.contour", xlab=labx, ylab=laby, col=fcol, show.points=show.points, pch=".", extra.points=cbind(bfx, bfy), xlim=limx, ylim=limy, nbins=bins, ci.levels=breaks)
        
    } else if (type == "contour") {
        fcol <- if (missing(col)) "black" else col
        
        .ci2d(datax, datay, show="contour", xlab=labx, ylab=laby, col=fcol, show.points=show.points, extra.points=cbind(bfx, bfy), xlim=limx, ylim=limy, nbins=bins,  ci.levels=breaks, pch=pch)
        
    } else if (type == "image") {
        fcol <- if (missing(col)) "black" else col
        
        .ci2d(datax, datay, show="image", xlab=labx, ylab=laby, col=fcol, show.points=show.points, extra.points=cbind(bfx, bfy),  xlim=limx, ylim=ylim, nbins=bins,  ci.levels=breaks)
        
    } else if (type == "all.histograms") {
        stop("Not yet implemented.")
    }
    title(main=ifelse(missing(main), paste(labx, "vs", laby), main))

    invisible()
}

plot.integration <- function(int, what=c('a', 'ecc'), legend=TRUE) {
    oldpar <- par(no.readonly=TRUE)
    on.exit(suppressWarnings(par(oldpar)))

    par(mfrow=c(length(what), 1), mar=c(4, 4, 1, 8))
    par(systemic.par)
    all <- do.call("rbind", int$els)
    times <- (int$times-int$times[1])/(YEAR/DAY)
    for (el in what) {
        ymin <- min(all[, el], na.rm=TRUE)
        ymax <- max(all[, el], na.rm=TRUE)
        
        xlab <- sprintf("Time - Epoch (years)")
        ylab <- .elements.labels[[el]]
        
        for (i in 1:int$nplanets) {
            if (i == 1)
                plot(times, int$els[[i]][, el], col=i, 
                     xlim=c(0, max(times)), ylim=c(ymin, ymax), xlab=xlab, ylab=ylab, type="l")
            else
                lines(times, int$els[[i]][, el], col=i)
        }
        
        legend("topright", inset=c(-0.35,0), col=1:int$nplanets, legend=sprintf("Planet %d", 1:int$nplanets), xpd=TRUE, lty=rep(1, int$nplanets), box.col="white")
    }

}
