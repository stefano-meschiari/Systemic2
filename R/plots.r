systemic.par <- list()
systemic.par$font.lab <- 2
systemic.par$tck <- 0.02

systemic.palette <- systemic.theme.tomorrow
systemic.palette.face <- systemic.theme.tomorrow.face

palette(systemic.palette)
par(systemic.par)

plot.kernel <- function(k, type = "rv", wrap=NA, plot.residuals=TRUE, transiting.planet = NA, transiting.per = NA, xlim = NULL, ylim=NULL, which.planets=1:k$nplanets,
                        breaks=NA, plot.gaussian=TRUE, density=FALSE, pch=21, lwd=2, layout=TRUE, separate.sets=TRUE, xlab=NULL, ylab=NULL, col=0, yshift=0, ...) {
    .check_kernel(k)
    par(systemic.par)
    if (is.nan(k$epoch)) {
        stop("No epoch set")
    }
    
    if (type == "rv") {
        rows <- if (plot.residuals) 2 else 1
        if (layout) {
            par(mfrow=c(1,1))
            par(mfrow=c(rows, 1), mar=c(4.1, 5.1, 2.1, 2.1))
        }
        data <- kdata(k)
        data <- data[data[, FLAG] == RV, ]
        rvsamples <- getOption("systemic.rvsamples", 5000)
        trange <- k$trange
        if (is.nan(trange[1])) {
            trange <- c(1, 1000)
        }

        xlab <- if (!is.null(xlab)) xlab else 'Time [JD]'
        ylab <- if (!is.null(ylab)) ylab else 'Radial velocity [m/s]'
        
        sl <- K_integrateRange(k$h, trange[1], trange[2], rvsamples, NULL, k$last.error.code)
        m <- .gsl_matrix_to_R(ok_get_rvs(sl, rvsamples))
        ok_free_systems(sl, rvsamples)

        ylim <- if (is.null(ylim)) c(min(data[,SVAL], m[,VAL]), max(data[, SVAL], m[,VAL]))
        ylim <- ylim + yshift
        if (! is.na(wrap)) {
            if (wrap == T) wrap <- k[1, 'period']
            data[, TIME] <- data[, TIME] %% wrap
            m[, TIME] <- m[, TIME] %% wrap
            m <- m[order(m[, TIME]), ]
        }

        suppressWarnings(plotCI(data[,TIME], data[, SVAL]+yshift, data[, ERR], sfrac=0, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=data[,SET]+2+col, pch=pch, gap=0, pt.bg=systemic.palette.face[data[,SET]+2+col], ...))
        lines(m[,TIME], m[,VAL] + yshift, lwd=lwd)
        axis(3, labels=FALSE)
        axis(4, labels=FALSE)
        
        if (plot.residuals) {
            suppressWarnings(plotCI(data[,TIME], data[,SVAL] - data[, PRED] + yshift, data[,ERR], xlab="Time [JD]", ylab="Residuals [m/s]", ylim=ylim, xlim=xlim, col=data[,SET]+2+col, pch=pch, sfrac=0, gap = 0, pt.bg=systemic.palette.face[data[,SET]+2+col], ...))
            axis(3, labels=FALSE)
            axis(4, labels=FALSE)
        }
        
    } else if (type == "allrv") {
        stopifnot(k$nplanets > 0)
        np <- k$nplanets
        rows <- if (plot.residuals) length(which.planets)+1 else length(which.planets)
        if (layout)
            par(mfrow=c(rows, 1))

        fpars <- list(...)
        xlab <- if (!is.null(fpars$xlab)) fpars$xlab else 'Time [JD]'
        ylab <- if (!is.null(fpars$ylab)) fpars$ylab else 'RV, Planet %d [m/s]'
        
        
        k <- kclone(k)
        rvsamples <- getOption("systemic.rvsamples", 5000)
        
        kcalculate(k)
        data <- kdata(k)
        data <- data[data[, FLAG] == RV, ]		
        trange <- c(min(data[,TIME]), max(data[,TIME]))
        
        ret <- list()
        for (i in 1:k$nplanets)
            krange(k, i, 'mass') <- c(0, NaN)
        
        for (i in which.planets) {
            masses <- k[, 'mass']

            
            k[i, 'mass'] <- 0

            kcalculate(k)
            data_i <- kdata(k)
            data_i <- data_i[data_i[, FLAG] == RV, ]			
            k[,'mass'] <- masses
            k[-i, 'mass'] <- 0
						print(k)
            print(i)
            sl <- K_integrateRange(k$h, trange[1], trange[2], rvsamples, NULL, k$last.error.code)
            m <- .gsl_matrix_to_R(ok_get_rvs(sl, rvsamples))
            ok_free_systems(sl, rvsamples)
            
            ylim <- if (is.null(ylim)) c(min(data_i[, SVAL] - data_i[, PRED], m[,VAL]), max(data_i[, SVAL] - data_i[, PRED], m[,VAL]))
            
            if (! is.na(wrap)) {
                data_i[, TIME] <- data_i[, TIME] %% k[i, 'period']
                m[, TIME] <- m[, TIME] %% k[i, 'period']
                m <- m[order(m[, TIME]), ]
            }
            xlim <- c(min(data_i[, TIME]), max(data_i[, TIME]))
            plotCI(data_i[,TIME], data_i[, SVAL] - data_i[,PRED], data_i[, ERR], sfrac=0, xlab=xlab, ylab=sprintf(ylab, i), ylim=ylim, col=data_i[,SET]+2, xlim=xlim, pch=pch, gap = 0, pt.bg=systemic.palette.face[data[,SET]+2])
            lines(m[,TIME], m[,VAL], xlim=xlim, lwd=lwd)
            
            axis(3, labels=FALSE)
            axis(4, labels=FALSE)
            k[,'mass'] <- masses
        }
        
        if (plot.residuals) {
            if (!is.na(wrap)) {
                data[,TIME] <- data[,TIME] %% k[which.planets[length(which.planets)], 'period']
            }
            plotCI(data[,TIME], data[,SVAL] - data[, PRED], data[,ERR], xlab="Time [JD]", ylab="Residuals [m/s]", sfrac=0, col=data[,SET]+2, gap = 0, pch=pch, pt.bg=systemic.palette.face[data[,SET]+2])
            
            axis(3, labels=FALSE)
            axis(4, labels=FALSE)
        }
        

    } else if (type == "periodogram") {
        if (layout)
            par(mfrow=c(2, 1), mar=c(4, 4, 1, 1))
        
        p <- kperiodogram(k, samples=getOption("systemic.psamples", 3e4), pmin=getOption("systemic.pmin", 0.5), pmax=getOption("systemic.pmax", 2e4))
        
        pr <- kperiodogram(k, per_type="res", samples=getOption("systemic.psamples", 3e4), pmin=getOption("systemic.pmin", 0.5), pmax=getOption("systemic.pmax", 2e4))
        
        plot(p, overplot.window = TRUE, lwd=lwd)
        plot(pr, overplot.window = TRUE, lwd=lwd)

    } else if (type == "residuals") {
        if (layout)
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
        if (layout)
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
            
            suppressWarnings(plotCI(kd[, TIME], omc, kd[, ERR], xlab="Time [d]", ylab="O-C [d]", pch=pch, sfrac=0, gap=0))
            
            lines(tsamp, (-fit$coefficients[1] - fit$coefficients[2] * idx2) + ret$transits[[pl]], col="red", lwd=lwd)
            
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
        xlim <- if (is.null(xlim)) c(min(nr), max(nr)) else xlim
        ylim <- c(0, -1e10)
        
        h <- list()
        cat("# Median, mad and results of K-S test compared to a unit gaussian\n")
        m <- matrix(nrow=k$nsets, ncol=4)
        nsets <- k$nsets
        if (!separate.sets)
            nsets <- 1
        
        for (set in 1:nsets) {
            if (separate.sets)
                ns <- nr[kd[, SET]==(set-1)]
            else
                ns <- nr
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
                    plot(c(hs$breaks, hs$breaks[length(hs$breaks)]), c(0, hs$density, 0), col=b+1, xlim=xlim, ylim=ylim, type="S", lwd=lwd, xlab="Normalized residuals", ylab="Density")
                else
                    lines(c(hs$breaks, hs$breaks[length(hs$breaks)]), c(0, hs$density, 0), col=b+1, xlim=xlim, ylim=ylim, type="S", lwd=lwd)
            } else {
                if (! add)
                    plot(hs$x, hs$y, xlim=xlim, ylim=ylim, lwd=lwd, col=b+1, xlab="Normalized residuals", ylab="Density", type="l", ...)
                else
                    lines(hs$x, hs$y, xlim=xlim, ylim=ylim, lwd=lwd, col=b+1, xlab="Normalized residuals", ylab="Density", ...)
            }
            add <- TRUE
            b <- b+1
        }
        x <- seq(min(xlim), max(xlim), length.out=40)
        lines(x, dnorm(x), lty="dotted")
        
        legend("topright", col=1:nsets, legend=sapply(k$datanames, basename), lty=rep(1, k$nsets), box.col="white")
        colnames(m) <- c("median", "mad", "p.value", "statistic")
        rownames(m) <- sapply(k$datanames, basename)
        print(m)
        cat("# Specify the breaks parameter to use a different number of bins\n")
    } else if (type == "orbits") {
        plot.orbit(k, ...)
    }
    
    axis(3, labels=FALSE)
}

plot.orbit <- function(k, planet=-1, nplanets=k$nplanets, samples=1000, samples.alpha=0.1, samples.lwd=1, xlim=NA, ylim=NA, add=FALSE, plot.pericenter=TRUE, plot.planet=TRUE, lwd=2, col='black', best.col='red', xlab="", ylab="", plot.scale=TRUE, axes=FALSE, ...) {
    oldpar <- par('pty')
    par(systemic.par)

    par(pty="s")
    on.exit(suppressWarnings(par(oldpar)))
    
    if (("kernel" %in% class(k)) || (".orbit" %in% class(k))) {
        if ("kernel" %in% class(k))
            if (k$element.type != ASTROCENTRIC)
                stop("Currently only plots orbits for ASTROCENTRIC orbital elements.")

        # Eliminates weird "pattern"
        r <- runif(1, max=2*pi)
        t <- seq(r, 2*pi+r, length.out=250)
        
        if (planet == -1)
            planet = 1:nplanets
        if (is.na(xlim) || is.na(ylim)) {
            xlim <- c(-max(k[, 'a'] * (1+k[, 'ecc'])), max(k[, 'a'] * (1+k[, 'ecc'])))
            ylim <- xlim
        }
        if (!add)
            plot(0, 0, pch=19, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, axes=axes)

        for (i in planet) {
            f <- k[i, 'trueanomaly'] * pi/180
            pom <- k[i, 'lop'] * pi/180
            e <- k[i, 'ecc']
            a <- k[i, 'a']
            r <- a*(1-e^2)/(1+e*cos(t-pom))
            q <- a*(1-e)
            rf <- a*(1-e^2)/(1+e*cos(f))
            lines(r*cos(t), r*sin(t), type='l', lwd=lwd, col=col)
            if (plot.pericenter && e > 0)
                lines(c(0, q*cos(pom)), c(0, q*sin(pom)), lty="dotted", col=col)
            if (plot.planet)
                points(rf*cos(f+pom), rf*sin(f+pom), pch=19, col=col)
        }

        if (plot.scale) {
            span <- par('usr')
            w <- diff(range(span))
            da <- 2^ceil(log2(0.1*w))
            x <- 0.9 * span[2]
            y <- 0.9 * span[3]
            lines(c(x-da, x),
                  c(y, y), col='black', lwd=2)
            text(x-0.5*da, y-0.05*w, adj=c(0.5, 0), labels=paste(da, 'AU'), cex=par('cex.lab'))
        }
        
    } else if ("error.est" %in% class(k)) {
        if (!is.null(k$element.type) && k$element.type != ASTROCENTRIC)
            stop("Currently only plots orbits for ASTROCENTRIC orbital elements.");
        if (k$size < samples)
            warning("There are fewer elements than samples.")
        
        samples = min(samples, k$size)
        samples.indices = sample(1:samples, samples)

        if (planet == -1)
            planet = 1:k$nplanets
        
        lim <- Reduce(max, sapply(planet, function(i) k[[i]][,'a'] * (1+k[[i]][, 'ecc'])))
        if (is.na(xlim) || is.na(ylim)) {
            xlim <- c(-lim, lim)
            ylim <- xlim
        }
        
        plot(0, 0, pch=19, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, axes=axes)

        obj <- matrix(nrow=length(planet), ncol=4)
        colnames(obj) <- c('a', 'ecc', 'lop', 'trueanomaly')
        class(obj) <- c('matrix', '.orbit')
        samples.col <- rgb(0, 0, 0, samples.alpha)
        
        for (j in samples.indices) {
            for (e in colnames(obj))
                obj[, e] <- sapply(planet, function(i) k[[i]][j, e])

            plot.orbit(obj, xlim=lim, ylim=lim, planet=planet, nplanets=k$nplanets,
                       lwd=samples.lwd, col=samples.col, xlab=xlab, ylab=ylab, add=TRUE, plot.pericenter=FALSE, plot.planet=FALSE, plot.scale=FALSE)
        }

        obj <- k$fit.els
        class(obj) <- c('matrix', '.orbit')
        plot.orbit(obj, xlim=lim, ylim=lim, planet=planet, nplanets=k$nplanets,
                   lwd=lwd, col=best.col, xlab=xlab, ylab=ylab, add=TRUE)
        
        
    } else {
        stop("First argument should be an object of class kernel, or of class error.est")
    }
}


plot.periodogram <- function(p, overplot.window = F, what = 'power', plot.fap = TRUE, xlim, ylim, xlab, ylab, show.resampled=FALSE, ...) {
    par(systemic.par)

    if (!is.null(attr(p, 'resampled')) && !(show.resampled)) {
        cat(sprintf("# Note: this periodogram samples a lot of frequencies [%d]\n", nrow(p)))
        cat(sprintf("# To avoid large files and slow PDFs, specify show.resampled=TRUE as a parameter to only plot important frequencies.\n"))
    } else if (show.resampled) {
        cat(sprintf("# Using resampled periodogram with %d frequencies instead of %d.\n",
                    nrow(attr(p, 'resampled')), nrow(p)))
        p <- attr(p, 'resampled')
    }
    
    ymax = max(p[, what])
    if (overplot.window)
        ymax = max(ymax, p[, 'window'])
    if (missing(ylim))
        ylim = c(0, ymax)
    if (missing(xlim))
        xlim = c(min(p[, 'period']), max(p[, 'period']))
    if (missing(xlab))
        xlab = 'Period [d]'
    if (missing(ylab))
        ylab = 'Normalized power'
    
    plot(p[, 'period'], p[, what], type="l", log="x", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
    
    if (overplot.window) {
        lines(p[, 'period'], p[, 'window'], col="red")
    }
    
    if (plot.fap && sum(p[, 'fap'] < 0.1) > 3) {
        
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
.ci2d <- function (x, y = NULL, nbins = 400, method = c("bkde2D", "hist2d"), 
                   bandwidth, factor = 1, ci.levels = c(0.5, 0.75, 0.9, 0.95, 0.975), show = c("filled.contour", "contour", "image", "none"), col = topo.colors(length(breaks) - 1), show.points = FALSE, pch = par("pch"), points.col = "black", xlab, ylab, xlim, ylim, extra.points=NULL, lwd=1, add=FALSE, ...) 
{
    show <- match.arg(show)
    method <- match.arg(method)
    breaks <- unique(c(0, ci.levels, 1))
    if (missing(xlab)) 
        xlab <- if (missing(x)) 
            ""
        else deparse(substitute(x))
    if (missing(ylab)) 
        ylab <- if (missing(y)) 
            ""
        else deparse(substitute(y))
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
                col = col, xlim=xlim, ylim=ylim, add=add, lwd=lwd)
        if (show.points) {
            points(x[, 1], x[, 2], pch = pch, col = points.col)
            contour(h2d$x, h2d$y, h2d$cumDensity, levels = tmpBreaks, 
                labels = tmpBreaks, xlab = xlab, ylab = ylab, nlevels = length(tmpBreaks), 
                col = col, xlim=xlim, ylim=ylim, add=TRUE, lwd=lwd)
            
        }
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


plot.error.est <- function(e, type="histogram", px=list(1, "period"), py=NULL, dev.factor = 5, planet, xlab, ylab, main, xlim, ylim, pch=16, col=ifelse(type=="histogram", 'white', 'black'), show.points = FALSE, points.col=20, breaks=c(0.5, 0.75, 0.9, 0.95, 0.99), cut.outliers = 12, scatter.bins = 8, add=FALSE, bf.color='red', subset=1:e$length, ...) {
    par(systemic.par)
    
    if (!is.null(px)) {
        if (is.character(px[[2]]) && px[1] != 'par')
            px[[2]] <- which(.allelements == px[[2]])
        else
            px[[2]] <- which(.params == px[[2]])
    }
    if (!is.null(py)) {
        if (is.character(py[[2]]) && py[1] != 'par')
            py[[2]] <- which(.allelements == py[[2]])
        else
            py[[2]] <- which(.params == py[[2]])
    }

    x <- px
    y <- py

    pars <- list(...)
    lwd <- if (is.null(pars$lwd)) 1 else pars$lwd
    
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

    if (! is.null(x)) {
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
    }
    if (! is.null(y)) {
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
    }

    if (!missing(xlab)) {
        labx <- xlab
    }
    if (!missing(ylab)) {
        laby <- ylab
    }

    datax <- datax[subset]
    datay <- datay[subset]
    
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
        if (length(scatter.bins) == 1)
            scatter.bins <- c(scatter.bins, scatter.bins)
        bins <- c(as.integer(diff(range(datax))/devx * scatter.bins), as.integer(diff(range(datay))/devy * scatter.bins))
    }
    
    
    if (type == "histogram") {
        a <- hist(datax, freq=FALSE, xlab=labx, main=ifelse(missing(main), labx, main), col=col, ...)
        points(c(bfx), c(0.), pch=19,  col=bf.color, ...)
        return(invisible())
    } else if (type == "scatter") {
        if (!add) {
            plot(datax, datay, xlab=labx, ylab=laby, xlim=limx, ylim=limy, pch=pch, col=col,  ...)
        }
        else 
            points(datax, datay, pch=pch, col=col,  ...)
        points(c(bfx), c(bfy), pch=19, col=bf.color, ...)
        
    } else if (type == "smoothScatter") {
        smoothScatter(datax, datay, nrpoints=0, xlab=labx, ylab=laby, col=col, xlim=limx, ylim=limy, pch=pch, ...)
        points(c(bfx), c(bfy), pch=19, col=bf.color)
        
        axis(3, labels=FALSE)
        axis(4, labels=FALSE)
    } else if (type == "filled.contour") {
        
        fcol <- if (missing(col)) .grayscale.pal else col

        .ci2d(datax, datay, show="filled.contour", xlab=labx, ylab=laby, col=fcol, show.points=show.points, pch=".", extra.points=cbind(bfx, bfy), xlim=limx, ylim=limy, nbins=bins, ci.levels=breaks)
        
    } else if (type == "contour") {
        fcol <- if (missing(col)) "black" else col
        
        .ci2d(datax, datay, show="contour", xlab=labx, ylab=laby, col=fcol, show.points=show.points, extra.points=cbind(bfx, bfy), xlim=limx, ylim=limy, nbins=bins,  ci.levels=breaks, pch=pch, points.col=points.col, add=add, lwd=lwd)
        
    } else if (type == "image") {
        fcol <- if (missing(col)) "black" else col
        
        .ci2d(datax, datay, show="image", xlab=labx, ylab=laby, col=fcol, show.points=show.points, extra.points=cbind(bfx, bfy),  xlim=limx, ylim=ylim, nbins=bins,  ci.levels=breaks)
        
    } else if (type == "all.histograms") {
        stop("Not yet implemented.")
    } else {
        stop(sprintf("Type '%s' not recognized.", type))
    }
    if (!add)
        title(main=ifelse(missing(main), paste(labx, "vs", laby), main))
    
    invisible()
}

plot.integration <- function(int, what=c('a', 'ecc'), legend=TRUE, xlab="Time - Epoch (years)", ylab='%s', ...) {
    par(mfrow=c(length(what), 1), mar=c(4, 4, 1, 8))
    par(systemic.par)

    all <- do.call("rbind", int$els)
    times <- (int$times-int$times[1])/(YEAR/DAY)
    for (el in what) {
        ymin <- min(all[, el], na.rm=TRUE)
        ymax <- max(all[, el], na.rm=TRUE)
        
        .ylab <- sprintf(ylab, .elements.labels[[el]])
        
        for (i in 1:int$nplanets) {
            if (i == 1)
                plot(times, int$els[[i]][, el], col=i+1, 
                     xlim=c(0, max(times)), ylim=c(ymin, ymax), xlab=xlab, ylab=.ylab, type="l", ...)
            else
                lines(times, int$els[[i]][, el], col=i+1, ...)
        }
        
        legend("topright", inset=c(-0.35,0), col=(1:int$nplanets)+1, legend=sprintf("Planet %d", 1:int$nplanets), xpd=TRUE, lty=rep(1, int$nplanets), box.col="white")
    }
    

}
