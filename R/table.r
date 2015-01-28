require('stringr')
require('xtable')



nformat <- function(n, err=0, digits=2, fmt="%s \\pm{} %s") {
    if (err != 0) {
        e10 <- floor(log10(abs(err)))
        if (e10 > 0) {
            err <- round(err)
            n <- round(n)
        } else {
            fmt2 <- sprintf("%%.%df", -e10)            
            err <- sprintf(fmt2, round(err, -e10))
            n <- sprintf(fmt2, round(n, -e10))
        }
        return(sprintf(fmt, n, err))
    } else {
        if (floor(n) == n)
            n <- sprintf("%d", n)
        else
            n <- sprintf(sprintf("%%.%df", digits), n)

        return(n)
    }
}

ktable <- function(k, what=c('period', 'mass', 'ma', 'ecc', 'lop', 'k', 'a', 'tperi', 'mstar', 'rms', 'jitter', 'epoch', 'ndata', 'trange', 'pars.minimized'), labels=systemic.names, units=systemic.units, star.names=NULL, default.format="%.2f", default.nf="%s \\pm{} %s") {
    systemic.names <- labels
    systemic.units <- units
    
    if (class(k) == 'kernel')
        k <- list(k)
    if (is.null(star.names))
        star.names <- rep('', length(k))
    if (length(star.names) != length(k))
        stop("There are more kernels than star names")
    
    df <- data.frame()
    
    if ('pars.minimized' %in% what) {
        idx <- which(what == 'pars.minimized')
        parnames <- lapply(k, function(k) {
            p <- kflags(k)$par
            p <- p[p == bitOr(ACTIVE, MINIMIZE)]
            p <- p[!is.na(systemic.names[names(p)])]
            return(names(p))
        })
        
        if (length(parnames[[1]]) > 0) {
            what <- c(what[-idx], unlist(parnames))
        } else
            what <- what[-idx]
    }
    
    planet.labels <- unlist(lapply(1:length(k), function(i) str_join(star.names[i], letters[1+1:k[[i]]$nplanets])))
    
    row.labels <- sapply(what, function(n) sprintf("%s %s", systemic.names[n], systemic.units[n]))
    
    df <- matrix('', nrow=length(row.labels), ncol=length(planet.labels))
    colnames(df) <- planet.labels
    rownames(df) <- row.labels

    col <- 1
    for (i in 1:length(k)) {
        kk <- k[[i]]
        if (is.null(kk$errors))
            cat(sprintf("[Warning: Errors were not calculated for kernel #%d]\n", i))

        for (pl in 1:kk$nplanets) {
                for (j in 1:length(what)) {
                    if (systemic.type[what[j]] == ELEMENT) {
                        if (!is.null(kk$errors)) 
                            df[j, col] <- nformat(kk$errors$stats[[pl]][what[j], 'median'],
                                              kk$errors$stats[[pl]][what[j], 'mad'], fmt=default.nf)
                        else
                            df[j, col] <- nformat(kk[pl, what[j]])
                        
                    }

                    if (pl > 1)
                        next
                    
                    if (systemic.type[what[j]] == PARAMETER) {
                        if (!is.null(kk$errors))
                            df[j, col] <- nformat(kk$errors$params.stats[what[j], 'median'],
                                              kk$errors$params.stats[what[j], 'mad'])
                        else
                            df[j, col] <- nformat(kk['par', what[j]])
                        
                    } else if (systemic.type[what[j]] == PROPERTY) {
                        
                        p <- kget(kk, what[j])
                        if (length(p) > 1)
                            p <- str_join(sapply(p, nformat), collapse=' - ')
                        else
                            p <- nformat(p)
                        df[j, col] <- p
                    }
                }
            col <- col+1
        }        
    }

    class(df) <- c('systemic.table', 'matrix')
    
    return(df)
}

.systemic.table.display <- function(df) {
    class(df) <- "matrix"
    df <- gsub("(\\\\pm\\{\\})", "+-", df)
    df <- gsub("(\\\\times)", "x", df)
    rownames(df) <- gsub("[\\{\\}]", "", rownames(df))
    
    return(df)
}

print.systemic.table <- function(df, type="text", file=stdout()) {
    if (class(file) == "character") {
        file <- file(file)
        on.exit(close(file))
    }
        
    if (type == "text") {
        sink(file)
        cat("\n")
        print(noquote(.systemic.table.display(df)))
        cat("\n")
        sink()
    } else if (type == "latex") {

    }
}

plot.systemic.table <- function(df, ...) {
    textplot(.systemic.table.display(df), ...)
}


options(xtable.sanitize.text.function=function(x) x)
