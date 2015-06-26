require('stringr')
require('xtable')

dev.copy2png <- function(file, ...) {
  current.device <- dev.cur()
  on.exit(dev.set(current.device))
  dev.copy(png, file = file,  bg = "white", ...)
  dev.off()
}

systemic.units.latex <- c(period='[days]', mass='[$\\mass_\\s{jup}$]', ma='[deg]', ecc='',
                         lop='[deg]', inc='[deg]', node='[deg]',
                         a='[AU]', k='[$\\mathrm{m s}^{-1}$]', tperi='[JD]',
                         rv.trend='[$\\mathrm{m s}^{-1}$]', rv.trend.quadratic='[$\\mathrm{m s}^{-1}$^2]',
                         mstar = '[$\\mass_\\odot$]', chi2='', jitter='[$\\mathrm{m s}^{-1}$]', rms='[$\\mathrm{m s}^{-1}$]',
                         epoch = '[JD]', ndata='', trange='[JD]', loglik='', chi2nr='')
for (i in 1:10) {
  systemic.units.latex[sprintf('data.noise%d', i)] <- '[$\\mathrm{m s}^{-1}$]'
  systemic.units.latex[sprintf('data%d', i)] <- '[$\\mathrm{m s}^{-1}$]'
}

nformat <- function(n, err=0, digits=3, fmt="%s [%s]") {
  if (is.na(n))
    return('NaN')
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

kval <- function(k, what, idx='all', digits=3) {
  if (idx == 'all')
    return(sapply(1:k$nplanets, function(i) kval(k, what, i)))
  else if (idx == 'min')
    return(kval(k, what, which.min(k[,what])))
  else if (idx == 'max')
    return(kval(k, what, which.max(k[,what])))
  else {
    if (!is.null(k$errors)) {
      if (what != 'par')
        return(nformat(k$errors$stats[[idx]][what, 'median'],
                       k$errors$stats[[idx]][what, 'mad'],
                       format='\\ensuremath{%s \\pm %s}'))
    } else {
      return(nformat(k[idx, what], digits=digits))
    }
  }
}

.kis.active <- function(k, a, b) {
  return(bitAnd(kflag(k, a, b), ACTIVE) == ACTIVE)
}

.linesep <- '        '

.namedata <- function(k, s) {
  if (str_detect(s, 'dataset [[:digit:]]')) {
    n <- as.numeric(str_match(s, 'dataset ([[:digit:]])')[,2])
    return(str_replace(s, 'dataset ([[:digit:]])', str_c('dataset ', k$datanames[n])))
  } else
    return(s)
}

ktable <- function(k, what=c('period', 'mass', 'ma', 'ecc', 'lop', 'k', 'a', 'tperi', '-', sprintf('data.noise%d', 1:9), sprintf('data%d', 1:9), '-', 'mstar', 'chi2nr', 'chi2', 'loglik', 'rms', 'jitter', 'epoch', 'ndata', 'trange'), labels=systemic.names, units=if (!latex) systemic.units else systemic.units.latex, caption=NULL, default.format="%.2f", default.nf="%s [%s]", latex=FALSE, sep.length=15) {
  
  systemic.names <- labels  
  systemic.units <- units

  
  
  stopifnot(class(k) == 'kernel')    
  df <- data.frame()
  
  if ('pars.minimized' %in% what) {
    idx <- which(what == 'pars.minimized')
    
    p <- kflags(k)$par
    p <- p[p == bitOr(ACTIVE, MINIMIZE)]
    p <- p[!is.na(systemic.names[names(p)])]
    parnames <- names(p)
    
    if (length(parnames[[1]]) > 0) {
      what <- c(what[-idx], unlist(parnames))
    } else
      what <- what[-idx]
  }
  what <- Filter(function(what) {
    if (what == '-')
      return(TRUE)
    else if (systemic.type[what] == ELEMENT) {
      if (what %in% .elements) {
        return(any(sapply(1:k$nplanets, function(i) {
          return(.kis.active(k, i, what))
        })))
      } else {
        return(TRUE)
      }
      
    } else if (systemic.type[what] == PARAMETER)
      return(.kis.active(k, 'par', what))
    else
      return(TRUE)
  }, what)
  
  planet.labels <- str_c(k$starname, letters[1+1:k$nplanets])
  
  row.labels <- sapply(what, function(n) {
    if (n == '-')
      return('-')
    else
      sprintf("%s %s", .namedata(k, systemic.names[n]), systemic.units[n])
  })

  
  df <- matrix('', nrow=length(row.labels), ncol=length(planet.labels))
  row.labels[row.labels == '-'] <- sapply(seq_along(which(row.labels == '-')), function(i) str_c(.linesep, str_rep(" ", i)))

  colnames(df) <- planet.labels
  rownames(df) <- row.labels
  
  col <- 1

  
  if (is.null(k$errors))
    warning(sprintf("[Warning: Errors were not calculated for kernel #%d]\n", i))

  errors <- k$errors
  if (!is.null(errors)) {
    if (errors$nplanets != k$nplanets || errors$nsets != k$nsets) {
      warning("Errors not applicable for current kernel (different # of planets or different data)")
      errors <- NULL
    }
  }
  
  for (pl in 1:k$nplanets) {
    for (j in 1:length(what)) {

      if (what[j] == '-') {
        df[j, col] <- .linesep
        next
      }
      
      if (systemic.type[what[j]] == ELEMENT) {

        if (!is.null(errors)) {
          df[j, col] <- nformat(errors$stats[[pl]][what[j], 'median'],
                                errors$stats[[pl]][what[j], 'mad'], fmt=default.nf)
        } else
          df[j, col] <- nformat(k[pl, what[j]])
        
      }
      
      if (pl > 1)
        next
      
      if (systemic.type[what[j]] == PARAMETER) {
        if (!is.null(errors))
          df[j, col] <- nformat(errors$params.stats[what[j], 'median'],
                                errors$params.stats[what[j], 'mad'])
        else
          df[j, col] <- nformat(k['par', what[j]])
        
      } else if (systemic.type[what[j]] == PROPERTY) {
        
        p <- kget(k, what[j])
        if (what[j] == 'trange') {
          p <- paste(nformat(p[1]), nformat(p[2]-p[1]), sep='+')
        } else {
          if (length(p) > 1)
            p <- str_c(sapply(p, nformat), collapse=' - ')
          else
            p <- nformat(p)
        }
        df[j, col] <- p
      }
    }
    col <- col+1
  }        

  class(df) <- c('systemic.table', 'matrix')
  attr(df, 'latex') <- latex
  attr(df, 'caption') <- caption
  
  return(df)
}

.systemic.table.display <- function(df) {
  class(df) <- "matrix"
  df <- gsub("(\\\\pm\\{\\})", "+-", df)
  df <- gsub("(\\\\times)", "x", df)
  rownames(df) <- gsub("[\\{\\}]", "", rownames(df))
  attr(df, 'latex') <- NULL
  attr(df, 'caption') <- NULL
  
  return(df)
}

print.systemic.table <- function(df, file=stdout(), type=NA, font.size="normalsize") {
  if (any(class(file) == "character")) {
    file <- file(file) 
    on.exit(close(file))
  }
  
  if (is.na(type)) {
    if (!attr(df, 'latex'))
      type = 'text'
    else
      type = 'latex'
  }
  if (type == "text") {
    sink(file)
    cat("\n")
    print(noquote(.systemic.table.display(df)))
    cat("\n")
    if (file != stdout())
      sink()
    return(invisible())
  } else if (type == "latex") {
    s <- paste(toLatex(xtable(df, caption=attr(df, 'caption'))), collapse='\n')
    s <- str_replace_all(s, str_c(.linesep, '.+'), "\\\\hline")
    s <- str_replace_all(s, "\\\\centering", paste0("\\\\centering\\\\", font.size))
    cat(s, file=file)
    return(invisible(s))
  }
}

plot.systemic.table <- function(df, ...) {
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  textplot(.systemic.table.display(df), ...)
}

.append.figure <- function(file, width, height) {
  dev.copy2pdf(file, width, height)
}

latex.report <- function(k, where=stop("Specify a folder where to save the file."), trials=1e4, integration=1e4, report.file=str_c(getOption('systemic.path'), '/report.tex')) {
  default.par <- par(no.readonly=TRUE)
  cur <- getwd()
  on.exit({ setwd(cur); par(default.par) })
  setwd(where)

  file.copy(report.file, getwd(), overwrite=FALSE)
  print(ktable(k, latex=TRUE, caption='Best fit'), file='table.tex', font.size='small')

  
  if (is.null(k$per_all)) {
    cat("Calculating periodogram...\n")
    kperiodogram.boot(k, trials=1e4)        
  }
  par(default.par)
  plot(k$per_all, show.resampled=TRUE)
  dev.copy2pdf(file='periodogram.pdf', width=8, height=6)
  dev.copy2png(file='periodogram.png', units='in', width=8, height=6, res=300)

  if (is.null(k$per_res)) {
    cat("Calculating periodogram of residuals...\n")
    k$per_res <- kperiodogram.boot(k, 'res', trials=trials)
  }
  par(default.par)
  plot(k$per_res, show.resampled=TRUE)
  dev.copy2pdf(file="residuals.pdf", width=8, height=6)
  dev.copy2png(file='residuals.png', units='in', width=8, height=6, res=300)
  par(default.par)
  cat(file='residuals.tex', paste(toLatex(xtable(k$per_res, caption='Residuals periodogram')), collapse='\n'))

  if (is.null(k$integration)) {
    cat("Integrating for 10,000 years...\n")
    k$integration <- kintegrate(k, integration * 365.25, RK89)
    saveRDS(k$integration, 'integration.rds')
  }

  par(default.par)
  plot(k$integration, show.resampled=TRUE)
  dev.copy2pdf(file="integration.pdf", width=9, height=12)
  par(default.par)

  par(default.par)
  plot(k)
  dev.copy2pdf(file="fit.pdf", width=8, height=12)

  par(default.par)
  plot(k, plot.rvline=FALSE)
  dev.copy2pdf(file="data.pdf", width=8, height=6)
  
  par(default.par)
  plot(k, 'allrv', wrap=TRUE, plot.residuals=FALSE)
  dev.copy2pdf(file="fit_wrap.pdf", width=4 * (floor(k$nplanets/4)+1), height=3 * min(k$nplanets, 4))
  par(default.par)

  if (is.null(k$errors)) {
    par(default.par)
    par(mar=c(0, 0, 0, 0))
    plot.orbit(k)
    dev.copy2pdf(file='orbit.pdf', width=8, height=8)
    par(default.par)
  } else {
    par(default.par)
    par(mar=c(0, 0, 0, 0))
    plot.orbit(k$errors)
    dev.copy2pdf(file='orbit.pdf', width=8, height=8)
    par(default.par)
  }

  system("pdflatex report.tex")
  ksave(k, 'fit.k')

}

options(xtable.sanitize.text.function=function(x) x)

xtable.periodogram <- function(p, ...) {
  p <- attr(p, 'peaks')
  p <- p[1:min(nrow(p), 10), c(1, 2, 3, 6)]
  display <- rep('f', ncol(p)+1)
  display[4] <- 'e'
  display[1] <- 'd'
  return(xtable(p, display=display, ...))
}

adrange <- function(p) abs(diff(range(p)))

