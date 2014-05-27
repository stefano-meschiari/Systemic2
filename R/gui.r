source("systemic.r")

if (!exists(".gui")) {
    stop("This file is not meant to be used by a regular R session. Start the Systemic GUI instead.")
}

options(systemic.url="http://www.stefanom.org/")
options(help_type="html")
options(max.print=200)
options(systemic.autosave=FALSE)
options(digits=10)

.gui.version <- SYSTEMIC.VERSION
.gui.env <- new.env()
.gui.silent <- FALSE
.gui.no.undo <- FALSE
.gui.debug.printevents <- FALSE
.gui.editor <- 'mate'
.gui.search.path <- c(0)
exec <- tryCatch

.gui.REQUEST_UPDATE <- K_GUI_RESERVED_1
.gui.os <- Sys.info()['sysname']
.gui.inline.plotting <- FALSE
.gui.inline.width <- 600
.gui.inline.height <- 600

if (.gui.os == "Linux") {
	options(device='x11')
} else {
	options(device='quartz')
}

.gui.event <- function(event, name = NA, size1 = NA, size2 = NA) {
	if (.gui.debug.printevents)
		cat(sprintf("%s,%s\n", event, name))
		
	if (is.na(name))
		cat(sprintf("@@@@`%s\n", event))
	else
		cat(sprintf("@@@@`%s`%s\n", event, name))
	
	if (! is.na(size1) && ! is.na(size2)) {
		cat(sprintf("%d %d\n", size1, size2))
	} else if (! is.na(size1)) {
		cat(sprintf("%d\n", size1))
	}
}

.gui.matrix <- function(m, free = F) {
	if (is.null(m) || is.nullptr(m)) {
		cat(-1, "\n")
	} else {
		if (class(m) == "matrix") {
			m <- .R_to_gsl_matrix(m)
			free <- TRUE
		}
		cat(ok_matrix_rows(m), " ", ok_matrix_cols(m), "\n")
		if (! .gui.silent) 
			ok_save_matrix_bin(m, NULL)
		
		if (free) {
			gsl_matrix_free(m)
		}
	}
}

.gui.matrix.int <- function(m) {

	if (is.null(m) || is.nullptr(m)) {
		cat(-1, "\n")
	} else {
		cat(ok_matrix_rows(m), " ", ok_matrix_cols(m), "\n")
		if (! .gui.silent) 	
			ok_fprintf_matrix_int(m, NULL, "%d ")
	}
}

.gui.vector.int <- function(m) {

	if (is.null(m) || is.nullptr(m)) {
		cat(-1, "\n")
	} else {
		cat(ok_vector_len(m), "\n")
		if (! .gui.silent) 
			ok_fprintf_vector_int(m, NULL, "%d ")
	}
}

.gui.vector <- function(m) {

	if (is.null(m) || is.nullptr(m)) {
		cat(-1, "\n")
	} else {
		cat(ok_vector_len(m), "\n")
		if (! .gui.silent) 
			ok_fprintf_vector(m, NULL, "%e ")
	}
}


.gui.nameof <- function(k) {
	names = c()
	
	for (v in ls(envir=globalenv())) {
			k2 <- get(v, envir=globalenv())
			
			if (class(k2)[1] == "kernel" && identical(k[['h']], k2[['h']])) {
				names = c(names, v)
			}
	}
	if (length(names) == 0)
		return(NULL)
	else
		return(names[1])
}

.gui.write.error <- function(name, obj) {
	con <- file(paste(.gui.path, "_err_", name, sep=""), open="w")
	cat(file=con, sprintf("type\t%s\n", obj$type))
	cat(file=con, sprintf("desc\t%s\n", obj$desc))	
	cat(file=con, sprintf("nplanets\t%s\n", obj$nplanets))
	cat(file=con, sprintf("nsets\t%s\n", obj$nsets))	
	close(con)
}

.gui.write.env <- function(k) {
	con <- file(paste(.gui.path, "_vars", sep=""), open="w")

	for (v in ls(envir=globalenv())) {
			k2 <- get(v, envir=globalenv())
			cat(sprintf("%s\t%s\n", v, class(k2)), file=con)
			
			if (class(k2)[1] == "error.est") {
				.gui.write.error(v, k2)
			} else if (class(k2)[1] == "kernel" && ! is.null(k2$errors)) {
				.gui.write.error(sprintf("%s$errors", v), k2$errors)
				cat(sprintf("%s$errors\t%s\n", v, class(k2$errors)), file=con)
			}
	}
	
	close(con)
}

.gui.autosave.dir <- paste(getOption("systemic.dir", "./"), "/autosave/", sep="")
.gui.last.autosave <- Sys.time()

if (! file.exists(.gui.autosave.dir)) {
    dir.create(.gui.autosave.dir, showWarnings=FALSE)
}
.gui.autosave <- function() {
    if (getOption('systemic.autosave') && (as.numeric(difftime(Sys.time(), .gui.last.autosave), unit="secs") > 60)) {
        .gui.last.autosave <<- Sys.time()
        fn <- paste(.gui.autosave.dir, "session_", .gui.session, sep="")
        .gui.event("#save_session", fn)
        .gui.event("#hint", "Autosaved.")
    }
}

.gui.undo.list = c()
.gui.add.undo <- function(name, k) {
	if (is.null(.gui.undo.list[[name]])) {
		.gui.undo.list[[name]] <<- c(new=kclone(k), old=kclone(k))	
	} else {

		# Try to avoid setting a new undo point when nothing
		# has changed
		old <- .gui.undo.list[[name]]$new
		
		if (isTRUE(all.equal(kallels(old, keep.first=T), kallels(k,  keep.first=T))) && isTRUE(all.equal(kpars(old), kpars(k))) && isTRUE(all.equal(kdata(old), kdata(k)))) {
		} else {
			.gui.undo.list[[name]]$old <<- .gui.undo.list[[name]]$new
			.gui.undo.list[[name]]$new <<- kclone(k)
		}
	}
}



.gui.last.str <- "?"


.gui.check.stop <- function() {
	return(file.exists(paste(.gui.path, "_stop", sep="")))	
}

gui.progress <- function(cur=0, max=100, val=0, job="") {
    .gui.event("progress", sprintf("%d|%d|%s|%s", cur, max, val, job))
    if (.gui.check.stop())
        stop("User stopped the current job.")
}

.gui.progress.create <- function(k) {
    k[['.progress']] <- function(cur, max, state, str) {
        if (is.null(k$silent)) {
            if (! is.nullptr(state)) {
                state <- as.struct(state, "ok_kernel")
                .gui.event("progress", sprintf("%d|%d|%e|%s", cur, max, K_getChi2(state), .job))
            } else {
                if (str == "")
                    str <- .gui.last.str
                else
                    .gui.last.str <<- str
                .gui.event("progress", sprintf("%d|%d|%s|%s", cur, max, str, .job))
            }
        }
       
        
        if (.gui.check.stop()) {
            cat("Stop requested, please wait...\n", file=stderr())
            return(K_PROGRESS_STOP)
        } else {
            if (!is.null(k[['progress']])) {
                ret <- k[['progress']](k, list(cur=cur, max=max, job=.job, str=str))
                if (is.numeric(ret))
                    return(ret)
            }
            return(K_PROGRESS_CONTINUE)
        }
    }
    
    k[['.progress.cb']] <- new.callback("iipZ)i", k[['.progress']])
    return(k[['.progress.cb']])
}





.gui.scratch_model <- NULL

.gui.periodogram.tol <- double(1)
.gui.periodogram.tol[1] <- 1e-5
.gui.rvsignal.tol <- double(1)
.gui.rvsignal.tol[1] <- 1e-4

.gui.update <- function(k, name = NA, calculate=TRUE) {
	if (class(k) == "character" && k == "all") {
		for (v in ls(envir=globalenv())) {
			k2 <- get(v, envir=globalenv())
			if (class(k2)[1] == "kernel") {
				.gui.update(k2, name = v, calculate=TRUE)
			}
		}
		return(invisible())
	}
	flush(stderr())
	
	if (is.na(name)) {
		name <- .gui.nameof(k)
	}
	
	if (is.null(name) || name == "*tmp*") {
		.gui.mark(k, .gui.REQUEST_UPDATE)
		return()
	}
	.gui.unmark(k, .gui.REQUEST_UPDATE)
	
	if ((! .gui.no.undo) && calculate)
		.gui.add.undo(name, k)
	K_setProgress(k$h, .gui.progress.create(k))

	ndata <- k$ndata
	if (calculate) {
		if (ndata > 0) {
			c <- K_getCompiled(k$h)
			.gui.event("data", name, ndata, DATA_SIZE)

			if (! .gui.silent) 
				ok_save_buf_bin(c, NULL, ndata, DATA_SIZE)
			flush(stdout())
			
			if (k$nrvs > 0) {
				p <- kperiodogram(k, samples=getOption("systemic.psamples", 50000), pmin=getOption("systemic.pmin", 0.5), pmax=getOption("systemic.pmax", 2e4), .keep.h=T)
				m <- ok_resample_curve(attr(p, "h"), 0, 1, 1, 10000,
                               2000, .gui.periodogram.tol, 5, TRUE)

				.gui.event("periodogram", name)
				.gui.matrix(m, free=T)
			}
			
		} else {
			.gui.event("data", name, -1)
			.gui.event("periodogram", name, -1)
		}
		rvdata <- kdata(k)
		rvdata <- rvdata[rvdata[, FLAG] == RV, ]
		
		if ((! is.nan(k$epoch)) && nrow(rvdata) > 1) {
			
		

			trange <- c(min(rvdata[,1]), max(rvdata[,1]))
			rvsamples <- getOption("systemic.rvsamples", 30000)
			a <- integer(1)
      dt <- k$dt
      K_setIntDt(k$h, min(k$dt, 0.05 * k[, 'period']))

			rvsignal <- K_integrateStellarVelocity(k$h, trange[1], trange[2], rvsamples, NULL, a)
      K_setIntDt(k$h, dt)
      
			if (a[1] != K_INTEGRATION_SUCCESS)
				print(.integration.errors[a[1]])
			if (is.nullptr(rvsignal)) {
				print("Error during integration")
				.gui.event("rvcurve", name)
				.gui.matrix(NULL)
			} else { 
				m <- ok_resample_curve(rvsignal, 0, 1, 1, 5000,
                     500, .gui.rvsignal.tol, 5, FALSE)
				.gui.event("rvcurve", name)
				.gui.matrix(m, free=T)
				gsl_matrix_free(rvsignal)
			}


			if (k$nrvs > 0) {
				p <- kperiodogram(k, per_type="res", samples=getOption("systemic.samples", 50000), pmin=getOption("systemic.pmin", 0.5), pmax=getOption("systemic.pmax", 2e4), .keep.h = TRUE)
				m <- ok_resample_curve(attr(p, "h"), 0, 1, 1., 10000,
                               2000, .gui.periodogram.tol, 5, TRUE)
        
				.gui.event("per_res", name)
				.gui.matrix(m, free=T)
			} else {
				.gui.event("per_res", name, -1)
			}
		} else {
			.gui.event("rvcurve", name, -1)
			.gui.event("per_res", name, -1)			
		}
		
		.gui.event("update", name)
	}	

	.gui.stats(k, name, calculate = calculate)
	
}


detach('.systemic.functions')
.kupdate <- .systemic.functions$kupdate
.systemic.functions$kupdate <- function(k, calculate=TRUE) {
	K_setProgress(k$h, .gui.progress.create(k))	
	.kupdate(k, calculate=calculate)
	.gui.update(k, calculate=calculate)
	invisible()
}


# Used to hide matrix attributes when printing
.systemic.functions$print.matrix <- function(m, ...) {
	at <- attributes(m)
	attr(m, 'h') <- NULL
	
	print.default(m, ...)
	attributes(m) <- at
}


attach(.systemic.functions)

.gui.mark <- function(k, flag) {
	K_setFlags(k$h, bitOr(K_getFlags(k$h), as.integer(flag)))
}

.gui.unmark <- function(k, flag) {
	K_setFlags(k$h, bitAnd(K_getFlags(k$h), bitFlip(as.integer(flag))))
}

.gui.marked <- function(k, flag) {
	return(bitAnd(K_getFlags(k$h), as.integer(flag)) == as.integer(flag))
}

.gui.stats <- function(k, name = NULL, calculate = TRUE) {
	if (is.null(name))
		name <- .gui.nameof(k)
	if (is.null(name))
		return()
	
	if (calculate) {
		nsets <- k$nsets
		
		con <- file(paste(.gui.path, name, "_stats", sep=""), open="w")
		
		for (field in names(.kdollartable)) {
			cat(sprintf("%s = %.15e\n", field, .kdollartable[[field]](k$h)), file=con)
		}
		for (field in ls(envir=k)) {
			if (class(k[[field]]) == "numeric")
				cat(sprintf("%s = %.15e\n", field, k[[field]]), file=con)
		}
		
		if (nsets > 0)
			for (i in 1:nsets) {
				cat(sprintf("%d!%s = %.15e\n", i, k$datanames[[i]], K_getPar(k$h, i-1)),
				file=con)
			}
	
		for (i in 1:PARAMS_SIZE) {
			cat(sprintf("%d = %.15e\n", i, K_getPar(k$h, i-1)),
				file=con)
		}
		close(con)
	}
		
	.gui.event("stats", name)
	.gui.event("state", name)
	if (! .gui.silent) {
		.gui.matrix(K_getAllElements(k$h), free = T)
		.gui.matrix.int(K_getElementFlags(k$h))
		.gui.vector.int(K_getParFlags(k$h))
		.gui.matrix(K_getElementSteps(k$h))
		.gui.vector(K_getParSteps(k$h))
	}
	
}

. <<- knew()

.gui.before <- function(currentk) {
	tryCatch(. <<- get(currentk, envir=globalenv()), error=function(e) {})
	.job <<- ""
	.gui.stopped = FALSE
	unlink(paste(.gui.path, "_stop", sep=""))
	.gui.event("busy")
	
	
	flush(stdin())
	flush(stdout())
	flush(stderr())
	
	if (.gui.inline.plotting) {
		fn <- paste(.gui.path, "inline_plot.png", sep="")
		if (file.exists(fn)) {
			unlink(fn)
		}
	}
	
	if (length(.gui.hooks[['before']]) > 0)
		for (fun in .gui.hooks[['before']])
			fun()
}

.gui.old = c()

.gui.after <- function() {
	tryCatch({	
		.gui.new = c()
		
		for (v in ls(envir=globalenv())) {
			k2 <- get(v, envir=globalenv())
			
			
			if (class(k2)[1] == "kernel") {
				if (! is.null(k2[['.new']])) {
					k2[['.new']] <- NULL
					k2[['auto']] <- TRUE
					kupdate(k2)
				}
			
				if (.gui.marked(k2, .gui.REQUEST_UPDATE)) {
					.gui.update(k2, name=v)
				}
				if (! is.null(.gui.old[[v]]) && ! identical(.gui.old[[v]][['h']], k2[['h']])) {
					.gui.update(k2, name=v)
				}
				.gui.new[[v]] = k2
			}

		}
		.gui.old <<- .gui.new
		
		if (.gui.inline.plotting) {
			fn <- paste(.gui.path, "inline_plot.png", sep="")			
			tryCatch(dev.print(device=png, width=.gui.inline.width, height=.gui.inline.height, file=fn), error=function(e) {})
			if (file.exists(fn)) {
				.gui.event("#plot", fn)
				dev.off()
			}
		}

    .gui.autosave()
    
		if (length(.gui.hooks[['after']]) > 0)
			for (fun in .gui.hooks[['after']])
				fun()
	}, finally={
		.gui.write.env()
		.gui.event("ready")				
		flush(stdin())
		flush(stderr())
		flush(stdout())
	})
}

.gui.phasedrvs.tol <- double(1)
.gui.phasedrvs.tol[1] <- 1e-3

.gui.phasedrvs <- function(k, name) {
	
	np <- k$nplanets

	if (np < 1 || is.nan(k$epoch) || k$nrvs < 2 || k$int.method != 0)
		return(invisible())
	
	k <- kclone(k)
	k$auto <- FALSE
	ret <- c()

	trange <- k$trange
	masses <- k[,'mass']
	rvsamples <- getOption("systemic.rvsamples", 10000)
	.gui.event("phasedrvs", name)
	
	a <- c()
	
	for (i in 1:np) {
		k[i, 'mass'] <- 0
		
		kcalculate(k)

		k[,'mass'] <- masses
		k[-i, 'mass'] <- 0
		
		sl <- K_integrateRange(k$h, trange[1], trange[2], rvsamples, NULL, k$last.error.code)
		if (is.nullptr(sl)) {
			.gui.matrix(NULL, free=TRUE)	
			.gui.matrix(NULL, free=TRUE)
		} else {
			rvsignal <- ok_get_rvs(sl, rvsamples)
			rvdata <- .buf_to_R(K_compileData(k$h), k$ndata, DATA_SIZE)
			rvdata <- rvdata[rvdata[, FLAG] == RV, ]
			.gui.matrix(rvdata, free=TRUE)	
			.gui.matrix(rvsignal)
			gsl_matrix_free(rvsignal)
			ok_free_systems(sl, rvsamples)
		}
		k[,'mass'] <- masses
	}

	invisible()

}

.gui.add.planet <- function(k, P, M, algo, np, circular = FALSE) {
	.mute(k)
	k[] <- c(period = P, mass = M)
	.flags <- kflags(k)

	kminimize1d(k, np, 'ma', algo=algo)
	kdeselect(k)
	kselect(k, np, 'all')
	
	if (circular) {
		kdeselect(k, np, 'ecc')
		kdeselect(k, np, 'lop')	
	}
	
	kminimize(k, algo=algo)
	kflags(k) <- .flags	
	if (circular) {
		kdeactivate(k, np, 'ecc')
		kdeactivate(k, np, 'lop')		
	}
	.unmute(k)
	kupdate(k)
}

.gui.env$undo <- function(k) {
	if (class(k) == "kernel") {
		k <- .gui.nameof(k)
	}
	
	if (is.null(k) || is.null(.gui.undo.list[[k]]))
		stop("No undo available")
	else {
		assign(k, .gui.undo.list[[k]]$old, envir=globalenv())
		kcalculate(get(k, envir=globalenv()))
		kupdate(get(k, envir=globalenv()))
	}
}

.gui.startup <- function() {
	if (file.exists("startup.r"))
		source("startup.r")
	if (getOption("systemic.installed", FALSE))
		.gui.event("installed")
  .gui.event("#pid", Sys.getpid())
}

.gui.hint <- function(str, func="") {
	.gui.event("#hint", gsub("function", func, str))
}

.gui.input <- function(str) {
	.gui.event("#input", str)
}

.gui.input.atcursor <- function(str) {
    .gui.event("#insert_input", str)
}

.gui.title <- function(str) {
    .gui.event("#frame_title", str)
}

.gui.complete <- function(partial, quotedPartial) {
    
    completions <- apropos(paste("^", quotedPartial, sep=""), ignore.case=FALSE)
    if (length(completions) == 0)
        .gui.hint(paste("No completion for ", partial, sep=""))
    else if (length(completions) == 1)
        .gui.input.atcursor(substr(completions[1], nchar(partial)+1, nchar(completions[1])))
    else {
        prefix <- longest.prefix(completions)
        print(prefix)
        .gui.input.atcursor(substr(prefix, nchar(partial)+1, nchar(prefix)))
        print(completions)
    }       
}

.gui.show <- function(str) {
	.gui.event("#switchSidebar", str)
}

.gui.flash <- function(str) {
	.gui.event("#highlight", str)
}

.gui.select <- function(str) {
	.gui.event("selkernel", str)
}

edit.script <- function(path) {
	.gui.event("edit", path)
}

.cats <- function() {
	.gui.event("catz")
}

attach(.gui.env)

acknowledgments <- function() {
	cat(readLines(paste(getOption("systemic.dir", "./"), "/acknowledgments.txt", sep="")), sep="\n")
  edit.script("License.txt")
}

tutorial <- function(which="rv") {
	.tutorial.fn <<- paste("doc/tutorial/tutorial_", which, ".txt", sep="")
	.tutorial.fn <<- paste(getOption("systemic.dir", "./"), "/", .tutorial.fn, sep="")
	source(paste(getOption("systemic.dir", "./"), "/doc/tutorial/tutorial.r", sep=""))
}

quickstart <- function() {
    edit.script("doc/quickstart.txt")
}

.gui.hooks <- list()
.gui.hooks[['before']] <- list()
.gui.hooks[['after']] <- list()

hook <- function(stage, fun) {
	if (! stage %in% c("before", "after"))
		stop("stage should be one of 'before' or 'after'")
	.gui.hooks[[stage]] <<- c(.gui.hooks[[stage]], fun)
}

systemic.plot.theme <- function(name) {
    col2hex <- function(c) {
        c <- col2rgb(c, alpha=TRUE)
        return(rgb(c['red',], c['green',], c['blue',], c['alpha',], maxColorValue=255))
    }
    systemic.palette <<- get(paste('systemic.theme.', name, sep=''))
    systemic.palette.face <<- get(paste('systemic.theme.', name, '.face', sep=''))
    cat("Theme set to ", name, "\n")
    for (i in 2:length(systemic.palette)) {
        
        .gui.event(paste('#palette', i, sep=""), col2hex(systemic.palette.face[i]))
    }
    .gui.event("#updated_palette")
}
    

.gui.path = ".temp"

options(error=function() {})
dev.off()

.gui.event("start")
.gui.startup()

if (getOption("systemic.installed", FALSE)) {
	.gui.event("clearsession")
	.gui.event("restart")	
}
