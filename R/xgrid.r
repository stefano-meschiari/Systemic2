.random.string <- function(len) {
	return(paste(sample(letters, len, replace=TRUE), collapse=""))
}

.xgrid.submit <- function(cmd, path, job, n, host, password) {
	fn <- paste(job, n, ".sh", sep="")
	cat(file=paste(path, fn, sep=""), sprintf("#!/bin/bash\n%s\n", cmd))
	system(sprintf("chmod +x %s", paste(path, fn, sep="")))
	
	cmd <- sprintf("cd \"%s\"; xgrid -h %s -p %s -in '.' -out '.' -job submit %s &",
		path, host, password, fn)
		
	ret <- system(cmd, intern=TRUE)

	if (grepl("error", paste(ret, collapse=""))) {
		stop(paste(ret, collapse=" "))
	}
	
	ret <- strsplit(sub(";", "", ret[[2]], fixed=TRUE), '=', fixed=TRUE)
	return(as.numeric(ret[[1]][2]))
}

.xgrid.wait <- function(id, host, password, dir) {
	ret <- system(sprintf("xgrid -h %s -p \"%s\" -job wait -id %s", host, password, id), intern=TRUE)
	pasted <- paste(ret, collapse="")
	if (grepl("Canceled", pasted)) 
		return("Canceled")
	else if (grepl("Failed", pasted))
		return("Failed")
	else {
		cmd <- sprintf("cd '%s'; xgrid -h %s -p \"%s\" -out '.' -job results -id %s", dir, host, password, id)

		ret <- system(cmd, intern=TRUE)
		
		return("Finished")
	}
}

kbootstrap.xgrid <- function(k, xgrid.host=getOption("systemic.xgrid.host", NA), xgrid.password=getOption("systemic.xgrid.password", NA), algo = NA, trials.each = 500, jobs = 10, warmup = 0, min_iter = 2000, plot = FALSE, print = FALSE, save=NA, systemic.dir=getOption("systemic.dir", "."),  job.name="bootstrap", save.script=NA) {
	.check_kernel(k)
	stopifnot(k$ndata > 0)	
	.job <<- "Bootstrap on Xgrid"
	
	if (is.na(xgrid.host) || is.na(xgrid.password))
		stop("Specify xgrid.host and xgrid.password\nUse option(systemic.xgrid.host='hostname', systemic.xgrid.password='pwd')")
		
	if (is.na(algo))
		algo <- k$min.method

	if (! file.exists("utils") || !file.exists("utils/systemic_cli")) {
		stop("Could not find utils/systemic_cli; specify the directory where systemic is installed using the argument systemic.dir")
	}
	
	dirs <- c()
	ids <- c()
	for (job in 1:jobs) {
		dir <- paste(tempdir(), "/", job, "_", .random.string(20), "/", sep='')
		dir.create(dir)
		system(paste("cp -R utils/* ", dir))
		ksave(k, paste(dir, "kernel.k", sep=""))
		
		id <- .xgrid.submit(paste("./systemic_cli -file kernel.k -action load -trials ", trials.each, " -algo ", algo, " -iters ", min_iter, " -save kl.kl -action bootstrap", sep=""), dir, job.name, job, xgrid.host, xgrid.password)
		dirs[job] <- dir
		ids[job] <- id
		
		cat(paste("Job ", job, " submitted, id = ", id, "\n", sep=""))
	}
	
	kls <- NULL
	success <- c()
	for (job in 1:jobs) {
		ret <- .xgrid.wait(ids[job], xgrid.host, xgrid.password, dirs[job])
		cat(paste("Job ", job, ", id = ", ids[job], ", status = ", ret, "\n", sep=""))
		success[job] <- FALSE
		
		for (iter in 1:10) {
			if (file.exists(paste(dirs[job], "kl.kl", sep=""))) {
				fn <- paste(dirs[job], "kl.kl", sep="")
				fid <- fopen(fn, "r")
				if (!is.nullptr(fid)) {
					if (is.null(kls))
						kls <- KL_load(fid, 0)
					else
						KL_append(kls, KL_load(fid, 0))
					fclose(fid)
				}
				success[job] <- TRUE
				fclose(fid)
				break;
			} 
		}
		
		if (! success[job]){
			cat(sprintf("Could not retrieve results for job %d [%d], reading others\n", job, ids[job]))
		}		
	}
	
	if (is.null(kls)) {
		stop("No jobs returned successfully")
	}
		
	a <- .klnew(kls, k, type="bootstrap.xgrid", desc=sprintf("algo = %d, trials = %d", algo, trials.each*jobs))
	
	a$jobs <- jobs
	a$ids <- ids
	a$dirs <- dirs
	
	if (plot)
		plot(a)
	if (print)
		print(a)
	if (!is.na(save))
		save(a, file=save)
		
	return(a)
}