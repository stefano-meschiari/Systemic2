.next <- 0
.tutorial <- c()

.wait <- function() {
	cat("# Press RETURN to continue the tutorial\n")
	readLines(n = 1)
}

STOP <- function() {
	.gui.flash("clear")
	cat("The tutorial is over.")
	invisible()
}

NEXT <- function() {
	if (.next == 0) {
		.next <<- 1
		if (! file.exists(.tutorial.fn))
			stop(paste("Could not find ", .tutorial.fn))
		
		.tutorial <<- readLines(.tutorial.fn)
	}
	.gui.event("clearsession")
	.gui.input("")
	while (.next <= length(.tutorial)) {
		line <- .tutorial[[.next]]
		.next <<- .next + 1
		if (is.null(line))
			break()
		if (substr(line, 1, 1) == "#") {
			cat("\n> ", substr(line, 2, nchar(line)), "\n")
			.gui.hint("Type STOP() to stop the tutorial")
			.gui.event("scrollup")
			return(invisible())
		} else if (substr(line, 1, 1) == "@")
			eval(parse(text = substr(line, 2, nchar(line))))
		else
			cat(line, "\n")
	} 
	
	cat("# The tutorial is over.\n")
	.gui.flash("clear")	
	invisible()
}

NEXT()