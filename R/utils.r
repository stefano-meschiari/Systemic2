write.fmatrix <- function(m, file="", col.names=TRUE, format="%18.10e", sformat="%18s") {
	f <- file(file, "w")

	if (col.names && !is.null(colnames(m)))
		writeLines(Reduce(paste, sprintf(sformat, colnames(m))), con=f)

	for (r in 1:nrow(m)) 
		writeLines(Reduce(paste, sprintf(format, m[r, ])), con=f)
	
	close(f)
}