write.f <- function(m, file="", col.names=TRUE, format="%18.10e", sformat="%18s", comments=NULL) {
    f <- file(file, "w")

    if (!isOpen(file))
        stop(sprintf("Could not open %s", file))

    if (!is.null(comments)) 
        writeLines(sprintf("# %s", comments))
    
    if (is.matrix(m)) {
        if (col.names && !is.null(colnames(m)))
            writeLines(Reduce(paste, sprintf(sformat, colnames(m))), con=f)

        for (r in 1:nrow(m)) 
            writeLines(Reduce(paste, sprintf(format, m[r, ])), con=f)
    } else if (is.vector(m)) {
        writeLines(sprintf(format, m), con=f)        
    } else {
        stop(sprintf("I don't know how to write out object of class %s; try to use the write, dump or save functions instead"))
    }
    close(f)
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
