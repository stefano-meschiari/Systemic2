karea <- function(x1, y1, x2, y2, xn, yn) {
    x1 <- x1-xn; y1 <- y1 - yn
    x2 <- x2-xn; y2 <- y2 - yn

    return(abs(x1*y2-x2*y1))
}

indexes <- c()
if (!exists('knew')) {
    setwd('..')
    source('R/systemic.r', chdir=TRUE)
}

# Re-sample curve based on a triangle area criterion.
# Adapted from http://ariel.chronotext.org/dd/defigueiredo93adaptive.pdf
kreduce.curve <- function(m, a = 1, b = nrow(m), tol=1e-5 * (max(m[, 1]) - min(m[, 1])) *  (max(m[, 2]) - min(m[, 2]))) {
    if (a == b || b == a + 1)
        return 
    n = trunc(runif(1, 0.45, 0.55) * (b-a) + a)
    if (n < 1 || n > nrow(m))
        stop(sprintf("%d %d %d", a, b, n))
    
    if (karea(m[a, 1], m[a, 2], m[b, 1], m[b, 2], m[n, 1], m[n, 2]) > tol) {
        kreduce.curve(m, a, n)
        kreduce.curve(m, n, b)
    } else {
        kadd.index(n)
    }

    
    return(invisible())
}

k <- knew()
kload.datafile(k, 'datafiles/51peg_B06L.sys')

k[] <- c(4)
k[] <- c(10, ecc=0.9)
#b <- krvcurve(k, times=seq(k$epoch, k$epoch + 40, length.out=5000))
b <- kperiodogram(k)
b[, 1] <- log10(b[, 1])
par(mfrow=c(2, 1))
plot(b[, 1], b[, 2], type='l')
print(diff(range(b[,1])) * diff(range(b[, 2])))
print(kreduce.curve(b))
indexes <<- sort(indexes)
print(length(indexes))
plot(b[indexes, 1], b[indexes, 2], type='o')
