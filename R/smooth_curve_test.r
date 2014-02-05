karea <- function(x1, y1, x2, y2, xn, yn) {
    x1 <- x1-xn; y1 <- y1 - yn
    x2 <- x2-xn; y2 <- y2 - yn

    return(abs(x1*y2-x2*y1))
}

indexes <- c()
if (!exists('knew')) {
    if (! file.exists('R/systemic.r'))
        setwd('..')
    source('R/systemic.r', chdir=TRUE)
}


# Re-sample curve based on a triangle area criterion.
# Adapted from http://ariel.chronotext.org/dd/defigueiredo93adaptive.pdf
kreduce.curve <- function(m, a = 1, b = nrow(m), tol=1e-5 * (max(m[, 1]) - min(m[, 1])) *  (max(m[, 2]) - min(m[, 2]))) {
 
    if (a == b || b == a + 1)
        return
    n = NULL

    if (is.null(n))
        n = trunc(runif(1, 0.45, 0.55) * (b-a) + a)
    if (n < 1 || n > nrow(m))
        stop(sprintf("%d %d %d", a, b, n))
    
    if (karea(m[a, 1], m[a, 2], m[b, 1], m[b, 2], m[n, 1], m[n, 2]) > tol) {
        kreduce.curve(m, a, n)
        kreduce.curve(m, n, b)
    } else {
        indexes <<- c(indexes, n)
    }

    
    return(invisible())
}

find.peaks <- function(m, max) {
    peaks <- c()
    for (i in 2:(nrow(m)-1)) {
        if ((abs(m[i, 2]) > abs(m[i-1, 2])) && (abs(m[i, 2]) > abs(m[i+1, 2])))
            peaks <- c(peaks, i)
    }
    return(peaks)
}

k <- knew()
kload.datafile(k, 'private/datafiles/HD115617.sys')

k[] <- c(4)
k[] <- c(10, ecc=0.9)
#b <- krvcurve(k, times=seq(k$epoch, k$epoch + 40, length.out=5000))
b <- kperiodogram(k, samples=40000)


par(mfrow=c(2, 1))
plot(b[, 1], b[, 2], type='l', xlab="Power", ylab="Log Period [d]", xlim=c(10^-0.3, 10^4), log='x')
#title("Original periodogram (40000 samples)")

print(diff(range(b[,1])) * diff(range(b[, 2])))
peaks <- find.peaks(b, 3.*sd(b[, 2]))
peaks <- c(1, peaks)
peaks <- c(peaks, nrow(b))
#print(peaks)


b2 <- .R_to_gsl_matrix(b)

t <- seq(1e-6, 1e-3, length.out=100)

v <- double(1)
v[1] <- 1e-3

ret <- .gsl_matrix_to_R(ok_resample_curve(b2, 0, 1, 0.1, 800, 50, v, 10, TRUE))
print(v[1])
print(nrow(ret))
plot(ret[, 1], ret[, 2], type='l', xlab="Power", ylab="Log Period [d]", xlim=c(10^-0.3, 10^4), log='x')
