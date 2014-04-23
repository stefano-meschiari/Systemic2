
kde <- function(k, n=50, nan.range = 1e3, check.stability=NA, dry.run=FALSE, log=c(PER, MASS), print=c(PER, MASS), F=1) {

    fitness <- function(k) {
        kcalculate(k)
        f <- k$chi2
        if (!is.na(check.stability)) {
            i <- kintegrate(k, times=seq(k$epoch, k$epoch + check.stability, length.out=100), int.method=RK89)

            for (j in 1:k$nplanets) {
                ma <- mad(i$els[[j]][,PER])
                me <- median(i$els[[j]][, PER])
                f <- f * if (is.na(ma/me)) 1 else ma/me
            }
        }
        return(f)
    }
    
    klist <- list()
    flags <- kflags(k)
    
    for (i in 1:n) {
        k2 <- kclone(k)

        for (jj in 1:nrow(flags$els))
            for (kk in 1:ncol(flags$els))
                if (flags$els[jj, kk] == 6) {
                    range <- krange(k, jj, kk)
                    if (is.nan(range[1]))
                        range[1] = -nan.range
                    if (is.nan(range[2]))
                        range[2] = -nan.range
                    if (kk %in% log)
                        k2[jj, kk] = 10^runif(1, log10(range[1]), log10(range[2]))
                    else
                        k2[jj, kk] = runif(1, range[1], range[2])
                }
        k2$fitness <- fitness(k2)
        print(i)
        klist[[i]] <- k2
        
    }

    
    if (dry.run)
        return(klist)

    
    
}
