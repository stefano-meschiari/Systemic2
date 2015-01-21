phases <- function(k, from, to, phases=c(0, 0.25, 0.5, 0.75), save=NULL, print=TRUE,  plIndex = 1) {
    .check_kernel(k)
    k1 <- kclone(k)

    P <- k[plIndex, 'period']
    M <- seq(0, 360, length.out=2e4)
    f <- c()
    for (m in M) {
        k1[1, 'ma'] <- m
        f <- c(f, k1[plIndex, 'trueanomaly'])
    }

    
    t <- M * P/360 + k[1, 'tperi']

    phases.deg <- phases * 360 + 90 - k[1, 'lop']
    phases.deg[phases.deg < 0] <- phases.deg[phases.deg < 0] + 360

    times <- c()
    for (phase.deg in phases.deg) {
        min <- which.min(abs(phase.deg-f))
        times <- c(times, t[min])
    }

    if (!is.null(save))
        save <- file(save, 'w')
    
    diff <- (from %% P) * P
    out <- list()

    if (print)
        cat(sprintf("# Between JD = %f and JD = %f:\n", from, to))
    
    for (i in 1:length(phases)) {
        phase <- phases[i]
        t.phase <- c()
        t <- times[i] + diff

        while (t <= to) {
            if (t >= from) {
                t.phase <- c(t.phase, t)
            }
            t <- t + P
        }
        
        out[[paste(phase)]] <- t.phase
        if (!is.null(save))
            cat(sprintf("# %f\n", phase), sprintf("%10.8e\n", t.phase), file=save)

        if (print) {
            cat(sprintf("# phase = %f\n", phase), sprintf("%10.5f\n", t.phase))
            cat("\n")
        }
    }
    if (!is.null(save))
        close(save)
    
    return(invisible(out))
}
