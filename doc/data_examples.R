# Use kdata(k) to access data associated with a kernel object.
# kdata returns a matrix where each row is a data point, and 
# each column represents, respectively:
# TIME [1]: the time stamp
# VAL [2]: the value associated with the data point (e.g. RV value)
# ERR [3]: the measured uncertainty (as reported by the data file)
# FLAG [7]: a flag indicating the type of data point (one of RV or TIMING)
# SVAL [8]: the value stored in VAL + any vertical shifts or trends (e.g.
# with the RV offset added)
# PRED [9]: the predicted value at that time stamp (e.g. the stellar RV
# response velocity) for the current fit values
# SET [10]: an index indicating the data-set containing the specific data
# point
# TDS_PLANET [4]: planet associated with a transit timing point (not
# used for RV data)
# TDS_FLAG [5]: primary transit or secondary eclipse? (not used for RV data)

# Data manipulation can be carried out using the R syntax, which is
# described here:
# http://cran.r-project.org/doc/manuals/R-intro.html#Arrays-and-matrices
# or in any elementary R book.

# A few examples are listed in the "Data" menu, and below:

# Normalized residuals
nr <- (kdata(k)[, PRED] - kdata(k)[, SVAL])/kdata(k)[, ERR]

# Print outliers that have normalized residuals > 5
print(nr[nr > 5, ])

# Remove data where normalized residuals are > 5
kdata(k) <- kdata(k)[nr < 5, ]

# Select a temporal subset of the data before JD = 2454047
subs <- kdata(k)[kdata(k)[, TIME] < 2454047, ]

# Select data from datasets 1 and 2, but not 0
d12 <- kdata(k)[kdata(k)[, SET] = 1 | kdata(k)[, SET] = 2, ]

# Add 2m/s error in quadrature for dataset # 1
d <- kdata(k, idx=1)
d[, ERR] <- sqrt(d[, ERR]^2 + 2^2)
kdata(k, idx=1) <- d