# This script shows how to add additional terms to the model (a "custom model")
# In this case, we want to add a quadratic term to represent curvature
# in the radial velocity.

# The custom model receives the kernel, a vector of measurement times and a vector
# of flags (specifying the types of data for each element of the times vector; e.g.
# RV or TIMING). It should return a vector of values of the same size, which will
# be added to the fit model.


# Choose parameter 21 as the quadratic coefficient; any parameters <= 20 are used 
# internally. Notice that a parameter to remove a linear trend is already part of the
# base model (P_TREND) and the offsets represent the constant part of the quadratic model,
# so we only need one coefficient
P_QUADRATIC <- 21

add.quadratic.model <- function(k) {
	k$custom.model <- function(k, times, flags) {
		# Prepare the return vector, initialize to 0
		ret <- rep(0., length(times))
		# For the datapoints that represent radial velocities, add the quadratic term
		ret[flags == T_RV] = k['par', P_QUADRATIC] * (times[flags == T_RV] - k$epoch)^2
		# Return the values
		return(ret)
	}
}