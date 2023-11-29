## mcmc diagnostics -- a convergence diagnostic based on a test for equality of the means of the first and last part of a Markov chain
## if the chain has converged, the two means are equal and Geweke's statistic has an asymptotically standard normal distribution.
z <- gewekediag(pdep_draws.dart)$z
qqnorm(z)
abline(0,1, col='red')

plot(pdep_draws.dart[,which.min(z)], type='l')
plot(pdep_draws.dart[,which.max(z)], type='l')

plot(dart.fit$sigma, type = 'l')
abline(v = 2000, col = 'red')
