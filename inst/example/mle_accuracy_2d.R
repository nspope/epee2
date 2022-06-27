#--------------------#
# compare accuracy of mle using loglikelihood approximation for 2d binary/binary problem
#--------------------#

#NATE: the unconstrained model is unidentifiable, and I haven't been able to figure out how to make it identifiable!

#devtools::install_github("nspope/epee2")
#devtools::install_github("nspope/hiscott2d")

library(ape)
library(epee2)

#simulate data
set.seed(2)
num_tips <- 500
num_traits <- 2

tree <- ape::reorder.phylo(ape::rtree(num_tips, br=function(n) runif(n, 0, 0.05)), "postorder") #must be in postorder ...

phylo_cov <- diag(num_traits) #rWishart(1, num_traits, diag(num_traits))[,,1]
error_cov <- diag(num_traits) #i think this needs to be fixed for identifiability?
ancestral <- c(-0.5, 0.5) #rnorm(num_traits)

covariance <- kronecker(ape::vcv(tree), phylo_cov) + kronecker(diag(num_tips), error_cov)
mean <- c(kronecker(rep(1,num_tips), ancestral))
trait_liability <- matrix(MASS::mvrnorm(1, mean, covariance), num_traits, num_tips)
trait_data <- ifelse(trait_liability >= 0, 1, 0)

#starting values
phylo_cov_start <- log_cholesky_parameterization(diag(num_traits))
error_cov_start <- log_cholesky_parameterization(diag(num_traits))
ancestral_start <- rep(0, num_traits)
start <- c(phylo_cov_start[lower.tri(phylo_cov_start, diag=TRUE)], error_cov_start[lower.tri(error_cov_start, diag=TRUE)], ancestral_start)

#functions for optimization (functions autodetect if trait is continuous/binary)
loss <- function(par)
{
  n_cov_par <- num_traits * (num_traits + 1) / 2
  lambda_idx <- 1:n_cov_par
  sigma_idx <- n_cov_par + 1:n_cov_par
  mean_idx <- 2*n_cov_par + 1:num_traits
  loglik <- epee_objective_fn(trait_data, tree, 
    Sigma_log_cholesky=par[sigma_idx],
    Lambda_log_cholesky=par[lambda_idx],
    mu=par[mean_idx],
    tol=1e-8, max_iter=100, verbose=FALSE)
  -loglik$loglikelihood
}

gradient <- function(par)
{
  n_cov_par <- num_traits * (num_traits + 1) / 2
  lambda_idx <- 1:n_cov_par
  sigma_idx <- n_cov_par + 1:n_cov_par
  mean_idx <- 2*n_cov_par + 1:num_traits
  loglik <- epee_objective_fn(trait_data, tree, 
    Sigma_log_cholesky=par[sigma_idx],
    Lambda_log_cholesky=par[lambda_idx],
    mu=par[mean_idx],
    tol=1e-8, max_iter=100, verbose=FALSE)
  -c(loglik$gradient_Lambda_log_cholesky, loglik$gradient_Sigma_log_cholesky, loglik$gradient_mu)
}

#double check that gradient is right -- any discrepancy should get smaller as "tol" gets smaller
numDeriv::grad(loss, start)
gradient(start)

#fit model
fit <- optim(start, fn=loss, gr=gradient, method="L-BFGS-B", hessian=TRUE)
mm <- matrix(0, 2, 2); mm[lower.tri(mm, diag=TRUE)] <- fit$par[1:3]; diag(mm) <- exp(diag(mm)); phylo_cov_est <- mm %*% t(mm)
mm <- matrix(0, 2, 2); mm[lower.tri(mm, diag=TRUE)] <- fit$par[4:6]; diag(mm) <- exp(diag(mm)); error_cov_est <- mm %*% t(mm)
ancestral_est <- fit$par[7:8]

##double check loss against other implementations
#original api
loss(fit$par)

#alternative api without log-Cholesky parameterization
alt_api <- ThresholdModel$new(tree$edge-1, tree$edge.length, trait_data, c(1, 1))
alt_api$loglikelihood(phylo_cov_est, error_cov_est, matrix(ancestral_est, num_traits, num_tips), 1e-8, 100)

#hiscott method (might need to up level to get accurate answer)
hiscott2d(trait_data, tree, Lambda=phylo_cov_est, Sigma=error_cov_est, mu=ancestral_est, level=10)

#this is not identifiable, can't proceed until I figure out what constraints are needed
