#--------------------#
# compare accuracy of mle using loglikelihood approximation for 1 binary trait
#--------------------#

#devtools::install_github("nspope/epee2")
#devtools::install_github("nspope/hiscott2d")

library(ape)
library(epee2)

set.seed(2)
num_tips <- 500
replicates <- 100

phylo_cov <- 1
error_cov <- 1 #fixed for identifiability
ancestral <- -0.5

sims_1d <- c()
for (i in 1:replicates)
{
  #simulate data
  tree <- ape::reorder.phylo(ape::rtree(num_tips), "postorder") #must be in postorder ...
  tree$edge.length <- tree$edge.length / max(node.depth.edgelength(tree)) #scale tree
  covariance <- ape::vcv(tree)*phylo_cov + diag(num_tips)*error_cov
  mean <- rep(ancestral,num_tips)
  trait_liability <- matrix(MASS::mvrnorm(1, mean, covariance), 1, num_tips)
  trait_data <- ifelse(trait_liability >= 0, 1, 0)
  
  #starting values
  phylo_cov_start <- log_cholesky_parameterization(matrix(1))
  error_cov_start <- log_cholesky_parameterization(matrix(1))
  ancestral_start <- rep(0, num_traits)
  start <- c(ancestral_start, phylo_cov_start, error_cov_start)
  
  #functions for optimization (functions autodetect if trait is continuous/binary)
  loss <- function(par)
  {
    loglik <- epee_objective_fn(trait_data, tree, 
      Sigma_log_cholesky=par[3],
      Lambda_log_cholesky=par[2],
      mu=par[1],
      tol=1e-8, max_iter=100, verbose=FALSE)
    -loglik$loglikelihood
  }
  gradient <- function(par)
  {
    loglik <- epee_objective_fn(trait_data, tree, 
      Sigma_log_cholesky=par[3],
      Lambda_log_cholesky=par[2],
      mu=par[1],
      tol=1e-8, max_iter=100, verbose=FALSE)
    loglik$gradient_Sigma_log_cholesky <- 0 #identifiability: fix error covariance
    -c(loglik$gradient_mu, loglik$gradient_Lambda_log_cholesky, loglik$gradient_Sigma_log_cholesky)
  }
  ##double check that gradient is right
  #numDeriv::grad(loss, start)
  #gradient(start)
  
  #fit model
  fit <- optim(start, fn=loss, gr=gradient, method="L-BFGS-B", hessian=TRUE)
  ancestral_est <- fit$par[1]
  phylo_cov_est <- exp(fit$par[2])^2 #log-cholesky so just log(sqrt(var)) for 1d case
  
  #double check loss against other implementations
  #need to be sure that degree ("level") is high enough for hiscott to be accurate
  loglik_epee <- -loss(fit$par)
  loglik_hiscott <- hiscott1d(trait_data, tree, Lambda=phylo_cov_est, Sigma=error_cov_est, mu=ancestral_est, level=50)#50 is slow but accurate

  #optimize using exact likelihood
  #this is the bottleneck, takes awhile
  loss_hiscott <- function(par)
  {
    -hiscott1d(trait_data, tree, Lambda=exp(par[2])^2, Sigma=error_cov, mu=par[1], level=20) #takes forever if level is high
  }
  fit_hiscott <- optim(fit$par, fn=loss_hiscott, method="L-BFGS-B", hessian=TRUE)
  ancestral_est_hiscott <- fit_hiscott$par[1]
  phylo_cov_est_hiscott <- exp(fit_hiscott$par[2])^2
  
  sims_1d <- rbind(sims_1d, 
    data.frame(num_tips=num_tips, 
               ancestral=ancestral, 
               phylo_cov=phylo_cov, 
               ancestral_est_epee=ancestral_est, 
               phylo_cov_est_epee=phylo_cov_est, 
               ancestral_est_hiscott=ancestral_est_hiscott, 
               phylo_cov_est_hiscott=phylo_cov_est_hiscott, 
               loglik_epee=loglik_epee,
               loglik_hiscott=loglik_hiscott)
  )
  cat(i, "\n")
}

#plots
library(ggplot2)
library(cowplot)

ggplot(sims_1d) + 
  geom_point(aes(x=ancestral_est_epee, y=phylo_cov_est_epee/(phylo_cov_est_epee + error_cov)), size=1) +
  annotate(geom="point", x=ancestral, y=phylo_cov/(phylo_cov + error_cov), color="red", size=3) +
  annotate(geom="text", x=ancestral, y=phylo_cov/(phylo_cov + error_cov), color="red", label="Truth", hjust=-0.25) +
  theme_cowplot() + ylim(0, 1) +
  xlab("Ancestral liability") + ylab("Phylogenetic signal") -> fig_mle_scatter

ggplot(sims_1d) + 
  geom_abline(intercept=0, slope=1) +
  geom_point(aes(y=ancestral_est_epee, x=ancestral_est_hiscott), size=1) +
  theme_cowplot() +
  ylab("Ancestral liability (approx MLE)") +
  xlab("Ancestral liability (exact MLE)") -> fig_acc_ancestral

ggplot(sims_1d) + 
  geom_abline(intercept=0, slope=1) +
  geom_point(aes(y=phylo_cov_est_epee/(phylo_cov_est_epee + error_cov), x=phylo_cov_est_hiscott/(phylo_cov_est_hiscott + error_cov)), size=1) +
  theme_cowplot() +
  ylab("Phylogenetic signal (approx MLE)") +
  xlab("Phylogenetic signal (exact MLE)") -> fig_acc_signal

ggplot(sims_1d) + 
  geom_abline(intercept=0, slope=1) +
  geom_point(aes(y=loglik_epee, x=loglik_hiscott), size=1) +
  theme_cowplot() +
  ylab("Loglikelihood (approx)") +
  xlab("Loglikelihood (exact)") -> fig_acc_loglik

png("1d_approximation_accuracy.png", height=4, width=16, units="in", res=300)
cowplot::plot_grid(fig_mle_scatter, fig_acc_ancestral, fig_acc_signal, fig_acc_loglik, ncol=4, align="h")
dev.off()
