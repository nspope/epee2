# Nate Pope 22 June 2022 : I can't get Hiscott method to work with mixed continuous/binary, so let's just work with 2 binary traits
# We can always show accuracy of the mixed case by comparing MLEs to true values across simulations rather than comparing likelihood to "gold standard"

#devtools::install_github("nspope/epee2")
#devtools::install_github("nspope/hiscott2d")

library(ape)
library(epee2)
library(hiscott2d)

#simulate data
set.seed(1024)
num_tips <- 25
tree <- ape::reorder.phylo(ape::rtree(num_tips), "postorder") #must be in postorder ...
trait_data <- matrix(rnorm(num_tips*2), 2, num_tips)
trait_data[] <- as.numeric(trait_data[] > 0)

#likelihood evaluations (functions autodetect trait type)
phylogenetic_covariance <- matrix(c(1.1, 0.4, 0.4, 0.8), 2, 2)
error_covariance <- matrix(c(0.5, -0.1, -0.1, 0.5), 2, 2)
mean <- c(-0.3, 0.3)
convergence_tol <- 1e-8
max_iterations <- 100

#approximate likelihood (will dump warning if doesn't converge)
epee_est <- epee2d(
  data = trait_data, 
  tree = tree, 
  Lambda = phylogenetic_covariance, 
  Sigma = error_covariance, 
  mu = mean, 
  tol = convergence_tol, 
  max_iter = max_iterations
)

#exact likelihood approximated through quadrature (hiscott)
# - the quadrature degree controls the accuracy of approximation -- higher is better but slower
# - so, we run across grid of quadrature degree to monitor convergence
# - it would be good to store timings, too, because hiscott gets suuupppeeer slow with high degree / many tips
# - unfortunately, it's going to be hard to use 2d hiscott for large trees b/c the necessary degree is very high ...
#   so we may just have to use 1d example to show accuracy for big trees.
quadrature_degree_grid <- seq(3,20,1)
hiscott_est <- sapply(quadrature_degree_grid, 
  function(level)
  hiscott2d(
    data = trait_data, 
    tree = tree, 
    Lambda = phylogenetic_covariance, 
    Sigma = error_covariance, 
    mu = mean,
    level = level
  )
)

#exact likelihood approximated by genz-bretz (only accurate for small problems, use to validate hiscott)
genz_est <- genz2d(
  data = trait_data, 
  tree = tree, 
  Lambda = phylogenetic_covariance, 
  Sigma = error_covariance, 
  mu = mean
)

#compare estimators
library(ggplot2)
library(cowplot)
df <- data.frame(degree=quadrature_degree_grid, hiscott=hiscott_est)
ggplot(df) + 
  geom_line(aes(x=quadrature_degree_grid, y=hiscott, color='hiscott')) +
  geom_hline(data=data.frame(y=epee_est), aes(yintercept=y, color='epee')) +
  geom_hline(data=data.frame(y=genz_est), aes(yintercept=y, color='genz')) +
  scale_color_manual("Algorithm", values=c('epee'='dodgerblue', 'genz'='black', 'hiscott'='coral')) +
  theme_cowplot() + 
  theme(plot.background=element_rect(fill="white")) + 
  xlab("Quadrature degree") + ylab("Loglikelihood")


