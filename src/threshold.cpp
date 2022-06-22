#include "epee.h"
#include "covariance.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// cleaner interface:
// simple wrapper around likelihood, gradient

struct ThresholdModel
{
  epee::Covariance covariance;

  ThresholdModel (arma::umat edge_matrix, arma::vec edge_length, arma::mat trait_data, arma::uvec trait_type)
    : covariance (build_covariance_matrix(edge_matrix, edge_length, trait_data, trait_type))
  {}

  epee::Covariance build_covariance_matrix (arma::umat edge_matrix, arma::vec edge_length, arma::mat trait_data, arma::uvec trait_type)
  {
    if (edge_matrix.n_rows != edge_length.n_elem)
    {
      Rcpp::stop("Must have edge length for each edge");
    }

    if (trait_data.n_rows != trait_type.n_elem)
    {
      Rcpp::stop("Must have type (0 for continuous, 1 for binary) for each row in trait matrix");
    }

    if (edge_matrix.min() > 0)
    {
      Rcpp::stop("Edge matrix must use 0-based indices");
    }

    epee::Covariance out (edge_matrix, edge_length, trait_data, trait_type);
    return out;
  }

  Rcpp::List loglikelihood (arma::mat phylogenetic_covariance, arma::mat error_covariance, arma::mat mean, const double tol = 1e-4, const unsigned max_iter = 100)
  {
    covariance.parameters(phylogenetic_covariance, error_covariance, mean);
    covariance.traverse_tree_marginal();

    epee::Epee ep_approximation (covariance);
    ep_approximation.options(tol, max_iter, true);
    ep_approximation.run(covariance.Z, covariance);

    return Rcpp::List::create(
        Rcpp::_["loglikelihood"] = ep_approximation.Z.tail(1),
        Rcpp::_["gradient_phylogenetic_covariance"] = ep_approximation.dL,
        Rcpp::_["gradient_error_covariance"] = ep_approximation.dS,
        Rcpp::_["gradient_mean"] = ep_approximation.dm,
        Rcpp::_["trajectory"] = ep_approximation.Z,
        Rcpp::_["converged"] = ep_approximation.converged);
  }
};

RCPP_EXPOSED_CLASS(ThresholdModel)

RCPP_MODULE(threshold) {
  using namespace Rcpp;

  class_<ThresholdModel>("ThresholdModel")
    .constructor<arma::umat, arma::vec, arma::mat, arma::uvec>()
    .method("loglikelihood", &ThresholdModel::loglikelihood)
    ;
}
