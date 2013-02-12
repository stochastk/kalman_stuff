library(inline)
library(Rcpp)
library(astsa)

#prior values
mu0 <- matrix(c(0,3,1), ncol = 1)
Sigma0 <- diag(c(.1,.1,1),3)
Phi <- diag(1,3)
cQ <- diag(c(.1,.1,1),3)
cR <- diag(c(.1,.1,1),3)
Q <- t(cQ)%*%(cQ)
R <- t(cR)%*%(cR)

#this is only to see if i can instantiate a Kalman object 
src <- 'arma::mat priMu = as<arma::mat>(mu0);
        arma::mat priSigma = as<arma::mat>(Sigma0);
        arma::mat priPhi = as<arma::mat>(Phi0);
        arma::mat priQ = as<arma::mat>(Q0);
        arma::mat priR = as<arma::mat>(R0);
        Kalman K (priMu, priSigma, priPhi, priQ, priR)
        return Rcpp::wrap(priMu);'
fun <- cxxfunction(signature(mu0="numeric", Sigma0 = "numeric", Phi0 = "numeric",
                             Q0 = "numeric", R0 = "numeric"), 
                   body=src, plugin="RcppArmadillo", 
                   includes = '#include "/home/taylor/Desktop/kalmanclass.txt"')

fun(mu0, Sigma0, Phi, Q, R)
