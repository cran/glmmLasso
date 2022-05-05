#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]] 

// [[Rcpp::export]]
SEXP  RcppEigenProd1(const Eigen::Map<Eigen::MatrixXd> Z_alles, 
                     const Eigen::Map<Eigen::SparseMatrix<double> > D, 
                     const Eigen::Map<Eigen::SparseMatrix<double> > SigmaInv,  
                     const Eigen::Map<Eigen::VectorXd> y, 
                     const Eigen::Map<Eigen::VectorXd> Mu){ 
  Eigen::MatrixXd ans = Z_alles.transpose() * D * SigmaInv * (y-Mu);
  return Rcpp::wrap(ans); 
}

// [[Rcpp::export]]
SEXP RcppEigenProd2(const Eigen::Map<Eigen::SparseMatrix<double> > D, 
                    const Eigen::Map<Eigen::SparseMatrix<double> > SigmaInv){ 
  Eigen::SparseMatrix<double> ans = D * (SigmaInv * D.transpose());
  return Rcpp::wrap(ans); 
}

// [[Rcpp::export]]
SEXP RcppEigenProd3(const Eigen::Map<Eigen::SparseMatrix<double> > W_inv_t, 
                    const Eigen::Map<Eigen::MatrixXd> Z_aktuell, 
                    const Eigen::Map<Eigen::MatrixXd> InvFisher2){ 
  Eigen::MatrixXd ans = W_inv_t * ( Z_aktuell * ( InvFisher2 * ( Z_aktuell.transpose() * W_inv_t.transpose())));
  return Rcpp::wrap(ans); 
}

// [[Rcpp::export]]
SEXP RcppEigenDiagSp(const Eigen::Map<Eigen::VectorXd> a){
  Eigen::SparseMatrix<double> ans = Eigen::SparseMatrix<double>(a.asDiagonal());
  return Rcpp::wrap(ans);
}

// [[Rcpp::export]]
SEXP RcppEigenSpChol(const Eigen::Map<Eigen::SparseMatrix<double> > W_opt){
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > cholesky(W_opt); // performs direct sparse LLT Cholesky factorization of A
  Eigen::SparseMatrix<double> ans = cholesky.matrixU();
  return Rcpp::wrap(ans); 
}

// [[Rcpp::export]]
SEXP RcppEigenSigmaInv(const Eigen::Map<Eigen::VectorXd> mu){
  Eigen::VectorXd oneVec;
  Eigen::MatrixXd sigma = mu * (oneVec.setOnes(mu.rows()) - mu).transpose();
  sigma.triangularView<Eigen::Lower>() = sigma.transpose().triangularView<Eigen::Lower>();
  int k=1;
  double det = sigma.determinant();
  
  while(det==0){
    
    Eigen::VectorXd kVec = Eigen::VectorXd::Constant(mu.rows(),pow(10,-7+k));
    Eigen::MatrixXd diagMatk = Eigen::SparseMatrix<double> (kVec.asDiagonal());
    sigma += diagMatk;
    det = sigma.determinant();
    k++;
    
    }
  Eigen::MatrixXd SigmaInv = sigma.inverse();
  return Rcpp::wrap(SigmaInv);
}

// [[Rcpp::export]]
SEXP RcppEigenInvMa(const Eigen::Map<Eigen::MatrixXd> Sigma){
  Eigen::MatrixXd SigmaInv = Sigma.inverse();
  return Rcpp::wrap(SigmaInv);
}

