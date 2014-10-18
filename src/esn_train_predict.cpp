
#include "utils.h"
#include "esn_train_predict.h"
#include <Eigen/Cholesky>

using namespace Rcpp;
using namespace RcppEigen;
using namespace Eigen;
using namespace std;
using Eigen::SelfAdjointView;


RcppExport SEXP train_esn(
  SEXP Yd,
  SEXP u,
  SEXP Win,
  SEXP W,
  SEXP Wback,
  SEXP tfRes,
  SEXP leakrate,
  SEXP lambda)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  using namespace std;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Rcpp::as;
    using Eigen::MappedSparseMatrix;
    using Eigen::SparseMatrix;
    typedef Eigen::MappedSparseMatrix<double> MSpMat;
    typedef Eigen::SparseMatrix<double> SpMat;
    typedef Map<MatrixXd> MapMatd;
    
    //assign constant variables
    const MatrixXd Y(as<MapMatd>(Yd));
    const MatrixXd UU(as<MapMatd>(u));
    const SpMat WWin(as<MSpMat>(Win));
    const SpMat WWback(as<MSpMat>(Wback));
    const SpMat WW(as<MSpMat>(W));
    const double lr(as<double>(leakrate));
    const double lam(as<double>(lambda));
    
    const int nrow_W = WW.rows();
    const int nobs = Y.rows();
    
    
    MatrixXd X(MatrixXd(1 + nrow_W, nobs - 1));
    VectorXd xx(nrow_W);
    MatrixXd WWout(MatrixXd(Y.cols(), nrow_W + 1));
    xx.fill(0);
    VectorXd y_add1(Y.cols() + 1);
    VectorXd uurow(UU.cols());


    for (int k = 0; k < nobs - 1; k++) {
      y_add1 << 1, Y.row(k);
      uurow = UU.row(k);

      xx = tanh((MatrixXd(WW * xx) + ( MatrixXd(WWin * uurow) + (MatrixXd(WWback * y_add1) )) )).array() * 
            (1 - lr) + xx.array() * lr;
      X.col(k) << 1, xx;
    }

    // + MatrixXd( (WWback * y_add1).addTo( WWin * UU.row(k)) ) 
    // + ( ( WWin * UU.row(k)).addTo(MatrixXd(WWback * y_add1) )) 

    //WWout = Y.bottomRows(nobs - 2).adjoint() * X.rightCols(X.cols() - 1).adjoint() * (XXt(X.rightCols(X.cols() - 1)).array() + 
    //        lam * MatrixXd::Identity(nrow_W + 1, nrow_W + 1).array()).matrix().llt().solve(MatrixXd::Identity(nrow_W + 1, nrow_W + 1));

    WWout = (XXt(X.rightCols(X.cols() - 1)).array() + 
            lam * MatrixXd::Identity(nrow_W + 1, nrow_W + 1).array()).matrix().llt().solve(X.rightCols(X.cols() - 1) * Y.bottomRows(nobs - 2)).adjoint();

    return (wrap(WWout));
    
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}


RcppExport SEXP predict_esn(
  SEXP ncolY,
  SEXP u,
  SEXP Win,
  SEXP W,
  SEXP Wback,
  SEXP Wout,
  SEXP tfRes,
  SEXP leakrate)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  using namespace std;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Rcpp::as;
    using Eigen::MappedSparseMatrix;
    using Eigen::SparseMatrix;
    typedef Eigen::MappedSparseMatrix<double> MSpMat;
    typedef Eigen::SparseMatrix<double> SpMat;
    typedef Map<MatrixXd> MapMatd;
    
    //assign constant variables
    const MatrixXd UU(as<MapMatd>(u));
    const MatrixXd WWout(as<MapMatd>(Wout));
    const SpMat WWin(as<MSpMat>(Win));
    const SpMat WWback(as<MSpMat>(Wback));    
    const SpMat WW(as<MSpMat>(W));
    const double lr(as<double>(leakrate));
    const int ncols(as<int>(ncolY));
    
    const int nrow_W = WW.rows();
    const int nobs = UU.rows();

    VectorXd xx(nrow_W);
    xx.fill(0);
    VectorXd xx1(nrow_W + 1);
    VectorXd y_add1(ncols + 1);
    VectorXd uurow(UU.cols());
    MatrixXd Yh(MatrixXd(UU.rows(), ncols));
    

    for (int k = 0; k < nobs - 1; k++) {
      y_add1 << 1, Yh.row(k);
      uurow = UU.row(k);

      xx = tanh((MatrixXd(WW * xx) + ( MatrixXd(WWin * uurow) + (MatrixXd(WWback * y_add1) )) )).array() * 
            (1 - lr) + xx.array() * lr;
      xx1 << 1, xx;
      Yh.row(k+1) = WWout * xx1;
    }


    return (wrap(Yh));
    
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}


