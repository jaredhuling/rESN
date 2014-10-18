

#include "utils.h"


//computes X'WX where W is diagonal (input w as vector)
MatrixXd XtWX(const MatrixXd& xx, const MatrixXd& ww) {
  const int n(xx.cols());
  MatrixXd AtWA(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(xx.adjoint() * ww.asDiagonal()));
  return (AtWA);
}

//computes X'X where 
MatrixXd XtX(const MatrixXd& xx) {
  const int n(xx.cols());
  MatrixXd AtA(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(xx.adjoint()));
  return (AtA);
}

//computes XX' where 
MatrixXd XXt(const MatrixXd& xx) {
  const int n(xx.rows());
  MatrixXd AAt(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(xx));
  return (AAt);
}


MatrixXd tanh(const MatrixXd& x) {
  const int n(x.rows());
  const int p(x.cols());
  MatrixXd expmat(MatrixXd(n, p));
  expmat = (2 * x.array()).array().exp();
  expmat = (expmat.array() - 1) / (expmat.array() + 1);
  return (expmat);
}



