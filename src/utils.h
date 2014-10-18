#ifndef _RESN_UTILS_H
#define _RESN_UTILS_H


#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector> 
#include <functional> 
#include <algorithm> 
#include <iostream>
#include <cmath>


using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::Lower;
using Eigen::Ref;

   
double threshold(double num);

//computes X'WX where W is diagonal (input w as vector)
MatrixXd XtWX(const MatrixXd& xx, const MatrixXd& ww);

//computes X'X
MatrixXd XtX(const MatrixXd& xx);

//computes XX'
MatrixXd XXt(const MatrixXd& xx);

//computes hyperbolic tangent
MatrixXd tanh(const MatrixXd& x);

#endif