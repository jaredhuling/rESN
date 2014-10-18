#ifndef _RESN_OEMFIT_H
#define _RESN_OEMFIT_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector> 
#include <functional> 
#include <algorithm> 
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <float.h>


RcppExport SEXP train_esn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

RcppExport SEXP predict_esn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


#endif


