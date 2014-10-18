
propackSVD <- function (X, neig = min(m, n), opts = list()) {
  if (is.matrix(X)) {
    m <- dim(X)[1]
    n <- dim(X)[2]
    storage.mode(X) <- "double"
  }
  else if (is.extmat(X)) {
    m <- extmat.nrow(X)
    n <- extmat.ncol(X)
  }
  else {
    stop("unsupported matrix type for SVD")
  }
  storage.mode(neig) <- "integer"
  storage.mode(opts) <- "list"
  .Call("propack_svd", X, neig, opts, PACKAGE = "svd")
}

check_esn <- function(object) {
  errors <- character()
  length_lr <- length(object@leak.rate)
  if (length_lr != 1) {
    msg <- paste("leak.rate is length ", length_lr, ".  Should be 1", sep = "")
    errors <- c(errors, msg)
  }
  
  length_lr <- length(object@spectral.radius)
  if (length_lr != 1) {
    msg <- paste("spectral.radius is length ", length_lr, ".  Should be 1", sep = "")
    errors <- c(errors, msg)
  }
  
  if (object@spectral.radius < 0) {
    msg <- paste("spectral.radius is negative. It should be positive", sep = "")
    errors <- c(errors, msg)
  }
  
  if (object@lambda < 0) {
    msg <- paste("lambda is negative. It should be positive", sep = "")
    errors <- c(errors, msg)
  }
  
  
  
  #length_name <- length(object@name)
  #if (length_name != 1) {
  #  msg <- paste("Name is length ", length_name, ".  Should be 1", sep = "")
  #  errors <- c(errors, msg)
  #}
  
  if (length(errors) == 0) TRUE else errors
}

setClassUnion("all.mat", c("matrix", "sparseMatrix", "dgCMatrix"))

setClass("esn", representation(tfRes = "function", 
                               tfReadout = "function", 
                               leak.rate = "numeric",
                               lambda = "numeric",
                               spectral.radius = "numeric",
                               u = "matrix",
                               Win = "all.mat",
                               W = "all.mat",
                               Wback = "all.mat",
                               Wout = "all.mat",
                               Yd = "matrix"), 
         prototype(Wout = matrix (0)),
         validity = check_esn)

initweights <- function(n1, n2, d, enforce_esp = TRUE, spectral.radius = 1.15) {
  require(Matrix)
  mask <- matrix(runif(n1 * n2), ncol = n2) < d
  W <- (matrix(runif(n1 * n2), ncol = n2) - 0.5) * mask
  if (n1 == n2 & enforce_esp) {
    rhoW <- abs(propackSVD(W, neig = 1L, opts = list(tol = 1e-7, maxiter = 1e4))$d)
    p <- spectral.radius / rhoW
    W <- p * W
  }
  as(W, "sparseMatrix")
}

#' Intialize an Echo State Network
#'
#' @param Y input matrix (can be multivariate) or vector of response. Each row is an observation in time, each column corresponds to a response
#' @param U covariate matrix. Each row is an observation in time. Each column is a covariate
#' @param n.neurons integer corresponding to the number of neurons in the reservoir (hidden layer)
#' @param density double between 0 and 1. Corresponds to the density of connections in the reservoir
#' @param back.density double between 0 and 1. Corresponds to the density of connections in the layer which connects the current 
#' response (in time) to the previous response. This is not always used in ESNs
#' @param leak.rate double between 0 and 1. Corresponds to the rate for leaky integration
#' @param lambda double corresponding to the tuning parameter for ridge regression. Higher values penalize more
#' @param spectral.radius positive double corresponding to the spectral radius of the reservoir
#' @return An S4 object with class "esn" 
#' @export
#' @examples
#' data(esn_data)
#' 
#' net <- newESN(example_data$Y.train.noisy, 
#'               example_data$U.train, 
#'               n.neurons = 1000, 
#'               density = 0.02, 
#'               back.density = 0.02, 
#'               leak.rate = 0.6, 
#'               lambda = 10)
#' 
#' net <- train(net)
#' 
#' ypred <- predict(net, u = example_data$U.test)
#' 
#' range <- 431:480
#' all <- c(example_data$Y.test.noisy[range,1], example_data$Y.true.test[range,1],
#'          ypred[range])
#' plot(example_data$Y.test.noisy[range,1], type = "l", lwd = 5, ylim = range(all), col = "green")
#' lines(ypred[range], col = "blue", lwd = 3)
#' lines(example_data$Y.test.signal[range,1], col = "black", lwd = 5)
#' 
#' legend("bottomleft", inset=0.01, legend=c("Noisy (observed) seq","ESN Prediction","True Signal"), 
#'        col=c("green", "blue", "black"), lwd = 3, cex = 0.75)
#' 
#' 
newESN <- function(Y, U, n.neurons = 50L, density = 0.5, 
                   back.density = 0.5,
                   leak.rate = 0.3, lambda = 1e-7,
                   spectral.radius = 1.15) {
  n.neurons <- as.integer(n.neurons)
  stopifnot(density >= 0.01 & density <= 1)
  net <- new("esn", tfRes = tanh, tfReadout = tanh, 
             leak.rate = leak.rate,
             lambda = lambda,
             u = U,
             Win = initweights(n.neurons, ncol(U), d = 1, enforce_esp = FALSE),
             W = initweights(n.neurons, n.neurons, density, spectral.radius = spectral.radius),
             Wback = initweights(n.neurons, 1L + ncol(Y), back.density, enforce_esp = FALSE),
             Yd = as.matrix(Y) )
}


setGeneric("trainR", function(net) 0)
setMethod("trainR", signature(net = "esn"), function(net) {
  X <- array(0, dim = c(1 + nrow(net@W), (nrow(net@Yd) - 1) ) )
  x <- array(0, dim = c(nrow(net@W),1) )
  for (k in 1:(nrow(net@Yd) - 1)) {
    x <- net@tfRes( net@W %*% x + net@Wback %*% c(1,net@Yd[k,]) + net@Win %*% net@u[k,] ) * 
      (1 - net@leak.rate) + net@leak.rate * x
    X[,k] <- c(1, as.numeric(x))
  }
  M <- X[,2:ncol(X)]
  #TT <- net@tfReadout(net@Yd[3:nrow(net@Yd),])
  TT <- net@Yd[3:nrow(net@Yd),]
  net@Wout <- t(TT) %*% t(M) %*% solve(tcrossprod(M) + net@lambda * diag(1 + nrow(net@W)))
  net
})

setGeneric("train", function(net) 0)
setMethod("train", signature(net = "esn"), function(net) {
  net@Wout <- .Call("train_esn", Yd = net@Yd, u = net@u, Win = net@Win, 
                    W = net@W, Wback = net@Wback, tfRes = "nothing",
                    leakrate = net@leak.rate, lambda = net@lambda, package = "rESN")
  net
})

setGeneric("predictR", function(net, u) 0)
setMethod("predictR", signature(net = "esn", u = "matrix"), function(net, u) {
  Yh <- array(0, dim = c((nrow(u)), (ncol(net@Yd)) ) )
  x <- array(0, dim = c(nrow(net@W),1) )
  for (k in 1:(nrow(u) - 1)) {
    x <- net@tfRes( net@W %*% x + net@Wback %*% c(1, Yh[k,]) + net@Win %*% u[k,] ) * 
      (1 - net@leak.rate) + net@leak.rate * x
    Yh[k+1,] <- net@Wout %*% c(1,as.numeric(x)) 
  }
  Yh
})

setGeneric("predict", function(net, u) 0)
setMethod("predict", signature(net = "esn", u = "matrix"), function(net, u) {
  Yh <- .Call("predict_esn", ncolY = ncol(net@Yd), u = u, Win = net@Win, 
              W = net@W, Wback = net@Wback, Wout = net@Wout, tfRes = "nothing",
              leakrate = net@leak.rate, lambda = net@lambda, package = "rESN")
  Yh
})


setMethod("predict", signature(net = "esn", u = "missing"), function(net) {
  Yh <- .Call("predict_esn", ncolY = ncol(net@Yd), u = net@u, Win = net@Win, 
              W = net@W, Wback = net@Wback, Wout = net@Wout, tfRes = "nothing",
              leakrate = net@leak.rate, lambda = net@lambda, package = "rESN")
  Yh
})



setMethod("predictR", signature(net = "esn", u = "missing"), function(net) {
  Yh <- array(0, dim = c((nrow(net@u)), (ncol(net@Yd)) ) )
  x <- array(0, dim = c(nrow(net@W),1) )
  for (k in 1:(nrow(net@Yd) - 1)) {
    x <- net@tfRes( net@W %*% x + net@Wback %*% c(1, Yh[k,]) + net@Win %*% net@u[k,] ) * 
      (1 - net@leak.rate) + net@leak.rate * x
    Yh[k+1,] <- net@Wout %*% c(1,as.numeric(x))
  }
  Yh
})
