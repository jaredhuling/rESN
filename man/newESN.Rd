\name{newESN}
\alias{newESN}
\title{Intialize an Echo State Network}
\usage{
newESN(Y, U, n.neurons = 50L, density = 0.5, back.density = 0.5,
  leak.rate = 0.3, lambda = 1e-07, spectral.radius = 1.15)
}
\arguments{
  \item{Y}{input matrix (can be multivariate) or vector of
  response. Each row is an observation in time, each column
  corresponds to a response}

  \item{U}{covariate matrix. Each row is an observation in
  time. Each column is a covariate}

  \item{n.neurons}{integer corresponding to the number of
  neurons in the reservoir (hidden layer)}

  \item{density}{double between 0 and 1. Corresponds to the
  density of connections in the reservoir}

  \item{back.density}{double between 0 and 1. Corresponds
  to the density of connections in the layer which connects
  the current response (in time) to the previous response.
  This is not always used in ESNs}

  \item{leak.rate}{double between 0 and 1. Corresponds to
  the rate for leaky integration}

  \item{lambda}{double corresponding to the tuning
  parameter for ridge regression. Higher values penalize
  more}

  \item{spectral.radius}{positive double corresponding to
  the spectral radius of the reservoirnewE}
}
\value{
An S4 object with class "esn"
}
\description{
Intialize an Echo State Network
}
\examples{
data(esn_data)

net <- newESN(example_data$Y.train.noisy,
              example_data$U.train,
              n.neurons = 1000,
              density = 0.02,
              back.density = 0.02,
              leak.rate = 0.6,
              lambda = 10)

net <- train(net)

ypred <- predict(net, u = example_data$U.test)

range <- 431:480
all <- c(example_data$Y.test.noisy[range,1], example_data$Y.true.test[range,1],
         ypred[range])
plot(example_data$Y.test.noisy[range,1], type = "l", lwd = 5, ylim = range(all), col = "green")
lines(ypred[range], col = "blue", lwd = 3)
lines(example_data$Y.test.signal[range,1], col = "black", lwd = 5)

legend("bottomleft", inset=0.01, legend=c("Noisy (observed) seq","ESN Prediction","True Signal"),
       col=c("green", "blue", "black"), lwd = 3, cex = 0.75)
}

