rESN
====

Echo state networks in R


## Installation

**rESN** is not available on CRAN, but can be installed using the R package **devtools**. **rESN** can be installed with the following R code:

```r
devtools::install_github("jaredhuling/rESN")
library(rESN)
```

## Example
U.train is the covariate matrix. Each row is an observation in time. Each column is a covariate.
Y.train.noisy is a response in that depends nonlinearly on U over time.
n.neurons defines the number of neurons in the "reservoir" (hidden layer).
density is the density of non-zero connections in the reservoir.
lambda is the ridge regression tuning parameter
```r
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


```