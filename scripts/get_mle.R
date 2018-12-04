libSize <- load("../MMUPHin/results/quick_and_easy/qe_biom.RData") %>%
  get %>%
  otu_table() %>%
  apply(2, sum)

estimate_expNormal <- function(x) {
  minLik <- optim(par = c(0, 1),
        fn = function(para) {
          sum(-dnorm(log(x), mean = para[1], sd = para[2], log = TRUE))
        })
  return(list(mu = minLik$par[1], sigma = minLik$par[2]))
}

estimate_exp <- function(x) {
  minLik <- optimize(f = function(para) {
    sum(-dexp(x, rate = para, log = TRUE))
  },
  interval = c(0, 10000)
  )
  return(list(rate = minLik$minimum))
}

estimate_gamma <- function(x) {
  minLik <- optim(par = c(1, 1),
                  fn = function(para) {
                    sum(-dgamma(x, shape = para[1], rate = para[2], log = TRUE))
                  })
  return(list(shape = minLik$par[1], rate = minLik$par[2]))
}

estimate_logGamma <- function(x) {
  minLik <- optim(par = c(1, 1),
                  fn = function(para) {
                    sum(-dgamma(log10(x), shape = para[1], rate = para[2], log = TRUE))
                  })
  return(list(shape = minLik$par[1], rate = minLik$par[2]))
}

para_expNormal <- estimate_expNormal(libSize)
para_exp <- list(rate = 1 / mean(libSize))
para_gamma <- estimate_gamma(libSize)
para_logGamma <- estimate_logGamma(libSize)
ggplot(data = tibble(`log10 Lib Size` = log10(libSize)),
       aes(x = `log10 Lib Size`)) +
  geom_density() +
  xlim(0, 8) +
  stat_function(fun = function(x) {ifelse(x < para_expNormal$mu + 3*para_expNormal$sigma,
                                          dnorm(x,
                                                mean = para_expNormal$mu,
                                                sd = para_expNormal$sigma) / pnorm(3),
                                          0)},
                aes(color = "log normal")) +
  stat_function(fun = function(x) {dexp(10^x, rate = para_exp$rate) * 10^x * log(10)},
                aes(color = "exponential")) +
  stat_function(fun = function(x) {dgamma(10^x,
                                          shape = para_gamma$shape,
                                          rate = para_gamma$rate) * 10^x * log(10)},
                aes(color = "gamma")) +
  stat_function(fun = function(x) {dgamma(x,
                                          shape = para_logGamma$shape,
                                          rate = para_logGamma$rate)},
                aes(color = "log gamma")) +
  theme_bw() +
  scale_color_manual("Distribution",
                     values = c("log normal" = "red",
                                "exponential" = "blue",
                                "gamma" = "orange",
                                "log gamma" = "green"))

