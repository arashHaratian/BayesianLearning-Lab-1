
# 1 -----------------------------------------------------------------------


s <- 22
n <- 70
f <- n - s
alpha0 <- beta0 <- 8


# 1.a
nDraws <- 10000

alpha_post <- alpha0+s
beta_post <- beta0+f

sample_post <- rbeta(nDraws, alpha_post, beta_post)

expected_mean <- alpha_post / (alpha_post + beta_post)
expected_var <- (alpha_post*beta_post) / ((alpha_post + beta_post)^2 * (alpha_post + beta_post + 1))

hist(sample_post, breaks = 30)
abline(v = expected_mean)
abline(v = expected_mean - sqrt(expected_var))
abline(v = expected_mean + sqrt(expected_var))


sample_means <- vector("numeric", nDraws)
sample_var <- vector("numeric", nDraws)
for(draw in 1:nDraws){
  samples <- rbeta(draw+1, alpha_post, beta_post)
  sample_means[draw] <- mean(samples)
  sample_var[draw] <- sd(samples)
}

plot(sample_means, type = "l")
abline(h = expected_mean, col = "blue")

plot(sample_var, type = "l")
abline(h = sqrt(expected_var), col = "blue")




# 1.b

exact_prob <- 1 - pbeta(0.3, alpha_post, beta_post)

estimated_prob <- sum(sample_post > 0.3) / nDraws

all.equal(exact_prob, estimated_prob)


# 1.c

odds <- sample_post / (1-sample_post)

hist(odds, freq = F)
lines(density(odds))

# 2 -----------------------------------------------------------------------

incomes <- c(33, 24, 48, 32, 55, 74, 23, 17)
mu <- 3.6
n <- length(incomes)
tau <- sum((log(incomes) - mu) ^ 2) / (n)
# tau <- var(log(incomes))


# 2.a
nDraws <- 10000

sample_post <- (tau * (n - 1)) / rchisq(nDraws, n)
hist(sample_post, breaks = 30)

# 2.b
gini <- 2 * pnorm(sqrt(sample_post)/ sqrt(2)) - 1

hist(gini, breaks = 30)

# 2.c

upper_CI <- quantile(gini, 0.975)
lower_CI <- quantile(gini, 0.025)

hist(gini, breaks = 30)
abline(v = upper_CI)
abline(v = lower_CI)

# 2.d

gini_density <- density(gini)
plot(gini_density)
gini_density <- data.frame("x" = gini_density$x, "y" = gini_density$y)
gini_density_sorted <- gini_density[order(gini_density$y, decreasing = T), ]
gini_density_sorted[["cumsum"]] <- cumsum(gini_density_sorted$y)

threshold <- 0.95 * sum(gini_density_sorted$y)
HDPI <- range(gini_density_sorted[gini_density_sorted[["cumsum"]] <= threshold, "x"])

hist(gini, breaks = 30)
abline(v = upper_CI, col = "blue")
abline(v = lower_CI, col = "blue")
abline(v = HDPI[1], col = "red")
abline(v = HDPI[2], col = "red")



# 3 -----------------------------------------------------------------------

wind_dir <- c(20, 314, 285, 40, 308, 314, 299, 296, 303, 326)
wind_dir_rad <- c(-2.79, 2.33, 1.83, -2.44, 2.23, 2.33, 2.07, 2.02, 2.14, 2.54)


lambda <- 0.5
mu <- 2.4
n <- length(wind_dir_rad)
# 3.a

posterior <- function(kapa, y = wind_dir_rad){
exp(kapa * (sum(cos(y - mu)) - lambda))/ (besselI(kapa, 0))^n 
}

kapas <- seq(0.05, 10, by = 0.05)

post_dist_y <- posterior(kapas)

post_dist_y <- post_dist_y/integrate(posterior, 0,10 )$value

plot(kapas, post_dist_y, type = "l")

# 3.b

mode_post <- kapas[which.max(post_dist_y)]
abline(v = mode_post)

