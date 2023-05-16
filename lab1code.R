library(ggplot2)


#1

s <- 22
n <- 70
f <- n-s
alpha0 <- 8
beta0 <- 8

#a
nDraws <- 10000
theta <- rbeta(nDraws,alpha0+s,beta0+f)
df <- data.frame(x=1:nDraws,y=theta, csmean=cumsum(theta)/1:10000)
csstd = vector(length=nDraws)
for(i in 1:nDraws){
  csstd[i] = sd(theta[1:i])
}
df["csstd"] = csstd
sstd_expected <- sqrt((alpha0 + s)*(beta0 + f)/((alpha0 + s + beta0 + f)^2 * (alpha0 + s + beta0 + f + 1)))
smean_expected <- (alpha0 + s)/(alpha0 + s + beta0 + f)
ggplot(df) +
  geom_histogram(aes(y),color = "white") +
  geom_vline(xintercept = smean_expected,color = "blue",size=1.5) +
  geom_vline(xintercept = smean_expected - sstd_expected,color = "red",size=1.1) +
  geom_vline(xintercept = smean_expected + sstd_expected,color = "red",size=1.1) +
  xlab("Theta")



ggplot(df[2:nDraws,]) +
  geom_line(aes(x=x,y=csstd)) +
  geom_hline(yintercept = sstd_expected,color = "red",size=0.4) +
  ylab("Sample Standard Deviation") +
  xlab("Number of Draws")


ggplot(df) +
  geom_line(aes(x=x,y=csmean)) +
  geom_hline(yintercept = smean_expected,color = "red",size=0.4) +
  ylab("Sample Mean") +
  xlab("Number of Draws")



#b

pr <- sum(theta>0.3)/length(theta)

pr_beta <- 1 - pbeta(0.3,alpha0+s,beta0+f)

#c

phi <- theta/(1-theta)
df["phi"] <- phi
ggplot(df) +
  geom_histogram(aes(phi, after_stat(density)),color = "white") +
  geom_density(aes(phi))



#2

#a
y <- c(33, 24, 48, 32, 55, 74, 23, 17)
mu <- 3.6
n <- 8
nDwaws <- 10000

tau_squared <- sum((log(y)-mu)^2)/n

posterior <- n*tau_squared / rchisq(nDwaws,n)

hist(posterior,breaks=30)



#b

gini <- 2*pnorm(sqrt(posterior/2),0,1)-1
hist(gini,breaks=30)



#c

interval <- quantile(gini,c(0.025,0.975))
df <- data.frame("y" = gini)
ggplot(df) +
  geom_histogram(aes(y),color = "white") +
  geom_vline(xintercept = interval[["2.5%"]],color = "red",size=1.5) +
  geom_vline(xintercept = interval[["97.5%"]],color = "red",size=1.1) +
  xlab("Gini")



#d
gini_density <- density(gini)
plot(gini_density)
gini_density <- data.frame("x" = gini_density$x, "y" = gini_density$y)
gini_density_sorted <- gini_density[order(gini_density$y, decreasing = T), ]
gini_density_sorted[["cumsum"]] <- cumsum(gini_density_sorted$y)

threshold <- 0.95 * sum(gini_density_sorted$y)
HDPI <- range(gini_density_sorted[gini_density_sorted[["cumsum"]] <= threshold, "x"])


ggplot(df) +
  geom_histogram(aes(y),color = "white") +
  geom_vline(xintercept = HDPI[1],color = "green",size=1) +
  geom_vline(xintercept = HDPI[2],color = "green",size=1) +
  xlab("Gini")


ggplot(df) +
  geom_histogram(aes(y),color = "white") +
  geom_vline(xintercept = interval[["2.5%"]],color = "red",size=1.1) +
  geom_vline(xintercept = interval[["97.5%"]],color = "red",size=1.1) +
  geom_vline(xintercept = HDPI[1],color = "green",size=1.1) +
  geom_vline(xintercept = HDPI[2],color = "green",size=1.1) +
  xlab("Gini")



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


