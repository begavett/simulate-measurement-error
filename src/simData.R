simData <- function(N = 100001, rel_x = .7225, true_beta = .6) {
  require(data.table)
  e <- function(N) qnorm(runif(N))
  a <- sqrt(rel_x) #0.85 # x' on x, a^2 is reliability of x
  b <- 0    # y on x'
  c <- true_beta  # y on x
  x <- e(N) 
  x_obs <- a*x + ((1-a^2)^.5)*e(N) # xprime (observed x)
  y <- c*x + b*x_obs + (1-c^2-b^2-2*a*b*c)^.5*e(N)
  df <- as.data.frame(t(rbind(x,x_obs,y)))
  return(as.data.table(df))
}