generate.mu <- function(DATA){
distDD <- dist(DATA)
distDD <- as.matrix(distDD)
mu <- numeric(nrow(distDD))
for(i in 1:nrow(distDD)){
  sD <- sort(distDD[i,])
  mu[i] <-sD[3]/sD[2]
}
return(mu)
}

mle.TwoNN <- function(X){
  n    <- length(X)
  logX <- log(X)
  n/sum(logX)
}

ID_linfit <- function(MU){
  N       <- length(MU)
  F_mui   <- (0:(N-1))/N
  Y       <-  -log(1-(F_mui))
  X       <- sort(log(MU))
  modlin1 <- lm(Y~X-1)
  c1      <- summary(modlin1)
  p1      <- qplot(X,Y)+theme_bw()+geom_abline(intercept = 0,slope = modlin1$coefficients,col=I("red"))+
    ylab("-log(1-(i/N))")+xlab(expression(log(mu)))
  return(list(coef=c1,plot=p1))
}



###########  compute Nq ############################
f_which.min <- function(vec, idx) sort(vec, index.return = TRUE)$ix[idx]

Nq_x2 <- function(q,dist_mat){
  n <- nrow(dist_mat)
  # negli q non Ã¨ contato i==j!
  # sort(c,index.return=T) !!!!!
  Nq <- t(apply(
    dist_mat, 1, function(x) { 
      index <- f_which.min(x,idx =2:(q+1) ) 
      l <- numeric(n)
      l[index] <- 1
      l
    }
  ) )
  return(Nq)
}



data.preprocessing <- function(DATA,q){
  
  mu      <- generate.mu(DATA)
  Nq      <- Nq_x2(q = q,dist_mat = as.matrix(dist(DATA)))

  return(list(mu=mu,Nq=Nq,q=q))  
  
}
