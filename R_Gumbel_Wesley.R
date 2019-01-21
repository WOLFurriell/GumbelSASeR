#-----------------------------------------------------------
dir <- "D:/Estatística/ESTATÍSTICA_COMPUTACIONAL_II/MLE/Gumbel_Latex/"
setwd(dir)
# Wesley Furriel RA: 61493
#-----------------------------------------------------------
library(MASS)
library(fitdistrplus)
#install.packages("reliaR",dependencies = T)
#Biblioteca da  Gumbel
library(reliaR)
rm(list = ls())

#-----------------------------------------------------------
# Funcao densidade de probabilidade (d) fdp
gumbeld <- function(x,mu,sigma){
  aux <- ((x - mu)/sigma)
  exp(-aux -exp(-aux)) / sigma
}
gumbeld(x=2,mu=1.5,sigma=2)
dgumbel(x=2,mu=1.5,sigma=2)

# Funcao distribuicao de prob. (p) dist. acumulada
gumbelp <- function(q,mu,sigma){
  exp(-exp((-q + mu) / sigma))
}
gumbelp(q = 2, mu = 1.5, sigma = 2)
pgumbel(q = 2, mu = 1.5, sigma = 2, lower.tail = TRUE, log.p = FALSE)

# Funcao quantil (q)
gumbelq <- function(p,mu,sigma){
  mu - sigma*log(-log(p))
}
gumbelq(p = 0.5, mu = 1.5, sigma = 2)
qgumbel(p = 0.5, mu = 1.5, sigma = 2, lower.tail = TRUE, log.p = FALSE)

#Funcao para valores aleatorios
gumbelr <- function(n, mu, sigma){
  U <- runif(n)
  gumbelq(U,mu,sigma)
}
gumbelr(n = 1, mu = 2, sigma = 3)
rgumbel(n = 1, mu = 2, sigma = 3)

#Verificando o ajuste
#-----------------------------------------------------------
va <- gumbelr(n = 1000, mu = 2, sigma = 1.5)
jpeg("graph0.jpeg",height = 6, width = 8, units = 'in', res = 800)
  hist(va,probability = T)
     curve(gumbeld(x, mu=2, sigma=1.5), add=TRUE,col = "blue", lwd = 2)
dev.off()
#-----------------------------------------------------------
# Grafico da densidade
x   <- seq(-4, 20, length=1000)
hx  <- gumbeld(x = x, mu = 0, sigma = 1)
hx1 <- gumbeld(x = x, mu = 0.5, sigma = 2)
hx2 <- gumbeld(x = x, mu = 1.5, sigma = 3)
hx3 <- gumbeld(x = x, mu = 3, sigma = 4)
jpeg("graph1.jpeg",height = 6, width = 8, units = 'in', res = 800)
plot(x, hx, type = "l", lty = 2, lwd = 2)
colors <- c("black","blue", "darkgreen", "red")
labels <- c(expression(mu ~ "=0"~~~~sigma~"=1"),
            expression(mu ~ "=0.5"~sigma~"=2"),
            expression(mu ~ "=1.5"~sigma~"=3"),
            expression(mu ~ "=3"~~~~sigma~"=4"))
lines(x, hx1, col = colors[2], lwd = 2)
lines(x, hx2, col = colors[3], lwd = 2)
lines(x, hx3, col = colors[4], lwd = 2)
legend("topright", inset=.05,
       labels, lwd = 2, lty=c(2, 1, 1, 1), col = colors)
dev.off()

#-----------------------------------------------------------
# Grafico da acumulada 
hx  <- gumbelp(q = x, mu = 0, sigma = 1)
hx1 <- gumbelp(q = x, mu = 0.5, sigma = 2)
hx2 <- gumbelp(q = x, mu = 1.5, sigma = 3)
hx3 <- gumbelp(q = x, mu = 3, sigma = 4)
jpeg("graph2.jpeg", height = 6, width = 8, units = 'in', res = 800)
  plot(x, hx, type = "l", lty = 2, lwd = 2)
      lines(x, hx1, col = colors[2], lwd = 2)
      lines(x, hx2, col = colors[3], lwd = 2)
      lines(x, hx3, col = colors[4], lwd = 2)
legend("bottomright", inset=.05,
       labels, lwd = 2, lty=c(2, 1, 1, 1), col = colors)
dev.off()
#----------------------------------------------------------
y   <- gumbelr(n = 200, mu = 2, sigma = 3)
fit <- fitdist(y, "gumbel", start=list(mu = 2, sigma = 3))
jpeg("graph3.jpeg",height = 6, width = 8, units = 'in', res = 800)
plot(fit)
dev.off()

#-----------------------------------------------------------
#Estimacao por MLE
mu    <- 1.5
sigma <- 2
theta <- c(mu,sigma)

set.seed(666)
x <- gumbelr(n=1000,mu=theta[1],sigma=theta[2])
n<-length(x)
logLEgumbel <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
  aux<-(-x+mu)/sigma
  mle<--n * log(sigma) + sum(log(exp(aux - exp(aux))))
  return(mle)
}
library(maxLik)
maxLik(logLik = logLEgumbel, start = c(1.3,2))

#-----------------------------------------------------------
#Metodo dos momentos
library(reliaR)
x<-rgumbel(n=500,mu=1.5,sigma=2)
MMgumbel <- function(theta,x) {
  mu     <- theta[1]
  sigma  <- theta[2]
  xbar   <- mu + 0.5772*sigma
  s      <- (pi^2/6)*sigma^2
  m1     <- xbar -x
  m2     <- s - (x - xbar)^2 
  f      <- cbind(m1,m2)
  return(f)
}
library(gmm)
  mm<-gmm(MMgumbel,x,c(mu=1.5,sigma=2));mm
  vcov(mm)
    muest<-mm$coefficients[1]
    sigmaest<-mm$coefficients[2]

#-----------------------------------------------------------
#QQplot
library(fastR)
qgumbel.plot <- function(x){qgumbel(x,muest,sigmaest)}
jpeg("momento_plot.jpeg",height = 6, width = 8, units = 'in',res=800)
xqqmath(x, distribution = qgumbel.plot, 
        xlab = "Quantil gumbel",
        ylab = "Quantil amostral")
dev.off()
#-----------------------------------------------------------
#Score
Uscore <- function(n, x, theta){
  mu <- theta[1]
  sigma <- theta[2]
  W <- -((x - mu)/sigma)
  Es <- matrix(nrow=2)  
  Es[1] <- n / sigma - sum(1 / sigma * exp(W))
  Es[2] <- (sum((mu-x) * exp(Z) + x) - (mu + sigma) * n) / sigma ^ 2
  return(Es)
}
#-----------------------------------------------------------
# Matriz Hessiana 
Hessian <- function(n, x, theta){
  mu <- theta[1]
  sigma <- theta[2]
  W <- -((x - mu)/sigma)
  Hs <- matrix(ncol=2,nrow=2)
  Hs[1,1] <- -1 / sigma ^ 2 * sum(exp(W))  
  Hs[1,2] <- (-n * sigma + sum(exp(W) * (sigma - x + mu))) / sigma ^ 3
  Hs[2,1] <- Hs[1,2]
  Hs[2,2] <- (2 * n * mu * sigma+n*sigma^2
              -sum((-x+mu)*(mu+2*sigma-x)*exp(W) + 2 * x * sigma)) / sigma ^ 4
  return(Hs)}
#-----------------------------------------------------------
#Iniciando a simulacao
par1<-1.5 #mu
par2<-2 #sigma
nrep<-1000
m<-seq(10,100,10)
results<-matrix(nrow=length(m),ncol=3)
mle<-matrix(ncol=2,nrow = nrep)
k<-1
for (n in m){
  set.seed(666)
  for (i in 1:nrep){
    x<- rgumbel(n=n,mu=par1,sigma=par2)
    mle[i,]<-fitdist(x,"gumbel",start=list(mu=par1,sigma=par2))$estimate}
  results[k,1]<-n
  results[k,2]<-mean(mle[,1]) - par1
  results[k,3]<-mean(mle[,2]) - par2
  k<-k+1
  cat(k,i,n,"\n")}
jpeg("graph4.jpeg",height = 6, width = 8, units = 'in',res=800)
par(mfrow = c(1,2)) 
plot(results[,1],results[,2],type="l",col="red",lwd=2)
plot(results[,1],results[,3],type="l",col="blue",lwd=2)
dev.off()

