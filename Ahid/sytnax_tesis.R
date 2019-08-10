#Initial value
library (polynom)
library (orthopolynom)
library (gaussquad)
library (gmp)
library (Rmpfr)
library (MASS)
library (MixedPoisson)
data0=read.csv('./Ahid/data_tesis.csv',header=TRUE)[,-c(1:3)]
data=data.frame((data0))
y1=as.matrix((data[,1]))
y2=as.matrix((data[,2]))
estdelta1=est.delta(y1+5)
tau0=1/(estdelta1$ll.delta.max)^2
n=nrow(data)
x=as.matrix(cbind(rep(1,n),(data[,-c(1,2)])))
t1=lambda_start(y1,x)
beta10=t1$beta
t2=lambda_start(y2,x)
beta20=t2$beta
beta10
beta20
tau0
b1 <- as.matrix(beta10)
b2 <- as.matrix(beta20)
beta1 <- cbind(b1,b2)
t <- 1/(abs(estdelta1$ll.delta.max))
Y1 <- as.matrix(y1)
Y2 <- as.matrix(y2)
Y <- cbind(Y1,Y2)
n <- nrow(Y1)
My <- matrix(0, nrow = n, ncol = 1)
vnorm <- 1
#Myi
for(i in 1:n){
  
  sumxb <- 0
  sumfact <- 0
  sumfact2 <- 0
  
  for(j in 1:2){
    sumxb1 <- sumxb + (exp(x[i,] %*% beta1[,j]))
    sumxb <- sumxb1
  }
  
  for (m in 0:(Y[i,1]+Y[i,2])){
    kurung <- ((2/t)*sqrt(1+(2*t*sumxb)))^(-m)
    sumfact1 <- sumfact +
      ((factorial(Y[i,1]+Y[i,2]+m)/((factorial(Y[i,1]+Y[i,2]-
                                                 m))*factorial(m)))*kurung)
    sumfact <- sumfact1
  }
  
  for (m in 0:abs((Y[i,1]+Y[i,2]-1))){
    kurung <- (((2/t)*sqrt(1+(2*t*sumxb)))^(-m))
    sumfact2.1 <- sumfact2 + ((factorial(abs(Y[i,1]+Y[i,2]-
                                               1)+m)/(factorial(abs(Y[i,1]+Y[i,2]-1)-m)*factorial(m)))*kurung)
    sumfact2 <- sumfact2.1
  }
  My[i,1] <- (1/sqrt(1+(2*t*sumxb)))*(sumfact/sumfact2)
}
while(vnorm >= 10^-3){
  #Proses Matriks Gt
  #mendapatkan deltab
  
  deltab1 <- 0
  for (i in 1:n){
    
    sumxb <- 0
    sumfact <- 0
    atas <- 0
    sumfact2 <- 0
    bawah <- 0
    
    for(j in 1:2){
      sumxb1 <- sumxb + (exp(x[i,] %*% beta1[,j]))
      sumxb <- sumxb1
    }
    
    for (m in 0:(Y[i,1]+Y[i,2])){
      kurung <- (((2/t)*sqrt(1+(2*t*sumxb)))^(-m))
      sumfact1 <- sumfact +
        ((factorial(Y[i,1]+Y[i,2]+m)/((factorial(Y[i,1]+Y[i,2]-
                                                   m))*factorial(m)))*kurung)
      sumfact <- sumfact1
    }
    
    atas <- ((exp(x[i,] %*% beta1[,1]))*sumfact)
    
    for (m in 0:abs((Y[i,1]+Y[i,2]-1))){
      kurung <- (((2/t)*sqrt(1+(2*t*sumxb)))^(-m))
      sumfact2.1 <- sumfact2 + ((factorial(abs(Y[i,1]+Y[i,2]-
                                                 1)+m)/(factorial(abs(Y[i,1]+Y[i,2]-1)-m)*factorial(m)))*kurung)
      sumfact2 <- sumfact2.1
    }
    
    bawah <- ((sqrt(1+(2*t*sumxb)))*sumfact2)
    delta1 <- deltab1+((Y[i,1]-(atas/bawah))*x[i,])
    deltab1 <- delta1
  }
  
  #mendapatkan deltab2
  deltab2 <- 0
  for (i in 1:n){
    
    sumxb <- 0
    sumfact <- 0
    atas <- 0
    sumfact2 <- 0
    bawah <- 0
    for(j in 1:2){
      sumxb1 <- sumxb + (exp(x[i,] %*% beta1[,j]))
      sumxb <- sumxb1
    }
    
    for (m in 0:(Y[i,1]+Y[i,2])){
      kurung <- (((2/t)*sqrt(1+(2*t*sumxb)))^(-m))
      sumfact1 <- sumfact +
        ((factorial(Y[i,1]+Y[i,2]+m)/((factorial(Y[i,1]+Y[i,2]-
                                                   m))*factorial(m)))*kurung)
      sumfact <- sumfact1
    }
    
    atas <- ((exp(x[i,] %*% beta1[,2]))*sumfact)
    
    for (m in 0:abs((Y[i,1]+Y[i,2]-1))){
      kurung <- (((2/t)*sqrt(1+(2*t*sumxb)))^(-m))
      sumfact2.1 <- sumfact2 + ((factorial(abs(Y[i,1]+Y[i,2]-
                                                 1)+m)/(factorial(abs(Y[i,1]+Y[i,2]-1)-m)*factorial(m)))*kurung)
      sumfact2 <- sumfact2.1
    }
    
    bawah <- ((sqrt(1+(2*t*sumxb)))*sumfact2)
    
    lol <- ((Y[i,2]-(atas/bawah))*x[i,])
    delta2 <- deltab2+lol
    deltab2 <- delta2
  }
  #mendapatkan deltaT
  deltaT <- 0
  for (i in 1:n){
    
    sumxb <- 0
    sumfact <- 0
    atas <- 0
    sumfact2 <- 0
    bawah <- 0
    sumy <- 0
    
    for(j in 1:2){
      sumxb1 <- sumxb + (exp(x[i,] %*% beta1[,j]))
      sumxb <- sumxb1
    }
    
    for (m in 0:(Y[i,1]+Y[i,2])){
      kurung <- (((2/t)*sqrt(1+(2*t*sumxb)))^(-m))
      sumfact1 <- sumfact +
        ((factorial(Y[i,1]+Y[i,2]+m)/((factorial(Y[i,1]+Y[i,2]-
                                                   m))*factorial(m)))*kurung)
      sumfact <- sumfact1
    }
    
    atas <- ((1+(t*sumxb))*sumfact)
    
    for (m in 0:abs((Y[i,1]+Y[i,2]-1))){
      kurung <- (((2/t)*sqrt(1+(2*t*sumxb)))^(-m))
      sumfact2.1 <- sumfact2 + ((factorial(abs(Y[i,1]+Y[i,2]-
                                                 1)+m)/(factorial(abs(Y[i,1]+Y[i,2]-1)-m)*factorial(m)))*kurung)
      sumfact2 <- sumfact2.1
    }
    
    bawah <- (sqrt((1+(2*t*sumxb)))*sumfact2)
    
    for(j in 1:2){
      sumy1 <- sumy + (Y[i,j]-(atas/bawah))
      sumy <- sumy1
    }
    
    deltaTfix <- deltaT + ((1+(t*sumy)))
    deltaT <- deltaTfix
  }
  
  delta1T <- (-(1/(t^2)))*deltaT
  
  #Penggabungan menjadi Matriks Gt
  Gt <- matrix(c(deltab1,deltab2,delta1T), nrow=13, ncol=1)
  
  #Mendapatkan matriks Ht
  #Mendapatkan diff2b1b1
  diff2b1 <- 0
  for (i in 1:n){
    sumy <- 0
    summiu <- 0
    
    for (j in 1:2){
      sumy1 <- sumy + Y[i,j]
      sumy <- sumy1
    }
    
    for (j in 1:2){
      summiu1 <- summiu+(exp(x[i,] %*% beta1[,j]))
      summiu <- summiu1
    }
    
    atas <- (exp(x[i,] %*% beta1[,1]))*(1+(t*My[i,1]*(1+2*sumy)))
    bawah <- 1+(2*t*summiu)
    kurung1 <- My[i,1] - (atas/bawah) + ((exp(x[i,] %*%
                                                beta1[,1]))*(My[i,1]^2))
    diff2b1.1 <- diff2b1 + ((as.vector(exp(x[i,] %*% beta1[,1]))) *
                              ((x[i,]) %*% t(x[i,])))*(as.vector(kurung1))
    diff2b1 <- diff2b1.1
  }
  
  diff2b1 <- -(diff2b1)
  
  #Mendapatkan diff2b1b2
  diff2b1b2 <- 0
  for (i in 1:n){
    
    sumy <- 0
    summiu <- 0
    
    for (j in 1:2){
      sumy1 <- sumy+Y[i,j]
      sumy <- sumy1
    }
    
    for (j in 1:2){
      summiu1 <- summiu+(exp(x[i,] %*% beta1[,j]))
      summiu <- summiu1
    }
    atas <- (1+(t*My[i,1]))*(1+2*sumy)
    bawah <- 1+(2*t*summiu)
    kurung1 <- (atas/bawah)-(My[i,1]^2)
    
    lol1 <- (as.vector((exp(x[i,] %*% beta1[,1])) * (exp(x[i,] %*%
                                                           beta1[,2]))) * ((x[i,]) %*% t(x[i,])))
    diff2b1b2.1 <- diff2b1b2+ (lol1 * as.vector(kurung1))
    diff2b1b2 <- diff2b1b2.1
  }
  
  #Mendapatkan diff2b1t
  diff2b1tfix <- 0
  
  for (i in 1:n){
    sumy <- 0
    summiu <- 0
    for (j in 1:2){
      sumy1 <- sumy+Y[i,j]
      sumy <- sumy1
    }
    
    for (j in 1:2){
      summiu1 <- summiu+(exp(x[i,] %*% beta1[,j]))
      summiu <- summiu1
    }
    
    atas <- (1+t*summiu)*(1+(t*My[i,1])*(1+2*sumy))
    bawah <- (1+2*t*summiu)
    kurung1 <- (t*My[i,1])-(atas/bawah)+((My[i,1]^2)*(1+t*summiu))
    
    diff2b1t.1 <- diff2b1tfix+((as.vector((exp(x[i,] %*% beta1[,1])))
                                * t(x[i,]))*as.vector(kurung1))
    diff2b1tfix <- diff2b1t.1
  }
  
  diff2b1t <- as.vector(diff2b1tfix)*(1/(t^2))
  
  #mendapatkan diff2b2b2
  diff2b2 <- 0
  
  for (i in 1:n){
    sumy <- 0
    summiu <- 0
    
    for (j in 1:2){
      sumy1 <- sumy+Y[i,j]
      sumy <- sumy1
    }
    
    for (j in 1:2){
      summiu1 <- summiu+(exp(x[i,] %*% beta1[,j]))
      summiu <- summiu1
    }
    
    atas <- (exp(x[i,] %*% beta1[,2]))*(1+(t*My[i,1])*(1+2*sumy))
    bawah <- 1+(2*t*summiu)
    kurung1 <- My[i,1] - (atas/bawah) + ((exp(x[i,] %*% 
                                                beta1[,2]))*(My[i,1]^2))
    
    diff2b2.1 <- diff2b2+((as.vector(exp(x[i,] %*% beta1[,2]))) *
                            ((x[i,]) %*% t(x[i,]))* as.vector(kurung1))
    diff2b2 <- diff2b2.1
  }
  
  diff2b2 <- -(diff2b2)
  
  #mendapatkan diff2b2t
  diff2b2tfix <- 0
  
  for (i in 1:n){
    sumy <- 0
    summiu <- 0
    
    for (j in 1:2){
      sumy1 <- sumy+Y[i,j]
      sumy <- sumy1
    }
    
    for (j in 1:2){
      summiu1 <- summiu+(exp(x[i,] %*% beta1[,j]))
      summiu <- summiu1
    }
    
    atas <- (1+t*summiu)*(1+t*My[i,1]*(1+2*sumy))
    bawah <- 1+(2*t*summiu)
    kurung1 <- (t*My[i,1])-(atas/bawah)+(My[i,1]^2*(1+t*summiu))
    
    diff2b2t.1 <- diff2b2tfix+((as.vector((exp(x[i,] %*% beta1[,2])))
                                * t(x[i,]))*as.vector(kurung1))
    diff2b2tfix <- diff2b2t.1
  }
  
  diff2b2t <- as.vector(diff2b2tfix)*(1/(t^2))
  
  #mendapatkan diff2t2
  
  diff2t2 <- 0
  
  for (i in 1:n){
    sumy <- 0
    sumy2 <- 0
    summiu <- 0
    for (j in 1:2){
      sumy1 <- sumy+Y[i,j]
      sumy <- sumy1
    }
    
    for (j in 1:2){
      summiu1 <- summiu+(exp(x[i,] %*% beta1[,j]))
      summiu <- summiu1
    }
    
    pers1 <- -((1+t*summiu)/((t^4)*(1+2*t*summiu)))
    pers2 <- ((t^2)*My[i,1]*summiu)-((1+t*summiu)*(1-
                                                     (My[i,1]^2)*(1+2*t*summiu)+(2*t*My[i,1]*sumy)))
    
    
    diff2t2.1 <- diff2t2 + (pers1*pers2-
                              ((My[i,1]*(2+2*t*summiu))/(t^3))+(2/(t^3))+((My[i,1]*summiu+sumy)/(t^2
                              )))
    diff2t2 <- diff2t2.1
  }
  
  #Matriks Ht
  col1 <- matrix(c(diff2b1,diff2b1b2,diff2b1t), nrow=6, ncol=13)
  col2 <- matrix(c(diff2b1b2,diff2b2,diff2b2t), nrow=6, ncol=13)
  col3 <- matrix(c(diff2b1t,diff2b2t,diff2t2), nrow=1, ncol=13)
  Ht <- matrix(rbind(col1,col2,col3), nrow=13, ncol=13)
  teta0 <- rbind(b1,b2,t)
  teta1 <- teta0 - (solve(Ht) %*% Gt)
  
  b1 <- as.matrix(teta1[1:6,1])
  b2 <- as.matrix(teta1[7:12,1])
  t <- as.matrix(teta1[13,1])
  
  if(t > 2){
    t <- 1/(abs(t))
  }
  
  beta1 <- cbind(b1,b2)
  ji <- teta1-teta0
  print(teta1)
  vnorm <- norm(ji, type = 'f')
  print(vnorm)
} 
