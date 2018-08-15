setwd('/home/joao_polo/Dropbox/Mestrado/Dissertacao/fbst-nparam')
source('funcoes-uteis.R') 
library(nortest)

#FBST 2 POP calculado por diferença de normais

#funcao que recebe como parametros os dados (pop x e y) ja tranformados pela funcao 
#phi barra e calcula a estatistica de teste do FBST para duas pop
testeFBST2 <- function(x, n, m, tau){
  sum( (n*m*x)^2 / (n*m+(n+m)*tau)   )
}

###################################################################################
############ ESSE MÉTODO NAO DÁ CERTO POR CAUSA DA SUPOSIÇÃO DO BRUNK  ############
###################################################################################
#funcao que calcula o evalor modificado
#dados - dados
#phi -base ortonormal
#r - tamanho da sequencia da base ortonormal
#peso - peso para cada elemento da sequencia
#simu - numero de observacoes que devem ser geradas da distribuicao de comparacao
ev_mod2 <- function(dados1, dados2, phi, r, peso, simu=10000){
  n <- length(dados1)
  m <- length(dados2)
  tau <- peso(1:r)
  x <- phi_barra(dados1, phi, r)
  y <- phi_barra(dados2, phi, r)
  w<-x-y
  aux <- testeFBST2(w, n, m, tau)
  #print(aux)

  q <- rep(NA, simu)
  for (i in 1:simu){
    er<-rnorm(r, mean=0, sd=1/sqrt((n+m)/(n*m+(n+m)*tau)))
    q[i] <- sum( (n*m*er) / sqrt((n+m)/(n*m+(n+m)*tau) )^2  )
  }
  
  1 - mean(q < aux)
}

##TESTE

todos<-rep(NA,100)
for (i in 1:100){
  data1<-rnorm(100)
  data2<-rnorm(100)
  todos[i]<-ev_mod2(data1, data2, phi, r =50, function(x) x^2)
  print(i)
}

hist(todos)
mean(todos<0.05)

#############################################################################

#############################################################################
### TESTANDO ENCONTRAR A DIST ATRAVES DE  BOOTSTRAP #########################
#############################################################################

#FBST 2 POP calculado por diferença de normais
#funcao que calcula o evalor modificado
#dados - dados
#phi -base ortonormal
#r - tamanho da sequencia da base ortonormal
#peso - peso para cada elemento da sequencia
#simu - numero de observacoes que devem ser geradas da distribuicao de comparacao
ev_mod2.2 <- function(dados1, dados2, phi, r, peso, simu=1000){
  n <- length(dados1)
  m <- length(dados2)
  tau <- peso(1:r)
  x <- phi_barra(dados1, phi, r)
  y <- phi_barra(dados2, phi, r)
  w<-x-y
  aux <- testeFBST2(x, n, m, tau)
  #print(aux)
  
  auxq <- rep(NA, simu)
  junto <- c(dados1,dados2)
  
  for (i in 1:simu){
    #amostra<-sample(junto, n+n)
    #amostra1 <- amostra[1:n]
    #amostra2 <- amostra[(n+1):(n+n)]
    amostra1 <- sample(junto, n, replace = T)
    amostra2 <- sample(junto, m, replace = T)
    x1 <- phi_barra(amostra1, phi, r)
    y1 <- phi_barra(amostra2, phi, r)
    w1<-x1-y1
    auxq[i] <- testeFBST2(w1, n, m, tau)
  }
  
  #hist(auxq)
  quantil <- sort(auxq)[round(simu*0.05)]
  
  #abline(v=1,col='red',lwd=2)
  #print(aux)
  aux < quantil
}

#usar essa estat de teste somente quando os tamanhos de amostra são iguais
testeIZBICKI<-function(x, n, tau){
  sum( ( (n^2)*(x^2) ) / 2*(n+tau^2) )
}

ev_mod2iz <- function(dados1, dados2, phi, r, peso, simu=1000){
  n <- length(dados1)
  m <- length(dados2)
  tau <- peso(1:r)
  x <- phi_barra(dados1, phi, r)
  y <- phi_barra(dados2, phi, r)
  w<-x-y
  aux <- testeIZBICKI(x, n, tau)
  #print(aux)
  
  auxq <- rep(NA, simu)
  junto <- c(dados1,dados2)
  
  for (i in 1:simu){
    #amostra<-sample(junto, n+n)
    #amostra1 <- amostra[1:n]
    #amostra2 <- amostra[(n+1):(n+n)]
    amostra1 <- sample(junto, n, replace = T)
    amostra2 <- sample(junto, m, replace = T)
    x1 <- phi_barra(amostra1, phi, r)
    y1 <- phi_barra(amostra2, phi, r)
    w1<-x1-y1
    auxq[i] <- testeIZBICKI(w1, n, tau)
  }
  
  #hist(auxq)
  quantil <- sort(auxq)[round(simu*0.05)]
  
  #abline(v=1,col='red',lwd=2)
  #print(aux)
  aux < quantil
}


#TESTAR PARA DIFERENTES DISTRIBUIÇÕES

medias<-seq(-2,2,by=0.2)
resultado<-rep(NA, length(medias))
resultadoks<-rep(NA, length(medias))
for (k in 1:length(medias)){
  testes.save<-rep(NA, 500)
  testes.ks<-rep(NA, 500)
  for (j in 1:500){
    data1<-rnorm(100, mean=medias[k])
    data2<-rnorm(100, mean=medias[k])
    testes.save[j]<-ev_mod2.2(data1, data2, phi, r =50, function(x) x^2, simu=100)
    testes.ks[j]<-ks.test(data1,data2)$'p.value'
    #print(j)
  }
  resultado[k] <- mean(testes.save)
  resultadoks[k] <- mean(testes.ks<0.05)
  #print(k)
}

plot(medias,resultado, type='l')
lines(medias, resultadoks, col='red')

#distribuicao exponecial
medias<-seq(1,5,by=0.2)
resultado<-rep(NA, length(medias))
resultadoks<-rep(NA, length(medias))
for (k in 1:length(medias)){
  testes.save<-rep(NA, 500)
  testes.ks<-rep(NA, 500)
  for (j in 1:500){
    data1<-rexp(100, rate=2)
    data2<-rexp(100, rate=medias[k])
    testes.save[j]<-ev_mod2.2(data1, data2, phi, r =50, function(x) x^2, simu=100)
    testes.ks[j]<-ks.test(data1,data2)$'p.value'
    #print(j)
  }
  resultado[k] <- mean(testes.save)
  resultadoks[k] <- mean(testes.ks<0.05)
  print(k)
}

plot(medias,resultado, type='l')
lines(medias, resultadoks, col='red')


#distribuicao multimodal (trigonometrica)
densidade<-function(x,k){
  abs(sin(pi*k*x)*pi/2)
}
gerar<-function(n,k){
  gerado<-rep(NA,n)
  i<-1
  while(is.na(gerado[n])){
    aux1<- runif(1)
    aux2<- runif(1,0,pi/2)
    
    if(densidade(aux1,k) < aux2){
      gerado[i]<-aux1
      i<-i+1
    }
  }
  gerado
}


medias<-seq(1,10,by=1)
resultado<-rep(NA, length(medias))
resultadoks<-rep(NA, length(medias))
for (k in 1:length(medias)){
  testes.save<-rep(NA, 500)
  testes.ks<-rep(NA, 500)
  for (j in 1:500){
    data1<-gerar(100, 5)
    data2<-gerar(100, medias[k])
    testes.save[j]<-ev_mod2.2(data1, data2, phi, r =50, function(x) x^2, simu=100)
    testes.ks[j]<-ks.test(data1,data2)$'p.value'
    #print(j)
  }
  resultado[k] <- mean(testes.save)
  resultadoks[k] <- mean(testes.ks<0.05)
  print(k)
}

plot(medias,resultado, type='l')
lines(medias, resultadoks, col='red')

#izbicki

medias<-seq(1,10,by=1)
resultado<-rep(NA, length(medias))
resultadoks<-rep(NA, length(medias))
for (k in 1:length(medias)){
  testes.save<-rep(NA, 500)
  testes.ks<-rep(NA, 500)
  for (j in 1:500){
    data1<-gerar(100, 5)
    data2<-gerar(100, medias[k])
    testes.save[j]<-ev_mod2iz(data1, data2, phi, r =50, function(x) x^2, simu=100)
    testes.ks[j]<-ks.test(data1,data2)$'p.value'
    #print(j)
  }
  resultado[k] <- mean(testes.save)
  resultadoks[k] <- mean(testes.ks<0.05)
  print(k)
}

plot(medias,resultado, type='l')
lines(medias, resultadoks, col='red')