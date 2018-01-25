setwd('/home/joao_polo/Dropbox/Mestrado/Dissertacao/fbst-nparam')
source('funcoes-uteis.R') 

#carregando pacotes de devem ser usados
library(nortest)
library(tidyverse)

#funcao que calcula o evalor não modificado
#dados - dados
#phi -base ortonormal
#r - tamanho da sequencia da base ortonormal
#peso - peso para cada elemento da sequencia
ev_nmod <- function(dados, phi, r, peso){
  n <- length(dados)
  tau <- peso^(1:r)
  x <- phi_barra(dados, phi, r)
  aux <- testeFBST(x, n, tau)

  1 - pchisq(aux, df=r)
}

#funcao que calcula o evalor modificado
#dados - dados
#phi -base ortonormal
#r - tamanho da sequencia da base ortonormal
#peso - peso para cada elemento da sequencia
#simu - numero de observacoes que devem ser geradas da distribuicao de comparacao
ev_mod <- function(dados, phi, r, peso, simu=1000){
  n <- length(dados)
  tau <- peso^(1:r)
  x <- phi_barra(dados, phi, r)
  aux <- testeFBST(x, n, tau)
  
  q <- rep(NA, simu)
  for (i in 1:simu){
    q[i] <- sum(((n*rnorm(r,sd=1/sqrt(n)))/sqrt(n+tau))^2)
  }
  
  1 - mean(q < aux)
}

##################################################################################
##### ESSE BLOCO COMPARA A DISTRIBUICAO DO EVALOR MODIFICADO E NAO MODIFICADO ####
##################################################################################

#H_0 verdadeira
evm0 <- rep(NA,1000)
evnm0 <- rep(NA,1000)
for (i in 1:1000){
  dados <- runif(10000)
  evm0[i] <- ev_mod(dados, phi, r = 50, peso = 5)
  evnm0[i] <- ev_nmod(dados, phi, r = 50, peso = 5)
}
par(mfrow=c(1,2))
hist(evm0)
hist(evnm0)

#H_0 falsa
evm1 <- rep(NA,1000)
evnm1 <- rep(NA,1000)
for (i in 1:1000){
  dados <- rnorm(10000, mean=5, sd=1)
  evm1[i] <- ev_mod(dados, phi, 50, 5)
  evnm1[i] <- ev_nmod(dados, phi, 50, 5)
}
hist(evm1)
hist(evnm1)

#################################################################################

#gerar dados e aplicar uma F nos dados e 
#verificar se os dados modificado tem distribuicao uniforme

#funcao que vai calcular o evalor nos dados transformados (padronizados)
#dados - dados
#phi -base ortonormal
#r - tamanho da sequencia da base ortonormal
#peso - peso para cada elemento da sequencia
#simu - numero de observacoes que devem ser geradas da distribuicao de comparacao
#pf - funcao F que deve ser usada para padronizar os dados
ev_mod_pf<-function(dados, phi, r = 50, peso = 5, simu=1000, pf){
  ev <- dados %>%  pf() %>% ev_mod(phi, r, peso, simu = simu)
  ev
}

#funcao que vai calcular o pvalor do ks nos dados transformados
#dados - dados
#pf - funcao F que deve ser usada para padronizar os dados
pvalor_ks <- function(dados, pf){
  test <- dados %>% pf() %>% ks.test(punif)
  test$'p.value'
}

#essa funcao retorna uma matriz com a primeira coluna de evalore e a outra de pvalor ks
evalor_vs_ks <- function(matriz_dados, phi, r = 50, peso = 5, pf){
  ev <- apply(matriz_dados, 2, ev_mod_pf, phi = phi, r = r, peso = peso, pf = pf)
  pv <- apply(matriz_dados, 2, pvalor_ks, pf = pf)
  cbind(ev, pv)
}

#funcao para ver a proporcao de rejeicao fixado um alpha
poder_ev_vs_ks <- function(matriz_dados, phi, r = 50, peso = 5, pf, alpha = 0.05){
  poder <- colMeans(evalor_vs_ks(dados, phi, 50, 5 , pf) < alpha)
  poder
}

#mesma funcao usando pipe
poder_ev_vs_ks1 <- function(matriz_dados, phi, r = 50, peso = 5, pf, alpha = 0.05){
  poder <- evalor_vs_ks(dados, phi, 50, 5 , pf) %>% 
          as.tibble() %>% 
          summarise(erro_ev = mean(ev < alpha), erro_ks = mean(pv < alpha))
  poder
}

#########################################################
############ COMPARACOES DE DIFERENTES CENARIOS #########
#########################################################

#H0 verdadeira
set.seed(87)

dados <- matrix(rnorm(100000,0,1),ncol=100)
pf<-function(x) pnorm (x , mean=0, sd=1)
poder_ev_vs_ks(dados, phi, pf= pf)
poder_ev_vs_ks1(dados, phi, pf= pf) # as vezes da um numero diferente 

dados <- matrix(rexp(100000, rate = 2),ncol=100)
pf<-function(x) pexp (x , rate = 2)
poder_ev_vs_ks(dados, phi, pf= pf)


#H0 falsa
dados <- matrix(rnorm(100000,0,1),ncol=100)
pf<-function(x) pnorm (x , mean=5, sd=1)
poder_ev_vs_ks(dados, phi, pf= pf)

dados <- matrix(rexp(100000, rate = 2),ncol=100)
pf<-function(x) pexp (x , rate = 1)
poder_ev_vs_ks(dados, phi, pf= pf)
 