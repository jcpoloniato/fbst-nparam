setwd('/home/joao_polo/Dropbox/Mestrado/Dissertacao/fbst-nparam')
source('funcoes-uteis.R') 

#carregando pacotes de devem ser usados
library(nortest)
library(tidyverse)

#funcao que calcula o evalor não modificado
#dados - dados
#phi -base ortonormal
#r - tamanho da sequencia da base ortonormal
#peso - funcao do peso para cada elemento da sequencia
ev_nmod <- function(dados, phi, r, peso){
  n <- length(dados)
  tau <- peso(1:r)
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
  tau <- peso(1:r)
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

###### H_0 verdadeira #######

# evm0 <- rep(NA,1000)
# evnm0 <- rep(NA,1000)
# for (i in 1:1000){
#   dados <- runif(10000)
#   evm0[i] <- ev_mod(dados, phi, r = 50, peso = 5)
#   evnm0[i] <- ev_nmod(dados, phi, r = 50, peso = 5)
# }
# par(mfrow=c(1,2))
# hist(evm0)
# hist(evnm0)

##### H_0 falsa #############

# evm1 <- rep(NA,1000)
# evnm1 <- rep(NA,1000)
# for (i in 1:1000){
#   dados <- rnorm(10000, mean=5, sd=1)
#   evm1[i] <- ev_mod(dados, phi, 50, 5)
#   evnm1[i] <- ev_nmod(dados, phi, 50, 5)
# }
# hist(evm1)
# hist(evnm1)

#################################################################################

#gerar dados e aplicar uma F nos dados e 
#verificar se os dados modificado tem distribuicao uniforme

#funcao que vai calcular o evalor nos dados transformados (padronizados)
#dados - dados
#phi -base ortonormal
#r - tamanho da sequencia da base ortonormal
#peso - funcao de r para cada elemento da sequencia
#simu - numero de observacoes que devem ser geradas da distribuicao de comparacao
#pf - funcao F que deve ser usada para padronizar os dados
ev_mod_pf<-function(dados, phi, r = 50, peso = function(r) 5^r, simu=1000, pf){
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

#essa funcao retorna uma matriz com a primeira coluna de evalor e a outra de pvalor ks
evalor_vs_ks <- function(matriz_dados, phi, r = 50, peso = function(r) 5^r, pf){
  ev <- apply(matriz_dados, 2, ev_mod_pf, phi = phi, r = r, peso = peso, pf = pf)
  pv <- apply(matriz_dados, 2, pvalor_ks, pf = pf)
  cbind(ev, pv)
}

#funcao para ver a proporcao de rejeicao fixado um alpha
poder_ev_vs_ks <- function(matriz_dados, phi, r = 50, peso = function(r) 5^r, pf, alpha = 0.05){
  poder <- colMeans(evalor_vs_ks(matriz_dados, phi, r, peso, pf) < alpha)
  poder
}

##mesma funcao usando pipe
# poder_ev_vs_ks1 <- function(matriz_dados, phi, r = 50, peso = 5, pf, alpha = 0.05){
#   poder <- evalor_vs_ks(dados, phi, r, peso, pf) %>% 
#           as.tibble() %>% 
#           summarise(erro_ev = mean(ev < alpha), erro_ks = mean(pv < alpha))
#   poder
# }

#########################################################
############ COMPARACOES DE DIFERENTES CENARIOS #########
#########################################################

#################################
#### H0 verdadeira ##############
#################################

dados <- matrix(rnorm(100000,0,1),ncol=1000)
pf<-function(x) pnorm (x , mean=0, sd=1)
poder_ev_vs_ks(dados, phi, pf= pf)

dados <- matrix(rexp(100000, rate = 2),ncol=100)
pf<-function(x) pexp (x , rate = 2)
poder_ev_vs_ks(dados, phi, pf= pf)

#################################
#### H0 falso ###################
#################################
dados <- matrix(rnorm(100000,0,1),ncol=100)
pf<-function(x) pnorm (x , mean=5, sd=1)
poder_ev_vs_ks(dados, phi, pf= pf)

dados <- matrix(rexp(100000, rate = 2),ncol=100)
pf<-function(x) pexp (x , rate = 1)
poder_ev_vs_ks(dados, phi, pf= pf)

#################################
## Obtendo a funcao poder #######
#################################

#distribuicao normal
medias<- seq(-1,1,by=0.01)
poderes<-matrix(NA,nrow=length(medias),ncol=2)
pf<-function(x) pnorm (x , mean=0, sd=1)

for (i in 1:length(medias)){
  dados <- matrix(rnorm(100000,medias[i],1),ncol=100)
  poderes[i,]<-poder_ev_vs_ks(dados, phi, pf= pf)
}

dat<-data.frame(medias=c(medias, medias),
                poder=c(poderes[,1], poderes[,2]), 
                teste=c(rep('ev', length(medias)), rep('pv', length(medias))))
saveRDS(dat, file='normal')

#distribuicao exponencial
medias<- seq(0.5,2,by=0.01)
poderes<-matrix(NA,nrow=length(medias),ncol=2)
pf<-function(x) pexp(x , rate=1)

for (i in 1:length(medias)){
  dados <- matrix(rexp(100000,rate=medias[i]),ncol=100)
  poderes[i,]<-poder_ev_vs_ks(dados, phi, pf= pf)
}

dat<-data.frame(medias=c(medias, medias),
                poder=c(poderes[,1], poderes[,2]), 
                teste=c(rep('ev', length(medias)), rep('pv', length(medias))))
saveRDS(dat, file='exponencial')

## funcao que gera mistura de distribuicoes normais

#pf - a distribuicao a ser testada nos dados
#probs - peso para cada distribuicao na mistura
#phi - base ortonormal
#medias - medias das normais
#variancias - variancias das normais
#n - tamanho das amostras
#replicas - quantidade de replicas das amostras
#O TAMANHO DO VETOR DE PROBS E VARIANCIAS DEVEM SER IGUAIS
mist_norm <- function (pf=function(x) pnorm (x , mean=0, sd=1), phi, probs, medias, variancias, n=1000 , replicas=100){
  m <- length(medias)
  k <- length(probs)
  nt <- n*replicas 
  poderes<-matrix(NA,nrow=m,ncol=2)
  
  for (i in 1:m){
    vetor<-rep(NA,nt)
    
    aux<-rmultinom(nt, size = 1, prob = probs)
    delta<-rep(NA,nt)
    j<-1
    while (j < (k+1)){
      delta[aux[j,]==1]<-j
      j <- j+1
    }

    if(k%%2 == 0){
      j<-1
      while(j < (k+1)){
        vetor[delta==j]<-rnorm(sum(delta==j),-medias[i]*j,variancias[j])
        j<-j+2
      }
      j<-2
      while(j < (k+1)){
        vetor[delta==j]<-rnorm(sum(delta==j),medias[i]*(j-1),variancias[j])
        j<-j+2
      }
    } else {
      vetor<-rnorm(nt, mean=0, sd=variancias[1])
      j<-2
      while(j < (k+1)){
        vetor[delta==j]<-rnorm(sum(delta==j),-medias[i]*(j-1),variancias[j])
        j<-j+2
      }
      j<-3
      while(j < (k+1)){
        vetor[delta==j]<-rnorm(sum(delta==j),medias[i]*(j-2),variancias[j])
        j<-j+2
      }
    }
    
    dados<-matrix(vetor,ncol=n,byrow=T)
    poderes[i,]<-poder_ev_vs_ks(dados, phi, pf= pf)
  }
  
  poderes
}

#simulando poder dos testes para diferentes normais multimodais

#2 modas
medias<-seq(0,2,by=0.01)
probs<-c(1/2,1/2)
variancias<-c(1,1)
test <- mist_norm(phi = phi, probs = probs, medias = medias, variancias = variancias)
dat<-data.frame(medias=c(medias, medias),
                poder=c(test[,1], test[,2]), 
                teste=c(rep('ev', length(medias)), rep('pv', length(medias))))
saveRDS(dat, file='normal2')

#3 modas
medias<-seq(0,2,by=0.01)
probs<-rep(1/3,3)
variancias<-c(1,1,1)
test <- mist_norm(phi = phi, probs = probs, medias = medias, variancias = variancias)
dat<-data.frame(medias=c(medias, medias),
                poder=c(test[,1], test[,2]), 
                teste=c(rep('ev', length(medias)), rep('pv', length(medias))))
saveRDS(dat, file='normal3')

#4 modas
medias<-seq(0,2,by=0.01)
probs<-rep(1/4,4)
variancias<-c(1,1,1,1)
test <- mist_norm(phi = phi, probs = probs, medias = medias, variancias = variancias)
dat<-data.frame(medias=c(medias, medias),
                poder=c(test[,1], test[,2]), 
                teste=c(rep('ev', length(medias)), rep('pv', length(medias))))
saveRDS(dat, file='normal4')

#5 modas
medias<-seq(0,2,by=0.01)
probs<-rep(1/5,5)
variancias<-c(1,1,1,1,1)
test <- mist_norm(phi = phi, probs = probs, medias = medias, variancias = variancias)
dat<-data.frame(medias=c(medias, medias),
                poder=c(test[,1], test[,2]), 
                teste=c(rep('ev', length(medias)), rep('pv', length(medias))))
saveRDS(dat, file='normal5')

#########
#distribuicao beta
#########

#funcao para testar alguma beta especifica variando shape1=shape2
tbeta <- function(pf, shape, n=100, replicas=1000){
  nt <- n*replicas
  poderes<-matrix(NA,nrow=length(shape),ncol=2)
  for (i in 1:length(shape)){
    dados<-matrix(rbeta(nt,shape1=shape[i],shape2=shape[i]),ncol=replicas,byrow=T)
    poderes[i,]<-poder_ev_vs_ks(dados, phi, pf= pf)
  }
  poderes
}

#testar beta 11
shape <- seq(0,3,by=0.01)
pf<-function(x) pbeta (x , shape1=1, shape2=1)
test <- tbeta(pf = pf, shape = shape)

dat<-data.frame(shape=c(shape, shape),
                poder=c(test[,1], test[,2]), 
                teste=c(rep('ev', length(shape)), rep('pv', length(shape))))
saveRDS(dat, file = 'beta 11')

#testar beta 1/2 1/2
shape <- seq(0,3,by=0.01)
pf<-function(x) pbeta (x , shape1=1/2, shape2=1/2)
test <- tbeta(pf = pf, shape = shape)

dat<-data.frame(shape=c(shape, shape),
                poder=c(test[,1], test[,2]), 
                teste=c(rep('ev', length(shape)), rep('pv', length(shape))))
saveRDS(dat, file = 'beta 0.5 0.5')

#############################
###### GRAFICOS #############
#############################

#
dat<-readRDS(file='normal')
ggplot(data = dat, aes(x=medias, y=poder, group=teste , colour=teste)) +
  geom_line() +
  scale_colour_discrete(labels=c("FBST", "KS"))+
  geom_hline(yintercept = 0.05, linetype = "dashed", col='green') +
  ggtitle("") +
  labs(y = "Poder do teste", x = "μ")+
  theme_bw() +
  theme(#axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #panel.background = element_blank(),
        #legend.position="top",
        legend.title=element_blank(),
        axis.text=element_text(size=15),
        legend.text = element_text(size=15),
        axis.title.x = element_text(size = 20, angle = 0),
        axis.title.y = element_text(size = 20, angle = 90))

dif <- dat$poder[1:201]-dat$poder[202:402]
dat<-data.frame(medias = dat$medias,
                diferenca = dif)

ggplot(data = dat, aes(x=medias, y=diferenca)) +
  geom_line()+
  geom_hline(yintercept = 0, col='green', linetype = "dashed") +
  labs(y = "Poder FBST - Poder KS", x = "μ")+  theme_bw() +
  ylim(c(-0.15, 0.6)) +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

#
dat<-readRDS(file='exponencial')
ggplot(data = dat, aes(x=medias, y=poder, group=teste , colour=teste)) +
  geom_line()+
  scale_colour_discrete(labels=c("FBST", "KS"))+
  geom_hline(yintercept = 0.05, linetype = "dashed", col='green') +
  ggtitle("") +
  labs(y = "Poder do teste", x = "λ")+
  theme_bw() +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

dif <- dat$poder[1:151]-dat$poder[152:302]
dat<-data.frame(medias = dat$medias,
                diferenca = dif)

ggplot(data = dat, aes(x=medias, y=diferenca)) +
  geom_line()+
  geom_hline(yintercept = 0, col='green', linetype = "dashed") +
  labs(y = "Poder FBST - Poder KS", x = "λ")+  theme_bw() +
  ylim(c(-0.15, 0.6)) +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

#
dat<-readRDS(file='normal2')
ggplot(data = dat, aes(x=medias, y=poder, group=teste , colour=teste)) +
  geom_line()+
  scale_colour_discrete(labels=c("FBST", "KS"))+
  geom_hline(yintercept = 0.05, linetype = "dashed", col='green') +
  ggtitle("") +
  labs(y = "Poder do teste", x = "μ")+
  theme_bw() +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

dif <- dat$poder[1:201]-dat$poder[202:402]
dat<-data.frame(medias = dat$medias,
                diferenca = dif)

ggplot(data = dat, aes(x=medias, y=diferenca)) +
  geom_line()+
  geom_hline(yintercept = 0, col='green', linetype = "dashed") +
  labs(y = "Poder FBST - Poder KS", x = "μ")+  theme_bw() +
  ylim(c(-0.15, 0.6)) +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

#
dat<-readRDS(file='normal3')
ggplot(data = dat, aes(x=medias, y=poder, group=teste , colour=teste)) +
  geom_line()+
  scale_colour_discrete(labels=c("FBST", "KS"))+
  geom_hline(yintercept = 0.05, linetype = "dashed", col='green') +
  ggtitle("") +
  labs(y = "Poder do teste", x = "μ")+
  theme_bw() +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

dif <- dat$poder[1:201]-dat$poder[202:402]
dat<-data.frame(medias = dat$medias,
                diferenca = dif)

ggplot(data = dat, aes(x=medias, y=diferenca)) +
  geom_line()+
  geom_hline(yintercept = 0, col='green', linetype = "dashed") +
  labs(y = "Poder FBST - Poder KS", x = "μ")+  theme_bw() +
  ylim(c(-0.15, 0.6)) +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

#
dat<-readRDS(file='normal4')
ggplot(data = dat, aes(x=medias, y=poder, group=teste , colour=teste)) +
  geom_line()+
  scale_colour_discrete(labels=c("FBST", "KS"))+
  geom_hline(yintercept = 0.05, linetype = "dashed", col='green') +
  ggtitle("") +
  labs(y = "Poder do teste", x = "μ")+
  theme_bw() +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

dif <- dat$poder[1:201]-dat$poder[202:402]
dat<-data.frame(medias = dat$medias,
                diferenca = dif)

ggplot(data = dat, aes(x=medias, y=diferenca)) +
  geom_line()+
  geom_hline(yintercept = 0, col='green', linetype = "dashed") +
  labs(y = "Poder FBST - Poder KS", x = "μ")+  theme_bw() +
  ylim(c(-0.15, 0.6)) +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

#
dat<-readRDS(file='normal5')
ggplot(data = dat, aes(x=medias, y=poder, group=teste , colour=teste)) +
  geom_line()+
  scale_colour_discrete(labels=c("FBST", "KS"))+
  geom_hline(yintercept = 0.05, linetype = "dashed", col='green') +
  ggtitle("") +
  labs(y = "Poder do teste", x = "μ")+
  theme_bw() +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

dif <- dat$poder[1:201]-dat$poder[202:402]
dat<-data.frame(medias = dat$medias,
                diferenca = dif)

ggplot(data = dat, aes(x=medias, y=diferenca)) +
  geom_line()+
  geom_hline(yintercept = 0, col='green', linetype = "dashed") +
  labs(y = "Poder FBST - Poder KS", x = "μ")+  theme_bw() +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

#
dat<-readRDS(file='beta 0.5 0.5')
ggplot(data = dat, aes(x=shape, y=poder, group=teste , colour=teste)) +
  geom_line()+
  scale_colour_discrete(labels=c("FBST", "KS"))+
  geom_hline(yintercept = 0.05, linetype = "dashed", col='green') +
  ggtitle("") +
  labs(y = "Poder do teste", x = "α = β")+
  theme_bw() +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

dif <- dat$poder[1:301]-dat$poder[302:602]
dat<-data.frame(shape = dat$shape,
                diferenca = dif)

ggplot(data = dat, aes(x=shape, y=diferenca)) +
  geom_line()+
  geom_hline(yintercept = 0, col='green', linetype = "dashed") +
  labs(y = "Poder FBST - Poder KS", x = "α = β")+  theme_bw() +
  ylim(c(-0.15, 0.6)) +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

#
dat<-readRDS(file='beta 11')
ggplot(data = dat, aes(x=shape, y=poder, group=teste , colour=teste)) +
  geom_line()+
  scale_colour_discrete(labels=c("FBST", "KS"))+
  geom_hline(yintercept = 0.05, linetype = "dashed", col='green') +
  ggtitle("") +
  labs(y = "Poder do teste", x = "α = β")+
  theme_bw() +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

dif <- dat$poder[1:301]-dat$poder[302:602]
dat<-data.frame(shape = dat$shape,
                diferenca = dif)

ggplot(data = dat, aes(x=shape, y=diferenca)) +
  geom_line()+
  geom_hline(yintercept = 0, col='green', linetype = "dashed") +
  labs(y = "Poder FBST - Poder KS", x = "α = β")+  theme_bw() +
  ylim(c(-0.15, 0.6)) +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

#DISTRIBUICAO TRIGONOMÉTRICA

densidade<-function(x,k){
  abs(sin(pi*k*x)*pi/2)
}

x<-seq(0,1,by=0.001)
fy<-densidade(x,k=5)
dens<-data.frame(x=x, fy=fy)
ggplot(data = dens, aes(x=x, y=fy))+
  geom_line()+
  ggtitle("") +
  labs(y = "f(x,k=5)", x = "")+
  theme_bw() +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

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

tseno<-function(pf = function(x) punif(x, 0, 1), phi, k, n = 1000, replicas = 1000, peso){
  nt <- n * replicas
  vetor <- gerar(n = nt, k = k)
  dados<-matrix(vetor, ncol=replicas, byrow=T)
  poder_ev_vs_ks(dados, phi, pf = pf, peso = peso)
}

#mudando as prioris

#n=100

k<-seq(1,20,by=1)
poderes<- matrix(NA, nrow=length(k), ncol=5)
for (i in 1:length(k)){
  poderes[i,1:2]<-tseno(phi=phi, k = k[i], n=100, peso=function(x) 5^x)
  poderes[i,3]<-tseno(phi=phi, k = k[i], n=100, peso=function(x) x^5)[1]
  poderes[i,4]<-tseno(phi=phi, k = k[i], n=100, peso=function(x) 2^x)[1]
  poderes[i,5]<-tseno(phi=phi, k = k[i], n=100, peso=function(x) x^2)[1]
  print(i)
}

dat<-data.frame(k=rep(k,5),
                poder=c(poderes[,1], poderes[,2], poderes[,3], 
                        poderes[,4], poderes[,5]), 
                teste=c(rep('FBST 5^r', length(k)), 
                        rep('KS', length(k)),
                        rep('FBST r^5', length(k)),
                        rep('FBST 2^r', length(k)),
                        rep('FBST r^2', length(k))
                        )
                )
saveRDS(dat, file = 'seno100')

#n=1000

k<-seq(1,20,by=1)
poderes<- matrix(NA, nrow=length(k), ncol=5)
for (i in 1:length(k)){
  poderes[i,1:2]<-tseno(phi=phi, k = k[i], n=1000, peso=function(x) 5^x)
  poderes[i,3]<-tseno(phi=phi, k = k[i], n=1000, peso=function(x) x^5)[1]
  poderes[i,4]<-tseno(phi=phi, k = k[i], n=1000, peso=function(x) 2^x)[1]
  poderes[i,5]<-tseno(phi=phi, k = k[i], n=1000, peso=function(x) x^2)[1]
  print(i)
}

dat<-data.frame(k=rep(k,5),
                poder=c(poderes[,1], poderes[,2], poderes[,3], 
                        poderes[,4], poderes[,5]), 
                teste=c(rep('FBST 5^r', length(k)), 
                        rep('KS', length(k)),
                        rep('FBST r^5', length(k)),
                        rep('FBST 2^r', length(k)),
                        rep('FBST r^2', length(k))
                )
)
saveRDS(dat, file = 'seno1000')

#graficos 

dat<-readRDS(file = 'seno100')
ggplot(data = dat, aes(x=k, y=poder, group=teste , colour=teste)) +
  geom_line()+
  #scale_colour_discrete(labels=c("FBST", "KS"))+
  geom_hline(yintercept = 0.05, linetype = "dashed", col='red') +
  ggtitle("") +
  labs(y = "Poder do teste", x = "k")+
  theme_bw() +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

dat<-readRDS(file = 'seno1000')
ggplot(data = dat, aes(x=k, y=poder, group=teste , colour=teste)) +
  geom_line()+
  #scale_colour_discrete(labels=c("FBST", "KS"))+
  geom_hline(yintercept = 0.05, linetype = "dashed", col='red') +
  ggtitle("") +
  labs(y = "Poder do teste", x = "k")+
  theme_bw() +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))