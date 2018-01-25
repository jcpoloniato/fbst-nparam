#funcoes uteis

#funcao de x para calcular a i-esimo funcao da base ortonormal (phi) de cossenos
phi <- function(i, x){
  sqrt(2) * cos(i * pi * x)
}

#funcao de x para calcular a i-esimo funcao da base ortonornal de  senos e cossenos
#PROBLEMA NAO ESTA ESTIMANDO BEM A DENSIDADE
# phi2 <- function(i, x){
#   return(ifelse(i%%2,sqrt(2)*sin((i+1)*pi*x),sqrt(2)*cos((i*pi*x))))
# } 

#funcao que transforma os dados usando a sequencia ortonormal de funcoes
#dados - dados a serem transformados
#phi - funcao da sequencia ortonormal
#r - tamanho da sequencia
phi_barra <- function(dados, phi, r){
  phi_barra <- rep(NA, r)
  for (i in 1:r){
    phi_barra[i] <- mean(phi(i, dados))
  } 
  phi_barra
}

#funcao que recebe como parametros os dados ja tranformados (phi barra)
#e calcula a estatistica de teste do FBST
testeFBST <- function(x, n, tau){
 sum(( n * x /sqrt(n + tau) ) ^2)
}
