#funcoes uteis

#funcao de x para calcular a i-esimo funcao da base ortonormal (phi) de cossenos
phi <- function(i, x){
  invisible( sqrt(2) * cos(i * pi * x) )
}

#funcao que transforma os dados usando a sequencia ortonormal de funcoes
#dados - dados a serem transformados
#phi - funcao da sequencia ortonormal
#r - tamanho da sequencia
phi_barra <- function(dados, phi, r){
  phi_barra <- rep(NA, r)
  for (i in 1:r){
    phi_barra[i] <- mean(phi(i, dados))
  } 
  invisible(phi_barra)
}

#funcao que recebe como parametros os dados ja tranformados (phi barra)
#e calcula a estatistica de teste do FBST
testeFBST <- function(x, n, tau){
  invisible( sum(( n* x /sqrt(n + tau) ) ^2) )
}