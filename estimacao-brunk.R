setwd('/home/joao_polo/Dropbox/Disciplinas/Mestrado/Dissertação/fbst-nparam')
source('funcoes-uteis.R') 

#dados = amostra da funcao f
#phi = base ortonormal
#chute = uma lista com os membros:
# f0 = chute de densidade que se ajusta a f real
# pf0 = acumulada da funcao chute
#r = quantidade de fun??es ortogonais a serem usadas
#d = peso da base de ordem (d^r)
#enp retorna uma lista com os membros: 
# f_chapeu: a densidade estimada 
# beta_chapeu: os coeficientes estimados da base 
enp <- function(dados, phi, chute, r, d, plot_out = FALSE){
  
  f0 <- chute$f0 #densidade da funcao chute
  pf0 <- chute$pf0 #acumulada da funcao chute
  
  n <- length(dados)
  z <- pf0(dados) #criar a variavel z para que a base (funcao phi) seja ortonormal
  beta_chapeu <- rep(NA, r)
  phi_barra <- rep(NA, r)
  
  #estimar os coeficientes betas que acompanham cada elemento da base
  for (i in 1:r){
    phi_barra[i] <- mean( phi(i, z) )
    beta_chapeu[i] <- (n/(n-1+d^i))*phi_barra[i]
  }
  
  #estimativa da funcao densidade de z (funcao p)
  p_chapeu <- function(x){
    out <- 1
    for (i in 1:r){
      out <- out + beta_chapeu[i] * phi(i, x)
    }
    return(out)
  }
  
  #estimativa da funcao densidade dos dados (f)
  f_chapeu <- function(x){
    out <- p_chapeu( pf0(x) )* f0(x)
    return(out)
  }
  
  if(plot_out == TRUE){
    #fazer um plot das funcoes densidades:
    # azul - densidade estimada pelo metodo do brunk
    # preto - densidade estimada pelo R
    plot(density(dados))
    y <- seq(min(dados), max(dados), by=0.01)
    lines(y, f_chapeu(y), col='blue')
  }

  return(list(f_chapeu=f_chapeu,
              beta_chapeu=beta_chapeu))
}

#log_lik = log verossimilhanca
#dados = dados
#chute0 = chute inicial para os parametros a serem estimados
#a funcao maximiza a log verossimlhanca e 
#retorna as estimativas de cada um dos parametros da densidade
max_log_lik <- function(log_lik, dados, chute0, lower=-Inf, upper=Inf){
  fit <- optim(par=chute0, fn=log_lik, dados=dados,
               control=list(fnscale=-1),method='L-BFGS-B',lower=lower,upper=upper)
  #  fit<-optim(par=fit$par,fn=log_lik,dados=dados,
  #             control=list(fnscale=-1),method='L-BFGS-B',lower=lower,upper=upper)
  if (fit$convergence==0) 
    theta_chapeu=fit$par
  else 
    print('Nao convergiu')
  return(theta_chapeu)
}

#funcoes que recebe como argumento os dados e retorna a densidade e a acumulada
#como parametros estimados pelo metodo da maxima verossilhanca
chute_norm_mv <- function(dados){
  lik_norm <- function(theta, dados){
    out <- sum(dnorm(dados, mean=theta[1], sd=theta[2], log=T))
    return(out)
  }
  m0 <- mean(dados)
  s0 <- sd(dados)
  theta <- max_log_lik(lik_norm, dados, c(m0,s0))
  m0 <- theta[1]
  s0 <- theta[2]
  return(list(f0=function(x) dnorm(x,mean=m0,sd=s0),
              pf0=function(x) pnorm(x,mean=m0,sd=s0)
  ))
}
chute_gama_mv <- function(dados){
  lik_gama <- function(theta, dados){
    out <- sum(dgamma(dados, shape = theta[1], scale = theta[2], log=T))
    return(out)
  }
  v0 <- sd(dados)^2
  m0 <- mean(dados)
  s0 <- v0/m0
  a0 <- m0/s0
  theta <- max_log_lik(lik_gama, dados, c(a0,s0), lower=c(0,0))
  a0 <- theta[1]
  s0 <- theta[2]
  return(list(pf0 = pf0 <- function(x) pgamma(x, shape=a0, scale=s0),
              f0 = f0 <- function(x) dgamma(x, shape=a0, scale=s0)
  ))
}
chute_exp_mv <- function(dados){
  lik_exp <- function(theta, dados){
    out <- sum(dexp(dados, rate=theta, log=T))
    return(out)
  }
  rate <- 1/mean(dados)
  theta <- max_log_lik(lik_exp, dados, rate, lower=0)
  rate <- theta
  
  return(list(  pf0 = pf0 <- function(x) pexp(x, rate=rate),
                f0 = f0 <- function(x) dexp(x, rate=rate)))
}
chute_beta_mv <- function(dados){
  lik_beta <- function(theta, dados){
    out <- sum(dbeta(dados, shape1 = theta[1], shape2 = theta[2], log=T))
    return(out)
  }
  v0 <- sd(dados)^2
  m0 <- mean(dados)
  a0 <- (-m0^3+m0^2-m0*v0)/v0
  b0 <- (m0^3-2*m0^2+m0*v0+m0-v0)/v0
  theta <- max_log_lik(lik_beta, dados, c(a0,b0), lower=c(0,0))
  a0 <- theta[1]
  b0 <- theta[2]
  return(list(pf0 = pf0 <- function(x) pbeta(x, shape1 = a0, shape2 = b0),
              f0 = f0 <- function(x) dbeta(x, shape1 = a0, shape2 = b0)))
}

#funcoes que recebe como argumento os dados e retorna a densidade e a acumulada
#como parametros estimados pelo metodo dos momentos
chute_norm_mm <- function(dados){
  m0 <- mean(dados)
  s0 <- sd(dados)
  return(list(f0 = f0 <- function(x) dnorm(x, mean = m0, sd = s0),
              pf0 = pf0 <- function(x) pnorm(x, mean = m0, sd = s0)
              ))
}
chute_gama_mm <- function(dados){
  v0 <- sd(dados)^2
  m0 <- mean(dados)
  s0 <- v0/m0
  a0 <- m0/s0
  return(list(pf0 = pf0 <- function(x) pgamma(x, shape = a0, scale = s0),
              f0 = f0 <- function(x) dgamma(x, shape = a0, scale = s0)
  ))
}
chute_exp_mm <- function(dados){
  rate <- 1/mean(dados)
  return(list(  pf0 = pf0 <- function(x) pexp(x, rate=rate),
                f0 = f0 <- function(x) dexp(x, rate=rate)))
}
chute_beta_mm <- function(dados){
  v0 <- sd(dados)^2
  m0 <- mean(dados)
  a0 <- (-m0^3+m0^2-m0*v0)/v0
  b0 <- (m0^3-2*m0^2+m0*v0+m0-v0)/v0
return(list( pf0 = pf0 <- function(x) pbeta(x, shape1 = a0, shape2=b0),
            f0 = f0 <- function(x) dbeta(x, shape1 = a0, shape2=b0)))
}

#Testes MV

dados <- rnorm(100, mean = 5, sd = 10)
est <- enp(dados, phi, chute_norm_mv(dados), 50, 5)

dados <- rbeta(10000, shape1 = 1.1, shape2 = 1)
est <- enp(dados, phi, chute_beta_mv(dados), 50, 5)

dados <- rgamma(10000, shape = 5, scale = 1/2)
est <- enp(dados, phi, chute_gama_mv(dados), 50, 5)

dados <- rexp(100000, 5)
est <- enp(dados, phi, chute_exp_mv(dados), 1, 5)

dados <- rexp(100000, 5)
est <- enp(dados, phi, chute_norm_mv(dados), 50, 5)

dados <- rexp(100000, 5)
est <- enp(dados, phi, chute_gama_mv(dados), 50, 5)

#Testes Momentos

dados <- rnorm(100, mean = 5, sd = 10)
est <- enp(dados, phi, chute_norm_mm(dados), 50, 5)

dados <- rnorm(100000, mean = 0, sd = 10)
est <- enp(dados, phi, chute_norm_mm(dados), 50, 5)

dados <- rnorm(100000, mean = 0, sd = 10)
chute <- list( pf0 = function(x) pnorm(x, mean = 5, sd = 1),
            f0 = function(x) dnorm(x, mean = 5, sd = 1))
est <- enp(dados, phi, chute, 5, 5)

dados <- rexp(100000, 5)
est <- enp(dados, phi, chute_norm_mm(dados), 50, 5)

dados <- rexp(100000, 5)
est <- enp(dados, phi, chute_exp_mm(dados), 1, 5)

dados <- rnorm(100000, mean = 2)
est <- enp(dados, phi, chute_exp_mm(dados), 50, 5)

dados <- rgamma(10000, 5)
est <- enp(dados, phi,chute_gama_mv(dados), 50, 5)

dados <- rgamma(10000, shape = 5, scale = 1/2)
est <- enp(dados, phi, chute_gama_mm(dados), 50, 5)