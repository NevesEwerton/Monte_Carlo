###########################################################################
#
#   NOME DO PROGRAMA: projeto_weibull.R
#
#   DISCIPLINA: MÉTODOS ESTATÍSTICOS COMPUTACIONAIS (UFPE)
#
#   PROFESSOR: FRANCISCO CRIBARI NETO
#
#   OBJETIVO: Pretendemos avaliar o desempenho dos estimadores de máxima
#	      		  verossimilhança da Distribuião Weibull. Porém, estes estima-
# 			      dores não tem forma fechada. Dito isso, é necessário apelar
#			        métodos de otimização não-linear. O método escolhido é o 
#			        desenvolvido por Broyden-Flecther-Goldfarb-Shannon, conheci-
# 			      do como BFGS. Para tal, utilizaremos a função optim da
#		          biblioteca padrão do R. Para avaliar tais estimadores,
#			        faremos uso do Método de Monte Carlo. Faremos dez mil
#			        repetições, com diferentes tamanhos de amostra, para avaliar
#			        o desempenho dos estimadores e a capacidade de convergência
#			        do método BFGS.
#
# 
#  AUTOR DO PROGRAMA: Ewerton Neves Cardoso
#
#  ÚLTIMA VERSÃO: 05/07/2022
#
############################################################################


weibull6.mc.ga.bfgs = function(nrep=10000, nobs=100, semente=2000,
                               shape1=4, shape2=5)
{
  library(moments)
  
  # função de log-verossimilhança
  logLikWeibull = function(theta){
    beta = abs(theta[1]); alpha = abs(theta[2])
    logLik = sum(dweibull(y, beta,alpha, log = TRUE))/nobs
    #logLik =  nobs*log(beta) - nobs*beta*log(alpha) + (beta-1)*(sum(log(y)) - (1/(alpha^beta))*(sum(y^beta)))
    return(logLik)
  }
  #funlça escore (gradiente)
  #scoreFnW = function(theta){
  #  beta = abs(theta[1]); alpha = abs(theta[2])
  #  cbind((beta/alpha) + beta*(alpha^(-(beta+1)))*(sum(y^beta)),
  #        (1/beta) - log(alpha) + sum(log(y))/nobs - (sum(y^beta))*(log(beta)/alpha^beta))[1,]
  #}
  
  # início da cronometragem
  tempo.inicio = Sys.time()
  
  # vetores para armanezar as estimativas
  emv_alpha = rep(0, nrep)
  emv_beta = rep(0, nrep)
  
  set.seed(2019) # semente do gerador
  contadorFalhas = 0 # contador de falhas
  
  # laço de Monte Carlo
  i = 1
  while(i <= nrep){
    
    #a=2^(beta)
    # amostra gerada
    y = rweibull(nobs, shape=4, scale=5)
    
    # chute inicial, estimador de momentos
    #ybar = mean(y)
    #ybaroneminus = 1 - ybar
    #yvar = crossprod(y - ybar)/nobs
    #secondterm = (ybar * ybaroneminus)/yvar
    #pmm = ybar * secondterm
    #qmm = ybaroneminus * secondterm
    
    shoot_beta = cbind(2)
    shoot_alpha = cbind(4)
    
    # maximização da log-verossimilhan
    ir = optim(c(shoot_beta,shoot_alpha), logLikWeibull, method="BFGS",
               control=list(fnscale=-1, reltol = 1e-12),hessian = T)
    
    # checagem de convergência
    if(ir$convergence == 0){
      emv_beta[i] = ir$par[1]
      emv_alpha[i] = ir$par[2]
      i = i + 1
    }
    else{
      contadorFalhas = contadorFalhas + 1
    }
  } # fim do laço de Monte Carlo
  # estimativas médias, vieses e vieses relativos
  beta_medio = mean(emv_beta)
  alpha_medio = mean(emv_alpha)
  beta_vies = beta_medio - shape1
  alpha_vies = alpha_medio - shape2
  beta_viesrel = (beta_vies/shape1)*100
  alpha_viesrel = (alpha_vies/shape2)*100
  
  # erros quadráticos médios
  beta_eqm = beta_vies^2 + var(emv_beta)
  alpha_eqm = alpha_vies^2 + var(emv_alpha)
  
  
  # Arredondando as estimativas de beta e alpha
  alpha_round = round(emv_alpha, 3)
  beta_round = round(emv_beta, 3)
  
  # 3 e 4 momentos 
  alpha_skew = skewness(emv_alpha)
  beta_skew = skewness(emv_beta)
  alpha_curt = kurtosis (emv_alpha)
  beta_curt = kurtosis (emv_beta)
  
  # Salvando as estimativas dos estimadores
  #write(alpha_round, file = "estimadores_BFGS250.text", ncolumns = 1, append = TRUE)
  
  #write(beta_round, file="estimador_beta_BFGS250.text", ncolumns =1, append = TRUE)
  
  
  # resultados armazenados em matriz
  mResultados = matrix(c(beta_medio, alpha_medio, beta_vies,
                         alpha_vies, beta_viesrel, alpha_viesrel,
                         beta_eqm, alpha_eqm,
                         beta_skew,alpha_skew,
                         beta_curt,alpha_curt), 6, 2, byrow=TRUE)
  
  rownames(mResultados) = c("média", "viés",
                            "viés rel.", "eqm","Assimetria","Curtose")
  
  colnames(mResultados) = c("EMV Beta", "EMV Alpha")
  
  # cálculo do tempo de execução
  tempo.fim = Sys.time()
  tempo.exec = tempo.fim - tempo.inicio
  
  # resultados colecionados em uma lista
  resultado = list(nobs=nobs, nrep=nrep, semente=semente, beta=shape1,
                   alpha=shape2, falhas=contadorFalhas, resultados=mResultados,
                   horario=tempo.inicio, tempoexec=tempo.exec)
  return(resultado)
}

weibull6.mc.ga.bfgs(nrep=10000, nobs=10, semente=2019, shape1=4, shape2=5)

