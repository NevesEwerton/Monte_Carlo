# Monte_Carlo

## OBJETIVO

O presente trabalho tem como objetivo avaliar o desempenho dos estimadores de máxima verossimilhança da Distribuição Weibull.
Entretanto, não há solução fechada para o estimadores de máxima-verossimilhança desta distribuição. Para tal, vamos fazer uso 
métodos de otimização não-linear. O método escolhido para o trabalho é o BFGS, desenvolvido por Broyden-Fletcher-Goldfarb-Shannon.
Esta ferramenta compõe a classe dos métodos de Quasi-Newton. A vantagem deste método é sua capacidade de convergência e a sua 
facilidade. Esta útima vem da não necessidade de usarmos as segundas derivadas. Dito isso, faremos uso apenas dos gradientes da
função log-verossimilhança. Isso diminui drasticamente o custo computacional do algoritmo. Pelo próprio nome, esta classe tem 
semelhanças com o método de Newton. Sua diferença mais óbvia é descrita anteriormente. Porém, de sua deferença vem sua semelhança também. Apesar de não calcular a matriz hessiana (segundas derivadas), os métodos de otimização Quasi-Newton tentam usar uma matriz de aproximação. Esta matriz se aproxima da hessiana a cada iteração. Para avaliar o desempenho dos estimadores via BFGS, vamos realizar simulações de Monte Carlo com cerca de dez mil repetições. O tamanho da amostra será aumentada progressivamente. 


