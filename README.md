# Monte_Carlo

## OBJETIVO

O presente trabalho tem como objetivo avaliar o desempenho dos estimadores de máxima verossimilhança da Distribuição Weibull.
Entretanto, não há solução fechada para o estimadores de máxima-verossimilhança desta distribuição. Para tal, vamos fazer uso 
métodos de otimização não-linear. O método escolhido para o trabalho é o BFGS, desenvolvido por Broyden-Fletcher-Goldfarb-Shannon.
Esta ferramenta compõe a classe dos métodos de quasi-newton. A vantagem deste método é sua capacidade de convergência e a sua 
facilidade. Esta útima vem da não necessidade de usarmos as segundas derivadas. Dito isso, faremos uso apenas dos gradientes da
função log-verossimilhança. Isso diminui drasticamente o custo computacional do algoritmo. 

## FUNÇÕES


