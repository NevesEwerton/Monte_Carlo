# Monte_Carlo

## OBJETIVO

O presente trabalho tem como objetivo avaliar o desempenho dos estimadores de máxima verossimilhança da Distribuição Weibull.
Entretanto, não há solução fechada para o estimadores de máxima-verossimilhança desta distribuição. Para tal, vamos fazer uso 
métodos de otimização não-linear. O método escolhido para o trabalho é o BFGS, desenvolvido por Broyden-Fletcher-Goldfarb-Shannon.
Esta ferramenta compõe a classe dos métodos de Quasi-Newton. A vantagem deste método é sua capacidade de convergência e a sua 
facilidade. Esta útima vem da não necessidade de usarmos as segundas derivadas. Dito isso, faremos uso apenas dos gradientes da
função log-verossimilhança. Isso diminui drasticamente o custo computacional do algoritmo. Pelo próprio nome, esta classe tem 
semelhanças com o método de Newton. Sua diferença mais óbvia é descrita anteriormente. Porém, de sua deferença vem sua semelhança também. Apesar de não calcular a matriz hessiana (segundas derivadas), os métodos de otimização Quasi-Newton tentam usar uma matriz de aproximação. Esta matriz se aproxima da hessiana a cada iteração. Para avaliar o desempenho dos estimadores via BFGS, vamos realizar simulações de Monte Carlo com cerca de dez mil repetições. O trabalho foi realizado em duas linguagens, Ox e R. A ideia é avaliar o desempenho das duas em rodar o programa. 

## FUNÇÕES

Aqui, vamos descrever as funções utilizadas para realizar a otimização do tipo BFGS e as respectivas bibliotecas de cada linguagem.

### Ox

A MaxBFGS maximiza a função, usando o método de Quasi-Newton desenvolvido por Broyden, Fletcher, Goldfarb, Shanno (BFGS).  A função usa as derivadas analíticas, ou as primeiras derivadas numéricas (no caso em que somente a função valor precisa ser usada: esta utiliza a função `Num1Derivative`). Usando as derivadas numéricas economiza tempo de programação, mas derivadas analíticas tendem a ser mais precisas e em uma faixa de parâmetros mais ampla. O processo de iteração não é afetado por essa escolha. 

Em Ox, a otimização do tipo BFGS é feita pela função MaxBFGS, da biblioteca `maximize`. Esta função tem os seguintes argumentos:

* `avP` -  vetor de parâmetros de dimensão px1. Primeiro, são colocados os chutes inciais. Depois da primeira iteração, as estimativas vão sendo armazenadas neste vetor;
* `adFunc` - Endereço na memória da função;
* `avScore` - Endereço dos gradientes. Se for diferente de zero, será uma matriz de dimensão px1 das primeias derivadas.
* `amHess` - É a matriz hessiana. Será sempre igual a zero para o caso da BFGS.

A função MaxBFGS retorna os seguintes valores:

* `MAX CONV`: Strong Convergence (o teste de convergência é feita usando a tolerancia \varepsilon = \varepsilon_{1}.
* `MAX_WEAK_CON`: Weak Convergence (no improvement in line search);
* `MAX_MAXIT`: No convergence (maximum no of iterations reached);
* `MAX_LINE_FAIL`: No Convergence (no improvment in line search);
* `MAX_FUNC_FAIL`: No Convergence (function evaluation failed);

Os valores padrões escolhidos para os níveis de tolerância são: \varepsilon_{1} = 10^-4 e \varepsilon_{2} = 5X10^-3.




