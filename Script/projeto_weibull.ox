/******************************************************************
*  NOME DO PROGRAMA: projeto_weibull.ox
*
*  DISCIPLINA: MÉTODOS ESTATÍSTICOS COMPUTACIONAIS (UFPE)
*
*  PROFESSOR: FRANCISCO CRIBARI NETO
*
*  OBJETIVO: Pretendemos avaliar o desempenho dos estimadores de máxima
*			 verossimilhança da Distribuião Weibull. Porém, estes estima-
*			 dores não tem forma fechada. Dito isso, é necessário apelar
*			 métodos de otimização não-linear. O método escolhido é o 
*			 desenvolvido por Broyden-Flecther-Goldfarb-Shannon, conheci-
*			 do como BFGS. Para tal, utilizaremos a função MaxBFGS da
*			 biblioteca maximize do Ox. Para avaliar tais estimadores,
*			 faremos uso do Método de Monte Carlo. Faremos dez mil
*			 repetições, com diferentes tamanhos de amostra, para avaliar
*			 o desempenho dos estimadores e a capacidade de convergência
*			 do método BFGS.
*
*
* AUTOR DO PROGRAMA: Ewerton Neves Cardoso
*
* ÚLTIMA VERSÃO: 05/07/2022
*******************************************************************/



// Setando as bibliotecas que serão usadas.
#include <oxstd.oxh>
#include<oxdraw.h>
#include <oxprob.oxh>
#import <maximize>

/*Variável global: vai receber a amostra aleatória gerada a partir de uma distribuição weibull
usando a função ranweibull. Lembre-se: essa variável vai ser reconhecida em todos os blocos do programa.
Usar com parcimônia.*/
static decl weibull; 


fLogLik(const vP, const adFunc, const avScore, const amHess)
{
	decl alpha = vP[0]; // parâmetro alpha (scale parameter)
	decl beta = vP[1]; // parâmetro beta   (shape parameter)
	decl N = rows(weibull); // Tamanho da amostra
	// decl vone = ones(1,N);
	
/* Função Log-Verossimialhança da Distribuição Weibull com dois parâmetros */
	adFunc[0]=N*log(beta)-N*beta*log(alpha)+(beta-1)*(sumc(log(weibull))) - (1/(alpha^beta))*(sumc(weibull.^beta));
	
	//adFunc[0]=N*log(beta)-N*beta*log(alpha)+(beta-1)*(vone*(log(weibull))) - (1/(alpha^beta))*(vone*(weibull.^beta));
	
 // Gradiente Analítico	da Distribuiçao Weibull com dois parâmetros
	if(avScore)
{
	(avScore[0])[0]=-N*(beta/alpha)+beta*(alpha^(-(beta+1)))*(sumc(weibull.^beta));
	(avScore[0])[1]=(N/beta)-N*log(alpha)+sumc(log(weibull)) - (sumc(weibull.^beta))*(log(beta)/(alpha^beta));

}

	if( isnan(adFunc[0]) || isdotinf(adFunc[0]) )
		return 0; // 0 significa falha
	else
		return 1; // 1 significa sucesso
}

main()
{
	decl vp,dfunc,alpha,beta,emv_alpha,emv_beta,
	media_alpha,media_beta,ir,i,hess_matriz,var,exectime,ep;

	decl N = 250; /*Tamanho da amostra na primeira simulação. Depois faremoso mesmo para
	diferentes tamanhos de amostra. Os valores utilizados serão N=10,20,50,100,150,200,250.*/
	decl Nrep = 10000; // Número de réplicas realizadas no laço de Monte Carlo
	decl Contador_Falhas = 0; // Contador de Falhas de Convergência 

	
	exectime = timer(); // inicia contador do tempo de execução


	ranseed(1987); // Semente
	
	/* Parâmetros verdadeiros */
	alpha = 5.0;
	beta = 4.0;
	decl a = double(2^(-beta));

	/* Vetores de zeros que armazenarão valores estimador */
	emv_alpha = zeros(Nrep, 1);
	emv_beta = zeros(Nrep, 1);

	/* Esta função irá limitar o número de iterações */
	MaxControl(50, -1); 

	
	// Início do laço de Monte Carlo
	for(i=0; i<Nrep; i++)
	{
		weibull = ranweibull(N,1,a,beta); // gerador de números aleatórios da distribuição weibull
	
		vp = <2.0; 2.0>; // chutes iniciais

		/* Maximização via BFGS */
		ir = MaxBFGS(fLogLik, &vp, &dfunc, 0, TRUE);
			if(ir == MAX_CONV || ir == MAX_WEAK_CONV)
			{
				emv_alpha[i] = vp[0];  // armezando os valores estimados de alpha
				emv_beta[i] = vp[1];   // armazenando os valores estimados de beta

				/*Calculando a média dos estimadores beta e alpha*/
				media_alpha = double(meanc(emv_alpha));
				media_beta = double(meanc(emv_beta));

				/*Calculando a matriz hessiana*/
				Num2Derivative(fLogLik, vp, &hess_matriz);
				var =  diagonal(invertsym(-hess_matriz))';
				ep = sqrt(diagonal(invertsym(-hess_matriz)))';
		
			}
			else
			{
				i--;
				Contador_Falhas++;
			}
	} // Fim do laço de Monte Carlo


	/*Plotando e salvando os gráficos dos estimadores beta e alpha*/ 
	//DrawDensity(0, emv_alpha', "EMV(alpha)", TRUE, TRUE);
	//SaveDrawWindow("plot_alpha5.pdf");

	//DrawDensity(0, emv_beta', "EMV(beta)", TRUE, TRUE);
	//SaveDrawWindow("plot_beta5.pdf");
	
	print("\nAutor do Programa: Ewerton Neves Cardoso");
	print("\nNome do arquivo: projeto_weibull.ox");
	print("\nGerador de Números Aleatórios:","ranweibull()");
	print("\nSemente: ", "1987");
	print("\nData de execução: ", date());
	print("\nHorário de execução: ",time());
	print("\n");
	print("\nDistribuição: Two-parameter-Weibull, W(shape,scale)");
	print("\nMétodo de estimação: Máxima Verossimilhança");
	print("\nTipo de otimização: Método de Broyden-Fletcher-Goldfarb-Shannon (BFGS)");
	print("\nTamanho da amostra: ", "%6d", N);
	print("\nNum Rep Monte Carlo: ", "%6d", Nrep);
	print("\nNum falhas: ", "%6d", Contador_Falhas);

	print("\n");

	/*Printando a média, o desvio padrão, o erro padrão assintótico, o número de observações, número de repetições
	da simulação de monte carlo  do parâmetro beta etc*/
	print("\nMomentos MLE para beta: ", "%10.4f", "%c", {"Média", " Desvio-padrão ",
	"Assimetria ", "Curtose "}, moments(emv_beta)[1:4][]');
//	print("\nTamanho da amostra: ", "%6d", N);
//	print("\nNum Rep Monte Carlo: ", "%6d", Nrep);
//	print("\nNum falhas: ", "%6d", Contador_Falhas);
	print("\nValor verdadeiro de beta: ", "%6.3f", beta);
	print("\nMédia EMVs de beta: ", "%6.3f", double(meanc(emv_beta)));
	print("\nVariância assintótica EMV de beta: ", "%6.3f", var[1]);
	print("\nErro Padrão Assintótico EMV de beta: ", "%6.3f", ep[1]);	
	print("\nViés de beta: ", "%6.3f", double(meanc(emv_beta)-beta));
	print("\nViés relativo de beta (%): ", "%6.3f", (double(meanc(emv_beta))-beta)/beta*100);
	print("\nEQM de beta: ", "%6.3f", double(varc(emv_beta)+ (meanc(emv_beta)-beta)^2));
	// limite inferior do intervalo de confiança para beta
	decl Linfbeta = double(media_beta - (1.96)*ep[1]);
    // limite inferior do intervalo de confiança para beta
	decl Lsupbeta = double(media_beta + (1.96)*ep[1]);

	print("\n");
	
	print("Limite inferior para o intervalo de confinca de beta","%10.3f\n",Linfbeta);
	print("Limite superior para o intervalo de confinca de beta","%10.3f\n",Lsupbeta);

	print("\n");
	

	print("\nMomentos MLE para alpha: ", "%10.4f", "%c", {"Média", " Desvio-padrão ",
	"Assimetria ", "Curtose"}, moments(emv_alpha)[1:4][]');
//	print("\nTamanho da amostra: ", "%6d", N);
//	print("\nNum Rep Monte Carlo: ", "%6d", Nrep);
//	print("\nNum falhas: ", "%6d", Contador_Falhas);
	print("\nValor verdadeiro de alpha: ", "%6.3f", alpha);
	print("\nMédia EMVs de alpha: ", "%6.3f", double(meanc(emv_alpha)));
	print("\nVariância assintótica EMV de alpha: ", "%6.3f", var[0]);
	print("\nErro Padrão Assintótico EMV de alpha: ", "%6.3f", ep[0]);	
	print("\nViés de alpha: ", "%6.3f", double(meanc(emv_alpha)-alpha));
	print("\nViés relativo de alpha (%): ", "%6.3f", (double(meanc(emv_alpha))-alpha)/alpha*100);
	print("\nEQM de alpha: ", "%6.3f", double(varc(emv_alpha)+ (meanc(emv_alpha)-alpha)^2));

	print("\n");
	
	// limite inferior do intervalo de confiança para alpha
	decl Linf = double(media_alpha - (1.96)*ep[0]);
    // limite superior do intervalo de confiançaa para alpha
	decl Lsup = double(media_alpha + (1.96)*ep[1]);
	print("Limite inferior para o intervalo de confinca de alpha","%10.3f \n",Linf);
	print("Limite superior para o intervalo de confinca de alpha","%10.3f\n",Lsup);
	//print("\n");
	
	println( "\nTempo de Execução: ", timespan(exectime) );

		// SALVANDO OS DADOS EM UM ARQUIVO TEXTO
	decl estimativa=fopen("Weibull_Estimadores.txt","w");
	fprintln(estimativa, "Alpha", emv_alpha,"Beta",emv_beta, "Media de Alpha",media_alpha,"Media de Beta",media_beta,
	"TEMPO DE EXECUCAO", timespan(exectime));

	
	}