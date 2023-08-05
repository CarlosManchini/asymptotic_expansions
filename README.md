## Expansões Assintóticas

As expansões assintóticas são essenciais para analisar e entender fenômenos complexos em diferentes campos da matemática, física, estatística, engenharia e ciências aplicadas. Elas são especialmente úteis em situações complexas  quando as abordagens tradicionais se tornam impraticáveis ou insuficientes pois a expressão analítica exata é difícil ou impossível de obter. Assim, elas permitem uma descrição **aproximada** de funções complexas em regimes extremos em torno de um ponto específico ou quando um parâmetro tende a valores muito grandes ou pequenos. 

No meio estatístico e na teoria das probabilidades, as expansões de Edgeworth e ponto de sela (Saddle point) são úteis em problemas relacionados a distribuições de médias amostrais e outras estatísticas em amostras grandes. Ou seja, obter aproximações assintóticas para a densidade de probabilidade de uma soma de variáveis aleatórias independentes e identicamente distribuídas.
Em especial, a expansão de Edgeworth melhora a precisão da aproximação de Laplace, que é uma aproximação de primeira ordem baseada no método de Laplace (o qual possui uma exemplo prático neste repositório em [Aplicação da aproximação de Laplace](https://github.com/CarlosManchini/asymptotic_expansions/blob/main/Aplica%C3%A7%C3%A3o%20da%20aproxima%C3%A7%C3%A3o%20de%20Laplace.pdf).


Esses métodos de aproximação nos fornecem a capacidade de reduzir a quantidade de variáveis aleatórias independentes de uma  estimativa e pode gerar resultados mais estáveis, eficientes, interpretáveis e teoricamente válidos. Isso é particularmente importante quando se trabalha com dados complexos e modelos estatísticos, onde a simplicidade e a precisão são considerações cruciais. 

a comparação com a distribuição normal pelo Teorema Central do Limite (TCL). 

um aplicativo desenvolvido em \texttt{Shiny} para explorar diferentes cenários em uma aplicação da aproximação de Edgeworth para a densidade da normal inversa.
#### Distribuição Normal Inversa

Para aplicação das expansões de Edgeworth e ponto de sela, irei utilizar a distribuição normal inversa para realizar a aproximação.


[Shiny App](https://ufsm.shinyapps.io/appig/ "Approximations of the inverse Gaussian distribution")
