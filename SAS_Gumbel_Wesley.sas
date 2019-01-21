proc delete data=_all_;
run;
goptions reset=all;
*Wesley oliveira Furriel RA: 61493;
*----------------------------------------------------------------------------------------------;
*Criando uma macro para plotar os gráficos ;
%macro graph(dist=, var=, dado=,nome=);
ods html style=harvest; 
options nodate nonumber;
options orientation=portrait;
ods pdf file="D:\Estatística\ESTATÍSTICA_COMPUTACIONAL_II\
MLE\Gumbel_Latex\&nome..pdf" startpage= yes;
legend1 position=(top left inside) mode=protect
          value=("mu=0.5  sigma=2")  label=none;
axis1 minor=none order=(-5 to 20 by 2);
symbol1 i=j v=none l=1 c=blue w=3;
proc gplot data=&dado;
	plot &dist*&var/overlay haxis=axis1 noframe autovref legend=legend1;
run;
ODS pdf close;
ods close;
%mend graph;

*----------------------------------------------------------------------------------------------;
* Criando as funções;
* Funcao densidade de probabilidade (d) fdp;
data dist1;
 do x = -5 to 20 by 0.1;
 	mu=0.5;
	sigma=2;
	dgumbel= exp((-x + mu) / sigma - exp((-x + mu) / sigma))/sigma;
    output;
 end;
run;
*Funcao distribuicao de prob. (p) dist. acumulada;
data dist2;
 do q = -5 to 20 by 0.1;
 	mu=0.5;
	sigma=2;
	pgumbel= exp(-exp((-q + mu) / sigma));
    output;
 end;
run;
%graph(dist=dgumbel, var=x, dado=dist1,nome=densidade);
%graph(dist=pgumbel, var=q, dado=dist2,nome=acumulada);
*Funcao quantil (q);
data dist3;
 do p = 0 to 1 by 0.01;
 	mu=0.5;
	sigma=2;
	qgumbel= mu - sigma*log(-log(p));
    output;
 end;
run;

*----------------------------------------------------------------------------------------------;
*Funcao para valores aleatorios;
data dist4(drop=mu sigma u i);
do i= 1 to 10;
call streaminit(666);/*Semente*/
	mu=0.5;
	sigma=2;
	u=rand("uniform");
	rgumbel=-log(-log(u))*sigma+mu; 
  output;
end;
run;
*Estimação dos Parametros;
proc nlp data=dist4 tech=nra;
	max l;
	parms mu=0.5,sigma=2;
	x=rgumbel;
	aux=(-x + mu) / sigma;
	l = -log(sigma)+(log(exp(aux-exp(aux))));
run;

*--------------------------------------------------------------------------------------------;
*Funções Gumbel com a proc IML;
proc delete data=_all_;
run;
goptions reset=all;

proc iml;
/*Funcao densidade de probabilidade (d) fdp*/
start dgumbel(x, mu, sigma);
	fdp= exp((-x + mu) / sigma - exp((-x + mu) / sigma))/sigma;
	return(fdp);
finish;
x = dgumbel(2,1.5,2); 
print(x);
quit;

*----------------------------------------------------------------------------------------------;
	/*Funcao distribuicao de prob. (p) dist. acumulada*/
	proc iml;
	start pgumbel(q,mu,sigma);
		cdf= exp(-exp((-q + mu) / sigma));;
		return(cdf);
	finish;
	q = pgumbel(2,1.5,2); 
	print(q);
	quit;
	/*Funcao quantil (q)*/
	proc iml;
	start qgumbel(p,mu,sigma);
		qf=mu - sigma*log(-log(p));
		return(qf);
	finish;
	p = qgumbel(0.5,1.5,2); 
	print(p);
	quit;

		/*Funcao para valores aleatorios*/
		proc iml;
		start rgumbel(n,mu,sigma);
			aux	= j(n,1);                
			call randgen(aux, "uniform");
			ran	= -log(-log(aux))*sigma+mu; 		
			return(ran);
		finish;
		call randseed(666);
		random_gumbel = rgumbel(10,1.5,2); 
		print(random_gumbel);
		quit;

*----------------------------------------------------------------------------------------;
*MLE Gumbel no SAS;
proc iml;
start MLEgumbel(param) global(x);
   	mu = param[1];
   	sigma = param[2];
   	n	= nrow(x);
	mle = -n#log(sigma) + sum(log(exp((-x+mu)/sigma - exp((-x+mu)/sigma))));
   	return (mle);
finish;

*----------------------------------------------------------------------------------------------;
/*Gerar valores aleatorios*/
start rgumbel(n,mu,sigma);
	aux	= j(n,1);                
	call randgen(aux,"uniform"); 
	ran	= -log(-log(aux))*sigma+mu; 		
	return(ran);
finish;
call randseed(666);
x = rgumbel(10000,1.5,2); 
print(x);
sup = {.   0,  /* Li -inf < mu; 0 < sigma */
       .   .}; /* Ls mu < inf; sigma < inf */
iniciais = {1.3  2};	/*Valores iniciais*/
opt = {1,2};   
call nlpnra(it, result, "MLEgumbel", iniciais, opt, sup);
print it result iniciais/*Mostrar Resultados*/;
quit;
/*Outra dorma de estimação*/
*Funcao para valores aleatorios;
data dist4(drop=mu sigma u i);
do i= 1 to 5000;
call streaminit(666);/*Semente*/
	mu=0.5;
	sigma=2;
	u=rand("uniform");
	rgumbel=-log(-log(u))*sigma+mu; 
  output;
end;
run;
proc nlp data=dist4 tech=nra;
	max l;
	parms mu=0.5,sigma=2;
	x=rgumbel;
	aux=(-x + mu) / sigma;
	l = -log(sigma)+(log(exp(aux-exp(aux))));
run;

*--------------------------------------------------------------------------------------------------;
proc delete data=_all_;
run;
/*SIMULAÇÃO*/
%let mu=1;
%let sigma=1.5;
%macro simul(stop=,step=,nsimulacao=);
 %do m = 10 %to &stop %by &step;	
data simul;
call streaminit(1212);
	do b = 1 to &nsimulacao;
		do i = 1 to &m;
		u=rand("uniform");
		x=&mu-log(-log(u))*&sigma;
		output;
		end;
	end;
run;

*----------------------------------------------------------------------------------------------;
proc nlp data=simul tech=nra noprint outest=estimacao(keep=_type_ mu sigma where=(_type_="PARMS"));
	max ll;
	parms mu=1,sigma=1.5;
	bounds sigma > 0;
	aux=(-x+mu)/sigma;
  	ll= -log(sigma)+(log(exp(aux-exp(aux)))); *EXPRESSAO SEM SOMATORIO;
	by b;
run;
data resultados;
set estimacao;
	vicio_mu=mu-&mu;
	vicio_sigma=sigma-&sigma;
run;
proc means data=resultados noprint;
	var _numeric_;
		output out = v&m mean=/autoname;
run;
%end;
%mend simul;
%simul(stop=150,step=10,nsimulacao=400);
/*Juntandos os data-sets e apresentandos os resultados*/
data resultados_simulacao;
	set v10 v20 v30 v40 v50 v60 v70 v80 v90 v100 v110 v120 v130 v140 v150;
	n+10;
		if last.mu_mean then n=10;
run;

*----------------------------------------------------------------------------------------------;
%macro graph(var=,nome=);
ods html style=harvest; 
options nodate nonumber;
options orientation=portrait;
ods pdf file="D:\Estatística\ESTATÍSTICA_COMPUTACIONAL_II\
MLE\Gumbel_Latex\&nome..pdf" startpage= yes;
axis1 minor=none order=(0 to 150 by 20);
symbol1 i=j v=none l=1 c=blue w=3;
	proc gplot data=Resultados_simulacao;
		plot &var*n/overlay haxis=axis1 noframe autovref;
	run;
		ODS pdf close;
		ods close;
	%mend graph;
%graph(var=vicio_mu_Mean,nome=simul_mu);
%graph(var=vicio_sigma_Mean,nome=simul_sigma);
