/*
 * 	Integral.cpp
 *
 *  Created on: 28 de jul de 2016
 *      Author: Ruan
 */
#include "Integral.hpp"
double E=1e-7;
#define pause cout<<"Aperte enter para continuar."<<endl; getchar();
#define pauseclear   pause system("clear || CLS");
double	Integral::getFunction(double * val)
{
	return function.Eval(val);
}
void  Trapezium::partitionateInterval(const bool &saveLog){
	double amp=(getUpperLimit()-getBottomLimit()),passo;
	int c,divisoes;
	if(amp<=1)
		divisoes=1;
	else
	if(amp<=10)
		divisoes=100000;
	else
	if(amp<=100)
		divisoes=1000000;
	else
	if(amp<=1000)
		divisoes=10000000;
	else
		divisoes=(int)(amp/0.001);

	passo=amp/divisoes;
	setError(divisoes*E);//Tratanto acumulo dos erros pela utilização da soma de integrais

	for(c=0;getBottomLimit()+(c+1)*passo<getUpperLimit();c++)
	{
		setResult(getResult()+partitionatedIntervalIntegralSolution(getBottomLimit()+c*passo,getBottomLimit()+(c+1)*passo));
		if (saveLog)
			file << setprecision(6) << "Intervalo: " << c*passo << " a " << (c+1)*passo << endl <<
				setprecision(6) << "Integral ate aqui: " << getResult() << endl;
	}

	if(c>0)
		c--;
	setResult(getResult()+partitionatedIntervalIntegralSolution(getBottomLimit()+c*passo,getUpperLimit()));
	if (saveLog)
		file<< setprecision(6) << "Intervalo: " << c*passo << " a " << getUpperLimit() << endl <<setprecision(6) << "Integral ate aqui: " << getResult() << endl;
}

double Trapezium::partitionatedIntervalIntegralSolution(double begin,double end,const bool &saveLog){
	double A,Aa,dx=end-begin;
	long long int c,i=1;//c - Contador e i = n° de intervalos
	double temp[] = { 0 };
	A=((getFunction(&begin)+getFunction(&end))*dx)/2;
	do{
		for(c=1,Aa=A,dx/=2,A=0,i<<=1; c<i && dx > E; c+=2){
			temp[0] = (c*dx)+begin;
			A+=getFunction(temp);
		}

		A*=dx;
		A+=(Aa/2);
	}while(fabs(A-Aa)>E);

	if (saveLog)
	file << setprecision(6) << "\n\tDivisao em " << i << " intervalos de " << dx <<
			" -- valor da integral = " << setprecision(6) << A << "\t" << endl;

	return A;
}

Integral::Integral(){
	setLimits(0,0);
	setResult(0);
	setError(0);
	file.open("log.txt");
}

Integral::~Integral(){
	file.close();
}

void Integral::setBottomLimit(const double &newBottomLimit){
	bottomLimit = newBottomLimit;
}

void Integral::setUpperLimit(const double &newUpperLimit){
	upperLimit = newUpperLimit;
}

void Integral::setLimits(const double &a, const double &b){
	if (a < b){
		setBottomLimit(a);
		setUpperLimit(b);
	}
	else {
		setBottomLimit(b);
		setUpperLimit(a);
	}
}

void Integral::setResult(const double &newResult){
	result = newResult;
}
void Integral::setError(const double &newError){
	error = newError;
}

double Integral::getBottomLimit(){
	return bottomLimit;
}

double Integral::getUpperLimit(){
	return upperLimit;
}

double Integral::getResult(){
	return result;
}
double	Integral::getError(){
	return error;
}
void	Integral::setFunction(string &f)
{
	sfunction = f;
}
string	Integral::getFunction()
{
	return sfunction;
}
void Integral::readFunction(string &f){
	setFunction(f);
	while (true){
		int res = function.Parse(f, "x");
		if (res < 0)
			break;
		std::cout << std::string(res+7, ' ') << "^\n"
					<< function.ErrorMsg() << "\n\n";
	}
}

void Trapezium::solveIntegration(const bool &saveLog)
{
	if (getUpperLimit() == getBottomLimit())
			setResult(0);
		else
			partitionateInterval(saveLog);
}

void Integral::showSolution(){
	file << setprecision(6)<< "A solucao da integral de f(x) = " << getFunction() <<
			" , no intervalo de " << getBottomLimit() <<
			" até " << setprecision(6) << getUpperLimit() << " eh: " << setprecision(6) <<
			getResult() << endl<<"O erro associado é de: "<<getError()<<endl;
	cout << setprecision(6)<< "A solucao da integral de f(x) = " << getFunction() <<
			" , no intervalo de " << getBottomLimit() <<
			" até " << setprecision(6) << getUpperLimit() << " eh: " << setprecision(6) <<
			getResult() << endl<<"O erro associado é de: "<<getError()<<endl;
}
void	GaussHermite::generateHermitePolinoms(int order)
{
	for(int c=getOrderPoli();c<order;c++)
	{
		int *aux = new int[(c+2)];
		double maior;

		for(int c1=0;c1<(c+2);c1++)
			aux[c1]=0;

		for(int c1=0;c1<=getOrderPoli();c1++)
			aux[c1+1]+=2*poli[c1];

		for(int c2=0;c2<=getOrderPoli2();c2++)
			aux[c2]-=2*(orderPoli)*poli2[c2];

		free(poli2);
		poli2=poli;
		poli=aux;
		setOrderPoli2(getOrderPoli());
		setOrderPoli(getOrderPoli()+1);

		for(int c1=0;c1<order-1;c1++)
			if(!c1 || fabs(poli[c1])>maior)
				maior=fabs(poli[c1]);

		R = (maior/poli[order])+1;


	}

}
double  GaussHermite::getHermitePolinom(bool Switch,double x)//Metodo dos parenteses encaixados de horner
{	int c;
	double value;

	if(Switch)
	{
		for(c=getOrderPoli(),value=poli[c];c>0;c--)
		{
			value*=x;
			value+=poli[c-1];
		}
	}
	else
	{
		for(c=getOrderPoli2(),value=poli2[c];c>0;c--)
		{
			value*=x;
			value+=poli2[c-1];
		}
	}
	return value;
}
double  GaussHermite::getHermitePolinomDerivative(double x)
{
	return (2*orderPoli*getHermitePolinom(false,x));
}
double	GaussHermite::newtonMethodPolinoms(double kick)
{	double R=kick+1,E1=10e-10,E2=10e-9;

	while((fabs(getHermitePolinom(true,kick))>E1) && (fabs(R-kick)>E2))
	{
			R=kick;
			kick-=(getHermitePolinom(true,kick)/getHermitePolinomDerivative(kick));
	}
	return  kick;

}
void	GaussHermite::allocRoots(int num)
{
	if(num>0)
	{
		if(roots)
			delete []roots;

		roots = new double[num];
		if(!roots)
		{
			cout<<"Erro ao alocar"<<endl;
		}

	}
	else
	{
		cout<<"Valor invalido para alocar"<<endl;
	}

}
void	GaussHermite::generateRoots()
{
	int i;
	double kick,passo=0.1;

	allocRoots(getOrderPoli());
	setNumRoots(0);
	if(getOrderPoli()&1)//Order é impar?
	{
		roots[0]=0;
		setNumRoots(getNumRoots()+1);
	}

	for(i=1;i*passo<=R && getNumRoots()<=getOrderPoli()/2;i++)
	{
		if((getHermitePolinom(true,i*passo)*getHermitePolinom(true,(i+1)*passo)) < 0)
		{	kick = (((2*i+1)*passo)/2);
			roots[getNumRoots()] = newtonMethodPolinoms(kick);
			setNumRoots(getNumRoots()+1);
		}

	}

	i = (getOrderPoli()&1)?1:0;

	for(;getNumRoots()<=getOrderPoli();i++)//espelha as raizes pela origem
	{
		roots[getNumRoots()]=(-roots[i]);
		setNumRoots(getNumRoots()+1);
	}

	setNumRoots(getNumRoots()-1);


}
GaussHermite::GaussHermite()
{
	poli = new int[1];
	poli2 = new int[1];
	roots = NULL;
	poli[0]=1;
	poli2[0]=0;
	R=0;
	orderPoli=0;
	orderPoli2=0;
	numRoots=0;
	file<<"Solução pela Quadratura de Gauss-Hermite:"<<endl;
}
GaussHermite::~GaussHermite()
{
	if(poli)
		free(poli);
	if(poli2)
		free(poli2);
}
int	GaussHermite::getOrderPoli()
{
	return orderPoli;
}
int	GaussHermite::getOrderPoli2()
{
	return orderPoli2;
}
int	GaussHermite::getNumRoots()
{
	return numRoots;
}
void	GaussHermite::setOrderPoli(int order)
{
	orderPoli=order;
}
void	GaussHermite::setOrderPoli2(int order)
{
	orderPoli2=order;
}
void	GaussHermite::setNumRoots(int num)
{
	numRoots = num;
}
void	GaussHermite::printPolinomsRoots()
{
	cout<<"Roots: ";
	for(int c=0; c<getNumRoots();c++)
		cout<<roots[c]<<" ";
	cout<<endl;
}
void	GaussHermite::printHermitePolinoms()
{
	int c;

	cout<<"P("<<orderPoli<<"): ";
	for(c=orderPoli; c>-1; c--)
		if(poli[c])
			cout<<poli[c]<<"*"<<"x^"<<c<<" ";

	cout<<endl;
	cout<<"P("<<orderPoli2<<"): ";
	for(c=orderPoli2; c>-1; c--)
			cout<<poli2[c]<<"*"<<"x^"<<c<<" ";

	cout<<endl;

}
void 	GaussHermite::solveIntegration(const bool &saveLog)
{	int c=2,fat=1;//O primeiro é pulado por ser trivial
	double newResult,oldResult=0,weight;
	newResult=getFunction(&oldResult)* 1.77245385090551;// 1.77245385090551 é o peso para P(1), que é constante e a raiz de P(1) é 0
	setError(E);

	while((fabs(newResult-oldResult)>E && c<10) || c==2)//Fatorial estoura lindamente!
	{
		if(saveLog)
			file<<"P("<<getOrderPoli()<<") gera a integral de valor:  "<<newResult<<endl;
		oldResult=newResult;
		newResult=0;
		generateHermitePolinoms(c);
		generateRoots();
		fat*=c;
		for(int c1=0,weight = ((2<<(c))*fat*1.77245385090551) ;c1<c; c1++)
			newResult+=((weight*getFunction(roots+c1))/(getHermitePolinomDerivative(roots[c1])*getHermitePolinomDerivative(roots[c1])));
		c++;
	}

	if(saveLog)
		file<<"P("<<getOrderPoli()<<") gera a integral de valor:  "<<newResult/2<<endl;
	setResult(newResult);
}
void	GaussHermite::showSolution()
{
	string aux="e^(-x^2)*";

	file << setprecision(6)<< "A solucao da integral de f(x) = " << aux+getFunction() <<
				" , em toda extensão da reta dos reais eh: " << setprecision(6) <<
				getResult() << endl<<"O erro associado é de: "<<getError()<<endl;
		cout << setprecision(6)<< "A solucao da integral de f(x) = " << aux+getFunction() <<
				" , em toda extensão da reta dos reais eh: " << setprecision(6) <<
								getResult() << endl<<"O erro associado é de: "<<getError()<<endl;
}
