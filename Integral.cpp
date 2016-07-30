/*
 * 	Integral.cpp
 *
 *  Created on: 28 de jul de 2016
 *      Author: Ruan
 */
#include "Integral.hpp"
double E=1e-7;
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
double Integral::getError(){
	return error;
}

void Integral::readFunction(string &f){
	sfunction = f;
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
	file << setprecision(6)<< "A solucao da integral de f(x) = " << sfunction <<
			" , no intervalo de " << getBottomLimit() <<
			" até " << setprecision(6) << getUpperLimit() << " eh: " << setprecision(6) <<
			getResult() << endl<<"O erro associado é de: "<<getError()<<endl;
	cout << setprecision(6)<< "A solucao da integral de f(x) = " << sfunction <<
			" , no intervalo de " << getBottomLimit() <<
			" até " << setprecision(6) << getUpperLimit() << " eh: " << setprecision(6) <<
			getResult() << endl<<"O erro associado é de: "<<getError()<<endl;
}
void	GaussHermite::generateHermitePolinom(int order)
{
	for(int c=0;c<order;c++)
	{
		int *aux = new int[(c+2)];

		for(int c1=0;c1<(c+2);c1++)
			aux[c1]=0;

		for(int c1=0;c1<ordemPoli;c1++)
		{
			aux[c1+1]+=2*poli[c1];
		}
		for(int c2=0;c2<ordemPoli2;c2++)
		{
			aux[c2]-=2*(ordemPoli-1)*poli2[c2];
		}

		free(poli2);
		poli2=poli;
		poli=aux;
		ordemPoli2=ordemPoli;
		ordemPoli++;

		/*//Debug
		for(int c1=(c+1);c1>-1;c1--)
			cout<<aux[c1]<<" ";
		cout<<endl;

		getchar();
		*/
	}

}
double  GaussHermite::getHermitePolinom(double x)
{
	double val;
	//usar parenteses encaixados...

	return val;
}

GaussHermite::GaussHermite()
{
	poli = new int[1];
	poli2 = new int[1];
	poli[0]=1;
	poli2[0]=0;
	ordemPoli=1;
	ordemPoli2=1;
}
GaussHermite::~GaussHermite()
{
	if(poli)
		free(poli);
	if(poli2)
		free(poli2);
}
void 	GaussHermite::solveIntegration(const bool &saveLog)
{

}
