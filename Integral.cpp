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

double	Integral::getFunction(const double &val){
	double t[] = { 0 };
	t[0] = val;
	return function.Eval(t);
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

void Integral::showSolution(){
	file << fixed << setprecision(7) << "A solucao da integral de f(x) = " << sfunction <<
			" , no intervalo de " << getBottomLimit() <<
			" até " << setprecision(7) << getUpperLimit() << " eh: " << setprecision(7) <<
			getResult() << endl<<"O erro associado é de: "<< setprecision(7) << getError()<<endl;
	cout << fixed << setprecision(7)<< "A solucao da integral de f(x) = " << sfunction <<
			" , no intervalo de " << getBottomLimit() <<
			" até " << setprecision(7) << getUpperLimit() << " eh: " << setprecision(7) <<
			getResult() << endl<<"O erro associado é de: "<< setprecision(7) << getError()<<endl;
}

void Trapezium::solveIntegration(const bool &saveLog)
{
	if (getUpperLimit() == getBottomLimit())
			setResult(0);
	else
		partitionateInterval(saveLog);
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
		file << setprecision(6) << "\n\tDivisão em " << i << " intervalos de " << dx <<
				" -- valor da integral = " << setprecision(6) << A << "\t" << endl;

	return A;
}

void 	Boole::solveIntegration(const bool &saveLog){
	if (getUpperLimit() == getBottomLimit())
			setResult(0);
	else {
		double aux;
		double n = 1;
		double divisoes = 4;	//Número de subintervalos
		double dx=getUpperLimit()-getBottomLimit();
		double c;
		do {
			aux = getResult();
			divisoes *=n;
			dx=(getUpperLimit()-getBottomLimit())/divisoes;

			setResult(7*getFunction(getUpperLimit())+7*getFunction(getBottomLimit()));

			for (c = 1;c < divisoes;c+=2)
				setResult(getResult()+32*getFunction((getBottomLimit()+(dx*c))));
			for (c = 2;c < divisoes;c+=4)
				setResult(getResult()+12*getFunction((getBottomLimit()+(dx*c))));
			for (c = 4;c < divisoes;c+=4)
				setResult(getResult()+14*getFunction((getBottomLimit()+(dx*c))));

			setResult((2.0/45.0)*dx*getResult());
			n++;
			if (saveLog) {
				file << endl << fixed << setprecision(0) <<
						"Divisão parcial em " << divisoes << " intervalos de " << setprecision(7) << dx <<
						"-- Valor da integral = " << setprecision(7) << getResult() << endl;
			}
		} while (fabs(getResult()-aux) > E);
		setError(E);

		if (saveLog){
			file << endl << fixed << setprecision(0) <<
					"Divisão em " << divisoes << " intervalos de " << setprecision(7) <<dx << '.' << endl;
			cout << endl << fixed << setprecision(0) <<
					"Divisão em " << divisoes << " intervalos de " << setprecision(7) << dx << '.' << endl;
		}
	}
}


void	GaussHermite::generateHermitePolinoms(int order)
{
	for(int c=getOrderPoli();c<order;c++)
	{
		int *aux = new int[(c+2)];

		for(int c1=0;c1<(c+2);c1++)
			aux[c1]=0;

		//Debug
		//cout<<"orderPoli: "<<getOrderPoli()<<" orderPoli2: "<<getOrderPoli2()<<endl;

		for(int c1=0;c1<=getOrderPoli();c1++)
		{
			aux[c1+1]+=2*poli[c1];
		}
		for(int c2=0;c2<=getOrderPoli2();c2++)
		{
			aux[c2]-=2*(orderPoli)*poli2[c2];
		}

		free(poli2);
		poli2=poli;
		poli=aux;
		setOrderPoli2(getOrderPoli());
		setOrderPoli(getOrderPoli()+1);

		/*//Debug
		for(int c1=(c+1);c1>-1;c1--)
			cout<<aux[c1]<<" ";
		cout<<endl;

		getchar();
		//*/
	}

}
double  GaussHermite::getHermitePolinom(bool Switch,double x)
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
	double kick;

	allocRoots(getOrderPoli());
	setNumRoots(0);

	for(i=0;i<getOrderPoli();i++)//Pelo pdf de edison seria getOrderPoli()-1, porém não faz sentido...
	{
		if(i==0)
			kick=sqrt(2.0*getOrderPoli()+1 - 1.85575*pow(2.0*getOrderPoli()+1,-0.16667));
		else
		if(i==1)
			kick-=1.14*pow(getOrderPoli(),0.426)/kick;
		else
		if(i==2)
			kick=1.86*kick-0.86*roots[0];
		else
		if(i==3)
			kick=1.91*kick-0.91*roots[1];
		else
			kick=2.0*kick-roots[i-2];

		roots[i]=newtonMethodPolinoms(kick);
		setNumRoots(getNumRoots()+1);

	}

}
GaussHermite::GaussHermite()
{
	poli = new int[1];
	poli2 = new int[1];
	roots = NULL;
	poli[0]=1;
	poli2[0]=0;
	orderPoli=0;
	orderPoli2=0;
	numRoots=0;
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
void 	GaussHermite::solveIntegration(const bool &saveLog){
	int c=2,fat=1;//O primeiro é pulado por ser trivial
	double newResult,oldResult=0,weight;
	newResult=getFunction(&oldResult)* 1.77245385090551;// 1.77245385090551 é o peso para P(1), que é constante e a raiz de P(1) é 0

	while(fabs(newResult-oldResult)>E && c>14)//Fatorial estoura lindamente!
	{
		oldResult=newResult;
		newResult=0;
		generateHermitePolinoms(c);
		printHermitePolinoms();
		generateRoots();
		printPolinomsRoots();
		fat*=c;
		for(int c1=0,weight = ((2<<(c+1))*fat*1.77245385090551) ;c1<c; c1++)
			newResult+=((weight*getFunction(roots+c1))/(getHermitePolinomDerivative(roots[c1])*getHermitePolinomDerivative(roots[c1])));

		pause;
		c++;


	}
}
