/*
 * 	Integral.cpp
 *
 *  Created on: 28 de jul de 2016
 *      Author: Ruan
 */
#include "Integral.hpp"

void  Integral::partitionateInterval(const bool &saveLog){
	double amp=(upperLimit-bottomLimit),passo;
	int c;
	if(amp<10)
		passo=amp/100000;
	else
	if(amp<100)
		passo=amp/1000000;
	else
	if(amp<1000)
		passo=amp/10000000;
	else
		passo=0.1;

	for(c=0;(c+1)*passo<upperLimit;c++)
	{
		setResult(getResult()+partitionatedIntervalIntegralSolution(c*passo,(c+1)*passo));
		if (saveLog)
			file << setprecision(6) << "Intervalo: " << c*passo << " a " << (c+1)*passo << endl <<
				setprecision(6) << "Integral ate aqui: " << getResult() << endl;
	}
	c--;
	setResult(getResult()+partitionatedIntervalIntegralSolution(c*passo,upperLimit));
	if (saveLog)
		file << setprecision(6) << "Intervalo: " << c*passo << " a " << upperLimit << endl <<
				setprecision(6) << "Integral ate aqui: " << getResult() << endl;
}

double Integral::partitionatedIntervalIntegralSolution(double begin,double end,const bool &saveLog){
	double A,Aa,dx=end-begin,E=1e-7;
	long long int c,i=1;//c - Contador e i = n° de intervalos
	double temp[] = { 0 };
	dx=end-begin;
	A=((function.Eval(&begin)+function.Eval(&end))*dx)/2;
	do{
		for(c=1,Aa=A,dx/=2,A=0,i<<=1; c<i && dx > E; c+=2){
			temp[0] = (c*dx)+begin;
			A+=function.Eval(temp);
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

void Integral::solveWithTrapeziumRule(const bool &saveLog){
	if (getUpperLimit() == getBottomLimit())
		setResult(0);
	else if((getUpperLimit()-getBottomLimit())>1)
		partitionateInterval(saveLog);
	else
		setResult(partitionatedIntervalIntegralSolution(getBottomLimit(),getUpperLimit(),saveLog));

}

void Integral::solveIntegration(){

}

void Integral::showSolution(){
	file << setprecision(6)<< "A solucao da integral de f(x) = " << sfunction <<
			" no intervalo de " << getBottomLimit() <<
			" até " << setprecision(6) << getUpperLimit() << " eh: " << setprecision(6) <<
			getResult() << endl;
	cout << setprecision(6)<< "A solucao da integral de f(x) = " << sfunction <<
			" no intervalo de " << getBottomLimit() <<
			" até " << setprecision(6) << getUpperLimit() << " eh: " << setprecision(6) <<
			getResult() << endl;
}
