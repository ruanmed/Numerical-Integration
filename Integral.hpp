/*
 * Integral.hpp
 *
 *  Created on: 28 de jul de 2016
 *      Author: Ruan
 */

#ifndef INTEGRAL_HPP_
#define INTEGRAL_HPP_

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "fparser.hh"
using namespace std;



class Integral{
	private:
		FunctionParser function;
		string 	sfunction;
		double 	bottomLimit;
		double 	upperLimit;
		double 	result;
		double 	error;
	public:
		ofstream file;
		Integral();
		virtual	~Integral();
		void 	setBottomLimit(const double &newBottomLimit);
		void 	setUpperLimit(const double &newUpperLimit);
		void 	setLimits(const double &a, const double &b);
		void 	setResult(const double &newResult);
		void 	setError(const double &newError);
		double 	getBottomLimit();
		double 	getUpperLimit();
		double 	getResult();
		double 	getError();
		void	setFunction(string &f);
		string	getFunction();
		void 	readFunction(string &f);
		double	getFunction(double * val);
		virtual void 	solveIntegration(const bool &saveLog) = 0;
		void 	showSolution();
};
class Trapezium: public Integral
{
	private:
		void  	partitionateInterval(const bool &saveLog = 0);
		double 	partitionatedIntervalIntegralSolution(double begin,double end,const bool &saveLog = 0);
	public:
		void 	solveIntegration(const bool &saveLog);
};
class GaussHermite: public Integral
{
	private:
		int	*poli;
		int *poli2;
		double *roots;
		double R;
		int orderPoli;
		int orderPoli2;
		int numRoots;

		void	generateHermitePolinoms(int order);
		double  getHermitePolinom(bool Switch,double x);
		double  getHermitePolinomDerivative(double x);
		double	newtonMethodPolinoms(double kick);
		void	allocRoots(int num);
		void	generateRoots();
	public:
		GaussHermite();
		~GaussHermite();
		int	getOrderPoli();
		int	getOrderPoli2();
		int	getNumRoots();
		void	setOrderPoli(int order);
		void	setOrderPoli2(int order);
		void	setNumRoots(int num);
		void	printPolinomsRoots();
		void	printHermitePolinoms();
		void 	solveIntegration(const bool &saveLog);
		void	showSolution();
};



#endif /* INTEGRAL_HPP_ */
