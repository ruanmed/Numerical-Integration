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
		ofstream file;
		string 	sfunction;
		double 	bottomLimit;
		double 	upperLimit;
		double 	result;
		double 	error;
		void  	partitionateInterval(const bool &saveLog = 0);
		double 	partitionatedIntervalIntegralSolution(double begin,double end,const bool &saveLog = 0);
		void	generateHermitePolinom();
		double  getHermitePolinom(double x);

	public:
		Integral();
		~Integral();
		void 	setBottomLimit(const double &newBottomLimit);
		void 	setUpperLimit(const double &newUpperLimit);
		void 	setLimits(const double &a, const double &b);
		void 	setResult(const double &newResult);
		void 	setError(const double &newError);
		double 	getBottomLimit();
		double 	getUpperLimit();
		double 	getResult();
		double 	getError();
		void 	readFunction(string &f);
		void	solveWithGaussHermite(const bool &saveLog = 0);
		void 	solveWithTrapeziumRule(const bool &saveLog = 0);
		void 	solveIntegration();
		void 	showSolution();
};



#endif /* INTEGRAL_HPP_ */
