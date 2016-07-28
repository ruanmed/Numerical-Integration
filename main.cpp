/*
 * Regra do Trapézio para cálculo numérico de integrais
 * Universidade Federal do Vale do São Francisco
 * Trabalho de Cálculo Numérico - Integração Numérica
 * Professor: Edson Leite Araújo
 */
#include "Integral.hpp"

#define pause cout<<"Aperte enter para continuar."<<endl; getchar();
#define pauseclear   pause system("clear || CLS");

using namespace std;

int main(){
	double a,b;
	bool op;
	string arquivo;
	string ler;
	Integral integral;
	
	cout << "Integração Numérica usando Regra do Trapézio" << endl;
	cout << "Digite a função para integração abaixo." << endl;
	cout << "f(x) = ";
	getline(cin, ler);
	if (cin.fail())
		return 0;
	integral.readFunction(ler);

	cout << "Digite os limites de integração [a,b] abaixo." << endl;
	cout << "a: ";
	cin >> a;
	cout << endl;
	cout << "b: ";
	cin >> b;
	cout << endl;
	cout << "Deseja salvar um registro das divisões do programa?" <<
			"(Isso afetará bastante o desempenho do programa)" << endl <<
			"[1] Sim [0] Não" << endl;
	cin >> op;
	if (op)
		cout << endl << "O registro de divisões do programa será salvo no arquivo log.txt"<<endl;
	pauseclear;
	
	integral.setLimits(a,b);

	cout << "Processando..." << endl;

	integral.solveWithTrapeziumRule(op);
	integral.showSolution();

	pauseclear;
	return 0;
}
