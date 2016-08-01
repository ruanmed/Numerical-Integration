/*
 * Universidade Federal do Vale do São Francisco
 * Trabalho de Cálculo Numérico - Integração Numérica
 * Professor: Edson Leite Araújo
 * Alunos:
 */
#include "Integral.hpp"
#define pause cout<<"Aperte enter para continuar."<<endl; getchar();
#define pauseclear   pause system("clear || CLS");

using namespace std;
void Problema_1();
void Problema_2();
void Problema_3();

int main(){
	int seletor;

	cout << "Problema 1: Integração utilizando a Regra do Trapézio" << endl;
	cout << "Problema 2: Integração utilizando a Regra de Boole (Regra de Simpson generalizada e exata para poliônimos de graur menor ou igual a 5)" << endl;
	cout << "Problema 3: Integração utilizando a Quadratura de Gauss-Hermite" << endl;
	cout << "Digite o número do problema que deseja resolver: " ;

	cin >> seletor;
	pauseclear;

	switch(seletor){
		case 1:
			Problema_1();
			break;
		case 2:
			Problema_2();
			break;
		case 3:
			Problema_3();
			break;
		default:
			cout<<"Opção inválida."<<endl;
			break;
	}
	pauseclear;

	return 0;
}

void Problema_1()
{
		double a,b;
		bool op;
		string ler;
		Trapezium integral;

		cout << "Integração Numérica usando Regra do Trapézio" << endl;
		cout << "Digite a função para integração abaixo." << endl;
		cout << "f(x) = ";
		getline(cin, ler);
		if (cin.fail())
			return;
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

		integral.solveIntegration(op);
		integral.showSolution();
		return;
}

void Problema_2(){
	double a,b;
	bool op;
	string arquivo;
	string ler;
	Boole integral;

	cout << "Integração Numérica usando Regra de Boole (Generalização da Regra de Simpson exata para poliônimos de graur menor ou igual a 5)" << endl;
	cout << "Digite a função para integração abaixo." << endl;
	cout << "f(x) = ";
	getline(cin, ler);
	if (cin.fail())
		return;
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
	cout << "Isso pode demorar um pouco para funções muito irregurales em intervalos grandes." << endl;

	integral.solveIntegration(op);
	integral.showSolution();

	return;
}

void Problema_3()
{
	GaussHermite integral;
	string ler;//Ler a função sem e^(-x^2) e deixar isso explicito//concaternar com e^(-x^2) depois..
	bool op;
	cout<<"Integração Numérica usando a quadratura de Gauss-Hermite"<<endl;
	cout<<"Digite uma função que multiplicada por e^(-x^2) no intervalo de menos infinito a mais infinito tem integral convergente."<<endl;
	cout << "f(x) = ";

	getline(cin, ler);
	if (cin.fail())
		return;

	integral.readFunction(ler);

	cout << "Deseja salvar um registro da expansão dos polinômios do programa?" <<"(Isso afetará bastante o desempenho do programa)" << endl <<"[1] Sim [0] Não" << endl;
	cin >> op;
	getchar();
	if (op)
		cout << endl << "O registro da expansão dos polinômios do programa será salvo no arquivo log.txt"<<endl;
	pauseclear;

	cout << "Processando..." << endl;
	cout << "Isso pode demorar um pouco se você selecionou para salvar o registro das expansões." << endl;

	integral.solveIntegration(op);
	pauseclear;
	integral.showSolution();


	return;
}
