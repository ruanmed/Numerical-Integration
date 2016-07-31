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
int Questao_1();
int Questao_3();

int main(){
	int seletor;

	cout<<"Questao 1: Metodo dos Trapezios"<<endl;
	cout<<"Questao 2: Metodo de Simpson"<<endl;
	cout<<"Questao 3: Quadratura de Gauss-Hermite"<<endl;
	cout<<"Digite a questao que deseja utilizar: ";
	cin>>seletor;
	pauseclear;

	switch(seletor)
	{
		case 1:
			Questao_1();
			pauseclear;
			break;
		case 3:
			Questao_3();
			pauseclear;
			break;

		default:
			cout<<"Opção invalida"<<endl;
			pauseclear;
			break;
	}
}

int Questao_1()
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

		integral.solveIntegration(op);
		integral.showSolution();

		return 0;
}
int Questao_3()
{
	GaussHermite integral;
	string ler;//Ler a função sem e^(-x^2) e deixar isso explicito
	cout<<"Integracao numerica usando a quadratura de Gauss-Hermite"<<endl;
	cout<<"Digite uma função que multiplicada por e^(-x^2) no intervalo de menos infinito a mais infinito tem  integral convergente"<<endl;
	cout << "f(x) = ";
	getline(cin, ler);
	if (cin.fail())
		return 0;
	integral.readFunction(ler);





	integral.solveIntegration(true);
	integral.showSolution();
	pause;


}
