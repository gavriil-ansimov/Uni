#include "gauss.h"
#include "mngs.h"
#include "mnps.h"
#include <iostream>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

string f = "2x^2+(3+0.1*1)y^2+(4+0.1*1)z^2+xy-yz+xz+x-2y+3z+1";
double eps = 10e-6;

int main() {
	vector<double> min = gauss();
	vector<double> x{ 0, 0, 0 };
	setlocale(LC_ALL, "Russian");
	cout << "Номер студента в группе: 1" << '\n';
	cout << "Функция:\nf=" << f;
	cout << "\nМатрица А:" << '\n';
	for (int i = 0; i < 3; i++) {
		cout << "(\t";
		for (int j = 0; j < 3; j++) {
			cout << matrix[i][j] << '\t';
		}
		cout << ')' << '\n';
	}
	cout << "Вектор В:" << '\n';
	cout << "(\t";
	for (int i = 0; i < 3; i++) {
		cout << -matrix[i][3] << '\t';
	}
	cout << ")\n";
	cout << "Начальный вектор Х:" << '\n';
	cout << "(\t";
	for (int i = 0; i < 3; i++) {
		cout << x[i] << '\t';
	}
	cout << ")\n";
	vector <double> gr = mngs(eps, x);
	vector <double> kor = mnps(eps, x);
	cout << "МНГС:\nПолученный вектор Х:\n";
	cout << "(\t";
	for (int i = 0; i < 3; i++) {
		cout << gr[i] << '\t';
	}
	cout << ")\n";
	cout << "Значение функции:\n" << fun(gr);
	cout << "\nАбсолютная погрешность: " << abs(fun(min) - fun(gr));
	cout << "\nМНПС:\nПолученный вектор Х:\n";
	cout << "(\t";
	for (int i = 0; i < 3; i++) {
		cout << kor[i] << '\t';
	}
	cout << ")\n";
	;
	cout << "Значение функции:\n" << fun(kor);
	cout << "\nАбсолютная погрешность: " << abs(fun(min) - fun(kor));
	cout << '\n' << min[0] << ' ' << min[1] << ' ' << min[2];
}