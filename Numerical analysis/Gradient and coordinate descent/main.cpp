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
	cout << "����� �������� � ������: 1" << '\n';
	cout << "�������:\nf=" << f;
	cout << "\n������� �:" << '\n';
	for (int i = 0; i < 3; i++) {
		cout << "(\t";
		for (int j = 0; j < 3; j++) {
			cout << matrix[i][j] << '\t';
		}
		cout << ')' << '\n';
	}
	cout << "������ �:" << '\n';
	cout << "(\t";
	for (int i = 0; i < 3; i++) {
		cout << -matrix[i][3] << '\t';
	}
	cout << ")\n";
	cout << "��������� ������ �:" << '\n';
	cout << "(\t";
	for (int i = 0; i < 3; i++) {
		cout << x[i] << '\t';
	}
	cout << ")\n";
	vector <double> gr = mngs(eps, x);
	vector <double> kor = mnps(eps, x);
	cout << "����:\n���������� ������ �:\n";
	cout << "(\t";
	for (int i = 0; i < 3; i++) {
		cout << gr[i] << '\t';
	}
	cout << ")\n";
	cout << "�������� �������:\n" << fun(gr);
	cout << "\n���������� �����������: " << abs(fun(min) - fun(gr));
	cout << "\n����:\n���������� ������ �:\n";
	cout << "(\t";
	for (int i = 0; i < 3; i++) {
		cout << kor[i] << '\t';
	}
	cout << ")\n";
	;
	cout << "�������� �������:\n" << fun(kor);
	cout << "\n���������� �����������: " << abs(fun(min) - fun(kor));
	cout << '\n' << min[0] << ' ' << min[1] << ' ' << min[2];
}