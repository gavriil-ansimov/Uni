#ifndef MNPS
#define MNPS
#include "gauss.h"
#include<iostream>
#include <vector>
using namespace std;

vector <double> mnps(double eps, vector <double> x) {
	double step = 0, k = 0, stop = 1;
	vector <double> q{ 0, 0, 0 };
	while (stop > eps) {
		for (int i = 0; i < 3; i++) {
			q = { 0, 0, 0 };
			q[i] = 1;
			step = (-1) * (strvec(q, multAplusB(x)) / strvec(q, multA(q)));
			x[i] = x[i] + step * q[i];
			k++;
		}
		//cout << fun(x);
		stop = abs(sqrt(strvec(multAplusB(x), multAplusB(x)))) / d;
	}
	cout << "iter mnps: " << k << '\n';
	//cout << fun(x);
	return x;
}













#endif
