#ifndef MNGS
#define MNGS
#include<iostream>
#include <vector>
using namespace std;

vector <double> mngs(double eps, vector<double> x) {
	double step = 0, x1 = 0, x2 = 0, x3 = 0, k = 0, stop = 1;
	while (stop > eps) {
		step = (-1) * (strvec(multAplusB(x), multAplusB(x)) / strvec(multAplusB(x), multA(multAplusB(x))));
		x1 = x[0] + step * gradx(x);
		x2 = x[1] + step * grady(x);
		x3 = x[2] + step * gradz(x);
		x[0] = x1;
		x[1] = x2;
		x[2] = x3;
		step = 0;
		stop = abs(sqrt(strvec(multAplusB(x), multAplusB(x)))) / d;
		//k++;
		//for (int i = 0; i < 3; i++) 
		//	cout << x[i] << ' ';
		// cout << '\n';
		//cout << "step" << step << "F = " << fun(x) << '\n';
	}
	cout << "delta: " << stop << '\n';
	return x;
}


#endif