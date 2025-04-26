#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<string>
using namespace std;

double fun(double x) {
	return -sin(x) + x - 0.25;
}

double lag(int n, vector<double>& x, double x1) {
	double res = 0, l = 1;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i != j) {
				l *= (x1 - x[j]) / (x[i] - x[j]);
			}
		}
		res += fun(x[i]) * l;
		l = 1;
	}
	return res;
}

void x_eq(double a, double b, int n, vector<double>& x) {
	double step = (b - a) / (n - 1);
	//cout << step << ' ';
	for (int i = 0; i < n; i++) {
		x[i] = a + i * step;
	}
}

void x_opt(double a, double b, int n, vector<double>& x) {
	double x_i = 0;
	for (int i = 0; i != n; i++) {
		x_i = 0.5 * ((b - a) * cos((2 * i + 1) * acos(-1.0) / (2 * n + 2)) + b + a);
		x[i] = x_i;
	}
}

double Newton(int n, vector<double>& x, vector<double>& y, double x1) {
	double res = y[0];
	double F = 0;
	for (int i = 1; i < n; ++i) {
		for (int j = 0; j <= i; ++j) {
			double den = 1;
			for (int k = 0; k <= i; ++k)
				if (k != j)
					den *= (x[j] - x[k]);
			F += y[j] / den;
		}
		for (int k = 0; k < i; ++k)
			F *= (x1 - x[k]);
		res += F;
		F = 0;
	}
	return res;
}

int main() {
	double a = 1, b = 10, r_Ln = 0, r_Lopt = 0, r_n = 0, r_nop = 0;
	double t_fun = 0, t_lag = 0, t_lag1 = 0;
	int n = 100;
	double max = 0;
	cout << "Lagrange:\n";
	cout << "N nodes: " << "Max err: \t" << "Max err opt:\n";
	for (int n = 5; n <= 10; n += 1) {
		vector<double> x(n);
		vector<double> x1(n);
		ofstream fout("r_eq_L" + to_string(n) + ".txt");
		ofstream fout1("r_opt_L" + to_string(n) + ".txt");
		x_eq(a, b, n, x);
		x_opt(a, b, n, x1);
		for (double i = 1; i < 10; i += 0.1) {
			t_fun = fun(i);
			t_lag = lag(n, x, i);
			t_lag1 = lag(n, x1, i);

			if (abs(t_fun - t_lag) > r_Ln) {
				r_Ln = abs(t_fun - t_lag);
			}
			if (abs(t_fun - t_lag1) > r_Lopt) {
				r_Lopt = abs(t_fun - t_lag1);
			}
			fout << i << '\t' << abs(t_fun - t_lag) << '\n';
			fout1 << i << '\t' << abs(t_fun - t_lag1) << '\n';
		}
		fout.close();
		fout1.close();

		cout << n << '\t' << r_Ln << '\t' << r_Lopt << '\n';
		r_Ln = 0;
		r_Lopt = 0;
	}
	cout << "Newton:\n";
	cout << "N nodes: " << "Max err: \t" << "Max err opt:\n";
	for (int n = 10; n != 100; n += 10) {
		vector<double> x(n);
		vector<double> x1(n);
		x_eq(a, b, n, x);
		x_opt(a, b, n, x1);
		vector<double> y;
		vector<double> y1;
		for (auto el : x) {
			y.push_back(fun(el));
		}
		for (auto el : x1) {
			y1.push_back(fun(el));
		}
		for (double i = 1.5; i < 9.5; i += 0.1) {
			t_fun = fun(i);
			t_lag = Newton(n, x, y, i);
			t_lag1 = Newton(n, x1, y1, i);
			if (abs(t_fun - t_lag) > r_n) {
				r_n = abs(t_fun - t_lag);
			}
			if (abs(t_fun - t_lag1) > r_nop) {
				r_nop = abs(t_fun - t_lag1);
			}
		}
		cout << n << "\t" << r_n << "\t" << r_nop << '\n';
		r_n = 0;
		r_nop = 0;
	}
	return 0;
}