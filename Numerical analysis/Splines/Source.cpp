#include<iostream>
#include<cmath>
#include<vector>
#include<algorithm>
#include<fstream>
#include<string>
using namespace std;

double fun(double x) {
	return x - sin(x) - 0.25;
}

vector<double> gauss(vector<vector<double>>& matrix) {
	int N = matrix[0].size() - 1;
	double tmp = 0;
	vector<double> xx(N);
	short int i, j, k;
	vector<double> res;
	for (i = 0; i < N; i++)
	{
		tmp = matrix[i][i];
		for (j = N; j >= i; j--)
			matrix[i][j] /= tmp;
		for (j = i + 1; j < N; j++)
		{
			tmp = matrix[j][i];
			for (k = N; k >= i; k--)
				matrix[j][k] -= tmp * matrix[i][k];
		}
	}
	xx[N - 1] = matrix[N - 1][N];
	for (i = N - 2; i >= 0; i--)
	{
		xx[i] = matrix[i][N];
		for (j = i + 1; j < N; j++) xx[i] -= matrix[i][j] * xx[j];
	}
	for (i = 0; i < N; i++)
		res.push_back(xx[i]);
	return res;
}

void x_eq(double a, double b, int n, vector<double>& x) {
	double step = (b - a) / (n - 1);
	for (int i = 0; i < n; i++) {
		x.push_back(a + i * step);
	}
}

void x_opt(double a, double b, int n, vector<double>& x) {
	double x_i = 0;
	for (int i = 0; i != n; i++) {
		x_i = 0.5 * ((b - a) * cos((2 * i + 1) * acos(-1.0) / (2 * n + 2)) + b + a);
		x.push_back(x_i);
	}
	reverse(x.begin(), x.end());
}

double linSpline(vector<double>& x, vector<double>& y, double x_) {
	int k = -1;
	int n = x.size();
	for (int i = 0; i < n - 1; i++) {
		if (x_ >= x[i] && x_ <= x[i + 1]) {
			k = i;
			break;
		}
	}
	vector<vector<double>> slae(2, vector<double>(3));
	for (int i = 0; i < 2; i++) {
		slae[i][0] = x[k + i];
		slae[i][1] = 1;
		slae[i][2] = y[k + i];
	}
	vector<double> a = gauss(slae);
	return a[0] * x_ + a[1];
}

struct quadSpline {
	quadSpline(vector<double>& x, vector<double>& y) {
		this->x = x;
		this->y = y;
		this->n = x.size();
		int k = 0;
		vector<vector<double>> slae(3 * (n - 1), vector<double>(3 * n - 2));
		for (int i = 2; i < 3 * (n - 1); i += 3)
		{
			slae[i - 2][i - 2] = x[k] * x[k];
			slae[i - 2][i - 1] = x[k];
			slae[i - 2][i] = 1;
			slae[i - 1][i - 2] = x[k + 1] * x[k + 1];
			slae[i - 1][i - 1] = x[k + 1];
			slae[i - 1][i] = 1;
			slae[i][i - 2] = 2 * x[k + 1];
			slae[i][i - 1] = 1;
			if (i != 3 * (n - 1) - 1) {
				slae[i][i + 1] = -2 * x[k + 1];
				slae[i][i + 2] = -1;
			}
			slae[i - 2][3*n - 3] = y[k];
			slae[i - 1][3*n - 3] = y[k+1];
			slae[i][3*n - 3] = 0;
			k++;
		}
		double temp = 0;
		for (int i = 0; i < 3 * (n - 1); i++)
		{
			if (slae[i][i] == 0 && i != (3 * (n - 1) - 1))
			{
				for (int j = i + 1; j < 3 * (n - 1); j++)
				{
					if (slae[j][i] != 0)
					{
						for (int m = 0; m < 3 * (n - 1); m++)
						{
							temp = slae[j][m];
							slae[j][m] = slae[i][m];
							slae[i][m] = temp;
						}
						temp = slae[i][3*n-3];
						slae[i][3*n-3] = slae[j][3*n - 3];
						slae[j][3*n-3] = temp;
						break;
					}
				}
			}
			else if (slae[i][i] == 0 && i == (3 * (n - 1) - 1))
			{
				for (int j = i - 1; j >= 0; j--)
				{
					if (slae[j][j] != 0)
					{
						for (int m = 0; m < 3 * (n - 1); m++)
						{
							temp = slae[j][m];
							slae[j][m] = slae[i][m];
							slae[i][m] = temp;
						}
						temp = slae[i][n-1];
						slae[i][3*n-3] = slae[j][3*n - 3];
						slae[j][3*n-3] = temp;
						break;
					}
				}
			}
		}
		this->M = slae;
		this->a = gauss(slae);
	}
	double get(double x_) {
		int k = -1;
		for (int i = 0; i < n - 1; i++) {
			if (x_ >= x[i] && x_ <= x[i + 1]) {
				k = i;
				break;
			}
		}
		return a[3*k] * x_ * x_ + a[3*k + 1] * x_ + a[3*k+2];
	}
	void printa() {
		for (auto i : a) {
			cout << i << '\n';
		}
	}
	void printM() {
		for (int i = 0; i < 3 * (n - 1); i++) {
			for (int j = 0; j < 3 * (n - 1) + 1; j++) {
				cout << M[i][j] << '\t';
			}
			cout << '\n';
		}
	}
private:
	vector<double> x;
	vector<double> y;
	vector<double> a;
	vector<vector<double>> M;
	int n;
};

struct cubSpline {
	cubSpline(vector<double>& x, vector<double>& y) {
		this->n = x.size();
		this->x = x;
		this->y = y;
		int k = -1;
		double res = 0;
		vector<double> h(n - 1);
		for (int i = 0; i < n - 1; i++) {
			h[i] = x[i + 1] - x[i];
		}
		vector<double> gamma(n - 1);
		for (int i = 1; i < n - 1; i++) {
			gamma[i] = 6 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
		}
		vector<vector<double>> H(n - 2, vector<double>(n - 1));
		for (int i = 0; i < n - 2; i++) {
			for (int j = 0; j < n - 2; j++) {
				if (i == j) {
					H[i][j] = 2 * (h[i] + h[i + 1]);
				}
				else if (i - j == 1) {
					H[i][j] = h[i];
				}
				else if (j - i == 1) {
					H[i][j] = h[j];
				}
				else {
					H[i][j] = 0;
				}
			}
		}
		for (int i = 0; i < n - 2; i++) {
			H[i][n - 2] = gamma[i + 1];
		}
		vector<double> y_sec(n - 2);
		y_sec = gauss(H);
		y_sec.insert(y_sec.begin(), 0);
		y_sec.push_back(0);
		vector<double> y_f(n - 1);
		for (int i = 0; i < n - 1; i++) {
			y_f[i] = (y[i + 1] - y[i]) / h[i] - y_sec[i + 1] * h[i] / 6 - y_sec[i] * h[i] / 3;
		}
		this->y_f = y_f;
		this->y_sec = y_sec;
		this->h = h;
	}
	double get(double x_) {
		int k = -1;
		for (int i = 0; i < n - 1; i++) {
			if (x_ >= x[i] && x_ <= x[i + 1]) {
				k = i;
				break;
			}
		}
		//cout << k << '\n';
		return y[k] + y_f[k] * (x_ - x[k]) + y_sec[k] * pow((x_ - x[k]), 2) / 2 + (y_sec[k + 1] - y_sec[k]) * pow((x_ - x[k]), 3) / (6 * h[k]);
	}
	void printx() {
		for (auto i : x) {
			cout << i << '\n';
		}
	}
private:
	int n;
	vector<double> x;
	vector<double> y;
	vector<double> h;
	vector<double> y_f;
	vector<double> y_sec;
};

int main() {
	double a = 1, b = 10, r_eq = 0, r_opt = 0;
	vector<double> x;
	vector<double> x1;
	vector<double> y;
	vector<double> y1;
	cout << "Linear spline:\n";
	cout << "N nodes:\tErr eq:\t\tErr opt:\n";
	for (int n = 10; n <= 100; n += 10) {
		x_eq(a, b, n, x);
		x_opt(a, b, n, x1);
		for (int i = 0; i < x.size(); i++) {
			y.push_back(fun(x[i]));
			y1.push_back(fun(x1[i]));
		}
		for (double i = x[0]; i <= x[n - 1]; i += 0.1) {
			if (abs(fun(i) - linSpline(x, y, i)) > r_eq) {
				r_eq = abs(fun(i) - linSpline(x, y, i));
			}
		}
		for (double i = x1[0]; i <= x1[n - 1]; i += 0.1) {
			if (abs(fun(i) - linSpline(x1, y1, i)) > r_opt) {
				r_opt = abs(fun(i) - linSpline(x1, y1, i));
			}
		}
		cout << n << "\t\t" << r_eq << '\t' << r_opt << '\n';
		r_eq = 0;
		r_opt = 0;
		x.clear();
		x1.clear();
		y.clear();
		y1.clear();
	}
	cout << "Quadratic spline:\n";
	cout << "N nodes:\tErr eq:\t\tErr opt:\n";
	/*
	x = { 0, 2, 3 };
	y = { 1, 3, 2 };
	quadSpline quad_eq(x, y);
	cout << quad_eq.get(1);
	*/
	for (int n = 10; n <= 100; n += 10) {
		x_eq(a, b, n, x);
		x_opt(a, b, n, x1);
		for (int i = 0; i < x.size(); i++) {
			y.push_back(fun(x[i]));
			y1.push_back(fun(x1[i]));
		}
		quadSpline quad_eq(x, y);
		quadSpline quad_opt(x1, y1);
		for (double i = x[0]; i <= x[n - 1]; i += 0.1) {
			if (abs(fun(i) - quad_eq.get(i)) > r_eq) {
				r_eq = abs(fun(i) - quad_eq.get(i));
			}
		}
		for (double i = x1[0]; i <= x1[n - 1]; i += 0.1) {
			if (abs(fun(i) - quad_opt.get(i)) > r_opt) {
				r_opt = abs(fun(i) - quad_opt.get(i));
			}

		}
		cout << n << "\t\t" << r_eq << '\t' << r_opt << '\n';
		r_eq = 0;
		r_opt = 0;
		x.clear();
		x1.clear();
		y.clear();
		y1.clear();

	}
	cout << "Cubic spline:\n";
	cout << "N nodes:\tErr eq:\t\tErr opt:\n";
	for (int n = 10; n <= 100; n += 10) {
		x_eq(a, b, n, x);
		x_opt(a, b, n, x1);
		for (int i = 0; i < x.size(); i++) {
			y.push_back(fun(x[i]));
			y1.push_back(fun(x1[i]));
		}
		cubSpline s_eq(x, y);
		cubSpline s_opt(x1, y1);
		ofstream fout3("r_eq_s" + to_string(n) + ".txt");
		for (double i = x[0]; i <= x[n - 1]; i += 0.1) {
			if (abs(fun(i) - s_eq.get(i)) > r_eq) {
				r_eq = abs(fun(i) - s_eq.get(i));
			}
			fout3 << i << '\t' << abs(fun(i) - s_eq.get(i)) << '\n';

		}
		fout3.close();
		ofstream fout4("r_opt_s" + to_string(n) + ".txt");
		for (double i = x1[0]; i <= x1[n-1]; i += 0.1) {
			if (abs(fun(i) - s_opt.get(i)) > r_opt) {
				r_opt = abs(fun(i) - s_opt.get(i));
			}
			fout4 << i << '\t' << abs(fun(i) - s_opt.get(i)) << '\n';
			
		}
		fout4.close();
		cout << n << "\t\t" << r_eq << '\t' << r_opt << '\n';
		r_eq = 0;
		r_opt = 0;
		x.clear();
		x1.clear();
		y.clear();
		y1.clear();
	}
}