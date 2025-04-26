#include<iostream>
#include<vector>
#include<cmath>
#include<string>
#include<fstream>
using namespace std;

double fun(double x) {
	return x * x * cos(x);
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


void createData(double a, double b, int n, vector<double>& x, vector<double>& y) {
	double x_i = 0;		
	for (int i = 0; i != n; i++) {
		double eps = 0.001 * (rand() % 100) + 0.01;
		double step = (b - a) / (n - 1);
		x_i = a + i * step;
		x.push_back(x_i);
		x.push_back(x_i);
		x.push_back(x_i);
		y.push_back(fun(x_i));
		y.push_back(fun(x_i) + eps);
		y.push_back(fun(x_i) - eps);
	}
}

double mnkOrt(vector<double>& x, vector<double>& y, int n) {
	vector<double> a(n);
	vector<vector<double>> q(n, vector<double>(x.size()));
	vector<vector<double>> data(n, vector<double>(22));
	vector<double> xn(22);
	int p = 0;
	for (double i = -1; i <= 1; i += 0.1) {
		xn[p] = i;
		p++;
	}
	p = 0;
	double sum = 0;
	for (auto i : x) {
		sum += i;
	}
	for (int j = 0; j < x.size(); j++) {
		q[0][j] = 1;
		q[1][j] = x[j] - sum / (x.size());

	}
	for (int j = 0; j < 21; j++) {
		data[0][j] = 1;
		data[1][j] = xn[j] - sum / (x.size());
	}
	if (n > 2) {
		for (int i = 2; i < n; i++) {
			double a_n = 0, a_d = 0;
			double b_n = 0, b_d = 0;
			for (int j = 0; j < x.size(); j++) {
				a_n += x[j] * q[i - 1][j] * q[i - 1][j];
				a_d += q[i - 1][j] * q[i - 1][j];
				b_n += x[j] * q[i - 1][j] * q[i - 2][j];
				b_d += q[i - 2][j] * q[i - 2][j];
			}
			a_n /= a_d;
			b_n /= b_d;
			for (int j = 0; j < x.size(); j++) {
				q[i][j] = x[j] * q[i - 1][j] - a_n * q[i - 1][j] - b_n * q[i - 2][j];
			}
		}
	}
	if (n > 2) {
		for (int i = 2; i < n; i++) {
			double a_n = 0, a_d = 0;
			double b_n = 0, b_d = 0;
			for (int j = 0; j < 21; j++) {
				a_n += xn[j] * data[i - 1][j] * data[i - 1][j];
				a_d += data[i - 1][j] * data[i - 1][j];
				b_n += xn[j] * data[i - 1][j] * data[i - 2][j];
				b_d += data[i - 2][j] * data[i - 2][j];
			}
			a_n /= a_d;
			b_n /= b_d;
			for (int j = 0; j < 21; j++) {
				data[i][j] = xn[j] * data[i - 1][j] - a_n * data[i - 1][j] - b_n * data[i - 2][j];
			}
		}
	}
	for (int i = 0; i < n; i++) {
		double Eqf = 0, Eqf_d = 0;
		for (int j = 0; j < x.size(); j++) {
			Eqf += q[i][j] * y[j];
			Eqf_d += q[i][j] * q[i][j];
		}
		Eqf /= Eqf_d;
		a[i] = Eqf;
	}
	vector<double> mnk(x.size());

	vector<double> out(22);
	ofstream fout("output_ort" + to_string(n - 1) + ".txt");
	ofstream fout1("EXPort" + to_string(n - 1) + ".txt");
	for (int i = 0; i < x.size(); i++) {
		fout1 << x[i] << '\t' << y[i] << '\n';
	}
	for (int i = 0; i < x.size(); i++) {
		double mnk_t = 0;
		for (int j = 0; j < n; j++) {
			mnk_t += a[j] * q[j][i];
		}
		mnk[i] = mnk_t;
	}
	for (int i = 0; i < 21; i++) {
		double mnk_t = 0;
		for (int j = 0; j < n; j++) {
			mnk_t += a[j] * data[j][i];
		}
		out[i] = mnk_t;
		fout << xn[i] << '\t' << out[i] << '\n';
	}
	double err_sum = 0;
	for (int i = 0; i < x.size(); i++) {
		double err = y[i] - mnk[i];
		err_sum += err * err;
	}
	return err_sum;
}

double mnkNorm(vector<double>& x, vector<double>& y, int n) {
	vector<double> sum(2 * n - 1);
	vector<double> sum_f(n);
	sum[0] = x.size();
	vector<double> xn(101);
	int p = 0;
	for (double i = -1; i <= 1; i += 0.02) {
		xn[p] = i;
		p++;
	}
	for (int i = 1; i < 2 * n - 1; i++) {
		for (int j = 0; j < x.size(); j++) {
			sum[i] += pow(x[j], i);
		}
	}
	for (int i = 0; i <  n; i++) {
		for (int j = 0; j < x.size(); j++) {
			sum_f[i] += pow(x[j], i) * y[j];
		}
	}
	ofstream fout1("normEXP" + to_string(n - 1) + ".txt");
	for (int i = 0; i < x.size(); i++) {
		fout1 << x[i] << '\t' << y[i] << '\n';
	}
	vector<vector<double>> matrix(n, vector<double>(n + 1));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			matrix[i][j] = sum[i + j];
		}
		matrix[i][n] = sum_f[i];
	}
	vector<double> a(n);

	a = gauss(matrix);
	vector<double> mnk(x.size());
	for (int i = 0; i < x.size(); i++) {
		double mnk_t = 0;
		for (int j = 0; j < n; j++) {
			mnk_t += a[j] * pow(x[i], j);
		}
		mnk[i] = mnk_t;
	}
	ofstream fout("NORM" + to_string(n) + ".txt");
	vector<double> out(101);
	for (int i = 0; i < 101; i++) {
		double mnk_t = 0;
		for (int j = 0; j < n; j++) {
			mnk_t += a[j] * pow(xn[i], j);
		}
		out[i] = mnk_t;
		fout << xn[i] << '\t' << out[i] << '\n';

	}
	double err_sum = 0;
	for (int i = 0; i < x.size(); i++) {
		double err = y[i] - mnk[i];
		err_sum += err * err;
	}
	return err_sum;
}

int main() {
	vector<double> x;
	vector<double> y;
	cout << "Pol degree\tSum of errors(norm)\tSumof errors(ort)\n";
	createData(-1, 1, 11, x, y);
	for (int n = 2; n <= 5; n++) {
		cout << n-1 << '\t' << mnkNorm(x, y, n) << "\t" << mnkOrt(x, y, n) << '\n';;

	}
	return 0;
}