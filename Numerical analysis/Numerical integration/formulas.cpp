#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include<fstream>
using namespace std;

double pi = 2 * acos(0.0);
double a = 0.1;//1.5;
double b = 2.3;//3.3;
double alpha = 1. / 5;//1. / 3;
double betha = 0;
double val = 7.258002984374563;
double val_p = 7.077031437995776;
double eps = 10e-6;
//gauss slae
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
//Cardano formulas
vector<double> crdn(vector<double>& a) {
    vector<double> res(3);
    double b = a[0];
    double c = a[1];
    double d = a[2];
    double y = b / 3;
    double q = (2 * (b*b*b) / 27 - (b * c / 3) + d) / 2;
    double p = (3 * c - b * b) / 9;
    double r = sqrt(abs(p)) * (signbit(q) ? -1 : 1);
    //cout << q << p << r << endl;
    if (p < 0 && (pow(q, 2) + pow(p, 3) <= 0)) {
        double phi = acos(q / pow(r, 3));
        res = { -2 * r * cos(phi / 3) - y, 2 * r * cos(pi / 3 - phi / 3) - y, 2 * r * cos(pi / 3 + phi / 3) - y };
        return res;
    }
    cout << "BAD ROOTS %)";
    return {0, 0, 0};
}
//f(x)
double fun(double x) {
    return 2.5 * cos(2 * x) * exp(2 * x / 3) + 4 * sin(3.5 * x) * exp(-3 * x) + 3 * x;
}
/*
double fun(double x) {
    return 2 * cos(2.5 * x) * exp(x / 3) + 
        4 * sin(3.5 * x) * exp(-3 * x) + x;
}
*/
//p(x)
double p(double a, double b, double x) {
    return pow((x - a), -alpha) * pow((b - x), -betha);
}
//nodes
void x_eq(double a, double b, int n, vector<double>& x) {
    double step = (b - a) / (n - 1);
    //cout << step << ' ';
    for (int i = 0; i < n; i++) {
        x[i] = a + i * step;
    }
}
//mid rct
double mid_rct(double a, double b) {
    return (b - a) * fun((a + b) / 2);
}
//left rct
double lft_rct(double a, double b) {
    return (b - a) * fun(a);
}
//trp
double trp(double a, double b) {
    return (b - a) / 2 * (fun(a) + fun(b));
}
//Simpson
double smp(double a, double b) {
    return (b - a) / 6 * (fun(a) + 4 * fun((a + b) / 2) + fun(b));
}
//Newton-Cotes
double nwt_cts(int n) {
    vector<double> z_(n + 1);
    x_eq(a, b, n + 1, z_);
    double res = 0;
    for (int i = 0; i < z_.size() - 1; i++) {
        double z0 = z_[i];
        double z1 = z_[i + 1];
        double z = (z0 + z1) / 2;
        double mu0 = (pow((z1 - a), (1 - alpha)) -
            pow((z0 - a), (1 - alpha))) / (1 - alpha);
        double mu1 = (pow((z1 - a), (2 - alpha)) -
            pow((z0 - a), (2 - alpha))) / (2 - alpha) + a * mu0;
        double mu2 = (pow((z1 - a), (3 - alpha)) -
            pow((z0 - a), (3 - alpha))) / (3 - alpha) + 2 * a * mu1 - a * a * mu0;
        double a_1 = (mu2 - mu1 * (z + z1) + mu0 * z * z1) /
            ((z - z0) * (z1 - z0));
        double a_2 = -(mu2 - mu1 * (z0 + z1) + mu0 * z0 * z1) /
            ((z - z0) * (z1 - z));
        double a_3 = (mu2 - mu1 * (z + z0) + mu0 * z * z0) /
            ((z1 - z) * (z1 - z0));
        res+= a_1 * fun(z0) + a_2 * fun(z) + a_3 * fun(z1);
    }
    return res;
}
//Gauss
double gs_qdr(int n) {
    vector<double> z_(n + 1);
    x_eq(a, b, n + 1, z_);
    double res = 0;
    for (int i = 0; i < z_.size() - 1; i++) {
        double z0 = z_[i];
        double z1 = z_[i + 1];
        double mu0 = (pow((z1 - a), (1 - alpha)) -
            pow((z0 - a), (1 - alpha))) / (1 - alpha);
        double mu1 = (pow((z1 - a), (2 - alpha)) -
            pow((z0 - a), (2 - alpha))) / (2 - alpha) + a * mu0;
        double mu2 = (pow((z1 - a), (3 - alpha)) -
            pow((z0 - a), (3 - alpha))) / (3 - alpha) + 2 * a * mu1 - a * a * mu0;
        double mu3 = (pow((z1 - a), (4 - alpha)) - pow((z0 - a), (4 - alpha))) /
            (4 - alpha) + 3 * a * mu2 - 3 * a * a * mu1 + pow(a, 3) * mu0;
        double mu4 = (pow((z1 - a), (5 - alpha)) - pow((z0 - a), (5 - alpha))) /
            (5 - alpha) + 4 * a * mu3 - 6 * a * a * mu2 + 4 * pow(a, 3) * mu1 - pow(a, 4) * mu0;
        double mu5 = (pow((z1 - a), (6 - alpha)) - pow((z0 - a), (6 - alpha))) /
            (6 - alpha) + 5 * a * mu4 - 10 * a * a * mu3 + 10 * pow(a, 3) * mu2 - 5 * pow(a, 4) * mu1 + pow(a, 5) * mu0;
        vector<vector<double>> slae(3, vector<double>(4));
        slae[0] = { mu0, mu1, mu2, -mu3 };
        slae[1] = { mu1, mu2, mu3, -mu4 };
        slae[2] = { mu2, mu3, mu4, -mu5 };

        vector<double> a(3);
        a = gauss(slae);
        reverse(a.begin(), a.end());
        //find roots
        vector<double> roots(3);
        roots = crdn(a);
        sort(roots.begin(), roots.end());
        slae[0] = { 1, 1, 1, mu0 };
        slae[1] = { roots[0], roots[1], roots[2], mu1 };
        slae[2] = { pow(roots[0], 2), pow(roots[1], 2), pow(roots[2], 2), mu2 };
        a = gauss(slae);
        res += a[0] * fun(roots[0]) + a[1] * fun(roots[1]) + a[2] * fun(roots[2]);
    }
    return res;
}
int main()
{   

    ofstream fout("data.txt");
    fout << "n\tmid\tleft\ttrp\tsimpson\tnew-cot\tgauss\n";
    for (int n = 2; n <= 100; n += 1) {
        double sum_mid = 0, sum_lft = 0, sum_trp = 0, sum_sim = 0, sum_nc = 0, sum_gs = 0;
        vector<double> x(n);
        x_eq(a, b, n, x);
        for (int i = 0; i < x.size() - 1; i++) {
            sum_mid += mid_rct(x[i], x[i + 1]);
            sum_lft += lft_rct(x[i], x[i + 1]);
            sum_trp += trp(x[i], x[i + 1]);
            sum_sim += smp(x[i], x[i + 1]);
            //sum_nc += nwt_cts(x[i], x[i + 1]);
            //sum_gs += gs_qdr(x[i], x[i + 1]);
        }
        sum_nc = nwt_cts(n);
        sum_gs = gs_qdr(n);
        cout << sum_mid << '\n';
        //fout << n << '\t' << abs(val - sum_mid) << '\t' << abs(val - sum_lft) << '\t' << abs(val - sum_trp) << '\t' << abs(val - sum_sim) << '\t' << abs(val_p - sum_nc) << '\t' << abs(val_p - sum_gs) << '\n';
    }
    
    cout << "Newton-Cotes formula\nValue\tStep\tError\n";
    double d = b - a;
    vector<double> s_h = {-nwt_cts(1), -nwt_cts(2)};

    for (int r = 2; r < 10; r++) {
        s_h.push_back(-nwt_cts(pow(2, r)));
        int size = s_h.size();
        double m = -log(abs((-s_h[size - 1] + s_h[size - 2]) 
            / (-s_h[size - 2] + s_h[size - 3]))) / log(2);
        cout << "m = " << m << endl;
        vector<vector<double>> h_m(r + 1, vector<double>(r + 2));
        for (int i = 0; i < r + 1; i++) {
            for (int j = 0; j < r + 2; j++) {
                if (j == r + 1) {
                    h_m[i][j] = s_h[i];
                }
                else if (j != r) {
                    h_m[i][j] = pow((d / pow(2, i)), m + j);
                }
                else {
                    h_m[i][j] = -1;
                }
            }
        }
        vector<double> c_m(r + 1);
        c_m = gauss(h_m);
        if (abs(c_m[r] + s_h[size - 1]) < eps) {
            cout << -s_h[size - 1] << '\t' << d/pow(2, r) << '\t' << abs(c_m[r] + s_h[size - 1]) << '\t' << r << '\n';
            break;
        }
    }

    cout << "Gaussian formula\nValue\tStep\n";
    s_h = { -gs_qdr(1), -gs_qdr(2) };
    for (int r = 2; r < 10; r++) {
        s_h.push_back(-gs_qdr(pow(2, r)));
        int size = s_h.size();
        double m = -log(abs((-s_h[size - 1] + s_h[size - 2])
            / (-s_h[size - 2] + s_h[size - 3]))) / log(2);
        cout << "m = " << m << endl;
        vector<vector<double>> h_m(r + 1, vector<double>(r + 2));
        for (int i = 0; i < r + 1; i++) {
            for (int j = 0; j < r + 2; j++) {
                if (j == r + 1) {
                    h_m[i][j] = s_h[i];
                }
                else if (j != r) {
                    h_m[i][j] = pow((d / pow(2, i)), m + j);
                }
                else {
                    h_m[i][j] = -1;
                }
            }
        }
        vector<double> c_m(r + 1);
        c_m = gauss(h_m);
        cout << "CM " << c_m[r] << '\n';
        if (abs(c_m[r] + s_h[size - 1]) < eps) {
            cout <<  -s_h[size - 1] << '\t' << d / pow(2, r) << '\t' << abs(-s_h[size - 1] - val_p) << '\n';
            break;
        }
    }

    s_h = { -nwt_cts(1), -nwt_cts(2), -nwt_cts(4) };
    int size = s_h.size();
    double m = -log(abs((-s_h[size - 1] + s_h[size - 2])
        / (-s_h[size - 2] + s_h[size - 3]))) / log(2);
    vector<vector<double>> h_m(3, vector<double>(4));
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            if (j == 3) {
                h_m[i][j] = s_h[i];
            }
            else if (j != 2) {
                h_m[i][j] = pow((d / pow(2, i)), m + j);
            }
            else {
                h_m[i][j] = -1;
            }
        }
    }
    vector<double> c_m(3);
    c_m = gauss(h_m);
    double h_opt = d * pow((eps / abs(s_h[size - 3] + c_m[2])), 1 / m);
    cout << "h opt = " << h_opt << '\n';
    cout << "Newton-Cotes formula with optimal step\nValue\tStep\tError\n";
    s_h = { -nwt_cts(floor(d/h_opt)), -nwt_cts(floor(2*d/(h_opt))) };
    for (int r = 2; r < 10; r++) {
        s_h.push_back(-nwt_cts(floor(pow(2, r)*d/(h_opt))));
        int size = s_h.size();
        double m = -log(abs((-s_h[size - 1] + s_h[size - 2])
            / (-s_h[size - 2] + s_h[size - 3]))) / log(2);
        //cout << "m = " << m << endl;
        vector<vector<double>> h_m(r + 1, vector<double>(r + 2));
        for (int i = 0; i < r + 1; i++) {
            for (int j = 0; j < r + 2; j++) {
                if (j == r + 1) {
                    h_m[i][j] = s_h[i];
                }
                else if (j != r) {
                    h_m[i][j] = pow((d / pow(2, i)), m + j);
                }
                else {
                    h_m[i][j] = -1;
                }
            }
        }
        vector<double> c_m(r + 1);
        c_m = gauss(h_m);
        cout << "CM " << c_m[r] + s_h[size - 1] << '\n';
        if (abs(c_m[r] + s_h[size - 1]) < eps) {
            cout << -s_h[size - 1] << '\t' << d / floor(pow(2, r) * d / (h_opt)) << '\t' << abs(c_m[r] + s_h[size - 1]) << '\n';
            break;
        }
    }
}
