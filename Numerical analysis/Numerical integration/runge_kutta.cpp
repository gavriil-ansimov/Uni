#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

double pi = 2 * acos(0.0);
double xi = 0.05;
double A = 1. / 30;
double B = 1. / 15;
double eps = pow(10, -5), ro = pow(10, -6);
double x0 = 0, x1 = pi;
double a = xi, b2 = 1.0 / (2 * xi), b1 = 1 - b2;
//dy1/dx = Ay2
//dy2/dx = -By1
//y1(0) = B*pi; y2(0) = A*pi
//x_k = pi; yi(x_k) = ?

double f1(double y) {
    return A * y;
}

double f2(double y) {
    return -B * y;
}

vector<double> c_n(double x, double y1, double y2) {
    vector<double> res(2);
    res[0] = (2 * y1 * cos(x / (15 * sqrt(2))) - sqrt(2) * y2 * sin(x / (15 * sqrt(2))))/2;
    res[1] = y2 * cos(x / (15 * sqrt(2))) + sqrt(2) * y1 * sin(x / (15 * sqrt(2)));
    return res;
}

vector<double> val_n(double x, double c1, double c2) {
    vector<double> res(2);
    res[0] = c2 * sin(x / (15 * sqrt(2))) / sqrt(2) + c1 * cos(x / (15 * sqrt(2)));
    res[1] = c2 * cos(x / (15 * sqrt(2))) - sqrt(2) * c1 * sin(x / (15 * sqrt(2)));
    return res;
}

vector<double> val(double x) {
    vector<double> res(2);
    res[0] = pi * sin(x / (15 * sqrt(2))) / (30 * sqrt(2)) + pi * cos(x / (15 * sqrt(2))) / 15;
    res[1] = pi * cos(x / (15 * sqrt(2))) / 30 - pi * sqrt(2) * sin(x / (15 * sqrt(2))) / 15;
    return res;
}

vector<double> sol2(double y_10, double y_20, double step) {
    int cnt = 0;
    vector<vector<double>> k(2, vector<double>(2));
    double y1 = y_10, y2 = y_20;
    double x = x0;
    double h = step;
    string name = "2x2_h.txt";
    ofstream fout(name);
    double norm = 0;
    vector<double> ref(2);
    while (x <= x1) {
        
        k[0][0] = step * f1(y2);
        k[1][0] = step * f2(y1);
        k[0][1] = step * f1(y2 + a * k[0][0]);
        k[1][1] = step * f2(y1 + a * k[1][0]);
        y1 = y1 + b1 * k[0][0] + b2 * k[0][1];
        y2 = y2 + b1 * k[1][0] + b2 * k[1][1];
        x += step;
        cnt++;
        ref = val(x);
        norm = sqrt(pow(y1 - ref[0], 2) + pow(y2 - ref[1], 2));
        
        fout << x << '\t' << norm << '\n';
        
        //cout << y1 << '\n';
    }
    fout << pi;
    vector<double> res = { y1, y2 };
    return res;
}

vector<double> sol3(double y_10, double y_20, double step) {
    int cnt = 0;
    vector<vector<double>> k(2, vector<double>(3));
    double y1 = y_10, y2 = y_20, y0 = 0, y00 = 0;
    double x = x0;
    string name = "3x3_h.txt";
    ofstream fout(name);
    vector<double> ref(2);
    double norm = 0;
    while (x < x1) {
        
        k[0][0] = step * f1(y2);
        k[0][1] = step * f1(y2 + 1./3 * k[0][0]);
        k[0][2] = step * f1(y2 +2./3 * k[0][1]);
        k[1][0] = step * f2(y1);
        k[1][1] = step * f2(y1 + 1./3 * k[1][0]);
        k[1][2] = step * f2(y1 +2./3 * k[1][1]);
        y00 = y1;
        y0 = y2;
        y1 = y1 + (1. / 4) * (k[0][0] + 3 * k[0][2]);
        y2 = y2 + (1. / 4) * (k[1][0] + 3 * k[1][2]);
        x += step;
        ref = val(x);
        norm = sqrt(pow(y1 - ref[0], 2) + pow(y2 - ref[1], 2));
        fout << x << '\t' << norm << '\n';
        
        cnt++;
        //cout << y1 << '\n';
    }
    vector<double> res = { y1, y2};
    return res;
}

vector<double>sch2(vector<double> y, double step) {
    vector<vector<double>> k(2, vector<double>(2));
    k[0][0] = step * f1(y[1]);
    k[1][0] = step * f2(y[0]);
    k[0][1] = step * f1(y[1] + a * k[0][0]);
    k[1][1] = step * f2(y[0] + a * k[1][0]);
    double y1 = y[0] + b1 * k[0][0] + b2 * k[0][1];
    double y2 = y[1] + b1 * k[1][0] + b2 * k[1][1];
    return { y1, y2 };
}

vector<double>sch3(vector<double> y, double step) {
    vector<vector<double>> k(2, vector<double>(3));
    k[0][0] = step * f1(y[1]);
    k[0][1] = step * f1(y[1] + 0.5 * k[0][0]);
    k[0][2] = step * f1(y[1] - k[0][0] + 2 * k[0][1]);
    k[1][0] = step * f2(y[0]);
    k[1][1] = step * f2(y[0] + 0.5 * k[1][0]);
    k[1][2] = step * f2(y[0] - k[1][0] + 2 * k[1][1]);
    double y1 = y[0] + 1. / 6 * (k[0][0] + 4 * k[0][1] + k[0][2]);
    double y2 = y[1] + 1. / 6 * (k[1][0] + 4 * k[1][1] + k[1][2]);
    return { y1, y2 };
}

vector<double> auto_step(double x, double y1, double y2, double h, vector<double>& x_, vector<double>& y, vector<double>& error) {
    vector<double> y_s(2);
    vector<double> y_2s(2);
    double step = h;
    
    y_s = sch2({ y1, y2 }, step);
    y_2s = sch2({ y1, y2 }, step / 2);
    y_2s = sch2({ y_2s[0], y_2s[1]}, step / 2);
    double err = max(abs((y_2s[0] - y_s[0]) / 0.75), abs((y_2s[1] - y_s[1]) / 0.75));
    if (err > ro * 4) {
        return auto_step(x, y1, y2, step / 2, x_, y, error);
    }
    else if (err > ro && err <= ro * 4) {
        
        //double norm = sqrt(pow(y_2s[0] - ref[0], 2) + pow(y_2s[1] - ref[1], 2));
        return { x + step, step / 2, y_2s[0], y_2s[1], err};
    }
    else if (err >= ro / 8 && err <= ro) {

        //double norm = abs(y_s[0] - ref[0]);
        //double norm = sqrt(pow(y_s[0] - ref[0], 2) + pow(y_s[1] - ref[1], 2));
        return { x + step, step, y_s[0], y_s[1], err};
    }
    else {
        
        //double norm = abs(y_s[0] - ref[0]);
        //double norm = sqrt(pow(y_s[0] - ref[0], 2) + pow(y_s[1] - ref[1], 2));
        
        return { x + step, 2 * step, y_s[0], y_s[1], err};
    }
}

vector<double> auto_step3(double x, double y1, double y2, double step, vector<double>& x_, vector<double>& y) {
    vector<double> y_s(2);
    vector<double> y_2s(2);

    y_s = sch3({ y1, y2 }, step);
    y_2s = sch3({ y1, y2 }, step / 2);
    y_2s = sch3({ y_2s[0], y_2s[1] }, step / 2);
    double err = max(abs((y_2s[0] - y_s[0]) / 0.875), abs((y_2s[1] - y_s[1]) / 0.875));
    if (err > ro * 8) {
        return auto_step3(x, y1, y2, step / 2, x_, y);
    }
    else if (err > ro) {
        return { x + step, step / 2, y_2s[0], y_2s[1], err };
    }
    else if (err >= ro / 16 && err <= ro) {
        return { x + step, step, y_2s[0], y_2s[1], err };
    }
    else {
        return { x + step, 2 * step, y_2s[0], y_2s[1], err };
    }
}

int main()
{
    double const y_10 = B * pi;
    double const y_20 = A * pi;
    //1.1
    double step2x2 = 0, step3x3 = 0;
    vector<double> value = val(pi);
    double step = 0.1;
    vector<double> res = sol2(y_10, y_20, step);
    cout << "1.1\n";
    cout << "True errors: " << abs(res[0] - value[0]) << ' ' << abs(res[1] - value[1]) << '\n';
    //1.2
    //start step
    double delta = pow(1.0 / pi, 3) + pow(sqrt(pow(f1(y_20), 2) + pow(f2(y_10), 2)), 3);
    double h = pow(eps / delta, 1.0 / 3);
    double h1 = pow(ro / delta, 1.0 / 3);
    double R = 1;
    vector<double>y_1;
    vector<double>y_2;
    step = h;
    int k = 0;
    while (R >= eps) {
        k++;
        y_1 = sol2(y_10, y_20, step);
        y_2 = sol2(y_10, y_20, step / 2);
        R = sqrt(pow((y_2[0] - y_1[0]) / 3, 2) + pow((y_2[1] - y_1[1]) / 3, 2));
        step /= 2;
        //cout << R << '\n';
    }
    step2x2 = 2 * step;
    cout << "1.2\n";
    cout << "True errors: " << ' ' << abs(y_2[0] - value[0]) << ' ' << abs(y_2[1] - value[1]) << '\n';
    cout << "Step = " << step2x2 << '\n';
    //2.1
    step = h1;
    double x = x0;
    double y1 = y_10, y2 = y_20;
    res = { 0 };
    int cnt = 0;
    vector<double> x_, errors, true_err;
    double err = 0, err1 = 0;
    while (x < x1) {
        step = min(h1, pi - x);
        res = auto_step(x, y1, y2, step, x_, errors, true_err);
        vector<double> c = c_n(x, y1, y2);
        x = res[0];
        step = res[1];
        y1 = res[2];
        y2 = res[3];
        vector<double> y = val(x);
        vector<double> ref = val_n(x, c[0], c[1]);
        double norm = max(abs(y1 - ref[0]), abs(y2 - ref[1]));
        err = res[4];
        err1 = norm;
        cnt++;
        x_.push_back(x);
        errors.push_back(err);
        true_err.push_back(err1);
    }

    cout << "2.1\n";
    cout << "True errors: " << ' ' << abs(y1 - value[0]) << ' ' << abs(y2 - value[1]) << '\n';
    /*
    x_.push_back(x);
    errors.push_back(err);
    true_err.push_back(err1);
    */
    ofstream fout("auto2x2.txt");
    
    for (int i = 0; i < x_.size(); i++) {
        fout << x_[i] << '\t' << true_err[i]/errors[i] << '\n';
    }
    //3.1

    cout << "3.1\n";
    R = 1;
    delta = pow(1.0 / pi, 4) + pow(sqrt(pow(f1(y_20), 2) + pow(f2(y_10), 2)), 4);
    h = pow(eps / delta, 1.0 / 4);
    h1 = pow(ro / delta, 1.0 / 4);
    step = h;
    k = 0;

    while (R >= eps) {
        k++;
        y_1 = sol3(y_10, y_20, step);
        y_2 = sol3(y_10, y_20, step / 2);
        //R = max(abs((y_2[0] - y_1[0]) / 7), abs((y_2[1] - y_1[1]) / 7));
        R = sqrt(pow((y_2[0] - y_1[0]) / 7, 2) + pow((y_2[1] - y_1[1]) / 7, 2));
        step /= 2;
    }
    step3x3 = step * 2;
    cout << "True errors: " << abs(y_2[0]-value[0]) << ' ' << abs(y_2[1] - value[1]) << '\n';
    cout << "Step = " << step3x3 << '\n';
    step = h1;
    x = x0;
    y1 = y_10, y2 = y_20;
    res = { 0, 0, 0, 0 };
    x_.clear();
    errors.clear();
    true_err.clear();
    cnt = 0;
    while (x < x1) {
        step = min(h1, pi - x);
        res = auto_step3(x, y1, y2, step, x_, errors);
        vector<double> c = c_n(x, y1, y2);
        x = res[0];
        step = res[1];
        y1 = res[2];
        y2 = res[3];
        vector<double> y = val(x);
        vector<double> ref = val_n(x, c[0], c[1]);
        double norm = max(abs(y1 - ref[0]), abs(y2 - ref[1]));
        err = res[4];
        err1 = norm;
        cnt++;
        x_.push_back(x);
        errors.push_back(err);
        true_err.push_back(err1);
        
    }

    ofstream fout2("auto3x3.txt");
    for (int i = 0; i < x_.size(); i++) {
        fout2 << x_[i] << '\t' << true_err[i]/errors[i] << '\n';
    }
    cout << "Auto step\n";
    cout << "True errors: " << abs(y1 - value[0]) << ' ' << abs(y2 - value[1]) << '\n';
    res = sol2(y_10, y_20, step2x2);
    res = sol3(y_10, y_20, step3x3);
    
}


