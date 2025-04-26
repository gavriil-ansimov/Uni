#ifndef GAUSS
#define GAUSS

#include<iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#define N    3
#define N1    N+1
using namespace std;
float matrix1[N][N1] = { {4, 1, 1, -1},
                        {1, 6.2, -1, 2},
                        {1, -1, 8.2, -3} };


float matrix[3][4] = { {4, 1, 1, -1},
                        {1, 6.2, -1, 2},
                        {1, -1, 8.2, -3} };

float epsilon = 0.001;

double delta() {
    double del1, del2, del3;
    del1 = abs(matrix[0][0]) - abs(matrix[0][1]) - abs(matrix[0][2]);
    del2 = abs(matrix[1][1]) - abs(matrix[1][0]) - abs(matrix[1][2]);
    del3 = abs(matrix[2][2]) - abs(matrix[2][1]) - abs(matrix[2][0]);
    if (del1 > 0 && del2 > 0 && del3 > 0) {
        return min(del1, min(del2, del3));
    } else {
        return 0;
    }
}

double d = delta();

double fun(vector<double> x) {
    return 2 * x[0] * x[0] + 6.2 * x[1] * x[1] + 8.2 * x[2] * x[2] + x[0] * x[1] - x[1] * x[2] + x[0] * x[2] + x[0] - 2 * x[1] + 3 * x[2] + 1;
}




double gradx(vector<double> x) {
    return 4 * x[0] + x[1] + x[2] + 1;
}

double grady(vector<double> x) {
    return 6.2 * x[1] + x[0] - x[2] - 2;
}

double gradz(vector<double> x) {
    return 8.2 * x[2] + x[0] - x[1] + 3;
}

vector<double> multAplusB(vector<double> x) {
    vector<double> res{ 0, 0, 0 };
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            res[i] += matrix[i][j] * x[j];
        }
    }
    for (int i = 0; i < 3; i++) {
        res[i] -= matrix[i][3];
    }
    return res;
}

vector<double> multA(vector<double> x) {
    vector<double> res(3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            res[i] += matrix[i][j] * x[j];
        }
    }
    return res;
}

double strvec(vector<double> str, vector<double> vec) {
    return str[0] * vec[0] + str[1] * vec[1] + str[2] * vec[2];
}

vector<double> gauss()
{
    double tmp, xx[N1];
    short int i, j, k;
    vector<double> res;
    /*Метод Гаусса*/
    /*прямой ход*/
    for (i = 0; i < N; i++)
    {
        tmp = matrix1[i][i];
        for (j = N; j >= i; j--)
            matrix1[i][j] /= tmp;
        for (j = i + 1; j < N; j++)
        {
            tmp = matrix1[j][i];
            for (k = N; k >= i; k--)
                matrix1[j][k] -= tmp * matrix1[i][k];
        }
    }
    /*обратный ход*/
    xx[N - 1] = matrix1[N - 1][N];
    for (i = N - 2; i >= 0; i--)
    {
        xx[i] = matrix1[i][N];
        for (j = i + 1; j < N; j++) xx[i] -= matrix1[i][j] * xx[j];
    }
    for (i = 0; i < N; i++)
        res.push_back(xx[i]);
    return res;
}

#endif