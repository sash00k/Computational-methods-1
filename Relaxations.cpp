#define _CRT_SECURE_nO_WARNINGS
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace std;
void yOut(double* y, int n);
void xOut(double* x, int n);
void sysOut(double** arr, int n);
void relax(double** arr, double* x, double* y, int n);
void gauss(double** arr, double* x, double* y, int n);

int main() {
    setlocale(LC_ALL, "Russian");

    int n = 3;

    double** arr = new double* [n];
    double* y = new double[n];
    double* x = new double[n];
    for (int i = 0; i < n; i++) arr[i] = new double[n];

    //читаем данные
    /*
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) cin >> arr[i][j];
        //cin >> y[i];
    }
    */
    /*
    arr[0][0] = 6.25;
    arr[0][1] = -10;
    arr[0][2] = 0.5;
    arr[1][0] = -1;
    arr[1][1] = 5;
    arr[1][2] = 2.12;
    arr[2][0] = 0.5;
    arr[2][1] = 2.12;
    arr[2][2] = 3.6;
    y[0] = 7.5;
    y[1] = -8.68;
    y[2] = -0.24;
    */

    //ввод с экспонентой и альфой 
    
    double alpha = 1;
    double* tmp = new double[n];

    //с квадратом
    for (int i = 0; i < n; i++)
    {
        x[i] = i + 1;
        tmp[i] = x[i];
        for (int j = 0; j < n; j++)
        {
            arr[i][j] = exp(-pow(abs(i - j), 2) * alpha);
        }
    }
    for (int i = 0; i < n; i++)
    {
        y[i] = 0;
        for (int j = 0; j < n; j++)
        {
            y[i] += arr[i][j] * x[j];
        }
    }
    

    //без квадрата
    /*
    for (int i = 0; i < n; i++)
    {
        x[i] = i + 1;
        tmp[i] = x[i];
        for (int j = 0; j < n; j++)
        {
            arr[i][j] = exp(-abs(i - j) * alpha);
        }
    }
    for (int i = 0; i < n; i++)
    {
        y[i] = 0;
        for (int j = 0; j < n; j++)
        {
            y[i] += arr[i][j] * x[j];
        }
    }
    */

    sysOut(arr, n);
    yOut(y, n);
    relax(arr, x, y, n);
    xOut(x, n);
    cout << "--------\n";
    gauss(arr, x, y, n);
    xOut(x, n);

    //погрешность для метода с заданными иксами
    
    double pogr = 0;
    for (int i = 0; i < n; i++) pogr += abs(tmp[i] - x[i]) / abs(x[i]);
    pogr = pogr / n;
    cout << "alpha: " << alpha << "\n";
    cout << "относительная погрешность: (" << pogr * 100 << ")%\n\n";
    

    return 0;
}

void sysOut(double** arr, int n) {
    cout << "Левая часть: \n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (abs(arr[i][j]) < 1e-10) cout << "0.0 ";
            else cout << arr[i][j] << setprecision(4) << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}

void xOut(double* x, int n) {
    cout << "Решение: \n";
    for (int i = 0; i < n; i++)
    {
        cout << "x[" << i << "]=" << x[i] << setprecision(10) << "\n";
    }
    cout << "\n";
}

void yOut(double* y, int n) {
    cout << "Свободные члены: \n";
    for (int i = 0; i < n; i++)
    {
        cout << "y[" << i << "]=" << y[i] << setprecision(4) << "\n";
    }
    cout << "\n";
}

void relax(double** arr, double* x, double* y, int n) {
    double* xn = new double[n];
    for (int i = 0; i < n; i++) xn[i] = 0;
    double w = 1, eps = 0.001, delta = 0;
    int counter = 0;
    do
    {
        for (int j = 0; j < n; j++)
        {
            double sum1 = 0;
            for (int k = j + 1; k < n; k++) 
                sum1 += (arr[j][k] / arr[j][j]) * xn[k];

            double sum2 = 0;
            for (int k = 0; k < j; k++)
                sum2 += (arr[j][k] / arr[j][j]) * x[k];

            x[j] = (1 - w) * xn[j] - w * sum1 + w * y[j] / arr[j][j] - w * sum2;
        }
        counter++;
        delta = 0;
        // delta по значениям матрицы
        for (int i = 0; i < n; i++)
        {
            double tmp = 0;
            for (int j = 0; j < n; j++)
                tmp += x[j] * arr[i][j];
            if (delta < abs(y[i] - tmp))
                delta = abs(y[i] - tmp);
            xn[i] = x[i];
        }

        // delta по соседним итерациям
        /*
        for (int i = 0; i < n; i++)
        {
            if (delta < abs(x[i] - xn[i]))
                delta = abs(x[i] - xn[i]);
            xn[i] = x[i];
        }
        */
    } while (delta > eps);
    cout << counter << " итераций выполнено для получения точности " << eps << /*" (между итерациями)"<<*/"\n\n";
}

void gauss(double** arr, double* x, double* y, int n) {
    const double eps = 0.00001;  // точность определени¤ нул¤
    int i, j, k = 0;
    int index_max, det_coef = 1;

    // поиск строки с максимальным a[i][k]
    double max = abs(arr[k][k]);
    index_max = k;
    for (i = 0; i < n; i++)
    {
        if (abs(arr[i][k]) > max)
        {
            max = abs(arr[i][k]);
            index_max = i;
        }
    }

    // проверка на нулевой столбец
    if (max < eps)
    {
        cout << "Единственное решение получить невозможно из-за нулевого столбца " << k << " матрицы A" << endl;
        exit(1);
    }

    // перестановка строк
    if (k != index_max) {
        for (j = 0; j < n; j++)
        {
            //в матрице
            double temp1 = arr[k][j];
            arr[k][j] = arr[index_max][j];
            arr[index_max][j] = temp1;
        }
        //в правой части
        double temp3 = y[k];
        y[k] = y[index_max];
        y[index_max] = temp3;
    }

    for (k = 0; k < n; k++)
    {

        //занул¤ем вычитанием столбцы, кроме arr[k][k]
        for (i = 0;i < n;i++)
        {
            if (i == k) continue;
            double temp = (double)arr[i][k] / arr[k][k];
            for (j = 0; j < n; j++)
            {
                arr[i][j] -= arr[k][j] * (double)temp;
            }
            y[i] -= y[k] * (double)temp;
        }
    }

    //ищем иксы в уже диагональный матрице
    for (i = 0;i < n;i++)
        x[i] = y[i] / arr[i][i];
}