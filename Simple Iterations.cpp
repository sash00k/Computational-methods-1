
#define _CRT_SECURE_nO_WARnInGS
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace std;
void iterations(double** arr, double* x, double* y, int n);
void gauss(double** arr, double* x, double* y, int n);
void yOut(double* y, int n);
void xOut(double* x, int n);
void sysOut(double** arr, int n);

void main() {
    setlocale(LC_ALL, "Russian");

    int n = 3;

    double** arr = new double* [n];
    for (int i = 0; i < n; i++) arr[i] = new double[n];
    double* y = new double[n];
    double* x = new double[n];
    for (int i = 0; i < n; i++) x[i] = 0;


    //������ ������
    /*
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) cin >> arr[i][j];
        //cin >> y[i];
    }
    */

    arr[0][0] = 6.25;
    arr[0][1] = -1;
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

    double* tmp = new double[n];
    for (int i = 0; i < n; i++) tmp[i] = x[i];

    //���� � ����������� � ������ 
    /*
    double alpha = 1;
    double* tmp = new double[n];

    //� ���������
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
    */

    //��� ��������
    /*�
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

    xOut(x, n);
    iterations(arr, x, y, n);
    sysOut(arr, n);
    yOut(y, n);
    xOut(x, n);

    //����������� ��� ������ � ��������� ������
    /*
    double pogr = 0;
    for (int i = 0; i < n; i++) pogr += abs(tmp[i] - x[i]) / abs(x[i]);
    pogr = pogr / n;
    cout << "alpha: " << alpha << "\n";
    cout << "������������� �����������: (" << pogr * 100 << ")%\n\n";
    */
}

void sysOut(double** arr, int n) {
    cout << "����� �����: \n";
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
    cout << "�������: \n";
    for (int i = 0; i < n; i++)
    {
        cout << "x[" << i << "]=" << x[i] << setprecision(4) << "\n";
    }
    cout << "\n";
}

void yOut(double* y, int n) {
    cout << "��������� �����: \n";
    for (int i = 0; i < n; i++)
    {
        cout << "y[" << i << "]=" << y[i] << setprecision(4) << "\n";
    }
    cout << "\n";
}

void gauss(double** arr, double* x, double* y, int n) {
    const double eps = 0.00001;  // �������� ���������� ���
    int i, j, k = 0;
    int index_max, det_coef = 1;

    // ����� ������ � ������������ a[i][k]
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

    // �������� �� ������� �������
    if (max < eps)
    {
        cout << "������������ ������� �������� ���������� ��-�� �������� ������� " << k << " ������� A" << endl;
        exit(1);
    }

    // ������������ �����
    if (k != index_max) {
        for (j = 0; j < n; j++)
        {
            //� �������
            double temp1 = arr[k][j];
            arr[k][j] = arr[index_max][j];
            arr[index_max][j] = temp1;
        }
        //� ������ �����
        double temp3 = y[k];
        y[k] = y[index_max];
        y[index_max] = temp3;
    }

    for (k = 0; k < n; k++)
    {

        //������� ���������� �������, ����� arr[k][k]
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

    //���� ���� � ��� ������������ �������
    for (i = 0;i < n;i++)
        x[i] = y[i] / arr[i][i];
}

void iterations(double** arr, double* x, double* y, int n) {

    double norma, c = 0, s = 0;
    double* xn = new double[n];
    double eps = 0.0001;

    double** r = new double* [n];
    for (int i = 0; i < n; i++) r[i] = new double[n];
    double* yr = new double[n];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) r[i][j] = 0;
    }

    do
    {
        norma = 0.0;
        for (int i = 0; i < n; i++)
        {
            xn[i] = -y[i];
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                    xn[i] += arr[i][j] * x[j];
            }
            xn[i] /= (-arr[i][i]);
        }
        for (int i = 0; i < n; i++)
        {
            if (fabs(x[i] - xn[i]) > norma)
                norma = fabs(x[i] - xn[i]);
            x[i] = xn[i];
        }
    } while (norma > eps);
}