#include <iostream>
#include <cmath>
#include <algorithm>

#include <stdio.h>
#include <locale>
#include <iostream>

using namespace std;

double khaletsky(double** arr, double** l1, double** l2, double* y, double* x, int n);
double gauss(double** arr, double* x, double* y, int n);

int main() {
	setlocale(LC_ALL, "Russian");

	int n = 3;
	int i, j, k;
	double** arr = new double* [n];
	double** l1 = new double* [n];
	double** l2 = new double* [n];
	double* x = new double[n];
	double* y = new double[n];

	//test
	for (i = 0;i < n;i++)
	{
		l1[i] = new double[n];
		l2[i] = new double[n];
		arr[i] = new double[n];
		for (j = 0; j < n; j++)
		{
			l1[i][j] = 0;
			l2[i][j] = 0;
		}
	}
	arr[0][0] = 81.0;
	arr[0][1] = -45.0;
	arr[0][2] = 45.0;
	arr[1][0] = -45.0;
	arr[1][1] = 50.0;
	arr[1][2] = -15.0;
	arr[2][0] = 45.0;
	arr[2][1] = -15.0;
	arr[2][2] = 38.0;
	y[0] = 531.0;
	y[1] = -460.0;
	y[2] = 193.0;

	//вывод исходной системы
	cout << "Исходная система: \n";
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			cout << arr[i][j] << " ";
		cout << "= " << y[i] << "\n";
	}
	cout << "\n";

	khaletsky(arr, l1, l2, y, x, n);

	//вывод ответа
	cout << "Решение: \n";
	for (i = 0; i < n; i++)
	{
		cout << "x[" << i << "]=" << x[i] << "\n";
	}
	cout << "\n";

	return 0;
}

double khaletsky(double** arr, double** l1, double** l2, double* y, double* x, int n)
{
	int i, j, k;
	double* tmp_x = new double[n];
	for (k = 0; k < n; k++)
	{
		double sum1 = 0;
		//считаем элементы на диагонали
		for (i = 0; i < k; i++)
			sum1 += pow(l1[k][i], 2);
		l1[k][k] = sqrt(arr[k][k]-sum1);
		l2[k][k] = l1[k][k];

		//считаем поддиагональные элементы
		for (i = 0; i < n; i++)
		{
			double sum2 = 0;
			for (j = 0; j < k; j++)
			{
				sum2 += l1[i][j] * l1[k][j];
			}
			l1[i][k] = (double)(arr[i][k] - sum2) / l1[k][k];
			l2[k][i] = l1[i][k];
		}

	}

	//вывод разложения
	cout << "L1: \n";
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			cout << l1[i][j] << " ";
		cout << "\n";
	}
	cout << "\n";
	cout << "L2: \n";
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			cout << l2[i][j] << " ";
		cout << "\n";
	}
	cout << "\n";

	gauss(l1, tmp_x, y, n);
	gauss(l2, x, tmp_x, n);

	return 0;
}

double gauss(double** arr, double* x, double* y, int n) {
	const double eps = 0.00001;  // точность определени¤ нул¤
	int i, j, k;
	int index_max, det_coef = 1;
	for (k = 0; k < n; k++)
	{
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

	return 0;
}