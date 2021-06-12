#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

double gauss_det(double** arr, double** rev, double* x, double* y, int n);

int main() {
	setlocale(LC_ALL, "Russian");

    int n = 4;
    int i, j, k;
    double** arr = new double* [n];
    double** rev = new double* [n];
    double* x = new double[n];
    double* y = new double[n];

    //создание единичной матрицы для rev
    for (i = 0;i < n;i++)
    {
        rev[i] = new double[n];
        for (j = 0;j < n;j++)
        {
            if (i==j) rev[i][j] = 1;
            else rev[i][j] = 0;
        }
    }

    //шаблоны исходной системы
    /*

    

    //2
    for (i = 0;i < n;i++)
    {
        arr[i] = new double[n];
        for (j = 0; j < n; j++)
            arr[i][j] = 1;
        y[i] = 1;
    }


    //3
    for (i = 0;i < n;i++)
    {
        arr[i] = new double[n];
        for (j = 0; j < n; j++)
            arr[i][j] = 0;
        y[i] = 1;
        arr[i][n - i - 1] = 1;
    }

    //4
    for (i = 0;i < n;i++)
    {
        arr[i] = new double[n];
        for (j = 0; j < n; j++)
            arr[i][j] = pow(2,-abs(i-j));
        y[i] = 1;
    }

    */  

    //1
    for (i = 0;i < n;i++)
    {
        arr[i] = new double[n];
        for (j = 0; j < n; j++)
            arr[i][j] = (double)min(i + 1, j + 1) / max(i + 1, j + 1);
        y[i] = 1;
    }
    
    //вывод исходной системы
    cout << "Исходная система: \n";
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
            cout << arr[i][j] << " ";
        cout << "= " << y[i] << "\n";
    }
    cout << "\n";

    //вывод опеределителя
    cout << "Определитель: \n";
    cout << gauss_det(arr, rev, x, y, n) << "\n\n";

    //вывод диагонального вида матрицы
    cout << "Диагональный вид: \n";
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
            cout << arr[i][j] << " ";
        cout << "\n";
    }
    cout << "\n";

    //вывод обратной матрицы
    cout << "Обратная матрица: \n";
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
            cout << rev[i][j] << " ";
        cout << "\n";
    }
    cout << "\n";

    return 0;
}

double gauss_det(double** arr, double** rev, double* x, double* y, int n) {
    const double eps = 0.00001;  // точность определения нуля
    double det = 1;
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
                //в единичной матрице
                double temp2 = rev[k][j];
                rev[k][j] = rev[index_max][j];
                rev[index_max][j] = temp2;
            }
            //в правой части
            double temp3 = y[k];
            y[k] = y[index_max];
            y[index_max] = temp3;
            //меняем знак определителя
            det_coef *= -1;
        }

        //зануляем вычитанием столбцы, кроме arr[k][k]
        for (i = 0;i < n;i++)
        {
            if (i == k) continue;
            double temp = (double)arr[i][k] / arr[k][k];
            for (j = 0; j < n; j++)
            {
                arr[i][j] -= arr[k][j] * (double)temp;
                rev[i][j] -= rev[k][j] * (double)temp;
            }
            y[i] -= y[k] * (double)temp;
        }
    }

    //ищем иксы в уже диагональный матрице
    for (i = 0;i < n;i++)
        x[i] = y[i] / arr[i][i];

    //ищем значение rev матрицы при нормирвоании arr
    for (int i = 0; i < n; i++)
    {
        double temp = arr[i][i];
        for (j = 0; j < n; j++)
            rev[i][j] /= (double)temp;
    }

    //найдем определитель
    for (i = 0;i < n;i++)
        det *= (double)arr[i][i];
    det *= det_coef;

    return det;
}

/*
 cout << "Исходные иксы: \n";
 for (int i = 0; i < n; i++)
 {
     x[i] = i;
     cout << "x[" << i << "]=" << x[i] << '\n';
 }
 cout << "\n";

 for (int i = 0; i < n;i++)
     temp[i] = x[i];



 cout << "Игреки через иксы: \n";
 for (int i = 0; i < n; i++)
 {
     y[i] = 0;
     for (int j = 0; j < n;j++)
     {
         y[i] += x[j] * arr[i][j];
     }
     cout << "y[" << i << "]=" << y[i] << '\n';
 }
 cout << "\n";

//---------------------------------------------------//

  double pogr = 0;
    for (int i = 0; i < n; i++)
    {
        pogr += abs(temp[i] - x[i]);
    }
    pogr /= (double)n;
    cout << "относитальная погрешность = "<<pogr<<" \n";
*/