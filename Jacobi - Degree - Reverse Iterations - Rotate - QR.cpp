#define _CRT_SECURE_NO_WRARNINGS
//#include <stdio.h>
#include <iostream>
#include <cmath>
#include <iomanip> 

using namespace std;

void jacobi(double** arr, int n);
void degree_method(double** arr, int n);
void reverse_iterations(double** arr, int n);
void qr(double** arr, int n);
void rotate(double** arr, double* x, int n);
void inversion(double** arr, int N);
void xOut(double* x, int n);
double norma_(double* x, int n);
void arrOut(double** arr, int n);
void arrIn(double** arr, int n)  {
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) cin >> arr[i][j];
}

int main(){
	setlocale(LC_ALL, "Russian");

    int n = 3;
    double** arr = new double* [n];
    for (int i = 0; i < n; i++) arr[i] = new double[n];

    //arrIn(arr, n);

    arr[0][0] = 5;
    arr[0][1] = 1;
    arr[0][2] = 2;
    arr[1][0] = 1;
    arr[1][1] = 4;
    arr[1][2] = 1;
    arr[2][0] = 2;
    arr[2][1] = 1;
    arr[2][2] = 3;

    /*for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            arr[i][j] = 0;
            if (i == j) arr[i][j] = 5 * i;
        }
    }
    arr[0][0] = 1;
    arr[0][2] = arr[2][0] = 0.1;
    arr[1][3] = arr[3][1] = -1.2;
    arr[2][4] = arr[4][2] = 2;
    arr[3][5] = arr[5][3] = -2.7;
    arr[4][6] = arr[6][4] = 3.4;
    arr[5][7] = arr[7][5] = 4.9;
    arr[6][8] = arr[8][6] = -5;
    arr[7][9] = arr[9][7] = 6.8;*/


    /*arr[0][0] = 2;
    arr[0][1] = 1;
    arr[1][0] = 1;
    arr[1][1] = 3;*/

    arrOut(arr, n);
    jacobi(arr, n);
    //degree_method(arr, n);
    reverse_iterations(arr, n);
    qr(arr, n);

    //ввод с экспонентой и альфой 
    /*
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
    */

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
    *

    //погрешность для метода с заданными иксами
    /*
    double pogr = 0;
    for (int i = 0; i < n; i++) pogr += abs(tmp[i] - x[i]) / abs(x[i]);
    pogr = pogr / n;
    cout << "alpha: " << alpha << "\n";
    cout << "относительная погрешность: (" << pogr * 100 << ")%\n\n";
    */

	return 0;
};

void jacobi(double** arr, int n) {

    double* owns = new double[n];
    double** vectors = new double* [n];
    double** vectors_tmp = new double* [n];
    double** h = new double* [n];
    double** h_t = new double* [n];
    double** a_tmp = new double* [n];
    double** a = new double* [n];
    
    for (int i = 0; i < n; i++)
    {
        a[i] = new double[n];
        a_tmp[i] = new double[n];
        h[i] = new double[n];
        h_t[i] = new double[n];
        vectors[i] = new double[n];
        vectors_tmp[i] = new double[n];
        for (int j = 0; j < n; j++) a[i][j] = arr[i][j];
    }

    int i_max, j_max;
    int counter = 0;
    double a_max = 0;
    double phi;
    while (1)
    {
        /*cout << "-------- Итерация " << counter << " --------\n";
        cout << "a(" << counter << "):\n";
        arrOut(a, n);*/
        //ищем максимальный внедиагональный элемент
        a_max = 0;
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                if (abs(a_max) < abs(a[i][j]))
                {
                    a_max = a[i][j];
                    i_max = i;
                    j_max = j;
                }
            }
        }
        /*cout << "max el = " << a[i_max][j_max] <<"\n\n";*/
        //если он меньше eps завершаем метод
        if (abs(a[i_max][j_max]) < 0.00001)
        {
            for (int i = 0; i < n; i++) owns[i] = a[i][i];
            
            //умножим исходную матрицу на матрицу собственных векторов
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    h[i][j] = 0;
                    for (int k = 0; k < n; k++) h[i][j] += arr[i][k] * vectors[k][j];
                }
            }
            
            cout << "Метод Якоби (" << counter + 1 << " итераций):\n";
            for (int j = 0; j < n;j++)
            {
                cout << "Для вектора " << j+1 << " с собств. числом "<< owns[j] <<" ( ";
                for (int i = 0; i < n;i++)
                {
                    cout << vectors[i][j] << " ";
                }
                cout << ")\n";

                cout << " - Вектор, умноженый на матрицу: ";
                for (int i = 0; i < n;i++)
                {
                    cout << h[i][j] << " ";
                }
                cout << "\n";
                cout << " - Затем деленный на собств число: ";
                for (int i = 0; i < n;i++)
                {
                    cout << h[i][j]/owns[j] << " ";
                }
                cout << "\n";
                cout << "\n";
            }
            cout << "\n";
            break;
        }

        //найдем угол поворота и матрицу h
        phi = 0.5 * atan(2 * a_max / (a[i_max][i_max] - a[j_max][j_max]));
        /*cout << "phi = " << phi << "\n\n";*/
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i != j) h[i][j] = h_t[i][j] = 0;
                else h[i][j] = h_t[i][j] = 1;
            }
        }
        h[i_max][i_max] = h_t[i_max][i_max] = h[j_max][j_max] = h_t[j_max][j_max] = cos(phi);
        h[i_max][j_max] = h_t[j_max][i_max] = -sin(phi);
        h_t[i_max][j_max] = h[j_max][i_max] = sin(phi);
        /*cout << "h("<<counter<<"):\n";
        arrOut(h, n);*/

        //умножим матрицы a = h_t * a * h = h_t * a_tmp;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                a_tmp[i][j] = 0;
                for (int k = 0; k < n; k++) a_tmp[i][j] += a[i][k] * h[k][j];
            }
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                a[i][j] = 0;
                for (int k = 0; k < n; k++) a[i][j] += h_t[i][k] * a_tmp[k][j];
            }
        }

        //считаем матрицы векторов как произведение vectors = h0 * h1 * ... * hn = vectors * hn
        if (counter == 0)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++) vectors[i][j] = h[i][j];
            }
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++) vectors_tmp[i][j] = vectors[i][j];
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    vectors[i][j] = 0;
                    for (int k = 0; k < n; k++) vectors[i][j] += vectors_tmp[i][k] * h[k][j];
                }
            }
        }

        /*cout << "vectors(" << counter << "):\n";
        arrOut(vectors, n);*/

        counter++;
    }
}

void degree_method(double** arr, int n) {

    /*Нормально работает только для симметричных матриц. В принципе может использоваться 
    для матриц общего вида, но если у такой матрицы собственные числа будут комплексными
    (в большинстве случаев) сходиться не будет. кроме того, метод находит максимальное по 
    модулю число, поэтому в случае отрицательных собственных чисел может сойтись не к тому 
    числу. В случае вещественных симметричных матриц собственные числа положительные 
    вещественные, то есть таких проблем нет и метод всегда сходится. Но скорость сходимости 
    зависит от отношения максимального собственного числа и второго по величине. В случае 
    если они будут почти равны,сходимость будет медленной и может не хватить итераций.*/
    double own = 0, own_n, delta;
    int counter = 0;
    double* x = new double[n];
    double* x_n = new double[n];
    for (int i = 0; i < n; i++) x[i] = 1/sqrt((double)n);

    do 
    {
        for (int i = 0; i < n; i++) x_n[i] = 0;

        //посчитали xn для следующей итерации
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                x_n[i] += arr[i][j] * x[j];
            }
        }
        /* double norma = 0;
        for (int i = 0; i < n; i++) norma += pow(x_n[i], 2);
        norma = sqrt(norma);
        for (int i = 0; i < n; i++) x_n[i] = x_n[i] / norma;*/

        //считаем own для этой итерации
        double numerator = 0, denomerator = 0;
        for (int i = 0; i < n; i++)
        {
            numerator += x_n[i] * x[i];
            denomerator += x[i] * x[i];
        }
        own_n = numerator / denomerator;

        //считаем приближеник к own и переопределяем x, own
        delta = abs(own - own_n);
        own = own_n;
        for (int i = 0; i < n; i++) x[i] = x_n[i];

        counter++;

    } while (delta > 0.00001);

    cout << "Степенной метод ("<< counter <<" итераций):\n";
    cout << "Наибольшее собственно число = " << own << "\n\n";
}

void inversion(double** arr, int N)
{
    double temp;

    double** E = new double* [N];

    for (int i = 0; i < N; i++)
        E[i] = new double[N];

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            E[i][j] = 0.0;

            if (i == j)
                E[i][j] = 1.0;
        }

    for (int k = 0; k < N; k++)
    {
        temp = arr[k][k];

        for (int j = 0; j < N; j++)
        {
            arr[k][j] /= temp;
            E[k][j] /= temp;
        }

        for (int i = k + 1; i < N; i++)
        {
            temp = arr[i][k];

            for (int j = 0; j < N; j++)
            {
                arr[i][j] -= arr[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int k = N - 1; k > 0; k--)
    {
        for (int i = k - 1; i >= 0; i--)
        {
            temp = arr[i][k];

            for (int j = 0; j < N; j++)
            {
                arr[i][j] -= arr[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            arr[i][j] = E[i][j];

    for (int i = 0; i < N; i++)
        delete[] E[i];

    delete[] E;
}

void reverse_iterations(double** arr, int n) {

    double** a = new double* [n];
    for (int i = 0; i < n; i++)
    {
        a[i] = new double[n];
        for (int j = 0; j < n; j++) a[i][j] = arr[i][j];
    }

    int counter = 0;
    double* x = new double[n];
    double* x_n = new double[n];

    double own = 0, own_n;
    double norma, delta;

    for (int i = 0; i < n; i++) x[i] = 1 / sqrt((double)n); //n-размерность матрицы
    for (int i = 0; i < n; i++) a[i][i] -= own;
    inversion(a, n);

    do
    {
        for (int i = 0; i < n; i++) x_n[i] = 0;

        //посчитали x_n для следующей итерации
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                x_n[i] += a[i][j] * x[j]; //x(n+1)=(A-own*E)^(-1)*x(n)
            }
        }

        norma = norma_(x_n, n);
        //for (int i = 0; i < n; i++) x_n[i] = x_n[i] / norma;

        //считаем own для этой итерации (расчет собственного числа по методы 3.187\3.188)
        double numerator = 0, denomerator = 0;
        for (int i = 0; i < n; i++)
        {
            numerator += x[i] * x_n[i];
            denomerator += x_n[i] * x_n[i];
        }
        own_n = numerator / denomerator;//приближенное собственное значение

        delta = abs(own - own_n);
        //cout << counter << ": " << own << "\n";
        own = own_n;
        for (int i = 0; i < n; i++) x[i] = x_n[i];

        counter++;

    } while (delta > 0.00001);
    cout << "Метод обратных итераций (" << counter << " итераций):\n";
    cout << "Собсвтенное число с наименьшим модулем = " << own << "\n";
    cout << "Отвечающий ему вектор ( ";

    norma = norma_(x, n);
    for (int i = 0; i < n; i++) cout << x[i] << " ";
    cout << ")\n\n";

}

double norma_(double* x, int n) {
    double norma = 0;
    for (int i = 0;i < n;i++) norma += pow(x[i], 2);
    norma = sqrt(norma);
    return norma;
}

void rotate(double** in, double** arr, double* x, double* y, int n) {
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) arr[i][j] = in[i][j];
    }

    double c, s;
    for (int b = 0; b < n - 1; b++)
    {
        for (int k = b + 1; k < n; k++)
        {
            c = arr[b][b] / pow(pow(arr[b][b], 2) + pow(arr[k][b], 2), 0.5);
            s = arr[k][b] / pow(pow(arr[b][b], 2) + pow(arr[k][b], 2), 0.5);

            //временные копии строк
            double* tmp_str_0 = new double[n];
            double* tmp_str_k = new double[n];
            for (int i = b;i < n;i++)
            {
                tmp_str_0[i] = arr[b][i];
                tmp_str_k[i] = arr[k][i];
            }
            double tmp_y_0 = y[b];
            double tmp_y_k = y[k];

            for (int i = b;i < n;i++)
            {
                arr[b][i] = c * tmp_str_0[i] + s * tmp_str_k[i];
                y[b] = c * tmp_y_0 + s * tmp_y_k;

                arr[k][i] = -s * tmp_str_0[i] + c * tmp_str_k[i];
                y[k] = -s * tmp_y_0 + c * tmp_y_k;
            }
        }
    }
    for (int k = n - 1; k > -1; k--)
    {
        double sum = 0;
        for (int i = k + 1; i < n; i++)
        {
            sum += arr[k][i] * x[i];
        }
        x[k] = (y[k] - sum) / arr[k][k];
    }
}

void qr(double** arr, int n) {
    double** a = new double* [n];
    double** L = new double* [n];
    double** U = new double* [n];
    double** tmp1 = new double* [n];
    double** vectors = new double* [n];
    double** vectors_prev = new double* [n];
    int counter = 0;
    double delta = 0;
    for (int i = 0; i < n; i++)
    {
        a[i] = new double[n];
        L[i] = new double[n];
        U[i] = new double[n];
        tmp1[i] = new double[n];
        vectors[i] = new double[n];
        vectors_prev[i] = new double[n];
        for (int j = 0; j < n; j++) a[i][j] = arr[i][j];
    }
    double* kos1 = new double[n];
    double* kos2 = new double[n];
    for (int i = 0; i < n; i++) kos1[i] = kos2[i] = 1;
    do {
        rotate(a, L, kos1, kos2, n);
        inversion(L, n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++) U[i][j]=0;
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++) 
            {
                for (int k = 0; k < n; k++) U[i][j] += a[i][k] * L[k][j];
            }
        }


        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++) tmp1[i][j] = 0;
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                for (int k = 0; k < n; k++) tmp1[i][j] += a[i][k] * U[k][j];
            }
        }
        inversion(U, n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++) a[i][j] = 0;
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                for (int k = 0; k < n; k++) a[i][j] += U[i][k] * tmp1[k][j];
            }
        }
        if (counter == 0)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++) vectors_prev[i][j] = vectors[i][j] = U[i][j];
            }
        }
        else 
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++) vectors[i][j] = 0;
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    for (int k = 0; k < n; k++) vectors[i][j] += U[i][k] * vectors_prev[k][j];
                }
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++) vectors_prev[i][j] = vectors[i][j];
            }
        }
        delta = 0;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++) 
            {
                if (i != j && U[i][j] > delta) delta = U[i][j];
            }
        }
        counter++;
    } while (delta > 0.00001);
    cout << "Метод QR-разложения (" << counter + 1 << " итераций ):\n";
    for (int j = 0; j < n; j++)
    {
        cout << "Вектору ( ";
        for (int i = 0; i < n; i++)
        {
            if (abs(vectors[i][j]) < 1e-10) cout << "0.0 ";
            else cout << vectors[i][j] << setprecision(4) << " ";
        }
        cout << ") соответствует собсв. число " << a[j][j] << "\n";
    }
    cout << "\n";
    /*arrOut(vectors, n);
    cout << "\n";
    arrOut(a, n);*/

}

//void rotate(double** arr, double* x, int n) {
//
//    double** a = new double* [n];
//    double** U = new double* [n];
//    double** u_n = new double* [n];
//    double** u_tmp = new double* [n];
//
//    double** U1 = new double* [n];
//    double** u_n1 = new double* [n];
//    double** u_tmp1 = new double* [n];
//    for (int i = 0; i < n; i++)
//    {
//        a[i] = new double[n];
//        U[i] = new double[n];
//        u_n[i] = new double[n];
//        u_tmp[i] = new double[n];
//        /*U1[i] = new double[n];
//        u_n1[i] = new double[n];
//        u_tmp1[i] = new double[n];*/
//        for (int j = 0; j < n; j++) 
//        {
//            a[i][j] = arr[i][j];
//            if (i == j) U[i][j] = 1;
//            else U[i][j] = 0;
//        }
//    }
//    //true_inv(U, U1, n);
//
//    double c, s;
//    for (int b = 0; b < n - 1; b++)
//    {
//        for (int k = b + 1; k < n; k++)
//        {
//            c = a[b][b] / pow(pow(a[b][b], 2) + pow(a[k][b], 2), 0.5);
//            s = a[k][b] / pow(pow(a[b][b], 2) + pow(a[k][b], 2), 0.5);
//
//            //u_n = ...
//            for (int i = 0; i < n; i++)
//            {
//                for (int j = 0; j < n; j++) {
//                    if (i == j) u_n[i][j] = 1;
//                    else u_n[i][j] = 0;
//                }
//            }
//            u_n[b][b] = c;
//            u_n[k][k] = c;
//            u_n[k][b] = -s;
//            u_n[b][k] = s;
//
//            //true_inv(u_n, u_n1, n);
//
//            ////u_tmp1 = U1
//            //for (int i = 0; i < n; i++)
//            //{
//            //    for (int j = 0; j < n; j++) {
//            //        u_tmp1[i][j] = U1[i][j];
//            //        U1[i][j] = 0;
//            //    }
//            //}
//            ////U1 = u_tmp1 * u_n1
//            //for (int i = 0; i < n; i++)
//            //{
//            //    for (int j = 0; j < n; j++) {
//            //        for (int k = 0; k < n;k++) U1[i][j] += u_tmp1[i][k] * u_n1[k][j];
//            //    }
//            //}
//
//            //u_tmp = U
//            for (int i = 0; i < n; i++)
//            {
//                for (int j = 0; j < n; j++) {
//                    u_tmp[i][j] = U[i][j];
//                    U[i][j] = 0;
//                }
//            }
//            //U = u_n * u_tmp
//            for (int i = 0; i < n; i++)
//            {
//                for (int j = 0; j < n; j++) {
//                    for (int k = 0; k < n;k++) U[i][j] += u_n[i][k] * u_tmp[k][j];
//                }
//            }
//
//            //временные копии строк
//            double* tmp_str_0 = new double[n];
//            double* tmp_str_k = new double[n];
//            for (int i = b;i < n;i++)
//            {
//                tmp_str_0[i] = a[b][i];
//                tmp_str_k[i] = a[k][i];
//            }
//
//            for (int i = b;i < n;i++)
//            {
//                a[b][i] = c * tmp_str_0[i] + s * tmp_str_k[i];
//
//                a[k][i] = -s * tmp_str_0[i] + c * tmp_str_k[i];
//            }
//        }
//    }
//    for (int k = n - 1; k > -1; k--)
//    {
//        double sum = 0;
//        for (int i = k + 1; i < n; i++)
//        {
//            sum += a[k][i] * x[i];
//        }
//    }
//    for (int i = 0; i < n; i++)
//    {
//        for (int j = 0; j < n; j++) {
//            arr[i][j] = U[i][j];
//        }
//    }
//    /*cout << "\nQ:\n";
//    arrOut(U, n);
//    cout << "\nR:\n";
//    arrOut(a, n);*/
//
//    //double** tmp2 = new double* [n]; // = A' * U1
//    //double** tmp1 = new double* [n]; // = U * (A' * U1)
//    //for (int i = 0; i < n; i++)
//    //{
//    //    tmp1[i] = new double[n];
//    //    tmp2[i] = new double[n];
//    //    for (int j = 0; j < n; j++) 
//    //    {
//    //        tmp1[i][j] = 0;
//    //        tmp2[i][j] = 0;
//    //    }
//    //}
//    //for (int i = 0; i < n; i++)
//    //{
//    //    for (int j = 0; j < n; j++) {
//    //        for (int k = 0; k < n;k++) tmp2[i][j] += a[i][k] * U1[k][j];
//    //    }
//    //}
//    //for (int i = 0; i < n; i++)
//    //{
//    //    for (int j = 0; j < n; j++) {
//    //        for (int k = 0; k < n;k++) tmp1[i][j] += U[i][k] * tmp2[k][j];
//    //    }
//    //}
//
//    //cout << "\nQ*R:\n";
//    //arrOut(tmp1, n);
//
//}

void arrOut(double** arr, int n) {
    //cout << "Матрица: \n";
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
