#include <iostream>
#include <cmath>

using namespace std;

double gauss(double** arr, double* x, double* y, int n, double** rev)
{
    const double eps = 0.00001;  // �������� ����������� ����
    double det = 1;
    int max, index_max, coef = 1;
    for (int k = 0; k < n; k++)
    {
        // ����� ������ � ������������ a[i][k]
        max = abs(arr[k][k]);
        index_max = k;
        for (int i = 0; i < n; i++)
        {
            if (abs(arr[i][k]) > max)
            {
                max = abs(arr[i][k]);
                index_max = i;
            }
        }
        //�������� ��������� �������
        for (int i = 0;i < n;i++)
        {
            rev[i] = new double[n];
            rev[i][i] = 1;
            for (int j = 0;j < n && j != i;j++) rev[i][j] = 0;
        }
        // �������� �� ������� �������
        if (max < eps)
        {
            cout << "������������ ������� �������� ���������� ��-�� �������� ������� " << index_max << " ������� A" << endl;
            exit(1);
        }
        // ������������ �����
        if (k != index_max) {
            for (int j = 0; j < n; j++)
            {
                //� �������
                double temp1 = arr[k][j];
                arr[k][j] = arr[index_max][j];
                arr[index_max][j] = temp1;
                //� ��������� �������
                double temp2 = rev[k][j];
                rev[k][j] = rev[index_max][j];
                rev[index_max][j] = temp2;
            }
            //� ������ �����
            double temp3 = y[k];
            y[k] = y[index_max];
            y[index_max] = temp3;
            //������ ���� ������������
            coef *= -1;
        }
        //�������� �������, ����� arr[k][k]
        for (int i = 0;i < n;i++)
        {
            if (i == k) continue;
            double temp = (double)arr[i][k] / arr[k][k];
            for (int j = 0; j < n; j++)
            {
                arr[i][j] -= arr[k][j] * (double)temp;
                rev[i][j] -= rev[k][j] * (double)temp;
            }
            y[i] -= y[k] * (double)temp;
        }
        //���� ���� � ��� ������������ �������
        for (int i = 0;i < n;i++)
        {
            x[i] = y[i] / arr[i][i];
        }
        //��������� ����� �������, ������� �� �������� �����
        for (int i = 0; i < n; i++)
        {
            double temp = arr[i][i];
            for (int j = 0; j < n; j++)
            {
                rev[i][j] /= temp;
            }
        }
    }
    //������ ������������
    for (int i = 0;i < n;i++)
    {
        det *= (double)arr[i][i];
    }
    det *= coef;
    return det;
}


int main()
{
    setlocale(LC_ALL, "Russian");

    int n = 3;
    double** arr = new double* [n];
    double** rev = new double* [n];
    double* x = new double[n];
    double* y = new double[n];

    //��������
    for (int i = 0;i < n;i++)
    {
        arr[i] = new double[n];
        y[i] = -1;
        for (int j = 0;j < n;j++)
        {
            if (i < j)
            {
                arr[i][j] = -1;
            }
            else
            {
                if (i == j)
                    arr[i][j] = 1;
                else
                    arr[i][j] = 0;
            }
        }
        if (i==n-1) y[n - 1] = 1;
    }

    //������ ����
    arr[0][2] = (double)(-1) / 3;
    arr[1][1] = 2;
    arr[1][2] = (double)(2) / 7;
    arr[2][2] = 7;
    arr[1][0] = 0.5;
    //������ ����

    //����� �������� �������
    cout << "�������� �������: \n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << arr[i][j] << " ";
        }
        cout << "= " << y[i] << "\n";
    }
    cout << "\n";

    double det = gauss(arr, x, y, n, rev);

    //����� ������������ �������
    cout << "������������ �������: \n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << arr[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
    
    //����� ������� �����
    cout << "�������� �������: \n";
    for (int i = 0; i < n; i++)
        cout << "x[" << i << "]=" << x[i] << ' ';
    cout << "\n\n";

    //����� ������������
    cout << "������������ ����� " << det << endl;
    cout << "\n";

    //�������� �������
    cout << "�������� �������: \n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << rev[i][j] << " ";
        }
        cout << "\n";
    }

    return 0;
}

    //����
   /*
   double alpha = 0.1;
   for (int i = 0;i < n;i++)
   {
       arr[i] = new double[n];
       for (int j = 0;j < n;j++)
       {
           y[i] = 1;
           arr[i][j] = exp(-alpha * (i - j) * (i - j));
       }
   }
   */

   //������ ����������
   /*for (int i = 0; i < n; i++)
   {
       arr[i] = new double[n];
       for (int j = 0; j < n; j++)
       {
           cout << "arr[" << i << "][" << j << "]= ";
           cin >> arr[i][j];
       }
   }

   cout << "\n";
   for (int i = 0; i < n; i++)
   {
       cout << "y[" << i << "]= ";
       cin >> y[i];
   }*/

   //�� ������ ��������
     /*
     for (int i = 0;i < n;i++)
     {
         arr[i] = new double[n];
         for (int j = 0;j < n;j++)
         {
             arr[i][j] = 1 / (i + j + 1);
         }
     }

     y[0] = 137 / 60;
     y[1] = 29 / 20;
     y[2] = 153 / 140;
     y[3] = 743 / 840;
     y[4] = 1879 / 2520;
     */