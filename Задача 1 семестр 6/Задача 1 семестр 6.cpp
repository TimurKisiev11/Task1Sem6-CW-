#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <windows.h>
#include <iomanip>
#include <cstdlib>
using namespace std;

void cop(double** a, double** b, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            a[i][j] = b[i][j];
}

double** get_null_matrix(int n)
{
    n = n + 1;
    double** A = new double* [2 * n];
    for (int i = 0; i <= n; i++)
    {
        A[i] = new double[2 * n];
        for (int j = 0; j <= n; j++)
            A[i][j] = 0.0;
    }
    return A;
}
void show(double** A, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << fixed << setprecision(2) << " \t" << A[i][j] << " \t";
        }
        cout << endl;
    }
}

void show_vect(double* v, int n)
{
    for (int i = 0; i < n; i++)
    {

        cout << fixed << setprecision(2) << " \t" << v[i] << " \t";
        cout << endl;
    }
}

void sysout(double** a, double* y, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << "(" << a[i][j] << ")" << " * x " << j;
            if (j < n)
                cout << " + ";
        }
        cout << " = " << y[i] << endl;
    }
    return;
}

//коэффициенты диф.ур.
double P(double x)
{
    return (-1) / (x - 3);
}
double Q(double x)
{
    return (1 + x / 2);
}
double R(double x)
{
    return exp(x / 2);
}
double F(double x)
{
    return (2 - x);
}
//для решения системы с матрицей "a" вектором правых частей "y" и порядком "n"
double* gauss2(double** a, double* y, int n)
{
    double* x, max;
    int k, index;
    const double eps = 0.00001;  // точность
    x = new double[n];
    k = 0;
    for (k = 0; k < n; k++)
    {
        // Поиск строки с максимальным a[i][k]
        max = fabs(a[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++)
        {
            if (fabs(a[i][k]) > max)
            {
                max = fabs(a[i][k]);
                index = i;
            }
        }
        // Перестановка строк
        if (max < eps)
        {
            // нет ненулевых диагональных элементов
            cout << "Решение получить невозможно из-за нулевого столбца ";
            cout << index << " матрицы A" << endl;
            return 0;
        }
        for (int j = 0; j < n; j++)
        {
            double temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;
        }
        double temp = y[k];
        y[k] = y[index];
        y[index] = temp;
        // Нормализация уравнений
        for (int i = k; i < n; i++)
        {
            double temp = a[i][k];
            if (fabs(temp) < eps) continue; // для нулевого коэффициента пропустить
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] / temp;
            y[i] = y[i] / temp;
            if (i == k)  continue; // уравнение не вычитать само из себя
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] - a[k][j];
            y[i] = y[i] - y[k];
        }
    }
    // обратная подстановка
    for (k = n - 1; k >= 0; k--)
    {
        x[k] = y[k];
        for (int i = 0; i < k; i++)
            y[i] = y[i] - a[i][k] * x[k];
    }
    return x;
}

//Многочлены Лежандра
double Legendre(int n, double x)
{
    if (n == 0)
        return 1;
    else if (n == 1)
        return x;
    else
        return (Legendre(n - 1, x) * x * (2 * n - 1) / (n)) - (Legendre(n - 2, x) * (n - 1) / (n));
}

//Локализация корней многочленов Лежандра
double* loc_root(double a, double b, double N, int n)
{
    double* mas = new double[100];
    mas[0] = n;
    double h = (b - a) / N;
    double tmp = 0;
    int j = 2;
    while (a + h <= b)
    {
        if ((Legendre(n, a) * Legendre(n, a + h)) <= 0)
        {
            mas[j] = a;
            mas[j + 1] = a + h;
            j += 2;
        }
        a = a + h;
    }
    if (mas[0] > 0)
    {
        //cout << fixed << setprecision(0) << "количество корней ортоганального многочлена на отрезке (a,b) = " << mas[0] << endl;
        return mas;
    }
    else
    {
        //cout << "нет корней на данном отрезке" << endl;
    }
    return 0;
}

//Коэффициенты для формулы Гаусса
double* koeff(double* z, int n)
{
    double* A = new double[100];
    for (int i = 0; i < n; i++)
    {
        A[i] = (2 * (1 - pow(z[i], 2)) / (pow(n, 2) * Legendre(n - 1, z[i]) * Legendre(n - 1, z[i])));
    }
    return A;
}

//Поиск корней многочленов Лежандра
double Secant(double a, double b, double EPS, int n)
{
    double fapprox, sapprox, x;
    int counter = 0;

    fapprox = a;
    sapprox = b;

    x = sapprox - (Legendre(n, sapprox) / (Legendre(n, sapprox) - Legendre(n, fapprox)) * (sapprox - fapprox));
    counter++;

    while (abs(x - sapprox) >= EPS)
    {
        fapprox = sapprox;
        sapprox = x;
        x = sapprox - (Legendre(n, sapprox) / (Legendre(n, sapprox) - Legendre(n, fapprox)) * (sapprox - fapprox));
        counter++;
    }

    return x;
}

//строит массив значений системы многочленов Якоби в данной точке (узле для вычисление интеграла по Гауссу)
double* Jacob(int k, int n, double x)
{
    double* P_x = new double[100];
    int a, b, c, d = 0;
    P_x[0] = 1;
    P_x[1] = (1 + k) * x;
    for (int i = 0; i < n - 2; i++)
    {
        a = (i + k + 2);
        b = (2 * i + 2 * k + 3);
        c = (i + k + 1);
        d = (i + 2 * k + 2) * (i + 2);
        //cout<<a<<" ---- "<<d<<" ---- " << b << " --- "<<c<<endl;
        P_x[i + 2] = (a * b * x * P_x[i + 1] - a * c * P_x[i]) / d;
    }
    /*
    cout << endl;
    cout << "Массив зачений системы многочленов Якоби n = " << n << " в точке = " << x << endl;
    for (int i = 0; i <= n; i++)
    {
        cout <<i<<"   "<< P_x[i] << " "<<endl;
    }
    */
    return P_x;
}
//строит массив значений координатоной системы многочленов Якоби в данной точке (узле для вычисление интеграла по Гауссу)
double* Jacob_coord_sys(int k, int n, double x)
{
    double* phi = Jacob(k, n, x);
    for (int i = 0; i < n; i++)
    {
        phi[i] *= (1 - x * x);
    }
    /*
    cout << endl;
    cout << "Массив зачений координатой системы многочленов Якоби n = " << n << " в точке = " << x << endl;
    for (int i = 0; i <= n; i++)
    {
        cout << i << "   " << phi[i] << " " << endl;
    }
    cout << endl;
    */
    return phi;
}

//строит массив значений производных координатоной системы многочленов Якоби в данной точке (узле для вычисление интеграла по Гауссу)
double* Jacob_Dir_coord_sys(int k, int n, double x)
{
    double* P_x = Jacob(k - 1, n + 1, x);
    double* phi_dir = new double[2 * n];
    phi_dir[0] = -2 * x;
    phi_dir[1] = 2 - 6 * x * x;
    for (int i = 2; i < n; i++)
    {
        phi_dir[i] = -2 * (i + 1) * P_x[i + 1];
    }
    /*
    cout << endl;
    cout << "Массив зачений производных координатой системы многочленов Якоби n = " << n << " в точке = " << x << endl;
    for (int i = 0; i < n; i++)
    {
        cout << i << "   " << phi_dir[i] << " " << endl;
    }
    */
    return phi_dir;
}

//строит массив значений вторых производных координатоной системы многочленов Якоби в данной точке (узле для вычисление интеграла по Гауссу)
double* Jacob_Sec_Dir_coord_sys(int k, int n, double x)
{
    double* P_x = Jacob(k, n, x);
    double* phi_sec_dir = new double[100];
    phi_sec_dir[0] = -2;
    phi_sec_dir[1] = -12 * x;
    for (int i = 2; i < n; i++)
    {
        phi_sec_dir[i] = -1 * (i + 1) * (i + 2 * k) * P_x[i];
    }
    /*
    cout << endl;
    cout << "Массив зачений вторых производных координатой системы многочленов Якоби n = " << n << " в точке = " << x << endl;
    for (int i = 0; i < n; i++)
    {
        cout << i << "   " << phi_sec_dir[i] << " " << endl;
    }
    */
    return phi_sec_dir;
}
/*
получает матрицу со значениями координатной системы n узлах для формулы Гаусса
на ее основе строит значение подыитегрального выражения в данном узле для кажого коэффиента
для вычисления коэффиентов для матрицы системы a[i,j] (интегралов от (Lw_j,w_i) ) и правой части системы
находит решение системы модифицированным методом Гаусса*/
double* Galerkin(int a, int b, int n)
{
    double** phi = get_null_matrix(n);
    double** phi_dir = get_null_matrix(n);
    double** phi_sec_dir = get_null_matrix(n);
    double* z = new double[2 * n];           //массив для узлов для формулы Гаусса
    int val = n;                             // количество узлов для формулы Гаусса
    double* MASS = loc_root(-1, 1, 100, val);//локализуем узлы многочлена Лежандра
    double counter = MASS[0];
    int cnt = 0;
    for (int j = 2; j <= 2 * counter + 1; j += 2)
    {
        z[cnt] = Secant(MASS[j], MASS[j + 1], 0.0000000001, val);//находим узлы (корни многочлена Лежандра)
        cnt++;
    }
    if (val % 2 > 0)
    {
        int m = val / 2;
        z[m] = 0.0;
    }
    double* A = koeff(z, val); //находим коэффициенты формулы гаусса
// cout << "Узлы и коэффициенты для формулы Гаусса: " << endl;
// cout << "      Узел      " << " <-> " << "    Коэффициент    " << endl;
    for (int k = 0; k < cnt; k++)
    {
        z[k] = (b - a) / 2 * z[k] + (b + a) / 2; //узлы для формулы гаусса
        A[k] = (b - a) / 2 * A[k];
// cout << fixed << setprecision(2) << z[k] << "            <->            " << A[k] << endl;
    }
//получим массивы значений координатной системы и ее производной в узлах для формулы гаусса и запишем их в матрицы, чтобы потом из них собирать значения подыитегральной функции.
    double intg = 0.0;
//cout << "Найдем массивы значений функций из координатой системы и их производных в узлах формулы Гаусса: " << endl;
    for (int k = 0; k < cnt; k++)
    {
        phi[k] = Jacob_coord_sys(1, n, z[k]);
        phi_dir[k] = Jacob_Dir_coord_sys(1, n, z[k]);
        phi_sec_dir[k] = Jacob_Sec_Dir_coord_sys(1, n, z[k]);
    }
/*
    cout << "Соберем их в матрицы, где i-я строка - значения производной функций на i-ом узле формулы Гаусса:" << endl;
    show(phi, n);
    cout << "------------------------------------------------------------------------------------------------" << endl;
    show(phi_dir, n);
    cout << "------------------------------------------------------------------------------------------------" << endl;
    show(phi_sec_dir, n);
    cout << "------------------------------------------------------------------------------------------------" << endl;

    cout << "Используя эти матрицы найдем коэффициенты исходной системы как интегралы, используя узлы и коэффициенты формулы Гаусса:" << endl;
*/
    double** AA = get_null_matrix(n);
    double* d = new double[2 * n];
    for (int i = 0; i < n; i++)
    {
        d[i] = 0.0;
        for (int j = 0; j < n; j++)
        {
            AA[i][j] = 0.0;
            for (int k = 0; k < cnt; k++)
                AA[i][j] += A[k] * (P(z[k]) * phi_sec_dir[k][j] * phi[k][i] + Q(z[k]) * phi_dir[k][j] * phi[k][i] + R(z[k]) * phi[k][j] * phi[k][i]);
        }
        for (int k = 0; k < cnt; k++)
        {
            d[i] += A[k] * F(z[k]) * phi[k][i];
        }
    }
    double* res = gauss2(AA, d, n);
/*
    cout << "------------------------------------------------------------------------------------------------" << endl;
    show(AA, n);
    cout << "------------------------------------------------------------------------------------------------" << endl;
    cout << "Правая часть системы:" << endl;
    show_vect(d, n);
    cout << "Система выглядит так:" << endl;
    sysout(AA, d, n);

    cout << "------------------------------------------------------------------------------------------------" << endl;
    cout << "Решение системы (искомые коэффициенты):" << endl;
    show_vect(res, n);
*/
    return res;
}

//Находит искомые коээфициенты как решение системы. Компоненты матрицы это значение L(phi) в узлах коллокации
double* Kolloc(int a, int b, int n)
{
    double** phi = get_null_matrix(n);
    double** phi_dir = get_null_matrix(n);
    double** phi_sec_dir = get_null_matrix(n);
    double* t = new double[2 * n];//узлы коллокации
    double* d = new double[2 * n];//вектор правых частей
    double c = 0.0;
    double bb = 0.0;
    double f = 0.0;
    for (int k = 1; k <= n; k++)
    {
        c = (2 * k - 1);
        bb = 2 * n;
        f = c / bb;
        t[k - 1] = cos(f * M_PI);
        //cout << k - 1 <<" --- "<<f * M_PI<< " --- " << t[k - 1] << " --- узк" << endl;
    }
//cout << "Найдем массивы значений функций из координатой системы и их производных в узлах Коллокации: " << endl;
    for (int k = 0; k < n; k++)
    {
        phi[k] = Jacob_coord_sys(1, n, t[k]);
        phi_dir[k] = Jacob_Dir_coord_sys(1, n, t[k]);
        phi_sec_dir[k] = Jacob_Sec_Dir_coord_sys(1, n, t[k]);
    }
/*
    cout << endl;
    cout << "Соберем их в матрицы, где i-я строка - значения производной функций на i-ом узле формулы Гаусса:" << endl;
    cout << "------------------------------------------------------------------------------------------------" << endl;
    show(phi, n);
    cout << "------------------------------------------------------------------------------------------------" << endl;
    show(phi_dir, n);
    cout << "------------------------------------------------------------------------------------------------" << endl;
    show(phi_sec_dir, n);
    cout << "------------------------------------------------------------------------------------------------" << endl;
    cout << "Используя эти матрицы найдем коэффициенты исходной системы:" << endl;
*/
    double** AA = get_null_matrix(n);
    double* res = new double[2 * n];
    for (int i = 0; i < n; i++)
    {
        res[i] = 0.0;
        d[i] = F(t[i]); // находим правую часть системы как f(t_k)
        for (int j = 0; j < n; j++)
        {
            AA[i][j] = P(t[i]) * phi_sec_dir[i][j] + Q(t[i]) * phi_dir[i][j] + R(t[i]) * phi[i][j];//находим коэффиценты системы как значения L(phi(t_k))
        }
    }
    res = gauss2(AA, d, n);//решаем систему
/*
    cout << "Матрица системы" << endl;
    show(AA, n);
    cout << "Решение системы (искомые коэффициенты):" << endl;
    show_vect(test, n);
*/
    return res;
}

//Собирает по коэффициентам и значениям функций координатной системы значение искомой функции
double func_value(double* koeffs, int n, double x)
{
    double value = 0.0;
    double* phi = Jacob_coord_sys(1, n, x); //массив значений координатной системы в точке x
    for (int k = 0; k < n; k++)
    {
        value += koeffs[k] * phi[k];
    }
    return value;
}
int main()
{
    system("chcp 1251");
    system("cls");
    cout << "Задача 1" << endl;
    cout << "Вариант 3" << endl;
    cout << endl;
    cout << "Уравнение вида: (-1)/(x - 3)u'' + (1 + x/2)u' - exp(x/2)u = 2 - x" << endl;
    cout << endl;
    cout << "Граничные условия: u(-1)=u(1)=0" << endl;
    cout << endl;
    int n = 0;
    double a = -1.0;
    double b = 1.0;
    double NN = 20.0;
    cout << "Введите порядок координатой системы: -> ";
    cin >> n;
    cout << "Введите частоту сетки значений решения: -> ";
    cin >> NN;
    double h = (b - a) / NN;
    cout << "Таблица значений решения. Метод Галеркина." << endl;
    double* test11 = Galerkin(-1, 1, n);
    double* test12 = Galerkin(-1, 1, n + 1);
    double* test13 = Galerkin(-1, 1, n + 2);
    cout << "|---------------|---------------|----------------|----------------|\n";
    cout << "|      x_i      |      n = "<<n<<"    |      n = " << n+1 << "     |      n = " << n +2 << "     | \n";
    cout << "|---------------|---------------|----------------|----------------|\n";
    for (double i = a; i <= b; i += h)
    {
        if (i <= 0)
            cout << fixed << setprecision(6) << "|   " << i << "   |   " << func_value(test11, n, i) << "   |    " << func_value(test12, n + 1, i) << "   |   " << func_value(test13, n + 2, i) << "    |" << endl;
        else
            cout << fixed << setprecision(6) << "|    " << i << "   |   " << func_value(test11, n, i) << "   |    " << func_value(test12, n + 1, i) << "   |   " << func_value(test13, n + 2, i) << "    |" << endl;
    }
    cout << "|---------------|---------------|----------------|----------------|\n";
    double* test2 = Kolloc(-1, 1, n);
    double* test3 = Kolloc(-1, 1, n + 1);
    double* test4 = Kolloc(-1, 1, n + 2);
    cout << "Таблица значений решения. Метод Коллокации." << endl;
    cout << "|---------------|---------------|----------------|----------------|\n";
    cout << "|      x_i      |      n = " << n << "    |      n = " << n + 1 << "     |      n = " << n + 2 << "     | \n";
    cout << "|---------------|---------------|----------------|----------------|\n";
    for (double i = a; i <= b; i += h)
    {
        if (i <= 0)
            cout << fixed << setprecision(6) << "|   " << i << "   |   " << func_value(test2, n, i) << "   |    " << func_value(test3, n + 1, i) << "   |   " << func_value(test4, n + 2, i) << "    |" << endl;
        else
            cout << fixed << setprecision(6) << "|    " << i << "   |   " << func_value(test2, n, i) << "   |    " << func_value(test3, n + 1, i) << "   |   " << func_value(test4, n + 2, i) << "    |" << endl;
    }
    cout << "|---------------|---------------|----------------|----------------|\n";
    cout << "Таблица разностей значений двух методов." << endl;
    cout << "|---------------|---------------|----------------|----------------|\n";
    cout << "|      x_i      |      n = " << n << "    |      n = " << n + 1 << "     |      n = " << n + 2 << "     | \n";
    cout << "|---------------|---------------|----------------|----------------|\n";
    for (double i = a; i <= b; i += h)
    {
        if (i <= 0)
            cout << fixed << setprecision(6) << "|   " << i << "   |   " << fabs(func_value(test11, n, i) -func_value(test2, n, i)) << "    |    " << fabs(func_value(test12, n + 1, i) - func_value(test3, n + 1, i)) << "    |   " << fabs(func_value(test13, n + 2, i) - func_value(test4, n + 2, i)) << "     |" << endl;
        else
            cout << fixed << setprecision(6) << "|    " << i << "   |   " << fabs(func_value(test11, n, i) - func_value(test2, n, i)) << "    |    " << fabs(func_value(test12, n + 1, i) - func_value(test3, n + 1, i)) << "    |   " << fabs(func_value(test13, n + 2, i) - func_value(test4, n + 2, i)) << "     |" << endl;
    }
    cout << "|---------------|---------------|----------------|----------------|\n";
    return 0;
}
