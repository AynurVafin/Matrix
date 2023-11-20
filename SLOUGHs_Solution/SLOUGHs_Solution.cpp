#include <iostream>
#include<cmath>
#include<vector>

using namespace std;

void mult_matrix(vector<vector<double> > A, vector<double> B, vector<double> &result) {
    for (int i = 0; i < A.size(); i++) {
        double ans = 0;
        for (int k = 0; k < B.size(); k++) {
            ans += A[i][k] * B[k];
        }
        result[i] = ans;
    }
}

void inverse_matrix(vector<vector<double> > &A, int n) {
    vector<vector<double> > B(n, vector<double> (n,0));
    for (int i = 0; i < n; ++i) {
        B[i][i] = 1;
    }
    for (int i = 0; i < n; ++i) {
        double tr =  A[i][i];
        for (int j = 0; j < n; ++j) {
            A[i][j] /= tr;
            B[i][j] /= tr;
        }
        for (int k = 0; k < n; ++k) {
            if(k == i) continue;
            double x = A[k][i] / A[i][i];
            for (int j = 0; j < n; ++j) {
                A[k][j] -= x*A[i][j];
                B[k][j] -= x*B[i][j];
            }
        }
    }
    A = B;
}


void determinant(vector<vector<double> > A, double *result) {
    if (A.size() == 1) {
        *result = A[0][0];
    } else if (A.size() == 2) {
        *result = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    } else {
        *result = 0;
        for (int i = 0; i < A.size(); i++) {
            vector<vector<double> > B(A.size() - 1, vector<double>(A.size() - 1, 0));
            for (int r = 1; r < A.size(); r++) {
                int count = 0;
                for (int c = 0; c < A.size(); c++) {
                    if (c != i) {
                        B[r - 1][c - count] = A[r][c];
                    } else {
                        count = 1;
                    }
                }
            }
            double res;
            determinant(B, &res);
            *result += (pow(-1, i % 2) * res * A[0][i]);
        }
    }
}



void clear_terminal() {
    cout << "\033[2J\033[1;1H";
}



void Gaus_matr(vector<vector<double> > &Work, int n, int m) {
    for (int i = 0; i < n - 1; ++i) {
        for (int k = i+1; k < n; ++k) {
            double x = Work[k][i] / Work[i][i];
            for (int j = 0; j < m; ++j) {
                Work[k][j] -= x*Work[i][j];
            }
        }
    }
}

bool check_prog(vector<vector<double> > Matrix, int n) {
    bool flag = true;
    for (int i = 0; i < n; ++i) {
        if(i == 0) {
            int kol = 0;
            for (int j = i+2; j < n; ++j) {
                if(Matrix[i][j]){
                    kol= 1;
                    break;
                }
            }
            if(Matrix[i][i]==0 || kol) {
                flag = false;
                break;
            }
        } else if (i < n - 1) {
            int kol = 0;
            for (int j = 0; j < n; ++j) {
                if(Matrix[i][j])kol++;
            }
            kol =kol - (Matrix[i][i-1] != 0) - (Matrix[i][i+1] !=0);
            if(Matrix[i][i] == 0 || kol!=1){
                flag = false;
                break;
            }
        } else {
            int kol =0;
            for (int j = 0; j < n -2; ++j) {
                if(Matrix[i][j]) {
                    kol = 1;
                    break;
                }
            }
            if(Matrix[i][i] == 0 || kol) {
                flag = false;
                break;
            }
        }
    }
    return flag;
}

void met_prog(vector<vector<double> > Matrix_cpy,int n,int m) {
    vector<double> x(n, 0);
    vector<vector<double> > Matrix(n, vector<double> (m - 1, 0));
    vector<double> B(n, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m - 1; ++j) {
            Matrix[i][j] = Matrix_cpy[i][j];
        }
        B[i] = Matrix_cpy[i][m-1];
    }
    if(check_prog(Matrix, n)) {
        vector<double> a(n,0), b(n, 0);
        double y;
        for (int i = 0; i < n; ++i) {
            if(i == 0) {
                y = Matrix[i][i];
                a[i] = -Matrix[i][i+1]/y;
                b[i] =B[i]/y;
            }else if(i == n-1) {
                y = Matrix[i][i] + Matrix[i][i-1]*a[i-1];
                b[i] = (B[i] - Matrix[i][i-1] *b[i-1])/y;
            } else {
                y = Matrix[i][i] + Matrix[i][i-1]*a[i-1];
                a[i] = -Matrix[i][i+1]/y;
                b[i] = (B[i] - Matrix[i][i-1] *b[i-1])/y;
            }
        }
        for (int i = n-1; i >=0 ; i--) {
            if(i == n-1) x[i] = b[i];
            else {
                x[i] = a[i]*x[i+1]+b[i];
            }
        }
        clear_terminal();
        cout << "Решение Вашего СЛАУ по методу прогонки:\n\n";
        for (int i = 0; i < n; ++i) {
            cout << x[i] << " ";
        }
        cout << "\n\n\n Нажми на клавишу \"Enter\" для продолжения...";
        getchar();
    } else {
        clear_terminal();
        cout << "Этот метод является модификацией метода Гаусса для частного случая"
                " разреженных систем – системы уравнений с трехдиагональной матрицей.\n\n";
        cout << "В данном случае это не трехдиагональная матрица, воспользуйтесь Методом Гауса";
        cout << "\n\nДля продолжения нажмити на Enter";
        getchar();
    }
}

void met_kram(vector<vector<double> > Matrix_cpy,int n,int m) {
    vector<double> x(n, 0);
    vector<vector<double> > Matrix(n, vector<double> (m - 1, 0));
    vector<double> B(n, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m - 1; ++j) {
            Matrix[i][j] = Matrix_cpy[i][j];
        }
        B[i] = Matrix_cpy[i][m-1];
    }

    double res = 0;
    determinant(Matrix, &res);
    for (int i = 0; i < n; ++i) {
        vector<vector<double> > Cpy(n, vector<double> (m - 1, 0));
        Cpy = Matrix;
        double res_cpy = 0;
        for (int j = 0; j < n; ++j) {
            Cpy[j][i] = B[j];
        }

        determinant(Cpy, &res_cpy);

        x[i] = res_cpy/res;
    }

    clear_terminal();
    cout << "Решение Вашего СЛАУ по методу Крамера:\n\n";
    for (int i = 0; i < n; ++i) {
        cout << x[i] << " ";
    }
    cout << "\n\n\n Нажми на клавишу \"Enter\" для продолжения...";
    getchar();
}

void met_matr(vector<vector<double> > Matrix_cpy,int n,int m) {
    vector<vector<double> > Matrix(n, vector<double> (m - 1));
    vector<double> B(n, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m - 1; ++j) {
            Matrix[i][j] = Matrix_cpy[i][j];
        }
        B[i] = Matrix_cpy[i][m-1];
    }
    vector<double>x(n, 0);
    inverse_matrix(Matrix, n);
    mult_matrix(Matrix, B, x);
    clear_terminal();
    cout << "Решение Вашего СЛАУ по методу обратной матрицы:\n\n";
    for (int i = 0; i < n; ++i) {
        cout << x[i] << " ";
    }
    cout << "\n\n\n Нажми на клавишу \"Enter\" для продолжения...";
    getchar();
}

void met_gaus(vector<vector<double> > Matrix_cpy,int n,int m) {
    Gaus_matr(Matrix_cpy,n , m);
    vector<double> x(n, 0);
    for (int i = n-1; i >=0 ; i--) {
        double ans = Matrix_cpy[i][m - 1];
        for (int j = n-1; j > i ; j--) {
            ans -= Matrix_cpy[i][j] * x[j];
        }
        x[i] = ans/Matrix_cpy[i][i];
    }
    clear_terminal();
    cout << "Решение Вашего СЛАУ по методу Гауса:\n\n";
    for (int i = 0; i < n; ++i) {
        cout << x[i] << " ";
    }
    cout << "\n\n\n Нажми на клавишу \"Enter\" для продолжения...";
    getchar();
}




int check_ans(vector<vector<double> > Matrix_cpy,int n,int m) {
    int flag = 0;
    int rA = n, rA1 = n;
    Gaus_matr(Matrix_cpy,n , m);
    for (int i = n-1; i >= 0; i--) {
        if(!Matrix_cpy[i][m-1]) {
            rA1--;
        } else {break;}
    }

    for (int i = n-1; i >= 0; i--) {
        int prov = 1;
        for (int j = 0; j < m - 1; ++j) {
            if(Matrix_cpy[i][j]) {
                prov = 0;
                break;
            }
        }
        if (prov) rA--;
        else break;
    }

    if(rA == rA1 && rA  < m - 1){
        flag = 2;
    } else if(rA == rA1 && rA == m - 1){
        flag = 1;
    }
    return flag;
}

int main() {
    char a;
    cout << "Здравствуйте! Это программа для решения СЛАУ\n\n\n Нажми на клавишу \"Enter\" для продолжения...";
    getchar();
    clear_terminal();
    int n, m;
    cout << "Введите n, m, где n - количество уравнений в системе, m - количество коэффициентов при неизвестных\n";
    int kol_flag;
    while ((kol_flag = scanf("%d %d%c", &n, &m, &a)) != 3 || a != '\n' || n <= 0 || m <= 0) {
        fflush(stdin);
        clear_terminal();
        cout << "Вы ввели неверно или некорректно. Повторите ещё раз\n\n";
        cout << "Введите n, m, где n - количество уравнений в системе, m - количество коэффициентов при неизвестных\n";
    }
    vector<vector<double> > Matrix_SLAY(n, vector<double>(m + 1)), Matrix_cpy(n, vector<double>(m + 1));
    cout << "Введите в виде матрицы коэффициенты при неизвестных и свободных коэффициентов.\n\n";
    cout << "(Сначала идут m коэффициентов при неизвестных уравнения, а потом свободный коэффициент)\n\n";
    vector<vector<double> > Matrix(n, vector<double>(m));
    m++;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double x;
            int kol;
            while ((kol = scanf("%lf%c", &x, &a)) != 2 || (a != '\n' && a != ' ')) {
                fflush(stdin);
                cout << "Вы ввели неверно или некорректно. Строчку: " << i + 1 << " Столбец: " << j + 1;
                cout << "\n\nПовторите ещё раз\n\n";
            }
            Matrix_SLAY[i][j] = x;
            Matrix_cpy[i][j] = x;
            if (j != m - 1) {
                Matrix[i][j] = x;
            }
        }
    }
    //Проверка на нет решения
    int flag = check_ans(Matrix_cpy, n, m);
    double res;
    determinant(Matrix, &res);
    if (!res || n < m -1) {
        if (flag == 2) {
            cout << "Данное СЛАУ имеет бесконечно число решений. \n\n"
                    "Для решения еще одного СЛАУ, перезапустите программу!\n\n";
        } else {
            cout << "Данное СЛАУ не имеет решений. \n\nДля решения еще одного СЛАУ, перезапустите программу!\n\n";
        }
    } else {
        while (true) {
            bool flag_exit = 0;
            clear_terminal();
            cout << "Выберите один из методов решения СЛАУ и подтвердите клавишей \"Enter\":\n\n";
            cout << "Метод прогонки - нажми на 1 и подтверди \"Enter\"" << endl;
            cout << "Метод Крамера - нажми на 2 и подтверди \"Enter\"" << endl;
            cout << "Метод обратной матрицы - нажми на 3 и подтверди \"Enter\"" << endl;
            cout << "Метод Гауса - нажми на 4 и подтверди \"Enter\"\n\n" << endl;
            cout << "Выход - нажми на 0 и подтверди \"Enter\"\n\n" << endl;
            fflush(stdin);
            a = getchar();
            char b;
            b = getchar();
            fflush(stdin);
            while ((b != ' ' && b != '\n') || a < '0' || a > '4') {


                cout << "Вы ввели неверно или некорректно. Повторите ещё раз\n\n";
                a = getchar();
                if (a == '\n') continue;
                b = getchar();
                fflush(stdin);
            }
            switch (a) {
                case '0':
                    flag_exit = 1;
                    break;
                case '1':
                    met_prog(Matrix_cpy, n, m);
                    break;// some method
                case '2':
                    met_kram(Matrix_cpy, n, m);
                    break;
                case '3':
                    met_matr(Matrix_cpy, n, m);
                    break;
                case '4':
                    met_gaus(Matrix_SLAY, n, m);
                    break;
                default:
                    break;
            }
            if (flag_exit) break;
        }
        clear_terminal();
    }
    cout << "Пока!";
}