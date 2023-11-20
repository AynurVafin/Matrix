#include"s21_matrix.h"

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
    int flag = 0;
    if (A->columns <= 0 || A->rows <= 0) flag = 1;
    else if (A->columns != A->rows) flag = 2;
    else {
        s21_remove_matrix(result);
        s21_create_matrix(A->rows, A->columns, result);
        double res_der;
        flag = s21_determinant(A, &res_der);
        if(res_der) {
            matrix_t res_prom, res;
            flag = s21_calc_complements(A, &res_prom);
            flag = s21_transpose(&res_prom, &res);
            flag = s21_mult_number(&res, res_der, result);
            s21_remove_matrix(&res_prom);
            s21_remove_matrix(&res);
        } else {
            flag = 2;
        }
    }
    return flag;
}