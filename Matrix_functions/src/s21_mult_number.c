#include"s21_matrix.h"

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
    int flag = 0;
    s21_remove_matrix(result);
    flag = s21_create_matrix(A->rows, A->columns, result);
    *result = *A;
    for (int i = 0; i < result->rows; i++)
    {
        for (int j = 0; j < result->columns; j++)
        {
            result->matrix[i][j]*=number;
        }
        
    }
    return flag;
}