#include"s21_matrix.h"

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int flag = 0;
     if(A->columns <= 0 || A->rows <= 0 || B->columns<=0 || B->rows <=0)
    flag = 1;
    else if (A->rows != B->rows || A->columns != B->columns) {
        flag = 2;
    } else {
        s21_remove_matrix(result);
        flag = s21_create_matrix(A->rows,A->columns, result);
        if (!flag)
        for (int i = 0; i < A->rows; i++)
        {
            for (int j = 0; j < A->columns; j++)
            {
                result->matrix[i][j] = A->matrix[i][j]+B->matrix[i][j];
            }
            
        }
    }
    return flag;
}