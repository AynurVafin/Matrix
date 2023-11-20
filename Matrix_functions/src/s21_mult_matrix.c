#include"s21_matrix.h"

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int flag = 0;
    if(A->columns <= 0 || A->rows <= 0 || B->columns<=0 || B->rows <=0)
    flag = 1;
    else if (A->columns != B->rows) {
        flag = 2;
    } else {
        s21_remove_matrix(result);
        flag = s21_create_matrix(A->rows, B->columns, result);
        for (int i = 0; i < A->rows; i++)
        {
            for (int  j= 0; j < B->columns; j++)
            {
                double ans = 0;
               for(int k = 0; k < A->columns; k++) {
                ans+= A->matrix[i][k] * B->matrix[k][j];
               }
               result->matrix[i][j] = ans;
            }
            
        }
    }
    return flag;
}