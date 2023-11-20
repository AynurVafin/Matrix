#include"s21_matrix.h"

int s21_determinant(matrix_t *A, double *result) {
    int flag = 0;
    if(A->columns <= 0 || A->rows <= 0)
    flag = 1;
    else if (A->columns != A->rows) {
        flag = 2;
    } else {
        if (A->columns == 1) {
            *result = A->matrix[0][0];
        } else  if (A->columns == 2){
            *result = A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
        } else {
            matrix_t B;
            *result = 0;
            for (int i = 0; i < A->columns; i++)
            {
               flag = s21_create_matrix(A->rows - 1, A->columns -1, &B);
                for (int r = 1; r < A->rows; r++)
                {
                    int count = 0;
                    for (int c = 0; c < A->columns; c++)
                    {
                        if(c != i) {
                            B.matrix[r-1][c - count] = A->matrix[r][c];
                        } else {
                            count = 1;
                        }
                    }
                }
                double res;
                flag = s21_determinant(&B, &res);
                *result += (pow(-1, i%2) * res * A->matrix[0][i]);
               s21_remove_matrix(&B);
            }
        }
    }
    return flag;
}