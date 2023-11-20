#include"s21_matrix.h"

int s21_calc_complements(matrix_t *A, matrix_t *result) {
    int flag = 0;
    if (A->rows <= 0 || A->columns <= 0) flag = 1;
    else if(A->rows != A->columns) flag = 2;
    else {
        s21_remove_matrix(result);
        s21_create_matrix(A->rows, A->columns, result);
        for (int i = 0; i < A->rows; i++)
        {
            for (int j = 0; j < A->columns; j++)
            {
                matrix_t B;
                s21_create_matrix(A->rows - 1, A->columns - 1, &B);
                int count_r = 0;
                for (int r = 0; r < A->rows; r++)
                {
                    if(r != i) {
                        int count_c = 0;
                    for (int c = 0; c < A->columns; c++)
                    {
                     if(c != j) {
                        B.matrix[r-count_r][c-count_c] = A->matrix[i][j];
                     }   else {
                        count_c = 1;
                     } 
                    }
                    } else {
                    count_r = 1;
                    }
                }
                double res;
                flag =  s21_determinant(&B, &res);
                result->matrix[i][j] = pow(-1, (i+j)%2) * res;
                s21_remove_matrix(&B);
            }
        }
        
    }
    return flag;
}