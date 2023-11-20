#include"s21_matrix.h"

void s21_remove_matrix(matrix_t *A) {
    if (A != NULL ){//
    int n = A->rows;
    for (int i = 0; i < n; i++)
    {
        free(A->matrix[i]);
    }
    //free(A->matrix);
    free(A);
    }
}