#include"s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
    int flag = 0;
    if (rows < 1 || columns < 1) flag = 1;
    else {
        matrix_t cpy_result;
        cpy_result.rows = rows;
        cpy_result.columns = columns;
        cpy_result.matrix = (double**)malloc(sizeof(double) * rows);
        for (int i = 0; i < rows; i++)
        {
           cpy_result.matrix[i] = (double*)malloc(sizeof(double)*columns);
        }
        *result = cpy_result;
    }
    return flag;
}