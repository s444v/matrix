#ifndef S21_MATRIX_H
#define S21_MATRIX_H

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

#define SUCCESS 1
#define FAILURE 0
#define EPSILON 1e-8

int s21_error(matrix_t *A);

int s21_create_matrix(int rows, int columns,
                      matrix_t *result);  // выделение памяти

void s21_remove_matrix(matrix_t *A);  // очистка

int s21_eq_matrix(matrix_t *A, matrix_t *B);  // сравнение

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);

int s21_mult_number(matrix_t *A, double number, matrix_t *result);

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);

int s21_transpose(matrix_t *A, matrix_t *result);

int s21_inverse_matrix(matrix_t *A, matrix_t *result);

int s21_determinant(matrix_t *A, double *result);

int s21_calc_complements(matrix_t *A, matrix_t *result);

void scan_matrix(matrix_t *A, matrix_t *B, int row, int col);
double det(matrix_t *A);

#endif