#include "s21_matrix.h"

#include <math.h>
#include <stdlib.h>

int s21_error(matrix_t *A) {
  int flag = 0;
  if (A == NULL || A->matrix == NULL) flag = 1;
  return flag;
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int flag = 1;
  if (rows > 0 && columns > 0) {
    flag = 0;
    result->rows = rows;
    result->columns = columns;
    result->matrix =
        malloc(rows * columns * sizeof(double) + rows * sizeof(double *));
    double *ptr = (double *)(result->matrix + rows);
    for (int i = 0; i < rows; i++) {
      result->matrix[i] = ptr + columns * i;
    }
    for (int i = 0; i < result->rows && !flag; i++) {
      for (int j = 0; j < result->columns; j++) {
        result->matrix[i][j] = 0;
      }
    }
    flag = s21_error(result);
  }
  return flag;
}

void s21_remove_matrix(matrix_t *A) {
  free(A->matrix);
  A->rows = 0;
  A->columns = 0;
  A->matrix = NULL;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int flag = 1;
  if (A->rows != B->rows || A->columns != B->columns) {
    flag = 0;
  }
  for (int i = 0; i < A->rows && flag; i++) {
    for (int j = 0; j < A->columns && flag; j++) {
      if (fabs(A->matrix[i][j] - B->matrix[i][j]) > EPSILON) {
        flag = 0;
      }
    }
  }
  return flag;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int flag = 0;
  flag = s21_create_matrix(A->rows, A->columns, result);
  if (s21_error(A) || s21_error(B)) {
    flag = 1;
  }
  if (A->rows != B->rows || A->columns != B->columns) {
    flag = 2;
  }
  for (int i = 0; i < A->rows && flag == 0; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    }
  }
  return flag;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int flag = 0;
  flag = s21_create_matrix(A->rows, A->columns, result);
  if (s21_error(A) || s21_error(B)) {
    flag = 1;
  }
  if (A->rows != B->rows || A->columns != B->columns) {
    flag = 2;
  }
  for (int i = 0; i < A->rows && flag == 0; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
    }
  }
  return flag;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int flag = 0;
  flag = s21_create_matrix(A->rows, A->columns, result);
  if (s21_error(A)) {
    flag = 1;
  }
  for (int i = 0; i < A->rows && flag == 0; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] * number;
    }
  }
  return flag;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int flag = 0;
  flag = s21_create_matrix(A->rows, B->columns, result);
  if (s21_error(A) || s21_error(B)) {
    flag = 1;
  }
  if (A->rows != B->columns || A->columns != B->rows) {
    flag = 2;
  }
  for (int i = 0; i < result->rows && flag == 0; i++) {
    for (int j = 0; j < result->columns; j++) {
      for (int k = 0; k < A->columns; k++) {
        result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
      }
    }
  }
  return flag;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int flag = 0;
  if (s21_error(A)) {
    flag = 1;
  } else {
    flag = s21_create_matrix(A->columns, A->rows, result);
    for (int i = 0; i < result->rows && flag == 0; i++) {
      for (int j = 0; j < result->columns; j++) {
        result->matrix[i][j] = A->matrix[j][i];
      }
    }
  }
  return flag;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int flag = 0;
  if (A == NULL || A->matrix == NULL || A->rows <= 0 || A->columns <= 0) {
    flag = 1;
  } else if (A->columns != A->rows) {
    flag = 2;
  } else if (A->columns == 1 && A->rows == 1) {
    s21_create_matrix(A->rows, A->columns, result);
    result->matrix[0][0] = 1;
  } else {
    int pow = 1;
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < result->rows; i++) {
      for (int j = 0; j < result->columns; j++) {
        if ((i + j) % 2 != 0)
          pow = -1;
        else
          pow = 1;
        matrix_t B = {0};
        s21_create_matrix(A->rows - 1, A->columns - 1, &B);
        scan_matrix(A, &B, i, j);
        result->matrix[i][j] = det(&B) * pow;
        s21_remove_matrix(&B);
      }
    }
  }
  return flag;
}

int s21_determinant(matrix_t *A, double *result) {
  int flag = 0;
  if (A == NULL || A->matrix == NULL || A->rows <= 0 || A->columns <= 0) {
    flag = 1;
  } else if (A->columns != A->rows) {
    flag = 2;
  } else {
    *result = det(A);
  }
  return flag;
}
double det(matrix_t *A) {
  double pow = 1;
  double deter = 0;
  if (A->rows == 1 && A->columns == 1) {
    deter = A->matrix[0][0];
  } else if (A->rows == 2 && A->columns == 2) {
    deter =
        A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
  } else {
    matrix_t B = {0};
    s21_create_matrix(A->rows - 1, A->columns - 1, &B);
    for (int i = 0; i < A->columns; i++) {
      scan_matrix(A, &B, 0, i);
      deter += pow * A->matrix[0][i] * det(&B);
      pow *= (-1);
    }
    s21_remove_matrix(&B);
  }
  return deter;
}

void scan_matrix(matrix_t *A, matrix_t *B, int row, int col) {
  int row_shift = 0;
  int col_shift = 0;
  for (int i = 0; i < A->rows - 1; i++) {
    if (i == row) {
      row_shift = 1;
    }
    for (int j = 0; j < A->columns - 1; j++) {
      if (j == col) {
        col_shift = 1;
      }
      B->matrix[i][j] = A->matrix[i + row_shift][j + col_shift];
    }
    col_shift = 0;
  }
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int flag = 0;
  double determinant = 0;
  if (A->columns != A->rows) {
    flag = 2;
  } else if (A->rows <= 0 || A->columns <= 0) {
    flag = 1;
  }
  if (flag == 0) {
    flag = s21_determinant(A, &determinant);
  }
  if (fabs(determinant) < 1e-6 && flag == 0) {
    flag = 2;
  } else if (flag == 0) {
    matrix_t B;
    s21_calc_complements(A, &B);
    s21_transpose(&B, result);
    s21_remove_matrix(&B);
    for (int i = 0; i < result->rows; i++) {
      for (int j = 0; j < result->columns; j++) {
        result->matrix[i][j] /= determinant;
      }
    }
  }
  return flag;
}