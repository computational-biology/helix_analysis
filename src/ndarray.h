//
// Created by parthajit on 15/8/20.
//

#ifndef CONSHELIX_NDARRAY_H
#define CONSHELIX_NDARRAY_H

#include <stdio.h>
#include <stdlib.h>


void matrixc_fprintf(FILE* fp, char** M, long row, long col);

char** matrixc_create(long row, long col);

void matrixc_setval(char** M, int row, int col, char val);
void matrixc_free(char** M, long row, long col);




void matrixl_fprintf(FILE* fp, long** M, long row, long col);

long** matrixl_create(long row, long col);

void matrixl_setval(long** M, long row, long col, long val);




void matrixl_free(long** M, long row, long col);







void matrixi_fprintf(FILE* fp, int** M, int row, int col);

int** matrixi_create(int row, int col);

void matrixi_setval(int** M, int row, int col, int val);

void matrixi_free(int** M, int row, int col);

#endif //CONSHELIX_NDARRAY_H
