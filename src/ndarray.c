//
// Created by parthajit on 15/8/20.
//


#include "ndarray.h"

void matrixc_fprintf(FILE* fp, char** M, long row, long col){
    for(long i=0; i<row; ++i){
        for(long j=0; j<col; ++j){
            fprintf(fp,"%c ", M[i][j]);
        }
        fprintf(fp, "\n");
    }
    return;
}

char** matrixc_create(long row, long col){
    char** M = (char**) malloc(row * sizeof(char*)); // Energy matrix
    for(long i=0; i<row; ++i){
        M[i] = (char*) malloc(col * sizeof(char));
    }
    return M;
}

void matrixc_setval(char** M, int row, int col, char val){
    for(long i=0; i<row; ++i){
        for(long j=0; j<col; ++j){
            M[i][j] = val;
        }
    }
}
void matrixc_free(char** M, long row, long col){
    for(long i=0; i<row; ++i){
        free(M[i]);
    }
    free(M);
    return;
}

void matrixi_fprintf(FILE* fp, int** M, int row, int col){
    for(int i=0; i<row; ++i){
        for(int j=0; j<col; ++j){
            fprintf(fp,"%2d ", M[i][j]);
        }
        fprintf(fp, "\n");
    }
    return;
}

int** matrixi_create(int row, int col){
    int** M = (int**) malloc(row * sizeof(int*)); // Energy matrix
    for(int i=0; i<row; ++i){
        M[i] = (int*) malloc(col * sizeof(int));
    }
    return M;
}

void matrixi_setval(int** M, int row, int col, int val){
    for(int i=0; i<row; ++i){
        for(int j=0; j<col; ++j){
            M[i][j] = val;
        }
    }
}

void matrixi_free(int** M, int row, int col){
    for(int i=0; i<row; ++i){
        free(M[i]);
    }
    free(M);
    return;
}








void matrixl_fprintf(FILE* fp, long** M, long row, long col){
    for(long i=0; i<row; ++i){
        for(long j=0; j<col; ++j){
            fprintf(fp,"%2ld ", M[i][j]);
        }
        fprintf(fp, "\n");
    }
    return;
}

long** matrixl_create(long row, long col){
    long** M = (long**) malloc(row * sizeof(long*)); // Energy matrix
    for(long i=0; i<row; ++i){
        M[i] = (long*) malloc(col * sizeof(long));
    }
    return M;
}

void matrixl_setval(long** M, long row, long col, long val){
    for(long i=0; i<row; ++i){
        for(long j=0; j<col; ++j){
            M[i][j] = val;
        }
    }
}




void matrixl_free(long** M, long row, long col){
    for(long i=0; i<row; ++i){
        free(M[i]);
    }
    free(M);
    return;
}


