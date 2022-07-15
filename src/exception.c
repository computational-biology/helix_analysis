//
// Created by parthajit on 15/8/20.
//

#include "exception.h"

void assert_else(int assert_expr, char* else_msg){
    if(assert_expr == 0) return;
    fprintf(stderr,"Error... %s\n",else_msg);
    exit(EXIT_FAILURE);
    return;
}
