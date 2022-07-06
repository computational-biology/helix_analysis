//
// Created by parthajit on 15/8/20.
//

#ifndef CONSHELIX_CHAR_EDITDIST_H
#define CONSHELIX_CHAR_EDITDIST_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "exception.h"
#include "ndarray.h"





void char_Needleman_Wunsch_seq_align(char* s, char* t, char** s_aligned, char** t_aligned);




void char_smith_waterman_seq_align(char* s, char* t, char** s_aligned, char** t_aligned);



















#endif //CONSHELIX_CHAR_EDITDIST_H
