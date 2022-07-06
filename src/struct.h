//
// Created by parthajit on 13/7/20.
//

#ifndef CPPMET_SECSTRUCT_H
#define CPPMET_SECSTRUCT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>



struct fasta{

    char* secseq;
    char* priseq;
    char chain[200][6];
    int num_chain;
    int num_res;
    };

    void fasta_init(struct fasta* self, char* file, int size);
    void fasta_free(struct fasta* self);





#endif //CPPMET_SECSTRUCT_H
