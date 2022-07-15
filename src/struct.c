//
// Created by parthajit on 13/7/20.
//


#include "struct.h"

void fasta_init(struct fasta* self, char *file, int size) {
    self->num_res = size;
    if(self->num_res <= 0) return;
    self->secseq = (char*) malloc ((self->num_res+1) * sizeof(char));
    self->num_chain = 0;
    FILE* fpdat = fopen(file, "r");
    assert( fpdat != NULL);
    char line[1024];
    char sep[]="\t :\n";
    char* token;
    long int index = 0;
    while(fgets(line, 1024, fpdat) != NULL){
        if(line[0] == '>'){
            token = strtok(line, sep);
            token = strtok(NULL, sep);//second one is the chain info
            strcpy(self->chain[self->num_chain], token);
            self->num_chain++;
        }else{
            token = strtok(line, sep);
            int len = strlen(token);
            assert(len>0);
            strcpy(self->secseq+index, token);
            index = index + len;
        }
    }
    if( self->num_res != index ){    /* Exception Handling */ 
	  fprintf(stderr, "Error in function %s()... count mismatch in sequence reading.\n", __func__);
	  exit(EXIT_FAILURE);
    }

    fclose(fpdat);
}


void fasta_free(struct fasta* self) {
    if(self->num_res > 0){
        free(self->secseq);
        free(self->priseq);
    }
    //printf("Sec free called\n");
}
