/*
 * =====================================================================================
 *
 *       Filename:  editdist.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  Monday 23 May 2022 05:36:54  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#include "char_editdist.h"


#define MIN3(a, b, c)  ((a) < (b) ? ((a) < (c)? (a):(c)): ((b) < (c) ? (b) : (c)))
#define MAX3(a, b, c)  ((a) > (b) ? ((a) > (c)? (a):(c)): ((b) > (c) ? (b) : (c)))

void char_NW_sim(long NW[][128]){
    for(int i=0; i<127; ++i){
        for(int j=0; j<127; ++j){
            if(i==j){
                NW[i][j] = 10;
            }else{
                NW[i][j] = -10;
            }
        }
    }
    /*
       NW['A']['A'] = NW['A']['A'] = 10;
       NW['A']['C'] = NW['C']['A'] = -5;
       NW['A']['U'] = NW['U']['A'] = -6;
       NW['A']['G'] = NW['G']['A'] = -1;

       NW['C']['G'] = NW['G']['C'] = -5;
       NW['C']['C'] = NW['C']['C'] =  8;
       NW['C']['U'] = NW['U']['C'] = -1;

       NW['G']['G'] = NW['G']['G'] = 10;
       NW['G']['U'] = NW['U']['G'] = -3;

       NW['U']['U'] = NW['U']['U'] =  5;
       */
}


void char_NW_kernel(const char* s, long sn, const char* t, long tn, long NWSim[][128], char* s_aligned, char* t_aligned){

    long penalty = -2;

    //long ** F = (long**)malloc((sn+1) * sizeof(long*));


    long** F = matrixl_create(sn+1, tn+1);


    /*init the 1st row and 1st col */
    F[0][0] = 0;
    for(long i=1; i<=sn; ++i){
        F[i][0] = F[i-1][0] + penalty;
    }
    for(long j=1; j<=tn; ++j){
        F[0][j] = F[0][j-1] + penalty;
    }

    /*populate F matrix*/
    for(long i=1; i<=sn; ++i){
        for(long j=1; j<=tn; ++j){
            long match = F[i-1][j-1] + NWSim[s[i-1]][t[j-1]];
            long del = F[i-1][j] + penalty;
            long ins = F[i][j-1] + penalty;
            F[i][j] = MAX3(match, del, ins);
        }
    }
    //matrixl_fprintf(stdout, F, sn+1, tn+1);
    char* A = (char*) malloc( (sn + tn + 1) * sizeof(char));
    char* B = (char*) malloc( (sn + tn + 1) * sizeof(char));
    long indx=0;
    long i = sn;
    long j = tn;

    while( i > 0 || j > 0){
        if(i > 0 && j > 0 && F[i][j] == F[i-1][j-1]+NWSim[s[i-1]][t[j-1]]){
            A[indx] = s[i-1];
            B[indx] = t[j-1];
            i--;
            j--;
            indx++;
        }else if( i > 0 && F[i][j] == F[i-1][j]+penalty){
            A[indx] = s[i-1];
            B[indx] = '-';
            i--;
            indx++;
        }else{
            A[indx] = '-';
            B[indx] = t[j-1];
            j--;
            indx++;
        }
    }
    //char* a_align = (char*) malloc((indx+1) * sizeof(char));
    //char* b_align = (char*) malloc((indx+1) * sizeof(char));
    s_aligned[indx] = '\0';
    t_aligned[indx] = '\0';
    for(long i=0; i<indx; ++i){
        s_aligned[indx-1-i] = A[i];
        t_aligned[indx-1-i] = B[i];
    }
    free(A);
    free(B);
    matrixl_free(F, sn+1, tn+1);

    return;
}




void char_Needleman_Wunsch_seq_align(char* s, char* t, char** s_aligned, char** t_aligned){
    long NWSim[127][128];
    char_NW_sim(NWSim);


    long sn = strlen(s);
    long tn = strlen(t);
    long size;

    char* s_align = (char*) malloc((sn+tn+1) * sizeof(char));
    char* t_align = (char*) malloc((sn+tn+1) * sizeof(char));

    char_NW_kernel(s, sn, t, tn, NWSim, s_align, t_align);
    *s_aligned = s_align;
    *t_aligned = t_align;
    return;
}

void char_Needleman_Wunsch_revscore(char* A, long na, char* B, long nb, long NWSim[][128], long* NW_score, long* scoreTmp){
    long* row1 = NW_score;//(long*) malloc( (nb+1) * sizeof(long));
    long* row2 = scoreTmp; //(long*) malloc( (nb+1) * sizeof(long));
    printf("nb=%ld\n",nb);
    long penalty = -2;
    row1[nb] = 0;
    for(long j=nb-1; j>=0; --j){
        printf("%c ",B[j]);
        row1[j] = row1[j+1] + penalty;
    }
    printf("\n");
    for(long i=na; i>0; --i){
        row2[nb] = row1[nb] + penalty;
        for(long j=nb; j>0; --j){
            //printf("%ld,%ld\n",i, j);
            long inscost = row2[j+1] + penalty;
            long delcost = row1[j]   + penalty;
            long subcost = row1[j+1] + NWSim[A[i-1]][B[j-1]];
            row2[j] = MAX3(inscost, delcost, subcost);
            printf("%ld,%ld,%ld, %c, %c,\n",i, j, row2[j], A[i-1],B[j-1]);
        }
        for(int j=0; j<=nb; ++j){
            printf("%ld,.,  ", row1[j]);
            row1[j] = row2[j];
        }
        printf("\n");
    }
    for(int j=0; j<=nb; ++j){
        printf("%ld,, ", row1[j]);
    }
    printf("\n");
    //free(row2);
    return;
}
void char_Needleman_Wunsch_score(char* A, int na, char* B, int nb, long NWSim[][128], long* NW_score, long* scoreTmp){
    long* row1 = NW_score;//(long*) malloc( (nb+1) * sizeof(long));
    long* row2 = scoreTmp;//(long*) malloc( (nb+1) * sizeof(long));
    long penalty = -2;
    /*printf("values are\n");
      for(long a=0;a<na; ++a){
      printf("%c",A[a]);
      }
      printf("\n");
      for(long a=0;a<nb; ++a){
      printf("%c",B[a]);
      }
      printf("\n");*/
    row1[0] = 0;
    for(long j=1; j<=nb; ++j){
        row1[j] = row1[j-1] + penalty;
    }
    for(long i = 1; i<= na; ++i){
        row2[0] = row1[0] + penalty;
        for(long j=1; j<=nb; ++j){
            long inscost = row2[j-1] + penalty;
            long delcost = row1[j]   + penalty;
            long subcost = row1[j-1] + NWSim[A[i-1]][B[j-1]];
            row2[j] = MAX3(inscost, delcost, subcost);
        }
        for(int j=0; j<=nb; ++j){
            printf("%ld,, ", row1[j]);
            row1[j] = row2[j];
        }
        printf("\n");
    }

    for(int j=0; j<=nb; ++j){
        printf("%ld,, ", row1[j]);
    }
    printf("\n---\n");
    return;
}
void char_sw_sim_populate(long NW[][128]){
           NW['A']['A'] = NW['A']['A'] = 10;
           NW['A']['C'] = NW['C']['A'] = -5;
           NW['A']['U'] = NW['U']['A'] = -6;
           NW['A']['G'] = NW['G']['A'] = -1;

           NW['C']['G'] = NW['G']['C'] = -5;
           NW['C']['C'] = NW['C']['C'] =  8;
           NW['C']['U'] = NW['U']['C'] = -1;

           NW['G']['G'] = NW['G']['G'] = 10;
           NW['G']['U'] = NW['U']['G'] = -3;

           NW['U']['U'] = NW['U']['U'] =  5;
}
long char_sw_sim_score(long NW[][128],char a, char b){
	if(a == b) return 3;
	return -3;
    return NW[a][b];
}

long char_gap_penalty(long gap, long extension, long opening){
	return extension * gap + opening ;
}




void char_smith_waterman_seq_align(char* s, char* t, char** s_aligned, char** t_aligned){
    long SW[128][128];
    char_sw_sim_populate(SW);
    long sn = strlen(s);
    long tn = strlen(t);
    long** H = matrixl_create(sn+1, tn+1); // Scoring matrix
    for(long i=0; i< sn+1; ++i){
        H[i][0] = 0;
    }
    for(long j=0; j<tn+1; ++j){
        H[0][j] = 0;
    }
    long diag, left, top;
    long score;
    long maxval = 0;
    long imax = 0;
    long jmax = 0;
    for(long i=1; i<sn+1; ++i){
        for(long j=1; j<tn+1; ++j){
            diag = H[i-1][j-1] + char_sw_sim_score(SW, s[i-1], t[j-1]);
            top = 0;
	    for(long k=1; k<=i; ++k){
		    score = H[i-k][j] - char_gap_penalty(k, 2, 0);
		    if(score > top){
			    top = score;
		    }
	    }
	    left = 0;
	    for(long l=1; l<=j; ++l){
		    score = H[i][j-l] - char_gap_penalty(l, 2, 0);
		    if(score > left){
			    left = score;
		    }
	    }
	    H[i][j] = MAX3(diag, top, left); // No need to check with zero. that is adjusted in top and left
	    if(H[i][j] > maxval){
		    maxval = H[i][j];
		    imax = i;
		    jmax = j;
	    }
	}
    }
    assert_else(maxval == 0,"Not a single character machct found in Smith-Waterman");
    char* tmp_s_align = (char*) malloc((sn+tn+1) * sizeof(char)); 
    char* tmp_t_align = (char*) malloc((sn+tn+1) * sizeof(char)); 
    long count = 0;
    long i = imax + 1;
    long j = jmax + 1;
    while(1){
	    double largest = MAX3(H[i-1][j-1], H[i][j-1], H[i-1][j]);
	    if(largest == 0){
		    tmp_s_align[count] = '\0';
		    tmp_t_align[count] = '\0';
		    break;
	    }
	    
	    if(largest == H[i-1][j-1]){
		    i = i-1;
		    j = j-1;
		    tmp_s_align[count] = s[i-1];
		    tmp_t_align[count] = t[j-1];

	    }else if(largest == H[i][j-1]){
		    j = j-1;
		    tmp_s_align[count] = '-';
		    tmp_t_align[count] = t[j-1];
	    }else{
		    i = i-1;
		    tmp_s_align[count] = s[i-1];
		    tmp_t_align[count] = '-';
	    }

	    count ++;
    }
    printf("\n%s\n%s\n", tmp_s_align,tmp_t_align);







    for(long i=0; i< sn+1; ++i){
	    for(long j=0; j<tn+1; ++j){
		    printf("%4ld", H[i][j]);
	    }
	    printf("\n");
    }

}


