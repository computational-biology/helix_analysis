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

#include "editdist.h"


#define MIN3(a, b, c)  ((a) < (b) ? ((a) < (c)? (a):(c)): ((b) < (c) ? (b) : (c)))
#define MAX3(a, b, c)  ((a) > (b) ? ((a) > (c)? (a):(c)): ((b) > (c) ? (b) : (c)))


int nw_simget(int a, int b){
      return (a == b)? 20 : -10;
}
void NW_sim(int NW[][128]){
    for(int i=0; i<127; ++i){
        for(int j=0; j<127; ++j){
            if(i==j){
                NW[i][j] = 2;
            }else{
                NW[i][j] = -1;
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


void NW_kernel(const int* s, int sn, const int* t, int tn, int* s_aligned, int* sz_sa, int* t_aligned, int* sz_ta){

    int penalty = -2;

    //int ** F = (int**)malloc((sn+1) * sizeof(int*));


    int** F = matrixi_create(sn+1, tn+1);


    /*init the 1st row and 1st col */
    F[0][0] = 0;
    for(int i=1; i<=sn; ++i){
        F[i][0] = F[i-1][0] + penalty;
    }
    for(int j=1; j<=tn; ++j){
        F[0][j] = F[0][j-1] + penalty;
    }

    /*populate F matrix*/
    for(int i=1; i<=sn; ++i){
        for(int j=1; j<=tn; ++j){
//            int match = F[i-1][j-1] + NWSim[s[i-1]][t[j-1]];
            int match = F[i-1][j-1] + nw_simget(s[i-1], t[j-1]);
            int del = F[i-1][j] + penalty;
            int ins = F[i][j-1] + penalty;
            F[i][j] = MAX3(match, del, ins);
        }
    }
    //matrixl_fprintf(stdout, F, sn+1, tn+1);
    int* A = (int*) malloc( (sn + tn + 1) * sizeof(int));
    int* B = (int*) malloc( (sn + tn + 1) * sizeof(int));
    int indx=0;
    int i = sn;
    int j = tn;

    while( i > 0 || j > 0){
//        if(i > 0 && j > 0 && F[i][j] == F[i-1][j-1]+NWSim[s[i-1]][t[j-1]]){
        if(i > 0 && j > 0 && F[i][j] == F[i-1][j-1]+nw_simget(s[i-1], t[j-1])){
            A[indx] = s[i-1];
            B[indx] = t[j-1];
            i--;
            j--;
            indx++;
        }else if( i > 0 && F[i][j] == F[i-1][j]+penalty){
            A[indx] = s[i-1];
            B[indx] = -1; // -1 means blank '-'
            i--;
            indx++;
        }else{
            A[indx] = -1;
            B[indx] = t[j-1];
            j--;
            indx++;
        }
    }
    //int* a_align = (int*) malloc((indx+1) * sizeof(int));
    //int* b_align = (int*) malloc((indx+1) * sizeof(int));
//    s_aligned[indx] = '\0';
//    t_aligned[indx] = '\0';
    *sz_sa = indx;
    *sz_ta = indx;
    for(int i=0; i<indx; ++i){
        s_aligned[indx-1-i] = A[i];
        t_aligned[indx-1-i] = B[i];
    }
    free(A);
    free(B);
    matrixi_free(F, sn+1, tn+1);

    return;
}




void needleman_wunsch_seq_align(int* s, int sz_s, int* t, int sz_t, int** s_aligned, int* sz_sa, int** t_aligned, int* sz_ta){
//    int NWSim[127][128];
//    NW_sim(NWSim);


    int sn = sz_s;
    int tn = sz_t; 
    int size;

    int* s_align = (int*) malloc((sn+tn+1) * sizeof(int));
    int* t_align = (int*) malloc((sn+tn+1) * sizeof(int));

    NW_kernel(s, sn, t, tn,  s_align, sz_sa, t_align, sz_ta);
    *s_aligned = s_align;
    *t_aligned = t_align;
    return;
}

//void needleman_wunsch_revscore(int* A, int na, int* B, int nb, int NWSim[][128], int* NW_score, int* scoreTmp){
//    int* row1 = NW_score;//(int*) malloc( (nb+1) * sizeof(int));
//    int* row2 = scoreTmp; //(int*) malloc( (nb+1) * sizeof(int));
//    printf("nb=%ld\n",nb);
//    int penalty = -2;
//    row1[nb] = 0;
//    for(int j=nb-1; j>=0; --j){
//        printf("%c ",B[j]);
//        row1[j] = row1[j+1] + penalty;
//    }
//    printf("\n");
//    for(int i=na; i>0; --i){
//        row2[nb] = row1[nb] + penalty;
//        for(int j=nb; j>0; --j){
//            //printf("%ld,%ld\n",i, j);
//            int inscost = row2[j+1] + penalty;
//            int delcost = row1[j]   + penalty;
//            int subcost = row1[j+1] + NWSim[A[i-1]][B[j-1]];
//            row2[j] = MAX3(inscost, delcost, subcost);
//            printf("%ld,%ld,%ld, %c, %c,\n",i, j, row2[j], A[i-1],B[j-1]);
//        }
//        for(int j=0; j<=nb; ++j){
//            printf("%ld,.,  ", row1[j]);
//            row1[j] = row2[j];
//        }
//        printf("\n");
//    }
//    for(int j=0; j<=nb; ++j){
//        printf("%ld,, ", row1[j]);
//    }
//    printf("\n");
//    //free(row2);
//    return;
//}
//void Needleman_Wunsch_score(int* A, int na, int* B, int nb, int NWSim[][128], int* NW_score, int* scoreTmp){
//    int* row1 = NW_score;//(int*) malloc( (nb+1) * sizeof(int));
//    int* row2 = scoreTmp;//(int*) malloc( (nb+1) * sizeof(int));
//    int penalty = -2;
//    /*printf("values are\n");
//      for(int a=0;a<na; ++a){
//      printf("%c",A[a]);
//      }
//      printf("\n");
//      for(int a=0;a<nb; ++a){
//      printf("%c",B[a]);
//      }
//      printf("\n");*/
//    row1[0] = 0;
//    for(int j=1; j<=nb; ++j){
//        row1[j] = row1[j-1] + penalty;
//    }
//    for(int i = 1; i<= na; ++i){
//        row2[0] = row1[0] + penalty;
//        for(int j=1; j<=nb; ++j){
//            int inscost = row2[j-1] + penalty;
//            int delcost = row1[j]   + penalty;
//            int subcost = row1[j-1] + NWSim[A[i-1]][B[j-1]];
//            row2[j] = MAX3(inscost, delcost, subcost);
//        }
//        for(int j=0; j<=nb; ++j){
//            printf("%ld,, ", row1[j]);
//            row1[j] = row2[j];
//        }
//        printf("\n");
//    }
//
//    for(int j=0; j<=nb; ++j){
//        printf("%ld,, ", row1[j]);
//    }
//    printf("\n---\n");
//    return;
//}
//void sw_sim_populate(int NW[][128]){
//           NW['A']['A'] = NW['A']['A'] = 10;
//           NW['A']['C'] = NW['C']['A'] = -5;
//           NW['A']['U'] = NW['U']['A'] = -6;
//           NW['A']['G'] = NW['G']['A'] = -1;
//
//           NW['C']['G'] = NW['G']['C'] = -5;
//           NW['C']['C'] = NW['C']['C'] =  8;
//           NW['C']['U'] = NW['U']['C'] = -1;
//
//           NW['G']['G'] = NW['G']['G'] = 10;
//           NW['G']['U'] = NW['U']['G'] = -3;
//
//           NW['U']['U'] = NW['U']['U'] =  5;
//}
//int sw_sim_score(int NW[][128],char a, char b){
//	if(a == b) return 3;
//	return -3;
//    return NW[a][b];
//}
//
//int gap_penalty(int gap, int extension, int opening){
//	return extension * gap + opening ;
//}
//
//
//
//
//void smith_waterman_seq_align(char* s, char* t, char** s_aligned, char** t_aligned){
//    int SW[128][128];
//    sw_sim_populate(SW);
//    int sn = strlen(s);
//    int tn = strlen(t);
//    int** H = matrixl_create(sn+1, tn+1); // Scoring matrix
//    for(int i=0; i< sn+1; ++i){
//        H[i][0] = 0;
//    }
//    for(int j=0; j<tn+1; ++j){
//        H[0][j] = 0;
//    }
//    int diag, left, top;
//    int score;
//    int maxval = 0;
//    int imax = 0;
//    int jmax = 0;
//    for(int i=1; i<sn+1; ++i){
//        for(int j=1; j<tn+1; ++j){
//            diag = H[i-1][j-1] + sw_sim_score(SW, s[i-1], t[j-1]);
//            top = 0;
//	    for(int k=1; k<=i; ++k){
//		    score = H[i-k][j] - gap_penalty(k, 2, 0);
//		    if(score > top){
//			    top = score;
//		    }
//	    }
//	    left = 0;
//	    for(int l=1; l<=j; ++l){
//		    score = H[i][j-l] - gap_penalty(l, 2, 0);
//		    if(score > left){
//			    left = score;
//		    }
//	    }
//	    H[i][j] = MAX3(diag, top, left); // No need to check with zero. that is adjusted in top and left
//	    if(H[i][j] > maxval){
//		    maxval = H[i][j];
//		    imax = i;
//		    jmax = j;
//	    }
//	}
//    }
//    assert_else(maxval == 0,"Not a single character machct found in Smith-Waterman");
//    char* tmp_s_align = (char*) malloc((sn+tn+1) * sizeof(char)); 
//    char* tmp_t_align = (char*) malloc((sn+tn+1) * sizeof(char)); 
//    int count = 0;
//    int i = imax + 1;
//    int j = jmax + 1;
//    while(1){
//	    double largest = MAX3(H[i-1][j-1], H[i][j-1], H[i-1][j]);
//	    if(largest == 0){
//		    tmp_s_align[count] = '\0';
//		    tmp_t_align[count] = '\0';
//		    break;
//	    }
//	    
//	    if(largest == H[i-1][j-1]){
//		    i = i-1;
//		    j = j-1;
//		    tmp_s_align[count] = s[i-1];
//		    tmp_t_align[count] = t[j-1];
//
//	    }else if(largest == H[i][j-1]){
//		    j = j-1;
//		    tmp_s_align[count] = '-';
//		    tmp_t_align[count] = t[j-1];
//	    }else{
//		    i = i-1;
//		    tmp_s_align[count] = s[i-1];
//		    tmp_t_align[count] = '-';
//	    }
//
//	    count ++;
//    }
//    printf("\n%s\n%s\n", tmp_s_align,tmp_t_align);
//
//
//
//
//
//
//
//    for(int i=0; i< sn+1; ++i){
//	    for(int j=0; j<tn+1; ++j){
//		    printf("%4ld", H[i][j]);
//	    }
//	    printf("\n");
//    }
//
//}
//
//
//void _HB_compute(char* s, int sbeg, int send, char* t, int tbeg, int tend, char* revs, char* revt, int NWSim[][128], char* s_align, char* t_align, int* scoreL, int* scoreR, int* scoreTmp){
//    int sn = (send+1 - sbeg);
//    int tn = (tend+1 - tbeg);
//    char* A = s_align;
//    char* B = t_align;
//    int indx=0;
//    if(sn == 0){
//        printf("Hi:%ld, %ld\n",send, sbeg);
//        for(int j=tbeg+0; j<tn; ++j){
//            A[j] = '-';
//            B[j] = t[j];
//        }
//        A[tn] = '\0';
//        B[tn] = '\0';
//    }else if(tn == 0){
//        printf("Ho:%ld, %ld\n",tend, tbeg);
//        for(int i=sbeg+0; i<sn; ++i){
//            A[i] = s[i];
//            B[i] = '-';
//        }
//        A[sn] = '\0';
//        B[sn] = '\0';
//    }else if(sn == 1 || tn == 1){
//        NW_kernel(s+sbeg, sn, t+tbeg, tn, NWSim, s_align, t_align);
//    }else{
//        //int smid = (sbeg + send) / 2;
//        int smid = sn/2;
//        smid = sbeg + smid -1;
//        int len = smid - sbeg + 1;
//        Needleman_Wunsch_score(s+sbeg, len, t+tbeg, tn, NWSim, scoreL, scoreTmp);
//        len = send - smid;
//
//        Needleman_Wunsch_score(revs+sn-send-1, len, revt+tn-tend-1, tn, NWSim, scoreR, scoreTmp);
//        int max=scoreL[0] + scoreR[tn];
//        int tmid=-1;
//        int currval;
//        for(int i=1; i<=tn; ++i){
//            currval = scoreL[i] + scoreR[tn-i];
//            if(currval > max){
//                max = currval;
//                tmid = i-1;
//            }
//        }
//
//        tmid = tbeg + tmid;
//
//        printf("Best s:");
//        for(int i=sbeg; i<=smid; ++i){
//            printf("%c",s[i]);
//        }
//        printf(", ");
//        for(int i=smid+1; i<=send; ++i){
//            printf("%c",s[i]);
//        }
//        printf("\n");
//
//        printf("Best t:");
//        for(int i=tbeg; i<=tmid; ++i){
//            printf("%c",t[i]);
//        }
//        printf(", ");
//        for(int i=tmid+1; i<=tend; ++i){
//            printf("%c",t[i]);
//        }
//        printf("\n");
//
//        char* ls_align = (char*) malloc(((smid-sbeg+1)+(tmid-tbeg+1) +10) * sizeof(char));
//        char* lt_align = (char*) malloc(((smid-sbeg+1)+(tmid-tbeg+1) +10) * sizeof(char));
//
//        _HB_compute(s, sbeg, smid, t, tbeg, tmid, revs, revt, NWSim, ls_align, lt_align, scoreL, scoreR, scoreTmp);
//
//        char* rs_align = (char*) malloc(((send - smid)+(tend - tmid) +10) * sizeof(char));
//        char* rt_align = (char*) malloc(((send - smid)+(tend - tmid) +10) * sizeof(char));
//
//        _HB_compute(s, smid+1, send, t, tmid+1, tend, revs, revt, NWSim, rs_align, rt_align, scoreL, scoreR, scoreTmp);
//
//        strcpy(s_align,ls_align);
//        strcat(s_align,rs_align);
//
//
//        strcpy(t_align,lt_align);
//        strcat(t_align,rt_align);
//        free(ls_align);
//        free(lt_align);
//        free(rs_align);
//        free(rt_align);
//
//    }
//    return;
//}
//
//
//
//
//void hirschberg_seq_align(char* s, char* t, char** s_aligned, char** t_aligned){
//    int NWSim[128][128];
//    NW_sim(NWSim);
//    int sn = strlen(s);
//    int tn = strlen(t);
//
//    char* revs = (char*) malloc( (sn + 1) * sizeof(char));
//    char* revt = (char*) malloc( (tn + 1) * sizeof(char));
//
//    for(int i=0; i<sn; ++i){
//        revs[sn-1 - i] = s[i];
//    }
//    revs[sn] = '\0';
//
//    for(int i=0; i<tn; ++i){
//        revt[tn-1 - i] = t[i];
//    }
//    revt[tn] = '\0';
//
//    char* s1 = (char*) malloc((sn+tn+1) * sizeof(char));
//    char* t1 = (char*) malloc((sn+tn+1) * sizeof(char));
//
//    int* scoreL = (int*) malloc((tn+10) * sizeof(int));
//    int* scoreR = (int*) malloc((tn+10) * sizeof(int));
//    int* scoreTmp = (int*) malloc((tn+10) * sizeof(int));
//    printf("w1: %ld, %ld\n",sn, tn);
//    _HB_compute(s, 0, sn-1, t, 0, tn-1, revs, revt, NWSim, s1, t1, scoreL, scoreR, scoreTmp);
//
//    char* s_align = (char*) malloc((strlen(s1)+1) * sizeof(char));
//    char* t_align = (char*) malloc((strlen(t1)+1) * sizeof(char));
//    strcpy(s_align, s1);
//    strcpy(t_align, t1);
//    *s_aligned = s_align;
//    *t_aligned = t_align;
//
//    printf("here6\n");
//    free(scoreL);
//    printf("here7\n");
//    free(scoreR);
//    printf("here8\n");
//    free(scoreTmp);
//    printf("here9\n");
//
//    free(s1);
//    free(t1);
//
//    free(revs);
//    free(revt);
//    //printf("Try:%s\nTry:%s\n",s1,t1);
//}
//
//
//
//
//
//
//
//
//
//
//int WFins(char a){
//    switch(a){
//        case 'A':
//        case 'C':
//        case 'G':
//        case 'U':
//        case 'T': return 1;
//        default:
//            return 1;
//    }
//}
//
//int WFdel(char a){
//    switch(a){
//        case 'A':
//        case 'C':
//        case 'G':
//        case 'U':
//        case 'T': return 1;
//        default:
//            return 1;
//    }
//}
//
//int WFsub(char a, char b){
//    if(a == 'A' && b == 'G'){
//        return 1;
//    }else{
//        return 1;
//    }
//}
//
//
//int edit_dist(const char* source, const char* target, char** s_aligned, char** t_aligned){
//    int sn = strlen(source);
//    int tn = strlen(target);
//
//    int ** DP = matrixl_create(sn+1, tn+1);
//    //char ** trace = matrixc_create(sn+1, tn+1);
//    for(int i=0; i<=sn; ++i){
//        DP[i][0] = i;
//    }
//    for(int j=0; j<=tn; ++j){
//        DP[0][j] = j;
//    }
//    for(int i=1; i<=sn; ++i){
//        for(int j=1; j<=tn; ++j){
//
//            if(source[i-1] == target[j-1]){
//                DP[i][j] = DP[i-1][j-1];
//            }else{
//                int min = MIN3(DP[i-1][j]+1, DP[i][j-1]+1,DP[i-1][j-1]+1);
//                DP[i][j] = min;
//            }
//        }
//    }
//    matrixl_fprintf(stdout, DP,sn+1, tn+1);
//
//    char* A = (char*) malloc((sn+tn+1)* sizeof(char));
//    char* B = (char*) malloc((sn+tn+1)* sizeof(char));
//    int indx = 0;
//    int i = sn;
//    int j = tn;
//    while(i>0 || j>0){
//        if( i>0 && j >0 && (source[i-1] == target[j-1] || DP[i][j] == DP[i-1][j-1]+1)){
//            A[indx] = source[i-1];
//            B[indx] = target[j-1];
//            i--;
//            j--;
//            indx++;
//        }else if(i>0 && DP[i][j] == DP[i-1][j] + 1){
//            A[indx] = source[i-1];
//            B[indx] = '-';
//            i--;
//            indx++;
//        }else{
//            A[indx] = '-';
//            B[indx] = target[j-1];
//            j--;
//            indx++;
//        }
//    }
//
//    char* a_align = (char*) malloc((indx+1) * sizeof(char));
//    char* b_align = (char*) malloc((indx+1) * sizeof(char));
//
//    a_align[indx] = '\0';
//    b_align[indx] = '\0';
//    for(int i=0; i<indx; ++i){
//        a_align[indx-1-i] = A[i];
//        b_align[indx-1-i] = B[i];
//    }
//    *s_aligned = a_align;
//    *t_aligned = b_align;
//    int dist= DP[sn][tn];
//
//    free(A);
//    free(B);
//    matrixl_free(DP, sn+1, tn+1);
//    return dist;
//}
//
//
//
//int Wagner_Fischer_weighted_edit_dist(const char* source, const char* target, char* action_seq, char* action_char){
//    //waterman book.
//    //dolittle
//    //1n32, 1aw4
//    int sn = strlen(source);
//    int tn = strlen(target);
//
//    int ** DP = matrixl_create(sn+1, tn+1);
//    char ** trace = matrixc_create(sn+1, tn+1);
//    for(int i=0; i<=sn; ++i){
//        DP[i][0] = i*WFdel(' ');
//        trace[i][0] = '|';
//    }
//    for(int j=0; j<=tn; ++j){
//        DP[0][j] = j*WFins(' ');
//        trace[0][j] = '-';
//    }
//    trace[0][0] = '.';
//    int  min;
//    int  mi, mj;
//    int subcost = 0;
//    for(int i=1; i<=sn; ++i){
//        for(int j=1; j<=tn; ++j){
//            if(source[i-1] == target[j-1]){
//                DP[i][j] = DP[i-1][j-1];
//                trace[i][j] = '\\';
//            }else{
//                int CD = WFdel(source[i-1]);
//                int CI = WFins(target[j-1]);
//                int CS = WFsub(source[i-1], target[j-1]);
//                min = MIN3(DP[i-1][j]+CD, DP[i][j-1]+CI,DP[i-1][j-1]+CS);
//                DP[i][j] = min;
//                if(min == DP[i-1][j]+CD){
//                    trace[i][j] = '|';
//                }else if(min == DP[i][j-1]+CI){
//                    trace[i][j] = '-';
//                }else{
//                    trace[i][j] = '\\';
//                }
//            }
//        }
//    }
//    matrixl_fprintf(stdout, DP,sn+1, tn+1);
//    matrixc_fprintf(stdout, trace, sn+1, tn+1);
//    if((action_seq == NULL && action_char != NULL) || (action_seq != NULL && action_char == NULL)){
//        fprintf(stderr, "Error... one of action_seq and action_char is not null and the other is null\n");
//        exit(3);
//    }
//    if(action_seq == NULL && action_char == NULL){
//        int dst= DP[sn][tn];
//        matrixl_free(DP, sn+1, tn+1);
//        matrixc_free(trace, sn+1, tn+1);
//        return dst;
//    }
//
//    int maxlen = sn>tn? sn:tn;
//    char* editact = (char*) malloc((sn+tn+1)* sizeof(char));//actions
//    char* editseq = (char*) malloc((sn+tn+1)* sizeof(char));
//    int indx = 0;
//    int i = sn;
//    int j = tn;
//    while(i>0 || j>0){
//        if(trace[i][j] == '-'){
//            editact[indx] = 'I';
//            editseq[indx] = target[j-1];
//            indx++;
//            j--;
//        }else if(trace[i][j] == '|'){
//            editact[indx] = 'D';
//            editseq[indx] = source[i-1];
//            indx++;
//            i--;
//        }else{
//            if(source[i-1] == target[j-1]){
//                editact[indx] = '.';
//                editseq[indx] = target[j-1];
//                indx++;
//            }else{
//                editact[indx] = 'S';
//                editseq[indx] = target[j-1];
//                indx++;
//            }
//            i--;
//            j--;
//        }
//    }
//
//    action_seq  = (char*) malloc((indx+1) * sizeof(char));
//    action_char = (char*) malloc((indx+1) * sizeof(char));
//    action_seq[indx] = '\0';
//    action_char[indx] = '\0';
//    for(int i=0; i<indx; ++i){
//        action_seq[indx-1-i] = editact[i];
//        action_char[indx-1-i] = editseq[i];
//    }
//
//    while(i>0){
//        printf("here 1... Program terminates\n");
//        exit(1);
//        editact[indx] = 'D';
//        editseq[indx] = source[i-1];
//        indx--;
//        i--;
//    }
//    while(j>0){
//        printf("here 2... Program terminates\n");
//        exit(2);
//        editact[indx] = 'I';
//        editseq[indx] = target[j-1];
//        indx--;
//        j--;
//    }
//    int dist= DP[sn][tn];
//    free(editseq);
//    free(editact);
//    matrixl_free(DP, sn+1, tn+1);
//    matrixc_free(trace, sn+1, tn+1);
//    printf("%s\n%s\n%s\n%s\n",source, action_seq, action_char, target);
//    return dist;
//}
//
//
//void opt_seq_align_Hirschberg(char* str1, int len1, char* str2, int len2){
//    //char* z = (char*) malloc( (maxlen+1) * sizeof(char));
//    //char* w = (char*) malloc( (maxlen+1) * sizeof(char));
//    ;
//}
//
//
//
///* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// *
// *  KMP_match_first(...)
// *
// *  This algorithm matches the first occurance of the pattern in the text.
// *  The return value is 0 or +ve value to indicate the match start index in the text array.
// *  The return value -1 indicates no match found.
// *
// *  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
//
//int KMP_match_first(char* pattern, int p_len, char* text, int t_len){
//
//    /* computation of lps(intest proper suffix) array starts here */
//    int* lps = (int*) malloc(p_len * sizeof(int));
//    int k = 0;
//    lps[0] = 0;
//    for(int i=1; i<p_len; ++i){
//        while(k>0 && pattern[k] != pattern[i]){
//            k = lps[k-1];
//        }
//        if(pattern[k] == pattern[i]){
//            k = k + 1;
//        }
//        lps[i] = k;
//    }
//
//    /* computation of lps(intest proper suffix) array ends here */
//    for(int i=0; i< p_len; ++i){
//        printf("%ld ",lps[i]);
//    }
//
//    int q=0;
//    int indx = -1;
//    for(int i=0; i<t_len; ++i){
//        while(q>0 && pattern[q] != text[i]){
//            q = lps[q-1];
//        }
//        if(pattern[q] == text[i]){
//            q = q + 1;
//        }
//        if(q == p_len){
//            indx = i-p_len+1;  // start of the match
//            break;
//        }
//    }
//    return indx;
//
//
//
//}
//
//
///* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// *
// *    lc_subseq(...)
// *
// *    This function computes the intest common subsequence of two strings. The Implementation  *
// *    is done based on Coremen et al. book on Algorithm. The Algorithm uses dynamic programming *
// *    and DP table for its execution.
// *
// *    INPUT: A is a string of size na. (array size is na+1 for accomodationof \0)
// *    INPUT: B is a string of size nb. (array size is nb+1 for accomodationof \0)
// *    OUTPUT: lcs_out is a pointer to char for storing the LCS. (size is min(na, nb) + 1
// *    OUTPUT: n is a pointer to int for storing the LCS size.
// *    OUTPUT: A1 is the modified A array with all non LCS positions as '.'. Some thing like
// *    if LCS is BAB, the A is like ..B.A..B.. (size is na+1)
// *    OUTPUT: B1 is same as A1 but for B.     (size is nb+1)
// *
// *    RETURN: It returns the length of LCS.
// *
// *    Note: In all the cases of lcs_out, A1 and B1 the memory SHOULD NOT be allocated by the user.
// *    The algorithm DOSE allocate all the memory for them. The user needs to call
// *    lc_subseq_free(...) when the task is done.
// *
// *    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
//
//void lc_subseq_free(char* lcs, char* a, char* b){
//    if(lcs != NULL)
//        free(lcs);
//    if(a != NULL)
//        free(a);
//    if(b != NULL)
//        free(b);
//}
//
//int lc_subseq(const char* A, const char* B, char** lcs, char** a, char** b){
//    const int na = strlen(A);
//    const int nb = strlen(B);
//    /* memory allocation starts*/
//    int** V = matrixl_create(na+1, nb+1);
//    char** U;
//    U = (char**) malloc((na+1) * sizeof(char*));
//    for(int i=0; i<na+1; ++i){
//        U[i] = (char*) malloc((nb+1) * sizeof(char));
//    }
//
//    /* momory allocation ends */
//
//
//    /* dynamic programming part starts */
//    for(int i=1; i<=na; ++i){
//        V[i][0] = 0;
//    }
//    for(int i=0; i<=nb; ++i){
//        V[0][i] = 0;
//    }
//    for(int i=1; i<=na; ++i){
//        for(int j=1; j<=nb; ++j){
//            if(A[i-1] == B[j-1]){
//                V[i][j] = V[i-1][j-1] + 1;
//                U[i][j] = '\\';
//            }else if(V[i-1][j] >= V[i][j-1]){
//                V[i][j] = V[i-1][j];
//                U[i][j] = '|';
//            }else{
//                V[i][j] = V[i][j-1];
//                U[i][j] = '-';
//            }
//        }
//    }
//    /* dynamic programming part ends */
//
//    /* output generationstarts */
//    int minsize = (na < nb) ? na : nb;
//    char* tmparr = (char*) malloc(minsize * sizeof(char));
//    int index = 0;
//    int i = na;
//    int j = nb;
//
//
//    char* A1 = (char*) malloc((na+1) * sizeof(char));
//    A1[i] = '\0';
//    char* B1 = (char*) malloc((nb+1) * sizeof(char));
//    B1[j] = '\0';
//    while( i != 0 && j != 0){
//        if( U[i][j] == '|'){
//            i = i - 1;
//            if(A1 != NULL){
//                A1[i] = '.';
//            }
//
//        }else if(U[i][j] == '-'){
//            j = j - 1;
//            if(B1 != NULL){
//                B1[j] = '.';
//            }
//        }else{
//            tmparr[index] = A[i-1]; // A[i-1] is the ith char
//            index++;
//            i = i - 1;
//            j = j - 1;
//            if(A1 != NULL){
//                A1[i] = A[i];
//            }
//            if(B1 != NULL){
//                B1[j] = B[j];
//            }
//        }
//    }
//    for(int k = i-1; k>=0; --k){
//        A1[k] = '.';
//    }
//    *a = A1;
//    for(int k = j-1; k>=0; --k){
//        B1[k] = '.';
//    }
//    *b = B1;
//
//
//    char* lcs_out = (char*) malloc((index+1) * sizeof(char));
//    for(int i=0; i<index; ++i){
//        lcs_out[index-1 - i] = tmparr[i];
//    }
//    lcs_out[index] = '\0';
//    *lcs = lcs_out;
//
//    /* output generation ends */
//
//
//    /* memory free part */
//    free(tmparr);
//
//    matrixl_free(V, na, nb);
//    for(int i=0; i< na+1; ++i){
//        free(U[i]);
//    }
//    free(U);
//    return index;
//}
