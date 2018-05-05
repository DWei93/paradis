/***************************************************************************
 *
 *      Module:       Math.C
 *      Description:  Contains various math functions
 *
 *      Includes:   
 *                  Arranege (int& n);
 *                  Arranege (int& n, int& m);
 *
 **************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>

#include "Math.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

void FatalMath(const char *format, ...) 
{
        char    msg[512];
        va_list args;

        va_start(args, format);
        vsnprintf(msg, sizeof(msg)-1, format, args);
        msg[sizeof(msg)-1] = 0;
        va_end(args);
        printf("Fatal - Math: %s\n", msg);

#ifdef PARALLEL
        MPI_Abort(MPI_COMM_WORLD, -1);
#endif
        exit(1);
}

int Arrange(int n){
    if(n <= 0)return(1);

    int a = 1, i;
    for(i=n; i>0; i--){
        a *= i;
    }
    return(a);
}

int Arrange(int n, int m){
    if(n <= 0 || m <= 0)return(1);

    int nSubm = n - m;
    return(Arrange(n)/Arrange(nSubm));
}

int Combine(int n, int m){
    if(n <= 0 || m <= 0)return(1);
    
    return(Arrange(n,m)/Arrange(m));
}

#if 0
void Perm(int deep, int n, int m, bool *isUsed, int **results){
    int i;
    if (deep > n-m)return;

    for(i=0; i<n; i++){
        if(!isUsed[i]){
            isUsed[i] = 1;
            results[deep] = i;
            Perm(deep+1, n, m, isUsed, results);
            isUsed[i] = 0;
        }
    }
}
#endif

void GetArrange(int n, int m, int **arra){
    bool        *isUsed;
    double      *results;

    isUsed = (bool *)calloc(n,sizeof(bool));
    results = (double *)calloc(Arrange(n,m), sizeof(double));
}

void GetCombination(int n, int m, int **comb){

    bool *isUsed;
    isUsed = (bool *)calloc(n,sizeof(bool));
    
    
}

double Volume(double **p, int num){
    int         i, numComb;
    double      volume;
    
    if(num < 3)FatalMath("Volume: the shape has %d edges",num);
    numComb = Combine(num, 3);

//    for(i=0; i<; i++){
        
    
//    }

    return(0);   
}
