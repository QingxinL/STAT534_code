//
//  main.cpp
//  Problem2
//
//  Created by Yuxuan Cheng on 5/11/17.
//  Copyright © 2017 Yuxuan Cheng. All rights reserved.
//

#include <iostream>

typedef struct
{
    size_t size1;
    size_t size2;
    size_t tda;
    double * data;
    gsl_block * block;
    int owner;
} gsl_matrix;

double marglik(gsl_matrix* data,int lenA,int* A)
{
    int n = data->size1;
    int p = data->size2;
    
    //creates a submatrix of matrix M
    int* IndRow = new int[n];
    for (int i=0; i<n; i++)
    {
        IndRow[i] = i;
    }
    int* firstColumn = new int[1];
    firstColumn[0] = 0;
    gsl_matrix* D1 = MakeSubmatrix(data, IndRow, n, firstColumn, 1);
    gsl_matrix* D_A = MakeSubmatrix(data, IndRow, n, A, lenA);
    
    gsl_matrix* D_At = transposematrix(D_A);
    
    gsl_matrix* M_A = gsl_matrix_alloc(lenA, lenA);
    matrixproduct(D_A, D_At, M_A);
    
    gsl_matrix* diagA = diagMatrix(lenA);
    gsl_matrix_add(M_A, diagA);
    printmatrix("M_A.txt",M_A);
    
    gsl_matrix* temp2 = gsl_matrix_alloc(1, 1);
    matrixproduct(transposematrix(D1), D1, temp2);
    
    last_term = 1 + gsl_matrix_get(temp2,0,0) - ;
    
    log_marglik = 
    
    gsl_matrix_free(D1);
    gsl_matrix_free(D_A);
    gsl_matrix_free(D_At);
    gsl_matrix_free(M_A);
    gsl_matrix_free(diagA);
    gsl_matrix_free(temp2);
    gsl_matrix_free();
    gsl_matrix_free();
    
    return ();
    
}


#include "matrices.h"
int main()
{
    int n = 158; //sample size
    int p = 51; //number of variables
    int i;
    int A[] = {2,5,10};//indices of the variables present in the regression
    int lenA = 3; //number of indices
    char datafilename[] = "erdata.txt";
    //allocate the data matrix
    gsl_matrix* data = gsl_matrix_alloc(n,p);
    //read the data
    FILE* datafile = fopen(datafilename,"r");
    if(NULL==datafile)
    {
        fprintf(stderr,"Cannot open data file [%s]\n",datafilename);
        return(0); }
    if(0!=gsl_matrix_fscanf(datafile,data))
    {
        fprintf(stderr,"File [%s] does not have the required format.\n",datafilename);
        return(0); }
    fclose(datafile);
    printf("Marginal likelihood of regression [1|%d",A[0]);
    for(i=1;i<lenA;i++)
    {
        printf(",%d",A[i]);
    }
    printf("] = %.5lf\n",marglik(data,lenA,A));
    //free memory
    gsl_matrix_free(data);
    return(1);
}

//Your task is to write the function

//If everything goes well, your program should run like this:
//stu5:~/534> ./matrices
//Marginal likelihood of regression [1|2,5,10] = -59.97893
//Points will be deducted if your version of the function uses numerical routines not defined in the library “GSL”. It is okay to use usual mathematical functions defined in “math.h”, e.g. “log” or “lgamma”.
