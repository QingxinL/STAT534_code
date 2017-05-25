//
//  main.cpp
//  Problem2
//
//  Created by Yuxuan Cheng on 5/11/17.
//  Copyright © 2017 Yuxuan Cheng. All rights reserved.
//


#include <iostream>
#include "matrices.h"

double marglik(gsl_matrix* data,int lenA, int* A)
{
    size_t n = data->size1;
    size_t p = data->size2;
    int i;
    
    //creates a submatrix of matrix M
    int* IndRow = new int[n];
    for (i=0; i<n; i++)
    {
        IndRow[i] = i;
    }
    int* firstColumn = new int[1];
    firstColumn[0] = 0;
    gsl_matrix* D1 = MakeSubmatrix(data, IndRow, n, firstColumn, 1);
    
    // sub 1 of the index select columns
    for (i = 0; i<lenA; i++)
        A[i]--;
    
    gsl_matrix* D_A = MakeSubmatrix(data, IndRow, n, A, lenA);
    
    gsl_matrix* D_At = transposematrix(D_A);
    
    gsl_matrix* M_A = gsl_matrix_alloc(lenA, lenA);
    matrixproduct(D_At, D_A, M_A);
    
    gsl_matrix* diagA = diagMatrix(lenA);
    //printmatrix("diag.txt", diagA);
    gsl_matrix_add(M_A, diagA);
    //printmatrix("M_A.txt",M_A);
    
    gsl_matrix* temp2 = gsl_matrix_alloc(1, 1);
    matrixproduct(transposematrix(D1), D1, temp2);
    
    gsl_matrix* temp3 = gsl_matrix_alloc(1, lenA);
    matrixproduct(transposematrix(D1), D_A, temp3);
    
    gsl_matrix* temp4 = gsl_matrix_alloc(1, lenA);
    matrixproduct(temp3, inverse(M_A), temp4);
    
    gsl_matrix* temp5 = gsl_matrix_alloc(1, n);
    matrixproduct(temp4, D_At, temp5);
    
    gsl_matrix* temp6 = gsl_matrix_alloc(1, 1);
    matrixproduct(temp5, D1, temp6);
    
    double last_term = 1 + gsl_matrix_get(temp2,0,0) - gsl_matrix_get(temp6,0,0);
    
    double log_marglik = lgamma((n+lenA+2)/2.0) - lgamma((lenA+2)/2.0) - 0.5 * logdet(M_A) - (n+lenA+2)/2.0 * log(last_term);
    
    gsl_matrix_free(D1);
    gsl_matrix_free(D_A);
    gsl_matrix_free(D_At);
    gsl_matrix_free(M_A);
    gsl_matrix_free(diagA);
    gsl_matrix_free(temp2);
    gsl_matrix_free(temp3);
    gsl_matrix_free(temp4);
    gsl_matrix_free(temp5);
    gsl_matrix_free(temp6);
    
    return (log_marglik);
    
}

// calculate the cofactor of m_(1,delj)
gsl_matrix* cofactor(gsl_matrix* m, size_t n, int delj)
{
    int i,j;
    gsl_matrix* cofMat = gsl_matrix_alloc(n-1,n-1);
    //left part
    for(i=1; i<n; i++)
    {
        for(j=0; j<delj; j++)
        {
            gsl_matrix_set(cofMat,i-1,j, gsl_matrix_get(m,i,j));
        }
    }
    
    // right part
    for(i=1; i<n; i++)
    {
        for(j=delj+1; j<n; j++)
        {
            gsl_matrix_set(cofMat,i-1,j-1, gsl_matrix_get(m,i,j));
        }
    }

    return (cofMat);

}

// det
double det(gsl_matrix* K)
{
    size_t n = K->size1;
    double detResult = 0;
    if (n==1){return (gsl_matrix_get(K,0,0));}
    if (n==2)
    {
        detResult = gsl_matrix_get(K,0,0)*gsl_matrix_get(K,1,1)-gsl_matrix_get(K,0,1)*gsl_matrix_get(K,1,0);
        return (detResult);
    }
    if (n>=3)
    {
        //double detSum = 0;
        for (int j=0; j<K->size2; j++)
        {
            gsl_matrix* cofMat = gsl_matrix_alloc(n-1, n-1);
            cofMat = cofactor(K, n, j);
            detResult += gsl_matrix_get(K,0,j)*pow(-1, j)*det(cofMat);
        }
    }
    
    return detResult;
}

// return the size of the matrix
int readLines(char datafilename[])
{
    FILE *pf = fopen(datafilename, "r"); //
    char buf[1000];
    int lineCnt = 0;
    if (!pf) // if success
        return -1;
    while (fgets(buf, 1000, pf)) // fgets
        lineCnt++; //
    fclose(pf);
    //printf("file line count = %d\n", lineCnt);
    return lineCnt;
}

// get the co variance matrix from m
gsl_matrix* co_variance(gsl_matrix* m)
{
    int i,j;
    size_t p = m->size2-1;
    size_t n = m->size1;
    
    gsl_matrix* coM = gsl_matrix_alloc(p, p);
    
    for(i=0;i<p;i++)
    {
        for(j=0;j<p;j++)
        {
            gsl_vector* v_i = gsl_vector_alloc(n);
            gsl_matrix_get_col(v_i, m, i);
            gsl_vector* v_j = gsl_vector_alloc(n);
            gsl_matrix_get_col(v_j, m, j);
            double coXY = gsl_stats_covariance (v_i->data, 1, v_j->data, 1, n);
            gsl_matrix_set(coM, i, j, coXY);
        }
    }	
    
    return (coM);
}

// cholesky decomposition
//gsl_linalg_cholesky_decomp1 (gsl_matrix * A);


#include "matrices.h"
int main()
{

    char datafilename[] = "erdata.txt";
    // get the size of the matriax
    //int lines = readLines(datafilename);
    int n = 158;
    int p = 51;
    // alloc the memory
    gsl_matrix* data = gsl_matrix_alloc(n, p);
    
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

    // calculate the deterimant of matrix
    //double detMat = det(data);
    //printf("the det of matrix is = %.5lf\n", detMat);
    
    
    // Problem 1 get the co variance matrix
    gsl_matrix* coVarMat = gsl_matrix_alloc(p-1, p-1);
    coVarMat = co_variance(data);
    printmatrix("coVarMat.txt", coVarMat);
    // FILE* output = fopen("coVarMat.txt","wb");
    // gsl_matrix_fprintf (output, coVarMat, "%f");
    // gsl_matrix_fwrite(output, coVarMat);
    //fclose(output);
    //free memory
    gsl_matrix_free(coVarMat);
    
    // Problem 2 draw independent samples
    
    
    
    gsl_matrix_free(data);
    
    
    return(1);
}

