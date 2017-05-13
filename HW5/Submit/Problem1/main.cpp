//
//  main.cpp
//  Hw5
//
//  Created by Yuxuan Cheng on 5/9/17.
//  Copyright © 2017 Yuxuan Cheng. All rights reserved.
//

#include <iostream>
#include <math.h>
#include "matrices.h"

double marglik(int n,int p,double** data,int lenA,int* A)
{
    double** D1 = allocmatrix(n,1);
    copymatrix(n, 1, data, D1);
    //D1 = selectData(n, data, 1, SelectA);
//    printmatrix("D1.txt", n, 1, D1);
    
    double** D_A = selectData(n, data, lenA, A);//allocmatrix(n, lenA);
    //copymatrix(n, , data, D_A);
    //D_A= selectData(n, data, lenA, A);
//    printmatrix("D_A.txt", n, lenA, D_A);

    double** D_At = transposematrix(n, lenA, D_A);
    // D_At = transposematrix(n, lenA, D_A);
//    printmatrix("D_At.txt", lenA, n, D_At);
    
    double** temp1 = allocmatrix(lenA, lenA);
    matrixproduct(lenA, n, lenA, D_At, D_A, temp1);
//    printmatrix("temp1.txt", lenA, lenA, temp1);
    
    double** M_A = allocmatrix(lenA, lenA);
    M_A = matrixplus(lenA, lenA, set_mat_identity2D(lenA), temp1);
//    printmatrix("M_A.txt", lenA, lenA, M_A);
    
    double** temp2 = allocmatrix(1, 1);
    matrixproduct(1, n, 1, transposematrix(n, 1, D1), D1, temp2);
//    printmatrix("temp2.txt", 1, 1, temp2);
    
    double** temp3 = allocmatrix(1, lenA);
    matrixproduct(1, n, lenA, transposematrix(n, 1, D1), D_A, temp3);
    
    double** temp4 = allocmatrix(1, lenA);
    copymatrix(lenA, lenA, M_A, temp1);
    inverse(lenA, temp1);
    matrixproduct(1, lenA, lenA, temp3, temp1, temp4);
    
    double** temp5 = allocmatrix(1, n);
    matrixproduct(1, lenA, n, temp4, D_At, temp5);
    
    
    double** temp6 = allocmatrix(1, 1);
    matrixproduct(1, n, 1, temp5, D1, temp6);
//    printmatrix("temp6.txt", 1, 1, temp6);
    
    double last_term = 1 + temp2[0][0] - temp6[0][0];
    double log_marglik = lgamma((n+lenA+2)/2.0) - lgamma((lenA+2)/2.0) - 0.5 * logdet(lenA, M_A) - (n+lenA+2)/2.0 * log(last_term);
    
    freematrix(1, D1);
    freematrix(n, D_A);
    freematrix(lenA, D_At);
    freematrix(lenA, M_A);
    freematrix(lenA, temp1);
    freematrix(1, temp2);
    freematrix(1, temp3);
    freematrix(1, temp4);
    freematrix(1, temp5);
    freematrix(1, temp6);
    return log_marglik;
}



int main()
{
    int n = 158; //sample size
    int p = 51; //number of variables
    int i;
    int A[] = {2,5,10};//indices of the variables present in the regression
    int lenA = 3; //number of indices
    char datafilename[] = "erdata.txt";
    //allocate the data matrix
    double** data = allocmatrix(n,p);
    //read the data
    readmatrix(datafilename,n,p,data);
    printf("Marginal likelihood of regression [1|%d",A[0]);
    for(i=1;i<lenA;i++)
    {
        printf(",%d",A[i]);
    }
    printf("] = %.5lf\n",marglik(n,p,data,lenA,A));
    //free memory
    freematrix(n,data);
    return(1);
}


