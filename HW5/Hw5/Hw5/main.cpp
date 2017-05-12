//
//  main.cpp
//  Hw5
//
//  Created by Yuxuan Cheng on 5/9/17.
//  Copyright © 2017 Yuxuan Cheng. All rights reserved.
//

#include <iostream>
#include <math.h>
/*
int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}
*/
#include "matrices.h"

double marglik(int n,int p,double** data,int lenA,int* A)
{
    double** D1 = allocmatrix(n,1);
    copymatrix(n, 1, data, D1);
    //D1 = selectData(n, data, 1, SelectA);
    
    double** D_A = allocmatrix(n, lenA);
    //copymatrix(n, , data, D_A);
    D_A= selectData(n, data, lenA, A);

    double** D_At = allocmatrix(lenA, n);
    D_At = transposematrix(n, lenA, D_A);
    
    double** M_A = allocmatrix(lenA, lenA);
    double** temp1 = allocmatrix(lenA, lenA);
    matrixproduct(lenA, n, lenA, D_At, D_A, temp1);
    M_A = matrixplus(lenA, lenA, set_mat_identity2D(lenA), temp1);
    
    double** temp2 = allocmatrix(1, 1);
    matrixproduct(1, n, 1, transposematrix(n, 1, D1), D1, temp2);
    
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
    
    double last_term = 1 + temp2[0][0] - temp6[0][0];
    
    double log_marglik = lgamma((n+lenA+2)/2.0) - lgamma((n+2)/2.0) - 0.5 * logdet(lenA, M_A) - (n+lenA+2)/2.0 * log(last_term);
    
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


