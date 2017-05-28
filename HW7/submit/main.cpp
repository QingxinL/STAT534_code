//
//  main.cpp
//  Problem1
//
//  Created by Yuxuan Cheng on 5/11/17.
//  Copyright © 2017 Yuxuan Cheng. All rights reserved.
//

#include <iostream>
#include "matrices.h"

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
    size_t p = m->size2;
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

// cholesky decomposition of matriax
gsl_matrix* choleskyComp(gsl_matrix* m)
{
    size_t n = m->size1;
    size_t p = m->size2;
    if (n!=p)
    {
        printf("matrix dimension dismatch\n");
        exit(1);
    }

    // get a copy of m
    gsl_matrix* copyM = gsl_matrix_alloc(n, n);
    if(GSL_SUCCESS!=gsl_matrix_memcpy(copyM, m))
    {
        printf("GSL failed to copy a matrix.\n");
        exit(1);
    }

    // set the upper trangular matrix 0
    for (int i=0; i<n; i++)
        for (int j=i+1; j<n; j++)
        {
            gsl_matrix_set(copyM, i, j, 0);
        }

    // calculate the cholesky decomposition
    gsl_linalg_cholesky_decomp(copyM); // for the cluster, but the GNU recommand gsl_linalg_cholesky_decomp1
    //gsl_linalg_cholesky_decomp1(copyM); // for local machine https://www.gnu.org/software/gsl/manual/html_node/Cholesky-Decomposition.html

    return (copyM);
}

// get the gaussian random matriax, *r is the random seed
gsl_matrix* gaussianRanM(int p, gsl_rng *r)
{
    gsl_matrix* GaussMat = gsl_matrix_alloc(p, 1);

    for (size_t i=0; i<p; i++)
    {
        double randTemp = gsl_ran_ugaussian(r);
        gsl_matrix_set(GaussMat, i, 0, randTemp);
    }

    return (GaussMat);
}

// draw the matriax
gsl_matrix* drawTheMat(gsl_matrix* cho, gsl_matrix* Z)
{
    size_t p = Z->size1;
    gsl_matrix* X = gsl_matrix_alloc(p, 1);
    matrixproduct(cho, Z, X);
    return (X);
}

// Get multiple independent samples and return the co variance matriax of sample matriax
gsl_matrix* multipleSample(gsl_matrix* choleskyMat, int sampleNum, int p)
{
    gsl_matrix* sampleMat = gsl_matrix_alloc(p, sampleNum);

    //Initialize the GSL generator
    gsl_rng *r = gsl_rng_alloc(gsl_rng_taus2);

    // Get the 10000 samples
    for (int j=0; j<sampleNum; j++)
    {
        gsl_matrix* randGaussMat = gaussianRanM(p, r);

        // Draw the Matriax
        gsl_matrix* drawMat = drawTheMat(choleskyMat, randGaussMat);
        gsl_vector* vec = gsl_vector_alloc(p);
        gsl_matrix_get_col (vec, drawMat, 0);
        gsl_matrix_set_col(sampleMat, j, vec);

        // free the memory
        gsl_vector_free(vec);
        gsl_matrix_free(randGaussMat);
        gsl_matrix_free(drawMat);
    }
    // printmatrix("10000SampleMat.txt", sampleMat);

    // transpose the sampleMat and then get the co variance matrix of sampleMat
    gsl_matrix* sampleMatT = transposematrix(sampleMat);
    gsl_matrix* covarMat = co_variance(sampleMatT);

    // free the memory
    gsl_rng_free (r);
    gsl_matrix_free(sampleMat);
    gsl_matrix_free(sampleMatT);

    return (covarMat);
}


int main()
{

    char datafilename[] = "erdata.txt";
    // get the size of the matriax
    //int lines = readLines(datafilename);
    int n = 158;
    int p = 51;

    // alloc the memory
    gsl_matrix* data = gsl_matrix_alloc(n, p);

    // read the data
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

    // Step 1 get the co variance matrix
    gsl_matrix* coVarMat = gsl_matrix_alloc(p, p);
    coVarMat = co_variance(data);
    printmatrix("matrix1.txt", coVarMat);

    // Step 2 get the cholesky decomposition of matriax
    gsl_matrix* choleskyMat = choleskyComp(coVarMat);
    //printmatrix("choleskyMat.txt", choleskyMat);

    // step 3
    gsl_matrix* mulSampleMat = multipleSample(choleskyMat, 10000, p);
    printmatrix("matrix1.3.txt.", mulSampleMat);

    //free memory
    gsl_matrix_free(mulSampleMat);
    gsl_matrix_free(coVarMat);
    gsl_matrix_free(choleskyMat);
    gsl_matrix_free(data);

    return(1);
}
