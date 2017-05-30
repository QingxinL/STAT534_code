/*
 * MatricesHW5.c
 */

#include<math.h>
#include"matrices.h"
#include<stdio.h>

/* subsetmatrix: function to extract predictor variables from the data matrix
 *
 * nrow: number of rows in both matrices
 * vars: columns of matrix to be extracted
 * npred: number of predictor variables to subset (i.e. length of vars)
 * */

void subsetDataMatrix(double** fulldata, double** subdata, int npred, int nrow, int* vars)
{
	int i,j;
	for(i = 0; i < npred; i++)
	{
		for(j = 0; j < nrow; j++)
		{
			subdata[j][i] = fulldata[j][(vars[i] - 1)]; //copy the predictor variables
		}
	}
	return;
}

double marglik(int n,int p, double** data, int lenA, int* A)
{
	double** D_A = allocmatrix(n, lenA);
	double** response = allocmatrix(n,1);

	int i,j,k;

	subsetDataMatrix(data, D_A, lenA, n, A);

	// fill the response "matrix"
	for(i = 0; i < n; i++)
	{
		response[i][0] = data[i][0];
	}

	double** D_A_transpose = transposematrix(n, lenA, D_A);
	double** response_transpose = transposematrix(n, 1, response);

	double gam_coef;
	double l_det;
	double matprod;
	double inner_d = 0.0;
	double total;

	//get the funky log gamma piece
	gam_coef = lgamma(((double)n + 2.0 + (double)lenA)/2.0) - lgamma(((double)lenA + 2.0)/2.0);

	double** M_A = allocmatrix(lenA, lenA);
	matrixproduct(lenA, n, lenA, D_A_transpose, D_A, M_A);

	for(j = 0; j < lenA;j ++) // Add the identity matrix
	{
		M_A[j][j] += 1;
	}

	l_det = logdet(lenA, M_A);

	inverse(lenA, M_A); // invert MA. MA will be destroyed, but we don't need it anymore.

	// get the inner product of the response
	for(k = 0; k < n; k++)
	{
		inner_d += response[k][0]*response[k][0];
	}

	// matrices for intermediate calculations
	double** work1 = allocmatrix(1,lenA);
	double** work2 = allocmatrix(1,lenA);
	double** work3 = allocmatrix(1, n);
	double** work4 = allocmatrix(1,1);

	// incrememntally compute t(D_1)%*%(D_A)%*%inv(M_A)%*%t(D_A)%*%D_1
	matrixproduct(1, n, lenA, response_transpose, D_A, work1);
	matrixproduct(1, lenA, lenA, work1, M_A, work2);
	matrixproduct(1, lenA, n, work2, D_A_transpose, work3);
	matrixproduct(1, n, 1, work3, response, work4);

	// technically work4 is a 1x1 matrix, we need to extract its element.
	matprod = work4[0][0];

	total = gam_coef - .5*l_det - ((n + lenA + 2.0)/2.0)*log(1.0 + inner_d - matprod);

	// delete matrices from memory
	freematrix(n, D_A);
	freematrix(n, response);
	freematrix(lenA, D_A_transpose);
	freematrix(1, response_transpose);
	freematrix(lenA, M_A);
	freematrix(1, work1);
	freematrix(1, work2);
	freematrix(1, work3);
	freematrix(1, work4);

	return(total);
}

int main()
{
	int n = 158;
	int p = 51;
	int i;

	int A[] = {2,5,10};
	int lenA = 3;
	char datafilename[] = "erdata.txt";

	double** data = allocmatrix(n,p);
	readmatrix(datafilename, n, p, data);

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

