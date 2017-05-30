/*
 * GSL version
 *
 */
#include"matrices.h"

void subset_gsl_matrix(gsl_matrix* fulldata, gsl_matrix* subdata,
		int n, int* A, int lenA);
gsl_matrix* create_MA(gsl_matrix* D_A);

void subset_gsl_matrix(gsl_matrix* fulldata, gsl_matrix* subdata,
		int n, int* A, int lenA)
{
	gsl_vector* tempvec = gsl_vector_alloc(n);
	int i;

	//set the columns of the matrix.
	for(i = 0; i < lenA; i++)
	{
		gsl_matrix_get_col(tempvec, fulldata, (A[i]-1));
		gsl_matrix_set_col(subdata, i, tempvec);
	}
	// clean the memory
	gsl_vector_free(tempvec);
	return;
}

gsl_matrix* create_MA(gsl_matrix* D_A)
{
	int j = 0;
	int len_A = D_A->size2;
	int n = D_A->size1;
	double x;

	gsl_matrix* D_A_transpose = gsl_matrix_alloc(len_A, n);
	gsl_matrix_transpose_memcpy(D_A_transpose, D_A);

	gsl_matrix* M_A = gsl_matrix_alloc(len_A, len_A);

	// Here, compute (X^T X), store in matrix MA. The function gsl_blas_dgemm
	// is a general matrix multiplication algorithm, and is several times faster than
	// the usual nested for-loop routine that you might use naively.

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, D_A_transpose, D_A, 0.0, M_A);
	gsl_matrix_free(D_A_transpose);

	for(j = 0; j < len_A;j ++) // Add the identity matrix
	{
		x = gsl_matrix_get(M_A, j, j) + 1.0;
		gsl_matrix_set(M_A, j , j , x);
	}

	return(M_A);
}

double compute_ml_matprod(gsl_vector* response, gsl_matrix* M_A, gsl_matrix* D_A)
{
	double matprod = 0.0;
	int lenA = D_A->size2;
	int n = D_A->size1;

	// allocate temporary vectors for intermediate calculations
	gsl_vector* tempvec = gsl_vector_calloc(lenA);
	gsl_vector* tempvec2 = gsl_vector_calloc(lenA);
	gsl_vector* tempvec3 = gsl_vector_calloc(n);
	gsl_matrix* choldecomp = gsl_matrix_alloc(lenA, lenA);
	gsl_matrix_memcpy(choldecomp, M_A);

	// product of t(D_A) %*% D_1, stored in tempvec
	gsl_blas_dgemv(CblasTrans, 1.0, D_A, response, 0.0, tempvec);

	// next two lines, tempvec2 stores inv(MA) %*% t(D_A) %*% D_1,
	// computed via the cholesky decomposition, since MA is pos-definite.
	gsl_linalg_cholesky_decomp(choldecomp);
	gsl_linalg_cholesky_solve(choldecomp, tempvec, tempvec2);

	// product of D_A %*% inv(MA) %*% t(D_A) %*% D_1, which is a vector of length
	// n, stored in tempvec3
	gsl_blas_dgemv(CblasNoTrans, 1.0, D_A, tempvec2, 0.0, tempvec3);

	// finally, t(D_1) %*% D_A %*% inv(MA) %*% t(D_A) %*% D_1
	gsl_blas_ddot(response, tempvec3, &matprod);

	gsl_matrix_free(choldecomp);
	gsl_vector_free(tempvec);
	gsl_vector_free(tempvec2);
	gsl_vector_free(tempvec3);

	return(matprod);
}

double marglik(gsl_matrix* data, int lenA, int* A)
{
	// declare variables for intermediate calculations
	double gam_coef = 0.0;
	double l_det = 0.0;
	double matprod = 0.0;
	double inner_d = 0.0;
	double total = 0.0;
	int n = data->size1;

	// extract submatrix D_A and create the matrix M_A
	gsl_matrix* D_A = gsl_matrix_alloc(n, lenA);
	subset_gsl_matrix(data, D_A, n, A, lenA);
	gsl_matrix* M_A = create_MA(D_A);

	// extract response vector
	gsl_vector* response = gsl_vector_alloc(n);
	gsl_matrix_get_col(response, data, 0);


	//get the log gamma piece
	gam_coef = lgamma((n + 2.0 + lenA)/2.0) - lgamma((lenA + 2.0)/2.0);
	
	// get the log determinant of MA
	l_det = logdet(M_A);

	// get the inner product of the response, ie compute (y^t y)
	gsl_blas_ddot(response, response, &inner_d);
	
	// compute the long matrix product
	matprod = compute_ml_matprod(response, M_A, D_A);
	
	total = gam_coef - .5*l_det - ((n + lenA + 2.0)/2.0)*log(1.0 + inner_d - matprod);

	// deallocate memory
	gsl_matrix_free(D_A);
	gsl_vector_free(response);
	gsl_matrix_free(M_A);

	return(total);
}

int main()
{
	int n = 158; //sample size
	int p = 51; //number of variables
	int i;
	int A[] = {2,5,10};//indices of the variables present in the regression
	int lenA = 3; //number of indices
	char datafilename[] = "erdata.txt";
	printf("HW 5 problem 2: GSL version \n");
	//allocate the data matrix
	gsl_matrix* data = gsl_matrix_alloc(n,p);

	//read the data
	FILE* datafile = fopen(datafilename,"r");
	if(NULL==datafile)
	{
		fprintf(stderr,"Cannot open data file [%s]\n",datafilename);
		return(0);
	}
	if(0!=gsl_matrix_fscanf(datafile,data))
	{
		fprintf(stderr,"File [%s] does not have the required format.\n",datafilename);
		return(0);
	}
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
