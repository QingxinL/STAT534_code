/*
 * determinant.c
 */

#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_vector.h>
#include<stdio.h>
#include<math.h>
#include<strings.h>

// Assume that we are deleting the first row and the "colout" column.
// The resulting matrix is stored in submatrix.
void getMinor(gsl_matrix* fullmatrix, gsl_matrix* submatrix, int colout)
{
	int i, j, n = fullmatrix->size1;
	double mat_el;
	int colreach = 0;

	for(j = 0; j < n - 1 ; j++)
	{
		if(j == colout)
		{
			/* set indicator that we've reached the column to skip,
			 * this triggers a change farther down in which column we
			 * access from fullmatrix.
			 */
			colreach = 1;
		}
		// copy the column of fullmatrix
		// into submatrix, omitting the first row element.
		for(i = 0; i < n - 1; i++)
		{

			mat_el = gsl_matrix_get(fullmatrix, i + 1, j + colreach);
			gsl_matrix_set(submatrix, i, j, mat_el);
		}
	}
	return;
}

// function to compute the determinant of a matrix, assuming that it's a square matrix
double getDeterminant(gsl_matrix* matrix)
{
	double det = 0;
	int n = matrix->size1;

	if(2 == n)
	{
		// formula for a 2x2 matrix
		det = gsl_matrix_get(matrix, 0, 0)*gsl_matrix_get(matrix,1,1) -
				gsl_matrix_get(matrix,1,0)*gsl_matrix_get(matrix,1,0);
		return(det);
	}

	if(1 == n)
	{
		det = gsl_matrix_get(matrix, 1, 1);
		return(det);
	}

	int i;
	double a, pow_neg_one = -1.0;
	gsl_matrix* submatrix = gsl_matrix_alloc(n - 1,n - 1);

	for(i = 0; i < n; i++)
	{
		a = gsl_matrix_get(matrix, 0, i);
		pow_neg_one *= -1.0;

		// if the a_ij element is zero, no need to compute the i'th sub-determinant so
		// we skip this iteration using the "continue" keyword. On the particular matrix
		// provided for this assignment, this saves a TREMENDOUS amount of computation.
		if(a == 0.0)
			continue;

		// omit the first row and i'th column, stored in tempSubMatrix,
		// and get the determinant of the submatrix
		getMinor(matrix, submatrix, i);
		det += a*pow_neg_one*getDeterminant(submatrix);
	}

	gsl_matrix_free(submatrix);
	return(det);
}


int	main()
{
	int n = 10;
	char filename[]	 = "mybandedmatrix.txt";
	double det = 0;
	FILE* dat = fopen(filename, "r");
	gsl_matrix* matrix= gsl_matrix_alloc(n,n);

	if(NULL == dat)
	{
		fprintf(stderr, "Cannot open data file [%s] \n", filename);
		return(0);
	}
	if(0 != gsl_matrix_fscanf(dat, matrix))
	{
		fprintf(stderr, "Failed to read matrix from file [%s] \n", filename);
		return(0);
	}

	fclose(dat);

	det = getDeterminant(matrix);

	printf("The determinant of your matrix is %.3lf \n", det);

	gsl_matrix_free(matrix);
	return(1);
}

