#include "matrices.h"

//allocates the memory for a matrix with 
//n rows and p columns
double ** allocmatrix(int n,int p)
{
	int i;
	double** m;
	
	m = new double*[n];
	for(i=0;i<n;i++)
	{
		m[i] = new double[p];
		memset(m[i],0,p*sizeof(double));
	}
	return(m);
}

//frees the memory for a matrix with n rows
void freematrix(int n,double** m)
{
	int i;
	
	for(i=0;i<n;i++)
	{
		delete[] m[i]; m[i] = NULL;
	}
	delete[] m; m = NULL;
	return;
}

//creates the copy of a matrix with n rows and p columns
void copymatrix(int n,int p,double** source,double** dest)
{
	int i,j;
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<p;j++)
		{
			dest[i][j] = source[i][j];
		}
	}
	return;
}

//reads from a file a matrix with n rows and p columns
void readmatrix(char* filename,int n,int p,double* m[])
{
	int i,j;
	double s;
	FILE* in = fopen(filename,"r");
	
	if(NULL==in)
	{
		printf("Cannot open input file [%s]\n",filename);
		exit(1);
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<p;j++)
		{
			fscanf(in,"%lf",&s);
			m[i][j] = s;
		}
	}
	fclose(in);
	return;
}

//prints the elements of a matrix in a file
void printmatrix(char* filename,int n,int p,double** m)
{
	int i,j;
	double s;
	FILE* out = fopen(filename,"w");
	
	if(NULL==out)
	{
		printf("Cannot open output file [%s]\n",filename);
		exit(1);
	}
	for(i=0;i<n;i++)
	{
		fprintf(out,"%.3lf",m[i][0]);
		for(j=1;j<p;j++)
		{
			fprintf(out,"\t%.3lf",m[i][j]);
		}
		fprintf(out,"\n");
	}
	fclose(out);
	return;
}

//creates the transpose of the matrix m
double** transposematrix(int n,int p,double** m)
{
	int i,j;
	
	double** tm = allocmatrix(p,n);
	
	for(i=0;i<p;i++)
	{
		for(j=0;j<n;j++)
		{
			tm[i][j] = m[j][i];
		}
	}	
	
	return(tm);
}

//calculates the dot (element by element) product of two matrices m1 and m2
//with n rows and p columns; the result is saved in m
void dotmatrixproduct(int n,int p,double** m1,double** m2,double** m)
{
	int i,j;
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<p;j++)
		{
			m[i][j] = m1[i][j]*m2[i][j];
		}
	}
	
	return;
}

//calculates the product of a nxp matrix m1 with a pxl matrix m2
//returns a nxl matrix m
void matrixproduct(int n,int p,int l,double** m1,double** m2,double** m)
{
	int i,j,k;
	double s;
	
	for(i=0;i<n;i++)
	{
		for(k=0;k<l;k++)
		{
			s = 0;
			for(j=0;j<p;j++)
			{
				s += m1[i][j]*m2[j][k];
			}
			m[i][k] = s;
		}
	}
	return;
}

void set_mat_identity(int p, double *A)
{
 int i;

 for(i = 0; i < p * p; i++) A[i] = 0;
 for(i = 0; i < p; i++) A[i * p + i] = 1;
 return;
}

//computes the inverse of a symmetric positive definite matrix
void inverse(int p,double** m)
{
  int i,j,k;
  double* m_copy = (double*)malloc((p * p) * sizeof(double));
  double* m_inv = (double*)malloc((p * p) * sizeof(double));

  k=0;
  for(i=0;i<p;i++)
  {
     for(j=0;j<p;j++)
     {
        m_copy[k] = m[i][j];
        k++;
     }
  }

  set_mat_identity(p, m_inv);

  //-----  Use LAPACK  -------
  if(0!=(k=clapack_dposv(CblasRowMajor, CblasUpper, p, p, m_copy, p, m_inv, p)))
  {
    fprintf(stderr,"Something was wrong with clapack_dposv [%d]\n",k);
     exit(1);
  }
  //--------------------------

  k=0;
  for(i=0;i<p;i++)
  {
     for(j=0;j<p;j++)
     {
        m[i][j] = m_inv[k];
        k++;
     }
  }  

  free(m_copy);
  free(m_inv);

  return;
}


//computes the log of the determinant of a symmetric positive definite matrix
double logdet(int p,double** m)
{
        //just take care of the 1x1 case
        if(1==p)
	{
	  return(log(m[0][0]));
	}

	int i,j;
	char jobvl = 'N';
	char jobvr = 'N';
	int lda = p;
	double wr[2*p];
	double wi[2*p];
	double vl[p][p];
	int ldvl = p*p;
	double vr[p][p];
	int ldvr = p*p;
	double work[p*p];
	int lwork = p*p;
	double a[p][p];
	int info = 1;
	
	for(i=0;i<p;i++)
	{
		for(j=0;j<p;j++)
		{
			a[i][j] = m[i][j];
		}
	}
	dgeev_(&jobvl,&jobvr,&p,(double*)a,&lda,(double*)wr,(double*)wi,(double*)vl, 
		  &ldvl,(double*)vr,&ldvr,(double*)work,&lwork,&info);

	if(0!=info)
	{
		printf("Smth wrong in the call of 'dgeev' error is [info = %d]\n",info);
		exit(1);
	}	   
	
	double logdet = 0;
	for(i=0;i<p;i++) logdet+=log(wr[i]);	
	return(logdet);
}

//this returns the columns specified by the set of indices A
//lenA is the number of indices of A
//data is an nxp matrix
//the returned submatrix is n x lenA
double** submatrix(int n,int p,double** data,int lenA,int* A)
{
  double** result = allocmatrix(n,lenA);

  int i,j;
  for(i=0;i<lenA;i++)
  {
    for(j=0;j<n;j++)
    {
      result[j][i] = data[j][A[i]-1];
    }
  }

  return(result);
}

//computes the marginal likelihood
double marglik(int n,int p,double** data,int lenA,int* A)
{
  double result = 0.0; //value is returned by the function
  int i;

  //nx1 submatrix corresponding with the response
  int response = 1;
  double** D1 = submatrix(n,p,data,1,&response);

  //n x lenA matrix correspoding with the regressors
  double** DA = submatrix(n,p,data,lenA,A);

  //get the tranposed of DA
  double** transpDA = transposematrix(n,lenA,DA); 

  //MA is a lenA x lenA matrix
  double** MA = allocmatrix(lenA,lenA);
  matrixproduct(lenA,n,lenA,transpDA,DA,MA);
  //add the identity matrix to MA
  //i.e., increment with 1 the diagonal elements
  for(i=0;i<lenA;i++)
  {
    MA[i][i] += 1;
  }

  result = lgamma((n+lenA+2.0)/2.0)-lgamma((lenA+2.0)/2.0)-0.5*logdet(lenA,MA);

  //calculate 1 + t(D1) %*% D1
  double s = 1.0;
  for(i=0;i<n;i++)
  {
    s += pow(D1[i][0],2);
  }
 
  //calculate t(DA) %*% D1
  //this is a lenA x 1 matrix
  double** transpDAD1 = allocmatrix(lenA,1);
  matrixproduct(lenA,n,1,transpDA,D1,transpDAD1);

  //calculate t(D1) %*% DA
  //just transpose transDAD1
  //this is a 1 x lenA matrix
  double** transpD1DA = transposematrix(lenA,1,transpDAD1);

  //get the inverse of MA
  inverse(lenA,MA);

  //need some place to store results
  double** temp = allocmatrix(lenA,1);
  matrixproduct(lenA,lenA,1,MA,transpDAD1,temp);

  //this is the final matrix multiplication
  for(i=0;i<lenA;i++)
  {
    s -= transpD1DA[0][i]*temp[i][0];
  }
  result -= 0.5*(n+lenA+2)*log(s);

  //clean memory
  freematrix(n,D1);
  freematrix(n,DA);
  freematrix(lenA,transpDA);
  freematrix(lenA,MA);
  freematrix(lenA,transpDAD1);
  freematrix(1,transpD1DA);
  freematrix(lenA,temp);

  return(result);
}
