/////////////////////////////////////////////////////////
// HOMEWORK 7 - COVARIANCE MATRICES AND DIRECT SAMPLING//
// FROM MULTIVARIATE NORMAL DISTRIBUTION               //
/////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_errno.h>

//prints line-by-line a GSL matrix in a file
void printGSLmatrix(FILE* stream,gsl_matrix* m)
{
  int i,j;

  for(i=0;i<m->size1;i++)
  {
    fprintf(stream,"%.5lf",gsl_matrix_get(m,i,0));
    for(j=1;j<m->size2;j++)
    {
      fprintf(stream,"\t%.5lf",gsl_matrix_get(m,i,j));
    }
    fprintf(stream,"\n");
  }

  return;
}

//this function makes the column sums
//of a matrix equal to zero
void Center(gsl_matrix* m)
{
  int j;
  double columnMean;
       
  for(j=0;j<m->size2;j++)
  {
    //obtain a column
    gsl_vector_view mycolumn = gsl_matrix_column(m,j);
    //calculate the mean of that column
    columnMean =  gsl_stats_mean(mycolumn.vector.data,
                                 mycolumn.vector.stride,
				 m->size1);
    //substract the column mean
    gsl_vector_add_constant(&(mycolumn.vector),-columnMean);
  }
  return;
}


//calculates U = t(X-mean(X)) * (X-mean(X))
void makeCovariance(gsl_matrix* covX,gsl_matrix* X)
{
	int i;
	double* s;
	
	//create a copy of the matrix X
	gsl_matrix* Xcopy = gsl_matrix_alloc(X->size1,X->size2);
	if(GSL_SUCCESS!=gsl_matrix_memcpy(Xcopy,X))
	{
		printf("GSL failed to copy a matrix.\n");
		exit(1);
	}
	Center(Xcopy);

	gsl_blas_dgemm (CblasTrans, CblasNoTrans,
			1.0, Xcopy, Xcopy,
			0.0, covX);

	gsl_matrix_scale(covX,1.0/((double)(X->size1-1)));

	gsl_matrix_free(Xcopy);
	return;
}

//creates the Cholesky decomposition of a matrix
gsl_matrix* makeCholesky(gsl_matrix* K)
{
	int i,j;
	
	gsl_matrix* Phi = gsl_matrix_alloc(K->size1,K->size1);
	if(GSL_SUCCESS!=gsl_matrix_memcpy(Phi,K))
	{
		printf("GSL failed to copy a matrix.\n");
		exit(1);
	}
	if(GSL_SUCCESS!=gsl_linalg_cholesky_decomp(Phi))
	{
		printf("GSL failed Cholesky decomposition.\n");
		exit(1);
	}
	for(i=0;i<Phi->size1;i++)
	{
		for(j=i+1;j<Phi->size2;j++)
		{
			gsl_matrix_set(Phi,i,j,0.0);
		}
	}
	
	return(Phi);
}

//samples from the multivariate normal distribution N(0,Sigma)
//the samples are saved in the matrix "Samples"
void randomMVN(gsl_rng* mystream,gsl_matrix* Samples,gsl_matrix* Sigma)
{
  int i;
  gsl_matrix* Psi = makeCholesky(Sigma);
  gsl_matrix* Z = gsl_matrix_alloc(Sigma->size1,1);
  gsl_matrix* X = gsl_matrix_alloc(Sigma->size1,1);

  for(int asample=0;asample<Samples->size1;asample++)
  {
    for(i=0;i<Sigma->size1;i++)
    {
      gsl_matrix_set(Z,i,0,gsl_ran_ugaussian(mystream));
    }
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
		    1.0, Psi, Z,
		    0.0, X);
    //record the sample we just generated
    for(i=0;i<Sigma->size1;i++)
    {
      gsl_matrix_set(Samples,asample,i,
		     gsl_matrix_get(X,i,0));
    }
  }

  //free memory
  gsl_matrix_free(Psi);
  gsl_matrix_free(Z);
  gsl_matrix_free(X);
  return;
}


int main()
{
  int i,j;
  FILE* file = NULL;
  
  clock_t start, stop;

  start = clock();

  int n = 158; //sample size
  int p = 51; //number of variables
  char datafilename[] = "erdata.txt"; //name of the data file

  //this is the number of samples we draw from 
  //the multivariate normal
  int NumberOfSamples = 1000000;

  //the name of the file where we save
  //the sample covariance matrix
  char covariancefilename[] = "samplecovariance.txt"; 

  //the name of the file where we save the estimate
  //of the sample covariance matrix
  char estimatefilename[] = "estimatedcovariance.txt";

  //allocate the data matrix
  gsl_matrix* data = gsl_matrix_alloc(n,p);
  //storage for the sample covariance matrix
  gsl_matrix* covMatrix = gsl_matrix_alloc(p,p);

  //storage for the samples we generate from the multivariate normal
  gsl_matrix* samplesFromMVN = gsl_matrix_alloc(NumberOfSamples,p);

  //our estimate of the sample covariance matrix will be stored here
  gsl_matrix* estimateMatrix = gsl_matrix_alloc(p,p);
  
  //open for reading the data file
  file = fopen(datafilename,"r");
  if(NULL==file)
  {
    fprintf(stderr,"Cannot open data file [%s].\n",
            datafilename);
    return(0);
  }

  //read the data
  if(GSL_SUCCESS!=gsl_matrix_fscanf(file,data))
  {
    fprintf(stderr,"Failed to read the data.\n");
    exit(1);
  }
  //close the data file
  fclose(file);

  //calculate the sample covariance matrix
  makeCovariance(covMatrix,data);
  //open for writing the file where we save
  //the sample covariance matrix
  file = fopen(covariancefilename,"w");
  if(NULL==file)
  {
    fprintf(stderr,"Cannot open data file [%s].\n",
            covariancefilename);
    return(0);
  } 
  printGSLmatrix(file,covMatrix);
  fclose(file);

  //initialize the GSL RNG
  const gsl_rng_type* T;
  gsl_rng* r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  //sample from the multivariate normal
  randomMVN(r,samplesFromMVN,covMatrix);

  //calculate the sample covariance matrix
  //associated with the samples we just generated
  makeCovariance(estimateMatrix,samplesFromMVN);
  //open for writing the file where we save
  //our estimate
  file = fopen(estimatefilename,"w");
  if(NULL==file)
  {
    fprintf(stderr,"Cannot open data file [%s].\n",
            estimatefilename);
    return(0);
  } 
  printGSLmatrix(file,estimateMatrix);
  fclose(file);

  //free memory
  gsl_rng_free(r);
  gsl_matrix_free(data);
  gsl_matrix_free(covMatrix);
  gsl_matrix_free(samplesFromMVN);
  gsl_matrix_free(estimateMatrix);

  stop = clock();
  fprintf(stderr,"Program took %.2lf seconds\n",
	  (double)(stop-start)/CLOCKS_PER_SEC);

  return(1);
}
