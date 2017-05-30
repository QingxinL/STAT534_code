#include "vectors.h"

//dynamically allocates a vector of double of length n
double* allocvector(int n)
{
	int i;
	
	//that's where the memory gets allocated
	double * v = new double[n];
	
	//we set the elements of our vector to zero
	
	/*
	for(i=0;i<n;i++)
	{
		v[i] = 0;
	}
	 */
	//the code above is slow
	//it's better and faster to use
	memset(v,0,n*sizeof(double));
	
	
	return(v);
}

//prints the elements of a vector of length n
void printvector(double* v,int n)
{
	int i;
	
	for(i=0;i<n;i++)
	{
		printf("element [%d] = %.3lf\n",i+1,v[i]);
	}
	
	return;
}

//frees the memory of this vector
void freevector(double* v)
{
	delete[] v; v = NULL; 
	return;
}

//calculates the dot product ofvectors v1 and v2 and saves the result in v
void vectorproduct(int n,double* v1,double* v2,double* v)
{
	int i;
	
	for(i=0;i<n;i++)
	{
		v[i] = v1[i]*v2[i];
	}
	
	return;
}