//////////////////////////
//PROBLEM 2 - HOMEWORK 7//
//////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

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

#include "vectors.h"

/////////////////////////////////////////////////////
//MPI MESSAGES TO THE SLAVES                       //
//the master process can ask the slaves to estimate//
//another pi=P(X=i) or to die.                     //
/////////////////////////////////////////////////////
#define GETPI 1
#define DIETAG 0

//////////////////////
//CALCULATES LOG(N!)//
//////////////////////
double logfact(double n)
{
	return(lgamma(n+1));
}

///////////////////////////////////////////////////
//CALCULATES THE LOG OF THE BINOMIAL COEFFICIENTS//
///////////////////////////////////////////////////
double logchoose(double n,double k)
{
	return(logfact(n)-logfact(k)-logfact(n-k));
}

////////////////////////////////////////////////
//BINARY SEARCH TO DETERMINE THE INDEX SAMPLED//
////////////////////////////////////////////////
int WeightedSampling(int maxk,double* weights,gsl_rng* mystream)
{
	double s;
	int i;
	int q1 = 0;
	int q2 = maxk;
	
	s = gsl_rng_uniform_pos(mystream);
	while(q1+1!=q2)
	{
		int q = (q1+q2)/2;
		if(weights[q]<s)
		{
			q1 = q;
		}
		else
		{
			q2 = q;
		}
	}
	return(q1);  
}

///////////////////////////////////////////////
//CALCULATE THE CUMULATIVE PROBABILITIES USED//
//IN DIRECT SAMPLING                         //
///////////////////////////////////////////////
double* NormalizeWeights(double* wOriginal,int nmax)
{
        int i;
	double* w = new double[nmax];
	double* cumw = new double[nmax+1];
	
	memcpy(w,wOriginal,nmax*sizeof(double));

	double maxw = w[0];
	for(i = 1; i < nmax; i++)
	{ 
		if(maxw < w[i])
		{ 
			maxw = w[i];
		} 
	}
	
	for(i = 0; i < nmax; i++)
	{ 
		w[i] -= maxw;
	}
	
	for(i = 0; i < nmax; i++)
	{
		w[i] = exp(w[i]);
	}
	
	cumw[0] = 0;
	for(i = 1; i <= nmax; i++)
	{
		cumw[i] = cumw[i-1]+w[i-1];
	}
	for(i = 1; i <= nmax; i++)
	{
		cumw[i] /= cumw[nmax];
	}
	delete[] w; w = NULL;
	return cumw;
}

void master(int n)
{
  int i;
  int rank;
  int ntasks;
  int jobsRunning;
  int work[1]; //used to transmit a work request to the slaves
  double workresults[3]; //used to receive results from the slaves
  MPI_Status status;	// MPI information
  
  //calculate the actual values of pi = P(X=i)
  double* probs = allocvector(n+1);		
  //calculate the log of binomial probabilities
  for(i=0;i<=n;i++)
  {
    probs[i] = logchoose(n,i);
  }
  for(i=0;i<=n;i++)
  {
    probs[i] = exp(probs[i]-n*log(2));
  }
  //this is where our estimates of pi are stored
  double* estimates = allocvector(n+1);
  //storage for the standard errors of the estimates
  double* estimatesSD = allocvector(n+1);

  // Find out how many slaves there are
  MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
  fprintf(stderr, "Total Number of processors = %d\n",ntasks);

  jobsRunning = 1;
  for(i=0;i<=n;i++)
  {
    work[0] = i; //work request to estimate pi=P(X=i)
    if(jobsRunning < ntasks) // Do we have an available processor?
    {
      // Send out a work request
      MPI_Send(&work, 	// the vector with the variable
	       1, 	// the size of the vector
	       MPI_INT,	// the type of the vector
               jobsRunning,	// the ID of the slave to use
               GETPI,	// tells the slave what to do
               MPI_COMM_WORLD); 

       printf("Master sends out work request [%d] to slave [%d]\n",
              work[0],jobsRunning);

       // Increase the # of processors in use
       jobsRunning++;
      }
      else // all the processors are in use!
      {
	 //RECEIVE RESULTS
         MPI_Recv(workresults,	// where to store the results
 		  3,		// the size of the vector
		  MPI_DOUBLE,	// the type of the vector
	 	  MPI_ANY_SOURCE,
		  MPI_ANY_TAG, 	
		  MPI_COMM_WORLD,
		  &status);
	 printf("Master has received the result of work request [%d] from slave [%d]\n",
                (int)workresults[0],status.MPI_SOURCE);

	 //STORE RESULTS
	 estimates[(int)workresults[0]] = workresults[1];
	 estimatesSD[(int)workresults[0]] = workresults[2];

	 printf("Master sends out work request [%d] to slave [%d]\n",
                work[0],status.MPI_SOURCE);
	 // Send out a new work order to the processors that just
         // returned
         MPI_Send(&work,
                  1,
                  MPI_INT,
                  status.MPI_SOURCE, // the slave that just returned
                  GETPI,
                  MPI_COMM_WORLD);
      }
  }

  /////////////////////////////////////////////////////////
  //COLLECT RESULTS FOR ALL THE OUTSTANDING WORK REQUESTS//
  /////////////////////////////////////////////////////////
  for(rank=1; rank<jobsRunning; rank++)
  {
    //RECEIVE RESULTS
    MPI_Recv(workresults,	// where to store the results
 	     3,		// the size of the vector
	     MPI_DOUBLE,	// the type of the vector
	     MPI_ANY_SOURCE,
	     MPI_ANY_TAG, 	
	     MPI_COMM_WORLD,
	     &status);
     printf("Master has received the result of work request [%d] from slave [%d]\n",
            (int)workresults[0],status.MPI_SOURCE);

     //STORE RESULTS
     estimates[(int)workresults[0]] = workresults[1];
     estimatesSD[(int)workresults[0]] = workresults[2];
  }

  printf("Master tells the slave to die\n");
  // Shut down the slave processes
  for(rank=1; rank<ntasks; rank++)
  {
    printf("Master is killing slave [%d]\n",rank);
    MPI_Send(0,
	     0,
             MPI_INT,
             rank,		// shutdown this particular node
             DIETAG,		// tell it to shutdown
	     MPI_COMM_WORLD);
  }

  ////////////////////
  //SHOW THE RESULTS//
  ////////////////////
  printf("\n\n");
  for(i=0;i<=n;i++)
  {
    printf("P[X = %d] = %.5lf is estimated to %.5lf with s.d. %.5lf\n",
	   i,probs[i],
	   estimates[i],estimatesSD[i]);
  }

  //free memory
  freevector(probs);
  freevector(estimates);
  freevector(estimatesSD);

  return;
}

void MHSampling(gsl_rng* mystream,
		int n,int NumberOfBatches,int NumberOfIterations,
		double* probs,
		gsl_vector* mysamples,
		int targetState)
{
  int i,j;
  double u;
  int currentState;
  int nextState;

  gsl_vector_set_zero(mysamples);
  for(i=0;i<NumberOfBatches;i++)
  {
    //start the chain at a randomly chose state
    currentState = gsl_rng_uniform_int(mystream,n+1);
    for(j=0;j<NumberOfIterations;j++)
    {
      //propose a new state from the uniform distribution
      //on the possible states {0,1,...,n}
      nextState = gsl_rng_uniform_int(mystream,n+1);
      if(probs[nextState]>probs[currentState])
      {
	//accept the move to the proposed state
	currentState = nextState;
      }
      else
      {
	u = gsl_rng_uniform_pos(mystream);
	if(log(u)<=probs[nextState]-probs[currentState])
	{
	  currentState = nextState;
	}
      }

      //record the occurence of the target state
      //this is the only state we keep track of
      if(targetState==currentState)
      {
	gsl_vector_set(mysamples,
		       i,
		       1+gsl_vector_get(mysamples,i));
      }
    }
  }

  gsl_vector_scale(mysamples,1.0/((double)NumberOfIterations));
  
  return;
}

void slave(int n,int NumberOfBatches,int NumberOfIterations,int slavename)
{
  int i,j,k;
  int work[1]; //used to receive a work request from the Master
  double workresults[3]; //used to transmit results back to the Master
  MPI_Status status; //used for MPI communication

  //we only store what happens with one value of the random variable X
  gsl_vector* mysamples = gsl_vector_alloc(NumberOfBatches);
  //true vector of probabilities
  double* probs = allocvector(n+1);		

  //calculate the log of binomial probabilities
  for(i=0;i<=n;i++)
  {
    probs[i] = logchoose(n,i);
  }

  //initialize the GSL RNG
  const gsl_rng_type* T;
  gsl_rng* r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  //VERY IMPORTANT: MAKE SURE THE SLAVES DO NOT GENERATE
  //THE SAME SEQUENCE OF RANDOM NUMBERS
  //each slave generates a different number of U(0,1)
  for(i=0;i<slavename*NumberOfIterations;i++)
  {
    gsl_rng_uniform(r);
  }

  //////////////////////////////////////////////////////
  //THE SLAVE LISTENS FOR INSTRUCTIONS FROM THE MASTER//
  //////////////////////////////////////////////////////
  int notDone = 1;
  while(notDone)
  {
     printf("Slave %d is waiting\n",slavename);
     MPI_Recv(&work,		// the inputs from the master
	      1,		// the size of the inputs
	      MPI_INT,		// the type of the inputs
              0,		// from the MASTER node (rank=0)
              MPI_ANY_TAG,	// any type of order is fine
              MPI_COMM_WORLD,
              &status);
      printf("Slave %d just received smth\n",slavename);

      // switch on the type of work request
      switch(status.MPI_TAG)
      {
         case GETPI:
           printf("Slave %d has received work request [%d]\n",
                  slavename,work[0]);
          
	   //sample
	   MHSampling(r,
		      n,NumberOfBatches,NumberOfIterations,
		      probs,
		      mysamples,
		      work[0]);

           workresults[0] = (double)work[0];
	   workresults[1] = gsl_stats_mean(mysamples->data,mysamples->stride,NumberOfBatches);
	   workresults[2] = gsl_stats_sd(mysamples->data,mysamples->stride,NumberOfBatches);

           // Send the results
           MPI_Send(&workresults,
                    3,
                    MPI_DOUBLE,
                    0,		// send it to the master
                    0,		// doesn't need a TAG
                    MPI_COMM_WORLD);

            printf("Slave %d finished processing work request [%d]\n",
                   slavename,work[0]);
            break;

         case DIETAG:
            printf("Slave %d was told to die\n",slavename);

         default:
            notDone = 0;
      }
  }

  //free memory
  gsl_vector_free(mysamples);
  gsl_rng_free(r);
  freevector(probs);
  return;
}

////////////////////////////////////////////////////////////////////
//The program should be called by specifying the number of batches//
//and the number of iterations per batch in the command line.     //
//As such, argc should be 4.                                      //
//For example, to complete Problem 2, you make the following call://
// mpirun -np 6 randbin 10 25 1000                                //
//in order to simulate from the discrete distribution associated  //
//with the binomial coefficients for n = 10 in 25 batches of      //
//1000 samples each. You requested MPI to create 6 copies of your //
//program, hence you get 1 Master and 5 Slaves                    //
////////////////////////////////////////////////////////////////////
int main(int argc,char** argv)
{
        int i,j,k;
        double s;
	int myrank;

	if(4!=argc)
	{
	  fprintf(stderr,"USAGE: mpirun -np <NumberOfProcesses> randbin <n> <Number of batches> <Number of samples per batch>\n");
	  return(0);
	}

	///////////////////////////
	// START THE MPI SESSION //
	///////////////////////////
	MPI_Init(&argc, &argv);

	/////////////////////////////////////
	// What is the ID for the process? //   
	/////////////////////////////////////
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

	int n = atoi(argv[1]);
	int NumberOfBatches = atoi(argv[2]);
	int NumberOfIterations = atoi(argv[3]);

	if(myrank==0)
	{
	  fprintf(stderr,
		"You asked for n=%d, %d batches and %d iterations per batch.\n",
		n,NumberOfBatches,NumberOfIterations);
	  master(n);
	}
	else
	{
	  slave(n,NumberOfBatches,NumberOfIterations,myrank);
	}

	// Finalize the MPI session
	MPI_Finalize();

	return(1);
}
