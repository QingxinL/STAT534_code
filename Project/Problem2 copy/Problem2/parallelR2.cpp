/*
 This program computes the R2 of several regressions in parallel.
 Compile the program using the makefile provided.
 
 Run the program using the command:

 mpirun -np 10 parallelR2 
*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <mpi.h>
#include <iomanip>

// For MPI communication
#define GETR2	1
#define DIETAG	0

// Used to determine MASTER or SLAVE
static int myrank;

// Global variables
int nobservations = 40;
int nvariables = 1000;
double* Y = NULL;
double** X = NULL;

double ssy = 0.0;	// used in R2 calculation

// Function Declarations
void NormStand();
void master();
void slave(int slavename);
double GetR2(int v);


int main(int argc, char* argv[])
{
   int i,j;
   FILE *yin, *xin;
   double tmp;

   ///////////////////////////
   // START THE MPI SESSION //
   ///////////////////////////
   MPI_Init(&argc, &argv);

   /////////////////////////////////////
   // What is the ID for the process? //   
   /////////////////////////////////////
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

   // Read in the data
   xin = fopen("X.txt", "r");
   yin = fopen("Y.txt", "r");
   if( (NULL == xin) || (NULL == yin))
   {
      printf("Cannot open data file!\n");
      MPI_Finalize();
      exit(1);
   }

   X = new double*[nobservations];
   Y = new double[nobservations];

   for(i=0; i<nobservations; i++)
   {
      X[i] = new double[nvariables];
      for(j=0; j<nvariables; j++)
      {
         fscanf(xin, "%lf", &tmp);
         X[i][j] = tmp;
      }
      fscanf(yin, "%lf", &tmp);
      Y[i] = tmp;
   }
   
   fclose(xin);
   fclose(yin);

   // Demean and standardize the data...
   NormStand();

   // Compute the ssy value
   // Used to calculate R2
   for(i=0; i<nobservations; i++) ssy += Y[i]*Y[i];

   // Branch off to master or slave function
   // Master has ID == 0, the slaves are then in order 1,2,3,...
   
   if(myrank==0)
   {
      master();
   }
   else
   {
      slave(myrank);
   }

   // clean memory
   for(i=0; i<nobservations; i++)
   {
      delete[] X[i]; X[i] = NULL;
   }
   delete[] Y; Y = NULL;
   delete[] X; X = NULL;

   // Finalize the MPI session
   MPI_Finalize();

   return(1);
}

void master()
{
   int var;		// to loop over the variables
   int rank;		// another looping variable
   int ntasks;		// the total number of slaves
   int jobsRunning;	// how many slaves we have working
   int work[1];		// information to send to the slaves
   double workresults[2]; // info received from the slaves
   FILE* fout;		// the output file
   MPI_Status status;	// MPI information

   fout = fopen("R2_values.txt","w");

   // Find out how many slaves there are
   MPI_Comm_size(MPI_COMM_WORLD,
		     &ntasks);

   fprintf(stdout, "Total Number of processors = %d\n",ntasks);

   // Now loop through the variables and compute the R2 values in
   // parallel
   jobsRunning = 1;

   for(var=0; var<nvariables; var++)
   {
      // This will tell the slave which variable to work on
      work[0] = var;

      if(jobsRunning < ntasks) // Do we have an available processor?
      {
         // Send out a work request
         MPI_Send(&work, 	// the vector with the variable
		  1, 		// the size of the vector
		  MPI_INT,	// the type of the vector
                  jobsRunning,	// the ID of the slave to use
                  GETR2,	// tells the slave what to do
                  MPI_COMM_WORLD); // send the request out to anyone
				   // who is available
         printf("Master sends out work request [%d] to slave [%d]\n",
                work[0],jobsRunning);

         // Increase the # of processors in use
         jobsRunning++;

      }
      else // all the processors are in use!
      {
         MPI_Recv(workresults,	// where to store the results
 		  2,		// the size of the vector
		  MPI_DOUBLE,	// the type of the vector
	 	  MPI_ANY_SOURCE,
		  MPI_ANY_TAG, 	
		  MPI_COMM_WORLD,
		  &status);     // lets us know which processor
				// returned these results

         printf("Master has received the result of work request [%d] from slave [%d]\n",
                (int) workresults[0],status.MPI_SOURCE);
 
         // Print out the results
         fprintf(fout, "%d\t%f\n", (int)workresults[0]+1, workresults[1]);

         printf("Master sends out work request [%d] to slave [%d]\n",
                work[0],status.MPI_SOURCE);

         // Send out a new work order to the processors that just
         // returned
         MPI_Send(&work,
                  1,
                  MPI_INT,
                  status.MPI_SOURCE, // the slave that just returned
                  GETR2,
                  MPI_COMM_WORLD); 
      } // using all the processors
   } // loop over all the variables


   ///////////////////////////////////////////////////////////////
   // NOTE: we still have some work requests out that need to be
   // collected. Collect those results now!
   ///////////////////////////////////////////////////////////////

   // loop over all the slaves
   for(rank=1; rank<jobsRunning; rank++)
   {
      MPI_Recv(workresults,
               2,
               MPI_DOUBLE,
               MPI_ANY_SOURCE,	// whoever is ready to report back
               MPI_ANY_TAG,
               MPI_COMM_WORLD,
               &status);

       printf("Master has received the result of work request [%d]\n",
                (int) workresults[0]);
 
      //save the results received
      fprintf(fout, "%d\t%f\n", (int)workresults[0]+1, workresults[1]);
   }

   printf("Tell the slave to die\n");

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

   printf("got to the end of Master code\n");

   fclose(fout);

   // return to the main function
   return;
  
}

void slave(int slavename)
{
   int work[1];			// the inputs from the master
   double workresults[2];	// the outputs for the master
   MPI_Status status;		// for MPI communication

   // the slave listens for instructions...
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
         case GETR2:
            // Get the R2 value for this variable
            // ...and save it in the results vector

           printf("Slave %d has received work request [%d]\n",
                  slavename,work[0]);
          
	    workresults[1] = GetR2(work[0]);

            // tell the master what variable you're returning

            workresults[0] = (double)work[0];

            // Send the results
            MPI_Send(&workresults,
                     2,
                     MPI_DOUBLE,
                     0,		// send it to the master
                     0,		// doesn't need a TAG
                     MPI_COMM_WORLD);

            printf("Slave %d finished processing work request [%d]\n",
                   slavename,work[0]);

            break;

         case DIETAG:
            printf("Slave %d was told to die\n",slavename);
            return;

         default:
            notDone = 0;
            printf("The slave code should never get here.\n");
            return;
      }
   }

   // No memory to clean up, so just return to the main function
   return;
}

// Data must have zero mean and unit variance
// This is only for 1 variable regressions without an intercept
double GetR2(int v)
{
   int	i;
   double tmp;

   tmp = 0.0;
   for(i=0; i<nobservations; i++)
   {
      tmp += (X[i][v]*Y[i]);
   }
   tmp = tmp*tmp / ((double)nobservations - 1.0);

   return(tmp / ssy);
}


void NormStand()
{
   int i, j;
   double tmp;

   for(j=0; j<nvariables; j++)
   {
      tmp = 0.0;
      for(i=0; i<nobservations; i++)
      {
         tmp += X[i][j];
      }
      for(i=0; i<nobservations; i++)
      {
         X[i][j] -= tmp/((double)nobservations);
      }
   }

   tmp = 0.0;
   for(i=0; i<nobservations; i++)
   {
      tmp += Y[i];
   }
   for(i=0; i<nobservations; i++) Y[i] -= tmp/((double)nobservations);

   // Now make the data have unit sample variance
   for(j=0; j<nvariables; j++)
   {
      tmp = 0.0;
      for(i=0; i<nobservations; i++)
      {
         tmp += X[i][j]*X[i][j];
      }

      tmp = sqrt(tmp / ((double)nobservations-1.0));

      for(i=0; i<nobservations; i++)
      {
         X[i][j] = X[i][j]/tmp;
      }
   }   

   // Do the same for Y
   tmp = 0.0;
   for(i=0; i<nobservations; i++)
   {
      tmp += Y[i]*Y[i];
   }
   tmp = sqrt( tmp / ((double)nobservations - 1.0));

   for(i=0; i<nobservations; i++)
   {
      Y[i] = Y[i]/tmp;
   }

   return;
}

