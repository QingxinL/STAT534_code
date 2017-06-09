/*
 This program use MCMC to sample the the P_i in parallel.
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
int nvariables = 10;
int numberOfGen = 25000;

// Function Declarations
void master();
void slave(int slavename);

double getPi(int n, int i);
int randInt(int a, int b);
double rand01();
double MCMC(int n, int i, int iter);

int main(int argc, char* argv[])
{
    // set the random seeds
    srand(unsigned(time(NULL)));
    
    ///////////////////////////
    // START THE MPI SESSION //
    ///////////////////////////
    MPI_Init(&argc, &argv);
    
    /////////////////////////////////////
    // What is the ID for the process? //
    /////////////////////////////////////
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    
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
    int work[2];		// information to send to the slaves  work[0] =i; work[1] = 25000;
    double workresults[2]; // info received from the slaves  workresults[0]=i; workresults[1] = p_e;
    FILE* fout;		// the output file
    MPI_Status status;	// MPI information
    double result[nvariables+1];    // save the results of [0, 1, 2,..., n]
    
    for (int i=0; i<=nvariables; i++)
        result[i] = 0.0;
    
    fout = fopen("output.txt","w");
    
    // Find out how many slaves there are
    MPI_Comm_size(MPI_COMM_WORLD,
                  &ntasks);
    
    fprintf(stdout, "Total Number of processors = %d\n",ntasks);
    
    // Now loop through the variables and compute the R2 values in
    // parallel
    jobsRunning = 1;
    work[1] = numberOfGen;
    
    for(var=0; var<=nvariables; var++)
    {
        // This will tell the slave which variable to work on
        work[0] = var;
        
        if(jobsRunning < ntasks) // Do we have an available processor?
        {
            // Send out a work request
            MPI_Send(&work, 	// the vector with the variable
                     2, 		// the size of the vector
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
            //fprintf(fout, "%d\t%f\n", (int)workresults[0]+1, workresults[1]);
            // save the result
            result[(int)workresults[0]] = workresults[1];
            
            printf("Master sends out work request [%d] to slave [%d]\n",
                   work[0],status.MPI_SOURCE);
            
            // Send out a new work order to the processors that just
            // returned
            MPI_Send(&work,
                     2,
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
        result[(int)workresults[0]] = workresults[1];
        //fprintf(fout, "%d\t%f\n", (int)workresults[0]+1, workresults[1]);
    }
    for (int i=0; i<=nvariables; i++)
        fprintf(fout, "%d\t%f\n", i, result[i]);
    
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
    int work[2];			// the inputs from the master
    double workresults[2];	// the outputs for the master
    MPI_Status status;		// for MPI communication
    
    // the slave listens for instructions...
    int notDone = 1;
    while(notDone)
    {
        printf("Slave %d is waiting\n",slavename);
        MPI_Recv(&work,		// the inputs from the master
                 2,		// the size of the inputs
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
                
                workresults[1] = MCMC(nvariables, work[0], work[1]);
                
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
                
            default:             notDone = 0;
                printf("The slave code should never get here.\n");
                return;
        }
    }
    
    // No memory to clean up, so just return to the main function
    return;
}

// This function return the Pi of n, i
double getPi(int n, int i)
{
    double Pi;
    Pi = (tgamma(n+1)/(tgamma(i+1)*tgamma(n-i+1)));
    return (Pi);
}

// return the random int belong [a, b]
int randInt(int a, int b)
{
    int out;
    out = (rand() % (b-a+1))+ a;
    return out;
}

// return the random number between [0,1]
double rand01()
{
    double out = rand()/double(RAND_MAX);
    return (out);
}

//Using the Metropolis-Hastings algorithm to sample the distribution
// iter=25000;
double MCMC(int n, int i, int iter)
{
    int sumI = 0;   //count the number that i_t = i
    double p_is, p_it;  //the prob
    int i_t;
    i_t = randInt(0, n); // for i_t(0)
    int i_t1 = 0;   // the i(t+1)
    int i_s = 0;    // the new random i*
    
    if (i_t==i) sumI++;  // consider the i_t(0)
    
    // generate iter samples
    for (int j=0; j<iter; j++)
    {
        i_s = randInt(0, n);
        p_is = getPi(n, i_s);   // get the P_i*
        p_it = getPi(n, i_t);   // get the P_it
        
        // MCMC set i(t+1) = i* with probability min(p_i*/p_it, 1)
        if (p_is>=p_it)
            i_t1 = i_s;
        else if (rand01() < (p_is/p_it))
            i_t1 = i_s;
        else
            i_t1 = i_t;
        
        // update the i_t
        i_t = i_t1;
        
        // count the sample i_t equal with i
        if (i_t==i) sumI++;
        
    }
    
    // return the estimated p_est
    double p_est = double(sumI)/iter;
    return (p_est);
}

