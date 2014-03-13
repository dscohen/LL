#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>
#include <math.h>

/***** Globals ******/
float **a; /* The coefficients */
float *x;  /* The unknowns */
float *b;  /* The constants */
float *curr;
float err = 10; /* The absolute relative error */
int num = 0;  /* number of unknowns */


#define COLUMNS num/comm_sz;

/****** Function declarations */
void check_matrix(); /* Check whether the matrix will converge */
void get_input();  /* Read input from file */

/********************************/



/* Function definitions: functions are ordered alphabetically ****/
/*****************************************************************/

/* 
   Conditions for convergence (diagonal dominance):
   1. diagonal element >= sum of all other elements of the row
   2. At least one diagonal element > sum of all other elements of the row
   */
void check_matrix()
{
  int bigger = 0; /* Set to 1 if at least one diag element > sum  */
  int i, j;
  float sum = 0;
  float aii = 0;

  for(i = 0; i < num; i++)
  {
    sum = 0;
    aii = fabs(a[i][i]);

    for(j = 0; j < num; j++)
      if( j != i)
        sum += fabs(a[i][j]);

    if( aii < sum)
    {
      printf("The matrix will not converge\n");
      exit(1);
    }

    if(aii > sum)
      bigger++;

  }

  if( !bigger )
  {
    printf("The matrix will not converge\n");
    exit(1);
  }
}


/******************************************************/
/* Read input from file */
void get_input(char filename[])
{
  FILE * fp;
  int i,j;  

  fp = fopen(filename, "r");
  if(!fp)
  {
    printf("Cannot open file %s\n", filename);
    exit(1);
  }

  fscanf(fp,"%d ",&num);
  fscanf(fp,"%f ",&err);

  /* Now, time to allocate the matrices and vectors */
  a = (float**)malloc(num * sizeof(float*));
  if( !a)
  {
    printf("Cannot allocate a!\n");
    exit(1);
  }

  for(i = 0; i < num; i++) 
  {
    a[i] = (float *)malloc(num * sizeof(float)); 
    if( !a[i])
    {
      printf("Cannot allocate a[%d]!\n",i);
      exit(1);
    }
  }

  x = (float *) malloc(num * sizeof(float));
  if( !x)
  {
    printf("Cannot allocate x!\n");
    exit(1);
  }

  curr = (float *) malloc(num * sizeof(float));
  if( !curr)
  {
    printf("Cannot allocate curr!\n");
    exit(1);
  }

  b = (float *) malloc(num * sizeof(float));
  if( !b)
  {
    printf("Cannot allocate b!\n");
    exit(1);
  }

  /* Now .. Filling the blanks */ 

  /* The initial values of Xs */
  for(i = 0; i < num; i++)
    fscanf(fp,"%f ", &x[i]);

  for(i = 0; i < num; i++)
  {
    for(j = 0; j < num; j++)
      fscanf(fp,"%f ",&a[i][j]);

    /* reading the b element */
    fscanf(fp,"%f ",&b[i]);
  }

  fclose(fp); 

}


/************************************************************/
/***** Globals ******/
/*float **a; [> The coefficients <]*/
/*float *x;  [> The unknowns <]*/
/*float *b;  [> The constants <]*/
/*float err; [> The absolute relative error <]*/
/*int num = 0;  [> number of unknowns <]*/


int main(int argc, char *argv[])
{
  int comm_sz;
  int my_rank;
  int root = 0;
  int irow, jrow, icol, index, real_row;
  int current_error = 1e5;

  int i;
  int nit = 0; /* number of iterations */
  int left_over = 0;
  int *displ, *sendcount, *Xsendcount;
  float *local_x, *local_A, *new_a;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  //restrict to each


  if(my_rank == root) {
    if( argc != 2)
    {
      printf("Usage: gsref filename\n");
      exit(1);
    }

    /* Read the input file and fill the global data structure above */ 
    get_input(argv[1]);


    /* Check for convergence condition */
    check_matrix();
    new_a = (float*)malloc(num*num*sizeof(float));
    int i, j = 0;
    for (i = 0; i < num; i++) {
      for (j = 0; j < num; j++) {
        new_a[i*num + j] = a[i][j];
      }
    }
  }


  //Send out updated globals
  //creating pounters
  int *n = &num;
  float *error = &err;
  //THIS PART IS FOR WHEN ALOT OF PROCESSORS ARE AVAILABLE
  /*if (comm_sz > num) {*/
    /*comm_sz = num;}*/

  MPI_Bcast(n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (my_rank != 0) {
    x = (float*) malloc(num * sizeof(float));
    b = (float*) malloc(num * sizeof(float));
  }
  local_x = (float*) malloc(num * sizeof(float));
  MPI_Bcast(error, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(x, num, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(b, num, MPI_FLOAT, 0, MPI_COMM_WORLD);

  //getting leftover rows for last process and making room for local partition of A
    left_over = num%comm_sz;
  if (my_rank == comm_sz-1) {
    local_A = (float*) malloc((num*num/comm_sz+num*left_over) * sizeof(float));
  } else {
    local_A = (float*) malloc((num*num/comm_sz) * sizeof(float));
  }

  //Make send counts and parameter for scatterv and gather
  sendcount = (int*)  calloc(comm_sz,sizeof(int));
  Xsendcount = (int*) malloc(comm_sz*sizeof(int));
  displ = (int*)  calloc(comm_sz,sizeof(int));
  int *Xdispl = (int*)  calloc(comm_sz, sizeof(int));
  for (i = 0; i < comm_sz; i++) {
    //how many data points each process is going to get
    Xdispl[i] = (num/comm_sz)*i;
    int displacement = i*num;
    displ[i] += displacement;
    sendcount[i] = num*(num/comm_sz);
    Xsendcount[i] = num/comm_sz;
  }
    sendcount[comm_sz-1] += num*left_over;
    Xsendcount[comm_sz-1] += left_over;
  //Make stuff for gatherV

  //sending out matrix A partitions
  MPI_Scatterv(new_a, sendcount, displ, MPI_FLOAT, local_A, num*(num/comm_sz)+num*left_over, MPI_FLOAT, 0, MPI_COMM_WORLD);


  float *pass_x;

  float temp_sum;
  int plus_extra;

  /*while (current_error > err) {*/
  int iterationsss = 0;
  if (my_rank == comm_sz-1) {
    plus_extra = left_over;
  } else {
    plus_extra = 0;
  }
  pass_x = malloc((num/comm_sz + plus_extra) * sizeof(float));

  while(iterationsss < 20){
    iterationsss++;
    //begin work loop for jacobi
    //
    int local_rows = num/comm_sz;
    int j;
    for(i = 0; i < local_rows+plus_extra; i++) {
    temp_sum = 0;
      //figure out real row in overall matrix
      //initialize sum of row at 0 for each new row
      temp_sum = 0;
      real_row = (num/comm_sz)*my_rank+i;
      
      for (j = 0; j < num; j++) {
        //if not aij (x to be solved)
        if (j != my_rank+i){
          temp_sum = local_A[i*num+j]*x[j] + temp_sum;
        }
      }

      //solve for x of that row for next iteration
      //see if I should change local or real x
      local_x[real_row] = (b[real_row] - temp_sum)/local_A[i*num + my_rank+i];
      x[real_row] = local_x[real_row];
      pass_x[i] = local_x[real_row];

      //calculate error for x, take max of error for all for's
    }

    MPI_Allgatherv(pass_x, num/comm_sz+plus_extra, MPI_FLOAT,x, Xsendcount, Xdispl, MPI_FLOAT, MPI_COMM_WORLD);
  }





  /* Writing to the stdout */
  /* Keep that same format */

  if (my_rank == 0){
    for( i = 0; i < num; i++)
      printf("%f\n",x[i]);
    printf("%d ::: total number of iterations: %d\n", my_rank, nit);
  }
  //FREE EVERYTHING
  MPI_Finalize();


  return 0;

}
