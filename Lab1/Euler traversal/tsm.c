#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

int logical_n = 1;
int n;
int test;
int **a;

void get_matrix(char filename[])
{
  int i,j;
  char temp[60];
  FILE * fptr;

  fptr = fopen(filename, "r");
  if(!fptr) {
    printf("Cannot open file %s\n", filename);
    exit(1);
  }
  while (fgets(temp, sizeof(temp), fptr)!=NULL){
    logical_n++;
  }
  rewind(fptr);

  a = (int**)malloc(logical_n*sizeof(int*));
  for(i = 0; i < logical_n; i++) 
  {
    a[i] = (int *)malloc(logical_n * sizeof(int)); 
    if( !a[i])
    {
      printf("Cannot allocate a[%d]!\n",i);
      exit(1);
    }
  }
  n = logical_n-1;
  for(i = 0; i < n; i++)
  {
    for(j = 0; j < n; j++){
      fscanf(fptr,"%d ", &a[i][j]);
    }
    /*fscanf(fptr,"%d ",&a[i][j]);*/
  }

  fclose(fptr); 
  for(j = 0; j < logical_n; j++) {
    a[n][j] = 0xFFFFFFF;
  }
  for(j = 0; j < logical_n; j++) {
    a[j][n] = 0xFFFFFFF;
  }
}

#define START       200 //starting population
#define TRUE        1 //starting population
#define FALSE       0 //starting population

int compute_length(int* g){
  int i = 0;
  int length = 0;
  for (i = 0; i < n-1; i++) {
    int depart = g[i];
    int arrive = g[i+1];
    length += a[depart][arrive];
  }
  return length;
}

void crossover(int *b1, int *b, int *c) {
  int *used = (int*)calloc(logical_n, sizeof(int));
  used[n] = -1;
  int f_a, present = FALSE;
  int valid, f_b = TRUE;
  int x = 0; int i = 0; int z = 0;
  int t1 = b1[x];
  int t2 = b1[x];
  int alt_t2, y, temp;
  //Find location of random town in each array
  while(z < n) {
    alt_t2 = 0;
    y = 0;
    //check to see if next city in b1 already exists
    while(TRUE){
      t1 = b1[x];
      if (used[t1] < 0) {
        /*printf("%d already exists,", t1);*/
        //if already n c, increase pointer in a.
        x = (x+1) % n;
        if (b1[x] == t1){
          /*printf("present bit has been set");*/
          present = TRUE;
          break;
        }
      } else {
        //town 1 isn't in c yet, then use that
        break;
      }
    }
    if (present)
      break;
    if(used[t1] < 0) {
      printf("adding repeat");
    }
    c[z++] = t1;
    used[t1] += -1;
    if(z >= n)
      break;
    //find equivalent city in 2nd chromosome
    y = rand() % n;
    //if not found (last city on tour, take random);
    for (i = 0; i < n-1; i++) {
      if (b[i] == t1){
        y = i;
        break;
      }
    }
    //find the path that has the shortest route
    if ((a[b1[x]][b1[x+1]] >= a[b[y]][b[y+1]])){
      t2 = b[y+1];
      //alternate city if t2 original exists already in new c
      alt_t2 = b1[x+1];
    } 
    //take t2 from b1 instead of b
    if ((a[b1[x]][b1[x+1]] < a[b[y]][b[y+1]]) && (b[y+1] != n)){
      t2 = b1[x+1]; alt_t2 = b[y+1];
      if (b[y+1] == n) {
        //would mean infinite distance,
        //again, biased but is ok for this purpose
        alt_t2 = rand() % n;
      }
    }
    //check to see if it already exists in new c
    if(used[t2] < 0){
      if(used[alt_t2] < 0){
        while(used[alt_t2] < 0) {
          alt_t2 = rand() % n;
        }
      }
      if (used[alt_t2] < 0)
        printf("adding repeat\n");
      c[z++] = alt_t2;
      used[alt_t2] += -1;
    } else {
      if (used[t2] < 0)
        printf("adding repeat\n");
      c[z++] = t2;
      used[t2] += -1;
    }
    //not present, so make alt town the next town
    x++;
  } 
  c[n] = n;
  free(used);
}


/*void greedy_swap(int *a, int *b) {*/
  /*int original_d = compute_length(a);*/
  /*int t1 = rand() % (n-1);*/
  /*int t2 = t1+1;*/
  /*int temp,i;*/
  /*for (i = 0; i < logical_n; i++) {*/
    /*b[i] = a[i];*/
  /*}*/
  /*for (i = 0; i < 3; i++) {*/
    /*temp = a[t1];*/
    /*a[t1] = a[t2];*/
    /*a[t2] = temp;*/
    /*if(compute_length(a) < original_d) {*/
      /*b[t1] = a[t1];*/
      /*b[t2] = a[t2];*/
      /*break;*/
    /*} else {*/
      /*b[t2] = b[t1];*/
      /*b[t1] = temp;*/
      /*t1 = rand() % (n-1);*/
      /*t2 = t1+1;*/
    /*}*/
  /*}*/
  /*int k;*/

/*}*/

int **gene;
int **child_gene;
int main(int argc, char *argv[])
{
  int i,j,r,hold, thread_count, original_d, t1, t2;
  int **elite;
  int not_converged = TRUE;
  /*thread_count = strtol(argv[1], NULL, 10);*/
  get_matrix(argv[1]);
  //create population 1 and 2
  elite = (int**)malloc(10 * sizeof(int*));
  for (i = 0; i < 10; i++) {
    elite[i] = (int*)malloc(logical_n*sizeof(int*));
  }
  gene = (int**)malloc(START * sizeof(int*));
  child_gene = (int**)malloc(START * sizeof(int*));
  for (i = 0; i < START; i++) {
    gene[i] = (int*)malloc(logical_n*sizeof(int*));
    child_gene[i] = (int*)malloc(START * sizeof(int*));
    for (j = 0; j < logical_n; j++) {
      gene[i][j] = j;
      child_gene[i][j] = j;
    }
    //This method introduces a slight bias, but still performs within error bounds
    //populates pop1 with randoms
    for (j = 0; j < n; j++) {
      r = rand() % n;
      hold = gene[i][r];
      gene[i][r] = gene[i][j];
      gene[i][j] = hold;
      child_gene[i][j] = hold;
      child_gene[i][r] = gene[i][r];
    }

  }
  //Begin genetic algorithm action
  int k,l;
  int cycle = 0;
  int jumping = -1;
  int previous = -1;
  int iter = 0;
  while(TRUE){
    int e_index = 1;
    cycle %= 2;
    //find parents via tournament
    int local_parent1, local_parent2;
    int max1, max2, temp, best_max = 100000;
    //gene is parent when cycle is even
    for (j = 20; j <= START; j+=20) {
      max1 = 100000; max2 = 100000;
      if(!cycle) {
        for (i = j-20; i < j; i++) {
          temp = compute_length(gene[i]);
          if (temp < max1) {
            local_parent1 = i;
            max1 = temp;
            if (temp < best_max) {
              e_index = 0;
              best_max = temp;
              for (k = 0; k <logical_n; k++) {
                elite[e_index][k] = gene[i][k];
              }
            }
            continue;
          }
          if (temp < max2) {
            local_parent2 = i;
            max2 = temp;
            continue;
          }
          int original_d = compute_length(gene[i]);
          int t1 = rand() % (n-1);
          int t2 = t1+1;
          int temp,i;
          for (k = 0; k < 2; k++) {

            for (l = 0; l < logical_n; l++) {
              child_gene[i] = gene[i];
            }
            temp = gene[i][t1];
            child_gene[i][t1] = gene[i][t2];
            child_gene[i][t2] = temp;
          }
        }
      } else {
        //child_gene is parent on odd cycles;
        for (i = j-20; i < j; i++) {
          temp = compute_length(child_gene[i]);
          if (temp < max1) {
            local_parent1 = i;
            max1 = temp;
            if (temp < best_max) {
              e_index = 0;
              best_max = temp;
              for (k = 0; k <logical_n; k++) {
                elite[e_index][k] = child_gene[i][k];
              }
            }
            continue;
          }
          if (temp < max2) {
            local_parent2 = i;
            max2 = temp;
            continue;
          }
          original_d = compute_length(gene[i]);
          t1 = rand() % (n-1);
          t2 = t1+1;
          for (k = 0; k < 2; k++) {

            for (l = 0; l < logical_n; l++) {
              gene[i] = child_gene[i];
            }
            temp = child_gene[i][t1];
            gene[i][t1] = child_gene[i][t2];
            gene[i][t2] = temp;
          }
        }
      }
    if (!cycle) {
      crossover(gene[local_parent1], gene[local_parent2], child_gene[local_parent2]);
    } else {
      crossover(child_gene[local_parent1], child_gene[local_parent2], gene[local_parent2]);
    }
    }
    /////with two champs found, create and breed child into next generation at random spot.
    e_index++;
    if (((abs(previous - (temp = compute_length(elite[0]))) < 1) || (jumping == previous)) && (iter > 3)) {
      printf("final solution = %d\n", temp);
      for (i = 0; i < n; i++) {
        printf("%d ", elite[0][i]);
      }
      printf("\n");
      break;
    } else {
      jumping = previous;
      previous = temp;
      printf("best result for iter %d = %d, check with %d\n", iter, temp, best_max);
    }
    if(!cycle){
      for (i = 0; i < logical_n; i++) {
        child_gene[0][i] = elite[0][i];
      }
    } else {
      for (i = 0; i < logical_n; i++) {
        gene[0][i] = elite[0][i];
      }
    }
    cycle++;
    iter++;
  }
free(gene);
free(child_gene);
free(a);
}

