#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

/* store bits as bytes */
typedef unsigned char bit_t;
/* pointer to individual = pointer to the first bit of it */
typedef bit_t *p_ind;
/* pointer to population = pointer to the first individual in it */
typedef p_ind *p_pop;

/*const int pop_size = 200;*/
#define pop_size 200
const int chromo_len = 100;
const double prob_cross = 0.9;
const double prob_mut = 0.5;
const int tournament_size = 3;
/*const int elitism_size = 2;*/
#define elitism_size 2

/* current population is stored in pop, new population is
   created in new_pop, pop and new_pop are then exchanged */
p_pop pop = NULL, new_pop = NULL;
/* fitness of each individual is stored here */
double rating[pop_size];
/* indices of first elitism_size best individuals */
int best[elitism_size+1];
/* generation counter */
int gen_count = 0;

/* useful wrappers around the fatal_error function */
#define fatal() fatal_error(__FILE__, __LINE__, NULL)
#define fatal2(msg) fatal_error(__FILE__, __LINE__, msg)

/* when called in case of a system error prints a meaningful error
   message and exits with return code 1 */
void fatal_error(char *filename, int linenum, char *error_msg)
{
  if (!error_msg)
    error_msg = strerror(errno);
  fprintf(stderr, "%s:%d: %s\n", filename, linenum, error_msg);
  exit(1);
}

/* returns a random integer in interval [0;max] */
int random_int(int max)
{
  return rand() % (max+1);
}

/* true if succeeded, false otherwise */
int random_trial(double prob)
{
  return ((double)rand()/RAND_MAX <= prob) ? 1 : 0;
}

p_pop allocate(p_pop pop)
{
  p_ind store = NULL;
  /* always allocate even number of individuals */
  int alloc_size = pop_size + ((pop_size-elitism_size) % 2);
  int i;
  if ((store = (p_ind)malloc(alloc_size*chromo_len*sizeof(bit_t))) == NULL)
    fatal();
  if ((pop = (p_pop)malloc(alloc_size*sizeof(p_ind))) == NULL)
    fatal();
  for (i = 0; i < pop_size; ++i)
    pop[i] = &store[i*chromo_len];
  return pop;
}

void deallocate(p_pop pop)
{
  free(pop[0]);
  free(pop);
}

void init()
{
  pop = allocate(pop);
  new_pop = allocate(new_pop);
  int i, j;
  /* all bits in all individuals are randomly generated */
  for (i = 0; i < pop_size; ++i)
    for (j = 0; j < chromo_len; ++j)
      pop[i][j] = random_int(1);
}

void output()
{
  printf("generation: %03d fitness: %f\n", gen_count, rating[best[0]]);
}

int job_done()
{
  return (gen_count >= 50) ? 1 : 0;
}

void update_best(int ind_idx)
{
  int i = elitism_size - 1;
  while (rating[ind_idx] > rating[best[i]] && i >= 0) {
    best[i+1] = best[i];
    --i;
  }
  best[i+1] = ind_idx;
}    

double fitness(p_ind ind)
{
  int i;
  int acc = 0;
  for (i = 0; i < chromo_len; ++i)
    acc += ind[i];
  return (double)acc;
}

void eval()
{
  int i, j;
  /* reset pointers to the best individuals in the population */
  for ( i = 0; i < elitism_size; ++i)
    best[i] = 0;
  /* evaluate fitness for all individuals and record
     elitism_size best of them */
  for (i = 0; i < pop_size; ++i) {
    rating[i] = fitness(pop[i]);
    update_best(i);
  }
}

/* selects an individual using the tournament selection scheme */
void select_ind(p_ind ind, p_pop pop)
{
  int best = random_int(pop_size-1);
  int candidate;
  int i;
  for (i = 0; i < tournament_size-1; ++i) {
    candidate = random_int(pop_size-1);
    if (rating[candidate] > rating[best])
      best = candidate;
  }
  for (i = 0; i < chromo_len; ++i)
    ind[i] = pop[best][i];
}

/* inserts individual ind into population pop at position pos*/
void insert(p_ind ind, p_pop pop, int position)
{

  int i;
  for (i = 0; i < chromo_len; ++i)
    pop[position][i] = ind[i];
}


/* the crossover operator, produces offspring in place of parents */
void cross(p_ind ind1, p_ind ind2)
{
  bit_t dummy[chromo_len];
  int point = random_int(chromo_len - 2);
  int i;
  for (i = 0; i < point; ++i)
    dummy[i] = ind1[i];
  for (i = 0; i < point; ++i)
    ind1[i] = ind2[i];
  for (i = 0; i < point; ++i)
    ind2[i] = dummy[i];
}

/* the mutation operator, modifies individual in place */
void mutate(p_ind ind)
{
  int bit = random_int(chromo_len-1);
  ind[bit] = ind[bit] ? 0 : 1;
}

/* generates new_pop from pop and exchanges pop and new_pop */
void renew_pop()
{
  /* elitism: copy elitism_size best individuals */
  int i;
  for (i = 0; i < elitism_size; ++i)
    insert(pop[best[i]], new_pop, i);

  /* generate the rest of population (pop_size - 2 individuals) */
  bit_t ind1[chromo_len], ind2[chromo_len];

  for (; i < pop_size; i += 2) {
    select_ind(ind1, pop);
    select_ind(ind2, pop);
    if (random_trial(prob_cross)) cross(ind1, ind2);
    if (random_trial(prob_mut)) mutate(ind1);
    if (random_trial(prob_mut)) mutate(ind2);
    insert(ind1, new_pop, i);
    insert(ind2, new_pop, i+1);
  }

  /* pop <-> new_pop */
  p_pop dummy;
  dummy = pop;
  pop = new_pop;
  new_pop = dummy;

  /* update generation counter */
  ++gen_count;
}

int main()
{
  init();
  eval();
  while (!job_done()) {
    renew_pop();
    eval();
    output();
  }
  deallocate(pop);
  deallocate(new_pop);
}