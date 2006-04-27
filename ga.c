#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>

/* store bits as bytes */
typedef bool bit_t;
/* pointer to individual = pointer to the first bit of it */
typedef bit_t *p_ind;
/* pointer to population = pointer to the first individual in it */
typedef p_ind *p_pop;

/*const int pop_size = 200;*/
#define pop_size 200
const size_t chromo_len = 100;
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
// Shuffle bit order?
bool shuffle = false;
// Mapping for bit number -> elemnt of chromosome.
size_t * bit_order = NULL;
// Pointer to chosen fitness function.
double (* fitness) (p_ind);

/* useful wrappers around the fatal_error function */
#define fatal() fatal_error(__FILE__, __LINE__, NULL)
#define fatal2(msg) fatal_error(__FILE__, __LINE__, msg)

/* when called in case of a system error prints a meaningful error
   message and exits with return code 1 */
void
fatal_error(char const *filename, int linenum, char const *error_msg)
{
  if (!error_msg)
    error_msg = strerror(errno);
  fprintf(stderr, "%s:%d: %s\n", filename, linenum, error_msg);
  exit(EXIT_FAILURE);
}


inline
bit_t
get_bit (p_ind const ind, size_t i)
{
  assert (i < (size_t)chromo_len);
  return ind[bit_order[i]];
}

inline
void
set_bit (p_ind ind, size_t i, bool val)
{
  assert (i < (size_t)chromo_len);
  ind[bit_order[i]] = val;
}


/* returns a random integer in interval [0;max] */
inline
int
random_int (int max)
{
  return nearbyint (((double)rand () / RAND_MAX) * max);
}

/* true if succeeded, false otherwise */
inline
bool
random_trial (double prob)
{
  return (double)rand()/RAND_MAX <= prob ? true : false;
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

void init(void)
{
  pop = allocate(pop);
  new_pop = allocate(new_pop);
  size_t i, j;
  /* all bits in all individuals are randomly generated */
  for (i = 0; i < pop_size; ++i)
    for (j = 0; j < chromo_len; ++j)
      pop[i][j] = random_int(1);
  // Shuffle bit order array.
  if ((bit_order = (size_t *)malloc (chromo_len * sizeof (size_t))) == NULL)
    fatal ();
  for (i = 0; i < chromo_len; ++i)
    {
      bit_order[i] = i;
    }
  if (shuffle)
    for (i = 0; i < chromo_len; ++i)
      {
        i = random_int (chromo_len);
        j = random_int (chromo_len);
        size_t tmp = bit_order[i];
        bit_order[i] = bit_order[j];
        bit_order[j] = tmp;
      }
}

void output(void)
{
  printf("generation: %03d fitness: %f\n", gen_count, rating[best[0]]);
}

bool job_done(void)
{
  return (gen_count >= 50) ? true : false;
}

void update_best(int ind_idx)
{
  int i = elitism_size - 1;
  while (rating[ind_idx] > rating[best[i]] && i >= 0)
    {
      best[i+1] = best[i];
      --i;
    }
  best[i+1] = ind_idx;
}


// OneMax function.
double 
fitness_onemax (p_ind ind)
{
  unsigned acc = 0;
  for (size_t i = 0; i < chromo_len; ++i)
    acc += get_bit (ind, i);
  return acc;
}


// Rosenbrock function.
double
fitness_rosenbrock (p_ind ind)
{
  int x1 = 0;
  int x2 = 0;
  double sum = 0;
  size_t const count = chromo_len / 12; // + ((chromo_len % 12) ? 1 : 0);
  for (size_t i = 0; i < count; ++i)
    {
      x2 = 0;
      size_t const base = i * 12;
      for (size_t j = 0; j <= 12 && base + j < chromo_len; ++j)
        x2 = x2 * 2 + get_bit (ind, base + j);
      // Assuming 2's complement code, convert bit pattern of negative number
      // into the number.
      if (x2 > 2047)
        x2 = -(4096 - x2);
      if (i >= 1)
        sum += 100 * pow (pow (x1, 2) - x2, 2) + pow (1 - x1, 2);
      x1 = x2;
    }
  // Minimalization, the smaller the sum is the bigger is resulting fitness.
  return -sum;
}

// F101 function.
double
fitness_f101 (p_ind ind)
{
  int x1 = 0;
  int x2 = 0;
  double sum = 0;
  //size_t const count = chromo_len / 10+ ((chromo_len % 10) ? 1 : 0);
  for (size_t i = 0; i <= 10; ++i)
    {
      x2 = 0;
      size_t const base = (i * 10) % chromo_len;;
      for (size_t j = 0; j <= 10 && base + j < chromo_len; ++j)
        x2 = x2 * 2 + get_bit (ind, base + j);
      // Assuming 2's complement code, convert bit pattern of negative number
      // into the number.
      if (x2 > 511)
        x2 = -(1024 - x2);
      if (i >= 1)
        sum += 
          -x1 * sin (sqrt (abs (x1 - x2 - 47.0)))
          -(x2 + 47.0) * sin (sqrt (abs (x2 + 47 + x1 / 2.0)));
      x1 = x2;
    }
  // Minimalization, thus the minus.
  return -sum;
}


void eval(void)
{
  int i;
  /* reset pointers to the best individuals in the population */
  for (i = 0; i < elitism_size; ++i)
    best[i] = 0;
  /* evaluate fitness for all individuals and record
     elitism_size best of them */
  for (i = 0; i < pop_size; ++i)
    {
      rating[i] = fitness(pop[i]);
      update_best(i);
    }
}

/* selects an individual using the tournament selection scheme */
void select_ind(p_ind ind, p_pop pop)
{
  int best = random_int(pop_size-1);
  int candidate;
  size_t i;
  for (i = 0; i < tournament_size-1; ++i)
    {
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
  for (size_t i = 0; i < chromo_len; ++i)
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
extern void mutate(p_ind ind)
{
  int bit = random_int(chromo_len-1);
  set_bit (ind, bit, get_bit (ind, bit) ? false : true);
}

/* generates new_pop from pop and exchanges pop and new_pop */
void renew_pop(void)
{
  /* elitism: copy elitism_size best individuals */
  int i;
  for (i = 0; i < elitism_size; ++i)
    insert(pop[best[i]], new_pop, i);

  /* generate the rest of population (pop_size - 2 individuals) */
  bit_t ind1[chromo_len], ind2[chromo_len];

  for (; i < pop_size; i += 2)
    {
      select_ind(ind1, pop);
      select_ind(ind2, pop);
      if (random_trial(prob_cross))
        cross(ind1, ind2);
      if (random_trial(prob_mut))
        mutate(ind1);
      if (random_trial(prob_mut))
        mutate(ind2);
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


void
usage (void)
{
  fprintf (stderr, 
           "Usage: ga -sffoh?\n"
           "\t-s\tshuffle chromosome\n"
           "\t-f\tF101 fitness function\n"
           "\t-r\tRosenbrock fitness function\n"
           "\t-o\tOneMax fitness function\n"
           "\t-h\tthis help\n"
           "\t-?\tthis help\n");
  exit (EXIT_FAILURE);
}


void
analyze_options (int argc, char * const argv[])
{
  static char const opts[] = "sfroh?";

  int ch;
  while ((ch = getopt(argc, argv, opts)) != -1)
    {
      switch (ch)
        {
        case 's':
          shuffle = true;
          break;
        case 'f':
          fitness = fitness_f101;
          break;
        case 'r':
          fitness = fitness_rosenbrock;
          break;
        case 'o':
          fitness = fitness_onemax;
          break;
        case 'h':
        case '?':
        default:
          usage ();
        }
    }
}


int
main (int argc, char * const argv[])
{
  analyze_options (argc, argv);
  init();
  eval();
  while (!job_done())
    {
      renew_pop();
      eval();
      output();
    }
  deallocate(pop);
  deallocate(new_pop);

  exit (EXIT_SUCCESS);
}
