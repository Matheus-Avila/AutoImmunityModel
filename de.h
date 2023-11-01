/*------Constants for rnd_uni()--------------------------------------------*/
#include <stdbool.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0 / IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1 / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

typedef double evalfn(double *, int );

int get_de_strategy(char *strat);
double *diffential_evolution(bool from_file, evalfn *evaluate, double *inibound_l, double *inibound_h, int genmax, double CR, double F, int NP, int D, int strategy, long seed, int num_days);

