#include "ag.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


static double rnd_uni(long *idum) {

    long j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    double temp;

    if(*idum <= 0) {
        if(-(*idum) < 1)
            *idum = 1;
        else
            *idum = -(*idum);
        idum2 = (*idum);
        for(j = NTAB + 7; j >= 0; j--) {
            k = (*idum) / IQ1;
            *idum = IA1 * (*idum - k * IQ1) - k * IR1;
            if(*idum < 0)
                *idum += IM1;
            if(j < NTAB)
                iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ1;
    *idum = IA1 * (*idum - k * IQ1) - k * IR1;
    if(*idum < 0)
        *idum += IM1;
    k = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
    if(idum2 < 0)
        idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if(iy < 1)
        iy += IMM1;
    if((temp = AM * iy) > RNMX)
        return RNMX;
    else
        return temp;
}

static struct individual * select_parents_by_rank(struct individual *pop, int NP, int D, long *seed) {

	double rand;
	int ind_pop;
	struct individual *parents = malloc(sizeof(struct individual)*2);
	parents[0].chromossome = (double*) malloc(sizeof(double)*D);
	parents[1].chromossome = (double*) malloc(sizeof(double)*D);

	ind_pop = (int)(NP * (1.5 - sqrt(2.25 - 2 * rnd_uni(seed))));

	if (ind_pop < 0) {
		ind_pop = 0;             
	}

	parents[0].fitness = pop[ind_pop].fitness;
	memcpy(parents[0].chromossome, pop[ind_pop].chromossome, sizeof(double)*D);

	ind_pop = (int)(NP * (1.5 - sqrt(2.25 - 2 * rnd_uni(seed))));

	if (ind_pop < 0) {
		ind_pop = 0;             
	}

	parents[1].fitness = pop[ind_pop].fitness;
	memcpy(parents[1].chromossome, pop[ind_pop].chromossome, sizeof(double)*D);

	return parents;

}


struct individual * blend_crossover(struct individual x, struct individual y, int D,  double prob_crossover, long *seed, double *lb, double *hb) {

	struct individual *offspring = malloc(sizeof(struct individual)*2);
	offspring[0].chromossome = (double*) calloc(sizeof(double), D);
	offspring[1].chromossome = (double*) calloc(sizeof(double), D);

	double prob = rnd_uni(seed);

	double a = 0.6, b = 0.7, di;
	double u;

	if (prob < prob_crossover) {

		for(int i = 0; i < D; i++) {

			di = fabs(x.chromossome[i] - y.chromossome[i]);

			if (x.chromossome[i] <= y.chromossome[i]) {

				double low = x.chromossome[i] - a*di;
				double high = y.chromossome[i] + b*di;

				if(low < lb[i]) low = lb[i];
				if(high > hb[i]) high = hb[i];

				u = low + rnd_uni(seed) * (high - low);
				offspring[0].chromossome[i] = u;

				u = low + rnd_uni(seed) * (high - low);
				offspring[1].chromossome[i] = u;
			}
			else {
				
				double low = y.chromossome[i] - b*di;
				double high = x.chromossome[i] + a*di;

				if(low < lb[i]) low = lb[i];
				if(high > hb[i]) high = hb[i];

				u = low + rnd_uni(seed) * (high - low);
				offspring[0].chromossome[i] = u;

				u = low + rnd_uni(seed) * (high - low);
				offspring[1].chromossome[i] = u;
			}
		}

		return offspring;
	}
	else {
		offspring[0].fitness = x.fitness;
		memcpy(offspring[0].chromossome, x.chromossome, sizeof(double)*D);

		offspring[1].fitness = y.fitness;
		memcpy(offspring[1].chromossome, y.chromossome, sizeof(double)*D);
		return offspring;
	}
}

int comp_ind(const void *vind1, const void* vind2) {
	double i1 = ((struct individual *)vind1)->fitness;
	double i2 = ((struct individual *)vind2)->fitness;
	
	if(i1 > i2) return 1;
	if(i1 < i2) return -1;
	
	return 0;

}

double delta(double t, double y, int T, long *seed) {

	double b = 5.0;

	return y*(1.0 - pow(rnd_uni(seed),pow(1.0-(t/T),b)));

}

void mutation(struct individual ind, int D, double MR, double *LB, double *UB, int t, int T, long *seed) {

	double ran;
	int i, d;
	int l = D;	

	for(i = 0; i < l; i++){
		double r = rnd_uni(seed);
		if(r < MR) {
			double p = rnd_uni(seed);
			if(p < 0.5) {
				ind.chromossome[i] = ind.chromossome[i] + delta(t, UB[i]-ind.chromossome[i], T, seed);
			}
			else {
				ind.chromossome[i] = ind.chromossome[i] - delta(t, ind.chromossome[i]-LB[i], T, seed);
			}
		}
	}
}	



struct individual *ag(evalfn *evaluate, double *inibound_l, double *inibound_h, int genmax, double CR, double MR, int NP, int D, long seed) {

    /*-----Initialize random number generator-----------------------------*/

    long rnd_uni_init = -(long)seed; /* initialization of rnd_uni() */

    struct individual best;
   	best.chromossome = malloc(sizeof(double)*D);

    struct individual *pop = (struct individual*)malloc(sizeof(struct individual)*NP);
    struct individual *new_pop = (struct individual*)malloc(sizeof(struct individual)*NP);

	for(int i = 0; i < NP; i++) {
		pop[i].chromossome = (double*)malloc(sizeof(double)*D);
    }

    for(int i = 0; i < NP; i++) {
        for(int j = 0; j < D; j++) {
            pop[i].chromossome[j] = inibound_l[j] + rnd_uni(&rnd_uni_init) * (inibound_h[j] - inibound_l[j]);
        }

        pop[i].fitness = evaluate(pop[i].chromossome); /* obj. funct. value */
    }

    /*=======================================================================*/
    /*=========Iteration loop================================================*/
    /*=======================================================================*/

    int gen = 0;

	while(gen < genmax) {

		qsort(pop, NP, sizeof(struct individual), comp_ind);
		best.fitness = pop[0].fitness;
		memcpy(best.chromossome, pop[0].chromossome, sizeof(double)*D);

		for(int i = 0; i < NP; i+=2) {
    		struct individual *offspring;

			struct individual *parents;
		   
			parents = select_parents_by_rank(pop, NP, D, &rnd_uni_init);

			if(parents[0].fitness < parents[1].fitness) {
				offspring = blend_crossover(parents[0], parents[1], D, CR, &rnd_uni_init, inibound_l, inibound_h);
			}
			else {
				offspring = blend_crossover(parents[1], parents[0], D, CR, &rnd_uni_init, inibound_l, inibound_h);
			}

			mutation(offspring[0], D, MR, inibound_l, inibound_h, gen, genmax, &rnd_uni_init);
			mutation(offspring[1], D, MR, inibound_l, inibound_l, gen, genmax, &rnd_uni_init);                              

			new_pop[i] = offspring[0];

			if((i+1) < NP) {
				new_pop[i+1] = offspring[1];
			}

		}

		new_pop[0] = best;

		#pragma omp parallel for
		for(int i = 0; i < NP; i++) {
			pop[i] = new_pop[i];
		}

		#pragma omp parallel for
		for(int i = 0; i < NP; i++) {
			pop[i].fitness = evaluate(pop[i].chromossome); /* obj. funct. value */
		}

		if(gen % 20 == 0) {
			printf("Gen %d Best: \n", gen);		
			for(int k = 0; k < D; k++) {
				printf("%lf ", best.chromossome[k]);
			}
			printf("\n");
			printf("Fit %lf: \n", best.fitness);		
		}

		gen++;

	}
    /*=======================================================================*/
    /*=========End of iteration loop=========================================*/
    /*=======================================================================*/

    return pop;
}

/*-----------End of main()------------------------------------------*/

