#include "de.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define STR_EQ(s1) strcmp(strat, s1) == 0

int get_de_strategy(char *strat) {

	if(STR_EQ("best1exp")) {
		return 1;
	}
	else if(STR_EQ("rand1exp")) {
		return 2;
	}
	else if(STR_EQ("randtobest1exp")) {
		return 3;
	}
	else if(STR_EQ("best2exp")) {
		return 4;
	}
	else if(STR_EQ("rand2exp")) {
		return 5;
	}

	else if(STR_EQ("best1bin")) {
		return 6;
	}
	else if(STR_EQ("rand1bin")) {
		return 7;
	}
	else if(STR_EQ("randtobest1bin")) {
		return 8;
	}
	else if(STR_EQ("best2bin")) {
		return 9;
	}
	else if(STR_EQ("rand2bin")) {
		return 10;
	}

	return 11;

}

double **load_pop(const char *filename, int pop_size, int num_par) {

	double **filedata = NULL;

	filedata = (double**)malloc(sizeof(double*)*pop_size);

	for (int i = 0; i < pop_size; i++) {
		filedata[i] = malloc(sizeof(double)*num_par);
	}


    size_t len = 0;

    FILE *fp;
    fp = fopen(filename, "r");

    if(fp == NULL) {
        fprintf(stderr, "Error reading file %s\n", filename);
        return NULL;
    }

    double data;
	char ch;
	for (int i = 0; i < pop_size; i++) {
		for(int j = 0; j < num_par; j++) {
			fscanf(fp, "%lf", &data);
			filedata[i][j] = data;
		}

		while((ch = fgetc(fp)) != '\n');
		ch=fgetc(fp);
    }

    fclose(fp);

    return filedata;
}

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

double *diffential_evolution(bool from_file, evalfn *evaluate, double *inibound_l, double *inibound_h, int genmax, double CR, double F, int NP, int D, int strategy, long seed, int num_days) {

    /*-----Initialize random number generator-----------------------------*/

    long rnd_uni_init = -(long)seed; /* initialization of rnd_uni() */

    double *cost = (double*)malloc(sizeof(double)*NP);

    double *best = malloc(sizeof(double)*D);

    double *bestit = (double*)malloc(sizeof(double)*D);
    double *tmp = (double*)malloc(sizeof(double)*D);

    double ***pold, ***pnew, ***pswap;

  	double **c = (double **)malloc(sizeof(double*)*NP);
 
    double **d = (double**)malloc(sizeof(double*)*NP);
    for(int i = 0; i < NP; i++) {
        d[i] = malloc(sizeof(double)*D);
    }

    double *original_ind = malloc(sizeof(double)*D);

    /*------Initialization------------------------------------------------*/
    /*------Right now this part is kept fairly simple and just generates--*/
    /*------random numbers in the range [-initfac, +initfac]. You might---*/
    /*------want to extend the init part such that you can initialize-----*/
    /*------each parameter separately.------------------------------------*/

	if(from_file) {
		c = load_pop("population.txt", NP, D);
	}

	else {
		c = (double **)malloc(sizeof(double*)*NP);
		for(int i = 0; i < NP; i++) {
			c[i] = malloc(sizeof(double)*D);
		}
	}


    for(int i = 0; i < NP; i++) {
        for(int j = 0; j < D; j++) {
			if(from_file) {
				c[i][j] = inibound_l[j] + c[i][j] * (inibound_h[j] - inibound_l[j]);
			}
			else {
				c[i][j] = inibound_l[j] + rnd_uni(&rnd_uni_init) * (inibound_h[j] - inibound_l[j]);
			}
        }
	
        cost[i] = evaluate(c[i], num_days); /* obj. funct. value */
    }


    double cmin = cost[0];
    double trial_cost;
    int imin = 0;

    for(int i = 1; i < NP; i++) {
        if(cost[i] < cmin) {
            cmin = cost[i];
            imin = i;
        }
    }

    memcpy(best, c[imin], D * sizeof(double));
    memcpy(bestit, c[imin], D * sizeof(double));

    pold = &c; /* old population (generation G)   */
    pnew = &d; /* new population (generation G+1) */

    int r1, r2, r3, r4, r5; /* placeholders for random indexes    */

    /*=======================================================================*/
    /*=========Iteration loop================================================*/
    /*=======================================================================*/

    int gen = 0;
    int n, L;
    double cmean, cvar;

    while(gen < genmax) {
        gen++;
        imin = 0;

        for(int i = 0; i < NP; i++) {

            do { /* Endless loop for NP < 2 !!!     */
                r1 = (int)(rnd_uni(&rnd_uni_init) * NP);
            } while(r1 == i);

            do { /* Endless loop for NP < 3 !!!     */
                r2 = (int)(rnd_uni(&rnd_uni_init) * NP);
            } while((r2 == i) || (r2 == r1));

            do { /* Endless loop for NP < 4 !!!     */
                r3 = (int)(rnd_uni(&rnd_uni_init) * NP);
            } while((r3 == i) || (r3 == r1) || (r3 == r2));

            do { /* Endless loop for NP < 5 !!!     */
                r4 = (int)(rnd_uni(&rnd_uni_init) * NP);
            } while((r4 == i) || (r4 == r1) || (r4 == r2) || (r4 == r3));

            do { /* Endless loop for NP < 6 !!!     */
                r5 = (int)(rnd_uni(&rnd_uni_init) * NP);
            } while((r5 == i) || (r5 == r1) || (r5 == r2) || (r5 == r3) || (r5 == r4));

            /*=======Choice of strategy===============================================================*/
            /*=======We have tried to come up with a sensible naming-convention: DE/x/y/z=============*/
            /*=======DE :  stands for Differential Evolution==========================================*/
            /*=======x  :  a string which denotes the vector to be perturbed==========================*/
            /*=======y  :  number of difference vectors taken for perturbation of x===================*/
            /*=======z  :  crossover method (exp = exponential, bin = binomial)=======================*/
            /*                                                                                        */
            /*=======There are some simple rules which are worth following:===========================*/
            /*=======1)  F is usually between 0.5 and 1 (in rare cases > 1)===========================*/
            /*=======2)  CR is between 0 and 1 with 0., 0.3, 0.7 and 1. being worth to be tried first=*/
            /*=======3)  To start off NP = 10*D is a reasonable choice. Increase NP if misconvergence=*/
            /*           happens.                                                                     */
            /*=======4)  If you increase NP, F usually has to be decreased============================*/
            /*=======5)  When the DE/best... schemes fail DE/rand... usually works and vice versa=====*/

            /*=======EXPONENTIAL CROSSOVER============================================================*/

            /*-------DE/best/1/exp--------------------------------------------------------------------*/
            /*-------Our oldest strategy but still not bad. However, we have found several------------*/
            /*-------optimization problems where misconvergence occurs.-------------------------------*/

            memcpy(original_ind, (*pold)[i], D * sizeof(double));
            memcpy(tmp, (*pold)[i], D * sizeof(double));

            if(strategy == 1) {
                n = (int)(rnd_uni(&rnd_uni_init) * D);
                L = 0;
                do {
                    tmp[n] = bestit[n] + F * ((*pold)[r2][n] - (*pold)[r3][n]);
                    n = (n + 1) % D;
                    L++;
                } while((rnd_uni(&rnd_uni_init) < CR) && (L < D));

            }
            /*-------DE/rand/1/exp-------------------------------------------------------------------*/
            /*-------This is one of my favourite strategies. It works especially well when the-------*/
            /*-------"bestit[]"-schemes experience misconvergence. Try e.g. F=0.7 and CR=0.5---------*/
            /*-------as a first guess.---------------------------------------------------------------*/
            else if(strategy == 2) /* strategy DE1 in the techreport */
            {
                n = (int)(rnd_uni(&rnd_uni_init) * D);
                L = 0;
                do {
                    tmp[n] = (*pold)[r1][n] + F * ((*pold)[r2][n] - (*pold)[r3][n]);
                    n = (n + 1) % D;
                    L++;
                } while((rnd_uni(&rnd_uni_init) < CR) && (L < D));


            }
            /*-------DE/rand-to-best/1/exp-----------------------------------------------------------*/
            /*-------This strategy seems to be one of the best strategies. Try F=0.85 and CR=1.------*/
            /*-------If you get misconvergence try to increase NP. If this doesn't help you----------*/
            /*-------should play around with all three control variables.----------------------------*/
            else if(strategy == 3) /* similiar to DE2 but generally better */
            {
                n = (int)(rnd_uni(&rnd_uni_init) * D);
                L = 0;
                do {
                    tmp[n] = tmp[n] + F * (bestit[n] - tmp[n]) + F * ((*pold)[r1][n] - (*pold)[r2][n]);
                    n = (n + 1) % D;
                    L++;
                } while((rnd_uni(&rnd_uni_init) < CR) && (L < D));
            }
            /*-------DE/best/2/exp is another powerful strategy worth trying--------------------------*/
            else if(strategy == 4) {
                n = (int)(rnd_uni(&rnd_uni_init) * D);
                L = 0;
                do {
                    tmp[n] = bestit[n] + ((*pold)[r1][n] + (*pold)[r2][n] - (*pold)[r3][n] - (*pold)[r4][n]) * F;
                    n = (n + 1) % D;
                    L++;
                } while((rnd_uni(&rnd_uni_init) < CR) && (L < D));
            }
            /*-------DE/rand/2/exp seems to be a robust optimizer for many functions-------------------*/
            else if(strategy == 5) {
                n = (int)(rnd_uni(&rnd_uni_init) * D);
                L = 0;
                do {
                    tmp[n] = (*pold)[r5][n] + ((*pold)[r1][n] + (*pold)[r2][n] - (*pold)[r3][n] - (*pold)[r4][n]) * F;
                    n = (n + 1) % D;
                    L++;
                } while((rnd_uni(&rnd_uni_init) < CR) && (L < D));
            }

            /*=======Essentially same strategies but BINOMIAL CROSSOVER===============================*/

            /*-------DE/best/1/bin--------------------------------------------------------------------*/
            else if(strategy == 6) {
                n = (int)(rnd_uni(&rnd_uni_init) * D);
                for(L = 0; L < D; L++) /* perform D binomial trials */
                {
                    if((rnd_uni(&rnd_uni_init) < CR) || L == (D - 1)) /* change at least one parameter */
                    {
                        tmp[n] = bestit[n] + F * ((*pold)[r2][n] - (*pold)[r3][n]);
                    }
                    n = (n + 1) % D;
                }

            }

            /*-------DE/rand/1/bin-------------------------------------------------------------------*/
            else if(strategy == 7) {
                n = (int)(rnd_uni(&rnd_uni_init) * D);
                for(L = 0; L < D; L++) /* perform D binomial trials */
                {
                    if((rnd_uni(&rnd_uni_init) < CR) || L == (D - 1)) /* change at least one parameter */
                    {
                        tmp[n] = (*pold)[r1][n] + F * ((*pold)[r2][n] - (*pold)[r3][n]);
                    }
                    n = (n + 1) % D;
                }
            }
            /*-------DE/rand-to-best/1/bin-----------------------------------------------------------*/
            else if(strategy == 8) {
                n = (int)(rnd_uni(&rnd_uni_init) * D);
                for(L = 0; L < D; L++) /* perform D binomial trials */
                {
                    if((rnd_uni(&rnd_uni_init) < CR) || L == (D - 1)) /* change at least one parameter */
                    {
                        tmp[n] = tmp[n] + F * (bestit[n] - tmp[n]) + F * ((*pold)[r1][n] - (*pold)[r2][n]);
                    }
                    n = (n + 1) % D;
                }
            }
            /*-------DE/best/2/bin--------------------------------------------------------------------*/
            else if(strategy == 9) {
                n = (int)(rnd_uni(&rnd_uni_init) * D);
                for(L = 0; L < D; L++) /* perform D binomial trials */
                {
                    if((rnd_uni(&rnd_uni_init) < CR) || L == (D - 1)) /* change at least one parameter */
                    {
                        tmp[n] = bestit[n] + ((*pold)[r1][n] + (*pold)[r2][n] - (*pold)[r3][n] - (*pold)[r4][n]) * F;
                    }
                    n = (n + 1) % D;
                }
            }
            /*-------DE/rand/2/bin--------------------------------------------------------------------*/
            else {
                n = (int)(rnd_uni(&rnd_uni_init) * D);
                for(L = 0; L < D; L++) /* perform D binomial trials */
                {
                    if((rnd_uni(&rnd_uni_init) < CR) || L == (D - 1)) /* change at least one parameter */
                    {
                        tmp[n] =
                            (*pold)[r5][n] + ((*pold)[r1][n] + (*pold)[r2][n] - (*pold)[r3][n] - (*pold)[r4][n]) * F;
                    }
                    n = (n + 1) % D;
                }
            }

			//Check bounds after mutation
            for(int j = 0; j < D; j++) { //----and bounce back----------------------------------------
                if(tmp[j] < inibound_l[j]) {
                    tmp[j] = inibound_l[j] + rnd_uni(&rnd_uni_init) * (original_ind[j] - inibound_l[j]);
                }
                if(tmp[j] > inibound_h[j]) {
                    tmp[j] = inibound_h[j] + rnd_uni(&rnd_uni_init) * (original_ind[j] - inibound_h[j]);
                }
            }

            /*=======Trial mutation now in tmp[]. Test how good this choice really was.==================*/

            trial_cost = evaluate(tmp, num_days); /* Evaluate new vector in tmp[] */

            if(trial_cost <= cost[i]) /* improved objective function value ? */
            {
                cost[i] = trial_cost;
                memcpy((*pnew)[i], tmp, D * sizeof(double));
                if(trial_cost < cmin)  /* Was this a new minimum? */
                {                      /* if so...*/
                    cmin = trial_cost; /* reset cmin to new low...*/
                    imin = i;
                    memcpy(best, tmp, D * sizeof(double));
                }
            } else {
                memcpy((*pnew)[i], (*pold)[i], D * sizeof(double)); /* replace target with old value */
            }
        } /* End mutation loop through pop. */

        memcpy(bestit, best, D * sizeof(double)); /* Save best population member of current iteration */

        /* swap population arrays. New generation becomes old one */

        pswap = pold;
        pold = pnew;
        pnew = pswap;

        /*----Compute the energy variance (just for monitoring purposes)-----------*/

        cmean = 0.; /* compute the mean value first */
        for(int j = 0; j < NP; j++) {
            cmean += cost[j];
        }
        cmean = cmean / NP;

        cvar = 0.; /* now the variance              */
        for(int j = 0; j < NP; j++) {
            cvar += (cost[j] - cmean) * (cost[j] - cmean);
        }
        cvar = cvar / (NP - 1);

         if(gen % 100 == 0) {
            fprintf(stderr, "\n\n\n Best-so-far cost funct. value=%-15.10g\n", cmin);

            for(int j = 0; j < D; j++) {
                fprintf(stderr, "\n best[%d]=%-15.10g", j, best[j]);
            }
            fprintf(stderr, "\n\n Generation=%d ", gen);
            fprintf(stderr, "\n NP=%d    F=%-4.2g    CR=%-4.2g   cost-variance=%-10.5g\n", NP, F, CR, cvar);
         }
		 
    }
    /*=======================================================================*/
    /*=========End of iteration loop=========================================*/
    /*=======================================================================*/

    return best;
}

/*-----------End of main()------------------------------------------*/

