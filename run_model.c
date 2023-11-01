#include "de.h"
#include "ini_parse.h"
#include "sirmodel.h"
#include <cvode/cvode.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"

#include "sds/sds.h"

#define NEQ 4

realtype *load_data(const char *filename);

realtype c0, c1;

realtype *confirmed = NULL;
realtype *infected = NULL;
realtype *deaths = NULL;
realtype *recovery = NULL;
realtype *mask = NULL;
realtype **best_de = NULL;

configuration *config;

int find_nearest(int num_days, realtype value) {

    realtype min = DBL_MAX;
    int min_index = 0;
    realtype aux;

    for(int i = 0; i < num_days; i++) {
        aux = fabs((realtype)i - value);
        if(aux < min) {
            min = aux;
            min_index = i;
        }
    }
    return min_index;
}

realtype *load_data(const char *filename) {

    printf("Loading %s\n", filename);

    realtype *filedata = NULL;

    size_t len = 0;

    FILE *fp;
    realtype data;

    fp = fopen(filename, "r");

    if(fp == NULL) {
        fprintf(stderr, "Error reading file %s\n", filename);
        return NULL;
    }

    char *line = NULL;
    while((getline(&line, &len, fp)) != -1) {
        if(line) {
            data = strtod(line, NULL);
            arrput(filedata, data);
        }
    }

    free(line);
    fclose(fp);

    return filedata;
}

static int check_flag(void *flagvalue, const char *funcname, int opt) {

    int *errflag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if(opt == 0 && flagvalue == NULL) {

        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return (1);
    }

    /* Check if flag < 0 */
    else if(opt == 1) {
        errflag = (int *)flagvalue;
        if(*errflag < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
            return (1);
        }
    }

    /* Check if function returned NULL pointer - no memory allocated */
    else if(opt == 2 && flagvalue == NULL) {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return (1);
    }

    return (0);
}

void exp_mmq(realtype *dia, realtype *casos, realtype *c0, realtype *c1, int n) {

    N_Vector phi_0 = N_VNew_Serial(n);

    for(int i = 0; i < n; i++) {
        NV_Ith_S(phi_0, i) = 1.0;
    }

    N_Vector phi_1 = N_VNew_Serial(n);
    for(int i = 0; i < n; i++) {
        NV_Ith_S(phi_1, i) = dia[i];
    }

    N_Vector F = N_VNew_Serial(n);
    for(int i = 0; i < n; i++) {
        NV_Ith_S(F, i) = log(casos[i]);
    }

    SUNMatrix A = SUNDenseMatrix(2, 2);

    SM_ELEMENT_D(A, 0, 0) = N_VDotProd_Serial(phi_0, phi_0);
    SM_ELEMENT_D(A, 0, 1) = N_VDotProd_Serial(phi_0, phi_1);
    SM_ELEMENT_D(A, 1, 0) = N_VDotProd_Serial(phi_1, phi_0);
    SM_ELEMENT_D(A, 1, 1) = N_VDotProd_Serial(phi_1, phi_1);

    N_Vector b = N_VNew_Serial(2);
    NV_Ith_S(b, 0) = N_VDotProd_Serial(F, phi_0);
    NV_Ith_S(b, 1) = N_VDotProd_Serial(F, phi_1);

    SUNLinearSolver LS = SUNDenseLinearSolver(b, A);
    SUNLinSolSetup_Dense(LS, A);

    N_Vector x = N_VNew_Serial(2);

    SUNLinSolSolve_Dense(LS, A, x, b, 1e-16);

    N_Vector y = N_VNew_Serial(2);
    SUNMatMatvec_Dense(A, x, y);

    *c0 = exp(NV_Ith_S(x, 0));
    *c1 = NV_Ith_S(x, 1);

    N_VDestroy(phi_0);
    N_VDestroy(phi_1);
    N_VDestroy(b);
    N_VDestroy(x);
    N_VDestroy(F);
    N_VDestroy(y);
    SUNMatDestroy(A);
    SUNLinSolFree(LS);
}

void solve_ode(N_Vector y, int num_days, realtype **S, realtype **I, realtype **R, realtype **D, struct SIR_params params, struct SIR_params params2) {

    void *cvode_mem = NULL;
    int flag;

    // Set up solver
    cvode_mem = CVodeCreate(CV_BDF);

    if(cvode_mem == 0) {
        fprintf(stderr, "Error in CVodeMalloc: could not allocate\n");
        return;
    }

    flag = CVodeInit(cvode_mem, SIR_PP_V4, 0, y);
    if(check_flag(&flag, "CVodeInit", 1))
        return;

    flag = CVodeSStolerances(cvode_mem, 1.49012e-9, 1.49012e-9);
    if(check_flag(&flag, "CVodeSStolerances", 1))
        return;

    /* Provide RHS flag as user data which can be access in user provided routines */
    flag = CVodeSetUserData(cvode_mem, (void *)&params);
    if(check_flag(&flag, "CVodeSetUserData", 1))
        return;

    // Create dense SUNMatrix for use in linear solver
    SUNMatrix A = SUNDenseMatrix(NEQ, NEQ);
    if(check_flag((void *)A, "SUNDenseMatrix", 0))
        return;

    // Create dense linear solver for use by CVode
    SUNLinearSolver LS = SUNLinSol_Dense(y, A);
    if(check_flag((void *)LS, "SUNLinSol_Dense", 0))
        return;

    // Attach the linear solver and matrix to CVode by calling CVodeSetLinearSolver
    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if(check_flag((void *)&flag, "CVodeSetLinearSolver", 1))
        return;

    int iout = 1;
    int retval;
    realtype t;

    arrput(*S, NV_Ith_S(y, 0));
    arrput(*I, NV_Ith_S(y, 1));
    arrput(*R, NV_Ith_S(y, 2));
    arrput(*D, NV_Ith_S(y, 3));

    while(iout < config->stage_two_start) {

        retval = CVode(cvode_mem, iout, y, &t, CV_NORMAL);

        if(retval == CV_SUCCESS) {
			arrput(*S, NV_Ith_S(y, 0));
			arrput(*I, NV_Ith_S(y, 1));
			arrput(*R, NV_Ith_S(y, 2));
			arrput(*D, NV_Ith_S(y, 3));
            iout++;
        }
    }

	if(config->run_two_parts) {

		realtype not1 = config->theta_l;
		realtype not2 = config->theta_h;

		N_Vector P0 = N_VNew_Serial(NEQ);
		NV_Ith_S(P0, 0) = config->total_population - ((NV_Ith_S(y,1)*not1)/not2 + (NV_Ith_S(y,2)*not1)/not2 + NV_Ith_S(y,3));
		NV_Ith_S(P0, 1) = (NV_Ith_S(y,1)*not1)/not2;
		NV_Ith_S(P0, 2) = (NV_Ith_S(y,2)*not1)/not2;
		NV_Ith_S(P0, 3) = NV_Ith_S(y,3);


		flag = CVodeReInit(cvode_mem, t, P0);
		if(check_flag(&flag, "CVodeInit", 1))
			return;

		//Provide RHS flag as user data which can be access in user provided routines 
		flag = CVodeSetUserData(cvode_mem, (void *)&params2);
		if(check_flag(&flag, "CVodeSetUserData", 1))
			return;

		while(iout < config->num_days) {

			retval = CVode(cvode_mem, iout, y, &t, CV_NORMAL);

			if(retval == CV_SUCCESS) {
				arrput(*S, NV_Ith_S(y, 0));
				arrput(*I, NV_Ith_S(y, 1));
				arrput(*R, NV_Ith_S(y, 2));
				arrput(*D, NV_Ith_S(y, 3));
				iout++;
			}
		}

	}

	// Free the linear solver memory
	SUNLinSolFree(LS);
	SUNMatDestroy(A);
	CVodeFree(&cvode_mem);

}

bool print_error = false;

#define MAX_NORM 1
#define L1_NORM 2

int selected_norm = L1_NORM;

static double model(realtype *x, realtype *x2, int num_days, bool run_two) {

	int first_day = config->first_day;
	realtype total_population = config->total_population;

	realtype not = x[5];
	realtype not2;
	realtype tau_1 = x[6];
	realtype tau_2 = x[7];
	realtype tau_3 = x[8];
	realtype m = x[9];

	N_Vector P0 = N_VNew_Serial(NEQ);
	NV_Ith_S(P0, 0) = total_population - infected[first_day] / not;
	NV_Ith_S(P0, 1) = infected[first_day] / not;
	NV_Ith_S(P0, 2) = recovery[first_day] / not;
	NV_Ith_S(P0, 3) = deaths[first_day];

	realtype r1 = 1.0 / (tau_2 + tau_1);
	realtype r2 = 1.0 / (tau_3 + tau_1);

	struct SIR_params params;
	struct SIR_params params2;

	params.amax = x[0] / total_population;
	params.ra = x[1];
	params.ti = x[2];
	params.delta = x[3];
	params.m = m;
	params.e = x[4];
	params.r1 = r1;
	params.r2 = r2;
	params.first_day = first_day;
	params.first_day_fit = config->first_day_fit;
	params.c0 = c0;
	params.c1 = c1;

	if(config->run_two_parts) {
		not2 = x2[5];
		tau_1 = x2[6];
		tau_2 = x2[7];
		tau_3 = x2[8];
		m = x2[9];

		r1 = 1.0 / (tau_2 + tau_1);
		r2 = 1.0 / (tau_3 + tau_1);

		params2.amax = x2[0] / total_population;
		params2.ra = x2[1];
		params2.ti = config->stage_two_start + x2[2];
		params2.delta = x2[3];
		params2.m = m;
		params2.e = x2[4];
		params2.r1 = r1;
		params2.r2 = r2;
	}

	realtype *S = NULL;
	realtype *I = NULL;
	realtype *R = NULL;
	realtype *D = NULL;

	solve_ode(P0, num_days, &S, &I, &R, &D, params, params2);

	N_Vector In = N_VNew_Serial(arrlen(I));
	N_Vector Rn = N_VNew_Serial(arrlen(R));

	N_Vector I_CV = N_VMake_Serial(arrlen(I), I);
	N_Vector R_CV = N_VMake_Serial(arrlen(R), R);
	N_Vector D_CV = N_VMake_Serial(arrlen(D), D);

	if(config->use_delay) {

		for(int i = 0; i < arrlen(I); i++) {
			NV_Ith_S(In, i) = 0.0;
		}

		for(int i = 0; i < arrlen(R); i++) {
			NV_Ith_S(Rn, i) = 0.0;
		}

		int ii = find_nearest(num_days, (1 - not) * tau_1);
		// int ii = find_nearest(num_days + opt_last, tau_1);

		if((ii - (1 - not) * tau_1) < 0.0)
			ii = ii + 1;

		for(int i = 0; i < ii; i++) {
			NV_Ith_S(In, i) = infected[first_day];
		}

		for(int i = ii; i < arrlen(I); i++) {
			int j = find_nearest(num_days, i - (1 - not) * tau_1);
			// int j = find_nearest(num_days + opt_last, i-tau_1);
			NV_Ith_S(In, i) = not * I[j];
			NV_Ith_S(Rn, i) = not * R[j];
		}
	} else {
		if(config->run_two_parts) {

			realtype n = not;
			for(int i = 0; i < arrlen(I); i++) {
				if(i >= config->stage_two_start) {
					n = not2;
				}
				NV_Ith_S(In, i) = n * NV_Ith_S(I_CV, i);
				NV_Ith_S(Rn, i) = n * NV_Ith_S(R_CV, i);

			}
		}

		else {
			N_VScale_Serial(not, I_CV, In);
			N_VScale_Serial(not, R_CV, Rn);
		}
	}

	// C = In+Rn+D;
	N_Vector C = N_VNew_Serial(arrlen(I));
	N_VLinearSum_Serial(RCONST(1.0), Rn, RCONST(1.0), In, C);
	N_VLinearSum_Serial(RCONST(1.0), D_CV, RCONST(1.0), C, C);

	N_Vector c = N_VMake_Serial(num_days, confirmed + first_day);
	N_Vector cminusC = N_VNew_Serial(num_days);
	N_VLinearSum_Serial(RCONST(1.0), c, RCONST(-1.0), C, cminusC);

	realtype n1, n2;

	if(selected_norm == MAX_NORM) {
		n1 = N_VMaxNorm_Serial(cminusC);
		n2 = N_VMaxNorm_Serial(c);
	} else {
		n1 = N_VL1Norm_Serial(cminusC);
		n2 = N_VL1Norm_Serial(c);
	}

	realtype erro_conf = n1 / n2;

	int n = arrlen(infected);

	N_Vector infected_CV = N_VMake_Serial(n, infected);
	N_Vector mask_CV = N_VMake_Serial(n, mask);
	N_Vector mask_CV_shifted = N_VMake_Serial(num_days, mask + first_day);

	N_Vector infected_aux = N_VNew_Serial(n);
	N_VProd_Serial(infected_CV, mask_CV, infected_aux);

	N_Vector In_aux = N_VNew_Serial(num_days);
	N_VProd_Serial(In, mask_CV_shifted, In_aux);

	N_Vector c2 = N_VNew_Serial(num_days);

	N_Vector infected_aux_shifted = N_VMake_Serial(num_days, NV_DATA_S(infected_aux) + first_day);
	N_VLinearSum_Serial(RCONST(1.0), infected_aux_shifted, RCONST(-1.0), In_aux, c2);

	if(selected_norm == MAX_NORM) {
		n1 = N_VMaxNorm_Serial(c2);
		n2 = N_VMaxNorm_Serial(infected_aux_shifted);
	} else {
		n1 = N_VL1Norm_Serial(c2);
		n2 = N_VL1Norm_Serial(infected_aux_shifted);
	}

	realtype erro_inf = n1 / n2;

	N_Vector deaths_CV_shifted = N_VMake_Serial(num_days, deaths + first_day);
	N_VLinearSum_Serial(RCONST(1.0), deaths_CV_shifted, RCONST(-1.0), D_CV, c2);

	if(selected_norm == MAX_NORM) {
		n1 = N_VMaxNorm_Serial(c2);
		n2 = N_VMaxNorm_Serial(deaths_CV_shifted);
	} else {
		n1 = N_VL1Norm_Serial(c2);
		n2 = N_VL1Norm_Serial(deaths_CV_shifted);
	}

	realtype erro_deaths = n1 / n2;

	N_Vector recovery_CV = N_VMake_Serial(n, recovery);
	N_Vector recovery_aux = N_VNew_Serial(n);
	N_VProd_Serial(recovery_CV, mask_CV, recovery_aux);

	N_Vector Rn_aux = N_VNew_Serial(num_days);
	N_Vector Rn_shifted = N_VMake_Serial(num_days, NV_DATA_S(Rn));
	N_VProd_Serial(Rn_shifted, mask_CV_shifted, Rn_aux);

	N_Vector recovery_aux_shifted = N_VMake_Serial(num_days, NV_DATA_S(recovery_aux) + first_day);

	N_VLinearSum_Serial(RCONST(1.0), recovery_aux_shifted, RCONST(-1.0), Rn_aux, c2);

	if(selected_norm == MAX_NORM) {
		n1 = N_VMaxNorm_Serial(c2);
		n2 = N_VMaxNorm_Serial(recovery_aux_shifted);
	} else {
		n1 = N_VL1Norm_Serial(c2);
		n2 = N_VL1Norm_Serial(recovery_aux_shifted);
	}

	realtype erro_recovery = n1 / n2;
	if(isnan(erro_deaths) || isinf(erro_deaths)) erro_deaths = 0.0;
	if(isnan(erro_conf) || isinf(erro_conf)) erro_conf = 0.0;
	if(isnan(erro_recovery) || isinf(erro_recovery)) erro_recovery = 0.0;
	if(isnan(erro_inf) || isinf(erro_inf)) erro_inf = 0.0;

	realtype error = config->w_d*erro_deaths + config->w_c*erro_conf + config->w_i*erro_inf + config->w_r*erro_recovery;
	
	if(print_error) {
		printf("Error confirmed %lf\n", erro_conf);
		printf("Error active %lf\n",    erro_inf);
		printf("Error deaths %lf\n",    erro_deaths);
		printf("Error recovery %lf\n",  erro_recovery);

		sds fname =
			sdscatprintf(sdsempty(), "%s/I_%s_f_%d_n_%d.txt", config->out_dir, config->place, first_day, num_days);

		FILE *f = fopen(fname, "w");
		N_VPrintFile_Serial(I_CV, f);
		fclose(f);

		sdsfree(fname);

		fname = sdscatprintf(sdsempty(), "%s/In_%s_f_%d_n_%d.txt", config->out_dir, config->place, first_day, num_days);

		f = fopen(fname, "w");
		N_VPrintFile_Serial(In, f);
		fclose(f);

		sdsfree(fname);
		fname = sdscatprintf(sdsempty(), "%s/D_%s_f_%d_n_%d.txt", config->out_dir, config->place, first_day, num_days);

		f = fopen(fname, "w");
		N_VPrintFile_Serial(D_CV, f);
		fclose(f);

		sdsfree(fname);
		fname = sdscatprintf(sdsempty(), "%s/Rn_%s_f_%d_n_%d.txt", config->out_dir, config->place, first_day, num_days);

		f = fopen(fname, "w");
		N_VPrintFile_Serial(Rn, f);
		fclose(f);

		sdsfree(fname);
		fname = sdscatprintf(sdsempty(), "%s/C_%s_f_%d_n_%d.txt", config->out_dir, config->place, first_day, num_days);

		f = fopen(fname, "w");
		N_VPrintFile_Serial(C, f);

		sdsfree(fname);
		fclose(f);
	}

	N_VDestroy(P0);
	N_VDestroy(In);
	N_VDestroy(Rn);
	N_VDestroy(I_CV);
	N_VDestroy(R_CV);
	N_VDestroy(D_CV);
	N_VDestroy(C);
	N_VDestroy(c);
	N_VDestroy(c2);
	N_VDestroy(cminusC);
	N_VDestroy(infected_CV);
	N_VDestroy(mask_CV);
	N_VDestroy(mask_CV_shifted);
	N_VDestroy(infected_aux);
	N_VDestroy(In_aux);
	N_VDestroy(infected_aux_shifted);
	N_VDestroy(deaths_CV_shifted);
	N_VDestroy(recovery_CV);
	N_VDestroy(recovery_aux);
	N_VDestroy(recovery_aux_shifted);
	N_VDestroy(Rn_aux);

	arrfree(S);
	arrfree(I);
	arrfree(R);
	arrfree(D);

	return error;
}

int main(int argc, char **argv) {

    config = new_configuration();

    if(argc != 2) {
        fprintf(stderr, "Usage %s config_file.ini\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    if(ini_parse(argv[1], parse_configuration_file, config) < 0) {
        printf("Can't load '%s'\n", argv[1]);
        return 1;
    }

    if(config->out_dir == NULL) {
        config->out_dir = strdup("./");
    } else {
        int status;
        status = mkdir(config->out_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if(status == -1) {
            if(errno != EEXIST) {
                fprintf(stderr, "Error creating output dir %s\n", config->out_dir);
                perror("Error was:");
            }
        }
    }

    confirmed = load_data(config->confirmed_path);
    infected = load_data(config->active_path);
    deaths = load_data(config->deaths_path);
    recovery = load_data(config->recovered_path);

    if(config->mask_path) {
        mask = load_data(config->mask_path);
    } else {
        printf("Making mask\n");
        for(int i = 0; i < arrlen(confirmed); i++) {
            arrput(mask, 1.0);
        }
    }


    if(config->num_days <= 0) {
        config->num_days = arrlen(infected) - config->first_day;
	}
    
	if(config->num_days_fit <= 0) {
        config->num_days_fit = arrlen(infected) - config->first_day_fit;
    } 

	if(config->first_day + config->num_days > arrlen(infected)) {
		fprintf(stderr, "first_day + num_days is greater than the data available. Reajusting from %d to %ld: \n", config->num_days, arrlen(infected) - config->first_day);
		config->num_days = arrlen(infected) - config->first_day;
	}

	if(config->num_days <= 0) {
		fprintf(stderr, "Invalid num_days parameter: %d\n", config->num_days);
	}

    if(config->use_ft) {
        realtype *casos = load_data(config->ft_data_path);

        realtype *dia = NULL;

        for(int i = config->first_day_fit; i < config->first_day_fit + config->num_days_fit; i++) {
            arrput(dia, i);
        }

        exp_mmq(dia, casos + config->first_day_fit, &c0, &c1, config->num_days_fit);
        arrfree(dia);
        arrfree(casos);
        printf("\nc0 %lf, c1 %lf\n", c0, c1);
    }

    if(!config->use_ft) {
        config->e_l = 0.0;
        config->e_h = 0.0;
    }

	printf("---------------------------------------------\n");
	printf("RUN INFO:\n");

		printf("\nPlace %s, Population %ld, First day %d, num_days: %d, use delay %d, use f(t) %d\n\n", config->place, 
				(long int)config->total_population, config->first_day, config->num_days, config->use_delay, config->use_ft);


    printf("Error function\n error = %lf*error_confirmed + %lf*error_active + %lf*error_recovered + %lf*error_deaths\n\n", config->w_c, config->w_i, config->w_r, config->w_d);

    int num_par = 10;

	realtype x[] = {config->amax_l, config->r_l, config->ti_l, config->delta_l, config->e_l, config->theta_l, config->tau1_l, config->tau2_l, config->tau3_l, config->m_l};
	realtype x2[] = {config->amax_h, config->r_h, config->ti_h, config->delta_h, config->e_h, config->theta_h, config->tau1_h, config->tau2_h, config->tau3_h, config->m_h};

    char *names[10] = {"b", "r", "ti", "delta", "e", "not", "tau_1", "tau_2", "tau_3", "m"};
	if(config->run_two_parts) {
		printf("P1 params:\n");
	}
	else {
		config->stage_two_start = config->num_days;
	}

    for(int i = 0; i < num_par; i++) {
        printf("%s: %e\n", names[i], x[i]);
    }
	
	printf("==============================\n");

	if(config->run_two_parts) {
		printf("P2 params:\n");
		for(int i = 0; i < num_par; i++) {
			printf("%s: %e\n", names[i], x2[i]);
		}
	}

    printf("\n");
    print_error = true;
    model(x, x2, config->num_days, config->run_two_parts);

    sds command = sdsempty();

    command = sdscatprintf(command, "python ./python_scripts/plot.py %s %s %s %s %s %d %d %s", config->place,
                           config->confirmed_path, config->active_path, config->deaths_path, config->recovered_path,
                           config->first_day, config->num_days, config->out_dir);
    printf("%s\n", command);
    system(command);

    sdsfree(command);

    return (0);
}
