#include <nvector/nvector_serial.h>
#include <math.h>

struct SIR_params {
    realtype amax, ra, ra2, ti, delta, delta2, delta3, m, e, r1, r2, first_day, first_day_fit, c0, c1;
};

static int SIR_PP_V5(realtype t, N_Vector y, N_Vector ydot, void *f_data) {

	realtype S = NV_Ith_S(y, 0);
    realtype I = NV_Ith_S(y, 1);
    realtype R = NV_Ith_S(y, 2);
    realtype D = NV_Ith_S(y, 3);

    struct SIR_params *params = (struct SIR_params *)f_data;

    realtype t_exp = t + fabs(params->first_day - params->first_day_fit);
	realtype ft = params->e * params->c0 * exp(params->c1 * t_exp);

    realtype amax = params->amax;
    realtype ra = params->ra;
    realtype ra2 = params->ra2;
    realtype ti = params->ti;
    realtype delta = params->delta;
    realtype delta2 = params->delta2;
    realtype delta3 = params->delta3;
    realtype m = params->m;
    realtype r1 = params->r1;
    realtype r2 = params->r2;

    realtype alfa_t;

	realtype ti2 = ti + delta + delta2;

    if(t < ti) {
        alfa_t = 1.0;
    } else if(t >= ti && t <= ti + delta) {
    	alfa_t = ( (1.0 - ra) / (-delta) ) * (t - ti) + 1;
	}
	else if( t >= ti + delta && t <= ti2) {
        alfa_t = ra;
    } 
	else if(t >= ti2  && t <= ti2 + delta3) {
    	alfa_t = ( (ra2 - ra) / (delta3) ) * (t - ti2) + ra;
    }
	else {
		alfa_t = ra2;
	}

	realtype alfa = amax*alfa_t;

    realtype beta1 = m * r1;
    realtype beta2 = (1.0 - m) * r2;

    NV_Ith_S(ydot, 0) = -alfa * S * I;
    NV_Ith_S(ydot, 1) = alfa * S * I + ft - beta1 * I - beta2 * I;
    NV_Ith_S(ydot, 2) = beta2 * I;
    NV_Ith_S(ydot, 3) = beta1 * I;

    return 0;
}

static int SIR_PP_V4(realtype t, N_Vector y, N_Vector ydot, void *f_data) {

	realtype S = NV_Ith_S(y, 0);
    realtype I = NV_Ith_S(y, 1);
    realtype R = NV_Ith_S(y, 2);
    realtype D = NV_Ith_S(y, 3);

    struct SIR_params *params = (struct SIR_params *)f_data;

    realtype t_exp = t + fabs(params->first_day - params->first_day_fit);
	realtype ft = params->e * params->c0 * exp(params->c1 * t_exp);

	if(params->e == 0.0) {
		ft = 0.0;
	}

    realtype amax = params->amax;
    realtype r = params->ra;
    realtype ti = params->ti;
    realtype delta = params->delta;
    realtype m = params->m;
    realtype r1 = params->r1;
    realtype r2 = params->r2;

    realtype alfa_t;


    if(t < ti) {
        alfa_t = 1.0;
    } else if(t >= ti && t <= ti + delta) {
    	alfa_t = ( (1.0 - r) / (-delta) ) * (t - ti) + 1;
    } else {
        alfa_t = r;
    }

	realtype alfa = amax*alfa_t;

    realtype beta1 = m * r1;
    realtype beta2 = (1.0 - m) * r2;

	realtype alfaSI = alfa*S*I;
	realtype beta1I = beta1*I;
	realtype beta2I = beta2*I;

    NV_Ith_S(ydot, 0) = -alfaSI;
    NV_Ith_S(ydot, 1) = alfaSI + ft - beta1I - beta2I;
    NV_Ith_S(ydot, 2) = beta2I;
    NV_Ith_S(ydot, 3) = beta1I;

    return 0;
}
