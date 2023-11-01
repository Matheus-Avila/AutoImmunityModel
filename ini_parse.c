#include "ini_parse.h"
#include <string.h>
#include <stdlib.h>

#define IS_TRUE(str) (strcmp((str), "true") == 0 || strcmp((str), "yes") == 0 || strcmp((str), "1") == 0)
#define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0

configuration *new_configuration() {
	configuration *config = calloc(1, sizeof(configuration));

	config->run_two_parts = false;

	config->use_delay = false;
	config->use_ft = false;
	config->n_gens = 1000;
	config->pop_size = 200;

    config->amax_l = 0;
    config->amax_h = 10;

    config->r_l = 0.001;
    config->r_h = 1;

	config->ti_l = -1;
	config->ti_h = -1;

	config->delta_l = -1;
	config->delta_h = -1;

    config->e_l = 0.0;
    config->e_h = 1e-3;

    config->theta_l = 0.001;
    config->theta_h = 0.4;

    config->tau1_l = 2;
    config->tau1_h = 14;

    config->tau2_l = 6;
    config->tau2_h = 22;

    config->tau3_l = 7;
    config->tau3_h = 17;

    config->m_l = 0.01;
    config->m_h = 0.034;

	config->uq_max_error = 0.15;
	config->w_d = 1.0;
	config->w_r = 1.0;
	config->w_c = 1.0;
	config->w_i = 1.0;

	config->de_strat = strdup("best1bin");
	config->de_cr = 0.7;
	config->de_f = 0.5;
	config->de_use_latin_hypercube = true;

	return config;
}

int parse_configuration_file(void* user, const char* section, const char* name, const char* value) {

	configuration* pconfig = (configuration*)user;

	if (MATCH("DE", "num_gens")) {
		pconfig->n_gens = strtol(value, NULL, 10);
	}	
	else if (MATCH("DE", "pop_size")) {
		pconfig->pop_size = strtol(value, NULL, 10);
	}
	else if (MATCH("DE", "F")) {
		pconfig->de_f = strtod(value, NULL);
	}
	else if (MATCH("DE", "CR")) {
		pconfig->de_cr = strtod(value, NULL);
	}
	else if (MATCH("DE", "strategy")) {
		free(pconfig->de_strat);
		pconfig->de_strat = strdup(value);
	}
	else if (MATCH("config", "num_days")) {
		pconfig->num_days = strtol(value, NULL, 10);
	} 
	else if (MATCH("config", "first_day")) {
		pconfig->first_day = strtol(value, NULL, 10);
	} 
	else if (MATCH("config", "first_day_fit_ft")) {
		pconfig->first_day_fit = strtol(value, NULL, 10);
	} 
	else if (MATCH("config", "num_days_fit_ft")) {
		pconfig->num_days_fit = strtol(value, NULL, 10);
	}
	else if (MATCH("config", "population")) {
		pconfig->total_population = strtod(value, NULL);
	}
	else if (MATCH("config", "weight_deaths")) {
		pconfig->w_d = strtod(value, NULL);
	}
	else if (MATCH("config", "weight_confirmed")) {
		pconfig->w_c = strtod(value, NULL);
	}
	else if (MATCH("config", "weight_active")) {
		pconfig->w_i = strtod(value, NULL);
	}
	else if (MATCH("config", "weight_recovered")) {
		pconfig->w_r = strtod(value, NULL);
	}
	else if (MATCH("config", "uq_max_error")) {
		pconfig->uq_max_error = strtod(value, NULL);
	}
	else if (MATCH("config", "use_ft")) {
		pconfig->use_ft = IS_TRUE(value);
	} 
	else if (MATCH("DE", "use_latin_hypercube")) {
		pconfig->de_use_latin_hypercube = IS_TRUE(value);
	}
	else if (MATCH("config", "use_delay")) {
		pconfig->use_delay = IS_TRUE(value);
	}
	else if (MATCH("config", "place")) {
		pconfig->place = strdup(value);
	}
	else if (MATCH("config", "output_dir")) {
		pconfig->out_dir = strdup(value);
	}
	else if (MATCH("config", "confirmed_path")) {
		pconfig->confirmed_path = strdup(value);
	}
	else if (MATCH("config", "active_path")) {
		pconfig->active_path = strdup(value);
	}
	else if (MATCH("config", "deaths_path")) {
		pconfig->deaths_path = strdup(value);
	}
	else if (MATCH("config", "recovered_path")) {
		pconfig->recovered_path = strdup(value);
	}
	else if (MATCH("config", "ft_data_path_path")) {
		pconfig->ft_data_path = strdup(value);
	}
	else if (MATCH("config", "mask_path")) {
		pconfig->mask_path = strdup(value);
	}
	else if (MATCH("model_params", "b_min")) { 
		pconfig->amax_l = strtod(value, NULL);
	}
	else if (MATCH("model_params", "b_max")) { 
		pconfig->amax_h = strtod(value, NULL);
	}
	else if (MATCH("model_params", "r_min")) { 
		pconfig->r_l = strtod(value, NULL);
	}
	else if (MATCH("model_params", "r_max")) { 
		pconfig->r_h = strtod(value, NULL);
	}
	else if (MATCH("model_params", "ti_min")) { 
		pconfig->ti_l = strtod(value, NULL);
	}
	else if (MATCH("model_params", "ti_max")) { 
		pconfig->ti_h = strtod(value, NULL);
	}
	else if (MATCH("model_params", "delta_min")) { 
		pconfig->delta_l = strtod(value, NULL);
	}
	else if (MATCH("model_params", "delta_max")) { 
		pconfig->delta_h = strtod(value, NULL);
	}
	else if (MATCH("model_params", "e_min")) { 
		pconfig->e_l = strtod(value, NULL);
	}
	else if (MATCH("model_params", "e_max")) { 
		pconfig->e_h = strtod(value, NULL);
	}
	else if (MATCH("model_params", "theta_min")) { 
		pconfig->theta_l = strtod(value, NULL);
	}
	else if (MATCH("model_params", "theta_max")) { 
		pconfig->theta_h = strtod(value, NULL);
	}
	else if (MATCH("model_params", "tau1_min")) { 
		pconfig->tau1_l = strtod(value, NULL);
	}
	else if (MATCH("model_params", "tau1_max")) { 
		pconfig->tau1_h = strtod(value, NULL);
	}
	else if (MATCH("model_params", "tau2_min")) { 
		pconfig->tau2_l = strtod(value, NULL);
	}
	else if (MATCH("model_params", "tau2_max")) { 
		pconfig->tau2_h = strtod(value, NULL);
	}
	else if (MATCH("model_params", "tau3_min")) { 
		pconfig->tau3_l = strtod(value, NULL);
	}
	else if (MATCH("model_params", "tau3_max")) { 
		pconfig->tau3_h = strtod(value, NULL);
	}
	else if (MATCH("model_params", "m_min")) { 
		pconfig->m_l = strtod(value, NULL);
	}
	else if (MATCH("model_params", "m_max")) { 
		pconfig->m_h = strtod(value, NULL);
	}

	//RUN CONFIG
	else if (MATCH("run_model_params", "b")) { 
		pconfig->amax_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params", "r")) { 
		pconfig->r_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params", "ti")) { 
		pconfig->ti_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params", "delta")) { 
		pconfig->delta_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params", "e")) { 
		pconfig->e_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params", "theta")) { 
		pconfig->theta_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params", "tau1")) { 
		pconfig->tau1_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params", "tau2")) { 
		pconfig->tau2_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params", "tau3")) { 
		pconfig->tau3_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params", "m")) { 
		pconfig->m_l = strtod(value, NULL);
	}


	//RUN CONFIG P1
	else if (MATCH("run_model_params1", "b")) { 
		pconfig->amax_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params1", "r")) { 
		pconfig->r_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params1", "ti")) { 
		pconfig->ti_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params1", "delta")) { 
		pconfig->delta_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params1", "e")) { 
		pconfig->e_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params1", "theta")) { 
		pconfig->theta_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params1", "tau1")) { 
		pconfig->tau1_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params1", "tau2")) { 
		pconfig->tau2_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params1", "tau3")) { 
		pconfig->tau3_l = strtod(value, NULL);
	}
	else if (MATCH("run_model_params1", "m")) { 
		pconfig->m_l = strtod(value, NULL);
	}

	//RUN CONFIG P2
	else if (MATCH("run_model_params2", "b")) { 
		pconfig->amax_h = strtod(value, NULL);
	}
	else if (MATCH("run_model_params2", "r")) { 
		pconfig->r_h = strtod(value, NULL);
	}
	else if (MATCH("run_model_params2", "ti")) { 
		pconfig->ti_h = strtod(value, NULL);
	}
	else if (MATCH("run_model_params2", "delta")) { 
		pconfig->delta_h = strtod(value, NULL);
	}
	else if (MATCH("run_model_params2", "e")) { 
		pconfig->e_h = strtod(value, NULL);
	}
	else if (MATCH("run_model_params2", "theta")) { 
		pconfig->theta_h = strtod(value, NULL);
	}
	else if (MATCH("run_model_params2", "tau1")) { 
		pconfig->tau1_h = strtod(value, NULL);
	}
	else if (MATCH("run_model_params2", "tau2")) { 
		pconfig->tau2_h = strtod(value, NULL);
	}
	else if (MATCH("run_model_params2", "tau3")) { 
		pconfig->tau3_h = strtod(value, NULL);
	}
	else if (MATCH("run_model_params2", "m")) { 
		pconfig->m_h = strtod(value, NULL);
	}

	else if (MATCH("config", "run_two_parts")) { 
		pconfig->run_two_parts = IS_TRUE(value);
	}
	else if (MATCH("config", "stage_two_start")) { 
		pconfig->stage_two_start = strtol(value, NULL, 10);
	}
	else {
		return 0;  /* unknown section/name, error */
	}
	return 1;
}
