#ifndef __INI_PARSE_H
#define __INI_PARSE_H 

#include <stdbool.h>
#include "ini.h"

typedef struct {
	double total_population;
	int first_day;
	int num_days;
	int stage_two_start;
	int first_day_fit;
	int num_days_fit;
	char *de_strat;
	char *place;
	char *confirmed_path;
    char *active_path;
    char *deaths_path;
    char *recovered_path;
    char *ft_data_path;
    char *mask_path;
	char *out_dir;
	bool use_delay;
	bool use_ft;
	bool de_use_latin_hypercube;
	bool run_two_parts;
	int n_gens;
	int pop_size;
	double w_d;
	double w_c;
	double w_i;
	double w_r;
	double uq_max_error;
	double de_f;
	double de_cr;
	

    // amax_r = b*r
    double amax_l;
    double amax_h;

    double r_l;
    double r_h;

	double ti_l;
	double ti_h;

	double delta_l;
	double delta_h;

    double e_l;
    double e_h;

    double theta_l;
    double theta_h;

    double tau1_l;
    double tau1_h;

    double tau2_l;
    double tau2_h;

    double tau3_l;
    double tau3_h;

    double m_l;
    double m_h;
	

} configuration;

configuration *new_configuration();
int parse_configuration_file(void* user, const char* section, const char* name, const char* value);

#endif /* __INI_PARSE_H */
