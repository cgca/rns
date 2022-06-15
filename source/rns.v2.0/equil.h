
void make_grid( double s_gp[SDIV+1], 
                double mu[MDIV+1]);                        

void load_eos( char eos_file[], 
               double log_e_tab[201], 
               double log_p_tab[201], 
               double log_h_tab[201],
               double log_n0_tab[201], 
               int *n_tab);

double e_of_rho0(double rho0, double Gamma_P);

double e_at_p(double pp, 
              double log_e_tab[201], 
              double log_p_tab[201],
              int    n_tab, 
              int    *n_nearest_pt,
              char eos_type[],
              double Gamma_P);

double p_at_e(double ee, 
              double log_p_tab[201], 
              double log_e_tab[201],
              int    n_tab, 
              int    *n_nearest_pt);

double p_at_h(double hh, 
              double log_p_tab[201], 
              double log_h_tab[201],
              int    n_tab, 
              int    *n_nearest_pt);

double h_at_p(double pp, 
              double log_h_tab[201], 
              double log_p_tab[201],
              int    n_tab, 
              int    *n_nearest_pt);

double n0_at_e(double ee, 
               double log_n0_tab[201], 
               double log_e_tab[201],
               int    n_tab, 
               int    *n_nearest_pt);

void make_center(
	       char eos_file[], 
               double log_e_tab[201], 
               double log_p_tab[201], 
               double log_h_tab[201],
               double log_n0_tab[201], 
               int n_tab,                 
	       char eos_type[],
	       double Gamma_P, 
	       double e_center,
	       double *p_center, 
	       double *h_center);




void mass_radius(
		 double s_gp[SDIV+1],
		 double mu[MDIV+1],
		 double log_e_tab[201], 
		 double log_p_tab[201], 
		 double log_h_tab[201],
		 double log_n0_tab[201], 
		 int n_tab,                 
		 char eos_type[],
		 double Gamma_P, 
		 double **rho,
		 double **gama,
		 double **alpha,
		 double **omega,
		 double **energy,
		 double **pressure,
		 double **enthalpy,
		 double **velocity_sq,
                 double r_ratio,
		 double e_surface,
                 double r_e,
                 double Omega,
                 double *Mass, 
		 double *Mass_0,
		 double *ang_mom,
                 double *R_e,
		 double *v_plus,
		 double *v_minus,
		 double *Omega_K);

double dm_dr_is(double r_is, 
                double r, 
                double m, 
                double p, 
                double e_center, 
                double p_surface,
                double log_e_tab[SDIV+1],
                double log_p_tab[SDIV+1],
                int    n_tab,
                int    *n_nearest_pt,
                char   eos_type[],
                double Gamma_P);

double dp_dr_is(double r_is, 
                double r, 
                double m, 
                double p,
                double e_center, 
                double p_surface,
                double log_e_tab[SDIV+1],
                double log_p_tab[SDIV+1],
                int    n_tab,
                int    *n_nearest_pt,
                char eos_type[],
                double Gamma_P);

double dr_dr_is(double r_is, double r, double m);

void TOV(int    i_check, 
               char   eos_type[],
               double e_center,
               double p_center,
               double p_surface,
               double e_surface,
               double Gamma_P, 
               double log_e_tab[201],
               double log_p_tab[201],
               double log_h_tab[201],
               int    n_tab,
               double r_is_gp[RDIV+1], 
               double lambda_gp[RDIV+1], 
               double nu_gp[RDIV+1], 
               double *r_is_final, 
               double *r_final, 
               double *m_final);

void sphere(double s_gp[SDIV+1], 
	    double log_e_tab[201], 
	    double log_p_tab[201], 
	    double log_h_tab[201],
	    double log_n0_tab[201], 
	    int n_tab,                 
	    char eos_type[],
	    double Gamma_P, 
	    double e_center,
	    double p_center, 
	    double h_center,
	    double p_surface,
	    double e_surface,
	    double **rho,
	    double **gama,
	    double **alpha,
	    double **omega,
	    double *r_e);


void spin(double s_gp[SDIV+1],
	  double mu[MDIV+1],
	  double log_e_tab[201], 
	  double log_p_tab[201], 
	  double log_h_tab[201],
	  double log_n0_tab[201], 
	  int n_tab,                 
	  char eos_type[],
	  double Gamma_P, 
	  double h_center,
	  double enthalpy_min,
	  double **rho,
	  double **gama,
	  double **alpha,
	  double **omega,
	  double **energy,
	  double **pressure,
	  double **enthalpy,
	  double **velocity_sq,
	  int    a_check, 
	  double accuracy,
	  double cf,
	  double r_ratio,
	  double *r_e_new,
	  double *Omega) ;
 
