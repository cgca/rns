/************************************************************************** 
*                            MAIN.C                                       * 
*                                                                         *
*       This is a sample program which uses the routines 
*       written by Nikolaos Stergioulas. 
*
* This sample program is meant merely to serve as an 
* example of how to use the functions! 
*
* In this example, the user specifies the star's equation of state,
* and the central energy density. The program computes models with
* increasing angular velocity until the star is spinning with 
* the same angular velocity as a particle orbiting the star
* at its equator. The last star computed will be spinning with
* the maximum allowed angular velocity for a star with the given
* central energy density. 
*
* Note that there is no guarantee that any of the intermediate
* stars are stable to radial perturbations. Also there is no
* guarantee that given any energy density, there will be a 
* stable rotating solution. As a rule of thumb, the highest 
* energy density you should use is the value which gives the
* maximum mass nonrotating star. 
*
* It would be a good idea to read some of the papers on rapidly
* rotating neutron stars (such as given on the rns homepage)
* before embarking on a study of rotating neutron stars
* using these routines. 
* 
* For example, a reasonable model would use the file eosA,
* central energy density = 10^{15} g/cm^3 
* To specify these parameters, run the executable:
 
kepler -f eosA -e 1e15 

* This program (compiled on the "standard" setting) 
* requires about 2.7 MBytes of RAM and took about 2 minutes to run.
*                                                                         *
**************************************************************************/
#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "consts.h"
#include "nrutil.h"
#include "equil.h"

/***************************************************************************
 * printing routine  
 ***************************************************************************/
									    
void print(double r_ratio,
	   double e_center, double Mass, double Mass_0, double R_e,
	   double Omega, double Omega_K, double J
	   )
{

  double I_45;
  
  if( Omega == 0.0) I_45 = 0.0;
  else I_45 = J/(Omega*1.0e45);

  printf(
	 "%4.3f \t%4.1f \t%4.3f \t%4.3f \t%4.3f \t%4.1f \t%5.1f \t%4.2f \t%4.3f \n",
	 r_ratio,
	 e_center,
	 Mass/MSUN,
	 Mass_0/MSUN,
	 R_e/1.0e5,
	 Omega,
	 Omega_K,
	 I_45,
	 ( C*J/(G*Mass*Mass)))
    ;
}

/***************************************************************************
 * printing routine for polytropic stars
 ***************************************************************************/

void printpoly(double r_ratio,
	   double e_center, double Mass, double Mass_0, double R_e,
	   double Omega, double Omega_K, double J
	   )
{

  double I;
  
  if( Omega == 0.0) I = 0.0;
  else I = J/(Omega);

  printf(
	 "%4.3f \t%4.3f \t%4.3f \t%4.3f \t%4.3f \t%4.3f \t%5.3f \t%4.2f \t%4.3f \n",
	 r_ratio,
	 e_center,
	 Mass,
	 Mass_0,
	 R_e,
	 Omega,
	 Omega_K,
	 I,
	 J/(Mass_0*Mass_0))
    ;
}


/*************************************************************************/
/* Main program.                                                         */
/*************************************************************************/
int main(int argc,                    /* Number of command line arguments */ 
         char **argv)                 /* Command line arguments */
{


 /* EQUILIBRIUM VARIABLES */

 int    n_tab;                     /* Number of points in EOS file */
       

 int i,
     j;

 double log_e_tab[201],               /* energy density/c^2 in tabulated EOS */
   log_p_tab[201],               /* pressure in tabulated EOS */
   log_h_tab[201],               /* enthalpy in EOS file */
   log_n0_tab[201],              /* number density in EOS file */  
   e_center,                     /* central en. density */
   p_center,                     /* central pressure */
   h_center,                     /* central enthalpy */
   n_P,                          /* Polytropic index N */
   Gamma_P,                      /* Gamma for polytropic EOS */     
   r_ratio,                      /* axis ratio */
   s_gp[SDIV+1],                 /* s grid points */
   mu[MDIV+1],                   /* \mu grid points */
   **rho,                          /* potential \rho */ 
   **gama,                         /* potential \gamma */ 
   **omega,                        /* potential \omega */ 
   **alpha,                        /* potential \alpha */ 
   Mass,                         /* Gravitational mass */
   e_surface,                    /* surface en. density */ 
   p_surface,                    /* surface pressure */
   enthalpy_min,                 /* minimum enthalpy in EOS */
   **energy,                       /* energy density \epsilon */
   **pressure,                     /* pressure */ 
   **enthalpy,                     /* enthalpy */
   **velocity_sq,                  /* square of velocity */ 
   Mass_0,                         /* Baryon Mass */
   Omega,			   /* Angular Velocity */
   J,				/* Angular Momentum */
   R_e,                          /* Circumferential radius at equator */
   *v_plus,			/* vel. of co-rot. particle wrt ZAMO */
   *v_minus,			/* vel. of counter-rot. ... */
   Omega_K, /* Keplerian velocity of particle orbiting at equator */
   r_e                         /* coord. radius at equator 	*/
   ;

 double
   cf=1, /* convergence factor */
   accuracy;

 /* Definitions used to find the interval containing the 
    correct spin frequency */

 double dr,
   diff_Omega,
   old_diff_Omega,
   a_check;
   
 /* Definitions used in the Ridder zero finding method */
 
 float ans, fh,fl,fm,fnew,sroot,xh,xl,xm,xnew,xacc;


  char eos_file[80] = "no EOS file specified";   /* EOS file name */
  char eos_type[80] = "tab";                    /* EOS type (poly or tab) */
 

  /* READ IN THE COMMAND LINE OPTIONS */

  for(i=1;i<argc;i++) 
    if(argv[i][0]=='-'){
      switch(argv[i][1]){

      case 'q':
	/* CHOOSE THE EOS TYPE: EITHER "tab" or "poly"
	   (default is tab) */
	sscanf(argv[i+1],"%s",eos_type);
	break;
               
      case 'N':
	/* IF A POLYTROPIC EOS WAS CHOSEN, CHOOSE THE 
	   POLYTROPIC INDEX "N" */
	sscanf(argv[i+1],"%lf",&n_P);
	Gamma_P=1.0+1.0/n_P;
	break;               

      case 'f':
	/* IF A TABULATED EOS WAS CHOSEN, CHOOSE THE
	   NAME OF THE FILE */
	sscanf(argv[i+1],"%s",eos_file);
	break;

      case 'e':
	/* CHOOSE THE CENTRAL ENERGY DENSITY OF THE 
	   NEUTRON STAR (IN g/cm^3) */
	sscanf(argv[i+1],"%lf",&e_center);
	if(strcmp(eos_type,"tab")==0)
	  e_center *= C*C*KSCALE;
	break;
		
      case 'h': 
	fprintf(stderr,"\n");
	fprintf(stderr,"Quick help:\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"  -q EOS type (tab)\n"); 
	fprintf(stderr,"     tab : tabulated \n");
	fprintf(stderr,"     poly : analytic polytropic \n");           
	fprintf(stderr,"  -N polytropic index (P=K*e^(1+1/N))\n");  
	fprintf(stderr,"  -f EOS file \n");
	fprintf(stderr,"  -e central energy density in gr/cm^3\n");
	fprintf(stderr,"  -h this menu\n");
	fprintf(stderr,"\n");
	exit(1);
	break;
     
      }
    }

  /* PRINT THE HEADER */

  if(strcmp(eos_type,"tab")==0)
    printf("%s,  MDIVxSDIV=%dx%d\n",eos_file,MDIV,SDIV);
  else printf("polytrope with N=%f, MDIVxSDIV=%dx%d\n",n_P,MDIV,SDIV);
  printf("ratio\te_15\tM\tM_0\tr_star\tspin\tOmega_K\tI\tJ/M^2\n");
  if(strcmp(eos_type,"tab")==0)
    printf("\tg/cm^3\tsun\tsun\tkm\ts-1\ts-1\tg cm^2\t\n");
  printf("\n");


  /* LOAD TABULATED EOS */ 

  if(strcmp(eos_type,"tab")==0) 
    load_eos( eos_file, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, 
                                                    &n_tab );
  /* SET UP GRID */

  make_grid(s_gp, mu);

  /* ALLLOCATE MEMORY */

  rho = dmatrix(1,SDIV,1,MDIV);
  gama = dmatrix(1,SDIV,1,MDIV);
  alpha = dmatrix(1,SDIV,1,MDIV);
  omega = dmatrix(1,SDIV,1,MDIV);

  energy = dmatrix(1,SDIV,1,MDIV);
  pressure = dmatrix(1,SDIV,1,MDIV);
  enthalpy = dmatrix(1,SDIV,1,MDIV);
  velocity_sq = dmatrix(1,SDIV,1,MDIV);

  v_plus = dvector(1,SDIV);
  v_minus = dvector(1,SDIV);

  /* set program defaults */
  cf=1.0;
  accuracy=1e-5;    
  xacc = 1e-4;  
 
  if(strcmp(eos_type,"tab")==0) {
    e_surface=7.8*C*C*KSCALE;
    p_surface=1.01e8*KSCALE;
    enthalpy_min=1.0/(C*C);
  }
  else{
    e_surface=0.0;
    p_surface=0.0;
    enthalpy_min=0.0;
  }

  /* CALCULATE THE PRESSURE AND ENTHALPY AT THE CENTRE OF THE STAR*/

  make_center(eos_file, log_e_tab, log_p_tab, 
	      log_h_tab, log_n0_tab, n_tab,eos_type, Gamma_P, 
	      e_center, &p_center, &h_center);

  /* COMPUTE A SPHERICAL STAR AS A FIRST GUESS FOR THE ROTATING STAR */

  sphere( s_gp, log_e_tab, log_p_tab, 
	  log_h_tab, log_n0_tab, n_tab, eos_type, Gamma_P, 
	  e_center, p_center, h_center, p_surface,e_surface,
	  rho, gama, alpha, omega, &r_e);

  r_ratio = 1.0; 

  /* THE PROCEDURE SPIN() WILL COMPUTE THE METRIC OF A STAR WITH
     GIVEN OBLATENESS. THE OBLATENESS IS SPECIFIED BY GIVING 
     THE RATIO OF THE LENGTH OF THE AXIS CONNECTING THE CENTRE OF THE STAR 
     TO ONE OF THE POLES TO THE RADIUS OF THE STAR'S EQUATOR. 
     THIS RATIO IS NAMED r_ratio.
     WHEN r_ratio = 1.0, THE STAR IS SPHERICAL */

  spin(s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, 
       n_tab, eos_type, Gamma_P, 
       h_center, enthalpy_min,
       rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq,
       a_check, accuracy, cf,
       r_ratio, &r_e, &Omega);
  
  /* THE METRIC FUNCTIONS ARE STORED IN THE FUNCTIONS 
     alpha, rho, gama, omega (see user's manual for the definition
     of the metric */


  /* COMPUTE THE VALUES OF VARIOUS EQUILIBRIUM QUANTITIES, SUCH
     AS MASS (Mass), RADIUS (R_e), BARYON MASS(Mass_0), 
     ANGULAR MOMENTUM (J), 
     KEPLERIAN ANGULAR VELOCITY OF PARTICLE ORBITING AT 
     THE EQUATOR,
     VELOCITIES OF CO-ROTATING PARTICLES (v_plus),
     AND COUNTER-ROTATING PARTICLES (v_minus) */

  mass_radius( s_gp, mu, log_e_tab, log_p_tab, 
	       log_h_tab, log_n0_tab, n_tab, eos_type, Gamma_P, 
	       rho, gama, alpha, omega, 
	       energy, pressure, enthalpy, velocity_sq,
	       r_ratio, e_surface, r_e, Omega,
	       &Mass, &Mass_0, &J, &R_e, v_plus, v_minus, &Omega_K);


  /* PRINT OUT INFORMATION ABOUT THE STELLAR MODEL */

  if(strcmp(eos_type,"tab")==0)
    print(r_ratio,e_center, Mass, Mass_0, R_e, Omega, Omega_K, J);
  else 
    printpoly(r_ratio,e_center, Mass, Mass_0, R_e, Omega, Omega_K, J);

  dr=0.05;
 

  /* THIS LOOP STARTS WITH A NON-ROTATING STAR AND INCREASES
     THE STAR'S OBLATENESS (BY DECREASING R_RATIO) AND 
     THEN CALCULATES THE STAR'S ANGULAR VELOCITY. ONCE THE
     COMPUTED VALUE OF ANGULAR VELOCITY IS LARGER THAN 
     THE ANGULAR VELOCITY OF A PARTICLE ORBITING THE STAR
     AT THE EQUATOR, (Omega_K), THE LOOP STOPS */
    
  diff_Omega = Omega_K - Omega;
  old_diff_Omega = diff_Omega;
  
  while( diff_Omega >0){
    /* Find the interval of r_ratio where the star has the
       correct angular velocity	*/
    r_ratio -= dr;
	
    /* Compute the star with the specified value of r_ratio	*/

    spin(s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, 
	 n_tab, eos_type, Gamma_P, 
	 h_center, enthalpy_min,
	 rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq,
	 a_check, accuracy, cf,
	 r_ratio, &r_e, &Omega);
   
    mass_radius( s_gp, mu, log_e_tab, log_p_tab, 
		 log_h_tab, log_n0_tab, n_tab, eos_type, Gamma_P, 
		 rho, gama, alpha, omega, 
		 energy, pressure, enthalpy, velocity_sq,
		 r_ratio, e_surface, r_e, Omega,
		 &Mass, &Mass_0, &J, &R_e, v_plus, v_minus, &Omega_K);

    if(strcmp(eos_type,"tab")==0)
      print(r_ratio,e_center, Mass, Mass_0, R_e, Omega, Omega_K, J);
    else 
      printpoly(r_ratio,e_center, Mass, Mass_0, R_e, Omega, Omega_K, J);
 



   
    old_diff_Omega = diff_Omega;
    diff_Omega = Omega_K - Omega;
   
  } 

  /* The correct star lies between r_ratio and r_ratio + dr */
  xl = r_ratio;
  xh = r_ratio + dr;
  fl = diff_Omega;
  fh = old_diff_Omega;

  /* Use Ridder's method to find the correct star (Taken from 
     Numerical Recipes)	*/


  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    ans=-1.11e30;
    for (j=1;j<=60;j++) {
      xm=0.5*(xl+xh);
      r_ratio = xm;

      spin(s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, 
	   n_tab, eos_type, Gamma_P, 
	   h_center, enthalpy_min,
	   rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq,
	   a_check, accuracy, cf,
	   r_ratio, &r_e, &Omega);

      mass_radius( s_gp, mu, log_e_tab, log_p_tab, 
		   log_h_tab, log_n0_tab, n_tab, eos_type, Gamma_P, 
		   rho, gama, alpha, omega, 
		   energy, pressure, enthalpy, velocity_sq,
		   r_ratio, e_surface, r_e, Omega,
		   &Mass, &Mass_0, &J, &R_e, v_plus, v_minus, &Omega_K);

      fm= Omega_K - Omega;
      sroot=sqrt(fm*fm-fl*fh);
      if (sroot == 0.0) {
	r_ratio = ans;
	break;
      }
       
      xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/sroot);
      if (fabs(xnew-ans) <= xacc) {
	r_ratio = ans;
	break;
      }
      ans=xnew;
      r_ratio = ans;

      spin(s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, 
	   n_tab, eos_type, Gamma_P, 
	   h_center, enthalpy_min,
	   rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq,
	   a_check, accuracy, cf,
	   r_ratio, &r_e, &Omega);

      mass_radius( s_gp, mu, log_e_tab, log_p_tab, 
		   log_h_tab, log_n0_tab, n_tab, eos_type, Gamma_P, 
		   rho, gama, alpha, omega, 
		   energy, pressure, enthalpy, velocity_sq,
		   r_ratio, e_surface, r_e, Omega,
		   &Mass, &Mass_0, &J, &R_e, v_plus, v_minus, &Omega_K);

      if(strcmp(eos_type,"tab")==0)
	print(r_ratio,e_center, Mass, Mass_0, R_e, Omega, Omega_K, J);
      else 
	printpoly(r_ratio,e_center, Mass, Mass_0, R_e, Omega, Omega_K, J);
 

      fnew =  Omega_K - Omega;
      if (fnew == 0.0){
	r_ratio = ans;
	break;
      }
       
      if (SIGN(fm,fnew) != fm) {
	xl=xm;
	fl=fm;
	xh=ans;
	fh=fnew;
      } else if (SIGN(fl,fnew) != fl) {
	xh=ans;
	fh=fnew;
      } else if (SIGN(fh,fnew) != fh) {
	xl=ans;
	fl=fnew;
      } else nrerror("never get here.");
      if (fabs(xh-xl) <= xacc){
	r_ratio = ans;
	break;
      }
    }
  }
  else {
    if (fh == 0.0){
      r_ratio +=dr;
    }
    nrerror("root must be bracketed in zriddr.");
  }
 
  /* THE RIDDER ZERO-FINDING ROUTINE HAS FOUND THE VALUE
     OF R_RATIO WHICH GIVES THE DESIRED STAR. */

  return 0;

}









