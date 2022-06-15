/************************************************************************** 
*                            RNS.C                                        * 
*                                                                         *
*  Relativistic models of rapidly rotating compact stars,                 *
*                                                                         * 
*  using tabulated or polytropic equations of state.                      *
*                                                                         *
*                                                                         *
*  Author:  Nikolaos Stergioulas                                          *
*                                                                         *   
*  Current Address:                                                       *
*                                                                         *
*           Department of Physics                                         *
*           University of Wisconsin-Milwaukee                             *    
*           PO Box 413, Milwaukee, WI 53201, USA                          *  
*                                                                         *
*           E-mail: niksterg@csd.uwm.edu, or                              *
*                   niksterg@pauli.phys.uwm.edu                           *
*                                                                         *
*  Version: 1.1                                                           *  
*                                                                         *
*  Date:    June, 1995                                                    * 
*                                                                         *
*  Usage:   Type rns -h for quick help. Read manual for instructions.     *
*                                                                         *
**************************************************************************/
 
#include <stdio.h>
#include <string.h> 
#include <math.h>

char version[30] = "1.1b, November 1995";


#define SMAX 0.9999                /* maximum value of s-coordinate */  
#define DS (SMAX/(SDIV-1.0))         /* spacing in s-direction */
#define DM (1.0/(MDIV-1.0))          /* spacing in mu direction */ 
#define RDIV 900                     /* grid point in RK integration */ 
#define C 2.9979e10                  /* speed of light in vacuum */
#define G 6.6732e-8                  /* gravitational constant */ 
#define KAPPA 1.346790806509621e+13  /* square of length scale = 1e-15*C*C/G */
#define KSCALE 1.112668301525780e-36 /* KAPPA*G/(C*C*C*C) */  
#define MSUN 1.987e33                /* Mass of Sun */
#define PI 3.1415926535              /* what else */
#define SQ(x) ((x)*(x))              /* square macro */

/*
#define MB 1.659e-24    used in v1.0 
*/

#define MB 1.66e-24

#define RMIN 1.0e-15


char file_name[80] = "no EOS file specified", 
     eos_type[30] = "tab"; 

int n_tab,                           /* number of tabulated EOS points */
    n_nearest=1,                     /* nearest grid point, used in interp. */ 
    a_check=0,                       /* indicates if iteration diverged */ 
    print_dif,                       /* 0 if not print dif */  
    print_option;                    /* select print out */ 


double energy[SDIV+1][MDIV+1],       /* energy density \epsilon */
       mu[MDIV+1],                   /* grid points in mu-direction */
       s_gp[SDIV+1],                 /* grid points in s-direction */
       rho[SDIV+1][MDIV+1],          /* potential \rho */ 
       gama[SDIV+1][MDIV+1],         /* potential \gamma */ 
       omega[SDIV+1][MDIV+1],        /* potential \omega */ 
       alpha[SDIV+1][MDIV+1],        /* potential \alpha */ 
       rho_guess[SDIV+1][MDIV+1],    /* guess for rho */ 
       gama_guess[SDIV+1][MDIV+1],   /* guess for gamma */
       omega_guess[SDIV+1][MDIV+1],  /* guess for alpha */
       alpha_guess[SDIV+1][MDIV+1],  /* guess for omega */
       pressure[SDIV+1][MDIV+1],     /* pressure */ 
       enthalpy[SDIV+1][MDIV+1],     /* enthalpy */
       velocity_sq[SDIV+1][MDIV+1],  /* proper velocity squared */
       gama_mu_1[SDIV+1],            /* gama at \mu=1 */
       gama_mu_0[SDIV+1],            /* gama at \mu=0 */
       rho_mu_1[SDIV+1],             /* rho at \mu=1 */
       rho_mu_0[SDIV+1],             /* rho at \mu=0 */
       omega_mu_0[SDIV+1],           /* omega at \mu=0 */
       r_p,                          /* radius at pole */      
       s_p,                          /* s-coordinate at pole */
       s_e=0.5,                      /* s-coordinate at equator */
       gama_pole_h,                  /* gama^hat at pole */  
       gama_center_h,                /* gama^hat at center */
       gama_equator_h,               /* gama^hat at equator */
       rho_pole_h,                   /* rho^hat at pole */ 
       rho_center_h,                 /* rho^hat at center */
       rho_equator_h,                /* rho^hat at equator */ 
       omega_equator_h,              /* omega^hat at equator */
       Omega_h,                      /* angular velocity \Omega^hat */
       Omega_K,                      /* Kepler frequency */  
       e_center,                     /* central energy density */
       p_center,                     /* central pressure */ 
       h_center,                     /* central enthalpy */
       r_ratio,                      /* ratio of radii r_p/r_e*/ 
       enthalpy_min,                 /* min. enthalpy in h file */
       log_e_tab[201],               /* rho points in tabulated EOS */
       log_p_tab[201],               /* p points in tabulated EOS */
       log_h_tab[201],               /* h points in EOS file */
       log_n0_tab[201],              /* number density in EOS file */  
       f_rho[SDIV+1][LMAX+1][SDIV+1],/* f_rho(s,n,s') */
       f_gama[SDIV+1][LMAX+1][SDIV+1],/* f_gama(s,n,s') */
       f_omega[SDIV+1][LMAX+1][SDIV+1],/* f_omega(s,n,s') */
       P_2n[MDIV+1][LMAX+1],         /* Legendre polynomial P_2n(mu,n) */  
       P1_2n_1[MDIV+1][LMAX+1],      /* ass. Leg. polynomial P^1_2n-1(mu,n) */ 
       r_e,                          /* radius at equator */ 
       r_e_guess,
       R_e,                          /* circumferential radius */
       error,                        /* error in routine ratint (not used) */
       velocity_equator,              
       Mass,
       Mass_0,
       Mass_p,
       J,
       T,
       I,
       W,
       Z_p,
       Z_f,
       Z_b,
       Omega,
       r_final=0.0,                 /* used in guess */
       m_final=0.0,
       r_is_final=0.0,
       r_gp[RDIV+1],
       r_is_gp[RDIV+1],
       m_gp[RDIV+1],
       lambda_gp[RDIV+1],
       e_d_gp[RDIV+1],   
       nu_gp[RDIV+1],
       k_rescale,
       gama_s[SDIV+1],
       rho_s[SDIV+1],
       da_dm_s[MDIV+1],         
       da_dm[SDIV+1][MDIV+1],        /* derivative alpha,mu */
       d_temp,                      /* temporary storage of derivative */
       accuracy,                    /* accuracy in r_e */
       D1_rho[LMAX+1][SDIV+1],  /* integrated term over m in eqn for rho */
       D1_gama[LMAX+1][SDIV+1], /* integrated term over m in eqn for gama */
       D1_omega[LMAX+1][SDIV+1],/* integ. term over m in eqn for omega */
       D2_rho[SDIV+1][LMAX+1],  /* integrated term over s in eqn for rho */
       D2_gama[SDIV+1][LMAX+1], /* integrated term over s in eqn for gama */
       D2_omega[SDIV+1][LMAX+1],/* integ. term over s in eqn for omega */
       S_gama[SDIV+1][MDIV+1],  /* source term in eqn for gama */
       S_rho[SDIV+1][MDIV+1],   /* source term in eqn for rho */
       S_omega[SDIV+1][MDIV+1], /* source term in eqn for omega */
       h_plus,
       h_minus,
       vel_plus,
       vel_minus,
       s_1_s[SDIV+1],
       one_s[SDIV+1],
       one_m2[MDIV+1],
       theta[MDIV+1],
       sin_theta[MDIV+1],
       d_Omega,
       diff_omega,
       sign,
       dr,
       omega_error,
       h_error,
       M_0const,
       J_const,
       M_0_error,
       M_error,
       J_error,
       dgds[SDIV+1][MDIV+1],
       dgdm[SDIV+1][MDIV+1],
       gama_center_old,
       rho_center_old,
       rho_old,
       gama_old,
       omega_old,
       cf,
       Omega_const,
       fix_error,
       p_surface,
       e_surface,
       n_P,
       Gamma_P,
       rho0_center,
       eccentricity,
       grv2,
       grv2_new,
       grv3,
       e_match,
       p_match,
       h_match,
       n0_match,
       de_pt,
       e_cl,
       Omega_plus,
       Omega_p,
       u_phi;

       
/*******************************************************************/
/* Create computational mesh.                                      */
/* Points in the mu-irection are stored in the array mu[i].        */
/* Points in the s-direction are stored in the array s_pg[j].      */
/*******************************************************************/
void make_grid(void)                                    
{ 
  int m,
      s;          /* counter */
   
 for(s=1;s<=SDIV;s++) {
    s_gp[s]=SMAX*(1.0*s-1.0)/(SDIV-1);
    s_1_s[s]=s_gp[s]*(1.0-s_gp[s]);
    one_s[s]=1.0-s_gp[s];
  }

 for(m=1;m<=MDIV;m++) {
    mu[m] = (1.0*m-1.0)/(MDIV-1);
    one_m2[m] = 1.0 - mu[m]*mu[m];
    theta[m] = acos(mu[m]);
    sin_theta[m] = sqrt(one_m2[m]);
 }

}

/*************************************************************************/
/* Load EOS file. Format:  mass_density   pressure   (in CGS units)      */ 
/*************************************************************************/
void load_eos( char eos_file[])
{
 int i;                    /* counter */

 double p,                 /* pressure */
        rho,               /* density */
        h,                 /* enthalpy */
        n0;                /* number density */    
      
 FILE *f_eos;              /* pointer to eos_file */
 

    /* OPEN FILE TO READ */

    if((f_eos=fopen(eos_file,"r")) == NULL ) {    
       printf("cannot open file:  %s\n",eos_file); 
       exit(0);
    }

 
    /* READ NUMBER OF TABULATED POINTS */

    fscanf(f_eos,"%d\n",&n_tab);


    /* READ EOS, H, N0 AND MAKE THEM DIMENSIONLESS (EXCEPT N0) */
 
    for(i=1;i<=n_tab;i++) {  
       fscanf(f_eos,"%lf %lf %lf %lf\n",&rho,&p,&h,&n0) ;
       log_e_tab[i]=log10(rho*C*C*KSCALE);       /* multiply by C^2 to get */ 
       log_p_tab[i]=log10(p*KSCALE);             /* energy density. */
       log_h_tab[i]=log10(h/(C*C));        
       log_n0_tab[i]=log10(n0);                  /* STILL IN CGS ! */
    }
}


/*******************************************************************/
double e_of_rho0(double rho0)
{
 return(pow(rho0,Gamma_P)/(Gamma_P-1.0)+rho0);
}
 
#define MAXIT 30
/*******************************************************************/
double rtsec(double (*func)(double), double x1, double x2, double xacc,double ee)
{
 int j;
 double fl,f,dx,swap, xl,rts;
 
 fl=(*func)(x1)-ee;
 f=(*func)(x2)-ee;

 if(fabs(fl)<fabs(f)) {
   rts=x1;
   xl=x2;
   swap=fl;
   fl=f;
   f=swap;
 } else {
         xl=x1;
         rts=x2;
        }

 
 for(j=1;j<=MAXIT;j++) {
    dx=(xl-rts)*f/(f-fl);
    xl=rts;
    fl=f;
    rts += dx;
    f=(*func)(rts)-ee;
    if(fabs(dx)<xacc||f==0.0) return rts;
  }
 
 printf("Maximum number of iterations exceeded in rtsec");  
 return 0.0;
}

/***************************************************************************/
/* Routine that locates nearest grid point for a given value.              */
/* Adapted from Numerical Recipes.                                         */
/***************************************************************************/
void hunt(double xx[], int n, double x, int *jlo)
{ 
	int jm,jhi,inc,ascnd;

	ascnd=(xx[n] > xx[1]);
	if (*jlo <= 0 || *jlo > n) {
		*jlo=0;
		jhi=n+1;
	} else {
		inc=1;
		if (x >= xx[*jlo] == ascnd) {
			if (*jlo == n) return;
			jhi=(*jlo)+1;
			while (x >= xx[jhi] == ascnd) {
				*jlo=jhi;
				inc += inc;
				jhi=(*jlo)+inc;
				if (jhi > n) {
					jhi=n+1;
					break;
				}
			}
		} else {
			if (*jlo == 1) {
				*jlo=0;
				return;
			}
			jhi=(*jlo);
			*jlo -= 1;
			while (x < xx[*jlo] == ascnd) {
				jhi=(*jlo);
				inc += inc;
				*jlo=jhi-inc;
				if (*jlo < 1) {
					*jlo=0;
					break;
				}
			}
		}
	}
	while (jhi-(*jlo) != 1) {
		jm=(jhi+(*jlo)) >> 1;
		if (x > xx[jm] == ascnd)
			*jlo=jm;
		else
			jhi=jm;
	}
}

/***************************************************************************/
/* Macros used in interp.                                                  */
/***************************************************************************/

static int imaxarg1, imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1)>(imaxarg2) ?\
   (imaxarg1):(imaxarg2))

static int iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1)<(iminarg2) ?\
   (iminarg1):(iminarg2))

/*************************************************************************/
/* Driver for the interpolation routine. First we find the tab. point    */
/* nearest to xb, then we interpolate using four points around xb.       */  
/*************************************************************************/
double interp(double xp[], double yp[], int np ,double xb)
{ 
 int k,        /* index of 1st point */
     m=4;      /* degree of interpolation */ 
 
 double y;     /* intermediate value */

 hunt(xp,np,xb,&n_nearest);

 k=IMIN(IMAX(n_nearest-(m-1)/2,1),np+1-m);

 if( xb==xp[k] ||  xb==xp[k+1] || xb==xp[k+2] || xb==xp[k+3]) xb += 1e-12;

 y= (xb-xp[k+1])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k]/
        ((xp[k]-xp[k+1])*(xp[k]-xp[k+2])*(xp[k]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k+1]/
       ((xp[k+1]-xp[k])*(xp[k+1]-xp[k+2])*(xp[k+1]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+3])*yp[k+2]/
       ((xp[k+2]-xp[k])*(xp[k+2]-xp[k+1])*(xp[k+2]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+2])*yp[k+3]/
       ((xp[k+3]-xp[k])*(xp[k+3]-xp[k+1])*(xp[k+3]-xp[k+2]));

 return (y);
}


/*************************************************************************/
/* Driver for the interpolation routine. Four point interpolation at a   */
/* given offset the index of the first point k.                          */
/*************************************************************************/
double interp_4_k(double xp[], double yp[], int np, double xb, int k)
{  
 double y;     /* intermediate value */

 if( xb==xp[k] ||  xb==xp[k+1] || xb==xp[k+2] || xb==xp[k+3]) xb += 1e-14;

 y= (xb-xp[k+1])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k]/
        ((xp[k]-xp[k+1])*(xp[k]-xp[k+2])*(xp[k]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k+1]/
       ((xp[k+1]-xp[k])*(xp[k+1]-xp[k+2])*(xp[k+1]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+3])*yp[k+2]/
       ((xp[k+2]-xp[k])*(xp[k+2]-xp[k+1])*(xp[k+2]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+2])*yp[k+3]/
       ((xp[k+3]-xp[k])*(xp[k+3]-xp[k+1])*(xp[k+3]-xp[k+2]));

 return (y);
}

/********************************************************************/
/* Integrate f[mu] from m-1 to m.                                   */
/********************************************************************/
double int_z(double f[MDIV+1],int m) 
{
 double x[9];

 x[1]=f[m-1];  
 x[2]=interp(mu,f,MDIV,mu[m-1]+DM/7.0);
 x[3]=interp(mu,f,MDIV,mu[m-1]+2.0*DM/7.0);
 x[4]=interp(mu,f,MDIV,mu[m-1]+3.0*DM/7.0);
 x[5]=interp(mu,f,MDIV,mu[m-1]+4.0*DM/7.0);
 x[6]=interp(mu,f,MDIV,mu[m-1]+5.0*DM/7.0);
 x[7]=interp(mu,f,MDIV,mu[m-1]+6.0*DM/7.0);
 x[8]=f[m];

 return( (DM/17280.0)*(751*x[1]+3577*x[2]+1323*x[3]+2989*x[4]+2989*x[5]
          +1323*x[6]+3577*x[7]+751*x[8]));
}


/*******************************************************************/
double e_at_p(double pp)
{
 if(strcmp(eos_type,"tab")==0) {

   if( strcmp(file_name,"eosCL_LOW")==0 && pp > p_match )
     return pp+e_match+de_pt-p_match;
   else
     return pow(10.0,interp(log_p_tab,log_e_tab,n_tab,log10(pp)));
 
 }else
   return pp/(Gamma_P-1.0) + pow(pp,1.0/Gamma_P); 
}

/*******************************************************************/
double p_at_e(double ee)
{
  if( strcmp(file_name,"eosCL_LOW")==0 && ee > e_match ) {

    if(ee<=e_cl)
      return p_match; 
    else
     return ee-e_match-de_pt+p_match;

  } else
     return pow(10.0,interp(log_e_tab,log_p_tab,n_tab,log10(ee)));
} 

/*******************************************************************/
double p_at_h(double hh)
{ 
  if( strcmp(file_name,"eosCL_LOW")==0 && hh > h_match )
     return 0.5*( (e_match+de_pt+p_match)*exp(2.0*(hh-h_match)) 
            + p_match-e_match-de_pt);
  else 
     return pow(10.0,interp(log_h_tab,log_p_tab,n_tab,log10(hh)));
}

/*******************************************************************/
double h_at_p(double pp)
{ 
  if( strcmp(file_name,"eosCL_LOW")==0 && pp > p_match )
     return h_match + 0.5*log( (2.0*pp+e_match+de_pt-p_match)/
                               (e_match+de_pt+p_match) );
  else 
     return pow(10.0,interp(log_p_tab,log_h_tab,n_tab,log10(pp)));
}

/*******************************************************************/
double n0_at_e(double ee)
{ 
  if( strcmp(file_name,"eosCL_LOW")==0 && ee > e_match ) {

    if(ee<=e_cl)
     return ( (ee+p_match)/(MB*C*C*KSCALE) )*exp(-h_match);
    else         
     return (n0_match+ (de_pt/(MB*C*C*KSCALE))*exp(-h_match) )  
            *sqrt(1.0 + 2.0*(ee-e_match-de_pt)/(e_match+de_pt+p_match));

  }else 
 
 return pow(10.0,interp(log_e_tab,log_n0_tab,n_tab,log10(ee)));
}


/*******************************************************************/
/* Returns the derivative w.r.t. s of an array f[SDIV+1].          */ 
/*******************************************************************/
double s_deriv(double f[SDIV+1],int s)
{   
 switch(s) { 
            case 1    : d_temp=(f[s+1]-f[s])/DS;
                        break;

            case SDIV   : d_temp=(f[s]-f[s-1])/DS;
                          break;
      
            default     : d_temp=(f[s+1]-f[s-1])/(2.0*DS);
                          break; 
 }
 
 return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. mu of an array f[MDIV+1].         */ 
/*******************************************************************/
double m_deriv(double f[MDIV+1], int m)
{
 switch(m) {  
            case 1    : d_temp=(f[m+1]-f[m])/DM;
                        break;

            case MDIV   : d_temp=(f[m]-f[m-1])/DM;
                          break;
      
            default     : d_temp=(f[m+1]-f[m-1])/(2.0*DM);
                          break; 
 }
 
 return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. s                                 */ 
/*******************************************************************/
double deriv_s(double f[SDIV+1][MDIV+1],int s, int m)
{ 
 switch(s) { 
            case 1    : d_temp=(f[s+1][m]-f[s][m])/DS;
                        break;

            case SDIV   : d_temp=(f[s][m]-f[s-1][m])/DS;
                          break;
      
            default     : d_temp=(f[s+1][m]-f[s-1][m])/(2.0*DS);
                          break; 
 } 
 return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. mu                                */ 
/*******************************************************************/
double deriv_m(double f[SDIV+1][MDIV+1],int s, int m)
{
 switch(m) { 
            case 1    : d_temp=(f[s][m+1]-f[s][m])/DM;
                        break;

            case MDIV   : d_temp=(f[s][m]-f[s][m-1])/DM;
                          break;
      
            default     : d_temp=(f[s][m+1]-f[s][m-1])/(2.0*DM);
                          break; 
 } 
 return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. s and mu                          */ 
/*******************************************************************/
double deriv_sm(double f[SDIV+1][MDIV+1],int s, int m)
{
 switch(s) {
     case 1 : if(m==1) {   
               d_temp=(f[s+1][m+1]-f[s][m+1]-f[s+1][m]+f[s][m])/(DM*DS);
              }else{
                if(m==MDIV) {
                 d_temp=(f[s+1][m]-f[s][m]-f[s+1][m-1]+f[s][m-1])/(DM*DS);
                }else{         
                   d_temp=(f[s+1][m+1]-f[s+1][m-1]-f[s][m+1]+f[s][m-1])/
                                                                (2.0*DM*DS);
                }
              }
              break;

     case SDIV : if(m==1) {   
               d_temp=(f[s][m+1]-f[s][m]-f[s-1][m+1]+f[s-1][m])/(DM*DS);
              }else{
                if(m==MDIV) {
                 d_temp=(f[s][m]-f[s-1][m]-f[s][m-1]+f[s-1][m-1])/(DM*DS);
                }else{         
                   d_temp=(f[s][m+1]-f[s][m-1]-f[s-1][m+1]+f[s-1][m-1])/
                                                                (2.0*DM*DS);
                }
             }
             break;
  
     default : if(m==1) {   
               d_temp=(f[s+1][m+1]-f[s-1][m+1]-f[s+1][m]+f[s-1][m])/(2.0*DM*DS);
              }else{
                if(m==MDIV) {
                 d_temp=(f[s+1][m]-f[s-1][m]-f[s+1][m-1]+f[s-1][m-1])/
                                                                (2.0*DM*DS);
                }else{         
                  d_temp=(f[s+1][m+1]-f[s-1][m+1]-f[s+1][m-1]+f[s-1][m-1])/
                                                                (4.0*DM*DS);
                }
             }
             break;
     }

  return d_temp;

}




/*******************************************************************/
/* Returns the Legendre polynomial of degree n, evaluated at x.    */
/*******************************************************************/
double legendre( int n, double x )                      /* checked */
{
  int i;           /* counter */

  double p,        /* Legendre polynomial of order n */
         p_1,      /*    "         "      "    "   n-1*/
         p_2;      /*    "         "      "    "   n-2 */

  p_2=1.0;
  p_1=x;

 if(n>=2) { 
  for(i=2;i<=n;i++){
     p=(x*(2.0*i-1.0)*p_1 - (i-1.0)*p_2)/i;
     p_2=p_1;
     p_1=p;
  }
  return p;
 } else { 
    if (n==1) return p_1;
      else return p_2;
   }
}

/*******************************************************************/
/* Returns the associated Legendre polynomial P_l^m(x).            */
/* Adapted from numerical recipes.                                 */
/*******************************************************************/
double plgndr(int l, int m, double x)
{
	double fact,pll,pmm,pmmp1,somx2;
	int i,ll;

	if (m < 0 || m > l || fabs(x) > 1.0)
		printf("Bad arguments in routine PLGNDR");
	pmm=1.0;
	if (m > 0) {
		somx2=sqrt((1.0-x)*(1.0+x));
		fact=1.0;
		for (i=1;i<=m;i++) {
			pmm *= -fact*somx2;
			fact += 2.0;
		}
	}
	if (l == m)
		return pmm;
	else {
		pmmp1=x*(2*m+1)*pmm;
		if (l == (m+1))
			return pmmp1;
		else {
			for (ll=(m+2);ll<=l;ll++) {
				pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
				pmm=pmmp1;
				pmmp1=pll;
			}
			return pll;
		}
	}
}
 
 
/********************************************************************/
/* Since the grid points are fixed, we can compute the functions    */
/* f_rho, f_gama, f_omega, P_2n, and P1_2n_1 once at the beginning. */
/********************************************************************/
void comp_f_P(void)                                      /* checked */
{ 
  int n,                 /* counter */
      k,                 /* counter for s' */
      j,                 /* counter for s */
      i,                 /* counter for P_2n*/
      s_temp;

  double sj, sk, sj1, sk1;

   if(SMAX==1.0) s_temp=SDIV-1;
     else
         s_temp=SDIV;


   j=1;
 
      n=0; 
      for(k=2;k<=SDIV;k++) {
         sk=s_gp[k];
         f_rho[j][n][k]=1.0/(sk*(1.0-sk));
      }

      n=1;
      for(k=2;k<=SDIV;k++) {
         sk=s_gp[k];
         sk1=1.0-sk;
         
         f_rho[j][n][k]=0.0;
         f_gama[j][n][k]=1.0/(sk*sk1);
         f_omega[j][n][k]=1.0/(sk*sk1);
      }

      for(n=2;n<=LMAX;n++)
         for(k=1;k<=SDIV;k++) {
            f_rho[j][n][k]=0.0;
            f_gama[j][n][k]=0.0;
            f_omega[j][n][k]=0.0;
         }


   k=1;

      n=0;
      for(j=1;j<=SDIV;j++)
         f_rho[j][n][k]=0.0;

      for(j=1;j<=SDIV;j++)
         for(n=1;n<=LMAX;n++) {
            f_rho[j][n][k]=0.0;
            f_gama[j][n][k]=0.0;
            f_omega[j][n][k]=0.0;
         }


   n=0;
   for(j=2;j<=SDIV;j++)
      for(k=2;k<=SDIV;k++) {

         if(SMAX==1.0 && (k==SDIV || j==SDIV))        
           f_rho[j][n][k] = 0.0;
         else {      
               sk=s_gp[k];
               sj=s_gp[j];
               sk1=1.0-sk;
               sj1=1.0-sj;

               if(k<j) 
                 f_rho[j][n][k] = pow(sj1/sj,2.0*n+1.0)*
                                  pow(sk,2.0*n)/pow(sk1,2.0*n+2.0);
               else     
                 f_rho[j][n][k] = pow(sj/sj1,2.0*n)*
                                  pow(sk1,2.0*n-1.0)/pow(sk,2.0*n+1.0);
         }
      }

   for(j=2;j<=SDIV;j++)
      for(n=1;n<=LMAX;n++)
         for(k=2;k<=SDIV;k++) {

            if(SMAX==1.0 && (k==SDIV || j==SDIV)) {
              f_rho[j][n][k] = 0.0;
              f_gama[j][n][k] = 0.0;
              f_omega[j][n][k] = 0.0;
            }
            else{      
                 sk=s_gp[k];
                 sj=s_gp[j];
                 sk1=1.0-sk;
                 sj1=1.0-sj;

                 if(k<j) {   
                          f_rho[j][n][k] = pow(sj1/sj,2.0*n+1.0)*
                                  pow(sk,2.0*n)/pow(sk1,2.0*n+2.0);
                          f_gama[j][n][k] = pow(sj1/sj,2.0*n)*
                                pow(sk,2.0*n-1.0)/pow(sk1,2.0*n+1.0);
                          f_omega[j][n][k] = pow(sj1/sj,2.0*n+1.0)*
                                  pow(sk,2.0*n)/pow(sk1,2.0*n+2.0);
                 }else {     
                        f_rho[j][n][k] = pow(sj/sj1,2.0*n)*
                               pow(sk1,2.0*n-1.0)/pow(sk,2.0*n+1.0);
                        f_gama[j][n][k] = pow(sj/sj1,2.0*n-2.0)*
                                pow(sk1,2.0*n-3.0)/pow(sk,2.0*n-1.0);
                        f_omega[j][n][k] = pow(sj/sj1,2.0*n-2.0)*
                                pow(sk1,2.0*n-3.0)/pow(sk,2.0*n-1.0);
                 }
	    }
	  }


   n=0;
   for(i=1;i<=MDIV;i++)
      P_2n[i][n]=legendre(2*n,mu[i]);

   for(i=1;i<=MDIV;i++)
     for(n=1;n<=LMAX;n++) {
      P_2n[i][n]=legendre(2*n,mu[i]);
      P1_2n_1[i][n] = plgndr(2*n-1 ,1,mu[i]);
    }

}

/***************************************************************/
void make_center(double e_center)
{
 double rho0_center;

 n_nearest=n_tab/2; 

 if(strcmp(eos_type,"tab")==0) {
   p_center=p_at_e(e_center);
   h_center=h_at_p(p_center);
 }
 else {
       rho0_center=rtsec(e_of_rho0,0.0,e_center,1e-16,e_center);
       p_center=pow(rho0_center,Gamma_P);
       h_center=log((e_center+p_center)/rho0_center);
 }
}


/********************************************************************/
/* Compute Omega and Omega_K.                                       */
/********************************************************************/
void comp_omega(void)
{
 int s,
     m;

 double d_o_e[SDIV+1],
        d_g_e[SDIV+1],
        d_r_e[SDIV+1],
        doe,
        dge, 
        dre,
        vek,
        rho_equator,               /* rho at equator */
        omega_equator;             /* omega at equator */

    if(strcmp(eos_type,"tab")==0)
      Omega=Omega_h*C/(r_e*sqrt(KAPPA));
    else
      Omega=Omega_h*C/r_e;

/* Kepler angular velocity */

   for(s=1;s<=SDIV;s++) {
      rho_mu_0[s]=rho[s][1];                     
      omega_mu_0[s]=omega[s][1];                 
   }

   n_nearest= SDIV/2;
   rho_equator=interp(s_gp,rho_mu_0,SDIV,s_e);   

   if(r_ratio==1.0) 
        omega_equator=0.0;
   else 
        omega_equator=interp(s_gp,omega_mu_0,SDIV,s_e);

   for(s=1;s<=SDIV;s++) { 
     d_o_e[s]=deriv_s(omega,s,1);
     d_g_e[s]=deriv_s(gama,s,1);
     d_r_e[s]=deriv_s(rho,s,1);
   }

   doe=interp(s_gp,d_o_e,SDIV,s_e);
   dge=interp(s_gp,d_g_e,SDIV,s_e);
   dre=interp(s_gp,d_r_e,SDIV,s_e);

  vek=(doe/(8.0+dge-dre))*r_e*exp(-rho_equator) + sqrt(((dge+dre)/(8.0+dge
        -dre)) + pow((doe/(8.0+dge-dre))*r_e*exp(-rho_equator),2.0));

  
  if(strcmp(eos_type,"tab")==0)
    Omega_K=(C/sqrt(KAPPA))*(omega_equator+vek*exp(rho_equator)/r_e);
  else 
    Omega_K = omega_equator + vek*exp(rho_equator)/r_e;

 for(s=1;s<=SDIV;s++)
    for(m=1;m<=MDIV;m++) {
       gama_guess[s][m]=gama[s][m];
       rho_guess[s][m]=rho[s][m];
       alpha_guess[s][m]=alpha[s][m];
       omega_guess[s][m]=omega[s][m];
    }

 r_e_guess=r_e;

}
  
/********************************************************************/
/* Compute rest mass and angular momentum.                          */
/********************************************************************/
void comp_M_J(void)
{
 int s,
     m,
     s_temp;

 double D_m_0[SDIV+1],             /* int. quantity for M_0 */ 
        D_m[SDIV+1],               /* int. quantity for M */ 
        D_J[SDIV+1],               /* int. quantity for J */
        rho_0[SDIV+1][MDIV+1],     /* rest mass density */
        rho_mu_0[SDIV+1],
        omega_mu_0[SDIV+1],
        rho_equator,
        omega_equator,                      
        d_o_e[SDIV+1],
        d_g_e[SDIV+1],
        d_r_e[SDIV+1],
        doe,
        dge,
        dre,
        vek;



/* Angular velocity */

    if(r_ratio==1.0) 
        Omega=0.0;
    else {
          if(strcmp(eos_type,"tab")==0)
            Omega=Omega_h*C/(r_e*sqrt(KAPPA));
          else
            Omega=Omega_h*C/r_e;
         }   
/* Kepler angular velocity */

   for(s=1;s<=SDIV;s++) {
      rho_mu_0[s]=rho[s][1];                     
      omega_mu_0[s]=omega[s][1];                 
   }

   n_nearest= SDIV/2;
   rho_equator=interp(s_gp,rho_mu_0,SDIV,s_e);   

   if(r_ratio==1.0)  
        omega_equator=0.0; 
   else 
        omega_equator=interp(s_gp,omega_mu_0,SDIV,s_e);
 
   for(s=1;s<=SDIV;s++) { 
     d_o_e[s]=deriv_s(omega,s,1);
     d_g_e[s]=deriv_s(gama,s,1);
     d_r_e[s]=deriv_s(rho,s,1);
    }

   doe=interp(s_gp,d_o_e,SDIV,s_e);
   dge=interp(s_gp,d_g_e,SDIV,s_e);
   dre=interp(s_gp,d_r_e,SDIV,s_e);
  
  vek=(doe/(8.0+dge-dre))*r_e*exp(-rho_equator) + sqrt(((dge+dre)/(8.0+dge
        -dre)) + pow((doe/(8.0+dge-dre))*r_e*exp(-rho_equator),2.0));

  if(strcmp(eos_type,"tab")==0)
    Omega_K=(C/sqrt(KAPPA))*(omega_equator+vek*exp(rho_equator)/r_e);
  else
    Omega_K=omega_equator+vek*exp(rho_equator)/r_e;

 /* Rest mass and angular momentum */

 Mass_0=0.0;
 Mass=0.0;
 J=0.0;



   if(strcmp(eos_type,"tab")==0) {
     n_nearest=n_tab/2;
     for(s=1;s<=SDIV;s++)
        for(m=1;m<=MDIV;m++) {
           if(energy[s][m]>e_surface)
             rho_0[s][m]=n0_at_e(energy[s][m])*MB*C*C*KSCALE;
           else
             rho_0[s][m]=0.0;
        } 
    }else {
           for(s=1;s<=SDIV;s++)
              for(m=1;m<=MDIV;m++)
                 rho_0[s][m]=(energy[s][m]+pressure[s][m])*exp(-enthalpy[s][m]);
    }


    if(SMAX==1.0) s_temp=SDIV-1;
    else  
        s_temp=SDIV;

    for(s=1;s<=s_temp;s++) {
       D_m[s]=0.0;           /* initialize */
       D_m_0[s]=0.0;
       D_J[s]=0.0;

    for(m=1;m<=MDIV-2;m+=2) {
     D_m[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+gama[s][m])*
              (((energy[s][m]+pressure[s][m])/(1.0-velocity_sq[s][m]))*
              (1.0+velocity_sq[s][m]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m]*mu[m])*r_e*omega[s][m]*
              exp(-rho[s][m])) + 2.0*pressure[s][m])

            + 4.0*exp(2.0*alpha[s][m+1]+gama[s][m+1])*
              (((energy[s][m+1]+pressure[s][m+1])/(1.0-velocity_sq[s][m+1]))*
              (1.0+velocity_sq[s][m+1]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m+1])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+1]*mu[m+1])*r_e*omega[s][m+1]*
              exp(-rho[s][m+1])) + 2.0*pressure[s][m+1]) 

            + exp(2.0*alpha[s][m+2]+gama[s][m+2])*
              (((energy[s][m+2]+pressure[s][m+2])/(1.0-velocity_sq[s][m+2]))*
              (1.0+velocity_sq[s][m+2]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m+2])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+2]*mu[m+2])*r_e*omega[s][m+2]*
              exp(-rho[s][m+2])) + 2.0*pressure[s][m+2]));    

       D_m_0[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+(gama[s][m]
                   -rho[s][m])/2.0)*rho_0[s][m]/sqrt(1.0-velocity_sq[s][m])
 
                   + 4.0* exp(2.0*alpha[s][m+1]+(gama[s][m+1]
                   -rho[s][m+1])/2.0)*rho_0[s][m+1]/sqrt(1.0
                   -velocity_sq[s][m+1])
         
                   + exp(2.0*alpha[s][m+2]+(gama[s][m+2]
                   -rho[s][m+2])/2.0)*rho_0[s][m+2]
                   /sqrt(1.0-velocity_sq[s][m+2])); 
 
        D_J[s] += (1.0/(3.0*(MDIV-1)))*( sin_theta[m]*
                  exp(2.0*alpha[s][m]+gama[s][m]-rho[s][m])*(energy[s][m]
                  +pressure[s][m])*sqrt(velocity_sq[s][m])/(1.0
                  -velocity_sq[s][m])
  
                  +4.0*sqrt(1.0-mu[m+1]*mu[m+1])*
                  exp(2.0*alpha[s][m+1]+gama[s][m+1]-rho[s][m+1])
                  *(energy[s][m+1]+pressure[s][m+1])*sqrt(velocity_sq[s][m+1])/
                  (1.0-velocity_sq[s][m+1])

                  +sqrt(1.0-mu[m+2]*mu[m+2])*
                  exp(2.0*alpha[s][m+2]+gama[s][m+2]-rho[s][m+2])
                  *(energy[s][m+2]+pressure[s][m+2])*sqrt(velocity_sq[s][m+2])/
                  (1.0-velocity_sq[s][m+2]));
    }
 }

   if(SMAX==1.0) {
        D_m[SDIV]=0.0;
        D_m_0[SDIV]=0.0;
        D_J[SDIV]=0.0;
    }

 for(s=1;s<=SDIV-4;s+=2) { 
     Mass += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
          D_m[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_m[s+1]
          +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_m[s+2]);

    Mass_0 += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
              D_m_0[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_m_0[s+1]
              +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_m_0[s+2]);
 
    J += (SMAX/(3.0*(SDIV-1)))*((pow(s_gp[s],3.0)/pow(1.0-s_gp[s],5.0))*
         D_J[s]+ 4.0*(pow(s_gp[s+1],3.0)/pow(1.0-s_gp[s+1],5.0))*
         D_J[s+1] + (pow(s_gp[s+2],3.0)/pow(1.0-s_gp[s+2],5.0))*
         D_J[s+2]);
 }

 if(strcmp(eos_type,"tab")==0) {
   Mass *= 4*PI*sqrt(KAPPA)*C*C*pow(r_e,3.0)/G; 
   Mass_0 *= 4*PI*sqrt(KAPPA)*C*C*pow(r_e,3.0)/G;
   J *= 4*PI*KAPPA*C*C*C*pow(r_e,4.0)/G;
 }
 else {
   Mass *= 4*PI*pow(r_e,3.0); 
   Mass_0 *= 4*PI*pow(r_e,3.0);
   J *= 4*PI*pow(r_e,4.0);
 }

 for(s=1;s<=SDIV;s++)
    for(m=1;m<=MDIV;m++) {
       gama_guess[s][m]=gama[s][m];
       rho_guess[s][m]=rho[s][m];
       alpha_guess[s][m]=alpha[s][m];
       omega_guess[s][m]=omega[s][m];
    }

 r_e_guess=r_e;

}

/********************************************************************/
/* Compute various quantities.                                      */
/********************************************************************/
void comp(void)
{
 int s_surf_1,
     s_temp,
     s,
     m,
     i;                /* counter */
 
 double velocity[SDIV+1][MDIV+1],  /* velocity */
        gama_pole,                 /* gama at pole */  
        gama_equator,              /* gama at equator */
        rho_pole,                  /* rho at pole */ 
        rho_equator,               /* rho at equator */
        omega_equator,             /* omega at equator */
        D_m[SDIV+1],               /* int. quantity for M */
        D_m_0[SDIV+1],             /* int. quantity for M_0 */ 
        D_m_p[SDIV+1],             /* int. quantity for M_p */
        D_J[SDIV+1],               /* int. quantity for J */
        rho_0[SDIV+1][MDIV+1],     /* rest mass density */
        d_o_e[SDIV+1],
        d_g_e[SDIV+1],
        d_r_e[SDIV+1],
        d_v_e[SDIV+1],
        doe,
        dge, 
        dre,
        dve,
        vek,
        d_gama_s,
        d_rho_s,
        d_omega_s,
        d_gama_m,
        d_rho_m,
        d_omega_m,
        d_alpha_s,
        d_alpha_m,
        v_plus[SDIV+1],
        v_minus[SDIV+1],
        sqrt_v,
        s1,
        s_1,
        s_plus,
        s_minus,
        r_plus,
        gama_plus,
        rho_plus,
        gama_minus,
        rho_minus,
        B_st_p[SDIV+1],
        B_st_m[SDIV+1],
        d_v_plus_s,
        d_v_minus_s,
        gama_m[SDIV+1],
        rho_m[SDIV+1],
        alpha_m[SDIV+1],
        enthalpy_m[SDIV+1],
        r_surface[MDIV+1],
        s_surface[MDIV+1],
        dpi_dm,
        d_r_s_dm,
        d_gama_s_dm,
        d_rho_s_dm,
        srt_1_m2,
        f_mu[MDIV+1],
         virial1,
         virial2,
         virial3,
        S_virial1[SDIV+1][MDIV+1],
        S_virial2[SDIV+1][MDIV+1],
        S_virial3[SDIV+1][MDIV+1],
        u_t,
        u_t_plus,
        omega_plus,
        vel_p,
        D_virial1[SDIV+1],
        D_virial2[SDIV+1],
        m1,
        temp1,
        temp2,
        temp3,
        temp4,
        temp5,
        temp6,
        temp_x[5],
        temp_y1[5],
        temp_y2[5],
        S_ad1[SDIV+1][4],
        S_ad2[SDIV+1][4],
        muCh[MDIV+1],
        t_rho[MDIV+1],
        t_alpha[MDIV+1],
        t_rho_s[MDIV+1],
        t_gama_s[MDIV+1],
        t_gama_m[MDIV+1],
        t_rho_m[MDIV+1],
        t_omega_s[MDIV+1],
        t_omega_m[MDIV+1],
        t_pressure[MDIV+1],
        t_energy[MDIV+1],
        t_v2[MDIV+1],
        rhoCh,
        alphaCh,
        rho_sCh,
        gama_sCh,
        gama_mCh,
        rho_mCh,
        omega_sCh,
        omega_mCh,
        pressureCh,
        energyCh,
        v2Ch,
        Sv1Ch[SDIV+1][MDIV+1],
        Sv2Ch[SDIV+1][MDIV+1],
        Dv1Ch[SDIV+1],
        Dv2Ch[SDIV+1],
        gama_surface[MDIV+1],
        rho_surface[MDIV+1],
        alpha_surface[MDIV+1], 
        pi_bar[MDIV+1],
        z_emb[MDIV+1],
        grv2_spherical,
        B_st_p_surface,
        B_st_m_surface,
        B_st_p_out[SDIV-(SDIV/2)+2+1],
        B_st_m_out[SDIV-(SDIV/2)+2+1],
        s_gp_out[SDIV-(SDIV/2)+2+1];



   r_p= r_ratio*r_e;                              /* radius at pole */
   s_p= r_p/(r_p+r_e);                            /* s-coordinate at pole */
   s_e=0.5;

   for(s=1;s<=SDIV;s++)
    for(m=1;m<=MDIV;m++)
     velocity[s][m]=sqrt(velocity_sq[s][m]);
 
   for(s=1;s<=SDIV;s++) {
      gama_mu_1[s]=gama[s][MDIV];                
      gama_mu_0[s]=gama[s][1];                   
      rho_mu_0[s]=rho[s][1];                     
      rho_mu_1[s]=rho[s][MDIV];                  
      omega_mu_0[s]=omega[s][1];                 
   }

   n_nearest= SDIV/2;
   gama_pole=interp(s_gp,gama_mu_1,SDIV,s_p);    
   gama_equator=interp(s_gp,gama_mu_0,SDIV,s_e); 
   rho_pole=interp(s_gp,rho_mu_1,SDIV,s_p);      
   rho_equator=interp(s_gp,rho_mu_0,SDIV,s_e);   

   if(r_ratio==1.0) {
    velocity_equator=0.0;
    omega_equator=0.0;
   }
   else {
    omega_equator=interp(s_gp,omega_mu_0,SDIV,s_e);
    velocity_equator=sqrt(1.0-exp(gama_pole+rho_pole-gama_equator
                                                              -rho_equator));
   }


   /* Circumferential radius */

   if(strcmp(eos_type,"tab")==0)
     R_e = sqrt(KAPPA)*r_e*exp((gama_equator-rho_equator)/2.0);
   else
     R_e = r_e*exp((gama_equator-rho_equator)/2.0);


   /* Masses and angular momentum */

   Mass=0.0;              /* initialize */
   Mass_0=0.0;
   Mass_p=0.0;
   J=0.0;



   if(strcmp(eos_type,"tab")==0) {
     n_nearest=n_tab/2;
     for(s=1;s<=SDIV;s++)
        for(m=1;m<=MDIV;m++) {
           if(energy[s][m]>e_surface)
             rho_0[s][m]=n0_at_e(energy[s][m])*MB*C*C*KSCALE;
           else
             rho_0[s][m]=0.0;
        } 
    }else {
           for(s=1;s<=SDIV;s++)
              for(m=1;m<=MDIV;m++)
                 rho_0[s][m]=(energy[s][m]+pressure[s][m])*exp(-enthalpy[s][m]);
    }


    if(SMAX==1.0) s_temp=SDIV-1;
    else  
        s_temp=SDIV;

    for(s=1;s<=s_temp;s++) {
       D_m[s]=0.0;           /* initialize */
       D_m_0[s]=0.0;
       D_m_p[s]=0.0;
       D_J[s]=0.0;

    for(m=1;m<=MDIV-2;m+=2) {
     D_m[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+gama[s][m])*
              (((energy[s][m]+pressure[s][m])/(1.0-velocity_sq[s][m]))*
              (1.0+velocity_sq[s][m]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m]*mu[m])*r_e*omega[s][m]*
              exp(-rho[s][m])) + 2.0*pressure[s][m])

            + 4.0*exp(2.0*alpha[s][m+1]+gama[s][m+1])*
              (((energy[s][m+1]+pressure[s][m+1])/(1.0-velocity_sq[s][m+1]))*
              (1.0+velocity_sq[s][m+1]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m+1])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+1]*mu[m+1])*r_e*omega[s][m+1]*
              exp(-rho[s][m+1])) + 2.0*pressure[s][m+1]) 

            + exp(2.0*alpha[s][m+2]+gama[s][m+2])*
              (((energy[s][m+2]+pressure[s][m+2])/(1.0-velocity_sq[s][m+2]))*
              (1.0+velocity_sq[s][m+2]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m+2])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+2]*mu[m+2])*r_e*omega[s][m+2]*
              exp(-rho[s][m+2])) + 2.0*pressure[s][m+2]));    
 
     D_m_0[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+(gama[s][m]
              -rho[s][m])/2.0)*rho_0[s][m]/sqrt(1.0-velocity_sq[s][m])

             + 4.0* exp(2.0*alpha[s][m+1]+(gama[s][m+1]
             -rho[s][m+1])/2.0)*rho_0[s][m+1]/sqrt(1.0-velocity_sq[s][m+1])
         
             + exp(2.0*alpha[s][m+2]+(gama[s][m+2]
             -rho[s][m+2])/2.0)*rho_0[s][m+2]/sqrt(1.0-velocity_sq[s][m+2])); 
 
     D_m_p[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+(gama[s][m]
              -rho[s][m])/2.0)*energy[s][m]/sqrt(1.0-velocity_sq[s][m])
 
              + 4.0* exp(2.0*alpha[s][m+1]+(gama[s][m+1]
              -rho[s][m+1])/2.0)*energy[s][m+1]/sqrt(1.0-velocity_sq[s][m+1])
          
             + exp(2.0*alpha[s][m+2]+(gama[s][m+2]
             -rho[s][m+2])/2.0)*energy[s][m+2]/sqrt(1.0-velocity_sq[s][m+2])); 
 
     D_J[s] += (1.0/(3.0*(MDIV-1)))*( sqrt(1.0-mu[m]*mu[m])*
              exp(2.0*alpha[s][m]+gama[s][m]-rho[s][m])*(energy[s][m]
              +pressure[s][m])*sqrt(velocity_sq[s][m])/(1.0-velocity_sq[s][m])
  
              +4.0*sqrt(1.0-mu[m+1]*mu[m+1])*
              exp(2.0*alpha[s][m+1]+gama[s][m+1]-rho[s][m+1])*(energy[s][m+1]
              +pressure[s][m+1])*sqrt(velocity_sq[s][m+1])/
              (1.0-velocity_sq[s][m+1])

              + sqrt(1.0-mu[m+2]*mu[m+2])*
              exp(2.0*alpha[s][m+2]+gama[s][m+2]-rho[s][m+2])*(energy[s][m+2]
              +pressure[s][m+2])*sqrt(velocity_sq[s][m+2])/
              (1.0-velocity_sq[s][m+2]));
 
    }
   }

   if(SMAX==1.0) {
        D_m[SDIV]=0.0;
        D_m_0[SDIV]=0.0;
        D_m_p[SDIV]=0.0;
        D_J[SDIV]=0.0;
    }


    for(s=1;s<=SDIV-4;s+=2) { 

     Mass += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
          D_m[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_m[s+1]
          +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_m[s+2]);
 
     Mass_0 += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
          D_m_0[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_m_0[s+1]
          +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_m_0[s+2]);
 
     Mass_p += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
          D_m_p[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_m_p[s+1]
          +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_m_p[s+2]);

     J += (SMAX/(3.0*(SDIV-1)))*((pow(s_gp[s],3.0)/pow(1.0-s_gp[s],5.0))*
          D_J[s]+ 4.0*(pow(s_gp[s+1],3.0)/pow(1.0-s_gp[s+1],5.0))*
          D_J[s+1] + (pow(s_gp[s+2],3.0)/pow(1.0-s_gp[s+2],5.0))*
          D_J[s+2]);
    }
   
    if(strcmp(eos_type,"tab")==0) {
      Mass *= 4*PI*sqrt(KAPPA)*C*C*pow(r_e,3.0)/G;
      Mass_0 *= 4*PI*sqrt(KAPPA)*C*C*pow(r_e,3.0)/G;
      Mass_p *= 4*PI*sqrt(KAPPA)*C*C*pow(r_e,3.0)/G;
    }
    else {
          Mass *= 4*PI*pow(r_e,3.0);
          Mass_0 *= 4*PI*pow(r_e,3.0);
          Mass_p *= 4*PI*pow(r_e,3.0);
    }
 
    if(r_ratio==1.0) {
       J=0.0; 
       Omega=0.0;
    }
    else {    
          if(strcmp(eos_type,"tab")==0) {
            Omega=Omega_h*C/(r_e*sqrt(KAPPA));
            J *= 4*PI*KAPPA*C*C*C*pow(r_e,4.0)/G;
          }
          else {
                Omega=Omega_h/r_e;
                J *= 4*PI*pow(r_e,4.0);
	  }
    }

    T=0.5*J*Omega;

    if(r_ratio==1.0) I=0.0;
    else
     I=J/Omega;
 

    if(strcmp(eos_type,"tab")==0)     
      W = Mass_p*C*C - Mass*C*C + T;
    else
      W = Mass_p - Mass + T;


   /* Redshifts */

   Z_p=exp(-0.5*(gama_pole+rho_pole))-1.0;

   Z_b=sqrt((1.0+velocity_equator)/(1.0-velocity_equator))*(exp(-0.5*(
      gama_equator+rho_equator))/(1.0-omega_equator*r_e
       *exp(-rho_equator)))-1.0;

   Z_f= sqrt((1.0-velocity_equator)/(1.0+velocity_equator))*(exp(-0.5*(
      gama_equator+rho_equator))/(1.0+omega_equator*r_e
       *exp(-rho_equator)))-1.0;


   /* Kepler angular velocity */

   for(s=1;s<=SDIV;s++) { 
     d_o_e[s]=deriv_s(omega,s,1);
     d_g_e[s]=deriv_s(gama,s,1);
     d_r_e[s]=deriv_s(rho,s,1);
     d_v_e[s]=deriv_s(velocity,s,1);
   }

   doe=interp(s_gp,d_o_e,SDIV,s_e);
   dge=interp(s_gp,d_g_e,SDIV,s_e);
   dre=interp(s_gp,d_r_e,SDIV,s_e);
   dve=interp(s_gp,d_v_e,SDIV,s_e);

  vek=(doe/(8.0+dge-dre))*r_e*exp(-rho_equator) + sqrt(((dge+dre)/(8.0+dge
        -dre)) + pow((doe/(8.0+dge-dre))*r_e*exp(-rho_equator),2.0));

  if(strcmp(eos_type,"tab")==0) 
    Omega_K=(C/sqrt(KAPPA))*(omega_equator+vek*exp(rho_equator)/r_e);
  else
    Omega_K= omega_equator + vek*exp(rho_equator)/r_e;



  /* Embedding */
 
  for(m=1;m<=MDIV;m++) { 

    /* Find last s grid point inside star. Construct arrays for the metric
     | potentials and the enthalpy along each direction \mu.
    */

    for(s=1;s<=SDIV;s++) {
       if(energy[s][m]>0) s_surf_1=s;       
       gama_m[s]=gama[s][m];                 
       rho_m[s]=rho[s][m];                    
       alpha_m[s]=alpha[s][m];         
       enthalpy_m[s]=enthalpy[s][m];
    }
 
    /* If the positive enthalpy region outside the star is at a distance  
     | greater than two gridpoints, interpolate using two points on either
     | side of the surface, else, if only one gridpoint has negative      
     | enthalpy, use linear interpolation, else use s_surf_1 as the       
     | surface. 
    */
        
    if(enthalpy_m[s_surf_1+2]<0.0) 
          s_surface[m]=interp_4_k(enthalpy_m,s_gp,SDIV,0.0,s_surf_1-1);
    else {     
          if(enthalpy_m[s_surf_1+1]<0.0) 
              s_surface[m]=s_gp[s_surf_1]-DS*enthalpy_m[s_surf_1]/
                        (enthalpy_m[s_surf_1+1]-enthalpy_m[s_surf_1]);
          else 
              s_surface[m]=s_gp[s_surf_1];
    }

    if(m==1)
      s_surface[m]=0.5;

    /* Construct arrays for the coordinate radius and the metric potentials
     | on the surface.
    */
 
    r_surface[m]=r_e*(s_surface[m]/(1.0-s_surface[m]));

    gama_surface[m]=interp(s_gp,gama_m,SDIV,s_surface[m]);
    rho_surface[m]=interp(s_gp,rho_m,SDIV,s_surface[m]);
    alpha_surface[m]=interp(s_gp,alpha_m,SDIV,s_surface[m]);
 }

 for(m=1;m<=MDIV;m++) {
    d_gama_s_dm=m_deriv(gama_surface,m);
    d_rho_s_dm=m_deriv(rho_surface,m);
    d_r_s_dm=m_deriv(r_surface,m);

    if(m==MDIV) 
          f_mu[m]=0.0;
    else {
          dpi_dm=exp(0.5*(gama_surface[m]-rho_surface[m]))*(0.5*r_surface[m]*
                (d_gama_s_dm-d_rho_s_dm)*sin_theta[m] + d_r_s_dm*sin_theta[m] -
                r_surface[m]*mu[m]/sin_theta[m]);
          f_mu[m]=sqrt(exp(2.0*alpha_surface[m])*(pow(d_r_s_dm,2.0)
                +pow(r_surface[m],2.0)/(1.0-mu[m]*mu[m])) - pow(dpi_dm,2.0));
         }
 }

   pi_bar[1]=R_e;
   z_emb[1]=0.0;                         /* integrate f_mu */
   for(m=2;m<=MDIV;m++) {                
      z_emb[m] = z_emb[m-1] + sqrt(KAPPA)*int_z(f_mu,m);    
      pi_bar[m]=sqrt(KAPPA)*exp((gama_surface[m]-rho_surface[m])/2.0)
                *r_surface[m]*sin_theta[m];
      if(pi_bar[m]>pi_bar[m-1] && m>=2)
 
        pi_bar[m]=pi_bar[m-1];

   }  

   eccentricity=sqrt(1.0-SQ(z_emb[MDIV]/R_e));



   /* Last stable circular orbit */


    for(s=1;s<=SDIV-1;s++) {
    s1= s_gp[s]*(1.0-s_gp[s]);
    s_1=1.0-s_gp[s];
        
    d_gama_s=deriv_s(gama,s,1);
    d_rho_s=deriv_s(rho,s,1);
    d_omega_s=deriv_s(omega,s,1);

    sqrt_v= exp(-2.0*rho[s][1])*r_e*r_e*pow(s_gp[s],4.0)*pow(d_omega_s,2.0) 
            + 2*s1*(d_gama_s+d_rho_s)+s1*s1*(d_gama_s*d_gama_s-d_rho_s*d_rho_s);

    if(sqrt_v>0.0) sqrt_v= sqrt(sqrt_v);
     else {
      sqrt_v=0.0;
      if(s_gp[s]>=s_e) printf("velocity imaginary at s=%3.2e\n",s_gp[s]);
     }

    v_plus[s]=(exp(-rho[s][1])*r_e*s_gp[s]*s_gp[s]*d_omega_s + sqrt_v)/
              (2.0+s1*(d_gama_s-d_rho_s));

    v_minus[s]=(exp(-rho[s][1])*r_e*s_gp[s]*s_gp[s]*d_omega_s - sqrt_v)/
               (2.0+s1*(d_gama_s-d_rho_s));
  }

    v_plus[SDIV]=0.0;
    v_minus[SDIV]=0.0;




  for(s=1;s<=SDIV;s++) {
    s1= s_gp[s]*(1.0-s_gp[s]);
        
    d_gama_s=deriv_s(gama,s,1);
    d_rho_s=deriv_s(rho,s,1);
    d_omega_s=deriv_s(omega,s,1);

    d_v_plus_s=s_deriv(v_plus,s);
    d_v_minus_s=s_deriv(v_minus,s);
 
    B_st_p[s]=v_plus[s]*(r_e*s_gp[s]*s_gp[s]*d_omega_s*exp(-rho[s][1])
 
              +s1*d_v_plus_s) + 0.5*s1*(d_gama_s+d_rho_s)

              - pow(v_plus[s],4.0)*(1.0+0.5*s1*(d_gama_s-d_rho_s));
 
    B_st_m[s]=v_minus[s]*(r_e*s_gp[s]*s_gp[s]*d_omega_s*exp(-rho[s][1])
 
              +s1*d_v_minus_s) + 0.5*s1*(d_gama_s+d_rho_s)

              - pow(v_minus[s],4.0)*(1.0+0.5*s1*(d_gama_s-d_rho_s));
  }  


  B_st_p_surface=interp(s_gp,B_st_p,SDIV,s_surface[1]);
  B_st_m_surface=interp(s_gp,B_st_m,SDIV,s_surface[1]);

  if(B_st_p_surface>0.0) {
    h_plus= 0.0;
    vel_plus= vek;
    Omega_plus= Omega_K;
  }
  else{
       for(i=1;i<=SDIV-(SDIV/2)+2;i++) {
          B_st_p_out[i]=B_st_p[(SDIV/2)-2+i];
          s_gp_out[i]=s_gp[(SDIV/2)-2+i];
       }

       n_nearest=SDIV/4;
       s_plus=interp(B_st_p_out, s_gp_out, SDIV-(SDIV/2)+2, 0.0);
       r_plus = r_e*s_plus/(1.0-s_plus);
       gama_plus=interp(s_gp,gama_mu_0,SDIV,s_plus);
       rho_plus=interp(s_gp,rho_mu_0,SDIV,s_plus);
       omega_plus=interp(s_gp,omega_mu_0,SDIV,s_plus);
       vel_plus=interp(s_gp,v_plus,SDIV,s_plus); 
 
       if(strcmp(eos_type,"tab")==0) {
         h_plus=sqrt(KAPPA)*(r_plus*exp((gama_plus-rho_plus)/2.0)
                               -r_e*exp((gama_equator-rho_equator)/2.0));
         Omega_plus=(C/sqrt(KAPPA))*(omega_plus+vel_plus*exp(rho_plus)/r_e);
       }
       else {
             h_plus=(r_plus*exp((gama_plus-rho_plus)/2.0)
                              -r_e*exp((gama_equator-rho_equator)/2.0));
             Omega_plus = (omega_plus+vel_plus*exp(rho_plus)/r_e);
       }

  }


  if(B_st_m_surface>0.0) {
    h_minus= 0.0;
    vel_minus= vek;
  }
  else{
       for(i=1;i<=SDIV-(SDIV/2)+2;i++) {
          B_st_m_out[i]=B_st_m[(SDIV/2)-2+i];
          s_gp_out[i]=s_gp[(SDIV/2)-2+i];
       }

       n_nearest=SDIV/4;
       s_minus=interp(B_st_m_out, s_gp_out, SDIV-(SDIV/2)+2, 0.0);
       gama_minus=interp(s_gp,gama_mu_0,SDIV,s_minus);
       rho_minus=interp(s_gp,rho_mu_0,SDIV,s_minus);
       vel_minus=interp(s_gp,v_plus,SDIV,s_minus); 
 
       if(strcmp(eos_type,"tab")==0) 
         h_minus=sqrt(KAPPA)*r_e*((s_minus/(1.0-s_minus))*exp((gama_minus
                 -rho_minus)/2.0) - exp((gama_equator-rho_equator)/2.0));      
       else 
             h_minus=r_e*((s_minus/(1.0-s_minus))*exp((gama_minus
                     -rho_minus)/2.0) - exp((gama_equator-rho_equator)/2.0));

  }

  if(h_plus!= 0.0) {
  
    u_phi= vel_plus*r_plus*exp((gama_plus-rho_plus)/2.0)/
                                                    sqrt(1.0 - SQ(vel_plus));
    vel_p = u_phi/sqrt(SQ(u_phi)+SQ(r_e)*exp(gama_equator-rho_equator));

    Omega_p = omega_equator + (vel_p/r_e)*exp(rho_equator);

    if(strcmp(eos_type,"tab")==0) 
         Omega_p *= C/sqrt(KAPPA);
  } 
  else {
        u_phi= vek*r_e*exp((gama_equator-rho_equator)/2.0)/
                                                sqrt(1.0 - SQ(vek));
        Omega_p= Omega_K;
  }  

  if(strcmp(eos_type,"tab")==0) {
    u_phi *= C*C*sqrt(KAPPA)/(G*MSUN);
  }
 



   /* Virial theorem */


   /* GRV2 spherical */

   virial1=0.0;
   virial2=0.0; 

   if(r_ratio==1.0) {


   if(SMAX==1.0) s_temp=SDIV-1;
     else  
         s_temp=SDIV;

   for(s=1;s<=s_temp;s++) {
      d_gama_s=deriv_s(gama,s,1);
      d_rho_s=deriv_s(rho,s,1);
 
      S_virial1[s][1] = 8*SQ(PI*r_e)*s_gp[s]*pressure[s][1]
                        *exp(2.0*alpha[s][1])/pow(1.0-s_gp[s],3.0);

      S_virial2[s][1] = PI*s_1_s[s]*SQ(d_gama_s+d_rho_s)/4.0;
   }
 
 
   if(SMAX==1.0) { 
     S_virial1[SDIV][1] =0.0;
 
     d_gama_s=deriv_s(gama,SDIV,1);
     d_rho_s=deriv_s(rho,SDIV,1);
     S_virial2[SDIV][1] =PI*s_1_s[s]*SQ(d_gama_s+d_rho_s)/4.0;
   }

   for(s=1;s<=SDIV-2;s+=2) {
       virial1 += (DS/3.0)*(S_virial1[s][1] + 4.0*S_virial1[s+1][1]
                                            + S_virial1[s+2][1]);
       virial2 += (DS/3.0)*(S_virial2[s][1] + 4.0*S_virial2[s+1][1]
                                            + S_virial2[s+2][1]);
   }


   grv2_spherical=fabs(1.0-virial1/virial2); 
/*
   printf("\n");
   printf("virial1=%6.5e, virial2=%6.5e e\n",virial1,virial2);  
   printf("grv2_spherical=%6.5e\n",grv2_spherical); 
*/
   }
 

   /* GRV2 */

   virial1=0.0;
   virial2=0.0; 
   virial3=0.0;

   for(i=1;i<=4;i++)
      temp_x[i]=mu[MDIV+1-i];

   temp6=cos(0.5*(theta[MDIV]+theta[MDIV-2]));


   /* Set equal to zero at center */

   for(m=1;m<=MDIV;m++) {
      S_virial1[1][m];
      S_virial2[1][m];
   }

   for(m=1;m<=3;m++) {   
      S_ad1[1][m];
      S_ad2[1][m];
   }


   for(s=2;s<=SDIV;s++) {
      for(m=1;m<=MDIV;m++) {
         s1=s_1_s[s];
         m1=one_m2[m]; 
         d_gama_s=deriv_s(gama,s,m);
         d_rho_s=deriv_s(rho,s,m);
         d_gama_m=deriv_m(gama,s,m);
         d_rho_m=deriv_m(rho,s,m);
         d_omega_s=deriv_s(omega,s,m);
         d_omega_m=deriv_m(omega,s,m);
         if(m==MDIV || m==MDIV-1) { 
           S_virial1[s][m]=0.0;
           S_virial3[s][m]=0.0;
         }
         else {
               if(SMAX==1 && s==SDIV) S_virial1[s][m]=0.0;
               else 
                   S_virial1[s][m] = 16*PI*SQ(r_e)*(s_gp[s]/pow(1.0-s_gp[s],3.0))
                                     *( pressure[s][m]+(energy[s][m]+
                                     pressure[s][m])
                                     *(velocity_sq[s][m]/(1.0-velocity_sq[s][m])))
                                     *exp(2.0*alpha[s][m])/sqrt(m1);

               S_virial3[s][m] = (s1/(2.0*sqrt(m1)))*SQ(d_gama_s+d_rho_s);
	 }

         if(SMAX==1.0 && s==SDIV) S_virial2[s][m]=0.0;
         else 
             S_virial2[s][m] = 2.0*( sqrt(m1)*SQ(d_gama_m+d_rho_m)/s1 
                               -(3.0*SQ(r_e)
                               *s_gp[s]/(4.0*pow(1.0-s_gp[s],3.0)))*sqrt(m1)
                               *exp(-2.0*rho[s][m])
                               *( SQ(s1*d_omega_s)+m1*SQ(d_omega_m) ) );
       }

   if(SMAX==1.0 && s==SDIV) {
      for(m=1;m<=4;m++) { 
         if(m!=2 || m!=4) {
           S_ad1[s][m]=0.0;
           S_ad2[s][m]=0.0;
         }
      }
      S_ad1[s][2]=0.0;
      S_ad2[s][2]=0.0;
   }
   else {

      for(m=1;m<=4;m++) { 
         temp_y1[m]= 16*PI*SQ(r_e)*(s_gp[s]/pow(1.0-s_gp[s],3.0))*( 
                     pressure[s][MDIV+1-m]+(energy[s][MDIV+1-m]
                     +pressure[s][MDIV+1-m])*(velocity_sq[s][MDIV+1-m]/
                     (1.0-velocity_sq[s][MDIV+1-m])))*exp(2.0*alpha[s][MDIV+1-m]);
 
         temp_y2[m]= 0.5*s1*SQ(deriv_s(gama,s,MDIV+1-m)+deriv_s(rho,s,MDIV+1-m));

         if(m!=2 || m!=4) {
           S_ad1[s][m]=temp_y1[m];
           S_ad2[s][m]=temp_y2[m];
         }            
      }  

      temp1=interp(temp_x,temp_y1,4,temp6);
      temp2=interp(temp_x,temp_y2,4,temp6);

      S_ad1[s][2]=temp1;
      S_ad2[s][2]=temp2;
    }

  }


   for(s=1;s<=SDIV;s++) {
      for(m=1;m<=MDIV-4;m+=2) {
         virial1 += (DM/3.0)*(S_virial1[s][m]+4.0*S_virial1[s][m+1]
                                  +S_virial1[s][m+2]);
         virial2 += (DM/3.0)*(S_virial2[s][m]+4.0*S_virial2[s][m+1]
                                  +S_virial2[s][m+2]);    
         virial3 += (DM/3.0)*(S_virial3[s][m]+4.0*S_virial3[s][m+1]
                                  +S_virial3[s][m+2]);    

      }

      virial1 += ((theta[MDIV-2]-theta[MDIV])/6.0)*(S_ad1[s][1]
                            +4.0*S_ad1[s][2]+S_ad1[s][3]);
      virial2 += (DM/3.0)*(S_virial2[s][MDIV-2]+4.0*S_virial2[s][MDIV-1]
                                  +S_virial2[s][MDIV]);    
      virial3 += ((theta[MDIV-2]-theta[MDIV])/6.0)*(S_ad2[s][1]
                            +4.0*S_ad2[s][2]+S_ad2[s][3]);

      virial2 += virial3;

      D_virial1[s]=virial1;
      D_virial2[s]=virial2;
      virial1=0.0;
      virial2=0.0;
      virial3=0.0;
   }

   for(s=1;s<=SDIV-2;s+=2) {
      virial1 += (DS/3.0)*(D_virial1[s]+4.0*D_virial1[s+1]
                                  +D_virial1[s+2]);
      virial2 += (DS/3.0)*(D_virial2[s]+4.0*D_virial2[s+1]
                                  +D_virial2[s+2]);
   }

   grv2=fabs(1.0-virial1/virial2); 
/*
   printf("\n");
   printf("virial1=%6.5e, virial2=%6.5e\n",virial1,virial2);  
   printf("grv2=%6.5e\n",grv2);
*/


   /* GRV2 GAUSS-CHEBYSHEV */

   virial1=0.0;
   virial2=0.0;

   for(i=1;i<=MDIV;i++) muCh[i]=cos((1.0*i-0.5)/((2.0*MDIV-1.0))*PI);

   for(s=1;s<=s_temp;s++) {
      for(m=1;m<=MDIV;m++) {
         t_rho[m]=rho[s][m];
         t_alpha[m]=alpha[s][m];
         t_rho_s[m]=deriv_s(rho,s,m);
         t_gama_s[m]=deriv_s(gama,s,m);
         t_gama_m[m]=deriv_m(gama,s,m);
         t_rho_m[m]=deriv_m(rho,s,m);
         t_omega_s[m]=deriv_s(omega,s,m);
         t_omega_m[m]=deriv_m(omega,s,m);
         t_pressure[m]=pressure[s][m];
         t_energy[m]=energy[s][m];
         t_v2[m]=velocity_sq[s][m];
      }

      for(m=1;m<=MDIV;m++) {
         rhoCh=interp(mu,t_rho,MDIV,muCh[m]);
         alphaCh=interp(mu,t_alpha,MDIV,muCh[m]);
         rho_sCh=interp(mu,t_rho_s,MDIV,muCh[m]);
         gama_sCh=interp(mu,t_gama_s,MDIV,muCh[m]);
         gama_mCh=interp(mu,t_gama_m,MDIV,muCh[m]);
         rho_mCh=interp(mu,t_rho_m,MDIV,muCh[m]);
         omega_sCh=interp(mu,t_omega_s,MDIV,muCh[m]);
         omega_mCh=interp(mu,t_omega_m,MDIV,muCh[m]);
         pressureCh=interp(mu,t_pressure,MDIV,muCh[m]);
         energyCh=interp(mu,t_energy,MDIV,muCh[m]);
         v2Ch=interp(mu,t_v2,MDIV,muCh[m]);
         s1=s_1_s[s];
         m1=1.0-SQ(muCh[m]); 

         if(s==1) temp1=0.0;
           else
               temp1=1.0/s1;

         Sv1Ch[s][m]= 8*PI*SQ(r_e)*(s_gp[s]/pow(1.0-s_gp[s],3.0))
                      *( pressureCh+(energyCh+pressureCh)*v2Ch/(1.0-v2Ch))
                      *exp(2.0*alphaCh);
       
         Sv2Ch[s][m]= s1*SQ(gama_sCh+rho_sCh)/4.0+m1*SQ(gama_mCh+rho_mCh)
                      *temp1-(3.0/4.0)*SQ(r_e)*(s_gp[s]/
                      pow(1.0-s_gp[s],3.0))*m1
                      *exp(-2.0*rhoCh)*(SQ(s1*omega_sCh)+m1*SQ(omega_mCh));	 
      }


      if(SMAX==1.0) {
        for(m=1;m<=MDIV;m++) {    
           Sv1Ch[SDIV][m]=0.0;
           Sv2Ch[SDIV][m]=0.0;
        }
      }

      Dv1Ch[s]=0.0;
      Dv2Ch[s]=0.0;

      for(m=1;m<=MDIV;m++) {
         Dv1Ch[s] += (2.0*PI/(2.0*MDIV-1.0))*Sv1Ch[s][m];
         Dv2Ch[s] += (2.0*PI/(2.0*MDIV-1.0))*Sv2Ch[s][m];
      }
   }

   for(s=1;s<=SDIV-2;s+=2) {
      virial1 += (DS/3.0)*(Dv1Ch[s]+4.0*Dv1Ch[s+1]
                                  +Dv1Ch[s+2]);
      virial2 += (DS/3.0)*(Dv2Ch[s]+4.0*Dv2Ch[s+1]
                                  +Dv2Ch[s+2]);
   }

   grv2_new=fabs(1.0-virial1/virial2); 

/*
   printf("\n");
   printf("virial1=%6.5e, virial2=%6.5e\n",virial1,virial2);  
   printf("grv2_new=%6.5e\n",grv2_new);
*/


   /* GRV2  BY PARTS */

   virial1=0.0;
   virial2=0.0; 

   /* Set equal to zero at center */
/*
   for(m=1;m<=MDIV;m++) {
      S_virial1[1][m];
      S_virial2[1][m];
      Sv1Ch[1][m];
      Sv2Ch[1][m];
   }

   for(s=2;s<=SDIV;s++) {
      for(m=1;m<=MDIV;m++) {
         s1=s_gp[s]*(1.0-s_gp[s]);
         m1=1.0-SQ(mu[m]); 
         d_gama_s=deriv_s(gama,s,m);
         d_rho_s=deriv_s(rho,s,m);
         d_gama_m=deriv_m(gama,s,m);
         d_rho_m=deriv_m(rho,s,m);
         d_omega_s=deriv_s(omega,s,m);
         d_omega_m=deriv_m(omega,s,m);
 
         S_virial1[s][m] = SQ(r_e)*(s_gp[s]/pow(1.0-s_gp[s],3.0))
                           *( pressure[s][m] + (energy[s][m]+pressure[s][m])
                           *(velocity_sq[s][m]/(1.0-velocity_sq[s][m])))
                           *exp(2.0*alpha[s][m]);
 
         S_virial2[s][m] = (s1/4.0)*SQ(d_gama_s+d_rho_s)+
                            m1*SQ(d_gama_m+d_rho_m)/s1 - (3.0*SQ(r_e)*s_gp[s]/
                           (4.0*pow(1.0-s_gp[s],3.0)))*m1*exp(-2.0*rho[s][m])
                           *( SQ(s1*d_omega_s)+m1*SQ(d_omega_m) ) ;
      }
   } 

   for(s=2;s<=SDIV;s++) 
      for(m=1;m<=MDIV;m++) {
         Sv1Ch[s][m]=deriv_m(S_virial1,s,m);
         Sv2Ch[s][m]=deriv_m(S_virial2,s,m);
      }

   D_virial1[1]=0.0;
   D_virial2[1]=0.0;

   for(s=2;s<=SDIV;s++) {
     s1=s_gp[s]*(1.0-s_gp[s]); 
     d_gama_s=deriv_s(gama,s,MDIV);
     d_rho_s=deriv_s(rho,s,MDIV);

     for(m=1;m<=MDIV-4;m+=2) {
         virial1 += (DM/3.0)*(asin(mu[m])*Sv1Ch[s][m]
                              +4.0*asin(mu[m+1])*Sv1Ch[s][m+1]
                              +asin(mu[m+2])*Sv1Ch[s][m+2]);

         virial2 += (DM/3.0)*(asin(mu[m])*Sv2Ch[s][m]
                              +4.0*asin(mu[m+1])*Sv2Ch[s][m+1]
                              +asin(mu[m+2])*Sv2Ch[s][m+2]);
      }

      virial1 *= (-1.0);
      virial1 += (PI/2.0)*SQ(r_e)*(s_gp[s]/pow(1.0-s_gp[s],3.0))
                 *pressure[s][MDIV]*exp(2.0*alpha[s][MDIV]);       

      virial2 *= (-1.0);
      virial2 += (PI/8.0)*s1*SQ(d_gama_s+d_rho_s);              
              
      D_virial1[s]=virial1;
      D_virial2[s]=virial2;
      virial1=0.0;
      virial2=0.0;
   }

   for(s=1;s<=SDIV-2;s+=2) {
      virial1 += (DS/3.0)*(D_virial1[s]+4.0*D_virial1[s+1]
                                  +D_virial1[s+2]);
      virial2 += (DS/3.0)*(D_virial2[s]+4.0*D_virial2[s+1]
                                  +D_virial2[s+2]);
   }

   virial1 *= 16.0*PI;
   virial2 *= 2.0;

   grv2_new=fabs(1.0-virial1/virial2); 

   printf("\n");
   printf("virial1=%6.5e, virial2=%6.5e\n",virial1,virial2);  
   printf("grv2 b.p. =%6.5e\n",grv2_new);
*/


   /* GRV3 */

   virial1=0.0;
   virial2=0.0;
 
   for(s=1;s<=SDIV;s++)
      for(m=1;m<=MDIV;m++) {
         s1=s_1_s[s];
         m1=one_m2[m];
         d_gama_s=deriv_s(gama,s,m);
         d_gama_m=deriv_m(gama,s,m);
         d_rho_s=deriv_s(rho,s,m);
         d_rho_m=deriv_m(rho,s,m);
         d_omega_s=deriv_s(omega,s,m);
         d_omega_m=deriv_m(omega,s,m);
         d_alpha_s=deriv_s(alpha,s,m);
         d_alpha_m=deriv_m(alpha,s,m);


         if(SMAX==1.0 && s==SDIV) {
           S_virial1[s][m] = 0.0;

/* CORRECT THIS ! */

           S_virial2[s][m] = 0.5*( SQ(s_gp[s]*(d_gama_s+d_rho_s))- d_alpha_s*
                              ( SQ(s_gp[s])*(d_gama_s-d_rho_s)))*
                              exp(0.5*(gama[s][m]-rho[s][m]))*r_e;
         }
         else{
  
         S_virial1[s][m] = 8.0*PI*( 3.0*pressure[s][m]+(energy[s][m]+
                           pressure[s][m])
                           *velocity_sq[s][m]/(1.0-velocity_sq[s][m]))*
                           exp(2.0*alpha[s][m]+0.5*(gama[s][m]-rho[s][m]))
                           *SQ(r_e)*r_e*SQ(s_gp[s]/SQ(1.0-s_gp[s]));

         S_virial2[s][m] = 0.5*( SQ(s1*(d_gama_s+d_rho_s)) + m1*SQ(d_gama_m+
                           d_rho_m)
                           - d_alpha_s*( SQ(s1)*(d_gama_s-d_rho_s) + 2.0*s1)
                           + d_alpha_m*( m1*(-d_gama_m+d_rho_m) + 2.0*mu[m]) 
                           + 2.0*exp(2.0*alpha[s][m]-gama[s][m]+rho[s][m])
                           *(s1*d_alpha_s - mu[m]*d_alpha_m)
                           + 0.5*(1.0-exp(2.0*alpha[s][m]-gama[s][m]+rho[s][m]))
                           *(s1*(d_gama_s-d_rho_s)-mu[m]*(d_gama_m-d_rho_m))
                           -1.5*exp(-2.0*rho[s][m])*SQ(r_e*s_gp[s]/(1.0-s_gp[s]
                           ))*m1*( SQ(s1*d_omega_s)+m1*SQ(d_omega_m)))
                           *exp(0.5*(gama[s][m]-rho[s][m]))*r_e/SQ(1.0-s_gp[s]);
         }
      }

   for(s=1;s<=SDIV;s++) {
      for(m=1;m<=MDIV-2;m+=2) {
         virial1 += (DM/3.0)*(S_virial1[s][m]+4.0*S_virial1[s][m+1]
                                  +S_virial1[s][m+2]);
         virial2 += (DM/3.0)*(S_virial2[s][m]+4.0*S_virial2[s][m+1]
                                  +S_virial2[s][m+2]);
      }     
      D_virial1[s]=virial1;
      D_virial2[s]=virial2;
      virial1=0.0;
      virial2=0.0;
   }

   for(s=1;s<=SDIV-2;s+=2) {
        virial1 += (DS/3.0)*(D_virial1[s]+4.0*D_virial1[s+1]
                                  +D_virial1[s+2]);
        virial2 += (DS/3.0)*(D_virial2[s]+4.0*D_virial2[s+1]
                                  +D_virial2[s+2]);
   }

   grv3=fabs(1.0-virial1/virial2);

/*
   printf("\n");
   printf("virial1=%6.5e, virial2=%6.5e\n",virial1,virial2);  
   printf("grv3=%6.5e\n",grv3);
   printf("\n");
*/
   /* Prepare next guess */

   for(s=1;s<=SDIV;s++)
      for(m=1;m<=MDIV;m++) {
         gama_guess[s][m]=gama[s][m];
         rho_guess[s][m]=rho[s][m];
         alpha_guess[s][m]=alpha[s][m];
         omega_guess[s][m]=omega[s][m];
      }

   r_e_guess=r_e;

}



/***************************************************************************/
void print_header(void)
{ 

 if(print_option==2 || print_option==3) {

 if(strcmp(eos_type,"tab")==0) {

 printf("\n");
 printf("-------------------------------------------------------------------------------\n");   
 printf("%s,  MDIVxSDIV=%dx%d,  accuracy=%1.0e\n",file_name,MDIV,SDIV,accuracy);

 printf("-------------------------------------------------------------------------------\n");   
 printf("   rho_c        M         M_0         R         Omega      Omega_p       T/W           C*J/GM_s^2      I        h_plus       h_minus       Z_p      Z_b        Z_f       omega_c/Omega   r_e     r_ratio     Omega_pa    Omega+    u_phi\n");

 printf("-------------------------------------------------------------------------------\n");   
}
 else {

 printf("\n");
 printf("-------------------------------------------------------------------------------\n");   
 printf(" N=%6.5e,  MDIVxSDIV=%dx%d,  accuracy=%1.0e\n",n_P,MDIV,SDIV,accuracy);

 printf("-------------------------------------------------------------------------------\n");   
 printf("   rho_c        M         M_0         R         Omega      Omega_p       T/W            J          I        h_plus       h_minus       Z_p      Z_b        Z_f       omega_c/Omega   r_e     r_ratio   Omega_pa   Omega+   u_phi \n");

 printf("-------------------------------------------------------------------------------\n");   

 }
}

} 


/***************************************************************************/
void print_table(void)
{ 
 int s,m;

 double om_over_Om,
        si;

  if(r_ratio!=1.0) om_over_Om=omega[1][1]*r_e/Omega_h;
    else
  om_over_Om=0.0;


  switch(print_option) {

  case 1:

  if(strcmp(eos_type,"tab")==0) {
       printf("\n");
       printf("    %s,  MDIVxSDIV=%dx%d,  accuracy=%1.0e\n",file_name,MDIV,SDIV,
              accuracy);
       printf("  ------------------------------------------\n");
       printf("  %6.5e  e_c           (10^15 gr/cm^3)\n",e_center);
       printf("  %6.5e  M             (M_sun)\n",Mass/MSUN);
       printf("  %6.5e  M_0           (M_sun)\n",Mass_0/MSUN);
       printf("  %6.5e  R_e           (km)\n",R_e/1.0e5);
       printf("  %6.5e  Omega         (10^4 s^-1)\n",Omega/1.0e4);
       printf("  %6.5e  Omega_p       (10^4 s^-1)\n",Omega_K/1.0e4);
       printf("  %6.5e  T/W\n",T/W);
       printf("  %6.5e  cJ/GM_sun^2\n",J*C/(G*MSUN*MSUN));
  if(r_ratio!=1.0)
       printf("  %6.5e  I             (10^45 gr cm^2)\n",I/1.0e45);
  else
       printf("      ---      I             (10^45 gr cm^2)\n",I/1.0e45);
  if(h_plus==-99.9) 
       printf("      ---      h+            (km)\n");
  else       
       printf("  %6.5e  h+            (km)\n",h_plus/1.0e5);
  if(h_minus==-99.9) 
       printf("      ---      h-            (km)\n");
  else       
       printf("  %6.5e  h-            (km)\n",h_minus/1.0e5);
  }



  else {

       printf("\n");
       printf("  N=%6.5e,  MDIVxSDIV=%dx%d,  accuracy=%1.0e\n",n_P,MDIV,SDIV,
              accuracy);
       printf("  ------------------------------------------\n");
       printf("  %6.5e  e_c       \n",e_center);
       printf("  %6.5e  M         \n",Mass);
       printf("  %6.5e  M_0       \n",Mass_0);
       printf("  %6.5e  R_e       \n",R_e);
       printf("  %6.5e  Omega     \n",Omega);
       printf("  %6.5e  Omega_p   \n",Omega_K);
       printf("  %6.5e  T/W       \n",T/W);
       printf("  %6.5e  J         \n",J);
  if(r_ratio!=1.0)
       printf("  %6.5e  I         \n",I);
  else
       printf("      ---      I   \n");
  if(h_plus==-99.9) 
       printf("      ---      h+  \n");
  else       
       printf("  %6.5e  h+        \n",h_plus);
  if(h_minus==-99.9) 
       printf("      ---      h-  \n");
  else       
       printf("  %6.5e  h-        \n",h_minus);
  }

       printf("  %6.5e  Z_p       \n",Z_p);
  if(Z_f<=0.0)
       printf(" %6.5e  Z_f        \n",Z_f);
  else
       printf("  %6.5e  Z_f       \n",Z_f); 
       printf("  %6.5e  Z_b       \n",Z_b);
       printf("  %6.5e  omega_c/Omega\n",om_over_Om);
  if(strcmp(eos_type,"tab")==0)
       printf("  %6.5e  r_e           (km)\n",r_e*sqrt(KAPPA)/1.0e5);      
  else
       printf("  %6.5e  r_e        \n",r_e);      
       printf("  %6.5e  r_p/r_e\n",r_ratio);
       printf("  ------------------------------------------\n");
       printf("\n");
       break;

  case 2:

  if(strcmp(eos_type,"tab")==0) {
  
    printf("%6.5e %6.5e %6.5e %6.5e %6.5e %6.5e",e_center*1e15,Mass/MSUN,
           Mass_0/MSUN,R_e/1.0e5,Omega/1.0e4,Omega_K/1e4);

if(r_ratio!=1.0) {
  if(h_plus==-99.9 && h_minus==-99.9) 
    printf(" %6.5e %6.5e %6.5e     ---         ---    ",T/W,C*J/(G*MSUN*MSUN),
                                                         I/1e45);
  else{
    if(h_plus==-99.9) 
     printf(" %6.5e %6.5e %6.5e     ---     %6.5e",T/W,C*J/(G*MSUN*MSUN),
                                                          I/1e45,h_minus/1.0e5);
    else{
         if(h_minus==-99.9)
           printf(" %6.5e %6.5e %6.5e %6.5e    ---    ",T/W,C*J/(G*MSUN*MSUN),
                   I/1e45,h_plus/1.0e5);  
         else
            printf(" %6.5e %6.5e %6.5e %6.5e %6.5e",T/W,C*J/(G*MSUN*MSUN),
                   I/1e45,h_plus/1.0e5,h_minus/1.0e5);  
    }
  }
}else {
  if(h_plus==-99.9 && h_minus==-99.9) 
    printf(" %6.5e %6.5e    ---       ---         ---    ",T/W,C*J/(G*MSUN*MSUN));
  else{
    if(h_plus==-99.9) 
     printf(" %6.5e %6.5e      ---        ---     %6.5e",T/W,C*J/(G*MSUN*MSUN),
                                                                   h_minus/1.0e5);
    else{
         if(h_minus==-99.9)
           printf(" %6.5e %6.5e    ---    %6.5e    ---    ",T/W,C*J/(G*MSUN*MSUN),
                                                                   h_plus/1.0e5);
  
         else
            printf(" %6.5e %6.5e    ---     %6.5e %6.5e",T/W,C*J/(G*MSUN*MSUN),
                                                     h_plus/1.0e5,h_minus/1.0e5); 
    }
  }
}
    printf(" %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e\n"
           ,Z_p,Z_b,Z_f,om_over_Om,r_e,r_ratio,Omega_p/1e4,Omega_plus/1e4,u_phi);

 }
 else {

  
    printf("%6.5e %6.5e %6.5e %6.5e %6.5e %6.5e",e_center,Mass,
           Mass_0,R_e,Omega,Omega_K);
 
if(r_ratio!=1.0) {
  if(h_plus==-99.9 && h_minus==-99.9) 
    printf(" %6.5e %6.5e %6.5e     ---         ---    ",T/W,J,I);
  else{
    if(h_plus==-99.9) 
     printf(" %6.5e %6.5e %6.5e     ---     %6.5e",T/W,J,I,h_minus);
    else{
         if(h_minus==-99.9)
           printf(" %6.5e %6.5e %6.5e %6.5e    ---    ",T/W,J,I,h_plus);  
         else
            printf(" %6.5e %6.5e %6.5e %6.5e %6.5e",T/W,J,I,h_plus,h_minus);  
    }
  }
}else {
  if(h_plus==-99.9 && h_minus==-99.9) 
    printf(" %6.5e %6.5e    ---       ---         ---    ",T/W,J,I);
  else{
    if(h_plus==-99.9) 
     printf(" %6.5e %6.5e      ---        ---     %6.5e",T/W,J,I,h_minus);
     else{
         if(h_minus==-99.9)
           printf(" %6.5e %6.5e    ---    %6.5e    ---    ",T/W,J,I,h_plus);  
         else
            printf(" %6.5e %6.5e    ---     %6.5e %6.5e",T/W,J,I,h_plus,
                                                                h_minus);  
    }
  }
}
    printf(" %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e\n"
           ,Z_p,Z_b,Z_f,om_over_Om,r_e,r_ratio,Omega_p,Omega_plus,u_phi);

 }
    break;
  


  case 3:

  if(strcmp(eos_type,"tab")==0) {

    printf("%6.5e %6.5e %6.5e %6.5e %6.5e %6.5e",e_center*1e15,Mass/MSUN,
           Mass_0/MSUN,R_e/1.0e5,Omega/1.0e4,Omega_K/1e4);

  if(h_plus==-99.9 && h_minus==-99.9) 
    printf(" %6.5e %6.5e %6.5e     ---         ---    ",T/W,C*J/(G*MSUN*MSUN),
          I/1e45);
  else{
    if(h_plus==-99.9) 
     printf(" %6.5e %6.5e %6.5e     ---     %6.5e",T/W,C*J/(G*MSUN*MSUN),
          I/1e45,h_minus/1.0e5);
    else{
         if(h_minus==-99.9)
           printf(" %6.5e %6.5e %6.5e %6.5e    ---    ",T/W,C*J/(G*MSUN*MSUN),
                   I/1e45,h_plus/1.0e5);  
         else
            printf(" %6.5e %6.5e %6.5e %6.5e %6.5e",T/W,C*J/(G*MSUN*MSUN),
                   I/1e45,h_plus/1.0e5,h_minus/1.0e5);  
    }
  }

    printf(" %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e\n"
           ,Z_p,Z_b,Z_f,om_over_Om,r_e,r_ratio);

  printf("-------------------------------------------------------------------------------\n");   
  printf("r/(r+r_e)  cos(theta)     rho        gamma      alpha      omega     pressure\n");
  printf("-------------------------------------------------------------------------------\n");   

  for(s=1;s<=SDIV;s++)
     for(m=1;m<=MDIV;m++)
        printf("%5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e\n",s_gp[s],mu[m],
               rho[s][m],gama[s][m],alpha[s][m],omega[s][m]*C/(sqrt(KAPPA)*1e4),
               pressure[s][m]/KSCALE);
  }
  else {


    printf("%6.5e %6.5e %6.5e %6.5e %6.5e %6.5e",e_center,Mass,
           Mass_0,R_e,Omega,Omega_K);

  if(h_plus==-99.9 && h_minus==-99.9) 
    printf(" %6.5e %6.5e %6.5e     ---         ---    ",T/W,J,I);
  else{
    if(h_plus==-99.9) 
     printf(" %6.5e %6.5e %6.5e     ---     %6.5e",T/W,I,J,h_minus);
    else{
         if(h_minus==-99.9)
           printf(" %6.5e %6.5e %6.5e %6.5e    ---    ",T/W,J,I,h_plus);  
         else
            printf(" %6.5e %6.5e %6.5e %6.5e %6.5e",T/W,J,I,h_plus,h_minus);  
    }
  }

    printf(" %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e\n"
           ,Z_p,Z_b,Z_f,om_over_Om,r_e,r_ratio);

  printf("-------------------------------------------------------------------------------\n");   
  printf("r/(r+r_e)  cos(theta)     rho        gamma      alpha      omega     pressure\n");
  printf("-------------------------------------------------------------------------------\n");   

  for(s=1;s<=SDIV;s++)
     for(m=1;m<=MDIV;m++)
        printf("%5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e\n",s_gp[s],mu[m],
               rho[s][m],gama[s][m],alpha[s][m],omega[s][m],pressure[s][m]);
  }

  break;

 case 4:

 if(r_ratio==1.0) {

  if(strcmp(eos_type,"tab")==0) {
       printf("  %s\n",file_name);
       printf("\n");
       printf("  ---------------------\n");
       printf("  %8.7e  e_c           (10^15 gr/cm^3)\n",e_center);
       printf("  %8.7e  Omega         (10^4 s^-1)\n",0.0);
       printf("  %8.7e  Omega_p       (10^4 s^-1)\n",0.0);
       printf("  %8.7e  M             (M_sun)\n",Mass/MSUN);
       printf("  %8.7e  M_0           (M_sun)\n",Mass_0/MSUN);
       printf("  %8.7e  R_eq          (km)\n",R_e/1.0e5);
       printf("  %8.7e  cJ/GM_sun^2       \n",0.0);
       printf("  ---            I         \n");
       printf("  %8.7e  T/W\n",0.0);
       printf("  %8.7e  Z_p\n",Z_p);
       printf("  %8.7e  Z_f\n",Z_f);
       printf("  %8.7e  Z_b\n",Z_b);
       printf("  %8.7e  r_p/r_eq\n",r_ratio);
       printf("  %8.7e  e\n",0.0);
       printf("  %8.7e  GRV2\n",grv2_new);
       printf("  %8.7e  GRV3\n",grv3);
       printf("  ---------------------\n");
       printf("\n");
       printf("  %8.7e  SMAX\n",SMAX);
       printf("          %dx%d         SDIVxMDIV\n",SDIV,MDIV);
       printf("          %d            LMAX\n",LMAX);
       printf("          %1.0e         accuracy\n",accuracy);
   }
   else{
       printf("  %8.7e       %8.7e    n, gamma\n",n_P,Gamma_P);
       printf("\n");
       printf("  ---------------------\n");
       printf("  %8.7e  e_c           \n",e_center);
       printf("  %8.7e  Omega         \n",0.0);
       printf("  %8.7e  Omega_p       \n",0.0);
       printf("  %8.7e  M             \n",Mass);
       printf("  %8.7e  M_0           \n",Mass_0);
       printf("  %8.7e  R_eq          \n",R_e);
       printf("  %8.7e  J             \n",J);
       printf("  ---            I     \n");
       printf("  %8.7e  T/W\n",0.0);
       printf("  %8.7e  Z_p\n",Z_p);
       printf("  %8.7e  Z_f\n",Z_f);
       printf("  %8.7e  Z_b\n",Z_b);
       printf("  %8.7e  r_p/r_eq\n",r_ratio);
       printf("  %8.7e  e\n",0.0);
       printf("  %8.7e  GRV2\n",grv2_new);
       printf("  %8.7e  GRV3\n",grv3);
       printf("  ---------------------\n");
       printf("\n");
       printf("  %8.7e  SMAX\n",SMAX);
       printf("          %dx%d         SDIVxMDIV\n",SDIV,MDIV);
       printf("          %d            LMAX\n",LMAX);
       printf("          %1.0e         accuracy\n",accuracy);
     }
 }   
 else {

     if(strcmp(eos_type,"tab")==0) {
       printf("  %s\n",file_name);
       printf("\n");
       printf("  ---------------------\n");
       printf("  %8.7e  e_c           (10^15 gr/cm^3)\n",e_center);
       printf("  %8.7e  Omega         (10^4 s^-1)\n",Omega/1.0e4);
       printf("  %8.7e  Omega_p       (10^4 s^-1)\n",Omega_K/1.0e4);
       printf("  %8.7e  M             (M_sun)\n",Mass/MSUN);
       printf("  %8.7e  M_0           (M_sun)\n",Mass_0/MSUN);
       printf("  %8.7e  R_eq          (km)\n",R_e/1.0e5);
       printf("  %8.7e  cJ/GM_sun^2   \n",J*C/(G*MSUN*MSUN));
       printf("  %8.7e  I             (10^45 gr cm^2)\n",I/1.0e45);
       printf("  %8.7e  T/W\n",T/W);
       printf("  %8.7e  Z_p\n",Z_p);
       printf("  %8.7e  Z_f\n",Z_f);
       printf("  %8.7e  Z_b\n",Z_b);
       printf("  %8.7e  r_p/r_eq\n",r_ratio);
       printf("  %8.7e  e\n",eccentricity);
       printf("  %8.7e  GRV2\n",grv2_new);
       printf("  %8.7e  GRV3\n",grv3);
       printf("  ---------------------\n");
       printf("\n");
       printf("  %8.7e  SMAX\n",SMAX);
       printf("          %dx%d        SDIVxMDIV\n",SDIV,MDIV);
       printf("          %d             LMAX\n",LMAX);
       printf("          %1.0e          accuracy\n",accuracy);
     }
     else{
       printf("  %8.7e       %8.7e    n, gamma\n",n_P,Gamma_P);
       printf("\n");
       printf("  ---------------------\n");
       printf("  %8.7e  e_c               \n",e_center);
       printf("  %8.7e  Omega             \n",Omega);
       printf("  %8.7e  Omega_p           \n",Omega_K);
       printf("  %8.7e  M                 \n",Mass);
       printf("  %8.7e  M_0               \n",Mass_0);
       printf("  %8.7e  R_eq              \n",R_e);
       printf("  %8.7e  J                 \n",J);
       printf("  %8.7e  I                 \n",I);
       printf("  %8.7e  T/W               \n",T/W);
       printf("  %8.7e  Z_p               \n",Z_p);
       printf("  %8.7e  Z_f               \n",Z_f);
       printf("  %8.7e  Z_b               \n",Z_b);
       printf("  %8.7e  r_p/r_eq          \n",r_ratio);
       printf("  %8.7e  e                 \n",eccentricity);
       printf("  %8.7e  GRV2              \n",grv2_new);
       printf("  %8.7e  GRV3              \n",grv3);
       printf("  ---------------------\n");
       printf("\n");
       printf("  %8.7e  SMAX\n",SMAX);
       printf("          %dx%d        SDIVxMDIV\n",SDIV,MDIV);
       printf("          %d             LMAX\n",LMAX);
       printf("          %1.0e          accuracy\n",accuracy);
     }
   }

 break;

 case 6:
        for(s=1;s<=SDIV/2;s++)
           for(m=1;m<=MDIV;m++) {
              si=s_gp[s]/(1.0-s_gp[s]);
              printf("%6.5e %6.5e %6.5e %6.5e\n",si,mu[m],pressure[s][m],
                                                 energy[s][m]);
	    }
 break;

 }
}   

/**************************************************************************/
double dm_dr_is(double r_is, double r, double m, double p)
{
 double dmdr,
        e_d;

 if(p<p_surface) 
    e_d=0.0;
 else 
    e_d=e_at_p(p);
 
 if(r_is<RMIN) 
    dmdr=4.0*PI*e_center*r*r*(1.0+4.0*PI*e_center*r*r/3.0);
 else
    dmdr=4.0*PI*e_d*r*r*r*sqrt(1.0-2.0*m/r)/r_is;
 
return dmdr;
}
 
/**************************************************************************/
double dp_dr_is(double r_is, double r, double m, double p)

{ double dpdr,
         e_d; 

  if(p<p_surface) e_d=0.0;
  else        
   e_d=e_at_p(p);
  
  if(r_is<RMIN) dpdr = -4.0*PI*(e_center+p)*(e_center+3.0*p)*r*(1.0
                     +4.0*e_center*r*r/3.0)/3.0;

  else 
   dpdr = -(e_d+p)*(m+4.0*PI*r*r*r*p)/(r*r_is*sqrt(1.0-2.0*m/r));

 return dpdr;
}

/**************************************************************************/
double dr_dr_is(double r_is, double r, double m)
{
 double drdris;

 if(r_is<RMIN) drdris=1.0;
  else
   drdris=(r/r_is)*sqrt(1.0-2.0*m/r);

 return drdris;
}

/************************************************************************/
void integrate(int i_check)
{
  int i=2;

  double r,                           /* radius */
         r_is,                        /* isotropic radial coordinate */
         m,                           /* mass   */
         p,                           /* pressure */
         e_d,                         /* density */
         r_is_est,                    /* estimate on final isotr. radius */ 
         dr_is_save,                  /* r_is saving interval */  
         r_is_check,                  /*                      */    
         nu_s,
         hh,
         h,                           /* stepsize during integration */
         a1,a2,a3,a4,b1,b2,b3,b4,     /* coeff. in Runge-Kutta equations */
         c1,c2,c3,c4,
         rho0;   

    if(i_check==1) {
      if(strcmp(eos_type,"tab")==0)
        r_is_est=1.5e6/sqrt(KAPPA);
      else
        r_is_est=2.0*sqrt(Gamma_P/(4.0*PI*(Gamma_P-1.0)))*
                            pow(e_center,(Gamma_P-2.0)/2.0);

      h=r_is_est/100;     
    }
    else {
          r_is_est=r_is_final;
          h=r_is_est/10000;     
	 }

    dr_is_save=r_is_final/RDIV;
    r_is_check=dr_is_save;
    r_is=0.0;                            /* initial isotropic radius */
    r=0.0;                               /* initial radius */
    m=0.0;                               /* initial mass */ 
    p=p_center;                          /* initial pressure */ 

    r_is_gp[1]=0.0;
    r_gp[1]=0.0;
    m_gp[1]=0.0;
    lambda_gp[1]=0.0;
    e_d_gp[1]=e_center; 

    while (p>=p_surface) { 

      e_d=e_at_p(p);

     if((i_check==3) && (r_is>r_is_check) && (i<=RDIV)) {
      r_is_gp[i]=r_is;
      r_gp[i]=r;
      m_gp[i]=m;
      e_d_gp[i]=e_d; 
      i++;   
      r_is_check += dr_is_save;
     }    
       
     r_is_final=r_is;
     r_final=r;
     m_final=m;
 
     a1=dr_dr_is(r_is,r,m);
     b1=dm_dr_is(r_is,r,m,p);
     c1=dp_dr_is(r_is,r,m,p);

 
     a2=dr_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0);
     b2=dm_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0, p+h*c1/2.0);
     c2=dp_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0, p+h*c1/2.0);

     a3=dr_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0);
     b3=dm_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0, p+h*c2/2.0);
     c3=dp_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0, p+h*c2/2.0);

     a4=dr_dr_is(r_is+h, r+h*a3, m+h*b3);
     b4=dm_dr_is(r_is+h, r+h*a3, m+h*b3, p+h*c3);
     c4=dp_dr_is(r_is+h, r+h*a3, m+h*b3, p+h*c3);

     
     r += (h/6.0)*(a1+2*a2+2*a3+a4);
     m += (h/6.0)*(b1+2*b2+2*b3+b4);
     p += (h/6.0)*(c1+2*c2+2*c3+c4);

     r_is += h;
    }

    r_is_gp[RDIV]=r_is_final;
    r_gp[RDIV]=r_final;
    m_gp[RDIV]=m_final;

/* Rescale r_is and compute lambda */

   if(i_check==3) {
      k_rescale=0.5*(r_final/r_is_final)*(1.0-m_final/r_final+
                sqrt(1.0-2.0*m_final/r_final));
 
      r_is_final *= k_rescale;
 
      nu_s = log((1.0-m_final/(2.0*r_is_final))/(1.0+m_final/
               (2.0*r_is_final)));

      for(i=1;i<=RDIV;i++) {
         r_is_gp[i] *= k_rescale;
 
         if(i==1) lambda_gp[1]= log(1.0/k_rescale);
           else
               lambda_gp[i]=log(r_gp[i]/r_is_gp[i]); 

         if(e_d_gp[i]<e_surface) 
 
           hh=0.0;
 
         else { if(strcmp(eos_type,"tab")==0) {
                  p=p_at_e(e_d_gp[i]);
                  hh=h_at_p(p);
                }
                else { 
                      rho0=rtsec(e_of_rho0,0.0,e_d_gp[i],1e-12,e_d_gp[i]);
                      p=pow(rho0,Gamma_P);
                      hh=log((e_d_gp[i]+p)/rho0);
                }
              }
 
         nu_gp[i]=nu_s-hh;
      }
      nu_gp[RDIV]=nu_s;

   }
}


/*************************************************************************/
void guess(void)
{
 int i,
     s,
     m,
     s_temp;

 double r_is_s,
        lambda_s,
        nu_s,
        gama_eq,
        rho_eq;

 n_nearest=n_tab/2;

 integrate(1);
 integrate(2);
 integrate(3);


 if(SMAX==1.0) s_temp=SDIV-1;
 else 
     s_temp=SDIV;


 n_nearest=RDIV/2;

 for(s=1;s<=s_temp;s++) {
    r_is_s=r_is_final*(s_gp[s]/(1.0-s_gp[s]));

    if(r_is_s<r_is_final) {
      lambda_s=interp(r_is_gp,lambda_gp,RDIV,r_is_s);
      nu_s=interp(r_is_gp,nu_gp,RDIV,r_is_s);
    }
    else {
      lambda_s=2.0*log(1.0+m_final/(2.0*r_is_s));
      nu_s=log((1.0-m_final/(2.0*r_is_s))/(1.0+m_final/(2*r_is_s)));
    }

    gama[s][1]=nu_s+lambda_s;
    rho[s][1]=nu_s-lambda_s;

    for(m=1;m<=MDIV;m++) {
        gama_guess[s][m]=gama[s][1],        
        rho_guess[s][m]=rho[s][1],
        alpha_guess[s][m]=(gama[s][1]-rho[s][1])/2.0;
        omega_guess[s][m]=0.0; 
    }
 
    gama_mu_0[s]=gama[s][1];                   /* gama at \mu=0 */
    rho_mu_0[s]=rho[s][1];                     /* rho at \mu=0 */
 }

 if(SMAX==1.0) {
    for(m=1;m<=MDIV;m++) {
        gama_guess[SDIV][m]=0.0;
        rho_guess[SDIV][m]=0.0;
        alpha_guess[SDIV][m]=0.0;
        omega_guess[SDIV][m]=0.0; 
    }
 
    gama_mu_0[SDIV]=0.0;
    rho_mu_0[SDIV]=0.0;
   
  }
   

   n_nearest=SDIV/2;
   gama_eq=interp(s_gp,gama_mu_0,SDIV,s_e);      /* gama at equator */
   rho_eq=interp(s_gp,rho_mu_0,SDIV,s_e);        /* rho at equator */
 
   r_e_guess= r_final*exp(0.5*(rho_eq-gama_eq)); 

   R_e=r_final*sqrt(KAPPA);
   Mass=m_final*sqrt(KAPPA)*(C*C/G);
   Z_p=exp(-0.5*(gama_eq+rho_eq))-1.0;
}


/*************************************************************************/
/* Main iteration cycle.                                                 */
/*************************************************************************/
void iterate(double r_ratio)
{
 int m,                      /* counter */
     s,                      /* counter */
     n,                      /* counter */
     k,                      /* counter */
     n_of_it=0,
     s_temp;
 
double   s_term=0.0,          /* term in sum */       
	 sum_rho=0.0,         /* intermediate sum in eqn for rho */
	 sum_gama=0.0,        /* intermediate sum in eqn for gama */
	 sum_omega=0.0,       /* intermediate sum in eqn for omega */
         r_e_old,             /* equatorial radius in previus cycle */
   	 dif=1.0,             /* difference | r_e_old -r_e | */
         d_gama_s,            /* derivative of gama w.r.t. s */
         d_gama_m,            /* derivative of gama w.r.t. m */
         d_rho_s,             /* derivative of rho w.r.t. s */
         d_rho_m,             /* derivative of rho w.r.t. m */
         d_omega_s,           /* derivative of omega w.r.t. s */
         d_omega_m,           /* derivative of omega w.r.t. m */
         d_gama_ss,           /* 2nd derivative of gama w.r.t. s */
         d_gama_mm,           /* 2nd derivative of gama w.r.t. m */
         d_gama_sm,           /* derivative of gama w.r.t. m and s */
         temp1,                /* temporary term in da_dm */ 
         temp2, 
         temp3,
         temp4,
         temp5,
         temp6,
         temp7,
         temp8,
         m1,                  
         s1,
         s2,
         ea,
         R_emb_eq,
         R_emb_pole,
         rsm,
         gsm,
         esm,
         psm,
         v2sm,
         mum,
         omsm,
         sgp,
         s_1,
         e_gsm,
         e_rsm, 
         rho0sm;

 
 if(SMAX==1.0) s_temp=SDIV-1;
 else
     s_temp=SDIV;


 for(s=1;s<=SDIV;s++)
    for(m=1;m<=MDIV;m++) {
       gama[s][m]=gama_guess[s][m];
       rho[s][m]=rho_guess[s][m];
       alpha[s][m]=alpha_guess[s][m];
       omega[s][m]=omega_guess[s][m];
    }

 r_e=r_e_guess;

 gama_center_old=0.0;
 rho_center_old=0.0;


 while(dif> accuracy || n_of_it<2) { 

 if(print_dif!=0)
      printf("%4.3e\n",dif);

 
      /* Rescale potentials and construct arrays with the potentials along
       | the equatorial and polar directions.
      */        

      for(s=1;s<=SDIV;s++) {
         for(m=1;m<=MDIV;m++) {
            rho[s][m] /= SQ(r_e);
            gama[s][m] /= SQ(r_e); 
            alpha[s][m] /= SQ(r_e);
            omega[s][m] *= r_e;
         }
         rho_mu_0[s]=rho[s][1];     
         gama_mu_0[s]=gama[s][1];   
         omega_mu_0[s]=omega[s][1]; 
         rho_mu_1[s]=rho[s][MDIV];  
         gama_mu_1[s]=gama[s][MDIV];
      }
 
      /* Compute new r_e. */ 

      r_e_old=r_e;
      r_p=r_ratio*r_e;                          
      s_p=r_p/(r_p+r_e);                        
  
      n_nearest= SDIV/2;
      gama_pole_h=interp(s_gp,gama_mu_1,SDIV,s_p); 
      gama_equator_h=interp(s_gp,gama_mu_0,SDIV,s_e);
      gama_center_h=gama[1][1];                    
  
      rho_pole_h=interp(s_gp,rho_mu_1,SDIV,s_p);   
      rho_equator_h=interp(s_gp,rho_mu_0,SDIV,s_e);
      rho_center_h=rho[1][1];                      
 
      r_e=sqrt(2*h_center/(gama_pole_h+rho_pole_h-gama_center_h-rho_center_h));

      /* Compute angular velocity Omega. */
 
      if(r_ratio==1.0) {
        Omega_h=0.0;
        omega_equator_h=0.0;
      } 
      else {
            omega_equator_h=interp(s_gp,omega_mu_0,SDIV,s_e);
 
            Omega_h = omega_equator_h + exp(SQ(r_e)*rho_equator_h)*sqrt(1.0
                      -exp(SQ(r_e)*(gama_pole_h+rho_pole_h-gama_equator_h
                      -rho_equator_h)));
      }
 
      /* Compute velocity, energy density and pressure. */
 
      n_nearest=n_tab/2; 

      for(s=1;s<=SDIV;s++) {
         sgp=s_gp[s];

         for(m=1;m<=MDIV;m++) {
            rsm=rho[s][m];
            
            if(r_ratio==1.0 || s > (SDIV/2+2) ) 
                velocity_sq[s][m]=0.0;
            else 
                velocity_sq[s][m]=SQ((Omega_h-omega[s][m])*(sgp/(1.0-sgp))
                                  *sin_theta[m]*exp(-rsm*SQ(r_e)));

            if(velocity_sq[s][m]>=1.0) 
              velocity_sq[s][m]=0.0;

            enthalpy[s][m]=enthalpy_min + 0.5*(SQ(r_e)*(gama_pole_h+rho_pole_h
                           -gama[s][m]-rsm)-log(1.0-velocity_sq[s][m]));
  
            if((enthalpy[s][m]<=enthalpy_min) || (sgp>s_e)) {
                  pressure[s][m]=0.0;
                  energy[s][m]=0.0; 
	    }
            else { if(strcmp(eos_type,"tab")==0) {
                     pressure[s][m]=p_at_h(enthalpy[s][m]);
                     energy[s][m]=e_at_p(pressure[s][m]);
	           }
                   else {
                         rho0sm=pow(((Gamma_P-1.0)/Gamma_P)
                                *(exp(enthalpy[s][m])-1.0),1.0/(Gamma_P-1.0));
 
                         pressure[s][m]=pow(rho0sm,Gamma_P);

                         energy[s][m]=pressure[s][m]/(Gamma_P-1.0)+rho0sm;
                   }
            }  

            /* Rescale back metric potentials (except omega) */

            rho[s][m] *= SQ(r_e);
            gama[s][m] *= SQ(r_e);
            alpha[s][m] *= SQ(r_e);
	 }
      }

      /* EVALUATION OF SOURCE TERMS */

      if(SMAX==1.0) {
        s=SDIV;
        for(m=1;m<=MDIV;m++) {
           S_rho[s][m] = 0.0;
           S_gama[s][m]=0.0;
           S_omega[s][m]=0.0;
        }
      }

      for(s=1;s<=s_temp;s++)
         for(m=1;m<=MDIV;m++) {
            rsm=rho[s][m];
            gsm=gama[s][m];
            omsm=omega[s][m];
            esm=energy[s][m];
            psm=pressure[s][m];
            e_gsm=exp(0.5*gsm);
            e_rsm=exp(-rsm);
            v2sm=velocity_sq[s][m];
            mum=mu[m];            
            m1=1.0-SQ(mum);
            sgp=s_gp[s];
            s_1=1.0-sgp;
            s1=sgp*s_1;
            s2=SQ(sgp/s_1);  

            ea=16.0*PI*exp(2.0*alpha[s][m])*SQ(r_e);
 
            if(s==1) {
              d_gama_s=0.0;
              d_gama_m=0.0;
              d_rho_s=0.0;
              d_rho_m=0.0;
              d_omega_s=0.0;
              d_omega_m=0.0;
            }else{
                 d_gama_s=deriv_s(gama,s,m);
                 d_gama_m=deriv_m(gama,s,m);
                 d_rho_s=deriv_s(rho,s,m);
                 d_rho_m=deriv_m(rho,s,m);
                 d_omega_s=deriv_s(omega,s,m);
                 d_omega_m=deriv_m(omega,s,m);
	     }

            S_rho[s][m] = e_gsm*(0.5*ea*(esm + psm)*s2*(1.0+v2sm)/(1.0-v2sm)
  
                          + s2*m1*SQ(e_rsm)*(SQ(s1*d_omega_s) + m1*SQ(d_omega_m))
                         
                          + s1*d_gama_s - mum*d_gama_m + 0.5*rsm*(ea*psm*s2  
 
                          - s1*d_gama_s*(0.5*s1*d_gama_s+1.0) 
 
                          - d_gama_m*(0.5*m1*d_gama_m-mum)));

            S_gama[s][m] = e_gsm*(ea*psm*s2 + 0.5*gsm*(ea*psm*s2 - 0.5*SQ(s1

                           *d_gama_s) - 0.5*m1*SQ(d_gama_m)));

            S_omega[s][m]=e_gsm*e_rsm*( -ea*(Omega_h-omsm)*(esm+psm)

                          *s2/(1.0-v2sm) + omsm*( -0.5*ea*(((1.0+v2sm)*esm + 2.0 
      
                          *v2sm*psm)/(1.0-v2sm))*s2 - s1*(2*d_rho_s+0.5*d_gama_s)

                          + mum*(2*d_rho_m+0.5*d_gama_m) + 0.25*SQ(s1)*(4

                          *SQ(d_rho_s)-SQ(d_gama_s)) + 0.25*m1*(4*SQ(d_rho_m)

                          - SQ(d_gama_m)) - m1*SQ(e_rsm)*(SQ(SQ(sgp)*d_omega_s) 

                          + s2*m1*SQ(d_omega_m))));
	 }

      
      /* ANGULAR INTEGRATION */
  
      n=0;
      for(k=1;k<=SDIV;k++) {      

         for(m=1;m<=MDIV-2;m+=2) {
               sum_rho += (DM/3.0)*(P_2n[m][n]*S_rho[k][m]
                          + 4.0*P_2n[m+1][n]*S_rho[k][m+1] 
                          + P_2n[m+2][n]*S_rho[k][m+2]);
	 }

         D1_rho[n][k]=sum_rho;
         D1_gama[n][k]=0.0;
         D1_omega[n][k]=0.0;
         sum_rho=0.0;

      }

      for(n=1;n<=LMAX;n++)
         for(k=1;k<=SDIV;k++) {      
            for(m=1;m<=MDIV-2;m+=2) {

               sum_rho += (DM/3.0)*(P_2n[m][n]*S_rho[k][m]
                          + 4.0*P_2n[m+1][n]*S_rho[k][m+1] 
                          + P_2n[m+2][n]*S_rho[k][m+2]);
                       
               sum_gama += (DM/3.0)*(sin((2.0*n-1.0)*theta[m])*S_gama[k][m]
                           +4.0*sin((2.0*n-1.0)*theta[m+1])*S_gama[k][m+1]
                           +sin((2.0*n-1.0)*theta[m+2])*S_gama[k][m+2]);
  
               sum_omega += (DM/3.0)*(sin_theta[m]*P1_2n_1[m][n]*S_omega[k][m]
                            +4.0*sin_theta[m+1]*P1_2n_1[m+1][n]*S_omega[k][m+1]
                            +sin_theta[m+2]*P1_2n_1[m+2][n]*S_omega[k][m+2]);
	    }
            D1_rho[n][k]=sum_rho;
            D1_gama[n][k]=sum_gama;
            D1_omega[n][k]=sum_omega;
            sum_rho=0.0;
            sum_gama=0.0;
            sum_omega=0.0;
	}


      /* RADIAL INTEGRATION */

      n=0;
      for(s=1;s<=SDIV;s++) {
            for(k=1;k<=SDIV-2;k+=2) { 
               sum_rho += (DS/3.0)*( f_rho[s][n][k]*D1_rho[n][k] 
                          + 4.0*f_rho[s][n][k+1]*D1_rho[n][k+1]
                          + f_rho[s][n][k+2]*D1_rho[n][k+2]);
 	    }
	    D2_rho[s][n]=sum_rho;
	    D2_gama[s][n]=0.0;
	    D2_omega[s][n]=0.0;
            sum_rho=0.0;
	 }


      for(s=1;s<=SDIV;s++)
         for(n=1;n<=LMAX;n++) {
            for(k=1;k<=SDIV-2;k+=2) { 
               sum_rho += (DS/3.0)*( f_rho[s][n][k]*D1_rho[n][k] 
                          + 4.0*f_rho[s][n][k+1]*D1_rho[n][k+1]
                          + f_rho[s][n][k+2]*D1_rho[n][k+2]);
 
               sum_gama += (DS/3.0)*( f_gama[s][n][k]*D1_gama[n][k] 
                           + 4.0*f_gama[s][n][k+1]*D1_gama[n][k+1]
                           + f_gama[s][n][k+2]*D1_gama[n][k+2]);
 
               sum_omega += (DS/3.0)*( f_omega[s][n][k]*D1_omega[n][k] 
                            + 4.0*f_omega[s][n][k+1]*D1_omega[n][k+1]
                            + f_omega[s][n][k+2]*D1_omega[n][k+2]);
	    }
	    D2_rho[s][n]=sum_rho;
	    D2_gama[s][n]=sum_gama;
	    D2_omega[s][n]=sum_omega;
            sum_rho=0.0;
            sum_gama=0.0;
            sum_omega=0.0;
	 }


      /* SUMMATION OF COEFFICIENTS */

      for(s=1;s<=SDIV;s++) 
         for(m=1;m<=MDIV;m++) {

            gsm=gama[s][m];
            rsm=rho[s][m];
            omsm=omega[s][m];             
            e_gsm=exp(-0.5*gsm);
            e_rsm=exp(rsm);
            temp1=sin_theta[m];

            sum_rho += -e_gsm*P_2n[m][0]*D2_rho[s][0]; 

            for(n=1;n<=LMAX;n++) {

               sum_rho += -e_gsm*P_2n[m][n]*D2_rho[s][n]; 

               if(m==MDIV) {             
                 sum_omega += 0.5*e_rsm*e_gsm*D2_omega[s][n]; 
                 sum_gama += -(2.0/PI)*e_gsm*D2_gama[s][n];   
	       }
               else { 
                     sum_omega += -e_rsm*e_gsm*(P1_2n_1[m][n]/(2.0*n*(2.0*n-1.0)
                                  *temp1))*D2_omega[s][n];  
                     sum_gama += -(2.0/PI)*e_gsm*(sin((2.0*n-1.0)*theta[m])
                                 /((2.0*n-1.0)*temp1))*D2_gama[s][n];   
	       }
	    }
	   
            rho[s][m]=rsm + cf*(sum_rho-rsm);
            gama[s][m]=gsm + cf*(sum_gama-gsm);
            omega[s][m]=omsm + cf*(sum_omega-omsm);

            sum_omega=0.0;
            sum_rho=0.0;
            sum_gama=0.0; 
	  }


      /* CHECK FOR DIVERGENCE */

      if(fabs(omega[2][1])>100.0 || fabs(rho[2][1])>100.0 
         || fabs(gama[2][1])>300.0) {
         a_check=200; 
         break;
      }


      /* TREAT SPHERICAL CASE */

      if(r_ratio==1.0) {
        for(s=1;s<=SDIV;s++)
           for(m=1;m<=MDIV;m++) {
              rho[s][m]=rho[s][1];
              gama[s][m]=gama[s][1];
              omega[s][m]=0.0;          
	   }
      }


      /* TREAT INFINITY WHEN SMAX=1.0 */

      if(SMAX==1.0) {
         for(m=1;m<=MDIV;m++) {
            rho[SDIV][m]=0.0;
            gama[SDIV][m]=0.0;
            omega[SDIV][m]=0.0;
	 }
      } 

      
      /* COMPUTE FIRST ORDER DERIVATIVES OF GAMA */ 
 
      for(s=1;s<=SDIV;s++)
         for(m=1;m<=MDIV;m++) {
            dgds[s][m]=deriv_s(gama,s,m);
            dgdm[s][m]=deriv_m(gama,s,m);
	 }



      /* ALPHA */
 
      if(r_ratio==1.0) {
        for(s=1;s<=SDIV;s++)
           for(m=1;m<=MDIV;m++)
              da_dm[s][m]=0.0;
      } 
      else {
            for(s=2;s<=s_temp;s++)
               for(m=1;m<=MDIV;m++) {

                  da_dm[1][m]=0.0; 
       
                  sgp=s_gp[s];
                  s1=sgp*(1.0-sgp);
                  mum=mu[m]; 
                  m1=1.0-SQ(mum);
          
                  d_gama_s=dgds[s][m];
                  d_gama_m=dgdm[s][m];
                  d_rho_s=deriv_s(rho,s,m);
                  d_rho_m=deriv_m(rho,s,m);
                  d_omega_s=deriv_s(omega,s,m);
                  d_omega_m=deriv_m(omega,s,m);
                  d_gama_ss=s1*deriv_s(dgds,s,m)+(1.0-2.0*sgp)*d_gama_s;
                  d_gama_mm=m1*deriv_m(dgdm,s,m)-2.0*mum*d_gama_m;  
                  d_gama_sm=deriv_sm(gama,s,m);

           temp1=2.0*SQ(sgp)*(sgp/(1.0-sgp))*m1*d_omega_s*d_omega_m

                *(1.0+s1*d_gama_s) - (SQ(SQ(sgp)*d_omega_s) - 
 
                SQ(sgp*d_omega_m/(1.0-sgp))*m1)*(-mum+m1*d_gama_m); 
  
           temp2=1.0/(m1 *SQ(1.0+s1*d_gama_s) + SQ(-mum+m1*d_gama_m));

           temp3=s1*d_gama_ss + SQ(s1*d_gama_s);
  
           temp4=d_gama_m*(-mum+m1*d_gama_m);
   
           temp5=(SQ(s1*(d_rho_s+d_gama_s)) - m1*SQ(d_rho_m+d_gama_m))

                 *(-mum+m1*d_gama_m);

           temp6=s1*m1*(0.5*(d_rho_s+d_gama_s)* (d_rho_m+d_gama_m) 
  
                + d_gama_sm + d_gama_s*d_gama_m)*(1.0 + s1*d_gama_s); 

           temp7=s1*mum*d_gama_s*(1.0+s1*d_gama_s);

           temp8=m1*exp(-2*rho[s][m]);
 
          da_dm[s][m] = -0.5*(d_rho_m+d_gama_m) - temp2*(0.5*(temp3 - 

            d_gama_mm - temp4)*(-mum+m1*d_gama_m) + 0.25*temp5 

            - temp6 +temp7 + 0.25*temp8*temp1);	 
       }
   }

      for(s=1;s<=s_temp;s++) {
         alpha[s][1]=0.0;
         for(m=1;m<=MDIV-1;m++) 
            alpha[s][m+1]=alpha[s][m]+0.5*DM*(da_dm[s][m+1]+
                          da_dm[s][m]);
      } 
 
      for(s=1;s<=s_temp;s++)          
         for(m=1;m<=MDIV;m++) {     
            alpha[s][m] += -alpha[s][MDIV]+0.5*(gama[s][MDIV]-rho[s][MDIV]);

            if(alpha[s][m]>=300.0) {
              a_check=200; 
              break;
            }

            omega[s][m] /= r_e;
         } 

      if(SMAX==1.0) {
         for(m=1;m<=MDIV;m++)      
            alpha[SDIV][m] = 0.0;
      }

      if(a_check==200)
        break;

      dif=fabs(r_e_old-r_e)/r_e;
      n_of_it++;


 }   /* end while */


 for(s=1;s<=SDIV;s++)
    for(m=1;m<=MDIV;m++) {
       gama_guess[s][m]=gama[s][m];
       rho_guess[s][m]=rho[s][m];
       alpha_guess[s][m]=alpha[s][m];
       omega_guess[s][m]=omega[s][m];
    }
 
 r_e_guess=r_e;

} 

/*************************************************************************/
/* Compute m/s model for current e_center.                               */
/*************************************************************************/
void ms_model(void)
{ 
    /* First model is guess */
  
    make_center(e_center);
    guess();             
    d_Omega=1.0;           
    sign=1.0;
    dr=0.1;            /* 0.2 */
    r_ratio = 1.0;

    diff_omega=1.0;
 
    /* Decrease r_p. Whenever Omega_K-Omega changes sign, or iteration does
     | not converge (a_check=200), cut the step dr in half and reverse 
     | its direction.
    */

    while((diff_omega>omega_error || d_Omega<0.0) && r_ratio<=1.0) {          
         if(d_Omega*sign<0.0) { 
           sign=d_Omega;        
           dr /= (-2.0);              
         } 
         r_ratio -= dr; 
         a_check=0.0;
         iterate(r_ratio);
         if(a_check==200) {
                d_Omega=-1.0;

	 }
         else {
                comp_omega();
                d_Omega=Omega_K-Omega;
                diff_omega=fabs(Omega-Omega_K)/Omega_K;
         }

         if(r_ratio>=1.0)
           printf("r_ratio>=1.0 !\n");
       }
    comp();
    print_table();
}


/*************************************************************************/
/* Compute h+ = 0 model for current e_center.                            */
/*************************************************************************/
void h_model(void)
{ 
 double d_h,
        diff_h;

    /* First model is guess */
  
    make_center(e_center);
    guess();             
    d_h=1.0;           
    sign=1.0;
    dr=0.1;            /* 0.2 */
    r_ratio = 1.0;

    diff_h=1.0;
 
    /* Decrease r_p. Whenever Omega_K-Omega changes sign, or iteration does
     | not converge (a_check=200), cut the step dr in half and reverse 
     | its direction.
    */

    while((diff_h>h_error || d_h<0.0) && r_ratio<=1.0) {          
         if(d_h*sign<0.0) { 
           sign=d_h;        
           dr /= (-2.0);              
         } 
         r_ratio -= dr;

         if(r_ratio >1.0) { 
           sign= -d_h;        
           dr /= (-2.0);              
           r_ratio -= 3.0*dr;
         }

         a_check=0.0;
         iterate(r_ratio);

         if(a_check==200) {
                d_h=-1.0;

	 }
         else {
                comp();
                if(h_plus>0.0) {
                     d_h=h_plus/1e5;
                     diff_h = h_plus/1e5;
                } 
                else {
                      d_h= -1.0;
                      diff_h = -1.0;
                }  
         }
/*
         if(r_ratio>1.0)           
           printf("r_ratio>1.0 !\n");
*/
       }
    comp();
    print_table();
}

/*************************************************************************/
/* Compute intermediate model for given Mass_0 and e_center.             */
/*************************************************************************/
void m0_model(double M_0)
{
 double dr,
        diff_M_0,
        d_ratio_M_0;

 dr=0.1;
 r_ratio=1.0-dr;

 /* Compute first rotating model */

 make_center(e_center);
 guess();
 a_check=0.0;
 iterate(r_ratio);
 if(a_check==200) {
   diff_M_0=-1.0;
   sign=-1.0;
} 
 else {   
       comp_M_J();
       diff_M_0=M_0-Mass_0;
       sign=diff_M_0;
       d_ratio_M_0=fabs(diff_M_0)/M_0;
 }

 /* If rest mass is greater than desired, reverse direction and cut stepsize
  | in half.
 */
 
 if(diff_M_0<0.0) {
   dr /= (-2.0);
 }

 while(d_ratio_M_0>M_0_error && r_ratio<=1.0) {
      if(diff_M_0*sign<0.0) {
        sign=diff_M_0;
        dr /= (-2.0);
      }
      r_ratio -= dr;
      a_check=0.0;
      iterate(r_ratio);
      if(a_check==200) {
            diff_M_0 = -1.0;
      } 
      else { 
            comp_M_J();      
            if(Omega>Omega_K)
                  diff_M_0=-1.0;
            else {
                  diff_M_0=M_0-Mass_0;
                  d_ratio_M_0=fabs(diff_M_0)/M_0;
            }     
      }
 }
 comp();
 print_table();
}
 

/*************************************************************************/
/* Compute intermediate model for given Mass and e_center.               */
/*************************************************************************/
void m_model(double M_fix)
{
 double dr,
        diff_M,
        d_ratio_M;

 dr=0.1;
 r_ratio=1.0-dr;

 /* Compute first rotating model */

 make_center(e_center);
 guess();
 a_check=0.0;
 iterate(r_ratio);
 if(a_check==200) {
   diff_M=-1.0;
   sign=-1.0;
} 
 else { 
       comp_M_J();
       diff_M=M_fix-Mass;
       sign=diff_M;
       d_ratio_M=fabs(diff_M)/M_fix;
 }

 /* If mass is greater than desired, reverse direction and cut stepsize
  | in half.
 */
 
 if(diff_M<0.0) {
   dr /= (-2.0);
 }

 while(d_ratio_M>M_error && r_ratio<=1.0) {
      if(diff_M*sign<0.0) {
        sign=diff_M;
        dr /= (-2.0);
      }
      r_ratio -= dr;
      a_check=0.0;
      iterate(r_ratio);
      if(a_check==200) {
            diff_M = -1.0;
      }
      else { 
            comp_M_J();      
            if(Omega>Omega_K) {
                  diff_M=-1.0;
            }else {
                  diff_M=M_fix-Mass;
                  d_ratio_M=fabs(diff_M)/M_fix;
            }     
      }
 } 
 comp();
 print_table();
}


/*************************************************************************/
/* Compute intermediate model for given Omega and e_center.              */
/*************************************************************************/
void omega_model(double Omega_const)
{
 double dr,
        diff_Omega,
        d_ratio_Omega;

 dr=0.1;
 r_ratio=1.0-dr;

 /* Compute first rotating model */

 make_center(e_center);
 guess();
 a_check=0.0;
 iterate(r_ratio);
 if(a_check==200) {
   diff_Omega=-1.0;
   sign=-1.0;
} 
 else { 
       comp_omega();
       diff_Omega=Omega_const-Omega;
       sign=diff_Omega;
       d_ratio_Omega=fabs(diff_Omega)/Omega_const;
 }

 /* If Omega is greater than desired, reverse direction and cut stepsize
  | in half.
 */
 
 if(diff_Omega<0.0) {
   dr /= (-2.0);
 }

 while(d_ratio_Omega>omega_error && r_ratio<=1.0) {
      if(diff_Omega*sign<0.0) {
        sign=diff_Omega;
        dr /= (-2.0);
      }
      r_ratio -= dr;
      a_check=0.0;
      iterate(r_ratio);
      if(a_check==200) {
            diff_Omega = -1.0;
      }
      else { 
            comp_omega();      
            if(Omega>Omega_K) {
                  diff_Omega=-1.0;
            }else {
                  diff_Omega=Omega_const-Omega;
                  d_ratio_Omega=fabs(diff_Omega)/Omega_const;
            }     
      }
 } 
 comp();
 print_table();
}


/*************************************************************************/
/* Compute model for given J and e_center.                               */
/*************************************************************************/
void J_model(double J_const)
{
 double dr,
        diff_J,
        d_ratio_J;

 dr=0.1;
 r_ratio=1.0-dr;

 /* Compute first rotating model */

 make_center(e_center);
 guess();
 a_check=0.0;
 iterate(r_ratio);
 if(a_check==200) {
   diff_J=-1.0;
   sign=-1.0;
} 
 else { 
       comp_M_J();
       diff_J=J_const-J;
       sign=diff_J;
       d_ratio_J=fabs(diff_J)/J_const;
 }

 /* If J is greater than desired, reverse direction and cut stepsize
  | in half.
 */
 
 if(diff_J<0.0) {
   dr /= (-2.0);
 }

 while(d_ratio_J>J_error && r_ratio<=1.0) {
      if(diff_J*sign<0.0) {
        sign=diff_J;
        dr /= (-2.0);
      }
      r_ratio -= dr;
      a_check=0.0;
      iterate(r_ratio);
      if(a_check==200) {
            diff_J = -1.0;
      }
      else { 
            comp_M_J();      
            if(Omega>Omega_K) {
                  diff_J=-1.0;
            }else {
                  diff_J=J_const-J;
                  d_ratio_J=fabs(diff_J)/J_const;
            }     
      }
 } 
 comp();
 print_table();
}


/*************************************************************************/
/* Main program.                                                         */
/*************************************************************************/
main(int argc, char **argv)
{
 int i,
     j,
     s,
     m,
     n_of_models,
     task_option,
     t_check=0;

 double e_min,
        e_max,
        a,
        e_stat_1,
        e_stat_2,
        e_stat_int,
        r_ratio_last,
        M_fix;

 
 char task[30] = "ref";


   /* DEFAULT POLYTROPIC EOS PARAMETERS */

   n_P=1.0;

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

   for(i=1;i<argc;i++) 
     if(argv[i][0]=='-'){
       switch(argv[i][1]){


       case 'q':
                 sscanf(argv[i+1],"%s",eos_type);
			break;               
       case 'N':
                 sscanf(argv[i+1],"%lf",&n_P);
			break;               

       case 'f':
                 sscanf(argv[i+1],"%s",file_name);
    	         break;

       case 'v': fprintf(stderr,"rns version %s\n",version);
                 exit(1);
                 break;

       case 'h': 
                 fprintf(stderr,"\n");
                 fprintf(stderr,"Quick help for rns:\n");
                 fprintf(stderr,"\n");
                 fprintf(stderr,"  -q EOS type \n"); 
                 fprintf(stderr,"     tab : tabulated \n");
                 fprintf(stderr,"     poly : analytic polytropic \n");           
                 fprintf(stderr,"  -N polytropic index (P=K*e^(1+1/N))\n");      
                 fprintf(stderr,"  -f EOS file \n");
	         fprintf(stderr,"  -e central energy density in gr/cm^3\n");
	         fprintf(stderr,"  -l central energy density of last model\n");
	         fprintf(stderr,"  -r axes ratio\n");
	         fprintf(stderr,"  -m mass in solar masses\n");
	         fprintf(stderr,"  -z rest mass in solar masses\n");
	         fprintf(stderr,"  -o angular velocity in 10^4 s^-1\n");
	         fprintf(stderr,"  -j angular momentum in G*M_SUN^2/C\n");
	         fprintf(stderr,"  -t task to be performed:\n");
	         fprintf(stderr,"     model : requires -e and -r\n");
	         fprintf(stderr,"     gmass : requires -e and -m\n");
	         fprintf(stderr,"     rmass : requires -e and -z\n");
	         fprintf(stderr,"     omega : requires -e and -o\n");
	         fprintf(stderr,"     jmoment : requires -e and -j\n");
	         fprintf(stderr,"     static : requires -e\n");
	         fprintf(stderr,"     kepler : requires -e\n");
	         fprintf(stderr,"     test\n");                             
                 fprintf(stderr,"  -p printing option   (1)\n");
	         fprintf(stderr,"  -c relaxation factor   (1.0)\n");
          fprintf(stderr,"  -d if 0, do not monitor convergence   (1)\n");
	         fprintf(stderr,"  -a accuracy in convergence   (1e-5)\n");
          fprintf(stderr,"  -b accuracy in fixing M, Omega, etc.   (1e-4)\n");
	         fprintf(stderr,"  -n number of models, if sequence   (1)\n");
	         fprintf(stderr,"  -v print version number\n"); 
	         fprintf(stderr,"  -h this menu\n");
                 fprintf(stderr,"\n");
                 exit(1);
	         break;
     }
   }  


   make_grid();
   
   if(strcmp(eos_type,"tab")==0) {
     load_eos(file_name);

     if(strcmp(file_name,"eosCL_LOW")==0) {
       n0_match=1.0e38;
       e_cl=0.0; 
       de_pt=0.0;
     }

     e_min=2.66e15*C*C*KSCALE;
     e_max=1e16*C*C*KSCALE;
     r_ratio=0.75;
     M_fix=1.4*MSUN;
     M_0const=1.4*MSUN;
     Omega_const=1e-5;
     J_const=0.0;
   }
   else {
     e_min=0.34;
     e_max=0.5;
     r_ratio=0.58447265625;
     M_fix=0.1;
     M_0const=0.1;
     Omega_const=1e-5;
     J_const=0.0;
   }

   print_option=1;  
   cf=1.0;
   print_dif=1;
   accuracy=1e-5;    
   fix_error=1e-4;   
   n_of_models=1;
 
   comp_f_P();


   for(i=1;i<argc;i++) 
      if(argv[i][0]=='-'){
	switch(argv[i][1]){

              case 'e':
			sscanf(argv[i+1],"%lf",&e_min);
                        if(strcmp(eos_type,"tab")==0)
                          e_min *= C*C*KSCALE;
			break;
 
              case 'l':
			sscanf(argv[i+1],"%lf",&e_max);
                        if(strcmp(eos_type,"tab")==0)
                          e_max *= C*C*KSCALE;
			break;

              case 'r':
			sscanf(argv[i+1],"%lf",&r_ratio);
					break;
              case 'm':
			sscanf(argv[i+1],"%lf",&M_fix);
                        M_fix *= MSUN;
			break;
              case 'z':
			sscanf(argv[i+1],"%lf",&M_0const);
                        M_0const *= MSUN;
			break;
              case 'o':
			sscanf(argv[i+1],"%lf",&Omega_const);
                        Omega_const *= 1.0e4;
			break;
              case 'j':
			sscanf(argv[i+1],"%lf",&J_const);
                        J_const *= (G*MSUN*MSUN/C);
			break;
              case 'p':
			sscanf(argv[i+1],"%d",&print_option);
					break;
              case 'c':
			sscanf(argv[i+1],"%lf",&cf);
					break;
              case 'd':
			sscanf(argv[i+1],"%d",&print_dif);
					break;
              case 'a':
			sscanf(argv[i+1],"%lf",&accuracy);
					break;
              case 'b':
			sscanf(argv[i+1],"%lf",&fix_error);
					break;
              case 'n':
			sscanf(argv[i+1],"%d",&n_of_models);
					break;
              case 'u':
			sscanf(argv[i+1],"%lf",&n0_match);
                        n0_match *= 1e39;
		        break;

              case 'y':
			sscanf(argv[i+1],"%lf",&e_cl);
                        e_cl *= C*C*KSCALE;
		        break;
                   
        } 
      }  


 Gamma_P=1.0+1.0/n_P;

 if(strcmp(file_name,"eosCL_LOW")==0) {
   n_nearest=n_tab/2;
   e_match=pow(10.0,interp(log_n0_tab,log_e_tab,n_tab,log10(n0_match)));
   p_match=pow(10.0,interp(log_n0_tab,log_p_tab,n_tab,log10(n0_match)));
   h_match=pow(10.0,interp(log_n0_tab,log_h_tab,n_tab,log10(n0_match)));
   
   if(e_cl != 0.0) de_pt = e_cl - e_match;   
 }

 print_header();

 for(i=1;i<argc;i++) 
    if(argv[i][0]=='-'){
      switch(argv[i][1]){          
      case 't':
   	       sscanf(argv[i+1],"%s",task);
               task_option=task[0];
  
               switch(task_option) {
	        case 't' :                                    
                            e_center=e_min;
                            make_center(e_center);
                            guess();
                            iterate(r_ratio);
                            comp();
                            print_table();
                          
                	    break;

	        case 'm' :
                            if(n_of_models==1)
                               a=1.0;
                            else
                               a=pow(e_max/e_min,1.0/(n_of_models-1.0));
 
                            for(j=1;j<=n_of_models;j++) {
                               e_center=pow(a,1.0*j-1.0)*e_min;  
                               make_center(e_center);
                               guess();

                               if(r_ratio<0.8) 
                                 iterate(0.8);

                               if(r_ratio<0.6) 
                                 iterate(0.6);

                               iterate(r_ratio);
                               comp();
                               print_table();
			    }
                               break;               

	        case 's' :
                            r_ratio=1.0;

                            if(n_of_models==1)
                               a=1.0;
                            else
                               a=pow(e_max/e_min,1.0/(n_of_models-1.0));
 
                            for(j=1;j<=n_of_models;j++) {
                               e_center=pow(a,1.0*j-1.0)*e_min;  
                               make_center(e_center);
                               guess();
                               iterate(r_ratio);
                               comp();
                               print_table();
                            }
                            
                            break;
 
                case 'k' : 
                           omega_error=fix_error;

                           if(n_of_models==1)
                              a=1.0;
                           else
                              a=pow(e_max/e_min,1.0/(n_of_models-1.0));
 
                           for(j=1;j<=n_of_models;j++) {
                              e_center=pow(a,1.0*j-1.0)*e_min;  
                              ms_model();
                           }
                
                           break;

                case 'g' : 
                           M_error=fix_error;

                           if(n_of_models==1)
                              a=1.0;
                           else
                              a=pow(e_max/e_min,1.0/(n_of_models-1.0));
 
                           for(j=1;j<=n_of_models;j++) {
                              e_center=pow(a,1.0*j-1.0)*e_min;  
                              m_model(M_fix);
                           }
                
                           break;


                case 'r' : 
                           M_0_error=fix_error;

                           if(n_of_models==1)
                              a=1.0;
                           else
                              a=pow(e_max/e_min,1.0/(n_of_models-1.0));
 
                           for(j=1;j<=n_of_models;j++) {
                              e_center=pow(a,1.0*j-1.0)*e_min;  
                              m0_model(M_0const);
                           }
                
                           break;

 
                case 'o' : 
                           omega_error=fix_error;

                           if(n_of_models==1)
                              a=1.0;
                           else
                              a=pow(e_max/e_min,1.0/(n_of_models-1.0));
 
                           for(j=1;j<=n_of_models;j++) {
                              e_center=pow(a,1.0*j-1.0)*e_min;  
                              omega_model(Omega_const);
                           }

                           break;

                case 'j' : 
                           J_error=fix_error;

                           if(n_of_models==1)
                              a=1.0;
                           else
                              a=pow(e_max/e_min,1.0/(n_of_models-1.0));
 
                           for(j=1;j<=n_of_models;j++) {
                              e_center=pow(a,1.0*j-1.0)*e_min;  
                              J_model(J_const);
                           }

                           break;

	        case 'n' :                                    
                            e_center=e_min;
                            make_center(e_center);
                            guess();
                            printf("M=%8.7e  R=%8.7e  Z=%8.7e\n",Mass/MSUN,
                                                                R_e/1e5,Z_p);   
          
                      	    break;

                case 'h' : 
                           h_error=fix_error;

                           if(n_of_models==1)
                              a=1.0;
                           else
                              a=pow(e_max/e_min,1.0/(n_of_models-1.0));
 
                           for(j=1;j<=n_of_models;j++) {
                              e_center=pow(a,1.0*j-1.0)*e_min;  
                              h_model();
                           }
                
                           break;
               }
               
               t_check=1;  
               break;
        }
      }  

   if(t_check==0) {
     fprintf(stderr,"No task selected ! Use -t option.\n");
     exit(1);
   }

}
