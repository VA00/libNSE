#include "../nse.h"
#include "../config.h"
#include "nse_interp.h"
#include "../data/nse_tbl.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if HAVE_LIBGSL
 #include <gsl/gsl_errno.h>
 #include <gsl/gsl_spline.h>
 #include <gsl/gsl_interp.h>
#endif

#define kT_to_T9 11.6045

  const double G0[niso]=  
  
  #include "../data/G0.dat"
  
  ;

  /* mass numbers */
  const double A[niso]=

  #include "../data/A.dat"

  ;

  const double Z[niso]=

  #include "../data/Z.dat"

  ;

  const double N[niso]=

  #include "../data/N.dat"

  ;

  const double Q[niso]=

  #include "../data/Q.dat"

  ;

  const double Ye_tbl[N_Ye_nse] = 
    #include "../data/Ye.dat"
  ;
  const double GkT[niso][N_kT_nse]={
    #include "../data/GkT.dat"
  };
  const double kT_tbl[N_kT_nse] = 
    #include "../data/kT.dat"
  ;
  
  const double lg_rho_tbl[N_lg_rho_nse] = 
    #include "../data/lg_rho.dat"
  ;

#include "../data/np_table.c"

inline double
trilinear_interp_unit_cube
(double X, double Y, double Z, double A000, double A100, double A101, double A001, double A010, double A110, double A111, double A011)
{

  return A000*(1.0-X)*(1.0-Y)*(1.0-Z) + A100*X*(1.0-Y)*(1.0-Z) + A010*(1.0-X)*Y*(1.0-Z) + A001*(1.0-X)*(1.0-Y)*Z 
       + A101* X*(1.0 - Y)*Z + A011* (1.0 - X)*Y*Z + A110* X* Y* (1.0 - Z) + A111* X*Y*Z;
}

double NSE(double T9, double rho, double Ye, int Z, int N)
{
  /* table consist of integer enumerator of subsequent nuclides,
     this is not really needed, in real application you
     probably will use this enumerator directly.
     Anyway, it is nice, natural and portable to select nuclides by Z and N ! 
  */
  const int nuclide_tbl[Z_max_nse][N_max_nse]=
    #include "../data/N_Z_tbl.dat"
  ;
  int inuc=-1;
  
  if(Z>=Z_max_nse) return 0.0;
  if(N>=N_max_nse) return 0.0;

  inuc = nuclide_tbl[Z][N];
  /* if given nuclide is not included in NSE
     calculations (inuc==-1), we return abundance equal to 0.0 */
  if(inuc==-1) return 0.0;

  return NSE_enum(T9, rho, Ye, inuc);

}





double NSE_enum_trilinear(double T9, double rho, double Ye, int inuc)
{
  /* temp.  dependent part of partition function*/
  double GT_A, GT_B, GT_C, GT_D, GT_E, GT_F, GT_G, GT_H;//8 corners of cube for trilinear interpolation
  double kT, lg_rho;
  /* trilinear interpolator variables */
  int i,j,k;
  double kT_lower,kT_upper,kT_A,kT_B, kT_C, kT_D, kT_E,kT_F, kT_G, kT_H;
  double lg_rho_lower, lg_rho_upper, lg_rho_A, lg_rho_B, lg_rho_C, lg_rho_D, lg_rho_E, lg_rho_F, lg_rho_G, lg_rho_H;
  double Ye_lower, Ye_upper;
  double n_A, n_B, n_C, n_D, n_E, n_F, n_G, n_H;
  double p_A, p_B, p_C, p_D, p_E, p_F, p_G, p_H; 
  double X_A, X_B, X_C, X_D, X_E, X_F, X_G, X_H; 
  double X_kT, X_lg_rho, X_Ye;
  double X_interp;
  
  /*
  For non-equidistant grid you have to handle all three vars
  like Ye. To speed-up selection of the point, you might
  try to use analytical estimate (limes inferior) for
  Ye(k) function. 
  TODO: do this!
  */

  kT = T9/kT_to_T9;
  lg_rho = log10(rho);
  
  /*
  i = (int) ceil( (kT-kT_MIN_nse)/delta_kT_nse);
  j = (int) ceil( (lg_rho-lg_rho_MIN_nse)/delta_lg_rho_nse);
  
  if( i<1 ) i=1; if( i>=N_kT_nse     ) i=N_kT_nse-1;
  if( j<1 ) j=1; if( j>=N_lg_rho_nse ) j=N_lg_rho_nse-1;
  */
  
  i=0;
  while( (kT_tbl[i]<kT) && (i<N_kT_nse-1) ) i++;
  kT_lower = kT_tbl[i-1]; kT_upper = kT_tbl[i];
  
  j=0;
  while( (lg_rho_tbl[j]<lg_rho) && (j<N_lg_rho_nse-1) ) j++;
  lg_rho_lower = lg_rho_tbl[j-1]; lg_rho_upper = lg_rho_tbl[j];
  
  k=0;
  while( (Ye_tbl[k]<Ye) && (k<N_Ye_nse-1) ) k++;
  Ye_lower=Ye_tbl[k-1]; Ye_upper=Ye_tbl[k];
  

  /* 8 corners, starting with upper right are: 
  A, B, C, D (first Ye plane), 
  E, F, G, H (second Ye plane) */
  kT_A = kT_lower;
  kT_B = kT_upper;
  kT_C = kT_upper;
  kT_D = kT_lower;
  kT_E = kT_lower;
  kT_F = kT_upper;
  kT_G = kT_upper;
  kT_H = kT_lower;
  

  
  #if !HAVE_LIBGSL
  GT_A = G0[inuc] + G_linear(kT_A,inuc);
  GT_B = G0[inuc] + G_linear(kT_B,inuc);
  #endif
  #if HAVE_LIBGSL
  GT_A = G0[inuc] + G_gsl_linear(kT_A,inuc);
  GT_B = G0[inuc] + G_gsl_linear(kT_B,inuc);
  #endif
  
  GT_C = GT_B;
  GT_D = GT_A;
  GT_E = GT_A;
  GT_F = GT_B;
  GT_G = GT_C;
  GT_H = GT_D;
  
  
  lg_rho_A = lg_rho_lower;
  lg_rho_B = lg_rho_lower;
  lg_rho_C = lg_rho_upper;
  lg_rho_D = lg_rho_upper;
  lg_rho_E = lg_rho_lower;
  lg_rho_F = lg_rho_lower;
  lg_rho_G = lg_rho_upper;
  lg_rho_H = lg_rho_upper;
  

  
  n_A = n[k-1][i-1][j-1];
  n_B = n[k-1][i][j-1]; 
  n_C = n[k-1][i][j];
  n_D = n[k-1][i-1][j];
  n_E = n[k][i-1][j-1];  
  n_F = n[k][i][j-1];
  n_G = n[k][i][j];
  n_H = n[k][i-1][j];
  
  
  p_A = p[k-1][i-1][j-1];
  p_B = p[k-1][i][j-1];
  p_C = p[k-1][i][j];
  p_D = p[k-1][i-1][j];
  p_E = p[k][i-1][j-1];
  p_F = p[k][i][j-1];
  p_G = p[k][i][j];
  p_H = p[k][i-1][j];
  
    
  X_A = X_NSE(kT_A, lg_rho_A, n_A, p_A, N[inuc], Z[inuc], A[inuc], GT_A, Q[inuc]);  
  X_B = X_NSE(kT_B, lg_rho_B, n_B, p_B, N[inuc], Z[inuc], A[inuc], GT_B, Q[inuc]);  
  X_C = X_NSE(kT_C, lg_rho_C, n_C, p_C, N[inuc], Z[inuc], A[inuc], GT_C, Q[inuc]);  
  X_D = X_NSE(kT_D, lg_rho_D, n_D, p_D, N[inuc], Z[inuc], A[inuc], GT_D, Q[inuc]);  
  X_E = X_NSE(kT_E, lg_rho_E, n_E, p_E, N[inuc], Z[inuc], A[inuc], GT_E, Q[inuc]);  
  X_F = X_NSE(kT_F, lg_rho_F, n_F, p_F, N[inuc], Z[inuc], A[inuc], GT_F, Q[inuc]);  
  X_G = X_NSE(kT_G, lg_rho_G, n_G, p_G, N[inuc], Z[inuc], A[inuc], GT_G, Q[inuc]);  
  X_H = X_NSE(kT_H, lg_rho_H, n_H, p_H, N[inuc], Z[inuc], A[inuc], GT_H, Q[inuc]);  

  //rescallig kT, lg_rho, Ye cell to unit cube 

  
  X_kT     = (kT-kT_lower)/(kT_upper-kT_lower);
  X_lg_rho = (lg_rho-lg_rho_lower)/(lg_rho_upper-lg_rho_lower);
  X_Ye     = (Ye-Ye_lower)/(Ye_upper-Ye_lower);

  
  X_interp = trilinear_interp_unit_cube(X_kT,X_Ye,X_lg_rho,
          X_A,X_B,X_C,X_D,X_E,X_F,X_G,X_H);

           return X_interp;
}


/* T9 - temperature in Kelvins divided by 1.0E09,
  rho - density in [g/cm^3] (WARNING! previously in [kg/m^3], tables STILL in [kg/m^3] !!!
  Ye - electron lepton number fraction ; number of ,,electrons'' (electrons minus positrons)  
       divided by number of baryons
  inuc - nuclide number inuc
*/

double NSE_enum(double T9, double rho, double Ye, int inuc)
 {
   
  if(inuc>niso) return 0.0;
 
  return NSE_enum_trilinear(T9, rho, Ye, inuc); 

 } 

double X_NSE(double kT, double lg_rho, double Xn, double Xp, 
             double N, double Z, double A, double G, double Q)
{
  double lambda, lambda_3, factor, log_X;
  const double mp=1.672E-27;
  double rho;
  
  if(A==1.0)
  {
    if(Z==1) return Xp;
    if(N==1) return Xn;
  }  

  rho = pow(10.0, lg_rho);

  lambda = 1.6150947165848602E-14/sqrt(kT);

  lambda_3 = lambda*lambda*lambda;

  factor =  1000.0*rho/(2.0*mp) * lambda_3; // rho initially was in kg/m^3, now in g/cc 2010-07-29
 
  log_X =log(0.5*G)
     +
     (A - 1.0)*log(factor)
     +
     2.5*log(A)
     +
     log(Xn)*N
     +
     log(Xp)*Z
     +
     Q/kT;


  return exp(log_X);


} 

/* calculates temperature dependent partition function for inuc-th nuclei based on tabulated
 Sum[exp(-E_i/kT] factor (GkT.dat) NOTE: WITHOUT ground state G=2*J0+1*/
#if 0  

Unknown bug in function G_linear below forced to remove it !!!
A.Odrzywolek, 2015-09-07

double G_linear(double kT, int inuc){

  double kT_1,kT_2, G_1, G_2;
  int i;

  if(kT<=kT_MIN_nse) return 0.0;
  if(kT>=kT_MAX_nse) return GkT[inuc][N_kT_nse-1];

  if(GkT[inuc][N_kT_nse-1]<0.000001) return 0.0;

  /* select upper right grid point */
  i=0;
  while( (kT<kT_tbl[i]) && (i<N_kT_nse-1) ) i++;

  /* neighbouring values */

  G_1 =  GkT[inuc][i-1];
  G_2 =  GkT[inuc][i];

  kT_1 = kT_tbl[i-1];
  kT_2 = kT_tbl[i];


  return G_1 + (G_2-G_1)/(kT_2-kT_1)*(kT-kT_1);
  

}
#endif

/* NOTE: function below do not include G0 ! It is added in NSE_enum_trilinear */
#if HAVE_LIBGSL
double G_gsl_linear(double kT, int inuc){

  double G;

  
 
    
  if(kT<=kT_MIN_nse) return 0.0;
  if(kT>=kT_MAX_nse) return GkT[inuc][N_kT_nse-1];
  
  gsl_interp_accel *acc     = gsl_interp_accel_alloc ();
  gsl_interp *lin     = gsl_interp_alloc (gsl_interp_linear, N_kT_nse);

  gsl_interp_init (lin, kT_tbl, GkT[inuc], N_kT_nse);

  G = gsl_interp_eval (lin, kT_tbl, GkT[inuc], kT, acc);

  gsl_interp_free (lin);
  gsl_interp_accel_free (acc);

  return G;


}
#endif


/* Fortran wrapper function. Usage:
call nse(X,T9,rho,Ye,inuc)
*/
void nse_(double *X, double *T9, double *rho, double *Ye, int *inuc) {
  
  double Xtmp;
  
  Xtmp = NSE_enum((*T9),(*rho),(*Ye),(*inuc));

  (*X)=Xtmp;
}
