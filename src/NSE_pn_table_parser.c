/*
Written by A. Odrzywolek
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "../data/nse_tbl.h"

int main( int argc, char** argv)
{
  int i,j,k,inuc;
  double tbl_n[N_Ye_nse][N_kT_nse][N_lg_rho_nse];
  double tbl_p[N_Ye_nse][N_kT_nse][N_lg_rho_nse];
  FILE *datafile_Xn,*datafile_Xp;
  double T, lg_rho, Ye, Xn, Xp, kT, Y_e;
  char string[810],FILENAMEn[2000],FILENAMEp[2000];
  double Ye_tbl[N_Ye_nse] = 

  #include "../data/Ye.dat"
  
  ;

  for(k=0;k<N_Ye_nse;k++){

    Ye = Ye_tbl[k];
    sprintf(FILENAMEn,"../data/X_n_%.4lf.dat",Ye);
    sprintf(FILENAMEp,"../data/X_p_%.4lf.dat",Ye);
    //printf("%s\n",FILENAMEp);

    datafile_Xn = fopen(FILENAMEn,"r");
    datafile_Xp = fopen(FILENAMEp,"r");
    if(datafile_Xn==NULL) { fprintf(stderr,"NSE data not found\n");exit(1);}
    if(datafile_Xp==NULL) { fprintf(stderr, "NSE data not found\n");exit(1);}


    for(i=0;i<N_kT_nse;i++){
      for(j=0;j<N_lg_rho_nse;j++){
        fgets(string ,810, datafile_Xn);
        sscanf(string,"%lf %lf %lf", &T, &lg_rho, &Xn);
        tbl_n[k][i][j]=Xn; 

        fgets(string ,810, datafile_Xp);
        sscanf(string,"%lf %lf %lf", &T, &lg_rho, &Xp);
        tbl_p[k][i][j]=Xp; 
      }
    }

    fclose(datafile_Xn);   
    fclose(datafile_Xp);   
    
  }

/*
This part generates
static initialized C code
*/
#if 1

  printf("static const double n[%d][%d][%d]={\n",N_Ye_nse,N_kT_nse,N_lg_rho_nse);

  for(k=0;k<N_Ye_nse;k++){
    printf("{");
    for(i=0;i<N_kT_nse;i++){
      printf("{");
      for(j=0;j<N_lg_rho_nse;j++){
        printf("%.16e",tbl_n[k][i][j]); 
        if(j!=N_lg_rho_nse-1) printf(", ");
      }
      printf("}");
      if(i!=N_kT_nse-1) printf(", \n");
    }
    printf("}");
    if(k!=N_Ye_nse-1) printf(", ");
  }
  printf("\n};\n\n");  


  printf("static const double p[%d][%d][%d]={\n",N_Ye_nse,N_kT_nse,N_lg_rho_nse);

  for(k=0;k<N_Ye_nse;k++){
    printf("{");
    for(i=0;i<N_kT_nse;i++){
      printf("{");
      for(j=0;j<N_lg_rho_nse;j++){
        printf("%.16e",tbl_p[k][i][j]); 
        if(j!=N_lg_rho_nse-1) printf(", ");
      }
      printf("}");
      if(i!=N_kT_nse-1) printf(", \n");
    }
    printf("}");
    if(k!=N_Ye_nse-1) printf(", ");
  }
  printf("\n};\n\n");  

#endif

/*
This part generates
LateX table
*/
#if 0



// for(k=0;k<N_Ye_nse;k=k+1){
  for(k=7;k<27;k=k+5){
    for(i=0;i<N_kT_nse;i=i+1){
      for(j=0;j<N_lg_rho_nse;j=j+1){
        kT    = kT_MIN_nse+i*delta_kT_nse;
        lg_rho = lg_rho_MIN_nse + j*delta_lg_rho_nse;
        Y_e=Ye_tbl[k];
        printf("%.2lf&\t%.0lf&\t%.3lf&\t%.16e&\t%.16e\\\\\n",
        kT, lg_rho, Y_e,
        tbl_p[k][i][j],tbl_n[k][i][j]); 
      }
    }
  }




#endif


#if 0


/* FULL TABLE */
  for(k=0;k<N_Ye_nse;k=k+1){
    for(i=0;i<N_kT_nse;i=i+1){
      for(j=0;j<N_lg_rho_nse;j=j+1){
        kT    = kT_MIN_nse+i*delta_kT_nse;
        lg_rho = lg_rho_MIN_nse + j*delta_lg_rho_nse;
        Y_e=Ye_tbl[k];
        printf("%.2lf\t%.0lf\t%.3lf\t%.16e\t%.16e\n",
        kT, lg_rho, Y_e,
        tbl_p[k][i][j],tbl_n[k][i][j]); 
      }
    }
  }




#endif


/*
Partition function txt  table generator


for(inuc=0;inuc<niso;inuc++)
 {
  for(i=0;i<N_kT_nse;i++) printf("%.6lf\t",GkT[inuc][i]);
  printf("\n");
 } 

printf("\n");

*/





return 0;
}


