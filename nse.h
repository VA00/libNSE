/* Interpolate NSE abundance at given triad:
(T9, rho [g/cc] , Ye), for subsequent nuclei defined in 
data/isotopes.dat starting with index ZERO (this is C, not Fortran)
*/
double NSE_enum(double, double, double, int     );


/* Interpolate NSE abundance at given triad:
(T9, rho [g/cc] , Ye), for nuclei with given Z, N
*/
double NSE(double, double, double, int, int);

void nse_(double *, double *, double *, double *, int *  );
