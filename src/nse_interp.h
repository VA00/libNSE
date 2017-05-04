double G_linear(double, int);

#if HAVE_LIBGSL
double G_gsl_linear(double, int);
#endif

double X_NSE(double, double, double, double, 
             double, double, double, double, double);


double NSE_enum_trilinear(double, double, double, int);


