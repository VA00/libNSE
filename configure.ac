#    -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.65)
AC_INIT([libnse],[0.14],[andrzej.odrzywolek@uj.edu.pl],[],[2024-11-12])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror gnu])
AC_CONFIG_MACRO_DIR([m4])
AM_PROG_AR
        
AC_CONFIG_FILES([                
Makefile
src/Makefile
])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
LT_INIT


# Checks for libraries.
AC_CHECK_LIB([m], [pow, sqrt, log10, exp])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt sscanf])

# Checks for header files.

AC_CHECK_HEADERS([stdio.h stdlib.h string.h unistd.h math.h gsl/gsl_errno.h gsl/gsl_spline.h gsl/gsl_interp.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_HEADER_STDBOOL


AC_OUTPUT
