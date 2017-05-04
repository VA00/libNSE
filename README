Library provides two functions via nse.h header and libnse library:


nse(double T9, double rho, double Ye, int Z, int N) 

where:



T9 - temperature in 1E9 K (Kelvins) 
rho - density in [g/cc]  (grams per cubic centimeter)
Ye  - electron fraction of the matter (ratio of lepton number to baryon number) 

Z - number of protons in the nuclei
N - number of neutrons in the nuclei

(Z,N) pair uniquely identify nuclei, e.g.

Z=2, N=2 -> Helium4


If for some values Z and N nuclei do not exist, or is not included in our ensemble,
function returns 0.0 . In "normal" operation mode nse() returns equilibrium fraction X(N,Z) of given nuclei.
Fractions are normalized to 1.




If you have MATHEMATICA follow the steps to generate custom tabulations:

1)  Open notebook NSE_pn_table_generator.nb
2)  select number of isotopes from Wolfram database:
    default niso=800 with He2 deleted, but Li3 retained
3)  select range of temperatures in MeV, logarithm of densities in kg/m^3
    and electron fraction Ye:
    default: kT = (0.2,1.0) with 0.02 MeV step (41 entries)
             lg_rho = (5 , 15) with 0.2 step (51 entries)
             Ye = (0.40, 0.50) with 0.01 step (11 entries)
4)  select precision:
    default: $MachinePrecision

NOTE: this could take approx. 10 hours on modern PC, pre-calculated table
is included as files ... TODO
NOTE: generating nuclear partition function directly from
excited state spins and energies is very time-consuming, especially
at high software-emulated precisions, by default stored pre-calculated
tables are used
    
5) run notebook <<Evaluation/Evaluate Notebook>>     

6) notebook will generate include files for interpolation procedures
(skip this if you are not interested in using NSE abundances in your
C/Fortran/ etc. code):

a)
    A.dat -> mass numbers
    N.dat -> neutron numbers
    Z.dat -> proton numbers
    Q.dat -> binding energies (in MeV, for entire nuclei)
    G0.dat -> 2*J+1, where J is ground state spin
    GkT.dat -> temperature dependent partition function (NOTE: without G0 ! )
    N_Z_tbl.dat -> nuclide chart on N-Z plane, -1 means nuclei has not been included, 
    other number
    position in original tables, for use in NSE_enum function (see code nse.c)

b)    
    set of files:   np_<<niso>>_<<Ye>>.dat -> neutron and proton abundances
    with (kT,lg_rho) grid stored as arbitrary precision and rational numbers
    to avoid possible precision loss
    NOTE: these tables are not readable in usual programing languages (C, Fortran)
    as both rational and arbitrary precision types do not exist in these 
    languages, use converted files (see below)

c)  machine readable neutron and proton abundance tables:
      X_n_<<Ye>>.dat
      X_p_<<Ye>>.dat
    these files are filled with double precision floats, you can easily generate
    higher (or lower) precision files if you want this
d)  nse_tbl.h include file with parametrs of the generated p,n abundance tables

e)  two files with positions of the FFN nuclei in our NSE calculations (required
    for neutrino emissivities under NSE)        
  
7) Compile manually and run "./NSE_pn_table_parser >np_table.c" command to generate staticaly included at compile time in-code 
     proton and neutron abundance tables, included in file nse.c (function NSE_enum) 
8) compile test_nse.c and enjoy
