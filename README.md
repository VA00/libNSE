# libnse
A library for interpolating abundances of nuclear species in a state of nuclear statistical equilibrium (NSE).

Library provides the function calculating the abundance of a specified nuclide defined by Z nad N:
```c
nse(double T9, double rho, double Ye, int Z, int N)
```
where:
* T9 - temperature in 1E9 K (Kelvins)
* rho - density in [g/cc]  (grams per cubic centimeter)
* Ye  - electron fraction of the matter (ratio of lepton number to baryon number)
* Z - number of protons in the nuclei
* N - number of neutrons in the nuclei

### Prerequisites
* Linux operating system
* GCC
* autoconf
* libtool
* Mathematica [link](http://www.wolfram.com/mathematica/)

### Installation & Usage
In your console go to libNSE folder and type:
```sh
autoreconf -vif
./configure
mathematica src/NSE_pn_table_generator.nb
```
In the new console window type:
```sh
cd path_to_your_library/libNSE/src
gcc NSE_pn_table_parser.c -o NSE_pn_parser
./NSE_pn_table_parser > ../data/np_table.c
cd ../
make
```
