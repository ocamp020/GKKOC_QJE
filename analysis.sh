#! /bin/tcsh
module load gcc/4.8.5
cd RA_Guvenen/GKKC_Code/
make clean
rm -f log_timing_tax_reform.txt
make GKK_Main > log_timing_tax_reform.txt