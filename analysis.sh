#! /bin/tcsh
module load gcc/5.2
cd RA_Guvenen/GKKC_Code/
make clean
rm -f log_analysis_sh_opt_tax_w_ss.txt
make GKK_Main > log_analysis_sh_opt_tax_w_ss.txt