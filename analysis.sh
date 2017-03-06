#! /bin/tcsh
module load gcc/5.2
cd RA_Guvenen/GKKC_Code/
make clean
rm -f log_analysis_sh_Hsieh_Klenow_Exp.txt
make GKK_Main > log_analysis_sh_Hsieh_Klenow_Exp.txt