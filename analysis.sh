#! /bin/tcsh
module load gcc/4.8.5
cd RA_Guvenen/GKKC_Code/
make clean
rm -f log_Transition.txt
make GKK_Main > log_Transition.txt