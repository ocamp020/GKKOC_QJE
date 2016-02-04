# To run the global optimization algorithm
#
# 1. Run "source compile.sh" 
# 2. Edit config.txt with the appropriate options
# 3. Follow Arun's instructions on config.txtR
#    For a cold start using Nelder Mead I use
#    ./GlobalSearch 0 config.txt a
#    Change to a for b for BOBYQ
#
# Set up environment variables for Mac and Intel Fortan 11.0
# change accordingly to the folder where teh ifortvars.sh is
source /opt/intel/bin/ifortvars.sh intel64
# Remove old libraries
rm *.o *.mod
# Compile
ifort -m64 -g -debug all -traceback -fp-stack-check -heap-arrays nrtype.f90 stateControl.f90 genericParams.f90 simplex.f90 global.f90 Globals.f90 random.f90 prctile.f90 utilities.f90 objective.f90 minimize.f90 GlobalSearch.f90 -o GlobalSearch -llapack -lblas 
