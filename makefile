
GKK_Main: GKK_Main.a
GKK_Main_Server: GKK_Main_Server.a
GKK_Opt_Taxes: GKK_Opt_Taxes.a
GKK_Wealth_Tax: Sergio.a
simple: Sergio_Simple.a
.PHONY: clean

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# Folders to place objects, modules and executables
Folder       = ./Compiled_Files_Baseline
GO_Folder    = ./Global_Opt
Objects_Main    = $(Folder)/NRTYPE.o $(Folder)/NRUTIL.o $(Folder)/Toolbox.o \
                  $(Folder)/parameters.o $(Folder)/global.o $(Folder)/programfunctions.o $(Folder)/Simulation_Module.o \
                  $(Folder)/GKK_Stats.o
Objects_Opt_Tax = $(Folder)/Opt_Tax_Parameters.o $(Folder)/Opt_Tax_Functions.o
Objects_GO      = $(Folder)/GKK_Calibration.o \
                  $(Folder)/stateControl.o $(Folder)/genericParams.o $(Folder)/utilities.o \
                  $(Folder)/simplex.o $(Folder)/objective.o $(Folder)/minimize.o

Flags    =  -fbounds-check
	#-mcmodel=large 
	#-fbounds-check -finline-functions
omp_flag = -fopenmp
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# Compile all the modules in the main folder. Modules are saved as .f90 files. M
$(Folder)/%.o: %.F90
	gfortran -O3 -J$(Folder) -c $< -o $@

$(Folder)/%.o: %.f90
	gfortran -O3 -J$(Folder) -c $< -o $@

$(Folder)/programfunctions.o: programfunctions.f90
	gfortran $(omp_flag) -O3 $(Flags) -J$(Folder) -c programfunctions.f90 -o $(Folder)/programfunctions.o

$(Folder)/Simulation_Module.o: Simulation_Module.f90
	gfortran $(omp_flag) -O3 $(Flags) -J$(Folder) -c Simulation_Module.f90 -o $(Folder)/Simulation_Module.o

$(Folder)/GKK_Stats.o: GKK_Stats.f90
	gfortran $(omp_flag) -O3 $(Flags) -J$(Folder) -c GKK_Stats.f90 -o $(Folder)/GKK_Stats.o

#??Compile and execute programs
Sergio_Simple.a: GKK_simple.f90 $(Folder)/NRTYPE.o $(Folder)/NRUTIL.o
	gfortran -I$(Folder) GKK_simple.f90 $(Folder)/NRUTIL.o $(Folder)/NRTYPE.o -o $(Folder)/Sergio_Simple.a
	$(Folder)/Sergio_Simple.a

Sergio.a: GKK_Wealth_Tax_Sergio.f90 $(Folder)/NRTYPE.o $(Folder)/NRUTIL.o $(Folder)/Toolbox.o
	gfortran -I$(Folder) GKK_Wealth_Tax_Sergio.f90 $(Folder)/NRUTIL.o $(Folder)/NRTYPE.o $(Folder)/Toolbox.o -o $(Folder)/Sergio.a
	$(Folder)/Sergio.a

GKK_Main.a: GKK_Main.f90 $(Objects_Main) $(Objects_Opt_Tax)
	gfortran $(omp_flag) $(Flags) -O3 -I$(Folder) GKK_Main.f90 $(Objects_Main) $(Objects_Opt_Tax) -o $(Folder)/GKK_Main.a
	time $(Folder)/GKK_Main.a

GKK_Main_pg.a: profiler.f90 $(Objects_Main)
	gfortran -pg $(omp_flag) -O3 -I$(Folder) profiler.f90 $(Objects_Main) -o $(Folder)/GKK_Main_pg.a
	time $(Folder)/GKK_Main_pg.a
	gprof $(Folder)/GKK_Main_pg.a

GKK_Main_Server.a: GKK_Main.f90 $(Objects_Main)
	gfortran -I$(Folder) GKK_Main.f90 $(Objects_Main) -o $(Folder)/GKK_Main.a

GKK_Opt_Taxes.a: GKK_Optimal_Taxes.f90 $(Objects_Main) $(Objects_Opt_Tax)
	gfortran $(omp_flag) -O3 -I$(Folder) GKK_Optimal_Taxes.f90 $(Objects_Main) $(Objects_Opt_Tax) -o $(Folder)/GKK_Opt_Taxes.a  
	time $(Folder)/GKK_Opt_Taxes.a

CE_program.a: Consumption_Equivalent.f90 $(Objects_Main)
	gfortran $(omp_flag) -I$(Folder) Consumption_Equivalent.f90 $(Objects_Main) -o $(Folder)/CE_program.a
	time $(Folder)/CE_program.a

GKK_Simul.a: GKK_Simul.f90 NRTYPE.o $(Objects_Main)
	gfortran -I$(Folder) GKK_Simul.f90 $(Objects_Main) -o $(Folder)/GKK_Simul.a
	$(Folder)//GKK_Simul.a

Aux.a: Aux_File.f90 $(Objects_Main)
	gfortran -I$(Folder) Aux_File.f90 $(Objects_Main) -o $(Folder)/Aux.a
	$(Folder)//Aux.a


#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# Compile all the modules for Global Optimisation. Modules are saved as .f90 files. 
$(Folder)/GKK_Calibration.o: GKK_Calibration.f90
	gfortran -J$(Folder) -c GKK_Calibration.f90 -o $(Folder)/GKK_Calibration.o

$(Folder)/%.o: $(GO_Folder)/%.f90
	gfortran -J$(Folder) -c $< -o $@

GO_Calibration.a: $(GO_Folder)/GlobalSearch.f90 $(Objects_Main) $(Objects_GO)
	gfortran -I$(Folder) $(GO_Folder)/GlobalSearch.f90 $(Objects_Main) $(Objects_GO) -o $(Folder)/GO_Calibration.a
	# cd Calibration
	# ../Compiled_Files/GO_Calibration.a 0 ../Global_Opt/config.txt b

GKK_Calibration.a: $(GO_Folder)/GlobalSearch.f90 $(Objects_Main) $(Objects_GO)
	gfortran -I$(Folder) $(GO_Folder)/GlobalSearch.f90 $(Objects_Main) $(Objects_GO) -o $(Folder)/GKK_Calibration.a
	cd Calibration
	rm -f log_calibration.txt
	nohup ../Compiled_Files/GKK_Calibration.a 0 ../Global_Opt/config.txt b >log_calibration.txt \
	& nohup ../Compiled_Files/GKK_Calibration.a 1 ../Global_Opt/config.txt b \
	& nohup ../Compiled_Files/GKK_Calibration.a 1 ../Global_Opt/config.txt b \
	& nohup ../Compiled_Files/GKK_Calibration.a 1 ../Global_Opt/config.txt b \
	& nohup ../Compiled_Files/GKK_Calibration.a 1 ../Global_Opt/config.txt b \
	& nohup ../Compiled_Files/GKK_Calibration.a 1 ../Global_Opt/config.txt b \



#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

GKK_Calibration_Loop.a: GKK_Calibration_Loop.f90 $(Objects_Main) $(Objects_Opt_Tax)
	gfortran $(omp_flag) -O3 -I$(Folder) GKK_Calibration_Loop.f90 $(Objects_Main) $(Objects_Opt_Tax) -o $(Folder)/GKK_Calibration_Loop.a

Calibration_Loop.a: Calibration_Loop.f90 $(Folder)/NRTYPE.o
	gfortran $(omp_flag) -O3 -I$(Folder) Calibration_Loop.f90 $(Folder)/NRTYPE.o -o $(Folder)/Calibration_Loop.a
	rm -f log_ind*
	time $(Folder)/Calibration_Loop.a


#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

Simulation_Labor_Income.a: Simulation_Labor_Income.f90 $(Objects_Main)
	gfortran $(omp_flag) -O3 -I$(Folder) Simulation_Labor_Income.f90 $(Objects_Main) -o $(Folder)/Simulation_Labor_Income.a
	time $(Folder)/Simulation_Labor_Income.a

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

clean: 
	cd $(Folder) ; pwd ; rm -f *.o *.mod *.a ; cd ..
