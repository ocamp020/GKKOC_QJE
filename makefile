
GKK_Main: GKK_Main.a
GKK_Main_Server: GKK_Main_Server.a
GKK_Opt_Taxes: GKK_Opt_Taxes.a
GKK_Wealth_Tax: Sergio.a
simple: Sergio_Simple.a

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# Folders to place objects, modules and executables
Folder       = ./Compiled_Files
GO_Folder    = ./Global_Opt
Objects_Main    = $(Folder)/NRTYPE.o $(Folder)/NRUTIL.o $(Folder)/Toolbox.o $(Folder)/parameters.o $(Folder)/global.o $(Folder)/programfunctions.o
Objects_Opt_Tax = $(Folder)/Opt_Tax_Parameters.o $(Folder)/Opt_Tax_Functions.o
Objects_GO      = $(Folder)/GKK_Calibration.o \
                  $(Folder)/stateControl.o $(Folder)/genericParams.o $(Folder)/utilities.o \
                  $(Folder)/simplex.o $(Folder)/objective.o $(Folder)/minimize.o
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# Compile all the modules in the main folder. Modules are saved as .f90 files. 
$(Folder)/%.o: %.F90
	gfortran -J$(Folder) -c $< -o $@

$(Folder)/%.o: %.f90
	gfortran -J$(Folder) -c $< -o $@

# Compile and execute programs
Sergio_Simple.a: GKK_simple.f95 $(Folder)/NRTYPE.o $(Folder)/NRUTIL.o
	gfortran -I$(Folder) GKK_simple.f95 $(Folder)/NRUTIL.o $(Folder)/NRTYPE.o -o $(Folder)/Sergio_Simple.a
	$(Folder)/Sergio_Simple.a

Sergio.a: GKK_Wealth_Tax_Sergio.f95 $(Folder)/NRTYPE.o $(Folder)/NRUTIL.o $(Folder)/Toolbox.o
	gfortran -I$(Folder) GKK_Wealth_Tax_Sergio.f95 $(Folder)/NRUTIL.o $(Folder)/NRTYPE.o $(Folder)/Toolbox.o -o $(Folder)/Sergio.a
	$(Folder)/Sergio.a

GKK_Main.a: GKK_Main.f95 $(Objects_Main)
	gfortran -I$(Folder) GKK_Main.f95 $(Objects_Main) -o $(Folder)/GKK_Main.a
	$(Folder)/GKK_Main.a

GKK_Main_Server.a: GKK_Main.f95 $(Objects_Main)
	gfortran -I$(Folder) GKK_Main.f95 $(Objects_Main) -o $(Folder)/GKK_Main.a

GKK_Opt_Taxes.a: GKK_Optimal_Taxes.f95 $(Objects_Main) $(Objects_Opt_Tax)
	gfortran -I$(Folder) GKK_Optimal_Taxes.f95 $(Objects_Main) $(Objects_Opt_Tax) -o $(Folder)/GKK_Opt_Taxes.a  

CE_program.a: Consumption_Equivalent.f95 $(Objects_Main)
	gfortran -I$(Folder) Consumption_Equivalent.f95 $(Objects_Main) -o $(Folder)/CE_program.a

GKK_Simul.a: GKK_Simul.f95 NRTYPE.o $(Objects_Main)
	gfortran -I$(Folder) GKK_Simul.f95 $(Objects_Main) -o $(Folder)/GKK_Simul.a
	$(Folder)//GKK_Simul.a

Aux.a: Aux_File.f95 $(Objects_Main)
	gfortran -I$(Folder) Aux_File.f95 $(Objects_Main) -o $(Folder)/Aux.a
	$(Folder)//Aux.a


#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# Compile all the modules for Global Optimisation. Modules are saved as .f90 files. 
$(Folder)/GKK_Calibration.o: GKK_Calibration.f95
	gfortran -J$(Folder) -c GKK_Calibration.f95 -o $(Folder)/GKK_Calibration.o

$(Folder)/%.o: $(GO_Folder)/%.f90
	gfortran -J$(Folder) -c $< -o $@

GO_Calibration.a: $(GO_Folder)/GlobalSearch.f95 $(Objects_Main) $(Objects_GO)
	gfortran -I$(Folder) $(GO_Folder)/GlobalSearch.f95 $(Objects_Main) $(Objects_GO) -o $(Folder)/GO_Calibration.a
	# cd Calibration
	# ../Compiled_Files/GO_Calibration.a 0 ../Global_Opt/config.txt b


