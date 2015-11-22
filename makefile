
GKK_Main: GKK_Main.a
GKK_Main_Server: GKK_Main_Server.a
GKK_Opt_Taxes: GKK_Opt_Taxes.a
GKK_Wealth_Tax: Sergio.a
simple: Sergio_Simple.a

# Folders to place objects, modules and executables
Folder       = ./Compiled_Files
Objects_Main = $(Folder)/NRTYPE.o $(Folder)/NRUTIL.o $(Folder)/Toolbox.o $(Folder)/parameters.o $(Folder)/global.o $(Folder)/programfunctions.o
Objects_Opt_Tax = $(Folder)/Opt_Tax_Parameters.o $(Folder)/Opt_Tax_Functions.o

# The following lines compile all the modules. Modules are saved as .f90 or .F90 files. 
$(Folder)/%.o: %.F90
	gfortran -J$(Folder) -c $< -o $@

$(Folder)/%.o: %.f90
	gfortran -J$(Folder) -c $< -o $@

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