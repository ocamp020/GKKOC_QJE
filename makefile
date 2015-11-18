
GKK_Main: GKK_Main.a
GKK_Main_Server: GKK_Main_Server.a
GKK_Opt_Taxes: GKK_Opt_Taxes.a
GKK_Wealth_Tax: Sergio.a
simple: Sergio_Simple.a
simple_2: Sergio_Simple_2.a

# The following lines compile all the modules. Modules are saved as .f90 or .F90 files. 
%.o: %.F90
	gfortran -c $<

%.o: %.f90
	gfortran -c $<

Sergio_Simple_2.a: GKK_simple_V2.f95 NRTYPE.o NRUTIL.o
	gfortran GKK_simple_V2.f95 NRUTIL.o NRTYPE.o -o Sergio_Simple_2.a
	./Sergio_Simple_2.a

Sergio_Simple.a: GKK_simple.f95 NRTYPE.o NRUTIL.o
	gfortran GKK_simple.f95 NRUTIL.o NRTYPE.o -o Sergio_Simple.a
	./Sergio_Simple.a

Sergio.a: GKK_Wealth_Tax_Sergio.f95 NRTYPE.o NRUTIL.o Toolbox.o
	gfortran GKK_Wealth_Tax_Sergio.f95 NRUTIL.o NRTYPE.o Toolbox.o -o Sergio.a
	./Sergio.a

GKK_Main.a: GKK_Main.f95 NRUTIL.o NRTYPE.o Toolbox.o parameters.o global.o programfunctions.o
	gfortran GKK_Main.f95 NRUTIL.o NRTYPE.o Toolbox.o parameters.o global.o programfunctions.o -o GKK_Main.a
	./GKK_Main.a

GKK_Main_Server.a: GKK_Main.f95 NRUTIL.o NRTYPE.o Toolbox.o parameters.o global.o programfunctions.o
	gfortran GKK_Main.f95 NRUTIL.o NRTYPE.o Toolbox.o parameters.o global.o programfunctions.o -o GKK_Main_Server.a

GKK_Opt_Taxes.a: GKK_Optimal_Taxes.f95 NRUTIL.o NRTYPE.o Toolbox.o Opt_Tax_Parameters.o Opt_Tax_Functions.o  parameters.o global.o programfunctions.o 
	gfortran GKK_Optimal_Taxes.f95 NRUTIL.o NRTYPE.o Toolbox.o Opt_Tax_Parameters.o Opt_Tax_Functions.o parameters.o global.o programfunctions.o -o GKK_Opt_Taxes.a  

CE_program.a: Consumption_Equivalent.f95 NRUTIL.o NRTYPE.o Toolbox.o parameters.o global.o programfunctions.o
	gfortran Consumption_Equivalent.f95 NRUTIL.o NRTYPE.o Toolbox.o parameters.o global.o programfunctions.o -o CE_program.a

GKK_Simul.a: GKK_Simul.f95 NRUTIL.o NRTYPE.o Toolbox.o parameters.o global.o programfunctions.o
	gfortran GKK_Simul.f95 NRUTIL.o NRTYPE.o Toolbox.o parameters.o global.o programfunctions.o -o GKK_Simul.a
	./GKK_Simul.a

Aux.a: Aux_File.f95 NRUTIL.o NRTYPE.o Toolbox.o parameters.o global.o programfunctions.o
	gfortran Aux_File.f95 NRUTIL.o NRTYPE.o Toolbox.o parameters.o global.o programfunctions.o -o Aux.a
	./Aux.a