# Use it or Lose it: Efficiency and Redistributional Effects of Wealth Taxation
#### **Fatih Guvenen & Gueorgui Kambourov & Burhan Kuruscu & Sergio Ocampo & Daphne Chen**

#### **Replication Files**

#### **Contact:** socampod@uwo.ca

<br/>
<br/>

---
### **Overview**

This repository contains code that implements the model of Guvenen, Kambourov, Krusucu, Ocampo & Chen (QJE). 
The code is in Fortran and compiled with GCC gfortran. 
Graphs and tables are generated in Matlab. 
The paper is located in the main folder or by [clicking here](https://github.com/ocamp020/GKKOC_QJE/blob/Model_2.1/GKKOC_2022_QJE.pdf).
A computational appendix is also located in the main folder or by [clicking here](https://github.com/ocamp020/GKKOC_QJE/blob/Model_2.1/GKKOC_Computational_Appendix.pdf).

The repository is organized by branches, with the main branch (Model_2.1) corresponding to the baseline model. 
Other branches correspond to extensions and robustness exercises. 
See below for the description of the main files and the form of the output. 

The Fortran code generates output files and simulated cross-sections of agents that are later used for tables and figures. 
Most tables are generated directly as part of the Fortran output, but some are generated in Matlab. 
All figures are generated in Matlab using output from Fortran. 
See the description of Graphs and Tables below for more details. 


---
### **Navigating the Folders**

The main folder contains this readme file and a copy of the paper and the computational appendix.

The main folder also contains all Fortran code (described below). 
The main files are the makefile and GKK_Main.f90 that control the execution of the code. 

The Graphs folder contains the Matlab code used to generate the figures and tables (described below).
It also contains the auxiliary data required for certain figures.


---
### **Reproducig Baseline Results**


The execution of the code is controlled from GKK_Main.f90. 
There is a set a logical flags defined at the beginning of the code that define what output is to be generated. 

Flags:
    * Tax_Reform: Defines if the tax reform experiment is being performed 
        * compute_bench: Defines if the benchmark model must be solved. If .false. the code loads the result from a previous run. 
        * compute_exp: Defines if the tax reform must be solved. If .false. the code loads results from a previous run. 
    * Tax_Reform_tktw: Defines if the tax reform experiment that adds wealth taxes on top of the baseline tax system is being performed.
		* budget_flag: Defines whether the additional revenue is being rebated by lowering labor income taxes
    * Opt_tax: Defines if the optimal tax experiments are being performed
        * Opt_Tax_KW: Defines whether OWT or OKIT is being performed. true=tau_K, false=tau_W.
    * Opt_Threshold: Defines if the optimal wealth tax with threshold is being performed. 
    * Transition_Tax_Reform: Defines if the transition for the tax reform is being performed 
        * Transition_OT: Defines if the transition for the optimal tax is being performed.
        * budget_balance: Defined whether we are balancing the budget while computing the transition. 
        * balance_tau_L: Defines which tax is being used to balance the budget. true=tau_L, false=tau_K or tau_W depending on Opt_Tax_KW.
    * Simul_Switch: Defines if we are simulating the model after computing the solution. The simulation provides a cross-section of 20M agents.



The parameters are declared in parameters.f90 and global.f90. 
The file programfunctions.f90 contains all the functions used for the solution of the model.
The file GKK_Stats.f90 produces the outptu files from the model solution (including welfare gains and the output.txt file used for tables).
The files Opt_Tax_Functions.f90 and Opt_Tax_Parameters.f90 contain additional functions and parameters used for the optimal tax experiments. 
The file Simulation_Module.f90 contains the routines used for simulating the model after computing the solution. 

All the code requires using the following toolbox: Toolbox.f90, NRTYPE.f90, NRUTIL.f90. 

The execution of the code uses the makefile to ensure that all modules are loaded. 
Use "make GKK_Main" to run the code after specifying all the flags. 
Before running the code for the first time create a folder called "Compiled_Files_Baseline" where the executable files are stored.
All code is compiled using gfortran. 

The code creates folders to store the results. 
The main folder is called Model_2.1. 
Inside it there will be a simulation folder (if Simul_Switch=.true.), a Bench_Files folder for the storage of the policy functions, value functions and distribution. 
Additional folders will be created for storing the results of the Tax_Reform and Opt_Tax_W, and so on. 
See the graphs and tables file for a description of the files containing the information presented in the paper. 


---
### **Figures and Tables**

The Figures_Tables folder contains a Matlab file that generates all the figures of the paper from the outptu of the model. 
The file also generates the tables that are not automatically generated by the Fortran code. 
Running this file requires running the baseline model as well as the relevant extensions. 
The folder also contains the necessary files for U.S. wealth data and intergenerational correlation in Norway. 

---
### **Extensions**

The repository is organized by branches. 
Each branch corresponds to a version of the model (an extension or robustness).
The main (default) branch is Model_2.1 and corresponds to the baseline model. 
When switching to other branches the code updates to run an alternative version of the model and save the results in its corresponding folder. 
Folders have the same name as the branch. 

Branches and Expensions:

    * Low-Inequality Calibration: Match_Return_Lambda
    * Awesome State Income Shocks: Super_Aiyagari  
    * Credit Spread (6%): Model_2.1_RW_Debt_1.5
    * Credit Spread (10%): Model_2.1_RW_Debt_2.0
    * Public Firms: Model_2.1_IPO and Model_2.1_IPO_High_Entry
    * Corporate Sector: Model_2.1_Corp_CD_Large_Corp
    * Pure Rents Model: Model_2.1_CKK_Markups
    * Non-linear OKIT: Model_2.1_NLKT
    * Looser Constraints (Debt/GDP=2.5): Model_2.1_Debt_3
    * Constant Constraints: Model_2.1_Constant_Vartheta
    * Higher Markups (mu=0.8): Model_2.1_mu80
    * Constant Productivity (z_ih=z_i): Model_2.1_Constant_Z
    * Individuals start in the normal lane and can move to high lane productivity: Model_2.1_X_Transition_Low_NB

