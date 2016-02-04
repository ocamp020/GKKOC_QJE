-----------------------------------
0. Contents
-----------------------------------
1 - About the program
2 - Executing the program
3 - Description of source files
4 - Description of text files
5 - Description of .dat files
6 - Specifying the objective function
7 - Troubleshooting

-----------------------------------
1. About the program
-----------------------------------
This code is based on the work by Fatih Guvenen and Tony Smith. Great
care was taken to make it as compliant with Fortran 90 as possible, but
there are a couple of invocations to Fortran 95 intrinsics (specifically
to timing variables such as SYSTEM_CLOCK)

-----------------------------------
2. Executing the program
-----------------------------------
To execute the program, run 
	./GlobalSearch <-1|0|1|2> configfile <a|b|d>
For help, run
	./GlobalSearch

 0 = cold start - The first invocation of the program, that will set up all the parallel helper files.
             Should always be the parameter when this is the first attempt at solving the problem.
 1 = warm start - after one process is already running, all helper programs should be invoked using
             a warm start.
 2 = update points - update the sobol point parameters over which to search, but assumes everything
             else in the config file has not been changed.
             
Note, amoeba and bobyq search have been implemented to date.

-----------------------------------
3. Description of source files
-----------------------------------
These files are specific for the generic search
GlobalSearch.f90  - the main driver program for the search.
genericParams.f90 - the parameters that the generic search program needs. Note that
					we do not put function specific parameters in this file
minimize.f90      - this module contains the code for minimization. 
nrtype.f90        - basic types used in all functions.
simplex.f90       - open source code that obtains an m-dimensional simplex centered
                    on the origin. Used for amoeba search
stateControl.f90  - module that manages the genericSearch states using file I/O.
utilities.f90     - implementation of sobol and other helper functions.
---------------------------------
These are all specific to the value function being solved.
global.f90        - global parameters for the specific value function being solved.
objective.f90     - the specific objective function being solved. Require the following functions
                    to be defined: objFun, dfovec, initial0 

-----------------------------------
4. Description of text files
-----------------------------------
config.txt        - the configuration file for execution
readme.txt        - this file
logfile.txt       - the log of the program execution

-----------------------------------
5. Description of .dat files
-----------------------------------
initTerm.dat      - the instance which is the main "driver" program
internalConfig.dat- a copy of the config.txt file, but for use by parallel instances
lastParam.dat     - the list of sobol point to be used with minimization
lastSobol.dat     - the list of sobol point at which the objective function is to be evaluated
searchResults.dat - the minimized objective function and associated parameters
searchStart.dat   - the starting point (and sobol point) for the associated line in
                    searchResults.dat (i.e. line 3 here is the starting point for line 3
                    in searchResults.dat)
seq.dat           - the number of the last concurrent instance started
sobol.dat         - the list of sobol points
sobolFnVal.dat    - the value of the objective function at each sobol point
state.dat         - the current state of the program
x_starts.dat      - sorted values of sobolFnVal, used by minimization routine to converge
                    to global minimum                    
stateTime.dat     - the time when the state file was last altered. This is used to check if
                    processes are still alive.                    

-----------------------------------
6. Specifying the objective function
-----------------------------------
The objective function should be specified in the file "otherStuff.f90". Currently two objective
functions need to be specified, one for each minimization routine.

A. amoeba
The amoeba routine takes a function of the following form:

			FUNCTION objFunc(theta)
			    use genericParams
			    use nrtype
                implicit none
			    REAL(DP),DIMENSION(p_nx),INTENT(IN) :: theta
			    REAL(DP) :: objFunc
			END FUNCTION objFunc
			
B. bobyq
The bobyq routine requires a function named dfovec. From the comments of the bobyq code,

            SUBROUTINE dfovec(n, mv, x, v_err)
     
It must provide the values of the vector function v_err(x) : R^n to R^{mv}
at the variables X(1),X(2),...,X(N), which are generated automatically in
a way that satisfies the bounds given in XL and XU.

-----------------------------------
7. Troubleshooting
-----------------------------------
1. XXXXXX file not found
Sometimes this error occurs if the objective function being solved runs too quickly (multiple
times per second). Essentially, the file system is not able to keep up with the state check calls,
resulting in the system being unable to find files that exist. Slowing down the objective function
execution (to once per second, for example) should solve this problem.
Note: This error has been seen on Windows systems, but not on linux systems.