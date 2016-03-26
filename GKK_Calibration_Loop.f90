

! This program solves for the wealth taxation model of 
! Guvenen, Kuruscu, Kambourov and Chen

! Three procedures are included
	! Solve for benchmark economy with capital income taxes
	! Solve for a tax reform to wealth taxes (Experiment)
	! Solve for optimal taxes:
		! Optimal capital income taxes and linear labor taxes to balance budget
		! Optimal wealth taxes and linear labor taxes to balance budget	

! The code allows for:
	! Progressive labor income taxation
	! Threshold on wealth for wealth taxation
	! Non-Separable utility

! When utility is separable it is:
	! U(c,h) = log(c) + phi*log(1-h)

! When utility is non-separable it is:
	! U(c,h) = (c^(gamma)(1-l)^(1-gamma))^(1-sigma) / (1-sigma)


!========================================================================================
!========================================================================================
!========================================================================================

PROGRAM GKK_Calibration_Loop
	USE parameters
	USE GLOBAL
	use programfunctions
	use Toolbox
	use omp_lib

	IMPLICIT NONE
	! Variables to measure running time
		REAL(DP) :: start_time, finish_time, Loop_Par(5)
	! Compute benchmark or load results
		logical  :: compute_bench, compute_exp, Opt_Tax, Opt_Tax_KW, Tax_Reform, Simul_Switch, Calibration_Switch
	! Auxiliary variable for writing file
		character(4)   :: string_theta
		character(100) :: folder_aux

	! Resutls Folder
		Result_Folder = './Calibration_Loop/'


	! Capital Market
		theta = 1.50_dp
	! Threshold 
		Threshold_Factor = 0.00_dp 

	! Switch for solving benchmark or just reading resutls
		! If compute_bench==.true. then just read resutls
		! If compute_bench==.false. then solve for benchmark and store results
		compute_bench = .true.

	! Switch for separable and non-separable utility
		! If NSU_Switch==.true. then do non-separable utility
		! If NSU_Switch==.false. then do separable utility
		NSU_Switch = .true.

	! Switch for log utility 
		! If Log_Switch==.true. then utility is log
		! If Log_Switch==.false. then utility is not log
		Log_Switch = .false.

	! Switch for labor taxes
		! If Progressive_Tax_Switch==.true. then use progressive taxes
		! If Progressive_Tax_Switch==.false. then use linear taxes
		Progressive_Tax_Switch = .false.

	! Set Parameters
		Params =[0.962_dp, 0.0_dp, 0.50_dp, 0.387_dp, 0.29_dp, 0.4494_dp] ! alpha=0.4, zgrid 11, m5, alpha=0.4, dep005, mu=090, K/Y=3, Top1PVa=0.36
		

		beta   	= params(1)
		mu_z   	= params(2) ! this is just shifting the z grids. it is zero now.
		rho_z  	= params(3) 
		sigma_z_eps      = params(4)
		sigma_lambda_eps = params(5)
		gamma  	= params(6)
		
		sigma  	= 4.0_dp
		phi    	= (1.0_dp-gamma)/gamma

		x_hi	= 3.00_dp
		x_lo	= 1.0_dp
		a_x 	= 0.10_dp
		b_x 	= 0.00_dp

	! Set Parameters from Loop
		OPEN (UNIT=3, FILE=trim(Result_Folder)//'Loop_Par'  , STATUS='old', ACTION='read')
		READ (UNIT=3,  FMT=*), Loop_Par
		CLOSE(unit=3)

		beta 		= Loop_Par(1)
		sigma_z_eps = Loop_Par(2)
		x_hi 		= Loop_Par(3)
		a_x  		= Loop_Par(4)
		b_x  		= Loop_Par(5)


		if (Log_Switch.eqv..true.) then
				sigma = 1.0_dp
			if (NSU_Switch.eqv..false.) then 
				gamma = 1.0_dp 
			endif 
		endif

	! Taxes
	! Wealth tax: minimum wealth tax to consider and increments for balancing budget
		tauWmin_bt=0.00_DP
		tauWinc_bt=0.000_DP ! Minimum tax below threshold and increments
		tauWmin_at=0.010_DP
		tauWinc_at=0.005_DP ! Minimum tax above threshold and increments
	! Consumption tax
		tauC=0.075_DP
	! Set Labor Tax Regime
		if (Progressive_Tax_Switch.eqv..true.) then 
			tauPL = 0.185_DP
			psi   = 0.776_DP  
		else 
			tauPL = 0.0_DP
	 		psi   = 0.776_DP  	
	 	endif 
	
		print*, "Results are stored in directory: ", Result_Folder
		print*,'na=',na,'update_period=',update_period
		print*, "NSU_Switch=",NSU_Switch,'sigma=',sigma,'gamma=',gamma,'phi',phi
		print*,'Labor Taxes: tauPl=',tauPl,'psi',psi
		print*, 'Borrowing Constraint: Theta=',theta
		print*, 'm tauchen for zgrid is ',mtauchen_z,'nz=',nz, 'amax=',amax, 'totpop=', totpop
		print*, 'x_hi', x_hi, 'a_x', a_x, 'b_x', b_x

	! Start timining of  the process
		call cpu_time(start_time) 

	! Set initia lvalues of R, Wage, Ebar to find equilibrium
		! ------- DO NOT REMOVE THE LINES BELOW
		R     =  0.05_dp
		P     =  4.906133597851297E-002 
		wage  =  1.97429920063330 
		Ebar  =  1.82928004963637  
		Ebar_bench = Ebar
		 

		! ------- DO NOT REMOVE THE LINES ABOVE


	! Call routines
	
	!====================================================================================================
	PRINT*,''
	Print*,'--------------- SOLVING BENCHMARK WITH BEST PARAMETERS -----------------'
	PRINT*,''
	print*,'CAPITAL TAX ECONOMY'

	! Benchmark economy
		solving_bench=1

	! Set taxes for benchmark economy
		tauK = 0.25_DP
		tauW_bt = 0.00_DP
		tauW_at = 0.00_DP
		Y_a_threshold = 0.00_DP 

	! Solve for the model and compute stats
		CALL INITIALIZE
		CALL FIND_DBN_EQ
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction)
		CALL Firm_Value
		CALL COMPUTE_STATS

	! Save Results 
	OPEN(UNIT=3, FILE=trim(Result_Folder)//'Calibration_Loop_Restuls.txt', STATUS='old', POSITION='append') 
	WRITE(unit=3, FMT=*) ' '
	WRITE(unit=3, FMT=*) beta,sigma_z_eps,x_hi,a_x,b_x,Wealth_Output,Std_Log_Earnings_25_60,meanhours_25_60,MeanReturn, &
		& FW_top_x_share(4), FW_top_x_share(3)
	CLOSE(unit=3)
		
end Program GKK_Calibration_Loop