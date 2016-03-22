
MODULE GKK_Calibration

Contains
!========================================================================================
!========================================================================================
!========================================================================================

Function Moments_Objective(par_in)
	USE parameters
	USE GLOBAL
	use programfunctions
	use Toolbox
	use omp_lib

	IMPLICIT NONE
	! Parameters
		REAL(DP), INTENT(IN) :: par_in(5)
	! Objective Moments
		REAL(DP) :: Obj_Moments(6)
		Real(DP) :: Moments_Objective

	! Resutls Folder
		Result_Folder = './'
		! call execute_command_line( 'mkdir -p ' // trim(Result_Folder) )
		! call system( 'mkdir -p ' // trim(Result_Folder) )

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
		!Params =[ 0.9436, 0.00, 0.50, 0.70444445, 0.34, 0.4494 ] ! tauL=0.224, tauC=0.075 calibration
		beta   			 = par_in(1)
		a_x   			 = par_in(2)
		sigma_z_eps      = par_in(3)
		sigma_lambda_eps = par_in(4)
		gamma  			 = par_in(5)

	! Utility Parameters
		sigma  = 4.0_dp
		phi    = (1.0_dp-gamma)/gamma

	! Taxes
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
 	! Capital and wealth taxes
 		tauK = 0.25_DP
		tauW_bt = 0.00_DP
		tauW_at = 0.00_DP
		Y_a_threshold = 0.00_DP   	


	!====================================================================================================

	! Benchmark economy
		solving_bench=1
		
	! Set initia lvalues of R, Wage, Ebar to find equilibrium
		rr    =  4.906133597851297E-002 
		wage  =  1.97429920063330 
		Ebar  =  1.82928004963637  
		Ebar_bench = Ebar

	! Solve for the model and compute stats
		CALL INITIALIZE
		CALL FIND_DBN_EQ
		CALL GOVNT_BUDGET 
		CALL Firm_Value
		CALL COMPUTE_MOMENTS

	!====================================================================================================
	! Print moments (square difference)
		Obj_Moments(1) = (Wealth_Output/3.00_dp 		 -1.0_dp)**2.0_dp 
		Obj_Moments(2) = (prct1_wealth/0.34_dp 			 -1.0_dp)**2.0_dp 
		Obj_Moments(3) = (prct10_wealth/0.71_dp          -1.0_dp)**2.0_dp 
		Obj_Moments(4) = (Std_Log_Earnings_25_60/0.80_dp -1.0_dp)**2.0_dp 
		Obj_Moments(5) = (meanhours_25_60/0.40_dp		 -1.0_dp)**2.0_dp 
		Obj_Moments(6) = (MeanReturn/0.069_DP 			 -1.0_dp)**2.0_dp 

		SSE_Moments = (1.0-Wealth_Output/3.0_DP)**2.0_DP  + (1.0_DP-prct1_wealth/0.34_DP)**2.0_DP  + (prct10_wealth-0.71_DP)**2.0_DP &
                   & + (1.0_DP-Std_Log_Earnings_25_60 / 0.8_DP)**2.0_DP + (1.0_DP-meanhours_25_60/0.4_DP)**2.0_DP &
                   & + (1.0_DP-MeanReturn/0.069_DP)**2.0_DP

		Moments_Objective = SSE_Moments 

END Function Moments_Objective


!========================================================================================
!========================================================================================
!========================================================================================

END Module GKK_Calibration