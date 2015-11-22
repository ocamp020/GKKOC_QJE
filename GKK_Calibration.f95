
MODULE GKK_Calibration

Contains
!========================================================================================
!========================================================================================
!========================================================================================

Function Moments_Objective(theta)
	USE parameters
	USE GLOBAL
	use programfunctions
	use Toolbox

	IMPLICIT NONE
	! Parameters
		REAL(DP), INTENT(IN) :: theta(6)
	! Objective Moments
		REAL(DP) :: Obj_Moments(6)
		Real(DP) :: Moments_Objective

	! Resutls Folder
		Result_Folder = './GKK_Files/'
		! call execute_command_line( 'mkdir -p ' // trim(Result_Folder) )
		call system( 'mkdir -p ' // trim(Result_Folder) )

	! Switch for separable and non-separable utility
		! If Utility_Switch==1 then do non-separable utility
		! If Utility_Switch==0 then do separable utility
		Utility_Switch = 1	

	! Set Parameters
		!Params =[ 0.9436, 0.00, 0.50, 0.70444445, 0.34, 0.4494 ] ! tauL=0.224, tauC=0.075 calibration
		beta   = theta(1)
		mu_z   = theta(2) 
		rho_z  = theta(3) 
		sigma_z_eps      = theta(4)
		sigma_lambda_eps = theta(5)
		gamma  = theta(6)

		sigma  = 4.0_dp
		phi    = (1.0_dp-gamma)/gamma

	! Taxes
	! Consumption tax
		tauC=0.075_DP
	! Set Labor Tax Regime
		!tauPL=0.185_DP
		!psi=0.77_DP  
 		tauPL=0.0_DP
 		psi=0.776_DP
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
		CALL COMPUTE_STATS

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