

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

PROGRAM main
	USE parameters
	USE GLOBAL
	use programfunctions
	use GKK_Stats
	use Simulation_Module
	use Toolbox
	use omp_lib

	IMPLICIT NONE
	! Variables to measure running time
		REAL(DP) :: start_time, finish_time
	! Compute benchmark or load results
		logical  :: compute_bench, compute_exp, Opt_Tax, Opt_Tax_KW, Tax_Reform, Simul_Switch, Calibration_Switch
		logical  :: Opt_Threshold, Opt_Tau_C, Opt_Tau_CX, Opt_Tax_K_and_W, Tax_Reform_KW
		logical  :: compute_exp_pf, Fixed_PF, Fixed_PF_interp, Fixed_PF_prices, compute_exp_fixed_prices_and_taxes
		logical  :: compute_exp_prices, Fixed_W, Fixed_P, Fixed_R , Tax_Reform_Decomposition
		logical  :: Transition_Tax_Reform, Transition_OT, budget_balance, balance_tau_L
		logical  :: Tax_Reform_tau_C, compute_exp_tau_c, Opt_Tax_KW_TR, Tax_Reform_tktw, budget_flag
	! Auxiliary variable for writing file
		character(4)   :: string_theta
		character(100) :: folder_aux

	! Allocate Variables
	call Allocate_Variables

	! Threshold 
		Threshold_Factor = 0.00_dp 

	! Switch for solving benchmark or just reading resutls
		Calibration_Switch = .false.
		! If compute_bench==.true. then just read resutls
		! If compute_bench==.false. then solve for benchmark and store results
		Tax_Reform    = .true.
			compute_bench = .false.
			compute_exp   = .false.
			compute_exp_pf= .false.
				Fixed_PF        = .true.
				Fixed_PF_interp = .true.
				Fixed_PF_prices = .false.
			compute_exp_prices    = .false.
				Fixed_W = .true. 
				Fixed_P = .true.
				Fixed_R = .true.
				
		Tax_Reform_tktw = .false.
			budget_flag = .false.

		Opt_Tax       = .false.
			Opt_Tax_KW    = .false. ! true=tau_K, false=tau_W

		Opt_Threshold = .false.

		
		Transition_Tax_Reform = .false.
		Transition_OT = .false.
			budget_balance = .true.
			balance_tau_L  = .true. ! true=tau_L, false=tau_K or tau_W depending on Opt_Tax_KW
			Opt_Tax_KW_TR  = .true. ! true=tau_K, false=tau_W
		
		Simul_Switch  = .false.



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
		! Calibration to book value and rho_z=0.1
		Params =[ 0.9475_dp, 0.00_dp, 0.1_dp, 0.072_dp , 0.305_dp, 0.46_dp ] 
		! Calibration to book value and rho_z=0.2
		Params =[ 0.9472_dp, 0.00_dp, 0.2_dp, 0.070_dp , 0.307_dp, 0.46_dp ] 
		! Calibration to book value and rho_z=0.3
		Params =[ 0.947_dp, 0.00_dp, 0.3_dp, 0.0669_dp , 0.307_dp, 0.46_dp ] 
		! Calibration to book value and rho_z=0.1 and x_hi=10
		Params =[ 0.9473_dp, 0.00_dp, 0.1_dp, 0.0352_dp , 0.307_dp, 0.46_dp ] 
		


		! Corporate Sector
			A_C    = 0.0_dp ! 
			! A_C    = 0.9590_dp ! for Corp model ! 0.9409_dp (value without estate tax)

		if (A_C.eq.0.0_dp) then
		
		! Debt/Output = 1.5, lambda = 1.5, no bequest fee
			! Main Parameters 
				beta   	= 0.9593_dp ! 0.9404_dp (Value without estate tax)! 0.9475_dp (value in old benchmark) ! params(1) !
				sigma_z_eps      = 0.277_dp ! 0.0867_dp (Value without estate tax) ! 0.072_dp (value in old benchmark) ! params(4) !
				sigma_lambda_eps = 0.309_dp ! 0.309_dp (Value without estate tax) ! 0.305_dp (value in old benchmark) ! params(5)
				gamma  	= 0.4450_dp ! 0.4580_dp (Value without estate tax) ! 0.46_dp (value in old benchmark) !  params(6) ! 
				sigma  	= 4.0_dp
				x_hi	= 1.50_dp

			! Bequeset parameters chi_bq*(bq+bq_0)^(1-sigma)
				bq_0   = 00.30_dp ! Level shift 00.30_dp (value without estate tax)
				chi_u  = 00.20_dp ! Scaling 03.55_dp (value without estate tax)
				chi_bq = chi_u*(1.0_dp-tau_bq) ! Auxiliary parameter for FOC and EGM

			! Bequest parameters for Model_2.1_no_bq
				! bq_0   = 0.0_dp 
				! chi_u  = 0.0_dp 
				! chi_bq = 0.0_dp 

			! Capital Market
				do zi=1,nz
				theta(zi)    = 1.00_dp+(2.80_dp-1.00_dp)/(nz-1)*(real(zi,8)-1.0_dp)
				enddo

			! No bequest fees
				bq_fee = 0.00_dp

		else

		! Corporate SEctor

		! Main Parameters 
			beta   	= 0.9581_dp ! 0.9404_dp (Value without estate tax)! 0.9475_dp (value in old benchmark) ! params(1) !
			sigma_z_eps      = 0.09333_dp ! 0.0867_dp (Value without estate tax) ! 0.072_dp (value in old benchmark) ! params(4) !
			sigma_lambda_eps = 0.314_dp ! 0.309_dp (Value without estate tax) ! 0.305_dp (value in old benchmark) ! params(5)
			gamma  	=  0.4400_dp ! 0.4580_dp (Value without estate tax) ! 0.46_dp (value in old benchmark) !  params(6) ! 
			sigma  	= 4.0_dp
			x_hi	= 5.00_dp
		
		! Bequeset parameters chi_bq*(bq+bq_0)^(1-sigma)
			bq_0   = 00.30_dp ! Level shift 00.30_dp (value without estate tax)
			chi_u  = 00.10_dp ! Scaling 03.55_dp (value without estate tax)
			chi_bq = chi_u*(1.0_dp-tau_bq) ! Auxiliary parameter for FOC and EGM


		! Capital Market
			do zi=1,nz
			theta(zi)    = 1.00_dp+(2.50_dp-1.00_dp)/(nz-1)*(real(zi,8)-1.0_dp)
			enddo

		endif 

		! Other parameters 
		rho_z  	= 0.1_dp ! params(3)
		phi    	= (1.0_dp-gamma)/gamma
		mu_z   	= params(2) ! this is just shifting the z grids. it is zero now.
		Params =[beta, mu_z, rho_z, sigma_z_eps, sigma_lambda_eps, gamma] 


		x_lo	= 1.00_dp
		x_0     = 0.00_dp
		a_x 	= 0.10_dp
		b_x 	= 0.00_dp

		if (Log_Switch.eqv..true.) then
				sigma = 1.0_dp
			if (NSU_Switch.eqv..false.) then 
				gamma = 1.0_dp 
			endif 
		endif

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

 	! Debt
 		! Start government debt at zero
 		Debt_SS = 0.0_DP

	! Resutls Folder
	if (A_C.eq.0.0_dp) then 
 		if ((Progressive_Tax_Switch.eqv..false.).and.(NSU_Switch.eqv..true.)) then 
			Result_Folder = './Revision/Model_2.1_Local/' 
			! Result_Folder = './Revision/Model_2.1_no_bq/'
		else if ((Progressive_Tax_Switch.eqv..true.).and.(NSU_Switch.eqv..true.)) then 
			Result_Folder = './Revision/Model_2.1_PT/' 
		else if ((Progressive_Tax_Switch.eqv..false.).and.(NSU_Switch.eqv..false.)) then 
			Result_Folder = './Revision/Model_2.1_SU/' 
		else if ((Progressive_Tax_Switch.eqv..true.).and.(NSU_Switch.eqv..false.)) then 
			Result_Folder = './Revision/Model_2.1_PT_SU/' 
		end if
	else 
 		if ((Progressive_Tax_Switch.eqv..false.).and.(NSU_Switch.eqv..true.)) then 
			Result_Folder = './Revision/Model_2.1_Corp/' 
		else if ((Progressive_Tax_Switch.eqv..true.).and.(NSU_Switch.eqv..true.)) then 
			Result_Folder = './Revision/Model_2.1_Corp_PT/' 
		else if ((Progressive_Tax_Switch.eqv..false.).and.(NSU_Switch.eqv..false.)) then 
			Result_Folder = './Revision/Model_2.1_Corp_SU/' 
		else if ((Progressive_Tax_Switch.eqv..true.).and.(NSU_Switch.eqv..false.)) then 
			Result_Folder = './Revision/Model_2.1_Corp_PT_SU/' 
		end if
	endif 

		

		! call execute_command_line( 'mkdir -p ' // trim(Result_Folder) )
		call system( 'mkdir -p ' // trim(Result_Folder) )
		print*,' '; print*,'----------------------------------------------------------------------------------'
		print*, "Results are stored in directory: ", Result_Folder
		print*,'	na=',na,'update_period=',update_period,"NSU_Switch=",NSU_Switch
		print '(A,F7.4,X,A,F7.4,X,A,F7.4)','	Utility parameters:  sigma=',sigma,'gamma=',gamma,'phi',phi
		print '(A,F7.4,X,A,F7.4)','	Labor Taxes: tauPl=',tauPl,'psi',psi
		print '(A,F7.3,X,F7.3,X,F7.3,X,F7.3,X,F7.3,X,F7.3,X,F7.3,X,F7.3,X,F7.3)', '	Borrowing Constraint: Theta=',theta
		print '(A,F7.3,X,A,I3,X,A,F11.2,X,A,I10)', '	m tauchen for zgrid is ',mtauchen_z,'nz=',nz, 'amax=',amax, 'totpop=', totpop
		print '(A,F7.4,X,A,F7.4,X,A,F7.4,X,A,F7.4)', '	x_hi', x_hi, 'a_x', a_x, 'b_x', b_x, 'sigmaz', sigma_z_eps
		print*,'----------------------------------------------------------------------------------'; print*,' '


	! Set parameters to be used in all simulations economy
		OPEN(UNIT=3, FILE=trim(Result_Folder)//'params', STATUS='replace')
		WRITE(unit=3, FMT=*) params
		CLOSE(unit=3)


	! Start timining of  the process
		call cpu_time(start_time) 

	! Set initia lvalues of R, Wage, Ebar to find equilibrium
		! ------- DO NOT REMOVE THE LINES BELOW
		R     =  0.05_dp
		P     =  4.906133597851297E-002_dp
		wage  =  1.97429920063330 
		Ebar  =  1.82928004963637  
		Ebar_bench = Ebar
		 

		! ------- DO NOT REMOVE THE LINES ABOVE


	! Call routines
		! Calibration
		! if (Calibration_Switch) then 
		! 	CALL CALIBRATION_TRIALS
		! endif 

		! Tax Reform experiment
		if (Tax_Reform) then 

			call Solve_Benchmark(compute_bench,Simul_Switch)
			Bench_Simul_Folder = trim(Result_Folder)//'Simul/'
			
			if ((compute_exp_pf).and.(compute_exp_prices.eqv..false.)) then 
				if     ((Fixed_PF).and.(Fixed_PF_interp.eqv..false.).and.(Fixed_PF_prices.eqv..false.)) then
					Result_Folder = trim(Result_Folder)//'Exp_Policy_Functions/'
					call system( 'mkdir -p ' // trim(Result_Folder) )
					call Solve_Experiment_Fixed_Policy_Functions(Fixed_PF,Simul_Switch)
				elseif ((Fixed_PF.eqv..false.).and.(Fixed_PF_interp).and.(Fixed_PF_prices.eqv..false.)) then
					! If using benchmark prices
						Result_Folder = trim(Result_Folder)//'Exp_Policy_Functions_Interp/'
						call system( 'mkdir -p ' // trim(Result_Folder) )
					! If using tax reform prices
						! call Solve_Experiment(.false.,.false.)
						! 	Ebar   = EBAR_bench
						! 	P      = P_exp
						! 	R	   = R_exp 
						! 	wage   = wage_exp
						! 	DBN1   = DBN_bench
						! 	Cons   = Cons_bench        
						! 	Hours  = Hours_bench
						! 	Aprime = Aprime_bench
						! Result_Folder = trim(Result_Folder)//'Exp_Policy_Functions_Interp_Prices_Exp/'
						! call system( 'mkdir -p ' // trim(Result_Folder) )
					call Solve_Experiment_Fixed_PF_Interp(Fixed_PF_interp,Simul_Switch)
				elseif ((Fixed_PF.eqv..false.).and.(Fixed_PF_interp.eqv..false.).and.(Fixed_PF_prices)) then
					Result_Folder = trim(Result_Folder)//'Exp_Policy_Functions_Prices/'
					call system( 'mkdir -p ' // trim(Result_Folder) )
					call Solve_Experiment_Fixed_PF_Prices(Fixed_PF_prices,Simul_Switch)
				endif 
			elseif ((compute_exp_pf.eqv..false.).and.(compute_exp_prices)) then 
					if (Fixed_W.and.(Fixed_P.eqv..false.).and.(Fixed_R.eqv..false.)) Result_Folder = trim(Result_Folder)//'Exp_Fixed_Prices_W/'
					if ((Fixed_W.eqv..false.).and.Fixed_P.and.(Fixed_R.eqv..false.)) Result_Folder = trim(Result_Folder)//'Exp_Fixed_Prices_P/'
					if ((Fixed_W.eqv..false.).and.(Fixed_P.eqv..false.).and.Fixed_R) Result_Folder = trim(Result_Folder)//'Exp_Fixed_Prices_R/'
					if (Fixed_W.and.Fixed_P.and.(Fixed_R.eqv..false.)) Result_Folder = trim(Result_Folder)//'Exp_Fixed_Prices_WP/'
					if (Fixed_W.and.(Fixed_P.eqv..false.).and.Fixed_R) Result_Folder = trim(Result_Folder)//'Exp_Fixed_Prices_WR/'
					if ((Fixed_W.eqv..false.).and.Fixed_P.and.Fixed_R) Result_Folder = trim(Result_Folder)//'Exp_Fixed_Prices_PR/'
					if (Fixed_W.and.Fixed_P.and.Fixed_R) Result_Folder = trim(Result_Folder)//'Exp_Fixed_Prices_WPR/'
					call system( 'mkdir -p ' // trim(Result_Folder) )
					call Solve_Experiment_Fixed_Prices(compute_exp_prices,Simul_Switch,Fixed_W,Fixed_P,Fixed_R)
			else
				if (KeepSSatBench .eq. 1) then 
				Result_Folder = trim(Result_Folder)//'Tax_Reform/'
				else 
				folder_aux = Result_Folder
				Result_Folder = trim(folder_aux)//'Tax_Reform/'
				CALL Write_Experimental_Results(.false.)
				Result_Folder = trim(folder_aux)//'Tax_Reform_SS/'
				endif
				call system( 'mkdir -p ' // trim(Result_Folder) )
				call Solve_Experiment(compute_exp,Simul_Switch)
			endif 

			compute_bench = .false.
		endif 

		if (Tax_Reform_tktw) then 
			call Solve_Experiment_tktw(budget_flag)
		endif 
	
		! Optimal Tax
		if (Opt_Tax) then 
			call Solve_Benchmark(compute_bench,Simul_Switch)

			folder_aux = Result_Folder
			if (Opt_Tax_KW) then 
				if (KeepSSatBench .eq. 1) then 
				Result_Folder = trim(folder_aux)//'Opt_Tax_K/'
				else 
				Result_Folder = trim(folder_aux)//'Opt_Tax_K_SS/'
				endif 
			else 
				if (KeepSSatBench .eq. 1) then 
				Result_Folder = trim(folder_aux)//'Opt_Tax_W/'
				else 
				Result_Folder = trim(folder_aux)//'Opt_Tax_W_SS/'
				endif 
			endif
			call system( 'mkdir -p ' // trim(Result_Folder) )

			
			call Solve_Opt_Tax(Opt_Tax_KW,Simul_Switch)


		endif 

		if (Opt_Threshold) then 
			

			! Solve for optmal threshold 
				call Solve_Benchmark(compute_bench,Simul_Switch)
				call Solve_Opt_Threshold

			! Optimal taxes when threshold = 2 
				! print*, "Solving for optimal wealth tax with threshold factor of 2"
				! ! Solve Benchmark 
				! call Solve_Benchmark(compute_bench,Simul_Switch)
				! ! Set up folder and flags 
				! Opt_Tax_KW    = .false. 
				! Simul_Switch  = .false. 
				! folder_aux    = Result_Folder
				! Result_Folder = trim(folder_aux)//'Opt_Tax_W_Threshold_2/'
				! call system( 'mkdir -p ' // trim(Result_Folder) )
				! ! Set up threshold 
				! Threshold_Factor = 2.0_dp 
				! print*, ' Threshold_Factor=',Threshold_Factor
				! ! Solve for optimal taxes 
				! call Solve_Opt_Tax(Opt_Tax_KW,Simul_Switch)

			! Optimal taxes when threshold = 1
				! print*, "Solving for optimal wealth tax with threshold factor of 1"
				! ! Solve Benchmark 
				! call Solve_Benchmark(compute_bench,Simul_Switch)
				! ! Set up folder and flags 
				! Opt_Tax_KW    = .false. 
				! Simul_Switch  = .false. 
				! folder_aux    = Result_Folder
				! Result_Folder = trim(folder_aux)//'Opt_Tax_W_Threshold_1/'
				! call system( 'mkdir -p ' // trim(Result_Folder) )
				! ! Set up threshold 
				! Threshold_Factor = 1.0_dp 
				! print*, ' Threshold_Factor=',Threshold_Factor
				! ! Solve for optimal taxes 
				! call Solve_Opt_Tax(Opt_Tax_KW,Simul_Switch)

		endif 

		
		if (Transition_Tax_Reform) then
			call Solve_Transition_Tax_Reform(budget_balance)
		endif

		if (Transition_OT) then
			call Solve_Transition_Opt_Taxes(Opt_Tax_KW_TR,budget_balance,balance_tau_L)
		endif


	call cpu_time(finish_time)
	print*,'Total time =',finish_time-start_time


END PROGRAM main


!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Benchmark(compute_bench,Simul_Switch)
	use parameters
	use global 
	use programfunctions
	use GKK_Stats
	use Simulation_Module
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: compute_bench, Simul_Switch

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
		psi_bench   = psi

	! Solve for the model and compute stats
	print*,"	Initializing program"
		CALL INITIALIZE
		
	if (compute_bench) then
		print*,"	Reading initial conditions from file"
		! CALL Write_Benchmark_Results(.false.)
		print*,"	Computing equilibrium distribution"
		CALL FIND_DBN_EQ
		print*,"	Computing government spending"
		CALL GOVNT_BUDGET(.true.)
		print*,"	Computing Value Function"
		! CALL COMPUTE_VALUE_FUNCTION_SPLINE
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
		print*,"	Computing Firm Value Function"
		CALL Firm_Value
		print*, " 	Computing After Tax Income"
		CALL Compute_After_Tax_Income
		print*,"	Saving results in text files to be read later"
		CALL Write_Benchmark_Results(compute_bench)
	else
		print*,"	Reading benchmark results from files"
		CALL Write_Benchmark_Results(compute_bench)
		CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
		K_mat  = K_Matrix(R,P)
		Pr_mat = Profit_Matrix(R,P)
		CALL ComputeLaborUnits(EBAR,wage)
		CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
			YGRID_aux = YGRID; MBGRID_aux = MBGRID ; RetY_lambda_e_aux = RetY_lambda_e;
			! tauK = 0.0_dp 
			! call Find_TauW_Threshold(DBN1,W_bench)  
			! print*,' ' 
			! print*,'W_Bench=',W_bench
			! STOP
		CALL GOVNT_BUDGET(.true.)
	end if 

	! Aggregate variables in benchmark economy
		GBAR_bench  = GBAR
		QBAR_bench  = QBAR 
		NBAR_bench  = NBAR 
		Ebar_bench  = EBAR
		P_bench     = P
		R_bench     = R
		wage_bench  = wage
		Y_bench     = YBAR
		tauK_bench  = tauK
		tauPL_bench = tauPL
		psi_bench   = psi
		DBN_bench   = DBN1
		tauw_bt_bench = tauW_bt
		tauw_at_bench = tauW_at
		Y_a_threshold_bench = Y_a_threshold

		DBN_bench           = DBN1
		ValueFunction_bench = ValueFunction
		Bq_Value_bench 		= Bq_Value
		Cons_bench          = Cons           
		Hours_bench         = Hours
		Aprime_bench        = Aprime 
		V_Pr_bench          = V_Pr
		V_Pr_nb_bench       = V_Pr_nb 

		SSC_Payments_bench  = SSC_Payments

		YBAR_C_bench = YBAR_C
		L_C_bench 	 = L_C
		K_C_bench    = K_C
		YBAR_P_bench = YBAR_P
		L_P_bench 	 = L_P
		K_P_bench    = K_P

		print*,"	Computing satitics"
		CALL COMPUTE_STATS
		print*,"	Writing variables"
		CALL WRITE_VARIABLES(1)
		if (Simul_Switch) then 
			print*,"	Simulation"
			CALL SIMULATION(solving_bench)
			! CALL Simulation_Life_Cycle_Patterns(solving_bench)
			CALL Simulation_Life_Cycle_Asset_Return_Panel(solving_bench)
		endif

		print*,' '
		print*,'-------------------------------------------------------------------------'
		write(*,*) " Benchmark variables"
		print 12345," 	GBAR=",GBAR,"EBAR=",EBAR,"NBAR=",NBAR,"QBAR=",QBAR,"P=",P,"wage=",wage,'R(%)=',100.0_dp*R
		print*,' '
		print*,' Targets:'
		print 12345, &
			& " 	A/Y=",MeanWealth/YBAR,'BQ/A=',100.0_dp*Bequest_Wealth/MeanWealth ,'BQ/Inc=',Bq_Inc(3,1),&
			& 'Top_1%=',100.0_dp*prct1_wealth,'L_C/N=',100.0_dp*L_C/NBAR,&
			& 'stdEarn=',Std_Log_Earnings_25_60,'N',meanhours_25_60,'D/Y',External_Debt_GDP,&
			& 'Luxury',(EBAR_data/(EBAR_bench*0.727853584919652_dp))*bq_0
		12345 format (A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F12.1)
		print*,'-------------------------------------------------------------------------'
		print*,' ';print*,' ';print*,' '

		! Deallocate variables
		if (compute_bench.eqv..false.) then 
		deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )
		endif 


		! print*,"	Efficiency Computation"
		! CALL Hsieh_Klenow_Efficiency(solving_bench)

		

end Subroutine Solve_Benchmark


!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Experiment(compute_exp,Simul_Switch)
	use parameters
	use global 
	use programfunctions
	use GKK_Stats
	use Simulation_Module
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: compute_exp, Simul_Switch

	!====================================================================================================
	PRINT*,''
	Print*,'--------------- SOLVING EXPERIMENT WITH BEST PARAMETERS -----------------'
	PRINT*,''
	print*,'Wealth Tax Economy'
	
	! Experiment economy
		solving_bench=0
	! Set capital taxes to zero
		tauK = 0.0_DP
	! Set Y_a_threshold
		write(*,*) "Y_a threshold is set to a proportion of the mean wealth under current distribution"
		!Y_a_threshold = 0.0_dp ! big_p   !8.1812138704441200
		call Find_TauW_Threshold(DBN_bench,W_bench)  
		Y_a_threshold = Threshold_Factor*Ebar_bench !0.75_dp
		Wealth_factor = Y_a_threshold/W_bench

		! Wealth tax: minimum wealth tax to consider and increments for balancing budget
		tauWmin_bt=0.00_DP
		tauWinc_bt=0.000_DP ! Minimum tax below threshold and increments
		tauWmin_at=0.010_DP
		tauWinc_at=0.0005_DP ! Minimum tax above threshold and increments
		if (KeepSSatBench .eq. 0) then
		tauWmin_at = tauW_at
		tauWinc_at = 0.0005_dp
		endif 

	if (compute_exp) then 
		! Find wealth taxes that balances budget
		print*, "	Computing Wealth Tax to balance the budget"
			! Set initial value for G in experimental economy and for wealth taxes
			GBAR_exp = 0.0_DP
			tauW_bt  = tauWmin_bt
			tauW_at  = tauWmin_at
			tauWindx = 0.0_DP
			! Solve for the model increasing wealth taxes until revenue is enough to finance G_benchamark
			DO WHILE (GBAR_exp .lt. GBAR_bench)
				! Set old G and new value of tauW
				GBAR_exp_old = GBAR_exp
				tauW_bt = tauWmin_bt + tauWindx * tauWinc_bt
				tauW_at = tauWmin_at + tauWindx * tauWinc_at
				
				! Solve the model
				CALL FIND_DBN_EQ
				! Solve Government Budget
				CALL GOVNT_BUDGET(.true.)
				! Get new G
				GBAR_exp = GBAR 
				! Iteratioins  
				tauWindx = tauWindx + 1.0_DP   
				print '(A,F7.3,X,X,A,F7.3)', "Bracketing GBAR: tauW_bt=", tauW_bt*100.0_dp, "And tauW_at=", tauW_at*100.0_dp
				print '(A,F7.3,X,X,A,F7.3)', "Current tauW Threshold=", Y_a_threshold, "Share above threshold=", Threshold_Share
				print '(A,F10.6,X,X,A,F10.6)','GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
			ENDDO

			! Set tauW as weighted average of point in  the grid to balance budget more precisely
				tauW_up_bt  = tauW_bt
				tauW_low_bt = tauW_bt  -  tauWinc_bt
				tauW_up_at  = tauW_at
				tauW_low_at = tauW_at  -  tauWinc_at
				if (KeepSSatBench .eq. 1) then 
				tauW_bt     = tauW_low_bt + tauWinc_bt * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)  
				tauW_at     = tauW_low_at + tauWinc_at * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
				else
				tauW_bt     = (tauW_low_bt + tauW_up_bt)*0.5_dp
				tauW_at     = (tauW_low_at + tauW_up_at)*0.5_dp
				endif
				print*,''
				print*,'GBAR bracketed by taxes:'
				print '(A,F7.3,A,F7.3,A,F7.3,A)',&
					&'tauW_low_bt =', tauW_low_bt*100.0_dp, '%  tauW_up_bt=', tauW_up_bt*100.0_dp, '%  tauW_bt=', tauW_bt*100.0_dp, "%"
				print '(A,F7.3,A,F7.3,A,F7.3,A)',&
					&'tauW_low_at =', tauW_low_at*100.0_dp, '%  tauW_up_at=', tauW_up_at*100.0_dp, '%  tauW_at=', tauW_at*100.0_dp, "%"
				print*,''

			! Solve (again) experimental economy
				CALL FIND_DBN_EQ
				CALL GOVNT_BUDGET(.true.)

			! Find tauW that exactly balances the budget (up to precisioin 0.1) using bisection
				GBAR_exp = GBAR
				print*,"Gbar at midpoint of bracket and GBAR at benchmark"
				print '(A,F10.6,X,X,A,F10.6)','GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
				print*,''
				print*,'Bisection for TauW:'
				DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.001 ) ! as long as the difference is greater than 0.1% continue
				    if (GBAR_exp .gt. GBAR_bench ) then
				        tauW_up_bt  = tauW_bt 
				        tauW_up_at  = tauW_at 
				    else
				        tauW_low_bt = tauW_bt
				        tauW_low_at = tauW_at
				    endif
				    tauW_bt = (tauW_low_bt + tauW_up_bt)/2.0_DP
				    tauW_at = (tauW_low_at + tauW_up_at)/2.0_DP
				    CALL FIND_DBN_EQ
				    CALL GOVNT_BUDGET(.true.)
				    GBAR_exp = GBAR
				    print '(A,F7.3,A,F7.3,A,F7.3,A)',&
					&'tauW_low_bt =', tauW_low_bt*100.0_dp, '%  tauW_up_bt=', tauW_up_bt*100.0_dp, '%  tauW_bt=', tauW_bt*100.0_dp, "%"
				print '(A,F7.3,A,F7.3,A,F7.3,A)',&
					&'tauW_low_at =', tauW_low_at*100.0_dp, '%  tauW_up_at=', tauW_up_at*100.0_dp, '%  tauW_at=', tauW_at*100.0_dp, "%"
					print '(A,F7.3,X,X,A,F7.3)', "Current tauW Threshold=", Y_a_threshold, "Share above threshold=", Threshold_Share
					print '(A,F10.6,X,X,A,F10.6)','GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
					print*, ' '
				ENDDO

		! Compute value function and store policy functions, value function and distribution in file
		! CALL COMPUTE_VALUE_FUNCTION_SPLINE 
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
		CALL Firm_Value
		print*, " 	Computing After Tax Income"
		CALL Compute_After_Tax_Income

	endif 
	
	CALL Write_Experimental_Results(compute_exp)
	CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
	K_mat  = K_Matrix(R,P)
	Pr_mat = Profit_Matrix(R,P)
	CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
	CALL ComputeLaborUnits(EBAR,wage)
	CALL GOVNT_BUDGET(.true.)

	! Aggregate variable in experimental economy
		GBAR_exp  = GBAR
		QBAR_exp  = QBAR 
		NBAR_exp  = NBAR  
		Y_exp 	  = YBAR
		Ebar_exp  = EBAR
		P_exp     = P
		R_exp	  = R
		wage_exp  = wage
		tauK_exp  = tauK
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

		ValueFunction_exp = ValueFunction
		Bq_Value_exp      = Bq_Value
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime
		V_Pr_exp          = V_Pr 
		V_Pr_nb_exp  	  = V_Pr_nb

		YBAR_C_exp = YBAR_C
		L_C_exp    = L_C
		K_C_exp    = K_C
		YBAR_P_exp = YBAR_P
		L_P_exp    = L_P
		K_P_exp    = K_P

	! Compute moments
	CALL COMPUTE_STATS
	
	! Compute welfare gain between economies
	CALL COMPUTE_WELFARE_GAIN
	CALL COMPUTE_WELFARE_DECOMPOSITION

	! Write experimental results in output.txt
	CALL WRITE_VARIABLES(0)

	! Simulation
	! if ((Simul_Switch)) then 
	!  	print*,"	Experiment Simulation"
	! 	CALL SIMULATION(solving_bench)
	! endif
	! Call Simulation_Life_Cycle_Patterns(solving_bench)
	Call Simulation_Life_Cycle_Asset_Return_Panel(solving_bench)

	print*,'---------------------------'
	print*,''
	print*,'Output Gain Prct=', 100.0_DP*(Y_exp/Y_bench-1.0) 
	print*,''
	print*,'---------------------------'

	print*," "
	print*,"Wealth_factor=",Wealth_factor
	print*," "


	! Deallocate variables
		if (allocated(YGRID_t)) then 
			deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )
		endif

	! print*,"	Efficiency Computation"
	! 	CALL Hsieh_Klenow_Efficiency(solving_bench)
	print*, ' End of Solve_Experiment'

end Subroutine Solve_Experiment

!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Experiment_tktw(budget_flag)
	use parameters
	use global 
	use programfunctions
	use GKK_Stats
	use Simulation_Module
	use Toolbox
	use Opt_Tax_Parameters
	use Opt_Tax_Functions
	use omp_lib
	implicit none 
	logical, intent(in) :: budget_flag
	character(100) :: folder_aux

	print*,' '; print*,'---------------------------------------------------------------------'
	print*,'Solve Experiment tk and tw'
	print*,'---------------------------------------------------------------------'; print*,' '

	! Save base folder
		folder_aux = Result_Folder

	! Load Benchmark Variables
		call Solve_Benchmark(.false.,.false.)
		! Change flag 
		solving_bench=0

	! Change taxes
		tauK = tauK_bench
		tauW_at = 0.0303366758418838_dp ! 0.0119357175267747_dp ! 0.02_dp 

	! Change Folder 
		if (budget_flag) then 
			Result_Folder = trim(folder_aux)//'Tax_Reform_tktw/Budget_Balance_303/'
			print*,' '; print*,'---------------------------------------------------------------------'
			print*,'Adding wealth taxes balancing the budget with labor income taxes'
		else 
			Result_Folder = trim(folder_aux)//'Tax_Reform_tktw/No_Budget_Balance_303/'
			print*,' '; print*,'---------------------------------------------------------------------'
			print*,'Adding wealth taxes without balancing the budget'
		endif 
		call system( 'mkdir -p ' // trim(Result_Folder) )

	! Balance the budget
	if (budget_flag) then 
		GBAR_exp = 0.0_dp 
		DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.01 ) ! as long as the difference is greater than 0.01% continue
		    CALL FIND_DBN_EQ
		    CALL GOVNT_BUDGET_OPT
		    GBAR_exp = GBAR    
		ENDDO
	endif 

	! Solve Model 
		print*,"	Computing equilibrium distribution"
		CALL FIND_DBN_EQ
		print*,"	Computing government spending"
		CALL GOVNT_BUDGET(.true.)
		print*,"	Computing Value Function"
		! CALL COMPUTE_VALUE_FUNCTION_SPLINE
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
		print*,"	Computing Firm Value Function"
		CALL Firm_Value
		print*, " 	Computing After Tax Income"
		CALL Compute_After_Tax_Income
		print*,"	Saving results in text files to be read later"
		CALL Write_Experimental_Results(.true.)

	! Aggregate variable in experimental economy
		GBAR_exp  = GBAR
		QBAR_exp  = QBAR 
		NBAR_exp  = NBAR  
		Y_exp 	  = YBAR
		Ebar_exp  = EBAR
		P_exp     = P
		R_exp	  = R
		wage_exp  = wage
		tauK_exp  = tauK
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

		ValueFunction_exp = ValueFunction
		Bq_Value_exp      = Bq_Value
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime
		V_Pr_exp          = V_Pr 
		V_Pr_nb_exp  	  = V_Pr_nb

		YBAR_C_exp = YBAR_C
		L_C_exp    = L_C
		K_C_exp    = K_C
		YBAR_P_exp = YBAR_P
		L_P_exp    = L_P
		K_P_exp    = K_P

	! Compute moments
		CALL COMPUTE_STATS
	
	! Compute welfare gain between economies
		CALL COMPUTE_WELFARE_GAIN
		CALL COMPUTE_WELFARE_DECOMPOSITION

	! Write experimental results in output.txt
		CALL WRITE_VARIABLES(0)

	! Deallocate variables
		if (allocated(YGRID_t)) then 
			deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )
		endif


end Subroutine Solve_Experiment_tktw


!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Experiment_Fixed_Policy_Functions(compute_exp_pf,Simul_Switch)
	use parameters
	use global 
	use programfunctions
	use GKK_Stats
	use Simulation_Module
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: compute_exp_pf, Simul_Switch
	integer :: aa, age1, a1, z1, lambda1, e1, x1

	!====================================================================================================
	PRINT*,''
	Print*,'--------------- SOLVING EXPERIMENT WITH BEST PARAMETERS -----------------'
	PRINT*,''
	print*,'Wealth Tax Economy'
	
	! Experiment economy
		solving_bench=0
	! Set capital taxes to zero
		tauK = 0.0_DP
		tauWmin_at= 0.005_DP
		tauWinc_at= 0.005_DP
	! Set Y_a_threshold
		write(*,*) "Y_a threshold is set to a proportion of the mean wealth under current distribution"
		!Y_a_threshold = 0.0_dp ! big_p   !8.1812138704441200
		call Find_TauW_Threshold(DBN_bench,W_bench)  
		Y_a_threshold = Threshold_Factor*Ebar_bench !0.75_dp
		Wealth_factor = Y_a_threshold/W_bench

	if (compute_exp_pf) then 
		! Find wealth taxes that balances budget
		print*, "	Computing Wealth Tax to balance the budget"
			! Set initial value for G in experimental economy and for wealth taxes
			GBAR_exp = 0.0_DP
			tauW_bt  = tauWmin_bt
			tauW_at  = tauWmin_at
			tauWindx = 0.0_DP
			! Solve for the model increasing wealth taxes until revenue is enough to finance G_benchamark
			! OPEN(UNIT=19, FILE=trim(Result_Folder)//'Laffer_Curve.txt', STATUS='replace')
			! DO aa = 1,101
			DO WHILE (GBAR_exp .lt. GBAR_bench)
				! Set old G and new value of tauW
				GBAR_exp_old = GBAR_exp
				tauW_bt = tauWmin_bt + tauWindx * tauWinc_bt
				tauW_at = tauWmin_at + tauWindx * tauWinc_at
				! Solve the model
				CALL FIND_DBN_EQ_PF
				CALL GOVNT_BUDGET(.true.)

				! Get new G
				GBAR_exp = GBAR 
				! Iteratioins  
				tauWindx = tauWindx + 1.0_DP   
				write(*,*) "Bracketing GBAR: tauW_bt=", tauW_bt*100, "And tauW_at=", tauW_at*100
				print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
				print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
				! WRITE(UNIT=19, FMT=*) 'tau_W',tauW_at,'GBAR',GBAR_exp,'Assets',sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid )
			ENDDO
			! CLOSE(UNIT=19)
			! print*,' Program stoped by Sergio'
			! STOP

			! Set tauW as weighted average of point in  the grid to balance budget more precisely
				tauW_up_bt  = tauW_bt
				tauW_low_bt = tauW_bt  -  tauWinc_bt
				tauW_bt     = tauW_low_bt + tauWinc_bt * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
				tauW_up_at  = tauW_at
				tauW_low_at = tauW_at  -  tauWinc_at  
				tauW_at     = tauW_low_at + tauWinc_at * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
				print*,''
				print*,'GBAR bracketed by taxes:'
				print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
				print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
				print*,''

			! Solve (again) experimental economy
				CALL FIND_DBN_EQ_PF
				CALL GOVNT_BUDGET(.true.)

			! Find tauW that exactly balances the budget (up to precisioin 0.1) using bisection
				GBAR_exp = GBAR
				print*,"Gbar at midpoint of bracket and GBAR at benchmark"
				print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
				print*,''
				print*,'Bisection for TauW:'
				DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.001 ) ! as long as the difference is greater than 0.1% continue
				    if (GBAR_exp .gt. GBAR_bench ) then
				        tauW_up_bt  = tauW_bt 
				        tauW_up_at  = tauW_at 
				    else
				        tauW_low_bt = tauW_bt
				        tauW_low_at = tauW_at
				    endif
				    tauW_bt = (tauW_low_bt + tauW_up_bt)/2.0_DP
				    tauW_at = (tauW_low_at + tauW_up_at)/2.0_DP
				    CALL FIND_DBN_EQ_PF
				    CALL GOVNT_BUDGET(.true.)
				    GBAR_exp = GBAR
				    print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
					print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
					print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
					print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
				ENDDO

	    	! Compute aggregates with current distribution
	        QBAR =0.0
	        NBAR =0.0
	        DO x1=1,nx
	        DO age1=1,MaxAge
	        DO z1=1,nz
	        DO a1=1,na
	        DO lambda1=1,nlambda
	        DO e1=1, ne
	             QBAR= QBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * ( xz_grid(x1,z1) * K_mat(a1,z1,x1) )**mu
	             NBAR= NBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * eff_un(age1, lambda1, e1) * Hours(age1, a1, z1, lambda1,e1,x1)
	        ENDDO
	        ENDDO
	        ENDDO
	        ENDDO    
	        ENDDO    
	        ENDDO    
	    
	        QBAR = ( QBAR)**(1.0_DP/mu)                
	        YBAR = QBAR ** alpha * NBAR **(1.0_DP-alpha)
	        Ebar = wage  * NBAR  * sum(pop)/sum(pop(1:RetAge-1))

		! Compute value function and store policy functions, value function and distribution in file
		! CALL COMPUTE_VALUE_FUNCTION_SPLINE 
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
		CALL Firm_Value

	endif 

	
	CALL Write_Experimental_Results(compute_exp_pf)
	CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
	K_mat  = K_Matrix(R,P)
	Pr_mat = Profit_Matrix(R,P)
	CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
	CALL ComputeLaborUnits(EBAR,wage)
	CALL GOVNT_BUDGET(.true.)


	! Aggregate variable in experimental economy
		GBAR_exp  = GBAR
		QBAR_exp  = QBAR 
		NBAR_exp  = NBAR  
		Y_exp 	  = YBAR
		Ebar_exp  = EBAR
		P_exp     = P
		R_exp	  = R
		wage_exp  = wage
		tauK_exp  = tauK
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

		ValueFunction_exp = ValueFunction
		Bq_Value_exp 	  = Bq_Value
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime
		V_Pr_exp          = V_Pr 
		V_Pr_nb_exp  	  = V_Pr_nb

	! Compute moments
	CALL COMPUTE_STATS
	
	! Compute welfare gain between economies
	CALL COMPUTE_WELFARE_GAIN

	! Write experimental results in output.txt
	CALL WRITE_VARIABLES(0)
	if ((Simul_Switch)) then 
	 	print*,"	Experiment Simulation"
		CALL SIMULATION(solving_bench)
	endif


	print*,'---------------------------'
	print*,''
	print*,'Output Gain Prct=', 100.0_DP*(Y_exp/Y_bench-1.0) 
	print*,''
	print*,'---------------------------'

	print*," "
	print*,"Wealth_factor=",Wealth_factor
	print*," "


	! Deallocate variables
		deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )

	print*,"	Efficiency Computation"
		CALL Hsieh_Klenow_Efficiency(solving_bench)


end Subroutine Solve_Experiment_Fixed_Policy_Functions


!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Experiment_Fixed_PF_Interp (compute_exp_pf_interp,Simul_Switch)
	use parameters
	use global 
	use programfunctions
	use GKK_Stats
	use Simulation_Module
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: compute_exp_pf_interp, Simul_Switch
	integer :: aa, age1, a1, z1, lambda1, e1, x1
	REAL(DP), DIMENSION(na,nz,nx) :: YGRID_bench

	!====================================================================================================
	PRINT*,''
	Print*,'--------------- Get YGRID from benchmark -----------------'
	PRINT*,''
	CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
	CALL FORM_Y_MB_GRID(YGRID_bench,MBGRID,YGRID_t,MBGRID_t)
	! deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )

	!====================================================================================================
	PRINT*,''
	Print*,'--------------- SOLVING EXPERIMENT WITH BEST PARAMETERS -----------------'
	PRINT*,''
	print*,'Wealth Tax Economy'
	
	! Experiment economy
		solving_bench=0
	! Set capital taxes to zero
		tauK = 0.0_DP
		tauWmin_at= 0.010_DP
		tauWinc_at= 0.005_DP
	! Set Y_a_threshold
		write(*,*) "Y_a threshold is set to a proportion of the mean wealth under current distribution"
		!Y_a_threshold = 0.0_dp ! big_p   !8.1812138704441200
		call Find_TauW_Threshold(DBN_bench,W_bench)  
		Y_a_threshold = Threshold_Factor*Ebar_bench !0.75_dp
		Wealth_factor = Y_a_threshold/W_bench

	if (compute_exp_pf_interp) then 
		! Find wealth taxes that balances budget
		print*, "	Computing Wealth Tax to balance the budget"
			! Set initial value for G in experimental economy and for wealth taxes
			GBAR_exp = 0.0_DP
			tauW_bt  = tauWmin_bt
			tauW_at  = tauWmin_at
			tauWindx = 0.0_DP
			! Solve for the model increasing wealth taxes until revenue is enough to finance G_benchamark
			! OPEN(UNIT=19, FILE=trim(Result_Folder)//'Laffer_Curve.txt', STATUS='replace')
			! DO aa = 1,41
			DO WHILE (GBAR_exp .lt. GBAR_bench)
				! Set old G and new value of tauW
				GBAR_exp_old = GBAR_exp
				tauW_bt = tauWmin_bt + tauWindx * tauWinc_bt
				tauW_at = tauWmin_at + tauWindx * tauWinc_at
				! Solve the model
				CALL FIND_DBN_EQ_PF_Interp(YGRID_bench)
				CALL GOVNT_BUDGET(.true.)

				! Get new G
				GBAR_exp = GBAR 
				! Iteratioins  
				tauWindx = tauWindx + 1.0_DP   
				write(*,*) "Bracketing GBAR: tauW_bt=", tauW_bt*100, "And tauW_at=", tauW_at*100
				print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
				print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
				! WRITE(UNIT=19, FMT=*) 'tau_W',tauW_at,'GBAR',GBAR_exp,'Assets',sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid )
			ENDDO
			! CLOSE(UNIT=19)
			! print*,' Program stoped by Sergio'
			! STOP

			! Set tauW as weighted average of point in  the grid to balance budget more precisely
				tauW_up_bt  = tauW_bt
				tauW_low_bt = tauW_bt  -  tauWinc_bt
				tauW_bt     = tauW_low_bt + tauWinc_bt * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
				tauW_up_at  = tauW_at
				tauW_low_at = tauW_at  -  tauWinc_at  
				tauW_at     = tauW_low_at + tauWinc_at * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
				print*,''
				print*,'GBAR bracketed by taxes:'
				print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
				print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
				print*,''

			! Solve (again) experimental economy
				CALL FIND_DBN_EQ_PF_Interp(YGRID_bench)
				CALL GOVNT_BUDGET(.true.)

			! Find tauW that exactly balances the budget (up to precisioin 0.1) using bisection
				GBAR_exp = GBAR
				print*,"Gbar at midpoint of bracket and GBAR at benchmark"
				print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
				print*,''
				print*,'Bisection for TauW:'
				DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.001 ) ! as long as the difference is greater than 0.1% continue
				    if (GBAR_exp .gt. GBAR_bench ) then
				        tauW_up_bt  = tauW_bt 
				        tauW_up_at  = tauW_at 
				    else
				        tauW_low_bt = tauW_bt
				        tauW_low_at = tauW_at
				    endif
				    tauW_bt = (tauW_low_bt + tauW_up_bt)/2.0_DP
				    tauW_at = (tauW_low_at + tauW_up_at)/2.0_DP
				    CALL FIND_DBN_EQ_PF_Interp(YGRID_bench)
				    CALL GOVNT_BUDGET(.true.)
				    GBAR_exp = GBAR
				    print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
					print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
					print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
					print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
				ENDDO

	    	! Compute aggregates with current distribution
	        QBAR =0.0
	        NBAR =0.0
	        DO x1=1,nx
	        DO age1=1,MaxAge
	        DO z1=1,nz
	        DO a1=1,na
	        DO lambda1=1,nlambda
	        DO e1=1, ne
	             QBAR= QBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * ( xz_grid(x1,z1) * K_mat(a1,z1,x1) )**mu
	             NBAR= NBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * eff_un(age1, lambda1, e1) * Hours(age1, a1, z1, lambda1,e1,x1)
	        ENDDO
	        ENDDO
	        ENDDO
	        ENDDO    
	        ENDDO    
	        ENDDO    
	    
	        QBAR = ( QBAR)**(1.0_DP/mu)                
	        YBAR = QBAR ** alpha * NBAR **(1.0_DP-alpha)
	        Ebar = wage  * NBAR  * sum(pop)/sum(pop(1:RetAge-1))

		! Compute value function and store policy functions, value function and distribution in file
		! CALL COMPUTE_VALUE_FUNCTION_SPLINE 
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
		CALL Firm_Value

	endif 


	CALL Write_Experimental_Results(compute_exp_pf_interp)
	K_mat  = K_Matrix(R,P)
	Pr_mat = Profit_Matrix(R,P)
	CALL ComputeLaborUnits(EBAR,wage)
	CALL GOVNT_BUDGET(.true.)
	
	! Aggregate variable in experimental economy
		GBAR_exp  = GBAR
		QBAR_exp  = QBAR 
		NBAR_exp  = NBAR  
		Y_exp 	  = YBAR
		Ebar_exp  = EBAR
		P_exp     = P
		R_exp	  = R
		wage_exp  = wage
		tauK_exp  = tauK
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

		ValueFunction_exp = ValueFunction
		Bq_Value_exp 	  = Bq_Value
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime
		V_Pr_exp          = V_Pr 
		V_Pr_nb_exp  	  = V_Pr_nb

	! Compute moments
	! CALL COMPUTE_STATS
	
	! Compute welfare gain between economies
	CALL COMPUTE_WELFARE_GAIN

	! Write experimental results in output.txt
	CALL WRITE_VARIABLES(0)
	if ((Simul_Switch)) then 
	 	print*,"	Experiment Simulation"
		CALL SIMULATION(solving_bench)
	endif


	print*,'---------------------------'
	print*,''
	print*,'Output Gain Prct=', 100.0_DP*(Y_exp/Y_bench-1.0) 
	print*,''
	print*,'---------------------------'

	print*," "
	print*,"Wealth_factor=",Wealth_factor
	print*," "


	! Deallocate variables
		deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )

	! print*,"	Efficiency Computation"
	! 	CALL Hsieh_Klenow_Efficiency(solving_bench)


end Subroutine Solve_Experiment_Fixed_PF_Interp


!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Experiment_Fixed_PF_Prices(compute_exp_prices,Simul_Switch)
	use parameters
	use global 
	use programfunctions
	use GKK_Stats
	use Simulation_Module
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: compute_exp_prices, Simul_Switch

	!====================================================================================================
	PRINT*,''
	Print*,'--------------- SOLVING EXPERIMENT WITH BEST PARAMETERS -----------------'
	PRINT*,''
	print*,'Wealth Tax Economy'
	
	! Experiment economy
		solving_bench=0
	! Set capital taxes to zero
		tauK = 0.0_DP
	! Set Y_a_threshold
		write(*,*) "Y_a threshold is set to a proportion of the mean wealth under current distribution"
		!Y_a_threshold = 0.0_dp ! big_p   !8.1812138704441200
		call Find_TauW_Threshold(DBN_bench,W_bench)  
		Y_a_threshold = Threshold_Factor*Ebar_bench !0.75_dp
		Wealth_factor = Y_a_threshold/W_bench

	! if (compute_exp_prices) then 
	! 	! Find wealth taxes that balances budget
	! 	print*, "	Computing Wealth Tax to balance the budget"
	! 		! Set initial value for G in experimental economy and for wealth taxes
	! 		GBAR_exp = 0.0_DP
	! 		tauW_bt  = tauWmin_bt
	! 		tauW_at  = tauWmin_at
	! 		tauWindx = 0.0_DP
	! 		! Solve for the model increasing wealth taxes until revenue is enough to finance G_benchamark
	! 		DO WHILE (GBAR_exp .lt. GBAR_bench)
	! 			! Set old G and new value of tauW
	! 			GBAR_exp_old = GBAR_exp
	! 			tauW_bt = tauWmin_bt + tauWindx * tauWinc_bt
	! 			tauW_at = tauWmin_at + tauWindx * tauWinc_at
	! 			! Solve the model
	! 			CALL FIND_DBN_EQ_PF_Prices
	! 			CALL GOVNT_BUDGET(.true.)

	! 			! Get new G
	! 			GBAR_exp = GBAR 
	! 			! Iteratioins  
	! 			tauWindx = tauWindx + 1.0_DP   
	! 			write(*,*) "Bracketing GBAR: tauW_bt=", tauW_bt*100, "And tauW_at=", tauW_at*100
	! 			print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
	! 			print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
	! 		ENDDO

	! 		! Set tauW as weighted average of point in  the grid to balance budget more precisely
	! 			tauW_up_bt  = tauW_bt
	! 			tauW_low_bt = tauW_bt  -  tauWinc_bt
	! 			tauW_bt     = tauW_low_bt + tauWinc_bt * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
	! 			tauW_up_at  = tauW_at
	! 			tauW_low_at = tauW_at  -  tauWinc_at  
	! 			tauW_at     = tauW_low_at + tauWinc_at * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
	! 			print*,''
	! 			print*,'GBAR bracketed by taxes:'
	! 			print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
	! 			print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
	! 			print*,''

	! 		! Solve (again) experimental economy
	! 			CALL FIND_DBN_EQ_PF_Prices
	! 			CALL GOVNT_BUDGET(.true.)

	! 		! Find tauW that exactly balances the budget (up to precisioin 0.1) using bisection
	! 			GBAR_exp = GBAR
	! 			print*,"Gbar at midpoint of bracket and GBAR at benchmark"
	! 			print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
	! 			print*,''
	! 			print*,'Bisection for TauW:'
	! 			DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.01 ) ! as long as the difference is greater than 0.1% continue
	! 			    if (GBAR_exp .gt. GBAR_bench ) then
	! 			        tauW_up_bt  = tauW_bt 
	! 			        tauW_up_at  = tauW_at 
	! 			    else
	! 			        tauW_low_bt = tauW_bt
	! 			        tauW_low_at = tauW_at
	! 			    endif
	! 			    tauW_bt = (tauW_low_bt + tauW_up_bt)/2.0_DP
	! 			    tauW_at = (tauW_low_at + tauW_up_at)/2.0_DP
	! 			    CALL FIND_DBN_EQ_PF_Prices
	! 			    CALL GOVNT_BUDGET(.true.)
	! 			    GBAR_exp = GBAR
	! 			    print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
	! 				print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
	! 				print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
	! 				print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
	! 			ENDDO

	! 	! Compute value function and store policy functions, value function and distribution in file
	! 	! CALL COMPUTE_VALUE_FUNCTION_SPLINE 
	! 	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
	! 	CALL Firm_Value

	! endif 
	
	! CALL Write_Experimental_Results(compute_exp_prices)
	CALL Write_Experimental_Results(.false.)
	CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
	K_mat  = K_Matrix(R,P)
	Pr_mat = Profit_Matrix(R,P)
	CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
	CALL ComputeLaborUnits(EBAR,wage)
	CALL GOVNT_BUDGET(.true.)


	! Aggregate variable in experimental economy
		GBAR_exp  = GBAR
		QBAR_exp  = QBAR 
		NBAR_exp  = NBAR  
		Y_exp 	  = YBAR
		Ebar_exp  = EBAR
		P_exp     = P
		R_exp	  = R
		wage_exp  = wage
		tauK_exp  = tauK
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

		ValueFunction_exp = ValueFunction
		Bq_Value_exp	  = Bq_Value
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime
		V_Pr_exp          = V_Pr 
		V_Pr_nb_exp  	  = V_Pr_nb

	! Compute moments
	CALL COMPUTE_STATS
	
	! Compute welfare gain between economies
	CALL COMPUTE_WELFARE_GAIN

	! Write experimental results in output.txt
	CALL WRITE_VARIABLES(0)
	if ((Simul_Switch)) then 
	 	print*,"	Experiment Simulation"
		CALL SIMULATION(solving_bench)
	endif
	Call Simulation_Life_Cycle_Patterns(solving_bench)

	print*,'---------------------------'
	print*,''
	print*,'Output Gain Prct=', 100.0_DP*(Y_exp/Y_bench-1.0) 
	print*,''
	print*,'---------------------------'

	print*," "
	print*,"Wealth_factor=",Wealth_factor
	print*," "


	! Deallocate variables
		deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )

	! print*,"	Efficiency Computation"
	! 	CALL Hsieh_Klenow_Efficiency(solving_bench)


end Subroutine Solve_Experiment_Fixed_PF_Prices

!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Experiment_Fixed_Prices(compute_exp_prices,Simul_Switch,Fixed_W,Fixed_P,Fixed_R)
	use parameters
	use global 
	use programfunctions
	use GKK_Stats
	use Simulation_Module
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: compute_exp_prices, Simul_Switch, Fixed_W, Fixed_P, Fixed_R
	real(dp) :: GBAR_low, GBAR_up

	!====================================================================================================
	PRINT*,''
	Print*,'--------------- SOLVING EXPERIMENT WITH BEST PARAMETERS -----------------'
	PRINT*,''
	print*,'Wealth Tax Economy'
	
	! Experiment economy
		solving_bench=0
	! Set capital taxes to zero
		tauK = 0.0_DP
		if (Fixed_P) then 
		tauWmin_at= 0.0050_DP
		tauWinc_at= 0.005_DP
		endif 
	! Set Y_a_threshold
		write(*,*) "Y_a threshold is set to a proportion of the mean wealth under current distribution"
		!Y_a_threshold = 0.0_dp ! big_p   !8.1812138704441200
		call Find_TauW_Threshold(DBN_bench,W_bench)  
		Y_a_threshold = Threshold_Factor*Ebar_bench !0.75_dp
		Wealth_factor = Y_a_threshold/W_bench

	! if (compute_exp_prices) then 
	! 	! Find wealth taxes that balances budget
	! 	print*, "	Computing Wealth Tax to balance the budget"
	! 		! Set initial value for G in experimental economy and for wealth taxes
	! 		GBAR_exp = 0.0_DP
	! 		tauW_bt  = tauWmin_bt
	! 		tauW_at  = tauWmin_at
	! 		tauWindx = 0.0_DP
	! 		! Solve for the model increasing wealth taxes until revenue is enough to finance G_benchamark
	! 		DO WHILE (GBAR_exp .lt. GBAR_bench)
	! 			! Set old G and new value of tauW
	! 			GBAR_exp_old = GBAR_exp
	! 			tauW_bt = tauWmin_bt + tauWindx * tauWinc_bt
	! 			tauW_at = tauWmin_at + tauWindx * tauWinc_at
	! 			! Solve the model
	! 			CALL FIND_DBN_EQ_Prices(Fixed_W,Fixed_P,Fixed_R)
	! 			CALL GOVNT_BUDGET(.true.)

	! 			! Get new G
	! 			GBAR_exp = GBAR 
	! 			! Iteratioins  
	! 			tauWindx = tauWindx + 1.0_DP   
	! 			write(*,*) "Bracketing GBAR: tauW_bt=", tauW_bt*100, "And tauW_at=", tauW_at*100
	! 			print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
	! 			print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
	! 			! Set GBAR for low and high point in grid
	! 			GBAR_low     = GBAR_exp_old 
	! 			GBAR_up      = GBAR_exp 
	! 		ENDDO

	! 		! Set tauW as weighted average of point in  the grid to balance budget more precisely
	! 			tauW_up_bt  = tauW_bt
	! 			tauW_low_bt = tauW_bt  -  tauWinc_bt
	! 			tauW_bt     = tauW_low_bt + tauWinc_bt * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
	! 			tauW_up_at  = tauW_at
	! 			tauW_low_at = tauW_at  -  tauWinc_at  
	! 			tauW_at     = tauW_low_at + tauWinc_at * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
	! 			print*,''
	! 			print*,'GBAR bracketed by taxes:'
	! 			print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
	! 			print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
	! 			print*,''

	! 		! Solve (again) experimental economy
	! 			CALL FIND_DBN_EQ_Prices(Fixed_W,Fixed_P,Fixed_R)
	! 			CALL GOVNT_BUDGET(.true.)

	! 		! Find tauW that exactly balances the budget (up to precisioin 0.1) using bisection
	! 			GBAR_exp = GBAR
	! 			print*,"Gbar at midpoint of bracket and GBAR at benchmark"
	! 			print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
	! 			print*,''
	! 			print*,'Bisection for TauW:'
	! 			DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.001 ) ! as long as the difference is greater than 0.1% continue
	! 			    if (GBAR_exp .gt. GBAR_bench ) then
	! 			        tauW_up_bt  = tauW_bt 
	! 			        tauW_up_at  = tauW_at 
	! 			        GBAR_up     = GBAR_exp 
	! 			    else
	! 			        tauW_low_bt = tauW_bt
	! 			        tauW_low_at = tauW_at
	! 			        GBAR_low    = GBAR_exp 
	! 			    endif
	! 			    tauW_bt = tauW_low_bt + (tauW_up_bt-tauW_low_bt)*(GBAR_bench-GBAR_low)/(GBAR_up-GBAR_low)
	! 			    tauW_at = tauW_low_at + (tauW_up_at-tauW_low_at)*(GBAR_bench-GBAR_low)/(GBAR_up-GBAR_low)
	! 			    CALL FIND_DBN_EQ_Prices(Fixed_W,Fixed_P,Fixed_R)
	! 			    CALL GOVNT_BUDGET(.true.)
	! 			    GBAR_exp = GBAR
	! 			    print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
	! 				print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
	! 				print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
	! 				print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
	! 			ENDDO

	! 	! Compute value function and store policy functions, value function and distribution in file
	! 	! CALL COMPUTE_VALUE_FUNCTION_SPLINE 
	! 	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
	! 	CALL Firm_Value

	! endif 
	
	! CALL Write_Experimental_Results(compute_exp_prices)
 	CALL Write_Experimental_Results(.false.)
	CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
	K_mat  = K_Matrix(R,P)
	Pr_mat = Profit_Matrix(R,P)
	CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
	CALL ComputeLaborUnits(EBAR,wage)
	CALL GOVNT_BUDGET(.true.)


	! Aggregate variable in experimental economy
		GBAR_exp  = GBAR
		QBAR_exp  = QBAR 
		NBAR_exp  = NBAR  
		Y_exp 	  = YBAR
		Ebar_exp  = EBAR
		P_exp     = P
		R_exp	  = R
		wage_exp  = wage
		tauK_exp  = tauK
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

		ValueFunction_exp = ValueFunction
		Bq_Value_exp	  = Bq_Value
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime
		V_Pr_exp          = V_Pr 
		V_Pr_nb_exp  	  = V_Pr_nb

	! Compute moments
	CALL COMPUTE_STATS
	
	! Compute welfare gain between economies
	CALL COMPUTE_WELFARE_GAIN

	! Write experimental results in output.txt
	CALL WRITE_VARIABLES(0)
	if ((Simul_Switch)) then 
	 	print*,"	Experiment Simulation"
		CALL SIMULATION(solving_bench)
	endif
	Call Simulation_Life_Cycle_Patterns(solving_bench)


	print*,'---------------------------'
	print*,''
	print*,'Output Gain Prct=', 100.0_DP*(Y_exp/Y_bench-1.0) 
	print*,''
	print*,'---------------------------'

	print*," "
	print*,"Wealth_factor=",Wealth_factor
	print*," "


	! Deallocate variables
		deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )

	! print*,"	Efficiency Computation"
	! 	CALL Hsieh_Klenow_Efficiency(solving_bench)


end Subroutine Solve_Experiment_Fixed_Prices



!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Experiment_Fixed_Prices_and_Taxes
	use parameters
	use global 
	use programfunctions
	use GKK_Stats
	use Simulation_Module
	use Toolbox
	use omp_lib
	implicit none 
	REAL(DP), DIMENSION(na,nz,nx) :: YGRID_exp
	integer :: aa, age1, a1, z1, lambda1, e1, x1

	!====================================================================================================
	! Get Benchmark Values
		call Solve_Benchmark(.false.,.false.)

	!====================================================================================================
	! Get Experiment Values
		call Solve_Experiment(.false.,.false.)
		! Set YGRID for interpolation
			CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
			CALL FORM_Y_MB_GRID(YGRID_exp,MBGRID,YGRID_t,MBGRID_t)
			CALL ComputeLaborUnits(EBAR,wage)
			! deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )

	!====================================================================================================
	! Change folder		
		Result_Folder = trim(Result_Folder)//'Exp_PF_Interp_fixed_price_and_tax/'
		call system( 'mkdir -p ' // trim(Result_Folder) )			


	!====================================================================================================
	! Choose values for prices and policy functions
		P_exp        = P_bench
		R_exp	     = R_bench
		wage_exp     = wage_bench
		tauK_exp     = tauK_bench
		tauPL_exp    = tauPL_bench
		psi_exp      = psi_bench
		tauw_bt_exp  = tauW_bt_bench
		tauw_at_exp  = tauW_at_bench 
		Y_a_threshold_exp = Y_a_threshold_bench 

		P 			 = P_bench 
		R 			 = R_bench 
		wage 		 = wage_bench 
		tauK 		 = tauK_bench 
		tauPL 		 = tauPL_bench 
		tauw_bt      = tauw_bt_bench 
		tauw_at      = tauw_at_bench 

		Cons_bench   = Cons_exp 
		Hours_bench  = Hours_exp
		Aprime_bench = Aprime_exp


	!====================================================================================================
	! Solve Equilibrium with policy function interpolation
		CALL FIND_DBN_EQ_PF_Interp(YGRID_exp)
		CALL GOVNT_BUDGET(.true.)


    	! Compute aggregates with current distribution
	        QBAR =0.0
	        NBAR =0.0
	        DO x1=1,nx
	        DO age1=1,MaxAge
	        DO z1=1,nz
	        DO a1=1,na
	        DO lambda1=1,nlambda
	        DO e1=1, ne
	             QBAR= QBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * ( xz_grid(x1,z1) * K_mat(a1,z1,x1) )**mu
	             NBAR= NBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * eff_un(age1, lambda1, e1) * Hours(age1, a1, z1, lambda1,e1,x1)
	        ENDDO
	        ENDDO
	        ENDDO
	        ENDDO    
	        ENDDO    
	        ENDDO    
	    
	        QBAR = ( QBAR)**(1.0_DP/mu)                
	        YBAR = QBAR ** alpha * NBAR **(1.0_DP-alpha)
	        Ebar = wage  * NBAR  * sum(pop)/sum(pop(1:RetAge-1))

		! Compute value function and store policy functions, value function and distribution in file
			! CALL COMPUTE_VALUE_FUNCTION_SPLINE 
			CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
			CALL Firm_Value


	!====================================================================================================
	! Compute stats and tables
		CALL Write_Experimental_Results(.true.)
		! CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
		K_mat  = K_Matrix(R,P)
		Pr_mat = Profit_Matrix(R,P)
		! CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
		! CALL ComputeLaborUnits(EBAR,wage)


		! Aggregate variable in experimental economy
		GBAR_exp  = GBAR
		QBAR_exp  = QBAR 
		NBAR_exp  = NBAR  
		Y_exp 	  = YBAR
		Ebar_exp  = EBAR
		P_exp     = P
		R_exp	  = R
		wage_exp  = wage
		tauK_exp  = tauK
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

		ValueFunction_exp = ValueFunction
		Bq_Value_exp	  = Bq_Value
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime
		V_Pr_exp          = V_Pr 
		V_Pr_nb_exp  	  = V_Pr_nb

		! Compute moments
		CALL COMPUTE_STATS
		! Compute welfare gain between economies
		CALL COMPUTE_WELFARE_GAIN


	! Deallocate variables
		deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )

end Subroutine Solve_Experiment_Fixed_Prices_and_Taxes




!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Tax_Reform_Decomposition
	use parameters
	use global 
	use programfunctions
	use GKK_Stats
	use Simulation_Module
	use Toolbox
	use omp_lib
	implicit none 
	REAL(DP), DIMENSION(na,nz,nx) :: YGRID_exp
	REAL(DP) :: P_aux, R_aux, wage_aux
	character(100) :: aux_folder

	aux_folder = Result_Folder

	!====================================================================================================
	! Get Benchmark Values
		call Solve_Benchmark(.false.,.false.)
		! Set YGRID for interpolation
			CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
			CALL FORM_Y_MB_GRID(YGRID_exp,MBGRID,YGRID_t,MBGRID_t)
			CALL ComputeLaborUnits(EBAR,wage)
			deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )

	!====================================================================================================
	! Get Experiment Values
		Result_Folder = trim(aux_folder)//'Tax_Reform/'
		call Solve_Experiment(.false.,.false.)
			P_aux = P_exp
			R_aux = R_exp 
			wage_aux = wage_exp

	!====================================================================================================
	! Change only taxes
	!====================================================================================================


	!====================================================================================================
	! Change folder		
		Result_Folder = trim(aux_folder)//'Tax_Reform_Decomposition/Only_Taxes/'
		call system( 'mkdir -p ' // trim(Result_Folder) )			

	!====================================================================================================
	! Choose values for prices and policy functions
		P_exp        = P_bench
		R_exp	     = R_bench
		wage_exp     = wage_bench
		
		P 			 = P_bench
		R 			 = R_bench 
		wage 		 = wage_bench 
		tauK 		 = tauK_exp
		tauPL 		 = tauPL_exp
		tauw_bt      = tauw_bt_exp
		tauw_at      = tauw_at_exp

		CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)

	! Solve Interpolated Economy
	call Solve_Interpolated_Economy(YGRID_exp)

	!====================================================================================================
	! Change taxes and prices
	!====================================================================================================


	!====================================================================================================
	! Change folder		
		Result_Folder = trim(aux_folder)//'Tax_Reform_Decomposition/Taxes_and_Prices/'
		call system( 'mkdir -p ' // trim(Result_Folder) )			

	!====================================================================================================
	! Choose values for prices and policy functions
		P_exp        = P_aux
		R_exp	     = R_aux
		wage_exp     = wage_aux
		
		P 			 = P_exp
		R 			 = R_exp 
		wage 		 = wage_exp

		CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)

	! Solve Interpolated Economy
	call Solve_Interpolated_Economy(YGRID_exp)

end Subroutine Solve_Tax_Reform_Decomposition


Subroutine Solve_Interpolated_Economy(YGRID_exp)
	use parameters
	use global 
	use programfunctions
	use GKK_Stats
	use Simulation_Module
	use Toolbox
	use omp_lib
	implicit none 
	REAL(DP), DIMENSION(na,nz,nx) :: YGRID_exp
	INTEGER  :: aa, age1, a1, z1, lambda1, e1, x1

	!====================================================================================================
	! Solve Equilibrium with policy function interpolation
	
		CALL FIND_DBN_EQ_PF_Interp(YGRID_exp)
		CALL GOVNT_BUDGET(.true.)

    	! Compute aggregates with current distribution
	        QBAR =0.0
	        NBAR =0.0
	        DO x1=1,nx
	        DO age1=1,MaxAge
	        DO z1=1,nz
	        DO a1=1,na
	        DO lambda1=1,nlambda
	        DO e1=1, ne
	             QBAR= QBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * ( xz_grid(x1,z1) * K_mat(a1,z1,x1) )**mu
	             NBAR= NBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * eff_un(age1, lambda1, e1) * Hours(age1, a1, z1, lambda1,e1,x1)
	        ENDDO
	        ENDDO
	        ENDDO
	        ENDDO    
	        ENDDO    
	        ENDDO    
	    
	        QBAR = ( QBAR)**(1.0_DP/mu)                
	        YBAR = QBAR ** alpha * NBAR **(1.0_DP-alpha)
	        Ebar = wage  * NBAR  * sum(pop)/sum(pop(1:RetAge-1))

		! Compute value function and store policy functions, value function and distribution in file
			! CALL COMPUTE_VALUE_FUNCTION_SPLINE 
			CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
			CALL Firm_Value



	!====================================================================================================
	! Compute stats and tables
		CALL Write_Experimental_Results(.true.)
		! CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
		K_mat  = K_Matrix(R,P)
		Pr_mat = Profit_Matrix(R,P)
		! CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
		! CALL ComputeLaborUnits(EBAR,wage)


		! Aggregate variable in experimental economy
		GBAR_exp  = GBAR
		QBAR_exp  = QBAR 
		NBAR_exp  = NBAR  
		Y_exp 	  = YBAR
		Ebar_exp  = EBAR
		P_exp     = P
		R_exp	  = R
		wage_exp  = wage
		tauK_exp  = tauK
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

		ValueFunction_exp = ValueFunction
		Bq_Value_exp	  = Bq_Value
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime
		V_Pr_exp          = V_Pr 
		V_Pr_nb_exp  	  = V_Pr_nb

		! Compute moments
		CALL COMPUTE_STATS
	
		! Compute welfare gain between economies
		CALL COMPUTE_WELFARE_GAIN

	! Deallocate variables
		! deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )


end Subroutine Solve_Interpolated_Economy


!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Opt_Tax(Opt_Tax_KW,Simul_Switch)
	use parameters 
	use global
	use Opt_Tax_Parameters
	use Opt_Tax_Functions
	use programfunctions
	use GKK_Stats
	use Simulation_Module
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: Opt_Tax_KW,Simul_Switch
	integer  :: tau_grid_min, tau_grid_max, tau_grid_step
	logical  :: read_results, load_seed
	real(dp) :: MeanCons_bench

	! Set benchmark consumption
		MeanCons_bench = MeanCons

	! Set flag for reading results or computing optimal taxes
		read_results = .false.
		load_seed    = .false.


	if (read_results.eqv..false.) then 
 	!====================================================================================================
	PRINT*,''
	Print*,'--------------- SOLVING OPTIMAL TAXES -----------------'
	PRINT*,''
	
	! Experiment economy
		solving_bench=0
	
	! Set initial taxes for finding optimal ones
		tauK     = 0.0_DP
		tauW_at  = 0.0_DP
		Opt_TauK = 0.0_DP
		Opt_TauW = 0.0_DP
		maxbrentvaluet=-10000.0_DP
	
	if (Opt_Tax_KW) then
		print*,''
		print*,'--------------- OPTIMAL CAPITAL TAXES -----------------'
		print*,''
    	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_k_2.txt', STATUS='replace')
    	
    	tau_grid_min  =  05
    	tau_grid_max  =  25
    	tau_grid_step =  1

    	! Set low psi
    	psi = 0.711398150665184
	else
		print*,''
		print*,'--------------- OPTIMAL WEALTH TAXES -----------------'
		print*,''
    	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_w.txt', STATUS='replace')
    	
    	tau_grid_min  = 27
    	tau_grid_max  = 38
    	tau_grid_step = 1

    	! Set Y_a_threshold
			call Find_TauW_Threshold(DBN_bench,W_bench)  
			Y_a_threshold = Threshold_Factor*Ebar_bench !0.75_dp
			Wealth_factor = Y_a_threshold/W_bench

		! Set high psi
		psi = 0.80_dp
	endif 
    	WRITE(UNIT=77, FMT=*) 'tauK ', 'tauW_at ', 'psi ', 'GBAR_K/Tax_Rev_bench ', &
		      & 'KBAR ','QBAR ','TFP ','NBAR ','YBAR ','Y_Growth ', 'CBAR ','C_Growth ', 'wage ','R ', &
		      & 'Wealth_Output ', 'prct1_wealth ' , 'prct10_wealth ', 'Std_Log_Earnings ', 'mean_hours ', &
	      	  & 'GBAR ', 'GBAR_K ', 'GBAR_W ', 'GBAR_L ', 'GBAR_C ','Tot_Cap_Inc ', &
	      	  & 'Av_Util_Pop ', 'Av_Util_NB ', 'brentvaluet '
      	CLOSE (unit=77) 

	! Load results form file for re-starts of the code
	if (load_seed) then
		CALL Write_Experimental_Results(.false.)
    	! psi = 1.0_dp-0.02811_dp
	endif 

	print*,'	Optimal Tax Loop'
	do tauindx = tau_grid_min,tau_grid_max,tau_grid_step
		if (Opt_Tax_KW) then 
			tauK        = real(tauindx,8)/100_DP
            brentvaluet = - EQ_WELFARE_GIVEN_TauK(tauK)
		else 
			tauw_at     = real(tauindx,8)/1000_DP
            brentvaluet = - EQ_WELFARE_GIVEN_TauW(tauW_at)
		endif 

		! Allocate variables
		CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)

       ! Aggregate variable in experimental economy
			GBAR_exp  = GBAR
			QBAR_exp  = QBAR 
			NBAR_exp  = NBAR  
			Y_exp 	  = YBAR
			Ebar_exp  = EBAR
			P_exp     = P
			R_exp	  = R
			wage_exp  = wage
			tauK_exp  = tauK
			tauPL_exp = tauPL
			psi_exp   = psi
			DBN_exp   = DBN1
			tauw_bt_exp = tauW_bt
			tauw_at_exp = tauW_at
			Y_a_threshold_exp = Y_a_threshold

			ValueFunction_exp = ValueFunction
			Bq_Value_exp	  = Bq_Value
			Cons_exp          = Cons           
			Hours_exp         = Hours
			Aprime_exp        = Aprime 

		! Compute moments
			CALL COMPUTE_STATS
			CALL GOVNT_BUDGET(.true.)
			
		! Compute welfare gain between economies
			CALL COMPUTE_WELFARE_GAIN

		! Write experimental results in output.txt
			CALL WRITE_VARIABLES(0)

		! Update optimal tax
	    if (brentvaluet .gt. maxbrentvaluet) then
	        maxbrentvaluet = brentvaluet
			OPT_tauK = tauK
			OPT_tauW = tauW_at
			OPT_psi  = psi

			! Save variables 
	    	Call Write_Experimental_Results(.true.)
		endif

		! Print Results 
			print*,' ';print*,'------------------------------------------------------------------------------'
		    print '(A,F7.3,X,A,F7.4,X,A,F7.3,X,A,F7.3,X,A,F7.3,X,A,F12.6,X,A,F12.6)', &
		    	  & 'tauK=', 100.0_dp*tauK,'tauW=', 100.0_dp*tauW_at,'tauL',100.0_dp*(1.0_dp-psi),&
		    	  & 'Y/Y_bench=',100.0_dp*(YBAR/Y_bench-1), 'CE2=',Av_Util_NB, 'Brentvalue=',brentvaluet, & 
		    	  & 'Av. Util=', sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
	      	print*,'------------------------------------------------------------------------------';print*,' '

	      	if (Opt_Tax_KW) then 
	      	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_k_2.txt', STATUS='old', POSITION='append')
	      	else 
	      	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_w.txt', STATUS='old', POSITION='append')
	      	endif 
		    WRITE  (UNIT=77, FMT=*) tauK, tauW_at, psi, (GBAR_K+GBAR_W)/(GBAR_bench+SSC_Payments_bench), & 
			      &  MeanWealth, QBAR, QBAR/MeanWealth,NBAR, &
			      &  YBAR, 100.0_DP*(Y_exp/Y_bench-1.0),MeanCons,100.0_DP*(MeanCons/MeanCons_bench-1.0), &
			      &  wage, R,  &
			      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, &
		      	  & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Tot_Cap_Inc, Av_Util_Pop, Av_Util_NB, brentvaluet
	      	CLOSE (unit=77) 

	enddo 
	print*,' ';print*,'------------------------------------------------------------------------------'
	print*,'	End of Optimal Tax Loop'
		tauK    = OPT_tauK
		tauW_at = OPT_tauW
		psi     = OPT_psi
		! Load variables 
	    	Call Write_Experimental_Results(.false.)
	print '(A,F7.3,X,A,F7.3,X,A,F7.3)', &
		&"	Current Optimal Taxes:  tau_K=", 100.0_dp*tauK,"tau_W=",100.0_dp*tauW_at,"tau_L=", 100.0_dp*(1.0_dp-psi)
	print*,'------------------------------------------------------------------------------';print*,' '


	! Search for optimal tax
	print*,'	Optimal Tax Search'
	if (Opt_Tax_KW) then 
		call Find_Opt_Tax(Opt_Tax_KW,Opt_TauK,Opt_TauK-0.01_dp,Opt_TauK+0.01_dp) 
		tauK    = Opt_tauK
	else
		call Find_Opt_Tax(Opt_Tax_KW,Opt_TauW,Opt_TauW-0.001_dp,Opt_TauW+0.001_dp)
		tauW_at = Opt_tauW
	endif 
		OPT_psi  = psi

	print*,' ';print*,'------------------------------------------------------------------------------'
	print*,'	End of Search for Optimal Tax'
	print '(A,F7.3,X,A,F7.3,X,A,F7.3)', &
		&"	Optimal Taxes:  tau_K=", 100.0_dp*tauK,"tau_W=",100.0_dp*tauW_at,"tau_L=", 100.0_dp*(1.0_dp-psi)
	print*,'------------------------------------------------------------------------------';print*,' '


	! Solve model with optimal taxes 
	print*,'	Solving model with optimal taxes'
		CALL FIND_DBN_EQ
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
		CALL Firm_Value

	! Allocate variables
		CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)	

	! Aggregate variable in experimental economy
		GBAR_exp  = GBAR
		QBAR_exp  = QBAR 
		NBAR_exp  = NBAR  
		Y_exp 	  = YBAR
		Ebar_exp  = EBAR
		P_exp     = P
		R_exp	  = R
		wage_exp  = wage
		tauK_exp  = tauK
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

		ValueFunction_exp = ValueFunction
		Bq_Value_exp	  = Bq_Value
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime 

	! Compute moments
		CALL COMPUTE_STATS
		CALL GOVNT_BUDGET(.true.)
		
	! Compute welfare gain between economies
		CALL COMPUTE_WELFARE_GAIN
		CALL COMPUTE_WELFARE_DECOMPOSITION

	! Write experimental results in output.txt
		CALL WRITE_VARIABLES(0)

	! Save files
		CALL Write_Experimental_Results(.true.)

	! Print resutls 
		if (Opt_Tax_KW) then 
	 	OPEN (UNIT=77, FILE=trim(Result_Folder)//'stat_opt_tau_k.txt', STATUS='replace')
	 	else
	 	OPEN (UNIT=77, FILE=trim(Result_Folder)//'stat_opt_tau_w.txt', STATUS='replace')
	 	endif 
		WRITE(UNIT=77, FMT=*) tauK, tauW_at, psi, (GBAR_K+GBAR_W)/(GBAR_bench +SSC_Payments_bench ), & 
		      &  MeanWealth, QBAR, QBAR/MeanWealth,NBAR, &
		      &  YBAR, 100.0_DP*(Y_exp/Y_bench-1.0),MeanCons,100.0_DP*(MeanCons/MeanCons_bench-1.0), &
		      &  wage, R, &
		      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, &
		  	  & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Tot_Cap_Inc, Av_Util_Pop, Av_Util_NB, brentvaluet
		CLOSE (UNIT=77)

	else ! (Read_Results.eqv..true.)

	! Read results from file 
		CALL Write_Experimental_Results(.false.)
		CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
		K_mat  = K_Matrix(R,P)
		Pr_mat = Profit_Matrix(R,P)
		CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
		CALL ComputeLaborUnits(EBAR,wage)
		CALL GOVNT_BUDGET(.true.)
		CALL Compute_After_Tax_Income
		! CALL Write_Experimental_Results(.true.)

	! Aggregate variable in experimental economy
		GBAR_exp  = GBAR
		QBAR_exp  = QBAR 
		NBAR_exp  = NBAR  
		Y_exp 	  = YBAR
		Ebar_exp  = EBAR
		P_exp     = P
		R_exp	  = R
		wage_exp  = wage
		tauK_exp  = tauK
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

		ValueFunction_exp = ValueFunction
		Bq_Value_exp	  = Bq_Value
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime 

	! Compute moments
		CALL COMPUTE_STATS
	
	! Compute welfare gain between economies
		CALL COMPUTE_WELFARE_GAIN
		CALL COMPUTE_WELFARE_DECOMPOSITION

	! Write experimental results in output.txt
		CALL WRITE_VARIABLES(0)

	endif 
	
	

	! ! if (((theta.eq.1.50_dp)).and.(Threshold_Factor.eq.0.0_dp).and.(Simul_Switch)) then 
	! !  	print*,"	Optimal Tax Simulation"
	! ! 	CALL SIMULATION(solving_bench)
	! ! endif
	
	! CALL SIMULATION(0)

end Subroutine Solve_Opt_Tax


!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Opt_Threshold
	use parameters 
	use global
	use Opt_Tax_Parameters
	use Opt_Tax_Functions
	use Simulation_Module
	use programfunctions
	use GKK_Stats
	use Toolbox
	use omp_lib
	implicit none 
	real(DP) :: OPT_Threshold, OPT_psi_th, OPT_tauW_th, maxbrentvaluet_th 
	INTEGER  :: Threshold_ind
	integer  :: tau_grid_min, tau_grid_max, tau_grid_step
	logical  :: read_results, load_seed, only_grid
	character(100) :: folder_aux, OTW_Folder

	! Save base folder
		folder_aux = Result_Folder
		OTW_Folder = trim(folder_aux)//'Opt_Tax_W/'

	! Experiment economy
		solving_bench=0

	! Set flag for reading results or computing optimal taxes
		read_results = .false.
		load_seed    = .false.
		only_grid    = .true.


	if (load_seed.eqv..false.) then 
	! Load Optimal Tax Variables (for first start)
		Result_Folder = trim(folder_aux)//'Opt_Tax_W/'
		
		! Load variables
		CALL Write_Experimental_Results(.false.)
	endif

	! Set Results Folder
		Result_Folder = trim(folder_aux)//'Opt_Tax_W_Threshold_high/'
		folder_aux    = Result_Folder
		call system( 'mkdir -p ' // trim(Result_Folder) )



	if (read_results.eqv..false.) then 
 	!====================================================================================================
	PRINT*,''
	Print*,'--------------- SOLVING OPTIMAL Threshold -----------------'
	PRINT*,''

	! Set initial taxes for finding optimal ones
		tauK     = 0.0_DP
		tauW_at  = 0.0_DP
		Opt_TauK = 0.0_DP
		Opt_TauW = 0.0_DP
		Opt_psi  = 0.0_DP
 		maxbrentvaluet=-10000.0_DP

 		OPT_psi_th  = 0.0_dp 
 		OPT_tauW_th = 0.0_dp
 		maxbrentvaluet_th = -10000.0_DP

	! Set Grid for Optimal Tax
    	tau_grid_min  = 25
    	tau_grid_max  = 45
    	tau_grid_step = 1

	! Load results form file for re-starts of the code
	if (load_seed) then
		CALL Write_Experimental_Results(.false.)

		! ! Adjust psi if needed 
    	! psi = 1.0_dp-0.02811_dp

		! ! Optimal value of taxes and threshold if needed 
		! OPT_tauW      = 0.0395_dp
		! OPT_psi       = 0.875502136392213_dp
		! OPT_Threshold = 0.2_dp 
		! OPT_psi_th    = OPT_psi 
		! OPT_tauW_th   = OPT_tauW_th
	endif 
	
	OPEN(UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_w_threshold_2.txt', STATUS='replace')
	WRITE(UNIT=77, FMT=*) 'Threshold_Factor ', 'tauK ', 'tauW_at ', 'psi ', 'GBAR_K/Tax_Rev_bench ', &
		      & 'MeanWealth ','QBAR ','NBAR ','YBAR ','Y_Growth ', 'wage ', &
		      & 'Av_Util_NB ', 'CE2_NB ', 'CE2_Pop ', &
		      & 'Wealth_Output ', 'prct1_wealth ' , 'prct10_wealth ', 'Std_Log_Earnings ', 'mean_hours ', &
	      	  & 'GBAR ', 'GBAR_K ', 'GBAR_W ', 'GBAR_L ', 'GBAR_C ', 'Av_Util_Pop ', 'Av_Util_NB ', 'brentvaluet ','Threshold_Share'
	CLOSE(unit=77) 


	print*, '	Threshold Loop'	
	DO Threshold_ind = 120,160,20

		! Load Optimal Wealth Tax Results to start iteration fresh
		Result_Folder = OTW_Folder
		CALL Write_Experimental_Results(.false.)
		Result_Folder = folder_aux

		! Compute Threshold Factor for each iteration 
		Threshold_Factor = real(Threshold_ind,8)/100.0_dp
		print*, ' Threshold_Factor=',Threshold_Factor
    	! Set Y_a_threshold
			call Find_TauW_Threshold(DBN_bench,W_bench)  
			Y_a_threshold = Threshold_Factor*Ebar_bench !0.75_dp
			Wealth_factor = Y_a_threshold/W_bench

		! Initialize maxbrentvaluet
			maxbrentvaluet = -10000.0_DP

	
		PRINT*,''
		Print*,'--------------- OPTIMAL WEALTH TAXES - Threshold -----------------'
		PRINT*,''  	
	    DO tauindx = tau_grid_min,tau_grid_max,tau_grid_step
            tauw_at     = real(tauindx,8)/1000_DP
            brentvaluet = - EQ_WELFARE_GIVEN_TauW(tauW_at)

            ! Allocate variables
			CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)

            ! Aggregate variable in experimental economy
				GBAR_exp  = GBAR
				QBAR_exp  = QBAR 
				NBAR_exp  = NBAR  
				Y_exp 	  = YBAR
				Ebar_exp  = EBAR
				P_exp     = P
				R_exp	  = R
				wage_exp  = wage
				tauK_exp  = tauK
				tauPL_exp = tauPL
				psi_exp   = psi
				DBN_exp   = DBN1
				tauw_bt_exp = tauW_bt
				tauw_at_exp = tauW_at
				Y_a_threshold_exp = Y_a_threshold

				ValueFunction_exp = ValueFunction
				Bq_Value_exp	  = Bq_Value
				Cons_exp          = Cons           
				Hours_exp         = Hours
				Aprime_exp        = Aprime 

			! Compute moments
				CALL COMPUTE_STATS
				CALL GOVNT_BUDGET(.false.)
				CALL Compute_After_Tax_Income
				
			! Compute welfare gain between economies
				CALL COMPUTE_WELFARE_GAIN

			! Write experimental results in output.txt
				CALL WRITE_VARIABLES(0)

		    if (Av_Util_NB .gt. maxbrentvaluet) then
		        maxbrentvaluet = Av_Util_NB
		        OPT_tauK = tauK
				OPT_tauW = tauW_at
				OPT_psi  = psi
			endif

			! Print Results 
			print*,' ';print*,'------------------------------------------------------------------------------'
		    print '(A,F7.3,X,A,F7.4,X,A,F7.3,X,A,F7.3,X,A,F7.3,X,A,F12.6,X,A,F12.6)', &
		    	  & 'Threshold=', 100.0_dp*Threshold_Factor,'tauW=', 100.0_dp*tauW_at,'tauL',100.0_dp*(1.0_dp-psi),&
		    	  & 'Y/Y_bench=',100.0_dp*(YBAR/Y_bench-1), 'CE2=',Av_Util_NB, 'Brentvalue=',brentvaluet, & 
		    	  & 'Av. Util=', sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
	      	print*,'------------------------------------------------------------------------------';print*,' '

		      
		    OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_w_threshold_2.txt', STATUS='old', POSITION='append')
		    WRITE  (UNIT=77, FMT=*) Threshold_Factor, tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), &
		      &  MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0), &
		      &  wage, sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)),  &
			!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)), &
		      &100*( (sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)) /&
		      &sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))) &
		      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
		      &100*( (sum(ValueFunction(:,:,:,:,:,:)*DBN1(:,:,:,:,:,:))/sum(DBN1(:,:,:,:,:,:)) /&
		      &sum(ValueFunction_bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:,:))) &
		      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
		      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, &
	      	  & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_Pop, Av_Util_NB, brentvaluet, Threshold_Share, Tot_Cap_Inc
      	  	CLOSE (unit=77)
      	  	CALL Write_Experimental_Results(.true.)

   			! ! Deallocate variables
			! deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )

	    ENDDO

    	print*,' ';print*,'------------------------------------------------------------------------------'
		print'(A,F7.3)','	End of Optimal Tax Loop - Threshold',100.0_dp*Threshold_Factor
			tauK    = OPT_tauK
			tauW_at = OPT_tauW
			psi     = OPT_psi
		print '(A,F7.3,X,A,F7.3,X,A,F7.3)', &
			&"	Current Optimal Taxes:  tau_K=", 100.0_dp*tauK,"tau_W=",100.0_dp*tauW_at,"tau_L=", 100.0_dp*(1.0_dp-psi)
		print*,'------------------------------------------------------------------------------';print*,' '


		if (only_grid) then

		! Solve model with optimal taxes
			CALL FIND_DBN_EQ
			CALL GOVNT_BUDGET(.true.)
			CALL Compute_After_Tax_Income

		! Compute value function and store policy functions, value function and distribution in file
			CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
			CALL Firm_Value
			CALL Write_Experimental_Results(.true.)

		else 

		! Search for optimal tax
		print*,'	Optimal Tax Search'
		call Find_Opt_Tax(.false.,Opt_TauW,Opt_TauW-0.001_dp,Opt_TauW+0.001_dp)
		tauW_at = Opt_tauW
		OPT_psi  = psi

		print*,' ';print*,'------------------------------------------------------------------------------'
		print'(A,F7.3)','	End of Search for Optimal Tax - Threshold',100.0_dp*Threshold_Factor
		print '(A,F7.3,X,A,F7.3,X,A,F7.3)', &
			&"	Optimal Taxes:  tau_K=", 100.0_dp*tauK,"tau_W=",100.0_dp*tauW_at,"tau_L=", 100.0_dp*(1.0_dp-psi)
		print*,'------------------------------------------------------------------------------';print*,' '

		endif 


		! Aggregate variable in experimental economy
			GBAR_exp  = GBAR
			QBAR_exp  = QBAR 
			NBAR_exp  = NBAR  
			Y_exp 	  = YBAR
			Ebar_exp  = EBAR
			P_exp     = P
			R_exp	  = R
			wage_exp  = wage
			tauK_exp  = tauK
			tauPL_exp = tauPL
			psi_exp   = psi
			DBN_exp   = DBN1
			tauw_bt_exp = tauW_bt
			tauw_at_exp = tauW_at
			Y_a_threshold_exp = Y_a_threshold

			ValueFunction_exp = ValueFunction
			Bq_Value_exp	  = Bq_Value
			Cons_exp          = Cons           
			Hours_exp         = Hours
			Aprime_exp        = Aprime 

		! Compute moments
			CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
			CALL COMPUTE_STATS
			CALL GOVNT_BUDGET(.true.)
			CALL Compute_After_Tax_Income
			
		! Compute welfare gain between economies
			CALL COMPUTE_WELFARE_GAIN
			CALL COMPUTE_WELFARE_DECOMPOSITION

		! Write experimental results in output.txt
			CALL WRITE_VARIABLES(0)

		! Update optimal threshold
		brentvaluet = sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:));
	    if (Av_Util_NB .gt. maxbrentvaluet_th) then
	        maxbrentvaluet_th = Av_Util_NB
			OPT_tauW_th 	  = tauW_at
			OPT_psi_th  	  = psi
			OPT_Threshold 	  = Threshold_Factor
		endif

		! Print Results 
			print*,' ';print*,'------------------------------------------------------------------------------'
		    print '(A,F7.3,X,A,F7.4,X,A,F7.3,X,A,F7.3,X,A,F7.3,X,A,F12.6,X,A,F12.6)', &
		    	  & 'Threshold=', 100.0_dp*Threshold_Factor,'tauW=', 100.0_dp*tauW_at,'tauL',100.0_dp*(1.0_dp-psi),&
		    	  & 'Y/Y_bench=',100.0_dp*(YBAR/Y_bench-1), 'CE2=',Av_Util_NB, 'Brentvalue=',brentvaluet, & 
		    	  & 'Av. Util=', sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
	      	print*,'------------------------------------------------------------------------------';print*,' '


	    OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_w_threshold_2.txt', STATUS='old', POSITION='append') 
	    WRITE  (UNIT=77, FMT=*) Threshold_Factor, tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), &
	      &  MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0), &
	      &  wage, sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)),  &
		!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)), &
	      &100*( (sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)) /&
	      &sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))) &
	      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
	      &100*( (sum(ValueFunction(:,:,:,:,:,:)*DBN1(:,:,:,:,:,:))/sum(DBN1(:,:,:,:,:,:)) /&
	      &sum(ValueFunction_bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:,:))) &
	      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
	      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, &
      	  & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_Pop, Av_Util_NB, brentvaluet, Threshold_Share
  	  	CLOSE (unit=77)

  	  	! Save results 
  	  	CALL Write_Experimental_Results(.true.)

    ENDDO ! Threshold loop

    	print*,' ';print*,'------------------------------------------------------------------------------'
		print*,'	End of Optimal Threshold Loop'
			tauW_at = OPT_tauW_th
			psi 	= OPT_psi_th
			Threshold_Factor = OPT_Threshold
		print '(A,F7.3,X,A,F7.3,X,A,F7.3)', &
			&"	Current Optimal Taxes:  Threshold=", 100.0_dp*Threshold_Factor,"tau_W=",100.0_dp*tauW_at,"tau_L=", 100.0_dp*(1.0_dp-psi)
		print*,'------------------------------------------------------------------------------';print*,' '

	
	! Set Y_a_threshold
		call Find_TauW_Threshold(DBN_bench,W_bench)  
		Y_a_threshold = Threshold_Factor*Ebar_bench !0.75_dp
		Wealth_factor = Y_a_threshold/W_bench

	! Solve model with optimal taxes
		CALL FIND_DBN_EQ
		CALL GOVNT_BUDGET(.true.)
		CALL Compute_After_Tax_Income

	! Compute value function and store policy functions, value function and distribution in file
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
		CALL Firm_Value
		CALL Write_Experimental_Results(.true.)
	
	! Aggregate variable in experimental economy
		GBAR_exp  = GBAR
		QBAR_exp  = QBAR 
		NBAR_exp  = NBAR  
		Y_exp 	  = YBAR
		Ebar_exp  = EBAR
		P_exp     = P
		R_exp	  = R
		wage_exp  = wage
		tauK_exp  = tauK
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

		ValueFunction_exp = ValueFunction
		Bq_Value_exp	  = Bq_Value
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime 

	! Compute moments
		CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
		CALL COMPUTE_STATS
	
	! Compute welfare gain between economies
		CALL COMPUTE_WELFARE_GAIN
		CALL COMPUTE_WELFARE_DECOMPOSITION

	! Write experimental results in output.txt
		CALL WRITE_VARIABLES(0)

	! CALL SIMULATION(0)

	OPEN  (UNIT=77, FILE=trim(Result_Folder)//'stat_opt_tau_w_threshold.txt', STATUS='replace')
    WRITE  (UNIT=77, FMT=*) Threshold_Factor, tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), &
      &  MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0), &
      &  wage, sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)),  &
	!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)), &
      &100*( (sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)) /&
      &sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))) &
      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
      &100*( (sum(ValueFunction(:,:,:,:,:,:)*DBN1(:,:,:,:,:,:))/sum(DBN1(:,:,:,:,:,:)) /&
      &sum(ValueFunction_bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:,:))) &
      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, &
  	  & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_Pop, Av_Util_NB, brentvaluet, Threshold_Share, Tot_Cap_Inc
  	CLOSE (unit=77)



	else ! Load results from file 

		! Set threshold manually 
			Threshold_Factor = 0.000_dp

		! Read results from file 
			CALL Write_Experimental_Results(.false.)
			CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
			K_mat  = K_Matrix(R,P)
			Pr_mat = Profit_Matrix(R,P)
			CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
			CALL ComputeLaborUnits(EBAR,wage)
			CALL GOVNT_BUDGET(.true.)
			CALL Compute_After_Tax_Income
			! CALL Write_Experimental_Results(.true.)

		! Aggregate variable in experimental economy
			GBAR_exp  = GBAR
			QBAR_exp  = QBAR 
			NBAR_exp  = NBAR  
			Y_exp 	  = YBAR
			Ebar_exp  = EBAR
			P_exp     = P
			R_exp	  = R
			wage_exp  = wage
			tauK_exp  = tauK
			tauPL_exp = tauPL
			psi_exp   = psi
			DBN_exp   = DBN1
			tauw_bt_exp = tauW_bt
			tauw_at_exp = tauW_at
			Y_a_threshold_exp = Y_a_threshold

			ValueFunction_exp = ValueFunction
			Bq_Value_exp	  = Bq_Value
			Cons_exp          = Cons           
			Hours_exp         = Hours
			Aprime_exp        = Aprime 

		! Compute moments
			CALL COMPUTE_STATS
		
		! Compute welfare gain between economies
			CALL COMPUTE_WELFARE_GAIN
			CALL COMPUTE_WELFARE_DECOMPOSITION

		! Write experimental results in output.txt
			CALL WRITE_VARIABLES(0)

	endif 
	

end Subroutine Solve_Opt_Threshold



!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Transition_Tax_Reform(budget_balance)
	use parameters
	use global 
	use programfunctions
	use GKK_Stats
	use Simulation_Module
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: budget_balance
	real(dp)            :: tauW_at_0, tauW_bt_0
	character(100)      :: folder_aux
	

	! Set step for increments
	tauWinc_bt=0.000_DP
	tauWinc_at=0.001_DP

	! Load Benchmark Variables
	call Solve_Benchmark(.false.,.false.)

	! Load Tax Reform Variables
	folder_aux = Result_Folder
	Result_Folder = trim(folder_aux)//'Tax_Reform/'
	call Solve_Experiment(.false.,.false.)

	! Base level for taxes (current experimental value)
	tauw_bt_0 = tauW_bt_exp 
	tauw_at_0 = tauW_at_exp 

	Debt_Absorption = 1.0_dp
	Use_Transition_Seed = .false.

	if (budget_balance) then 

		! Set Results Folder
		Result_Folder = trim(Result_Folder)//'Transition_Tax_Reform_BB_cfm/'
		call system( 'mkdir -p ' // trim(Result_Folder) )

	if (.true.) then 

		! Find the Distribution and Policy Functions Along Transition Path
		! This is done for the tax reform steady state
		! call Find_DBN_Transition 
			CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
			GBAR_exp = GBAR_bench - 1 ; Debt_tr = 0 ; 

		! Find Taxes that balance the budget 
		print*,' '
		print*,'---------------------------------------------------'
		print*,' 	Balancing the Budget'
		print*,'---------------------------------------------------'
			! Solve for the model increasing wealth taxes until revenue is enough to finance G_benchamark
			tauWindx = 4.0_DP
			Debt_tr(T+1)  = 1.0_DP
			DO WHILE (GBAR_exp .lt. (GBAR_bench+R_exp*Debt_tr(T+1)))
				! Set old G and new value of tauW
				GBAR_exp_old = GBAR_exp
				tauW_bt = tauw_bt_0 + tauWindx * tauWinc_bt
				tauW_at = tauw_at_0 + tauWindx * tauWinc_at
				print*, 'Bracketing Iteration',tauWindx,'tauW_bt=',tauW_bt*100,"tauW_at=",tauW_at*100
				! Solve for New Steady State
				deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )
				CALL FIND_DBN_EQ
					GBAR_exp  = GBAR
					QBAR_exp  = QBAR 
					NBAR_exp  = NBAR  
					Y_exp 	  = YBAR
					Ebar_exp  = EBAR
					P_exp     = P
					R_exp	  = R
					wage_exp  = wage
					tauK_exp  = tauK
					tauPL_exp = tauPL
					psi_exp   = psi
					DBN_exp   = DBN1
					tauw_bt_exp = tauW_bt
					tauw_at_exp = tauW_at
					Y_a_threshold_exp = Y_a_threshold
					Cons_exp          = Cons           
					Hours_exp         = Hours
					Aprime_exp        = Aprime
				! Find the Distribution and Policy Functions Along Transition Path
				call Find_DBN_Transition 
				! Get new G
				GBAR_exp = GBAR_tr(T+1) 
				! Iteratioins  
				tauWindx = tauWindx + 1.0_DP  
				print*,' ' 
				print*,' ' 
				write(*,*) "Bracketing GBAR: tauW_bt=", tauW_bt*100, "And tauW_at=", tauW_at*100
				print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
				print*,'GBAR_exp =', GBAR_exp,'GBAR_bench+R*Debt=',GBAR_bench+R_exp*Debt_tr(T+1),'Debt',Debt_tr(T+1)
			ENDDO

			! Set tauW as weighted average of point in  the grid to balance budget more precisely
				tauW_up_bt  = tauW_bt
				tauW_low_bt = tauW_bt  -  tauWinc_bt
				tauW_bt     = tauW_low_bt + tauWinc_bt * 0.5_dp ! GBAR_bench+R_exp*Debt_tr(T+1) - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
				tauW_up_at  = tauW_at
				tauW_low_at = tauW_at  -  tauWinc_at  
				tauW_at     = tauW_low_at + tauWinc_at * 0.5_dp ! (GBAR_bench+R_exp*Debt_tr(T+1) - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
				print*,''
				print*,'GBAR bracketed by taxes:'
				print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_bt=', tauW_bt*100, "%", '% tauW_up_bt=', tauW_up_bt*100
				print*,'tauW_low_at =', tauW_low_at*100, '% tauW_at=', tauW_at*100, "%", '% tauW_up_at=', tauW_up_at*100
				print*,''

			! Solve (again) experimental economy
				! Solve for New Steady State
				deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )
				CALL FIND_DBN_EQ
					GBAR_exp  = GBAR
					QBAR_exp  = QBAR 
					NBAR_exp  = NBAR  
					Y_exp 	  = YBAR
					Ebar_exp  = EBAR
					P_exp     = P
					R_exp	  = R
					wage_exp  = wage
					tauK_exp  = tauK
					tauPL_exp = tauPL
					psi_exp   = psi
					DBN_exp   = DBN1
					tauw_bt_exp = tauW_bt
					tauw_at_exp = tauW_at
					Y_a_threshold_exp = Y_a_threshold
					Cons_exp          = Cons           
					Hours_exp         = Hours
					Aprime_exp        = Aprime
				! Find the Distribution and Policy Functions Along Transition Path
				call Find_DBN_Transition 
				! Get new G
				GBAR_exp = GBAR_tr(T+1) 

			! Find tauW that exactly balances the budget (up to precisioin 0.1) using bisection
				print*,"Gbar at midpoint of bracket"
				print*,'GBAR_exp =', GBAR_exp,'GBAR_bench+R*Debt=',GBAR_bench+R_exp*Debt_tr(T+1)
				print*,''
				print*,'Bisection for TauW:'
				DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/(GBAR_bench+R_exp*Debt_tr(T+1)))) .gt. 0.01 ) ! as long as the difference is greater than 0.1% continue
				    if (GBAR_exp .gt. GBAR_bench+R_exp*Debt_tr(T+1) ) then
				        tauW_up_bt  = tauW_bt 
				        tauW_up_at  = tauW_at 
				    else
				        tauW_low_bt = tauW_bt
				        tauW_low_at = tauW_at
				    endif
				    tauW_bt = (tauW_low_bt + tauW_up_bt)/2.0_DP
				    tauW_at = (tauW_low_at + tauW_up_at)/2.0_DP
					! Solve for New Steady State
					deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )
					CALL FIND_DBN_EQ
						GBAR_exp  = GBAR
						QBAR_exp  = QBAR 
						NBAR_exp  = NBAR  
						Y_exp 	  = YBAR
						Ebar_exp  = EBAR
						P_exp     = P
						R_exp	  = R
						wage_exp  = wage
						tauK_exp  = tauK
						tauPL_exp = tauPL
						psi_exp   = psi
						DBN_exp   = DBN1
						tauw_bt_exp = tauW_bt
						tauw_at_exp = tauW_at
						Y_a_threshold_exp = Y_a_threshold
						Cons_exp          = Cons           
						Hours_exp         = Hours
						Aprime_exp        = Aprime
					! Find the Distribution and Policy Functions Along Transition Path
					call Find_DBN_Transition 
					! Get new G
					GBAR_exp = GBAR_tr(T+1) 
					! Print Results 
				    print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
					print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
					print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
					print*,'GBAR_exp =', GBAR_exp,'GBAR_bench+R*Debt=',GBAR_bench+R_exp*Debt_tr(T+1),'Debt',Debt_tr(T+1)
				ENDDO

	else 

		print*,' '
		print*,'---------------------------------------------------'
		print*,' 	Computing Steady State at desired tax level'
		print*,'---------------------------------------------------'
		! Read Tax
			! OPEN (UNIT=4,  FILE=trim(Result_Folder)//'tauW_at_tr', STATUS='old', ACTION='read')
			! READ (UNIT=4,  FMT=*) tauW_at
			! CLOSE(unit=4)
			! R=   1.9480499900981853E-002 
			! P=  0.13847606093758086 
			! print*,' '
			! print*,'	Wealth taxes =',tauW_at*100,'%'

		! Solve for New Steady State
			CALL FIND_DBN_EQ
			CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
			CALL Write_Experimental_Results(.true.)
			CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
			K_mat  = K_Matrix(R,P)
			Pr_mat = Profit_Matrix(R,P)
			CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
			CALL ComputeLaborUnits(EBAR,wage)
			CALL GOVNT_BUDGET(.true.)

			! Aggregate variable in experimental economy
				GBAR_exp  = GBAR
				QBAR_exp  = QBAR 
				NBAR_exp  = NBAR  
				Y_exp 	  = YBAR
				Ebar_exp  = EBAR
				P_exp     = P
				R_exp	  = R
				wage_exp  = wage
				tauK_exp  = tauK
				tauPL_exp = tauPL
				psi_exp   = psi
				DBN_exp   = DBN1
				tauw_bt_exp = tauW_bt
				tauw_at_exp = tauW_at
				Y_a_threshold_exp = Y_a_threshold

				ValueFunction_exp = ValueFunction
				Bq_Value_exp	  = Bq_Value
				Cons_exp          = Cons           
				Hours_exp         = Hours
				Aprime_exp        = Aprime
				V_Pr_exp          = V_Pr 
				V_Pr_nb_exp  	  = V_Pr_nb

			! Compute moments
			CALL COMPUTE_STATS
			
			! Compute welfare and output gain between economies
				CALL COMPUTE_WELFARE_GAIN

				print*,'---------------------------'
				print*,'SS Output Gain Prct=', 100.0_DP*(Y_exp/Y_bench-1.0) 
				print*,'---------------------------'

			! Write experimental results in output.txt
				CALL WRITE_VARIABLES(0)

			! Deallocate variables
				deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )


		! Find the Distribution and Policy Functions Along Transition Path
			call Find_DBN_Transition 
	endif

	else ! If budget isn't balanced: Transition between steady states

		! Set Results Folder
		Result_Folder = trim(folder_aux)//'Transition_Tax_Reform/'
		call system( 'mkdir -p ' // trim(Result_Folder) )

		! Find the Distribution and Policy Functions Along Transition Path
		! This is done for the tax reform steady state
		call Find_DBN_Transition 

	endif 




	! Compute Value Functions for Cohorts Alive at Time of Policy Change
	call COMPUTE_VALUE_FUNCTION_TRANSITION


	! Compute Welfare Gain
	call COMPUTE_WELFARE_GAIN_TRANSITION


End Subroutine Solve_Transition_Tax_Reform


!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Transition_Opt_Taxes(Opt_Tax_KW,budget_balance,balance_tau_L)
	use parameters
	use global 
	use programfunctions
	use GKK_Stats
	use Simulation_Module
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: budget_balance, Opt_Tax_KW, balance_tau_L
	real(dp) :: psi_0, tauK_0, tauW_0
	integer  :: Debt_Absorption_iter 
	character(100) :: folder_aux
	logical  :: read_results

	! Set flag for reading results or computing optimal taxes
		read_results = .false.

	! Save base folder
		folder_aux = Result_Folder

	! Load Benchmark Variables
		call Solve_Benchmark(.false.,.false.)
		! Change flag 
		solving_bench=0

	! Load Optimal Tax Variables
		! Change to optimal tax folder 
		if (Opt_Tax_KW) then 
		Result_Folder = trim(folder_aux)//'Opt_Tax_K/'
		call system( 'mkdir -p ' // trim(Result_Folder) )
		else 
		Result_Folder = trim(folder_aux)//'Opt_Tax_W/'
		call system( 'mkdir -p ' // trim(Result_Folder) )
		endif 
		
		! Load variables
		CALL Write_Experimental_Results(.false.)
		
		! Compute auxiliary variables
		CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
		K_mat  = K_Matrix(R,P)
		Pr_mat = Profit_Matrix(R,P)
		CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
		CALL ComputeLaborUnits(EBAR,wage)
		CALL GOVNT_BUDGET(.false.)
		deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )

		! Aggregate variable in experimental economy
		GBAR_exp  = GBAR
		QBAR_exp  = QBAR 
		NBAR_exp  = NBAR  
		Y_exp 	  = YBAR
		Ebar_exp  = EBAR
		P_exp     = P
		R_exp	  = R
		wage_exp  = wage
		tauK_exp  = tauK
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

		ValueFunction_exp = ValueFunction
		Bq_Value_exp	  = Bq_Value
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime 

	! Set reference value for psi, tau_K and tau_W
		psi_0  = 1.0_dp - 0.35_dp
			! OTW with tauL set to 1.0_dp - 0.16989829280240965_dp
			! OTK with tauL set to 1.0_dp-0.3740_dp
			! If not using tauL set to psi
		tauK_0 = tauK 
		tauW_0 = tauW_at

	Use_Transition_Seed = .true.

		
	if (budget_balance) then 

		! Set Results Folder
			if (balance_tau_L) then 
				if (Opt_Tax_KW) then 
				Result_Folder = trim(folder_aux)//'Transition_OTK_tL/'
				else 
				Result_Folder = trim(folder_aux)//'Transition_OTW_tL/'
				endif 
			elseif (Opt_Tax_KW) then 
				Result_Folder = trim(folder_aux)//'Transition_OTK_tK/'
			else
				Result_Folder = trim(folder_aux)//'Transition_OTW_tW/'
			endif 
			call system( 'mkdir -p ' // trim(Result_Folder) )

		! Find the Distribution and Policy Functions Along Transition Path
		! This is done for the tax reform steady state
		! call Find_DBN_Transition 
			! If previous line is commented you need this:
			CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)

		! Find Taxes that balance the budget 
		if (balance_tau_L) then
		print*,' '
		print*,'---------------------------------------------------'
		print*,' 	Balancing the Budget with Labor Taxes'
		print*,'---------------------------------------------------'
		else 
		print*,' '
		print*,'---------------------------------------------------'
		print*,' 	Balancing the Budget with Capital/Wealth Taxes'
		print*,'---------------------------------------------------'
		endif 

		if (read_results.eqv..false.) then 
		! Solve for the optimal tax for iterative loops of Debt_Absorption

		DO Debt_Absorption_iter=6,10,2

			! Set Debt_Absorption
			Debt_Absorption = real(Debt_Absorption_iter,8)/10.0_dp

			! Solve for the model increasing wealth taxes until revenue is enough to finance G_benchamark
			BB_tax_ind = 0.0_DP ! Originally 1.0_DP
			BB_tax_chg = 0.002_DP ! Originally 0.0002_DP
			Debt_tr  = 1.0_DP
			DO WHILE (GBAR_exp .lt. (GBAR_bench+R_exp*Debt_tr(T+1)))
				! Set old G and new value of tauW
				GBAR_exp_old = GBAR_exp
				print*,' '
				if (balance_tau_L) then 
				psi = psi_0 - BB_tax_ind*BB_tax_chg ! Decreasing psi increases labor taxes
				print*, 'Bracketing Iteration',BB_tax_ind,'Increment',BB_tax_chg,'tau_L=',(1.0_dp-psi)*100
				elseif (Opt_Tax_KW) then 
				tauK = tauK_0 + BB_tax_ind * BB_tax_chg 
				print*, 'Bracketing Iteration',BB_tax_ind,'Increment',BB_tax_chg,"tauK=",tauK*100
				else
				tauW_at = tauW_0 + BB_tax_ind * BB_tax_chg 
				print*, 'Bracketing Iteration',BB_tax_ind,'Increment',BB_tax_chg,"tauW_at=",tauW_at*100
				endif 
				print*,' '
				

				! Find the Distribution and Policy Functions Along Transition Path
				call Find_DBN_Transition 
				! Get new G
				GBAR_exp = GBAR_tr(T+1) 
				! Iteratioins  
				BB_tax_ind = BB_tax_ind + 1.0_DP  
				print*,' ' 
				print*,'Bracketing GBAR: tau_L=', (1.0_dp-psi)*100,'tauK=',100*tauK,'tauW=',100*tauW_at
				print*,'GBAR_exp =', GBAR_exp,'GBAR_bench+R*Debt=',GBAR_bench+R_exp*Debt_tr(T+1),'Debt',Debt_tr(T+1)
				print*,'	Debt Absortion',Debt_Absorption_iter,Debt_Absorption
				print*,' ' 
			ENDDO

				print*,''
				print*,'GBAR bracketed by taxes:'
			if (balance_tau_L) then 
			! Set the upper bound of psi as the second to last iteration
				BB_tax_up  = psi + BB_tax_chg
			! Set the lower bound of psi as the last iteration
				BB_tax_low = psi
				print*,'tauL_low =',(1.0-BB_tax_up)*100,'tauL=',(1.0-psi-0.5_dp*BB_tax_chg)*100,'tauL_up=',(1.0-BB_tax_low)*100
			elseif (Opt_Tax_KW) then 
			! Set the upper bound of tau_K as the last iteration
				BB_tax_up  = tauK
			! Set the lower bound of tau_K as the second to last iteration
				BB_tax_low = tauK - BB_tax_chg
				print*,'tauK_low =',BB_tax_low*100,'tauK=',(tauK-0.5_dp*BB_tax_chg)*100,'tauK_up=',BB_tax_up*100
			else
			! Set the upper bound of tau_W as the last iteration
				BB_tax_up  = tauW_at 
			! Set the lower bound of tau_W as the second to last iteration
				BB_tax_low = tauW_at - BB_tax_chg
				print*,'tauW_low =',BB_tax_low*100,'tauW=',(tauW_at-0.5_dp*BB_tax_chg)*100,'tauW_up=',BB_tax_up*100
			endif 
				print*,''

			! Find psi that exactly balances the budget (up to precisioin 0.1%) using bisection
				print*,'Bisection for Taxes:'
				DO WHILE ((  abs(100.0_DP*(1.0_DP-GBAR_exp/(GBAR_bench+R_exp*Debt_tr(T+1)))) .gt. 0.05 )&
						&.and.(abs(BB_tax_up-BB_tax_low).gt.1.5E-05_DP)) ! as long as the difference is greater than 0.1% continue
			    	
					if (balance_tau_L) then 
					    if (GBAR_exp .gt. GBAR_bench+R_exp*Debt_tr(T+1) ) then
					        BB_tax_low  = psi ! If there is a surplus don't decrease psi (increase tau_L). Set a floor.
					    else
					        BB_tax_up   = psi ! If there is a deficit don't increase psi (decrease tau_L). Set a ceiling.
					    endif
    				    ! Set new tax
						    psi = (BB_tax_low + BB_tax_up)/2.0_DP
					elseif (Opt_Tax_KW) then 
					    if (GBAR_exp .gt. GBAR_bench+R_exp*Debt_tr(T+1) ) then
					        BB_tax_up   = tauK ! If there is a surplus don't increase tau_K. Set a ceiling.
					    else
					        BB_tax_low  = tauK ! If there is a deficit don't decrease tau_K. Set a floor.
					    endif
    				    ! Set new tax
						    tauK = (BB_tax_low + BB_tax_up)/2.0_DP
					else
					    if (GBAR_exp .gt. GBAR_bench+R_exp*Debt_tr(T+1) ) then
					        BB_tax_up   = tauW_at ! If there is a surplus don't increase tau_W. Set a ceiling.
					    else
					        BB_tax_low  = tauW_at ! If there is a deficit don't decrease tau_W. Set a floor.
					    endif
    				    ! Set new tax
						    tauW_at = (BB_tax_low + BB_tax_up)/2.0_DP
					endif 


					! Find the Distribution and Policy Functions Along Transition Path
					call Find_DBN_Transition 
					! Get new G
					GBAR_exp = GBAR_tr(T+1) 
					! Print Results 
					print*,' '
					if (balance_tau_L) then 
					print*,'tax_L_low =', (1.0_dp-BB_tax_up)*100, 'tau_L_up=', (1.0_dp-BB_tax_low)*100, 'tau_L=', (1.0_dp-psi)*100
					elseif (Opt_Tax_KW) then 
					print*,'tax_K_low =', BB_tax_low*100, 'tau_W_up=', BB_tax_up*100, 'tau_K=', tauK*100
					else
					print*,'tax_W_low =', BB_tax_low*100, 'tau_W_up=', BB_tax_up*100, 'tau_W=', tauW_at*100
					endif 
					print*,'GBAR_exp =', GBAR_exp,'GBAR_bench+R*Debt=',GBAR_bench+R_exp*Debt_tr(T+1),'Debt',Debt_tr(T+1)
					print*,'	Debt Absortion',Debt_Absorption_iter,Debt_Absorption
				ENDDO

			! Set reference value of taxes to current solution
				psi_0  = psi
				tauK_0 = tauK 
				tauW_0 = tauW_at


		ENDDO

		else 
		! Read results from main folder 

			Use_Transition_Seed = .true.
			Debt_Absorption     = 1.0_dp

			! Read taxes 
				print*, ' ' ; print*, 'Loading taxes from file'
				OPEN (UNIT=1,  FILE=trim(Result_Folder)//'tauC_tr'   , STATUS='old', ACTION='read')
				OPEN (UNIT=2,  FILE=trim(Result_Folder)//'tauK_tr'	 , STATUS='old', ACTION='read')
				OPEN (UNIT=3,  FILE=trim(Result_Folder)//'tauL_tr'	 , STATUS='old', ACTION='read')
				OPEN (UNIT=4,  FILE=trim(Result_Folder)//'tauW_at_tr', STATUS='old', ACTION='read')
				READ (UNIT=1,  FMT=*) tauC
				READ (UNIT=2,  FMT=*) tauK
				READ (UNIT=3,  FMT=*) psi
				READ (UNIT=4,  FMT=*) tauW_at
				CLOSE (unit=1); CLOSE (unit=2); CLOSE (unit=3); CLOSE (unit=4);
				print*, '	Reading completed'; print*, ' '

			! Adjust psi (file saves tauL = 1 - psi)
				psi = 1-psi

			! Find the Distribution and Policy Functions Along Transition Path
				call Find_DBN_Transition 


		endif

	else

		! Set Results Folder for no Budget Balance
			if (Opt_Tax_KW) then 
				Result_Folder = trim(folder_aux)//'Transition_OTK_no_BB/'
			else 
				Result_Folder = trim(folder_aux)//'Transition_OTW_no_BB/'
			endif 
			call system( 'mkdir -p ' // trim(Result_Folder) )

		! Set absorption to zero
			Debt_Absorption = 0.0_dp 

		! Find the Distribution and Policy Functions Along Transition Path
			call Find_DBN_Transition 

	endif


	! Compute Value Functions for Cohorts Alive at Time of Policy Change
		call COMPUTE_VALUE_FUNCTION_TRANSITION

	! Compute Welfare Gain
		call COMPUTE_WELFARE_GAIN_TRANSITION


End Subroutine Solve_Transition_Opt_Taxes

