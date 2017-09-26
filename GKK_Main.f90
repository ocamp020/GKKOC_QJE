

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
	use Simulation_Module
	use Toolbox
	use omp_lib

	IMPLICIT NONE
	! Variables to measure running time
		REAL(DP) :: start_time, finish_time
	! Compute benchmark or load results
		logical  :: compute_bench, compute_exp, Opt_Tax, Opt_Tax_KW, Tax_Reform, Simul_Switch, Calibration_Switch
		logical  :: Opt_Threshold, Opt_Tau_C, Opt_Tau_CX, Opt_Tax_K_and_W, Tax_Reform_KW
		logical  :: compute_exp_pf, Fixed_PF, Fixed_PF_interp, Fixed_PF_prices
		logical  :: compute_exp_prices, Fixed_W, Fixed_P, Fixed_R 
	! Auxiliary variable for writing file
		character(4)   :: string_theta
		character(100) :: folder_aux

	! Allocate Variables
	call Allocate_Variables

	! Capital Market
		theta_folder = 1.50_dp
		do zi=1,nz
		theta(zi)    = 1.00_dp+(2.50_dp-1.00_dp)/(nz-1)*(real(zi,8)-1.0_dp)
		enddo
	! Threshold 
		Threshold_Factor = 0.00_dp 

	! Switch for solving benchmark or just reading resutls
		Calibration_Switch = .false.
		! If compute_bench==.true. then just read resutls
		! If compute_bench==.false. then solve for benchmark and store results
		Tax_Reform    = .false.
			compute_bench = .false.
			compute_exp   = .false.
			compute_exp_pf= .false.
				Fixed_PF        = .false.
				Fixed_PF_interp = .false.
				Fixed_PF_prices = .true.
			compute_exp_prices    = .false.
				Fixed_W = .true. 
				Fixed_P = .true.
				Fixed_R = .true.
		Opt_Tax       = .true.
			Opt_Tax_KW    = .false. ! true=tau_K false=tau_W
		Opt_Tax_K_and_W = .false.
		Tax_Reform_KW   = .false.
		Opt_Threshold = .false.
		Opt_Tau_C = .false.
		Opt_Tau_CX = .false.
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
		
		beta   	= 0.9475_dp! 0.95_dp ! params(1) !
		mu_z   	= params(2) ! this is just shifting the z grids. it is zero now.
		rho_z  	= 0.1_dp ! params(3) 
		sigma_z_eps      =  0.072_dp !0.115_dp ! params(4) ! 0.01_dp ! ! 
		sigma_lambda_eps = 0.305_dp ! params(5)
		gamma  	=  0.46_dp !  0.471_dp ! params(6) ! 
		Params =[beta, mu_z, rho_z, sigma_z_eps, sigma_lambda_eps, gamma] 
		
		sigma  	= 4.0_dp
		phi    	= (1.0_dp-gamma)/gamma

		x_hi	= 5.00_dp
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

	! Resutls Folder
		write(Result_Folder,'(f4.2)') Threshold_Factor
		write(string_theta,'(f4.2)')  theta_folder

		if ((Progressive_Tax_Switch.eqv..false.).and.(NSU_Switch.eqv..true.)) then 
			Result_Folder = './NSU_ZS_LT_Results/Theta_'//trim(string_theta)//'/Factor_'//trim(Result_Folder)//'/'
		else if ((Progressive_Tax_Switch.eqv..true.).and.(NSU_Switch.eqv..true.)) then 
			Result_Folder = './NSU_ZS_PT_Results/Theta_'//trim(string_theta)//'/Factor_'//trim(Result_Folder)//'/'
		else if ((Progressive_Tax_Switch.eqv..false.).and.(NSU_Switch.eqv..false.)) then 
			Result_Folder = './SU_ZS_LT_Results/Theta_'//trim(string_theta)//'/Factor_'//trim(Result_Folder)//'/'
		else if ((Progressive_Tax_Switch.eqv..true.).and.(NSU_Switch.eqv..false.)) then 
			Result_Folder = './SU_ZS_PT_Results/Theta_'//trim(string_theta)//'/Factor_'//trim(Result_Folder)//'/'
		end if

		Result_Folder = trim(Result_Folder)//'Model_1.2_bv/' 

		! call execute_command_line( 'mkdir -p ' // trim(Result_Folder) )
		call system( 'mkdir -p ' // trim(Result_Folder) )
		print*, "Results are stored in directory: ", Result_Folder
		print*,'na=',na,'update_period=',update_period
		print*, "NSU_Switch=",NSU_Switch,'sigma=',sigma,'gamma=',gamma,'phi',phi
		print*,'Labor Taxes: tauPl=',tauPl,'psi',psi
		print*, 'Borrowing Constraint: Theta=',theta
		print*, 'm tauchen for zgrid is ',mtauchen_z,'nz=',nz, 'amax=',amax, 'totpop=', totpop
		print*, 'x_hi', x_hi, 'a_x', a_x, 'b_x', b_x, 'sigmaz', sigma_z_eps


	! Set parameters to be used in all simulations economy
		OPEN(UNIT=3, FILE=trim(Result_Folder)//'params', STATUS='replace')
		WRITE(unit=3, FMT=*) params
		CLOSE(unit=3)


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
					Result_Folder = trim(Result_Folder)//'Exp_Policy_Functions_Interp/'
					call system( 'mkdir -p ' // trim(Result_Folder) )
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
				call Solve_Experiment(compute_exp,Simul_Switch)
			endif 

			! Result_Folder = trim(Result_Folder)//'Tau_C_Experiment/'
			! call system( 'mkdir -p ' // trim(Result_Folder) )
			! call Solve_Experiment_tauC(compute_exp,Simul_Switch)

			compute_bench = .false.
		endif 

		! Optimal Tax
		if (Opt_Tax) then 
			call Solve_Benchmark(compute_bench,Simul_Switch)

			folder_aux = Result_Folder
			if (Opt_Tax_KW) then 
				Result_Folder = trim(folder_aux)//'Opt_Tax_K/'
			else 
				Result_Folder = trim(folder_aux)//'Opt_Tax_W/'
			endif
			call system( 'mkdir -p ' // trim(Result_Folder) )

			
			! call Solve_Opt_Tax(Opt_Tax_KW,Simul_Switch)
			
			CALL Write_Experimental_Results(.false.)
			CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
			K_mat  = K_Matrix(R,P)
			Pr_mat = Profit_Matrix(R,P)
			CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
			CALL ComputeLaborUnits(EBAR,wage)
			CALL GOVNT_BUDGET

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
				Cons_exp          = Cons           
				Hours_exp         = Hours
				Aprime_exp        = Aprime
				V_Pr_exp          = V_Pr 
				V_Pr_nb_exp  	  = V_Pr_nb

			! Compute moments
			CALL COMPUTE_STATS
			
			! Compute welfare gain between economies
			CALL COMPUTE_WELFARE_GAIN


		endif 

		if (Opt_Threshold) then 
			call Solve_Benchmark(compute_bench,Simul_Switch)

			folder_aux = Result_Folder
			Result_Folder = trim(folder_aux)//'Opt_Tax_W_Threshold/'
			call system( 'mkdir -p ' // trim(Result_Folder) )

			
			call Solve_Opt_Threshold
		endif 

		if (Opt_Tau_C) then 
			call Solve_Benchmark(compute_bench,Simul_Switch)

			folder_aux = Result_Folder
			if (Opt_Tax_KW) then 
				Result_Folder = trim(folder_aux)//'Opt_Tax_K_Tau_C_aux/'
			else 
				Result_Folder = trim(folder_aux)//'Opt_Tax_W_Tau_C/'
			endif
			call system( 'mkdir -p ' // trim(Result_Folder) )

			
			call Solve_Opt_Tau_C(Opt_Tax_KW)
			
			! print*, Result_Folder
			! CALL Write_Experimental_Results(.false.)
			! CALL SIMULATION(0)

		endif
		

		if (Opt_Tau_CX) then 
			call Solve_Benchmark(compute_bench,Simul_Switch)

			folder_aux = Result_Folder
			if (Opt_Tax_KW) then 
				Result_Folder = trim(folder_aux)//'Opt_Tax_K_Tau_C_No_Labor_Tax/'
			else 
				Result_Folder = trim(folder_aux)//'Opt_Tax_W_Tau_C_No_Labor_Tax/'
			endif
			call system( 'mkdir -p ' // trim(Result_Folder) )

			
			call Solve_Opt_Tau_CX(Opt_Tax_KW)
			
			! print*, Result_Folder
			! CALL Write_Experimental_Results(.false.)
			! CALL SIMULATION(0)

		endif

		if (Opt_Tax_K_and_W) then 
			call Solve_Benchmark(compute_bench,Simul_Switch)

			folder_aux = Result_Folder
			Result_Folder = trim(folder_aux)//'Opt_Tax_KW/'
			call system( 'mkdir -p ' // trim(Result_Folder) )

			
			call Solve_Opt_Tax_K_and_W(Simul_Switch)
			
			! print*, Result_Folder
			! CALL Write_Experimental_Results(.false.)
			! CALL SIMULATION(0)


		endif

		if (Tax_Reform_KW) then 
			call Find_Capital_and_Wealth_Tax(compute_exp,Simul_Switch)
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
		! CALL Write_Benchmark_Results(.false.)
	if (compute_bench) then
		print*,"	Computing equilibrium distribution"
		CALL FIND_DBN_EQ
		print*,"	Computing government spending"
		CALL GOVNT_BUDGET
		print*,"	Computing Value Function"
		! CALL COMPUTE_VALUE_FUNCTION_SPLINE
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction)
		print*,"	Computing Firm Value Function"
		CALL Firm_Value
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
			! tauK = 0.0_dp 
			! call Find_TauW_Threshold(DBN1,W_bench)  
			! print*,' ' 
			! print*,'W_Bench=',W_bench
			! STOP
		CALL GOVNT_BUDGET
	end if 
		print*,"	Computing satitics"
		CALL COMPUTE_STATS
		print*,"	Writing variables"
		CALL WRITE_VARIABLES(1)
		if (Simul_Switch) then 
			print*,"	Simulation"
			CALL SIMULATION(solving_bench)
		endif
		

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
		Cons_bench          = Cons           
		Hours_bench         = Hours
		Aprime_bench        = Aprime 
		V_Pr_bench          = V_Pr
		V_Pr_nb_bench       = V_Pr_nb 

		SSC_Payments_bench  = SSC_Payments

		write(*,*) "Benchmark variables"
		write(*,*) "GBAR=",GBAR,"EBAR=",EBAR,"NBAR=",NBAR,"QBAR=",QBAR,"P=",P,"wage=",wage,'R=',R

		! Deallocate variables
		if (compute_bench.eqv..false.) then 
		deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )
		endif 

		! Call Simulation_Life_Cycle_Patterns(solving_bench)


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
				CALL GOVNT_BUDGET

				! Get new G
				GBAR_exp = GBAR 
				! Iteratioins  
				tauWindx = tauWindx + 1.0_DP   
				write(*,*) "Bracketing GBAR: tauW_bt=", tauW_bt*100, "And tauW_at=", tauW_at*100
				print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
				print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
			ENDDO

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
				CALL FIND_DBN_EQ
				CALL GOVNT_BUDGET

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
				    CALL FIND_DBN_EQ
				    CALL GOVNT_BUDGET
				    GBAR_exp = GBAR
				    print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
					print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
					print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
					print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
				ENDDO

		! Compute value function and store policy functions, value function and distribution in file
		! CALL COMPUTE_VALUE_FUNCTION_SPLINE 
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction)
		CALL Firm_Value

	endif 
	
	CALL Write_Experimental_Results(compute_exp)
	CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
	K_mat  = K_Matrix(R,P)
	Pr_mat = Profit_Matrix(R,P)
	CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
	CALL ComputeLaborUnits(EBAR,wage)
	CALL GOVNT_BUDGET


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

	print*,"	Efficiency Computation"
		CALL Hsieh_Klenow_Efficiency(solving_bench)


end Subroutine Solve_Experiment



!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Find_Capital_and_Wealth_Tax(compute_exp,Simul_Switch)
	use parameters
	use global 
	use programfunctions
	use Simulation_Module
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: compute_exp, Simul_Switch
	real(dp) :: brentvaluet, Opt_TauK, Opt_tauW, CE2_NB, max_CE2_NB
	integer  :: tauindx 
	character(100) :: folder_aux

	! Solve Benchmark
	call Solve_Benchmark(.false.,.false.)

	! Set Folder
	folder_aux = Result_Folder
	Result_Folder = trim(folder_aux)//'Tax_Reform_KW/'
	call system( 'mkdir -p ' // trim(Result_Folder) )

	! Set initial wealth tax 
	tauW_bt = 0.0665_dp!0.01_dp 
	tauW_at = 0.0665_dp!0.01_dp 

	! Load current solution
	CALL Write_Experimental_Results(.false.)

	! Grid for tauK 
	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_k_9.txt', STATUS='replace')
	CLOSE (unit=77) 
	max_CE2_NB = -10000.0_dp 
	do tauindx=-202,-240,-2
    	tauK = real(tauindx,8)/100_DP
    	print*, 'Capital Tax Grid: tauK=',tauK
    	CE2_NB = Tax_Reform_Welfare(tauK)
	    if (CE2_NB .gt. max_CE2_NB) then
	        max_CE2_NB = CE2_NB
			Opt_tauK = tauK
			Opt_tauW = tauW_at 
		endif

		OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_k_9.txt', STATUS='old', POSITION='append')
		WRITE  (UNIT=77, FMT=*) tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), & 
			      &  MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0), &
			      &  wage, sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)), &
			      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, &
		      	  & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_Pop, Av_Util_NB
      	CLOSE (unit=77) 
	enddo 

	! Find optimal capital income tax (subsidy)
	! brentvaluet = brent(Opt_TauK-0.05_dp,Opt_TauK,Opt_TauK+0.05_dp,Tax_Reform_Welfare , brent_tol, Opt_TauK) 

	! Set Optimal Taxes
	tauK = Opt_TauK
	tauW_at = Opt_tauW
	CE2_NB = Tax_Reform_Welfare(Opt_TauK)

	! Save Results and Simulation
	CALL Write_Experimental_Results(.true.)
	if ((Simul_Switch)) then 
	 	print*,"	Experiment Simulation"
		CALL SIMULATION(solving_bench)
	endif
	Call Simulation_Life_Cycle_Patterns(solving_bench)


	! Deallocate variables
		deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )

	! print*,"	Efficiency Computation"
	! 	CALL Hsieh_Klenow_Efficiency(solving_bench)

end Subroutine Find_Capital_and_Wealth_Tax


!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Experiment_tauC(compute_exp,Simul_Switch)
	use parameters
	use global 
	use programfunctions
	use Simulation_Module
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: compute_exp, Simul_Switch
	real(dp)            :: tauCmin=0.075_dp, tauCindx=0.0_dp, tauCinc=0.05_dp, tauC_up, tauC_low

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

	if (compute_exp) then 
		! Find wealth taxes that balances budget
		print*, "	Computing Wealth Tax to balance the budget"
			! Set initial value for G in experimental economy and for wealth taxes
			GBAR_exp = 0.0_DP
			tauC  = tauCmin
			tauCindx = 0.0_DP
			! Solve for the model increasing wealth taxes until revenue is enough to finance G_benchamark
			DO WHILE (GBAR_exp .lt. GBAR_bench)
				! Set old G and new value of tauW
				GBAR_exp_old = GBAR_exp
				tauC = tauCmin + tauCindx * tauCinc
				! Solve the model
				CALL FIND_DBN_EQ
				CALL GOVNT_BUDGET

				! Get new G
				GBAR_exp = GBAR 
				! Iteratioins  
				tauCindx = tauCindx + 1.0_DP   
				write(*,*) "Bracketing GBAR: tauC=", tauC*100
				print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
			ENDDO

			! Set tauW as weighted average of point in  the grid to balance budget more precisely
				tauC_up  = tauC
				tauC_low = tauC  -  tauCinc
				tauC     = tauC_low + tauCinc * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
				print*,''
				print*,'GBAR bracketed by taxes:'
				print*,'tauC_low =', tauC_low*100, '% tauC_up=', tauC_up*100, '% tauC=', tauC*100, "%"
				print*,''

			! Solve (again) experimental economy
				CALL FIND_DBN_EQ
				CALL GOVNT_BUDGET

			! Find tauW that exactly balances the budget (up to precisioin 0.1) using bisection
				GBAR_exp = GBAR
				print*,"Gbar at midpoint of bracket and GBAR at benchmark"
				print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
				print*,''
				print*,'Bisection for TauW:'
				DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.001 ) ! as long as the difference is greater than 0.1% continue
				    if (GBAR_exp .gt. GBAR_bench ) then
				        tauC_up  = tauC
				    else
				        tauC_low = tauC
				    endif
				    tauC = (tauC_low + tauC_up)/2.0_DP
				    CALL FIND_DBN_EQ
				    CALL GOVNT_BUDGET
				    GBAR_exp = GBAR
				    print*,'tauC_low =', tauC_low*100, '% tauC_up=', tauC_up*100, '% tauC=', tauC*100, "%"
					print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
				ENDDO

		! Compute value function and store policy functions, value function and distribution in file
		! CALL COMPUTE_VALUE_FUNCTION_SPLINE 
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction)
		CALL Firm_Value

	endif 
	
	CALL Write_Experimental_Results(compute_exp)
	CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
	K_mat  = K_Matrix(R,P)
	Pr_mat = Profit_Matrix(R,P)
	CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
	CALL ComputeLaborUnits(EBAR,wage)
	CALL GOVNT_BUDGET


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
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime
		V_Pr_exp          = V_Pr 
		V_Pr_nb_exp  	  = V_Pr_nb

	! Compute moments
	CALL COMPUTE_STATS
	
	! Compute welfare gain between economies
	CALL COMPUTE_WELFARE_GAIN


	print*,'---------------------------'
	print*,''
	print*,'Output Gain Prct=', 100.0_DP*(Y_exp/Y_bench-1.0) 
	print*,''
	print*,'---------------------------'

	print*," "
	print*,"Wealth_factor=",Wealth_factor
	print*," "

	! Write experimental results in output.txt
	CALL WRITE_VARIABLES(0)
	if ((Simul_Switch)) then 
	 	print*,"	Experiment Simulation"
		CALL SIMULATION(solving_bench)
	endif


	! Deallocate variables
		deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )

	print*,"	Efficiency Computation"
		CALL Hsieh_Klenow_Efficiency(solving_bench)


end Subroutine Solve_Experiment_tauC

!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Experiment_Fixed_Policy_Functions(compute_exp_pf,Simul_Switch)
	use parameters
	use global 
	use programfunctions
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
				CALL GOVNT_BUDGET

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
				CALL GOVNT_BUDGET

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
				    CALL GOVNT_BUDGET
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
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction)
		CALL Firm_Value

	endif 

	
	CALL Write_Experimental_Results(compute_exp_pf)
	CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
	K_mat  = K_Matrix(R,P)
	Pr_mat = Profit_Matrix(R,P)
	CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
	CALL ComputeLaborUnits(EBAR,wage)
	CALL GOVNT_BUDGET


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

Subroutine Solve_Experiment_Fixed_PF_Interp(compute_exp_pf_interp,Simul_Switch)
	use parameters
	use global 
	use programfunctions
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
				CALL GOVNT_BUDGET

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
				CALL GOVNT_BUDGET

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
				    CALL GOVNT_BUDGET
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
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction)
		CALL Firm_Value

	endif 

	
	CALL Write_Experimental_Results(compute_exp_pf_interp)
	CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
	CALL ComputeLaborUnits(EBAR,wage)
	CALL GOVNT_BUDGET


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
	! 			CALL GOVNT_BUDGET

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
	! 			CALL GOVNT_BUDGET

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
	! 			    CALL GOVNT_BUDGET
	! 			    GBAR_exp = GBAR
	! 			    print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
	! 				print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
	! 				print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
	! 				print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
	! 			ENDDO

	! 	! Compute value function and store policy functions, value function and distribution in file
	! 	! CALL COMPUTE_VALUE_FUNCTION_SPLINE 
	! 	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction)
	! 	CALL Firm_Value

	! endif 
	
	! CALL Write_Experimental_Results(compute_exp_prices)
	CALL Write_Experimental_Results(.false.)
	CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
	K_mat  = K_Matrix(R,P)
	Pr_mat = Profit_Matrix(R,P)
	CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
	CALL ComputeLaborUnits(EBAR,wage)
	CALL GOVNT_BUDGET


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
	! 			CALL GOVNT_BUDGET

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
	! 			CALL GOVNT_BUDGET

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
	! 			    CALL GOVNT_BUDGET
	! 			    GBAR_exp = GBAR
	! 			    print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
	! 				print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
	! 				print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
	! 				print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
	! 			ENDDO

	! 	! Compute value function and store policy functions, value function and distribution in file
	! 	! CALL COMPUTE_VALUE_FUNCTION_SPLINE 
	! 	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction)
	! 	CALL Firm_Value

	! endif 
	
	! CALL Write_Experimental_Results(compute_exp_prices)
 	CALL Write_Experimental_Results(.false.)
	CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
	K_mat  = K_Matrix(R,P)
	Pr_mat = Profit_Matrix(R,P)
	CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
	CALL ComputeLaborUnits(EBAR,wage)
	CALL GOVNT_BUDGET


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

Subroutine Solve_Opt_Tax(Opt_Tax_KW,Simul_Switch)
	use parameters 
	use global
	use Opt_Tax_Parameters
	use Opt_Tax_Functions
	use programfunctions
	use Simulation_Module
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: Opt_Tax_KW,Simul_Switch


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
	
	print*,'Optimal Tax Loop'
	If (Opt_Tax_KW) then
		PRINT*,''
		Print*,'--------------- OPTIMAL CAPITAL TAXES -----------------'
		PRINT*,''
    	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_k_no_tax_3.txt', STATUS='replace')
    	CLOSE (unit=77) 

	    DO tauindx=0,25,1
	    	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_k_no_tax_3.txt', STATUS='old', POSITION='append')
            
            tauK        = real(tauindx,8)/100_DP
            brentvaluet = - EQ_WELFARE_GIVEN_TauK(tauK)

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
				Cons_exp          = Cons           
				Hours_exp         = Hours
				Aprime_exp        = Aprime 

				! Compute moments
				CALL COMPUTE_STATS
				CALL GOVNT_BUDGET
				
				! Compute welfare gain between economies
				CALL COMPUTE_WELFARE_GAIN

				! Write experimental results in output.txt
				CALL WRITE_VARIABLES(0)

		    if (brentvaluet .gt. maxbrentvaluet) then
		        maxbrentvaluet = brentvaluet
				OPT_tauK = tauK
				OPT_psi  = psi
			endif

			! Print Results 
		    print*, 'tauK=', tauK, 'YBAR=', YBAR, & 
		    	  & 'Av. Util=', sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
		      
		    WRITE  (UNIT=77, FMT=*) tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), & 
			      &  MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0), &
			      &  wage, sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)),  &
				!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)), &
			      &100*( (sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)) /&
			      &sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))) &
			      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
			      &100*( (sum(ValueFunction(:,:,:,:,:,:)*DBN1(:,:,:,:,:,:))/sum(DBN1(:,:,:,:,:,:)) /&
			      &sum(ValueFunction_bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:,:))) &
			      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
			      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, &
		      	  & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_Pop, Av_Util_NB, brentvaluet
	      	CLOSE (unit=77) 
	    	Call Write_Experimental_Results(.true.)
	    ENDDO 

	    tauK = OPT_tauK
		psi  = OPT_psi

	 	OPEN (UNIT=77, FILE=trim(Result_Folder)//'stat_opt_tau_k.txt', STATUS='replace')

		tauK = OPT_tauK
		psi  = OPT_psi
		call Find_Opt_Tax(Opt_Tax_KW,Opt_TauK,Opt_TauK-0.02_dp,Opt_TauK+0.02_dp) 

		tauK     = OPT_tauK
		OPT_psi  = psi

		! Compute moments
		CALL COMPUTE_STATS

		print*, "Optimal tau_K=", tauK, "Optimal psi=", psi

		WRITE  (UNIT=77, FMT=*) tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), & 
		    &  MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0), &
		    &  wage, sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)),  &
		!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)), &
	 	    &100*( (sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)) /&
	 	    &sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))) &
	 	    &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
	 	    &100*( (sum(ValueFunction(:,:,:,:,:,:)*DBN1(:,:,:,:,:,:))/sum(DBN1(:,:,:,:,:,:)) /&
	 	    &sum(ValueFunction_bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:,:))) &
	 	    &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
	 	    & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, &
  		    & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_Pop, Av_Util_NB


		CLOSE (UNIT=77)
		Call Write_Experimental_Results(.true.)

	else
		PRINT*,''
		Print*,'--------------- OPTIMAL WEALTH TAXES -----------------'
		PRINT*,''
    	! Set Y_a_threshold
			call Find_TauW_Threshold(DBN_bench,W_bench)  
			Y_a_threshold = Threshold_Factor*Ebar_bench !0.75_dp
			Wealth_factor = Y_a_threshold/W_bench
    	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_w_no_tax_2.txt', STATUS='replace')
    	CLOSE (unit=77) 

    	! CALL Write_Experimental_Results(.false.)

	    DO tauindx=40,60,2
	    	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_w_no_tax_2.txt', STATUS='old', POSITION='append')

            tauw_at     = real(tauindx,8)/1000_DP
            brentvaluet = - EQ_WELFARE_GIVEN_TauW(tauW_at)

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
				Cons_exp          = Cons           
				Hours_exp         = Hours
				Aprime_exp        = Aprime 

				! Compute moments
				CALL COMPUTE_STATS
				CALL GOVNT_BUDGET
				
				! Compute welfare gain between economies
				CALL COMPUTE_WELFARE_GAIN

				! Write experimental results in output.txt
				CALL WRITE_VARIABLES(0)

		    if (brentvaluet .gt. maxbrentvaluet) then
		        maxbrentvaluet = brentvaluet
				OPT_tauW = tauW_at
				OPT_psi  = psi
			endif

			! Print Results 
		    print*, 'tauW=', tauW_at, 'YBAR=', YBAR, & 
		    	  & 'Av. Util=', sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
		      
		    WRITE  (UNIT=77, FMT=*) tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), & 
			      &  MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0), &
			      &  wage, sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)),  &
				!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)), &
			      &100*( (sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)) /&
			      &sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))) &
			      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
			      &100*( (sum(ValueFunction(:,:,:,:,:,:)*DBN1(:,:,:,:,:,:))/sum(DBN1(:,:,:,:,:,:)) /&
			      &sum(ValueFunction_bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:,:))) &
			      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
			      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, &
		      	  & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_Pop, Av_Util_NB, brentvaluet

      	  	CLOSE (unit=77)
		    Call Write_Experimental_Results(.true.)
	    ENDDO 


	    OPEN (UNIT=77, FILE=trim(Result_Folder)//'stat_opt_tau_w.txt', STATUS='replace')

	    ! opt_psi = 0.860830826876844_dp 
	    ! Opt_TauW = 0.031_dp 
		tauW_at = OPT_tauW
		psi     = OPT_psi
		call Find_Opt_Tax(Opt_Tax_KW,Opt_TauW,Opt_TauW-0.002_dp,Opt_TauW+0.002_dp)

		! tauW_at = 0.025_dp
		! call Find_Opt_Tax(Opt_Tax_KW,Opt_TauW,0.023_dp,0.028_dp)

		tauW_at = OPT_tauW
		OPT_psi = psi

		! Compute moments
		CALL COMPUTE_STATS

		print*, "Optimal tau_W=", tauW_at, "Optimal psi=", psi
		
		WRITE  (UNIT=77, FMT=*) tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), & 
	      &  MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0), &
	      &  wage, sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)),  &
		!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)), &
	      &100*( (sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)) /&
	      &sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))) &
	      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
	      &100*( (sum(ValueFunction(:,:,:,:,:,:)*DBN1(:,:,:,:,:,:))/sum(DBN1(:,:,:,:,:,:)) /&
	      &sum(ValueFunction_bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:,:))) &
	      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
	      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, &
      	  & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_Pop, Av_Util_NB


		CLOSE (UNIT=77)
	endif 

	CALL FIND_DBN_EQ
	CALL GOVNT_BUDGET

	! Compute value function and store policy functions, value function and distribution in file
	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction)
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
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime 

	! Compute moments
	CALL COMPUTE_STATS
	
	! Compute welfare gain between economies
	CALL COMPUTE_WELFARE_GAIN

	! Write experimental results in output.txt
	CALL WRITE_VARIABLES(0)
	! if (((theta.eq.1.50_dp)).and.(Threshold_Factor.eq.0.0_dp).and.(Simul_Switch)) then 
	!  	print*,"	Optimal Tax Simulation"
	! 	CALL SIMULATION(solving_bench)
	! endif
	CALL SIMULATION(0)

end Subroutine Solve_Opt_Tax


!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Opt_Tax_K_and_W(Simul_Switch)
	use parameters 
	use global
	use Opt_Tax_Parameters
	use Opt_Tax_Functions
	use programfunctions
	use Simulation_Module
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: Simul_Switch
	INTEGER  :: tauindx_w
	Logical  :: Opt_Tax_KW


 	!====================================================================================================
	PRINT*,''
	Print*,'--------------- SOLVING OPTIMAL TAXES - Joint K and W -----------------'
	PRINT*,''
	
	! Experiment economy
		solving_bench=0
	
	! Set initial taxes for finding optimal ones
		tauK     = 0.0_DP
		tauW_at  = 0.0_DP
		Opt_TauK = 0.0_DP
		Opt_TauW = 0.0_DP
		maxbrentvaluet=-10000.0_DP
	
	print*,'Optimal Tax Loop'
    	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_kw.txt', STATUS='replace')
    	CLOSE (unit=77) 
    	! CALL Write_Experimental_Results(.false.)

    	! Set Y_a_threshold
			call Find_TauW_Threshold(DBN_bench,W_bench)  
			Y_a_threshold = Threshold_Factor*Ebar_bench !0.75_dp
			Wealth_factor = Y_a_threshold/W_bench

		Opt_Tax_KW = .true.

		DO tauindx_w=-00,20,5
			tauw_at     = real(tauindx_w,8)/1000_DP

	    DO tauindx=-40,-20,5

	    	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_kw.txt', STATUS='old', POSITION='append')
            
            tauK        = real(tauindx,8)/100_DP

            print*, 'Optimal Tax Loop', 'tauW',tauW_at,'tauK=', tauK

            brentvaluet = - EQ_WELFARE_GIVEN_TauK(tauK)

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
				Cons_exp          = Cons           
				Hours_exp         = Hours
				Aprime_exp        = Aprime 

				! Compute moments
				CALL COMPUTE_STATS
				CALL GOVNT_BUDGET
				
				! Compute welfare gain between economies
				CALL COMPUTE_WELFARE_GAIN

				! Write experimental results in output.txt
				CALL WRITE_VARIABLES(0)

		    if (brentvaluet .gt. maxbrentvaluet) then
		        maxbrentvaluet = brentvaluet
				OPT_tauK = tauK
				OPT_psi  = psi
				OPT_tauW = tauW_at
			endif

			! Print Results 
		    print*, 'tauW',tauW_at,'tauK=', tauK, 'YBAR=', YBAR, & 
		    	  & 'Av. Util=', sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
		      
		    WRITE  (UNIT=77, FMT=*) tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), & 
			      &  MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0), &
			      &  wage, sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)),  &
				!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)), &
			      &100*( (sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)) /&
			      &sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))) &
			      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
			      &100*( (sum(ValueFunction(:,:,:,:,:,:)*DBN1(:,:,:,:,:,:))/sum(DBN1(:,:,:,:,:,:)) /&
			      &sum(ValueFunction_bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:,:))) &
			      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
			      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, &
		      	  & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_Pop, Av_Util_NB
	      	CLOSE (unit=77) 
	    	Call Write_Experimental_Results(.true.)
	    ENDDO 
	    ENDDO



	    OPEN (UNIT=77, FILE=trim(Result_Folder)//'stat_opt_tau_kw.txt', STATUS='replace')

		tauw_at = OPT_tauW
		tauK = OPT_tauK
		psi  = OPT_psi
		call Find_Opt_Tax(Opt_Tax_KW,Opt_TauK,Opt_TauK-0.05_dp,Opt_TauK+0.05_dp) 

		tauK     = OPT_tauK
		OPT_psi  = psi

		! Compute moments
		CALL COMPUTE_STATS

		print*, "Optimal tau_W=",tauW_at, "Optimal tau_K=", tauK, "Optimal psi=", psi

		WRITE  (UNIT=77, FMT=*) tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), & 
	      &  MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0), &
	      &  wage, sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)),  &
		!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)), &
	      &100*( (sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)) /&
	      &sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))) &
	      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
	      &100*( (sum(ValueFunction(:,:,:,:,:,:)*DBN1(:,:,:,:,:,:))/sum(DBN1(:,:,:,:,:,:)) /&
	      &sum(ValueFunction_bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:,:))) &
	      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
	      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, &
      	  & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_Pop, Av_Util_NB


		CLOSE (UNIT=77)
		Call Write_Experimental_Results(.true.)

	print*,'Optimal Tax Loop Finished'

	CALL FIND_DBN_EQ
	CALL GOVNT_BUDGET

	! Compute value function and store policy functions, value function and distribution in file
	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction)
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
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime 

	! Compute moments
	CALL COMPUTE_STATS
	
	! Compute welfare gain between economies
	CALL COMPUTE_WELFARE_GAIN

	! Write experimental results in output.txt
	CALL WRITE_VARIABLES(0)
	! if (((theta.eq.1.50_dp)).and.(Threshold_Factor.eq.0.0_dp).and.(Simul_Switch)) then 
	!  	print*,"	Optimal Tax Simulation"
	! 	CALL SIMULATION(solving_bench)
	! endif
	CALL SIMULATION(0)

end Subroutine Solve_Opt_Tax_K_and_W

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
	use Toolbox
	use omp_lib
	implicit none 
	real(DP) :: OPT_Threshold
	INTEGER  :: Threshold_ind


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
	
	print*,'Optimal Tax Loop'
	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_w_threshold_2.txt', STATUS='replace')
		DO Threshold_ind = 0,3

		Threshold_Factor = real(Threshold_ind,8)/4.0_dp
		print*, ' Threshold_Factor=',Threshold_Factor
	
		PRINT*,''
		Print*,'--------------- OPTIMAL WEALTH TAXES - Threshold -----------------'
		PRINT*,''
    	! Set Y_a_threshold
			call Find_TauW_Threshold(DBN_bench,W_bench)  
			Y_a_threshold = Threshold_Factor*Ebar_bench !0.75_dp
			Wealth_factor = Y_a_threshold/W_bench
    	
	    DO tauindx=25,50
            tauw_at     = real(tauindx,8)/1000_DP
            brentvaluet = - EQ_WELFARE_GIVEN_TauW(tauW_at)

            ! Compute moments
			CALL COMPUTE_STATS

		    if (brentvaluet .gt. maxbrentvaluet) then
		        maxbrentvaluet = brentvaluet
				OPT_tauW = tauW_at
				OPT_psi  = psi
				OPT_Threshold = Threshold_Factor
			endif

			! Print Results 
		    print*, 'tauW=', tauW_at, 'YBAR=', YBAR, & 
		    	  & 'Av. Util=', sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
		      
		    WRITE  (UNIT=77, FMT=*) tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), &
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
		      & Threshold_Factor
	    ENDDO 

	    ENDDO

    CLOSE (unit=77)


	! Evaluate optimal point in grid
		tauW_at = OPT_tauW
		psi 	= OPT_psi
		Threshold_Factor = OPT_Threshold
			! Set Y_a_threshold
			call Find_TauW_Threshold(DBN_bench,W_bench)  
			Y_a_threshold = Threshold_Factor*Ebar_bench !0.75_dp
			Wealth_factor = Y_a_threshold/W_bench


	CALL FIND_DBN_EQ
	CALL GOVNT_BUDGET

	! Compute value function and store policy functions, value function and distribution in file
	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction)
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
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime 

	! Compute moments
	CALL COMPUTE_STATS
	
	! Compute welfare gain between economies
	CALL COMPUTE_WELFARE_GAIN

	! Write experimental results in output.txt
	CALL WRITE_VARIABLES(0)

	CALL SIMULATION(0)

end Subroutine Solve_Opt_Threshold



!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Opt_Tau_C(Opt_Tax_KW)
	use parameters 
	use global
	use Opt_Tax_Parameters
	use Opt_Tax_Functions
	use Simulation_Module
	use programfunctions
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: Opt_Tax_KW
	real(DP) :: OPT_tauC
	INTEGER  :: tauC_ind


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
	
	print*,'Optimal Tax Loop'
	
		

	If (Opt_Tax_KW) then
		PRINT*,''
		Print*,'--------------- OPTIMAL CAPITAL TAXES - Consumption Taxes -----------------'
		PRINT*,''
    	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_k_cons_tax_2.txt', STATUS='replace')
    	CLOSE (unit=77) 
    	CALL Write_Experimental_Results(.false.)
    	! psi = 1.9

    	DO tauC_ind = 11,15,1

    		OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_k_cons_tax_2.txt', STATUS='old', POSITION='append')

			tauC = real(tauC_ind,8)/10.0_dp
			print*, ' '
			print*, ' Consumption Taxes=',tauC
			print*, ' '

			! psi = 1.4_dp

		    DO tauindx=-50,-50,10
	            tauK        = real(tauindx,8)/100_DP
	            brentvaluet = - EQ_WELFARE_GIVEN_TauK(tauK)

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
				Cons_exp          = Cons           
				Hours_exp         = Hours
				Aprime_exp        = Aprime 

				! Compute moments
				CALL COMPUTE_STATS
				CALL GOVNT_BUDGET
				
				! Compute welfare gain between economies
				CALL COMPUTE_WELFARE_GAIN

				! Write experimental results in output.txt
				CALL WRITE_VARIABLES(0)

			    if (brentvaluet .gt. maxbrentvaluet) then
			        maxbrentvaluet = brentvaluet
					OPT_tauK = tauK
					OPT_psi  = psi
					OPT_tauC = tauC
				endif

				! Print Results 
			    print*, 'tauK=', tauK, 'YBAR=', YBAR, & 
			    	  & 'Av. Util=', sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
			      
			    WRITE  (UNIT=77, FMT=*) tauC, tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), & 
			      &  MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0), &
			      &  wage, sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)),  &
				!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)), &
			      &100*( (sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)) /&
			      &sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))) &
			      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
			      &100*( (sum(ValueFunction(:,:,:,:,:,:)*DBN1(:,:,:,:,:,:))/sum(DBN1(:,:,:,:,:,:)) /&
			      &sum(ValueFunction_bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:,:))) &
			      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
			      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, &
		      	  & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_NB
		    ENDDO 

	    	CLOSE (unit=77) 
	    	Call Write_Experimental_Results(.true.)
	    ENDDO


	    OPEN (UNIT=77, FILE=trim(Result_Folder)//'stat_opt_tau_k_cons_tax_2.txt', STATUS='replace')

		tauK = OPT_tauK
		psi  = OPT_psi
		tauC = OPT_tauC
		! call Find_Opt_Tax(Opt_Tax_KW,Opt_TauK,Opt_TauK-0.05_dp,Opt_TauK+0.05_dp) 

		tauK     = OPT_tauK
		OPT_psi  = psi

		! Compute moments
		CALL COMPUTE_STATS

		print*, "Optimal tau_K=", tauK, "Optimal psi=", psi, 'Optimal tauC=',tauC

		WRITE  (UNIT=77, FMT=*) tauC, tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), & 
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
		      & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_NB


		CLOSE (UNIT=77)
		Call Write_Experimental_Results(.true.)


	else 
		PRINT*,''
		Print*,'--------------- OPTIMAL WEALTH TAXES - Consumption Taxes -----------------'
		PRINT*,''
    	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_w_cons_tax_8.txt', STATUS='replace')
    	CLOSE (unit=77) 

    	CALL Write_Experimental_Results(.false.)
    	! psi = 1.8

    	DO tauC_ind = 10,15,1

    		OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_w_cons_tax_8.txt', STATUS='old', POSITION='append')

			tauC = real(tauC_ind,8)/10.0_dp
			print*, ' '
			print*, ' Consumption Taxes=',tauC
			print*, ' '

			! psi = 0.776_dp
			! psi = 1.50_dp 

		    DO tauindx=-03,03,1
	            tauw_at     = real(tauindx,8)/1000_DP
	            brentvaluet = - EQ_WELFARE_GIVEN_TauW(tauW_at)

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
				Cons_exp          = Cons           
				Hours_exp         = Hours
				Aprime_exp        = Aprime 

				! Compute moments
				CALL COMPUTE_STATS
				CALL GOVNT_BUDGET
				
				! Compute welfare gain between economies
				CALL COMPUTE_WELFARE_GAIN

				! Write experimental results in output.txt
				CALL WRITE_VARIABLES(0)

			    if (brentvaluet .gt. maxbrentvaluet) then
			        maxbrentvaluet = brentvaluet
					OPT_tauW = tauW_at
					OPT_psi  = psi
					OPT_tauC = tauC
				endif

				! Print Results 
			    print*, 'tauW=', tauW_at, 'YBAR=', YBAR, & 
			    	  & 'Av. Util=', sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)), &
			    	  & 'tauC=', tauC,'psi=',psi
			      
			    WRITE  (UNIT=77, FMT=*) tauC, tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), &
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
		      & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_NB
		    ENDDO 
		    CLOSE (unit=77)
		    Call Write_Experimental_Results(.true.)

	    ENDDO


	    OPEN (UNIT=77, FILE=trim(Result_Folder)//'stat_opt_tau_w_cons_tax_8.txt', STATUS='replace')

		tauW_at = OPT_tauW
		psi 	= OPT_psi
		tauC    = OPT_tauC
		! call Find_Opt_Tax(Opt_Tax_KW,Opt_TauW,Opt_TauW-0.001_dp,Opt_TauW+0.001_dp) 

		tauW_at = OPT_tauW
		OPT_psi  = psi

		! Compute moments
		CALL COMPUTE_STATS

		print*, "Optimal tau_W=", tauW_at, "Optimal psi=", psi
		
		WRITE  (UNIT=77, FMT=*) tauC, tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), & 
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
		      & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_NB


		CLOSE (UNIT=77)
		Call Write_Experimental_Results(.true.)
	endif 


	CALL FIND_DBN_EQ
	CALL GOVNT_BUDGET

	! Compute value function and store policy functions, value function and distribution in file
	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction)
	CALL Firm_Value
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
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime 

	! Compute moments
	CALL COMPUTE_STATS
	
	! Compute welfare gain between economies
	CALL COMPUTE_WELFARE_GAIN

	! Write experimental results in output.txt
	CALL WRITE_VARIABLES(0)

	print*,"	Efficiency Computation"
		CALL Hsieh_Klenow_Efficiency(solving_bench)

	CALL SIMULATION(0)

    

end Subroutine Solve_Opt_Tau_C


!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Opt_Tau_CX(Opt_Tax_KW)
	use parameters 
	use global
	use Opt_Tax_Parameters
	use Opt_Tax_Functions
	use Simulation_Module
	use programfunctions
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: Opt_Tax_KW
	real(DP) :: OPT_tauC
	INTEGER  :: tauC_ind


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
	
	print*,'Optimal Tax Loop'
	
		

	If (Opt_Tax_KW) then
		PRINT*,''
		Print*,'--------------- OPTIMAL CAPITAL TAXES - Consumption Taxes -----------------'
		PRINT*,''
    	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_c_k.txt', STATUS='replace')
    	CLOSE (unit=77) 
    	! CALL Write_Experimental_Results(.false.)

    	! Start tauK at some level that guarantees enough revenue
    	tauK = 0.25_dp 
    	tauC = 0.20_dp ! Initial Value for consumption taxes
    	
    	DO tauC_ind = 0,2,1

    		OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_c_k.txt', STATUS='old', POSITION='append')

			psi = 0.90_dp + 0.05_dp*real(tauC_ind,8)
			print*, ' '
			print*, ' Consumption Taxes=',tauC
			print*, ' Labor Taxes=',psi
			print*, ' '

            brentvaluet = - EQ_WELFARE_GIVEN_TauC(tauC,Opt_Tax_KW)

            
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
			Cons_exp          = Cons           
			Hours_exp         = Hours
			Aprime_exp        = Aprime 

			! Compute moments
			CALL COMPUTE_STATS
			CALL GOVNT_BUDGET
			
			! Compute welfare gain between economies
			CALL COMPUTE_WELFARE_GAIN

			! Write experimental results in output.txt
			CALL WRITE_VARIABLES(0)

		    if (brentvaluet .gt. maxbrentvaluet) then
		        maxbrentvaluet = brentvaluet
				OPT_tauK = tauK
				OPT_psi  = psi
				OPT_tauC = tauC
			endif

			! Print Results 
		    print*,'';print*, 'tauC',tauC,'tauK=', tauK, 'YBAR=', YBAR, & 
		    	  & 'Av. Util=', sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
		      
		    WRITE  (UNIT=77, FMT=*) tauC, tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), & 
		      &  MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0), &
		      &  wage, sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)),  &
			!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)), &
		      &100*( (sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)) /&
		      &sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))) &
		      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
		      &100*( (sum(ValueFunction(:,:,:,:,:,:)*DBN1(:,:,:,:,:,:))/sum(DBN1(:,:,:,:,:,:)) /&
		      &sum(ValueFunction_bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:,:))) &
		      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
		      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, &
		      & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_NB

	    	CLOSE (unit=77) 
	    	Call Write_Experimental_Results(.true.)
	    ENDDO

    	DO tauC_ind = 4,20,1

    		OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_c_k.txt', STATUS='old', POSITION='append')

			tauC = real(tauC_ind,8)/20.0_dp
			print*, ' '
			print*, ' Consumption Taxes=',tauC
			print*, ' '

            brentvaluet = - EQ_WELFARE_GIVEN_TauC(tauC,Opt_Tax_KW)

            
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
			Cons_exp          = Cons           
			Hours_exp         = Hours
			Aprime_exp        = Aprime 

			! Compute moments
			CALL COMPUTE_STATS
			CALL GOVNT_BUDGET
			
			! Compute welfare gain between economies
			CALL COMPUTE_WELFARE_GAIN

			! Write experimental results in output.txt
			CALL WRITE_VARIABLES(0)

		    if (brentvaluet .gt. maxbrentvaluet) then
		        maxbrentvaluet = brentvaluet
				OPT_tauK = tauK
				OPT_psi  = psi
				OPT_tauC = tauC
			endif

			! Print Results 
		    print*,'';print*, 'tauC',tauC,'tauK=', tauK, 'YBAR=', YBAR, & 
		    	  & 'Av. Util=', sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
		      
		    WRITE  (UNIT=77, FMT=*) tauC, tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), & 
		      &  MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0), &
		      &  wage, sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)),  &
			!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)), &
		      &100*( (sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)) /&
		      &sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))) &
		      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
		      &100*( (sum(ValueFunction(:,:,:,:,:,:)*DBN1(:,:,:,:,:,:))/sum(DBN1(:,:,:,:,:,:)) /&
		      &sum(ValueFunction_bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:,:))) &
		      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
		      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, &
		      & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_NB

	    	CLOSE (unit=77) 
	    	Call Write_Experimental_Results(.true.)
	    ENDDO


	    OPEN (UNIT=77, FILE=trim(Result_Folder)//'stat_opt_tau_c_k.txt', STATUS='replace')

		tauK = OPT_tauK
		psi  = OPT_psi
		tauC = OPT_tauC

		! Compute moments
		CALL COMPUTE_STATS

		print*, "Optimal tau_K=", tauK, "Optimal psi=", psi, 'Optimal tauC=',tauC

		WRITE  (UNIT=77, FMT=*) tauC, tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), & 
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
		      & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_NB


		CLOSE (UNIT=77)
		Call Write_Experimental_Results(.true.)


	else 
		PRINT*,''
		Print*,'--------------- OPTIMAL WEALTH TAXES - Consumption Taxes -----------------'
		PRINT*,''
    	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_c_w.txt', STATUS='replace')
    	CLOSE (unit=77) 

    	! Start tauW at some level that guarantees enough revenue
    	tauW_at = 0.03_dp

    	! CALL Write_Experimental_Results(.false.)

    	DO tauC_ind = 1,8,1

    		OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_c_w.txt', STATUS='old', POSITION='append')

			tauC = real(tauC_ind,8)/10.0_dp
			print*, ' '
			print*, ' Consumption Taxes=',tauC
			print*, ' '

            brentvaluet = - EQ_WELFARE_GIVEN_TauC(tauC,Opt_Tax_KW)

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
			Cons_exp          = Cons           
			Hours_exp         = Hours
			Aprime_exp        = Aprime 

			! Compute moments
			CALL COMPUTE_STATS
			CALL GOVNT_BUDGET
			
			! Compute welfare gain between economies
			CALL COMPUTE_WELFARE_GAIN

			! Write experimental results in output.txt
			CALL WRITE_VARIABLES(0)

		    if (brentvaluet .gt. maxbrentvaluet) then
		        maxbrentvaluet = brentvaluet
				OPT_tauW = tauW_at
				OPT_psi  = psi
				OPT_tauC = tauC
			endif

			! Print Results 
		    print*,''; print*, 'tauC=', tauC,'psi=',psi,'tauW=', tauW_at, 'YBAR=', YBAR, & 
		    	  & 'Av. Util=', sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
		      
		    WRITE  (UNIT=77, FMT=*) tauC, tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), &
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
		      & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_NB

		    CLOSE (unit=77)
		    Call Write_Experimental_Results(.true.)

	    ENDDO


	    OPEN (UNIT=77, FILE=trim(Result_Folder)//'stat_opt_tau_c_w.txt', STATUS='replace')

		tauW_at = OPT_tauW
		psi 	= OPT_psi
		tauC    = OPT_tauC

		! Compute moments
		CALL COMPUTE_STATS

		print*, "Optimal tau_C=",tauC,"Optimal tau_W=", tauW_at, "Optimal psi=", psi
		
		WRITE  (UNIT=77, FMT=*) tauC, tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), & 
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
		      & GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, Av_Util_NB


		CLOSE (UNIT=77)
		Call Write_Experimental_Results(.true.)
	endif 


	CALL FIND_DBN_EQ
	CALL GOVNT_BUDGET

	! Compute value function and store policy functions, value function and distribution in file
	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction)
	CALL Firm_Value
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
		Cons_exp          = Cons           
		Hours_exp         = Hours
		Aprime_exp        = Aprime 

	! Compute moments
	CALL COMPUTE_STATS
	
	! Compute welfare gain between economies
	CALL COMPUTE_WELFARE_GAIN

	! Write experimental results in output.txt
	CALL WRITE_VARIABLES(0)

	print*,"	Efficiency Computation"
		CALL Hsieh_Klenow_Efficiency(solving_bench)

	CALL SIMULATION(0)

    

end Subroutine Solve_Opt_Tau_CX


!========================================================================================
!========================================================================================
!========================================================================================


! SUBROUTINE CALIBRATION_TRIALS
! 	use parameters
! 	use global 
! 	use programfunctions
! 	use Toolbox
! 	use omp_lib

! 	IMPLICIT NONE
! 	real(DP)::  beta_L, beta_H, sigmaz_L, sigmaz_H, x_hi_L, x_hi_H
! 	integer :: parindx1,  parindx2, parindx3, parindx4, parindx5, parindx6, n_beta, n_sigmaz, n_x_hi
! 	real(DP), dimension(6):: paramsL, paramsH
! 	real(DP), dimension(3) :: opt_par

! 	!$ call omp_set_num_threads(3)

! 	print*,'SOLVING CALIBRATION'
! 	print*,'-----------------------------------------------------------------------'


! 	! CALIBRATION STARTS

! 	beta_L=0.952_dp
! 	beta_H=0.957_dp

! 	sigmaz_L = 0.11_dp
! 	sigmaz_H = 0.15_dp

! 	x_hi_L  = 3.0_dp 
! 	x_hi_H  = 3.0_dp


! 	n_beta   = 3
! 	n_sigmaz = 3
! 	n_x_hi   = 1

! 	Min_SSE_Moments=1000000.0_DP

! 	DO parindx1=1,n_beta
! 	DO parindx6=1,n_x_hi
! 	DO parindx5=1,n_sigmaz

! 	    beta  		= beta_L  	+ real(parindx1-1,8) *(beta_H-beta_L)/max(real(n_beta-1,8),1.0_DP)
! 	    sigma_z_eps = sigmaz_L 	+ real(parindx5-1,8)*(sigmaz_H-sigmaz_L) / max(real(n_sigmaz-1,8),1.0_DP)
! 	    x_hi 		= x_hi_L 	+ real(parindx6-1,8)*(x_hi_H-x_hi_L) / max(real(n_x_hi-1,8),1.0_DP)

! 	    print*, ' '
! 	    print*,'parameters=', beta, sigma_z_eps, x_hi
	        
	    
! 		! Benchmark economy
! 			solving_bench=1

! 		! Set taxes for benchmark economy
! 			tauK = 0.25_DP
! 			tauW_bt = 0.00_DP
! 			tauW_at = 0.00_DP
! 			Y_a_threshold = 0.00_DP 

! 	    CALL INITIALIZE
! 	    CALL FIND_DBN_EQ
! 	    CALL Firm_Value
! 	    CALL COMPUTE_MOMENTS
! 	    SSE_Moments = (1.0-Wealth_Output/3.0_DP)**2.0_DP  + (1.0_DP-FW_top_x_share(4)/0.357_DP)**2.0_DP  & !+ (FW_top_x_share(3)-0.75_DP)**2.0_DP &
! 	                   & + (1.0_DP-Std_Log_Earnings_25_60 / 0.8_DP)**2.0_DP + (1.0_DP-meanhours_25_60/0.4_DP)**2.0_DP !&
! 	                   !& + (1.0_DP-MeanReturn/0.069_DP)**2.0_DP
	                   
! 	!    SSE_Moments = (1.0-Wealth_Output/3.0_DP)**2.0_DP + (1.0_DP-MeanReturn/0.069_DP)**2.0_DP
! 		print*, ' '
! 		print*,'parameters=',beta, sigma_z_eps, x_hi,'SSE_Moments =',SSE_Moments
! 		print*,'Wealth_Output', Wealth_Output,'Top 1%',FW_top_x_share(4),'Std_Log_E',Std_Log_Earnings_25_60,'meanhours',meanhours_25_60

! 	    IF (SSE_Moments .lt. Min_SSE_Moments ) THEN
! 	        Min_SSE_Moments =SSE_Moments
! 	        !params= [ beta, mu_z, rho_z, sigma_z_eps, sigma_lambda_eps, gamma ]
! 	        opt_par = [beta, sigma_z_eps, x_hi]
! 	        Min_Moments = [  Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, MeanReturn  ]
! 	    ENDIF
! 	    !CALL WRITE_TO_FILE
	    
! 	ENDDO
! 	ENDDO
! 	ENDDO

! 	print*, opt_par
! 	print*, Min_Moments

! 	! beta=params(1)
! 	! mu_z=params(2)
! 	! rho_z=params(3)
! 	! sigma_z_eps =params(4)
! 	! sigma_lambda_eps = params(5)
! 	! gamma= params(6)

! END SUBROUTINE CALIBRATION_TRIALS

!====================================================================


