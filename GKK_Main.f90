

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
	use Toolbox
	use omp_lib

	IMPLICIT NONE
	! Variables to measure running time
		REAL(DP) :: start_time, finish_time
	! Compute benchmark or load results
		logical  :: compute_bench, compute_exp, Opt_Tax, Opt_Tax_KW, Tax_Reform, Simul_Switch, Calibration_Switch
	! Auxiliary variable for writing file
		character(4)   :: string_theta
		character(100) :: folder_aux

	! Capital Market
		theta = 1.50_dp
	! Threshold 
		Threshold_Factor = 0.00_dp 

	! Switch for solving benchmark or just reading resutls
		Calibration_Switch = .false.
		! If compute_bench==.true. then just read resutls
		! If compute_bench==.false. then solve for benchmark and store results
		Tax_Reform    = .true.
			compute_bench = .true.
			compute_exp   = .false.
		Opt_Tax       = .false.
			Opt_Tax_KW    = .false. ! true=tau_K false=tau_W
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
		Params =[0.962_dp, 0.0_dp, 0.50_dp, 0.387_dp, 0.29_dp, 0.4494_dp] ! alpha=0.4, zgrid 11, m5, alpha=0.4, dep005, mu=090, K/Y=3, Top1PVa=0.36
		

		beta   	= params(1)
		mu_z   	= params(2) ! this is just shifting the z grids. it is zero now.
		rho_z  	= params(3) 
		sigma_z_eps      = 0.125_dp ! params(4)
		sigma_lambda_eps = params(5)
		gamma  	= params(6)
		
		sigma  	= 4.0_dp
		phi    	= (1.0_dp-gamma)/gamma

		x_hi	= 3.0_dp
		x_lo	= 1.0_dp
		a_x 	= 0.10_dp
		b_x 	= 0.0_dp

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
		write(string_theta,'(f4.2)')  theta

		if ((Progressive_Tax_Switch.eqv..false.).and.(NSU_Switch.eqv..true.)) then 
			Result_Folder = './NSU_ZS_LT_Results/Theta_'//trim(string_theta)//'/Factor_'//trim(Result_Folder)//'/'
		else if ((Progressive_Tax_Switch.eqv..true.).and.(NSU_Switch.eqv..true.)) then 
			Result_Folder = './NSU_ZS_PT_Results/Theta_'//trim(string_theta)//'/Factor_'//trim(Result_Folder)//'/'
		else if ((Progressive_Tax_Switch.eqv..false.).and.(NSU_Switch.eqv..false.)) then 
			Result_Folder = './SU_ZS_LT_Results/Theta_'//trim(string_theta)//'/Factor_'//trim(Result_Folder)//'/'
		else if ((Progressive_Tax_Switch.eqv..true.).and.(NSU_Switch.eqv..false.)) then 
			Result_Folder = './SU_ZS_PT_Results/Theta_'//trim(string_theta)//'/Factor_'//trim(Result_Folder)//'/'
		end if

		Result_Folder = trim(Result_Folder)//'x_hi_30_zeps_125_z11/' 

		! call execute_command_line( 'mkdir -p ' // trim(Result_Folder) )
		call system( 'mkdir -p ' // trim(Result_Folder) )
		print*, "Results are stored in directory: ", Result_Folder
		print*,'na=',na,'update_period=',update_period
		print*, "NSU_Switch=",NSU_Switch,'sigma=',sigma,'gamma=',gamma,'phi',phi
		print*,'Labor Taxes: tauPl=',tauPl,'psi',psi
		print*, 'Borrowing Constraint: Theta=',theta
		print*, 'm tauchen for zgrid is ',mtauchen_z,'nz=',nz, 'amax=',amax, 'totpop=', totpop
		print*, 'x_hi', x_hi, 'a_x', a_x, 'b_x', b_x


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
		if (Calibration_Switch) then 
			CALL CALIBRATION_TRIALS
		endif 

		! Tax Reform experiment
		if (Tax_Reform) then 
			call Solve_Benchmark(compute_bench,Simul_Switch)
			!call Solve_Experiment(compute_exp,Simul_Switch)

			compute_bench = .false.
		endif 

		! Optimal Tax
		if (Opt_Tax) then 
			folder_aux = Result_Folder
			if (Opt_Tax_KW) then 
				Result_Folder = trim(folder_aux)//'Opt_Tax_K/'
			else 
				Result_Folder = trim(folder_aux)//'Opt_Tax_W/'
			endif
			call system( 'mkdir -p ' // trim(Result_Folder) )

			call Solve_Benchmark(compute_bench,Simul_Switch)
			call Solve_Opt_Tax(Opt_Tax_KW,Simul_Switch)
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

	! Solve for the model and compute stats
	print*,"	Initializing program"
		CALL INITIALIZE
	if (compute_bench) then
		print*,"	Computing equilibrium distribution"
		CALL FIND_DBN_EQ
		print*,"	Computing government spending"
		CALL GOVNT_BUDGET
		print*,"	Computing Value Function"
		! CALL COMPUTE_VALUE_FUNCTION_SPLINE
		CALL COMPUTE_VALUE_FUNCTION_LINEAR 
		print*,"	Computing Firm Value Function"
		CALL Firm_Value
		print*,"	Saving results in text files to be read later"
		CALL Write_Benchmark_Results(compute_bench)
	else
		print*,"	Reading benchmark results from files"
		CALL Write_Benchmark_Results(compute_bench)
	end if 

		print*,"	Computing satitics"
		CALL COMPUTE_STATS
		print*,"	Writing variables"
		CALL WRITE_VARIABLES(1)
		if (((theta.eq.1.0_dp).or.(theta.eq.1.50_dp)).and.(Threshold_Factor.eq.0.00_dp).and.(Simul_Switch)) then 
			print*,"	Simulation"
			CALL SIMULATION(solving_bench)
		endif
		!!CALL SIMULATION_TOP(solving_bench)

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

		write(*,*) "Benchmark variables"
		write(*,*) "GBAR=",GBAR,"EBAR=",EBAR,"NBAR=",NBAR,"QBAR=",QBAR,"P=",P,"wage=",wage,'R=',R

end Subroutine Solve_Benchmark


!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Experiment(compute_exp,Simul_Switch)
	use parameters
	use global 
	use programfunctions
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
		CALL COMPUTE_VALUE_FUNCTION_LINEAR
		CALL Firm_Value

	endif 
	
	CALL Write_Experimental_Results(compute_exp)

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

	! Compute moments
	CALL COMPUTE_STATS
	
	! Compute welfare gain between economies
	CALL COMPUTE_WELFARE_GAIN

	! Write experimental results in output.txt
	CALL WRITE_VARIABLES(0)
	! if (((theta.eq.1.0_dp).or.(theta.eq.1.50_dp)).and.(Threshold_Factor.eq.0.0_dp).and.(Simul_Switch)) then 
	!  	print*,"	Experiment Simulation"
	! 	CALL SIMULATION(solving_bench)
	! endif


	print*,'---------------------------'
	print*,''
	print*,'Output Gain Prct=', 100.0_DP*(Y_exp/Y_bench-1.0) 
	print*,''
	print*,'---------------------------'

	print*," "
	print*,"Wealth_factor=",Wealth_factor
	print*," "

end Subroutine Solve_Experiment

!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Opt_Tax(Opt_Tax_KW,Simul_Switch)
	use parameters 
	use global
	use Opt_Tax_Parameters
	use Opt_Tax_Functions
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: Opt_Tax_KW,Simul_Switch

	SSC_Payments_bench = SSC_Payments

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
    	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_k.txt', STATUS='replace')
	    DO tauindx=0,40
            tauK        = real(tauindx,8)/100_DP
            brentvaluet = - EQ_WELFARE_GIVEN_TauK(tauK)

            ! Compute moments
			CALL COMPUTE_STATS

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
		      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60
	    ENDDO 

	    CLOSE (unit=77) 


	    OPEN (UNIT=77, FILE=trim(Result_Folder)//'stat_opt_tau_k.txt', STATUS='replace')

		tauK = OPT_tauK
		psi  = OPT_psi
		call Find_Opt_Tax(Opt_Tax_KW,Opt_TauK,Opt_TauK-0.01_dp,Opt_TauK+0.01_dp) 

		tauK     = OPT_tauK
		OPT_psi  = psi

		! Compute moments
		CALL COMPUTE_STATS

		print*, "Optimal tau_K=", tauK, "Optimal psi=", psi

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
		      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60


		CLOSE (UNIT=77)

	else
    	! Set Y_a_threshold
			call Find_TauW_Threshold(DBN_bench,W_bench)  
			Y_a_threshold = Threshold_Factor*Ebar_bench !0.75_dp
			Wealth_factor = Y_a_threshold/W_bench
    	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Stats_by_tau_w.txt', STATUS='replace')
	    DO tauindx=0,40
            tauw_at     = real(tauindx,8)/1000_DP
            brentvaluet = - EQ_WELFARE_GIVEN_TauW(tauW_at)

            ! Compute moments
			CALL COMPUTE_STATS

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
			!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)), &
		      &100*( (sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)) /&
		      &sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))) &
		      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
		      &100*( (sum(ValueFunction(:,:,:,:,:,:)*DBN1(:,:,:,:,:,:))/sum(DBN1(:,:,:,:,:,:)) /&
		      &sum(ValueFunction_bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:,:))) &
		      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
		      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60
	    ENDDO 

	    CLOSE (unit=77)

	    OPEN (UNIT=77, FILE=trim(Result_Folder)//'stat_opt_tau_w.txt', STATUS='replace')

		tauW_at = OPT_tauW
		psi     = OPT_psi
		call Find_Opt_Tax(Opt_Tax_KW,Opt_TauW,Opt_TauW-0.001_dp,Opt_TauW+0.001_dp)

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
			!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)), &
		      &100*( (sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:)) /&
		      &sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))) &
		      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
		      &100*( (sum(ValueFunction(:,:,:,:,:,:)*DBN1(:,:,:,:,:,:))/sum(DBN1(:,:,:,:,:,:)) /&
		      &sum(ValueFunction_bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:,:))) &
		      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
		      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60


		CLOSE (UNIT=77)
	endif 

	CALL FIND_DBN_EQ
	CALL GOVNT_BUDGET

	! Compute value function and store policy functions, value function and distribution in file
	CALL COMPUTE_VALUE_FUNCTION_SPLINE 
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

end Subroutine Solve_Opt_Tax

!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE CALIBRATION_TRIALS
	use parameters
	use global 
	use programfunctions
	use Toolbox
	use omp_lib

	IMPLICIT NONE
	real(DP)::  betaL, betaH, sigmalambL,sigmalambH, sigmazL, sigmazH, rhozL, rhozH, gammaL, gammaH, muzL, muzH
	integer :: parindx1,  parindx2, parindx3, parindx4, parindx5, parindx6, nbeta, ngamma, nsigmalambda, nrhoz, nsigmaz, nmuz
	real(DP), dimension(6):: paramsL, paramsH

	print*,'SOLVING CALIBRATION'
	print*,'-----------------------------------------------------------------------'



	! CALIBRATIONS WITH POSITIVE DEPRECIATION

	ParamsL = [0.973, 0.0,  0.50, 0.45,  0.34,  0.4494]      
	ParamsH = [0.977, 0.0,  0.50, 0.64,  0.34,  0.4494]     ! dep. rat=0.04, mu=0.9 calibration, 7z,m3, vartheta1.5


	ParamsL= [0.935, 0.0,  0.50, 0.46,  0.34,  0.4494]       
	ParamsH= [0.937, 0.0,  0.50, 0.48,  0.34,  0.4494]       ! mu=0.95, lower sigma_z, vartheta1.5

	ParamsL = [0.968, 0.0,  0.50, 0.46,  0.34,  0.4494]      
	ParamsH = [0.979, 0.0,  0.50, 0.49,  0.34,  0.4494]      ! dep. rate=0.04, mu=0.95, calibration, 7z, m=3, vartheta=1.5

	ParamsL = [0.968, 0.0,  0.50, 0.41,  0.34,  0.4494]      
	ParamsH = [0.979, 0.0,  0.50, 0.43,  0.34,  0.4494]     ! dep. rat=0.04, mu=0.95, calibration, 9z,m4, vartheta1.5

	ParamsL = [0.973, 0.0,  0.50, 0.45,  0.34,  0.4494]      
	ParamsH = [0.975, 0.0,  0.50, 0.64,  0.34,  0.4494]     ! dep. rat=0.04, mu=0.9 calibration, 9 z, m=4, vartheta1.5

	ParamsL = [0.964, 0.0,  0.50, 0.38,  0.34,  0.4494]      
	ParamsH = [0.966, 0.0,  0.50, 0.42,  0.34,  0.4494]      ! dep. rate=0.04, mu=0.98, calibration, 7z, m=3, vartheta=1.5

	ParamsL = [0.967, 0.0,  0.50, 0.40,  0.34,  0.4494]      
	ParamsH = [0.968, 0.0,  0.50, 0.42,  0.34,  0.4494]     ! dep. rate=0.04, mu=0.95, calibration, 11 z, m = 5, vartheta1.5

	ParamsL = [0.968, 0.0,  0.50, 0.42,  0.34,  0.4494]      
	ParamsH = [0.969, 0.0,  0.50, 0.48,  0.34,  0.4494]     ! dep. rate=0.04, mu=0.93, calibration, 11 z, m = 5, vartheta1.5

	ParamsL= [0.935, 0.0,  0.50, 0.45,  0.34,  0.4494]     ! mu=0.95 calibration, targetting 0.34, 0.69, vartheta1.5
	ParamsH= [0.94, 0.0,  0.50, 0.5,  0.34,  0.4494]       ! mu=0.95 calibration, targetting 0.34, 0.69, vartheta1.5

	paramsL=  [0.969,  0.00, 0.50,  0.46,  0.34, 0.4494]       ! zgrid 19, m5, dep004, mu=093
	paramsH=  [0.969,  0.00, 0.50,  0.48,  0.34, 0.4494]       ! zgrid 19, m5, dep004, mu=093

	paramsL=  [0.9675,  0.00, 0.50,  0.41,  0.34, 0.4494]       ! zgrid 19, 17, m5, dep004, mu=095
	paramsH=  [0.9675,  0.00, 0.50,  0.43,  0.34, 0.4494]       ! zgrid 19, 17, m5, dep004, mu=095

	ParamsL = [0.95,  0.0,  0.50, 0.38,  0.29,  0.4494] ! zgrid 11, m5, alpha=0.4, dep005, mu=090, K/Y=3, Top1PVa=0.36 
	ParamsH = [0.955, 0.0,  0.50, 0.40,  0.29,  0.4494] ! phi_B=1

	ParamsL = [0.945, 0.0,  0.50, 0.38, 0.29, 0.4494] ! zgrid 11, m5, alpha=0.4, dep005, mu=090, K/Y=3, Top1PVa=0.36 
	ParamsH = [0.95,  0.0,  0.50, 0.40, 0.29, 0.4494] ! phi_B=1.5


	! CALIBRATION STARTS

	betaL=paramsL(1)
	betaH=paramsH(1)

	muzL = paramsL(2)
	muzH = paramsH(2)

	rhozL = paramsL(3)
	rhozH = paramsH(3)

	sigmazL = paramsL(4)
	sigmazH = paramsH(4)


	sigmalambL= paramsL(5)
	sigmalambH= paramsH(5)

	gammaL  = paramsL(6)
	gammaH  = paramsH(6)

	nbeta =2
	nmuz  =1
	nsigmaz=3
	nrhoz=1
	nsigmalambda=1
	ngamma=1

	Min_SSE_Moments=1000.0_DP

	DO parindx3=1,nsigmalambda
	DO parindx2=1,ngamma
	DO parindx4=1,nrhoz
	DO parindx1=1,nbeta
	DO parindx6=1,nmuz
	DO parindx5=1,nsigmaz

	    beta  = betaL  + real(parindx1-1,8) *(betaH-betaL)/max(real(nbeta-1,8),1.0_DP)
	    gamma = gammaL + real(parindx2-1,8) *(gammaH-gammaL)/max(real(ngamma-1,8),1.0_DP)
	    sigma_lambda_eps = sigmalambL + real(parindx3-1,8)*(sigmalambH -sigmalambL) / max(real(nsigmalambda-1,8),1.0_DP)
	    rho_z= rhozL   +  real(parindx4-1,8)*(rhozH-rhozL) / max(real(nrhoz-1,8),1.0_DP)
	    sigma_z_eps = sigmazL +  real(parindx5-1,8)*(sigmazH-sigmazL) / max(real(nsigmaz-1,8),1.0_DP)
	    mu_z = muzL +  real(parindx6-1,8)*(muzH-muzL) / max(real(nmuz-1,8),1.0_DP)

	    print*,'parameters=',beta, mu_z, rho_z,sigma_z_eps,sigma_lambda_eps, gamma
	        
	    solving_bench=1
	    CALL INITIALIZE
	    CALL FIND_DBN_EQ
	    CALL Firm_Value
	    CALL COMPUTE_MOMENTS
	    SSE_Moments = (1.0-Wealth_Output/3.0_DP)**2.0_DP  + (1.0_DP-FW_top_x_share(4)/0.357_DP)**2.0_DP  & !+ (FW_top_x_share(3)-0.75_DP)**2.0_DP &
	                   & + (1.0_DP-Std_Log_Earnings_25_60 / 0.8_DP)**2.0_DP + (1.0_DP-meanhours_25_60/0.4_DP)**2.0_DP !&
	                   !& + (1.0_DP-MeanReturn/0.069_DP)**2.0_DP
	                   
	!    SSE_Moments = (1.0-Wealth_Output/3.0_DP)**2.0_DP + (1.0_DP-MeanReturn/0.069_DP)**2.0_DP

	!    print*,'parameters=',beta, mu_z, rho_z,sigma_z_eps,sigma_lambda_eps, gamma,'SSE_Moments =',SSE_Moments

	    IF (SSE_Moments .lt. Min_SSE_Moments ) THEN
	        Min_SSE_Moments =SSE_Moments
	        params= [ beta, mu_z, rho_z, sigma_z_eps, sigma_lambda_eps, gamma ]
	        Min_Moments = [  Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, MeanReturn  ]
	    ENDIF
	    !CALL WRITE_TO_FILE
	    
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	print*, params
	print*, Min_Moments

	beta=params(1)
	mu_z=params(2)
	rho_z=params(3)
	sigma_z_eps =params(4)
	sigma_lambda_eps = params(5)
	gamma= params(6)

END SUBROUTINE CALIBRATION_TRIALS

!====================================================================


