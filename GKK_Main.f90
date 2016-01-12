

! This program solves for the wealth taxation model of 
! Guvenen, Kuruscu and Kambourov

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
		logical  :: compute_bench, compute_exp, Opt_Tax, Opt_Tax_KW, Tax_Reform
	! Auxiliary variable for writing file
		character(4)   :: string_theta
		character(100) :: folder_aux

	! Capital Market
		theta = 2.50_dp
	! Threshold 
		Threshold_Factor = 4.00_dp 

	! Switch for solving benchmark or just reading resutls
		! If compute_bench==.true. then just read resutls
		! If compute_bench==.false. then solve for benchmark and store results
		Tax_Reform    = .false.
			compute_bench = .false.
			compute_exp   = .false.
		Opt_Tax       = .true.
			Opt_Tax_KW    = .false. ! true=tau_K false=tau_W


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
		if (theta.eq.1.0_dp) then 
			Params =[0.9415_dp,  0.00_dp,  0.50_dp,  0.65_dp,  0.34_dp,  0.4494_dp] ! tauL=0.224, tauC=0.075 calibration
		else if (theta.eq.1.50_dp) then 
			if (mu.eq.0.85_dp) then 
			Params= [0.945_dp, 0.00_dp, 0.50_dp, 0.7889_dp, 0.34_dp, 0.4494_dp] ! mu=0.85 calibration, targetting 0.34, 0.69, vartheta1.5 
			else 
			Params= [0.9412_dp, 0.0_dp, 0.50_dp, 0.640_dp, 0.34_dp, 0.4494_dp] ! mu=0.9 calibration, targetting 0.34, 0.69, vartheta1.5
			endif 
		else if (theta.eq.1.60_dp) then 
			Params= [0.9409_dp, 0.0_dp, 0.50_dp, 0.640_dp, 0.34_dp, 0.4494_dp] ! mu=0.9 calibration, targetting 0.34, 0.69, vartheta1.6
		else if (theta.eq.2.00_dp) then 
			Params= [0.9405_dp, 0.0_dp, 0.50_dp, 0.639_dp, 0.34_dp, 0.4494_dp] ! mu=0.9 calibration, targetting 0.34, 0.69, vartheta2
		else if (theta.eq.2.50_dp) then 
			Params= [0.9400_dp, 0.0_dp, 0.50_dp, 0.639_dp, 0.34_dp, 0.4494_dp] ! mu=0.9 calibration, targetting 0.34, 0.69, vartheta2.5
		else
			print*, "No parameters for this theta, changing to default parameters (theta=1)"
			Params =[0.9436_dp, 0.0_dp, 0.50_dp, 0.70444445_dp, 0.34_dp, 0.4494_dp] ! tauL=0.224, tauC=0.075 calibration
		endif 
		beta   = params(1)
		mu_z   = params(2) ! this is just shifting the z grids. it is zero now.
		rho_z  = params(3) 
		sigma_z_eps      =params(4)
		sigma_lambda_eps = params(5)
		gamma  = params(6)
		
		sigma  = 4.0_dp
		phi    = (1.0_dp-gamma)/gamma

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
		tauWmin_at=0.014_DP
		tauWinc_at=0.002_DP ! Minimum tax above threshold and increments
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
			Result_Folder = './NSU_F_LT_Results/Theta_'//trim(string_theta)//'/Factor_'//trim(Result_Folder)//'/'
		else if ((Progressive_Tax_Switch.eqv..true.).and.(NSU_Switch.eqv..true.)) then 
			Result_Folder = './NSU_F_PT_Results/Theta_'//trim(string_theta)//'/Factor_'//trim(Result_Folder)//'/'
		else if ((Progressive_Tax_Switch.eqv..false.).and.(NSU_Switch.eqv..false.)) then 
			Result_Folder = './SU_F_LT_Results/Theta_'//trim(string_theta)//'/Factor_'//trim(Result_Folder)//'/'
		else if ((Progressive_Tax_Switch.eqv..true.).and.(NSU_Switch.eqv..false.)) then 
			Result_Folder = './SU_F_PT_Results/Theta_'//trim(string_theta)//'/Factor_'//trim(Result_Folder)//'/'
		end if 

		! call execute_command_line( 'mkdir -p ' // trim(Result_Folder) )
		call system( 'mkdir -p ' // trim(Result_Folder) )
		print*, "Results are stored in directory: ", Result_Folder
		print*,'na=',na,'update_period=',update_period
		print*, "NSU_Switch=",NSU_Switch,'sigma=',sigma,'gamma=',gamma,'phi',phi
		print*,'Labor Taxes: tauPl=',tauPl,'psi',psi
		print*, 'Borrowing Constraint: Theta=',theta


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
		! Tax Reform experiment
		if (Tax_Reform) then 
			call Solve_Benchmark(compute_bench)
			call Solve_Experiment(compute_exp)

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

			call Solve_Benchmark(compute_bench)
			call Solve_Opt_Tax(Opt_Tax_KW)
		endif 



	call cpu_time(finish_time)
	print*,'Total time =',finish_time-start_time


END PROGRAM main


!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Benchmark(compute_bench)
	use parameters
	use global 
	use programfunctions
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: compute_bench

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
		CALL COMPUTE_VALUE_FUNCTION_SPLINE 
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
		if (((theta.eq.1.0_dp).or.(theta.eq.1.50_dp)).and.(compute_bench)) then 
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

		write(*,*) "Benchmark variables"
		write(*,*) "GBAR=",GBAR,"EBAR=",EBAR,"NBAR=",NBAR,"QBAR=",QBAR,"P=",P,"wage=",wage,'R=',R

end Subroutine Solve_Benchmark


!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Experiment(compute_exp)
	use parameters
	use global 
	use programfunctions
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: compute_exp

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
		CALL COMPUTE_VALUE_FUNCTION_SPLINE 

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

	! Compute moments
	CALL COMPUTE_STATS
	
	! Compute welfare gain between economies
	CALL COMPUTE_WELFARE_GAIN

	! Write experimental results in output.txt
	CALL WRITE_VARIABLES(0)
	if (((theta.eq.1.0_dp).or.(theta.eq.1.50_dp)).and.(Threshold_Factor.eq.0.0_dp).and.(compute_exp)) then 
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

end Subroutine Solve_Experiment

!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Solve_Opt_Tax(Opt_Tax_KW)
	use parameters 
	use global
	use Opt_Tax_Parameters
	use Opt_Tax_Functions
	use Toolbox
	use omp_lib
	implicit none 
	logical, intent(in) :: Opt_Tax_KW

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
		if ((theta.eq.1.50_dp).and.(Threshold_Factor.eq.0.00_dp)) then 
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
		    	  & 'Av. Util=', sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))
		      
		    WRITE  (UNIT=77, FMT=*) tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), & 
		      &  MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0), &
		      &  wage, sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)),  &
			!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)), &
		      &100*( (sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)) /&
		      &sum(ValueFunction_bench(1,:,:,:,:)*DBN_bench(1,:,:,:,:))/sum(DBN_bench(1,:,:,:,:))) &
		      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
		      &100*( (sum(ValueFunction(:,:,:,:,:)*DBN1(:,:,:,:,:))/sum(DBN1(:,:,:,:,:)) /&
		      &sum(ValueFunction_bench(:,:,:,:,:)*DBN_bench(:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:))) &
		      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
		      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60
	    ENDDO 

	    CLOSE (unit=77) 


	    OPEN (UNIT=77, FILE=trim(Result_Folder)//'stat_opt_tau_k.txt', STATUS='replace')

		tauK = OPT_tauK
		psi  = OPT_psi
		call Find_Opt_Tax(Opt_Tax_KW,Opt_TauK,Opt_TauK-0.01_dp,Opt_TauK+0.01_dp) 

		else 
			OPEN (UNIT=77, FILE=trim(Result_Folder)//'stat_opt_tau_k.txt', STATUS='replace')

			call Find_Opt_Tax(Opt_Tax_KW,Opt_TauK,0.01_dp,0.10_dp) 
		endif 


		tauK = OPT_tauK
		psi  = OPT_psi

		! Compute moments
		CALL COMPUTE_STATS

		print*, "Optimal tau_K=", tauK, "Optimal psi=", psi

		WRITE  (UNIT=77, FMT=*) tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), & 
		      &  MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0), &
		      &  wage, sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)),  &
			!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)), &
		      &100*( (sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)) /&
		      &sum(ValueFunction_bench(1,:,:,:,:)*DBN_bench(1,:,:,:,:))/sum(DBN_bench(1,:,:,:,:))) &
		      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
		      &100*( (sum(ValueFunction(:,:,:,:,:)*DBN1(:,:,:,:,:))/sum(DBN1(:,:,:,:,:)) /&
		      &sum(ValueFunction_bench(:,:,:,:,:)*DBN_bench(:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:))) &
		      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
		      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60


		CLOSE (UNIT=77)

	else
    	! Set Y_a_threshold
			call Find_TauW_Threshold(DBN_bench,W_bench)  
			Y_a_threshold = Threshold_Factor*Ebar_bench !0.75_dp
			Wealth_factor = Y_a_threshold/W_bench
		if ((theta.eq.1.50_dp).and.(Threshold_Factor.eq.0.00_dp)) then 
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
		    	  & 'Av. Util=', sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))
		      
		    WRITE  (UNIT=77, FMT=*) tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), &
		      &  MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0), &
		      &  wage, sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)),  &
			!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)), &
		      &100*( (sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)) /&
		      &sum(ValueFunction_bench(1,:,:,:,:)*DBN_bench(1,:,:,:,:))/sum(DBN_bench(1,:,:,:,:))) &
		      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
		      &100*( (sum(ValueFunction(:,:,:,:,:)*DBN1(:,:,:,:,:))/sum(DBN1(:,:,:,:,:)) /&
		      &sum(ValueFunction_bench(:,:,:,:,:)*DBN_bench(:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:))) &
		      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
		      & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60
	    ENDDO 

	    CLOSE (unit=77)

	    OPEN (UNIT=77, FILE=trim(Result_Folder)//'stat_opt_tau_w.txt', STATUS='replace')

		tauW_at = OPT_tauW
		psi     = OPT_psi
		call Find_Opt_Tax(Opt_Tax_KW,Opt_TauW,Opt_TauW-0.001_dp,Opt_TauW+0.001_dp)

		else 
			OPEN (UNIT=77, FILE=trim(Result_Folder)//'stat_opt_tau_w.txt', STATUS='replace')

			call Find_Opt_Tax(Opt_Tax_KW,Opt_TauW,0.01_dp,0.05_dp)
		endif 

		tauW_at = OPT_tauW
		psi     = OPT_psi

		! Compute moments
		CALL COMPUTE_STATS

		print*, "Optimal tau_W=", tauW_at, "Optimal psi=", psi
		
		WRITE  (UNIT=77, FMT=*) tauK, tauW_at, psi, GBAR_K/(GBAR_bench +SSC_Payments_bench ), & 
		      &  MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0), &
		      &  wage, sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)),  &
			!      & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)), &
		      &100*( (sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)) /&
		      &sum(ValueFunction_bench(1,:,:,:,:)*DBN_bench(1,:,:,:,:))/sum(DBN_bench(1,:,:,:,:))) &
		      &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP ) , &
		      &100*( (sum(ValueFunction(:,:,:,:,:)*DBN1(:,:,:,:,:))/sum(DBN1(:,:,:,:,:)) /&
		      &sum(ValueFunction_bench(:,:,:,:,:)*DBN_bench(:,:,:,:,:))/sum(DBN_bench(:,:,:,:,:))) &
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
	if (((theta.eq.1.50_dp)).and.(Threshold_Factor.eq.0.0_dp)) then 
	 	print*,"	Optimal Tax Simulation"
		CALL SIMULATION(solving_bench)
	endif

end Subroutine Solve_Opt_Tax



