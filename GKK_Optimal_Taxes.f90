

! This program solves for the wealth taxation model of 
! Guvenen, Kuruscu and Kambourov
! It then computes welfare over a grid of taxes to determine the optimal ones.

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

PROGRAM Optimal_Taxes
	use parameters 
	use global
	use Opt_Tax_Parameters
	use Opt_Tax_Functions
	use Toolbox
	use omp_lib

	IMPLICIT NONE
	! Variables to measure running time
		REAL(DP) :: start_time, finish_time
	! Compute benchmark or load results
		logical  :: read_write_bench
	! Auxiliary variable for writing file
		character(4) :: string_theta

	! Set type of optimal taxe 1->TauK 0->TauW
		opt_tax_switch = 1

	! Switch for solving benchmark or just reading resutls
		! If read_write_bench==.true. then just read resutls
		! If read_write_bench==.false. then solve for benchmark and store results
		read_write_bench = .true.	

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
			Params =[0.9436_dp, 0.0_dp, 0.50_dp, 0.70444445_dp, 0.34_dp, 0.4494_dp] ! tauL=0.224, tauC=0.075 calibration
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
		Threshold_Factor = 2.00_dp 
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

		if (opt_tax_switch.eq.1) then 
			Result_Folder = trim(Result_Folder)//'Opt_Tax_K/'
		else 
			Result_Folder = trim(Result_Folder)//'Opt_Tax_W/'
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
	if (read_write_bench.eqv..false.) then
		print*,"	Computing equilibrium distribution"
		CALL FIND_DBN_EQ
		print*,"	Computing government spending"
		CALL GOVNT_BUDGET
		print*,"	Computing Value Function"
		CALL COMPUTE_VALUE_FUNCTION_SPLINE 
		print*,"	Saving results in text files to be read later"
		CALL Write_Benchmark_Results(read_write_bench)
	else
		print*,"	Reading benchmark results from files"
		CALL Write_Benchmark_Results(read_write_bench)
	end if 

		print*,"	Computing satitics"
		CALL COMPUTE_STATS
		print*,"	Writing variables"
		CALL WRITE_VARIABLES(1)

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
	If ( opt_tax_switch .eq. 1 ) then
		
		OPEN (UNIT=77, FILE=trim(Result_Folder)//'stat_all_tau_k', STATUS='replace')
	    
	    !brentvaluet = brent(0.00_DP, 0.1_DP , 0.4_DP, EQ_WELFARE_GIVEN_TauK, brent_tol, Opt_TauK)  
	    call Find_Opt_Tax(opt_tax_switch,Opt_TauK)
	    tauK = Opt_TauK

	    brentvaluet = - EQ_WELFARE_GIVEN_TauK(tauK)

	    ! Print results
            print*, ' '
            print*, 'Values optimal capital taxes'
            print*, tauK, tauPL, psi, GBAR_K, MeanWealth, QBAR, NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0) , &
              & CE_NEWBORN, sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))

            WRITE  (UNIT=77, FMT=*) tauK, tauW_at, tauPL, psi, GBAR_K, MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0) , &
              & CE_NEWBORN, sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)), &
              & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60	    


		! 	    DO tauindx=0,50

		!             tauK=real(tauindx,8)/100_DP
		!             brentvaluet = - EQ_WELFARE_GIVEN_TauK(tauK)

		!             if (brentvaluet .gt. maxbrentvaluet) then
		!                 maxbrentvaluet = brentvaluet
		!                 Opt_TauK=tauK
		!             endif

		!             ! Print results
		!             print*, ' '
		!             print*, 'Iteration',tauindx
		!             print*, tauK, tauPL, psi, GBAR_K, MeanWealth, QBAR, NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0) , &
		!               & CE_NEWBORN, sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))

		!             WRITE  (UNIT=77, FMT=*) tauK, tauW_at, tauPL, psi, GBAR_K, MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0) , &
		!               & CE_NEWBORN, sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)), &
		!               & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60	    
		! 	    ENDDO                
	else

		OPEN (UNIT=77, FILE=trim(Result_Folder)//'stat_all_tau_w', STATUS='replace')

		! Set Y_a_threshold
		call Find_TauW_Threshold(DBN_bench,W_bench)  
		Y_a_threshold = Threshold_Factor*Ebar_bench !0.75_dp
		Wealth_factor = Y_a_threshold/W_bench

	    !brentvaluet = brent(0.00_DP, 0.016_DP , 0.05_DP, EQ_WELFARE_GIVEN_TauW, brent_tol, Opt_TauW)
	    call Find_Opt_Tax(opt_tax_switch,Opt_TauW)
	    tauW_at = Opt_TauW

	    brentvaluet = -EQ_WELFARE_GIVEN_TauW(tauW_at)

	    ! Print results
            print*, ' '
            print*, 'Values optimal wealth taxes'
            print*, tauW_at, tauPL, psi, GBAR_K, MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0) , &
              & CE_NEWBORN, sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))

            WRITE  (UNIT=77, FMT=*) tauK, tauW_at, tauPL, psi, GBAR_K, MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0) , &
              & CE_NEWBORN, sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)), &
              & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60

		! 	    DO tauindx=0,50
		!             tauW_at=real(tauindx,8)/1000_DP
		!             brentvaluet = -EQ_WELFARE_GIVEN_TauW(tauW_at)
		            
		!             if (brentvaluet .gt. maxbrentvaluet) then
		!                 maxbrentvaluet = brentvaluet
		!                 Opt_TauW=tauW_at
		!             endif

		!             ! Print results
		!             print*, ' '
		!             print*, 'Iteration',tauindx
		!             print*, tauW_at, tauPL, psi, GBAR_K, MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0) , &
		!               & CE_NEWBORN, sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))

		!             WRITE  (UNIT=77, FMT=*) tauK, tauW_at, tauPL, psi, GBAR_K, MeanWealth, QBAR,NBAR, YBAR, 100.0_DP*(Y_exp/Y_bench-1.0) , &
		!               & CE_NEWBORN, sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:)), &
		!               & Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60
		! 	    ENDDO                
	endif
	close (unit=77)

	tauK    = Opt_TauK
	tauW_at = Opt_TauW

	CALL FIND_DBN_EQ
	CALL GOVNT_BUDGET

	! Compute value function and store policy functions, value function and distribution in file
	CALL COMPUTE_VALUE_FUNCTION_SPLINE 
	CALL Write_Experimental_Results()
	
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


	print*,'---------------------------'
	print*,''
	print*,'Output Gain Prct=', 100.0_DP*(Y_exp/Y_bench-1.0) 
	print*,''
	print*,'---------------------------'

	call cpu_time(finish_time)
	print*,'Total time =',finish_time-start_time


END PROGRAM Optimal_Taxes


!========================================================================================
!========================================================================================
!========================================================================================

