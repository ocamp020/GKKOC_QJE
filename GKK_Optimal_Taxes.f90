

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

	IMPLICIT NONE
	! Variables to measure running time
		REAL(DP) :: start_time, finish_time
	! Compute benchmark or load results
		INTEGER  :: read_write_bench

	! Set type of optimal taxe 1->TauK 0->TauW
		opt_tax_switch = 0	

	! Set Parameters 
		Params =[ 0.9436, 0.00, 0.50, 0.70444445, 0.34, 0.4494 ] ! tauL=0.224, tauC=0.075 calibration
		beta   = params(1)
		mu_z   = params(2) ! this is just shifting the z grids. it is zero now.
		rho_z  = params(3) 
		sigma_z_eps      =params(4)
		sigma_lambda_eps = params(5)
		gamma  = params(6)
		sigma  = 4.0_dp

	! Taxes
	! Wealth tax: minimum wealth tax to consider and increments for balancing budget
		tauWmin_bt=0.00_DP
		tauWinc_bt=0.000_DP ! Minimum tax below threshold and increments
		tauWmin_at=0.012_DP
		tauWinc_at=0.002_DP ! Minimum tax above threshold and increments
		Threshold_Factor = 0.00_dp 
	! Consumption tax
		tauC=0.075_DP
	! Set Labor Tax Regime
		!tauPL=0.185_DP
		!psi=0.77_DP  
 		tauPL=0.0_DP
 		psi=0.776_DP  	

	! Resutls Folder
		write(Result_Folder,'(f4.2)') Threshold_Factor

		if ((TauPL.eq.0.0_dp).and.(sigma.ne.1.0_dp)) then 
			Result_Folder = './NSU_LT_Results/Factor_'//trim(Result_Folder)//'/'
		else if ((TauPL.ne.0.0_dp).and.(sigma.ne.1.0_dp)) then 
			Result_Folder = './NSU_PT_Results/Factor_'//trim(Result_Folder)//'/'
		else if ((TauPL.eq.0.0_dp).and.(sigma.eq.1.0_dp)) then 
			Result_Folder = './SU_LT_Results/Factor_'//trim(Result_Folder)//'/'
		else if ((TauPL.ne.0.0_dp).and.(sigma.eq.1.0_dp)) then 
			Result_Folder = './SU_PT_Results/Factor_'//trim(Result_Folder)//'/'
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

	! Set parameters to be used in all simulations economy
		OPEN(UNIT=3, FILE=trim(Result_Folder)//'params', STATUS='replace')
		WRITE(unit=3, FMT=*) params
		CLOSE(unit=3)


	! Start timining of  the process
		call cpu_time(start_time) 

	! Set initia lvalues of R, Wage, Ebar to find equilibrium
		! ------- DO NOT REMOVE THE LINES BELOW

		rr    =  4.906133597851297E-002 
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

	! Solve for the model and compute stats
	read_write_bench = 1
	print*,"	Initializing program"
		CALL INITIALIZE
	if (read_write_bench.eq.0) then
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
		rr_bench    = rr
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
		write(*,*) "GBAR=",GBAR,"EBAR=",EBAR,"NBAR=",NBAR,"QBAR=",QBAR,"rr=",rr,"wage=",wage

	!====================================================================================================
	PRINT*,''
	Print*,'--------------- SOLVING OPTIMAL TAXES -----------------'
	PRINT*,''
	
	! Experiment economy
		solving_bench=0

	! Set Y_a_threshold
		call Find_TauW_Threshold(DBN_bench,W_bench)  
		Y_a_threshold = Threshold_Factor*Ebar_bench !0.75_dp
		Wealth_factor = Y_a_threshold/W_bench
	
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
	
	! AGgregate variable in experimental economy
		GBAR_exp  = GBAR
		QBAR_exp  = QBAR 
		NBAR_exp  = NBAR  
		Y_exp 	  = YBAR
		Ebar_exp  = EBAR
		rr_exp    = rr
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

