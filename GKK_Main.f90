

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

	IMPLICIT NONE
	! Variables to measure running time
		REAL(DP) :: start_time, finish_time
	! Compute benchmark or load results
		logical  :: read_write_bench
	! Auxiliary variable for writing file
		character(4) :: string_theta

	! Switch for solving benchmark or just reading resutls
		! If read_write_bench==.true. then just read resutls
		! If read_write_bench==.false. then solve for benchmark and store results
		read_write_bench = .false.

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
		Params =[ 0.9436, 0.00, 0.50, 0.70444445, 0.34, 0.4494 ] ! tauL=0.224, tauC=0.075 calibration
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
		Threshold_Factor = 0.00_dp 
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

	print*," "
	print*,"Wealth_factor=",Wealth_factor
	print*," "

	call cpu_time(finish_time)
	print*,'Total time =',finish_time-start_time


END PROGRAM main


!========================================================================================
!========================================================================================
!========================================================================================

