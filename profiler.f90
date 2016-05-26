

!========================================================================================
!========================================================================================
!========================================================================================

PROGRAM profiler 
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
	! Counters
		integer :: zi 

	! Capital Market
		theta_folder = 1.50_dp
		do zi=1,nz
		theta(zi)        = 1.00_dp+(2.50_dp-1.00_dp)/(nz-1)*(real(zi,8)-1.0_dp)
		enddo

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

		! New Exponential Shock that only affects high Z and includes Z-varying theta 2.5
		! beta 		sigmaz 		x_hi 	rho_z 	gamma
		! 0.9485_dp 0.0665_dp  	5.00_dp 0.1_dp 	0.470_dp

		beta   	= 0.9485_dp ! params(1) !
		mu_z   	= params(2) ! this is just shifting the z grids. it is zero now.
		rho_z  	= 0.1_dp ! params(3) 
		sigma_z_eps      =  0.0665_dp ! params(4) ! 0.01_dp ! ! 
		sigma_lambda_eps = params(5)
		gamma  	=  0.470_dp !  0.465_dp ! params(6) ! 
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
		tauWmin=0.010_DP
		tauWinc=0.005_DP ! Minimum tax above threshold and increments
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
		write(string_theta,'(f4.2)')  theta(nz)

		if ((Progressive_Tax_Switch.eqv..false.).and.(NSU_Switch.eqv..true.)) then 
			Result_Folder = './NSU_ZS_LT_Results/Theta_'//trim(string_theta)//'/No_Threshold/'
		else if ((Progressive_Tax_Switch.eqv..true.).and.(NSU_Switch.eqv..true.)) then 
			Result_Folder = './NSU_ZS_PT_Results/Theta_'//trim(string_theta)//'/No_Threshold/'
		else if ((Progressive_Tax_Switch.eqv..false.).and.(NSU_Switch.eqv..false.)) then 
			Result_Folder = './SU_ZS_LT_Results/Theta_'//trim(string_theta)//'/No_Threshold/'
		else if ((Progressive_Tax_Switch.eqv..true.).and.(NSU_Switch.eqv..false.)) then 
			Result_Folder = './SU_ZS_PT_Results/Theta_'//trim(string_theta)//'/No_Threshold/'
		end if

		Result_Folder = trim(Result_Folder)//'Model_1.3/' 

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


			CALL INITIALIZE
		! Compute labor units 
			CALL ComputeLaborUnits(Ebar, wage) 
		! Solve for policy and value functions 
			CALL EGM_RETIREMENT_WORKING_PERIOD 



	call cpu_time(finish_time)
	print*,'Total time =',finish_time-start_time


END PROGRAM profiler