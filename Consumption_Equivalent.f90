Program Consumption_Equivalent
	use parameters
	use global
	use programfunctions
	use Toolbox
	Implicit None
	character(100) :: Bench_Folder, Result_Folder_aux
	character(4)   :: string_theta
	character(100) :: folder_aux
	integer        :: i
	! Compute benchmark or load results
		logical  :: compute_bench, compute_exp, Opt_Tax, Opt_Tax_KW, Tax_Reform, Simul_Switch, Calibration_Switch
	real(dp), dimension(MaxAge,na,nz,nlambda,ne,nx) :: CE_total, CE_c, CE_cl, CE_cd, CE_h, CE_hl, CE_hd
	real(dp), dimension(MaxAge,na,nz,nlambda,ne,nx) :: CE_nb_cl, CE_nb_cd, CE_nb_hl, CE_nb_hd
	real(dp), dimension(MaxAge,na,nz,nlambda,ne,nx) :: Value_aux
	real(dp) :: CE2_nb_total, CE2_nb_c, CE2_nb_cl, CE2_nb_cd, CE2_nb_h, CE2_nb_hl, CE2_nb_hd
	real(dp) :: CE2_pop_total, CE2_pop_c, CE2_pop_cl, CE2_pop_cd, CE2_pop_h, CE2_pop_hl, CE2_pop_hd
	real(dp) :: C_bench, C_exp, H_bench, H_exp, C_NB_bench, C_NB_exp, H_NB_bench, H_NB_exp
	real(dp) :: CE_total_bench, CE_c_bench, CE_cl_bench, CE_cd_bench, CE_h_bench, CE_hl_bench, CE_hd_bench
	real(dp) :: CE_total_NB_bench, CE_c_NB_bench, CE_cl_NB_bench, CE_cd_NB_bench, CE_h_NB_bench, CE_hl_NB_bench, CE_hd_NB_bench
	real(dp) :: CE_total_exp, CE_c_exp, CE_cl_exp, CE_cd_exp, CE_h_exp, CE_hl_exp, CE_hd_exp
	real(dp) :: CE_total_NB_exp, CE_c_NB_exp, CE_cl_NB_exp, CE_cd_NB_exp, CE_h_NB_exp, CE_hl_NB_exp, CE_hd_NB_exp

	print*, 'CE Program starting'
	! Switch for solving benchmark or just reading resutls
		Calibration_Switch = .false.
		! If compute_bench==.true. then just read resutls
		! If compute_bench==.false. then solve for benchmark and store results
		Tax_Reform    = .false.
			compute_bench = .false.
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

	! Capital Market
		theta_folder = 1.50_dp
		do zi=1,nz
		theta(zi)        = 1.00_dp+(2.50_dp-1.00_dp)/(nz-1)*(real(zi,8)-1.0_dp)
		enddo
	! Threshold 
		Threshold_Factor = 0.00_dp 

	! Set Parameters 
		! Present value calibration
		Params =[ 0.9485_dp, 0.00_dp, 0.1_dp, 0.0665_dp, 0.29_dp , 0.470_dp] ! tauL=0.224, tauC=0.075 calibration
		! Book value calibration
		Params =[ 0.9475_dp, 0.00_dp, 0.1_dp, 0.072_dp , 0.305_dp, 0.46_dp ] ! tauL=0.224, tauC=0.075 calibration
		

		beta   = params(1)
		mu_z   = params(2) ! this is just shifting the z grids. it is zero now.
		rho_z  = params(3) 
		sigma_z_eps      = params(4)
		sigma_lambda_eps = params(5)
		gamma  = params(6)
		sigma  = 4.0_dp

		x_hi	= 5.00_dp
		x_lo	= 1.00_dp
		x_0     = 0.00_dp
		a_x 	= 0.10_dp
		b_x 	= 0.00_dp

	! Taxes
	! Wealth tax: minimum wealth tax to consider and increments for balancing budget
		TauW_bt          = 0.00_dp
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
		write(string_theta,'(f4.2)')  theta_folder
		
		Result_Folder = './NSU_ZS_LT_Results/Theta_'//trim(string_theta)//'/Factor_'//trim(Result_Folder)//'/'
		Result_Folder = trim(Result_Folder)//'Model_1.2_bv/' 
		Result_Folder_aux = Result_Folder

		! call execute_command_line( 'mkdir -p ' // trim(Result_Folder) )
		call system( 'mkdir -p ' // trim(Result_Folder) )
		call system( 'mkdir -p ' // trim(Result_Folder)//'CE_Files/' )
		print*, "Results are stored in directory: ", Result_Folder

	! Bench_Folder
		Bench_Folder = trim(Result_Folder)//'Bench_Files/'

	! Initialize program and load functions
		print*,"	Initializing program"
		CALL INITIALIZE

!====================================================================================================
	PRINT*,''
	Print*,'Loading benchmark'
	PRINT*,''

	! Set taxes for benchmark economy
		tauK = 0.25_DP
		tauW_at = 0.00_DP
		Y_a_threshold = 0.00_DP 

	print*,"	Reading benchmark results from files"
		CALL Write_Benchmark_Results(.false.)

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

print*, ' '
print*, ' Test 1'
print*, ' '

		DBN_bench           = DBN1
		print*, ' Test 2.1'
		print*, size(ValueFunction), size(ValueFunction_bench)
		print*, shape(ValueFunction),' ', shape(ValueFunction_bench)
		ValueFunction_bench = ValueFunction
	print*, ' Test 2.2'
		Cons_bench          = Cons           
		print*, ' Test 2.3'
		Hours_bench         = Hours
		print*, ' Test 2.4'
		Aprime_bench        = Aprime 
		print*, ' Test 2.5'
		V_Pr_bench          = V_Pr
		print*, ' Test 2.6'
		V_Pr_nb_bench       = V_Pr_nb 
		print*, ' Test 2.7'


		CALL ComputeLaborUnits(EBAR, wage) 
print*, ' '
print*, ' Test 3'
print*, ' '


! !====================================================================================================
! do i=1,3

! 	if (i.eq.1) then 
! 		! Tax Reform
! 		Result_Folder = Result_Folder 
! 	elseif (i.eq.2) then 
! 		! Optimal tax: Capital
! 		Result_Folder = trim(Result_Folder_aux)//'Opt_Tax_K/' 
! 	elseif (i.eq.3) then 
! 		! Optimal tax: Wealth
! 		Result_Folder = trim(Result_Folder_aux)//'Opt_Tax_W/' 
! 	endif 

! !====================================================================================================
! 	PRINT*,''
! 	Print*,'Loading experiment'
! 	PRINT*,''

! 	print*,"	Reading benchmark results from files"
! 		CALL Write_Experimental_Results(.false.)
! 	! Aggregate variable in experimental economy
! 		GBAR_exp  = GBAR
! 		QBAR_exp  = QBAR 
! 		NBAR_exp  = NBAR  
! 		Y_exp 	  = YBAR
! 		Ebar_exp  = EBAR
! 		P_exp     = P
! 		R_exp	  = R
! 		wage_exp  = wage
! 		tauK_exp  = tauK
! 		tauPL_exp = tauPL
! 		psi_exp   = psi
! 		DBN_exp   = DBN1
! 		tauw_bt_exp = tauW_bt
! 		tauw_at_exp = tauW_at
! 		Y_a_threshold_exp = Y_a_threshold

! 		ValueFunction_exp = ValueFunction
! 		Cons_exp          = Cons           
! 		Hours_exp         = Hours
! 		Aprime_exp        = Aprime
! 		V_Pr_exp          = V_Pr 
! 		V_Pr_nb_exp  	  = V_Pr_nb

! !====================================================================================================
! 	PRINT*,''
! 	Print*,'Total Consumption and Hours '
! 	PRINT*,''

! 	C_bench = sum(Cons_bench*DBN_bench)
! 	C_exp   = sum(Cons_exp*DBN_exp)
! 	C_NB_bench = sum(Cons_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
! 	C_NB_exp   = sum(Cons_exp(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:))/sum(DBN_exp(1,:,:,:,:,:))

! 	H_bench = sum(Hours_bench*DBN_bench)
! 	H_exp   = sum(Hours_exp*DBN_exp)
! 	H_NB_bench = sum(Hours_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
! 	H_NB_exp   = sum(Hours_exp(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:))/sum(DBN_exp(1,:,:,:,:,:))

! !====================================================================================================
! 	PRINT*,''
! 	Print*,'Consumption Equivalent - Total'
! 	PRINT*,''

! 	CE_total = 100.0_dp*((ValueFunction_exp/ValueFunction_bench)**(1.0_dp/((1.0_dp-sigma)*gamma)) - 1.0_dp ) ;

! 	CE2_nb_total = 100.0_dp*&
! 				& ((sum(ValueFunction_exp(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:))&
! 				&	/sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:)))&
! 				& **(1.0_dp/((1.0_dp-sigma)*gamma)) - 1.0_dp ) ;

! 	CE2_pop_total = 100.0_dp*((sum(ValueFunction_exp*DBN_exp)/sum(ValueFunction_bench*DBN_bench))&
! 				& **(1.0_dp/((1.0_dp-sigma)*gamma)) - 1.0_dp ) ;

! !====================================================================================================
! 	PRINT*,''
! 	Print*,'Consumption Equivalent - Consumption'
! 	PRINT*,''

! 	Cons    = Cons_exp 
! 	Aprime  = Aprime_exp
! 	Hours   = Hours_bench 

! 	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,Value_aux)

! 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	CE_c  = 100.0_dp*((Value_aux/ValueFunction_bench)**(1.0_dp/((1.0_dp-sigma)*gamma)) - 1.0_dp ) ;

! 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	CE_cl = 100_dp*(C_exp/C_bench - 1.0_dp)

! 	CE_cd = 100.0_dp*((Value_aux/ValueFunction_bench)**(1.0_dp/((1.0_dp-sigma)*gamma)) * C_bench/C_exp - 1.0_dp ) ;

! 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	CE_nb_cl = 100_dp*(C_NB_exp/C_NB_bench - 1.0_dp)

! 	CE_nb_cd = 100.0_dp*((Value_aux/ValueFunction_bench)**(1.0_dp/((1.0_dp-sigma)*gamma)) * C_NB_bench/C_NB_exp - 1.0_dp ) ;

! 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	CE2_nb_c = 100.0_dp*((sum(Value_aux(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/&
! 		&                  sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:)))&
! 		&                  **(1.0_dp/((1.0_dp-sigma)*gamma)) - 1.0_dp ) ;

! 	CE2_nb_cl = 100_dp*(C_NB_exp/C_NB_bench - 1.0_dp)

! 	CE2_nb_cd = 100.0_dp*((sum(Value_aux(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/&
! 		&				   sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:)))&
! 		& 					**(1.0_dp/((1.0_dp-sigma)*gamma)) * C_NB_bench/C_NB_exp - 1.0_dp ) ;

! 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	CE2_pop_c  = 100.0_dp*((sum(Value_aux*DBN_bench)/sum(ValueFunction_bench*DBN_bench))&
! 		& **(1.0_dp/((1.0_dp-sigma)*gamma)) - 1.0_dp ) ;

! 	CE2_pop_cl = 100_dp*(C_exp/C_bench - 1.0_dp)

! 	CE2_pop_cd = 100.0_dp*((sum(Value_aux*DBN_exp)/sum(ValueFunction_bench*DBN_bench))&
! 		& **(1.0_dp/((1.0_dp-sigma)*gamma)) * C_bench/C_exp - 1.0_dp ) ;

! !====================================================================================================
! 	PRINT*,''
! 	Print*,'Consumption Equivalent - Hours'
! 	PRINT*,''


! 	Cons    = Cons_exp
! 	Aprime  = Aprime_exp
! 	Hours   = Hours_bench 
	
! 	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,Value_aux)

! 	Cons    = Cons_bench 
! 	Aprime  = Aprime_bench
! 	Hours   = Hours_exp
! 	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction)

! 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	CE_h  = 100.0_dp*((ValueFunction_exp/Value_aux)**(1.0_dp/((1.0_dp-sigma)*gamma)) - 1.0_dp ) ;

! 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	CE_hl = 100_dp*(((1.0_dp-H_exp)/(1.0_dp-H_bench))**((1.0_dp-gamma)/gamma) - 1.0_dp)

! 	CE_hd = 100.0_dp*((ValueFunction/ValueFunction_bench)**(1.0_dp/((1.0_dp-sigma)*gamma)) &
! 	        &  * ((1.0_dp-H_bench)/(1.0_dp-H_exp))**((1.0_dp-gamma)/gamma) - 1.0_dp ) ;

! 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	CE_nb_hl = 100_dp*(((1.0_dp-H_NB_exp)/(1.0_dp-H_NB_bench))**((1.0_dp-gamma)/gamma) - 1.0_dp)

! 	CE_nb_hd = 100.0_dp*((ValueFunction/ValueFunction_bench)**(1.0_dp/((1.0_dp-sigma)*gamma)) &
! 	        &  * ((1.0_dp-H_NB_bench)/(1.0_dp-H_NB_exp))**((1.0_dp-gamma)/gamma) - 1.0_dp ) ;

! 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	CE2_nb_h  = 100.0_dp*((sum(ValueFunction_exp(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:))/&
! 			&				sum(Value_aux(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:)))&
! 			& **(1.0_dp/((1.0_dp-sigma)*gamma)) - 1.0_dp ) ;

! 	CE2_nb_hl = 100_dp*(((1.0_dp-H_NB_exp)/(1.0_dp-H_NB_bench))**((1.0_dp-gamma)/gamma) - 1.0_dp)

! 	CE2_nb_hd = 100.0_dp*((sum(ValueFunction(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:))/&
! 			& sum(ValueFunction_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:)))&
! 			& **(1.0_dp/((1.0_dp-sigma)*gamma)) &
! 	        &  * ((1.0_dp-H_NB_bench)/(1.0_dp-H_NB_exp))**((1.0_dp-gamma)/gamma) - 1.0_dp ) ;

! 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	CE2_pop_h  = 100.0_dp*((sum(ValueFunction_exp*DBN_exp)/sum(Value_aux*DBN_bench))&
! 			& **(1.0_dp/((1.0_dp-sigma)*gamma)) - 1.0_dp ) ;

! 	CE2_pop_hl = 100_dp*(((1.0_dp-H_exp)/(1.0_dp-H_bench))**((1.0_dp-gamma)/gamma) - 1.0_dp)

! 	CE2_pop_hd = 100.0_dp*((sum(ValueFunction*DBN_exp)/sum(ValueFunction_bench*DBN_bench))&
! 			& **(1.0_dp/((1.0_dp-sigma)*gamma)) &
! 	        &  * ((1.0_dp-H_bench)/(1.0_dp-H_exp))**((1.0_dp-gamma)/gamma) - 1.0_dp ) ;

! !====================================================================================================
! 	PRINT*,''
! 	Print*,'Report Output'
! 	PRINT*,''

! 	CE_total_bench 	= sum(CE_total*DBN_bench )
! 	CE_c_bench     	= sum(CE_c*DBN_bench )
! 	CE_cl_bench    	= sum(CE_cl*DBN_bench )
! 	CE_cd_bench 	= sum(CE_cd*DBN_bench )
! 	CE_h_bench 		= sum(CE_h*DBN_bench )
! 	CE_hl_bench 	= sum(CE_hl*DBN_bench )
! 	CE_hd_bench 	= sum(CE_hd*DBN_bench )

! 	CE_total_NB_bench 	= sum(CE_total(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
! 	CE_c_NB_bench 		= sum(CE_c(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
! 	CE_cl_NB_bench 		= sum(CE_nb_cl(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
! 	CE_cd_NB_bench 		= sum(CE_nb_cd(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
! 	CE_h_NB_bench 		= sum(CE_h(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
! 	CE_hl_NB_bench 		= sum(CE_nb_hl(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
! 	CE_hd_NB_bench 		= sum(CE_nb_hd(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))

! 	CE_total_exp 	= sum(CE_total*DBN_exp )
! 	CE_c_exp 		= sum(CE_c*DBN_exp )
! 	CE_cl_exp 		= sum(CE_cl*DBN_exp )
! 	CE_cd_exp 		= sum(CE_cd*DBN_exp )
! 	CE_h_exp 		= sum(CE_h*DBN_exp )
! 	CE_hl_exp 		= sum(CE_hl*DBN_exp )
! 	CE_hd_exp 		= sum(CE_hd*DBN_exp )

! 	CE_total_NB_exp 	= sum(CE_total(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:))/sum(DBN_exp(1,:,:,:,:,:))
! 	CE_c_NB_exp 		= sum(CE_c(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:))/sum(DBN_exp(1,:,:,:,:,:))
! 	CE_cl_NB_exp 		= sum(CE_nb_cl(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:))/sum(DBN_exp(1,:,:,:,:,:))
! 	CE_cd_NB_exp 		= sum(CE_nb_cd(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:))/sum(DBN_exp(1,:,:,:,:,:))
! 	CE_h_NB_exp 		= sum(CE_h(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:))/sum(DBN_exp(1,:,:,:,:,:))
! 	CE_hl_NB_exp 		= sum(CE_nb_hl(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:))/sum(DBN_exp(1,:,:,:,:,:))
! 	CE_hd_NB_exp 		= sum(CE_nb_hd(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:))/sum(DBN_exp(1,:,:,:,:,:))


! 	if (i.eq.1) then 
! 		! Tax Reform
! 		OPEN  (UNIT=1,  FILE=trim(Result_Folder_aux)//'CE_output.txt'  , STATUS='replace')
! 		WRITE (UNIT=1,  FMT=*) ' '
! 		WRITE (UNIT=1,  FMT=*) 'Consumption Equivalent Welfare - Tax Reform'
! 		WRITE (UNIT=1,  FMT=*) ' '
! 	elseif (i.eq.2) then 
! 		! Optimal tax: Capital
! 		OPEN  (UNIT=1,  FILE=trim(Result_Folder_aux)//'CE_output_otk.txt'  , STATUS='replace')
! 		WRITE (UNIT=1,  FMT=*) ' '
! 		WRITE (UNIT=1,  FMT=*) 'Consumption Equivalent Welfare - Optimal Capital Taxes'
! 		WRITE (UNIT=1,  FMT=*) ' '
! 	elseif (i.eq.3) then 
! 		! Optimal tax: Wealth
! 		OPEN  (UNIT=1,  FILE=trim(Result_Folder_aux)//'CE_output_otw.txt'  , STATUS='replace')
! 		WRITE (UNIT=1,  FMT=*) ' '
! 		WRITE (UNIT=1,  FMT=*) 'Consumption Equivalent Welfare - Optimal Wealth Taxes'
! 		WRITE (UNIT=1,  FMT=*) ' '
! 	endif 

! 	WRITE (UNIT=1,  FMT=*) 'Benchmark - Aggregate'
! 	WRITE (UNIT=1,  FMT=*) 'CE',CE_total_bench, 'CE_c',CE_c_bench,'CE_h',CE_h_bench
! 	WRITE (UNIT=1,  FMT=*) 'CE_c',CE_c_bench,'CE_cl',CE_cl_bench,'CE_cd',CE_cd_bench
! 	WRITE (UNIT=1,  FMT=*) 'CE_h',CE_h_bench,'CE_hl',CE_hl_bench,'CE_hd',CE_hd_bench
! 	WRITE (UNIT=1,  FMT=*)' '

! 	WRITE (UNIT=1,  FMT=*)'Benchmark - NewBorn'
! 	WRITE (UNIT=1,  FMT=*) 'CE',CE_total_NB_bench, 'CE_c',CE_c_NB_bench,'CE_h',CE_h_NB_bench
! 	WRITE (UNIT=1,  FMT=*) 'CE_c',CE_c_NB_bench,'CE_cl',CE_cl_NB_bench,'CE_cd',CE_cd_NB_bench
! 	WRITE (UNIT=1,  FMT=*) 'CE_h',CE_h_NB_bench,'CE_hl',CE_hl_NB_bench,'CE_hd',CE_hd_NB_bench
! 	WRITE (UNIT=1,  FMT=*)' '

! 	WRITE (UNIT=1,  FMT=*)'Experimental - Aggregate'
! 	WRITE (UNIT=1,  FMT=*) 'CE',CE_total_exp, 'CE_c',CE_c_exp,'CE_h',CE_h_exp
! 	WRITE (UNIT=1,  FMT=*) 'CE_c',CE_c_exp,'CE_cl',CE_cl_exp,'CE_cd',CE_cd_exp
! 	WRITE (UNIT=1,  FMT=*) 'CE_h',CE_h_exp,'CE_hl',CE_hl_exp,'CE_hd',CE_hd_exp
! 	WRITE (UNIT=1,  FMT=*)' '

! 	WRITE (UNIT=1,  FMT=*)'Experimental - NewBorn'
! 	WRITE (UNIT=1,  FMT=*) 'CE',CE_total_NB_exp, 'CE_c',CE_c_NB_exp,'CE_h',CE_h_NB_exp
! 	WRITE (UNIT=1,  FMT=*) 'CE_c',CE_c_NB_exp,'CE_cl',CE_cl_NB_exp,'CE_cd',CE_cd_NB_exp
! 	WRITE (UNIT=1,  FMT=*) 'CE_h',CE_h_NB_exp,'CE_hl',CE_hl_NB_exp,'CE_hd',CE_hd_NB_exp
! 	WRITE (UNIT=1,  FMT=*)' '

! 	WRITE (UNIT=1,  FMT=*) 'CE1 '	, 'NB '				, 'Pop '			, 'Test' 	
! 	WRITE (UNIT=1,  FMT=*) 'CE1'	, CE_total_NB_bench , CE_total_bench 	, 100*((1+CE_c_NB_bench/100)*(1+CE_h_NB_bench/100)-1)
! 	WRITE (UNIT=1,  FMT=*) 'CE1_c'	, CE_c_NB_bench 	, CE_c_bench 		, 100*((1+CE_cl_NB_bench/100)*(1+CE_cd_NB_bench/100)-1)		
! 	WRITE (UNIT=1,  FMT=*) 'CE1_cl'	, CE_cl_NB_bench	, CE_cl_bench
! 	WRITE (UNIT=1,  FMT=*) 'CE_cd'	, CE_cd_NB_bench 	, CE_cd_bench
! 	WRITE (UNIT=1,  FMT=*) 'CE1_h'	, CE_h_NB_bench 	, CE_h_bench 		, 100*((1+CE_hl_NB_bench/100)*(1+CE_hd_NB_bench/100)-1)		
! 	WRITE (UNIT=1,  FMT=*) 'CE1_hl'	, CE_hl_NB_bench 	, CE_hl_bench
! 	WRITE (UNIT=1,  FMT=*) 'CE1_hd'	, CE_hd_NB_bench 	, CE_hd_bench 
! 	WRITE (UNIT=1,  FMT=*)' '

! 	WRITE (UNIT=1,  FMT=*) 'CE2 '	, 'NB '			, 'Pop '		, 'Test' 	
! 	WRITE (UNIT=1,  FMT=*) 'CE2'	, CE2_nb_total 	, CE2_pop_total , 100*((1+CE2_nb_c/100)*(1+CE2_nb_h/100)-1)
! 	WRITE (UNIT=1,  FMT=*) 'CE2_c'	, CE2_nb_c 		, CE2_pop_c 	, 100*((1+CE2_nb_cl/100)*(1+CE2_nb_cd/100)-1)		
! 	WRITE (UNIT=1,  FMT=*) 'CE2_cl'	, CE2_nb_cl		, CE2_pop_cl
! 	WRITE (UNIT=1,  FMT=*) 'CE2_cd'	, CE2_nb_cd 	, CE2_pop_cd
! 	WRITE (UNIT=1,  FMT=*) 'CE2_h'	, CE2_nb_h 		, CE2_pop_h 	, 100*((1+CE2_nb_hl/100)*(1+CE2_nb_hd/100)-1)
! 	WRITE (UNIT=1,  FMT=*) 'CE2_hl'	, CE2_nb_hl 	, CE2_pop_hl
! 	WRITE (UNIT=1,  FMT=*) 'CE2_hd'	, CE2_nb_hd 	, CE2_pop_hd 
! 	WRITE (UNIT=1,  FMT=*)' '

! 	CLOSE (unit=1)

! 	! OPEN  (UNIT=1,  FILE=trim(Result_Folder)//'CE_Files/CE_total'  , STATUS='replace')
! 	! OPEN  (UNIT=2,  FILE=trim(Result_Folder)//'CE_Files/CE_c'  , STATUS='replace')
! 	! OPEN  (UNIT=3,  FILE=trim(Result_Folder)//'CE_Files/CE_cl'  , STATUS='replace')
! 	! OPEN  (UNIT=4,  FILE=trim(Result_Folder)//'CE_Files/CE_cd'  , STATUS='replace')
! 	! OPEN  (UNIT=5,  FILE=trim(Result_Folder)//'CE_Files/CE_h'  , STATUS='replace')
! 	! OPEN  (UNIT=6,  FILE=trim(Result_Folder)//'CE_Files/CE_hl'  , STATUS='replace')
! 	! OPEN  (UNIT=7,  FILE=trim(Result_Folder)//'CE_Files/CE_hd'  , STATUS='replace')
	
! 	! WRITE (UNIT=1,  FMT=*) CE_total
! 	! WRITE (UNIT=2,  FMT=*) CE_c
! 	! WRITE (UNIT=3,  FMT=*) CE_cl
! 	! WRITE (UNIT=4,  FMT=*) CE_cd
! 	! WRITE (UNIT=5,  FMT=*) CE_h
! 	! WRITE (UNIT=6,  FMT=*) CE_hl
! 	! WRITE (UNIT=7,  FMT=*) CE_hd
	
! 	! CLOSE (unit=1)
! 	! CLOSE (unit=2)
! 	! CLOSE (unit=3)
! 	! CLOSE (unit=4)
! 	! CLOSE (unit=5)
! 	! CLOSE (unit=6)
! 	! CLOSE (unit=7)

! 	if (i.eq.1) then 
! 		! Tax Reform
! 		print*, 'Consumption Equivalent Welfare - Tax Reform'
! 	elseif (i.eq.2) then 
! 		! Optimal tax: Capital
! 		print*, 'Consumption Equivalent Welfare - Optimal Capital Taxes'
! 	elseif (i.eq.3) then 
! 		! Optimal tax: Wealth
! 		print*, 'Consumption Equivalent Welfare - Optimal Wealth Taxes'
! 	endif 

! 	print*,' '
! 	print*, 'CE1 '	, 'NB '				, 'Pop '			, 'Test' 	
! 	print*, 'CE1'	, CE_total_NB_bench , CE_total_bench 	, 100*((1+CE_c_NB_bench/100)*(1+CE_h_NB_bench/100)-1)
! 	print*, 'CE1_c'	, CE_c_NB_bench 	, CE_c_bench 		, 100*((1+CE_cl_NB_bench/100)*(1+CE_cd_NB_bench/100)-1)		
! 	print*, 'CE1_cl', CE_cl_NB_bench	, CE_cl_bench
! 	print*, 'CE_cd'	, CE_cd_NB_bench 	, CE_cd_bench
! 	print*, 'CE1_h'	, CE_h_NB_bench 	, CE_h_bench 		, 100*((1+CE_hl_NB_bench/100)*(1+CE_hd_NB_bench/100)-1)		
! 	print*, 'CE1_hl', CE_hl_NB_bench 	, CE_hl_bench
! 	print*, 'CE1_hd', CE_hd_NB_bench 	, CE_hd_bench 
! 	print*,' '
! 	print*,  'CE2 '		, 'NB '			, 'Pop '		, 'Test '
! 	print*,  'CE2'		, CE2_nb_total 	, CE2_pop_total , 100*((1+CE2_nb_c/100)*(1+CE2_nb_h/100)-1)
! 	print*,  'CE2_c'	, CE2_nb_c 		, CE2_pop_c 	, 100*((1+CE2_nb_cl/100)*(1+CE2_nb_cd/100)-1)
! 	print*,  'CE2_cl'	, CE2_nb_cl		, CE2_pop_cl
! 	print*,  'CE2_cd'	, CE2_nb_cd 	, CE2_pop_cd
! 	print*,  'CE2_h'	, CE2_nb_h 		, CE2_pop_h 	, 100*((1+CE2_nb_hl/100)*(1+CE2_nb_hd/100)-1)
! 	print*,  'CE2_hl'	, CE2_nb_hl 	, CE2_pop_hl
! 	print*,  'CE2_hd'	, CE2_nb_hd 	, CE2_pop_hd 
! 	print*, ' '
	
! enddo 


end Program Consumption_Equivalent