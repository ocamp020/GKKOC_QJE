Program Consumption_Equivalent
	use parameters
	use global
	use programfunctions
	use Toolbox
	Implicit None
	character(100) :: Bench_Folder
	real(dp), dimension(MaxAge,na,nz,nlambda,ne) :: Cons_bench, Hours_bench, Aprime_bench, Value_bench
	real(dp), dimension(MaxAge,na,nz,nlambda,ne) :: Cons_exp, Hours_exp, Aprime_exp, Value_exp
	real(dp), dimension(MaxAge,na,nz,nlambda,ne) :: CE_total, CE_c, CE_cl, CE_cd, CE_h, CE_hl, CE_hd
	real(dp), dimension(MaxAge,na,nz,nlambda,ne) :: Value_aux
	real(dp) :: C_bench, C_exp, H_bench, H_exp
	real(dp) :: CE_total_bench, CE_c_bench, CE_cl_bench, CE_cd_bench, CE_h_bench, CE_hl_bench, CE_hd_bench
	real(dp) :: CE_total_NB_bench, CE_c_NB_bench, CE_cl_NB_bench, CE_cd_NB_bench, CE_h_NB_bench, CE_hl_NB_bench, CE_hd_NB_bench
	real(dp) :: CE_total_exp, CE_c_exp, CE_cl_exp, CE_cd_exp, CE_h_exp, CE_hl_exp, CE_hd_exp
	real(dp) :: CE_total_NB_exp, CE_c_NB_exp, CE_cl_NB_exp, CE_cd_NB_exp, CE_h_NB_exp, CE_hl_NB_exp, CE_hd_NB_exp


	! Set Parameters 
		Params =[ 0.9436, 0.00, 0.50, 0.70444445, 0.34, 0.4494 ] ! tauL=0.224, tauC=0.075 calibration
		beta   = params(1)
		mu_z   = params(2) ! this is just shifting the z grids. it is zero now.
		rho_z  = params(3) 
		sigma_z_eps      = params(4)
		sigma_lambda_eps = params(5)
		gamma  = params(6)
		sigma  = 4.0_dp

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

		if ((TauPL.eq.0.0_dp).and.(sigma.ne.1.0_dp)) then 
			Result_Folder = './NSU_LT_Results/Factor_'//trim(Result_Folder)//'/'
		else if ((TauPL.ne.0.0_dp).and.(sigma.ne.1.0_dp)) then 
			Result_Folder = './NSU_PT_Results/Factor_'//trim(Result_Folder)//'/'
		else if ((TauPL.eq.0.0_dp).and.(sigma.eq.1.0_dp)) then 
			Result_Folder = './SU_LT_Results/Factor_'//trim(Result_Folder)//'/'
		else if ((TauPL.ne.0.0_dp).and.(sigma.eq.1.0_dp)) then 
			Result_Folder = './SU_PT_Results/Factor_'//trim(Result_Folder)//'/'
		end if 
		print*, "Results are stored in directory: ", Result_Folder

	! Bench_Folder
		if ((TauPL.eq.0.0_dp).and.(sigma.ne.1.0_dp)) then 
			Bench_Folder = './NSU_LT_Results/Bench_Files/'
		else if ((TauPL.ne.0.0_dp).and.(sigma.ne.1.0_dp)) then 
			Bench_Folder = './NSU_PT_Results/Bench_Files/'
		else if ((TauPL.eq.0.0_dp).and.(sigma.eq.1.0_dp)) then 
			Bench_Folder = './SU_LT_Results/Bench_Files/'
		else if ((TauPL.ne.0.0_dp).and.(sigma.eq.1.0_dp)) then 
			Bench_Folder = './SU_PT_Results/Bench_Files/'
		end if 

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
		CALL Write_Benchmark_Results(1)

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
		psi_bench   = psi_bench
		DBN_bench   = DBN1
		tauw_at_bench = tauW_at
		Y_a_threshold_bench = Y_a_threshold

		Cons_bench   = Cons 
		Hours_bench  = Hours
		Aprime_bench = Aprime
		Value_bench  = ValueFunction

		CALL ComputeLaborUnits(EBAR, wage) 

!====================================================================================================
	PRINT*,''
	Print*,'Loading experiment'
	PRINT*,''

	tauK = 0.0_DP
	Y_a_threshold = Threshold_Factor*Ebar_bench 
	tauW_at = 0.017072675596579098_dp

	OPEN  (UNIT=1,  FILE=trim(Result_Folder)//'Exp_results_cons'  , STATUS='old', ACTION='read')
	OPEN  (UNIT=2,  FILE=trim(Result_Folder)//'Exp_results_aprime', STATUS='old', ACTION='read')
	OPEN  (UNIT=3,  FILE=trim(Result_Folder)//'Exp_results_hours' , STATUS='old', ACTION='read')
	OPEN  (UNIT=4,  FILE=trim(Result_Folder)//'Exp_results_value' , STATUS='old', ACTION='read')
	OPEN  (UNIT=5,  FILE=trim(Result_Folder)//'Exp_results_DBN'   , STATUS='old', ACTION='read')
	OPEN  (UNIT=60, FILE=trim(Result_Folder)//'Exp_results_GBAR'  , STATUS='old', ACTION='read')
	OPEN  (UNIT=7,  FILE=trim(Result_Folder)//'Exp_results_EBAR'  , STATUS='old', ACTION='read')
	OPEN  (UNIT=8,  FILE=trim(Result_Folder)//'Exp_results_NBAR'  , STATUS='old', ACTION='read')
	OPEN  (UNIT=9,  FILE=trim(Result_Folder)//'Exp_results_QBAR'  , STATUS='old', ACTION='read')
	OPEN  (UNIT=10, FILE=trim(Result_Folder)//'Exp_results_rr'    , STATUS='old', ACTION='read')
	OPEN  (UNIT=11, FILE=trim(Result_Folder)//'Exp_results_wage'  , STATUS='old', ACTION='read')
	OPEN  (UNIT=12, FILE=trim(Result_Folder)//'Exp_results_YBAR'  , STATUS='old', ACTION='read')

	READ (UNIT=1,  FMT=*), cons_exp
	READ (UNIT=2,  FMT=*), aprime_exp
	READ (UNIT=3,  FMT=*), hours_exp
	READ (UNIT=4,  FMT=*), Value_exp
	READ (UNIT=5,  FMT=*), DBN_exp
	READ (UNIT=60, FMT=*), GBAR_exp
	READ (UNIT=7,  FMT=*), EBAR_exp
	READ (UNIT=8,  FMT=*), NBAR_exp
	READ (UNIT=9,  FMT=*), QBAR_exp
	READ (UNIT=10, FMT=*), rr_exp
	READ (UNIT=11, FMT=*), wage_exp
	READ (UNIT=12, FMT=*), Y_exp

	CLOSE (unit=1)
	CLOSE (unit=2)
	CLOSE (unit=3)
	CLOSE (unit=4)
	CLOSE (unit=5)
	CLOSE (unit=60)
	CLOSE (unit=7)
	CLOSE (unit=8)
	CLOSE (unit=9)
	CLOSE (unit=10)
	CLOSE (unit=11)
	CLOSE (unit=12)

	tauK_exp  = tauK
	tauPL_exp = tauPL
	psi_exp   = psi
	tauw_at_exp = tauW_at
	Y_a_threshold_exp = Y_a_threshold

!====================================================================================================
	PRINT*,''
	Print*,'Total Consumption and Hours '
	PRINT*,''

	C_bench = sum(Cons_bench*DBN_bench)
	C_exp   = sum(Cons_exp*DBN_exp)

	H_bench = sum(Hours_bench*DBN_bench)
	H_exp   = sum(Hours_exp*DBN_exp)

!====================================================================================================
	PRINT*,''
	Print*,'Consumption Equivalent - Total'
	PRINT*,''

	CE_total = 100.0_dp*((Value_exp/Value_bench)**(1.0_dp/((1.0_dp-sigma)*gamma)) - 1.0_dp ) ;

!====================================================================================================
	PRINT*,''
	Print*,'Consumption Equivalent - Consumption'
	PRINT*,''

	Cons    = Cons_exp 
	Aprime  = Aprime_exp
	Hours   = Hours_bench 
	tauK    = tauK_exp
	tauW_at = tauW_at_exp 
	GBAR  	= GBAR_exp
	QBAR  	= QBAR_exp
	NBAR  	= NBAR_exp
	Ebar  	= EBAR_exp
	rr   	= rr_exp
	wage  	= wage_exp
	YBAR    = Y_exp
	Y_a_threshold = Y_a_threshold_exp

	CALL COMPUTE_VALUE_FUNCTION_SPLINE 
	Value_aux = ValueFunction

	CE_c  = 100.0_dp*((Value_aux/Value_bench)**(1.0_dp/((1.0_dp-sigma)*gamma)) - 1.0_dp ) ;

	CE_cl = 100_dp*(C_exp/C_bench - 1.0_dp)

	CE_cd = 100.0_dp*((Value_aux/Value_bench)**(1.0_dp/((1.0_dp-sigma)*gamma)) * C_bench/C_exp - 1.0_dp ) ;

!====================================================================================================
	PRINT*,''
	Print*,'Consumption Equivalent - Hours'
	PRINT*,''


	Cons    = Cons_exp
	Aprime  = Aprime_exp
	Hours   = Hours_bench 
	tauK    = tauK_bench
	tauW_at = tauW_at_bench 
	GBAR  	= GBAR_bench
	QBAR  	= QBAR_bench
	NBAR  	= NBAR_bench
	Ebar  	= EBAR_bench
	rr   	= rr_bench
	wage  	= wage_bench
	YBAR    = Y_bench
	Y_a_threshold = Y_a_threshold_bench

	CALL COMPUTE_VALUE_FUNCTION_SPLINE 
	Value_aux = ValueFunction

	Cons    = Cons_bench 
	Aprime  = Aprime_bench
	Hours   = Hours_exp
	tauK    = tauK_bench
	tauW_at = tauW_at_bench 
	GBAR  	= GBAR_bench
	QBAR  	= QBAR_bench
	NBAR  	= NBAR_bench
	Ebar  	= EBAR_bench
	rr   	= rr_bench
	wage  	= wage_bench
	YBAR    = Y_bench
	Y_a_threshold = Y_a_threshold_bench

	CALL COMPUTE_VALUE_FUNCTION_SPLINE

	CE_h  = 100.0_dp*((Value_exp/Value_aux)**(1.0_dp/((1.0_dp-sigma)*gamma)) - 1.0_dp ) ;

	CE_hl = 100_dp*(((1.0_dp-H_exp)/(1.0_dp-H_bench))**((1.0_dp-gamma)/gamma) - 1.0_dp)

	CE_hd = 100.0_dp*((ValueFunction/Value_bench)**(1.0_dp/((1.0_dp-sigma)*gamma)) &
	        &  * ((1.0_dp-H_bench)/(1.0_dp-H_exp))**((1.0_dp-gamma)/gamma) - 1.0_dp ) ;

!====================================================================================================
	PRINT*,''
	Print*,'Report Output'
	PRINT*,''

	CE_total_bench 	= sum(CE_total*DBN_bench )
	CE_c_bench     	= sum(CE_c*DBN_bench )
	CE_cl_bench    	= sum(CE_cl*DBN_bench )
	CE_cd_bench 	= sum(CE_cd*DBN_bench )
	CE_h_bench 		= sum(CE_h*DBN_bench )
	CE_hl_bench 	= sum(CE_hl*DBN_bench )
	CE_hd_bench 	= sum(CE_hd*DBN_bench )

	CE_total_NB_bench 	= sum(CE_total(1,:,:,:,:)*DBN_bench(1,:,:,:,:))/sum(DBN_bench(1,:,:,:,:))
	CE_c_NB_bench 		= sum(CE_c(1,:,:,:,:)*DBN_bench(1,:,:,:,:))/sum(DBN_bench(1,:,:,:,:))
	CE_cl_NB_bench 		= sum(CE_cl(1,:,:,:,:)*DBN_bench(1,:,:,:,:))/sum(DBN_bench(1,:,:,:,:))
	CE_cd_NB_bench 		= sum(CE_cd(1,:,:,:,:)*DBN_bench(1,:,:,:,:))/sum(DBN_bench(1,:,:,:,:))
	CE_h_NB_bench 		= sum(CE_h(1,:,:,:,:)*DBN_bench(1,:,:,:,:))/sum(DBN_bench(1,:,:,:,:))
	CE_hl_NB_bench 		= sum(CE_hl(1,:,:,:,:)*DBN_bench(1,:,:,:,:))/sum(DBN_bench(1,:,:,:,:))
	CE_hd_NB_bench 		= sum(CE_hd(1,:,:,:,:)*DBN_bench(1,:,:,:,:))/sum(DBN_bench(1,:,:,:,:))

	CE_total_exp 	= sum(CE_total*DBN_exp )
	CE_c_exp 		= sum(CE_c*DBN_exp )
	CE_cl_exp 		= sum(CE_cl*DBN_exp )
	CE_cd_exp 		= sum(CE_cd*DBN_exp )
	CE_h_exp 		= sum(CE_h*DBN_exp )
	CE_hl_exp 		= sum(CE_hl*DBN_exp )
	CE_hd_exp 		= sum(CE_hd*DBN_exp )

	CE_total_NB_exp 	= sum(CE_total(1,:,:,:,:)*DBN_exp(1,:,:,:,:))/sum(DBN_exp(1,:,:,:,:))
	CE_c_NB_exp 		= sum(CE_c(1,:,:,:,:)*DBN_exp(1,:,:,:,:))/sum(DBN_exp(1,:,:,:,:))
	CE_cl_NB_exp 		= sum(CE_cl(1,:,:,:,:)*DBN_exp(1,:,:,:,:))/sum(DBN_exp(1,:,:,:,:))
	CE_cd_NB_exp 		= sum(CE_cd(1,:,:,:,:)*DBN_exp(1,:,:,:,:))/sum(DBN_exp(1,:,:,:,:))
	CE_h_NB_exp 		= sum(CE_h(1,:,:,:,:)*DBN_exp(1,:,:,:,:))/sum(DBN_exp(1,:,:,:,:))
	CE_hl_NB_exp 		= sum(CE_hl(1,:,:,:,:)*DBN_exp(1,:,:,:,:))/sum(DBN_exp(1,:,:,:,:))
	CE_hd_NB_exp 		= sum(CE_hd(1,:,:,:,:)*DBN_exp(1,:,:,:,:))/sum(DBN_exp(1,:,:,:,:))

	OPEN  (UNIT=1,  FILE=trim(Result_Folder)//'CE_output.txt'  , STATUS='replace')
	WRITE (UNIT=1,  FMT=*) 'Benchmark - Aggregate'
	WRITE (UNIT=1,  FMT=*) 'CE',CE_total_bench, 'CE_c',CE_c_bench,'CE_h',CE_h_bench
	WRITE (UNIT=1,  FMT=*) 'CE_c',CE_c_bench,'CE_cl',CE_cl_bench,'CE_cd',CE_cd_bench
	WRITE (UNIT=1,  FMT=*) 'CE_h',CE_h_bench,'CE_hl',CE_hl_bench,'CE_hd',CE_hd_bench
	WRITE (UNIT=1,  FMT=*)' '

	WRITE (UNIT=1,  FMT=*)'Benchmark - NewBorn'
	WRITE (UNIT=1,  FMT=*) 'CE',CE_total_NB_bench, 'CE_c',CE_c_NB_bench,'CE_h',CE_h_NB_bench
	WRITE (UNIT=1,  FMT=*) 'CE_c',CE_c_NB_bench,'CE_cl',CE_cl_NB_bench,'CE_cd',CE_cd_NB_bench
	WRITE (UNIT=1,  FMT=*) 'CE_h',CE_h_NB_bench,'CE_hl',CE_hl_NB_bench,'CE_hd',CE_hd_NB_bench
	WRITE (UNIT=1,  FMT=*)' '

	WRITE (UNIT=1,  FMT=*)'Experimental - Aggregate'
	WRITE (UNIT=1,  FMT=*) 'CE',CE_total_exp, 'CE_c',CE_c_exp,'CE_h',CE_h_exp
	WRITE (UNIT=1,  FMT=*) 'CE_c',CE_c_exp,'CE_cl',CE_cl_exp,'CE_cd',CE_cd_exp
	WRITE (UNIT=1,  FMT=*) 'CE_h',CE_h_exp,'CE_hl',CE_hl_exp,'CE_hd',CE_hd_exp
	WRITE (UNIT=1,  FMT=*)' '

	WRITE (UNIT=1,  FMT=*)'Experimental - NewBorn'
	WRITE (UNIT=1,  FMT=*) 'CE',CE_total_NB_exp, 'CE_c',CE_c_NB_exp,'CE_h',CE_h_NB_exp
	WRITE (UNIT=1,  FMT=*) 'CE_c',CE_c_NB_exp,'CE_cl',CE_cl_NB_exp,'CE_cd',CE_cd_NB_exp
	WRITE (UNIT=1,  FMT=*) 'CE_h',CE_h_NB_exp,'CE_hl',CE_hl_NB_exp,'CE_hd',CE_hd_NB_exp
	WRITE (UNIT=1,  FMT=*)' '

	CLOSE (unit=1)

	OPEN  (UNIT=1,  FILE=trim(Result_Folder)//'CE_Files/CE_total'  , STATUS='replace')
	OPEN  (UNIT=2,  FILE=trim(Result_Folder)//'CE_Files/CE_c'  , STATUS='replace')
	OPEN  (UNIT=3,  FILE=trim(Result_Folder)//'CE_Files/CE_cl'  , STATUS='replace')
	OPEN  (UNIT=4,  FILE=trim(Result_Folder)//'CE_Files/CE_cd'  , STATUS='replace')
	OPEN  (UNIT=5,  FILE=trim(Result_Folder)//'CE_Files/CE_h'  , STATUS='replace')
	OPEN  (UNIT=6,  FILE=trim(Result_Folder)//'CE_Files/CE_hl'  , STATUS='replace')
	OPEN  (UNIT=7,  FILE=trim(Result_Folder)//'CE_Files/CE_hd'  , STATUS='replace')
	
	WRITE (UNIT=1,  FMT=*) CE_total
	WRITE (UNIT=2,  FMT=*) CE_c
	WRITE (UNIT=3,  FMT=*) CE_cl
	WRITE (UNIT=4,  FMT=*) CE_cd
	WRITE (UNIT=5,  FMT=*) CE_h
	WRITE (UNIT=6,  FMT=*) CE_hl
	WRITE (UNIT=7,  FMT=*) CE_hd
	
	CLOSE (unit=1)
	CLOSE (unit=2)
	CLOSE (unit=3)
	CLOSE (unit=4)
	CLOSE (unit=5)
	CLOSE (unit=6)
	CLOSE (unit=7)
	



end Program Consumption_Equivalent