MODULE global
    USE parameters

    ! Switch for non-separable vs separable utility
    logical :: NSU_Switch, Log_Switch

    ! Switch for progressive vs linear labor income taxes
    logical :: Progressive_Tax_Switch

    ! Foder to store results
    character(100) :: Result_Folder, Bench_Simul_Folder

    ! "params" determines the economy to be simulated. 
    	! It contains the values of: beta, rho_z, sigma_z_eps, sigma_lambda_eps and phi. In order.
    real(DP) , dimension(6)     :: params

    ! Population size (pop) and survival probability by age (survP)
    REAL(DP), DIMENSION(MaxAge) :: pop, survP          
	
    ! Labor efficiency shocks: 
    	! Lyfe cycle component (by age)
		REAL(DP), DIMENSION(RetAge) :: kappagrid

		! Transitory component: 
			! grid (egrid), invariant distribution (Ge), CDF of invariant distribution (cdf_Ge)
			REAL(DP), DIMENSION(ne)        :: egrid, Ge, cdf_Ge
			! transition matrix (pr_e), CDF of transition matrix (by row) (cdf_pr_e)
			REAL(DP), DIMENSION(ne,ne)     :: pr_e, cdf_pr_e
			! distribution or "e" by age. This is constructed with pr_e and the assumption that in the first period, everyone starts at median "e"
   			REAL(DP), DIMENSION(MaxAge,ne) :: Ge_byage, cdf_Ge_byage

		! Permanent component:
			! grid (lambdagrid), invariant distribution (Glambda), CDF of invariant distribution (cdf_Glambda)
            REAL(DP), DIMENSION(nlambda)    :: lambdagrid,  Glambda, cdf_Glambda
            ! transition matrix (pr_lambda), CDF of transition matrix (by row) (cdf_pr_lambda)
			REAL(DP), DIMENSION(nlambda,nlambda) :: pr_lambda, cdf_pr_lambda

	! Entrepreneurial ability
		! grid (zgrid), invariant distribution (Glz), CDF of invariant distribution (cdf_Gz)
		REAL(DP), DIMENSION(nz)    :: zgrid , Gz, cdf_Gz
		! transition matrix (pr_z), CDF of transition matrix (by row) (cdf_pr_z)
		REAL(DP), DIMENSION(nz,nz) :: pr_z, cdf_pr_z
		! grid (xgrid), invariant distribution (Gx), CDF of invariant distribution (cdf_Gx)
		REAL(DP), DIMENSION(nx)    :: xgrid 
		REAL(DP), DIMENSION(nx,nz,MaxAge) :: Gx, cdf_Gx
		REAl(DP), DIMENSION(nx,nz) :: xz_grid
		! transition matrix (pr_x), CDF of transition matrix (by row) (cdf_pr_x)
		REAL(DP), DIMENSION(nx,nx,nz,MaxAge) :: pr_x, cdf_pr_x

	! Retirement income 
	REAL(DP), DIMENSION(nlambda,ne) :: phi_lambda_e   ! phi_lambda_e is the income replacement ratio used to compute SS payments
	REAL(DP), DIMENSION(nlambda,ne) :: RetY_lambda_e, RetY_lambda_e_aux  ! Retirement income = phi_lambda_e*Ebar, the first 45 years are assigned zero. 

	! Labor efficiency units
		! eff_un(age,lambda,e) = kappa(age)*lambda*e
		! yh(age,lambda,e)     = Wage * eff_un
		REAL(DP), DIMENSION(MaxAge  , nlambda, ne) :: eff_un,  yh

	! Policy function and value function (defined on the exogenous grid)
    REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: Cons, Hours, Aprime, Income_AT
    REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: Cons_bench, Hours_bench, Aprime_bench, Cons_exp, Hours_exp, Aprime_exp 
    REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: ValueFunction, ValueFunction_bench, ValueFunction_exp
    REAL(DP), dimension(:,:,:,:,:,:), allocatable :: Cons_Eq_Welfare
	! Policy function and value function (defined on the adjusted grid for breakpoints)
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: Cons_t, Hours_t, Aprime_t
	!REAL(DP), DIMENSION(MaxAge,na+nz,nz,nlambda,ne) :: Cons_t, Hours_t, Aprime_t
 
 	! Aggregate variables
	 	! Benchmark values of Q, N, E, Wage, R, G, Y
	    REAL(DP) :: QBAR_bench, NBAR_bench, Ebar_bench, wage_bench, P_bench, GBAR_bench, Y_bench, W_bench, R_bench, SSC_Payments_bench
	    ! Experiment values of Q, N, E, Wage, R, G, Y
	    REAL(DP) :: QBAR_exp,   NBAR_exp,   Ebar_exp,   wage_exp,   P_exp,   GBAR_exp, GBAR_exp_old, Y_exp, R_exp
	    ! Values for aggregate variables (used when solving a given economy)
	    REAL(DP) :: P, Ebar , wage, NBAR, QBAR, YBAR, GBAR, R, SSC_Payments
	    ! Wealth tax threshold as proportion of mean benchmark wealth
	    REAL(DP) :: Wealth_factor, Threshold_Share

	! Asset, resources and marginal benefit of wealth grids
		! Note that these grids will change as we change tax rates and for each intermediate good price
		! For each tax rate we weill have different Y grids
		REAL(DP), DIMENSION(na)      :: agrid
	    REAL(DP), DIMENSION(fine_na) :: fine_agrid
	    REAL(DP), DIMENSION(na,nz,nx):: YGRID, MBGRID, YGRID_aux, MBGRID_aux
	    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: YGRID_t, MBGRID_t
	    REAL(DP), DIMENSION(:)    , ALLOCATABLE :: agrid_t
	    INTEGER                      :: na_t

    ! Capital markets
    	REAL(DP), DIMENSION(na,nz,nx) :: K_mat, Pr_mat, Wealth_mat
    	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: V_Pr, V_Pr_bench, V_Pr_exp, Firm_Wealth
    	REAL(DP), DIMENSION(na,nz,nlambda) :: V_Pr_nb, V_Pr_nb_bench, V_Pr_nb_exp
	
	! Values for taxes in benchmark and experiment
    REAL(DP) :: tauk_bench, tauPL_bench, psi_bench, tauw_bt_bench, tauw_at_bench, Y_a_threshold_bench 
    REAL(DP) :: tauk_exp,   tauPL_exp,   psi_exp,   tauw_bt_exp,   tauw_at_exp,   Y_a_threshold_exp

    ! Values for taxes (when solving an economy)
    REAL(DP) :: tauK
    REAL(DP) :: tauW_bt, tauW_at ! Wealth taxes below threshold and above threshold
    REAL(DP) :: Y_a_threshold = 0.0_dp ! Value of the threshold for change in tauW

	! Taxes
	! Wealth tax: minimum wealth tax to consider and increments for balancing budget
		REAL(DP) :: tauWmin_bt=0.00_DP, tauWinc_bt=0.000_DP ! Minimum tax below threshold and increments
		REAL(DP) :: tauWmin_at=0.012_DP, tauWinc_at=0.002_DP ! Minimum tax above threshold and increments
		REAL(DP) :: Threshold_Factor = 0.00_dp 
	! Consumption tax
		REAL(DP) :: tauC=0.075_DP
    ! Labor income tax: This is a progresive tax.
	! 1-psi controls the level of tax, and tauPL controls progressivity
		REAL(DP) :: tauPL, psi


    ! Auxiliary variables to find wealth tax that balances the budget in experiment economy
    REAL(DP) :: tauWindx, tauW_low_bt, tauW_up_bt, tauW_low_at, tauW_up_at

	! Counters for the age, and index of lamnbda, z, a and e
    INTEGER :: age, lambdai, zi, ai, ei, xi    

    ! Distribution of population by age, a, z, lambda, e
    REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne, nx) ::DBN1, DBN_bench, DBN_exp

    ! Stats and distribution in equilibrium
	    ! Distribution of assets
	    REAL(DP), DIMENSION(na)  :: pr_a_dbn, cdf_a_dbn
	    ! Mass of assets by a_grid (defined as a_grid(i)sum(DBN(:,i,:,:,:)))
	    REAL(DP), DIMENSION(na)  :: tot_a_by_grid, cdf_tot_a_by_grid
	    ! Percentiles of the asset distribution and mass of total assets by them
	    INTEGER , DIMENSION(100) ::  prctile_ai_ind
	    REAL(DP), DIMENSION(100) ::  cdf_tot_a_by_prctile, prctile_ai
    	! Other stats
	    REAL(DP) :: pop_25_60 , tothours_25_60, pop_pos_earn_25_60, tot_log_earnings_25_60, mean_log_earnings_25_60 
	    REAL(DP) :: meanhours_25_60, Var_Log_Earnings_25_60, Std_Log_Earnings_25_60, MeanWealth, Wealth_Output
	    REAL(DP) :: MeanReturn, StdReturn, VarReturn
	    REAL(DP), DIMENSION(nz) :: MeanReturn_by_z(nz), size_by_z(nz), Wealth_by_z(nz)
	    REAL(DP) :: prct1_wealth, prct10_wealth, prct20_wealth, prct40_wealth
	    REAL(DP) :: SSE_Moments, Min_SSE_Moments
	    REAL(DP) :: MeanCons, External_Debt_GDP
	    REAL(DP) :: Mean_Firm_Wealth, Firm_Wealth_by_z
	    ! Welfare measures
	    REAL(DP) :: Welfare_Gain_Pop_bench, Welfare_Gain_Pop_exp, Welfare_Gain_NB_bench, Welfare_Gain_NB_exp
	    REAL(DP) :: CE_NEWBORN, Av_Util_Pop, Av_Util_NB

    ! Objective moments (not currently in use)	
    	REAL(DP), DIMENSION(6)  ::  Min_Moments
    
    ! Switch to solve benchmark or experiment
    INTEGER :: solving_bench
    
    ! Seed for random number generator
    INTEGER :: newiseed

    ! Bequest Stats
	REAL(DP) :: Threshold_Share_bq, Bequest_Wealth, BQ_top_x(7), Mean_Bequest

	! top certain percentiles of the PV_a distribution
    REAL(DP) :: FW_top_x_share(6) 

    ! Frisch Elasticity
    REAL(DP) :: Frisch_Elasticity, Size_Frisch, Hours_Frisch, Frisch_Elasticity_2

    ! Extra tax information
    REAL(DP) ::  GBAR_K,  GBAR_W, GBAR_L, GBAR_C, Tot_Lab_Inc

    ! Transition Yo
    INTEGER :: ti 
    ! Prices and Quantities
    REAL(DP), DIMENSION(T+1) :: R_tr, P_tr, QBAR_tr, NBAR_tr, YBAR_tr, Wage_tr, EBAR_tr, K_tr, C_tr
    ! Government Budget
    REAL(DP), DIMENSION(T+1) :: GBAR_tr, GBAR_K_tr, GBAR_W_tr, GBAR_L_tr, GBAR_C_tr, SSC_Payments_tr, Tot_Lab_Inc_tr
    ! Policy function and value function (defined on the exogenous grid)
    REAL(DP), DIMENSION(:,:,:,:,:,:,:), allocatable :: Cons_tr, Hours_tr, Aprime_tr, DBN_tr, ValueFunction_tr
    ! Welfare Gain
    REAL(DP), DIMENSION(:,:,:,:,:,:)  , allocatable :: CE1_tr
    REAL(DP) :: CE1_nb_tr,  CE1_pop_tr, CE2_nb_tr,  CE2_pop_tr
	! Policy function and value function (defined on the adjusted grid for breakpoints)
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: Cons_t_pr, Hours_t_pr


Contains 

Subroutine Allocate_Variables
    !! Allocate variables
    allocate( Cons(               MaxAge,na,nz,nlambda,ne,nx) )
    allocate( Hours(              MaxAge,na,nz,nlambda,ne,nx) )
    allocate( Aprime(             MaxAge,na,nz,nlambda,ne,nx) )
    allocate( Income_AT(          MaxAge,na,nz,nlambda,ne,nx) )
    allocate( Cons_bench(         MaxAge,na,nz,nlambda,ne,nx) )
    allocate( Hours_bench(        MaxAge,na,nz,nlambda,ne,nx) )
    allocate( Aprime_bench(       MaxAge,na,nz,nlambda,ne,nx) )
    allocate( Cons_exp(           MaxAge,na,nz,nlambda,ne,nx) )
    allocate( Hours_exp(          MaxAge,na,nz,nlambda,ne,nx) )
    allocate( Aprime_exp(         MaxAge,na,nz,nlambda,ne,nx) )
    allocate( ValueFunction(      MaxAge,na,nz,nlambda,ne,nx) )
    allocate( ValueFunction_bench(MaxAge,na,nz,nlambda,ne,nx) )
    allocate( ValueFunction_exp(  MaxAge,na,nz,nlambda,ne,nx) )
    allocate( Cons_Eq_Welfare(    MaxAge,na,nz,nlambda,ne,nx) )
    allocate( V_Pr(               MaxAge,na,nz,nlambda,ne,nx) )
    allocate( V_Pr_bench(         MaxAge,na,nz,nlambda,ne,nx) )
    allocate( V_Pr_exp(           MaxAge,na,nz,nlambda,ne,nx) )
    allocate( Firm_Wealth(        MaxAge,na,nz,nlambda,ne,nx) )


    allocate( Cons_tr(            MaxAge,na,nz,nlambda,ne,nx,T+1) )
    allocate( Hours_tr(           MaxAge,na,nz,nlambda,ne,nx,T+1) )
    allocate( Aprime_tr(          MaxAge,na,nz,nlambda,ne,nx,T+1) )
    allocate( DBN_tr( 			  MaxAge,na,nz,nlambda,ne,nx,T+1) )
    allocate( ValueFunction_tr(   MaxAge,na,nz,nlambda,ne,nx,MaxAge) )
	allocate( Cons_t_pr(          MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Hours_t_pr(         MaxAge,na,nz,nlambda,ne,nx) )
    allocate( CE1_tr( 			  MaxAge,na,nz,nlambda,ne,nx) )
end Subroutine Allocate_Variables
   
END MODULE global