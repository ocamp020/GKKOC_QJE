
Module GKK_Stats
	use parameters
	use global
	use Toolbox
	use programfunctions
	    
	Contains




!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE COMPUTE_STATS()
	use omp_lib
	IMPLICIT NONE
	integer  :: prctile, group, i, j, age_group_counter, age2
	real(DP), dimension(nz)    :: Capital_by_z, DBN_Z, CDF_Z
	real(DP) :: MeanATReturn, StdATReturn, VarATReturn, MeanATReturn_by_z(nz), Mean_Capital
	real(DP) :: Std_k_Return,    Var_K_Return,    Mean_K_Return_by_z(nz)
	real(DP) :: Std_AT_K_Return, Var_AT_K_Return, Mean_AT_K_Return_by_z(nz)
	real(DP), dimension(max_age_category)    :: A_Age, Ap_Age, Y_Age, S_Age, S_Rate_A_Age, S_Rate_Y_Age 
	real(DP), dimension(max_age_category,nz) :: A_AZ , Ap_AZ , Y_AZ , S_AZ , S_Rate_A_AZ , S_Rate_Y_AZ  
	real(DP), dimension(3) 					 :: A_W  , Ap_W  , Y_W  , S_W  , S_Rate_A_W  , S_Rate_Y_W  
	real(DP), dimension(MaxAge) 			 :: constrained_firms_age, size_by_age
	real(DP) 	   :: Pr_mat_bench(na,nz,nx), FW_top_x(6), prctile_FW(6), prctile_bq(5)
	real(DP)  	   :: low_pct, high_pct, a, b, c, CCDF_c, c_low, c_high
	real(DP)       :: Frisch_Aux, Frisch_Aux_2
	real(DP)       :: K_Inc_aux, L_Inc_aux, K_Tax_aux, L_Tax_aux, K_Inc_bench, L_Inc_bench, K_Tax_bench, L_Tax_bench
	character(100) :: rowname
	integer        :: age_limit(max_age_category+1), draft_age_limit(draft_age_category+1)
	real(DP), dimension(MaxAge, nz) 		   :: size_by_age_z, leverage_age_z, constrained_firms_age_z 
	real(DP), dimension(draft_age_category,nz) :: size_draft_group_z, wealth_draft_group_z, capital_draft_group_z, &
		& Cons_draft_group_z, Hours_draft_group_z, Ap_draft_group_z, & 
		& K_Tax_draft_group_z, L_Tax_draft_group_z, K_Tax_Inc_draft_group_z, L_Tax_Inc_draft_group_z, &
		& T_Inc_draft_group_z, K_Inc_draft_group_z, L_Inc_draft_group_z, av_K_Inc_draft_group_z, av_L_Inc_draft_group_z, &
		&      Tax_Increase_tk_draft_group_z,      Tax_Increase_tl_draft_group_z,      Tax_Increase_draft_group_z, &
		& Tax_Rate_Increase_tk_draft_group_z, Tax_Rate_Increase_tl_draft_group_z, Tax_Rate_Increase_draft_group_z, & 
		& Inc_Increase_draft_group_z, K_Inc_Increase_draft_group_z, L_Inc_Increase_draft_group_z, &
		& Return_draft_group_z, Return_AT_draft_group_z, &
		& Entrepreneur_10_draft_group_z, Entrepreneur_50_draft_group_z
	real(DP), dimension(draft_age_category,draft_z_category) :: size_draft_group, &
		& wealth_draft_group,  av_wealth_draft_group, frac_wealth_draft_group, & 
		& capital_draft_group,  av_capital_draft_group, frac_capital_draft_group, &
		& cons_draft_group, hours_draft_group, Ap_draft_group,  av_Ap_draft_group, frac_Ap_draft_group, &
		& K_Tax_draft_group, L_Tax_draft_group, K_Tax_Inc_draft_group, L_Tax_Inc_draft_group, &
		& T_Inc_draft_group, K_Inc_draft_group, L_Inc_draft_group, av_K_Inc_draft_group, av_L_Inc_draft_group, &
		& frac_K_Tax_draft_group, frac_L_Tax_draft_group, frac_K_Inc_draft_group, frac_L_Inc_draft_group, &
		& Tax_Increase_tk_draft_group, Tax_Increase_tl_draft_group, Tax_Increase_draft_group, &
		& Tax_Rate_Increase_tk_draft_group, Tax_Rate_Increase_tl_draft_group, Tax_Rate_Increase_draft_group, &
		& Inc_Increase_draft_group, K_Inc_Increase_draft_group, L_Inc_Increase_draft_group, &
		& Return_draft_group, Return_AT_draft_group, &
		& Entrepreneur_10_draft_group, Entrepreneur_50_draft_group
	real(DP) :: DBN_az(na,nz)
	real(DP) :: Z_share_top_wealth(draft_age_category,nz), draft_group_share_top_wealth(draft_age_category,draft_z_category), &
			&	A_share_top_wealth(draft_age_category,nz), draft_group_wealth_share_top_wealth(draft_age_category,draft_z_category) 
	real(DP) :: DBN_azx(na,nz,nx), BT_Return(na,nz,nx), DBN_azx_vec(na*nz*nx), Return_vec(na*nz*nx)
	integer  :: ind_lo, ind_hi, prctile_ai_ind_age(14)
	real(DP) :: pct_graph_lim(14), ret_by_wealth(draft_age_category+1,13), pct_graph_wealth(draft_age_category+1,13)
	real(DP), dimension(:,:,:,:,:,:), allocatable :: DBN_bq, Total_Income ! , Firm_Output, Firm_Profit
	integer , dimension(:,:,:,:,:,:), allocatable :: constrained_firm_ind
	real(DP), dimension(:), allocatable :: DBN_vec, Firm_Wealth_vec, CDF_Firm_Wealth, BQ_vec, DBN_bq_vec, CDF_bq, Inc_vec


	allocate(DBN_vec(			size(DBN1)))
	allocate(Firm_Wealth_vec(	size(DBN1)))
	allocate(CDF_Firm_Wealth(	size(DBN1)))
	allocate(BQ_vec(			size(DBN1)))
	allocate(DBN_bq_vec(		size(DBN1)))
	allocate(CDF_bq(			size(DBN1)))
	allocate(Inc_vec(			size(DBN1)))
	! allocate(Firm_Output( MaxAge,na,nz,nlambda,ne,nx))
	! allocate(Firm_Profit( MaxAge,na,nz,nlambda,ne,nx))
	allocate(DBN_bq(       MaxAge,na,nz,nlambda,ne,nx))
	allocate(Total_Income( MaxAge,na,nz,nlambda,ne,nx))
	allocate(constrained_firm_ind(MaxAge,na,nz,nlambda,ne,nx))


	!$ call omp_set_num_threads(20)
	! $ print *, "OMP Test Message"

	! print*, ' '; print*,' Entering Compute_Stats'; print*, ' '; 
	
	! Age Brackets
		draft_age_limit = [0, 1, 15, 30, 45, MaxAge ] 

	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	! Size by age-z and age_group_z
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
		! print*, 'Test DBN_az'
		OPEN (UNIT=90, FILE=trim(Result_Folder)//'size_by_age_z.txt', STATUS='replace')   
		size_by_age_z =0.0_DP
		DO age=1,MaxAge 
		    DO zi=1,nz
		        size_by_age_z(age, zi) = sum(DBN1(age,:,zi,:,:,:))
		    ENDDO ! zi

		    WRITE  (UNIT=90, FMT=*)  size_by_age_z(age, :) 
		ENDDO
		CLOSE(unit=90)

		DBN_Z = sum(sum(sum(sum(sum(DBN1,6),5),4),2),1) 
		do zi=1,nz 
			CDF_Z(zi) = sum(DBN_Z(1:zi))
		enddo 
		
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	! Distribution of Assets
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
		! print*, 'Test assets'
		DO ai=1,na
		     pr_a_dbn(ai)          = sum(DBN1(:,ai,:,:,:,:)) 
		     cdf_a_dbn(ai)         = sum( pr_a_dbn(1:ai) )      
		     tot_a_by_grid(ai)     = sum(DBN1(:,ai,:,:,:,:) * agrid(ai) )
		     cdf_tot_a_by_grid(ai) = sum(tot_a_by_grid(1:ai))   
		!     print*, pr_a_dbn(ai), cdf_a_dbn(ai)
		ENDDO
		cdf_a_dbn = cdf_a_dbn + 1.0_DP - cdf_a_dbn(na)

		! Percentage of the population above wealth tax threshold
		Threshold_Share = 0.0_dp
		do ai=1,na 
			if (agrid(ai).gt.Y_a_threshold) then 
				Threshold_Share = Threshold_Share + pr_a_dbn(ai)
			end if 
		end do 

		
		! FIND THE ai THAT CORRESPONDS TO EACH PRCTILE OF WEALTH DBN & WEALTH HELD BY PEOPLE LOWER THAN THAT PRCTILE
		DO prctile=1,100
		    ai=1
		    DO while (cdf_a_dbn(ai) .lt. (REAL(prctile,8)/100.0_DP-0.000000000000001))
		        ai=ai+1
		    ENDDO
		    prctile_ai_ind(prctile) = ai
		    prctile_ai(prctile)     = agrid(ai)
		    ! print*,prctile, REAL(prctile,8)/100.0_DP,  ai
		    IF (ai .gt. 1) THEN
		        cdf_tot_a_by_prctile(prctile)  = cdf_tot_a_by_grid(ai-1) + (REAL(prctile,8)/100.0_DP - cdf_a_dbn(ai-1))*agrid(ai)
		    else
		        cdf_tot_a_by_prctile(prctile)  = (REAL(prctile,8)/100.0_DP )*agrid(ai)     
		    ENDIF
		ENDDO
		
		prct1_wealth  = 1.0_DP-cdf_tot_a_by_prctile(99)/cdf_tot_a_by_prctile(100)
		prct10_wealth = 1.0_DP-cdf_tot_a_by_prctile(90)/cdf_tot_a_by_prctile(100)
		prct20_wealth = 1.0_DP-cdf_tot_a_by_prctile(80)/cdf_tot_a_by_prctile(100)
		prct40_wealth = 1.0_DP-cdf_tot_a_by_prctile(60)/cdf_tot_a_by_prctile(100)

		! Z Composition of top wealth groups
		DBN_az = sum(sum(sum(sum(DBN1,6),5),4),1) 
		do zi=1,nz 
		Z_share_top_wealth(1,zi) = sum( DBN_az(:,zi) , (agrid.ge.prctile_ai(99)) )/sum( DBN_az , spread((agrid.ge.prctile_ai(99)),2,nz) )  
		Z_share_top_wealth(2,zi) = sum( DBN_az(:,zi) , (agrid.ge.prctile_ai(95)) )/sum( DBN_az , spread((agrid.ge.prctile_ai(95)),2,nz) )  
		Z_share_top_wealth(3,zi) = sum( DBN_az(:,zi) , (agrid.ge.prctile_ai(90)) )/sum( DBN_az , spread((agrid.ge.prctile_ai(90)),2,nz) )  
		Z_share_top_wealth(4,zi) = sum( DBN_az(:,zi) , (agrid.ge.prctile_ai(50)) )/sum( DBN_az , spread((agrid.ge.prctile_ai(50)),2,nz) )  
		Z_share_top_wealth(5,zi) = sum( DBN_az(:,zi) , (agrid.ge.prctile_ai(25)) )/sum( DBN_az , spread((agrid.ge.prctile_ai(25)),2,nz) )  

		A_share_top_wealth(1,zi) = sum( agrid*DBN_az(:,zi) , (agrid.ge.prctile_ai(99)) )&
										& /sum( agrid*sum(DBN_az,2) , agrid.ge.prctile_ai(99) )  
		A_share_top_wealth(2,zi) = sum( agrid*DBN_az(:,zi) , (agrid.ge.prctile_ai(95)) )&
										& /sum( agrid*sum(DBN_az,2) , agrid.ge.prctile_ai(95) )  
		A_share_top_wealth(3,zi) = sum( agrid*DBN_az(:,zi) , (agrid.ge.prctile_ai(90)) )&
										& /sum( agrid*sum(DBN_az,2) , agrid.ge.prctile_ai(90) )  
		A_share_top_wealth(4,zi) = sum( agrid*DBN_az(:,zi) , (agrid.ge.prctile_ai(50)) )&
										& /sum( agrid*sum(DBN_az,2) , agrid.ge.prctile_ai(50) )  
		A_share_top_wealth(5,zi) = sum( agrid*DBN_az(:,zi) , (agrid.ge.prctile_ai(25)) )&
										& /sum( agrid*sum(DBN_az,2) , agrid.ge.prctile_ai(25) )  
		enddo 

		! Composition by draft groups of top wealth groups
		draft_group_share_top_wealth        = Draft_Table(Z_share_top_wealth,DBN_z,.true.)
		draft_group_wealth_share_top_wealth = Draft_Table(A_share_top_wealth,DBN_z,.true.)

	    ! Write results in file
	    OPEN (UNIT=81, FILE=trim(Result_Folder)//'draft_group_share_top_wealth.txt', STATUS='replace') 
	    do age = 1,draft_age_category
		    WRITE  (UNIT=81, FMT=*)  draft_group_share_top_wealth(age,:)
		ENDDO
		close(unit=81)

		OPEN (UNIT=81, FILE=trim(Result_Folder)//'draft_group_wealth_share_top_wealth.txt', STATUS='replace') 
	    do age = 1,draft_age_category
		    WRITE  (UNIT=81, FMT=*)  draft_group_wealth_share_top_wealth(age,:)
		ENDDO
		close(unit=81)

		! Find top 99.9 and top 99.99
		! Top 99.9
		 	ai=1
		    DO while (cdf_a_dbn(ai) .lt. (99.9_dp/100.0_DP-0.000000000000001))
		        ai=ai+1
		    ENDDO
			prct999_wealth  = sum(tot_a_by_grid(ai:))/sum(tot_a_by_grid)
			
		! Top 99.99
		 	ai=1
		    DO while (cdf_a_dbn(ai) .lt. (99.99_dp/100.0_DP-0.000000000000001))
		        ai=ai+1
		    ENDDO
			prct9999_wealth = sum(tot_a_by_grid(ai:))/sum(tot_a_by_grid)

		! Print
		print*,' '
		print*,'Top Wealth Shares'
			print*,'	Top  0.01: ai',ai,'Wealth_Share',prct9999_wealth!,'Total_Wealth',sum(tot_a_by_grid)
			print*,'	Top  0.10: ai',ai,'Wealth_Share',prct999_wealth!,'Total_Wealth',sum(tot_a_by_grid)
			print*,'	Top  1.00: ai',prctile_ai_ind(99),'Wealth_Share',prct1_wealth
			print*,'	Top 10.00: ai',prctile_ai_ind(90),'Wealth_Share',prct10_wealth
		print*,' '

	
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	! Labor Earnings and Hours of working age population
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
		! print*, 'Test labor earnings'
		! COMPUTE AVERAGE HOURS FOR AGES 25-60 (5-40 IN THE MODEL) INCLUDING NON-WORKERS
		! COMPUTE VARIANCE OF LOG EARNINGS FOR 25-60 FOR THOSE WHO WORK MORE THAN 260 HOURS
		! WHICH CORRESPOND TO 0.055 IN THE MODEL
		pop_25_60        	   = 0.0_DP
		tothours_25_60         = 0.0_DP
		pop_pos_earn_25_60     = 0.0_DP
		tot_log_earnings_25_60 = 0.0_DP 
		Var_Log_Earnings_25_60 = 0.0_DP
		do xi=1,nx
		DO age=5,40
		DO ai=1,na
		DO zi=1,nz
		DO lambdai=1,nlambda
		DO ei=1,ne
			tothours_25_60 = tothours_25_60 + DBN1(age, ai, zi, lambdai, ei, xi)  * Hours(age, ai, zi, lambdai,ei,xi)
			pop_25_60      = pop_25_60 +  DBN1(age, ai, zi, lambdai, ei,xi)
			IF (Hours(age, ai, zi, lambdai, ei,xi) .ge. 0.055) THEN
			tot_log_earnings_25_60 = tot_log_earnings_25_60 + DBN1(age, ai, zi, lambdai, ei,xi)  &
			                 		& *  log( Y_h(Hours(age, ai, zi, lambdai, ei,xi),age,lambdai,ei,wage) )
			pop_pos_earn_25_60     = pop_pos_earn_25_60 +  DBN1(age, ai, zi, lambdai, ei,xi)
			ENDIF
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		meanhours_25_60         = tothours_25_60 / pop_25_60
		mean_log_earnings_25_60 = tot_log_earnings_25_60 / pop_pos_earn_25_60

		DO xi=1,nx
		DO age=5,40
		DO ai=1,na
		DO zi=1,nz
		DO lambdai=1,nlambda
		DO ei=1,ne
			IF (Hours(age, ai, zi, lambdai, ei,xi) .ge. 0.055) THEN
			    Var_Log_Earnings_25_60 =  Var_Log_Earnings_25_60 + DBN1(age, ai, zi, lambdai, ei,xi)  &
			                 			& * ( log( Y_h(Hours(age, ai, zi, lambdai, ei,xi),age,lambdai,ei,wage) ) &
			                 			& -   mean_log_earnings_25_60 ) ** 2.0_DP
			ENDIF
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		Var_Log_Earnings_25_60 = Var_Log_Earnings_25_60 / pop_pos_earn_25_60
		Std_Log_Earnings_25_60 = Var_Log_Earnings_25_60 ** 0.5_DP


	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	! Income, Wealth, Returns to Capital
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
		! print*, 'Test wealth'
		! Sources of income
		Pr_mat = Profit_Matrix(R,P)
		K_mat  = K_Matrix(R,P)
		do zi=1,nz
		do xi=1,nx 
		do ai=1,na 
			YGRID(ai,zi,xi) = Y_a(agrid(ai),zi,xi)
		enddo 
		enddo 
		enddo
		CALL ComputeLaborUnits(Ebar, wage) 

		MeanWealth 	 = 0.0_dp
		MeanATReturn = 0.0_DP
		MeanReturn 	 = 0.0_DP
		MeanCons  	 = 0.0_DP
		Mean_Capital = 0.0_DP
		
		MeanATReturn_by_z     = 0.0_DP
		MeanReturn_by_z       = 0.0_DP
		Mean_AT_K_Return_by_z = 0.0_DP
		Mean_K_Return_by_z    = 0.0_DP
		size_by_z         	  = 0.0_DP
		Wealth_by_z 	  	  = 0.0_DP
		Capital_by_z 	  	  = 0.0_DP
		DO xi=1,nx
		DO age=1,MaxAge
		DO zi=1,nz
		DO ai=1,na
		DO lambdai=1,nlambda
		DO ei=1, ne

			if (age.lt.RetAge) then
	    	Total_Income(age,ai,zi,lambdai,ei,xi) = R*agrid(ai)+Pr_mat(ai,zi,xi) + yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei,xi)
	    	else
	    	Total_Income(age,ai,zi,lambdai,ei,xi) = R*agrid(ai)+Pr_mat(ai,zi,xi) + RetY_lambda_e(lambdai,ei) 
	    	endif 

		    MeanWealth   = MeanWealth   + DBN1(age, ai, zi, lambdai, ei, xi)*agrid(ai)
		    Mean_Capital = Mean_Capital + DBN1(age, ai, zi, lambdai, ei, xi)*K_mat(ai,zi,xi)
		    MeanCons     = MeanCons     + DBN1(age, ai, zi, lambdai, ei, xi)*cons(age, ai, zi, lambdai, ei, xi)

		    size_by_z(zi)    = size_by_z(zi)    + DBN1(age, ai, zi, lambdai, ei, xi) 
		    Wealth_by_z(zi)  = Wealth_by_z(zi)  + DBN1(age, ai, zi, lambdai, ei, xi) * agrid(ai)
		    Capital_by_z(zi) = Capital_by_z(zi) + DBN1(age, ai, zi, lambdai, ei, xi) * K_mat(ai,zi,xi)

		    MeanReturn           = MeanReturn          + DBN1(age, ai, zi, lambdai, ei, xi) * &
		    							& (R*agrid(ai) + Pr_mat(ai,zi,xi))
		    MeanReturn_by_z(zi)  = MeanReturn_by_z(zi) + DBN1(age, ai, zi, lambdai, ei, xi) * &
		    							& (R*agrid(ai) + Pr_mat(ai,zi,xi)) 
		    
		    MeanATReturn           = MeanATReturn          + DBN1(age, ai, zi, lambdai, ei, xi) * (YGRID(ai,zi,xi)-agrid(ai))
		    MeanATReturn_by_z(zi)  = MeanATReturn_by_z(zi) + DBN1(age, ai, zi, lambdai, ei, xi) * (YGRID(ai,zi,xi)-agrid(ai))
		    

		    !MeanATReturn = MeanATReturn + DBN1(age, ai, zi, lambdai, ei) * (MBGRID(ai,zi)-1.0_DP)* agrid(ai)
		    !if (K_mat(ai,zi) .lt. (theta*agrid(ai)) ) then
		    !  MeanReturn   = MeanReturn   + DBN1(age, ai, zi, lambdai, ei) * R * agrid(ai)    
		    !else
		    !  MeanReturn = MeanReturn+ DBN1(age, ai, zi, lambdai, ei) * agrid(ai) * & 
		    !   			& ( R + (P*mu*((theta*zgrid(zi))**mu)*(agrid(ai))**(mu-1.0_DP)-(R+DepRate)*theta)) 
		    !endif      
		ENDDO
		ENDDO
		ENDDO
		ENDDO    
		ENDDO    
		ENDDO

			! Allocate MeanReturn to K
			Mean_K_Return_by_z    = MeanReturn_by_z
			Mean_AT_K_Return_by_z = MeanATReturn_by_z     
		Wealth_Output = MeanWealth/YBAR 
		MeanReturn    = MeanReturn/MeanWealth
		MeanATReturn  = MeanATReturn/MeanWealth
		MeanATReturn_by_z     = MeanATReturn_by_z / Wealth_by_z
	    MeanReturn_by_z       = MeanReturn_by_z   / Wealth_by_z
	    Mean_AT_K_Return_by_z = Mean_AT_K_Return_by_z / Capital_by_z
	    Mean_K_Return_by_z    = Mean_K_Return_by_z   / Capital_by_z


		VarATReturn = 0.0_DP
		VarReturn 	= 0.0_DP
		Var_AT_K_Return = 0.0_DP
		Var_K_Return 	= 0.0_DP
		DO xi=1,nx
		DO age=1,MaxAge
		DO zi=1,nz
		DO ai=1,na
		DO lambdai=1,nlambda
		DO ei=1, ne  

		    VarReturn    = VarReturn +  DBN1(age, ai, zi, lambdai, ei, xi) * agrid(ai)/MeanWealth * &
		    				& ((R*agrid(ai) + Pr_mat(ai,zi,xi))/agrid(ai)-MeanReturn)**2.0_dp

		    Var_K_Return = Var_K_Return +  DBN1(age, ai, zi, lambdai, ei, xi) * K_mat(ai,zi,xi)/MeanWealth * &
		    				& ((R*agrid(ai) + Pr_mat(ai,zi,xi))/K_mat(ai,zi,xi)-MeanATReturn)**2.0_dp

		    VarATReturn  = VarATReturn +  DBN1(age, ai, zi, lambdai, ei, xi) * agrid(ai)/MeanWealth * &
		    				& ((YGRID(ai,zi,xi)-agrid(ai))/agrid(ai)-MeanATReturn)**2.0_dp 

		    Var_AT_K_Return  = Var_AT_K_Return +  DBN1(age, ai, zi, lambdai, ei, xi) * K_mat(ai,zi,xi)/MeanWealth * &
		    				& ((YGRID(ai,zi,xi)-agrid(ai))/K_mat(ai,zi,xi)-MeanATReturn)**2.0_dp 
		    
		    !VarATReturn = VarATReturn + DBN1(age, ai, zi, lambdai, ei) * agrid(ai)/MeanWealth * &
		    !				& ((MBGRID(ai,zi)-1.0_DP)-MeanATReturn)**2.0_dp

		    !if (K_mat(ai,zi) .lt. (theta*agrid(ai)) ) then
			!    VarReturn = VarReturn   + DBN1(age, ai, zi, lambdai, ei) * agrid(ai)/MeanWealth * (R-MeanReturn)**2.0_dp
		   	!else
			!	VarReturn = VarReturn+ DBN1(age, ai, zi, lambdai, ei) * agrid(ai)/MeanWealth * & 
			!				& (( R + (P*mu*((theta*zgrid(zi))**mu)*(agrid(ai))**(mu-1.0_DP)-(R+DepRate)*theta)) -MeanReturn)**2.0_dp
		   	!endif  

		ENDDO
		ENDDO
		ENDDO
		ENDDO    
		ENDDO    
		ENDDO  
		StdATReturn     = VarATReturn**0.5_DP
		StdReturn       = VarReturn**0.5_DP
		Std_AT_K_Return = Var_AT_K_Return**0.5_DP
		Std_K_Return    = Var_K_Return**0.5_DP


	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	! Distribution of Returns (this recycles bequest variables)
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------ 
		DBN_azx = sum(sum(sum(DBN1,5),4),1) 
		do xi=1,nx 
		do zi=1,nz
		do ai=1,na 
			BT_Return(ai,zi,xi)    = 100.0_dp*(R+Pr_mat(ai,zi,xi)/agrid(ai))
		enddo 
		enddo 
		enddo 

		! Vectorizations
		DBN_azx_vec = reshape(DBN_azx   ,(/size(DBN_azx)/)); 
		Return_vec  = reshape(BT_Return ,(/size(DBN_azx)/)); 

		! Compute bequest by percentile (percentiles for counter CDF)
		prctile_bq = (/0.9_dp, 0.50_dp, 0.10_dp, 0.05_dp, 0.01_dp/)
		a = minval(Return_vec)
		b = maxval(Return_vec) 
		c = a
		do i=1,size(prctile_bq)
			a = c
			b = maxval(Return_vec)
			c = (a+b)/2.0_dp
			CCDF_c = sum(DBN_azx_vec,Return_vec>=c)
			!print*, ' '
			!print*, 'Percentile', prctile_bq(i)
			do while ((abs(CCDF_c-prctile_bq(i))>0.00001_dp).and.(b-a>1e-9))
				if (CCDF_c<prctile_bq(i)) then 
					b = c 
					c = (a+b)/2.0_dp
					CCDF_c = sum(DBN_azx_vec,Return_vec>=c)
				else 
					a = c 
					c = (a+b)/2.0_dp
					CCDF_c = sum(DBN_azx_vec,Return_vec>=c)
				endif
				! print*, 'a',a,'c',c,'b',b,'CCDF',CCDF_c,'obj',prctile_bq(i),'Error', abs(CCDF_c-prctile_bq(i))
			enddo 
			BQ_top_x(i) = c 
		enddo 
		if (solving_bench.eq.1) then
			OPEN(UNIT=11, FILE=trim(Result_Folder)//'Return_Pct_Bench.txt', STATUS='replace')
		else
			OPEN(UNIT=11, FILE=trim(Result_Folder)//'Return_Pct_Exp.txt', STATUS='replace')
		end if 
		print*,' '
		print*,'-----------------------------------------------------'
		WRITE(UNIT=11, FMT=*) 'Return Percentiles'
		WRITE(UNIT=11, FMT=*) 'Tax p10 p50 p90 p95 p99'
		WRITE(UNIT=11, FMT=*) 'Before_Tax',BQ_top_x
		if (solving_bench.eq.1) then 
		WRITE(UNIT=11, FMT=*) 'Before_Tax',BQ_top_x*(1.0_dp-tauK)
		else
		WRITE(UNIT=11, FMT=*) 'After_Tax',BQ_top_x-100.0_dp*tauW_at
		endif 
		CLOSE(UNIT=11)
		print*,' Return Percentiles'
		print '(A,X,X,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3)',&
			& ' 	p10',BQ_top_x(1),'p50',BQ_top_x(2),'p90',BQ_top_x(3),'p95',BQ_top_x(4),'p99',BQ_top_x(5)
		print*,'-----------------------------------------------------'; print*, ' '

	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	! Returns and Wealth 
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------ 
		! Compute average returns (unweighted) for bins of the wealth distribution
		! Bins chosen as 0-10%, 10-20%, 20-30%, ...80-90%, 90-95%, 95-99%, 99%+ (12 bins)
		pct_graph_lim = (/0.0_dp, 10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp, 70.0_dp,&
								& 80.0_dp, 90.0_dp, 95.0_dp, 99.0_dp, 99.90_dp, 100.0_dp/)

		do age=1,draft_age_category+1

		! Distribution of assets by age group
		if (age.le.draft_age_category) then
			! Select ages in age group 
			pr_a_dbn = sum(sum(sum(sum(sum(DBN1(draft_age_limit(age)+1:draft_age_limit(age+1),:,:,:,:,:),6),5),4),3),1) 
			! Distribution of returns for age group
			DBN_azx  = sum(sum(sum(DBN1(draft_age_limit(age)+1:draft_age_limit(age+1),:,:,:,:,:),5),4),1) 
		else 
			! All ages 
			pr_a_dbn = sum(sum(sum(sum(sum(DBN1,6),5),4),3),1) 
			! Distribution of returns for age group
			DBN_azx  = sum(sum(sum(DBN1,5),4),1) 
		endif
		 	pr_a_dbn = pr_a_dbn/sum(pr_a_dbn)
		 	DBN_azx  = DBN_azx/sum(DBN_azx)
		
		! CDF of assets 
		do ai=1,na
		    cdf_a_dbn(ai)         = sum( pr_a_dbn(1:ai) )      
		enddo 
		cdf_a_dbn = cdf_a_dbn + 1.0_DP - cdf_a_dbn(na)

		! Index of percentiles
		prctile_ai_ind_age(1) = 1 
		DO prctile=2,14
		    ai=1
		    DO while (cdf_a_dbn(ai) .lt. (pct_graph_lim(prctile)/100.0_DP-0.000000000000001))
		        ai=ai+1
		    ENDDO
		    prctile_ai_ind_age(prctile) = ai
		ENDDO


		! Average return by bin
		do i=1,13
			ind_lo = prctile_ai_ind_age(i  )
			ind_hi = prctile_ai_ind_age(i+1)

			pct_graph_wealth(age,i) = (EBAR_data/(EBAR_bench*0.727853584919652_dp))*agrid(prctile_ai_ind_age(i+1))
			ret_by_wealth(age,i)    = sum(BT_Return(ind_lo:ind_hi,:,:)*DBN_azx(ind_lo:ind_hi,:,:))/sum(DBN_azx(ind_lo:ind_hi,:,:))
		enddo 

		enddo 

		if (solving_bench.eq.1) then 
		OPEN (UNIT=81, FILE=trim(Result_Folder)//'Returns_by_Wealth_pct.txt', STATUS='replace') 
		OPEN (UNIT=82, FILE=trim(Result_Folder)//'Wealth_pct_by_age_group.txt', STATUS='replace') 
		else 
		OPEN (UNIT=81, FILE=trim(Result_Folder)//'Returns_by_Wealth_pct.txt', STATUS='replace') 
		OPEN (UNIT=82, FILE=trim(Result_Folder)//'Wealth_pct_by_age_group.txt', STATUS='replace') 
		endif 
			WRITE (UNIT=81, FMT=*)  'Returns by Percntile of Wealth'
			WRITE (UNIT=82, FMT=*)  'Percentile of Wealth by age group'
			WRITE (UNIT=81, FMT=*)  'p10 p20 p30 p40 p50 p60 p70 p80 p90 p95 p99 p99.9 p100'
			WRITE (UNIT=82, FMT=*)  'p10 p20 p30 p40 p50 p60 p70 p80 p90 p95 p99 p99.9 p100'
		do age=1,draft_age_category+1
			WRITE (UNIT=81, FMT=*)  ret_by_wealth(age,:) 
			WRITE (UNIT=82, FMT=*)  pct_graph_wealth(age,:) 
		enddo 
		close(unit=81); close(unit=82);

	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	! Distribution of bequest
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
		Bequest_Wealth=0.0_DP
		DO xi=1,nx
		DO zi=1,nz
		DO ai=1,na
		DO lambdai=1,nlambda
		DO ei=1, ne
		   Bequest_Wealth = Bequest_Wealth + DBN1(1, ai, zi, lambdai, ei, xi)*agrid(ai)/((1.0_dp-tau_bq)*(1.0_dp-bq_fee))
		ENDDO
		ENDDO
		ENDDO    
		ENDDO 
		ENDDO  

		! Distribution of bequest (matrix)	
		do ai=1,MaxAge
			DBN_bq(ai,:,:,:,:,:) = DBN1(ai,:,:,:,:,:)*(1.0_DP-survP(ai))
		enddo 
		DBN_bq = DBN_bq/sum(DBN_bq)

		! Vectorizations
		DBN_bq_vec        = reshape(DBN_bq      ,(/size(DBN1)/)); 
		BQ_vec            = reshape(Aprime      ,(/size(DBN1)/)); 
		Inc_vec 		  = reshape(Total_Income,(/size(DBN1)/)); 

		! Mean Bequest
		Mean_Bequest      = sum(BQ_vec*DBN_bq_vec)

		if (solving_bench.eq.1) then
			OPEN(UNIT=11, FILE=trim(Result_Folder)//'Bequest_Stats_Bench.txt', STATUS='replace')
		else
			OPEN(UNIT=11, FILE=trim(Result_Folder)//'Bequest_Stats_Exp.txt', STATUS='replace')
		end if 
		print*,' '
		print*,'-----------------------------------------------------'
		print*,' Bequest Stats'
		print*,' '
			WRITE(UNIT=11, FMT=*) ' '
			WRITE(UNIT=11, FMT=*) 'Bequest_Stats'
			WRITE(UNIT=11, FMT=*) 'Total_Bequest/Wealth= '		, Bequest_Wealth/MeanWealth 
			WRITE(UNIT=11, FMT=*) 'Mean_Bequest/Wealth= '		, Mean_Bequest/MeanWealth 
			WRITE(UNIT=11, FMT=*) 'Mean_Bequest/PV_Wealth= '	, Bequest_Wealth/Mean_Firm_Wealth 
			WRITE(UNIT=11, FMT=*) 'Bequests_Above_Threshold= '	, Threshold_Share_bq
			WRITE(UNIT=11, FMT=*) 'Bequest_Revenue/YBAR= '		, 0
			WRITE(UNIT=11, FMT=*) 'Prctile ','Bequest ','Bq/EBAR ','Bq/Inc 0.5%','Bq_Inc 1% ','Bq_Inc 2% '
		print '(A,F7.3)', ' 	Total_Bequest/Wealth= '		, 100.0_dp*Bequest_Wealth/MeanWealth 
		print '(A,F7.3)', ' 	Mean_Bequest/Wealth= '		, 100.0_dp*Mean_Bequest/MeanWealth 
		print*, ' 	Prctile  ','Bequest  ','Bq/EBAR  ','Bq/Inc 0.5%  ','Bq_Inc 1%  ','Bq_Inc 2%  '

		! Compute bequest by percentile (percentiles for counter CDF)
		prctile_bq = (/0.4_dp, 0.25_dp, 0.10_dp, 0.05_dp, 0.01_dp/)
		a = minval(BQ_vec)
		b = maxval(BQ_vec) 
		c = a
		do i=1,size(prctile_bq)
			a = c
			b = maxval(BQ_vec)
			c = (a+b)/2.0_dp
			CCDF_c = sum(DBN_bq_vec,BQ_vec>=c)
			!print*, ' '
			!print*, 'Percentile', prctile_bq(i)
			do while ((abs(CCDF_c-prctile_bq(i))>0.00001_dp).and.(b-a>1e-9))
				if (CCDF_c<prctile_bq(i)) then 
					b = c 
					c = (a+b)/2.0_dp
					CCDF_c = sum(DBN_bq_vec,BQ_vec>=c)
				else 
					a = c 
					c = (a+b)/2.0_dp
					CCDF_c = sum(DBN_bq_vec,BQ_vec>=c)
				endif
				! print*, 'a',a,'c',c,'b',b,'CCDF',CCDF_c,'obj',prctile_bq(i),'Error', abs(CCDF_c-prctile_bq(i))
			enddo 
			BQ_top_x(i) = c 
			do j=1,3
				! Get low end of range 	
				low_pct  = prctile_bq(i)+0.005_dp*(2**(j-1))
				if (low_pct<1.0_dp) then
					a = minval(BQ_vec)
					b = c
					c_low = c
					CCDF_c = sum(DBN_bq_vec,BQ_vec>=c)
					do while ((abs(CCDF_c-low_pct)>0.00001_dp).and.(b-a>1e-9))
						if (CCDF_c<low_pct) then 
							b = c_low
						else 
							a = c_low 
						endif
						c_low = (a+b)/2.0_dp
						CCDF_c = sum(DBN_bq_vec,BQ_vec>=c_low)
					enddo 
					! print*,'		Bisection results'
					! print*, '		a',a,'c',c_low,'b',b,'CCDF',CCDF_c,'Obj',low_pct,'Error', abs(CCDF_c-low_pct)
				else
					c_low = minval(BQ_vec)
				endif 
				! Get low end of range 	
				high_pct = prctile_bq(i)-0.005_dp*(2**(j-1))
				if (high_pct>0.0_dp) then
					a = c
					b = maxval(BQ_vec)
					c_high = c
					CCDF_c = sum(DBN_bq_vec,BQ_vec>=c)
					do while ((abs(CCDF_c-high_pct)>0.00001_dp).and.(b-a>1e-9))
						if (CCDF_c<high_pct) then 
							b = c_high
						else 
							a = c_high 
						endif
						c_high = (a+b)/2.0_dp
						CCDF_c = sum(DBN_bq_vec,BQ_vec>=c_high)
					enddo 
					! print*,'		Bisection results'
					! print*, '		a',a,'c',c_high,'b',b,'CCDF',CCDF_c,'Obj',high_pct,'Error', abs(CCDF_c-high_pct)
				else
					c_high = maxval(BQ_vec)
				endif 
				! print*, ' low_pct=',low_pct,'pct=',prctile_bq(i),'high_pct=',high_pct
				! print*, ' Test:','pct=',prctile_bq(i),'c_low=',c_low,'c=',c,'c_high=',c_high

				! Get Average Bequest/Income
				Bq_Inc(i,j) = sum( (BQ_vec/Inc_vec*DBN_bq_vec) , ((BQ_vec>=c_low).and.(BQ_vec<=c_high)) )&
							&  				/sum( (DBN_bq_vec) , ((BQ_vec>=c_low).and.(BQ_vec<=c_high)) )
				! print*, ' Test logical',count((BQ_vec>=c_low)),count((BQ_vec<=c_high)),count((BQ_vec>=c_low).and.(BQ_vec<=c_high))
				! print*, ' Test DBN',sum(DBN_bq_vec,(BQ_vec>=c_low)),sum(DBN_bq_vec,(BQ_vec<=c_high)),&
				! 		& sum(DBN_bq_vec,(BQ_vec>=c_low).and.(BQ_vec<=c_high))
				! print*, ' Test BQ',sum(BQ_vec/Inc_vec,(BQ_vec>=c_low)),sum(BQ_vec/Inc_vec,(BQ_vec<=c_high)),&
				! 		& sum(BQ_vec/Inc_vec,(BQ_vec>=c_low).and.(BQ_vec<=c_high))
				! print*, ' Test Sum', sum( (BQ_vec/Inc_vec*DBN_bq_vec) , ((BQ_vec>=c_low).and.(BQ_vec<=c_high)) )/&
				! 		& sum(DBN_bq_vec,((BQ_vec>=c_low).and.(BQ_vec<=c_high))), Bq_Inc(i,j)
			enddo 
			! Write down results 
			WRITE(UNIT=11, FMT=*) 100_dp*(1.0_dp-prctile_bq(i)),BQ_top_x(i),BQ_top_x(i)/EBAR_bench,Bq_Inc(i,:)
			print '(A,X,X,F7.3,X,X,F7.3,X,X,F7.3,X,X,F7.3,X,X,F7.3,X,X,F7.3)',&
				& ' 	', 100_dp*(1.0_dp-prctile_bq(i)),BQ_top_x(i),BQ_top_x(i)/EBAR_bench,Bq_Inc(i,:)
		enddo 
			CLOSE(UNIT=11)
			print*,'-----------------------------------------------------'; print*, ' '


	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
    ! Debt to GDP Ratio
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
	    ! print*, 'Test Debt'
	    External_Debt_GDP = 0.0_DP
		DO xi=1,nx
		DO zi=1,nz
		DO ai=1,na
			if (K_mat(ai,zi,xi).gt.agrid(ai)) then 
		    External_Debt_GDP = External_Debt_GDP + sum(DBN1(:, ai, zi, :, :,xi))*(K_mat(ai,zi,xi)-agrid(ai))
		    endif 
		ENDDO
		ENDDO
		ENDDO
		External_Debt_GDP = External_Debt_GDP / YBAR

	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	! Savings Rate
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
		! print*, 'Test Savings'
		group 	 = 1
		A_Age 	 = 0.0_dp
		A_AZ  	 = 0.0_dp 
		A_W  	 = 0.0_dp
		Ap_Age 	 = 0.0_dp
		Ap_AZ  	 = 0.0_dp 
		Ap_W  	 = 0.0_dp
		Y_Age 	 = 0.0_dp
		Y_AZ  	 = 0.0_dp 
		Y_W   	 = 0.0_dp
		DO age=1,MaxAge 

		    DO while (age.gt.age_limit(group+1))
		        group = group+1
		    ENDDO    
		 
		 	DO xi=1,nx
		    DO ai=1,na
	        DO zi=1,nz
	        DO lambdai=1,nlambda
	        DO ei=1,ne

	        	A_Age(group)      = A_Age(group)   + DBN1(age,ai,zi,lambdai,ei,xi)*agrid(ai)
	        	A_AZ(group,zi)    = A_AZ(group,zi) + DBN1(age,ai,zi,lambdai,ei,xi)*agrid(ai)

	        	Ap_Age(group)     = Ap_Age(group)   + DBN1(age,ai,zi,lambdai,ei,xi)*Aprime(age,ai,zi,lambdai,ei,xi)
	        	Ap_AZ(group,zi)   = Ap_AZ(group,zi) + DBN1(age,ai,zi,lambdai,ei,xi)*Aprime(age,ai,zi,lambdai,ei,xi)

	        	if (age.lt.RetAge) then
	        	Y_Age(group)      = Y_Age(group)   + DBN1(age,ai,zi,lambdai,ei,xi)*&
	        						& ( YGRID(ai,zi,xi)+ Y_h(Hours(age,ai,zi,lambdai,ei,xi),age,lambdai,ei,Wage))
	        	Y_AZ(group,zi)    = Y_AZ(group,zi) + DBN1(age,ai,zi,lambdai,ei,xi)*&
	        						& ( YGRID(ai,zi,xi)+ Y_h(Hours(age,ai,zi,lambdai,ei,xi),age,lambdai,ei,Wage))
	        	else 
	        	Y_Age(group)      = Y_Age(group)   + DBN1(age,ai,zi,lambdai,ei,xi)*( YGRID(ai,zi,xi)+ RetY_lambda_e(lambdai,ei) )
	        	Y_AZ(group,zi)    = Y_AZ(group,zi) + DBN1(age,ai,zi,lambdai,ei,xi)*( YGRID(ai,zi,xi)+ RetY_lambda_e(lambdai,ei) )
	        	endif

	        	if (ai.le.prctile_ai_ind(90)) then 
	        		A_W(1)    = A_W(1)    + DBN1(age,ai,zi,lambdai,ei,xi)*agrid(ai)
	        		Ap_W(1)   = Ap_W(1)   + DBN1(age,ai,zi,lambdai,ei,xi)*Aprime(age,ai,zi,lambdai,ei,xi)
	        		if (age.lt.RetAge) then 
	        		Y_W(1)    = Y_W(1)    + DBN1(age,ai,zi,lambdai,ei,xi)*&
	        					& (YGRID(ai,zi,xi)+ Y_h(Hours(age,ai,zi,lambdai,ei,xi),age,lambdai,ei,Wage))
	        		else 
	        		Y_W(1)    = Y_W(1)    + DBN1(age,ai,zi,lambdai,ei,xi)*(YGRID(ai,zi,xi)+ RetY_lambda_e(lambdai,ei))
	        		endif 
	        	else if  ((ai.gt.prctile_ai_ind(90)).and.(ai.le.prctile_ai_ind(99))) then
	        		A_W(2)    = A_W(2)    + DBN1(age,ai,zi,lambdai,ei,xi)*agrid(ai)
	        		Ap_W(2)   = Ap_W(2)   + DBN1(age,ai,zi,lambdai,ei,xi)*Aprime(age,ai,zi,lambdai,ei,xi)
	        		if (age.lt.RetAge) then 
	        		Y_W(2)    = Y_W(2)    + DBN1(age,ai,zi,lambdai,ei,xi)*&
	        					& (YGRID(ai,zi,xi)+ Y_h(Hours(age,ai,zi,lambdai,ei,xi),age,lambdai,ei,Wage))
	        		else 
	        		Y_W(2)    = Y_W(2)    + DBN1(age,ai,zi,lambdai,ei,xi)*&
	        					& (YGRID(ai,zi,xi)+ RetY_lambda_e(lambdai,ei))
	        		endif
	        	else 
	        		A_W(3)    = A_W(3)    + DBN1(age,ai,zi,lambdai,ei,xi)*agrid(ai)
	        		Ap_W(3)   = Ap_W(3)   + DBN1(age,ai,zi,lambdai,ei,xi)*Aprime(age,ai,zi,lambdai,ei,xi)
	        		if (age.lt.RetAge) then 
	        		Y_W(3)    = Y_W(3)    + DBN1(age,ai,zi,lambdai,ei,xi)*&
	        					&	(YGRID(ai,zi,xi)+ Y_h(Hours(age,ai,zi,lambdai,ei,xi),age,lambdai,ei,Wage))
	        		else 
	        		Y_W(3)    = Y_W(3)    + DBN1(age,ai,zi,lambdai,ei,xi)*(YGRID(ai,zi,xi)+ RetY_lambda_e(lambdai,ei))
	        		endif
	        	endif 
	        ENDDO
	        ENDDO
	        ENDDO
	        ENDDO
		    ENDDO
		ENDDO
		S_Rate_A_Age = (Ap_Age-A_Age)/A_Age
		S_Rate_A_AZ  = (Ap_AZ-A_AZ)/A_AZ
		S_Rate_A_W   = (Ap_W-A_W)/A_W
		S_Rate_Y_Age = (Ap_Age-A_Age)/Y_Age
		S_Rate_Y_AZ  = (Ap_AZ-A_AZ)/Y_AZ
		S_Rate_Y_W   = (Ap_W-A_W)/Y_W

	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	! Leverage Ratio and fraction of constrainted firms 
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
		! print*, 'Test Leverage'
		leverage_age_z = 0.0_dp 
		size_by_age_z  = 0.0_dp 
		constrained_firms_age_z = 0.0_dp
		constrained_firm_ind = 0
		do xi=1,nx
		do age = 1,MaxAge 
		do zi  = 1,nz 
	        DO lambdai=1,nlambda
	        DO ei=1,ne
	        DO ai=1,na
	        	size_by_age(age)       = size_by_age(age)       + DBN1(age,ai,zi,lambdai,ei,xi)
				size_by_age_z(age,zi)  = size_by_age_z(age,zi)  + DBN1(age,ai,zi,lambdai,ei,xi)
				leverage_age_z(age,zi) = leverage_age_z(age,zi) + DBN1(age,ai,zi,lambdai,ei,xi)*K_mat(ai,zi,xi)/agrid(ai)
				if (K_mat(ai,zi,xi).ge.(theta(zi)*agrid(ai))) then 
					constrained_firms_age(age)      = constrained_firms_age(age)      + DBN1(age,ai,zi,lambdai,ei,xi)
					constrained_firms_age_z(age,zi) = constrained_firms_age_z(age,zi) + DBN1(age,ai,zi,lambdai,ei,xi)
					constrained_firm_ind(age,ai,zi,lambdai,ei,xi) = 1
				endif 
				! Firm_Output(age,ai,zi,lambdai,ei,xi) = xz_grid(xi,zi)*K_mat(ai,zi,xi)
				! Firm_Profit(age,ai,zi,lambdai,ei,xi) = Pr_mat(ai,zi,xi)
			enddo
			enddo 
			enddo 
		enddo 
		enddo 
		enddo 
		leverage_age_z = leverage_age_z/size_by_age_z
		constrained_firms_age_z = constrained_firms_age_z/size_by_age_z 
		constrained_firms_age   = constrained_firms_age/size_by_age 


		if (solving_bench.eq.1) then
			OPEN(UNIT=11, FILE=trim(Result_Folder)//'Constrained_firms_stats.txt', STATUS='replace')
			WRITE(UNIT=11, FMT=*) ' '
			WRITE(UNIT=11, FMT=*) 'Z ','Const_firms_by_z: ','Const_firms_z_x1 ','Const_firms_z_x1 ','Opt_K_x_1 ','Opt_K_x2 '
			WRITE(UNIT=11, FMT=*) 'Benchmark'
		else
			OPEN(UNIT=11, FILE=trim(Result_Folder)//'Constrained_firms_stats_exp.txt', STATUS='replace') 
			WRITE(UNIT=11, FMT=*) ' '
			WRITE(UNIT=11, FMT=*) 'Z ','Const_firms_by_z: ','Const_firms_z_x1 ','Const_firms_z_x1 ','Opt_K_x_1 ','Opt_K_x2 '
			WRITE(UNIT=11, FMT=*) 'Tax_Reform'
		end if 
			WRITE(UNIT=11, FMT=*) ' '
		do zi=1,nz
			WRITE(UNIT=11, FMT=*) zi, & 
				100.0_dp*sum(constrained_firm_ind(:,:,zi,:,:,:)*DBN1(:,:,zi,:,:,:))/sum(DBN1(:,:,zi,:,:,:)), &
				100.0_dp*sum(constrained_firm_ind(:,:,zi,:,:,1)*DBN1(:,:,zi,:,:,1))/sum(DBN1(:,:,zi,:,:,1)), &
				100.0_dp*sum(constrained_firm_ind(:,:,zi,:,:,2)*DBN1(:,:,zi,:,:,2))/sum(DBN1(:,:,zi,:,:,2)), &
				(EBAR_data/(EBAR*0.727853584919652_dp))*(mu*P*xz_grid(1,zi)**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) , & 
				(EBAR_data/(EBAR*0.727853584919652_dp))*(mu*P*xz_grid(2,zi)**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu))  
		enddo 	
			WRITE(UNIT=11, FMT=*) 'Total', 100.0_dp*sum(constrained_firm_ind*DBN1)

			CLOSE(UNIT=11)

			! deallocate(Firm_Output,Firm_Profit)

	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	! Distribution of firm wealth
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
		! print*, 'Test Firm Wealth'
		do ai=1,na
			Firm_Wealth(:,ai,:,:,:,:) = V_Pr(:,ai,:,:,:,:) + (1.0_dp+R)*agrid(ai)
		enddo
		Mean_Firm_Wealth = sum(Firm_Wealth*DBN1)

		DBN_vec         = reshape(DBN1       ,(/size(DBN1)/))
		Firm_Wealth_vec = reshape(Firm_Wealth,(/size(DBN1)/))

		if (solving_bench.eq.1) then
			OPEN(UNIT=11, FILE=trim(Result_Folder)//'Firm_Wealth_Bench.txt', STATUS='replace')
		else
			OPEN(UNIT=11, FILE=trim(Result_Folder)//'Firm_Wealth_Exp.txt', STATUS='replace')
		end if 
			WRITE(UNIT=11, FMT=*) ' '
			WRITE(UNIT=11, FMT=*) 'Firm_Wealth_Stats'
			WRITE(UNIT=11, FMT=*) 'Mean_Firm_Wealth= ', Mean_Firm_Wealth
			WRITE(UNIT=11, FMT=*) 'Top_x% ','x_percentile ','wealth_share_above_x ', 'Counter_CDF'

		prctile_FW = (/0.40_DP, 0.20_dp, 0.10_dp, 0.01_dp, 0.001_dp, 0.0001_dp/)
		a = minval(Firm_Wealth_vec)
		b = maxval(Firm_Wealth_vec) 
		c = a
		do i=1,size(prctile_FW)
			a = c
			b = maxval(Firm_Wealth_vec)
			c = (a+b)/2.0_dp
			CCDF_c = sum(DBN_vec,Firm_Wealth_vec>=c)
			! print*, ' '
			!print*, 'Percentile', prctile_FW(i)
			do while ((abs(CCDF_c-prctile_FW(i))>0.0001_dp).and.(b-a>1e-8))
				if (CCDF_c<prctile_FW(i)) then 
					b = c 
					c = (a+b)/2.0_dp
					CCDF_c = sum(DBN_vec,Firm_Wealth_vec>=c)
				else 
					a = c 
					c = (a+b)/2.0_dp
					CCDF_c = sum(DBN_vec,Firm_Wealth_vec>=c)
				endif
				!print*, 'a',a,'c',c,'b',b,'CCDF',CCDF_c,'Error', CCDF_c-prctile_FW(i)
			enddo 
			FW_top_x(i)       = c 
			FW_top_x_share(i) = 100.0_dp*sum(Firm_Wealth_vec*DBN_vec,Firm_Wealth_vec>=c)/Mean_Firm_Wealth
			WRITE(UNIT=11, FMT=*) 100_dp*prctile_FW(i),FW_top_x(i),FW_top_x_share(i), CCDF_c
		enddo 

			CLOSE(UNIT=11)

	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	! Frisch Elasticity 
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
		! This is only for agents with positive hours worked
		Frisch_Elasticity = 0.0_dp
		Size_Frisch       = 0.0_dp 
		Hours_Frisch	  = 0.0_dp
		DO xi=1,nx
		DO ei=1, ne
		DO lambdai=1,nlambda
		DO zi=1,nz
		DO ai=1,na
		DO age=1,RetAge-1
		if (HOURS(age,ai,zi,lambdai,ei,xi).gt.0.001_dp) then 
		Size_Frisch = Size_Frisch + DBN1(age,ai,zi,lambdai,ei,xi)
		Frisch_Elasticity = Frisch_Elasticity + DBN1(age,ai,zi,lambdai,ei,xi)*(1.0_dp-tauPL)/ &
		& ( sigma/(1.0_dp-(1.0_dp-sigma)*gamma) * HOURS(age,ai,zi,lambdai,ei,xi)/(1-HOURS(age,ai,zi,lambdai,ei,xi)) - tauPL )
		Hours_Frisch = Hours_Frisch + DBN1(age,ai,zi,lambdai,ei,xi)*HOURS(age,ai,zi,lambdai,ei,xi)
		endif 
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		Frisch_Elasticity   = Frisch_Elasticity/Size_Frisch
		Hours_Frisch 	    = Hours_Frisch/Size_Frisch
		Frisch_Elasticity_2 = (1.0_dp-tauPL)/( sigma/(1.0_dp-(1.0_dp-sigma)*gamma) * Hours_Frisch/(1-Hours_Frisch) - tauPL )
		print*,' '; print*,'-----------------------------------------------------';
		print*,'	Frisch_Elasticity'
		print '(A,F7.3)',' Leisure Elasticity (LT) = ',Frisch_Elasticity
		print '(A,F7.3)',' Leisure Elasticity (PT) = ',Frisch_Elasticity_2
		print '(A,F7.3)','   Hours Elasticity (LT) = ',Hours_Frisch
		print '(A,F7.3)','    Share of pop H>0     = ',Size_Frisch
		print*,'-----------------------------------------------------'; print*, ' '


	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	! Draft Tables 1: Aggregates 
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
		size_draft_group_z    = 0.0_dp
		wealth_draft_group_z  = 0.0_dp; capital_draft_group_z = 0.0_dp ;
		Cons_draft_group_z    = 0.0_dp; Hours_draft_group_z = 0.0_dp   ; Ap_draft_group_z 	 = 0.0_dp ;
		do zi  = 1,nz
		do age = 1,draft_age_category
			size_draft_group_z(age,zi) = sum(DBN1(draft_age_limit(age)+1:draft_age_limit(age+1),:,zi,:,:,:))

	        do xi=1,nx
		    do ei=1,ne
		    do lambdai=1,nlambda
		    do ai=1,na
		    do age2=draft_age_limit(age)+1,draft_age_limit(age+1)
		    	! Wealth and Private Capital by group
	        	 wealth_draft_group_z(age,zi) =  wealth_draft_group_z(age,zi) + agrid(ai)      *DBN1(age2,ai,zi,lambdai,ei,xi)
	        	capital_draft_group_z(age,zi) = capital_draft_group_z(age,zi) + K_mat(ai,zi,xi)*DBN1(age2,ai,zi,lambdai,ei,xi)

	        	! Consumption, Labor and Savings by group
	        	 Cons_draft_group_z(age,zi) =  Cons_draft_group_z(age,zi) + &
	        	 								&  Cons(age2,ai,zi,lambdai,ei,xi)*DBN1(age2,ai,zi,lambdai,ei,xi)
	        	Hours_draft_group_z(age,zi) = Hours_draft_group_z(age,zi) + &
	        									& Hours(age2,ai,zi,lambdai,ei,xi)*DBN1(age2,ai,zi,lambdai,ei,xi)
	        	   Ap_draft_group_z(age,zi) =    Ap_draft_group_z(age,zi) + &
	        	   								&Aprime(age2,ai,zi,lambdai,ei,xi)*DBN1(age2,ai,zi,lambdai,ei,xi)
	    	enddo 
	    	enddo 
	    	enddo 
	    	enddo 
	    	enddo  
	    	! print*,'age',age,'zi',zi,'wealth',wealth_draft_group_z(age,zi),'size',size_draft_group_z(age,zi)
		enddo
		enddo 

		DBN_Z = sum(sum(sum(sum(sum(DBN1,6),5),4),2),1) 
		do zi=1,nz 
			CDF_Z(zi) = sum(DBN_Z(1:zi))
		enddo 
		! print*,' '
		! print*,'DBN_Z=',DBN_Z
		! print*,'CDF_Z=',CDF_Z 
		! print*,' '

		! Size of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		size_draft_group = Draft_Table(size_draft_group_z,DBN_z,.true.)

		! Wealth of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		wealth_draft_group = Draft_Table(wealth_draft_group_z,DBN_z,.true.)

		! Capital of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		capital_draft_group = Draft_Table(capital_draft_group_z,DBN_z,.true.)

		! Cons, Hours, Ap of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		 Cons_draft_group = Draft_Table(Cons_draft_group_z,DBN_z,.true.)*(EBAR_data/(EBAR_bench*0.727853584919652_dp))/size_draft_group
		Hours_draft_group = Draft_Table( Hours_draft_group_z,DBN_z,.true.)/size_draft_group
		   Ap_draft_group = Draft_Table(    Ap_draft_group_z,DBN_z,.true.)

		! Fix fractions
	    av_wealth_draft_group        = (EBAR_data/(EBAR_bench*0.727853584919652_dp))*wealth_draft_group/size_draft_group
	    av_capital_draft_group       = (EBAR_data/(EBAR_bench*0.727853584919652_dp))*capital_draft_group/size_draft_group
	    av_Ap_draft_group            = (EBAR_data/(EBAR_bench*0.727853584919652_dp))*Ap_draft_group/size_draft_group
	    frac_wealth_draft_group      = 100.0_dp*wealth_draft_group/sum(wealth_draft_group(:,1:6))
	    frac_capital_draft_group     = 100.0_dp*capital_draft_group/sum(capital_draft_group(:,1:6))
	    frac_Ap_draft_group      	 = 100.0_dp*Ap_draft_group/sum(Ap_draft_group(:,1:6))


	    ! Write results in file
	    OPEN (UNIT=81, FILE=trim(Result_Folder)//'draft_group_size.txt'			, STATUS='replace') 
	    OPEN (UNIT=83, FILE=trim(Result_Folder)//'draft_group_wealth.txt'		, STATUS='replace') 
	    OPEN (UNIT=84, FILE=trim(Result_Folder)//'draft_group_wealth_av.txt'	, STATUS='replace') 
	    OPEN (UNIT=85, FILE=trim(Result_Folder)//'draft_group_wealth_frac.txt'	, STATUS='replace') 
	    OPEN (UNIT=86, FILE=trim(Result_Folder)//'draft_group_capital.txt'		, STATUS='replace') 
	    OPEN (UNIT=87, FILE=trim(Result_Folder)//'draft_group_capital_av.txt'	, STATUS='replace') 
	    OPEN (UNIT=88, FILE=trim(Result_Folder)//'draft_group_capital_frac.txt'	, STATUS='replace') 
	    OPEN (UNIT=89, FILE=trim(Result_Folder)//'draft_group_cons.txt'			, STATUS='replace') 
	    OPEN (UNIT=91, FILE=trim(Result_Folder)//'draft_group_hours.txt'		, STATUS='replace') 
	    OPEN (UNIT=92, FILE=trim(Result_Folder)//'draft_group_savings.txt'		, STATUS='replace') 
	    OPEN (UNIT=93, FILE=trim(Result_Folder)//'draft_group_savings_av.txt'	, STATUS='replace') 
	    OPEN (UNIT=94, FILE=trim(Result_Folder)//'draft_group_savings_frac.txt'	, STATUS='replace') 
		do age = 1,draft_age_category
		    WRITE  (UNIT=81, FMT=*)  size_draft_group(age,:)

		    WRITE  (UNIT=83, FMT=*)  wealth_draft_group(age,:)
		    WRITE  (UNIT=84, FMT=*)  av_wealth_draft_group(age,:)
		    WRITE  (UNIT=85, FMT=*)  frac_wealth_draft_group(age,:)
		    WRITE  (UNIT=86, FMT=*)  capital_draft_group(age,:)
		    WRITE  (UNIT=87, FMT=*)  av_capital_draft_group(age,:)
		    WRITE  (UNIT=88, FMT=*)  frac_capital_draft_group(age,:)

		    WRITE  (UNIT=89, FMT=*)  Cons_draft_group(age,:)
		    WRITE  (UNIT=91, FMT=*)  Hours_draft_group(age,:)
		    WRITE  (UNIT=92, FMT=*)  Ap_draft_group(age,:)
		    WRITE  (UNIT=93, FMT=*)  av_Ap_draft_group(age,:)
		    WRITE  (UNIT=94, FMT=*)  frac_Ap_draft_group(age,:)

		ENDDO
		close(unit=81)
		close(unit=83); close(unit=84); close(unit=85); close(unit=86); close(unit=87); close(unit=88);
		close(unit=89); close(unit=91); close(unit=92); close(unit=93); close(unit=94);


	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	! Draft Tables 2: Income and Taxes 
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
		Pr_mat = Profit_Matrix(R,P)
		K_mat  = K_Matrix(R,P)
		CALL ComputeLaborUnits(EBAR,wage)

		K_Tax_draft_group_z			= 0.0_dp
		L_Tax_draft_group_z			= 0.0_dp
		T_Inc_draft_group_z			= 0.0_dp
		K_Inc_draft_group_z         = 0.0_dp
		L_Inc_draft_group_z         = 0.0_dp
		av_K_Inc_draft_group_z    	= 0.0_dp
		av_L_Inc_draft_group_z    	= 0.0_dp
		K_Tax_Inc_draft_group_z		= 0.0_dp
		L_Tax_Inc_draft_group_z		= 0.0_dp
		Return_draft_group_z 		= 0.0_dp
		Return_AT_draft_group_z     = 0.0_dp
		Entrepreneur_10_draft_group = 0.0_dp
		Entrepreneur_50_draft_group = 0.0_dp
		Entrepreneur_10				= 0.0_dp
		Entrepreneur_50				= 0.0_dp

		do zi  = 1,nz
		do age = 1,draft_age_category
	        do xi=1,nx
		    do ei=1,ne
		    do lambdai=1,nlambda
		    do ai=1,na
		    do age2=draft_age_limit(age)+1,draft_age_limit(age+1)

		    	! Income for each agent 
		    	if (age2.lt.RetAge) then
		    	L_Inc_aux   = yh(age2,lambdai,ei)*Hours(age2,ai,zi,lambdai,ei,xi)
		    	! if (solving_bench.ne.1) then 
		    	! print*, L_Inc_aux,yh(age2,lambdai,ei),Hours(age2,ai,zi,lambdai,ei,xi)
		    	! endif 
		    	else
		    	L_Inc_aux   = RetY_lambda_e(lambdai,ei) 
		    	endif 
		    	K_Inc_aux   = R*agrid(ai) + Pr_mat(ai,zi,xi)

		    	! Income by group (total, capital, labor)
	    		T_Inc_draft_group_z(age,zi) = T_Inc_draft_group_z(age,zi) + (K_Inc_aux + L_Inc_aux)*DBN1(age2,ai,zi,lambdai,ei,xi)
        		K_Inc_draft_group_z(age,zi) = K_Inc_draft_group_z(age,zi) + (K_Inc_aux 		  	  )*DBN1(age2,ai,zi,lambdai,ei,xi)
        		L_Inc_draft_group_z(age,zi) = L_Inc_draft_group_z(age,zi) + (L_Inc_aux 		  	  )*DBN1(age2,ai,zi,lambdai,ei,xi)

        		! Fraction of capital and labor income by agent (averaged)
        		av_K_Inc_draft_group_z(age,zi) = av_K_Inc_draft_group_z(age,zi) + & 
	        		& ( K_Inc_aux )/( K_Inc_aux + L_Inc_aux )*DBN1(age2,ai,zi,lambdai,ei,xi)

        		av_L_Inc_draft_group_z(age,zi) = av_L_Inc_draft_group_z(age,zi) + & 
	        		& ( L_Inc_aux )/( K_Inc_aux + L_Inc_aux )*DBN1(age2,ai,zi,lambdai,ei,xi)

        		! Tax by agent (capital and labor)
		    	K_Tax_aux  = K_Inc_aux - (YGRID(ai,zi,xi) - agrid(ai))
		    	if (age2.lt.RetAge) then
		    	L_Tax_aux  = L_Inc_aux - psi*(L_Inc_aux)**(1.0_DP-tauPL)
		    	else 
		    	L_Tax_aux  = 0.0_dp 
		    	endif 

		    	! Tax by group (capital and labor)
	        	K_Tax_draft_group_z(age,zi) = K_Tax_draft_group_z(age,zi) + & 
	        		& K_Tax_aux*DBN1(age2,ai,zi,lambdai,ei,xi)
	        	if (age2.lt.RetAge) then
        		L_Tax_draft_group_z(age,zi) = L_Tax_draft_group_z(age,zi) + & 
	        		& L_Tax_aux* DBN1(age2,ai,zi,lambdai,ei,xi)
	        	endif 

	        	! Tax to total income by agent (averaged)
        		K_Tax_Inc_draft_group_z(age,zi) = K_Tax_Inc_draft_group_z(age,zi) + & 
	        		& K_Tax_aux*DBN1(age2,ai,zi,lambdai,ei,xi)/( K_Inc_aux + L_Inc_aux )
	        	if (age2.lt.RetAge) then
        		L_Tax_Inc_draft_group_z(age,zi) = L_Tax_Inc_draft_group_z(age,zi) + & 
	        		& L_Tax_aux*DBN1(age2,ai,zi,lambdai,ei,xi)/( K_Inc_aux + L_Inc_aux )
	        	endif 

	        	! Return by agent (averaged)
	        	Return_draft_group_z(age,zi)    = Return_draft_group_z(age,zi) + &
	        		&                     K_Inc_aux/agrid(ai)*DBN1(age2,ai,zi,lambdai,ei,xi)
        		Return_AT_draft_group_z(age,zi) = Return_AT_draft_group_z(age,zi) + &
	        		& (YGRID(ai,zi,xi)-agrid(ai))/agrid(ai)*DBN1(age2,ai,zi,lambdai,ei,xi)

        		! Share of entrepreneurs

        		if ((R*min(agrid(ai),K_mat(ai,zi,xi))+Pr_mat(ai,zi,xi)/(K_Inc_aux + L_Inc_aux)).gt.0.10_dp) then 
        		Entrepreneur_10_draft_group_z(age,zi) = Entrepreneur_10_draft_group_z(age,zi) + DBN1(age2,ai,zi,lambdai,ei,xi)
        		Entrepreneur_10 = Entrepreneur_10 + DBN1(age2,ai,zi,lambdai,ei,xi)
        		endif 
        		if ((R*min(agrid(ai),K_mat(ai,zi,xi))+Pr_mat(ai,zi,xi)/(K_Inc_aux + L_Inc_aux)).gt.0.50_dp) then 
        		Entrepreneur_50_draft_group_z(age,zi) = Entrepreneur_50_draft_group_z(age,zi) + DBN1(age2,ai,zi,lambdai,ei,xi)
        		Entrepreneur_50 = Entrepreneur_50 + DBN1(age2,ai,zi,lambdai,ei,xi)
        		endif 

	    	enddo 
	    	enddo 
	    	enddo 
	    	enddo 
	    	enddo  
		enddo
		enddo

		! Total Capital Tax adjusted by productivity group
		K_Tax_draft_group = Draft_Table(K_Tax_draft_group_z,DBN_z,.true.)

		! Total Labor Tax adjusted by productivity group
		L_Tax_draft_group = Draft_Table(L_Tax_draft_group_z,DBN_z,.true.)

		! Total Pre-Tax Income adjusted by productivity group
		T_Inc_draft_group = Draft_Table(T_Inc_draft_group_z,DBN_z,.true.)

		! Total Pre-Tax Capital Income adjusted by productivity group
		K_Inc_draft_group = Draft_Table(K_Inc_draft_group_z,DBN_z,.true.)

		! Total Pre-Tax Labor Income adjusted by productivity group
		L_Inc_draft_group = Draft_Table(L_Inc_draft_group_z,DBN_z,.true.)


		! Get fraction of K or L income by group
		frac_K_Tax_draft_group = 100.0_dp*K_Tax_draft_group/sum(K_Tax_draft_group(:,1:6))
		frac_L_Tax_draft_group = 100.0_dp*L_Tax_draft_group/sum(L_Tax_draft_group(:,1:6))
		frac_K_Inc_draft_group = 100.0_dp*K_Inc_draft_group/sum(K_Inc_draft_group(:,1:6))
		frac_L_Inc_draft_group = 100.0_dp*L_Inc_draft_group/sum(L_Inc_draft_group(:,1:6))

		! Get ratios to total income in group
		K_Tax_draft_group = 100.0_dp*K_Tax_draft_group/T_Inc_draft_group
		L_Tax_draft_group = 100.0_dp*L_Tax_draft_group/T_Inc_draft_group
		K_Inc_draft_group = 100.0_dp*K_Inc_draft_group/T_Inc_draft_group
		L_Inc_draft_group = 100.0_dp*L_Inc_draft_group/T_Inc_draft_group

		! Divide by mass in group
		K_Tax_Inc_draft_group_z  = 100.0_dp*K_Tax_Inc_draft_group_z/size_draft_group_z
		L_Tax_Inc_draft_group_z  = 100.0_dp*L_Tax_Inc_draft_group_z/size_draft_group_z
		av_K_Inc_draft_group_z 	 = 100.0_dp*av_K_Inc_draft_group_z/size_draft_group_z 
		av_L_Inc_draft_group_z 	 = 100.0_dp*av_L_Inc_draft_group_z/size_draft_group_z 

		! Average return by productivity group
		Return_draft_group_z     = 100.0_dp*   Return_draft_group_z/size_draft_group_z
		Return_AT_draft_group_z  = 100.0_dp*Return_AT_draft_group_z/size_draft_group_z 

		! Average Capital Tax to income ratio adjusted by productivity group
		K_Tax_Inc_draft_group = Draft_Table(K_Tax_Inc_draft_group_z,DBN_z,.false.)

		! Average Labor Tax to income ratio adjusted by productivity group
		L_Tax_Inc_draft_group = Draft_Table(L_Tax_Inc_draft_group_z,DBN_z,.false.)

		! Capital Income Share Tax adjusted by productivity group
		av_K_Inc_draft_group = Draft_Table(av_K_Inc_draft_group_z,DBN_z,.false.)

		! Labor Income Share Tax adjusted by productivity group
		av_L_Inc_draft_group = Draft_Table(av_L_Inc_draft_group_z,DBN_z,.false.)

		! Average Return by productivity group
		   Return_draft_group = Draft_Table(   Return_draft_group_z,DBN_z,.false.)
		Return_AT_draft_group = Draft_Table(Return_AT_draft_group_z,DBN_z,.false.)

		! Frac. Entrepreneurs
		Entrepreneur_10_draft_group = Draft_Table(Entrepreneur_10_draft_group_z,DBN_z,.true.)
		Entrepreneur_50_draft_group = Draft_Table(Entrepreneur_50_draft_group_z,DBN_z,.true.)
			! Fix fractions
			Entrepreneur_10_draft_group = 100.0_dp*Entrepreneur_10_draft_group/size_draft_group
			Entrepreneur_50_draft_group = 100.0_dp*Entrepreneur_50_draft_group/size_draft_group

		OPEN (UNIT=80, FILE=trim(Result_Folder)//'draft_group_Tax_K.txt', STATUS='replace') 
	    OPEN (UNIT=81, FILE=trim(Result_Folder)//'draft_group_Tax_L.txt', STATUS='replace') 
	    OPEN (UNIT=82, FILE=trim(Result_Folder)//'draft_group_Tax_K_frac.txt', STATUS='replace') 
	    OPEN (UNIT=83, FILE=trim(Result_Folder)//'draft_group_Tax_L_frac.txt', STATUS='replace') 
	    OPEN (UNIT=84, FILE=trim(Result_Folder)//'draft_group_Inc.txt', STATUS='replace') 
	    OPEN (UNIT=85, FILE=trim(Result_Folder)//'draft_group_Tax_K_Inc.txt', STATUS='replace') 
	    OPEN (UNIT=86, FILE=trim(Result_Folder)//'draft_group_Tax_L_Inc.txt', STATUS='replace') 
	    OPEN (UNIT=87, FILE=trim(Result_Folder)//'draft_group_K_Inc.txt', STATUS='replace') 
	    OPEN (UNIT=88, FILE=trim(Result_Folder)//'draft_group_L_Inc.txt', STATUS='replace') 
	    OPEN (UNIT=89, FILE=trim(Result_Folder)//'draft_group_K_Inc_av.txt', STATUS='replace') 
	    OPEN (UNIT=90, FILE=trim(Result_Folder)//'draft_group_L_Inc_av.txt', STATUS='replace') 
	    OPEN (UNIT=91, FILE=trim(Result_Folder)//'draft_group_K_Inc_frac.txt', STATUS='replace') 
	    OPEN (UNIT=92, FILE=trim(Result_Folder)//'draft_group_L_Inc_frac.txt', STATUS='replace') 
	    OPEN (UNIT=93, FILE=trim(Result_Folder)//'draft_group_Return.txt', STATUS='replace') 
	    OPEN (UNIT=94, FILE=trim(Result_Folder)//'draft_group_Return_AT.txt', STATUS='replace') 
	    OPEN (UNIT=95, FILE=trim(Result_Folder)//'draft_group_Entrepreneur_10.txt', STATUS='replace') 
	    OPEN (UNIT=96, FILE=trim(Result_Folder)//'draft_group_Entrepreneur_50.txt', STATUS='replace') 
	    do age = 1,draft_age_category
		    WRITE  (UNIT=80, FMT=*)  K_Tax_draft_group(age,:)
		    WRITE  (UNIT=81, FMT=*)  L_Tax_draft_group(age,:)
		    WRITE  (UNIT=82, FMT=*)  frac_K_Tax_draft_group(age,:)
		    WRITE  (UNIT=83, FMT=*)  frac_L_Tax_draft_group(age,:)
		    WRITE  (UNIT=84, FMT=*)  T_Inc_draft_group(age,:)
		    WRITE  (UNIT=85, FMT=*)  K_Tax_Inc_draft_group(age,:)
		    WRITE  (UNIT=86, FMT=*)  L_Tax_Inc_draft_group(age,:)
		    WRITE  (UNIT=87, FMT=*)  K_Inc_draft_group(age,:)
		    WRITE  (UNIT=88, FMT=*)  L_Inc_draft_group(age,:)
		    WRITE  (UNIT=89, FMT=*)  av_K_Inc_draft_group(age,:)
		    WRITE  (UNIT=90, FMT=*)  av_L_Inc_draft_group(age,:)
		    WRITE  (UNIT=91, FMT=*)  frac_K_Inc_draft_group(age,:)
		    WRITE  (UNIT=92, FMT=*)  frac_L_Inc_draft_group(age,:)
		    WRITE  (UNIT=93, FMT=*)  Return_draft_group(age,:)
		    WRITE  (UNIT=94, FMT=*)  Return_AT_draft_group(age,:)
		    WRITE  (UNIT=95, FMT=*)  Entrepreneur_10_draft_group(age,:)
		    WRITE  (UNIT=96, FMT=*)  Entrepreneur_50_draft_group(age,:)
		ENDDO
		close(unit=80); close(unit=81); close(unit=83); close(unit=84); close(unit=85)
		close(unit=87); close(unit=88); close(unit=89); close(unit=90); close(unit=91); close(unit=92); 
		close(unit=93); close(unit=94); close(unit=95); close(unit=96);

		OPEN (UNIT=80, FILE=trim(Result_Folder)//'Entrepreneur_10_50.txt', STATUS='replace') 
		WRITE  (UNIT=80, FMT=*)  'Share of entrepreneurs '
		WRITE  (UNIT=80, FMT=*)  'Profits/Before_Tax_Income>10% ',100.0_dp*Entrepreneur_10
		WRITE  (UNIT=80, FMT=*)  'Profits/Before_Tax_Income>50% ',100.0_dp*Entrepreneur_50
		close(unit=80); 




	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	! Draft Tables 3: Income and Tax Changes
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	if (solving_bench.ne.1) then 

		! Set variables to benchmark values
		Pr_mat_bench = Profit_Matrix(R_bench,P_bench)

		Tax_Increase_tk_draft_group_z = 0.0_dp ; Tax_Increase_tl_draft_group_z = 0.0_dp ;
		Tax_Increase_draft_group_z = 0.0_dp ; Tax_Rate_Increase_draft_group_z = 0.0_dp ;
		Tax_Rate_Increase_tk_draft_group_z = 0.0_dp ; Tax_Rate_Increase_tl_draft_group_z = 0.0_dp ;
		Inc_Increase_draft_group_z = 0.0_dp  ; 
		K_Inc_Increase_draft_group_z = 0.0_dp ; L_Inc_Increase_draft_group_z = 0.0_dp;

		do zi  = 1,nz
		do age = 1,draft_age_category
	        do xi=1,nx
		    do ei=1,ne
		    do lambdai=1,nlambda
		    do ai=1,na
		    do age2=draft_age_limit(age)+1,draft_age_limit(age+1)

		    	! Income for each agent in benchmark
		    	if (age2.lt.RetAge) then
		    	L_Inc_aux     = yh(age2,lambdai,ei)*Hours(age2,ai,zi,lambdai,ei,xi)
		    	L_Inc_bench   = Wage_bench*eff_un(age2,lambdai,ei)*Hours_bench(age2,ai,zi,lambdai,ei,xi)
		    	! print*, L_Inc_aux,Wage_bench,eff_un(age2,lambdai,ei),Hours_bench(age2,ai,zi,lambdai,ei,xi)
		    	else
	    		L_Inc_aux   = RetY_lambda_e(lambdai,ei) 
		    		if (KeepSSatBench .eq. 1) then 
		    		L_Inc_bench = RetY_lambda_e(lambdai,ei) 
		    		else 
		    		L_Inc_bench = RetY_lambda_e(lambdai,ei)*EBAR_bench/EBAR
		    		endif 
		    	endif 
		    	K_Inc_aux   = R*agrid(ai) + Pr_mat(ai,zi,xi)
		    	K_Inc_bench = R_bench*agrid(ai) + Pr_mat_bench(ai,zi,xi)

		    	! print*, K_Inc_aux, K_Inc_bench, L_Inc_aux, L_Inc_bench

		    	! Tax by agent (capital and labor) in benchmark
		    	K_Tax_aux   = K_Inc_aux - (YGRID(ai,zi,xi) - agrid(ai))
		    	K_Tax_bench = tauK_bench*K_Inc_bench
		    	if (age2.lt.RetAge) then
		    	L_Tax_aux   = L_Inc_aux   - psi*(L_Inc_aux)**(1.0_DP-tauPL)
		    	L_Tax_bench = L_Inc_bench - psi_bench*(L_Inc_bench)**(1.0_DP-tauPL_bench)
		    	else 
		    	L_Tax_aux   = 0.0_dp
		    	L_Tax_bench = 0.0_dp 
		    	endif 

		    	! Compare income 
		    	if ((L_Inc_aux + K_Inc_aux).gt.(L_Inc_bench + K_Inc_bench)) then
		    	Inc_Increase_draft_group_z(age,zi) = Inc_Increase_draft_group_z(age,zi) + DBN_bench(age2,ai,zi,lambdai,ei,xi)
		    	endif 
		    	if (K_Inc_aux.gt.K_Inc_bench) then
		    	K_Inc_Increase_draft_group_z(age,zi) = K_Inc_Increase_draft_group_z(age,zi) + DBN_bench(age2,ai,zi,lambdai,ei,xi)
		    	endif 
		    	if (L_Inc_aux.gt.L_Inc_bench) then
		    	L_Inc_Increase_draft_group_z(age,zi) = L_Inc_Increase_draft_group_z(age,zi) + DBN_bench(age2,ai,zi,lambdai,ei,xi)
		    	endif 


		    	! Compare Capital taxes
		    	If ( K_Tax_aux .gt. tauK_bench*K_Inc_bench ) then
		    	Tax_Increase_tk_draft_group_z(age,zi) = Tax_Increase_tk_draft_group_z(age,zi) + DBN_bench(age2,ai,zi,lambdai,ei,xi)
		    	endif 

		    	! Compare Labor taxes
		    	If ( L_Tax_aux .gt. L_Tax_bench) then
		    	Tax_Increase_tl_draft_group_z(age,zi) = Tax_Increase_tl_draft_group_z(age,zi) + DBN_bench(age2,ai,zi,lambdai,ei,xi)
		    	endif 

		    	! Compare Total taxes
		    	If ( (K_Tax_aux+L_Tax_aux) .gt. ( tauK_bench*K_Inc_bench +L_Tax_bench) ) then
		    	Tax_Increase_draft_group_z(age,zi) = Tax_Increase_draft_group_z(age,zi) + DBN_bench(age2,ai,zi,lambdai,ei,xi)
		    	endif 

				! Compare Capital tax rate
		    	If (K_Tax_aux/K_Inc_aux .gt. tauK_bench)  then
		    	Tax_Rate_Increase_tk_draft_group_z(age,zi) = Tax_Rate_Increase_tk_draft_group_z(age,zi) + & 
		    		&  DBN_bench(age2,ai,zi,lambdai,ei,xi)
		    	endif 

		    	! Compare Labor tax rate
		    	if ((L_inc_aux.gt.0.0_dp).and.(L_inc_bench.gt.0.0_dp)) then
		    	If ( (L_Tax_aux/L_Inc_aux - L_Tax_bench/L_Inc_bench).gt.1.0E-07_dp ) then
		    	Tax_Rate_Increase_tl_draft_group_z(age,zi) = Tax_Rate_Increase_tl_draft_group_z(age,zi) + & 
		    		& DBN_bench(age2,ai,zi,lambdai,ei,xi)
		    	endif 
		    	endif 

		    	! Compare Total tax rates
		    	If ( ( ( K_Tax_aux+L_Tax_aux)/(K_Inc_aux+L_Inc_aux) &
		    	& - ( tauK_bench*K_Inc_bench + L_Tax_bench)/(K_Inc_bench + L_Inc_bench) ).gt.1.0E-07_dp ) then
		    	Tax_Rate_Increase_draft_group_z(age,zi) = Tax_Rate_Increase_draft_group_z(age,zi)+DBN_bench(age2,ai,zi,lambdai,ei,xi)
		    	endif 

		  !   	if (age2.lt.RetAge) then
		  !   	print*, 'test rates',(K_Tax_aux/K_Inc_aux),&
		  !   						& K_Tax_bench/K_Inc_bench, &
		  !   						& L_Tax_aux/L_Inc_aux, &
		  !   						& L_Tax_bench/L_Inc_bench, &
		  !   						& K_Inc_aux+L_Inc_aux,K_Inc_bench+L_Inc_bench 
				! endif
	    	enddo 
	    	enddo 
	    	enddo 
	    	enddo 
	    	enddo  
		enddo
		enddo

		! Frac. capital tax increase by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		Tax_Increase_tk_draft_group = Draft_Table(Tax_Increase_tk_draft_group_z,DBN_z,.true.)
			! Fix fractions
			Tax_Increase_tk_draft_group = 100.0_dp*Tax_Increase_tk_draft_group/size_draft_group

		! Frac. labor tax increase by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		Tax_Increase_tl_draft_group = Draft_Table(Tax_Increase_tl_draft_group_z,DBN_z,.true.)
			! Fix fractions
			Tax_Increase_tl_draft_group = 100.0_dp*Tax_Increase_tl_draft_group/size_draft_group

		! Frac. total tax increase by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		Tax_Increase_draft_group = Draft_Table(Tax_Increase_draft_group_z,DBN_z,.true.)
			! Fix fractions
			Tax_Increase_draft_group = 100.0_dp*Tax_Increase_draft_group/size_draft_group

		! Frac. capital tax rate increase by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		Tax_Rate_Increase_tk_draft_group = Draft_Table(Tax_Rate_Increase_tk_draft_group_z,DBN_z,.true.)
			! Fix fractions
			Tax_Rate_Increase_tk_draft_group = 100.0_dp*Tax_Rate_Increase_tk_draft_group/size_draft_group

		! Frac. labor tax rate increase by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		Tax_Rate_Increase_tl_draft_group = Draft_Table(Tax_Rate_Increase_tl_draft_group_z,DBN_z,.true.)
			! Fix fractions
			Tax_Rate_Increase_tl_draft_group = 100.0_dp*Tax_Rate_Increase_tl_draft_group/size_draft_group

		! Frac. total tax rate increase by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		Tax_Rate_Increase_draft_group = Draft_Table(Tax_Rate_Increase_draft_group_z,DBN_z,.true.)
			! Fix fractions
			Tax_Rate_Increase_draft_group = 100.0_dp*Tax_Rate_Increase_draft_group/size_draft_group

		! Frac. total income increase by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		Inc_Increase_draft_group   = Draft_Table(  Inc_Increase_draft_group_z,DBN_z,.true.)
		K_Inc_Increase_draft_group = Draft_Table(K_Inc_Increase_draft_group_z,DBN_z,.true.)
		L_Inc_Increase_draft_group = Draft_Table(L_Inc_Increase_draft_group_z,DBN_z,.true.)
			! Fix fractions
			  Inc_Increase_draft_group = 100.0_dp*  Inc_Increase_draft_group/size_draft_group
			K_Inc_Increase_draft_group = 100.0_dp*K_Inc_Increase_draft_group/size_draft_group
			L_Inc_Increase_draft_group = 100.0_dp*L_Inc_Increase_draft_group/size_draft_group

		OPEN (UNIT=80, FILE=trim(Result_Folder)//'draft_group_frac_Tax_K.txt', STATUS='replace') 
	    OPEN (UNIT=81, FILE=trim(Result_Folder)//'draft_group_frac_Tax_L.txt', STATUS='replace') 
	    OPEN (UNIT=82, FILE=trim(Result_Folder)//'draft_group_frac_Tax.txt', STATUS='replace') 
	    OPEN (UNIT=83, FILE=trim(Result_Folder)//'draft_group_frac_Tax_Rate_K.txt', STATUS='replace') 
	    OPEN (UNIT=84, FILE=trim(Result_Folder)//'draft_group_frac_Tax_Rate_L.txt', STATUS='replace') 
	    OPEN (UNIT=85, FILE=trim(Result_Folder)//'draft_group_frac_Tax_Rate.txt', STATUS='replace') 
	    OPEN (UNIT=86, FILE=trim(Result_Folder)//'draft_group_frac_Tot_Inc_Incr.txt', STATUS='replace') 
	    OPEN (UNIT=87, FILE=trim(Result_Folder)//'draft_group_frac_K_Inc_Incr.txt', STATUS='replace') 
	    OPEN (UNIT=88, FILE=trim(Result_Folder)//'draft_group_frac_L_Inc_Incr.txt', STATUS='replace') 
		do age = 1,draft_age_category
		    WRITE  (UNIT=80, FMT=*)  Tax_Increase_tk_draft_group(age,:)
		    WRITE  (UNIT=81, FMT=*)  Tax_Increase_tl_draft_group(age,:)
		    WRITE  (UNIT=82, FMT=*)  Tax_Increase_draft_group(age,:)
		    WRITE  (UNIT=83, FMT=*)  Tax_Rate_Increase_tk_draft_group(age,:)
		    WRITE  (UNIT=84, FMT=*)  Tax_Rate_Increase_tl_draft_group(age,:)
		    WRITE  (UNIT=85, FMT=*)  Tax_Rate_Increase_draft_group(age,:)
		    WRITE  (UNIT=86, FMT=*)  Inc_Increase_draft_group(age,:)
		    WRITE  (UNIT=87, FMT=*)  K_Inc_Increase_draft_group(age,:)
		    WRITE  (UNIT=88, FMT=*)  L_Inc_Increase_draft_group(age,:)
		ENDDO
		close(unit=80); close(unit=81); close(unit=82); close(unit=83); close(unit=84); close(unit=85)
		close(unit=86); close(unit=87); close(unit=88);
	endif 


	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	! Moments and Print Results
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	SSE_Moments = (1.0-Wealth_Output/3.0_DP)**2.0_DP  + (1.0_DP-prct1_wealth/0.34_DP)**2.0_DP  + (prct10_wealth-0.69_DP)**2.0_DP &
                   & + (1.0_DP-Std_Log_Earnings_25_60 / 0.8_DP)**2.0_DP + (1.0_DP-meanhours_25_60/0.4_DP)**2.0_DP &
                   & + (1.0_DP-MeanReturn/0.069_DP)**2.0_DP
	!print*,''
	!print*,"Current parameters"
	!print*,'beta',beta,'rho_z',rho_z,'sigma_z',sigma_z_eps,'sigma_lam',sigma_lambda_eps,'phi',phi
	print*,' '
	print*,'-----------------------------------------------------'
	print*,"Statistics"
	print*,' '
	print*,'	Debt/GDP',External_Debt_GDP,'A/GDP',Wealth_Output,'Top 1% A',prct1_wealth,'Top 10% A',prct10_wealth
	print*,'	STD Labor Earnings',Std_Log_Earnings_25_60,'Mean Labor (hours 25-60)',meanhours_25_60,'MeanReturn',MeanReturn
	print*,' '; print*,' Constrainted Firms and Demand for Capital'
	print*,' Z ','  Constrained_firms_by_z:     ',' Capital_high_shock ',' Capital_low_shock '
	do zi=1,nz
		print 12345, zi, & 
			100.0_dp*sum(constrained_firm_ind(:,:,zi,:,:,:)*DBN1(:,:,zi,:,:,:))/sum(DBN1(:,:,zi,:,:,:)), &
			100.0_dp*sum(constrained_firm_ind(:,:,zi,:,:,1)*DBN1(:,:,zi,:,:,1))/sum(DBN1(:,:,zi,:,:,1)), &
			100.0_dp*sum(constrained_firm_ind(:,:,zi,:,:,2)*DBN1(:,:,zi,:,:,2))/sum(DBN1(:,:,zi,:,:,2)), &
			(EBAR_data/(EBAR*0.727853584919652_dp))*(mu*P*xz_grid(1,zi)**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) , & 
			(EBAR_data/(EBAR*0.727853584919652_dp))*(mu*P*xz_grid(2,zi)**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu))  
	enddo 
	12345 format (I3,X,X,F7.2,X,X,F7.2,X,X,F7.2,X,X,E12.4,X,X,E12.4)
	print '(A,F7.3)', 'Total Constrained', 100.0_dp*sum(constrained_firm_ind*DBN1)
	print*, ' '
	print*,'Moments',SSE_Moments 
	print*,'-----------------------------------------------------'
	print*,''

	print*,' '
	print*,' Printing results into files'
	! Write in files some stats
	if (solving_bench.eq.1) then
		OPEN (UNIT=19, FILE=trim(Result_Folder)//'Asset_Stats_bench.txt', STATUS='replace')
		OPEN (UNIT=20, FILE=trim(Result_Folder)//'leverage_bench.txt', STATUS='replace')
	else
		OPEN (UNIT=19, FILE=trim(Result_Folder)//'Asset_Stats_exp.txt', STATUS='replace')
		OPEN (UNIT=20, FILE=trim(Result_Folder)//'leverage_exp.txt', STATUS='replace')
	end if 

		WRITE(UNIT=19, FMT=*) 'Stats on assets, return and savings'
		WRITE(UNIT=19, FMT=*) ' '
	! Aggregate Assets and Capital
		WRITE(UNIT=21, FMT=*) 'Assets and Capital'
		WRITE(UNIT=21, FMT=*) ' ', 'z1 ', 'z2 ', 'z3 ', 'z4 ', 'z5 ', 'z6 ', 'z7 ', 'Total'
		WRITE(UNIT=19, FMT=*) 'A', Wealth_by_z  , MeanWealth
		WRITE(UNIT=19, FMT=*) 'K', Capital_by_z , Mean_Capital 
		WRITE(UNIT=19, FMT=*) ' '
	! Return
		WRITE(UNIT=19, FMT=*) 'Return Stats'
		WRITE(UNIT=19, FMT=*) ' ', 'Assets', 'Capital'
		WRITE(UNIT=19, FMT=*) 'Mean_Return', MeanReturn, MeanReturn
		WRITE(UNIT=19, FMT=*) 'Std_Return', StdReturn, Std_K_Return
		do zi=1,nz
			write(rowname,*) zi 
			rowname = 'Mean_Return_z'//trim(rowname)
			WRITE(UNIT=19, FMT=*) trim(rowname), MeanReturn_by_z(zi), Mean_K_Return_by_z(zi)
		enddo 
		WRITE(UNIT=19, FMT=*) ' '
		WRITE(UNIT=19, FMT=*) 'Mean_Return_AT', MeanATReturn, MeanATReturn
		WRITE(UNIT=19, FMT=*) 'Std_Return_AT', StdATReturn, Std_AT_K_Return
		do zi=1,nz
			write(rowname,*) zi 
			rowname = 'Mean_Return_z'//trim(rowname)//'_AT'
			WRITE(UNIT=19, FMT=*) trim(rowname), MeanATReturn_by_z(zi), Mean_AT_K_Return_by_z(zi)
		enddo 
		WRITE(UNIT=19, FMT=*) ' '
	! Savings
		WRITE(UNIT=19, FMT=*) 'Saving Rates'
		WRITE(UNIT=19, FMT=*) '(Ap-A)/A'
		WRITE(UNIT=19, FMT=*) ' ', 'z1 ', 'z2 ', 'z3 ', 'z4 ', 'z5 ', 'z6 ', 'z7 ', 'Total'
		do group=1,max_age_category
			write(rowname,*) group  
			rowname = 'Age_Group_'//trim(rowname)
			WRITE(UNIT=19, FMT=*) rowname, S_Rate_A_AZ(group,:), S_Rate_A_Age(group)
		enddo 
		WRITE(UNIT=19, FMT=*) ' '
		WRITE(UNIT=19, FMT=*) '(Ap-A)/Y'
		WRITE(UNIT=19, FMT=*) ' ', 'z1 ', 'z2 ', 'z3 ', 'z4 ', 'z5 ', 'z6 ', 'z7 ', 'Total'
		do group=1,max_age_category
			write(rowname,*) group  
			rowname = 'Age_Group_'//trim(rowname)
			WRITE(UNIT=19, FMT=*) rowname, S_Rate_Y_AZ(group,:), S_Rate_Y_Age(group)
		enddo
		WRITE(UNIT=19, FMT=*) ' '
		WRITE(UNIT=19, FMT=*) ' ', '<90% ','90%<99% ','99%<'
		WRITE(UNIT=19, FMT=*) '(Ap-A)/A', S_Rate_A_W
		WRITE(UNIT=19, FMT=*) '(Ap-A)/Y', S_Rate_Y_W
		WRITE(UNIT=19, FMT=*) 'A', A_W
		WRITE(UNIT=19, FMT=*) 'Ap', Ap_W
		WRITE(UNIT=19, FMT=*) 'Y', Y_W
		WRITE(UNIT=19, FMT=*) ' '
		WRITE(UNIT=19, FMT=*) 'A'
		WRITE(UNIT=19, FMT=*) ' ', 'z1 ', 'z2 ', 'z3 ', 'z4 ', 'z5 ', 'z6 ', 'z7 ', 'Total'
		do group=1,max_age_category
			write(rowname,*) group  
			rowname = 'Age_Group_'//trim(rowname)
			WRITE(UNIT=19, FMT=*) rowname, A_AZ(group,:), A_Age(group)
		enddo 
		WRITE(UNIT=19, FMT=*) ' '
		WRITE(UNIT=19, FMT=*) 'Ap'
		WRITE(UNIT=19, FMT=*) ' ', 'z1 ', 'z2 ', 'z3 ', 'z4 ', 'z5 ', 'z6 ', 'z7 ', 'Total'
		do group=1,max_age_category
			write(rowname,*) group  
			rowname = 'Age_Group_'//trim(rowname)
			WRITE(UNIT=19, FMT=*) rowname, Ap_AZ(group,:), Ap_Age(group)
		enddo 
		WRITE(UNIT=19, FMT=*) ' '
		WRITE(UNIT=19, FMT=*) 'Y'
		WRITE(UNIT=19, FMT=*) ' ', 'z1 ', 'z2 ', 'z3 ', 'z4 ', 'z5 ', 'z6 ', 'z7 ', 'Total'
		do group=1,max_age_category
			write(rowname,*) group  
			rowname = 'Age_Group_'//trim(rowname)
			WRITE(UNIT=19, FMT=*) rowname, Y_AZ(group,:), Y_Age(group)
		enddo 
		WRITE(UNIT=19, FMT=*) ' '
	! Wealth percentiles
		WRITE(UNIT=19, FMT=*) ' '
		WRITE(UNIT=19, FMT=*) 'Percentiles of Asset Distribution'
		WRITE(UNIT=19, FMT=*) prctile_ai

	CLOSE(Unit=19)

	! Leverage and constrained firms 
		WRITE(UNIT=20, FMT=*) 'Leverage ','z1 ','z2 ','z3 ','z4 ','z5 ','z6 ','z7 ', ' ', 	&
							& 'Cons_Firm',' ',' ',											&	 
							& 'Cons_Firms ','z1 ','z2 ','z3 ','z4 ','z5 ','z6 ','z7 ', ' ', & 
							& 'Size_AZ ','z1 ','z2 ','z3 ','z4 ','z5 ','z6 ','z7 '
		do age=1,MaxAge 
			WRITE(UNIT=20, FMT=*) age, leverage_age_z(age,:), ' ',&
								& age, constrained_firms_age(age), ' ', age, constrained_firms_age_z(age,:), ' ', &
								& age, size_by_age_z(age,:)
		enddo 
	CLOSE(UNIT=20)

	! ! Save files of constrained index, output and profits
	! if (solving_bench.eq.1) then
	! 	OPEN(UNIT=1,  FILE=trim(Result_Folder)//'constrained_ind_bench'  , STATUS='replace')
	! 	OPEN(UNIT=2,  FILE=trim(Result_Folder)//'firm_output_bench'      , STATUS='replace')
	! 	OPEN(UNIT=3,  FILE=trim(Result_Folder)//'firm_profit_bench' 	    , STATUS='replace')
	! else 
	! 	OPEN(UNIT=1,  FILE=trim(Result_Folder)//'constrained_ind_exp'    , STATUS='replace')
	! 	OPEN(UNIT=2,  FILE=trim(Result_Folder)//'firm_output_exp'        , STATUS='replace')
	! 	OPEN(UNIT=3,  FILE=trim(Result_Folder)//'firm_profit_exp'  	    , STATUS='replace')
	! endif
	! 	WRITE(UNIT=1,FMT=*) constrained_firm_ind
	! 	WRITE(UNIT=2,FMT=*) Firm_Output
	! 	WRITE(UNIT=3,FMT=*) Firm_Profit
	! 	CLOSE(UNIT=1)
	! 	CLOSE(UNIT=2)
	! 	CLOSE(UNIT=3)


	print*, ' '; print*,' End of Compute_Stats'; print*, ' '


END SUBROUTINE COMPUTE_STATS

!========================================================================================
!========================================================================================
!========================================================================================

SUBROUTINE COMPUTE_WELFARE_GAIN()
	IMPLICIT NONE
	real(DP), dimension(MaxAge):: CumDiscountF
	! REAL(DP), dimension(MaxAge, na, nz, nlambda, ne) ::  Cons_Eq_Welfare ! ValueFunction_Bench, ValueFunction_Exp,
	REAL(DP), dimension(nz) ::  temp_ce_by_z
	REAL(DP), dimension(MaxAge, nz) :: frac_pos_welfare_by_age_z, size_pos_welfare_by_age_z, size_by_age_z
	INTEGER, dimension(max_age_category+1) :: age_limit
	INTEGER :: age_group_counter
	REAL(DP), dimension(max_age_category,nz) :: CE_by_agegroup_z, size_by_agegroup_z 
	REAL(DP), dimension(max_age_category,nz) :: size_pos_welfare_by_agegroup_z, frac_pos_welfare_by_agegroup_z  
	REAL(DP), dimension(draft_age_category,nz) :: CE_draft_group_z,  size_draft_group_z, frac_pos_welfare_draft_group_z
	REAL(DP), dimension(draft_age_category,draft_z_category) :: CE_draft_group,  size_draft_group, frac_pos_welfare_draft_group
	REAL(DP), dimension(draft_age_category,nz,nx) :: CE_draft_group_xz,  size_draft_group_xz, frac_pos_welfare_draft_group_xz
	REAL(DP), dimension(nz) 	:: DBN_Z, CDF_Z 
	REAL(DP), dimension(nz,nx) 	:: DBN_XZ
	INTEGER , dimension(draft_age_category+1) :: draft_age_limit
	INTEGER :: age2, z2
	REAL(DP):: K_Inc_aux, L_Inc_aux, cdf_xz, cdf_xz_low


	
	print*, ' Start of compute welfare gain'

	! Age Brackets
		age_limit       = [0, 5, 15, 25, 35, 45, 55, MaxAge ]
		draft_age_limit = [0, 1, 15, 30, 45, MaxAge ] 

	! Discount factor
		CumDiscountF(MaxAge)=1.0_DP
		DO age=MaxAge-1,1,-1
		    CumDiscountF(age)   = 1.0_DP + beta * survP(age) *CumDiscountF(age+1) 
		ENDDO
		!print*,CumDiscountF
		!PAUSE

	! Solve for the benchmark economy 
		solving_bench = 1
		tauK    = tauK_bench
		R       = R_bench
		P       = P_bench
		wage    = wage_bench
		Ebar    = Ebar_bench
		tauW_bt = tauW_bt_bench
		tauW_at = tauW_at_bench
		psi     = psi_bench
		tauPL   = tauPL_bench
		Y_a_threshold = Y_a_threshold_bench

		Cons   = Cons_bench
		Hours  = Hours_bench
		Aprime = Aprime_bench

		! print*,'BENCH: P=',P,'wage=',wage,'Ebar=',Ebar
		! CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
		! CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
		! CALL ComputeLaborUnits(Ebar, wage) 
		! CALL EGM_RETIREMENT_WORKING_PERIOD 

	! Compute the value function using interpolation and save it
		! CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction_Bench,Bq_Value_Bench)
		! CALL COMPUTE_VALUE_FUNCTION_SPLINE  
		! ValueFunction_Bench = ValueFunction


	! Profit Matrix
		Pr_mat = Profit_Matrix(R,P)

	! Get distribution by age and z
		DO age=1,MaxAge    
		    DO zi=1,nz
		        size_by_age_z(age, zi) = sum(DBN_bench(age,:,zi,:,:,:))
		    ENDDO ! zi
		ENDDO

		

	!=========================================== SOLVING EXP NEXT =================================
	! Solve the experimental economy  
		solving_bench = 0  
		tauK    = tauK_exp
		R       = R_exp
		P       = P_exp
		wage    = wage_exp
		Ebar    = Ebar_exp
		tauW_bt = tauW_bt_exp
		tauW_at = tauW_at_exp
		psi     = psi_exp
		tauPL   = tauPL_exp
		Y_a_threshold = Y_a_threshold_exp

		Cons   = Cons_exp
		Hours  = Hours_exp
		Aprime = Aprime_exp
		DBN1   = DBN_exp 

		!print*,' EXP: P=',P,'wage=',wage,'Ebar=',Ebar
		!CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
		!CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
		!CALL ComputeLaborUnits(Ebar, wage) 
		!CALL EGM_RETIREMENT_WORKING_PERIOD 

	! Compute the value function using interpolation and save it
		! CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction_Exp,Bq_Value_Exp)
		!CALL COMPUTE_VALUE_FUNCTION_SPLINE 
		! ValueFunction_Exp = ValueFunction


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Consumption Equivalent Welfare
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Compute Average Utility - CE2
		Av_Util_NB  =  100.0_dp*(( (sum(ValueFunction_exp(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:)) - &
	    					&  sum(Bq_Value_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))) / &
	                        &  sum((ValueFunction_Bench(1,:,:,:,:,:)-Bq_Value_bench(1,:,:,:,:,:))*DBN_bench(1,:,:,:,:,:))  ) &
	                        &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

		Av_Util_Pop = 100.0_dp*(( (sum(ValueFunction_exp*DBN1)-sum(Bq_Value_Bench*DBN_bench)) / &
	    					& sum((ValueFunction_Bench-Bq_Value_Bench)*DBN_bench)  ) &
	                        &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
	    print*,' '
	    print*,'---------------------------'
	    print*, ' CE 2 Computation'
	    print*, 'Av Utility NB (bench) =',sum(ValueFunction_Bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))  
	    print*, 'Av Utility NB (exp  ) =',sum(ValueFunction_exp(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))    
	    print*, 'CE2_NB =',Av_Util_NB
	    print*, 'CE2_Pop =', Av_Util_Pop
	    print*,'---------------------------'
	    print*,' '

	    ! Decomposition using distribution
		OPEN (UNIT=50, FILE=trim(Result_Folder)//'CE_alternative.txt', STATUS='replace') 
	    WRITE  (UNIT=50, FMT=*) 'CE2_NB =', 100.0_dp*(( (sum(ValueFunction_exp(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:)) - &
	    					&  sum(Bq_Value_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))) / &
	                        &  sum((ValueFunction_Bench(1,:,:,:,:,:)-Bq_Value_bench(1,:,:,:,:,:))*DBN_bench(1,:,:,:,:,:))  ) &
	                        &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
	    WRITE  (UNIT=50, FMT=*) 'CE2_Pop =', 100.0_dp*(( (sum(ValueFunction_exp*DBN1)-sum(Bq_Value_Bench*DBN_bench)) / &
	    					& sum((ValueFunction_Bench-Bq_Value_Bench)*DBN_bench)  ) &
	                        &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
	    WRITE  (UNIT=50, FMT=*) 'CE2_NB_bench =', &
	    					& 100.0_dp*(( sum((ValueFunction_exp(1,:,:,:,:,:)-Bq_Value_bench(1,:,:,:,:,:))*DBN_bench(1,:,:,:,:,:)) / &
	                        &  sum((ValueFunction_Bench(1,:,:,:,:,:)-Bq_Value_bench(1,:,:,:,:,:))*DBN_bench(1,:,:,:,:,:))  ) &
	                        &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
	    WRITE  (UNIT=50, FMT=*) 'CE2_Pop_bench =', 100.0_dp*(( sum((ValueFunction_exp-Bq_Value_Bench)*DBN_bench)/ &
	    					&  sum((ValueFunction_Bench-Bq_Value_Bench)*DBN_bench)) &
	                        &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
        WRITE  (UNIT=50, FMT=*) 'CE2_NB_exp =', &
        					& 100.0_dp*(( sum((ValueFunction_exp(1,:,:,:,:,:)-Bq_Value_bench(1,:,:,:,:,:))*DBN1(1,:,:,:,:,:)) / &
	                        &  sum((ValueFunction_Bench(1,:,:,:,:,:)-Bq_Value_bench(1,:,:,:,:,:))*DBN1(1,:,:,:,:,:))  ) &
	                        &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
	    WRITE  (UNIT=50, FMT=*) 'CE2_Pop_exp =', 100.0_dp*(( sum((ValueFunction_exp-Bq_Value_Bench)*DBN1) / &
	    					&  sum((ValueFunction_Bench-Bq_Value_Bench)*DBN1)  ) &
	                        &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
        WRITE  (UNIT=50, FMT=*) ' '

	    WRITE  (UNIT=50, FMT=*) 'Av Utility NB (bench) =',sum(ValueFunction_Bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))  
	    WRITE  (UNIT=50, FMT=*) 'Av Utility NB (exp  ) =',sum(ValueFunction_exp(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))    
	    WRITE  (UNIT=50, FMT=*) 'Av Utility  (bench)   =',sum(ValueFunction_Bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:))  
	    WRITE  (UNIT=50, FMT=*) 'Av Utility  (exp)     =',sum(ValueFunction_exp(:,:,:,:,:,:)*DBN1(:,:,:,:,:,:))
	    WRITE  (UNIT=50, FMT=*) 'Av BQ Utility NB (bench) =',sum(BQ_Value_Bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))  
	    WRITE  (UNIT=50, FMT=*) 'Av BQ Utility  (bench)   =',sum(BQ_Value_Bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:))  
	    close (unit=50)	




    ! Compute Individual Utility - CE1
		OPEN (UNIT=50, FILE=trim(Result_Folder)//'CE_NEWBORN.txt', STATUS='replace')  
		OPEN (UNIT=60, FILE=trim(Result_Folder)//'CE.txt', STATUS='replace')  
		OPEN (UNIT=70, FILE=trim(Result_Folder)//'CE_by_age.txt', STATUS='replace')  
		OPEN (UNIT=80, FILE=trim(Result_Folder)//'CE_by_age_z.txt', STATUS='replace')  

		DO age=1,MaxAge
			if (Log_Switch.eqv..true.) then 
		    	Cons_Eq_Welfare(age,:,:,:,:,:)= & 
		    		& exp((ValueFunction_exp(age,:,:,:,:,:)-ValueFunction_Bench(age,:,:,:,:,:))/CumDiscountF(age))-1.0_DP
		    else 
		    	Cons_Eq_Welfare(age,:,:,:,:,:)=((ValueFunction_exp(age,:,:,:,:,:)-Bq_Value_bench(age,:,:,:,:,:))/&
		    									& (ValueFunction_Bench(age,:,:,:,:,:)-Bq_Value_bench(age,:,:,:,:,:)) ) &
                                				&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP
		    end if 

		    WRITE  (UNIT=70, FMT=*) 100.0_dp*sum(Cons_Eq_Welfare(age,:,:,:,:,:)*DBN_bench(age,:,:,:,:,:))/sum(DBN_bench(age,:,:,:,:,:))
		    DO zi=1,nz
		         temp_ce_by_z(zi) = 100.0_dp*sum(Cons_Eq_Welfare(age,:,zi,:,:,:)*DBN_bench(age,:,zi,:,:,:))/sum(DBN_bench(age,:,zi,:,:,:))
		    ENDDO
		    WRITE  (UNIT=80, FMT=*) temp_ce_by_z
		    ! print*,'age=',age, temp_ce_by_z, ', mean:  ', &
		    !    & 100.0_dp*sum(Cons_Eq_Welfare(age,:,:,:,:,:)*DBN_bench(age,:,:,:,:,:))/sum(DBN_bench(age,:,:,:,:,:))
		ENDDO

		print*,'Test X for age=81 and z=8 and x=1'
		do age=1,na 
		print*,age, sum( ( ((ValueFunction_exp(81,age,8,3,3,1)-Bq_Value_bench(81,age,8,3,3,1))/&
		    									& (ValueFunction_Bench(81,age,8,3,3,1)-Bq_Value_bench(81,age,8,3,3,1)) ) &
                                				&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP )*DBN_bench(81,age,8,3,3,1) )& 
                                				!& /sum(DBN_bench(81,age,8,3,3,1)) ,& 
				& sum( ((ValueFunction_exp(81,age,8,3,3,1)-Bq_Value_bench(81,age,8,3,3,1)) )*DBN_bench(81,age,8,3,3,1) )& 
                                				!& /sum(DBN_bench(81,age,8,3,3,1)) , &
				& sum( ((ValueFunction_Bench(81,age,8,3,3,1)-Bq_Value_bench(81,age,8,3,3,1)) )*DBN_bench(81,age,8,3,3,1) )& 
                                				!& /sum(DBN_bench(81,age,8,3,3,1))
		enddo 

		CE_NEWBORN = 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))

		WRITE  (UNIT=50, FMT=*) 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
		WRITE  (UNIT=50, FMT=*) 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
		WRITE  (UNIT=60, FMT=*) 100.0_DP*sum(Cons_Eq_Welfare*DBN_bench)
		WRITE  (UNIT=60, FMT=*) 100.0_DP*sum(Cons_Eq_Welfare*DBN1)

		close (unit=50)
		close (unit=60)
		close (unit=70)
		close (unit=80)


		! CE by AGE-Z GROUP
		OPEN (UNIT=80, FILE=trim(Result_Folder)//'CE_by_AgeGroup_z.txt', STATUS='replace') 
		DO zi=1,nz
		    DO age_group_counter=1,max_age_category
		         CE_by_agegroup_z(age_group_counter,zi)= &
		            & 100.0_dp*sum(Cons_Eq_Welfare(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:,:)* &
		            &                         DBN_bench(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:,:))/&
		            &                sum( DBN_bench(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:,:))
		    ENDDO
		ENDDO
		DO age_group_counter=1,max_age_category
		    WRITE  (UNIT=80, FMT=*)  CE_by_agegroup_z(age_group_counter,:)
		ENDDO
		close (unit=80)


		! FRACTION POSITIVE WELFARE BY AGE-Z GROUP
		frac_pos_welfare 				= 0.0_DP
		size_pos_welfare_by_age_z 		= 0.0_DP
		size_by_agegroup_z 				= 0.0_DP
		size_pos_welfare_by_agegroup_z 	= 0.0_DP

		age_group_counter=1
		DO age=1,MaxAge 

		    DO while (age .gt. age_limit(age_group_counter+1) )
		        age_group_counter=age_group_counter+1
		    ENDDO    
		 
		 	DO xi=1,nx
		    DO ai=1,na
		    DO zi=1,nz
            DO lambdai=1,nlambda
            DO ei=1,ne
		                    
                If ( Cons_Eq_Welfare(age,ai,zi,lambdai,ei,xi) .ge. 0.0_DP) then
                    frac_pos_welfare = frac_pos_welfare + DBN_bench(age,ai,zi,lambdai,ei,xi)
                    size_pos_welfare_by_age_z(age,zi) = size_pos_welfare_by_age_z(age,zi) + DBN_bench(age,ai,zi,lambdai,ei,xi)
                    size_pos_welfare_by_agegroup_z(age_group_counter,zi) = size_pos_welfare_by_agegroup_z(age_group_counter,zi) &
                         &+  DBN_bench(age,ai,zi,lambdai,ei,xi)       
                endif 
                size_by_agegroup_z(age_group_counter,zi) = size_by_agegroup_z(age_group_counter,zi) + &
                         & DBN1(age,ai,zi,lambdai,ei,xi)       

	        ENDDO
	        ENDDO
	        ENDDO
		    ENDDO
		    ENDDO
		ENDDO

  
	OPEN (UNIT=60, FILE=trim(Result_Folder)//'frac_pos_welfare_by_agegroup_z.txt', STATUS='replace')  
	DO age_group_counter=1,max_age_category
	    WRITE (UNIT=60, FMT=*) size_pos_welfare_by_agegroup_z(age_group_counter,:)/ size_by_agegroup_z(age_group_counter,:)
	ENDDO
	CLOSE (UNIT=60);

	OPEN (UNIT=60, FILE=trim(Result_Folder)//'frac_pos_welfare.txt', STATUS='replace')  
	WRITE  (UNIT=60, FMT=*) 100.0_dp*frac_pos_welfare
	close (unit=60)

	OPEN (UNIT=60, FILE=trim(Result_Folder)//'frac_pos_welfare_by_age_z.txt', STATUS='replace')  
	DO age=1, MaxAge
	    WRITE  (UNIT=60, FMT=*) size_pos_welfare_by_age_z(age,:)/size_by_age_z(age,:)
	ENDDO
	close (UNIT=60)



	! Compute average welfare
		Welfare_Gain_Pop_bench = 100.0_DP*sum(Cons_Eq_Welfare*DBN_bench)
		Welfare_Gain_Pop_exp   = 100.0_DP*sum(Cons_Eq_Welfare*DBN1)
		Welfare_Gain_NB_bench  = 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
		Welfare_Gain_NB_exp    = 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))

	print*,''
	print*,'---------------------------'
	print*, ' CE 1 Computation'
	print*,' Pop_bench_DBN =',100.0_DP*sum(Cons_Eq_Welfare*DBN_bench),&
				& sum(ValueFunction_exp*DBN_bench),sum(ValueFunction_Bench*DBN_bench),sum(Bq_Value_bench*DBN_bench)
	print*,' Pop_exp_DBN   =',100.0_DP*sum(Cons_Eq_Welfare*DBN1),&
				& sum(ValueFunction_exp*DBN1),sum(ValueFunction_Bench*DBN1),sum(Bq_Value_bench*DBN1)
	print*,' NB_bench_DBN  =',&
	    & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
	print*,' NB_exp_dbn    =',&
	    & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
    print*,'Frac_Pos_Wel   =',100.0_dp*frac_pos_welfare
	print*,'---------------------------'
	print*,''


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Draft Tables
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	frac_pos_welfare_draft_group_z = 0.0_dp 
	size_draft_group_z             = 0.0_dp 
	CE_draft_group_z               = 0.0_dp 
	do zi  = 1,nz
	do age = 1,draft_age_category
		size_draft_group_z(age,zi) = &
			& sum(DBN_bench(draft_age_limit(age)+1:draft_age_limit(age+1),:,zi,:,:,:))

		CE_draft_group_z(age,zi) =  100.0_dp* &
			& sum(Cons_Eq_Welfare(draft_age_limit(age)+1:draft_age_limit(age+1),:,zi,:,:,:)* &
            &     DBN_bench(draft_age_limit(age)+1:draft_age_limit(age+1),:,zi,:,:,:))/&
            & sum(DBN_bench(draft_age_limit(age)+1:draft_age_limit(age+1),:,zi,:,:,:))

        do xi=1,nx
	    do ei=1,ne
	    do lambdai=1,nlambda
	    do ai=1,na
	    do age2=draft_age_limit(age)+1,draft_age_limit(age+1)
	    	If ( Cons_Eq_Welfare(age2,ai,zi,lambdai,ei,xi) .ge. 0.0_DP) then
        	frac_pos_welfare_draft_group_z(age,zi) = frac_pos_welfare_draft_group_z(age,zi) + DBN_bench(age2,ai,zi,lambdai,ei,xi)
        	endif 
    	enddo 
    	enddo 
    	enddo 
    	enddo 
    	enddo  
	enddo
	enddo 

	DBN_Z = sum(sum(sum(sum(sum(DBN_bench,6),5),4),2),1) 
	do zi=1,nz 
		CDF_Z(zi) = sum(DBN_Z(1:zi))
	enddo 
	! print*,' '
	! print*,'DBN_Z=',DBN_Z
	! print*,'CDF_Z=',CDF_Z 
	! print*,' '

	! Size of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
	size_draft_group = Draft_Table(size_draft_group_z,DBN_z,.true.)

	! CE of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
	CE_draft_group = Draft_Table(CE_draft_group_z,DBN_z,.false.)

	! Frac. pos. welfare by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
	frac_pos_welfare_draft_group = Draft_Table(frac_pos_welfare_draft_group_z,DBN_z,.true.)

	! Fix fractions
	frac_pos_welfare_draft_group = 100.0_dp*frac_pos_welfare_draft_group/size_draft_group
    
    OPEN (UNIT=80, FILE=trim(Result_Folder)//'draft_group_CE.txt', STATUS='replace') 
    OPEN (UNIT=82, FILE=trim(Result_Folder)//'draft_group_fpos_welfare.txt', STATUS='replace') 
    do age = 1,draft_age_category
	    WRITE  (UNIT=80, FMT=*)  CE_draft_group(age,:)
	    WRITE  (UNIT=82, FMT=*)  frac_pos_welfare_draft_group(age,:)
	ENDDO
	close(unit=80)
	close(unit=82)
	
	! print*,' '
	! print*,'CE_groups_z'
	! do age = 1,draft_age_category
	! 	print*, size_draft_group_z(age,:)
	! enddo 
	! print*,' '
	! print*,'CE_groups'
	! do age = 1,draft_age_category
	! 	print*, size_draft_group(age,:)
	! enddo 
	! print*,' '
	! print*,'Weights'
	! print*,(/ DBN_Z(1) , DBN_Z(2), DBN_Z(3), (0.40_dp-CDF_Z(3)) , DBN_Z(1)+DBN_Z(2)+DBN_Z(3)+(0.40_dp-CDF_Z(3)) /)/0.40_dp
	! print*,(/ (CDF_Z(4)-0.40_dp) , (0.80_dp-CDF_Z(4)) , (CDF_Z(4)-0.40_dp)+(0.80_dp-CDF_Z(4)) /)/0.40_dp
	! print*, 1.0_dp
	! print*,(/ (CDF_Z(5)-0.90_dp) , (0.99_dp-CDF_Z(5)) , (CDF_Z(5)-0.90_dp)+(0.99_dp-CDF_Z(5)) /)/0.09_dp
	! print*,(/ (CDF_Z(6)-0.99_dp) , (0.999_dp-CDF_Z(6)) , (CDF_Z(6)-0.99_dp)+(0.999_dp-CDF_Z(6)) /)/0.009_dp
	! print*,(/ (CDF_Z(7)-0.999_dp) ,  DBN_Z(8) , DBN_Z(9) , (CDF_Z(7)-0.999_dp)+DBN_Z(8)+DBN_Z(9) /)/0.001_dp
	! print*,(/ (CDF_Z(7)-0.999_dp) , (0.9999_dp-CDF_Z(7)) , (CDF_Z(7)-0.999_dp)+(0.9999_dp-CDF_Z(7)) /)/0.0009_dp
	! print*,(/ (CDF_Z(8)-0.9999_dp) , DBN_Z(9) , (CDF_Z(8)-0.9999_dp)+DBN_Z(9) /)/0.0001_dp


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Draft Tables by xz producitivity groups (cross-section of productivity, regardless of age)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	frac_pos_welfare_draft_group_xz = 0.0_dp 
	size_draft_group_xz = 0.0_dp
	DBN_XZ = sum(sum(sum(sum(DBN_bench,5),4),2),1) 

	do xi  = 1,nx
	do zi  = 1,nz
	do age = 1,draft_age_category
		size_draft_group_xz(age,zi,xi) = &
			& sum(DBN_bench(draft_age_limit(age)+1:draft_age_limit(age+1),:,zi,:,:,xi))

		CE_draft_group_xz(age,zi,xi) =  100.0_dp* &
			& sum(Cons_Eq_Welfare(draft_age_limit(age)+1:draft_age_limit(age+1),:,zi,:,:,xi)* &
            &     DBN_bench(draft_age_limit(age)+1:draft_age_limit(age+1),:,zi,:,:,xi))/&
            & sum(DBN_bench(draft_age_limit(age)+1:draft_age_limit(age+1),:,zi,:,:,xi))

	    do ei=1,ne
	    do lambdai=1,nlambda
	    do ai=1,na
	    do age2=draft_age_limit(age)+1,draft_age_limit(age+1)
	    	If ( Cons_Eq_Welfare(age2,ai,zi,lambdai,ei,xi) .ge. 0.0_DP) then
        	frac_pos_welfare_draft_group_xz(age,zi,xi) = frac_pos_welfare_draft_group_xz(age,zi,xi) + & 
        												&  DBN_bench(age2,ai,zi,lambdai,ei,xi)
        	endif 
    	enddo 
    	enddo 
    	enddo 
    	enddo 	
	enddo  
	enddo
	enddo 


	! CE of groups adjusting by x-z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
	CE_draft_group = Draft_Table_X(CE_draft_group_xz,DBN_xz,.false.)

	! Size of groups adjusting by x-z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
	size_draft_group = Draft_Table_X(size_draft_group_xz,DBN_xz,.true.)

	! Frac. pos. welfare by groups adjusting by x-z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
	frac_pos_welfare_draft_group = Draft_Table_X(frac_pos_welfare_draft_group_xz,DBN_xz,.true.)

	! Fix fractions
	frac_pos_welfare_draft_group = 100.0_dp*frac_pos_welfare_draft_group/size_draft_group

    OPEN (UNIT=80, FILE=trim(Result_Folder)//'draft_group_CE_xz.txt', STATUS='replace') 
    OPEN (UNIT=81, FILE=trim(Result_Folder)//'draft_group_size_xz.txt', STATUS='replace') 
    OPEN (UNIT=82, FILE=trim(Result_Folder)//'draft_group_fpos_welfare_xz.txt', STATUS='replace') 
	do age = 1,draft_age_category
	    WRITE  (UNIT=80, FMT=*)  CE_draft_group(age,:)
	    WRITE  (UNIT=81, FMT=*)  size_draft_group(age,:)
	    WRITE  (UNIT=82, FMT=*)  frac_pos_welfare_draft_group(age,:)
	ENDDO
	close(unit=80)
	close(unit=81)
	close(unit=82)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Draft Tables by xz producitivity groups by age-group
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Groups based on age group productivity
		! Age Group 1 
			! 00%-40%   : (z1,x1), (z2,x1), (z3,x1), (z4,x1)
			! 40%-80%   : (z4,x1), (z5,x1)
			! 80%-90%   : (z5,x1)
			! 90%-99%   : (z5,x1), (z6,x1)
			! 99%-99.9% : (z6,x1), (z7,x1)
			! 99.9%+    : (z7,x1), (z8,x1), (z9,x1)
			age = 1 
			DBN_XZ = sum(sum(sum(sum(DBN_bench(draft_age_limit(age)+1:draft_age_limit(age+1),:,:,:,:,:),5),4),2),1)
			DBN_XZ = DBN_XZ/sum(DBN_XZ) 

		! Welfare 
	    	cdf_xz = sum(DBN_XZ(1:3,1))
		CE_draft_group(age,1)   = ( DBN_XZ(1,1)*CE_draft_group_xz(age,1,1) + DBN_XZ(2,1)*CE_draft_group_xz(age,2,1) + &
								&   DBN_XZ(3,1)*CE_draft_group_xz(age,3,1) +   &
								&   (0.40_dp-cdf_xz)*CE_draft_group_xz(age,4,1) )/0.40_dp 
			cdf_xz = cdf_xz + DBN_XZ(4,1)
		CE_draft_group(age,2)   = ( (cdf_xz-0.40_dp)*CE_draft_group_xz(age,4,1) + (0.80_dp-cdf_xz)*CE_draft_group_xz(age,5,1) )/0.40_dp
			cdf_xz = cdf_xz + DBN_XZ(5,1)
		CE_draft_group(age,3)   = CE_draft_group_xz(age,5,1)
		CE_draft_group(age,4)   = ( (cdf_xz-0.80_dp)*CE_draft_group_xz(age,5,1) + (0.99_dp-cdf_xz)*CE_draft_group_xz(age,6,1) )/0.09_dp
			cdf_xz = cdf_xz + DBN_XZ(6,1)
		CE_draft_group(age,5)   = ( (cdf_xz-0.99_dp)*CE_draft_group_xz(age,6,1) + (0.999_dp-cdf_xz)*CE_draft_group_xz(age,7,1) )/0.009_dp
			cdf_xz = cdf_xz + DBN_XZ(7,1)
		CE_draft_group(age,6)   = ( (CDF_xz-0.999_dp)*CE_draft_group_xz(age,7,1) + &
								&  DBN_XZ(8,1)*CE_draft_group_xz(age,8,1) + DBN_XZ(9,1)*CE_draft_group_xz(age,9,1) )/0.001_dp
		! Size of Each group
	    	cdf_xz = sum(DBN_XZ(1:3,1))
		size_draft_group(age,1)   = size_draft_group_xz(age,1,1) + size_draft_group_xz(age,2,1) + size_draft_group_xz(age,3,1) +   &
								&   (0.40_dp-cdf_xz)*size_draft_group_xz(age,4,1)/DBN_XZ(4,1)
			cdf_xz = cdf_xz + DBN_XZ(4,1)
		size_draft_group(age,2)   = (cdf_xz-0.40_dp)*size_draft_group_xz(age,4,1)/DBN_XZ(4,1) + &
								&               (0.80_dp-cdf_xz)*size_draft_group_xz(age,5,1)/DBN_XZ(5,1)
			cdf_xz = cdf_xz + DBN_XZ(5,1)
		size_draft_group(age,3)   = 0.1_dp*size_draft_group_xz(age,5,1)/DBN_XZ(5,1)
		size_draft_group(age,4)   = (cdf_xz-0.80_dp)*size_draft_group_xz(age,5,1)/DBN_XZ(5,1) + &
								&               (0.99_dp-cdf_xz)*size_draft_group_xz(age,6,1)/DBN_XZ(6,1)
			cdf_xz = cdf_xz + DBN_XZ(6,1)
		size_draft_group(age,5)   = (cdf_xz-0.99_dp)*size_draft_group_xz(age,6,1)/DBN_XZ(6,1) + &
								&              (0.999_dp-cdf_xz)*size_draft_group_xz(age,7,1)/DBN_XZ(7,1)
			cdf_xz = cdf_xz + DBN_XZ(7,1)
		size_draft_group(age,6)   = (CDF_xz-0.999_dp)*size_draft_group_xz(age,7,1)/DBN_XZ(7,1) + &
								&  size_draft_group_xz(age,8,1) + size_draft_group_xz(age,9,1)


		! Fraction Positive Welfare
	    	cdf_xz = sum(DBN_XZ(1:3,1))
		frac_pos_welfare_draft_group(age,1)   = frac_pos_welfare_draft_group_xz(age,1,1) + &
								&   frac_pos_welfare_draft_group_xz(age,2,1) + &
								&   frac_pos_welfare_draft_group_xz(age,3,1) +   &
								&   (0.40_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,4,1)/DBN_XZ(4,1)
			cdf_xz = cdf_xz + DBN_XZ(4,1)
		frac_pos_welfare_draft_group(age,2)   = (cdf_xz-0.40_dp)*frac_pos_welfare_draft_group_xz(age,4,1)/DBN_XZ(4,1) + &
								&               (0.80_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,5,1)/DBN_XZ(5,1)
			cdf_xz = cdf_xz + DBN_XZ(5,1)
		frac_pos_welfare_draft_group(age,3)   = 0.1_dp*frac_pos_welfare_draft_group_xz(age,5,1)/DBN_XZ(5,1)
		frac_pos_welfare_draft_group(age,4)   = (cdf_xz-0.80_dp)*frac_pos_welfare_draft_group_xz(age,5,1)/DBN_XZ(5,1) + &
								&               (0.99_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,6,1)/DBN_XZ(6,1)
			cdf_xz = cdf_xz + DBN_XZ(6,1)
		frac_pos_welfare_draft_group(age,5)   = (cdf_xz-0.99_dp)*frac_pos_welfare_draft_group_xz(age,6,1)/DBN_XZ(6,1) + &
								&              (0.999_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,7,1)/DBN_XZ(7,1)
			cdf_xz = cdf_xz + DBN_XZ(7,1)
		frac_pos_welfare_draft_group(age,6)   = (CDF_xz-0.999_dp)*frac_pos_welfare_draft_group_xz(age,7,1)/DBN_XZ(7,1) + &
								&  frac_pos_welfare_draft_group_xz(age,8,1) + frac_pos_welfare_draft_group_xz(age,9,1)
			
		! Age Group 2
			! 00%-40%   : (z1,x3), (z2,x3), (z3,x3), (z4,x3), (z5,x3), (z6,x3), (z7,x3), (z8,x3), (z9,x3), (z1,x1), (z2,x1), (z3,x1)
			! 40%-80%   : (z3,x1), (z4,x1), (z5,x2)
			! 80%-90%   : (z5,x2), (z6,x2), (z7,x2), (z8,x2), (z9,x2), (z5,x1)
			! 90%-99%   : (z5,x1), (z6,x1)
			! 99%-99.9% : (z6,x1), (z7,x1)
			! 99.9%+    : (z7,x1), (z8,x1), (z9,x1)
			age = 2
			DBN_XZ = sum(sum(sum(sum(DBN_bench(draft_age_limit(age)+1:draft_age_limit(age+1),:,:,:,:,:),5),4),2),1) 
			DBN_XZ = DBN_XZ/sum(DBN_XZ)

		! Welfare
	    	cdf_xz = sum(DBN_XZ(:,3)) + sum(DBN_XZ(1:2,1))
		CE_draft_group(age,1)   = ( DBN_XZ(1,3)*CE_draft_group_xz(age,1,3) + DBN_XZ(2,3)*CE_draft_group_xz(age,2,3) + &
								& DBN_XZ(3,3)*CE_draft_group_xz(age,3,3) + DBN_XZ(4,3)*CE_draft_group_xz(age,4,3) + &
								& DBN_XZ(5,3)*CE_draft_group_xz(age,5,3) + DBN_XZ(6,3)*CE_draft_group_xz(age,6,3) + &
								& DBN_XZ(7,3)*CE_draft_group_xz(age,7,3) + DBN_XZ(8,3)*CE_draft_group_xz(age,8,3) + &
								& DBN_XZ(9,3)*CE_draft_group_xz(age,9,3) + &
								& DBN_XZ(1,1)*CE_draft_group_xz(age,1,1) + DBN_XZ(2,1)*CE_draft_group_xz(age,2,1) + &
								& (0.40_dp-cdf_xz)*CE_draft_group_xz(age,3,1) )/0.40_dp 
			cdf_xz = cdf_xz + DBN_XZ(3,1)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + DBN_XZ(4,1)
		CE_draft_group(age,2)   = ( (cdf_xz_low-0.40_dp)*CE_draft_group_xz(age,3,1) + &
								& DBN_XZ(4,1)*CE_draft_group_xz(age,4,1) + &
								& (0.80_dp-cdf_xz)*CE_draft_group_xz(age,5,2) )/0.40_dp
			cdf_xz = cdf_xz + DBN_XZ(5,2)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(6:,2))
		CE_draft_group(age,3)   = ( (cdf_xz-0.80_dp)*CE_draft_group_xz(age,5,2) + &
								& DBN_XZ(6,2)*CE_draft_group_xz(age,6,2) + DBN_XZ(7,2)*CE_draft_group_xz(age,7,2) + &
								& DBN_XZ(8,2)*CE_draft_group_xz(age,8,2) + DBN_XZ(9,2)*CE_draft_group_xz(age,9,2) + &
								& (0.90_dp-cdf_xz)*CE_draft_group_xz(age,5,1) )/0.10_dp
			cdf_xz = cdf_xz + DBN_XZ(5,1)
		CE_draft_group(age,4)   = ( (cdf_xz-0.90_dp)*CE_draft_group_xz(age,5,1) + (0.99_dp-cdf_xz)*CE_draft_group_xz(age,6,1) )/0.09_dp
			cdf_xz = cdf_xz + DBN_XZ(6,1)
		CE_draft_group(age,5)   = ( (cdf_xz-0.99_dp)*CE_draft_group_xz(age,6,1) + (0.999_dp-cdf_xz)*CE_draft_group_xz(age,7,1) )/0.009_dp
			cdf_xz = cdf_xz + DBN_XZ(7,1)
		CE_draft_group(age,6)   = ( (CDF_xz-0.999_dp)*CE_draft_group_xz(age,7,1) + &
								&  DBN_XZ(8,1)*CE_draft_group_xz(age,8,1) + DBN_XZ(9,1)*CE_draft_group_xz(age,9,1) )/0.001_dp

		! Size of Each group
	    	cdf_xz = sum(DBN_XZ(:,3)) + sum(DBN_XZ(1:2,1))
		size_draft_group(age,1)   = size_draft_group_xz(age,1,3) + size_draft_group_xz(age,2,3) + &
								& size_draft_group_xz(age,3,3) + size_draft_group_xz(age,4,3) + &
								& size_draft_group_xz(age,5,3) + size_draft_group_xz(age,6,3) + &
								& size_draft_group_xz(age,7,3) + size_draft_group_xz(age,8,3) + &
								& size_draft_group_xz(age,9,3) + &
								& size_draft_group_xz(age,1,1) + size_draft_group_xz(age,2,1) + &
								& (0.40_dp-cdf_xz)*size_draft_group_xz(age,3,1)/DBN_XZ(3,1)
			cdf_xz = cdf_xz + DBN_XZ(3,1)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + DBN_XZ(4,1)
		size_draft_group(age,2)   = (cdf_xz_low-0.40_dp)*size_draft_group_xz(age,3,1)/DBN_XZ(3,1) + &
								& size_draft_group_xz(age,4,1) + &
								& (0.80_dp-cdf_xz)*size_draft_group_xz(age,5,2)/DBN_XZ(5,2)
			cdf_xz = cdf_xz + DBN_XZ(5,2)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(6:,2))
		size_draft_group(age,3)   = (cdf_xz-0.80_dp)*size_draft_group_xz(age,5,2)/DBN_XZ(5,2) + &
								& size_draft_group_xz(age,6,2) + size_draft_group_xz(age,7,2) + &
								& size_draft_group_xz(age,8,2) + size_draft_group_xz(age,9,2) + &
								& (0.90_dp-cdf_xz)*size_draft_group_xz(age,5,1)/DBN_XZ(5,1)
			cdf_xz = cdf_xz + DBN_XZ(5,1)
		size_draft_group(age,4)   = (cdf_xz-0.90_dp)*size_draft_group_xz(age,5,1)/DBN_XZ(5,1) + &
								&   (0.99_dp-cdf_xz)*size_draft_group_xz(age,6,1)/DBN_XZ(6,1)
			cdf_xz = cdf_xz + DBN_XZ(6,1)
		size_draft_group(age,5)   = (cdf_xz-0.99_dp)*size_draft_group_xz(age,6,1)/DBN_XZ(6,1) + &
								&   (0.999_dp-cdf_xz)*size_draft_group_xz(age,7,1)/DBN_XZ(7,1)
			cdf_xz = cdf_xz + DBN_XZ(7,1)
		size_draft_group(age,6)   = (CDF_xz-0.999_dp)*size_draft_group_xz(age,7,1)/DBN_XZ(7,1) + &
								&  size_draft_group_xz(age,8,1) + size_draft_group_xz(age,9,1)

		! Fraction Positive Welfare
	    	cdf_xz = sum(DBN_XZ(:,3)) + sum(DBN_XZ(1:2,1))
		frac_pos_welfare_draft_group(age,1)   = frac_pos_welfare_draft_group_xz(age,1,3) + frac_pos_welfare_draft_group_xz(age,2,3) + &
								& frac_pos_welfare_draft_group_xz(age,3,3) + frac_pos_welfare_draft_group_xz(age,4,3) + &
								& frac_pos_welfare_draft_group_xz(age,5,3) + frac_pos_welfare_draft_group_xz(age,6,3) + &
								& frac_pos_welfare_draft_group_xz(age,7,3) + frac_pos_welfare_draft_group_xz(age,8,3) + &
								& frac_pos_welfare_draft_group_xz(age,9,3) + &
								& frac_pos_welfare_draft_group_xz(age,1,1) + frac_pos_welfare_draft_group_xz(age,2,1) + &
								& (0.40_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,3,1)/DBN_XZ(3,1)
			cdf_xz = cdf_xz + DBN_XZ(3,1)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + DBN_XZ(4,1)
		frac_pos_welfare_draft_group(age,2)   = (cdf_xz_low-0.40_dp)*frac_pos_welfare_draft_group_xz(age,3,1)/DBN_XZ(3,1) + &
								& frac_pos_welfare_draft_group_xz(age,4,1) + &
								& (0.80_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,5,2)/DBN_XZ(5,2)
			cdf_xz = cdf_xz + DBN_XZ(5,2)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(6:,2))
		frac_pos_welfare_draft_group(age,3)   = (cdf_xz-0.80_dp)*frac_pos_welfare_draft_group_xz(age,5,2)/DBN_XZ(5,2) + &
								& frac_pos_welfare_draft_group_xz(age,6,2) + frac_pos_welfare_draft_group_xz(age,7,2) + &
								& frac_pos_welfare_draft_group_xz(age,8,2) + frac_pos_welfare_draft_group_xz(age,9,2) + &
								& (0.90_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,5,1)/DBN_XZ(5,1)
			cdf_xz = cdf_xz + DBN_XZ(5,1)
		frac_pos_welfare_draft_group(age,4)   = (cdf_xz-0.90_dp)*frac_pos_welfare_draft_group_xz(age,5,1)/DBN_XZ(5,1) + &
								&               (0.99_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,6,1)/DBN_XZ(6,1)
			cdf_xz = cdf_xz + DBN_XZ(6,1)
		frac_pos_welfare_draft_group(age,5)   = (cdf_xz-0.99_dp)*frac_pos_welfare_draft_group_xz(age,6,1)/DBN_XZ(6,1) + &
								&				(0.999_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,7,1)/DBN_XZ(7,1)
			cdf_xz = cdf_xz + DBN_XZ(7,1)
		frac_pos_welfare_draft_group(age,6)   = (CDF_xz-0.999_dp)*frac_pos_welfare_draft_group_xz(age,7,1)/DBN_XZ(7,1) + &
								&  frac_pos_welfare_draft_group_xz(age,8,1) + frac_pos_welfare_draft_group_xz(age,9,1)

		! Age Group 3
			! 00%-40%   : (z1,x3), (z2,x3), (z3,x3), (z4,x3), (z5,x3)
			! 40%-80%   : (z5,x3), (z6,x3), (z7,x3), (z8,x3), (z9,x3), (z1,x1), (z2,x1), (z3,x1), (z4,x1)
			! 80%-90%   : (z4,x1), (z5,x2)
			! 90%-99%   : (z5,x2), (z6,x2), (z7,x2), (z8,x2), (z9,x2), (z5,x1), (z6,x1)
			! 99%-99.9% : (z6,x1), (z7,x1)
			! 99.9%+    : (z7,x1), (z8,x1), (z9,x1)
			age = 3
			DBN_XZ = sum(sum(sum(sum(DBN_bench(draft_age_limit(age)+1:draft_age_limit(age+1),:,:,:,:,:),5),4),2),1) 
			DBN_XZ = DBN_XZ/sum(DBN_XZ)

		! Welfare
	    	cdf_xz = sum(DBN_XZ(1:4,3))
		CE_draft_group(age,1)   = ( DBN_XZ(1,3)*CE_draft_group_xz(age,1,3) + DBN_XZ(2,3)*CE_draft_group_xz(age,2,3) + &
								& DBN_XZ(3,3)*CE_draft_group_xz(age,3,3) + DBN_XZ(4,3)*CE_draft_group_xz(age,4,3) + &
								& (0.40_dp-cdf_xz)*CE_draft_group_xz(age,5,3) )/0.40_dp 
			cdf_xz = cdf_xz + DBN_XZ(5,3)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(6:,3)) + sum(DBN_XZ(1:3,1))
		CE_draft_group(age,2)   = ( (cdf_xz_low-0.40_dp)*CE_draft_group_xz(age,5,3) + &
								& DBN_XZ(6,3)*CE_draft_group_xz(age,6,3) + &
								& DBN_XZ(7,3)*CE_draft_group_xz(age,7,3) + DBN_XZ(8,3)*CE_draft_group_xz(age,8,3) + &
								& DBN_XZ(9,3)*CE_draft_group_xz(age,9,3) + &
								& DBN_XZ(1,1)*CE_draft_group_xz(age,1,1) + DBN_XZ(2,1)*CE_draft_group_xz(age,2,1) + &
								& DBN_XZ(4,1)*CE_draft_group_xz(age,3,1) + &
								& (0.80_dp-cdf_xz)*CE_draft_group_xz(age,4,1) )/0.40_dp
			cdf_xz = cdf_xz + DBN_XZ(4,1)
		CE_draft_group(age,3)   = ( (cdf_xz-0.80_dp)*CE_draft_group_xz(age,4,1) + (0.90_dp-cdf_xz)*CE_draft_group_xz(age,5,2) )/0.10_dp
			cdf_xz = cdf_xz + DBN_XZ(5,2)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(6:,2)) + DBN_XZ(5,1)
		CE_draft_group(age,4)   = ( (cdf_xz-0.90_dp)*CE_draft_group_xz(age,5,2) + &
								& DBN_XZ(6,2)*CE_draft_group_xz(age,6,2) + DBN_XZ(7,2)*CE_draft_group_xz(age,7,2) + &
								& DBN_XZ(8,2)*CE_draft_group_xz(age,8,2) + DBN_XZ(9,2)*CE_draft_group_xz(age,9,2) + &
								& DBN_XZ(5,1)*CE_draft_group_xz(age,5,1) + &
								& (0.99_dp-cdf_xz)*CE_draft_group_xz(age,6,1) )/0.09_dp
			cdf_xz = cdf_xz + DBN_XZ(6,1)
		CE_draft_group(age,5)   = ( (cdf_xz-0.99_dp)*CE_draft_group_xz(age,6,1) + (0.999_dp-cdf_xz)*CE_draft_group_xz(age,7,1) )/0.009_dp
			cdf_xz = cdf_xz + DBN_XZ(7,1)
		CE_draft_group(age,6)   = ( (CDF_xz-0.999_dp)*CE_draft_group_xz(age,7,1) + &
								&  DBN_XZ(8,1)*CE_draft_group_xz(age,8,1) + DBN_XZ(9,1)*CE_draft_group_xz(age,9,1) )/0.001_dp

		! Size of Each group
	    	cdf_xz = sum(DBN_XZ(1:4,3))
		size_draft_group(age,1)   = size_draft_group_xz(age,1,3) + size_draft_group_xz(age,2,3) + &
								& size_draft_group_xz(age,3,3) + size_draft_group_xz(age,4,3) + &
								& (0.40_dp-cdf_xz)*size_draft_group_xz(age,5,3)/DBN_XZ(5,3)
			cdf_xz = cdf_xz + DBN_XZ(5,3)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(6:,3)) + sum(DBN_XZ(1:3,1))
		size_draft_group(age,2)   = (cdf_xz_low-0.40_dp)*size_draft_group_xz(age,5,3)/DBN_XZ(5,3) + &
								& size_draft_group_xz(age,6,3) + &
								& size_draft_group_xz(age,7,3) + size_draft_group_xz(age,8,3) + &
								& size_draft_group_xz(age,9,3) + &
								& size_draft_group_xz(age,1,1) + size_draft_group_xz(age,2,1) + &
								& size_draft_group_xz(age,3,1) + &
								& (0.80_dp-cdf_xz)*size_draft_group_xz(age,4,1)/DBN_XZ(4,1)
			cdf_xz = cdf_xz + DBN_XZ(4,1)
		size_draft_group(age,3)   = (cdf_xz-0.80_dp)*size_draft_group_xz(age,4,1)/DBN_XZ(4,1) + &
								&   (0.90_dp-cdf_xz)*size_draft_group_xz(age,5,2)/DBN_XZ(5,2)
			cdf_xz = cdf_xz + DBN_XZ(5,2)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(6:,2)) + DBN_XZ(5,1)
		size_draft_group(age,4)   = (cdf_xz-0.90_dp)*size_draft_group_xz(age,5,2)/DBN_XZ(5,2) + &
								& size_draft_group_xz(age,6,2) + size_draft_group_xz(age,7,2) + &
								& size_draft_group_xz(age,8,2) + size_draft_group_xz(age,9,2) + &
								& size_draft_group_xz(age,5,1) + &
								& (0.99_dp-cdf_xz)*size_draft_group_xz(age,6,1)/DBN_XZ(6,1)
			cdf_xz = cdf_xz + DBN_XZ(6,1)
		size_draft_group(age,5)   = (cdf_xz-0.99_dp)*size_draft_group_xz(age,6,1)/DBN_XZ(6,1) + &
								& (0.999_dp-cdf_xz)*size_draft_group_xz(age,7,1)/DBN_XZ(7,1)
			cdf_xz = cdf_xz + DBN_XZ(7,1)
		size_draft_group(age,6)   = (CDF_xz-0.999_dp)*size_draft_group_xz(age,7,1)/DBN_XZ(7,1) + &
								&  size_draft_group_xz(age,8,1) + size_draft_group_xz(age,9,1)


		! Fraction Positive Welfare
	    	cdf_xz = sum(DBN_XZ(1:4,3))
		frac_pos_welfare_draft_group(age,1)   = frac_pos_welfare_draft_group_xz(age,1,3) + frac_pos_welfare_draft_group_xz(age,2,3) + &
								& frac_pos_welfare_draft_group_xz(age,3,3) + frac_pos_welfare_draft_group_xz(age,4,3) + &
								& (0.40_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,5,3)/DBN_XZ(5,3)
			cdf_xz = cdf_xz + DBN_XZ(5,3)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(6:,3)) + sum(DBN_XZ(1:3,1))
		frac_pos_welfare_draft_group(age,2)   = (cdf_xz_low-0.40_dp)*frac_pos_welfare_draft_group_xz(age,5,3)/DBN_XZ(5,3) + &
								& frac_pos_welfare_draft_group_xz(age,6,3) + &
								& frac_pos_welfare_draft_group_xz(age,7,3) + frac_pos_welfare_draft_group_xz(age,8,3) + &
								& frac_pos_welfare_draft_group_xz(age,9,3) + &
								& frac_pos_welfare_draft_group_xz(age,1,1) + frac_pos_welfare_draft_group_xz(age,2,1) + &
								& frac_pos_welfare_draft_group_xz(age,3,1) + &
								& (0.80_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,4,1)/DBN_XZ(4,1)
			cdf_xz = cdf_xz + DBN_XZ(4,1)
		frac_pos_welfare_draft_group(age,3)   = (cdf_xz-0.80_dp)*frac_pos_welfare_draft_group_xz(age,4,1)/DBN_XZ(4,1) + &
								& (0.90_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,5,2)/DBN_XZ(5,2)
			cdf_xz = cdf_xz + DBN_XZ(5,2)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(6:,2)) + DBN_XZ(5,1)
		frac_pos_welfare_draft_group(age,4)   = (cdf_xz-0.90_dp)*frac_pos_welfare_draft_group_xz(age,5,2)/DBN_XZ(5,2) + &
								& frac_pos_welfare_draft_group_xz(age,6,2) + frac_pos_welfare_draft_group_xz(age,7,2) + &
								& frac_pos_welfare_draft_group_xz(age,8,2) + frac_pos_welfare_draft_group_xz(age,9,2) + &
								& frac_pos_welfare_draft_group_xz(age,5,1) + &
								& (0.99_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,6,1)/DBN_XZ(6,1)
			cdf_xz = cdf_xz + DBN_XZ(6,1)
		frac_pos_welfare_draft_group(age,5)   = (cdf_xz-0.99_dp)*frac_pos_welfare_draft_group_xz(age,6,1)/DBN_XZ(6,1) + &
								& (0.999_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,7,1)/DBN_XZ(7,1)
			cdf_xz = cdf_xz + DBN_XZ(7,1)
		frac_pos_welfare_draft_group(age,6)   = (CDF_xz-0.999_dp)*frac_pos_welfare_draft_group_xz(age,7,1)/DBN_XZ(7,1) + &
								&  frac_pos_welfare_draft_group_xz(age,8,1) + frac_pos_welfare_draft_group_xz(age,9,1)

		! Age Group 4
			! 00%-40%   : (z1,x3), (z2,x3), (z3,x3), (z4,x3)
			! 40%-80%   : (z4,x3), (z5,x3), (z6,x3), (z7,x3), (z8,x3), (z9,x3), (z1,x1), (z2,x1), (z3,x1), (z4,x1)
			! 80%-90%   : (z4,x1), (z5,x2)
			! 90%-99%   : (z5,x2), (z6,x2), (z7,x2), (z8,x2), (z9,x2), (z5,x1)
			! 99%-99.9% : (z5,x1), (z6,x1)
			! 99.9%+    : (z6,x1), (z7,x1), (z8,x1), (z9,x1)
			age = 4
			DBN_XZ = sum(sum(sum(sum(DBN_bench(draft_age_limit(age)+1:draft_age_limit(age+1),:,:,:,:,:),5),4),2),1)
			DBN_XZ = DBN_XZ/sum(DBN_XZ)

		! Welfare
	    	cdf_xz = sum(DBN_XZ(1:3,3))
		CE_draft_group(age,1)   = ( DBN_XZ(1,3)*CE_draft_group_xz(age,1,3) + DBN_XZ(2,3)*CE_draft_group_xz(age,2,3) + &
								& DBN_XZ(3,3)*CE_draft_group_xz(age,3,3) + &
								& (0.40_dp-cdf_xz)*CE_draft_group_xz(age,4,3) )/0.40_dp 
			cdf_xz = cdf_xz + DBN_XZ(4,3)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(5:,3)) + sum(DBN_XZ(1:3,1))
		CE_draft_group(age,2)   = ( (cdf_xz_low-0.40_dp)*CE_draft_group_xz(age,4,3) + &
								& DBN_XZ(5,3)*CE_draft_group_xz(age,5,3) + DBN_XZ(6,3)*CE_draft_group_xz(age,6,3) + &
								& DBN_XZ(7,3)*CE_draft_group_xz(age,7,3) + DBN_XZ(8,3)*CE_draft_group_xz(age,8,3) + &
								& DBN_XZ(9,3)*CE_draft_group_xz(age,9,3) + &
								& DBN_XZ(1,1)*CE_draft_group_xz(age,1,1) + DBN_XZ(2,1)*CE_draft_group_xz(age,2,1) + &
								& DBN_XZ(3,1)*CE_draft_group_xz(age,3,1) + &
								& (0.80_dp-cdf_xz)*CE_draft_group_xz(age,4,1) )/0.40_dp
			cdf_xz = cdf_xz + DBN_XZ(4,1)
		CE_draft_group(age,3)   = ( (cdf_xz-0.80_dp)*CE_draft_group_xz(age,4,1) + (0.90_dp-cdf_xz)*CE_draft_group_xz(age,5,2) )/0.10_dp
			cdf_xz = cdf_xz + DBN_XZ(5,2)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(6:,2)) 
		CE_draft_group(age,4)   = ( (cdf_xz-0.90_dp)*CE_draft_group_xz(age,5,2) + &
								& DBN_XZ(6,2)*CE_draft_group_xz(age,6,2) + DBN_XZ(7,2)*CE_draft_group_xz(age,7,2) + &
								& DBN_XZ(8,2)*CE_draft_group_xz(age,8,2) + DBN_XZ(9,2)*CE_draft_group_xz(age,9,2) + &
								& (0.99_dp-cdf_xz)*CE_draft_group_xz(age,5,1) )/0.09_dp
			cdf_xz = cdf_xz + DBN_XZ(5,1)
		CE_draft_group(age,5)   = ( (cdf_xz-0.99_dp)*CE_draft_group_xz(age,5,1) + (0.999_dp-cdf_xz)*CE_draft_group_xz(age,6,1) )/0.009_dp
			cdf_xz = cdf_xz + DBN_XZ(6,1)
		CE_draft_group(age,6)   = ( (CDF_xz-0.999_dp)*CE_draft_group_xz(age,6,1) + DBN_XZ(7,1)*CE_draft_group_xz(age,7,1) + &
								&  DBN_XZ(8,1)*CE_draft_group_xz(age,8,1) + DBN_XZ(9,1)*CE_draft_group_xz(age,9,1) )/0.001_dp

		! Size of Each group
	    	cdf_xz = sum(DBN_XZ(1:3,3))
		size_draft_group(age,1)   = size_draft_group_xz(age,1,3) + size_draft_group_xz(age,2,3) + &
								& size_draft_group_xz(age,3,3) + &
								& (0.40_dp-cdf_xz)*size_draft_group_xz(age,4,3)/DBN_XZ(4,3)
			cdf_xz = cdf_xz + DBN_XZ(4,3)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(5:,3)) + sum(DBN_XZ(1:3,1))
		size_draft_group(age,2)   = (cdf_xz_low-0.40_dp)*size_draft_group_xz(age,4,3)/DBN_XZ(4,3) + &
								& size_draft_group_xz(age,5,3) + size_draft_group_xz(age,6,3) + &
								& size_draft_group_xz(age,7,3) + size_draft_group_xz(age,8,3) + &
								& size_draft_group_xz(age,9,3) + &
								& size_draft_group_xz(age,1,1) + size_draft_group_xz(age,2,1) + &
								& size_draft_group_xz(age,3,1) + &
								& (0.80_dp-cdf_xz)*size_draft_group_xz(age,4,1)/DBN_XZ(4,1)
			cdf_xz = cdf_xz + DBN_XZ(4,1)
		size_draft_group(age,3)   = (cdf_xz-0.80_dp)*size_draft_group_xz(age,4,1)/DBN_XZ(4,1) + &
								& (0.90_dp-cdf_xz)*size_draft_group_xz(age,5,2)/DBN_XZ(5,2)
			cdf_xz = cdf_xz + DBN_XZ(5,2)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(6:,2)) 
		size_draft_group(age,4)   = (cdf_xz-0.90_dp)*size_draft_group_xz(age,5,2)/DBN_XZ(5,2) + &
								& size_draft_group_xz(age,6,2) + size_draft_group_xz(age,7,2) + &
								& size_draft_group_xz(age,8,2) + size_draft_group_xz(age,9,2) + &
								& (0.99_dp-cdf_xz)*size_draft_group_xz(age,5,1)/DBN_XZ(5,1)
			cdf_xz = cdf_xz + DBN_XZ(5,1)
		size_draft_group(age,5)   = (cdf_xz-0.99_dp)*size_draft_group_xz(age,5,1)/DBN_XZ(5,1) + &
								& (0.999_dp-cdf_xz)*size_draft_group_xz(age,6,1)/DBN_XZ(6,1)
			cdf_xz = cdf_xz + DBN_XZ(6,1)
		size_draft_group(age,6)   = (CDF_xz-0.999_dp)*size_draft_group_xz(age,6,1)/DBN_XZ(6,1) + size_draft_group_xz(age,7,1) + &
								&  size_draft_group_xz(age,8,1) + size_draft_group_xz(age,9,1)

		! Fraction Positive Welfare
	    	cdf_xz = sum(DBN_XZ(1:3,3))
		frac_pos_welfare_draft_group(age,1)   = frac_pos_welfare_draft_group_xz(age,1,3) + frac_pos_welfare_draft_group_xz(age,2,3) + &
								& frac_pos_welfare_draft_group_xz(age,3,3) + &
								& (0.40_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,4,3)/DBN_XZ(4,3)
			cdf_xz = cdf_xz + DBN_XZ(4,3)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(5:,3)) + sum(DBN_XZ(1:3,1))
		frac_pos_welfare_draft_group(age,2)   = (cdf_xz_low-0.40_dp)*frac_pos_welfare_draft_group_xz(age,4,3)/DBN_XZ(4,3) + &
								& frac_pos_welfare_draft_group_xz(age,5,3) + frac_pos_welfare_draft_group_xz(age,6,3) + &
								& frac_pos_welfare_draft_group_xz(age,7,3) + frac_pos_welfare_draft_group_xz(age,8,3) + &
								& frac_pos_welfare_draft_group_xz(age,9,3) + &
								& frac_pos_welfare_draft_group_xz(age,1,1) + frac_pos_welfare_draft_group_xz(age,2,1) + &
								& frac_pos_welfare_draft_group_xz(age,3,1) + &
								& (0.80_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,4,1)/DBN_XZ(4,1)
			cdf_xz = cdf_xz + DBN_XZ(4,1)
		frac_pos_welfare_draft_group(age,3)   = (cdf_xz-0.80_dp)*frac_pos_welfare_draft_group_xz(age,4,1)/DBN_XZ(4,1) + &
								&               (0.90_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,5,2)/DBN_XZ(5,2)
			cdf_xz = cdf_xz + DBN_XZ(5,2)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(6:,2)) 
		frac_pos_welfare_draft_group(age,4)   = (cdf_xz-0.90_dp)*frac_pos_welfare_draft_group_xz(age,5,2)/DBN_XZ(5,2) + &
								& frac_pos_welfare_draft_group_xz(age,6,2) + frac_pos_welfare_draft_group_xz(age,7,2) + &
								& frac_pos_welfare_draft_group_xz(age,8,2) + frac_pos_welfare_draft_group_xz(age,9,2) + &
								& (0.99_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,5,1)/DBN_XZ(5,1)
			cdf_xz = cdf_xz + DBN_XZ(5,1)
		frac_pos_welfare_draft_group(age,5)   = (cdf_xz-0.99_dp)*frac_pos_welfare_draft_group_xz(age,5,1)/DBN_XZ(5,1) + &
								&				(0.999_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,6,1)/DBN_XZ(6,1)
			cdf_xz = cdf_xz + DBN_XZ(6,1)
		frac_pos_welfare_draft_group(age,6)   = (CDF_xz-0.999_dp)*frac_pos_welfare_draft_group_xz(age,6,1)/DBN_XZ(6,1) + &
								&  frac_pos_welfare_draft_group_xz(age,7,1) + &
								&  frac_pos_welfare_draft_group_xz(age,8,1) + frac_pos_welfare_draft_group_xz(age,9,1)


		! Age Group 5
			! 00%-40%   : (z1,x3), (z2,x3), (z3,x3), (z4,x3)
			! 40%-80%   : (z4,x3), (z5,x3), (z6,x3)
			! 80%-90%   : (z6,x3), (z7,x3), (z8,x3), (z9,x3), (z1,x1), (z2,x1), (z3,x1), (z4,x1)
			! 90%-99%   : (z4,x1), (z5,x2), (z6,x2)
			! 99%-99.9% : (z6,x2), (z7,x2), (z8,x2), (z9,x2), (z5,x1)
			! 99.9%+    : (z5,x1), (z6,x1), (z7,x1), (z8,x1), (z9,x1)
			age = 5
			DBN_XZ = sum(sum(sum(sum(DBN_bench(draft_age_limit(age)+1:draft_age_limit(age+1),:,:,:,:,:),5),4),2),1) 
			DBN_XZ = DBN_XZ/sum(DBN_XZ)

		! Welfare
	    	cdf_xz = sum(DBN_XZ(1:3,3))
		CE_draft_group(age,1)   = ( DBN_XZ(1,3)*CE_draft_group_xz(age,1,3) + DBN_XZ(2,3)*CE_draft_group_xz(age,2,3) + &
								& DBN_XZ(3,3)*CE_draft_group_xz(age,3,3) + &
								& (0.40_dp-cdf_xz)*CE_draft_group_xz(age,4,3) )/0.40_dp 
			cdf_xz = cdf_xz + DBN_XZ(4,3)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + DBN_XZ(5,3)
		CE_draft_group(age,2)   = ( (cdf_xz_low-0.40_dp)*CE_draft_group_xz(age,4,3) + &
								& DBN_XZ(5,3)*CE_draft_group_xz(age,5,3) + &
								& (0.80_dp-cdf_xz)*CE_draft_group_xz(age,6,3) )/0.40_dp
			cdf_xz = cdf_xz + DBN_XZ(6,3)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(7:,3)) + sum(DBN_XZ(1:3,1)) 
		CE_draft_group(age,3)   = ( (cdf_xz_low-0.80_dp)*CE_draft_group_xz(age,6,3) + &
								& DBN_XZ(7,3)*CE_draft_group_xz(age,7,3) + DBN_XZ(8,3)*CE_draft_group_xz(age,8,3) + &
								& DBN_XZ(9,3)*CE_draft_group_xz(age,9,3) + &
								& DBN_XZ(1,1)*CE_draft_group_xz(age,1,1) + DBN_XZ(2,1)*CE_draft_group_xz(age,2,1) + &
								& DBN_XZ(3,1)*CE_draft_group_xz(age,3,1) + &
								& (0.90_dp-cdf_xz)*CE_draft_group_xz(age,4,1) )/0.10_dp
			cdf_xz = cdf_xz + DBN_XZ(4,1)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + DBN_XZ(5,2)
		CE_draft_group(age,4)   = ( (cdf_xz-0.90_dp)*CE_draft_group_xz(age,4,1) + &
								& DBN_XZ(5,2)*CE_draft_group_xz(age,5,2) + &
								& (0.99_dp-cdf_xz)*CE_draft_group_xz(age,6,2) )/0.09_dp
			cdf_xz = cdf_xz + DBN_XZ(6,2)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(7:,2))
		CE_draft_group(age,5)   = ( (cdf_xz-0.99_dp)*CE_draft_group_xz(age,6,2) + &
								& DBN_XZ(7,2)*CE_draft_group_xz(age,7,2) + DBN_XZ(8,2)*CE_draft_group_xz(age,8,2) + &
								& DBN_XZ(9,2)*CE_draft_group_xz(age,9,2) + &
								& (0.999_dp-cdf_xz)*CE_draft_group_xz(age,5,1) )/0.009_dp
			cdf_xz = cdf_xz + DBN_XZ(5,1)
		CE_draft_group(age,6)   = ( (CDF_xz-0.999_dp)*CE_draft_group_xz(age,5,1) + &
								& DBN_XZ(6,1)*CE_draft_group_xz(age,6,1) + DBN_XZ(7,1)*CE_draft_group_xz(age,7,1) + &
								& DBN_XZ(8,1)*CE_draft_group_xz(age,8,1) + DBN_XZ(9,1)*CE_draft_group_xz(age,9,1) )/0.001_dp

		! Size of Each group
	    	cdf_xz = sum(DBN_XZ(1:3,3))
		size_draft_group(age,1)   = size_draft_group_xz(age,1,3) + size_draft_group_xz(age,2,3) + &
								& size_draft_group_xz(age,3,3) + &
								& (0.40_dp-cdf_xz)*size_draft_group_xz(age,4,3)/DBN_XZ(4,3)
			cdf_xz = cdf_xz + DBN_XZ(4,3)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + DBN_XZ(5,3)
		size_draft_group(age,2)   = (cdf_xz_low-0.40_dp)*size_draft_group_xz(age,4,3)/DBN_XZ(4,3) + &
								& size_draft_group_xz(age,5,3) + &
								& (0.80_dp-cdf_xz)*size_draft_group_xz(age,6,3)/DBN_XZ(6,3)
			cdf_xz = cdf_xz + DBN_XZ(6,3)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(7:,3)) + sum(DBN_XZ(1:3,1)) 
		size_draft_group(age,3)   = (cdf_xz_low-0.80_dp)*size_draft_group_xz(age,6,3)/DBN_XZ(6,3) + &
								& size_draft_group_xz(age,7,3) + size_draft_group_xz(age,8,3) + &
								& size_draft_group_xz(age,9,3) + &
								& size_draft_group_xz(age,1,1) + size_draft_group_xz(age,2,1) + &
								& size_draft_group_xz(age,3,1) + &
								& (0.90_dp-cdf_xz)*size_draft_group_xz(age,4,1)/DBN_XZ(4,1)
			cdf_xz = cdf_xz + DBN_XZ(4,1)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + DBN_XZ(5,2)
		size_draft_group(age,4)   = (cdf_xz-0.90_dp)*size_draft_group_xz(age,4,1)/DBN_XZ(4,1) + &
								& size_draft_group_xz(age,5,2) + &
								& (0.99_dp-cdf_xz)*size_draft_group_xz(age,6,2)/DBN_XZ(6,2)
			cdf_xz = cdf_xz + DBN_XZ(6,2)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(7:,2))
		size_draft_group(age,5)   = (cdf_xz-0.99_dp)*size_draft_group_xz(age,6,2)/DBN_XZ(6,2) + &
								& size_draft_group_xz(age,7,2) + size_draft_group_xz(age,8,2) + &
								& size_draft_group_xz(age,9,2) + &
								& (0.999_dp-cdf_xz)*size_draft_group_xz(age,5,1)/DBN_XZ(5,1)
			cdf_xz = cdf_xz + DBN_XZ(5,1)
		size_draft_group(age,6)   = (CDF_xz-0.999_dp)*size_draft_group_xz(age,5,1)/DBN_XZ(5,1) + &
								& size_draft_group_xz(age,6,1) + size_draft_group_xz(age,7,1) + &
								& size_draft_group_xz(age,8,1) + size_draft_group_xz(age,9,1)

		! Fraction Positive Welfare
	    	cdf_xz = sum(DBN_XZ(1:3,3))
		frac_pos_welfare_draft_group(age,1)   = frac_pos_welfare_draft_group_xz(age,1,3) + frac_pos_welfare_draft_group_xz(age,2,3) + &
								& frac_pos_welfare_draft_group_xz(age,3,3) + &
								& (0.40_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,4,3)/DBN_XZ(4,3)
			cdf_xz = cdf_xz + DBN_XZ(4,3)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + DBN_XZ(5,3)
		frac_pos_welfare_draft_group(age,2)   = (cdf_xz_low-0.40_dp)*frac_pos_welfare_draft_group_xz(age,4,3)/DBN_XZ(4,3) + &
								& frac_pos_welfare_draft_group_xz(age,5,3) + &
								& (0.80_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,6,3)/DBN_XZ(6,3)
			cdf_xz = cdf_xz + DBN_XZ(6,3)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(7:,3)) + sum(DBN_XZ(1:3,1)) 
		frac_pos_welfare_draft_group(age,3)   = (cdf_xz_low-0.80_dp)*frac_pos_welfare_draft_group_xz(age,6,3)/DBN_XZ(6,3) + &
								& frac_pos_welfare_draft_group_xz(age,7,3) + frac_pos_welfare_draft_group_xz(age,8,3) + &
								& frac_pos_welfare_draft_group_xz(age,9,3) + &
								& frac_pos_welfare_draft_group_xz(age,1,1) + frac_pos_welfare_draft_group_xz(age,2,1) + &
								& frac_pos_welfare_draft_group_xz(age,3,1) + &
								& (0.90_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,4,1)/DBN_XZ(4,1)
			cdf_xz = cdf_xz + DBN_XZ(4,1)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + DBN_XZ(5,2)
		frac_pos_welfare_draft_group(age,4)   = (cdf_xz-0.90_dp)*frac_pos_welfare_draft_group_xz(age,4,1)/DBN_XZ(4,1) + &
								& frac_pos_welfare_draft_group_xz(age,5,2) + &
								& (0.99_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,6,2)/DBN_XZ(6,2)
			cdf_xz = cdf_xz + DBN_XZ(6,2)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(7:,2))
		frac_pos_welfare_draft_group(age,5)   = (cdf_xz-0.99_dp)*frac_pos_welfare_draft_group_xz(age,6,2)/DBN_XZ(6,2) + &
								& frac_pos_welfare_draft_group_xz(age,7,2) + frac_pos_welfare_draft_group_xz(age,8,2) + &
								& frac_pos_welfare_draft_group_xz(age,9,2) + &
								& (0.999_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(age,5,1)/DBN_XZ(5,1)
			cdf_xz = cdf_xz + DBN_XZ(5,1)
		frac_pos_welfare_draft_group(age,6)   = (CDF_xz-0.999_dp)*frac_pos_welfare_draft_group_xz(age,5,1)/DBN_XZ(5,1) + &
								& frac_pos_welfare_draft_group_xz(age,6,1) + frac_pos_welfare_draft_group_xz(age,7,1) + &
								& frac_pos_welfare_draft_group_xz(age,8,1) + frac_pos_welfare_draft_group_xz(age,9,1)

	! Fill in additional unused columns
		CE_draft_group(:,7)   = 0
		CE_draft_group(:,8)   = 0

		size_draft_group(:,7)   = 0
		size_draft_group(:,8)   = 0

		frac_pos_welfare_draft_group(:,7)   = 0
		frac_pos_welfare_draft_group(:,8)   = 0


	! Fix fractions
	frac_pos_welfare_draft_group = 100.0_dp*frac_pos_welfare_draft_group/size_draft_group

    OPEN (UNIT=80, FILE=trim(Result_Folder)//'draft_group_CE_axz.txt', STATUS='replace') 
    OPEN (UNIT=81, FILE=trim(Result_Folder)//'draft_group_size_axz.txt', STATUS='replace') 
    OPEN (UNIT=82, FILE=trim(Result_Folder)//'draft_group_fpos_welfare_axz.txt', STATUS='replace') 
	do age = 1,draft_age_category
	    WRITE  (UNIT=80, FMT=*)  CE_draft_group(age,:)
	    WRITE  (UNIT=81, FMT=*)  size_draft_group(age,:)
	    WRITE  (UNIT=82, FMT=*)  frac_pos_welfare_draft_group(age,:)
	ENDDO
	close(unit=80)
	close(unit=81)
	close(unit=82)


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! End of Compute Welfare Gain
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	print*, ' End of compute welfare gain'; print*, ' '

	! Deallocate policy functions on adjusted grid (so that they can be allocated later)
	if (allocated(YGRID_t)) then 
		deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )
	endif 
	

END SUBROUTINE  COMPUTE_WELFARE_GAIN


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE COMPUTE_WELFARE_GAIN_TRANSITION()
	IMPLICIT NONE
	real(DP), dimension(MaxAge):: CumDiscountF
	REAL(DP), dimension(draft_age_category,nz) :: CE_draft_group_z,  size_draft_group_z, frac_pos_welfare_draft_group_z
	REAL(DP), dimension(draft_age_category,draft_z_category) :: CE_draft_group,  size_draft_group, frac_pos_welfare_draft_group
	REAL(DP), dimension(nz) :: DBN_Z, CDF_Z 
	REAL(DP), dimension(nz,nx) :: DBN_XZ
	INTEGER , dimension(draft_age_category+1) :: draft_age_limit
	INTEGER :: age2, z2
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: Value_mat, Bq_Value_mat

	allocate( Value_mat(  	MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Bq_Value_mat( MaxAge,na,nz,nlambda,ne,nx) )

	print*,' '
	print*,'Computing Welfare Gain for Transition'
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Routine Set Up
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	print*,'	 Set up'
	! Age Brackets
		draft_age_limit = [0, 1, 15, 30, 45, MaxAge ] 

	! Discount factor
		CumDiscountF(MaxAge)=1.0_DP
		DO age=MaxAge-1,1,-1
		    CumDiscountF(age) = 1.0_DP + beta * survP(age) *CumDiscountF(age+1) 
		ENDDO

	! Select relevant values of value function
		DO age=1,MaxAge
			Value_mat(age,:,:,:,:,:)    = ValueFunction_tr(age,:,:,:,:,:,age)
			Bq_Value_mat(age,:,:,:,:,:) = BQ_Value_tr(age,:,:,:,:,:,age)
		ENDDO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Measuring Consumption Equivalent Welfare
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! Consumption Equivalent Welfare: CE 1
	print*,'	 Consumption Equivalent: CE 1'
		DO age=1,MaxAge
			if (Log_Switch.eqv..true.) then 
		    	CE1_tr(age,:,:,:,:,:)= & 
		    		& exp((Value_mat(age,:,:,:,:,:)-ValueFunction_Bench(age,:,:,:,:,:))/CumDiscountF(age))-1.0_DP
		    else 
		    	CE1_tr(age,:,:,:,:,:)=((Value_mat(age,:,:,:,:,:)-BQ_Value_mat(age,:,:,:,:,:))/&
    									& (ValueFunction_Bench(age,:,:,:,:,:)-BQ_Value_mat(age,:,:,:,:,:)) ) &
                        				&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP
		    end if 
		ENDDO

		! Aggregates 
		CE1_nb_tr  = 100.0_DP*sum(CE1_tr(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
		CE1_pop_tr = 100.0_DP*sum(CE1_tr*DBN_bench)/sum(DBN_bench)

	! Consumption Equivalent Welfare: CE 2
	print*,'	 Consumption Equivalent: CE 2'
		CE2_nb_tr  = 100.0_dp * (( sum((Value_mat(1,:,:,:,:,:)-BQ_Value_mat(1,:,:,:,:,:))*DBN_bench(1,:,:,:,:,:)) / &
				&               sum((ValueFunction_Bench(1,:,:,:,:,:)-BQ_Value_mat(1,:,:,:,:,:))*DBN_bench(1,:,:,:,:,:)) ) &
				& ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
		CE2_pop_tr = 100.0_dp*(( sum((Value_mat-BQ_Value_mat)*DBN_bench) / sum((ValueFunction_Bench-BQ_Value_mat)*DBN_bench)  ) &
		                                &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

	! Write Aggregate Results
		OPEN (UNIT=50, FILE=trim(Result_Folder)//'CE_Transition.txt', STATUS='replace') 
			WRITE  (UNIT=50, FMT=*) 'Measures of Consumption Equivalent Welfare for Transition'
			WRITE  (UNIT=50, FMT=*) 'Welfare measured for cohorts alive at time of policy change'
			WRITE  (UNIT=50, FMT=*) ' '
			WRITE  (UNIT=50, FMT=*) 'Aggregate Measures'
			WRITE  (UNIT=50, FMT=*) ' '
			WRITE  (UNIT=50, FMT=*) 'CE 1: Average of welfare gain across agents'
			WRITE  (UNIT=50, FMT=*) 'CE1_nb =',CE1_nb_tr
			WRITE  (UNIT=50, FMT=*) 'CE1_pop =',CE1_pop_tr
			WRITE  (UNIT=50, FMT=*) ' '
			WRITE  (UNIT=50, FMT=*) 'CE 2: Welfare gain average agent'
			WRITE  (UNIT=50, FMT=*) 'CE2_nb =',CE2_nb_tr
			WRITE  (UNIT=50, FMT=*) 'CE2_pop =',CE2_pop_tr
			WRITE  (UNIT=50, FMT=*) ' '
			WRITE  (UNIT=50, FMT=*) 'All unites in percentage points (already multiplied by 100)'
		close (unit=50)	

	print*,''
	print*,'---------------------------'
	print*, 'CE_1: Average of welfare gain across agents'
	print '(A,F7.3)', 'CE1_nb  =',CE1_nb_tr
	print '(A,F7.3)', 'CE1_pop =',CE1_pop_tr
	print*, ' '
	print*, 'CE_2: Welfare gain average agent'
	print '(A,F7.3)', 'CE2_nb  =',CE2_nb_tr
	print '(A,F7.3)', 'CE2_pop =',CE2_pop_tr
	print*,'---------------------------'
	print*,''


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Draft Tables
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	frac_pos_welfare_draft_group_z = 0.0_dp 

	do zi  = 1,nz
	do age = 1,draft_age_category
		size_draft_group_z(age,zi) = &
			& sum(DBN_bench(draft_age_limit(age)+1:draft_age_limit(age+1),:,zi,:,:,:))

		CE_draft_group_z(age,zi) =  100.0_dp* &
			& sum(CE1_tr(draft_age_limit(age)+1:draft_age_limit(age+1),:,zi,:,:,:)* &
            &     DBN_bench(draft_age_limit(age)+1:draft_age_limit(age+1),:,zi,:,:,:))/&
            & sum(DBN_bench(draft_age_limit(age)+1:draft_age_limit(age+1),:,zi,:,:,:))

        do xi=1,nx
	    do ei=1,ne
	    do lambdai=1,nlambda
	    do ai=1,na
	    do age2=draft_age_limit(age)+1,draft_age_limit(age+1)
	    	If ( CE1_tr(age2,ai,zi,lambdai,ei,xi) .ge. 0.0_DP) then
        	frac_pos_welfare_draft_group_z(age,zi) = frac_pos_welfare_draft_group_z(age,zi) + DBN_bench(age2,ai,zi,lambdai,ei,xi)
        	endif 
    	enddo 
    	enddo 
    	enddo 
    	enddo 
    	enddo  
	enddo
	enddo 

	DBN_Z = sum(sum(sum(sum(sum(DBN_bench,6),5),4),2),1) 

	! Size of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
	size_draft_group = Draft_Table(size_draft_group_z,DBN_z,.true.)

	! CE of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
	CE_draft_group = Draft_Table(CE_draft_group_z,DBN_z,.false.)

	! Frac. pos. welfare by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
	frac_pos_welfare_draft_group = Draft_Table(frac_pos_welfare_draft_group_z,DBN_z,.true.)
		! Fix fractions
		frac_pos_welfare_draft_group = 100.0_dp*frac_pos_welfare_draft_group/size_draft_group

    OPEN (UNIT=80, FILE=trim(Result_Folder)//'draft_group_CE_tr.txt', STATUS='replace') 
    OPEN (UNIT=81, FILE=trim(Result_Folder)//'draft_group_size_tr.txt', STATUS='replace') 
    OPEN (UNIT=82, FILE=trim(Result_Folder)//'draft_group_fpos_welfare_tr.txt', STATUS='replace') 
	do age = 1,draft_age_category
	    WRITE  (UNIT=80, FMT=*)  CE_draft_group(age,:)
	    WRITE  (UNIT=81, FMT=*)  size_draft_group(age,:)
	    WRITE  (UNIT=82, FMT=*)  frac_pos_welfare_draft_group(age,:)
	ENDDO
	close(unit=80); close(unit=81); close(unit=82)


END SUBROUTINE  COMPUTE_WELFARE_GAIN_TRANSITION


!========================================================================================
!========================================================================================
!========================================================================================

SUBROUTINE COMPUTE_WELFARE_DECOMPOSITION
	IMPLICIT NONE
	REAL(dp) :: size_nb
	REAL(dp) :: C_bench, C_exp, C_nb_bench, C_nb_exp, H_bench, H_exp, H_NB_bench, H_NB_exp, BQ_bench, BQ_exp, BQ_NB_bench, BQ_NB_exp
	REAL(dp) :: CE1_nb, CE1_nb_c, CE1_nb_cl, CE1_nb_cd, CE1_nb_h, CE1_nb_hl, CE1_nb_hd, CE1_nb_b, CE1_nb_bl, CE1_nb_bd
	REAL(dp) :: CE1_pop, CE1_pop_c, CE1_pop_cl, CE1_pop_cd, CE1_pop_h, CE1_pop_hl, CE1_pop_hd, CE1_pop_b, CE1_pop_bl, CE1_pop_bd
	REAL(dp) :: CE1_nb_ch, CE1_nb_chl, CE1_nb_chd, CE1_pop_ch, CE1_pop_chl, CE1_pop_chd
	REAL(dp) :: CE2_nb, CE2_nb_c, CE2_nb_cl, CE2_nb_cd, CE2_nb_h, CE2_nb_hl, CE2_nb_hd, CE2_nb_b, CE2_nb_bl, CE2_nb_bd
	REAL(dp) :: CE2_pop, CE2_pop_c, CE2_pop_cl, CE2_pop_cd, CE2_pop_h, CE2_pop_hl, CE2_pop_hd, CE2_pop_b, CE2_pop_bl, CE2_pop_bd
	REAL(dp) :: CE2_nb_ch, CE2_nb_chl, CE2_nb_chd, CE2_pop_ch, CE2_pop_chl, CE2_pop_chd
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: Value_aux, Bq_Value_aux
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: CE1_mat, CE1_c_mat, CE1_h_mat, CE1_b_mat, CE1_ch_mat

	allocate( Value_aux(      MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Bq_Value_aux(   MaxAge,na,nz,nlambda,ne,nx) )
	allocate( CE1_mat(   	  MaxAge,na,nz,nlambda,ne,nx) )
	allocate( CE1_c_mat(   	  MaxAge,na,nz,nlambda,ne,nx) )
	allocate( CE1_h_mat(   	  MaxAge,na,nz,nlambda,ne,nx) )
	allocate( CE1_b_mat(   	  MaxAge,na,nz,nlambda,ne,nx) )
	allocate( CE1_ch_mat(  	  MaxAge,na,nz,nlambda,ne,nx) )

	! Size of new borns 
	size_nb     = sum(DBN_bench(1,:,:,:,:,:))
	
	! Define benchmark and experiment variales
	C_bench 	= sum(Cons_bench*DBN_bench)
	C_exp   	= sum(Cons_exp*DBN_exp)
	C_nb_bench 	= sum(Cons_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
	C_nb_exp   	= sum(Cons_exp(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:))/sum(DBN_exp(1,:,:,:,:,:))

	H_bench 	= sum(Hours_bench*DBN_bench)
	H_exp   	= sum(Hours_exp*DBN_exp)
	H_NB_bench 	= sum(Hours_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
	H_NB_exp   	= sum(Hours_exp(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:))/sum(DBN_exp(1,:,:,:,:,:)) 

	BQ_bench 	= sum(Aprime_bench*DBN_bench)
	BQ_exp   	= sum(Aprime_exp*DBN_exp)
	BQ_NB_bench = sum(Aprime_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
	BQ_NB_exp   = sum(Aprime_exp(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:))/sum(DBN_exp(1,:,:,:,:,:)) 

	! Auxiliary Value Functions
		! Note: we are keeping Aprime at bench values, which affects bequest and expected values trhough 
		!		changes in states. An alternative is to use Aprime_exp and adjust net out the value of 
		!		bequests by subtracting Bq_Value_exp and adding Bq_value_bench from the auxiliary values
		! Change Consumption 
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons_exp,Hours_bench,Aprime_exp,ValueFunction,Bq_Value 	 )
			! Adjust Value to keep bequest value at benchmark level
			ValueFunction = ValueFunction - Bq_Value + Bq_Value_bench
		! Change Consumption and Leisure
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons_exp,Hours_exp  ,Aprime_exp,Value_aux    ,Bq_Value_aux)
			! Adjust Value to keep bequest value at benchmark level
			Value_aux = Value_aux - Bq_Value_aux + Bq_Value_bench

	! CE1 Measures 
		CE1_mat =((ValueFunction_exp  -Bq_Value_bench)/&
				& (ValueFunction_Bench-Bq_Value_bench) )**( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP

		! Decomposition: Consumption
			! Total
			CE1_c_mat  =((ValueFunction      -Bq_Value_bench)/&
					   & (ValueFunction_Bench-Bq_Value_bench) )**( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP
			! Level
			CE1_pop_cl = 100.0_dp*( C_exp/C_bench-1.0_dp )
			CE1_nb_cl  = 100.0_dp*( C_exp/C_bench-1.0_dp )! 100.0_dp*( C_nb_exp/C_nb_bench-1.0_dp )

		! Decomposition: Leisure
			! Total
			CE1_h_mat  =((Value_aux    -Bq_Value_bench)/&
					   & (ValueFunction-Bq_Value_bench) )**( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP
			! Level
			CE1_pop_hl = 100.0_dp*( ((1.0_dp-H_exp   )/(1.0_dp-H_bench   ))**((1.0_dp-gamma)/gamma) -1.0_dp )
			CE1_nb_hl  = 100.0_dp*( ((1.0_dp-H_exp   )/(1.0_dp-H_bench   ))**((1.0_dp-gamma)/gamma) -1.0_dp ) ! 100.0_dp*( ((1.0_dp-H_nb_exp)/(1.0_dp-H_nb_bench))**((1.0_dp-gamma)/gamma) -1.0_dp )

		! Decomposition: Bequests
			! Total
			CE1_b_mat  =((ValueFunction_exp-Bq_Value_bench)/&
					   & (Value_aux        -Bq_Value_bench) )**( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP
			! Level
			CE1_pop_bl = 100.0_dp*( BQ_exp/BQ_bench-1.0_dp )
			CE1_nb_bl  = 100.0_dp*( BQ_exp/BQ_bench-1.0_dp ) ! 100.0_dp*( BQ_nb_exp/BQ_nb_bench-1.0_dp )

		! Decomposition: Consumption and Leisure
			CE1_ch_mat = ((Value_aux    -Bq_Value_bench)/&
					   & (ValueFunction_Bench-Bq_Value_bench) )**( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP

			! Level
			CE1_pop_chl = 100.0_dp*( C_exp/C_bench*((1.0_dp-H_exp)/(1.0_dp-H_bench))**((1.0_dp-gamma)/gamma) -1.0_dp )
			CE1_nb_chl  = 100.0_dp*( C_exp/C_bench*((1.0_dp-H_exp)/(1.0_dp-H_bench))**((1.0_dp-gamma)/gamma) -1.0_dp ) ! 100.0_dp*( C_nb_exp/C_nb_bench*((1.0_dp-H_nb_exp)/(1.0_dp-H_nb_bench))**((1.0_dp-gamma)/gamma) -1.0_dp )

		! Aggregating CE1 for population or newborns 
			CE1_nb		= 100.0_dp*sum( CE1_mat(1,:,:,:,:,:)  *DBN_bench(1,:,:,:,:,:) )/size_nb
			CE1_nb_c	= 100.0_dp*sum( CE1_c_mat(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:) )/size_nb
			CE1_nb_cd	= 100.0_dp*( (CE1_nb_c/100.0_dp+1.0_dp)/(CE1_nb_cl/100.0_dp+1.0_dp) - 1.0_dp )
			CE1_nb_h	= 100.0_dp*sum( CE1_h_mat(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:) )/size_nb
			CE1_nb_hd	= 100.0_dp*( (CE1_nb_h/100.0_dp+1.0_dp)/(CE1_nb_hl/100.0_dp+1.0_dp) - 1.0_dp )
			CE1_nb_b	= 100.0_dp*sum( CE1_b_mat(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:) )/size_nb
			CE1_nb_bd	= 100.0_dp*( (CE1_nb_b/100.0_dp+1.0_dp)/(CE1_nb_bl/100.0_dp+1.0_dp) - 1.0_dp )
			CE1_nb_ch	= 100.0_dp*sum( CE1_ch_mat(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:) )/size_nb
			CE1_nb_chd	= 100.0_dp*( (CE1_nb_ch/100.0_dp+1.0_dp)/(CE1_nb_chl/100.0_dp+1.0_dp) - 1.0_dp )

			CE1_pop		= 100.0_dp*sum( CE1_mat  *DBN_bench )
			CE1_pop_c	= 100.0_dp*sum( CE1_c_mat*DBN_bench )
			CE1_pop_cd	= 100.0_dp*( (CE1_pop_c/100.0_dp+1.0_dp)/(CE1_pop_cl/100.0_dp+1.0_dp) - 1.0_dp )
			CE1_pop_h	= 100.0_dp*sum( CE1_h_mat*DBN_bench )
			CE1_pop_hd	= 100.0_dp*( (CE1_pop_h/100.0_dp+1.0_dp)/(CE1_pop_hl/100.0_dp+1.0_dp) - 1.0_dp )
			CE1_pop_b	= 100.0_dp*sum( CE1_b_mat*DBN_bench )
			CE1_pop_bd	= 100.0_dp*( (CE1_pop_b/100.0_dp+1.0_dp)/(CE1_pop_bl/100.0_dp+1.0_dp) - 1.0_dp )
			CE1_pop_ch	= 100.0_dp*sum( CE1_ch_mat*DBN_bench )
			CE1_pop_chd	= 100.0_dp*( (CE1_pop_ch/100.0_dp+1.0_dp)/(CE1_pop_chl/100.0_dp+1.0_dp) - 1.0_dp )



	! CE2 Measures 
		CE2_nb  =  100.0_dp*(( (sum(ValueFunction_exp(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:)) - &
	    					&  sum(Bq_Value_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))) / &
	                        &  sum((ValueFunction_Bench(1,:,:,:,:,:)-Bq_Value_bench(1,:,:,:,:,:))*DBN_bench(1,:,:,:,:,:))  ) &
	                        &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

		CE2_pop = 100.0_dp*(( (sum(ValueFunction_exp*DBN_exp)-sum(Bq_Value_Bench*DBN_bench)) / &
	    					& sum((ValueFunction_Bench-Bq_Value_Bench)*DBN_bench)  ) &
	                        &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

		! Decomposition: Consumption
			! Total (note that auxiliary value functions are weighted with benchmark distribution)
			CE2_nb_c   =  100.0_dp*(( (sum(ValueFunction(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:)) - &
	    						&  sum(Bq_Value_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))) / &
	                        	&  sum((ValueFunction_Bench(1,:,:,:,:,:)-Bq_Value_bench(1,:,:,:,:,:))*DBN_bench(1,:,:,:,:,:))  ) &
	                        	&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

			CE2_pop_c  = 100.0_dp*(( (sum(ValueFunction*DBN_exp)-sum(Bq_Value_Bench*DBN_bench)) / &
	    						& sum((ValueFunction_Bench-Bq_Value_Bench)*DBN_bench)  ) &
	                        	&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
			! Level
			CE2_nb_cl  = 100.0_dp*( C_exp/C_bench-1.0_dp ) ! 100.0_dp*( C_nb_exp/C_nb_bench-1.0_dp )
			CE2_pop_cl = 100.0_dp*( C_exp/C_bench-1.0_dp )
			! Distribution
			CE2_nb_cd  = 100.0_dp*( (CE2_nb_c/100.0_dp+1.0_dp)/(CE2_nb_cl/100.0_dp+1.0_dp) - 1.0_dp )
			CE2_pop_cd = 100.0_dp*( (CE2_pop_c/100.0_dp+1.0_dp)/(CE2_pop_cl/100.0_dp+1.0_dp) - 1.0_dp )

		! Decomposition: Leisure
			! Total (note that auxiliary value functions are weighted with benchmark distribution)
			CE2_nb_h   =  100.0_dp*(( (sum(Value_aux(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:)) - &
	    						&  sum(Bq_Value_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))) / &
	                        	&  (sum(ValueFunction(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:)) - &
	    						&  sum(Bq_Value_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:)))  ) &
	                        	&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

			CE2_pop_h  = 100.0_dp*(( (sum(Value_aux*DBN_exp)-sum(Bq_Value_Bench*DBN_bench)) / &
	    						& (sum(ValueFunction*DBN_exp)-sum(Bq_Value_Bench*DBN_bench))  ) &
	                        	&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
			! Level
			CE2_nb_hl  = 100.0_dp*( ((1.0_dp-H_exp   )/(1.0_dp-H_bench   ))**((1.0_dp-gamma)/gamma) -1.0_dp ) ! 100.0_dp*( ((1.0_dp-H_nb_exp)/(1.0_dp-H_nb_bench))**((1.0_dp-gamma)/gamma) -1.0_dp )
			CE2_pop_hl = 100.0_dp*( ((1.0_dp-H_exp   )/(1.0_dp-H_bench   ))**((1.0_dp-gamma)/gamma) -1.0_dp )
			! Distribution
			CE2_nb_hd  = 100.0_dp*( (CE2_nb_h/100.0_dp+1.0_dp)/(CE2_nb_hl/100.0_dp+1.0_dp) - 1.0_dp )
			CE2_pop_hd = 100.0_dp*( (CE2_pop_h/100.0_dp+1.0_dp)/(CE2_pop_hl/100.0_dp+1.0_dp) - 1.0_dp )

		! Decomposition: Bequests
			! Total (note that auxiliary value functions are weighted with experiment distribution)
			CE2_nb_b   =  100.0_dp*(( (sum(ValueFunction_exp(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:)) - &
	    						&  sum(Bq_Value_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))) / &
	                        	&  (sum(Value_aux(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:)) - &
	    						&  sum(Bq_Value_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:)))  ) &
	                        	&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

			CE2_pop_b  = 100.0_dp*(( (sum(ValueFunction_exp*DBN_exp)-sum(Bq_Value_Bench*DBN_bench)) / &
	    						& (sum(Value_aux*DBN_exp)-sum(Bq_Value_Bench*DBN_bench))  ) &
	                        	&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

			! Level
			CE2_nb_bl  = 100.0_dp*( BQ_exp/BQ_bench-1.0_dp ) ! 100.0_dp*( BQ_nb_exp/BQ_nb_bench-1.0_dp )
			CE2_pop_bl = 100.0_dp*( BQ_exp/BQ_bench-1.0_dp )
			! Distribution
			CE2_nb_bd  = 100.0_dp*( (CE2_nb_b/100.0_dp+1.0_dp)/(CE2_nb_bl/100.0_dp+1.0_dp) - 1.0_dp )
			CE2_pop_bd = 100.0_dp*( (CE2_pop_b/100.0_dp+1.0_dp)/(CE2_pop_bl/100.0_dp+1.0_dp) - 1.0_dp )

		! Decomposition: Consumption and Leisure
			! Total (note that auxiliary value functions are weighted with benchmark distribution)
			CE2_nb_ch   =  100.0_dp*(( (sum(Value_aux(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:)) - &
	    						&  sum(Bq_Value_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))) / &
	                        	&  sum((ValueFunction_Bench(1,:,:,:,:,:)-Bq_Value_bench(1,:,:,:,:,:))*DBN_bench(1,:,:,:,:,:))  ) &
	                        	&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

			CE2_pop_ch  = 100.0_dp*(( (sum(Value_aux*DBN_exp)-sum(Bq_Value_Bench*DBN_bench)) / &
	    						& sum((ValueFunction_Bench-Bq_Value_Bench)*DBN_bench)  ) &
	                        	&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
			! Level
			CE2_pop_chl = 100.0_dp*( C_exp/C_bench*((1.0_dp-H_exp)/(1.0_dp-H_bench))**((1.0_dp-gamma)/gamma) -1.0_dp )
			CE2_nb_chl  = 100.0_dp*( C_exp/C_bench*((1.0_dp-H_exp)/(1.0_dp-H_bench))**((1.0_dp-gamma)/gamma) -1.0_dp ) ! 100.0_dp*( C_nb_exp/C_nb_bench*((1.0_dp-H_nb_exp)/(1.0_dp-H_nb_bench))**((1.0_dp-gamma)/gamma) -1.0_dp )
			! Distribution
			CE2_nb_chd  = 100.0_dp*( (CE2_nb_ch/100.0_dp+1.0_dp)/(CE2_nb_chl/100.0_dp+1.0_dp) - 1.0_dp )
			CE2_pop_chd = 100.0_dp*( (CE2_pop_ch/100.0_dp+1.0_dp)/(CE2_pop_chl/100.0_dp+1.0_dp) - 1.0_dp )




	! Tables 
		OPEN  (UNIT=1,  FILE=trim(Result_Folder)//'CE_Decomposition.txt'  , STATUS='replace')
		WRITE (UNIT=1,  FMT=*) ' '
		WRITE (UNIT=1,  FMT=*) 'Decomposition: Consumption Equivalent Welfare'
		WRITE (UNIT=1,  FMT=*) ' '
		WRITE (UNIT=1,  FMT=*) 'CE1 '	, 'NB '		, 'Pop '	 , 'Test ' 	
		WRITE (UNIT=1,  FMT=*) 'CE1'	, CE1_nb 	, CE1_pop    , 100*((1+CE1_nb_c/100)*(1+CE1_nb_h/100)*(1+CE1_nb_b/100)-1)
		WRITE (UNIT=1,  FMT=*) 'CE1_c'	, CE1_nb_c 	, CE1_pop_c   
		WRITE (UNIT=1,  FMT=*) 'CE1_h'	, CE1_nb_h 	, CE1_pop_h  
		WRITE (UNIT=1,  FMT=*) 'CE1_b'	, CE1_nb_b 	, CE1_pop_b 
		WRITE (UNIT=1,  FMT=*) 'CE1_ch'	, CE1_nb_ch , CE1_pop_ch
		WRITE (UNIT=1,  FMT=*) '------'
		WRITE (UNIT=1,  FMT=*) 'CE1_c'	, CE1_nb_c 	, CE1_pop_c  , 100*((1+CE1_nb_cl/100)*(1+CE1_nb_cd/100)-1)				
		WRITE (UNIT=1,  FMT=*) 'CE1_cl'	, CE1_nb_cl	, CE1_pop_cl
		WRITE (UNIT=1,  FMT=*) 'CE1_cd'	, CE1_nb_cd , CE1_pop_cd
		WRITE (UNIT=1,  FMT=*) '------'
		WRITE (UNIT=1,  FMT=*) 'CE1_h'	, CE1_nb_h 	, CE1_pop_h  , 100*((1+CE1_nb_hl/100)*(1+CE1_nb_hd/100)-1)		
		WRITE (UNIT=1,  FMT=*) 'CE1_hl'	, CE1_nb_hl , CE1_pop_hl
		WRITE (UNIT=1,  FMT=*) 'CE1_hd'	, CE1_nb_hd , CE1_pop_hd 
		WRITE (UNIT=1,  FMT=*) '------'
		WRITE (UNIT=1,  FMT=*) 'CE1_b'	, CE1_nb_b 	, CE1_pop_b  , 100*((1+CE1_nb_bl/100)*(1+CE1_nb_bd/100)-1)		
		! WRITE (UNIT=1,  FMT=*) 'CE1_bl'	, CE1_nb_bl , CE1_pop_bl
		! WRITE (UNIT=1,  FMT=*) 'CE1_bd'	, CE1_nb_bd , CE1_pop_bd 
		WRITE (UNIT=1,  FMT=*) '------'
		WRITE (UNIT=1,  FMT=*) 'CE1_ch'	, CE1_nb_ch  , CE1_pop_ch  , 100*((1+CE1_nb_chl/100)*(1+CE1_nb_chd/100)-1)		
		WRITE (UNIT=1,  FMT=*) 'CE1_chl', CE1_nb_chl , CE1_pop_chl
		WRITE (UNIT=1,  FMT=*) 'CE1_chd', CE1_nb_chd , CE1_pop_chd
		WRITE (UNIT=1,  FMT=*) ' '
		WRITE (UNIT=1,  FMT=*) '-------------------------------'
		WRITE (UNIT=1,  FMT=*) ' '
		WRITE (UNIT=1,  FMT=*) 'CE2 '	, 'NB '		, 'Pop '	 , 'Test ' 	
		WRITE (UNIT=1,  FMT=*) 'CE2'	, CE2_nb 	, CE2_pop 	 , 100*((1+CE2_nb_c/100)*(1+CE2_nb_h/100)*(1+CE2_nb_b/100)-1)
		WRITE (UNIT=1,  FMT=*) 'CE2_c'	, CE2_nb_c 	, CE2_pop_c  , 100*((1+CE2_nb_cl/100)*(1+CE2_nb_cd/100)-1)		
		WRITE (UNIT=1,  FMT=*) 'CE2_h'	, CE2_nb_h 	, CE2_pop_h  , 100*((1+CE2_nb_hl/100)*(1+CE2_nb_hd/100)-1)
		WRITE (UNIT=1,  FMT=*) 'CE2_b'	, CE2_nb_b 	, CE2_pop_b  , 100*((1+CE2_nb_bl/100)*(1+CE2_nb_bd/100)-1)
		WRITE (UNIT=1,  FMT=*) 'CE2_ch'	, CE2_nb_ch , CE2_pop_ch
		WRITE (UNIT=1,  FMT=*) '------'
		WRITE (UNIT=1,  FMT=*) 'CE2_c'	, CE2_nb_c 	, CE2_pop_c  , 100*((1+CE2_nb_cl/100)*(1+CE2_nb_cd/100)-1)		
		WRITE (UNIT=1,  FMT=*) 'CE2_cl'	, CE2_nb_cl	, CE2_pop_cl
		WRITE (UNIT=1,  FMT=*) 'CE2_cd'	, CE2_nb_cd , CE2_pop_cd
		WRITE (UNIT=1,  FMT=*) '------'
		WRITE (UNIT=1,  FMT=*) 'CE2_h'	, CE2_nb_h 	, CE2_pop_h  , 100*((1+CE2_nb_hl/100)*(1+CE2_nb_hd/100)-1)
		WRITE (UNIT=1,  FMT=*) 'CE2_hl'	, CE2_nb_hl , CE2_pop_hl
		WRITE (UNIT=1,  FMT=*) 'CE2_hd'	, CE2_nb_hd , CE2_pop_hd 
		WRITE (UNIT=1,  FMT=*) '------'
		WRITE (UNIT=1,  FMT=*) 'CE2_b'	, CE2_nb_b 	, CE2_pop_b  , 100*((1+CE2_nb_bl/100)*(1+CE2_nb_bd/100)-1)		
		! WRITE (UNIT=1,  FMT=*) 'CE2_bl'	, CE2_nb_bl , CE2_pop_bl
		! WRITE (UNIT=1,  FMT=*) 'CE2_bd'	, CE2_nb_bd , CE2_pop_bd 
		WRITE (UNIT=1,  FMT=*) '------'
		WRITE (UNIT=1,  FMT=*) 'CE2_ch'	, CE2_nb_ch  , CE2_pop_ch  , 100*((1+CE2_nb_chl/100)*(1+CE2_nb_chd/100)-1)		
		WRITE (UNIT=1,  FMT=*) 'CE2_chl', CE2_nb_chl , CE2_pop_chl
		WRITE (UNIT=1,  FMT=*) 'CE2_chd', CE2_nb_chd , CE2_pop_chd
		WRITE (UNIT=1,  FMT=*) ' '
		CLOSE (unit=1)

		print*, ' '
		print*, '-----------------------------------------------------------------------------------'
		print*, ' '
		print*, 'Aggregates'
		print '(A,F7.3,F7.3,F6.2)', 'Consumption: ',C_bench,C_Exp,100.0_dp*(C_exp/C_bench-1.0_dp)
		print '(A,F7.3,F7.3,F6.2)', 'Hours      : ',H_bench,H_Exp,100.0_dp*(H_exp/H_bench-1.0_dp)
		print '(A,F7.3,F7.3,F6.2)', 'Leisure    : ',(1.0_dp-H_bench),(1.0_dp-H_Exp),100.0_dp*((1.0_dp-H_exp)/(1.0_dp-H_bench)-1.0_dp)
		print '(A,F7.3,F7.3,F6.2)', 'Bequest    : ',Bq_bench,Bq_Exp,100.0_dp*(Bq_exp/Bq_bench-1.0_dp)
		print*, ' '
		print*, 'Decomposition: Consumption Equivalent Welfare'
		print*, ' '
		print*, '	CE1 '	, 'NB '		, 'Pop '	 , 'Test ' 
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1'	, CE1_nb 	, CE1_pop   , 100*((1+CE1_nb_c/100)*(1+CE1_nb_h/100)*(1+CE1_nb_b/100)-1)
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_c'	, CE1_nb_c 	, CE1_pop_c   
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_h'	, CE1_nb_h 	, CE1_pop_h  
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_b'	, CE1_nb_b 	, CE1_pop_b  
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_ch'	, CE1_nb_ch , CE1_pop_ch  
		print*, '	------'
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_c'	, CE1_nb_c 	, CE1_pop_c  , 100*((1+CE1_nb_cl/100)*(1+CE1_nb_cd/100)-1)				
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_cl'	, CE1_nb_cl	, CE1_pop_cl
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_cd'	, CE1_nb_cd , CE1_pop_cd
		print*, '	------'
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_h'	, CE1_nb_h 	, CE1_pop_h  , 100*((1+CE1_nb_hl/100)*(1+CE1_nb_hd/100)-1)		
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_hl'	, CE1_nb_hl , CE1_pop_hl
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_hd'	, CE1_nb_hd , CE1_pop_hd 
		print*, '	------'
		print '(A,X,F6.3,X,F6.3,X,F6.3)', '	CE1_b'	, CE1_nb_b 	, CE1_pop_b  , 100*((1+CE1_nb_bl/100)*(1+CE1_nb_bd/100)-1)		
		! print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_bl'	, CE1_nb_bl , CE1_pop_bl
		! print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_bd'	, CE1_nb_bd , CE1_pop_bd 
		print*, '	------'
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_ch'	, CE1_nb_ch  , CE1_pop_ch  , 100*((1+CE1_nb_chl/100)*(1+CE1_nb_chd/100)-1)		
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_chl', CE1_nb_chl , CE1_pop_chl
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_chd', CE1_nb_chd , CE1_pop_chd 
		print*, ' '
		print*, '-------------------------------'
		print*, ' '
		print*, '	CE2 '	, 'NB '		, 'Pop '	 , 'Test ' 	
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2'	, CE2_nb 	, CE2_pop 	 , 100*((1+CE2_nb_c/100)*(1+CE2_nb_h/100)*(1+CE2_nb_b/100)-1)
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_c'	, CE2_nb_c 	, CE2_pop_c  , 100*((1+CE2_nb_cl/100)*(1+CE2_nb_cd/100)-1)		
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_h'	, CE2_nb_h 	, CE2_pop_h  , 100*((1+CE2_nb_hl/100)*(1+CE2_nb_hd/100)-1)
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_b'	, CE2_nb_b 	, CE2_pop_b  , 100*((1+CE2_nb_bl/100)*(1+CE2_nb_bd/100)-1)
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_ch'	, CE2_nb_ch , CE2_pop_ch , 100*((1+CE2_nb_chl/100)*(1+CE2_nb_chd/100)-1)
		print*, '	------'
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_c'	, CE2_nb_c 	, CE2_pop_c  , 100*((1+CE2_nb_cl/100)*(1+CE2_nb_cd/100)-1)		
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_cl'	, CE2_nb_cl	, CE2_pop_cl
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_cd'	, CE2_nb_cd , CE2_pop_cd
		print*, '	------'
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_h'	, CE2_nb_h 	, CE2_pop_h  , 100*((1+CE2_nb_hl/100)*(1+CE2_nb_hd/100)-1)
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_hl'	, CE2_nb_hl , CE2_pop_hl
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_hd'	, CE2_nb_hd , CE2_pop_hd 
		print*, '	------'
		print '(A,X,F6.3,X,F6.3,X,F6.3)', '	CE2_b'	, CE2_nb_b 	, CE2_pop_b  , 100*((1+CE2_nb_bl/100)*(1+CE2_nb_bd/100)-1)		
		! print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_bl'	, CE2_nb_bl , CE2_pop_bl
		! print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_bd'	, CE2_nb_bd , CE2_pop_bd 	
		print*, '	------'
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_ch'	, CE2_nb_ch  , CE2_pop_ch  , 100*((1+CE2_nb_chl/100)*(1+CE2_nb_chd/100)-1)
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_chl', CE2_nb_chl , CE2_pop_chl
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_chd', CE2_nb_chd , CE2_pop_chd 
		print*, ' '
		print*, '-----------------------------------------------------------------------------------'
		print*, ' '

!========================================================================================
!========================================================================================

	! Auxiliary Value Functions
		! Note: we are keeping Aprime at bench values, which affects bequest and expected values trhough 
		!		changes in states. An alternative is to use Aprime_exp and adjust net out the value of 
		!		bequests by subtracting Bq_Value_exp and adding Bq_value_bench from the auxiliary values
		! Change Consumption 
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons_bench,Hours_exp,Aprime_exp,ValueFunction,Bq_Value 	 )
			! Adjust Value to keep bequest value at benchmark level
			ValueFunction = ValueFunction - Bq_Value + Bq_Value_bench
		! Change Consumption and Leisure
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons_exp,Hours_exp  ,Aprime_exp,Value_aux    ,Bq_Value_aux)
			! Adjust Value to keep bequest value at benchmark level
			Value_aux = Value_aux - Bq_Value_aux + Bq_Value_bench

	! CE1 Measures 
		CE1_mat =((ValueFunction_exp  -Bq_Value_bench)/&
				& (ValueFunction_Bench-Bq_Value_bench) )**( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP

		! Decomposition: Leisure
			! Total
			CE1_h_mat  =((ValueFunction      -Bq_Value_bench)/&
					   & (ValueFunction_Bench-Bq_Value_bench) )**( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP
			! Level
			CE1_pop_hl = 100.0_dp*( ((1.0_dp-H_exp   )/(1.0_dp-H_bench   ))**((1.0_dp-gamma)/gamma) -1.0_dp )
			CE1_nb_hl  = 100.0_dp*( ((1.0_dp-H_exp   )/(1.0_dp-H_bench   ))**((1.0_dp-gamma)/gamma) -1.0_dp ) ! 100.0_dp*( ((1.0_dp-H_nb_exp)/(1.0_dp-H_nb_bench))**((1.0_dp-gamma)/gamma) -1.0_dp )


		! Decomposition: Consumption
			! Total
			CE1_c_mat  =((Value_aux      -Bq_Value_bench)/&
					   & (ValueFunction-Bq_Value_bench) )**( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP
			! Level
			CE1_pop_cl = 100.0_dp*( C_exp/C_bench-1.0_dp )
			CE1_nb_cl  = 100.0_dp*( C_exp/C_bench-1.0_dp )! 100.0_dp*( C_nb_exp/C_nb_bench-1.0_dp )

		! Decomposition: Bequests
			! Total
			CE1_b_mat  =((ValueFunction_exp-Bq_Value_bench)/&
					   & (Value_aux        -Bq_Value_bench) )**( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP
			! Level
			CE1_pop_bl = 100.0_dp*( BQ_exp/BQ_bench-1.0_dp )
			CE1_nb_bl  = 100.0_dp*( BQ_exp/BQ_bench-1.0_dp ) ! 100.0_dp*( BQ_nb_exp/BQ_nb_bench-1.0_dp )

		! Decomposition: Consumption and Leisure
			CE1_ch_mat = ((Value_aux    -Bq_Value_bench)/&
					   & (ValueFunction_Bench-Bq_Value_bench) )**( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP

			! Level
			CE1_pop_chl = 100.0_dp*( C_exp/C_bench*((1.0_dp-H_exp)/(1.0_dp-H_bench))**((1.0_dp-gamma)/gamma) -1.0_dp )
			CE1_nb_chl  = 100.0_dp*( C_exp/C_bench*((1.0_dp-H_exp)/(1.0_dp-H_bench))**((1.0_dp-gamma)/gamma) -1.0_dp ) ! 100.0_dp*( C_nb_exp/C_nb_bench*((1.0_dp-H_nb_exp)/(1.0_dp-H_nb_bench))**((1.0_dp-gamma)/gamma) -1.0_dp )

		! Aggregating CE1 for population or newborns 
			CE1_nb		= 100.0_dp*sum( CE1_mat(1,:,:,:,:,:)  *DBN_bench(1,:,:,:,:,:) )/size_nb
			CE1_nb_c	= 100.0_dp*sum( CE1_c_mat(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:) )/size_nb
			CE1_nb_cd	= 100.0_dp*( (CE1_nb_c/100.0_dp+1.0_dp)/(CE1_nb_cl/100.0_dp+1.0_dp) - 1.0_dp )
			CE1_nb_h	= 100.0_dp*sum( CE1_h_mat(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:) )/size_nb
			CE1_nb_hd	= 100.0_dp*( (CE1_nb_h/100.0_dp+1.0_dp)/(CE1_nb_hl/100.0_dp+1.0_dp) - 1.0_dp )
			CE1_nb_b	= 100.0_dp*sum( CE1_b_mat(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:) )/size_nb
			CE1_nb_bd	= 100.0_dp*( (CE1_nb_b/100.0_dp+1.0_dp)/(CE1_nb_bl/100.0_dp+1.0_dp) - 1.0_dp )
			CE1_nb_ch	= 100.0_dp*sum( CE1_ch_mat(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:) )/size_nb
			CE1_nb_chd	= 100.0_dp*( (CE1_nb_ch/100.0_dp+1.0_dp)/(CE1_nb_chl/100.0_dp+1.0_dp) - 1.0_dp )

			CE1_pop		= 100.0_dp*sum( CE1_mat  *DBN_bench )
			CE1_pop_c	= 100.0_dp*sum( CE1_c_mat*DBN_bench )
			CE1_pop_cd	= 100.0_dp*( (CE1_pop_c/100.0_dp+1.0_dp)/(CE1_pop_cl/100.0_dp+1.0_dp) - 1.0_dp )
			CE1_pop_h	= 100.0_dp*sum( CE1_h_mat*DBN_bench )
			CE1_pop_hd	= 100.0_dp*( (CE1_pop_h/100.0_dp+1.0_dp)/(CE1_pop_hl/100.0_dp+1.0_dp) - 1.0_dp )
			CE1_pop_b	= 100.0_dp*sum( CE1_b_mat*DBN_bench )
			CE1_pop_bd	= 100.0_dp*( (CE1_pop_b/100.0_dp+1.0_dp)/(CE1_pop_bl/100.0_dp+1.0_dp) - 1.0_dp )
			CE1_pop_ch	= 100.0_dp*sum( CE1_ch_mat*DBN_bench )
			CE1_pop_chd	= 100.0_dp*( (CE1_pop_ch/100.0_dp+1.0_dp)/(CE1_pop_chl/100.0_dp+1.0_dp) - 1.0_dp )



	! CE2 Measures 
		CE2_nb  =  100.0_dp*(( (sum(ValueFunction_exp(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:)) - &
	    					&  sum(Bq_Value_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))) / &
	                        &  sum((ValueFunction_Bench(1,:,:,:,:,:)-Bq_Value_bench(1,:,:,:,:,:))*DBN_bench(1,:,:,:,:,:))  ) &
	                        &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

		CE2_pop = 100.0_dp*(( (sum(ValueFunction_exp*DBN_exp)-sum(Bq_Value_Bench*DBN_bench)) / &
	    					& sum((ValueFunction_Bench-Bq_Value_Bench)*DBN_bench)  ) &
	                        &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

		! Decomposition: Leisure
			! Total (note that auxiliary value functions are weighted with benchmark distribution)
			CE2_nb_h   =  100.0_dp*(( (sum(ValueFunction(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:)) - &
	    						&  sum(Bq_Value_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))) / &
	                        	&  sum((ValueFunction_Bench(1,:,:,:,:,:)-Bq_Value_bench(1,:,:,:,:,:))*DBN_bench(1,:,:,:,:,:))  ) &
	                        	&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

			CE2_pop_h  = 100.0_dp*(( (sum(ValueFunction*DBN_exp)-sum(Bq_Value_Bench*DBN_bench)) / &
	    						& sum((ValueFunction_Bench-Bq_Value_Bench)*DBN_bench)  ) &
	                        	&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
			! Level
			CE2_nb_hl  = 100.0_dp*( ((1.0_dp-H_exp   )/(1.0_dp-H_bench   ))**((1.0_dp-gamma)/gamma) -1.0_dp ) ! 100.0_dp*( ((1.0_dp-H_nb_exp)/(1.0_dp-H_nb_bench))**((1.0_dp-gamma)/gamma) -1.0_dp )
			CE2_pop_hl = 100.0_dp*( ((1.0_dp-H_exp   )/(1.0_dp-H_bench   ))**((1.0_dp-gamma)/gamma) -1.0_dp )
			! Distribution
			CE2_nb_hd  = 100.0_dp*( (CE2_nb_h/100.0_dp+1.0_dp)/(CE2_nb_hl/100.0_dp+1.0_dp) - 1.0_dp )
			CE2_pop_hd = 100.0_dp*( (CE2_pop_h/100.0_dp+1.0_dp)/(CE2_pop_hl/100.0_dp+1.0_dp) - 1.0_dp )

		! Decomposition: Consumption
			! Total (note that auxiliary value functions are weighted with benchmark distribution)
			CE2_nb_c   =  100.0_dp*(( (sum(Value_aux(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:)) - &
	    						&  sum(Bq_Value_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))) / &
	                        	&  (sum(ValueFunction(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:)) - &
	    						&  sum(Bq_Value_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:)))  ) &
	                        	&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

			CE2_pop_c  = 100.0_dp*(( (sum(Value_aux*DBN_exp)-sum(Bq_Value_Bench*DBN_bench)) / &
	    						& (sum(ValueFunction*DBN_exp)-sum(Bq_Value_Bench*DBN_bench))  ) &
	                        	&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

			! Level
			CE2_nb_cl  = 100.0_dp*( C_exp/C_bench-1.0_dp ) ! 100.0_dp*( C_nb_exp/C_nb_bench-1.0_dp )
			CE2_pop_cl = 100.0_dp*( C_exp/C_bench-1.0_dp )
			! Distribution
			CE2_nb_cd  = 100.0_dp*( (CE2_nb_c/100.0_dp+1.0_dp)/(CE2_nb_cl/100.0_dp+1.0_dp) - 1.0_dp )
			CE2_pop_cd = 100.0_dp*( (CE2_pop_c/100.0_dp+1.0_dp)/(CE2_pop_cl/100.0_dp+1.0_dp) - 1.0_dp )



		! Decomposition: Bequests
			! Total (note that auxiliary value functions are weighted with experiment distribution)
			CE2_nb_b   =  100.0_dp*(( (sum(ValueFunction_exp(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:)) - &
	    						&  sum(Bq_Value_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))) / &
	                        	&  (sum(Value_aux(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:)) - &
	    						&  sum(Bq_Value_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:)))  ) &
	                        	&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

			CE2_pop_b  = 100.0_dp*(( (sum(ValueFunction_exp*DBN_exp)-sum(Bq_Value_Bench*DBN_bench)) / &
	    						& (sum(Value_aux*DBN_exp)-sum(Bq_Value_Bench*DBN_bench))  ) &
	                        	&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

			! Level
			CE2_nb_bl  = 100.0_dp*( BQ_exp/BQ_bench-1.0_dp ) ! 100.0_dp*( BQ_nb_exp/BQ_nb_bench-1.0_dp )
			CE2_pop_bl = 100.0_dp*( BQ_exp/BQ_bench-1.0_dp )
			! Distribution
			CE2_nb_bd  = 100.0_dp*( (CE2_nb_b/100.0_dp+1.0_dp)/(CE2_nb_bl/100.0_dp+1.0_dp) - 1.0_dp )
			CE2_pop_bd = 100.0_dp*( (CE2_pop_b/100.0_dp+1.0_dp)/(CE2_pop_bl/100.0_dp+1.0_dp) - 1.0_dp )

		! Decomposition: Consumption and Leisure
			! Total (note that auxiliary value functions are weighted with benchmark distribution)
			CE2_nb_ch   =  100.0_dp*(( (sum(Value_aux(1,:,:,:,:,:)*DBN_exp(1,:,:,:,:,:)) - &
	    						&  sum(Bq_Value_bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))) / &
	                        	&  sum((ValueFunction_Bench(1,:,:,:,:,:)-Bq_Value_bench(1,:,:,:,:,:))*DBN_bench(1,:,:,:,:,:))  ) &
	                        	&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

			CE2_pop_ch  = 100.0_dp*(( (sum(Value_aux*DBN_exp)-sum(Bq_Value_Bench*DBN_bench)) / &
	    						& sum((ValueFunction_Bench-Bq_Value_Bench)*DBN_bench)  ) &
	                        	&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
			! Level
			CE2_pop_chl = 100.0_dp*( C_exp/C_bench*((1.0_dp-H_exp)/(1.0_dp-H_bench))**((1.0_dp-gamma)/gamma) -1.0_dp )
			CE2_nb_chl  = 100.0_dp*( C_exp/C_bench*((1.0_dp-H_exp)/(1.0_dp-H_bench))**((1.0_dp-gamma)/gamma) -1.0_dp ) ! 100.0_dp*( C_nb_exp/C_nb_bench*((1.0_dp-H_nb_exp)/(1.0_dp-H_nb_bench))**((1.0_dp-gamma)/gamma) -1.0_dp )
			! Distribution
			CE2_nb_chd  = 100.0_dp*( (CE2_nb_ch/100.0_dp+1.0_dp)/(CE2_nb_chl/100.0_dp+1.0_dp) - 1.0_dp )
			CE2_pop_chd = 100.0_dp*( (CE2_pop_ch/100.0_dp+1.0_dp)/(CE2_pop_chl/100.0_dp+1.0_dp) - 1.0_dp )




	! Tables 
		OPEN  (UNIT=1,  FILE=trim(Result_Folder)//'CE_Decomposition_alt_order.txt'  , STATUS='replace')
		WRITE (UNIT=1,  FMT=*) ' '
		WRITE (UNIT=1,  FMT=*) 'Decomposition: Consumption Equivalent Welfare'
		WRITE (UNIT=1,  FMT=*) ' '
		WRITE (UNIT=1,  FMT=*) 'CE1 '	, 'NB '		, 'Pop '	 , 'Test ' 	
		WRITE (UNIT=1,  FMT=*) 'CE1'	, CE1_nb 	, CE1_pop    , 100*((1+CE1_nb_c/100)*(1+CE1_nb_h/100)*(1+CE1_nb_b/100)-1)
		WRITE (UNIT=1,  FMT=*) 'CE1_c'	, CE1_nb_c 	, CE1_pop_c   
		WRITE (UNIT=1,  FMT=*) 'CE1_h'	, CE1_nb_h 	, CE1_pop_h  
		WRITE (UNIT=1,  FMT=*) 'CE1_b'	, CE1_nb_b 	, CE1_pop_b 
		WRITE (UNIT=1,  FMT=*) 'CE1_ch'	, CE1_nb_ch , CE1_pop_ch
		WRITE (UNIT=1,  FMT=*) '------'
		WRITE (UNIT=1,  FMT=*) 'CE1_c'	, CE1_nb_c 	, CE1_pop_c  , 100*((1+CE1_nb_cl/100)*(1+CE1_nb_cd/100)-1)				
		WRITE (UNIT=1,  FMT=*) 'CE1_cl'	, CE1_nb_cl	, CE1_pop_cl
		WRITE (UNIT=1,  FMT=*) 'CE1_cd'	, CE1_nb_cd , CE1_pop_cd
		WRITE (UNIT=1,  FMT=*) '------'
		WRITE (UNIT=1,  FMT=*) 'CE1_h'	, CE1_nb_h 	, CE1_pop_h  , 100*((1+CE1_nb_hl/100)*(1+CE1_nb_hd/100)-1)		
		WRITE (UNIT=1,  FMT=*) 'CE1_hl'	, CE1_nb_hl , CE1_pop_hl
		WRITE (UNIT=1,  FMT=*) 'CE1_hd'	, CE1_nb_hd , CE1_pop_hd 
		WRITE (UNIT=1,  FMT=*) '------'
		WRITE (UNIT=1,  FMT=*) 'CE1_b'	, CE1_nb_b 	, CE1_pop_b  , 100*((1+CE1_nb_bl/100)*(1+CE1_nb_bd/100)-1)		
		! WRITE (UNIT=1,  FMT=*) 'CE1_bl'	, CE1_nb_bl , CE1_pop_bl
		! WRITE (UNIT=1,  FMT=*) 'CE1_bd'	, CE1_nb_bd , CE1_pop_bd 
		WRITE (UNIT=1,  FMT=*) '------'
		WRITE (UNIT=1,  FMT=*) 'CE1_ch'	, CE1_nb_ch  , CE1_pop_ch  , 100*((1+CE1_nb_chl/100)*(1+CE1_nb_chd/100)-1)		
		WRITE (UNIT=1,  FMT=*) 'CE1_chl', CE1_nb_chl , CE1_pop_chl
		WRITE (UNIT=1,  FMT=*) 'CE1_chd', CE1_nb_chd , CE1_pop_chd
		WRITE (UNIT=1,  FMT=*) ' '
		WRITE (UNIT=1,  FMT=*) '-------------------------------'
		WRITE (UNIT=1,  FMT=*) ' '
		WRITE (UNIT=1,  FMT=*) 'CE2 '	, 'NB '		, 'Pop '	 , 'Test ' 	
		WRITE (UNIT=1,  FMT=*) 'CE2'	, CE2_nb 	, CE2_pop 	 , 100*((1+CE2_nb_c/100)*(1+CE2_nb_h/100)*(1+CE2_nb_b/100)-1)
		WRITE (UNIT=1,  FMT=*) 'CE2_c'	, CE2_nb_c 	, CE2_pop_c  , 100*((1+CE2_nb_cl/100)*(1+CE2_nb_cd/100)-1)		
		WRITE (UNIT=1,  FMT=*) 'CE2_h'	, CE2_nb_h 	, CE2_pop_h  , 100*((1+CE2_nb_hl/100)*(1+CE2_nb_hd/100)-1)
		WRITE (UNIT=1,  FMT=*) 'CE2_b'	, CE2_nb_b 	, CE2_pop_b  , 100*((1+CE2_nb_bl/100)*(1+CE2_nb_bd/100)-1)
		WRITE (UNIT=1,  FMT=*) 'CE2_ch'	, CE2_nb_ch , CE2_pop_ch
		WRITE (UNIT=1,  FMT=*) '------'
		WRITE (UNIT=1,  FMT=*) 'CE2_c'	, CE2_nb_c 	, CE2_pop_c  , 100*((1+CE2_nb_cl/100)*(1+CE2_nb_cd/100)-1)		
		WRITE (UNIT=1,  FMT=*) 'CE2_cl'	, CE2_nb_cl	, CE2_pop_cl
		WRITE (UNIT=1,  FMT=*) 'CE2_cd'	, CE2_nb_cd , CE2_pop_cd
		WRITE (UNIT=1,  FMT=*) '------'
		WRITE (UNIT=1,  FMT=*) 'CE2_h'	, CE2_nb_h 	, CE2_pop_h  , 100*((1+CE2_nb_hl/100)*(1+CE2_nb_hd/100)-1)
		WRITE (UNIT=1,  FMT=*) 'CE2_hl'	, CE2_nb_hl , CE2_pop_hl
		WRITE (UNIT=1,  FMT=*) 'CE2_hd'	, CE2_nb_hd , CE2_pop_hd 
		WRITE (UNIT=1,  FMT=*) '------'
		WRITE (UNIT=1,  FMT=*) 'CE2_b'	, CE2_nb_b 	, CE2_pop_b  , 100*((1+CE2_nb_bl/100)*(1+CE2_nb_bd/100)-1)		
		! WRITE (UNIT=1,  FMT=*) 'CE2_bl'	, CE2_nb_bl , CE2_pop_bl
		! WRITE (UNIT=1,  FMT=*) 'CE2_bd'	, CE2_nb_bd , CE2_pop_bd 
		WRITE (UNIT=1,  FMT=*) '------'
		WRITE (UNIT=1,  FMT=*) 'CE2_ch'	, CE2_nb_ch  , CE2_pop_ch  , 100*((1+CE2_nb_chl/100)*(1+CE2_nb_chd/100)-1)		
		WRITE (UNIT=1,  FMT=*) 'CE2_chl', CE2_nb_chl , CE2_pop_chl
		WRITE (UNIT=1,  FMT=*) 'CE2_chd', CE2_nb_chd , CE2_pop_chd
		WRITE (UNIT=1,  FMT=*) ' '
		CLOSE (unit=1)

		print*, ' '
		print*, '-----------------------------------------------------------------------------------'
		print*, ' '
		print*, 'Decomposition: Consumption Equivalent Welfare'
		print*, ' 	Alternative Order'
		print*, ' '
		print*, '	CE1 '	, 'NB '		, 'Pop '	 , 'Test ' 
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1'	, CE1_nb 	, CE1_pop   , 100*((1+CE1_nb_c/100)*(1+CE1_nb_h/100)*(1+CE1_nb_b/100)-1)
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_c'	, CE1_nb_c 	, CE1_pop_c   
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_h'	, CE1_nb_h 	, CE1_pop_h  
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_b'	, CE1_nb_b 	, CE1_pop_b  
		print*, '	------'
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_c'	, CE1_nb_c 	, CE1_pop_c  , 100*((1+CE1_nb_cl/100)*(1+CE1_nb_cd/100)-1)				
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_cl'	, CE1_nb_cl	, CE1_pop_cl
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_cd'	, CE1_nb_cd , CE1_pop_cd
		print*, '	------'
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_h'	, CE1_nb_h 	, CE1_pop_h  , 100*((1+CE1_nb_hl/100)*(1+CE1_nb_hd/100)-1)		
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_hl'	, CE1_nb_hl , CE1_pop_hl
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE1_hd'	, CE1_nb_hd , CE1_pop_hd 
		print*, ' '
		print*, '-------------------------------'
		print*, ' '
		print*, '	CE2 '	, 'NB '		, 'Pop '	 , 'Test ' 	
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2'	, CE2_nb 	, CE2_pop 	 , 100*((1+CE2_nb_c/100)*(1+CE2_nb_h/100)*(1+CE2_nb_b/100)-1)
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_c'	, CE2_nb_c 	, CE2_pop_c  , 100*((1+CE2_nb_cl/100)*(1+CE2_nb_cd/100)-1)		
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_h'	, CE2_nb_h 	, CE2_pop_h  , 100*((1+CE2_nb_hl/100)*(1+CE2_nb_hd/100)-1)
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_b'	, CE2_nb_b 	, CE2_pop_b  , 100*((1+CE2_nb_bl/100)*(1+CE2_nb_bd/100)-1)
		print*, '	------'
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_c'	, CE2_nb_c 	, CE2_pop_c  , 100*((1+CE2_nb_cl/100)*(1+CE2_nb_cd/100)-1)		
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_cl'	, CE2_nb_cl	, CE2_pop_cl
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_cd'	, CE2_nb_cd , CE2_pop_cd
		print*, '	------'
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_h'	, CE2_nb_h 	, CE2_pop_h  , 100*((1+CE2_nb_hl/100)*(1+CE2_nb_hd/100)-1)
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_hl'	, CE2_nb_hl , CE2_pop_hl
		print '(A,X,F6.2,X,F6.2,X,F6.2)', '	CE2_hd'	, CE2_nb_hd , CE2_pop_hd 
		print*, ' '
		print*, '-----------------------------------------------------------------------------------'
		print*, ' '




END SUBROUTINE COMPUTE_WELFARE_DECOMPOSITION


!========================================================================================
!========================================================================================
!========================================================================================


Function Draft_Table(Table_az,DBN_z,Cum_flag)
	implicit none
	real(dp), dimension(draft_age_category,nz), intent(in)   :: Table_az
	real(dp), dimension(nz), intent(in) :: DBN_z
	real(dp), dimension(nz) 			:: CDF_z
	logical, intent(in)					:: Cum_flag
	integer  						    :: zi
	real(dp), dimension(draft_age_category,draft_z_category) :: Draft_Table


	! CDF
	do zi=1,nz 
		CDF_Z(zi) = sum(DBN_Z(1:zi))
	enddo 

	if (Cum_flag) then
	! Cumulative value of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		Draft_Table(:,1) = Table_az(:,1) + Table_az(:,2) + Table_az(:,3) & 
								&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*Table_az(:,4)
		Draft_Table(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*Table_az(:,4)+((0.80_dp-CDF_Z(4))/DBN_Z(5))*Table_az(:,5)
		Draft_Table(:,3) = (0.10_dp/DBN_Z(5))*Table_az(:,5)
		Draft_Table(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*Table_az(:,5) + &
							& ((0.99_dp-CDF_Z(5))/DBN_Z(6))*Table_az(:,6) 
		Draft_Table(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*Table_az(:,6) + & 
							& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*Table_az(:,7) 
		Draft_Table(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Table_az(:,7) + Table_az(:,8) + Table_az(:,9) 
		Draft_Table(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Table_az(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*Table_az(:,8)
		Draft_Table(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*Table_az(:,8) + Table_az(:,9)

	else
	! Average value of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		Draft_Table(:,1) = ( DBN_Z(1)*Table_az(:,1) + DBN_Z(2)*Table_az(:,2) + & 
								& DBN_Z(3)*Table_az(:,3) + (0.40_dp-CDF_Z(3))*Table_az(:,4) )/0.40_dp
		Draft_Table(:,2) = ( (CDF_Z(4)-0.40_dp)*Table_az(:,4) + (0.80_dp-CDF_Z(4))*Table_az(:,5) )/0.40_dp
		Draft_Table(:,3) = Table_az(:,5)
		Draft_Table(:,4) = ( (CDF_Z(5)-0.90_dp)*Table_az(:,5) + (0.99_dp-CDF_Z(5))*Table_az(:,6) )/0.09_dp
		Draft_Table(:,5) = ( (CDF_Z(6)-0.99_dp)*Table_az(:,6) + (0.999_dp-CDF_Z(6))*Table_az(:,7) )/0.009_dp
		Draft_Table(:,6) = ( (CDF_Z(7)-0.999_dp)*Table_az(:,7) + &
								&  DBN_Z(8)*Table_az(:,8) + DBN_Z(9)*Table_az(:,9) )/0.001_dp
		Draft_Table(:,7) = ( (CDF_Z(7)-0.999_dp)*Table_az(:,7) + (0.9999_dp-CDF_Z(7))*Table_az(:,8) )/0.0009_dp
		Draft_Table(:,8) = ( (CDF_Z(8)-0.9999_dp)*Table_az(:,8) + DBN_Z(9)*Table_az(:,9) )/0.0001_dp

	endif 

END Function Draft_Table


!========================================================================================
!========================================================================================
!========================================================================================


Function Draft_Table_X(Table_axz,DBN_xz,Cum_flag)
	implicit none
	real(dp), dimension(draft_age_category,nz,nx), intent(in)   :: Table_axz
	real(dp), dimension(nz,nx), intent(in) 	:: DBN_xz
	real(dp)				   				:: CDF_xz, cdf_xz_low
	logical, intent(in)						:: Cum_flag
	real(dp), dimension(draft_age_category,draft_z_category) :: Draft_Table_X

	! Groups based on cross sectional productivity, not age dependent 
	    ! % 0%-40%     of Current Productivity
	    !     % (z1,x3), (z2,x3), (z3,x3), (z4,x3), (z5,x3)
	    ! % 40%-80%    of Current Productivity 
	    !     % (z5,x3), (z6,x3), (z7,x3), (z8,x3), (z9,x3), (z1,x1), (z2,x1), (z3,x1), (z4,x1)
	    ! % 80%-90%    of Current Productivity
	    !     % (z4,x1), (z5,x2)
	    ! % 90%-99%    of Current Productivity 
	    !     % (z5,x2), (z6,x2), (z7,x2), (z8,x2), (z9,x2), (z5,x1), (z6,x1)
	    ! % 99%-99.9%  of Current Productivity 
	    !     % (z6,x1), (z7,x1)
	    ! % 99.9%-100% of Current Productivity 
	    !     % (z7,x1), (z8,x1), (z9,x1)

	if (Cum_flag) then
	! Cumulative value of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
    	cdf_xz = sum(DBN_XZ(1:4,3))
		Draft_Table_X(:,1)   = Table_axz(:,1,3) + Table_axz(:,2,3) + Table_axz(:,3,3) + Table_axz(:,4,3) + &
							& (0.40_dp-cdf_xz)*Table_axz(:,5,3)/DBN_XZ(5,3) 
			cdf_xz = cdf_xz + DBN_XZ(5,3)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(6:,3)) + sum(DBN_XZ(1:3,1))
		Draft_Table_X(:,2)   = (cdf_xz_low-0.40_dp)*Table_axz(:,5,3)/DBN_XZ(5,3) + &
						& Table_axz(:,6,3) + Table_axz(:,7,3) + Table_axz(:,8,3) + Table_axz(:,9,3) + &
						& Table_axz(:,1,1) + Table_axz(:,2,1) + Table_axz(:,3,1) + &
						& (0.80_dp-cdf_xz)*Table_axz(:,4,1)/DBN_XZ(4,1)
			cdf_xz = cdf_xz + DBN_XZ(4,1)
		Draft_Table_X(:,3)   = (cdf_xz-0.80_dp)*Table_axz(:,4,1)/DBN_XZ(4,1) + &
						& (0.90_dp-cdf_xz)*Table_axz(:,5,2)/DBN_XZ(5,2)
			cdf_xz = cdf_xz + DBN_XZ(5,2)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(6:,2)) + DBN_XZ(5,1)
		Draft_Table_X(:,4)   = (cdf_xz_low-0.90_dp)*Table_axz(:,5,2)/DBN_XZ(5,2) + & 
						& Table_axz(:,6,2) + Table_axz(:,7,2) + Table_axz(:,8,2) + Table_axz(:,9,2) + &
						& Table_axz(:,5,1) + (0.99_dp-cdf_xz)*Table_axz(:,6,1)/DBN_XZ(6,1)
			cdf_xz = cdf_xz + DBN_XZ(6,1)
		Draft_Table_X(:,5)   = (cdf_xz-0.99_dp)*Table_axz(:,6,1)/DBN_XZ(6,1) + &
						& (0.999_dp-cdf_xz)*Table_axz(:,7,1)/DBN_XZ(7,1)
			cdf_xz = cdf_xz + DBN_XZ(7,1)
		Draft_Table_X(:,6)   = (CDF_xz-0.999_dp)*Table_axz(:,7,1)/DBN_XZ(7,1) + Table_axz(:,8,1) + Table_axz(:,9,1)
		Draft_Table_X(:,7)   = 0
		Draft_Table_X(:,8)   = 0

	else
	! Average value of groups adjusting by xz group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
    	cdf_xz = sum(DBN_XZ(1:4,3))
		Draft_Table_X(:,1)   = ( DBN_XZ(1,3)*Table_axz(:,1,3) + DBN_XZ(2,3)*Table_axz(:,2,3) + &
						&  DBN_XZ(3,3)*Table_axz(:,3,3) + DBN_XZ(4,3)*Table_axz(:,4,3) + &
						& (0.40_dp-cdf_xz)*Table_axz(:,5,3) )/0.40_dp 
			cdf_xz = cdf_xz + DBN_XZ(5,3)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(6:,3)) + sum(DBN_XZ(1:3,1))
		Draft_Table_X(:,2)   = ( (cdf_xz_low-0.40_dp)*Table_axz(:,5,3) + &
						& DBN_XZ(6,3)*Table_axz(:,6,3) + DBN_XZ(7,3)*Table_axz(:,7,3) + &
						& DBN_XZ(8,3)*Table_axz(:,8,3) + DBN_XZ(9,3)*Table_axz(:,9,3) + &
						& DBN_XZ(1,1)*Table_axz(:,1,1) + DBN_XZ(2,1)*Table_axz(:,2,1) + &
						& DBN_XZ(3,1)*Table_axz(:,3,1) + (0.80_dp-cdf_xz)*Table_axz(:,4,1) )/0.40_dp
			cdf_xz = cdf_xz + DBN_XZ(4,1)
		Draft_Table_X(:,3)   = ( (cdf_xz-0.80_dp)*Table_axz(:,4,1) + (0.90_dp-cdf_xz)*Table_axz(:,5,2) )/0.10_dp
			cdf_xz = cdf_xz + DBN_XZ(5,2)
			cdf_xz_low = cdf_xz 
			cdf_xz = cdf_xz + sum(DBN_XZ(6:,2)) + DBN_XZ(5,1)
		Draft_Table_X(:,4)   = ( (cdf_xz_low-0.90_dp)*Table_axz(:,5,2) + & 
						& DBN_XZ(6,2)*Table_axz(:,6,2) + DBN_XZ(7,2)*Table_axz(:,7,2) + &
						& DBN_XZ(8,2)*Table_axz(:,8,2) + DBN_XZ(9,2)*Table_axz(:,9,2) + &
						& DBN_XZ(5,1)*Table_axz(:,5,1) + (0.99_dp-cdf_xz)*Table_axz(:,6,1) )/0.09_dp
			cdf_xz = cdf_xz + DBN_XZ(6,1)
		Draft_Table_X(:,5)   = ( (cdf_xz-0.99_dp)*Table_axz(:,6,1) + (0.999_dp-cdf_xz)*Table_axz(:,7,1) )/0.009_dp
			cdf_xz = cdf_xz + DBN_XZ(7,1)
		Draft_Table_X(:,6)   = ( (CDF_xz-0.999_dp)*Table_axz(:,7,1)+DBN_XZ(8,1)*Table_axz(:,8,1)+DBN_XZ(9,1)*Table_axz(:,9,1) )/0.001_dp
		Draft_Table_X(:,7)   = 0
		Draft_Table_X(:,8)   = 0

			! Adjustment for first age group
			    ! % 40%-80%    of Current Productivity 
			    !     % (z1,x1), (z2,x1), (z3,x1), (z4,x1)
		    	Draft_Table_X(1,2)   = ( DBN_XZ(1,1)*Table_axz(1,1,1) + DBN_XZ(2,1)*Table_axz(1,2,1) + &
							& DBN_XZ(3,1)*Table_axz(1,3,1) + DBN_XZ(4,1)*Table_axz(1,4,1) )/sum(DBN_XZ(1:4,1))
			    ! % 80%-90%    of Current Productivity
			    !     % (z4,x1)
			    Draft_Table_X(1,3)   = Table_axz(1,4,1)
			    ! % 90%-99%    of Current Productivity 
			    !     % (z5,x1), (z6,x1)
			    Draft_Table_X(1,4)   = ( DBN_XZ(5,1)*Table_axz(1,5,1) + DBN_XZ(6,1)*Table_axz(1,6,1) )/sum(DBN_XZ(5:6,1))


	endif 

END Function Draft_Table_X


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE COMPUTE_MOMENTS
	IMPLICIT NONE
	real(DP), dimension(size(DBN1)) :: DBN_vec, Firm_Wealth_vec, CDF_Firm_Wealth, BQ_vec, DBN_bq_vec, CDF_bq
	real(DP) :: FW_top_x(4),  prctile_FW(4), prctile_bq(7), a, b, c, CCDF_c
	real(DP) :: DBN_bq(MaxAge,na,nz,nlambda,ne)
	integer  :: i, prctile


	! COMPUTE AVERAGE HOURS FOR AGES 25-60 (5-40 IN THE MODEL) INCLUDING NON-WORKERS
		! COMPUTE VARIANCE OF LOG EARNINGS FOR 25-60 FOR THOSE WHO WORK MORE THAN 260 HOURS
		! WHICH CORRESPOND TO 0.055 IN THE MODEL
		pop_25_60        	   = 0.0_DP
		tothours_25_60         = 0.0_DP
		pop_pos_earn_25_60     = 0.0_DP
		tot_log_earnings_25_60 = 0.0_DP 
		Var_Log_Earnings_25_60 = 0.0_DP
		do xi=1,nx
		DO age=5,40
		DO ai=1,na
		DO zi=1,nz
		DO lambdai=1,nlambda
		DO ei=1,ne
			tothours_25_60 = tothours_25_60 + DBN1(age, ai, zi, lambdai, ei, xi)  * Hours(age, ai, zi, lambdai,ei,xi)
			pop_25_60      = pop_25_60 +  DBN1(age, ai, zi, lambdai, ei,xi)
			IF (Hours(age, ai, zi, lambdai, ei,xi) .ge. 0.055) THEN
			tot_log_earnings_25_60 = tot_log_earnings_25_60 + DBN1(age, ai, zi, lambdai, ei,xi)  &
			                 		& *  log( Y_h(Hours(age, ai, zi, lambdai, ei,xi),age,lambdai,ei,wage) )
			pop_pos_earn_25_60     = pop_pos_earn_25_60 +  DBN1(age, ai, zi, lambdai, ei,xi)
			ENDIF
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		meanhours_25_60         = tothours_25_60 / pop_25_60
		mean_log_earnings_25_60 = tot_log_earnings_25_60 / pop_pos_earn_25_60

		DO xi=1,nx
		DO age=5,40
		DO ai=1,na
		DO zi=1,nz
		DO lambdai=1,nlambda
		DO ei=1,ne
			IF (Hours(age, ai, zi, lambdai, ei,xi) .ge. 0.055) THEN
			    Var_Log_Earnings_25_60 =  Var_Log_Earnings_25_60 + DBN1(age, ai, zi, lambdai, ei,xi)  &
			                 			& * ( log( Y_h(Hours(age, ai, zi, lambdai, ei,xi),age,lambdai,ei,wage) ) &
			                 			& -   mean_log_earnings_25_60 ) ** 2.0_DP
			ENDIF
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		Var_Log_Earnings_25_60 = Var_Log_Earnings_25_60 / pop_pos_earn_25_60
		Std_Log_Earnings_25_60 = Var_Log_Earnings_25_60 ** 0.5_DP


	! Sources of income
		Pr_mat = Profit_Matrix(R,P)

		MeanWealth 	 = 0.0_DP
		MeanReturn 	 = 0.0_DP
		
		DO xi=1,nx
		DO age=1,MaxAge
		DO zi=1,nz
		DO ai=1,na
		DO lambdai=1,nlambda
		DO ei=1, ne
		    MeanWealth   = MeanWealth   + DBN1(age, ai, zi, lambdai, ei, xi)*agrid(ai)
		    MeanReturn           = MeanReturn          + DBN1(age, ai, zi, lambdai, ei, xi) * &
		    							& (R*agrid(ai) + Pr_mat(ai,zi,xi))
		    
		ENDDO
		ENDDO
		ENDDO
		ENDDO    
		ENDDO    
		ENDDO
		Wealth_Output = MeanWealth/YBAR 
		MeanReturn    = MeanReturn/MeanWealth
		



	! Distribution of Assets
		DO ai=1,na
		     pr_a_dbn(ai)          = sum(DBN1(:,ai,:,:,:,:)) 
		     cdf_a_dbn(ai)         = sum( pr_a_dbn(1:ai) )      
		     tot_a_by_grid(ai)     = sum(DBN1(:,ai,:,:,:,:) * agrid(ai) )
		     cdf_tot_a_by_grid(ai) = sum(tot_a_by_grid(1:ai))   
		!     print*, pr_a_dbn(ai), cdf_a_dbn(ai)
		ENDDO
		cdf_a_dbn = cdf_a_dbn + 1.0_DP - cdf_a_dbn(na)

		
		! FIND THE ai THAT CORRESPONDS TO EACH PRCTILE OF WEALTH DBN & WEALTH HELD BY PEOPLE LOWER THAN THAT PRCTILE
		DO prctile=1,100
		    ai=1
		    DO while (cdf_a_dbn(ai) .lt. (REAL(prctile,8)/100.0_DP-0.000000000000001))
		        ai=ai+1
		    ENDDO
		    prctile_ai_ind(prctile) = ai
		    prctile_ai(prctile)     = agrid(ai)
		    ! print*,prctile, REAL(prctile,8)/100.0_DP,  ai
		    IF (ai .gt. 1) THEN
		        cdf_tot_a_by_prctile(prctile)  = cdf_tot_a_by_grid(ai-1) + (REAL(prctile,8)/100.0_DP - cdf_a_dbn(ai-1))*agrid(ai) 
		    else
		        cdf_tot_a_by_prctile(prctile)  = (REAL(prctile,8)/100.0_DP )*agrid(ai)     
		    ENDIF
		ENDDO
		print*,''
		prct1_wealth  = 1.0_DP-cdf_tot_a_by_prctile(99)/cdf_tot_a_by_prctile(100)
		prct10_wealth = 1.0_DP-cdf_tot_a_by_prctile(90)/cdf_tot_a_by_prctile(100)
		prct20_wealth = 1.0_DP-cdf_tot_a_by_prctile(80)/cdf_tot_a_by_prctile(100)
		prct40_wealth = 1.0_DP-cdf_tot_a_by_prctile(60)/cdf_tot_a_by_prctile(100)

	! Distribution of firm wealth
		do ai=1,na
			Firm_Wealth(:,ai,:,:,:,:) = V_Pr(:,ai,:,:,:,:) + (1.0_dp+R)*agrid(ai)
		enddo
		Mean_Firm_Wealth = sum(Firm_Wealth*DBN1)

		DBN_vec         = reshape(DBN1       ,(/size(DBN1)/))
		Firm_Wealth_vec = reshape(Firm_Wealth,(/size(DBN1)/))


		prctile_FW = (/0.40_DP, 0.20_dp, 0.10_dp, 0.01_dp/)
		a = minval(Firm_Wealth_vec)
		b = maxval(Firm_Wealth_vec) 
		c = a
		do i=1,size(prctile_FW)
			a = c
			b = maxval(Firm_Wealth_vec)
			c = (a+b)/2.0_dp
			CCDF_c = sum(DBN_vec,Firm_Wealth_vec>=c)
			print*, ' '
			!print*, 'Percentile', prctile_FW(i)
			do while ((abs(CCDF_c-prctile_FW(i))>0.0001_dp).and.(b-a>1e-8))
				if (CCDF_c<prctile_FW(i)) then 
					b = c 
					c = (a+b)/2.0_dp
					CCDF_c = sum(DBN_vec,Firm_Wealth_vec>=c)
				else 
					a = c 
					c = (a+b)/2.0_dp
					CCDF_c = sum(DBN_vec,Firm_Wealth_vec>=c)
				endif
				!print*, 'a',a,'c',c,'b',b,'CCDF',CCDF_c,'Error', CCDF_c-prctile_FW(i)
			enddo 
			FW_top_x(i)       = c 
			FW_top_x_share(i) = 100.0_dp*sum(Firm_Wealth_vec*DBN_vec,Firm_Wealth_vec>=c)/Mean_Firm_Wealth
		enddo 

			CLOSE(UNIT=11)

	! Define SSE_Moments
		SSE_Moments = (1.0-Wealth_Output/3.0_DP)**2.0_DP  + (1.0_DP-prct1_wealth/0.34_DP)**2.0_DP  + (prct10_wealth-0.69_DP)**2.0_DP &
	                   & + (1.0_DP-Std_Log_Earnings_25_60 / 0.8_DP)**2.0_DP + (1.0_DP-meanhours_25_60/0.4_DP)**2.0_DP &
	                   & + (1.0_DP-MeanReturn/0.069_DP)**2.0_DP


END SUBROUTINE COMPUTE_MOMENTS


!========================================================================================
!========================================================================================
!========================================================================================

SUBROUTINE Hsieh_Klenow_Efficiency(bench_indx)
	IMPLICIT NONE
	integer, intent(in) :: bench_indx
	real(dp), dimension(na,nz,nx) :: K_Mat, TFPR_i=0.0_dp
	real(dp) :: TFP, TFP_star, TFPR_bar, size, K, theta_aux(nz),YBAR_aux,QBAR_aux,K_aux,NBAR_aux
	integer  :: i_a, i_z, i_x, i_theta

	size = 1.0_dp ! sum(DBN1(:,:,:,:,:,1:2))

	K_mat  = K_Matrix(R,P)

	TFPR_bar = 0.0_dp
	K 		 = 0.0_dp 
	do i_a = 1,na
	do i_z = 1,nz 
	do i_x = 1,2
		TFPR_i(i_a,i_z,i_x) = P * xz_grid(i_x,i_z)** mu * K_mat(i_a,i_z,i_x)**(mu-1.0_dp)
		TFPR_bar = TFPR_bar + sum(DBN1(:,i_a,i_z,:,:,i_x))/size * K_mat(i_a,i_z,i_x) / (alpha*QBAR**alpha*NBAR**(1.0_dp-alpha))
		K 		 = K    	+ sum(DBN1(:,i_a,i_z,:,:,i_x))      * K_mat(i_a,i_z,i_x) 
	enddo 
	enddo 
	enddo
	TFPR_bar = 1.0_dp / TFPR_bar

	TFP 	 = 0.0_dp
	TFP_star = 0.0_dp
	do i_a = 1,na
	do i_z = 1,nz 
	do i_x = 1,2
		TFP      = TFP      + sum(DBN1(:,i_a,i_z,:,:,i_x))/size*&
					&	( xz_grid(i_x,i_z) * TFPR_bar / TFPR_i(i_a,i_z,i_x) )**(mu/(1.0_dp-mu))
		TFP_star = TFP_star + sum(DBN1(:,i_a,i_z,:,:,i_x))/size*&
					&	( xz_grid(i_x,i_z)  								)**(mu/(1.0_dp-mu))
	enddo 
	enddo 
	enddo
	TFP 	 = TFP ** (alpha*(1.0_dp-mu)/mu)
	TFP_star = TFP_star ** (alpha*(1.0_dp-mu)/mu)

	! ! Compute output without distortions
	theta_aux = theta
	YBAR_aux  = YBAR
	QBAR_aux  = QBAR 
	NBAR_aux  = NBAR 
	K_aux     = K
	do i_theta = 1,10000,1
	theta     = 4.0_dp + real(i_theta,8)/10.0_dp; print*, ' '; print*, 'theta= ',theta(1); print*, ' '
	CALL FIND_DBN_EQ
	enddo 
	theta     = big_p ; print*, ' ';print*, 'theta= ',theta(1); print*, ' '
	CALL FIND_DBN_EQ
	theta     = theta_aux 
	K         = sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid )
		


	if (bench_indx.eq.1) then
	OPEN(UNIT=10, FILE=trim(Result_Folder)//'Hsieh_Klenow_Efficiency_bench.txt', STATUS='replace')
	else
	OPEN(UNIT=10, FILE=trim(Result_Folder)//'Hsieh_Klenow_Efficiency_exp.txt'  , STATUS='replace')
	endif 

	WRITE(UNIT=10, FMT=*) ' '
	WRITE(UNIT=10, FMT=*) 'Variable ','Distorted_Equilibrium ','Frictionless_Equilibrium ','Gain '
	WRITE(UNIT=10, FMT=*) 'TFP ', TFP , TFP_star, TFP_star/TFP 
	WRITE(UNIT=10, FMT=*) 'Y '  , YBAR_aux , YBAR, YBAR/YBAR_aux
	WRITE(UNIT=10, FMT=*) 'Q '  , QBAR_aux , QBAR, QBAR/QBAR_aux
	WRITE(UNIT=10, FMT=*) 'K '  , K_aux    , K   , K   /K_aux
	WRITE(UNIT=10, FMT=*) 'N '  , NBAR_aux , NBAR, NBAR/NBAR_aux
	WRITE(UNIT=10, FMT=*) 'Check'
	WRITE(UNIT=10, FMT=*) 'YBAR ','TFP*K^a*N^(1-a) ','Q ',' TFP*K ','MeanWealth ','K '
	WRITE(UNIT=10, FMT=*)  YBAR_aux , TFP*K_aux**alpha*NBAR_aux**(1.0_DP-alpha),QBAR_aux,TFP**(1.0_dp/alpha)*K_aux,MeanWealth,K_aux

	CLOSE(UNIT=10)

	if (bench_indx.eq.1) then
	OPEN(UNIT=10, FILE=trim(Result_Folder)//'Bench_Files/TFPR_bench', STATUS='replace')
	else
	OPEN(UNIT=10, FILE=trim(Result_Folder)//'Exp_Files/TFPR_exp'  , STATUS='replace')
	endif 

	WRITE(UNIT=10, FMT=*)  TFPR_i

	CLOSE(UNIT=10)

	print*, ' '
	print*, 'Variable ','Distorted_Equilibrium ','Frictionless_Equilibrium ','Gain '
	print*, 'TFP ', TFP , TFP_star, TFP_star/TFP 
	print*, 'Y '  , YBAR_aux , YBAR, YBAR/YBAR_aux
	print*, 'Q '  , QBAR_aux , QBAR, QBAR/QBAR_aux
	print*, 'K '  , K_aux    , K   , K   /K_aux
	print*, 'N '  , NBAR_aux , NBAR, NBAR/NBAR_aux
	print*, 'Check'
	print*, 'YBAR ','TFP*K^a*N^(1-a) ','Q ',' TFP*K ','MeanWealth ','K '
	print*,  YBAR_aux , TFP*K_aux**alpha*NBAR_aux**(1.0_DP-alpha),QBAR_aux,TFP**(1.0_dp/alpha)*K_aux,MeanWealth,K_aux



END SUBROUTINE Hsieh_Klenow_Efficiency

!========================================================================================
!========================================================================================
!========================================================================================

SUBROUTINE WRITE_VARIABLES(bench_indx)
	IMPLICIT NONE
	integer, intent(in) :: bench_indx
	integer :: prctile, status, zi

	print*, ' Writing summary of results in output.txt'
	if (bench_indx.eq.1) then
		OPEN (UNIT=19, FILE=trim(Result_Folder)//'output.txt', STATUS='replace') 
			WRITE(UNIT=19, FMT=*) "Parameters"
			WRITE(UNIT=19, FMT=*) params
			WRITE(UNIT=19, FMT=*) "NSU_Switch",NSU_Switch
 			WRITE(UNIT=19, FMT=*) "sigma",sigma,'gamma',gamma,'phi',phi,'beta',beta
			WRITE(UNIT=19, FMT=*) 'Theta',theta, 'mu',mu
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) "Results for benchmark economy"
			WRITE(UNIT=19, FMT=*) ' '
	else
		OPEN (UNIT=19, FILE=trim(Result_Folder)//'output_exp.txt', STATUS='replace') 
		! OPEN (UNIT=19, FILE=trim(Result_Folder)//'output.txt', STATUS='old', POSITION='append', iostat=status) 
		! 	if (status.ne.0) then 
		! 	OPEN (UNIT=19, FILE=trim(Result_Folder)//'output.txt', STATUS='replace') 
		! 	end if 
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) 'Wealth Taxes'
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) "Threshold_Factor"     	, Threshold_Factor
			WRITE(UNIT=19, FMT=*) "Wealth_Factor"		  	, Wealth_Factor
			WRITE(UNIT=19, FMT=*) "Threshold"			  	, Y_a_Threshold
			WRITE(UNIT=19, FMT=*) "Wealth_Tax_Above"	    , TauW_at
			WRITE(UNIT=19, FMT=*) 'Share_Above_Threshold'	, Threshold_Share
			WRITE(UNIT=19, FMT=*) 'Psi'						, Psi
			WRITE(UNIT=19, FMT=*) 'Tau_K'					, TauK
			WRITE(UNIT=19, FMT=*) 'Tau_C'					, TauC
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) "Results for experimental economy"
			WRITE(UNIT=19, FMT=*) ' '

	end if 
			WRITE(UNIT=19, FMT=*) 'Aggregate Variables'
			WRITE(UNIT=19, FMT=*) 'GBAR'	, GBAR
			WRITE(UNIT=19, FMT=*) 'KBAR'	, MeanWealth
			WRITE(UNIT=19, FMT=*) 'QBAR'	, QBAR
			WRITE(UNIT=19, FMT=*) 'NBAR'	, NBAR
			WRITE(UNIT=19, FMT=*) 'EBAR'	, EBAR
			WRITE(UNIT=19, FMT=*) 'YBAR'	, YBAR
			WRITE(UNIT=19, FMT=*) 'CBAR'    , MeanCons
			WRITE(UNIT=19, FMT=*) 'P'		, 100.0_dp*P
			WRITE(UNIT=19, FMT=*) 'wage'	, wage
			WRITE(UNIT=19, FMT=*) 'R'		, 100.0_dp*R
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) 'After_Tax_Prices'
			WRITE(UNIT=19, FMT=*) 'R_AT' 	, 100.0_dp*(-tauW_at+(1.0_dp-tauK)*R)
			WRITE(UNIT=19, FMT=*) 'wage_AT' , psi*wage
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) 'Moments:'
			WRITE(UNIT=19, FMT=*) 'Debt_Output'		  	, External_Debt_GDP
			WRITE(UNIT=19, FMT=*) 'Wealth_Output'	  	, Wealth_Output
			if (A_C.gt.0.0_dp) then 
			WRITE(UNIT=19, FMT=*) 'TFP_Q' 				, QBAR/K_P
			else 
			WRITE(UNIT=19, FMT=*) 'TFP_Q' 				, QBAR/MeanWealth
			endif 
			WRITE(UNIT=19, FMT=*) 'Wealth_Top_0.01%' 	, prct9999_wealth
			WRITE(UNIT=19, FMT=*) 'Wealth_Top_0.1%' 	, prct999_wealth
			WRITE(UNIT=19, FMT=*) 'Wealth_Top_1%' 		, prct1_wealth
			WRITE(UNIT=19, FMT=*) 'Wealth_Top_10%'		, prct10_wealth
			WRITE(UNIT=19, FMT=*) 'Wealth_Top_20%'		, prct20_wealth
			WRITE(UNIT=19, FMT=*) 'Wealth_Top_40%'		, prct40_wealth
			WRITE(UNIT=19, FMT=*) 'Mean_Labor_Earnings'	, mean_log_earnings_25_60
			WRITE(UNIT=19, FMT=*) 'STD_Labor_Earnings'	, Std_Log_Earnings_25_60
			WRITE(UNIT=19, FMT=*) 'Mean_Labor_25_60'	, meanhours_25_60
			WRITE(UNIT=19, FMT=*) 'Mean_Return'			, 100.0_dp*MeanReturn
			WRITE(UNIT=19, FMT=*) 'Std_Return'			, StdReturn
			do zi=1,nz
			WRITE(UNIT=19, FMT=*) 'Mean_Return_by_z'	, 100.0_dp*MeanReturn_by_z(zi)
			enddo 
			! WRITE(UNIT=19, FMT=*) ' '
			! WRITE(UNIT=19, FMT=*) 'Present_Value_Wealth'
			! WRITE(UNIT=19, FMT=*) "Mean_PV_Wealth"		    , Mean_Firm_Wealth
			! WRITE(UNIT=19, FMT=*) 'PV_Wealth_Top_0.01%' 	, FW_top_x_share(6)
			! WRITE(UNIT=19, FMT=*) 'PV_Wealth_Top_0.1%' 		, FW_top_x_share(5)
			! WRITE(UNIT=19, FMT=*) 'PV_Wealth_Top_1%' 		, FW_top_x_share(4)
			! WRITE(UNIT=19, FMT=*) 'PV_Wealth_Top_10%'		, FW_top_x_share(3)
			! WRITE(UNIT=19, FMT=*) 'PV_Wealth_Top_20%'		, FW_top_x_share(2)
			! WRITE(UNIT=19, FMT=*) 'PV_Wealth_Top_40%'		, FW_top_x_share(1)
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) 'Bequest'
			WRITE(UNIT=19, FMT=*) 'Total_Bequest_Wealth'	, Bequest_Wealth/MeanWealth 
			WRITE(UNIT=19, FMT=*) 'Mean_Bequest_Wealth'	    , Mean_Bequest/MeanWealth 
			WRITE(UNIT=19, FMT=*) 'BQ/Inc_for_90th_pct' 	, Bq_Inc(3,:)
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) 'Labor'
			WRITE(UNIT=19, FMT=*) 'Fraction_of_workers'		, Size_Frisch/sum(DBN1(1:RetAge-1,:,:,:,:,:))
			WRITE(UNIT=19, FMT=*) 'Av.Hours'				, Hours_Frisch
			WRITE(UNIT=19, FMT=*) 'Av.Frisch_Elasticity'   	, Frisch_Elasticity
			WRITE(UNIT=19, FMT=*) 'Frisch_Elasticity_Av'   	, Frisch_Elasticity_2
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) 'Taxes'
			WRITE(UNIT=19, FMT=*) 'Tax_Rev/GDP'				, (GBAR_K+GBAR_W+GBAR_L+GBAR_C+GBAR_BQ)/YBAR
			WRITE(UNIT=19, FMT=*) 'Capital_Tax/Total_Tax'	, (GBAR_K+GBAR_W)/(GBAR_K+GBAR_W+GBAR_L+GBAR_C+GBAR_BQ)
			WRITE(UNIT=19, FMT=*) 'Capital_Tax/_GDP'		, (GBAR_K+GBAR_W)/YBAR
			WRITE(UNIT=19, FMT=*) 'Labor_Tax/Total_Tax'		, GBAR_L/(GBAR_K+GBAR_W+GBAR_L+GBAR_C+GBAR_BQ)
			WRITE(UNIT=19, FMT=*) 'Labor_Tax/GDP'			, GBAR_L/YBAR
			WRITE(UNIT=19, FMT=*) 'Average_Labor_Tax'		, GBAR_L/Tot_Lab_Inc
			WRITE(UNIT=19, FMT=*) 'Cons_Tax/Total_Tax'		, GBAR_C/(GBAR_K+GBAR_W+GBAR_L+GBAR_C+GBAR_BQ)
			WRITE(UNIT=19, FMT=*) 'Cons_Tax/GDP'			, GBAR_C/YBAR
			WRITE(UNIT=19, FMT=*) 'Estate_Tax/Total_Tax'	, GBAR_BQ/(GBAR_K+GBAR_W+GBAR_L+GBAR_C+GBAR_BQ)
			WRITE(UNIT=19, FMT=*) 'Estate_Tax/GDP'			, GBAR_BQ/YBAR
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) 'Public_Corporate_Sector'
			WRITE(UNIT=19, FMT=*) 'Y_C/YBAR'				, 100.0_dp*YBAR_C/YBAR
			WRITE(UNIT=19, FMT=*) 'K_C/KBAR'				, 100.0_dp*K_C/MeanWealth
			WRITE(UNIT=19, FMT=*) 'L_C/NBAR'   				, 100.0_dp*L_C/NBAR
			WRITE(UNIT=19, FMT=*) 'Y_C'		   				, YBAR_C 
			WRITE(UNIT=19, FMT=*) 'L_C'		   				, L_C 
			WRITE(UNIT=19, FMT=*) 'K_C'		   				, K_C 
			WRITE(UNIT=19, FMT=*) 'Y_P'		   				, YBAR_P
			WRITE(UNIT=19, FMT=*) 'L_P'		   				, L_P 
			WRITE(UNIT=19, FMT=*) 'K_P'		   				, K_P 
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) 'Share_Entrepreneurs'
			WRITE(UNIT=19, FMT=*) 'Pr/BT_Income>10%'		, 100.0_dp*Entrepreneur_10
			WRITE(UNIT=19, FMT=*) 'Pr/BT_Income>50%'		, 100.0_dp*Entrepreneur_50
		CLOSE(Unit=19)
	if (bench_indx.ne.1) then
		OPEN (UNIT=19, FILE=trim(Result_Folder)//'output_exp.txt', STATUS='old', POSITION='append') 
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) 'Welfare and output gain'
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) "CE1_Pop(bench)" , Welfare_Gain_Pop_bench
			WRITE(UNIT=19, FMT=*) "CE1_Pop(exp)"   , Welfare_Gain_Pop_exp
			WRITE(UNIT=19, FMT=*) "CE1_NB(bench)"  , Welfare_Gain_NB_bench
			WRITE(UNIT=19, FMT=*) "CE1_NB(exp)"    , Welfare_Gain_NB_exp
			WRITE(UNIT=19, FMT=*) "Frac_pos_welfare(bench)" , 100.0_dp*frac_pos_welfare
			WRITE(UNIT=19, FMT=*) "CE2_Pop(exp)"   , Av_Util_Pop
			WRITE(UNIT=19, FMT=*) "CE2_NB(exp)"	   , Av_Util_NB
			WRITE(UNIT=19, FMT=*) ' '
		CLOSE(Unit=19)
	end if 

			
END SUBROUTINE WRITE_VARIABLES

!========================================================================================
!========================================================================================
!========================================================================================

SUBROUTINE Write_Benchmark_Results(Compute_bench)
	IMPLICIT NONE
	logical :: Compute_bench
	character(100) :: bench_folder

	bench_folder = trim(Result_Folder)//'Bench_Files/'
		call system( 'mkdir -p ' // trim(bench_folder) )
		print*,' '; print*, "Bench Files Folder:", bench_folder
	
	IF (Compute_bench) then 
		OPEN  (UNIT=1,  FILE=trim(bench_folder)//'cons'  , STATUS='replace')
		WRITE (UNIT=1,  FMT=*) cons
		CLOSE (unit=1)
		OPEN  (UNIT=2,  FILE=trim(bench_folder)//'aprime', STATUS='replace')
		WRITE (UNIT=2,  FMT=*) aprime
		CLOSE (unit=2)
		OPEN  (UNIT=3,  FILE=trim(bench_folder)//'hours' , STATUS='replace')
		WRITE (UNIT=3,  FMT=*) hours
		CLOSE (unit=3)
		OPEN  (UNIT=4,  FILE=trim(bench_folder)//'value' , STATUS='replace')
		WRITE (UNIT=4,  FMT=*) ValueFunction
		CLOSE (unit=4)
		OPEN  (UNIT=4,  FILE=trim(bench_folder)//'bq_value', STATUS='replace')
		WRITE (UNIT=4,  FMT=*) Bq_Value
		CLOSE (unit=4)
		OPEN  (UNIT=4,  FILE=trim(bench_folder)//'v_pr'  , STATUS='replace')
		WRITE (UNIT=4,  FMT=*) V_Pr
		CLOSE (unit=4)
		OPEN  (UNIT=4,  FILE=trim(bench_folder)//'v_pr_nb', STATUS='replace')
		WRITE (UNIT=4,  FMT=*) V_Pr_nb
		CLOSE (unit=4)
		OPEN  (UNIT=4,  FILE=trim(bench_folder)//'Income_AT', STATUS='replace')
		WRITE (UNIT=4,  FMT=*) Income_AT
		CLOSE (unit=4)

		OPEN  (UNIT=5,  FILE=trim(bench_folder)//'DBN'   , STATUS='replace')
		WRITE (UNIT=5,  FMT=*) DBN1 
		CLOSE (UNIT=5)
		OPEN  (UNIT=60,  FILE=trim(bench_folder)//'GBAR'  , STATUS='replace')
		WRITE (UNIT=60,  FMT=*) GBAR
		CLOSE (UNIT=60)
		OPEN  (UNIT=7,  FILE=trim(bench_folder)//'EBAR'  , STATUS='replace')
		WRITE (UNIT=7,  FMT=*) EBAR
		CLOSE (UNIT=7)
		OPEN  (UNIT=8,  FILE=trim(bench_folder)//'NBAR'  , STATUS='replace')
		WRITE (UNIT=8,  FMT=*) NBAR
		CLOSE (UNIT=8)
		OPEN  (UNIT=9,  FILE=trim(bench_folder)//'QBAR'  , STATUS='replace')
		WRITE (UNIT=9,  FMT=*) QBAR
		CLOSE (UNIT=9)
		OPEN  (UNIT=10, FILE=trim(bench_folder)//'P'    , STATUS='replace')
		WRITE (UNIT=10, FMT=*) P
		CLOSE (UNIT=10)
		OPEN  (UNIT=10, FILE=trim(bench_folder)//'R'    , STATUS='replace')
		WRITE (UNIT=10, FMT=*) R
		CLOSE (UNIT=10)
		OPEN  (UNIT=11, FILE=trim(bench_folder)//'wage'  , STATUS='replace')
		WRITE (UNIT=11, FMT=*) wage 
		CLOSE (UNIT=11)
		OPEN  (UNIT=12, FILE=trim(bench_folder)//'YBAR'  , STATUS='replace')
		WRITE (UNIT=12, FMT=*) YBAR
		CLOSE (UNIT=12)

		OPEN  (UNIT=12, FILE=trim(bench_folder)//'tauK'  , STATUS='replace')
		WRITE (UNIT=12, FMT=*) tauK
		CLOSE (UNIT=12)
		OPEN  (UNIT=12, FILE=trim(bench_folder)//'tauPL'  , STATUS='replace')
		WRITE (UNIT=12, FMT=*) tauPL
		CLOSE (UNIT=12)
		OPEN  (UNIT=12, FILE=trim(bench_folder)//'psi'  , STATUS='replace')
		WRITE (UNIT=12, FMT=*) psi
		CLOSE (UNIT=12)
		OPEN  (UNIT=12, FILE=trim(bench_folder)//'tauW_bt'  , STATUS='replace')
		WRITE (UNIT=12, FMT=*) tauW_bt
		CLOSE (UNIT=12)
		OPEN  (UNIT=12, FILE=trim(bench_folder)//'tauW_at'  , STATUS='replace')
		WRITE (UNIT=12, FMT=*) tauW_at
		CLOSE (UNIT=12)

		OPEN  (UNIT=12, FILE=trim(bench_folder)//'SSC_Payments'  , STATUS='replace')
		WRITE (UNIT=12, FMT=*) SSC_Payments
		CLOSE (UNIT=12)

		OPEN  (UNIT=12, FILE=trim(bench_folder)//'YBAR_C'  	, STATUS='replace')
		OPEN  (UNIT=13, FILE=trim(bench_folder)//'L_C'  	, STATUS='replace')
		OPEN  (UNIT=14, FILE=trim(bench_folder)//'K_C'  	, STATUS='replace')
		OPEN  (UNIT=15, FILE=trim(bench_folder)//'YBAR_P'  	, STATUS='replace')
		OPEN  (UNIT=16, FILE=trim(bench_folder)//'L_P'  	, STATUS='replace')
		OPEN  (UNIT=17, FILE=trim(bench_folder)//'K_P'  	, STATUS='replace')
		WRITE (UNIT=12, FMT=*) YBAR_C
		WRITE (UNIT=13, FMT=*) L_C
		WRITE (UNIT=14, FMT=*) K_C
		WRITE (UNIT=15, FMT=*) YBAR_P
		WRITE (UNIT=16, FMT=*) L_P
		WRITE (UNIT=17, FMT=*) K_P
		CLOSE (UNIT=12); CLOSE (UNIT=13); CLOSE (UNIT=14); CLOSE (UNIT=15); CLOSE (UNIT=16); CLOSE (UNIT=17); 

		print*, "Writing of benchmark results completed"; print*, ' '
	ELSE 
		OPEN (UNIT=1,  FILE=trim(bench_folder)//'cons'  , STATUS='old', ACTION='read')
		OPEN (UNIT=2,  FILE=trim(bench_folder)//'aprime', STATUS='old', ACTION='read')
		OPEN (UNIT=3,  FILE=trim(bench_folder)//'hours' , STATUS='old', ACTION='read')
		OPEN (UNIT=4,  FILE=trim(bench_folder)//'value' , STATUS='old', ACTION='read')
		OPEN (UNIT=70, FILE=trim(bench_folder)//'bq_value', STATUS='old', ACTION='read')
		OPEN (UNIT=5,  FILE=trim(bench_folder)//'DBN'   , STATUS='old', ACTION='read')
		OPEN (UNIT=60, FILE=trim(bench_folder)//'GBAR'  , STATUS='old', ACTION='read')
		OPEN (UNIT=7,  FILE=trim(bench_folder)//'EBAR'  , STATUS='old', ACTION='read')
		OPEN (UNIT=8,  FILE=trim(bench_folder)//'NBAR'  , STATUS='old', ACTION='read')
		OPEN (UNIT=9,  FILE=trim(bench_folder)//'QBAR'  , STATUS='old', ACTION='read')
		OPEN (UNIT=10, FILE=trim(bench_folder)//'P'     , STATUS='old', ACTION='read')
		OPEN (UNIT=11, FILE=trim(bench_folder)//'R'     , STATUS='old', ACTION='read')
		OPEN (UNIT=12, FILE=trim(bench_folder)//'wage'  , STATUS='old', ACTION='read')
		OPEN (UNIT=13, FILE=trim(bench_folder)//'YBAR'  , STATUS='old', ACTION='read')

		OPEN (UNIT=14, FILE=trim(bench_folder)//'tauK'    , STATUS='old', ACTION='read')
		OPEN (UNIT=15, FILE=trim(bench_folder)//'tauPL'   , STATUS='old', ACTION='read')
		OPEN (UNIT=16, FILE=trim(bench_folder)//'psi'  	  , STATUS='old', ACTION='read')
		OPEN (UNIT=17, FILE=trim(bench_folder)//'tauW_bt' , STATUS='old', ACTION='read')
		OPEN (UNIT=18, FILE=trim(bench_folder)//'tauW_at' , STATUS='old', ACTION='read')

		OPEN (UNIT=19, FILE=trim(bench_folder)//'v_pr'    , STATUS='old', ACTION='read')
		! OPEN (UNIT=20, FILE=trim(bench_folder)//'v_pr_nb' , STATUS='old', ACTION='read')		

		OPEN (UNIT=21, FILE=trim(bench_folder)//'SSC_Payments' , STATUS='old', ACTION='read')
		OPEN (UNIT=22, FILE=trim(bench_folder)//'Income_AT' , STATUS='old', ACTION='read')

		OPEN (UNIT=23, FILE=trim(bench_folder)//'YBAR_C' , STATUS='old', ACTION='read')
		OPEN (UNIT=24, FILE=trim(bench_folder)//'L_C'    , STATUS='old', ACTION='read')
		OPEN (UNIT=25, FILE=trim(bench_folder)//'K_C'    , STATUS='old', ACTION='read')
		OPEN (UNIT=26, FILE=trim(bench_folder)//'YBAR_P' , STATUS='old', ACTION='read')		
		OPEN (UNIT=27, FILE=trim(bench_folder)//'L_P'    , STATUS='old', ACTION='read')		
		OPEN (UNIT=28, FILE=trim(bench_folder)//'K_P'    , STATUS='old', ACTION='read')		

		READ (UNIT=1,  FMT=*) cons
		READ (UNIT=2,  FMT=*) aprime
		READ (UNIT=3,  FMT=*) hours
		READ (UNIT=4,  FMT=*) ValueFunction
		READ (UNIT=70, FMT=*) Bq_Value
		READ (UNIT=5,  FMT=*) DBN1 
		READ (UNIT=60, FMT=*) GBAR 
		READ (UNIT=7,  FMT=*) EBAR
		READ (UNIT=8,  FMT=*) NBAR
		READ (UNIT=9,  FMT=*) QBAR
		READ (UNIT=10, FMT=*) P
		READ (UNIT=11, FMT=*) R
		READ (UNIT=12, FMT=*) wage 
		READ (UNIT=13, FMT=*) YBAR

		READ (UNIT=14, FMT=*) tauK
		READ (UNIT=15, FMT=*) tauPL
		READ (UNIT=16, FMT=*) psi
		READ (UNIT=17, FMT=*) tauW_bt
		READ (UNIT=18, FMT=*) tauW_at

		READ (UNIT=19, FMT=*) V_Pr
		! READ (UNIT=20, FMT=*) V_Pr_nb
		READ (UNIT=21, FMT=*) SSC_Payments
		READ (UNIT=22, FMT=*) Income_AT

		READ (UNIT=23, FMT=*) YBAR_C
		READ (UNIT=24, FMT=*) L_C
		READ (UNIT=25, FMT=*) K_C
		READ (UNIT=26, FMT=*) YBAR_P
		READ (UNIT=27, FMT=*) L_P
		READ (UNIT=28, FMT=*) K_P

		CLOSE (unit=1); CLOSE (unit=2); CLOSE (unit=3); CLOSE (unit=4); CLOSE (unit=70); CLOSE (unit=5)
		CLOSE (unit=60); CLOSE (unit=7); CLOSE (unit=8); CLOSE (unit=9); CLOSE (unit=10)
		CLOSE (unit=11); CLOSE (unit=12); CLOSE (unit=13); CLOSE (unit=14); CLOSE (unit=15)
		CLOSE (unit=16); CLOSE (unit=17); CLOSE (unit=18); CLOSE (unit=19); !CLOSE (unit=20)
		CLOSE (unit=21); CLOSE (unit=22); 
		CLOSE (unit=23); CLOSE (unit=24); CLOSE (unit=25); CLOSE (unit=26); CLOSE (unit=27); CLOSE (unit=28);

		print*, "Reading of benchmark results completed"; print*, ' ';
	END IF 
END SUBROUTINE Write_Benchmark_Results


SUBROUTINE Write_Experimental_Results(compute_exp)
	IMPLICIT NONE
	logical, intent(in) :: compute_exp

	call system( 'mkdir -p ' // trim(Result_Folder) // 'Exp_Files/' )

	if (compute_exp) then 
		print*,' '; print*, "Writing experimental results in folder", trim(Result_Folder) // 'Exp_Files/'
		OPEN  (UNIT=1 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_cons'  , STATUS='replace')
		WRITE (UNIT=1 , FMT=*) cons
		OPEN  (UNIT=2 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_aprime', STATUS='replace')
		WRITE (UNIT=2 , FMT=*) aprime
		OPEN  (UNIT=3 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_hours' , STATUS='replace')
		WRITE (UNIT=3 , FMT=*) hours
		OPEN  (UNIT=4 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_value' , STATUS='replace')
		WRITE (UNIT=4 , FMT=*) ValueFunction
		OPEN  (UNIT=70, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_bq_value' , STATUS='replace')
		WRITE (UNIT=70, FMT=*) Bq_Value

		OPEN  (UNIT=5 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_DBN'   , STATUS='replace')
		WRITE (UNIT=5 , FMT=*) DBN1 
		OPEN  (UNIT=60, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_GBAR' , STATUS='replace')
		WRITE (UNIT=60, FMT=*) GBAR
		OPEN  (UNIT=7 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_EBAR'  , STATUS='replace')
		WRITE (UNIT=7 , FMT=*) EBAR
		OPEN  (UNIT=8 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_NBAR'  , STATUS='replace')
		WRITE (UNIT=8 , FMT=*) NBAR
		OPEN  (UNIT=9 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_QBAR'  , STATUS='replace')
		WRITE (UNIT=9 , FMT=*) QBAR
		OPEN  (UNIT=10, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_P'     , STATUS='replace')
		WRITE (UNIT=10, FMT=*) P
		OPEN  (UNIT=11, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_R'     , STATUS='replace')
		WRITE (UNIT=11, FMT=*) R
		OPEN  (UNIT=12, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_wage'  , STATUS='replace')
		WRITE (UNIT=12, FMT=*) wage 
		OPEN  (UNIT=13, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_YBAR'  , STATUS='replace')
		WRITE (UNIT=13, FMT=*) YBAR

		OPEN  (UNIT=14, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_tauK'   , STATUS='replace')
		WRITE (UNIT=14, FMT=*) tauK
		OPEN  (UNIT=15, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_psi'  	 , STATUS='replace')
		WRITE (UNIT=15, FMT=*) psi
		OPEN  (UNIT=16, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_tauPL'  , STATUS='replace')
		WRITE (UNIT=16, FMT=*) tauPL
		OPEN  (UNIT=17, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_tauW_bt', STATUS='replace')
		WRITE (UNIT=17, FMT=*) tauW_bt
		OPEN  (UNIT=18, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_tauW_at', STATUS='replace')
		WRITE (UNIT=18, FMT=*) tauW_at

		OPEN  (UNIT=19,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_v_pr'  , STATUS='replace')
		WRITE (UNIT=19,  FMT=*) V_Pr
		OPEN  (UNIT=20,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_v_pr_nb', STATUS='replace')
		WRITE (UNIT=20,  FMT=*) V_Pr_nb

		OPEN  (UNIT=21,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_Income_AT', STATUS='replace')
		WRITE (UNIT=21,  FMT=*) Income_AT

		OPEN  (UNIT=22,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_YBAR_C', STATUS='replace')
		WRITE (UNIT=22,  FMT=*) YBAR_C

		OPEN  (UNIT=23,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_L_C', STATUS='replace')
		WRITE (UNIT=23,  FMT=*) L_C

		OPEN  (UNIT=24,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_K_C', STATUS='replace')
		WRITE (UNIT=24,  FMT=*) K_C

		OPEN  (UNIT=25,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_YBAR_P', STATUS='replace')
		WRITE (UNIT=25,  FMT=*) YBAR_P

		OPEN  (UNIT=26,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_L_P', STATUS='replace')
		WRITE (UNIT=26,  FMT=*) L_P

		OPEN  (UNIT=27,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_K_P', STATUS='replace')
		WRITE (UNIT=27,  FMT=*) K_P

		print*, "Writing of experimental results completed"; print*, ' '

	else 
		print*,' '; print*, "Reading experimental results from folder", trim(Result_Folder) // 'Exp_Files/'
		OPEN (UNIT=1 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_cons'  	, STATUS='old', ACTION='read')
		OPEN (UNIT=2 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_aprime'	, STATUS='old', ACTION='read')
		OPEN (UNIT=3 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_hours' 	, STATUS='old', ACTION='read')
		OPEN (UNIT=4 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_value' 	, STATUS='old', ACTION='read')
		OPEN (UNIT=70, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_bq_value', STATUS='old', ACTION='read')
		OPEN (UNIT=5 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_DBN'   	, STATUS='old', ACTION='read')
		OPEN (UNIT=60, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_GBAR'  	, STATUS='old', ACTION='read')
		OPEN (UNIT=7 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_EBAR'  	, STATUS='old', ACTION='read')
		OPEN (UNIT=8 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_NBAR'  	, STATUS='old', ACTION='read')
		OPEN (UNIT=9 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_QBAR'  	, STATUS='old', ACTION='read')
		OPEN (UNIT=10, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_P'     	, STATUS='old', ACTION='read')
		OPEN (UNIT=11, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_R'     	, STATUS='old', ACTION='read')
		OPEN (UNIT=12, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_wage'  	, STATUS='old', ACTION='read')
		OPEN (UNIT=13, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_YBAR'  	, STATUS='old', ACTION='read')
		OPEN (UNIT=14, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_tauK'	, STATUS='old', ACTION='read')
		OPEN (UNIT=15, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_psi'	, STATUS='old', ACTION='read')
		OPEN (UNIT=16, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_tauPL'	, STATUS='old', ACTION='read')
		OPEN (UNIT=17, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_tauW_bt', STATUS='old', ACTION='read')
		OPEN (UNIT=18, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_tauW_at', STATUS='old', ACTION='read')
		OPEN (UNIT=19, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_v_pr'   , STATUS='old', ACTION='read')
		OPEN (UNIT=20, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_v_pr_nb', STATUS='old', ACTION='read')
		OPEN (UNIT=21, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_Income_AT', STATUS='old', ACTION='read')
		OPEN (UNIT=22, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_YBAR_C' , STATUS='old', ACTION='read')
		OPEN (UNIT=23, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_L_C' 	, STATUS='old', ACTION='read')
		OPEN (UNIT=24, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_K_C' 	, STATUS='old', ACTION='read')
		OPEN (UNIT=25, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_YBAR_P' , STATUS='old', ACTION='read')
		OPEN (UNIT=26, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_L_P' 	, STATUS='old', ACTION='read')
		OPEN (UNIT=27, FILE=trim(Result_Folder)//'Exp_Files/Exp_results_K_P' 	, STATUS='old', ACTION='read')

		READ (UNIT=1,  FMT=*) cons
		READ (UNIT=2,  FMT=*) aprime
		READ (UNIT=3,  FMT=*) hours
		READ (UNIT=4,  FMT=*) ValueFunction
		READ (UNIT=70, FMT=*) Bq_Value
		READ (UNIT=5,  FMT=*) DBN1 
		READ (UNIT=60, FMT=*) GBAR 
		READ (UNIT=7,  FMT=*) EBAR
		READ (UNIT=8,  FMT=*) NBAR
		READ (UNIT=9,  FMT=*) QBAR
		READ (UNIT=10, FMT=*) P
		READ (UNIT=11, FMT=*) R
		READ (UNIT=12, FMT=*) wage 
		READ (UNIT=13, FMT=*) YBAR
		READ (UNIT=14, FMT=*) tauK
		READ (UNIT=15, FMT=*) psi
		READ (UNIT=16, FMT=*) tauPL
		READ (UNIT=17, FMT=*) tauW_bt
		READ (UNIT=18, FMT=*) tauW_at
		READ (UNIT=19, FMT=*) V_Pr
		READ (UNIT=20, FMT=*) V_Pr_nb
		READ (UNIT=21, FMT=*) Income_AT
		READ (UNIT=22, FMT=*) YBAR_C 
		READ (UNIT=23, FMT=*) L_C 
		READ (UNIT=24, FMT=*) K_C 
		READ (UNIT=25, FMT=*) YBAR_P
		READ (UNIT=26, FMT=*) L_P 
		READ (UNIT=27, FMT=*) K_P 
		print*, "Reading of experimental results completed"; print*, ' '
	endif 

	CLOSE (unit=1); CLOSE (unit=2); CLOSE (unit=3); CLOSE (unit=4); CLOSE (unit=70);
	CLOSE (unit=5); CLOSE (unit=60); CLOSE (unit=7); CLOSE (unit=8); CLOSE (unit=9);
	CLOSE (unit=10); CLOSE (unit=11); CLOSE (unit=12); CLOSE (unit=13); CLOSE (unit=14);
	CLOSE (unit=15); CLOSE (unit=16); CLOSE (unit=17); CLOSE (unit=18); CLOSE (unit=19)
	CLOSE (unit=20); CLOSE (unit=21); CLOSE (unit=22); CLOSE (unit=23); CLOSE (unit=24); 
	CLOSE (unit=25); CLOSE (unit=26); CLOSE (unit=27); 

END SUBROUTINE Write_Experimental_Results


! !========================================================================================
! !========================================================================================
! !========================================================================================

Function Tax_Reform_Welfare(tk)
	implicit none 
	real(dp), intent(in)  :: tk 
	real(dp) :: Tax_Reform_Welfare
	real(dp) :: tau_diff, tau_tol

	!====================================================================================================
	PRINT*,''
	Print*,'--------------- SOLVING EXPERIMENT WITH BEST PARAMETERS -----------------'
	PRINT*,''
	print*,'Wealth Tax Economy','tauK',tk
	tau_tol = 0.0000001_dp 
	! Experiment economy
		solving_bench=0
	! Set capital taxes to tk
		tauK = tk
	! Set initial wealth taxes
		tauWmin_bt  = tauW_bt
		tauWmin_at  = tauW_at
		tauWinc_bt  = 0.002_dp
		tauWinc_at  = 0.002_dp
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
			! Solve the model
				CALL FIND_DBN_EQ
				CALL GOVNT_BUDGET(.true.)
			! Get new G
				GBAR_exp = GBAR 
			tauWindx = 0.0_DP

			if (GBAR_exp.lt.GBAR_bench) then 
				! Solve for the model increasing wealth taxes until revenue is enough to finance G_benchamark
				DO WHILE (GBAR_exp .lt. GBAR_bench)
					! Set old G and new value of tauW
					GBAR_exp_old = GBAR_exp
					tauW_bt = tauWmin_bt + tauWindx * tauWinc_bt
					tauW_at = tauWmin_at + tauWindx * tauWinc_at
					! Solve the model
					CALL FIND_DBN_EQ
					CALL GOVNT_BUDGET(.true.)

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
			elseif (GBAR_exp.gt.GBAR_bench) then 
				! Solve for the model increasing wealth taxes until revenue is enough to finance G_benchamark
				DO WHILE (GBAR_exp .gt. GBAR_bench)
					! Set old G and new value of tauW
					GBAR_exp_old = GBAR_exp
					tauW_bt = tauWmin_bt - tauWindx * tauWinc_bt
					tauW_at = tauWmin_at - tauWindx * tauWinc_at
					! Solve the model
					CALL FIND_DBN_EQ
					CALL GOVNT_BUDGET(.true.)

					! Get new G
					GBAR_exp = GBAR 
					! Iteratioins  
					tauWindx = tauWindx + 1.0_DP   
					write(*,*) "Bracketing GBAR: tauW_bt=", tauW_bt*100, "And tauW_at=", tauW_at*100
					print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
					print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
				ENDDO
				! Set tauW as weighted average of point in  the grid to balance budget more precisely
				tauW_up_bt  = tauW_bt  +  tauWinc_bt
				tauW_low_bt = tauW_bt  
				tauW_bt     = tauW_up_bt - tauWinc_bt * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
				tauW_up_at  = tauW_at  +  tauWinc_at  
				tauW_low_at = tauW_at
				tauW_at     = tauW_up_at - tauWinc_at * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
			endif 

				print*,''
				print*,'GBAR bracketed by taxes:'
				print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
				print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
				print*,''

			! Solve (again) experimental economy
				CALL FIND_DBN_EQ
				CALL GOVNT_BUDGET(.true.)

			! Find tauW that exactly balances the budget (up to precisioin 0.1) using bisection
				GBAR_exp = GBAR
				tau_diff = tauW_up_at - tauW_low_at
				print*,"Gbar at midpoint of bracket and GBAR at benchmark"
				print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
				print*,''
				print*,'Bisection for TauW:'
				DO WHILE ((  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.001 ).and.(tau_diff.gt.tau_tol))! as long as the difference is greater than 0.1% continue
				    if (GBAR_exp .gt. GBAR_bench ) then
				        tauW_up_bt  = tauW_bt 
				        tauW_up_at  = tauW_at 
				    else
				        tauW_low_bt = tauW_bt
				        tauW_low_at = tauW_at
				    endif
				    tau_diff = tauW_up_at - tauW_low_at
				    tauW_bt = (tauW_low_bt + tauW_up_bt)/2.0_DP
				    tauW_at = (tauW_low_at + tauW_up_at)/2.0_DP
				    CALL FIND_DBN_EQ
				    CALL GOVNT_BUDGET(.true.)
				    GBAR_exp = GBAR
				    print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
					print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
					print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
					print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
				ENDDO

		! Compute value function and store policy functions, value function and distribution in file
		! CALL COMPUTE_VALUE_FUNCTION_SPLINE 
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
		CALL Firm_Value

	
	! CALL Write_Experimental_Results(.true.)
	CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
	K_mat  = K_Matrix(R,P)
	Pr_mat = Profit_Matrix(R,P)
	CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
	CALL ComputeLaborUnits(EBAR,wage)
	CALL GOVNT_BUDGET(.true.)


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


	! Deallocate variables
	if (allocated(YGRID_t)) then 
		deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )
	endif 

	! Output variable
	Tax_Reform_Welfare = Av_Util_NB 

	print*, 'Tax Reform: tauK=',tauK,'tauW=',tauW_at,'Gain',Tax_Reform_Welfare


end Function Tax_Reform_Welfare

!========================================================================================
!========================================================================================
!========================================================================================



end module GKK_Stats
