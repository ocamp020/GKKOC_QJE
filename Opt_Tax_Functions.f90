
Module Opt_Tax_Functions
	use parameters 
	use global
	use Opt_Tax_Parameters
	use programfunctions
	use Toolbox
	use omp_lib
	    
	Contains

!================================================================================
SUBROUTINE Find_Opt_Tax(switch,opt_Tau,a,b)
	IMPLICIT NONE
	logical , intent(in)  :: switch
	real(dp), intent(in)  :: a,b
	real(dp), intent(out) :: opt_Tau
	real(dp)              :: brentvaluet

	if (switch) then 
		brentvaluet = brent( a, (a+b)/2.0_dp , b , EQ_WELFARE_GIVEN_TauK, brent_tol, Opt_Tau)  
	else 
		brentvaluet = brent( a, (a+b)/2.0_dp , b , EQ_WELFARE_GIVEN_TauW, brent_tol, Opt_Tau)
	end if 


END SUBROUTINE Find_Opt_Tax

!================================================================================

FUNCTION EQ_WELFARE_GIVEN_TauK(tauk_in)
	IMPLICIT NONE 
	real(DP), intent(in) :: tauk_in
	real(DP) ::EQ_WELFARE_GIVEN_TauK

	tauK    = tauk_in
	tauW_at = 0.0_DP

	GBAR_exp = 0.0_DP
	DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.003 ) ! as long as the difference is greater than 0.05% continue
	    CALL FIND_DBN_EQ
	    CALL GOVNT_BUDGET_OPT
	    GBAR_exp = GBAR    
	ENDDO

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
		eta_K_exp = eta_K
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

	! CALL COMPUTE_WELFARE_GAIN
	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
	CALL Firm_Value
    	EQ_WELFARE_GIVEN_TAUK = - sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))

	!CALL COMPUTE_STATS

	!
	!OPEN   (UNIT=3, FILE='psi', STATUS='replace')
	!print*,'Budget balancing psi=',psi
	!WRITE(unit=3, FMT=*) psi
	!CLOSE (unit=3)
	!
	!print*,'tauK=', tauK, ' CE_NEWBORN=', CE_NEWBORN, 'Av. Util=',sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))
	print*,' ';print*,'------------------------------------------------------------------------'
	print*, 'eta_K=',eta_K,'tauK=', tauK,'Av_Utility_NEWBORN=', -EQ_WELFARE_GIVEN_TAUK, 'Psi=', psi
	print*,'------------------------------------------------------------------------';print*,' '

END  FUNCTION EQ_WELFARE_GIVEN_TAUK

!================================================================================

FUNCTION EQ_WELFARE_GIVEN_TauW(tauW_in)
	IMPLICIT NONE 
	real(DP), intent(in) :: tauW_in
	real(DP) ::EQ_WELFARE_GIVEN_TauW

	tauK = 0.0_DP
	tauW_at =  tauW_in

	GBAR_exp = 0.0_DP
	DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.01 ) ! as long as the difference is greater than 0.05% continue
	    CALL FIND_DBN_EQ
	    CALL GOVNT_BUDGET_OPT
	    GBAR_exp = GBAR    
	ENDDO

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
		eta_K_exp = eta_K
		eta_K_exp = eta_K
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

	!CALL COMPUTE_WELFARE_GAIN
	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
	CALL Firm_Value
    	EQ_WELFARE_GIVEN_TauW = - sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
    !CALL COMPUTE_STATS

	!OPEN   (UNIT=3, FILE='psi', STATUS='replace')
	!print*,'Budget balancing psi=',psi
	!WRITE(unit=3, FMT=*) psi
	!CLOSE (unit=3)
	!
	!print*,'tauW=', tauW, ' CE_NEWBORN=', CE_NEWBORN, 'Av. Util=',sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))
	print*, ' '
	print*, 'EQ_WELFARE_GIVEN_TauW Completed'
	print*, 'tauW_at=', tauW_at,'Av_Utility_NEWBORN=', -EQ_WELFARE_GIVEN_TAUW, 'Psi=', psi
	print*, ' '

END  FUNCTION EQ_WELFARE_GIVEN_TauW

!================================================================================

FUNCTION EQ_WELFARE_GIVEN_TauC(tauC_in,Opt_Tax_KW)
	IMPLICIT NONE 
	real(DP), intent(in) :: tauC_in
	logical , intent(in) :: Opt_Tax_KW ! Find capital taxes if true, and wealth taxes if false
	real(DP) :: EQ_WELFARE_GIVEN_TauC
	real(DP) :: brentvaluet, tau_indicator(1)

	tauC = tauC_in 

	! Get capital and wealth taxes that balance the budget
	if (Opt_Tax_KW) then 
		tauW_at = 0.0_dp 
		tau_indicator = 1.0_dp 
		! Solve model at current taxes 
		print*,'tauC=',tauC,'tauK=',tauK
		CALL FIND_DBN_EQ
	    CALL GOVNT_BUDGET(.false.)
	    GBAR_exp = GBAR  
	    print*, ' '; print*, 'tauC=',tauC,'tauK=',tauK,'GBAR_exp=',GBAR_exp,'GBAR_bench=',GBAR_bench; print*,' ' 
	    if (GBAR_exp>GBAR_bench) then 
		! Bracket budget balancing taxes 
		do while (GBAR_exp>GBAR_bench) 
			tauK = tauK - 0.05_dp
			print*,'tauC=',tauC,'tauK=',tauK
			CALL FIND_DBN_EQ
		    CALL GOVNT_BUDGET(.false.)
		    GBAR_exp = GBAR 
		    print*, ' '; print*, 'tauC=',tauC,'tauK=',tauK,'GBAR_exp=',GBAR_exp,'GBAR_bench=',GBAR_bench; print*,' ' 
		enddo 
		! Find budget balancing taxes	
		brentvaluet = brent_p( tauK , tauK + 0.025 , tauK + 0.05_dp , diff_GBAR , brent_tol, tauK , tau_indicator ) 
		else 
		! Bracket budget balancing taxes 
		do while (GBAR_exp<GBAR_bench) 
			tauK = tauK + 0.05_dp
			print*,'tauC=',tauC,'tauK=',tauK
			CALL FIND_DBN_EQ
		    CALL GOVNT_BUDGET(.false.)
		    GBAR_exp = GBAR 
		    print*, ' '; print*, 'tauC=',tauC,'tauK=',tauK,'GBAR_exp=',GBAR_exp,'GBAR_bench=',GBAR_bench; print*,' ' 
		enddo 
		! Find budget balancing taxes	
		brentvaluet = brent_p( tauK-0.05_dp , tauK - 0.025 , tauK , diff_GBAR , brent_tol, tauK , tau_indicator ) 
		endif 
	else 
		tauK = 0.0_dp 
		tau_indicator = 0.0_dp 
		! Solve model at current taxes 
		CALL FIND_DBN_EQ
	    CALL GOVNT_BUDGET(.false.)
	    GBAR_exp = GBAR  
	    print*, ' '; print*, 'tauC=',tauC,'tauK=',tauK,'GBAR_exp=',GBAR_exp,'GBAR_bench=',GBAR_bench; print*,' ' 
	    if (GBAR_exp>GBAR_bench) then 
		! Bracket budget balancing taxes 
		do while (GBAR_exp>GBAR_bench) 
			tauW_at = tauW_at - 0.001_dp
			CALL FIND_DBN_EQ
		    CALL GOVNT_BUDGET(.false.)
		    GBAR_exp = GBAR  
		    print*, ' '; print*, 'tauC=',tauC,'tauW=',tauW_at,'GBAR_exp=',GBAR_exp,'GBAR_bench=',GBAR_bench; print*,' ' 
		enddo 
		! Find budget balancing taxes	
		brentvaluet = brent_p( tauW_at , tauW_at + 0.0005 , tauW_at + 0.001_dp , diff_GBAR , brent_tol, tauK , tau_indicator ) 
		else 
		! Bracket budget balancing taxes 
		do while (GBAR_exp>GBAR_bench) 
			tauW_at = tauW_at + 0.001_dp
			CALL FIND_DBN_EQ
		    CALL GOVNT_BUDGET(.false.)
		    GBAR_exp = GBAR  
		    print*, ' '; print*, 'tauC=',tauC,'tauW=',tauW_at,'GBAR_exp=',GBAR_exp,'GBAR_bench=',GBAR_bench; print*,' ' 
		enddo 
		! Find budget balancing taxes	
		brentvaluet = brent_p( tauW_at - 0.001_dp , tauW_at - 0.0005 , tauW_at , diff_GBAR , brent_tol, tauK , tau_indicator ) 
		endif 
	endif 


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
		eta_K_exp = eta_K
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

	!CALL COMPUTE_WELFARE_GAIN
	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
	CALL Firm_Value
    	EQ_WELFARE_GIVEN_TauC = - sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
    !CALL COMPUTE_STATS

	!OPEN   (UNIT=3, FILE='psi', STATUS='replace')
	!print*,'Budget balancing psi=',psi
	!WRITE(unit=3, FMT=*) psi
	!CLOSE (unit=3)
	!
	!print*,'tauW=', tauW, ' CE_NEWBORN=', CE_NEWBORN, 'Av. Util=',sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))
	print*, ' '
	print*, 'tauC',tauC,'tauK',tauK,'tauW_at=', tauW_at, 'Psi=', psi,'Av_Utility_NEWBORN=', -EQ_WELFARE_GIVEN_TauC
	print*, ' '

END  FUNCTION EQ_WELFARE_GIVEN_TauC


!================================================================================

FUNCTION EQ_WELFARE_GIVEN_PSI(tauK_in,tauW_in,psi_in)
	IMPLICIT NONE 
	real(DP), intent(in) :: tauK_in, tauW_in, psi_in
	real(DP) ::EQ_WELFARE_GIVEN_PSI

	tauK    = tauK_in
	tauW_at = tauW_in
	psi     = psi_in

	GBAR_exp = 0.0_DP
	DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.001 ) ! as long as the difference is greater than 0.05% continue
	    CALL FIND_DBN_EQ
	    CALL GOVNT_BUDGET_OPT_TauC
	    GBAR_exp = GBAR    
	ENDDO

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
		eta_K_exp = eta_K
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

	!CALL COMPUTE_WELFARE_GAIN
	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
	CALL Firm_Value
    	EQ_WELFARE_GIVEN_PSI = - sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
    !CALL COMPUTE_STATS

	!OPEN   (UNIT=3, FILE='psi', STATUS='replace')
	!print*,'Budget balancing psi=',psi
	!WRITE(unit=3, FMT=*) psi
	!CLOSE (unit=3)
	!
	!print*,'tauW=', tauW, ' CE_NEWBORN=', CE_NEWBORN, 'Av. Util=',sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))
	print*, 'tauK=', tauK,'tauW_at=', tauW_at,'Psi=', psi,'Av_Utility_NEWBORN=', -EQ_WELFARE_GIVEN_PSI,'tauC=',tauC 

END  FUNCTION EQ_WELFARE_GIVEN_PSI


!================================================================================

FUNCTION EQ_WELFARE_GIVEN_etaK(eta_K_in)
	IMPLICIT NONE 
	real(DP), intent(in) :: eta_k_in
	real(DP) :: EQ_WELFARE_GIVEN_etaK

	eta_K   = eta_k_in
	tauW_at = 0.0_DP

	GBAR_exp = 0.0_DP
	DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.003 ) ! as long as the difference is greater than 0.05% continue
	    CALL FIND_DBN_EQ
	    CALL GOVNT_BUDGET_OPT_tauK
	    GBAR_exp = GBAR    
	ENDDO

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
		eta_K_exp = eta_K
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

	! CALL COMPUTE_WELFARE_GAIN
	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction,Bq_Value)
	EQ_WELFARE_GIVEN_etaK = - sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))

	!CALL COMPUTE_STATS

	!
	!OPEN   (UNIT=3, FILE='psi', STATUS='replace')
	!print*,'Budget balancing psi=',psi
	!WRITE(unit=3, FMT=*) psi
	!CLOSE (unit=3)
	!
	!print*,'tauK=', tauK, ' CE_NEWBORN=', CE_NEWBORN, 'Av. Util=',sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))
	print*,' ';print*,'------------------------------------------------------------------------'
	print*, 'eta_K=',eta_K,'tauK=', tauK,'Av_Utility_NEWBORN=', -EQ_WELFARE_GIVEN_etaK
	print*,'------------------------------------------------------------------------';print*,' '

END  FUNCTION EQ_WELFARE_GIVEN_etaK


!====================================================================

SUBROUTINE GOVNT_BUDGET_OPT_tauK()
	IMPLICIT NONE

	real(DP) :: GBAR_W,  GBAR_L, GBAR_C, GBAR_NK, BT_EARNINGS , A_EARNINGS, SSC_Payments
	real(DP) :: new_tauK

	GBAR        = 0.0_DP
	GBAR_K 		= 0.0_DP
	GBAR_W		= 0.0_DP
	GBAR_C 		= 0.0_DP
	GBAR_L 		= 0.0_DP
	GBAR_NK 	= 0.0_DP
	GBAR_BQ     = 0.0_DP
	BT_EARNINGS = 0.0_DP
	A_EARNINGS 	= 0.0_DP

	DO xi=1,nx
	DO age=1, MaxAge
	DO ai=1,na
	DO zi=1,nz
	DO lambdai=1,nlambda
	DO ei=1,ne
		GBAR = GBAR + DBN1(age,ai,zi,lambdai,ei,xi) * ( &
	    	  &   ( R*agrid(ai) + Pr_mat(ai,zi,xi) - (1.0_DP-tauK)*( R*agrid(ai) + Pr_mat(ai,zi,xi) )**(1.0_dp-eta_K) ) &
	          & + ( agrid(ai) + (1.0_DP-tauK)*( R*agrid(ai) + Pr_mat(ai,zi,xi) )**(1.0_dp-eta_K)  ) - YGRID(ai,zi,xi)   &
	          & + yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei,xi)  									&
	          & - psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei,xi))**(1.0_DP-tauPL)  			&
	          & + tauC * cons(age, ai, zi, lambdai,ei,xi)  												&   
	          & + tau_bq*aprime(age,ai,zi,lambdai,ei,xi)*(1.0_DP-survP(age)) )

	    GBAR_L = GBAR_L  + DBN1(age,ai,zi,lambdai,ei,xi) * (  yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei,xi) &
	          &- psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei,xi))**(1.0_DP-tauPL) )

	    BT_EARNINGS = BT_EARNINGS + DBN1(age,ai,zi,lambdai,ei,xi) * (R*agrid(ai) + Pr_mat(ai,zi,xi)) 
	    
	    A_EARNINGS  =  A_EARNINGS + DBN1(age,ai,zi,lambdai,ei,xi) * (R*agrid(ai) + Pr_mat(ai,zi,xi))**(1.0_dp-eta_K)
	    
	    GBAR_K = GBAR_K +DBN1(age,ai,zi,lambdai,ei,xi) * (  &
	    	  &   ( R*agrid(ai) + Pr_mat(ai,zi,xi) - (1.0_DP-tauK)*( R*agrid(ai) + Pr_mat(ai,zi,xi) )**(1.0_dp-eta_K) ) )

	    GBAR_W = GBAR_W +DBN1(age,ai,zi,lambdai,ei,xi) * (  &
	    	  &   ( agrid(ai) + (1.0_DP-tauK)*( R*agrid(ai) + Pr_mat(ai,zi,xi) )**(1.0_dp-eta_K)  ) - YGRID(ai,zi,xi)    )

      	GBAR_C = GBAR_C +  DBN1(age,ai,zi,lambdai,ei,xi) * tauC * cons(age, ai, zi, lambdai,ei,xi)
	    
	    GBAR_BQ = GBAR_BQ +  DBN1(age,ai,zi,lambdai,ei,xi) * tau_bq * aprime(age,ai,zi,lambdai,ei,xi) * (1.0_dp-survP(age))
	   
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	GBAR_NK = GBAR - GBAR_K 

	SSC_Payments = 0.0_DP

	DO xi=1,nx
	DO age=RetAge, MaxAge
	DO ai=1,na
	DO zi=1,nz
	DO lambdai=1,nlambda
	DO ei=1,ne

	    SSC_Payments = SSC_Payments + DBN1(age,ai,zi,lambdai,ei,xi) * RetY_lambda_e(lambdai,ei) 
	    
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	GBAR = GBAR -  SSC_Payments


	print*,' '; print*,'-----------------------------------------------------------------------------'
	print*, "Government Budget - Revenues and taxes"
	print*,' '
	print '(A,F8.5,X,A,F8.5)','GBAR_bench=',GBAR_bench, 'GBAR=',GBAR
	print '(A,F7.4,X,A,F7.4,X,A,F7.4)','	SSC_Payments=', SSC_Payments, 'GBAR_L=',GBAR_L,'Av. Labor Tax=',100.0_dp*GBAR_L/Ebar
	print '(A,F7.4,X,A,F7.4,X,A,F7.4,X,A,F7.4)','	GBAR_K=', GBAR_K, "GBAR_W=", GBAR_W, 'GBAR_C=', GBAR_C, 'GBAR_BQ=', GBAR_BQ
	print '(A,F7.4,X,A,F7.4,X,A,F7.4)','	GBAR_NK=',GBAR_NK, 'BT_EARNINGS=',BT_EARNINGS,'A_EARNINGS=',A_EARNINGS 
	print '(A,F7.4,X,X,A,F7.4,X,X,A,F7.4,X,X,A,F7.4,X,X,A,F7.4,X,X,A,F7.4,X,X,A,F7.4)',&
			&'	Tau_L=',100.0_dp*(1.0_dp-psi),'Tau_K=', 100.0_dp*tauK, 'Tau_W=', 100.0_dp*tauW_at,&
			& 'Tau_C=', 100.0_dp*tauC, 'Tau_BQ=', 100.0_dp*tau_bq, "Threshold", Y_a_threshold, "eta_K", eta_K
	print*,'-----------------------------------------------------------------------------'
	print*, ' '

	! OBTAIN NEW tauK IF GOVETNMENT BUDGET DOES NOT BALANCE
	if (solving_bench .eq. 0) then
	    IF (  abs(100.0_DP*(1.0_DP-GBAR/GBAR_bench)) .gt. 0.01 ) THEN
	        !new_psi =  ( BT_EARNINGS - GBAR_bench -  SSC_Payments   + GBAR_NL ) / A_EARNINGS
	        ! new_tauK =  1.0_dp - (  GBAR_bench  - GBAR_NK ) / BT_EARNINGS
	        new_tauK =  1.0_dp - ( BT_EARNINGS + GBAR_NK - GBAR_bench - SSC_Payments) / A_EARNINGS
	        PRINT*,'New tauK=',new_tauK,'Old tauK=',tauK,'G_gap=',100.0_DP*(1.0_DP-GBAR/GBAR_bench)
	        print*,'test GBAR_K:',BT_EARNINGS-A_EARNINGS*(1.0_dp-tauK),GBAR_K
	        print*,'test GBAR  :',BT_EARNINGS-A_EARNINGS*(1.0_dp-tauK)+GBAR_NK,GBAR+SSC_Payments
	        tauK = 0.5_dp*new_tauK+0.5_dp*tauK
	    ENDIF
	endif     



END SUBROUTINE GOVNT_BUDGET_OPT_tauK


!====================================================================

SUBROUTINE GOVNT_BUDGET_OPT()
	IMPLICIT NONE

	real(DP) :: GBAR_W,  GBAR_L, GBAR_C, GBAR_NL, BT_EARNINGS , A_EARNINGS, SSC_Payments
	real(DP) :: new_psi

	GBAR        = 0.0_DP
	GBAR_K 		= 0.0_DP
	GBAR_W		= 0.0_DP
	GBAR_C 		= 0.0_DP
	GBAR_L 		= 0.0_DP
	GBAR_NL 	= 0.0_DP
	GBAR_BQ     = 0.0_DP
	BT_EARNINGS = 0.0_DP
	A_EARNINGS 	= 0.0_DP

	DO xi=1,nx
	DO age=1, MaxAge
	DO ai=1,na
	DO zi=1,nz
	DO lambdai=1,nlambda
	DO ei=1,ne
		GBAR = GBAR + DBN1(age,ai,zi,lambdai,ei,xi) * ( &
	    	  &   ( R*agrid(ai) + Pr_mat(ai,zi,xi) - (1.0_DP-tauK)*( R*agrid(ai) + Pr_mat(ai,zi,xi) )**(1.0_dp-eta_K) ) &  
	          & + ( agrid(ai) + (1.0_DP-tauK)*( R*agrid(ai) + Pr_mat(ai,zi,xi) )**(1.0_dp-eta_K)  ) - YGRID(ai,zi,xi)   &
	          & + yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei,xi)  									&
	          & - psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei,xi))**(1.0_DP-tauPL)  			&
	          & + tauC * cons(age, ai, zi, lambdai,ei,xi)  												&   
	          & + tau_bq*aprime(age,ai,zi,lambdai,ei,xi)*(1.0_DP-survP(age)) )

	    GBAR_NL = GBAR_NL + DBN1(age,ai,zi,lambdai,ei,xi) * ( &
	    	  &   ( R*agrid(ai) + Pr_mat(ai,zi,xi) - (1.0_DP-tauK)*( R*agrid(ai) + Pr_mat(ai,zi,xi) )**(1.0_dp-eta_K) ) &  	
	          & + ( agrid(ai) + (1.0_DP-tauK)*( R*agrid(ai) + Pr_mat(ai,zi,xi) )**(1.0_dp-eta_K)  ) - YGRID(ai,zi,xi)   &
	          & + tauC * cons(age, ai, zi, lambdai,ei,xi)   												&   
	          & + tau_bq*aprime(age,ai,zi,lambdai,ei,xi)*(1.0_DP-survP(age)) )

	    GBAR_L = GBAR_L  + DBN1(age,ai,zi,lambdai,ei,xi) * (  yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei,xi) &
	          &- psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei,xi))**(1.0_DP-tauPL) )

	    BT_EARNINGS = BT_EARNINGS + DBN1(age,ai,zi,lambdai,ei,xi) * yh(age,lambdai,ei)* Hours(age, ai, zi, lambdai,ei,xi) 
	    
	    A_EARNINGS  = A_EARNINGS  + DBN1(age,ai,zi,lambdai,ei,xi) *&
	                    & (yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei,xi))**(1.0_DP-tauPL)
	    
	    GBAR_K = GBAR_K +DBN1(age,ai,zi,lambdai,ei,xi) * (  &
	    	  &   ( R*agrid(ai) + Pr_mat(ai,zi,xi) - (1.0_DP-tauK)*( R*agrid(ai) + Pr_mat(ai,zi,xi) )**(1.0_dp-eta_K) ) &  	 
	          & + ( agrid(ai) + (1.0_DP-tauK)*( R*agrid(ai) + Pr_mat(ai,zi,xi) )**(1.0_dp-eta_K)  ) - YGRID(ai,zi,xi)    )

      	GBAR_C = GBAR_C +  DBN1(age,ai,zi,lambdai,ei,xi) * tauC * cons(age, ai, zi, lambdai,ei,xi)
	    
	    GBAR_BQ = GBAR_BQ +  DBN1(age,ai,zi,lambdai,ei,xi) * tau_bq * aprime(age,ai,zi,lambdai,ei,xi) * (1.0_dp-survP(age))
	   
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	SSC_Payments = 0.0_DP

	DO xi=1,nx
	DO age=RetAge, MaxAge
	DO ai=1,na
	DO zi=1,nz
	DO lambdai=1,nlambda
	DO ei=1,ne

	    SSC_Payments = SSC_Payments + DBN1(age,ai,zi,lambdai,ei,xi) * RetY_lambda_e(lambdai,ei) 
	    
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	GBAR = GBAR -  SSC_Payments


	print*,' '; print*,'-----------------------------------------------------------------------------'
	print*, "Government Budget - Revenues and taxes"
	print*,' '
	print '(A,F8.5,X,A,F8.5)','GBAR_bench=',GBAR_bench, 'GBAR=',GBAR
	print '(A,F7.4,X,A,F7.4,X,A,F7.4)','	SSC_Payments=', SSC_Payments, 'GBAR_L=',GBAR_L,'Av. Labor Tax=',100.0_dp*GBAR_L/Ebar
	print '(A,F7.4,X,A,F7.4,X,A,F7.4,X,A,F7.4)','	GBAR_K=', GBAR_K, "GBAR_W=", GBAR_W, 'GBAR_C=', GBAR_C, 'GBAR_BQ=', GBAR_BQ
	print '(A,F7.4,X,A,F7.4,X,A,F7.4)','	GBAR_NL=',GBAR_NL, 'BT_EARNINGS=',BT_EARNINGS,'A_EARNINGS=',A_EARNINGS 
	print '(A,F7.4,X,X,A,F7.4,X,X,A,F7.4,X,X,A,F7.4,X,X,A,F7.4,X,X,A,F7.4)',&
			&'	Tau_L=',100.0_dp*(1.0_dp-psi),'Tau_K=', 100.0_dp*tauK, 'Tau_W=', 100.0_dp*tauW_at,&
			& 'Tau_C=', 100.0_dp*tauC, 'Tau_BQ=', 100.0_dp*tau_bq, "Threshold", Y_a_threshold
	print*,'-----------------------------------------------------------------------------'
	print*, ' '

	! OBTAIN NEW PSI IF GOVETNMENT BUDGET DOES NOT BALANCE
	if (solving_bench .eq. 0) then
	    IF (  abs(100.0_DP*(1.0_DP-GBAR/GBAR_bench)) .gt. 0.001 ) THEN
	        !new_psi =  ( BT_EARNINGS - GBAR_bench -  SSC_Payments   + GBAR_NL ) / A_EARNINGS
	        new_psi = 1.0_dp - (  GBAR_bench + SSC_Payments - GBAR_NL ) / BT_EARNINGS
	        PRINT*,'NEW PSI=',new_psi,'Old Psi=',psi
	        psi = 0.5_dp*new_psi+0.5_dp*psi
	    ENDIF
	endif     



END SUBROUTINE GOVNT_BUDGET_OPT


!====================================================================

SUBROUTINE GOVNT_BUDGET_OPT_TauC()
	IMPLICIT NONE

	real(DP) :: GBAR_W,  GBAR_L, GBAR_C, GBAR_NL, BT_EARNINGS , A_EARNINGS, SSC_Payments, C_bar
	real(DP) :: new_tauC, aux_var 

	GBAR        = 0.0_DP
	GBAR_K 		= 0.0_DP
	GBAR_W		= 0.0_DP
	GBAR_C 		= 0.0_DP
	GBAR_L 		= 0.0_DP
	GBAR_NL 	= 0.0_DP
	GBAR_BQ     = 0.0_DP
	BT_EARNINGS = 0.0_DP
	A_EARNINGS 	= 0.0_DP
	C_bar       = 0.0_DP

	DO xi=1,nx
	DO age=1, MaxAge
	DO ai=1,na
	DO zi=1,nz
	DO lambdai=1,nlambda
	DO ei=1,ne
		GBAR = GBAR + DBN1(age,ai,zi,lambdai,ei,xi) * ( &
	    	  &   ( R*agrid(ai) + Pr_mat(ai,zi,xi) - (1.0_DP-tauK)*( R*agrid(ai) + Pr_mat(ai,zi,xi) )**(1.0_dp-eta_K) ) &  
	          & + ( agrid(ai) + (1.0_DP-tauK)*( R*agrid(ai) + Pr_mat(ai,zi,xi) )**(1.0_dp-eta_K)  ) - YGRID(ai,zi,xi)   &	
	          & + yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei,xi)  									&
	          & - psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei,xi))**(1.0_DP-tauPL)  			&
	          & + tauC * cons(age, ai, zi, lambdai,ei,xi)  										    & 
	          & + tau_bq*aprime(age,ai,zi,lambdai,ei,xi)*(1.0_DP-survP(age)) )

	    GBAR_NL = GBAR_NL + DBN1(age,ai,zi,lambdai,ei,xi) * ( &
	    	  &   ( R*agrid(ai) + Pr_mat(ai,zi,xi) - (1.0_DP-tauK)*( R*agrid(ai) + Pr_mat(ai,zi,xi) )**(1.0_dp-eta_K) ) &
	          & + ( agrid(ai) + (1.0_DP-tauK)*( R*agrid(ai) + Pr_mat(ai,zi,xi) )**(1.0_dp-eta_K)  ) - YGRID(ai,zi,xi)   & 
	          & + tauC * cons(age, ai, zi, lambdai,ei,xi)  &   
	          & + tau_bq*aprime(age,ai,zi,lambdai,ei,xi)*(1.0_DP-survP(age)) )         

	    GBAR_L = GBAR_L  + DBN1(age,ai,zi,lambdai,ei,xi) * (  yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei,xi) &
	          &- psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei,xi))**(1.0_DP-tauPL) )

	    BT_EARNINGS = BT_EARNINGS + DBN1(age,ai,zi,lambdai,ei,xi) * yh(age,lambdai,ei)* Hours(age, ai, zi, lambdai,ei,xi) 
	    
	    A_EARNINGS  = A_EARNINGS  + DBN1(age,ai,zi,lambdai,ei,xi) *&
	                    & (yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei,xi))**(1.0_DP-tauPL)
	    
	    GBAR_K = GBAR_K +DBN1(age,ai,zi,lambdai,ei,xi) * (  &
	    	  &   ( R*agrid(ai) + Pr_mat(ai,zi,xi) - (1.0_DP-tauK)*( R*agrid(ai) + Pr_mat(ai,zi,xi) )**(1.0_dp-eta_K) ) & 
	          & + ( agrid(ai) + (1.0_DP-tauK)*( R*agrid(ai) + Pr_mat(ai,zi,xi) )**(1.0_dp-eta_K)  ) - YGRID(ai,zi,xi) )

      	GBAR_C = GBAR_C +  DBN1(age,ai,zi,lambdai,ei,xi) * tauC * cons(age, ai, zi, lambdai,ei,xi)

      	C_bar = C_bar +  DBN1(age,ai,zi,lambdai,ei,xi) * cons(age, ai, zi, lambdai,ei,xi)

      	GBAR_BQ = GBAR_BQ +  DBN1(age,ai,zi,lambdai,ei,xi) * tau_bq * aprime(age,ai,zi,lambdai,ei,xi) * (1.0_dp-survP(age))
	    
	   
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	SSC_Payments = 0.0_DP

	DO xi=1,nx
	DO age=RetAge, MaxAge
	DO ai=1,na
	DO zi=1,nz
	DO lambdai=1,nlambda
	DO ei=1,ne

	    SSC_Payments = SSC_Payments + DBN1(age,ai,zi,lambdai,ei,xi) * RetY_lambda_e(lambdai,ei) 
	    
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	GBAR = GBAR -  SSC_Payments

	print*,' '
	print*,'Results from GOVNT_BUDGET'
	print*, 'GBAR_bench',GBAR_bench, 'GBAR=',GBAR, 'SSC_Payments=', SSC_Payments, 'GBAR_L=',GBAR_L,'Av. Labor Tax=',GBAR_L/Ebar 
	print*, 'GBAR_NL  =',GBAR_NL, 'BT_EARNINGS=',BT_EARNINGS,'A_EARNINGS=',A_EARNINGS 
	PRINT*,'PSI=',psi, 'tauK=', tauK, 'tauW_at=',tauW_at, 'tauC', tauC

	! OBTAIN NEW Tau_C IF GOVETNMENT BUDGET DOES NOT BALANCE
	! Formula exploits that consumption expenditure is constant
	if (solving_bench .eq. 0) then
	    IF (  abs(100.0_DP*(1.0_DP-GBAR/GBAR_bench)) .gt. 0.001 ) THEN
	    	aux_var  = ( GBAR_bench - (GBAR - GBAR_C) )/(C_bar+GBAR_C)
	    	new_tauC = aux_var/(1.0_dp-aux_var)
	        PRINT*,'NEW Tau_C=',new_tauC,'Old Tau_C=',tauC,'aux_var=',aux_var,'Rev=',GBAR - GBAR_C
	        tauC = new_tauC
	        ! tauC = 0.5_dp*new_tauC+0.5_dp*tauC
	    ENDIF
	endif     



END SUBROUTINE GOVNT_BUDGET_OPT_TauC


!================================================================================


Function diff_GBAR(tau_in,tau_indicator)
	implicit none 
	real(dp), intent(in) :: tau_in
	real(dp), dimension(:), intent(in) :: tau_indicator 
	real(dp) :: diff_GBAR

	! Set taxes
	if (tau_indicator(1).eq.1.0_dp) then 
		tauK = tau_in 
	else 
		tauW_at = tau_in
	endif 

	! Solve the model
	CALL FIND_DBN_EQ
    CALL GOVNT_BUDGET(.false.)
    GBAR_exp = GBAR  

    ! Get GBAR difference
    diff_GBAR = (GBAR_exp-gBAR_bench)**2 

    print*, ' '; print*, 'tauC=',tauC,'tau_in=',tau_in,'GBAR_exp=',GBAR_exp,'GBAR_bench=',GBAR_bench; print*,' '


end Function diff_GBAR


!================================================================================



end Module Opt_Tax_Functions