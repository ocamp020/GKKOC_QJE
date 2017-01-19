
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
	DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.001 ) ! as long as the difference is greater than 0.05% continue
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
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

	! CALL COMPUTE_WELFARE_GAIN
	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction)
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
	print*, 'tauK=', tauK,'Av_Utility_NEWBORN=', -EQ_WELFARE_GIVEN_TAUK, 'Psi=', psi

END  FUNCTION EQ_WELFARE_GIVEN_TAUK

!================================================================================

FUNCTION EQ_WELFARE_GIVEN_TauW(tauW_in)
	IMPLICIT NONE 
	real(DP), intent(in) :: tauW_in
	real(DP) ::EQ_WELFARE_GIVEN_TauW

	tauK = 0.0_DP
	tauW_at =  tauW_in

	GBAR_exp = 0.0_DP
	DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.001 ) ! as long as the difference is greater than 0.05% continue
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
		tauPL_exp = tauPL
		psi_exp   = psi
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

	!CALL COMPUTE_WELFARE_GAIN
	CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction)
	CALL Firm_Value
    	EQ_WELFARE_GIVEN_TauW = - sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
    !CALL COMPUTE_STATS

	!OPEN   (UNIT=3, FILE='psi', STATUS='replace')
	!print*,'Budget balancing psi=',psi
	!WRITE(unit=3, FMT=*) psi
	!CLOSE (unit=3)
	!
	!print*,'tauW=', tauW, ' CE_NEWBORN=', CE_NEWBORN, 'Av. Util=',sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))
	print*, 'tauW_at=', tauW_at,'Av_Utility_NEWBORN=', -EQ_WELFARE_GIVEN_TAUW, 'Psi=', psi

END  FUNCTION EQ_WELFARE_GIVEN_TauW


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
	BT_EARNINGS = 0.0_DP
	A_EARNINGS 	= 0.0_DP

	DO xi=1,nx
	DO age=1, MaxAge
	DO ai=1,na
	DO zi=1,nz
	DO lambdai=1,nlambda
	DO ei=1,ne
		GBAR = GBAR + DBN1(age,ai,zi,lambdai,ei,xi) * ( tauK*( R*agrid(ai) + Pr_mat(ai,zi,xi) )  	    &
	          & + ( agrid(ai) + ( R*agrid(ai) + Pr_mat(ai,zi,xi) ) *(1.0_DP-tauK)  ) - YGRID(ai,zi,xi) 	&	
	          & + yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei,xi)  									&
	          & - psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei,xi))**(1.0_DP-tauPL)  			&
	          & + tauC * cons(age, ai, zi, lambdai,ei,xi)  )   

	    GBAR_NL = GBAR_NL + DBN1(age,ai,zi,lambdai,ei,xi) * ( tauK*( R*agrid(ai) + Pr_mat(ai,zi,xi) )   &
	          & + ( agrid(ai) + ( R*agrid(ai) + Pr_mat(ai,zi,xi) ) *(1.0_DP-tauK)  ) - YGRID(ai,zi,xi)  & 
	          & + tauC * cons(age, ai, zi, lambdai,ei,xi)  )         

	    GBAR_L = GBAR_L  + DBN1(age,ai,zi,lambdai,ei,xi) * (  yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei,xi) &
	          &- psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei,xi))**(1.0_DP-tauPL) )

	    BT_EARNINGS = BT_EARNINGS + DBN1(age,ai,zi,lambdai,ei,xi) * yh(age,lambdai,ei)* Hours(age, ai, zi, lambdai,ei,xi) 
	    
	    A_EARNINGS  = A_EARNINGS  + DBN1(age,ai,zi,lambdai,ei,xi) *&
	                    & (yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei,xi))**(1.0_DP-tauPL)
	    
	    GBAR_K = GBAR_K +DBN1(age,ai,zi,lambdai,ei,xi) * (  tauK*( R*agrid(ai) + Pr_mat(ai,zi,xi) )   &
	          & + ( agrid(ai) + ( R*agrid(ai) + Pr_mat(ai,zi,xi) ) *(1.0_DP-tauK)  ) - YGRID(ai,zi,xi) )

      	GBAR_C = GBAR_C +  DBN1(age,ai,zi,lambdai,ei,xi) * tauC * cons(age, ai, zi, lambdai,ei,xi)
	    
	   
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



!================================================================================





end Module Opt_Tax_Functions