
Module Opt_Tax_Functions
	use parameters 
	use global
	use Opt_Tax_Parameters
	use programfunctions
	use Toolbox
	    
	Contains

!================================================================================
SUBROUTINE Find_Opt_Tax(switch,opt_Tau)
	IMPLICIT NONE
	integer , intent(in)  :: switch
	real(dp), intent(out) :: opt_Tau
	real(dp)              :: brentvaluet

	if (switch.eq.1) then 
		brentvaluet = brent(0.00_DP, 0.1_DP , 0.4_DP, EQ_WELFARE_GIVEN_TauK, brent_tol, Opt_Tau)  
	else 
		brentvaluet = brent(0.00_DP, 0.016_DP , 0.03_DP, EQ_WELFARE_GIVEN_TauW, brent_tol, Opt_Tau)
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
	DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.01 ) ! as long as the difference is greater than 0.05% continue
	    CALL FIND_DBN_EQ
	    CALL GOVNT_BUDGET_OPT
	    GBAR_exp = GBAR    
	ENDDO

	GBAR_exp = GBAR
	QBAR_exp = QBAR 
	NBAR_exp = NBAR  
	Y_exp 	 = YBAR
	Ebar_exp = EBAR
	rr_exp   = rr
	wage_exp = wage
	tauW_at_exp = tauW_at
	tauK_exp    = tauK
	tauPL_exp   = tauPL

	! CALL COMPUTE_WELFARE_GAIN
	CALL COMPUTE_VALUE_FUNCTION_SPLINE 
	EQ_WELFARE_GIVEN_TAUK = - sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))
	!CALL COMPUTE_STATS

	!
	!OPEN   (UNIT=3, FILE='psi', STATUS='replace')
	!print*,'Budget balancing psi=',psi
	!WRITE(unit=3, FMT=*) psi
	!CLOSE (unit=3)
	!
	!print*,'tauK=', tauK, ' CE_NEWBORN=', CE_NEWBORN, 'Av. Util=',sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))

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

	GBAR_exp = GBAR
	QBAR_exp = QBAR 
	NBAR_exp = NBAR  
	Y_exp    = YBAR
	Ebar_exp = EBAR
	rr_exp   = rr
	wage_exp = wage
	tauW_at_exp = tauW_at
	tauK_exp    = tauK
	tauPL_exp   = tauPL

	!CALL COMPUTE_WELFARE_GAIN
	CALL COMPUTE_VALUE_FUNCTION_SPLINE 
	EQ_WELFARE_GIVEN_TauW = - sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))
	!CALL COMPUTE_STATS

	!OPEN   (UNIT=3, FILE='psi', STATUS='replace')
	!print*,'Budget balancing psi=',psi
	!WRITE(unit=3, FMT=*) psi
	!CLOSE (unit=3)
	!
	!print*,'tauW=', tauW, ' CE_NEWBORN=', CE_NEWBORN, 'Av. Util=',sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))


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

	DO age=1, MaxAge
	DO ai=1,na
	DO zi=1,nz
	DO lambdai=1,nlambda
	DO ei=1,ne
	    GBAR_NL = GBAR_NL + DBN1(age,ai,zi,lambdai,ei) * ( tauK*( rr*(agrid(ai)*zgrid(zi))**mu-DepRate*agrid(ai) )  &
	          & + (agrid(ai) + ( rr * (zgrid(zi) * agrid(ai) )**mu - DepRate*agrid(ai) ) *(1.0_DP-tauK)  )		&
	          & - Y_a(agrid(ai),zgrid(zi)) 																		&	
	          & + tauC * cons(age, ai, zi, lambdai,ei)  )         

	    GBAR_L = GBAR_L  + DBN1(age,ai,zi,lambdai,ei) * (  yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei) &
	          &- psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei))**(1.0_DP-tauPL) )

	    BT_EARNINGS = BT_EARNINGS + DBN1(age,ai,zi,lambdai,ei) *     yh(age,lambdai,ei)* Hours(age, ai, zi, lambdai,ei) 
	    A_EARNINGS  = A_EARNINGS  + DBN1(age,ai,zi,lambdai,ei) *&
	                    & (yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei))**(1.0_DP-tauPL)
	    
	    GBAR_K = GBAR_K +DBN1(age,ai,zi,lambdai,ei) * ( tauK*( rr*(agrid(ai)*zgrid(zi))**mu-DepRate*agrid(ai) ) &
	            & + (agrid(ai) + ( rr * (zgrid(zi) * agrid(ai) )**mu - DepRate*agrid(ai) ) *(1.0_DP-tauK) )		&
	          	& - Y_a(agrid(ai),zgrid(zi)) )

      	GBAR_C = GBAR_C +  DBN1(age,ai,zi,lambdai,ei) * tauC * cons(age, ai, zi, lambdai,ei)
	    
	   
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	GBAR = GBAR_L + GBAR_NL


	SSC_Payments = 0.0_DP

	DO age=RetAge, MaxAge
	DO ai=1,na
	DO zi=1,nz
	DO lambdai=1,nlambda
	DO ei=1,ne

	    SSC_Payments = SSC_Payments + DBN1(age,ai,zi,lambdai,ei) * RetY_lambda_e(lambdai,ei) 
	    
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
	PRINT*,'PSI=',psi

	! OBTAIN NEW PSI IF GOVETNMENT BUDGET DOES NOT BALANCE
	if (solving_bench .eq. 0) then
	    IF (  abs(100.0_DP*(1.0_DP-GBAR/GBAR_bench)) .gt. 0.01 ) THEN
	        new_psi =  ( BT_EARNINGS - GBAR_bench -  SSC_Payments   + GBAR_NL ) / A_EARNINGS
	        PRINT*,'NEW PSI=',new_psi
	        psi = new_psi
	    ENDIF
	endif     



END SUBROUTINE GOVNT_BUDGET_OPT



!================================================================================





end Module Opt_Tax_Functions