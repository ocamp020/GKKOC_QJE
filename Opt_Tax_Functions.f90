
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
	tauW    = 0.0_DP

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
		tauw_exp  = tauW

	! CALL COMPUTE_WELFARE_GAIN
	CALL COMPUTE_VALUE_FUNCTION_SPLINE 
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
	tauW =  tauW_in

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
		tauw_exp  = tauW

	!CALL COMPUTE_WELFARE_GAIN
	CALL COMPUTE_VALUE_FUNCTION_SPLINE 
    	EQ_WELFARE_GIVEN_TauW = - sum(ValueFunction(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
    !CALL COMPUTE_STATS

	!OPEN   (UNIT=3, FILE='psi', STATUS='replace')
	!print*,'Budget balancing psi=',psi
	!WRITE(unit=3, FMT=*) psi
	!CLOSE (unit=3)
	!
	!print*,'tauW=', tauW, ' CE_NEWBORN=', CE_NEWBORN, 'Av. Util=',sum(ValueFunction(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))
	print*, 'tauW=', tauW,'Av_Utility_NEWBORN=', -EQ_WELFARE_GIVEN_TAUW, 'Psi=', psi

END  FUNCTION EQ_WELFARE_GIVEN_TauW


!====================================================================

SUBROUTINE GOVNT_BUDGET_OPT()
	IMPLICIT NONE

	real(DP) :: GBAR_W,  GBAR_L, GBAR_C, GBAR_NL, BT_EARNINGS , A_EARNINGS, SSC_Payments
	real(DP) :: new_psi, ap, kp, k_nb, Yp, y_nb, PrAprimelo, PrAprimehi
	integer  :: xi, age, yi, zi, lambdai, ei, xp_ind, ep_ind, zp_ind, lp_ind, tklo, tkhi

	GBAR        = 0.0_DP
	GBAR_K 		= 0.0_DP
	GBAR_W		= 0.0_DP
	GBAR_C 		= 0.0_DP
	GBAR_L 		= 0.0_DP
	GBAR_NL 	= 0.0_DP
	BT_EARNINGS = 0.0_DP
	A_EARNINGS 	= 0.0_DP



	DO xi=1,nx
	DO age=1, MaxAge-1
	DO yi=1,ny
	DO zi=1,nz
	DO lambdai=1,nlambda
	DO ei=1,ne
		ap = Aprime(age,yi,zi,lambdai,ei,xi)
		kp = Kprime(age,yi,zi,lambdai,ei,xi)

	do xp_ind = 1,nx
		Yp = Y_a(ap,kp,zi,xp_ind)
		if ( Yp .ge. ymax) then
		    tklo =ny-1
		else if (Yp .lt. ymin) then
		    tklo = 1
		else
		    tklo = ((Yp - ymin)/(ymax-ymin))**(1.0_DP/y_theta)*(ny-1)+1          
		end if            
		tkhi = tklo + 1        
		PrAprimelo = ( ygrid(tkhi) - Yp ) / ( ygrid(tkhi) -ygrid(tklo) )
		PrAprimehi = ( Yp - ygrid(tklo) ) / ( ygrid(tkhi) -ygrid(tklo) )        
		PrAprimelo = min(PrAprimelo, 1.0_DP); PrAprimelo = max(PrAprimelo, 0.0_DP)
		PrAprimehi = min(PrAprimehi, 1.0_DP); PrAprimehi = max(PrAprimehi, 0.0_DP)  

			GBAR_K = GBAR_K + pr_x(xi,xp_ind,zi,age)*survP(age)*DBN1(age,yi,zi,lambdai,ei,xi) * ( tauK* &
						&   ( R*ap + P*(xz_grid(xp_ind,zi)*kp)**mu - (R+DepRate)*kp ) + &
						&   (( ap + ( R*ap + P*(xz_grid(xp_ind,zi)*kp)**mu - (R+DepRate)*kp ) *(1.0_DP-tauK)  ) - Yp )  ) 


		if (age.lt.RetAge) then 
			GBAR_C = GBAR_C +  pr_x(xi,xp_ind,zi,age)*survP(age)*DBN1(age,yi,zi,lambdai,ei,xi) * tauC *&
			&  sum(pr_e(ei,:)*(PrAprimelo*Cons(age+1,tklo,zi,lambdai,:,xp_ind) + PrAprimehi*Cons(age+1,tkhi,zi,lambdai,:,xp_ind) ) ) 

			GBAR_L = GBAR_L  + pr_x(xi,xp_ind,zi,age)*survP(age)*DBN1(age,yi,zi,lambdai,ei,xi) * &
				&  sum(pr_e(ei,:)* ( PrAprimelo*(yh(age+1,lambdai,:)*Hours(age+1,tklo,zi,lambdai,:,xp_ind) &
           		&                     - psi*(yh(age+1,lambdai,:)*Hours(age+1,tklo,zi,lambdai,:,xp_ind))**(1.0_DP-tauPL) )  &
				& 				   + PrAprimehi*(yh(age+1,lambdai,:)*Hours(age+1,tkhi,zi,lambdai,:,xp_ind) &
            	&                     - psi*(yh(age+1,lambdai,:)*Hours(age+1,tkhi,zi,lambdai,:,xp_ind))**(1.0_DP-tauPL) )  ) )

            BT_EARNINGS = BT_EARNINGS +  pr_x(xi,xp_ind,zi,age)*survP(age)*DBN1(age,yi,zi,lambdai,ei,xi) * &
				&  sum(pr_e(ei,:)* ( PrAprimelo*(yh(age+1,lambdai,:)*Hours(age+1,tklo,zi,lambdai,:,xp_ind) )  &
				& 				   + PrAprimehi*(yh(age+1,lambdai,:)*Hours(age+1,tkhi,zi,lambdai,:,xp_ind) )  ) )

			A_EARNINGS  = A_EARNINGS  + pr_x(xi,xp_ind,zi,age)*survP(age)*DBN1(age,yi,zi,lambdai,ei,xi) * &
				&  sum(pr_e(ei,:)* ( PrAprimelo*((yh(age+1,lambdai,:)*Hours(age+1,tklo,zi,lambdai,:,xp_ind))**(1.0_DP-tauPL) )  &
				& 				   + PrAprimehi*((yh(age+1,lambdai,:)*Hours(age+1,tkhi,zi,lambdai,:,xp_ind))**(1.0_DP-tauPL) )  ) )

		else 
			GBAR_C = GBAR_C +  pr_x(xi,xp_ind,zi,age)*survP(age)*DBN1(age,yi,zi,lambdai,ei,xi) * tauC *&
					&  (PrAprimelo*Cons(age+1,tklo,zi,lambdai,ei,xp_ind) + PrAprimehi*Cons(age+1,tkhi,zi,lambdai,ei,xp_ind) )
		endif 
	enddo

	! New Borns
		do zp_ind=1,nz
		do lp_ind=1,nlambda
			k_nb = Linear_Int(agrid, opt_K_nb(:,zp_ind,lp_ind),na, Ap)    
			k_nb = max(theta(zp_ind)*ap , k_nb ) 
			
			do xp_ind=1,nx  
			
			Y_nb = Y_a(ap,k_nb,zp_ind,xp_ind)
			if ( Y_nb .ge. ymax) then
		    	tklo =ny-1
			else if (Y_nb .lt. ymin) then
			    tklo = 1
			else
			    tklo = ((Y_nb - ymin)/(ymax-ymin))**(1.0_DP/y_theta)*(ny-1)+1          
			end if            
			tkhi = tklo + 1        
			PrAprimelo = ( ygrid(tkhi) - Y_nb ) / ( ygrid(tkhi) -ygrid(tklo) )
			PrAprimehi = ( Y_nb - ygrid(tklo) ) / ( ygrid(tkhi) -ygrid(tklo) )        
			PrAprimelo = min(PrAprimelo, 1.0_DP); PrAprimelo = max(PrAprimelo, 0.0_DP)
			PrAprimehi = min(PrAprimehi, 1.0_DP); PrAprimehi = max(PrAprimehi, 0.0_DP)

			GBAR_K = GBAR_K + pr_x_nb(xp_ind,zp_ind)*pr_lambda(lambdai,lp_ind)*pr_z(zi,zp_ind)*(1.0_DP-survP(age))*&
						& DBN1(age,yi,zi,lambdai,ei,xi) * (tauK*( R*ap + P*(xz_grid(xp_ind,zp_ind)*k_nb)**mu - (R+DepRate)*k_nb ) + &
						& (( ap + ( R*ap + P*(xz_grid(xp_ind,zp_ind)*k_nb)**mu - (R+DepRate)*k_nb ) *(1.0_DP-tauK)  ) - Y_nb ) )

			GBAR_C = GBAR_C + pr_x_nb(xp_ind,zp_ind)*pr_lambda(lambdai,lp_ind)*pr_z(zi,zp_ind)*(1.0_DP-survP(age))*&
					&  DBN1(age,yi,zi,lambdai,ei,xi) * tauC *&
					&  (PrAprimelo*Cons(1,tklo,zp_ind,lp_ind,ne/2+1,xp_ind) + PrAprimehi*Cons(1,tkhi,zp_ind,lp_ind,ne/2+1,xp_ind) ) 

			GBAR_L = GBAR_L  + pr_x_nb(xp_ind,zp_ind)*pr_lambda(lambdai,lp_ind)*pr_z(zi,zp_ind)*(1.0_DP-survP(age))*&
				&  DBN1(age,yi,zi,lambdai,ei,xi) * &
				&  ( PrAprimelo*(yh(1,lp_ind,ne/2+1)*Hours(1,tklo,zp_ind,lp_ind,ne/2+1,xp_ind) &
           		&                 - psi*(yh(1,lp_ind,ne/2+1)*Hours(1,tklo,zp_ind,lp_ind,ne/2+1,xp_ind))**(1.0_DP-tauPL) )  &
				&  + PrAprimehi*(yh(1,lp_ind,ne/2+1)*Hours(age+1,tkhi,zi,lambdai,ne/2+1,xp_ind) &
            	&                 - psi*(yh(1,lp_ind,ne/2+1)*Hours(1,tkhi,zp_ind,lp_ind,ne/2+1,xp_ind))**(1.0_DP-tauPL) )  ) 

            BT_EARNINGS = BT_EARNINGS + pr_x_nb(xp_ind,zp_ind)*pr_lambda(lambdai,lp_ind)*pr_z(zi,zp_ind)*(1.0_DP-survP(age))*&
				&  DBN1(age,yi,zi,lambdai,ei,xi) * &
				&  ( PrAprimelo*(yh(1,lp_ind,ne/2+1)*Hours(1,tklo,zp_ind,lp_ind,ne/2+1,xp_ind)  )  &
				&  + PrAprimehi*(yh(1,lp_ind,ne/2+1)*Hours(age+1,tkhi,zi,lambdai,ne/2+1,xp_ind) )  ) 
	    
	   		A_EARNINGS  = A_EARNINGS  + pr_x_nb(xp_ind,zp_ind)*pr_lambda(lambdai,lp_ind)*pr_z(zi,zp_ind)*(1.0_DP-survP(age))*&
				&  DBN1(age,yi,zi,lambdai,ei,xi) * &
				&  ( PrAprimelo*((yh(1,lp_ind,ne/2+1)*Hours(1,tklo,zp_ind,lp_ind,ne/2+1,xp_ind))**(1.0_DP-tauPL) )  &
				&  + PrAprimehi*((yh(1,lp_ind,ne/2+1)*Hours(1,tkhi,zp_ind,lp_ind,ne/2+1,xp_ind))**(1.0_DP-tauPL) )  ) 
        	enddo 
        enddo 
		enddo  
	    
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	GBAR_NL = GBAR_K  + GBAR_C 
	GBAR    = GBAR_NL + GBAR_L

	SSC_Payments = 0.0_DP

	DO xi=1,nx
	DO age=RetAge, MaxAge
	DO yi=1,ny
	DO zi=1,nz
	DO lambdai=1,nlambda
	DO ei=1,ne

	    SSC_Payments = SSC_Payments + DBN1(age,yi,zi,lambdai,ei,xi) * RetY_lambda_e(lambdai,ei) 
	    
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
	PRINT*,'PSI=',psi, 'tauK=', tauK, 'tauW=',tauW

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