
Module programfunctions
	use parameters
	use global
	use Toolbox
	    
	Contains
		! Asset_Grid_Threshold
		
		! Y_h
		! MB_h
		! Y_a
		! MB_a
		! MB_a_bt
		! MB_a_at
		
		! FOC_R
		! FOC_WH
		! FOC_WH_NSU
		! FOC_H
		! FOC_H_NSU
		! FOC_HA

!========================================================================================
!========================================================================================

Subroutine Asset_Grid_Threshold(Y_a_threshold_in,agrid_t,na_t)
	IMPLICIT NONE
	real(dp), intent(in)       :: Y_a_threshold_in
	integer , intent(out)      :: na_t
	real(dp), dimension(:), allocatable, intent(out) :: agrid_t
	real(dp), dimension(nz*nx) :: a_aux
	integer                    :: a_ind
	integer , dimension(:), allocatable :: agrid_t_ind
	real(dp), dimension(:), allocatable :: par
	!real(dp), dimension(2)  :: par
	real(dp)                   :: max_wealth, K

	allocate( par(3) )
	a_ind = 0
	! If the threshold for wealth taxes is positive then agrid is adjusted
	if (Y_a_threshold_in.gt.0.0_dp) then 
 		na_t = na + 1
 		allocate( agrid_t(1:na_t) )
 		allocate( agrid_t_ind(1:na_t) )
 		agrid_t = [agrid,Y_a_threshold]
 		call Sort(na_t,agrid_t,agrid_t,agrid_t_ind)
 	else 
 		na_t    = na
 		allocate( agrid_t(1:na_t) )
 		agrid_t = agrid
	end if
	! print*, 'a_aux=',a_aux

	! Allocate variables that depend on na_t
		! Grids for Y and MB
		allocate( YGRID_t(na_t,nz,nx)  )
		allocate( MBGRID_t(na_t,nz,nx) )
		! Allocate size of policy function on adjusted grid:
		allocate( Cons_t(MaxAge,na_t,nz,nlambda,ne,nx) )
		allocate( Hours_t(MaxAge,na_t,nz,nlambda,ne,nx) )
		allocate( Aprime_t(MaxAge,na_t,nz,nlambda,ne,nx) )


		
end Subroutine Asset_Grid_Threshold

!========================================================================================
!========================================================================================
! Y_labor: Evaluate after tax labor income, takes into account state and age
!
! Usage: Y_h = Y_h(h_in,age_in,lambda_in,e_in,Wage_in)
!
! Input: h_in  	  , real(dp), value of hours
!		 age_in   , integer , age
!		 lambda_in, integer , index of lambda
!		 e_in     , integer , index of e
!		 Wage_in  , real(dp), Wage
!
! Output: Y_h , real(dp), After tax labor income
!
!
	FUNCTION Y_h(h_in,age_in,lambda_in,e_in,Wage_in)
		IMPLICIT NONE   
		real(DP), intent(in)  :: h_in, Wage_in
		integer , intent(in)  :: age_in, lambda_in, e_in
		real(DP)              :: Y_h
		
		Y_h = psi*( Wage_in*eff_un(age_in,lambda_in,e_in)*h_in)**(1.0_dp-tauPL)

	END  FUNCTION Y_h

	FUNCTION MB_h(h_in,age_in,lambda_in,e_in,Wage_in)
		IMPLICIT NONE   
		real(DP), intent(in)  :: h_in, Wage_in
		integer , intent(in)  :: age_in, lambda_in, e_in
		real(DP)              :: MB_h
		
		MB_h = (1.0_dp-tauPL)*psi*( Wage_in*eff_un(age_in,lambda_in,e_in))**(1.0_dp-tauPL) * h_in**(-tauPL)

	END  FUNCTION MB_h


!========================================================================================
!========================================================================================
! Y_a: Evaluate asset income (after tax wealth) at given asset level and entrepreneurial ability
!
! Usage: Y_a = Y_a(a_in,z_in)
!
! Input: a_in, real(dp), value of assets
!		 z_in, real(dp), value of entrepreneurial abillity (not necessariliy a grid point)
!
! Output: Y_a , real(dp), After tax wealth
!
	FUNCTION Y_a(a_in,z_in,x_in)
		IMPLICIT NONE   
		real(DP), intent(in)  :: a_in
		integer , intent(in)  :: z_in, x_in
		real(DP)              :: Y_a, K, Pr

		! Capital demand 
		K   = min( theta(z_in)*a_in , (mu*P*xz_grid(x_in,z_in)**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )
		! Profits 
		Pr  = P*(xz_grid(x_in,z_in)*K)**mu - (R+DepRate)*K
		! Before tax wealth
		Y_a = ( a_in +  ( Pr + R_z(z_in)*a_in ) *(1.0_DP-tauK) )

		! Compute after tax wealth according to threshold
		if (a_in.le.Y_a_threshold) then 
			Y_a = a_in* (1.0_dp-tauW_bt) + ( Pr + R_z(z_in)*a_in ) *(1.0_DP-tauK)  
		else
			Y_a = Y_a_threshold*(1.0_dp-tauW_bt) + &
				& (a_in - Y_a_threshold) * (1.0_dp-tauW_at) + ( Pr + R_z(z_in)*a_in ) *(1.0_DP-tauK)
		end if
	END  FUNCTION Y_a


!========================================================================================
!========================================================================================
! MB_a: Evaluate asset marginal benefit at given asset level and entrepreneurial ability
!
! Usage: MB_a = MB_a(a_in,z_in)
!
! Input: a_in, real(dp), value of assets
!		 z_in, real(dp), value of entrepreneurial abillity (not necessariliy a grid point)
!
! Output: MB_a , real(dp), marginal benefit from assets
!
	FUNCTION MB_a(a_in,z_in,x_in)
		IMPLICIT NONE   
		real(DP), intent(in)  :: a_in
		integer , intent(in)  :: x_in, z_in
		real(DP) 			  :: MB_a
		real(DP) :: K, Pr, Y_a, tauW

		! Capital demand 
		K   = min( theta(z_in)*a_in , (mu*P*xz_grid(x_in,z_in)**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )
		! Profits 
		Pr  = P*(xz_grid(x_in,z_in)*K)**mu - (R+DepRate)*K
		! Before tax wealth
		Y_a = ( a_in +  ( Pr + R_z(z_in)*a_in ) *(1.0_DP-tauK) )
		if (a_in.le.Y_a_threshold) then 
			tauW = tauW_bt 
		else
			tauW = tauW_at 
		end if

		! After tax marginal benefit of assets
		if (K.lt.theta(z_in)*a_in) then 
			MB_a = (1.0_dp*(1.0_dp-tauW) + R_z(z_in)*(1.0_dp-tauK))
		else 
			MB_a = (1.0_dp*(1.0_dp-tauW) + R_z(z_in)*(1.0_dp-tauK)) &
         	& + (P*mu*((theta(z_in)*xz_grid(x_in,z_in))**mu)*a_in**(mu-1.0_DP)-(R+DepRate)*theta(z_in))*(1.0_dp-tauK)
		endif 

		! print*,' '; print*,'------------------------------'
		! print*,'Inside MB_a(',a_in,z_in,x_in,')=',MB_a
		! print*,'R_z(z_in)=',R_z(z_in)
		! print*,'R_z=',R_z
		! print*,'K=',K,'Pr=',PR,'Y_a=',Y_a
		! print*,'P'
		! print*,'------------------------------'; print*,' '

	END  FUNCTION MB_a

	FUNCTION MB_a_at(a_in,z_in,x_in)
		IMPLICIT NONE   
		real(DP), intent(in)  :: a_in
		integer , intent(in)  :: x_in, z_in
		real(DP)			  :: MB_a_at, K

		! Capital demand 
		K   = min( theta(z_in)*a_in , (mu*P*xz_grid(x_in,z_in)**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )
		! Compute asset marginal benefit - subject to taxes
		if (K.lt.theta(z_in)*a_in) then 
			MB_a_at = (1.0_dp*(1.0_dp-tauW_at) + R_z(z_in)*(1.0_dp-tauK))
		else 
			MB_a_at = (1.0_dp*(1.0_dp-tauW_at) + R_z(z_in)*(1.0_dp-tauK)) &
         	& + (P*mu*((theta(z_in)*xz_grid(x_in,z_in))**mu)*a_in**(mu-1.0_DP)-(R+DepRate)*theta(z_in))*(1.0_dp-tauK)
		endif 

	END  FUNCTION MB_a_at

	FUNCTION MB_a_bt(a_in,z_in,x_in)
		IMPLICIT NONE   
		real(DP), intent(in)  :: a_in
		integer , intent(in)  :: x_in, z_in
		real(DP)             :: MB_a_bt, K

		! Capital demand 
		K   = min( theta(z_in)*a_in , (mu*P*xz_grid(x_in,z_in)**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )
		! Compute asset marginal benefit - subject to taxes
		if (K.lt.theta(z_in)*a_in) then 
			MB_a_bt = (1.0_dp*(1.0_dp-tauW_bt) + R_z(z_in)*(1.0_dp-tauK))
		else 
			MB_a_bt = (1.0_dp*(1.0_dp-tauW_bt) + R_z(z_in)*(1.0_dp-tauK)) &
         	& + (P*mu*((theta(z_in)*xz_grid(x_in,z_in))**mu)*a_in**(mu-1.0_DP)-(R+DepRate)*theta(z_in))*(1.0_dp-tauK)
		endif 

	END  FUNCTION MB_a_bt


!========================================================================================
!========================================================================================
! FOC_R: Evaluate square residual of Euler equation at current state and candidate savings a'
!		 The FOC is for the retirement period (hence the R subscript)
!
! Usage: res = FOC_R(aprimet)
!
! Input: aprimet   , real(dp), value of assets in next period
!
! Implicit inputs: ai     , integer, index of current level of assets
!				   zi     , integer, index of agent's permanent entreprenurial ability
!				   lambdai, integer, index of agent's permanent labor income component
!				   ei     , integer, index of agent's transitory labor income component at age Ret_Age
!				   age    , integer, agent's age
!
! Output: res , real(dp), residual of the Euler equatioin under retirement
!
	FUNCTION FOC_R(aprimet,state)
		IMPLICIT NONE   
		real(DP), intent(in)  :: aprimet
		real(DP), dimension(:), intent(in) :: state
		real(DP)              :: FOC_R
		real(DP)              :: MB_aprime(nx), yprime(nx), cprime(nx), euler_power
		integer               :: age_in, a_in, z_in, l_in, e_in, x_in, xp_ind 

		! Allocate state and get indeces
		age_in = int(state(1))
		a_in   = int(state(2))
		z_in   = int(state(3))
		l_in   = int(state(4))
		e_in   = int(state(5))
		x_in   = int(state(6))

		! Set the power used in the Euler equation for the retirement period
		if (NSU_Switch.eqv..true.) then 
			! Non-Separable Utility
			euler_power = (1.0_dp/((1.0_dp-sigma)*gamma-1.0_dp))
		else 
			! Separable Utility
			euler_power = (-1.0_dp/sigma)
		end if 

		do xp_ind=1,nx
		! Compute marginal benefit of next period at a' (given zi)
		MB_aprime(xp_ind) = MB_a(aprimet,z_in,xp_ind)
		 
		! Compute asset income of next period at a' (given zi)
		yprime(xp_ind)   = Y_a(aprimet,z_in,xp_ind)
		 
		! Compute consumption of next period given a' (given zi, lambdai and ei)
			! The value of c' comes from interpolating next period's consumption policy function
			! The policy function is defined over a grid of asset income Y, and is interpolated for y'
		cprime(xp_ind)   = Linear_Int(Ygrid_t(:,z_in,xp_ind), Cons_t(age_in+1,:,z_in,l_in,e_in,xp_ind),na_t, yprime(xp_ind))    
		enddo 

		! Evaluate square residual of Euler equation at current state (given by (ai,zi,lambdai,ei)) and savings given by a'
		FOC_R	= ( (YGRID_t(a_in,z_in,x_in)+RetY_lambda_e(l_in,e_in)-aprimet)    &
		    & - ( beta*survP(age_in) * sum( pr_x(x_in,:,z_in,age_in) * MB_aprime * cprime**(1.0_dp/euler_power) ) &
		    &   + (1.0_dp-survP(age_in))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)*chi_bq* &
		    &                  ((1.0_dp-tau_bq)*aprimet+bq_0)**(1.0_dp/euler_power) )  &
		    & **euler_power) ** 2.0_DP

	END  FUNCTION FOC_R

	FUNCTION FOC_R_Transition(aprimet,state)
		IMPLICIT NONE   
		real(DP), intent(in)  :: aprimet
		real(DP), dimension(:), intent(in) :: state
		real(DP)              :: FOC_R_Transition
		real(DP)              :: MB_aprime(nx), yprime(nx), cprime(nx), euler_power
		integer               :: age_in, a_in, z_in, l_in, e_in, x_in, t_in, xp_ind 

		! Allocate state and get indeces
		age_in = int(state(1))
		a_in   = int(state(2))
		z_in   = int(state(3))
		l_in   = int(state(4))
		e_in   = int(state(5))
		x_in   = int(state(6))
		t_in   = int(state(7))

		! Set the power used in the Euler equation for the retirement period
		if (NSU_Switch.eqv..true.) then 
			! Non-Separable Utility
			euler_power = (1.0_dp/((1.0_dp-sigma)*gamma-1.0_dp))
		else 
			! Separable Utility
			euler_power = (-1.0_dp/sigma)
		end if 

		! Set prices for next period
		R = R_tr(t_in+1); P = P_tr(t_in+1);

		do xp_ind=1,nx
		! Compute marginal benefit of next period at a' (given zi)
		MB_aprime(xp_ind) = MB_a(aprimet,z_in,xp_ind)
		 
		! Compute asset income of next period at a' (given zi)
		yprime(xp_ind)   = Y_a(aprimet,z_in,xp_ind)
		 
		! Compute consumption of next period given a' (given zi, lambdai and ei)
			! The value of c' comes from interpolating next period's consumption policy function
			! The policy function is defined over a grid of asset income Y, and is interpolated for y'
		cprime(xp_ind)   = Linear_Int(Ygrid_t(:,z_in,xp_ind), Cons_t_pr(age_in+1,:,z_in,l_in,e_in,xp_ind),na_t, yprime(xp_ind))    
		enddo 

		! Evaluate square residual of Euler equation at current state (given by (ai,zi,lambdai,ei)) and savings given by a'
		FOC_R_Transition	= ( (YGRID_t(a_in,z_in,x_in)+RetY_lambda_e(l_in,e_in)-aprimet)    &
		        & - ( beta*survP(age_in) * sum( pr_x(x_in,:,z_in,age_in) * MB_aprime * cprime**(1.0_dp/euler_power) ) & 
		        & + (1.0_dp-survP(age_in))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)*chi_bq* &
		        &				  ((1.0_dp-tau_bq)*aprimet+bq_0)**(1.0_dp/euler_power) ) &
		        & **euler_power) ** 2.0_DP

	END  FUNCTION FOC_R_Transition

!========================================================================================
!========================================================================================
!========================================================================================
!========================================================================================
! FOC_WH: Evaluate square residual of Euler equation at current state and candidate savings a'
!		  The FOC is for the working period (hence the W subscript)
!         Hours are found by solving the FOC for labor (hence the H subscript)
!
! Usage: res = FOC_WH(aprimet)
!
! Input: aprimet   , real(dp), value of assets in next period
!
! Implicit inputs: ai     , integer, index of current level of assets
!				   zi     , integer, index of agent's permanent entreprenurial ability
!				   lambdai, integer, index of agent's permanent labor income component
!				   ei     , integer, index of agent's transitory labor income component
!				   age    , integer, agent's age
!
! Output: res , real(dp), residual of the Euler equation
!
! Remarks: This function performs a optimization for labor. Labor taxes are taken into account
!
FUNCTION FOC_WH(aprimet,state)
	IMPLICIT NONE   
	real(DP), intent(in) :: aprimet
	real(DP), dimension(:), intent(in) :: state
	real(DP)             :: ctemp, ntemp, yprime(nx), MB_aprime(nx), FOC_WH, exp1overcprime, E_MU_cp(nx)
	REAL(DP)             :: brentvaluet, consin, H_min, c_foc, c_budget, cp, hp
	integer              :: ep_ind, xp_ind
	real(DP), dimension(ne) :: cprime, nprime, MU_cp
	REAL(DP), DIMENSION(7)  :: par_FOC
	integer                 :: age_in, a_in, z_in, l_in, e_in, x_in

	! Allocate state and get indeces
	age_in = int(state(1))
	a_in   = int(state(2))
	z_in   = int(state(3))
	l_in   = int(state(4))
	e_in   = int(state(5))
	x_in   = int(state(6))

	par_FOC(1:6) = state

	H_min = 0.000001 

	do xp_ind=1,nx 
	! Compute marginal benefit of next period at a' (given zi)
	MB_aprime(xp_ind) = MB_a(aprimet,z_in,x_in)
	! Compute asset income of next period at a' (given zi)
	yprime(xp_ind)    = Y_a(aprimet,z_in,x_in)
	enddo 

	! Set auxiliary variable for FOC_HA
	par_FOC(7) = aprimet
	! Solve for hours choice by solving the FOC for labor
		brentvaluet = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_HA, brent_tol, ntemp, par_FOC)  
			! 			if (NSU_Switch.eqv..true.) then
			! 				c_foc = (gamma/(1.0_dp-gamma))*(1.0_dp-H_min)*MB_h(H_min,age,lambdai,ei,wage)
			! 			else 
			! 				c_foc = (MB_h(H_min,age,lambdai,ei,wage)*(1.0_dp-H_min)**(gamma)/phi)**(1.0_dp/sigma)
			! 			end if 
			! 			c_budget = YGRID_t(ai,zi) - aprimet   
			! 			if (c_budget.ge.c_foc) then
			! 				ntemp = 0.0_dp
			! 			else
			! 				brentvaluet = brent(H_min, 0.4_DP, 0.99_DP, FOC_HA, brent_tol, ntemp)  
			! 			end if 
	! Current consumption given ntemp
		ctemp   = YGRID_t(a_in,z_in,x_in) + Y_h(ntemp,age_in,l_in,e_in,wage) - aprimet   

	if (NSU_Switch.eqv..true.) then
		if (Progressive_Tax_Switch) then 
			! Non-Separable Utility
			if (Log_Switch.eqv..true.) then
				! I have to evaluate the FOC in expectation over eindx prime given eindx
				! Compute c' for each value of e'
				do xp_ind=1,nx
				DO ep_ind=1,ne
				    cprime(ep_ind) = Linear_Int(Ygrid_t(:,z_in,xp_ind),Cons_t(age_in+1,:,z_in,l_in,ep_ind,xp_ind),&
				    						&	 na_t, yprime(xp_ind) )
				ENDDO
					! Compute the expected value of 1/c' conditional on current ei
					E_MU_cp(xp_ind) = SUM( pr_e(e_in,:) / cprime )
				enddo 
				

				! Evaluate the squared residual of the Euler equation for working period
				FOC_WH   = ( (1.0_dp/ctemp)- ( beta*survP(age_in) * sum( pr_x(x_in,:,z_in,age_in) * MB_aprime * E_MU_cp) &
					         &     +(1.0_dp-survP(age_in))*chi_bq/((1.0_dp-tau_bq)*aprimet+bq_0) ) ) **2.0_DP 
			else
				! I have to evaluate the FOC in expectation over eindx prime given eindx
				! Compute consumption and labor for eachvalue of eindx prime
				do xp_ind=1,nx 
				DO ep_ind=1,ne
				    cp     = Linear_Int(Ygrid_t(:,z_in,xp_ind),Cons_t(age_in+1,:,z_in,l_in,ep_ind,xp_ind), na_t, yprime(xp_ind))
					c_foc  = (gamma/(1.0_dp-gamma))*(1.0_dp-H_min)*MB_h(H_min,age_in,l_in,e_in,wage)
						if (cp.ge.c_foc) then
							hp = 0.0_dp
						else
							par_FOC(4) = ep_ind
							par_FOC(7) = cp 
							brentvaluet = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_H, brent_tol, hp, par_FOC ) 
						end if 
				    MU_cp(ep_ind) = cp*((1.0_dp-sigma)*gamma-1.0_dp) * (1.0_dp-hp)**((1.0_dp-sigma)*(1.0_dp-gamma))
				END DO
				! Compute the expected value of 1/c' conditional on current ei
				E_MU_cp(xp_ind) = SUM( pr_e(e_in,:) * MU_cp )
				enddo 

				! Evaluate the squared residual of the Euler equation for working period
				FOC_WH = ( ctemp**((1.0_dp-sigma)*gamma-1.0_dp) * (1.0_dp-ntemp)**((1.0_dp-sigma)*(1.0_dp-gamma)) & 
					         & - (beta*survP(age_in)* sum(pr_x(x_in,:,z_in,age_in)*MB_aprime*E_MU_cp) & 
					         & + (1.0_dp-survP(age_in))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
					         & *chi_bq*((1.0_dp-tau_bq)*aprimet+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ))**2.0_DP
			end if 
		else ! Linear Taxes 
			ntemp = max(0.0_DP , gamma - (1.0_DP-gamma)*(YGRID_t(a_in,z_in,x_in) - aprimet)/(psi*yh(age_in,l_in,e_in)) )
			ctemp = YGRID_t(a_in,z_in,x_in) + psi*yh(age_in,l_in,e_in) * ntemp - aprimet

			do xp_ind=1,nx
			DO ep_ind=1,ne
			    cprime(ep_ind) = Linear_Int(Ygrid(:,z_in,xp_ind), Cons_t(age_in+1,:,z_in,l_in,ep_ind,xp_ind),na,yprime(xp_ind))
			    nprime(ep_ind) = 1.0_DP-(1.0_DP-gamma)*cprime(ep_ind)/(gamma*psi*yh(age_in,l_in,e_in))
			ENDDO
				nprime = max(0.0_DP,nprime)
				E_MU_cp(xp_ind) = sum( pr_e(e_in,:) * (cprime**(gamma*(1.0_DP-sigma)-1.0_DP)) &
					& *((1.0_DP-nprime)**((1.0_DP-gamma)*(1.0_DP-sigma))))
			enddo 

			FOC_WH = ((ctemp**(gamma*(1.0_DP-sigma)-1))*((1.0_DP-ntemp)**((1.0_DP-gamma)*(1.0_DP-sigma))) &
					& - (beta*survP(age_in)* sum( pr_x(x_in,:,z_in,age_in)*MB_aprime*E_MU_cp ) & 
					&  + (1.0_dp-survP(age_in))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
					&  *chi_bq*((1.0_dp-tau_bq)*aprimet+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ) )**2.0_DP
		end if 
	else 
		! Separable Utility
			! I have to evaluate the FOC in expectation over eindx prime given eindx
			! Compute c' for each value of e'
			do xp_ind=1,nx 
			DO ep_ind=1,ne
				cprime(ep_ind) = Linear_Int(Ygrid_t(:,z_in,xp_ind),Cons_t(age+1,:,z_in,l_in,ep_ind,xp_ind),na_t,yprime(xp_ind))
			ENDDO
			! Compute the expected value of 1/c' conditional on current ei
			E_MU_cp(xp_ind) = SUM( pr_e(e_in,:) / cprime**sigma )
			enddo 

			! Evaluate the squared residual of the Euler equation for working period
			FOC_WH   = (ctemp - 1.0_dp/(beta*survP(age_in)*sum(pr_x(x_in,:,z_in,age_in)*MB_aprime*E_MU_cp) &
			&  + (1.0_dp-survP(age_in))*(1.0_dp+tauC)**(1.0_dp-sigma)*chi_bq/((1.0_dp-tau_bq)*aprimet+bq_0)**sigma)**(1.0_dp/sigma) )**2.0_DP 
	end if 

END  FUNCTION FOC_WH


FUNCTION FOC_WH_Transition(aprimet,state)
	IMPLICIT NONE   
	real(DP), intent(in) :: aprimet
	real(DP), dimension(:), intent(in) :: state
	real(DP)             :: ctemp, ntemp, yprime(nx), MB_aprime(nx), FOC_WH_Transition, exp1overcprime, E_MU_cp(nx)
	REAL(DP)             :: brentvaluet, consin, H_min, c_foc, c_budget, cp, hp
	integer              :: ep_ind, xp_ind
	real(DP), dimension(ne) :: cprime, nprime, MU_cp
	REAL(DP), DIMENSION(7)  :: par_FOC
	integer                 :: age_in, a_in, z_in, l_in, e_in, x_in, t_in

	! Allocate state and get indeces
	age_in = int(state(1))
	a_in   = int(state(2))
	z_in   = int(state(3))
	l_in   = int(state(4))
	e_in   = int(state(5))
	x_in   = int(state(6))
	t_in   = int(state(7))

	par_FOC(1:6) = state

	H_min = 0.000001 


	! Update prices for computing MB and yprime
		P = P_tr(t_in+1); R = R_tr(t_in+1) 
	do xp_ind=1,nx 
	! Compute marginal benefit of next period at a' (given zi)
	MB_aprime(xp_ind) = MB_a(aprimet,z_in,x_in)
	! Compute asset income of next period at a' (given zi)
	yprime(xp_ind)    = Y_a(aprimet,z_in,x_in)
	enddo 
	! Set prices to current period
		P = P_tr(t_in); R = R_tr(t_in) 

	! Set auxiliary variable for FOC_HA
	par_FOC(7) = aprimet
	! Solve for hours choice by solving the FOC for labor
		brentvaluet = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_HA, brent_tol, ntemp, par_FOC)  
			! 			if (NSU_Switch.eqv..true.) then
			! 				c_foc = (gamma/(1.0_dp-gamma))*(1.0_dp-H_min)*MB_h(H_min,age,lambdai,ei,wage)
			! 			else 
			! 				c_foc = (MB_h(H_min,age,lambdai,ei,wage)*(1.0_dp-H_min)**(gamma)/phi)**(1.0_dp/sigma)
			! 			end if 
			! 			c_budget = YGRID_t(ai,zi) - aprimet   
			! 			if (c_budget.ge.c_foc) then
			! 				ntemp = 0.0_dp
			! 			else
			! 				brentvaluet = brent(H_min, 0.4_DP, 0.99_DP, FOC_HA, brent_tol, ntemp)  
			! 			end if 
	! Current consumption given ntemp
		ctemp   = YGRID_t(a_in,z_in,x_in) + Y_h(ntemp,age_in,l_in,e_in,wage) - aprimet   

	if (NSU_Switch.eqv..true.) then
		if (Progressive_Tax_Switch) then 
			! Non-Separable Utility
			if (Log_Switch.eqv..true.) then
				! I have to evaluate the FOC in expectation over eindx prime given eindx
				! Compute c' for each value of e'
				do xp_ind=1,nx
				DO ep_ind=1,ne
				    cprime(ep_ind) = Linear_Int(Ygrid_t(:,z_in,xp_ind),Cons_t_pr(age_in+1,:,z_in,l_in,ep_ind,xp_ind),&
				    						&	 na_t, yprime(xp_ind) )
				ENDDO
					! Compute the expected value of 1/c' conditional on current ei
					E_MU_cp(xp_ind) = SUM( pr_e(e_in,:) / cprime )
				enddo 
				

				! Evaluate the squared residual of the Euler equation for working period
				FOC_WH_Transition   = &
					& ( (1.0_dp/ctemp)- (beta*survP(age_in) * sum( pr_x(x_in,:,z_in,age_in) * MB_aprime * E_MU_cp) & 
						& + (1.0_dp-survP(age_in))*chi_bq/((1.0_dp-tau_bq)*aprimet+bq_0) ) ) **2.0_DP 
			else
				! I have to evaluate the FOC in expectation over eindx prime given eindx
				! Compute consumption and labor for eachvalue of eindx prime
				do xp_ind=1,nx 
				DO ep_ind=1,ne
				    cp     = Linear_Int(Ygrid_t(:,z_in,xp_ind),Cons_t_pr(age_in+1,:,z_in,l_in,ep_ind,xp_ind), na_t, yprime(xp_ind))
					c_foc  = (gamma/(1.0_dp-gamma))*(1.0_dp-H_min)*MB_h(H_min,age_in,l_in,e_in,wage_tr(t_in))
						if (cp.ge.c_foc) then
							hp = 0.0_dp
						else
							par_FOC(4) = ep_ind
							par_FOC(7) = cp 
							brentvaluet = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_H, brent_tol, hp, par_FOC ) 
						end if 
				    MU_cp(ep_ind) = cp*((1.0_dp-sigma)*gamma-1.0_dp) * (1.0_dp-hp)**((1.0_dp-sigma)*(1.0_dp-gamma))
				END DO
				! Compute the expected value of 1/c' conditional on current ei
				E_MU_cp(xp_ind) = SUM( pr_e(e_in,:) * MU_cp )
				enddo 

				! Evaluate the squared residual of the Euler equation for working period
				FOC_WH_Transition = &
					& ( ctemp**((1.0_dp-sigma)*gamma-1.0_dp) * (1.0_dp-ntemp)**((1.0_dp-sigma)*(1.0_dp-gamma)) & 
					& - (beta*survP(age_in)* sum(pr_x(x_in,:,z_in,age_in)*MB_aprime*E_MU_cp) &
					& + (1.0_dp-survP(age_in))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
					&  *chi_bq*((1.0_dp-tau_bq)*aprimet+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ))**2.0_DP
			end if 
		else ! Linear Taxes 
			ntemp = max(0.0_DP , gamma - (1.0_DP-gamma)*(YGRID_t(a_in,z_in,x_in) - aprimet)/(psi*yh(age_in,l_in,e_in)) )
			ctemp = YGRID_t(a_in,z_in,x_in) + psi*yh(age_in,l_in,e_in) * ntemp - aprimet

			do xp_ind=1,nx
			DO ep_ind=1,ne
			    cprime(ep_ind) = &
			    	& Linear_Int(Ygrid(:,z_in,xp_ind), Cons_t_pr(age_in+1,:,z_in,l_in,ep_ind,xp_ind),na,yprime(xp_ind))
			    nprime(ep_ind) = 1.0_DP-(1.0_DP-gamma)*cprime(ep_ind)/(gamma*psi*yh(age_in,l_in,e_in))
			ENDDO
				nprime = max(0.0_DP,nprime)
				E_MU_cp(xp_ind) = sum( pr_e(e_in,:) * (cprime**(gamma*(1.0_DP-sigma)-1.0_DP)) &
					& *((1.0_DP-nprime)**((1.0_DP-gamma)*(1.0_DP-sigma))))
			enddo 

			FOC_WH_Transition = ((ctemp**(gamma*(1.0_DP-sigma)-1))*((1.0_DP-ntemp)**((1.0_DP-gamma)*(1.0_DP-sigma))) &
					& - (beta*survP(age_in)* sum( pr_x(x_in,:,z_in,age_in)*MB_aprime*E_MU_cp ) &
				    &  + (1.0_dp-survP(age_in))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
				    & *chi_bq*((1.0_dp-tau_bq)*aprimet+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ) )**2.0_DP
		end if 
	else 
		! Separable Utility
			! I have to evaluate the FOC in expectation over eindx prime given eindx
			! Compute c' for each value of e'
			do xp_ind=1,nx 
			DO ep_ind=1,ne
				cprime(ep_ind) = Linear_Int(Ygrid_t(:,z_in,xp_ind),Cons_t_pr(age+1,:,z_in,l_in,ep_ind,xp_ind),na_t,yprime(xp_ind))
			ENDDO
			! Compute the expected value of 1/c' conditional on current ei
			E_MU_cp(xp_ind) = SUM( pr_e(e_in,:) / cprime**sigma )
			enddo 

			! Evaluate the squared residual of the Euler equation for working period
			FOC_WH_Transition  = &
			& (ctemp - 1.0_dp/(beta*survP(age_in)*sum(pr_x(x_in,:,z_in,age_in)*MB_aprime*E_MU_cp) & 
			&  + (1.0_dp-survP(age_in))*(1.0_dp+tauC)**(1.0_dp-sigma)*chi_bq/ &
			& 			((1.0_dp-tau_bq)*aprimet+bq_0)**sigma )**(1.0_dp/sigma) )**2.0_DP
	end if 

END  FUNCTION FOC_WH_Transition


!========================================================================================
!========================================================================================
! FOC_H: Evaluate square residual of FOC for hours at current state and candidate hours h
!	     This function takes as given current consumption
!
! Usage: res = FOC_H(hoursin)
!
! Input: hoursin   , real(dp), value of assets in next period
!
! Implicit inputs: ai     , integer , index of current level of assets
!				   zi     , integer , index of agent's permanent entreprenurial ability
!				   lambdai, integer , index of agent's permanent labor income component
!				   ei     , integer , index of agent's transitory labor income component
!				   age    , integer , agent's age
!                  consin , real(dp), Current consumption 
!
! Output: res , real(dp), residual of the FOC for hours
!
! Remarks: Labor taxes are taken into account. But consumption is given and does not 
!          internalize changes in taxes. This is necessary when solving EGM.
!
	FUNCTION FOC_H(hoursin,par)
		IMPLICIT NONE   
		real(DP), intent(in) 	:: hoursin
		real(DP), dimension(:), intent(in) :: par
		real(DP)             	:: FOC_H
		real(DP)                :: cons_in
		integer                 :: age_in, a_in, z_in, l_in, e_in, x_in

		! Allocate state and get indeces
		age_in = int(par(1))
		a_in   = int(par(2))
		z_in   = int(par(3))
		l_in   = int(par(4))
		e_in   = int(par(5))
		x_in   = int(par(6))
		cons_in = par(7)

		if (NSU_Switch.eqv..true.) then 
			! Non-Separable Utility 
			FOC_H = ( cons_in - (gamma/(1.0_dp-gamma))*(1.0_dp-hoursin)*MB_h(hoursin,age_in,l_in,e_in,wage) )**2.0_DP 
		else 
			! Separable Utility 
			FOC_H = ( MB_h(hoursin,age_in,l_in,e_in,wage)*(1.0_dp-hoursin)**(gamma) - phi*cons_in**(sigma) )**2.0_DP 
		end if 

	END  FUNCTION FOC_H

	FUNCTION FOC_H_NSU(hoursin,par)
		IMPLICIT NONE   
		real(DP), intent(in) 	:: hoursin
		real(DP), dimension(:), intent(in) :: par
		real(DP)             	:: FOC_H_NSU, cons, E_MU_cp(nx)
		real(DP), dimension(ne) :: MU_cp
		real(DP)                :: MB_aprime(nx)
		integer                 :: age_in, a_in, z_in, l_in, e_in, x_in, xp_ind

		! Allocate state and get indeces
		age_in = int(par(1))
		a_in   = int(par(2))
		z_in   = int(par(3))
		l_in   = int(par(4))
		e_in   = int(par(5))
		x_in   = int(par(6))
		MB_aprime = par(7:6+nx)

		! Compute current consumption implied by hours and labor FOC
			cons  = (gamma/(1.0_dp-gamma))*(1.0_dp-hoursin)*MB_h(hoursin,age_in,l_in,e_in,wage)
		do xp_ind=1,nx
		! Compute marginal utility of consumption for next period at a' (for all values of e)
	    	MU_cp = Cons_t(age_in+1,a_in,z_in,l_in,:,xp_ind)**((1.0_dp-sigma)*gamma-1.0_dp) &
	    	        * (1.0_dp-Hours_t(age_in+1,a_in,z_in,l_in,:,xp_ind))**((1.0_dp-sigma)*(1.0_dp-gamma))
	    ! Compute expected value of marginal uitility of consumption
	    	E_MU_cp(xp_ind)  = SUM( pr_e(e_in,:) * MU_cp )
	    enddo 
		! Compute square residual of Euler FOC
			FOC_H_NSU = ( cons**((1.0_dp-sigma)*gamma-1.0_dp) * (1.0_dp-hoursin)**((1.0_dp-sigma)*(1.0_dp-gamma)) & 
			    & - (beta*survP(age_in)*sum(pr_x(x_in,:,z_in,age_in)*MB_aprime*E_MU_cp) &
			    &  + (1.0_dp-survP(age_in))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
			    & *chi_bq*((1.0_dp-tau_bq)*agrid_t(a_in)+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ))**2.0_DP

	END  FUNCTION FOC_H_NSU


	FUNCTION FOC_H_NSU_Transition(hoursin,par)
		IMPLICIT NONE   
		real(DP), intent(in) 	:: hoursin
		real(DP), dimension(:), intent(in) :: par
		real(DP)             	:: FOC_H_NSU_Transition, cons, E_MU_cp(nx)
		real(DP), dimension(ne) :: MU_cp
		real(DP)                :: MB_aprime(nx)
		integer                 :: age_in, a_in, z_in, l_in, e_in, x_in, t_in, xp_ind

		! Allocate state and get indeces
		age_in = int(par(1))
		a_in   = int(par(2))
		z_in   = int(par(3))
		l_in   = int(par(4))
		e_in   = int(par(5))
		x_in   = int(par(6))
		t_in   = int(par(7))
		MB_aprime = par(8:7+nx)

		! Compute current consumption implied by hours and labor FOC
			cons  = (gamma/(1.0_dp-gamma))*(1.0_dp-hoursin)*MB_h(hoursin,age_in,l_in,e_in,wage)
		do xp_ind=1,nx
		! Compute marginal utility of consumption for next period at a' (for all values of e)
	    	MU_cp = Cons_t_pr(age_in+1,a_in,z_in,l_in,:,xp_ind)**((1.0_dp-sigma)*gamma-1.0_dp) &
	    	        * (1.0_dp-Hours_t_pr(age_in+1,a_in,z_in,l_in,:,xp_ind))**((1.0_dp-sigma)*(1.0_dp-gamma))
	    ! Compute expected value of marginal uitility of consumption
	    	E_MU_cp(xp_ind)  = SUM( pr_e(e_in,:) * MU_cp )
	    enddo 
		! Compute square residual of Euler FOC
			FOC_H_NSU_Transition = &
				& ( cons**((1.0_dp-sigma)*gamma-1.0_dp) * (1.0_dp-hoursin)**((1.0_dp-sigma)*(1.0_dp-gamma)) & 
			    & - (beta*survP(age_in)*sum(pr_x(x_in,:,z_in,age_in)*MB_aprime*E_MU_cp) &
			    &  + (1.0_dp-survP(age_in))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
			    & *chi_bq*((1.0_dp-tau_bq)*agrid_t(a_in)+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ))**2.0_DP

	END  FUNCTION FOC_H_NSU_Transition

!========================================================================================
!========================================================================================
!========================================================================================
! FOC_HA: Evaluate square residual of FOC for hours at current state and candidate hours h
!	      This function takes as given the saving decision a'
!
! Usage: res = FOC_HA(hoursin)
!
! Input: hoursin   , real(dp), value of assets in next period
!
! Implicit inputs: ai     , integer , index of current level of assets
!				   zi     , integer , index of agent's permanent entreprenurial ability
!				   lambdai, integer , index of agent's permanent labor income component
!				   ei     , integer , index of agent's transitory labor income component
!				   age    , integer , agent's age
!                  ain    , real(dp), savings a'
!
! Output: res , real(dp), residual of the FOC for hours
!
! Remarks: Labor taxes are taken into account.
!
	FUNCTION FOC_HA(hoursin,par)
		IMPLICIT NONE   
		real(DP), intent(in) :: hoursin
		real(DP), dimension(:), intent(in) :: par
		real(DP)   			 :: FOC_HA, cons
		real(DP)             :: ap_in
		integer              :: age_in, a_in, z_in, l_in, e_in, x_in

		! Allocate state and get indeces
		age_in = int(par(1))
		a_in   = int(par(2))
		z_in   = int(par(3))
		l_in   = int(par(4))
		e_in   = int(par(5))
		x_in   = int(par(6))
		ap_in  = par(7)

		! Consumption given ain and hoursin
			cons   = YGRID_t(a_in,z_in,x_in)+  Y_h(hoursin,age_in,l_in,e_in,wage) - ap_in
		if (NSU_Switch.eqv..true.) then
			! Non-Separable Utility 
			FOC_HA = ( cons - (gamma/(1.0_dp-gamma))*(1.0_dp-hoursin)*MB_h(hoursin,age_in,l_in,e_in,wage) )**2.0_DP 
		else
			! Non-Separable Utility 
			FOC_HA = ( MB_h(hoursin,age_in,l_in,e_in,wage)*(1.0_dp-hoursin)**(gamma) - phi*cons**(sigma) )**2.0_DP 
		end if 
	END  FUNCTION FOC_HA


!========================================================================================
!========================================================================================
!========================================================================================

SUBROUTINE EGM_RETIREMENT_WORKING_PERIOD()
	use omp_lib
	IMPLICIT NONE
	REAL(DP) :: brentvalue, C_foc, H_min, euler_power, chi_aux
	INTEGER  :: tempai, sw 
	REAL(DP), DIMENSION(na_t+nz*nx+1)  	:: EndoCons, EndoYgrid, EndoHours
	INTEGER , DIMENSION(na_t+nz*nx+1)   :: sort_ind 
	REAL(DP), DIMENSION(na_t,nz,nx) :: Wealth_mat
	REAL(DP), DIMENSION(6)       	:: state_FOC
	REAL(DP), DIMENSION(7)       	:: par_FOC
	REAL(DP), DIMENSION(nx)       	:: MB_aprime_t
	integer  :: age, ai, zi, lambdai, ei, xi, xp_ind
	real(dp), dimension(na_t+nz*nx+1) :: EndoYgrid_sort

	!$ call omp_set_num_threads(nz)

	! Set a minimum value for labor to check in the FOC
		H_min = 0.000001_dp

	! Set the power used in the Euler equation for the retirement period
		if (NSU_Switch.eqv..true.) then 
			! Non-Separable Utility
				euler_power = (1.0_dp/((1.0_dp-sigma)*gamma-1.0_dp))
		else 
			! Separable Utility
				euler_power = (-1.0_dp/sigma)
		end if 

	! Set value of auxiliary chi parameter for last period's saving
		chi_aux = ((1.0_dp+tauC)*chi_bq)**(1.0_dp/(1-gamma*(1.0_dp-sigma)))/(1.0_dp+tauC)

	! Compute wealth given current R and P
		Wealth_mat = Wealth_Matrix_t(R,P)

		! print*, 'R=',R,'P=',P, 'W=',wage, 'na=', na, 'na_t=', na_t
	!========================================================================================
	!------RETIREMENT PERIOD-----------------------------------------------------------------

	! SET HOURS TO ZERO SO THAT RETIRED HAS ZERO HOURS
	Hours_t = 0.0_DP

	! Last period of life
	age=MaxAge
	!$omp parallel do private(lambdai,ei,ai,xi)
	DO zi=1,nz
	!$omp parallel do private(lambdai,ei,ai)
	DO xi=1,nx
    DO lambdai=1,nlambda
    DO ei=1,ne
    DO ai=1,na_t
    	if ((YGRID_t(ai,zi,xi) + RetY_lambda_e(lambdai,ei)).le.(bq_0/chi_aux) ) then 
        Cons_t(age,ai,zi,lambdai,ei,xi)   =  YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei) 
        else
        Cons_t(age,ai,zi,lambdai,ei,xi)   = (YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei)+bq_0)/(1.0_dp+chi_aux)
        endif
        Aprime_t(age,ai,zi,lambdai,ei,xi) = YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei)-Cons_t(age,ai,zi,lambdai,ei,xi)
	ENDDO ! ai
    ENDDO ! ei
	ENDDO ! lambdai
	ENDDO ! xi
	ENDDO ! zi
	
	
	! Rest of retirement
	DO age=MaxAge-1,RetAge,-1
	!$omp parallel do private(lambdai,ei,ai,xi,xp_ind,EndoCons,EndoYgrid,sw,sort_ind,tempai, &
	!$omp& state_FOC,par_FOC,MB_aprime_t,EndoYgrid_sort)
    DO zi=1,nz
    !$omp parallel do private(lambdai,ei,ai,xp_ind,EndoCons,EndoYgrid,sw,sort_ind,tempai, &
	!$omp& state_FOC,par_FOC,MB_aprime_t,EndoYgrid_sort)
    DO xi=1,nx
    DO lambdai=1,nlambda
    DO ei=1,ne
    	! Endogenous grid and consumption are initialized
	    EndoCons  = big_p 
	    EndoYgrid = big_p 
	    sw 		  = 0
    DO ai=1,na_t 
		if (abs(agrid_t(ai)-Y_a_threshold).lt.1e-8) then 
			sw 			  = sw+1	
    		
    		! Consumption on endogenous grid and implied asset income under tauW_bt
    		do xp_ind = 1,nx
				MB_aprime_t(xp_ind) = MB_a_bt(agrid_t(ai),zi,xp_ind)
    		enddo 
    		EndoCons(ai)  =  (beta*survP(age)* 	&
				& sum(pr_x(xi,:,zi,age)*MB_aprime_t*Cons_t(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power))&
				& + (1.0_dp-survP(age))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
		    	& * chi_bq*((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ) **euler_power
	        EndoYgrid(ai) = agrid_t(ai) +  EndoCons(ai) - RetY_lambda_e(lambdai,ei)
	        
	        ! Consumption on endogenous grid and implied asset income under tauW_at
	        do xp_ind = 1,nx 	
				MB_aprime_t(xp_ind) = MB_a_at(agrid_t(ai),zi,xp_ind)
    		enddo 
	        EndoCons(na_t+sw)  = (beta*survP(age)*	&
    			& sum(pr_x(xi,:,zi,age)*MB_aprime_t*Cons_t(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power)) &
    			& + (1.0_dp-survP(age))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
		    	& * chi_bq*((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ) **euler_power
	    	EndoYgrid(na_t+sw) = agrid_t(ai) +  EndoCons(na_t+sw) - RetY_lambda_e(lambdai,ei)

	  !   	print*, ' '
			! print*, ' Threshold test - Retirement'
			! print*, ' Current State', age, ai, zi, lambdai, ei, xi
			! print*, ' ', (pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold))
			! print*, ' ', (any((pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold)).lt.1e-8))
			! print*, ' ', MBGRID_t(ai,zi,:)
			! print*, ' ', MB_a_bt(agrid(ai),zi,xi)
			! print*, ' ', MB_a_at(agrid(ai),zi,xi)
			! print*, ' ', EndoCons(ai)
			! print*, ' ', EndoCons(na_t+sw)
			! print*, ' '
	    else 
	    	! Consumption on endogenous grid and implied asset income
	    	EndoCons(ai)  = (beta*survP(age)* 	&
	    				& sum(pr_x(xi,:,zi,age)*MBGRID_t(ai,zi,:)*Cons_t(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power)) &
	    				& + (1.0_dp-survP(age))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
				    	& * chi_bq*((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ) **euler_power
	        EndoYgrid(ai) = agrid_t(ai) +  EndoCons(ai) - RetY_lambda_e(lambdai,ei)
	    end if 


	 !    if (any(isnan(EndoCons))) then 
	 !    	print*,' '
		! 	print*,' '
		! 	print*,' '
		! 	print*,' Endo Consumption in retirement', age,ai,zi,lambdai,ei,xi
		! 	print*,' ',EndoCons(ai),(beta*survP(age)* 	&
	 !    				& sum(pr_x(xi,:,zi,age)*MBGRID_t(ai,zi,:)*Cons_t(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power)) ) **euler_power
		! 	print*,' MBGRID_t=',MBGRID_t(ai,zi,:)
		! 	print*,' cons(t+1)=',Cons_t(age+1,ai,zi,lambdai,ei,:)
		! 	print*,' Size(endoCons)=',size(EndoCons)
		! 	print*,' sw=',sw,'na_t=',na_t
		! 	print*,' Threshold_test=',(any((pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold)).lt.1e-8))
		! 	print*,' '
		! 	print*,' ',EndoCons
		! 	print*,' '
		! 	STOP
		! endif 
	ENDDO ! ai

		
	
	! Sort endogenous grid for interpolation
	call Sort(na_t+nx*nz+1,EndoYgrid,EndoYgrid_sort,sort_ind)
	! print*,' Yendo'
	! print*, EndoYgrid
	! print*,' Yendo_sort'
	! print*, EndoYgrid_sort
	! print*,' indices'
	! print*, sort_ind

	EndoYgrid = EndoYgrid_sort
	EndoCons = EndoCons(sort_ind)

	! print*, ' '
	! print*, ' isnan(endocons)', any(isnan(EndoCons))

	! Find  decision rules on exogenous grids
		! decision rules are obtained taking care of extrapolations
		tempai=1           
		DO WHILE ( YGRID_t(tempai,zi,xi) .lt. EndoYgrid(1) )
            tempai = tempai + 1
        ENDDO
	                
		DO ai=tempai,na_t              
		    ! CONSUMPTION ON EXOGENOUS GRIDS 
		    Cons_t(age,ai,zi,lambdai, ei,xi)  = Linear_Int(EndoYgrid(1:na_t+sw), EndoCons(1:na_t+sw),na_t+sw, YGRID_t(ai,zi,xi))
		    Aprime_t(age,ai,zi,lambdai,ei,xi) = YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei)-Cons_t(age, ai, zi, lambdai, ei,xi)
		    
		    If (Aprime_t(age,ai,zi,lambdai,ei,xi).lt.amin) then
		    	Aprime_t(age,ai,zi,lambdai,ei,xi) = amin
	            Cons_t(age,ai,zi,lambdai,ei,xi)   = YGRID_t(ai,zi,xi) + RetY_lambda_e(lambdai,ei) - Aprime_t(age,ai,zi,lambdai,ei,xi)
				IF (Cons_t(age,ai,zi,lambdai,ei,xi) .le. 0.0_DP)  THEN
				    print*,'r1: Cons(age, ai, zi, lambdai,ei)=',Cons_t(age, ai, zi, lambdai, ei,xi)
				ENDIF                   
	        endif 

	        ! ! !$omp critical 
	        ! if (isnan(Cons_t(age,ai,zi,lambdai,ei,xi))) then
	        ! 	print*,' '
	        ! 	print*,' '
	        ! 	print*,' '
	        ! 	print*,' isnan Consumption', age,ai,zi,lambdai,ei,xi
	        ! 	print*,' sw=',sw 
	        ! 	print*,' Yendo'
	        ! 	print*, EndoYgrid(1:na_t+sw)
	        ! 	print*,' Cendo'
	        ! 	print*, EndoCons(1:na_t+sw)
	        ! 	print*,' YGRID'
	        ! 	print*, YGRID_t(ai,zi,xi)
	        ! 	print*,' Linear Interpolation',Linear_Int(EndoYgrid(1:na_t+sw), EndoCons(1:na_t+sw),na_t+sw, YGRID_t(ai,zi,xi))
	        ! 	print*,' Aprime=',Aprime_t(age,ai,zi,lambdai,ei,xi)
	        ! 	call Sort(na_t+nx*nz+1,EndoYgrid,EndoYgrid_sort,sort_ind)
	        ! 	print*,' Yendo'
	        ! 	print*, EndoYgrid_sort
	        ! 	print*,' indices'
	        ! 	print*, sort_ind
	        ! 	print*, ' '
	        ! 	print*, ' ',minval(EndoYgrid),maxval(EndoYgrid)
	        ! 	print*, ' ',minval(EndoYgrid_sort),maxval(EndoYgrid_sort)
	        ! 	print*, ' The end'
	        ! 	STOP
	        ! endif 
	        ! ! !$omp end critical
		ENDDO ! ai  

        ai=1           
        DO WHILE ( YGRID_t(ai,zi,xi) .lt. EndoYgrid(1) )
			! Solve for a' directly by solving the Euler equation for retirement FOC_R
			state_FOC  = (/age,ai,zi,lambdai,ei,xi/)
			brentvalue = brent_p(min(amin,YGRID_t(ai,zi,xi)), (amin+YGRID_t(ai,zi,xi))/2.0_DP , &
			                & YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei) *0.95_DP, &
			                & FOC_R, brent_tol, Aprime_t(age, ai, zi, lambdai,ei,xi),state_FOC)
			
			Cons_t(age,ai,zi,lambdai,ei,xi) = YGRID_t(ai,zi,xi) + RetY_lambda_e(lambdai,ei) - Aprime_t(age,ai,zi,lambdai,ei,xi)
			
			IF (Cons_t(age, ai, zi, lambdai, ei,xi) .le. 0.0_DP)  THEN
			    print*,'r2: Cons(age, ai, zi, lambdai,ei,xi)=',Cons_t(age, ai, zi, lambdai,ei,xi)
			    print*,'Aprime(age, ai, zi, lambdai,ei,xi)=',Aprime_t(age, ai, zi, lambdai,ei,xi)
			    print*,'YGRID(ai,zi,xi)+RetY_lambda_e(lambdai,ei)=',YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei)
			    print*,'YGRID(ai,zi,xi)=',YGRID(ai,zi,xi),'EndoYgrid(1)=',EndoYgrid(1)
			    print*,'RetY_lambda_e(lambdai,ei)=',RetY_lambda_e(lambdai,ei)
			    print*,'lambdai=',lambdai
			ENDIF                   
			ai = ai + 1
		ENDDO  
	              
    ENDDO ! ei     
    ENDDO ! lambda
    ENDDO ! xi 
	ENDDO ! zi
    ENDDO !age

	!------RETIREMENT PERIOD ENDS------------------------------------------------------------
	!========================================================================================
	
	!========================================================================================
	!------Working Period Starts-------------------------------------------------------------

	DO age=RetAge-1,1,-1
	!$omp parallel do private(lambdai,ei,ai,xi,xp_ind,EndoCons,EndoHours,EndoYgrid,sw,sort_ind,tempai, &
	!$omp& C_foc,state_FOC,par_FOC,MB_aprime_t)
    DO zi=1,nz
    !$omp parallel do private(lambdai,ei,ai,xp_ind,EndoCons,EndoHours,EndoYgrid,sw,sort_ind,tempai, &
	!$omp& C_foc,state_FOC,par_FOC,MB_aprime_t)
    DO xi=1,nx
    DO lambdai=1,nlambda
    DO ei=1,ne	
    	! Endogenous grid and consumption are initialized
	    EndoCons  = big_p 
	    EndoYgrid = big_p 
	    EndoHours = big_p 
	    sw 		  = 0                    
        DO ai=1,na_t
        state_FOC  = (/age,ai,zi,lambdai,ei,xi/)
		if (abs(agrid_t(ai)-Y_a_threshold).lt.1e-8) then 
			sw 			  = sw+1	

	    	! Below threshold
    		do xp_ind = 1,nx 	
				MB_aprime_t(xp_ind) = MB_a_bt(agrid_t(ai),zi,xp_ind)
    		enddo 
			call EGM_Working_Period( MB_aprime_t , H_min , state_FOC , & 
			      & EndoCons(ai), EndoHours(ai) , EndoYgrid(ai)  )
			
			! Above threshold
			do xp_ind = 1,nx 	
    				MB_aprime_t(xp_ind) = MB_a_at(agrid_t(ai),zi,xp_ind)
    		enddo 
			call EGM_Working_Period( MB_aprime_t , H_min , state_FOC , & 
			      & EndoCons(na_t+sw), EndoHours(na_t+sw) , EndoYgrid(na_t+sw)  )

	    	!print*, ' '
	    	!print*, State_FOC 
	    	!print*, MB_a_bt(agrid_t(ai),zgrid(zi)), MB_a_at(agrid_t(ai),zgrid(zi))
	  !   	print*, ' '
			! print*, ' Threshold test - Working Period'
			! print*, ' Current State', age, ai, zi, lambdai, ei, xi
			! print*, ' ', (pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold))
			! print*, ' ', (any((pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold)).lt.1e-8))
			! print*, ' ', MBGRID_t(ai,zi,:)
			! print*, ' ', MB_a_bt(agrid(ai),zi,xi)
			! print*, ' ', MB_a_at(agrid(ai),zi,xi)
			! print*, ' ', MB_aprime_t
			! print*, ' ', EndoCons(ai)
			! print*, ' ', EndoCons(na_t+sw)
			! print*, ' '
		else 
			! Usual EGM
			call EGM_Working_Period( MBGRID_t(ai,zi,:) , H_min , state_FOC , & 
			      & EndoCons(ai), EndoHours(ai) , EndoYgrid(ai)  )
			! print*, ' '; print*, 'Test MB=',MBGRID_t(ai,zi,:); print*, ' '
		end if 

    	ENDDO ! ai

	    ! !$omp critical
		if (any(isnan(EndoCons))) then 
			print*, "isnan - Consumption endogenous"
			print*, age,zi,ai,lambdai,ei,xi
			print*, EndoCons
			print*, 'wage, R, R_C',wage, R, R_C 
			print*, 'R_z',R_z 
			STOP 
		end if 
		! !$omp end critical

    ! Sort endogenous grid for interpolation
	call Sort(na_t+nx*nz+1,EndoYgrid,EndoYgrid,sort_ind)
	EndoHours = EndoHours(sort_ind)
	EndoCons  = EndoCons(sort_ind)
! 	print*, ' '
! 	do ai=1,na_t+1
! 		print*, EndoHours(ai), EndoCons(ai), EndoYgrid(ai)
! 	end do
! 	print*, ' '

		if (any(isnan(EndoCons)).or.any(isnan(EndoHours)).or.any(isnan(EndoYgrid))) then 
			print*, "isnan - Consumption working 4"
			print*, age,lambdai,ai,zi,ei,xi
			STOP 
		end if 

    	! Find  decision rules on exogenous grids
        tempai=1           
        DO WHILE ( YGRID_t(tempai,zi,xi) .lt. EndoYgrid(1) )
              tempai = tempai + 1
        ENDDO
	                
		! decision rules are obtained taking care of extrapolations
		DO ai=tempai,na_t 
			! Interpolate for value of consumption in exogenous grid
			Cons_t(age,ai,zi,lambdai,ei,xi) = Linear_Int(EndoYgrid(1:na_t+sw), EndoCons(1:na_t+sw),na_t+sw, YGRID_t(ai,zi,xi))

			if (Progressive_Tax_Switch.eqv..true.) then
			! Check FOC for implied consumption 
				if (NSU_Switch.eqv..true.) then 
					! Non-Separable Utility
					C_foc = (gamma/(1.0_dp-gamma))*(1.0_dp-H_min)*MB_h(H_min,age,lambdai,ei,wage)
				else
					! Separable Utility
					C_foc = (MB_h(H_min,age,lambdai,ei,wage)*(1.0_dp-H_min)**(gamma)/phi)**(1.0_dp/sigma)
				end if 
			! Hours
				if (Cons_t(age,ai,zi,lambdai,ei,xi).ge.C_foc) then
					Hours_t(age,ai,zi,lambdai,ei,xi) = 0.0_dp
				else
					! Auxiliary variables for solving FOC for hours 
					par_FOC(1:6) = (/age,ai,zi,lambdai,ei,xi/) 
					par_FOC(7)   = Cons_t(age, ai, zi, lambdai,ei,xi)
					! FOC
					brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_H, brent_tol, Hours_t(age,ai,zi,lambdai,ei,xi) , par_FOC ) 
				end if
			else 
				if (NSU_Switch.eqv..true.) then 
					Hours_t(age,ai,zi,lambdai,ei,xi) = max( 0.0_dp , &
						&  1.0_DP - (1.0_DP-gamma)*Cons_t(age,ai,zi,lambdai,ei,xi)/(gamma*psi*yh(age,lambdai,ei)) )
				else 
					Hours_t(age,ai,zi,lambdai,ei,xi) = max( 0.0_DP , &
						&  1.0_DP - phi*Cons_t(age,ai,zi,lambdai,ei,xi)/(psi*yh(age, lambdai,ei)) )
				end if 
			end if 

			if (Hours_t(age,ai,zi,lambdai,ei,xi).ge.1.0_dp) then 
				print*, ' '
				print*, 'Hours highers than 1', age,ai,zi,lambdai,ei,xi
				print*, Cons_t(age,ai,zi,lambdai,ei,xi)
				STOP
			endif 


			! Savings 
				Aprime_t(age, ai, zi, lambdai,ei,xi) = YGRID_t(ai,zi,xi)  + Y_h(Hours_t(age, ai, zi, lambdai,ei,xi),age,lambdai,ei,wage)  & 
		                    					& - Cons_t(age, ai, zi, lambdai,ei,xi) 
		                    
		    If (Aprime_t(age, ai, zi, lambdai,ei,xi)  .lt. amin) then

		    	! print*, ' Aprime was below minimum!!!!'
            	Aprime_t(age, ai, zi, lambdai,ei,xi) = amin
		         
	           	! Compute hours
		        if (NSU_Switch.and.(Progressive_Tax_Switch.eqv..false.)) then 
		        	Hours_t(age, ai, zi, lambdai,ei,xi)  = max( 0.0_dp , &
		        			&  gamma - (1.0_DP-gamma) *( YGRID_t(ai,zi,xi) - Aprime_t(age, ai, zi, lambdai,ei,xi)) / (psi*yh(age, lambdai,ei)) )
                else               	        
                	!compute  hours using FOC_HA                              
		        	par_FOC(1:6) = (/age,ai,zi,lambdai,ei,xi/) 
		        	par_FOC(7)   = amin
		        	brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_HA, brent_tol, Hours_t(age, ai, zi, lambdai,ei,xi), par_FOC)
		        end if  

	            Cons_t(age,ai,zi,lambdai,ei,xi) = YGRID_t(ai,zi,xi)+  Y_h(Hours_t(age, ai, zi, lambdai,ei,xi),age,lambdai,ei,wage)  &
		                            		  & - Aprime_t(age, ai, zi, lambdai,ei,xi)      

				IF (Cons_t(age, ai, zi, lambdai,ei,xi) .le. 0.0_DP)  THEN
					print*,'w1: Cons(age, ai, zi, lambdai,ei,xi)=',Cons_t(age, ai, zi, lambdai,ei,xi)
					STOP
				ENDIF 
				if (Hours_t(age,ai,zi,lambdai,ei,xi).ge.1.0_dp) then 
					print*, ' '
					print*, 'w1 Hours highers than 1', age,ai,zi,lambdai,ei,xi
					print*, Cons_t(age,ai,zi,lambdai,ei,xi)
					STOP
				endif                   
		     endif  
		    !$omp critical
		    IF (Cons_t(age, ai, zi, lambdai,ei,xi) .le. 0.0_DP)  THEN
				print*,'w1: Cons(age, ai, zi, lambdai,ei,xi)=',Cons_t(age, ai, zi, lambdai,ei,xi)
				STOP
			ENDIF 
			if (Hours_t(age,ai,zi,lambdai,ei,xi).ge.1.0_dp) then 
				print*, ' '
				print*, 'w1 Hours highers than 1', age,ai,zi,lambdai,ei,xi
				print*, Cons_t(age,ai,zi,lambdai,ei,xi)
				STOP
			endif    
			!$omp end critical
		ENDDO ! ai   

		!$omp critical
		if (any(isnan(Cons_t(age,:,zi,lambdai,ei,xi)))) then 
			print*, "isnan - Consumption working 3"
			print*, age,lambdai,ai,zi,ei,xi
			print*, 'Cons'
			print*, Cons_t(age,:,zi,lambdai,ei,xi)
			print*, 'Hours'
			print*, Hours_t(age,:,zi,lambdai,ei,xi) 
			STOP 
		end if         
		!$omp end critical        

		ai=1           
        DO WHILE ( YGRID_t(ai,zi,xi) .lt. EndoYgrid(1) )
        	! print*, ' Extrapolation between YGRID and EndoYgrid!!!!'
	        ! Solve for the Euler equation directly
	        state_FOC  = (/age,ai,zi,lambdai,ei,xi/)
			brentvalue = brent_p( min(amin,YGRID_t(ai,zi,xi))   ,  (amin+YGRID_t(ai,zi,xi))/2.0_DP  ,  &
                             	& YGRID_t(ai,zi,xi)  + psi * ( yh(age, lambdai,ei)*0.95_DP )**(1.0_DP-tauPL)  ,  &
	                             & FOC_WH, brent_tol, Aprime_t(age, ai, zi, lambdai,ei,xi) , state_FOC )

			! Compute hours
	        if (NSU_Switch.and.(Progressive_Tax_Switch.eqv..false.)) then 
	        	Hours_t(age, ai, zi, lambdai,ei,xi)  = max( 0.0_dp , &
	        		&  gamma - (1.0_DP-gamma) *( YGRID_t(ai,zi,xi) - Aprime_t(age, ai, zi, lambdai,ei,xi)) / (psi*yh(age, lambdai,ei)) )
            else               	        
				!compute  hours using FOC_HA                              
				par_FOC(1:6) = (/age,ai,zi,lambdai,ei,xi/)
				par_FOC(7)   = Aprime_t(age, ai, zi, lambdai,ei,xi)
	            brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_HA, brent_tol, Hours_t(age, ai, zi, lambdai,ei,xi), par_FOC) 
 			end if 

            Cons_t(age, ai, zi, lambdai,ei,xi)=  YGRID_t(ai,zi,xi) + Y_h(Hours_t(age, ai, zi, lambdai,ei,xi),age,lambdai,ei,wage)  &
                           						 & - Aprime_t(age, ai, zi, lambdai,ei,xi)
           IF (Cons_t(age, ai, zi, lambdai,ei,xi) .le. 0.0_DP)  THEN
                print*,'w2:',age,zi,lambdai,ei,ai,xi, 'Cons(age, ai, zi, lambdai,ei,xi)=',Cons_t(age, ai, zi, lambdai,ei,xi), &
                    & 'Aprime(age, ai, zi, lambdai,ei,xi)=',Aprime_t(age, ai, zi, lambdai,ei,xi), &
                    & 'Hours(age, ai, zi, lambdai,ei,xi)=', Hours_t(age, ai, zi, lambdai,ei,xi), &
                    & 'yh(age, lambdai,ei)=', yh(age, lambdai,ei),'YGRID(ai,zi)=',YGRID_t(ai,zi,xi)
                !pause
                print*, "there is a problem in line 2063"
                STOP
            ENDIF 

            if (Hours_t(age,ai,zi,lambdai,ei,xi).ge.1.0_dp) then 
				print*, ' '
				print*, 'w2 Hours highers than 1', age,ai,zi,lambdai,ei,xi
				print*, Cons_t(age,ai,zi,lambdai,ei,xi)
				STOP
			endif 

            ai = ai + 1
        ENDDO  

	    if (any(isnan(Cons_t(age,:,zi,lambdai,ei,xi)))) then 
			print*, "isnan - Consumption working 2"
			print*, age,lambdai,ai,zi,ei,xi
			print*, 'Cons'
			print*, Cons_t(age,:,zi,lambdai,ei,xi)
			print*, 'Hours'
			print*, Hours_t(age,:,zi,lambdai,ei,xi) 
			STOP 
		end if 

	                 
    ENDDO !ei         
    ENDDO !lambdai
    ENDDO ! xi 
	ENDDO !zi
    ENDDO !age


	! Interpolate to get values of policy functions on agrid (note that policy functions are defined on agrid_t)
	
	if (Y_a_threshold.eq.0.0_dp) then 
		Cons   = Cons_t
		Hours  = Hours_t
		Aprime = Aprime_t
	else 
		DO xi=1,nx
		DO age=1,MaxAge
	    DO lambdai=1,nlambda
	    DO zi=1,nz
	    DO ei=1,ne	                
	    DO ai=1,na
			Cons(age,ai,zi,lambdai,ei,xi)   = Linear_Int(YGRID_t(:,zi,xi) , Cons_t(age,:,zi,lambdai,ei,xi)   , na_t , YGRID(ai,zi,xi))
	    	Hours(age,ai,zi,lambdai,ei,xi)  = Linear_Int(YGRID_t(:,zi,xi) , Hours_t(age,:,zi,lambdai,ei,xi)  , na_t , YGRID(ai,zi,xi))
	    	Aprime(age,ai,zi,lambdai,ei,xi) = Linear_Int(YGRID_t(:,zi,xi) , Aprime_t(age,:,zi,lambdai,ei,xi) , na_t , YGRID(ai,zi,xi))
		ENDDO !ai         
		ENDDO !ei         
	    ENDDO !zi
	    ENDDO !lambdai
		ENDDO !age
		ENDDO !xi
	end if 

	! Deallocate policy functions on adjusted grid (so that they can be allocated later)
	deallocate( YGRID_t  )
	deallocate( MBGRID_t ) 
	deallocate( Cons_t   )
	deallocate( Hours_t  )
	deallocate( Aprime_t )

	! Adjust consumption by taxes
	Cons = Cons/(1.0_DP+tauC)

	if (any(isnan(Cons))) then 
		print*, "isnan - Consumption"
		STOP 
	end if 
	if (any(isnan(Hours))) then 
		print*, "isnan - Hours"
		STOP 
	end if 
	if (any(isnan(Aprime))) then 
		print*, "isnan - Hours"
		STOP 
	end if 

! 	print*, 'Policy Functions'
! 	print*, ' YGRID ',' Cons ',' Hours ',' Aprime ',' A '
! 	do ai=1,na
! 		! print*, YGRID(ai,4), Cons(40,ai,4,3,3), Hours(40,ai,4,3,3), Aprime(40,ai,4,3,3), agrid(ai),&
! 		! 		& Y_h(Hours(40,ai,4,3,3),40,3,3,wage),&
! 		! 		& (1+tauC)*Cons(40,ai,4,3,3) + Aprime(40,ai,4,3,3) - YGRID(ai,4) -Y_h(Hours(40,ai,4,3,3),40,3,3,wage)
! 	end do 
! 	print*, ' '

END SUBROUTINE EGM_RETIREMENT_WORKING_PERIOD	


!========================================================================================
!========================================================================================
!========================================================================================
! EGM_Working_Period: Solves for the endogenous grid during the working period
!
! Usage: call EGM_Working_Period(MB_in,H_min,C_endo,H_endo,Y_endo)
!
! Input: MB_in   , real(dp), marginal benefit from assets in t+1
!        H_min   , real(dp), minimum value of hours
!
! Implicit inputs: ai     , integer , index of current level of assets
!				   zi     , integer , index of agent's permanent entreprenurial ability
!				   lambdai, integer , index of agent's permanent labor income component
!				   ei     , integer , index of agent's transitory labor income component
!				   age    , integer , agent's age
!
! Output: C_endo , real(dp), Endogenous value of consumption in t
! 		  H_endo , real(dp), Endogenous value of hours in t
! 		  Y_endo , real(dp), Endogenous value of Y grid in t
!
Subroutine EGM_Working_Period(MB_in,H_min,state_FOC,C_endo,H_endo,Y_endo)
	Implicit None 
	real(dp), intent(in)   	  :: MB_in(nx), H_min, state_FOC(6)
	real(dp), intent(out)  	  :: C_endo, H_endo, Y_endo
	real(dp)               	  :: C_euler, C_foc, brentvalue, E_MU_cp(nx), MB_a_vec(nx)
	REAL(DP), DIMENSION(6+nx) :: par_FOC
	integer                	  :: age, ai, zi, lambdai, ei, xi, xp_ind

	! Allocate state and get indeces
	age 	= int(state_FOC(1))
	ai   	= int(state_FOC(2))
	zi      = int(state_FOC(3))
	lambdai = int(state_FOC(4))
	ei      = int(state_FOC(5))
	xi      = int(state_FOC(6))

	par_FOC(1:6) = state_FOC

	if (Progressive_Tax_Switch.eqv..true.) then !Progressive labor taxes 
		if (NSU_Switch.eqv..true.) then 
			! Non-Separable Utility
			do xp_ind=1,nx
				E_MU_cp(xp_ind) = SUM(pr_e(ei,:) * &
			        & Cons_t(age+1,ai,zi,lambdai,:,xp_ind)**((1.0_dp-sigma)*gamma-1.0_dp)    * &
			        & (1.0_dp-Hours_t(age+1,ai,zi,lambdai,:,xp_ind))**((1.0_dp-sigma)*(1.0_dp-gamma))  )
			enddo

			C_euler = (beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_MU_cp) & 
					  & + (1.0_dp-survP(age))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma) &
					  &   *chi_bq*((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ) &
					  & **(1.0_dp/((1.0_dp-sigma)*gamma-1.0_dp))
			C_foc   = (gamma/(1.0_dp-gamma))*(1.0_dp-H_min)*MB_h(H_min,age,lambdai,ei,wage)

			if (C_euler.ge.C_foc) then
				C_endo  = C_euler 
				H_endo = 0.0_dp
			else 
				if (Log_Switch.eqv..true.) then
		  			! Auxiliary consumption variable for FOC_H        
				    par_FOC(7) = C_euler
				    ! Solution for hours from Euler or FOC  equation according to sigma 
				    brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_H, brent_tol, H_endo , par_FOC ) 
				    C_endo = C_euler
				else 
				    ! Set Marginal benefit of assets to the below threshold level
				    par_FOC(7:6+nx) = MB_in
				    ! Solution for hours from Euler or FOC  equation according to sigma 
				    brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_H_NSU, brent_tol, H_endo , par_FOC ) 
				    ! Implied consumption by hours from Labor FOC
				    C_endo = (gamma/(1.0_dp-gamma))*(1.0_dp-H_endo)*MB_h(H_endo,age,lambdai,ei,wage)
				end if 
			end if 

		else 
			! Separable Utility
			do xp_ind=1,nx
				E_MU_cp(xp_ind) = SUM(pr_e(ei,:) * Cons_t(age+1,ai,zi,lambdai,:,xp_ind)**(-sigma) )
			enddo
			C_endo = 1.0_dp/( (beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_mu_cp) &
					&   + (1.0_dp-survP(age))*(1.0_dp+tauC)**(1.0_dp-sigma)&
					&   *chi_bq/((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**sigma ) ) **(1.0_dp/sigma) 
			C_foc  = (MB_h(H_min,age,lambdai,ei,wage)*(1.0_dp-H_min)**(gamma)/phi)**(1.0_dp/sigma)

			if (C_endo.ge.C_foc) then
			  H_endo = 0.0_dp
			else 
			  ! Auxiliary consumption variable for FOC_H        
			  par_FOC(7) = C_endo
			  ! Solution for hours from Euler or FOC  equation according to sigma 
			  brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_H, brent_tol, H_endo , par_FOC ) 
			end if  

			end if 

	else ! Linear labor taxes 
		if (NSU_Switch.eqv..true.) then 
			! Non-Separable Utility
			do xp_ind=1,nx
				E_MU_cp(xp_ind) = sum( pr_e(ei,:) * (Cons_t(age+1,ai,zi,lambdai,:,xp_ind)**(gamma*(1.0_DP-sigma)-1.0_DP)) &
			    & *  ( (1.0_DP-Hours_t(age+1, ai, zi, lambdai,:,xp_ind))**((1.0_DP-gamma)*(1.0_DP-sigma))))
			enddo
			  C_endo = ((gamma*psi*yh(age, lambdai,ei)/(1.0_DP-gamma))**((1.0_DP-gamma)*(1.0_DP-sigma)) &
			    & *  (beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_MU_cp) &
			    & + (1.0_dp-survP(age))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
			    & *   chi_bq*((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ) )**(-1.0_DP/sigma)

			  H_endo = 1.0_DP - (1.0_DP-gamma)*C_endo/(gamma*psi*yh(age,lambdai,ei))   

			If (H_endo .lt. 0.0_DP) then
			    H_endo = 0.0_DP 
			    C_endo  = ( beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_MU_cp) &
			    		& + (1.0_dp-survP(age))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
			    		& *chi_bq*((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp))&
			    		& **(1.0_DP/(gamma*(1.0_DP-sigma)-1.0_DP))
			endif 

			! print*,' '
			! print*,' 	Inside EGM'
			! print*,' 	Consumption_t+1=',Cons_t(age+1,ai,zi,lambdai,:,:)
		else 
			! Separable Utility
			do xp_ind=1,nx
				E_MU_cp(xp_ind) = sum( pr_e(ei,:) * (Cons_t(age+1,ai,zi,lambdai,:,xp_ind)**(-sigma)) )
			enddo
			C_endo  = 1.0_DP/( beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_MU_cp) &
				&   + (1.0_dp-survP(age))*(1.0_dp+tauC)**(1.0_dp-sigma)&
				&   * chi_bq/((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**sigma ) **(1.0_dp/sigma) 

			H_endo = max(0.0_DP , 1.0_DP - (phi*C_endo**sigma/(psi*yh(age, lambdai,ei)))**(1.0_dp/gamma) )  

		end if 
	end if


	! Endogenous grid for asset income
	Y_endo = agrid_t(ai) + C_endo - Y_h(H_endo,age,lambdai,ei,wage)

	! 		!$omp critical
	! 		if ((zi.eq.1).and.(ai.eq.16)) then 
	! 		print*, "EGM Working Periods"
	! 		print*, C_endo,H_endo,Y_endo
	! 		print*, MB_in,state_FOC
	! 		print*, ' '
	! 		endif 
	! 		!$omp end critical
end Subroutine EGM_Working_Period

!========================================================================================
!========================================================================================
!========================================================================================
! Utility: Evaluate the utility at given consumption and hours
!
! Usage: U = Utility(c,h)
!
! Input: c   , real(dp), value of consumption
!		 h   , real(dp), value of hours
!
! Output: U , real(dp), utility
!
	FUNCTION Utility(c,h)
		IMPLICIT NONE   
		real(DP), intent(in) :: c,h
		real(DP)   			 :: Utility

		if (NSU_Switch.eqv..true.) then 
			if (Log_Switch.eqv..true.) then 
				Utility = log(c) + (1.0_dp-gamma)/gamma*log(1.0_dp-h) 
			else 
				Utility = ( c**gamma * (1.0_dp-h)**(1.0_dp-gamma) )**(1.0_dp-sigma) / (1.0_dp-sigma)
			end if 
		else 
			Utility = U_C(c) + U_H(h)
		end if 

	END FUNCTION Utility

	FUNCTION U_C(c)
		IMPLICIT NONE   
		real(DP), intent(in) :: c
		real(DP)   			 :: U_C

		if ((Log_Switch.eqv..true.).or.(sigma.eq.1.0_dp)) then
			U_C = log(c)
		else 
			U_C = c**(1.0_dp-sigma)/(1.0_dp-sigma)
		end if 

	END FUNCTION U_C

	FUNCTION U_H(h)
		IMPLICIT NONE   
		real(DP), intent(in) :: h
		real(DP)   			 :: U_H

		if ((Log_Switch.eqv..true.).or.(gamma.eq.1.0_dp)) then
			U_H = phi*log(1.0_dp-h)
		else 
			U_H = phi*(1.0_dp-h)**(1.0_dp-gamma)/(1.0_dp-gamma)
		end if 

	END FUNCTION U_H

	FUNCTION v_bq(a) 
		IMPLICIT NONE
		real(DP), intent(in) :: a
		real(DP)             :: v_bq

		if (chi_bq.gt.0.0_dp) then 
		v_bq = chi_u*((1.0_dp-tau_bq)*a+bq_0)**(gamma*(1.0_dp-sigma))/(1.0_dp-sigma) 
		else 
		v_bq = 0.0_dp 
		endif 


	END FUNCTION v_bq


!========================================================================================
!========================================================================================
!========================================================================================


! SUBROUTINE COMPUTE_VALUE_FUNCTION_SPLINE()
! 	IMPLICIT NONE
! 	REAL(DP), DIMENSION(na) :: ValueP1, ValueP2, ValueP, ExpValueP
! 	REAL(DP), DIMENSION(nx) :: ExpValueP_e, D_V_1, D_V_na
! 	INTEGER                 :: xp_ind

! 	! Announce method of interpolation
! 	! Interpolation is used to get value function at optimal a' (when computing expectations)
! 	!print*,'VALUE FUNCTION SPLINE'

! 	! Compute value function at grid nodes
! 	! Final age
! 	age=MaxAge
! 	DO xi=1,xi
! 	DO ai=1,na    
!     DO zi=1,nz
!     DO lambdai=1,nlambda          
! 	DO ei=1,ne
!       	ValueFunction(age, ai, zi, lambdai, ei,xi) = Utility(Cons(age, ai, zi, lambdai, ei,xi),Hours(age, ai, zi, lambdai, ei,xi))
! 		! print*,Cons(age, ai, zi, lambdai, ei),  ValueFunction(age, ai, zi, lambdai, ei) 
! 		! pause
! 	ENDDO ! ei          
!     ENDDO ! lambdai
!     ENDDO ! zi
! 	ENDDO ! ai
! 	ENDDO ! xi

! 	! Retirement Period
! 	DO xi=1,nx
! 	DO age=MaxAge-1,RetAge,-1
! 	DO zi=1,nz
! 	DO lambdai=1,nlambda          
! 	DO ei=1,ne   
! 		DO ai=1,na  
! 			ExpValueP(ai) = sum(ValueFunction(age+1, ai, zi, lambdai, ei,:) * pr_x(xi,:,zi,age))
!         ENDDO    

!   		if (NSU_Switch.eqv..true.) then
!             CALL spline( agrid, ValueFunction(age+1,  :, zi, lambdai, ei, xi) , na , &
!       		& sum(pr_x(xi,:,zi,age)*gamma*MBGRID(1 ,zi,:)*Cons(age+1, 1,zi,lambdai,ei,:)**((1.0_DP-sigma)*gamma-1.0_DP)/(1_DP+tauC)), &
!         	& sum(pr_x(xi,:,zi,age)*gamma*MBGRID(na,zi,:)*Cons(age+1,na,zi,lambdai,ei,:)**((1.0_DP-sigma)*gamma-1.0_DP)/(1_DP+tauC)), &
!         	& ValueP2)  
!         else 
!         	CALL spline( agrid, ValueFunction(age+1,  :, zi, lambdai, ei, xi) , na , &
!       			& sum(pr_x(xi,:,zi,age)* gamma*MBGRID(1,zi,:) /Cons(age+1,  1, zi, lambdai, ei, :) **(sigma)/(1_DP+tauC)), &
!         		& sum(pr_x(xi,:,zi,age)* gamma*MBGRID(na,zi,:)/Cons(age+1, na, zi, lambdai, ei, :)**(sigma)/(1_DP+tauC)) , ValueP2)  
!       	end if 
	                  
!         DO ai=1,na    
!             call splint( agrid, ValueFunction(age+1, :, zi, lambdai, ei, xi), &
!                     & ValueP2, na, Aprime(age,ai,zi,lambdai, ei, xi), ValueP(ai))  

!             ValueFunction(age,ai,zi,lambdai,ei,xi) = Utility(Cons(age,ai,zi,lambdai,ei,xi),Hours(age,ai,zi,lambdai,ei,xi)) &
!               										& + beta*survP(age)* ValueP(ai)      
!         ENDDO ! ai
	              
!     ENDDO ! ei          
! 	ENDDO ! lambdai
!     ENDDO ! zi
! 	ENDDO ! age
! 	ENDDO ! xi
   

! 	! Working Period
! 	DO xi=1,nx 
! 	DO age=RetAge-1,1,-1
!     DO zi=1,nz
!     DO lambdai=1,nlambda          
! 	DO ei=1,ne
!         DO ai=1,na
!         	DO xp_ind=1,nx 
!         		ExpValueP_e(xp_ind) = sum(ValueFunction(age+1, ai, zi, lambdai, :, xp_ind) * pr_e(ei,:))
!         	ENDDO    
!               ExpValueP(ai) = sum(ExpValueP_e * pr_x(xi,:,zi,age))
!         ENDDO

!         if (NSU_Switch.eqv..true.) then
! 	        DO xp_ind=1,nx 
! 	        	D_V_1(xp_ind)  = (gamma*MBGRID(1,zi,xp_ind)/(1.0_DP+tauC)) * sum(pr_e(ei,:)* &
! 			            & Cons(age+1, 1, zi, lambdai, :,xp_ind)**((1.0_DP-sigma)*gamma-1.0_DP) * &
! 			            & (1.0_DP-Hours(age+1,1,zi,lambdai,:,xp_ind))**((1.0_DP-gamma)*(1.0_DP-sigma)))
! 	        	D_V_na(xp_ind) = (gamma*MBGRID(na,zi,xp_ind)/(1.0_DP+tauC)) * sum(pr_e(ei,:)* &
! 			            & Cons(age+1, na, zi, lambdai, :, xp_ind)**((1.0_DP-sigma)*gamma-1.0_DP) * &
! 			            & (1.0_DP-Hours(age+1,na,zi,lambdai,:,xp_ind))**((1.0_DP-gamma)*(1.0_DP-sigma)))
! 	        ENDDO
! 	    else 
! 	    	DO xp_ind=1,nx 
! 	        	D_V_1(xp_ind)  = (gamma*MBGRID(1 ,zi,xp_ind)/(1.0_DP+tauC))*&
! 	        			& sum(pr_e(ei,:) / Cons(age+1,  1, zi, lambdai, :,xp_ind)**(sigma) )
! 	        	D_V_na(xp_ind) = (gamma*MBGRID(na,zi,xp_ind)/(1.0_DP+tauC))*&
! 	        			& sum(pr_e(ei,:) / Cons(age+1, na, zi, lambdai, :,xp_ind)**(sigma) )
! 	        ENDDO
! 	    endif 

!     	CALL spline( agrid, ExpValueP , na , sum(D_V_1  * pr_x(xi,:,zi,age)), sum(D_V_na * pr_x(xi,:,zi,age)), ValueP2)   

!         DO ai=1,na 
!         	call splint( agrid, ExpValueP, ValueP2, na, Aprime(age,ai,zi,lambdai, ei, xi), ValueP(ai))   
        
!         	ValueFunction(age, ai, zi, lambdai, ei, xi) = Utility(Cons(age,ai,zi,lambdai,ei,xi),Hours(age,ai,zi,lambdai,ei,xi)) &
!         	                                         & + beta*survP(age)*ValueP(ai)
! 		ENDDO ! ai
! 	ENDDO ! ei          
!     ENDDO ! lambdai
!     ENDDO ! zi
! 	ENDDO ! age
! 	ENDDO ! xi

! END SUBROUTINE COMPUTE_VALUE_FUNCTION_SPLINE 


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE COMPUTE_VALUE_FUNCTION_LINEAR(Cons_mat,Hours_mat,Aprime_mat,Value_mat,Bq_Value_mat)
	IMPLICIT NONE
	REAL(DP), DIMENSION(MaxAge,na,nz,nlambda,ne,nx), INTENT(in)  :: Cons_mat, Hours_mat, Aprime_mat
	REAL(DP), DIMENSION(MaxAge,na,nz,nlambda,ne,nx), INTENT(out) :: Value_mat, Bq_Value_mat
	INTEGER  :: tklo, tkhi, xp_ind, age, xi, ai, zi, lambdai, ei
	REAL(DP) :: PrAprimelo, PrAprimehi, E_MU_cp(nx), Bq_E_MU_cp(nx), aux

	! print*,'VALUE FUNCTION LINEAR'

	age=MaxAge
	DO xi=1,nx
	DO ai=1,na
    DO zi=1,nz
    DO lambdai=1,nlambda          
	DO ei=1,ne
      	   Value_mat(age, ai, zi, lambdai, ei, xi) = &
      				& Utility(Cons_mat(age, ai, zi, lambdai, ei, xi),Hours_mat(age, ai, zi, lambdai, ei, xi))
		Bq_Value_mat(age, ai, zi, lambdai, ei, xi) = v_bq(Aprime_mat(age, ai, zi, lambdai, ei, xi))
      	! Value_mat(age, ai, zi, lambdai, ei, xi) = ((Cons_mat(age,ai,zi,lambdai,ei,xi)**gamma) &
       !             & * (1.0_DP-Hours_mat(age,ai,zi,lambdai,ei,xi))**(1.0_DP-gamma))**(1.0_DP-sigma)/(1.0_DP-sigma)  
		! print*,Cons_mat(age, ai, zi, lambdai, ei),  Value_mat(age, ai, zi, lambdai, ei) 
		! pause
	ENDDO ! ei          
    ENDDO ! lambdai
    ENDDO ! zi
	ENDDO ! ai
	ENDDO ! xi

	aux = 0.0_dp
	! Retirement Period
	DO age=MaxAge-1,RetAge,-1
	DO xi=1,nx
	DO ai=1,na    
    DO zi=1,nz
    DO lambdai=1,nlambda          
	DO ei=1,ne
		if ( Aprime_mat(age,ai,zi,lambdai, ei,xi) .ge. amax) then
		    tklo =na-1
		else if (Aprime_mat(age,ai,zi,lambdai, ei, xi) .lt. amin) then
		    tklo = 1
		else
		    tklo = ((Aprime_mat(age,ai,zi,lambdai, ei, xi) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
		end if            
		tkhi = tklo + 1        
		PrAprimelo = ( agrid(tkhi) - Aprime_mat(age,ai,zi,lambdai, ei, xi) ) / ( agrid(tkhi) -agrid(tklo) )
		PrAprimehi = ( Aprime_mat(age,ai,zi,lambdai, ei, xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )        
		PrAprimelo = min(PrAprimelo, 1.0_DP)
		PrAprimelo = max(PrAprimelo, 0.0_DP)
		PrAprimehi = min(PrAprimehi, 1.0_DP)
		PrAprimehi = max(PrAprimehi, 0.0_DP)    

		   Value_mat(age, ai, zi, lambdai, ei, xi) = Utility(Cons_mat(age,ai,zi,lambdai,ei,xi),Hours_mat(age,ai,zi,lambdai,ei,xi)) &
			  & + beta*survP(age)* sum( pr_x(xi,:,zi,age)* (PrAprimelo*Value_mat(age+1, tklo, zi, lambdai, ei, :) &
			  & 				                     	 +  PrAprimehi*Value_mat(age+1, tkhi, zi, lambdai, ei, :)) ) 

		Bq_Value_mat(age, ai, zi, lambdai, ei, xi) = (1.0_dp-survP(age))*v_bq(Aprime_mat(age,ai,zi,lambdai,ei,xi)) &
			  & + beta*survP(age)* sum( pr_x(xi,:,zi,age)* (PrAprimelo*Bq_Value_mat(age+1, tklo, zi, lambdai, ei, :) &
			  & 				                     	 +  PrAprimehi*Bq_Value_mat(age+1, tkhi, zi, lambdai, ei, :)) ) 
        ! Value_mat(age, ai, zi, lambdai, ei, xi) = ((Cons_mat(age,ai,zi,lambdai,ei,xi)**gamma) &
        !               & * (1.0_DP-Hours_mat(age,ai,zi,lambdai,ei,xi))**(1.0_DP-gamma))**(1.0_DP-sigma)/(1.0_DP-sigma) &
        !               & + beta*survP(age)* sum(pr_x(xi,:,zi) * (PrAprimelo*Value_mat(age+1, tklo, zi, lambdai, ei, :)&
        !               & +                   PrAprimehi*Value_mat(age+1, tkhi, zi, lambdai, ei, :)) )
	ENDDO ! ei          
    ENDDO ! lambdai
    ENDDO ! zi
	ENDDO ! ai
	ENDDO ! age
	ENDDO ! xi
	!print*,Value_mat


	! Working Period
	DO age=RetAge-1,1,-1
	DO xi=1,nx
	DO ai=1,na    
    DO zi=1,nz
    DO lambdai=1,nlambda          
	DO ei=1,ne
		if ( Aprime_mat(age,ai,zi,lambdai,ei,xi) .ge. amax) then
		    tklo =na-1
	    elseif (Aprime_mat(age,ai,zi,lambdai,ei,xi) .lt. amin) then
	        tklo = 1
        else
            tklo = ((Aprime_mat(age,ai,zi,lambdai,ei,xi) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
		endif  

		tkhi = tklo + 1        
		PrAprimelo = ( agrid(tkhi) - Aprime_mat(age,ai,zi,lambdai, ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
		PrAprimehi = ( Aprime_mat(age,ai,zi,lambdai, ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )        
		PrAprimelo = min(PrAprimelo, 1.0_DP)
		PrAprimelo = max(PrAprimelo, 0.0_DP)
		PrAprimehi = min(PrAprimehi, 1.0_DP)
		PrAprimehi = max(PrAprimehi, 0.0_DP) 


		do xp_ind=1,nx
			E_MU_cp(xp_ind)    = sum( ( PrAprimelo *    Value_mat(age+1, tklo, zi, lambdai,:,xp_ind)  &
		   & 					   +    PrAprimehi *    Value_mat(age+1, tkhi, zi, lambdai,:,xp_ind)) * pr_e(ei,:))
			Bq_E_MU_cp(xp_ind) = sum( ( PrAprimelo * Bq_Value_mat(age+1, tklo, zi, lambdai,:,xp_ind)  &
		   & 					   +    PrAprimehi * Bq_Value_mat(age+1, tkhi, zi, lambdai,:,xp_ind)) * pr_e(ei,:))
		enddo    

		   Value_mat(age, ai, zi, lambdai, ei, xi) = Utility(Cons_mat(age,ai,zi,lambdai,ei,xi),Hours_mat(age,ai,zi,lambdai,ei,xi))  &
		   											& + beta*survP(age)* sum( pr_x(xi,:,zi,age)*E_mu_cp )
	   	Bq_Value_mat(age, ai, zi, lambdai, ei, xi) = (1.0_dp-survP(age))*v_bq(Aprime_mat(age,ai,zi,lambdai,ei,xi)) &
		   											& + beta*survP(age)* sum( pr_x(xi,:,zi,age)*Bq_E_mu_cp )
		! Value_mat(age, ai, zi, lambdai, ei) = ((Cons_mat(age,ai,zi,lambdai,ei,xi)**gamma) &
  !                      & * (1.0_DP-Hours_mat(age,ai,zi,lambdai,ei,xi))**(1.0_DP-gamma))**(1.0_DP-sigma)/(1.0_DP-sigma) &
  !                      & + beta*survP(age)* sum( pr_x(xi,:,zi)*E_mu_cp )
		! if ( Value_mat(age, ai, zi, lambdai, ei) .lt. (-100.0_DP) ) then
		!    print*,'Value_mat(age, ai, zi, lambdai, ei)=',Value_mat(age, ai, zi, lambdai, ei)
		! endif
	ENDDO ! ei          
    ENDDO ! lambdai
    ENDDO ! zi
	ENDDO ! ai
	ENDDO ! age
	ENDDO ! xi

	! Add bequest value to consumption/leisure value t get total value
	if (chi_bq.eq.0.0_dp) then 
		Bq_Value_mat = 0.0_dp  
	endif 
	Value_mat = Value_mat + Bq_Value_mat

END SUBROUTINE COMPUTE_VALUE_FUNCTION_LINEAR 


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE COMPUTE_VALUE_FUNCTION_TRANSITION
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: Cons_mat, Hours_mat, Aprime_mat, Value_mat, Bq_Value_mat
	INTEGER :: aux_T

	allocate( Cons_mat(    MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Hours_mat(   MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Aprime_mat(  MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Value_mat(   MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Bq_Value_mat(MaxAge,na,nz,nlambda,ne,nx) )

	print*,'Computing Value Function'
	print*,'	Value functions computed only for agents alive at the time of policy change'
	print*,' 	Last index corresponds to cohort (age at time of policy change), not transition period'
	do age=1,MaxAge
		! print*,'		Cohort of age',age,'at time of policy change'

		! Set policy functions for cohort of age "age" at time of policy change
			! Previous periods of life from benchmark
			if (age>1) then 
			Cons_mat(1:age-1,:,:,:,:,:)   = Cons_bench(1:age-1,:,:,:,:,:)
			Hours_mat(1:age-1,:,:,:,:,:)  = Hours_bench(1:age-1,:,:,:,:,:)
			Aprime_mat(1:age-1,:,:,:,:,:) = Aprime_bench(1:age-1,:,:,:,:,:)
			endif 
			! Remaining periods of life from transition matrices
			aux_T = max(MaxAge-age+1,T+1)
			do ti=1,aux_T
			Cons_mat(age+ti-1,:,:,:,:,:)   = Cons_tr(age+ti-1,:,:,:,:,:,ti)
			Hours_mat(age+ti-1,:,:,:,:,:)  = Hours_tr(age+ti-1,:,:,:,:,:,ti)
			Aprime_mat(age+ti-1,:,:,:,:,:) = Aprime_tr(age+ti-1,:,:,:,:,:,ti)
			enddo 
			! Supplement with second steady state if needed
			if (aux_T.lt.(MaxAge-age+1)) then
			Cons_mat(T+2:,:,:,:,:,:)   = Cons_exp(T+2:,:,:,:,:,:)
			Hours_mat(T+2:,:,:,:,:,:)  = Hours_exp(T+2:,:,:,:,:,:)
			Aprime_mat(T+2:,:,:,:,:,:) = Aprime_exp(T+2:,:,:,:,:,:)
			endif 
		! Obtain value function 
			call COMPUTE_VALUE_FUNCTION_LINEAR(Cons_mat,Hours_mat,Aprime_mat,Value_mat,Bq_Value_mat)

		! Save value function
			ValueFunction_tr(:,:,:,:,:,:,age) = Value_mat
		! Save bequest value function
			Bq_Value_tr(:,:,:,:,:,:,age) = Bq_Value_mat
	enddo 
	print*,'Computing Value Function Completed'

END SUBROUTINE COMPUTE_VALUE_FUNCTION_TRANSITION




!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE GOVNT_BUDGET(print_flag)
	IMPLICIT NONE
	LOGICAL, INTENT(IN) :: print_flag
	real(dp) :: test 

	! Initialize variables
	GBAR 		 = 0.0_DP
	GBAR_K 		 = 0.0_DP
	GBAR_W		 = 0.0_DP
	GBAR_L 		 = 0.0_DP
	GBAR_C       = 0.0_DP
	GBAR_BQ      = 0.0_DP
	SSC_Payments = 0.0_DP
	Tot_Lab_Inc  = 0.0_DP
	Tot_Cap_Inc  = 0.0_DP
	test         = 0.0_DP

	! Compute total expenditure = total revenue
		! For each state accumulate tax income weighted by equilibrium distribution
		! Taxes come from capital income, wealth, labor income and consumption
		! G_L is the tax revenue that comes exclusively from labor income
	DO xi=1,nx
	DO age=1, MaxAge
	DO ai=1,na
	DO zi=1,nz
	DO lambdai=1,nlambda
	DO ei=1,ne
	    GBAR = GBAR + DBN1(age,ai,zi,lambdai,ei,xi) * ( tauK*( R_z(zi)*agrid(ai) + Pr_mat(ai,zi,xi) )  	    &
	          & + ( agrid(ai) + ( R_z(zi)*agrid(ai) + Pr_mat(ai,zi,xi) ) *(1.0_DP-tauK)  ) - YGRID(ai,zi,xi)  &	
	          & + yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei,xi)  									&
	          & - psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei,xi))**(1.0_DP-tauPL)  			&
	          & + tauC  *  cons(age,ai,zi,lambdai,ei,xi)  										    & 
	          & + tau_bq*aprime(age,ai,zi,lambdai,ei,xi)*(1.0_DP-survP(age)) )


	    GBAR_L = GBAR_L  + DBN1(age,ai,zi,lambdai,ei,xi) * (  yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei,xi) &
	               &- psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei,xi))**(1.0_DP-tauPL) )

	    GBAR_K = GBAR_K + DBN1(age,ai,zi,lambdai,ei,xi) * tauK*( R_z(zi)*agrid(ai) + Pr_mat(ai,zi,xi) )

	    GBAR_W = GBAR_W + DBN1(age,ai,zi,lambdai,ei,xi) * &
	            & (( agrid(ai) + ( R_z(zi)*agrid(ai) + Pr_mat(ai,zi,xi) ) *(1.0_DP-tauK)  ) - YGRID(ai,zi,xi) )

      	GBAR_C = GBAR_C +  DBN1(age,ai,zi,lambdai,ei,xi) * tauC * cons(age,ai,zi,lambdai,ei,xi)

      	GBAR_BQ = GBAR_BQ +  DBN1(age,ai,zi,lambdai,ei,xi) * tau_bq * aprime(age,ai,zi,lambdai,ei,xi) * (1.0_dp-survP(age))

	    Tot_Lab_Inc = Tot_Lab_Inc + DBN1(age,ai,zi,lambdai,ei,xi) * yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei,xi)

	    Tot_Cap_Inc = Tot_Cap_Inc + DBN1(age,ai,zi,lambdai,ei,xi) * ( R_z(zi)*agrid(ai) + Pr_mat(ai,zi,xi) )
	    test = test + DBN1(age,ai,zi,lambdai,ei,xi) * ( R_z(zi)*agrid(ai) )
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	! Compute social security expenditure
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

	! Set government expenditures as total expenditures minus social security payments
	GBAR = GBAR -  SSC_Payments

	IF (print_flag) THEN
	! Print Results
	print*,' '; print*,'-----------------------------------------------------------------------------'
	print*, "Government Budget - Revenues and taxes"
	print*,' '
	print*,'	GBAR=',GBAR,'SSC_Payments=', SSC_Payments, 'GBAR_L=',GBAR_L,'Av. Labor Tax=', GBAR_L/Ebar 
	print*,'	GBAR_K=', GBAR_K, "GBAR_W=", GBAR_W, 'GBAR_C=', GBAR_C, 'GBAR_BQ=', GBAR_BQ
	print '(A,F7.4,X,X,A,F7.4,X,X,A,F7.4,X,X,A,F7.4,X,X,A,F7.4,X,X,A,F7.4)',&
			&'	Tau_K=', tauK, 'Tau_W_bt=', tauW_bt, 'Tau_W_at=', tauW_at,&
			& 'Tau_C=', tauC, 'Tau_BQ=', tau_bq, "Threshold", Y_a_threshold
	print*, ' '
	print '(A,F7.3)', ' 	Tax Revenue/GDP=      ', 100.0_dp*(GBAR_K+GBAR_W+GBAR_L+GBAR_C+GBAR_BQ)/YBAR
	print '(A,F7.3)', ' 	Capital_Tax/Total_Tax=', 100.0_dp*(GBAR_K+GBAR_W)/(GBAR_K+GBAR_W+GBAR_L+GBAR_C+GBAR_BQ)
	print '(A,F7.3)', ' 	Labor_Tax/Total_Tax  =', 100.0_dp*GBAR_L/(GBAR_K+GBAR_W+GBAR_L+GBAR_C+GBAR_BQ)
	print '(A,F7.3)', ' 	Estate_Tax/Total_Tax =', 100.0_dp*GBAR_BQ/(GBAR_K+GBAR_W+GBAR_L+GBAR_C+GBAR_BQ)
	print '(A,F7.3)', ' 	Capital_Tax/GDP      =', 100.0_dp*(GBAR_K+GBAR_W)/YBAR
	print '(A,F7.3)', ' 	Average Capital Tax  =', 100.0_dp*(GBAR_K+GBAR_W)/Tot_Cap_Inc
	print '(A,F7.3)', ' 	Labor_Tax/GDP        =', 100.0_dp*GBAR_L/YBAR
	print '(A,F7.3)', ' 	Average Labor Tax    =', 100.0_dp*GBAR_L/Tot_Lab_Inc
	print '(A,F7.3)', ' 	Estate_Tax/GDP       =', 100.0_dp*GBAR_BQ/YBAR
	print '(A,F7.3,X,X,A,F7.3)', ' 		Total Labor Income  =', Tot_Lab_Inc , 'EBAR=', EBAR
	print '(A,F7.3,X,X,A,F7.3)', ' 		Total Capital Income=', Tot_Cap_Inc , 'EBAR=', alpha*YBAR - DepRate*MeanWealth
	print*,'KBAR=',MeanWealth,'Integral=',sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid ),'K_C+K_P=',K_C+K_P
	print*,'K_C=',K_C,'Integral=',sum( sum(sum(sum(sum(sum(DBN1(:,:,z_C:,:,:,:),6),5),4),3),1)*agrid ), &
		   'K_P=',K_P,'Integral=',sum( sum(sum(sum(DBN1,5),4),1)*K_matrix(R,P) ) 
	print*,'R_Z=',100.0_dp*R_z
	print*,'Profit Income',sum( sum(sum(sum(DBN1,5),4),1)*Pr_mat ),alpha*YBAR_P-(R+DepRate)*K_P 
	print*,'Interest Income',test,R_C*K_C + R*K_P, alpha*YBAR_C-DepRate*K_C + R*K_P
	print*,'R_C=',R_C,alpha*A_C*K_C**(alpha-1)*L_C**(1-alpha)-DepRate
	print*,'alpha=',alpha,alpha_C
	print*,'-----------------------------------------------------------------------------'
	print*, ' '

	if (solving_bench.eq.1) then 
	OPEN(UNIT=11, FILE=trim(Result_Folder)//'Govnt_Budget_Bench.txt', STATUS='replace')
	else
	OPEN(UNIT=11, FILE=trim(Result_Folder)//'Govnt_Budget_Exp.txt', STATUS='replace')
	endif 
	WRITE(UNIT=11, FMT=*) ' '
	WRITE(UNIT=11, FMT=*) "Government Budget - Revenues and taxes"
	WRITE(UNIT=11, FMT=*) 'GBAR=',GBAR,'SSC_Payments=', SSC_Payments, 'GBAR_L=',GBAR_L,'Av.Labor_Tax=', GBAR_L/Ebar 
	WRITE(UNIT=11, FMT=*) 'GBAR_K=', GBAR_K, "GBAR_W=", GBAR_W, 'GBAR_C=', GBAR_C, 'GBAR_BQ=', GBAR_BQ
	WRITE(UNIT=11, FMT=*) 'Tau_K=', tauK, 'Tau_W=', tauW_at, 'Tau_C=', tauC, 'Tau_BQ=', tau_bq, "Threshold", Y_a_threshold
	WRITE(UNIT=11, FMT=*) 'Tax Revenue/GDP=      ', 100.0_dp*(GBAR_K+GBAR_W+GBAR_L+GBAR_C+GBAR_BQ)/YBAR
	WRITE(UNIT=11, FMT=*) 'Capital_Tax/Total_Tax=', 100.0_dp*(GBAR_K+GBAR_W)/(GBAR_K+GBAR_W+GBAR_L+GBAR_C+GBAR_BQ)
	WRITE(UNIT=11, FMT=*) 'Labor_Tax/Total_Tax  =', 100.0_dp*GBAR_L/(GBAR_K+GBAR_W+GBAR_L+GBAR_C+GBAR_BQ)
	WRITE(UNIT=11, FMT=*) 'Estate_Tax/Total_Tax =', 100.0_dp*GBAR_BQ/(GBAR_K+GBAR_W+GBAR_L+GBAR_C+GBAR_BQ)
	WRITE(UNIT=11, FMT=*) 'Capital_Tax/GDP      =', 100.0_dp*(GBAR_K+GBAR_W)/YBAR
	WRITE(UNIT=11, FMT=*) 'Average Capital Tax  =', 100.0_dp*(GBAR_K+GBAR_W)/Tot_Cap_Inc
	WRITE(UNIT=11, FMT=*) 'Labor_Tax/GDP        =', 100.0_dp*GBAR_L/YBAR
	WRITE(UNIT=11, FMT=*) 'Average Labor Tax    =', 100.0_dp*GBAR_L/Tot_Lab_Inc
	WRITE(UNIT=11, FMT=*) 'Estate_Tax/GDP       =', 100.0_dp*GBAR_BQ/YBAR
	WRITE(UNIT=11, FMT=*) 'Total_Labor_Income', Tot_Lab_Inc , 'EBAR', EBAR, 'Total_Capital_Income', Tot_Cap_Inc
	Close(UNIT=11)
	ENDIF  
END  SUBROUTINE GOVNT_BUDGET

!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE FIND_DBN_EQ()
	use omp_lib
	IMPLICIT NONE
	INTEGER  :: tklo, tkhi, age1, age2, z1, z2, a1, a2, lambda1, lambda2, e1, e2, DBN_iter, simutime, iter_indx, x1, x2
	REAL   	 :: DBN_dist, DBN_criteria
	REAL(DP) :: BBAR, Wealth, brent_value
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable ::  PrAprimelo, PrAprimehi, PrBqlo, PrBqhi, DBN2
	INTEGER , DIMENSION(:,:,:,:,:,:), allocatable ::  Aplo, Aphi, Bqlo, Bqhi
	! Timing
	real(kind=8)    :: t1, t2, elapsed_time
    integer(kind=8) :: tclock1, tclock2, clock_rate

	allocate( DBN2(         MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrAprimelo(   MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrAprimehi(   MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Aplo(         MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Aphi(         MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrBqlo(   	MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrBqhi(   	MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Bqlo(         MaxAge,na,nz,nlambda,ne,nx) )
	allocate( BQhi(         MaxAge,na,nz,nlambda,ne,nx) )

	!$ call omp_set_num_threads(nz)
	DBN_criteria = 1.00E-07_DP


	! ! Current aggregate values given QBAR and Wage and R
	! if (A_C.gt.0.0_dp) then 
	! 	L_P    = ( (1.0_dp-alpha)*Aprod/Wage )**(1.0_dp/alpha) * QBAR 
	! 	P      = alpha * QBAR**(alpha-mu) * L_P**(1.0_DP-alpha)
 !    	K_P    = sum( (sum(sum(sum(DBN1,5),4),1)) *(K_Matrix(R,P))) ! Note: DBN_azx  = sum(sum(sum(DBN1,5),4),1)
    	
 !    	Wealth = sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid )

 !    	K_C    = sum( sum(sum(sum(sum(sum(DBN1(:,:,z_C:,:,:,:),6),5),4),3),1) * agrid ) 
 !    	L_C    = ( (1.0_dp-alpha_C)*A_C/Wage )**(1.0_dp/alpha_C) * K_C
 !    	R_C    = alpha_C * A_C * (L_C/K_C)**(1.0_dp-alpha_C) - DepRate

	! 	YBAR_P = AProd * QBAR**alpha   * L_P**(1.0_DP-alpha  ) 
	! 	YBAR_C = A_C   * K_C **alpha_C * L_C**(1.0_DP-alpha_C) 
	! 	YBAR   = YBAR_P + YBAR_C

	! 	! print*,'initial values:',Wealth,K_P,K_C,L_C,L_P,YBAR,YBAR_P,YBAR_C,P,R,R_C
	! endif ! If A_C=0 then all aggregates must be provided


	! Solve the model at current aggregate values
		! Find the threshold for wealth taxes (a_bar)
			!call Find_TauW_Threshold(DBN1,Y_a_threshold)
		! Adjust grid to include breaking points
			CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
		! Compute labor units 
			CALL ComputeLaborUnits(Ebar, wage) 
		! Compute Capital demand and Profits by (a,z)
			K_mat  = K_Matrix(R,P)
			Pr_mat = Profit_Matrix(R,P)
		! Form YGRID for the capital income economy given interest rate "P"
			CALL FORM_Y_MB_GRID(YGRID,MBGRID,YGRID_t,MBGRID_t)
		! Solve for policy and value functions 
			CALL EGM_RETIREMENT_WORKING_PERIOD 

	! Discretize policy function for assets (a')
		! For each age and state vector bracket optimal a' between two grid points
		! When at that age and state the optimal decision is approximated by selecting one the grid points
		! The grid points are selected with probability proportional to their distance to the optimal a'
	DO age=1,MaxAge
	!$omp parallel do private(lambdai,ei,ai,xi,tklo,tkhi)
	DO zi=1,nz
	!$omp parallel do private(lambdai,ei,ai,tklo,tkhi)
	DO xi=1,nx
	DO ai=1,na
	DO lambdai=1,nlambda
	DO ei=1, ne
        if ( Aprime(age,ai,zi,lambdai,ei,xi) .ge. amax) then
            tklo =na-1
        elseif (Aprime(age,ai,zi,lambdai, ei,xi) .lt. amin) then
            tklo = 1
        else
            tklo = ((Aprime(age,ai,zi,lambdai, ei,xi) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
        endif
        tkhi = tklo + 1        
        Aplo(age,ai,zi,lambdai,ei,xi)  		= tklo
        Aphi(age,ai,zi,lambdai,ei,xi)  		= tkhi        
        PrAprimelo(age,ai,zi,lambdai,ei,xi) = ( agrid(tkhi) - Aprime(age,ai,zi,lambdai,ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
        PrAprimehi(age,ai,zi,lambdai,ei,xi) = ( Aprime(age,ai,zi,lambdai,ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )

        if ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) .ge. amax) then
            tklo =na-1
        elseif ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) .lt. amin) then
            tklo = 1
        else
            tklo = (((1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi)-amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
        endif
        tkhi = tklo + 1        
        Bqlo(age,ai,zi,lambdai,ei,xi)  	= tklo
        Bqhi(age,ai,zi,lambdai,ei,xi)  	= tkhi        
        PrBqlo(age,ai,zi,lambdai,ei,xi) = &
        	& ( agrid(tkhi) - (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
        PrBqhi(age,ai,zi,lambdai,ei,xi) = &
        	& ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	! Probablities are adjusted to lie in [0,1]
		PrAprimelo = min(PrAprimelo, 1.0_DP); PrAprimelo = max(PrAprimelo, 0.0_DP)
		PrAprimehi = min(PrAprimehi, 1.0_DP); PrAprimehi = max(PrAprimehi, 0.0_DP)

		PrBqlo = min(PrBqlo, 1.0_DP); PrBqlo = max(PrBqlo, 0.0_DP)
		PrBqhi = min(PrBqhi, 1.0_DP); PrBqhi = max(PrBqhi, 0.0_DP)

	! Compute distribution of assets by age and state
		! Distribution is obtained by iterating over an initial distribution using policy functions
	DBN_dist=1.0_DP
	simutime = 1
	iter_indx = 1
	print*,' '; print*,'-----------------------------------------------------------------------------'
	print*, 'Computing Equilibrium Distribution'; print*,' '
	DO WHILE ( ( DBN_dist .ge. DBN_criteria ) .and. ( simutime .le. 700 ) )
		! print*, 'DBN_dist=', DBN_dist

	    DBN2=0.0_DP

		! Everyone in MaxAge dies. Those who die, switch to z2, lambda2 and start at ne/2+1 and x=1
	    age1=MaxAge
	    !$omp parallel do reduction(+:DBN2) private(x1,a1,lambda1,e1,z2,lambda2)
	    DO z1=1,nz
	    DO x1=1,nx
	    DO a1=1,na
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN2(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqhi(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO    
	    ENDDO
	    ENDDO  
	    !$omp barrier  

		! retirees "e" stays the same for benefit retirement calculation purposes
		!$omp parallel do reduction(+:DBN2) private(x1,age1,a1,lambda1,e1,z2,lambda2,x2)
		DO z1=1,nz
	    DO x1=1,nx
	    DO age1=RetAge-1, MaxAge-1
	    DO a1=1,na
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        ! Those who die, switch to z2, lambda2 and start at ne/2+1
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN2(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqhi(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, and also e1 since they are retired
	        !e2=e1
	        DO x2=1,nx
		        DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) =  &
		        	& DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
		            & * survP(age1) * PrAprimelo(age1, a1, z1, lambda1, e1,x1) * pr_x(x1,x2,z1,age1)     
		        DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) =  &
		          	& DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
		            & * survP(age1) * PrAprimehi(age1, a1, z1, lambda1, e1,x1) * pr_x(x1,x2,z1,age1)  
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    !$omp barrier
	    
	    ! Working age agents
	    !$omp parallel do reduction(+:DBN2) private(x1,age1,a1,lambda1,e1,z2,lambda2,x2,e2)
	    DO z1=1,nz
	    DO x1=1,nx
	    DO age1=1,RetAge-2
	    DO a1=1,na
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        ! Those who die, switch to z2, lambda2 and start at ne/2+1
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN2(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqhi(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, but switch to e2 and x2
	        DO x2=1,nx
	        DO e2=1,ne
				DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) =  &
	          		& DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
	                & * survP(age1) * pr_e(e1,e2) * pr_x(x1,x2,z1,age1) * PrAprimelo(age1, a1, z1, lambda1, e1,x1)     
	            DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) =  &
	          		& DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
	                & * survP(age1) * pr_e(e1,e2) * pr_x(x1,x2,z1,age1) * PrAprimehi(age1, a1, z1, lambda1, e1,x1) 
	        ENDDO
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    !$omp barrier
	    
	    DBN2 = 0.7_dp*DBN1 + 0.3_dp*DBN2
	    DBN_dist = maxval(abs(DBN2-DBN1))
	    ! print*, DBN_dist
	    DBN1 = DBN2

	    ! Instead of adjusting policy functions to be consistent with current DBN at each iteration
	    ! The code iterates the distribution and only updates policy functions every "update_period"

	    ! Update of policy function with current aggregates
	    IF (iter_indx .ge. update_period) THEN

	    	! Compute aggregates with current distribution: New QBAR and Agg. Labor Supply (NBAR)
	        QBAR =0.0
	        NBAR =0.0
	        !$omp parallel do reduction(+:QBAR,NBAR) private(x1,age1,a1,lambda1,e1)
	        DO z1=1,nz
	        DO x1=1,nx
	        DO age1=1,MaxAge
	        DO a1=1,na
	        DO lambda1=1,nlambda
	        DO e1=1, ne
	             QBAR= QBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * ( xz_grid(x1,z1) * K_mat(a1,z1,x1) )**mu
	             NBAR= NBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * eff_un(age1, lambda1, e1) * Hours(age1, a1, z1, lambda1,e1,x1)
	        ENDDO
	        ENDDO
	        ENDDO
	        ENDDO    
	        ENDDO    
	        ENDDO 
	        !$omp barrier
	        QBAR = ( QBAR)**(1.0_DP/mu)


	        if (A_C.gt.0.0_dp) then 

	        	! Update Interest Rate (using old P)
	        	brent_value = brent(-0.1_DP,0.01_DP,10.0_DP,Agg_Debt_C, brent_tol,R)

		    	! Total wealth 
		    	Wealth = sum( sum(sum(sum(sum(sum(DBN1                ,6),5),4),3),1)*agrid )

	        	! Capital
        		K_C    = sum( sum(sum(sum(sum(sum(DBN1(:,:,z_C:,:,:,:),6),5),4),3),1)*agrid ) 
        		K_P    = Wealth - K_C - Debt_Absorption*Debt_SS 

	        	! Update Wage 
	        	if (alpha_C.eq.alpha) then 
	        		Wage = ((((1.0_dp-alpha)*Aprod)**(1.0_dp/alpha)*QBAR + ((1.0_dp-alpha)*A_C)**(1.0_dp/alpha)*K_C)/NBAR)**(alpha)
	        	else 
	        		print*, ' Solving labor market numerically'
	        		brent_value = brent(0.001_DP,Wage,10.0_DP,Labor_Market_Clearing, brent_tol,Wage)
	        		print*, ' New Wage=',Wage,'Error=',brent_value
	        	endif 

	        	! Update Other Corporate values 
		    	L_C    = ( (1.0_dp-alpha_C)*A_C/Wage )**(1.0_dp/alpha_C) * K_C
		    	R_C    = alpha_C * A_C * (L_C/K_C)**(1.0_dp-alpha_C) - DepRate
		    	

	    		! Check that R_C>R (if not solve the equilibrium with a single interest rate)
		    	if (R_C.lt.R) then 
		    		print*,' '; 
		    		print*,'Interest Rates Not in Order'
		    		print*,'R',100.0_dp*R,'R_C',100.0_dp*R_C,'brentvalue',brent_value

		        	! Private demand for capital
		        	K_P    = sum( (sum(sum(sum(DBN1,5),4),1)) *(K_matrix(R,P))) ! Note: DBN_azx  = sum(sum(sum(DBN1,5),4),1)
	        	
		        	if (Wealth.gt.(K_P+Debt_Absorption*Debt_SS)) then 
			        	! Corporate demand for capital
			        	K_C    = Wealth - K_P - Debt_Absorption*Debt_SS

			        	! Update Wage 
			        	if (alpha_C.eq.alpha) then 
			        		Wage = ((((1.0_dp-alpha)*Aprod)**(1.0_dp/alpha)*QBAR + ((1.0_dp-alpha)*A_C)**(1.0_dp/alpha)*K_C)/NBAR)**(alpha)
			        	else 
			        		print*, ' Solving labor market numerically'
			        		brent_value = brent(0.001_DP,Wage,10.0_DP,Labor_Market_Clearing, brent_tol,Wage)
			        		print*, ' New Wage=',Wage,'Error=',brent_value
			        	endif 

		        	else ! Capital Market did not clear

		        		print*, ' ';print*, ' 	Warning! Capital Market Did not Clear!';print*, ' ';
			        	! Modify prices and quantities 
			        	QBAR = 0.50_dp*QBAR 
			        	Wage = Wage 

		        	endif 

		        	! Update Interest Rates 
		        	R    = alpha_C * A_C * ( Wage/((1.0_dp-alpha_C)*A_C) )**(-(1.0_dp-alpha_C)/alpha_C) - DepRate
					R_C  = R 

		    	endif 


		    	! Update Output and Aggregates 
				L_P    = ( (1.0_dp-alpha)*Aprod/Wage )**(1.0_dp/alpha  ) * QBAR 
				L_C    = ( (1.0_dp-alpha_C)*A_C/Wage )**(1.0_dp/alpha_C) * K_C
				P      = alpha * QBAR**(alpha-mu) * L_P**(1.0_DP-alpha)
				Ebar   = wage  * NBAR  * sum(pop)/sum(pop(1:RetAge-1))
				YBAR_P = AProd * QBAR**alpha   * L_P**(1.0_DP-alpha  ) 
				YBAR_C = A_C   * K_C **alpha_C * L_C**(1.0_DP-alpha_C) 
				YBAR   = YBAR_P + YBAR_C

	        	
	        	! Update interest rate vector 
        		R_z(1:z_C-1) = R 
        		R_z(z_C:)    = R_C

				! print '(A,F9.3,F9.3,F9.3,X,X,A,F9.3,F9.3,F9.3)',&
				! 		& ' 				Corp. Sector Levels:', YBAR_C, K_C, L_C , &
				! 		& ' Ratios ', 100.0_dp*YBAR_C/YBAR, 100.0_dp*K_C/Wealth, 100.0_dp*L_C/NBAR

	        else 
	        	! Set Corporate values to zero
        		K_C    = 0.0_dp; L_C    = 0.0_dp; YBAR_C = 0.0_dp;  R_C = 0.0_dp 
        		K_P    = 0.0_dp; L_P    = 0.0_dp; YBAR_P = 0.0_dp
        		Wealth = sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid )

	        	! Solve for aggregates and clear capital market with R
		        P    = alpha* QBAR **(alpha-mu) * NBAR **(1.0_DP-alpha)
		        YBAR = QBAR ** alpha * NBAR **(1.0_DP-alpha)
		        wage = (1.0_DP-alpha)*QBAR **alpha * NBAR  **(-alpha)
		        Ebar = wage  * NBAR  * sum(pop)/sum(pop(1:RetAge-1))

		    	! Solve for new R 
		    	! R = zbrent(Agg_Debt,0.1_dp,1.00_dp,brent_tol) 
		    	if (sum(theta)/nz .gt. 1.0_DP) then
		    		P = min(P,1.0_dp)
		            brent_value = brent(-0.1_DP,0.01_DP,10.0_DP,Agg_Debt, brent_tol,R)
		        else
		            R = 0.0_DP
		        endif
		        	R_z = R

	        endif 

	    	!!
	    	print 12345, &
	    		& ' DBN_diff=', DBN_dist,'A=',sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid ),&
	    		& 'W=',wage,'R=',R,'R_C=',R_C,'P=',P,'Q=',QBAR, &
	    		& 'K_C/A=',100.0_dp*K_C/Wealth,'L_C/N=',100.0_dp*L_C/NBAR,'K_C=',K_C,'L_C=',L_C
    		print 12347,' 		K_C=',K_C,'K_P=',K_P,'A=',sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid ),'K_C+K_P=',K_P+K_C
			print 12347,' 		L_C=',L_C,'L_P=',L_P,'N=',NBAR,'L_C+L_P=',L_C+L_P
    		12345 format &
    		&(A,E12.5,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X)
    		12347 format &
    		&(A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X)
	    	!!

	    	! Solve the model at current aggregate values
				! Find the threshold for wealth taxes (a_bar)
					!call Find_TauW_Threshold(DBN1,Y_a_threshold)
				! Adjust grid to include breaking points
					CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
				! Compute labor units 
					CALL ComputeLaborUnits(Ebar, wage) 
				! Compute Capital demand and Profits by (a,z)
					K_mat  = K_Matrix(R,P)
					Pr_mat = Profit_Matrix(R,P)
				! Form YGRID for the capital income economy given interest rate "P"
					CALL FORM_Y_MB_GRID(YGRID,MBGRID,YGRID_t,MBGRID_t)
				! Solve for policy and value functions 
					CALL EGM_RETIREMENT_WORKING_PERIOD 
	        
				! Discretize policy function for assets (a')
					! For each age and state vector bracket optimal a' between two grid points
					! When at that age and state the optimal decision is approximated by selecting one the grid points
					! The grid points are selected with probability proportional to their distance to the optimal a'
				DO age=1,MaxAge
				!$omp parallel do private(lambdai,ei,ai,xi,tklo,tkhi)
				DO zi=1,nz
				!$omp parallel do private(lambdai,ei,ai,tklo,tkhi)
				DO xi=1,nx
				DO ai=1,na
				DO lambdai=1,nlambda
				DO ei=1, ne
			        if ( Aprime(age,ai,zi,lambdai,ei,xi) .ge. amax) then
			            tklo =na-1
			        elseif (Aprime(age,ai,zi,lambdai, ei,xi) .lt. amin) then
			            tklo = 1
			        else
			            tklo = ((Aprime(age,ai,zi,lambdai, ei,xi) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
			        endif
			        tkhi = tklo + 1        
			        Aplo(age,ai,zi,lambdai,ei,xi)  		= tklo
			        Aphi(age,ai,zi,lambdai,ei,xi)  		= tkhi        
			        PrAprimelo(age,ai,zi,lambdai,ei,xi) = &
			        	& ( agrid(tkhi) - Aprime(age,ai,zi,lambdai,ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
			        PrAprimehi(age,ai,zi,lambdai,ei,xi) = &
			        	& ( Aprime(age,ai,zi,lambdai,ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )

			        if ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) .ge. amax) then
			            tklo =na-1
			        elseif ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) .lt. amin) then
			            tklo = 1
			        else
			            tklo = (((1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi)-amin)/(amax-amin))&
			            	&**(1.0_DP/a_theta)*(na-1)+1          
			        endif
			        tkhi = tklo + 1        
			        Bqlo(age,ai,zi,lambdai,ei,xi)  	= tklo
			        Bqhi(age,ai,zi,lambdai,ei,xi)  	= tkhi        
			        PrBqlo(age,ai,zi,lambdai,ei,xi) = ( agrid(tkhi)-(1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) )&
			        									& / ( agrid(tkhi) -agrid(tklo) )
			        PrBqhi(age,ai,zi,lambdai,ei,xi) = ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi)-agrid(tklo) )&
			        									&  / ( agrid(tkhi) -agrid(tklo) )
				ENDDO
				ENDDO
				ENDDO
				ENDDO
				ENDDO
				ENDDO

				! Probablities are adjusted to lie in [0,1]
					PrAprimelo = min(PrAprimelo, 1.0_DP); PrAprimelo = max(PrAprimelo, 0.0_DP)
					PrAprimehi = min(PrAprimehi, 1.0_DP); PrAprimehi = max(PrAprimehi, 0.0_DP)

					PrBqlo = min(PrBqlo, 1.0_DP); PrBqlo = max(PrBqlo, 0.0_DP)
					PrBqhi = min(PrBqhi, 1.0_DP); PrBqhi = max(PrBqhi, 0.0_DP) 

		    ! Reset counter for next update of policy functions
	        iter_indx=0
	    ENDIF
	    
	    iter_indx = iter_indx + 1
	    simutime  = simutime +1 
	 
	ENDDO ! WHILE
	print*,' '
	print*,' 	Stationary Equilibrium Found: '
	print 12346, &
		& ' 	DBN_diff=', DBN_dist,'A=',sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid ),&
		& 'W=',wage,'R=',R,'R_C=',R_C,'P=',P,'Q=',QBAR, &
		& 'K_C/A=',100.0_dp*K_C/Wealth,'L_C/N=',100.0_dp*L_C/NBAR,'Y_C/Y=',100.0_dp*YBAR_C/YBAR,'Iter=',simutime
	print 12347,' 		K_C=',K_C,'K_P=',K_P,'A=',sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid ),'K_C+K_P=',K_P+K_C
	print 12347,' 		L_C=',L_C,'L_P=',L_P,'N=',NBAR,'L_C+L_P=',L_C+L_P
	12346 format &
	& (A,E12.5,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,I5)
	print*,' '; print*,'-----------------------------------------------------------------------------'
	print*,' '


	! Write
		!OPEN(UNIT=3, FILE='agrid', STATUS='replace')
		!WRITE(unit=3, FMT=*) agrid
		!CLOSE (unit=3)
		!
		!OPEN(UNIT=3, FILE='aprime', STATUS='replace')
		!WRITE(unit=3, FMT=*) Aprime(1, :, nz/2+1, nlambda/2+1, ne/2+1)
		!CLOSE (unit=3)
		!
		!OPEN(UNIT=3, FILE='aplo', STATUS='replace')
		!WRITE(unit=3, FMT=*) Aplo(1, :,nz/2+1, nlambda/2+1, ne/2+1)
		!CLOSE (unit=3)
		!
		!OPEN(UNIT=3, FILE='aphi', STATUS='replace')
		!WRITE(unit=3, FMT=*) Aphi(1, :, nz/2+1, nlambda/2+1, ne/2+1)
		!CLOSE (unit=3)
		!
		!OPEN(UNIT=3, FILE='Pr_aplo', STATUS='replace')
		!WRITE(unit=3, FMT=*) PrAprimelo(1, :, nz/2+1, nlambda/2+1, ne/2+1)
		!CLOSE (unit=3)
		!
		!OPEN(UNIT=3, FILE='pr_aphi', STATUS='replace')
		!WRITE(unit=3, FMT=*) PrAprimehi(1, :, nz/2+1, nlambda/2+1, ne/2+1)
		!CLOSE (unit=3)

END SUBROUTINE FIND_DBN_EQ


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE FIND_DBN_EQ_PF()
	use omp_lib
	IMPLICIT NONE
	INTEGER:: tklo, tkhi, age1, age2, z1, z2, a1, a2, lambda1, lambda2, e1, e2, DBN_iter, simutime, iter_indx, x1, x2
	REAL   :: DBN_dist, DBN_criteria
	real(dp)   ::BBAR, MeanWealth, brent_value
	REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne, nx) :: PrAprimelo, PrAprimehi, PrBqlo, PrBqhi, DBN2
	INTEGER,  DIMENSION(MaxAge, na, nz, nlambda, ne, nx) :: Aplo, Aphi, Bqlo, Bqhi
	real(dp), dimension(na,nz,nx) :: YGRID_old

	!$ call omp_set_num_threads(nz)
	DBN_criteria = 1.0E-07_DP
	DBN1 = DBN_bench

	! Solve the model at current aggregate values
		! Find the threshold for wealth taxes (a_bar)
			!call Find_TauW_Threshold(DBN1,Y_a_threshold)
		! Adjust grid to include breaking points
			CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
		! Compute labor units 
			CALL ComputeLaborUnits(Ebar, wage) 
		! Compute Capital demand and Profits by (a,z)
			K_mat  = K_Matrix(R,P)
			Pr_mat = Profit_Matrix(R,P)
		! Form YGRID for the capital income economy given interest rate "P"
			YGRID_old = YGRID 
			CALL FORM_Y_MB_GRID(YGRID,MBGRID,YGRID_t,MBGRID_t)
			! print*,' '; print*, 'YGRID';
			! print*, YGRID_old(5,4,1) , YGRID(5,4,1)
			! print*, YGRID_old(55,4,1) , YGRID(55,4,1)
			! print*, YGRID_old(105,4,1) , YGRID(105,4,1)
			! print*, YGRID_old(155,4,1) , YGRID(155,4,1)
		

	! Solve for policy and value functions 
			! Instead of EGM I only update the savings policy function
	! Discretize policy function for assets (a')
		! For each age and state vector bracket optimal a' between two grid points
		! When at that age and state the optimal decision is approximated by selecting one the grid points
		! The grid points are selected with probability proportional to their distance to the optimal a'
	DO age=1,MaxAge
	! !$omp parallel do private(lambdai,ei,ai,xi,tklo,tkhi)
	DO zi=1,nz
	DO xi=1,nx
	DO ai=1,na
	DO lambdai=1,nlambda
	DO ei=1, ne
		! Get Aprime
		if (age.lt.RetAge) then 
		Aprime(age, ai, zi, lambdai,ei,xi) = max( amin , min( amax , &
					& YGRID(ai,zi,xi)  + Y_h(Hours(age, ai, zi, lambdai,ei,xi),age,lambdai,ei,wage)  & 
		            & - (1.0_dp+tauC)*Cons(age, ai, zi, lambdai,ei,xi) ) )
        else 
        Aprime(age, ai, zi, lambdai,ei,xi) = max( amin , min( amax , &
					& YGRID(ai,zi,xi)  + RetY_lambda_e(lambdai,ei) - (1.0_dp+tauC)*Cons(age, ai, zi, lambdai,ei,xi) ))
        endif 
        ! Discretize Aprime
        if ( Aprime(age,ai,zi,lambdai,ei,xi) .ge. amax) then
            tklo =na-1
        elseif (Aprime(age,ai,zi,lambdai, ei,xi) .lt. amin) then
            tklo = 1
        else
            tklo = ((Aprime(age,ai,zi,lambdai, ei,xi) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
        endif
        tkhi = tklo + 1        
        Aplo(age,ai,zi,lambdai,ei,xi)  		= tklo
        Aphi(age,ai,zi,lambdai,ei,xi)  		= tkhi        
        PrAprimelo(age,ai,zi,lambdai,ei,xi) = ( agrid(tkhi) - Aprime(age,ai,zi,lambdai,ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
        PrAprimehi(age,ai,zi,lambdai,ei,xi) = ( Aprime(age,ai,zi,lambdai,ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )

        if ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) .ge. amax) then
            tklo =na-1
        elseif ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) .lt. amin) then
            tklo = 1
        else
            tklo = (((1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi)-amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
        endif
        tkhi = tklo + 1        
        Bqlo(age,ai,zi,lambdai,ei,xi)  	= tklo
        Bqhi(age,ai,zi,lambdai,ei,xi)  	= tkhi        
        PrBqlo(age,ai,zi,lambdai,ei,xi) = &
        	& ( agrid(tkhi) - (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
        PrBqhi(age,ai,zi,lambdai,ei,xi) = &
        	& ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	! Probablities are adjusted to lie in [0,1]
		PrAprimelo = min(PrAprimelo, 1.0_DP); PrAprimelo = max(PrAprimelo, 0.0_DP)
		PrAprimehi = min(PrAprimehi, 1.0_DP); PrAprimehi = max(PrAprimehi, 0.0_DP)

		PrBqlo = min(PrBqlo, 1.0_DP); PrBqlo = max(PrBqlo, 0.0_DP)
		PrBqhi = min(PrBqhi, 1.0_DP); PrBqhi = max(PrBqhi, 0.0_DP)

	! Compute distribution of assets by age and state
		! Distribution is obtained by iterating over an initial distribution using policy functions
	DBN_dist=1.0_DP
	simutime = 1
	iter_indx = 1
	!print*, 'Computing Equilibrium Distribution'
	DO WHILE ( ( DBN_dist .ge. DBN_criteria ) .and. ( simutime .le. 700 ) )
		! print*, 'DBN_dist=', DBN_dist

	    DBN2=0.0_DP

		! Everyone in MaxAge dies. Those who die, switch to z2, lambda2 and start at ne/2+1 and x=1
	    age1=MaxAge
	    DO x1=1,nx
	    DO z1=1,nz
	    DO a1=1,na
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        ! Those who die, switch to z2, lambda2 and start at ne/2+1
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN2(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqhi(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO    
	    ENDDO
	    ENDDO    

		! retirees "e" stays the same for benefit retirement calculation purposes
	    DO x1=1,nx
	    DO age1=RetAge-1, MaxAge-1
	    DO a1=1,na
	    DO z1=1,nz
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        ! Those who die, switch to z2, lambda2 and start at ne/2+1
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN2(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqhi(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, and also e1 since they are retired
	        !e2=e1
	        DO x2=1,nx
		        DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) =  &
		        	& DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
		            & * survP(age1) * PrAprimelo(age1, a1, z1, lambda1, e1,x1) * pr_x(x1,x2,z1,age1)     
		        DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) =  &
		          	& DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
		            & * survP(age1) * PrAprimehi(age1, a1, z1, lambda1, e1,x1) * pr_x(x1,x2,z1,age1)  
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    
	    ! Working age agents
	    DO x1=1,nx
	    DO age1=1,RetAge-2
	    DO a1=1,na
	    DO z1=1,nz
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        ! Those who die, switch to z2, lambda2 and start at ne/2+1
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN2(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqhi(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, but switch to e2 and x2
	        DO x2=1,nx
	        DO e2=1,ne
				DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) =  &
	          		& DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
	                & * survP(age1) * pr_e(e1,e2) * pr_x(x1,x2,z1,age1) * PrAprimelo(age1, a1, z1, lambda1, e1,x1)     
	            DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) =  &
	          		& DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
	                & * survP(age1) * pr_e(e1,e2) * pr_x(x1,x2,z1,age1) * PrAprimehi(age1, a1, z1, lambda1, e1,x1) 
	        ENDDO
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    
	    DBN2 = 0.8*DBN1 + 0.2*DBN2
	    DBN_dist = maxval(abs(DBN2-DBN1))
	    ! print*, DBN_dist
	    DBN1 = DBN2

	    ! Instead of adjusting policy functions to be consistent with current DBN at each iteration
	    ! The code iterates the distribution and only updates policy functions every "update_period"


	    	!!
	    	print*, 'DBN_diff=', DBN_dist, 'R=',R,'P=',P,'Assets',sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid )
	    	!!
	    simutime  = simutime +1 
	 
	ENDDO ! WHILE

    	! ! Compute aggregates with current distribution
	    !     QBAR =0.0
	    !     NBAR =0.0
	    !     DO x1=1,nx
	    !     DO age1=1,MaxAge
	    !     DO z1=1,nz
	    !     DO a1=1,na
	    !     DO lambda1=1,nlambda
	    !     DO e1=1, ne
	    !          QBAR= QBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * ( xz_grid(x1,z1) * K_mat(a1,z1,x1) )**mu
	    !          NBAR= NBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * eff_un(age1, lambda1, e1) * Hours(age1, a1, z1, lambda1,e1,x1)
	    !     ENDDO
	    !     ENDDO
	    !     ENDDO
	    !     ENDDO    
	    !     ENDDO    
	    !     ENDDO    
	    
	    !     QBAR = ( QBAR)**(1.0_DP/mu)                
	    !     YBAR = QBAR ** alpha * NBAR **(1.0_DP-alpha)
	    !     Ebar = wage  * NBAR  * sum(pop)/sum(pop(1:RetAge-1))

    	! Deallocate policy functions on adjusted grid (so that they can be allocated later)
			deallocate( YGRID_t  )
			deallocate( MBGRID_t ) 
			deallocate( Cons_t   )
			deallocate( Hours_t  )
			deallocate( Aprime_t )


END SUBROUTINE FIND_DBN_EQ_PF


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE FIND_DBN_EQ_PF_Interp(YGRID_bench)
	use omp_lib
	IMPLICIT NONE
	REAL(DP), DIMENSION(na,nz,nx), intent(in) :: YGRID_bench
	INTEGER:: tklo, tkhi, age1, age2, z1, z2, a1, a2, lambda1, lambda2, e1, e2, DBN_iter, simutime, iter_indx, x1, x2
	REAL   :: DBN_dist, DBN_criteria
	real(dp)   :: BBAR, MeanWealth, brent_value
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable ::  PrAprimelo, PrAprimehi, PrBqlo, PrBqhi, DBN2
	INTEGER , DIMENSION(:,:,:,:,:,:), allocatable ::  Aplo, Aphi, Bqlo, Bqhi

	allocate( DBN2(         MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrAprimelo(   MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrAprimehi(   MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Aplo(         MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Aphi(         MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrBqlo(   	MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrBqhi(   	MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Bqlo(         MaxAge,na,nz,nlambda,ne,nx) )
	allocate( BQhi(         MaxAge,na,nz,nlambda,ne,nx) )

	!$ call omp_set_num_threads(nz)
	DBN_criteria = 1.0E-07_DP
	DBN1 = DBN_bench

		! Form YGRID for the capital income economy given interest rate "P"
			CALL FORM_Y_MB_GRID(YGRID,MBGRID,YGRID_t,MBGRID_t)

		

	! Solve for policy and value functions 
			! Instead of EGM I update policy funtions by interpolating from the benchmark policy functions
	! Discretize policy function for assets (a')
		! For each age and state vector bracket optimal a' between two grid points
		! When at that age and state the optimal decision is approximated by selecting one the grid points
		! The grid points are selected with probability proportional to their distance to the optimal a'
	DO age=1,MaxAge
	!$omp parallel do private(lambdai,ei,ai,xi,tklo,tkhi)
	DO zi=1,nz
	DO xi=1,nx
	DO ai=1,na
	DO lambdai=1,nlambda
	DO ei=1, ne
		! Get Policy Functions
		Aprime(age,ai,zi,lambdai,ei,xi) = Linear_Int(YGRID_bench(:,zi,xi),Aprime_bench(age,:,zi,lambdai,ei,xi),na, YGRID(ai,zi,xi))
		Cons(age,ai,zi,lambdai,ei,xi)   = Linear_Int(YGRID_bench(:,zi,xi),Cons_bench(age,:,zi,lambdai,ei,xi)  ,na, YGRID(ai,zi,xi))
		Hours(age,ai,zi,lambdai,ei,xi)  = Linear_Int(YGRID_bench(:,zi,xi),Hours_bench(age,:,zi,lambdai,ei,xi) ,na, YGRID(ai,zi,xi)) 
		! ! Check
		! if (age.lt.RetAge) then 
		! print*, 'Residual',  YGRID(ai,zi,xi)  + Y_h(Hours(age, ai, zi, lambdai,ei,xi),age,lambdai,ei,wage)  & 
		! & - (1.0_dp+tauC)*Cons(age, ai, zi, lambdai,ei,xi) - Aprime(age, ai, zi, lambdai,ei,xi), &
		! & YGRID_bench(ai,zi,xi)  + Y_h(Hours_bench(age, ai, zi, lambdai,ei,xi),age,lambdai,ei,wage)  & 
		! & - (1.0_dp+tauC)*Cons_bench(age, ai, zi, lambdai,ei,xi) - Aprime_bench(age, ai, zi, lambdai,ei,xi),&
		! &'State', age,ai,zi,lambdai,ei,xi
  !       else 
  !       print*, 'Residual',  YGRID(ai,zi,xi)  + RetY_lambda_e(lambdai,ei)  & 
		! & - (1.0_dp+tauC)*Cons(age, ai, zi, lambdai,ei,xi) - Aprime(age, ai, zi, lambdai,ei,xi),&
		!  YGRID_bench(ai,zi,xi)  + RetY_lambda_e(lambdai,ei)  & 
		! & - (1.0_dp+tauC)*Cons_bench(age, ai, zi, lambdai,ei,xi) - Aprime_bench(age, ai, zi, lambdai,ei,xi),&
		! & 'State', age,ai,zi,lambdai,ei,xi
  !       endif 
        ! Discretize Aprime
        if ( Aprime(age,ai,zi,lambdai,ei,xi) .ge. amax) then
            tklo =na-1
        elseif (Aprime(age,ai,zi,lambdai, ei,xi) .lt. amin) then
            tklo = 1
        else
            tklo = ((Aprime(age,ai,zi,lambdai, ei,xi) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
        endif
        tkhi = tklo + 1        
        Aplo(age,ai,zi,lambdai,ei,xi)  		= tklo
        Aphi(age,ai,zi,lambdai,ei,xi)  		= tkhi        
        PrAprimelo(age,ai,zi,lambdai,ei,xi) = ( agrid(tkhi) - Aprime(age,ai,zi,lambdai,ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
        PrAprimehi(age,ai,zi,lambdai,ei,xi) = ( Aprime(age,ai,zi,lambdai,ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )

        if ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) .ge. amax) then
            tklo =na-1
        elseif ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) .lt. amin) then
            tklo = 1
        else
            tklo = (((1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi)-amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
        endif
        tkhi = tklo + 1        
        Bqlo(age,ai,zi,lambdai,ei,xi)  	= tklo
        Bqhi(age,ai,zi,lambdai,ei,xi)  	= tkhi        
        PrBqlo(age,ai,zi,lambdai,ei,xi) = &
        	& ( agrid(tkhi) - (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
        PrBqhi(age,ai,zi,lambdai,ei,xi) = &
        	& ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	! Probablities are adjusted to lie in [0,1]
		PrAprimelo = min(PrAprimelo, 1.0_DP); PrAprimelo = max(PrAprimelo, 0.0_DP)
		PrAprimehi = min(PrAprimehi, 1.0_DP); PrAprimehi = max(PrAprimehi, 0.0_DP)

		PrBqlo = min(PrBqlo, 1.0_DP); PrBqlo = max(PrBqlo, 0.0_DP)
		PrBqhi = min(PrBqhi, 1.0_DP); PrBqhi = max(PrBqhi, 0.0_DP)

	! Compute distribution of assets by age and state
		! Distribution is obtained by iterating over an initial distribution using policy functions
	DBN_dist=1.0_DP
	simutime = 1
	iter_indx = 1
	!print*, 'Computing Equilibrium Distribution'
	DO WHILE ( ( DBN_dist .ge. DBN_criteria ) .and. ( simutime .le. 700 ) )
		! print*, 'DBN_dist=', DBN_dist

	    DBN2=0.0_DP

		! Everyone in MaxAge dies. Those who die, switch to z2, lambda2 and start at ne/2+1 and x=1
	    age1=MaxAge
	    DO x1=1,nx
	    DO z1=1,nz
	    DO a1=1,na
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN2(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqhi(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO    
	    ENDDO
	    ENDDO    

		! retirees "e" stays the same for benefit retirement calculation purposes
	    DO x1=1,nx
	    DO age1=RetAge-1, MaxAge-1
	    DO a1=1,na
	    DO z1=1,nz
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        ! Those who die, switch to z2, lambda2 and start at ne/2+1
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN2(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqhi(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, and also e1 since they are retired
	        !e2=e1
	        DO x2=1,nx
		        DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) =  &
		        	& DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
		            & * survP(age1) * PrAprimelo(age1, a1, z1, lambda1, e1,x1) * pr_x(x1,x2,z1,age1)     
		        DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) =  &
		          	& DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
		            & * survP(age1) * PrAprimehi(age1, a1, z1, lambda1, e1,x1) * pr_x(x1,x2,z1,age1)  
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    
	    ! Working age agents
	    DO x1=1,nx
	    DO age1=1,RetAge-2
	    DO a1=1,na
	    DO z1=1,nz
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        ! Those who die, switch to z2, lambda2 and start at ne/2+1
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN2(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqhi(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, but switch to e2 and x2
	        DO x2=1,nx
	        DO e2=1,ne
				DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) =  &
	          		& DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
	                & * survP(age1) * pr_e(e1,e2) * pr_x(x1,x2,z1,age1) * PrAprimelo(age1, a1, z1, lambda1, e1,x1)     
	            DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) =  &
	          		& DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
	                & * survP(age1) * pr_e(e1,e2) * pr_x(x1,x2,z1,age1) * PrAprimehi(age1, a1, z1, lambda1, e1,x1) 
	        ENDDO
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    
	    DBN2 = 0.9*DBN1 + 0.1*DBN2
	    DBN_dist = maxval(abs(DBN2-DBN1))
	    ! print*, DBN_dist
	    DBN1 = DBN2

	    ! Instead of adjusting policy functions to be consistent with current DBN at each iteration
	    ! The code iterates the distribution and only updates policy functions every "update_period"


	    	!!
	    	print*, 'DBN_diff=', DBN_dist, 'R=',R,'P=',P,'Assets',sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid )
	    	!!
	    simutime  = simutime +1 
	 
	ENDDO ! WHILE

    	! ! Compute aggregates with current distribution
	    !     QBAR =0.0
	    !     NBAR =0.0
	    !     DO x1=1,nx
	    !     DO age1=1,MaxAge
	    !     DO z1=1,nz
	    !     DO a1=1,na
	    !     DO lambda1=1,nlambda
	    !     DO e1=1, ne
	    !          QBAR= QBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * ( xz_grid(x1,z1) * K_mat(a1,z1,x1) )**mu
	    !          NBAR= NBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * eff_un(age1, lambda1, e1) * Hours(age1, a1, z1, lambda1,e1,x1)
	    !     ENDDO
	    !     ENDDO
	    !     ENDDO
	    !     ENDDO    
	    !     ENDDO    
	    !     ENDDO    
	    
	    !     QBAR = ( QBAR)**(1.0_DP/mu)                
	    !     YBAR = QBAR ** alpha * NBAR **(1.0_DP-alpha)
	    !     Ebar = wage  * NBAR  * sum(pop)/sum(pop(1:RetAge-1))

    	! Deallocate policy functions on adjusted grid (so that they can be allocated later)
			! deallocate( YGRID_t  )
			! deallocate( MBGRID_t ) 
			! deallocate( Cons_t   )
			! deallocate( Hours_t  )
			! deallocate( Aprime_t )


END SUBROUTINE FIND_DBN_EQ_PF_Interp

!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE FIND_DBN_EQ_PF_Prices()
	use omp_lib
	IMPLICIT NONE
	INTEGER    :: tklo, tkhi, age1, age2, z1, z2, a1, a2, lambda1, lambda2, e1, e2, DBN_iter, simutime, iter_indx, x1, x2
	REAL       :: DBN_dist, DBN_criteria, tauW_aux
	real(dp)   :: BBAR, Wealth, brent_value
	REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne, nx) :: PrAprimelo, PrAprimehi, PrBqlo, PrBqhi, DBN2
	INTEGER,  DIMENSION(MaxAge, na, nz, nlambda, ne, nx) :: Aplo, Aphi, Bqlo, Bqhi
	real(dp), dimension(na,nz,nx) :: YGRID_bench
	real(dp), dimension(MaxAge, na, nz, nlambda, ne, nx) :: Aprime_aux, Cons_aux, Hours_aux 

	!$ call omp_set_num_threads(nz)
	DBN_criteria = 1.0E-07_DP
	! DBN1 = DBN_bench

	! Solve the model at current aggregate values
		! Find the threshold for wealth taxes (a_bar)
			!call Find_TauW_Threshold(DBN1,Y_a_threshold)
		! Adjust grid to include breaking points
			CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
		! Compute labor units 
			CALL ComputeLaborUnits(Ebar, wage) 
		! Compute Capital demand and Profits by (a,z)
			K_mat  = K_Matrix(R,P)
			Pr_mat = Profit_Matrix(R,P)
		! Compute policy functions and YGRID for benchmark economy
			tauW_aux = tauW_at
			tauW_at  = 0.00_dp 
			tauK     = 0.25_dp 
			CALL FORM_Y_MB_GRID(YGRID,MBGRID,YGRID_t,MBGRID_t)
			CALL EGM_RETIREMENT_WORKING_PERIOD 
			YGRID_bench = YGRID
			Aprime_aux  = Aprime 
			Cons_aux    = Cons 
			Hours_aux   = Hours 
			tauK     = 0.00_dp 
			tauW_at  = tauW_aux
			CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
			CALL FORM_Y_MB_GRID(YGRID,MBGRID,YGRID_t,MBGRID_t)
			deallocate(YGRID_t) ; deallocate(MBGRID_t)  
			deallocate(Cons_t)  ; deallocate(Hours_t)  ; deallocate( Aprime_t )			

	! Solve for policy and value functions 
			! Instead of EGM I update policy funtions by interpolating from the benchmark policy functions
	! Discretize policy function for assets (a')
		! For each age and state vector bracket optimal a' between two grid points
		! When at that age and state the optimal decision is approximated by selecting one the grid points
		! The grid points are selected with probability proportional to their distance to the optimal a'
	DO age=1,MaxAge
	!$omp parallel do private(lambdai,ei,ai,xi,tklo,tkhi)
	DO zi=1,nz
	DO xi=1,nx
	DO ai=1,na
	DO lambdai=1,nlambda
	DO ei=1, ne
		! Get Policy Functions
		Aprime(age,ai,zi,lambdai,ei,xi) = Linear_Int(YGRID_bench(:,zi,xi),Aprime_aux(age,:,zi,lambdai,ei,xi),na, YGRID(ai,zi,xi))
		Cons(age,ai,zi,lambdai,ei,xi)   = Linear_Int(YGRID_bench(:,zi,xi),Cons_aux(age,:,zi,lambdai,ei,xi)  ,na, YGRID(ai,zi,xi))
		Hours(age,ai,zi,lambdai,ei,xi)  = Linear_Int(YGRID_bench(:,zi,xi),Hours_aux(age,:,zi,lambdai,ei,xi) ,na, YGRID(ai,zi,xi)) 
		! ! Check
		! if (age.lt.RetAge) then 
		! print*, 'Residual',  YGRID(ai,zi,xi)  + Y_h(Hours(age, ai, zi, lambdai,ei,xi),age,lambdai,ei,wage)  & 
		! & - (1.0_dp+tauC)*Cons(age, ai, zi, lambdai,ei,xi) - Aprime(age, ai, zi, lambdai,ei,xi), &
		! & YGRID_bench(ai,zi,xi)  + Y_h(Hours_bench(age, ai, zi, lambdai,ei,xi),age,lambdai,ei,wage)  & 
		! & - (1.0_dp+tauC)*Cons_bench(age, ai, zi, lambdai,ei,xi) - Aprime_bench(age, ai, zi, lambdai,ei,xi),&
		! &'State', age,ai,zi,lambdai,ei,xi
  !       else 
  !       print*, 'Residual',  YGRID(ai,zi,xi)  + RetY_lambda_e(lambdai,ei)  & 
		! & - (1.0_dp+tauC)*Cons(age, ai, zi, lambdai,ei,xi) - Aprime(age, ai, zi, lambdai,ei,xi),&
		!  YGRID_bench(ai,zi,xi)  + RetY_lambda_e(lambdai,ei)  & 
		! & - (1.0_dp+tauC)*Cons_bench(age, ai, zi, lambdai,ei,xi) - Aprime_bench(age, ai, zi, lambdai,ei,xi),&
		! & 'State', age,ai,zi,lambdai,ei,xi
  !       endif 
        ! Discretize Aprime
        if ( Aprime(age,ai,zi,lambdai,ei,xi) .ge. amax) then
            tklo =na-1
        elseif (Aprime(age,ai,zi,lambdai, ei,xi) .lt. amin) then
            tklo = 1
        else
            tklo = ((Aprime(age,ai,zi,lambdai, ei,xi) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
        endif
        tkhi = tklo + 1        
        Aplo(age,ai,zi,lambdai,ei,xi)  		= tklo
        Aphi(age,ai,zi,lambdai,ei,xi)  		= tkhi        
        PrAprimelo(age,ai,zi,lambdai,ei,xi) = ( agrid(tkhi) - Aprime(age,ai,zi,lambdai,ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
        PrAprimehi(age,ai,zi,lambdai,ei,xi) = ( Aprime(age,ai,zi,lambdai,ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )

        if ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) .ge. amax) then
            tklo =na-1
        elseif ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) .lt. amin) then
            tklo = 1
        else
            tklo = (((1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi)-amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
        endif
        tkhi = tklo + 1        
        Bqlo(age,ai,zi,lambdai,ei,xi)  	= tklo
        Bqhi(age,ai,zi,lambdai,ei,xi)  	= tkhi        
        PrBqlo(age,ai,zi,lambdai,ei,xi) = &
        	& ( agrid(tkhi) - (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
        PrBqhi(age,ai,zi,lambdai,ei,xi) = &
        	& ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	! Probablities are adjusted to lie in [0,1]
		PrAprimelo = min(PrAprimelo, 1.0_DP); PrAprimelo = max(PrAprimelo, 0.0_DP)
		PrAprimehi = min(PrAprimehi, 1.0_DP); PrAprimehi = max(PrAprimehi, 0.0_DP)

		PrBqlo = min(PrBqlo, 1.0_DP); PrBqlo = max(PrBqlo, 0.0_DP)
		PrBqhi = min(PrBqhi, 1.0_DP); PrBqhi = max(PrBqhi, 0.0_DP)

	! Compute distribution of assets by age and state
		! Distribution is obtained by iterating over an initial distribution using policy functions
	DBN_dist=1.0_DP
	simutime = 1
	iter_indx = 1
	!print*, 'Computing Equilibrium Distribution'
	DO WHILE ( ( DBN_dist .ge. DBN_criteria ) .and. ( simutime .le. 700 ) )
		! print*, 'DBN_dist=', DBN_dist

	    DBN2=0.0_DP

		! Everyone in MaxAge dies. Those who die, switch to z2, lambda2 and start at ne/2+1 and x=1
	    age1=MaxAge
	    DO x1=1,nx
	    DO z1=1,nz
	    DO a1=1,na
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN2(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqhi(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO    
	    ENDDO
	    ENDDO    

		! retirees "e" stays the same for benefit retirement calculation purposes
	    DO x1=1,nx
	    DO age1=RetAge-1, MaxAge-1
	    DO a1=1,na
	    DO z1=1,nz
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        ! Those who die, switch to z2, lambda2 and start at ne/2+1
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN2(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqhi(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, and also e1 since they are retired
	        !e2=e1
	        DO x2=1,nx
		        DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) =  &
		        	& DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
		            & * survP(age1) * PrAprimelo(age1, a1, z1, lambda1, e1,x1) * pr_x(x1,x2,z1,age1)     
		        DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) =  &
		          	& DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
		            & * survP(age1) * PrAprimehi(age1, a1, z1, lambda1, e1,x1) * pr_x(x1,x2,z1,age1)  
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    
	    ! Working age agents
	    DO x1=1,nx
	    DO age1=1,RetAge-2
	    DO a1=1,na
	    DO z1=1,nz
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        ! Those who die, switch to z2, lambda2 and start at ne/2+1
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN2(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqhi(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, but switch to e2 and x2
	        DO x2=1,nx
	        DO e2=1,ne
				DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) =  &
	          		& DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
	                & * survP(age1) * pr_e(e1,e2) * pr_x(x1,x2,z1,age1) * PrAprimelo(age1, a1, z1, lambda1, e1,x1)     
	            DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) =  &
	          		& DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
	                & * survP(age1) * pr_e(e1,e2) * pr_x(x1,x2,z1,age1) * PrAprimehi(age1, a1, z1, lambda1, e1,x1) 
	        ENDDO
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    
	    DBN2 = 0.8*DBN1 + 0.2*DBN2
	    DBN_dist = maxval(abs(DBN2-DBN1))
	    ! print*, 'Iteration',simutime,'DBN_diff=', DBN_dist, 'R=',R,'P=',P
	    ! print*, DBN_dist
	    DBN1 = DBN2

	    ! Instead of adjusting policy functions to be consistent with current DBN at each iteration
	    ! The code iterates the distribution and only updates policy functions every "update_period"

	    ! Update of policy function with current aggregates
	    IF (iter_indx .ge. update_period) THEN

	    	! Compute aggregates with current distribution
	        QBAR =0.0
	        NBAR =0.0
	        DO x1=1,nx
	        DO age1=1,MaxAge
	        DO z1=1,nz
	        DO a1=1,na
	        DO lambda1=1,nlambda
	        DO e1=1, ne
	             QBAR= QBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * ( xz_grid(x1,z1) * K_mat(a1,z1,x1) )**mu
	             NBAR= NBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * eff_un(age1, lambda1, e1) * Hours(age1, a1, z1, lambda1,e1,x1)
	        ENDDO
	        ENDDO
	        ENDDO
	        ENDDO    
	        ENDDO    
	        ENDDO    
	        QBAR = ( QBAR)**(1.0_DP/mu) 

	        if (A_C.gt.0.0_dp) then 

	        	! Update Interest Rate (using old P)
	        	brent_value = brent(-0.1_DP,0.01_DP,10.0_DP,Agg_Debt_C, brent_tol,R)

		    	! Total wealth 
		    	Wealth = sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid )

	        	! Corporate Capital
        		K_C    = sum( sum(sum(sum(sum(sum(DBN1(:,:,z_C:,:,:,:),6),5),4),3),1) * agrid ) 
        		K_P    = Wealth - K_C

	        	! Update Wage 
	        	if (alpha_C.eq.alpha) then 
	        		Wage = ((((1.0_dp-alpha)*Aprod)**(1.0_dp/alpha)*QBAR + ((1.0_dp-alpha)*A_C)**(1.0_dp/alpha)*K_C)/NBAR)**(alpha)
	        	else 
	        		print*, ' Solving labor market numerically'
	        		brent_value = brent(0.001_DP,Wage,10.0_DP,Labor_Market_Clearing, brent_tol,Wage)
	        		print*, ' New Wage=',Wage,'Error=',brent_value
	        	endif 

	        	! Update Other Corporate values 
		    	L_C    = ( (1.0_dp-alpha_C)*A_C/Wage )**(1.0_dp/alpha_C) * K_C
		    	R_C    = alpha_C * A_C * (L_C/K_C)**(1.0_dp-alpha_C) - DepRate
		    	

	    		! Check that R_C>R (if not solve the equilibrium with a single interest rate)
		    	if (R_C.lt.R) then 
		        	! Private demand for capital
		        	K_P    = sum( (sum(sum(sum(DBN1,5),4),1)) *(K_matrix(R,P))) ! Note: DBN_azx  = sum(sum(sum(DBN1,5),4),1)
	        	
		        	if (Wealth.gt.(K_P+Debt_Absorption*Debt_SS)) then 
			        	! Corporate demand for capital
			        	K_C    = Wealth - K_P - Debt_Absorption*Debt_SS

			        	! Update Wage 
			        	if (alpha_C.eq.alpha) then 
			        		Wage = ((((1.0_dp-alpha)*Aprod)**(1.0_dp/alpha)*QBAR + ((1.0_dp-alpha)*A_C)**(1.0_dp/alpha)*K_C)/NBAR)**(alpha)
			        	else 
			        		print*, ' Solving labor market numerically'
			        		brent_value = brent(0.001_DP,Wage,10.0_DP,Labor_Market_Clearing, brent_tol,Wage)
			        		print*, ' New Wage=',Wage,'Error=',brent_value
			        	endif 

		        	else ! Capital Market did not clear

		        		print*, ' ';print*, ' 	Warning! Capital Market Did not Clear!';print*, ' ';
			        	! Modify prices and quantities 
			        	QBAR = 0.50_dp*QBAR 
			        	Wage = Wage 

		        	endif 

		        	! Update Interest Rates 
		        	R    = alpha_C * A_C * ( Wage/((1.0_dp-alpha_C)*A_C) )**(-(1.0_dp-alpha_C)/alpha_C) - DepRate
					R_C  = R 

		    	endif 


		    	! Update Output and Aggregates 
				L_P    = ( (1.0_dp-alpha)*Aprod/Wage )**(1.0_dp/alpha  ) * QBAR 
				L_C    = ( (1.0_dp-alpha_C)*A_C/Wage )**(1.0_dp/alpha_C) * K_C
				P      = alpha * QBAR**(alpha-mu) * L_P**(1.0_DP-alpha)
				Ebar   = wage  * NBAR  * sum(pop)/sum(pop(1:RetAge-1))
				YBAR_P = AProd * QBAR**alpha   * L_P**(1.0_DP-alpha  ) 
				YBAR_C = A_C   * K_C **alpha_C * L_C**(1.0_DP-alpha_C) 
				YBAR   = YBAR_P + YBAR_C

	        	
	        	! Update interest rate vector 
        		R_z(1:z_C-1) = R 
        		R_z(z_C:)    = R_C

				! print '(A,F9.3,F9.3,F9.3,X,X,A,F9.3,F9.3,F9.3)',&
				! 		& ' 				Corp. Sector Levels:', YBAR_C, K_C, L_C , &
				! 		& ' Ratios ', 100.0_dp*YBAR_C/YBAR, 100.0_dp*K_C/Wealth, 100.0_dp*L_C/NBAR

	        else 
	        	! Set Corporate values to zero
        		K_C    = 0.0_dp; L_C    = 0.0_dp; YBAR_C = 0.0_dp
        		K_P    = 0.0_dp; L_P    = 0.0_dp; YBAR_P = 0.0_dp

	        	! Solve for aggregates and clear capital market with R
		        P    = alpha* QBAR **(alpha-mu) * NBAR **(1.0_DP-alpha)
		        YBAR = QBAR ** alpha * NBAR **(1.0_DP-alpha)
		        wage = (1.0_DP-alpha)*QBAR **alpha * NBAR  **(-alpha)
		        Ebar = wage  * NBAR  * sum(pop)/sum(pop(1:RetAge-1))

		    	! Solve for new R 
		    	! R = zbrent(Agg_Debt,0.1_dp,1.00_dp,brent_tol) 
		    	if (sum(theta)/nz .gt. 1.0_DP) then
		    		P = min(P,1.0_dp)
		            brent_value = brent(-0.1_DP,0.01_DP,10.0_DP,Agg_Debt, brent_tol,R)
		        else
		            R = 0.0_DP
		        endif

	        endif 

	    	!!
	    	print 12345, &
	    		& ' DBN_diff=', DBN_dist,'A=',sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid ),&
	    		& 'W=',wage,'R=',R,'P=',P,'Q=',QBAR, &
	    		& 'K_C/A=',100.0_dp*K_C/Wealth,'L_C/N=',100.0_dp*L_C/NBAR,'K_C=',K_C,'L_C=',L_C
    		12345 format &
    		&(A,E12.5,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X)
	    	!!

	    	! Solve the model at current aggregate values
				! Find the threshold for wealth taxes (a_bar)
					!call Find_TauW_Threshold(DBN1,Y_a_threshold)
				! Adjust grid to include breaking points
					CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
				! Compute labor units 
					CALL ComputeLaborUnits(Ebar, wage) 
				! Compute Capital demand and Profits by (a,z)
					K_mat  = K_Matrix(R,P)
					Pr_mat = Profit_Matrix(R,P)
				! Solve for policy and value functions 
					! Policy functions are computed for the capital tax economy 
					tauW_aux = tauW_at
					tauW_at  = 0.00_dp 
					tauK     = 0.25_dp 
					CALL FORM_Y_MB_GRID(YGRID,MBGRID,YGRID_t,MBGRID_t)
					CALL EGM_RETIREMENT_WORKING_PERIOD 
					YGRID_bench = YGRID
					Aprime_aux  = Aprime 
					Cons_aux    = Cons 
					Hours_aux   = Hours 
					tauK     = 0.00_dp 
					tauW_at  = tauW_aux
					CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
					CALL FORM_Y_MB_GRID(YGRID,MBGRID,YGRID_t,MBGRID_t)
					deallocate(YGRID_t) ; deallocate(MBGRID_t)
					deallocate(Cons_t)  ; deallocate(Hours_t)  ; deallocate( Aprime_t )
	        
				! Solve for policy and value functions 
						! Instead of EGM I update policy funtions by interpolating from the benchmark policy functions
				! Discretize policy function for assets (a')
					! For each age and state vector bracket optimal a' between two grid points
					! When at that age and state the optimal decision is approximated by selecting one the grid points
					! The grid points are selected with probability proportional to their distance to the optimal a'
				DO age=1,MaxAge
				!$omp parallel do private(lambdai,ei,ai,xi,tklo,tkhi)
				DO zi=1,nz
				DO xi=1,nx
				DO ai=1,na
				DO lambdai=1,nlambda
				DO ei=1, ne
					! Get Policy Functions
					Aprime(age,ai,zi,lambdai,ei,xi) = Linear_Int(YGRID_bench(:,zi,xi),Aprime_aux(age,:,zi,lambdai,ei,xi),na, YGRID(ai,zi,xi))
					Cons(age,ai,zi,lambdai,ei,xi)   = Linear_Int(YGRID_bench(:,zi,xi),Cons_aux(age,:,zi,lambdai,ei,xi)  ,na, YGRID(ai,zi,xi))
					Hours(age,ai,zi,lambdai,ei,xi)  = Linear_Int(YGRID_bench(:,zi,xi),Hours_aux(age,:,zi,lambdai,ei,xi) ,na, YGRID(ai,zi,xi)) 
			        ! Discretize Aprime
			        if ( Aprime(age,ai,zi,lambdai,ei,xi) .ge. amax) then
			            tklo =na-1
			        elseif (Aprime(age,ai,zi,lambdai, ei,xi) .lt. amin) then
			            tklo = 1
			        else
			            tklo = ((Aprime(age,ai,zi,lambdai, ei,xi) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
			        endif
			        tkhi = tklo + 1        
			        Aplo(age,ai,zi,lambdai,ei,xi)  		= tklo
			        Aphi(age,ai,zi,lambdai,ei,xi)  		= tkhi        
			        PrAprimelo(age,ai,zi,lambdai,ei,xi) = &
			        	& ( agrid(tkhi) - Aprime(age,ai,zi,lambdai,ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
			        PrAprimehi(age,ai,zi,lambdai,ei,xi) = &
			        	& ( Aprime(age,ai,zi,lambdai,ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )

			        if ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) .ge. amax) then
			            tklo =na-1
			        elseif ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) .lt. amin) then
			            tklo = 1
			        else
			            tklo = (((1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi)-amin)/(amax-amin))&
			            	&**(1.0_DP/a_theta)*(na-1)+1          
			        endif
			        tkhi = tklo + 1        
			        Bqlo(age,ai,zi,lambdai,ei,xi)  	= tklo
			        Bqhi(age,ai,zi,lambdai,ei,xi)  	= tkhi        
			        PrBqlo(age,ai,zi,lambdai,ei,xi) = ( agrid(tkhi)-(1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) )&
			        									& / ( agrid(tkhi) -agrid(tklo) )
			        PrBqhi(age,ai,zi,lambdai,ei,xi) = ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi)-agrid(tklo) )&
			        									&  / ( agrid(tkhi) -agrid(tklo) )
				ENDDO
				ENDDO
				ENDDO
				ENDDO
				ENDDO
				ENDDO
	
	        ! Probablities are adjusted to lie in [0,1]
					PrAprimelo = min(PrAprimelo, 1.0_DP); PrAprimelo = max(PrAprimelo, 0.0_DP)
					PrAprimehi = min(PrAprimehi, 1.0_DP); PrAprimehi = max(PrAprimehi, 0.0_DP)

					PrBqlo = min(PrBqlo, 1.0_DP); PrBqlo = max(PrBqlo, 0.0_DP)
					PrBqhi = min(PrBqhi, 1.0_DP); PrBqhi = max(PrBqhi, 0.0_DP) 

		    ! Reset counter for next update of policy functions
	        iter_indx=0
	    ENDIF
	    
	    iter_indx = iter_indx + 1
	    simutime  = simutime +1 
	 
	ENDDO ! WHILE


END SUBROUTINE FIND_DBN_EQ_PF_Prices

!========================================================================================
!========================================================================================
!========================================================================================

SUBROUTINE FIND_DBN_EQ_Prices(Fixed_W,Fixed_P,Fixed_R)
	use omp_lib
	IMPLICIT NONE
	logical, intent(in) :: Fixed_W, Fixed_P, Fixed_R
	INTEGER:: tklo, tkhi, age1, age2, z1, z2, a1, a2, lambda1, lambda2, e1, e2, DBN_iter, simutime, iter_indx, x1, x2
	REAL   :: DBN_dist, DBN_criteria
	real(dp)   ::BBAR, Wealth, brent_value
	REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne, nx) :: PrAprimelo, PrAprimehi, PrBqlo, PrBqhi, DBN2
	INTEGER,  DIMENSION(MaxAge, na, nz, nlambda, ne, nx) :: Aplo, Aphi, Bqlo, Bqhi

	!$ call omp_set_num_threads(nz)
	DBN_criteria = 1.0E-07_DP
	DBN1 = DBN_bench

	! Solve the model at current aggregate values
		! Find the threshold for wealth taxes (a_bar)
			!call Find_TauW_Threshold(DBN1,Y_a_threshold)
		! Adjust grid to include breaking points
			CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
		! Compute labor units 
			CALL ComputeLaborUnits(Ebar, wage) 
		! Compute Capital demand and Profits by (a,z)
			K_mat  = K_Matrix(R,P)
			Pr_mat = Profit_Matrix(R,P)
		! Form YGRID for the capital income economy given interest rate "P"
			CALL FORM_Y_MB_GRID(YGRID,MBGRID,YGRID_t,MBGRID_t)
		! Solve for policy and value functions 
			CALL EGM_RETIREMENT_WORKING_PERIOD 

	! Discretize policy function for assets (a')
		! For each age and state vector bracket optimal a' between two grid points
		! When at that age and state the optimal decision is approximated by selecting one the grid points
		! The grid points are selected with probability proportional to their distance to the optimal a'
	DO age=1,MaxAge
	!$omp parallel do private(lambdai,ei,ai,xi,tklo,tkhi)
	DO zi=1,nz
	DO xi=1,nx
	DO ai=1,na
	DO lambdai=1,nlambda
	DO ei=1, ne
        if ( Aprime(age,ai,zi,lambdai,ei,xi) .ge. amax) then
            tklo =na-1
        elseif (Aprime(age,ai,zi,lambdai, ei,xi) .lt. amin) then
            tklo = 1
        else
            tklo = ((Aprime(age,ai,zi,lambdai, ei,xi) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
        endif
        tkhi = tklo + 1        
        Aplo(age,ai,zi,lambdai,ei,xi)  		= tklo
        Aphi(age,ai,zi,lambdai,ei,xi)  		= tkhi        
        PrAprimelo(age,ai,zi,lambdai,ei,xi) = ( agrid(tkhi) - Aprime(age,ai,zi,lambdai,ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
        PrAprimehi(age,ai,zi,lambdai,ei,xi) = ( Aprime(age,ai,zi,lambdai,ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )

        if ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) .ge. amax) then
            tklo =na-1
        elseif ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) .lt. amin) then
            tklo = 1
        else
            tklo = (((1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi)-amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
        endif
        tkhi = tklo + 1        
        Bqlo(age,ai,zi,lambdai,ei,xi)  	= tklo
        Bqhi(age,ai,zi,lambdai,ei,xi)  	= tkhi        
        PrBqlo(age,ai,zi,lambdai,ei,xi) = &
        	& ( agrid(tkhi) - (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
        PrBqhi(age,ai,zi,lambdai,ei,xi) = &
        	& ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	! Probablities are adjusted to lie in [0,1]
		PrAprimelo = min(PrAprimelo, 1.0_DP); PrAprimelo = max(PrAprimelo, 0.0_DP)
		PrAprimehi = min(PrAprimehi, 1.0_DP); PrAprimehi = max(PrAprimehi, 0.0_DP)

		PrBqlo = min(PrBqlo, 1.0_DP); PrBqlo = max(PrBqlo, 0.0_DP)
		PrBqhi = min(PrBqhi, 1.0_DP); PrBqhi = max(PrBqhi, 0.0_DP)

	! Compute distribution of assets by age and state
		! Distribution is obtained by iterating over an initial distribution using policy functions
	DBN_dist=1.0_DP
	simutime = 1
	iter_indx = 1
	!print*, 'Computing Equilibrium Distribution'
	DO WHILE ( ( DBN_dist .ge. DBN_criteria ) .and. ( simutime .le. 700 ) )
		! print*, 'DBN_dist=', DBN_dist

	    DBN2=0.0_DP

		! Everyone in MaxAge dies. Those who die, switch to z2, lambda2 and start at ne/2+1 and x=1
	    age1=MaxAge
	    DO x1=1,nx
	    DO z1=1,nz
	    DO a1=1,na
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN2(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqhi(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO    
	    ENDDO
	    ENDDO    

		! retirees "e" stays the same for benefit retirement calculation purposes
	    DO x1=1,nx
	    DO age1=RetAge-1, MaxAge-1
	    DO a1=1,na
	    DO z1=1,nz
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        ! Those who die, switch to z2, lambda2 and start at ne/2+1
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN2(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqhi(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, and also e1 since they are retired
	        !e2=e1
	        DO x2=1,nx
		        DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) =  &
		        	& DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
		            & * survP(age1) * PrAprimelo(age1, a1, z1, lambda1, e1,x1) * pr_x(x1,x2,z1,age1)     
		        DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) =  &
		          	& DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
		            & * survP(age1) * PrAprimehi(age1, a1, z1, lambda1, e1,x1) * pr_x(x1,x2,z1,age1)  
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    
	    ! Working age agents
	    DO x1=1,nx
	    DO age1=1,RetAge-2
	    DO a1=1,na
	    DO z1=1,nz
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        ! Those who die, switch to z2, lambda2 and start at ne/2+1
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqlo(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN2(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN2(1, Bqhi(age1,a1,z1,lambda1,e1,x1) ,z2,lambda2,ne/2+1,1) + DBN1(age1,a1,z1,lambda1,e1,x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, but switch to e2 and x2
	        DO x2=1,nx
	        DO e2=1,ne
				DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) =  &
	          		& DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
	                & * survP(age1) * pr_e(e1,e2) * pr_x(x1,x2,z1,age1) * PrAprimelo(age1, a1, z1, lambda1, e1,x1)     
	            DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) =  &
	          		& DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) + DBN1(age1, a1, z1, lambda1, e1,x1) &
	                & * survP(age1) * pr_e(e1,e2) * pr_x(x1,x2,z1,age1) * PrAprimehi(age1, a1, z1, lambda1, e1,x1) 
	        ENDDO
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    
	    DBN2 = 0.1*DBN1 + 0.9*DBN2
	    DBN_dist = maxval(abs(DBN2-DBN1))
	    ! print*, DBN_dist
	    DBN1 = DBN2

	    ! Instead of adjusting policy functions to be consistent with current DBN at each iteration
	    ! The code iterates the distribution and only updates policy functions every "update_period"

	    ! Update of policy function with current aggregates
	    IF (iter_indx .ge. update_period) THEN

	    	! Compute aggregates with current distribution
	        QBAR =0.0
	        NBAR =0.0
	        DO x1=1,nx
	        DO age1=1,MaxAge
	        DO z1=1,nz
	        DO a1=1,na
	        DO lambda1=1,nlambda
	        DO e1=1, ne
	             QBAR= QBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * ( xz_grid(x1,z1) * K_mat(a1,z1,x1) )**mu
	             NBAR= NBAR+ DBN1(age1, a1, z1, lambda1, e1, x1) * eff_un(age1, lambda1, e1) * Hours(age1, a1, z1, lambda1,e1,x1)
	        ENDDO
	        ENDDO
	        ENDDO
	        ENDDO    
	        ENDDO    
	        ENDDO    
	        QBAR = ( QBAR)**(1.0_DP/mu)                

	        if (A_C.gt.0.0_dp) then 

	        	! Update Interest Rate (using old P)
	        	brent_value = brent(-0.1_DP,0.01_DP,10.0_DP,Agg_Debt_C, brent_tol,R)

		    	! Total wealth 
		    	Wealth = sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid )

	        	! Corporate Capital
        		K_C    = sum( sum(sum(sum(sum(sum(DBN1(:,:,z_C:,:,:,:),6),5),4),3),1) * agrid ) 
        		K_P    = Wealth - K_C

	        	! Update Wage 
	        	if (alpha_C.eq.alpha) then 
	        		Wage = ((((1.0_dp-alpha)*Aprod)**(1.0_dp/alpha)*QBAR + ((1.0_dp-alpha)*A_C)**(1.0_dp/alpha)*K_C)/NBAR)**(alpha)
	        	else 
	        		print*, ' Solving labor market numerically'
	        		brent_value = brent(0.001_DP,Wage,10.0_DP,Labor_Market_Clearing, brent_tol,Wage)
	        		print*, ' New Wage=',Wage,'Error=',brent_value
	        	endif 

	        	! Update Other Corporate values 
		    	L_C    = ( (1.0_dp-alpha_C)*A_C/Wage )**(1.0_dp/alpha_C) * K_C
		    	R_C    = alpha_C * A_C * (L_C/K_C)**(1.0_dp-alpha_C) - DepRate
		    	

	    		! Check that R_C>R (if not solve the equilibrium with a single interest rate)
		    	if (R_C.lt.R) then 
		        	! Private demand for capital
		        	K_P    = sum( (sum(sum(sum(DBN1,5),4),1)) *(K_matrix(R,P))) ! Note: DBN_azx  = sum(sum(sum(DBN1,5),4),1)
	        	
		        	if (Wealth.gt.(K_P+Debt_Absorption*Debt_SS)) then 
			        	! Corporate demand for capital
			        	K_C    = Wealth - K_P - Debt_Absorption*Debt_SS

			        	! Update Wage 
			        	if (alpha_C.eq.alpha) then 
			        		Wage = ((((1.0_dp-alpha)*Aprod)**(1.0_dp/alpha)*QBAR + ((1.0_dp-alpha)*A_C)**(1.0_dp/alpha)*K_C)/NBAR)**(alpha)
			        	else 
			        		print*, ' Solving labor market numerically'
			        		brent_value = brent(0.001_DP,Wage,10.0_DP,Labor_Market_Clearing, brent_tol,Wage)
			        		print*, ' New Wage=',Wage,'Error=',brent_value
			        	endif 

		        	else ! Capital Market did not clear

		        		print*, ' ';print*, ' 	Warning! Capital Market Did not Clear!';print*, ' ';
			        	! Modify prices and quantities 
			        	QBAR = 0.50_dp*QBAR 
			        	Wage = Wage 

		        	endif 

		        	! Update Interest Rates 
		        	R    = alpha_C * A_C * ( Wage/((1.0_dp-alpha_C)*A_C) )**(-(1.0_dp-alpha_C)/alpha_C) - DepRate
					R_C  = R 

		    	endif 


		    	! Update Output and Aggregates 
				L_P    = ( (1.0_dp-alpha)*Aprod/Wage )**(1.0_dp/alpha  ) * QBAR 
				L_C    = ( (1.0_dp-alpha_C)*A_C/Wage )**(1.0_dp/alpha_C) * K_C
				P      = alpha * QBAR**(alpha-mu) * L_P**(1.0_DP-alpha)
				Ebar   = wage  * NBAR  * sum(pop)/sum(pop(1:RetAge-1))
				YBAR_P = AProd * QBAR**alpha   * L_P**(1.0_DP-alpha  ) 
				YBAR_C = A_C   * K_C **alpha_C * L_C**(1.0_DP-alpha_C) 
				YBAR   = YBAR_P + YBAR_C

	        	
	        	! Update interest rate vector 
        		R_z(1:z_C-1) = R 
        		R_z(z_C:)    = R_C

				! print '(A,F9.3,F9.3,F9.3,X,X,A,F9.3,F9.3,F9.3)',&
				! 		& ' 				Corp. Sector Levels:', YBAR_C, K_C, L_C , &
				! 		& ' Ratios ', 100.0_dp*YBAR_C/YBAR, 100.0_dp*K_C/Wealth, 100.0_dp*L_C/NBAR

	        else 
	        	! Set Corporate values to zero
        		K_C    = 0.0_dp; L_C    = 0.0_dp; YBAR_C = 0.0_dp
        		K_P    = 0.0_dp; L_P    = 0.0_dp; YBAR_P = 0.0_dp

		        YBAR = QBAR ** alpha * NBAR **(1.0_DP-alpha)

		        if (Fixed_W.eqv..false.) wage = (1.0_DP-alpha)*QBAR **alpha * NBAR  **(-alpha)
		        if (Fixed_P.eqv..false.) P    = alpha* QBAR **(alpha-mu) * NBAR **(1.0_DP-alpha)

		        Ebar = wage  * NBAR  * sum(pop)/sum(pop(1:RetAge-1))

		        if (Fixed_R.eqv..false.) then 
		    	! Solve for new R 
		    	! R = zbrent(Agg_Debt,0.1_dp,1.00_dp,brent_tol) 
		    	if (sum(theta)/nz .gt. 1.0_DP) then
		    		P = min(P,1.0_dp)
		           brent_value = brent(-0.1_DP,0.01_DP,10.0_DP,Agg_Debt, brent_tol,R)
	            else
	                R = 0.0_DP
		        endif
	        endif 


	        endif 

	    	!!
	    	print 12345, &
	    		& ' DBN_diff=', DBN_dist,'A=',sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid ),&
	    		& 'W=',wage,'R=',R,'P=',P,'Q=',QBAR, &
	    		& 'K_C/A=',100.0_dp*K_C/Wealth,'L_C/N=',100.0_dp*L_C/NBAR,'K_C=',K_C,'L_C=',L_C
    		12345 format &
    		&(A,E12.5,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X)
	    	!!

	    	! Solve the model at current aggregate values
				! Find the threshold for wealth taxes (a_bar)
					!call Find_TauW_Threshold(DBN1,Y_a_threshold)
				! Adjust grid to include breaking points
					CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
				! Compute labor units 
					CALL ComputeLaborUnits(Ebar, wage) 
				! Compute Capital demand and Profits by (a,z)
					K_mat  = K_Matrix(R,P)
					Pr_mat = Profit_Matrix(R,P)
				! Form YGRID for the capital income economy given interest rate "P"
					CALL FORM_Y_MB_GRID(YGRID,MBGRID,YGRID_t,MBGRID_t)
				! Solve for policy and value functions 
					CALL EGM_RETIREMENT_WORKING_PERIOD 
	        
				! Discretize policy function for assets (a')
					! For each age and state vector bracket optimal a' between two grid points
					! When at that age and state the optimal decision is approximated by selecting one the grid points
					! The grid points are selected with probability proportional to their distance to the optimal a'
				DO age=1,MaxAge
				!$omp parallel do private(lambdai,ei,ai,xi,tklo,tkhi)
				DO zi=1,nz
				DO xi=1,nx
				DO ai=1,na
				DO lambdai=1,nlambda
				DO ei=1, ne
			        if ( Aprime(age,ai,zi,lambdai,ei,xi) .ge. amax) then
			            tklo =na-1
			        elseif (Aprime(age,ai,zi,lambdai, ei,xi) .lt. amin) then
			            tklo = 1
			        else
			            tklo = ((Aprime(age,ai,zi,lambdai, ei,xi) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
			        endif            
			        tkhi = tklo + 1        
			        Aplo(age,ai,zi,lambdai,ei,xi)  		= tklo
			        Aphi(age,ai,zi,lambdai,ei,xi)  		= tkhi        
			        PrAprimelo(age,ai,zi,lambdai,ei,xi) = ( agrid(tkhi) - Aprime(age,ai,zi,lambdai,ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
			        PrAprimehi(age,ai,zi,lambdai,ei,xi) = ( Aprime(age,ai,zi,lambdai,ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )

			        if ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) .ge. amax) then
			            tklo =na-1
			        elseif ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) .lt. amin) then
			            tklo = 1
			        else
			            tklo = (((1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi)-amin)/(amax-amin))&
			            	&**(1.0_DP/a_theta)*(na-1)+1          
			        endif
			        tkhi = tklo + 1        
			        Bqlo(age,ai,zi,lambdai,ei,xi)  	= tklo
			        Bqhi(age,ai,zi,lambdai,ei,xi)  	= tkhi        
			        PrBqlo(age,ai,zi,lambdai,ei,xi) = ( agrid(tkhi)-(1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi) )&
			        									& / ( agrid(tkhi) -agrid(tklo) )
			        PrBqhi(age,ai,zi,lambdai,ei,xi) = ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime(age,ai,zi,lambdai,ei,xi)-agrid(tklo) )&
			        									&  / ( agrid(tkhi) -agrid(tklo) )
				ENDDO
				ENDDO
				ENDDO
				ENDDO
				ENDDO
				ENDDO
	
	        ! Probablities are adjusted to lie in [0,1]
					PrAprimelo = min(PrAprimelo, 1.0_DP); PrAprimelo = max(PrAprimelo, 0.0_DP)
					PrAprimehi = min(PrAprimehi, 1.0_DP); PrAprimehi = max(PrAprimehi, 0.0_DP)

					PrBqlo = min(PrBqlo, 1.0_DP); PrBqlo = max(PrBqlo, 0.0_DP)
					PrBqhi = min(PrBqhi, 1.0_DP); PrBqhi = max(PrBqhi, 0.0_DP)   

		    ! Reset counter for next update of policy functions
	        iter_indx=0
	    ENDIF
	    
	    iter_indx = iter_indx + 1
	    simutime  = simutime +1 
	 
	ENDDO ! WHILE

	print*,' '
	print*,' 	Stationary Equilibrium Found: '
	print 12346, &
		& ' 	DBN_diff=', DBN_dist,'A=',sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid ),&
		& 'W=',wage,'R=',R,'P=',P,'Q=',QBAR, &
		& 'K_C/A=',100.0_dp*K_C/Wealth,'L_C/N=',100.0_dp*L_C/NBAR,'Y_C/Y=',100.0_dp*YBAR_C/YBAR,'Iter=',simutime
	12346 format &
	& (A,E12.5,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,I5)
	print*,' '; print*,'-----------------------------------------------------------------------------'
	print*,' '



END SUBROUTINE FIND_DBN_EQ_Prices



!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE FIND_DBN_Transition()
	use omp_lib
	IMPLICIT NONE
	INTEGER    :: tklo, tkhi, age1, age2, z1, z2, a1, a2, lambda1, lambda2, e1, e2, DBN_iter, simutime, iter_indx, x1, x2
	REAL       :: DBN_dist, DBN_criteria, Q_dist, NW_dist, Price_criteria, Chg_criteria, Old_DBN_dist, Chg_dist, R_dist, Db_dist
	REAL(DP)   :: BBAR, Wealth, brent_value, K_bench, C_bench, K_exp, C_exp, R_old, Db_old, R_slope 
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: PrAprimelo, PrAprimehi, PrBqlo, PrBqhi
	INTEGER , DIMENSION(:,:,:,:,:,:), allocatable :: Aplo, Aphi, Bqlo, Bqhi
	REAL(DP), DIMENSION(T+1) :: QBAR2_tr, NBAR2_tr, wage2_tr, Wealth_Top_1_Tr, Wealth_Top_10_Tr, R2_tr
	INTEGER    :: prctile, ind_R 
	REAL(DP)   :: Dampen
	! Timing
	real(kind=8)    :: t1, t2, elapsed_time, t1_R, t2_R, elapsed_time_R
    integer(kind=8) :: tclock1, tclock2, clock_rate, tclock1_R, tclock2_R, clock_rate_R


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Set Up
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Allocate
	allocate( PrAprimelo(   MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrAprimehi(   MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Aplo(         MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Aphi(         MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrBqlo(   	MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrBqhi(   	MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Bqlo(         MaxAge,na,nz,nlambda,ne,nx) )
	allocate( BQhi(         MaxAge,na,nz,nlambda,ne,nx) )

	!$ call omp_set_num_threads(nz)
	DBN_criteria    = 1.0E-05_DP
	Price_criteria  = 1.0E-04_DP
	Chg_criteria    = 5.0E-06_DP
	ind_R   		= 3
	Dampen   	    = 0.70_dp


	! Wealth and consumption in benchmark and experiment
		K_bench = sum( sum(sum(sum(sum(sum(DBN_bench,6),5),4),3),1)*agrid )
		C_bench = sum( DBN_bench*Cons_bench )


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Initial guess for transition path variables
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (Use_Transition_Seed.eqv..false.) then 
		! Guess NBAR, QBAR and R as a linear combination of starting and end values
			if (A_C.gt.0.0_dp) then 
			wage_tr(1)   = wage_bench ; wage_tr(T+1) = wage_exp   ;
			else
			NBAR_tr(1)   = NBAR_bench ; NBAR_tr(T+1) = NBAR_exp   ;
			endif 
			QBAR_tr(1)   = QBAR_bench ; QBAR_tr(T+1) = QBAR_exp   ;
			R_tr(1)      = R_bench    ; R_tr(T+1)    = R_exp      ;
			do ti=2,T
				if (A_C.gt.0.0_dp) then 
				wage_tr(ti) = wage_tr(ti-1) + (wage_tr(T+1)-wage_tr(1))/T
				else
				NBAR_tr(ti) = NBAR_tr(ti-1) + (NBAR_tr(T+1)-NBAR_tr(1))/T
				endif 
				QBAR_tr(ti) = QBAR_tr(ti-1) + (QBAR_tr(T+1)-QBAR_tr(1))/T
				R_tr(ti)    = R_tr(ti-1) + (R_tr(T+1)-R_tr(1))/T
			enddo 
			Debt_SS = 0.0_dp
			Use_Transition_Seed = .true.
		else
		! Load Guess From Files
			print*, 'Loading initial variables from file'
			OPEN (UNIT=1,  FILE=trim(Result_Folder)//'QBAR_tr'   , STATUS='old', ACTION='read')
			OPEN (UNIT=3,  FILE=trim(Result_Folder)//'R_tr'		 , STATUS='old', ACTION='read')
			OPEN (UNIT=4,  FILE=trim(Result_Folder)//'Debt_SS' 	 , STATUS='old', ACTION='read')
			OPEN (UNIT=5,  FILE=trim(Result_Folder)//'DBN_SS' 	 , STATUS='old', ACTION='read')
			READ (UNIT=1,  FMT=*) QBAR_tr
			READ (UNIT=3,  FMT=*) R_tr
			READ (UNIT=4,  FMT=*) Debt_SS
			READ (UNIT=5,  FMT=*) DBN1
			CLOSE (unit=1); CLOSE (unit=3); CLOSE (unit=4); CLOSE (unit=5);

			if (A_C.gt.0.0_dp) then 
			OPEN (UNIT=2,  FILE=trim(Result_Folder)//'Wage_tr'	 , STATUS='old', ACTION='read')
			READ (UNIT=2,  FMT=*) wage_tr
			else
			OPEN (UNIT=2,  FILE=trim(Result_Folder)//'NBAR_tr'	 , STATUS='old', ACTION='read')
			READ (UNIT=2,  FMT=*) NBAR_tr
			endif 
			CLOSE (unit=2); 
			print*, 'Reading completed'
		endif 


		! Current aggregate values given QBAR and Wage
		if (A_C.gt.0.0_dp) then 
			L_P_tr    = ( (1.0_dp-alpha)*Aprod/Wage_tr )**(1.0_dp/alpha) * QBAR_tr 
			P_tr      = alpha * QBAR_tr**(alpha-mu) * L_P_tr**(1.0_DP-alpha)

			K_tr 	  = sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid )
	    	K_P_tr    = K_P
	    	K_C_tr    = K_tr - K_P_tr 
	    	L_C_tr    = ( (1.0_dp-alpha_C)*A_C/Wage_tr )**(1.0_dp/alpha_C) * K_C_tr
			YBAR_P_tr = AProd * QBAR_tr**alpha   * L_P_tr**(1.0_DP-alpha  ) 
			YBAR_C_tr = A_C   * K_C_tr **alpha_C * L_C_tr**(1.0_DP-alpha_C) 
			YBAR_tr   = YBAR_P_tr + YBAR_C_tr
			NBAR_tr   =    L_P_tr +    L_C_tr

			! print*,'initial values:',Wealth,K_P,K_C,L_C,L_P,YBAR,YBAR_P,YBAR_C,P,R
		else 
			! Choose YBAR, EBAR, P and Wage to be consistent
			P_tr    = alpha* QBAR_tr**(alpha-mu) * NBAR_tr**(1.0_DP-alpha)
	        YBAR_tr = QBAR_tr ** alpha * NBAR_tr **(1.0_DP-alpha)
	        wage_tr = (1.0_DP-alpha)*QBAR_tr**alpha * NBAR_tr**(-alpha)
        endif 

        Ebar_tr = wage_tr  * NBAR_tr  * sum(pop)/sum(pop(1:RetAge-1))
        DBN_exp = DBN1 ; DBN_tr(:,:,:,:,:,:,T+1) = DBN1 ; 

        	print*, "Test Prices"
	        print*, "P   ", P_bench, P_tr(1), P_tr(T), P_tr(T+1)
	        print*, "wage", wage_bench, wage_tr(1), wage_tr(T), wage_tr(T+1)
	        print*, "EBAR", EBAR_bench, EBAR_tr(1), EBAR_tr(T), EBAR_tr(T+1)
	        print*, "Tau_K", tauK,"Tau_W_at",tauW_at,"Tau_W_bt",tauW_bt,"tau_L",1.0_dp-psi

        ! Save initial guess of prices
        OPEN (UNIT=76, FILE=trim(Result_Folder)//'Transition_Distance.txt', STATUS='replace')
        OPEN (UNIT=77, FILE=trim(Result_Folder)//'Transition_NBAR.txt', STATUS='replace')
        OPEN (UNIT=78, FILE=trim(Result_Folder)//'Transition_QBAR.txt', STATUS='replace')
        OPEN (UNIT=79, FILE=trim(Result_Folder)//'Transition_R.txt', STATUS='replace')
        OPEN (UNIT=80, FILE=trim(Result_Folder)//'Transition_GBAR.txt', STATUS='replace')
        OPEN (UNIT=81, FILE=trim(Result_Folder)//'Transition_GBAR_K.txt', STATUS='replace')
        OPEN (UNIT=82, FILE=trim(Result_Folder)//'Transition_GBAR_W.txt', STATUS='replace')
        OPEN (UNIT=83, FILE=trim(Result_Folder)//'Transition_GBAR_L.txt', STATUS='replace')
        OPEN (UNIT=84, FILE=trim(Result_Folder)//'Transition_GBAR_C.txt', STATUS='replace')
        OPEN (UNIT=85, FILE=trim(Result_Folder)//'Transition_SSC.txt', STATUS='replace')
        OPEN (UNIT=86, FILE=trim(Result_Folder)//'Transition_K.txt', STATUS='replace')
        OPEN (UNIT=87, FILE=trim(Result_Folder)//'Transition_C.txt', STATUS='replace')
        OPEN (UNIT=88, FILE=trim(Result_Folder)//'Transition_Top1.txt', STATUS='replace')
        OPEN (UNIT=89, FILE=trim(Result_Folder)//'Transition_Top10.txt', STATUS='replace')
        OPEN (UNIT=90, FILE=trim(Result_Folder)//'Transition_Debt.txt', STATUS='replace')
        OPEN (UNIT=91, FILE=trim(Result_Folder)//'Transition_Wage.txt', STATUS='replace')
     		
     		WRITE(UNIT=76, FMT=*) 'Iteration, DBN_dist, Q_dist, NW_dist, R_dist, Db_dist',&
     								&' Q(T)/Q(SS), N(T)/N(SS), W(T)/W(SS), R(T)/R(SS), Db(T)/Db(SS)'
     		WRITE(UNIT=77, FMT=*) 'NBAR'
        	WRITE(UNIT=78, FMT=*) 'QBAR'
        	WRITE(UNIT=79, FMT=*) 'R'
        	WRITE(UNIT=80, FMT=*) 'GBAR'
        	WRITE(UNIT=81, FMT=*) 'GBAR_K'
        	WRITE(UNIT=82, FMT=*) 'GBAR_W'
        	WRITE(UNIT=83, FMT=*) 'GBAR_L'
        	WRITE(UNIT=84, FMT=*) 'GBAR_C'
        	WRITE(UNIT=85, FMT=*) 'SSC'
        	WRITE(UNIT=86, FMT=*) 'K'
        	WRITE(UNIT=87, FMT=*) 'C'
        	WRITE(UNIT=88, FMT=*) 'Wealth_Top_1'
        	WRITE(UNIT=89, FMT=*) 'Wealth_Top_10'
        	WRITE(UNIT=90, FMT=*) 'Debt'
        	WRITE(UNIT=91, FMT=*) 'Wage'

        	WRITE(UNIT=77, FMT=*) NBAR_bench, NBAR_tr, NBAR_exp
        	WRITE(UNIT=78, FMT=*) QBAR_bench, QBAR_tr, QBAR_exp
        	WRITE(UNIT=79, FMT=*) R_bench, R_tr, R_exp
        	WRITE(UNIT=91, FMT=*) wage_bench, wage_tr, wage_exp


    	CLOSE (unit=76); CLOSE (unit=77); CLOSE (unit=78); CLOSE (unit=79);
    	CLOSE (unit=80); CLOSE (unit=81); CLOSE (unit=82); CLOSE (unit=83); CLOSE (unit=84); CLOSE (unit=85);
    	CLOSE (unit=86); CLOSE (unit=87); CLOSE (unit=88); CLOSE (unit=89); CLOSE (unit=90); CLOSE (unit=91);

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! DBN Iteration 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	print*,' '
	print*,'---------------------------------------------------'
	print*,' 	Starting Transition'
	print*,'---------------------------------------------------'
	! Compute distribution of assets by age and state
		! Distribution is obtained by iterating over an initial distribution using policy functions
	DBN_dist     = 1.0_DP
	Q_dist       = 1.0_DP
	NW_dist      = 1.0_DP
	R_dist       = 1.0_DP
	Db_dist      = 1.0_DP
	Chg_dist     = 1.0_DP 
	Old_DBN_dist = 0.0_DP 
	simutime     = 1
	iter_indx    = 1
	!print*, 'Computing Equilibrium Distribution'
	DO WHILE ((DBN_dist.ge.DBN_criteria).and.(max(Q_dist,NW_dist,R_dist,Db_dist).ge.Price_criteria)&
			& .and.(simutime.le.7).and.(Chg_dist.ge.Chg_criteria) )
		! print*, 'DBN_dist=', DBN_dist

		! Start Q_dist, N_dist, R_dist, Db_dist
		Q_dist = 0.0 ; NW_dist = 0.0 ; R_dist = 0.0 ; Db_dist = 0.0 ;

		! Start Debt to zero
		Debt_tr = 0.0_dp

		! Set initial value for aggregate variables 
			R = R_tr(T+1) ; P = P_tr(T+1) ; wage = wage_tr(T+1) ; DBN1 = DBN_tr(:,:,:,:,:,:,T+1) ; 
		! Deallocate grids that include thresholds, they are reset in Find_DBN_EQ
			deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )
		! Solve for New Steady State
			print*,' '; print*,'-----------------------------------------------------------------------------'
			print*,' 	Finding SS: Debt_SS=',Debt_SS,'R_0=',R,'Q_0=',QBAR_exp
			CALL FIND_DBN_EQ
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
				Cons_exp          = Cons           
				Hours_exp         = Hours
				Aprime_exp        = Aprime
				Debt_exp  		  = Debt_SS

				YBAR_C_exp = YBAR_C
				L_C_exp    = L_C
				K_C_exp    = K_C
				YBAR_P_exp = YBAR_P
				L_P_exp    = L_P
				K_P_exp    = K_P
		! Wealth and consumption in benchmark and experiment
			K_exp   = sum( sum(sum(sum(sum(sum(DBN_exp  ,6),5),4),3),1)*agrid )
			C_exp   = sum( DBN_exp  *Cons_exp   )
			print*,' 	New SS: Debt_SS=',Debt_SS,'R_SS=',R,'Q_SS=',QBAR
			print*,'-----------------------------------------------------------------------------';print*,' '



	    ! Solve for policy functions by backwards induction and EGM
	    	! Output is policy functions for all times and all ages
	    	CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
	    	CALL EGM_Transition
	    	print*, ' '

	    ! First Period Starts at DBN_bench
	    ! print*, ' Set DBN_tr for first period to benchmark distribution'
	    DBN_tr(:,:,:,:,:,:,1) = DBN_bench



		! Initialize DBN_tr to zero for other periods
		! print*, ' Initializing remaining periods of DBN_tr to zero'
			! DO ti=2,T+1
			! 	! print*,' 	Transition Period',ti
		 	! 		DBN_tr(:,:,:,:,:,:,ti)=0.0_DP
		 	! ENDDO
	 	DBN_tr(:,:,:,:,:,:,2:) = 0.0_DP

    	print*,' '
		print*,' 	--------------------------------------'
		print*,' 	Starting DBN Forward Iteration'
		print*,' 	--------------------------------------'
	    ! Fill in other periods starting at DBN bench following policy functions
	    DO ti=1,T
	    	! print*,' 	Transition Period ',ti

    	! Initialize DBN1 
    		DBN1 = 0.0_DP
    		call system_clock(tclock1)  ! start wall timer
    		call cpu_time(t1)   		! start cpu timer

		! Discretize policy function for assets (a') for current period
			! For each age and state vector bracket optimal a' between two grid points
			! When at that age and state the optimal decision is approximated by selecting one the grid points
			! The grid points are selected with probability proportional to their distance to the optimal a'
		! print*,' Discretizing Policy Functions'
		DO age=1,MaxAge
		!$omp parallel do private(lambdai,ei,ai,xi,tklo,tkhi)
		DO zi=1,nz
		!$omp parallel do private(lambdai,ei,ai,tklo,tkhi)
		DO xi=1,nx
		DO ai=1,na
		DO lambdai=1,nlambda
		DO ei=1, ne
	        if ( Aprime_tr(age,ai,zi,lambdai,ei,xi,ti) .ge. amax) then
	            tklo =na-1
	        elseif (Aprime_tr(age,ai,zi,lambdai,ei,xi,ti) .lt. amin) then
	            tklo = 1
	        else
	            tklo = ((Aprime_tr(age,ai,zi,lambdai,ei,xi,ti) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
	        endif
	        tkhi = tklo + 1        
	        Aplo(age,ai,zi,lambdai,ei,xi)  	   = tklo
	        Aphi(age,ai,zi,lambdai,ei,xi)  	   = tkhi        
	        PrAprimelo(age,ai,zi,lambdai,ei,xi) = (agrid(tkhi) - Aprime_tr(age,ai,zi,lambdai,ei,xi,ti))/(agrid(tkhi)-agrid(tklo))
	        PrAprimehi(age,ai,zi,lambdai,ei,xi) = (Aprime_tr(age,ai,zi,lambdai,ei,xi,ti) - agrid(tklo))/(agrid(tkhi)-agrid(tklo))

	        if ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime_tr(age,ai,zi,lambdai,ei,xi,ti) .ge. amax) then
	            tklo =na-1
	        elseif ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime_tr(age,ai,zi,lambdai,ei,xi,ti) .lt. amin) then
	            tklo = 1
	        else
	            tklo = (((1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime_tr(age,ai,zi,lambdai,ei,xi,ti)-amin)/(amax-amin)) & 
	            			& **(1.0_DP/a_theta)*(na-1)+1        
	        endif
	        tkhi = tklo + 1        
	        Bqlo(age,ai,zi,lambdai,ei,xi)  	= tklo
	        Bqhi(age,ai,zi,lambdai,ei,xi)  	= tkhi        
	        PrBqlo(age,ai,zi,lambdai,ei,xi) = ( agrid(tkhi) - (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime_tr(age,ai,zi,lambdai,ei,xi,ti) )&
	        									& / ( agrid(tkhi) -agrid(tklo) )
	        PrBqhi(age,ai,zi,lambdai,ei,xi) = ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Aprime_tr(age,ai,zi,lambdai,ei,xi,ti) - agrid(tklo) )&
	        									& / ( agrid(tkhi) -agrid(tklo) )
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		! print*, ' Discretizing Policy Functions Completed'

		! Probablities are adjusted to lie in [0,1]
		PrAprimelo = min(PrAprimelo, 1.0_DP); PrAprimelo = max(PrAprimelo, 0.0_DP)
		PrAprimehi = min(PrAprimehi, 1.0_DP); PrAprimehi = max(PrAprimehi, 0.0_DP)

		PrBqlo = min(PrBqlo, 1.0_DP); PrBqlo = max(PrBqlo, 0.0_DP)
		PrBqhi = min(PrBqhi, 1.0_DP); PrBqhi = max(PrBqhi, 0.0_DP)


		! Everyone in MaxAge dies. Those who die, switch to z2, lambda2 and start at ne/2+1 and x=1
	    age1=MaxAge
	    ! $omp parallel do reduction(+:DBN1) private(x1,a1,lambda1,e1,z2,lambda2)
	    DO z1=1,nz
	    DO x1=1,nx
	    DO a1=1,na
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN1(1,Bqlo(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1 ,1)   =  &
	           		& DBN1(1,Bqlo(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1) + DBN_tr(age1,a1,z1,lambda1,e1,x1,ti) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN1(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1 ,1)   =  &
	           		& DBN1(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1) + DBN_tr(age1,a1,z1,lambda1,e1,x1,ti) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO    
	    ENDDO
	    ENDDO  
	    !$omp barrier  

		! retirees "e" stays the same for benefit retirement calculation purposes
		! $omp parallel do reduction(+:DBN1) private(x1,age1,a1,lambda1,e1,z2,lambda2,x2)
	    DO z1=1,nz
	    DO x1=1,nx
	    DO age1=RetAge-1, MaxAge-1
	    DO a1=1,na
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        ! Those who die, switch to z2, lambda2 and start at ne/2+1
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN1(1,Bqlo(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1 ,1)   =  &
	           		& DBN1(1,Bqlo(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1) + DBN_tr(age1,a1,z1,lambda1,e1,x1,ti) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN1(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1 ,1)   =  &
	           		& DBN1(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1) + DBN_tr(age1,a1,z1,lambda1,e1,x1,ti) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, and also e1 since they are retired
	        !e2=e1
	        DO x2=1,nx
		        DBN1(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2 ) =  &
		        	& DBN1(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2 ) + &
		        	& DBN_tr(age1, a1, z1, lambda1, e1,x1,ti) &
		            & * survP(age1) * PrAprimelo(age1,a1,z1,lambda1,e1,x1) * pr_x(x1,x2,z1,age1)     
		        DBN1(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2 ) =  &
		          	& DBN1(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2 ) + &
		          	& DBN_tr(age1, a1, z1, lambda1, e1,x1,ti) &
		            & * survP(age1) * PrAprimehi(age1,a1,z1,lambda1,e1,x1) * pr_x(x1,x2,z1,age1)  
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    !$omp barrier
	    
	    ! Working age agents
	    ! $omp parallel do reduction(+:DBN1) private(x1,age1,a1,lambda1,e1,z2,lambda2,x2,e2)
	    DO z1=1,nz
	    DO x1=1,nx
	    DO age1=1,RetAge-2
	    DO a1=1,na
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        ! Those who die, switch to z2, lambda2 and start at ne/2+1
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN1(1,Bqlo(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1 ,1)   =  &
	           		& DBN1(1,Bqlo(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1) + DBN_tr(age1,a1,z1,lambda1,e1,x1,ti) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN1(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1 ,1)   =  &
	           		& DBN1(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1) + DBN_tr(age1,a1,z1,lambda1,e1,x1,ti) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, but switch to e2 and x2
	        DO x2=1,nx
	        DO e2=1,ne
				DBN1(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2 ) =  &
	          		& DBN1(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2 ) + &
	          		& DBN_tr(age1, a1, z1, lambda1, e1,x1,ti) &
	                & * survP(age1) * pr_e(e1,e2) * pr_x(x1,x2,z1,age1) * PrAprimelo(age1,a1,z1,lambda1,e1,x1)
	            DBN1(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2 ) =  &
	          		& DBN1(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2 ) + &
	          		& DBN_tr(age1, a1, z1, lambda1, e1,x1,ti) &
	                & * survP(age1) * pr_e(e1,e2) * pr_x(x1,x2,z1,age1) * PrAprimehi(age1,a1,z1,lambda1,e1,x1) 
	        ENDDO
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    !$omp barrier

	    call cpu_time(t2)   ! end cpu timer
    	call system_clock(tclock2, clock_rate); elapsed_time = float(tclock2 - tclock1) / float(clock_rate)

	    ! Update distribution
	    	DBN_tr(:,:,:,:,:,:,ti+1) = DBN1

	    ! Set global variables to current period  (Time: ti)
	    	R = R_tr(ti) ; P = P_tr(ti) ; wage = wage_tr(ti) ;
	    	K_mat  = K_Matrix(R,P)
			Pr_mat = Profit_Matrix(R,P)
			CALL ComputeLaborUnits(Ebar_tr(ti), wage_tr(ti)) 
			CALL FORM_Y_MB_GRID(YGRID,MBGRID,YGRID_t,MBGRID_t)
			DBN1  = DBN_tr(:,:,:,:,:,:,ti)
			Cons  = Cons_tr(:,:,:,:,:,:,ti)
			Hours = Hours_tr(:,:,:,:,:,:,ti)
			    


	    ! Compute prices and aggregates for the current period (Time: ti)
	    ! print*,' Updating Prices and Quantities'
    		! Compute aggregates with current distribution (Time: ti)
	        QBAR2_tr(ti) =0.0
	        NBAR2_tr(ti) =0.0
	        DO x1=1,nx
	        DO age1=1,MaxAge
	        DO z1=1,nz
	        DO a1=1,na
	        DO lambda1=1,nlambda
	        DO e1=1, ne
	             QBAR2_tr(ti)= QBAR2_tr(ti)+ DBN_tr(age1,a1,z1,lambda1,e1,x1,ti) * ( xz_grid(x1,z1) * K_mat(a1,z1,x1) )**mu
	             NBAR2_tr(ti)= NBAR2_tr(ti)+ &
	             			& DBN_tr(age1,a1,z1,lambda1,e1,x1,ti) * eff_un(age1,lambda1,e1) * Hours_tr(age1,a1,z1,lambda1,e1,x1,ti)
	        ENDDO
	        ENDDO
	        ENDDO
	        ENDDO    
	        ENDDO    
	        ENDDO    
	        QBAR2_tr(ti) = ( QBAR2_tr(ti))**(1.0_DP/mu) 

	        ! Get Q_dist and NW_dist before dampening 
	        	Q_dist  	 = max( Q_dist,abs(QBAR2_tr(ti)/QBAR_tr(ti)-1))
	        	QBAR_tr(ti)  = Dampen*QBAR_tr(ti) + (1.0_dp-Dampen)*QBAR2_tr(ti)
	        	if (A_C.gt.0.0_dp) then
	        	NBAR_tr(ti)  = NBAR2_tr(ti)
	        	else
	        	NW_dist 	 = max(NW_dist,abs(NBAR2_tr(ti)/NBAR_tr(ti)-1))
	        	NBAR_tr(ti)  = Dampen*NBAR_tr(ti) + (1.0_dp-Dampen)*NBAR2_tr(ti)
	        	endif 

        	! Compute total assets and consumption
		        K_tr(ti) = sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid )
		        C_tr(ti) = sum( DBN1*Cons )

			! Update prices 		        
	        if (A_C.gt.0.0_dp) then 
	        	! Capital Market: Private demand 
	        	K_P_tr(ti)    = sum( (sum(sum(sum(DBN1,5),4),1)) *(K_mat)) ! Note: DBN_azx  = sum(sum(sum(DBN1,5),4),1)
	        	
	        	if (K_tr(ti).gt.(K_P_tr(ti)+Debt_Absorption*Debt_tr(ti-1))) then 
		        	! Corporate demand for capital
		        	K_C_tr(ti)   = K_tr(ti)- K_P_tr(ti)-Debt_Absorption*Debt_tr(ti-1)

		        	! Update Wage 
		        	if (alpha_C.eq.alpha) then 
		        		Wage2_tr(ti) = ((((1.0_dp-alpha)*Aprod)**(1.0_dp/alpha)*QBAR_tr(ti) + &
		        				& ((1.0_dp-alpha)*A_C)**(1.0_dp/alpha)*K_C_tr(ti))/NBAR_tr(ti))**(alpha)
		        	else 
		        		print*, ' Solving labor market numerically - Transition Period:', ti
		        		QBAR = QBAR_tr(ti); NBAR = NBAR_tr(ti); K_C = K_C_tr(ti); 
		        		brent_value = brent(0.001_DP,Wage_tr(ti),10.0_DP,Labor_Market_Clearing, brent_tol,Wage2_tr(ti))
		        		print*, ' New Wage=',Wage,'Error=',brent_value
		        	endif 
		        	NW_dist = max(NW_dist,abs(wage2_tr(ti)/wage_tr(ti)-1))
		        	wage_tr(ti)  = Dampen*wage_tr(ti) + (1.0_dp-Dampen)*wage2_tr(ti)

	        	else ! Capital Market did not clear

	        		print*, ' ';
	        		print*, ' 	Warning! Capital Market Did not Clear! Transition Period:',ti;
	        		print*, ' ';
		        	! Modify prices and quantities 
		        	QBAR_tr(ti) = 0.50_dp*QBAR_tr(ti) 
		        	Wage_tr(ti) = Wage_tr(ti) 

	        	endif 

	        	! Update other prices and quantities
        		L_P_tr(ti)  = ( (1.0_dp-alpha)*Aprod/Wage_tr(ti) )**(1.0_dp/alpha  ) * QBAR_tr(ti) 
        		L_C_tr(ti)  = ( (1.0_dp-alpha_C)*A_C/Wage_tr(ti) )**(1.0_dp/alpha_C) * K_C_tr(ti)
				
				P_tr(ti)    = alpha * QBAR_tr(ti)**(alpha-mu) * L_P_tr(ti)**(1.0_DP-alpha)
				R2_tr(ti)   = alpha_C * A_C * ( Wage_tr(ti)/((1.0_dp-alpha_C)*A_C) )**(-(1.0_dp-alpha_C)/alpha_C) - DepRate
					R_dist  = max(R_dist,abs(R2_tr(ti)/R_tr(ti)-1))
					R_tr(ti)= R2_tr(ti) ! Dampen*R_old + (1.0_dp-Dampen)*R2_tr(ti)
				
				Ebar_tr(ti)   = wage_tr(ti)  * NBAR_tr(ti)  * sum(pop)/sum(pop(1:RetAge-1))
				YBAR_P_tr(ti) = AProd * QBAR_tr(ti)**alpha   * L_P_tr(ti)**(1.0_DP-alpha  ) 
				YBAR_C_tr(ti) = A_C   * K_C_tr(ti) **alpha_C * L_C_tr(ti)**(1.0_DP-alpha_C) 
				YBAR_tr(ti)   = YBAR_P_tr(ti) + YBAR_C_tr(ti)
				! print '(A,F9.3,F9.3,F9.3,X,X,A,F9.3,F9.3,F9.3)',&
				! 		& ' 				Corp. Sector Levels:', YBAR_C_tr(ti), K_C_tr(ti), L_C_tr(ti) , &
				! 		& ' Ratios ', 100.0_dp*YBAR_C_tr(ti)/YBAR_tr(ti), 100.0_dp*K_C_tr(ti)/K_tr(ti), 100.0_dp*L_C_tr(ti)/NBAR_tr(ti)

	        else 
	        	! Set Corporate values to zero
        		K_C_tr(ti) = 0.0_dp; L_C_tr(ti) = 0.0_dp; YBAR_C_tr(ti) = 0.0_dp
        		K_P_tr(ti) = 0.0_dp; L_P_tr(ti) = 0.0_dp; YBAR_P_tr(ti) = 0.0_dp

	        	! Update other prices and quantities             
		        P_tr(ti)     = alpha* QBAR_tr(ti)**(alpha-mu) * NBAR_tr(ti)**(1.0_DP-alpha)
		        YBAR_tr(ti)  = QBAR_tr(ti)**alpha * NBAR_tr(ti)**(1.0_DP-alpha)
		        wage_tr(ti)  = (1.0_DP-alpha)*QBAR_tr(ti)**alpha * NBAR_tr(ti)**(-alpha)
		        Ebar_tr(ti)  = wage_tr(ti)  * NBAR_tr(ti) * sum(pop)/sum(pop(1:RetAge-1))

		    	! Solve for new R (that clears market under new guess for prices)
		    	! Update only every third period or in the first and last periods 
		    	if ((ind_R.eq.3)) then ! .or.(ti.eq.T)
			    		! Save old R for updating 
			    		R_old = R_tr(ti)
			    	if (sum(theta)/nz .gt. 1.0_DP) then
			    		! Set price 
			    		P = min(P_tr(ti),1.0_dp)
			    		! Solve for R using brent
			    			! call system_clock(tclock1_R)  ! start wall timer
	    					! call cpu_time(t1_R)   		! start cpu timer
			            brent_value = brent(-0.05_DP,R_old,0.15_DP,Agg_Debt_Tr,0.000001_DP,R2_tr(ti))
			            	! Usually brent_tol=0.00000001_DP
			             	! call cpu_time(t2_R)   ! end cpu timer
	    					! call system_clock(tclock2_R, clock_rate_R); elapsed_time_R = float(tclock2_R - tclock1_R) / float(clock_rate_R)
			            	print*, ' 	Solving for equilibrium interest rate (R)  -  Error=',brent_value,&
			            		& 'R_out=',R2_tr(ti),'Debt_Absorption=',Debt_Absorption
			            	! print*, '	Brent Time:		CPU time:',t2_R-t1_R,'sec		Elapsed time:',elapsed_time_R,'sec'
		            		print*, ' '
		            else
		                R2_tr(ti) = 0.0_DP
			        endif

			        ! Get R_dist before dampening 
		        	R_dist = max(R_dist,abs(R2_tr(ti)/R_tr(ti)-1))

		        	! Dampened Update of R
		        	R_tr(ti)  = Dampen*R_old + (1.0_dp-Dampen)*R2_tr(ti)

		        	if ((ti.gt.1).and.(ti.lt.T)) then 
		        		! Update extrapolation by interpolating between last update and this update
	    				R_tr(ti-2) = real(2,8)/3.0_dp*R_tr(ti-3)+real(1,8)/3.0_dp*R_tr(ti)
	    				R_tr(ti-1) = real(1,8)/3.0_dp*R_tr(ti-3)+real(2,8)/3.0_dp*R_tr(ti)
		        		! print*, '	R Interpolation'
		        		! print*, ' 	R1=',R_tr(ti-3),'R2=',R_tr(ti-2),'R3=',R_tr(ti-1),'R4=',R_tr(ti)
		        		! print*, ' '
	 

		        		! Update Slope for extrapolation 
		        		R_slope = (R_tr(ti)-R_tr(ti-3))/3.0_dp 

	        		!elseif (ti.eq.T) then 

	        		elseif (ti.eq.1) then 

	        			R_slope = 0.0_dp

		        	endif 

		        	! Update index
	        		ind_R = 1 

		        else 
		        	! Update by extrapolating from last update 
		    		R_tr(ti) = R_slope*ind_R+R_tr(ti-ind_R)
	        		! print*, '	R Interpolation'
	        		! print*, '		p1=',real(3-ind_R,8)/3.0_dp,'p2=',real(ind_R,8)/3.0_dp,'R1=',R_tr(ti-ind_R),'R2=',R_tr(ti+3-ind_R)
	        		! print*, ' '

		        	! Update index
		        	ind_R = ind_R+1
		        endif 


	        endif 


		    ! Compute government budget for the current preiod (Time: ti)
		    ! print*,' Calculating tax revenue'
	    	CALL GOVNT_BUDGET(.false.)
			    GBAR_tr(ti) 		= GBAR 
			    GBAR_K_tr(ti) 		= GBAR_K 
			    GBAR_W_tr(ti) 		= GBAR_W
			    GBAR_L_tr(ti) 		= GBAR_L 
			    GBAR_C_tr(ti) 		= GBAR_C 
			    SSC_Payments_tr(ti) = SSC_Payments
			    Tot_Lab_Inc_tr(ti) 	= Tot_Lab_Inc
			    if (ti.eq.1) then 
			    	Debt_tr(ti)     = (GBAR_bench-GBAR_tr(ti)) ! First period debt is equal to  the deficit
			    else
			    	Debt_tr(ti)     = Debt_tr(ti-1)*(1+R_tr(ti)) + (GBAR_bench-GBAR_tr(ti)) ! Debt is equal to  the deficit plus (compunded) previous debt
			    endif 



	        ! Compute top wealth concentration
	    		DO ai=1,na
					pr_a_dbn(ai)          = sum(DBN1(:,ai,:,:,:,:)) 
					cdf_a_dbn(ai)         = sum( pr_a_dbn(1:ai) )      
					tot_a_by_grid(ai)     = sum(DBN1(:,ai,:,:,:,:) * agrid(ai) )
					cdf_tot_a_by_grid(ai) = sum(tot_a_by_grid(1:ai))   
				ENDDO
				cdf_a_dbn = cdf_a_dbn + 1.0_DP - cdf_a_dbn(na)

		
				! FIND THE ai THAT CORRESPONDS TO EACH PRCTILE OF WEALTH DBN & WEALTH HELD BY PEOPLE LOWER THAN THAT PRCTILE
				DO prctile=90,100
				    ai=1
				    DO while (cdf_a_dbn(ai) .lt. (REAL(prctile,8)/100.0_DP-0.000000000000001))
				        ai=ai+1
				    ENDDO
				    prctile_ai_ind(prctile) = ai
				    prctile_ai(prctile)     = agrid(ai)

				    IF (ai .gt. 1) THEN
				        cdf_tot_a_by_prctile(prctile)  = cdf_tot_a_by_grid(ai-1) + (REAL(prctile,8)/100.0_DP - cdf_a_dbn(ai-1))*agrid(ai) 
				    else
				        cdf_tot_a_by_prctile(prctile)  = (REAL(prctile,8)/100.0_DP )*agrid(ai)     
				    ENDIF
				ENDDO
				! Store top wealth shares 
			        Wealth_Top_1_Tr(ti)  = 1.0_DP-cdf_tot_a_by_prctile(99)/cdf_tot_a_by_prctile(100)
			        Wealth_Top_10_Tr(ti) = 1.0_DP-cdf_tot_a_by_prctile(90)/cdf_tot_a_by_prctile(100)

	        ! print*, '	Prices and Dampening'
	        ! print*, '		Q_full=',QBAR2_tr(ti),'Q_tr=',QBAR_tr(ti)
	        ! print*, '		N_full=',NBAR2_tr(ti),'N_tr=',NBAR_tr(ti)
	        ! print*, '		R_full=',R2_tr(ti),'R_tr=',R_tr(ti)
	        ! print*, ' '

	        !!
	        print 12345, 't=',ti,'Deficit=',(GBAR_bench-GBAR_tr(ti)),'Debt=',Debt_tr(ti),&
	        	& 'Wealth=',K_tr(ti),'R=',100_dp*R_tr(ti),'Q=',QBAR_tr(ti),'W=',Wage_tr(ti),&
	        	& 'K_C/A=',100.0_dp*K_C_tr(ti)/K_tr(ti),'L_C/N=',100.0_dp*L_C_tr(ti)/NBAR_tr(ti)	    		
    		12345 format &
    		&(A,I3,X,X,A,F8.5,X,X,A,F8.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3)
	    	!!
	        ! print*, '		CPU time:',t2-t1,'sec		Elapsed time:',elapsed_time,'sec'
	        print*, ' '

	        if (Debt_Absorption*Debt_tr(ti).gt.K_tr(ti)) then 
	        	print*,' '
	        	print*,' 	--------------------------------------'
				print*,' 	Error in transition - Taxes are too low to sustain debt'
				print*,' 	Exiting Transition Loop'
				print*,' 	--------------------------------------'

				! Set GBAR and Debt to increase tax  
				GBAR_exp = 0; Debt_tr(T+1) = Debt_tr(ti)
				! Set flag to avoid using current solution as seed
				Use_Transition_Seed = .false.
				return 
			endif 
    	ENDDO ! Transition Time

		print*,' 	--------------------------------------'
		print*,' 	DBN Forward Iteration Completed'
		print*,' 	--------------------------------------'
		print*,' '


	    ! Set global variables to current period  (Time: ti)
	    	R = R_tr(T+1) ; P = P_tr(T+1) ; wage = wage_tr(T+1) ;
	    	K_mat  = K_Matrix(R,P)
			Pr_mat = Profit_Matrix(R,P)
			CALL ComputeLaborUnits(Ebar_tr(T+1), wage_tr(T+1)) 
			CALL FORM_Y_MB_GRID(YGRID,MBGRID,YGRID_t,MBGRID_t)
			DBN1  = DBN_tr(:,:,:,:,:,:,T+1)
			Cons  = Cons_tr(:,:,:,:,:,:,T+1)
			Hours = Hours_tr(:,:,:,:,:,:,T+1)

	    
	    ! Compute prices and aggregates for the current period (Time: T+1)
	    ! print*,' Updating Prices and Quantities'
    		! Compute aggregates with current distribution (Time: T+1)
	        QBAR2_tr(T+1) =0.0
	        NBAR2_tr(T+1) =0.0
	        DO x1=1,nx
	        DO age1=1,MaxAge
	        DO z1=1,nz
	        DO a1=1,na
	        DO lambda1=1,nlambda
	        DO e1=1, ne
	             QBAR2_tr(T+1)= QBAR2_tr(T+1)+ DBN_tr(age1,a1,z1,lambda1,e1,x1,T+1) * ( xz_grid(x1,z1) * K_mat(a1,z1,x1) )**mu
	             NBAR2_tr(T+1)= NBAR2_tr(T+1)+ &
	             			& DBN_tr(age1,a1,z1,lambda1,e1,x1,T+1) * eff_un(age1,lambda1,e1) * Hours_tr(age1,a1,z1,lambda1,e1,x1,T+1)
	        ENDDO
	        ENDDO
	        ENDDO
	        ENDDO    
	        ENDDO    
	        ENDDO    
	        QBAR2_tr(T+1) = ( QBAR2_tr(T+1))**(1.0_DP/mu) 

	        ! Get Q_dist and NW_dist before dampening 
	        	Q_dist  	 = max( Q_dist,abs(QBAR2_tr(T+1)/QBAR_tr(T+1)-1))
	        	QBAR_tr(T+1)  = Dampen*QBAR_tr(T+1) + (1.0_dp-Dampen)*QBAR2_tr(T+1)
	        	if (A_C.gt.0.0_dp) then
	        	NBAR_tr(T+1)  = NBAR2_tr(T+1)
	        	else
	        	NW_dist 	 = max(NW_dist,abs(NBAR2_tr(T+1)/NBAR_tr(T+1)-1))
	        	NBAR_tr(T+1)  = Dampen*NBAR_tr(T+1) + (1.0_dp-Dampen)*NBAR2_tr(T+1)
	        	endif 

        	! Compute total assets and consumption
		        K_tr(T+1) = sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid )
		        C_tr(T+1) = sum( DBN1*Cons )

			! Update prices 		        
	        if (A_C.gt.0.0_dp) then 
	        	! Capital Market: Private demand 
	        	K_P_tr(T+1)    = sum( (sum(sum(sum(DBN1,5),4),1)) *(K_mat)) ! Note: DBN_azx  = sum(sum(sum(DBN1,5),4),1)
	        	
	        	if (K_tr(T+1).gt.(K_P_tr(T+1)+Debt_Absorption*Debt_tr(T))) then 
		        	! Corporate demand for capital
		        	K_C_tr(T+1)   = K_tr(T+1)- K_P_tr(T+1)-Debt_Absorption*Debt_tr(T)

		        	! Update Wage 
		        	if (alpha_C.eq.alpha) then 
		        		Wage2_tr(T+1) = ((((1.0_dp-alpha)*Aprod)**(1.0_dp/alpha)*QBAR_tr(T+1) + &
		        				& ((1.0_dp-alpha)*A_C)**(1.0_dp/alpha)*K_C_tr(T+1))/NBAR_tr(T+1))**(alpha)
		        	else 
		        		print*, ' Solving labor market numerically - Transition Period:', ti
		        		QBAR = QBAR_tr(T+1); NBAR = NBAR_tr(T+1); K_C = K_C_tr(T+1); 
		        		brent_value = brent(0.001_DP,Wage_tr(T+1),10.0_DP,Labor_Market_Clearing, brent_tol,Wage2_tr(T+1))
		        		print*, ' New Wage=',Wage,'Error=',brent_value
		        	endif 
		        	NW_dist = max(NW_dist,abs(wage2_tr(T+1)/wage_tr(T+1)-1))
		        	wage_tr(T+1)  = Dampen*wage_tr(T+1) + (1.0_dp-Dampen)*wage2_tr(T+1)

	        	else ! Capital Market did not clear

	        		print*, ' ';
	        		print*, ' 	Warning! Capital Market Did not Clear! Transition Period:',ti;
	        		print*, ' ';
		        	! Modify prices and quantities 
		        	QBAR_tr(T+1) = 0.50_dp*QBAR_tr(T+1) 
		        	Wage_tr(T+1) = Wage_tr(T+1) 

	        	endif 

	        	! Update other prices and quantities
        		L_P_tr(T+1)  = ( (1.0_dp-alpha)*Aprod/Wage_tr(T+1) )**(1.0_dp/alpha  ) * QBAR_tr(T+1) 
        		L_C_tr(T+1)  = ( (1.0_dp-alpha_C)*A_C/Wage_tr(T+1) )**(1.0_dp/alpha_C) * K_C_tr(T+1)
				
				P_tr(T+1)    = alpha * QBAR_tr(T+1)**(alpha-mu) * L_P_tr(T+1)**(1.0_DP-alpha)
				R2_tr(T+1)   = alpha_C * A_C * ( Wage_tr(T+1)/((1.0_dp-alpha_C)*A_C) )**(-(1.0_dp-alpha_C)/alpha_C) - DepRate
					R_dist  = max(R_dist,abs(R2_tr(T+1)/R_tr(T+1)-1))
					R_tr(T+1)= R2_tr(T+1) ! Dampen*R_old + (1.0_dp-Dampen)*R2_tr(T+1)
				
				Ebar_tr(T+1)   = wage_tr(T+1)  * NBAR_tr(T+1)  * sum(pop)/sum(pop(1:RetAge-1))
				YBAR_P_tr(T+1) = AProd * QBAR_tr(T+1)**alpha   * L_P_tr(T+1)**(1.0_DP-alpha  ) 
				YBAR_C_tr(T+1) = A_C   * K_C_tr(T+1) **alpha_C * L_C_tr(T+1)**(1.0_DP-alpha_C) 
				YBAR_tr(T+1)   = YBAR_P_tr(T+1) + YBAR_C_tr(T+1)
				! print '(A,F9.3,F9.3,F9.3,X,X,A,F9.3,F9.3,F9.3)',&
				! 		& ' 				Corp. Sector Levels:', YBAR_C_tr(T+1), K_C_tr(T+1), L_C_tr(T+1) , &
				! 		& ' Ratios ', 100.0_dp*YBAR_C_tr(T+1)/YBAR_tr(T+1), 100.0_dp*K_C_tr(T+1)/K_tr(T+1), 100.0_dp*L_C_tr(T+1)/NBAR_tr(T+1)

	        else 
	        	! Set Corporate values to zero
        		K_C_tr(T+1) = 0.0_dp; L_C_tr(T+1) = 0.0_dp; YBAR_C_tr(T+1) = 0.0_dp
        		K_P_tr(T+1) = 0.0_dp; L_P_tr(T+1) = 0.0_dp; YBAR_P_tr(T+1) = 0.0_dp

	        	! Update other prices and quantities             
		        P_tr(T+1)     = alpha* QBAR_tr(T+1)**(alpha-mu) * NBAR_tr(T+1)**(1.0_DP-alpha)
		        YBAR_tr(T+1)  = QBAR_tr(T+1)**alpha * NBAR_tr(T+1)**(1.0_DP-alpha)
		        wage_tr(T+1)  = (1.0_DP-alpha)*QBAR_tr(T+1)**alpha * NBAR_tr(T+1)**(-alpha)
		        Ebar_tr(T+1)  = wage_tr(T+1)  * NBAR_tr(T+1) * sum(pop)/sum(pop(1:RetAge-1))

		    	! Solve for new R (that clears market under new guess for prices)
    				! Save old R for updating 
		    		R_old = R_tr(T+1)
		    	if (sum(theta)/nz .gt. 1.0_DP) then
		    		P = min(P_tr(T+1),1.0_dp)
		            brent_value = brent(-0.05_DP,R_old,0.15_DP,Agg_Debt_Tr,0.000001_DP,R2_tr(T+1))
		            	! Usually brent_tol=0.00000001_DP
		            	print*, ' 	Solving for equilibrium interest rate (R)  -  Error=',brent_value,&
		            		& 'R_out=',R2_tr(T+1),'Debt_Absorption=',Debt_Absorption
	            		print*, ' '
	            else
	                R2_tr(T+1) = 0.0_DP
		        endif

		        ! Get R_dist before dampening 
	        	R_dist = max(R_dist,abs(R2_tr(T+1)/R_tr(T+1)-1))

	        	! Dampened Update of R
	        	R_tr(T+1)  = Dampen*R_old + (1.0_dp-Dampen)*R2_tr(T+1)

	        endif 

	        ! Compute government budget for the current preiod (Time: T+1)
	    	! print*,' Calculating tax revenue'
	    	CALL GOVNT_BUDGET(.false.)
			    GBAR_tr(T+1) 		 = GBAR 
			    GBAR_K_tr(T+1) 		 = GBAR_K 
			    GBAR_W_tr(T+1) 		 = GBAR_W
			    GBAR_L_tr(T+1) 		 = GBAR_L 
			    GBAR_C_tr(T+1) 		 = GBAR_C 
			    SSC_Payments_tr(T+1) = SSC_Payments
			    Tot_Lab_Inc_tr(T+1)  = Tot_Lab_Inc
			    Debt_tr(T+1) 		 = Debt_tr(T)           
			    ! Note that Debt_tr is not computed for T+1 since that would increase interest payments

		    ! Get Db_dist before dampening 
	        	Db_dist = max(Db_dist,abs(Debt_SS/Debt_tr(T+1)-1))

	        ! Dampened Update of Db_ss
	        	Debt_SS  = 0.8*Debt_SS + 0.2*Debt_tr(T+1)	


	        print*,' '
	        print 12345, 't=',T+1,'Deficit=',(GBAR_bench-GBAR_tr(T+1)),'Debt=',Debt_tr(T+1),&
	        	& 'Wealth=',K_tr(T+1),'R=',100_dp*R_tr(T+1),'Q=',QBAR_tr(T+1),'W=',Wage_tr(T+1),&
	        	& 'K_C/A=',100.0_dp*K_C_tr(T+1)/K_tr(T+1),'L_C/N=',100.0_dp*L_C_tr(T+1)/NBAR_tr(T+1)	    		
	        print*,' '

	        ! Compute top wealth concentration
	    		DO ai=1,na
					pr_a_dbn(ai)          = sum(DBN1(:,ai,:,:,:,:)) 
					cdf_a_dbn(ai)         = sum( pr_a_dbn(1:ai) )      
					tot_a_by_grid(ai)     = sum(DBN1(:,ai,:,:,:,:) * agrid(ai) )
					cdf_tot_a_by_grid(ai) = sum(tot_a_by_grid(1:ai))   
				ENDDO
				cdf_a_dbn = cdf_a_dbn + 1.0_DP - cdf_a_dbn(na)

		
				! FIND THE ai THAT CORRESPONDS TO EACH PRCTILE OF WEALTH DBN & WEALTH HELD BY PEOPLE LOWER THAN THAT PRCTILE
				DO prctile=90,100
				    ai=1
				    DO while (cdf_a_dbn(ai) .lt. (REAL(prctile,8)/100.0_DP-0.000000000000001))
				        ai=ai+1
				    ENDDO
				    prctile_ai_ind(prctile) = ai
				    prctile_ai(prctile)     = agrid(ai)

				    IF (ai .gt. 1) THEN
				        cdf_tot_a_by_prctile(prctile)  = cdf_tot_a_by_grid(ai-1) + (REAL(prctile,8)/100.0_DP - cdf_a_dbn(ai-1))*agrid(ai) 
				    else
				        cdf_tot_a_by_prctile(prctile)  = (REAL(prctile,8)/100.0_DP )*agrid(ai)     
				    ENDIF
				ENDDO
				! Store top wealth shares 
			        Wealth_Top_1_Tr(T+1)  = 1.0_DP-cdf_tot_a_by_prctile(99)/cdf_tot_a_by_prctile(100)
			        Wealth_Top_10_Tr(T+1) = 1.0_DP-cdf_tot_a_by_prctile(90)/cdf_tot_a_by_prctile(100)
	    

	    	! Save Prices
	    	OPEN (UNIT=77, FILE=trim(Result_Folder)//'Transition_NBAR.txt', STATUS='old', POSITION='append')
	    	OPEN (UNIT=78, FILE=trim(Result_Folder)//'Transition_QBAR.txt', STATUS='old', POSITION='append')
	    	OPEN (UNIT=79, FILE=trim(Result_Folder)//'Transition_R.txt', STATUS='old', POSITION='append')
	    	OPEN (UNIT=80, FILE=trim(Result_Folder)//'Transition_GBAR.txt', STATUS='old', POSITION='append')
	    	OPEN (UNIT=81, FILE=trim(Result_Folder)//'Transition_GBAR_K.txt', STATUS='old', POSITION='append')
	    	OPEN (UNIT=82, FILE=trim(Result_Folder)//'Transition_GBAR_W.txt', STATUS='old', POSITION='append')
	    	OPEN (UNIT=83, FILE=trim(Result_Folder)//'Transition_GBAR_L.txt', STATUS='old', POSITION='append')
	    	OPEN (UNIT=84, FILE=trim(Result_Folder)//'Transition_GBAR_C.txt', STATUS='old', POSITION='append')
	    	OPEN (UNIT=85, FILE=trim(Result_Folder)//'Transition_SSC.txt', STATUS='old', POSITION='append')
	    	OPEN (UNIT=86, FILE=trim(Result_Folder)//'Transition_K.txt', STATUS='old', POSITION='append')
	    	OPEN (UNIT=87, FILE=trim(Result_Folder)//'Transition_C.txt', STATUS='old', POSITION='append')
	    	OPEN (UNIT=88, FILE=trim(Result_Folder)//'Transition_Top1.txt', STATUS='old', POSITION='append')
	    	OPEN (UNIT=89, FILE=trim(Result_Folder)//'Transition_Top10.txt', STATUS='old', POSITION='append')
	    	OPEN (UNIT=90, FILE=trim(Result_Folder)//'Transition_Debt.txt', STATUS='old', POSITION='append')
	    	OPEN (UNIT=91, FILE=trim(Result_Folder)//'Transition_Wage.txt', STATUS='old', POSITION='append')
	        	WRITE(UNIT=77, FMT=*) NBAR_bench, NBAR2_tr, NBAR_exp
	        	WRITE(UNIT=78, FMT=*) QBAR_bench, QBAR2_tr, QBAR_exp
	        	WRITE(UNIT=79, FMT=*) R_bench, R_tr, R_exp
	        	WRITE(UNIT=80, FMT=*) GBAR_bench, GBAR_tr, GBAR_exp
	        	WRITE(UNIT=81, FMT=*) GBAR_K_tr
	        	WRITE(UNIT=82, FMT=*) GBAR_W_tr
	        	WRITE(UNIT=83, FMT=*) GBAR_L_tr
	        	WRITE(UNIT=84, FMT=*) GBAR_C_tr
	        	WRITE(UNIT=85, FMT=*) SSC_Payments_tr
	        	WRITE(UNIT=86, FMT=*) K_bench, K_tr, K_exp
	        	WRITE(UNIT=87, FMT=*) C_bench, C_tr, C_exp
	        	WRITE(UNIT=88, FMT=*) Wealth_Top_1_Tr
	        	WRITE(UNIT=89, FMT=*) Wealth_Top_10_Tr
	        	WRITE(UNIT=90, FMT=*) 0, Debt_tr, Debt_SS
	        	WRITE(UNIT=90, FMT=*) wage_bench, wage_tr, wage_exp
	    	CLOSE (unit=77); CLOSE (unit=78); CLOSE (unit=79);
	    	CLOSE (unit=80); CLOSE (unit=81); CLOSE (unit=82); CLOSE (unit=83); CLOSE (unit=84); CLOSE (unit=85);
	    	CLOSE (unit=86); CLOSE (unit=87); CLOSE (unit=88); CLOSE (unit=89); CLOSE (unit=90); CLOSE (unit=91);


	    	! Write Variable Paths
	    	print*,' '
			print*,' 	--------------------------------------'
			print*,' 	Printing Variables'
			print*,' 	--------------------------------------'
			! OPEN  (UNIT=1,  FILE=trim(Result_Folder)//'Cons_tr'  , STATUS='replace')
			! OPEN  (UNIT=2,  FILE=trim(Result_Folder)//'Hours_tr'  , STATUS='replace')
			! OPEN  (UNIT=3,  FILE=trim(Result_Folder)//'Aprime_tr'  , STATUS='replace')
			OPEN  (UNIT=77,  FILE=trim(Result_Folder)//'GBAR_tr'   , STATUS='replace')
			OPEN  (UNIT=78,  FILE=trim(Result_Folder)//'Wage_tr'   , STATUS='replace')
			OPEN  (UNIT=79,  FILE=trim(Result_Folder)//'R_tr'      , STATUS='replace')
			OPEN  (UNIT=80,  FILE=trim(Result_Folder)//'P_tr'      , STATUS='replace')
			OPEN  (UNIT=81,  FILE=trim(Result_Folder)//'QBAR_tr'   , STATUS='replace')
			OPEN  (UNIT=82,  FILE=trim(Result_Folder)//'NBAR_tr'   , STATUS='replace')
			OPEN  (UNIT=83,  FILE=trim(Result_Folder)//'YBAR_tr'   , STATUS='replace')
			OPEN  (UNIT=84,  FILE=trim(Result_Folder)//'K_tr'      , STATUS='replace')
			OPEN  (UNIT=85,  FILE=trim(Result_Folder)//'C_tr'      , STATUS='replace')
			OPEN  (UNIT=86,  FILE=trim(Result_Folder)//'Debt_tr'   , STATUS='replace')
			OPEN  (UNIT=87,  FILE=trim(Result_Folder)//'tauW_at_tr', STATUS='replace')
			OPEN  (UNIT=88,  FILE=trim(Result_Folder)//'tauL_tr'   , STATUS='replace')
			OPEN  (UNIT=89,  FILE=trim(Result_Folder)//'tauK_tr'   , STATUS='replace')
			OPEN  (UNIT=90,  FILE=trim(Result_Folder)//'tauC_tr'   , STATUS='replace')
			OPEN  (UNIT=91,  FILE=trim(Result_Folder)//'Debt_SS'   , STATUS='replace')
			OPEN  (UNIT=92,  FILE=trim(Result_Folder)//'DBN_SS'   , STATUS='replace')
				! WRITE (UNIT=1,  FMT=*) Cons_tr
				! WRITE (UNIT=2,  FMT=*) Hours_tr
				! WRITE (UNIT=3,  FMT=*) Aprime_tr
				WRITE (UNIT=77,  FMT=*) GBAR_tr
				WRITE (UNIT=78,  FMT=*) Wage_tr
				WRITE (UNIT=79,  FMT=*) R_tr
				WRITE (UNIT=80,  FMT=*) P_tr
				WRITE (UNIT=81,  FMT=*) QBAR_tr
				WRITE (UNIT=82,  FMT=*) NBAR_tr
				WRITE (UNIT=83,  FMT=*) YBAR_tr
				WRITE (UNIT=84,  FMT=*) K_tr
				WRITE (UNIT=85,  FMT=*) C_tr
				WRITE (UNIT=86,  FMT=*) Debt_tr
				WRITE (UNIT=87,  FMT=*) tauW_at
				WRITE (UNIT=88,  FMT=*) 1.0_dp-psi
				WRITE (UNIT=89,  FMT=*) tauK
				WRITE (UNIT=90,  FMT=*) tauC
				WRITE (UNIT=91,  FMT=*) Debt_SS
				WRITE (UNIT=92,  FMT=*) DBN_exp
			! CLOSE (unit=1); CLOSE (unit=2); CLOSE (unit=3); 
			CLOSE (unit=77); CLOSE (unit=78); CLOSE (unit=79);
	    	CLOSE (unit=80); CLOSE (unit=81); CLOSE (unit=82); CLOSE (unit=83); CLOSE (unit=84); 
	    	CLOSE (unit=85); CLOSE (unit=86); CLOSE (unit=87); CLOSE (unit=88); CLOSE (unit=89); CLOSE (unit=90); 
	    	CLOSE (unit=91); CLOSE (unit=92); 
			print*,' 	--------------------------------------'
			print*,' 	Variable Printing Completed'
			print*,' 	--------------------------------------'
			print*,' '

		! Compare distance to tax reform distribution
		    DBN_dist = maxval(abs(DBN_tr(:,:,:,:,:,:,T+1)-DBN_exp))
		    Chg_dist = abs(DBN_dist-Old_DBN_dist)
		    Old_DBN_dist = DBN_dist

		    print*, ' '; print*,'-------------------------------------------------------------------'
		    print*, 'Iteration=',simutime
		    print '(A,E12.5,X,X,A,E12.5,X,X,A,E12.5,X,X,A,E12.5,X,X,A,E12.5,X,X,A,E12.5,X,X)', &
		    	&'	Distance: DBN=', DBN_dist,' Q=',Q_dist,' NW=',NW_dist,' R=',R_dist,' Db=',Db_dist,'Chg_dist=',Chg_dist
		    print '(A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X)', &
	    		&'	X(T)/X(SS): Q=',100.0_dp*(QBAR2_tr(T+1)/QBAR_exp-1),' N=',100.0_dp*(NBAR2_tr(T+1)/NBAR_exp-1),&
		    		' W=',100.0_dp*(Wage_tr(T+1)-wage_exp),' R=',100.0_dp*(R2_tr(T+1)-R_exp),' Db=',100.0_dp*(Debt_tr(T+1)/Debt_exp-1)
	    	print*, ' '
	    	print*, ' Government Budget:'
	    	print '(A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X)', &
	    		&'	tau_K=',100.0_dp*tauK,'tau_W',100.0_dp*tauW_at,'tau_L',100.0_dp*(1.0_dp-psi)
	    	print '(A,F7.3,X,X,A,F7.3)', '	Debt=',Debt_tr(T+1),'Debt/GDP=',Debt_tr(T+1)/YBAR_tr(T+1)
	    	print '(A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X)', &
	    		&'	GBAR_T+1=',GBAR_tr(T+1),'Expenditure',GBAR_bench+R_tr(T+1)*Debt_tr(T+1),&
	    					& 'Deficit=',GBAR_tr(T+1)-GBAR_bench-R_tr(T+1)*Debt_tr(T+1)
			print '(A,F7.3,X,X,A,F7.3,X,X,A,F7.3,X,X)', &
				&'	GBAR_exp=',GBAR_exp,'Expenditure',GBAR_bench+R_exp*Debt_tr(T+1),&
	    					& 'Deficit=',GBAR_exp-GBAR_bench-R_exp*Debt_tr(T+1)
			print '(A,F7.4)', '	Debt Absortion=',Debt_Absorption
			print*,'-------------------------------------------------------------------';print*, ' '

	    	OPEN (UNIT=76, FILE=trim(Result_Folder)//'Transition_Distance.txt', STATUS='old', POSITION='append')
	    	WRITE(UNIT=76, FMT=*) simutime,DBN_dist,Q_dist,NW_dist,R_dist,Db_dist,&
	    						&	100.0_dp*(QBAR2_tr(T+1)/QBAR_exp-1),100.0_dp*(NBAR2_tr(T+1)/NBAR_exp-1),&
	    						&	100.0_dp*(Wage_tr(T+1)/wage_exp-1),100.0_dp*(R2_tr(T+1)/R_exp-1),100.0_dp*(Debt_tr(T+1)/Debt_exp-1)
	    	CLOSE(UNIT=76)


	    	! Print Summary File
			OPEN  (UNIT=78,FILE=trim(Result_Folder)//'Transition_Summary.txt'   , STATUS='replace')
			WRITE (UNIT=78,FMT=*) 'Period Q N R Wage Y K C Y_C L_C K_C Y_P L_P K_P Debt ', &
					&'GBAR GBAR_K GBAR_W GBAR_L GBAR_C SSC  Wealth_Top_1 Wealth_Top_10'
			WRITE (UNIT=78,FMT=*) 'SS_1',QBAR_bench,NBAR_bench,R_bench,wage_bench & 
										& ,Y_bench,K_bench,C_bench &
										& ,YBAR_C_bench,L_C_bench,K_C_bench,YBAR_P_bench,L_P_bench,K_P_bench &
										& ,0,GBAR_bench
			do ti=1,T+1
			WRITE (UNIT=78,FMT=*) ti,QBAR_tr(ti),NBAR_tr(ti),R_tr(ti),Wage_tr(ti) & 
								& ,YBAR_tr(ti),K_tr(ti),C_tr(ti) &
								& ,YBAR_C_tr(ti),L_C_tr(ti),K_C_tr(ti),YBAR_P_tr(ti),L_P_tr(ti),K_P_tr(ti) &
								& ,Debt_tr(ti),GBAR_tr(ti) &
								& ,GBAR_K_tr(ti),GBAR_W_tr(ti),GBAR_L_tr(ti),GBAR_C_tr(ti),SSC_Payments_tr(ti) &
								& ,Wealth_Top_1_tr(ti),Wealth_Top_10_tr(ti)
			enddo 
			WRITE (UNIT=78,FMT=*) 'SS_2',QBAR_exp,NBAR_exp,R_exp,wage_exp & 
								&  ,Y_exp,K_exp,C_exp &
								&  ,YBAR_C_exp,L_C_exp,K_C_exp,YBAR_P_exp,L_P_exp,K_P_exp &
								&  ,Debt_exp,GBAR_exp,'99 99 99 99 99 99' &
								&  ,Wealth_Top_1_tr(T+1),Wealth_Top_10_tr(T+1)
			CLOSE (UNIT=78);


	    simutime  = simutime +1 
	 
	ENDDO ! WHILE

	print*,'---------------------------------------------------'
	print*,' 	Transition Completed'
	print*,'---------------------------------------------------'
	print*,' '

	



END SUBROUTINE FIND_DBN_Transition

!========================================================================================
!========================================================================================
!========================================================================================

SUBROUTINE EGM_Transition()
	use omp_lib
	IMPLICIT NONE
	REAL(DP) :: brentvalue, C_foc, H_min, euler_power, chi_aux
	INTEGER  :: tempai, sw 
	REAL(DP), DIMENSION(na_t+nz*nx+1)  	:: EndoCons, EndoYgrid, EndoHours
	INTEGER , DIMENSION(na_t+nz*nx+1)   :: sort_ind 
	REAL(DP), DIMENSION(na_t,nz,nx) :: Wealth_mat
	REAL(DP), DIMENSION(7)       	:: state_FOC
	REAL(DP), DIMENSION(7)       	:: par_FOC
	REAL(DP), DIMENSION(nx)       	:: MB_aprime_t
	integer  :: age, ai, zi, lambdai, ei, xi, xp_ind
	real(dp), dimension(na_t+nz*nx+1) :: EndoYgrid_sort

	print*,' '
	print*,' 	----------------------------'
	print*,' 	Starting EGM Transition'
	print*,' 	----------------------------'

	!$ call omp_set_num_threads(nz)

	! Set a minimum value for labor to check in the FOC
		H_min = 0.000001_dp

	! Set the power used in the Euler equation for the retirement period
		if (NSU_Switch.eqv..true.) then 
			! Non-Separable Utility
				euler_power = (1.0_dp/((1.0_dp-sigma)*gamma-1.0_dp))
		else 
			! Separable Utility
				euler_power = (-1.0_dp/sigma)
		end if 

	! Set value of auxiliary chi parameter for last period's saving
		chi_aux = ((1.0_dp+tauC)*chi_bq)**(1.0_dp/(1-gamma*(1.0_dp-sigma)))/(1.0_dp+tauC)

	! Policy functions set to tax reform in final period
		Cons_tr(:,:,:,:,:,:,T+1)   = Cons_exp*(1.0_DP+tauC) ; Cons_t_pr   = Cons_exp*(1.0_DP+tauC) ;
		Hours_tr(:,:,:,:,:,:,T+1)  = Hours_exp  			; Hours_t_pr  = Hours_exp  			   ;
		Aprime_tr(:,:,:,:,:,:,T+1) = Aprime_exp 			; 
		! print*,'Cons_tr(T+1)=',Cons_tr(81,:,5,3,3,1,T+1)
		! print*,"|Const_exp-Const_bench|=",maxval(abs(Cons_exp-Cons_bench))

	! Solve backwards for all transition periods
	do ti=T,1,-1
		! print*,' Solving EGM for transition period ',ti


	! Grids and auxiliary variables, must be solved per period 

		! Solve the model at current aggregate values
		! Compute labor units 
			CALL ComputeLaborUnits(Ebar_tr(ti), wage_tr(ti)) 
		! Compute Capital demand and Profits by (a,z)
			K_mat  = K_Matrix(R_tr(ti),P_tr(ti))
			Pr_mat = Profit_Matrix(R_tr(ti),P_tr(ti))
		! Compute wealth given current R and P
			Wealth_mat = Wealth_Matrix_t(R_tr(ti+1),P_tr(ti+1))
		! Form YGRID for the capital income economy given interest rate "P"
			CALL FORM_Y_MB_GRID_Transition(YGRID,MBGRID,YGRID_t,MBGRID_t,ti)
			! print*, 'YGRID_t(:,5,1)=',YGRID_t(:,5,1)
			! print*, 'RetY_lambda_e=',RetY_lambda_e
			! print*,"Period",ti,"|YGRID-YGRID_bench|=",maxval(abs(YGRID_t-YGRID_aux)),&
			! 	& "|RetY-RetY_bench|=",maxval(abs(RetY_lambda_e-RetY_lambda_e_aux))
			! print*,"Period",ti,"max(YGRID)=",maxval(YGRID_t),"max(YGRID_bench)=",maxval(YGRID_aux),&
			! 	& "max(RetY)=",maxval(RetY_lambda_e),"max(RetY_bench)=",minval(RetY_lambda_e_aux)

		! print*, 'R=',R,'P=',P, 'W=',wage, 'na=', na, 'na_t=', na_t
	!========================================================================================
	!------RETIREMENT PERIOD-----------------------------------------------------------------

	! SET HOURS TO ZERO SO THAT RETIRED HAS ZERO HOURS
	Hours_t = 0.0_DP

	! Last period of life
	age=MaxAge
		! print*,' 	Age=',age
	!$omp parallel do private(lambdai,ei,ai,xi)
	DO zi=1,nz
	!$omp parallel do private(lambdai,ei,ai)
	DO xi=1,nx
    DO lambdai=1,nlambda
    DO ei=1,ne
    DO ai=1,na_t
    	if ((YGRID_t(ai,zi,xi) + RetY_lambda_e(lambdai,ei)).le.(bq_0/chi_aux) ) then 
        Cons_t(age,ai,zi,lambdai,ei,xi)   =  YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei) 
        else
        Cons_t(age,ai,zi,lambdai,ei,xi)   = (YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei)+bq_0)/(1.0_dp+chi_aux)
        endif
        Aprime_t(age,ai,zi,lambdai,ei,xi) = YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei)-Cons_t(age,ai,zi,lambdai,ei,xi)
	ENDDO ! ai
    ENDDO ! ei
	ENDDO ! lambdai
	ENDDO ! xi
	ENDDO ! zi
	! print*,"Period",ti,"|Const_t-Const_bench|=",maxval(abs(Cons_t(MaxAge,:,:,:,:,:)-Cons_t_pr(MaxAge,:,:,:,:,:)))
	! print*,"Period",ti,"Test-Const_bench=",(YGRID_aux(50,5,2) + RetY_lambda_e_aux(3,3))/(1.0_DP+tauC)-Cons_bench(MaxAge,50,5,3,3,2)
	
	! Rest of retirement
	DO age=MaxAge-1,RetAge,-1
		! print*,' 	Age=',age
	!$omp parallel do private(lambdai,ei,ai,xi,xp_ind,EndoCons,EndoYgrid,sw,sort_ind,tempai, &
	!$omp& state_FOC,par_FOC,MB_aprime_t,EndoYgrid_sort)
    DO zi=1,nz
    !$omp parallel do private(lambdai,ei,ai,xp_ind,EndoCons,EndoYgrid,sw,sort_ind,tempai, &
	!$omp& state_FOC,par_FOC,MB_aprime_t,EndoYgrid_sort)
    DO xi=1,nx
    DO lambdai=1,nlambda
    DO ei=1,ne
    	! Endogenous grid and consumption are initialized
	    EndoCons  = big_p 
	    EndoYgrid = big_p 
	    sw 		  = 0
    DO ai=1,na_t 
		if (abs(agrid_t(ai)-Y_a_threshold).lt.1e-8) then 
			! print*, ' Threshold section - EGM Retirement'
			sw 			  = sw+1	

    		! Consumption on endogenous grid and implied asset income under tauW_bt
    		do xp_ind = 1,nx 	
				P = P_tr(ti+1); R = R_tr(ti+1) 
				MB_aprime_t(xp_ind) = MB_a_bt(agrid_t(ai),zi,xp_ind)
    		enddo 
    		EndoCons(ai)  =  (beta*survP(age)* 	&
    			& sum(pr_x(xi,:,zi,age)*MB_aprime_t*Cons_t_pr(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power)) &
				& + (1.0_dp-survP(age))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
		    	& *   chi_bq*((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ) **euler_power
	        EndoYgrid(ai) = agrid_t(ai) +  EndoCons(ai) - RetY_lambda_e(lambdai,ei)
	        ! Consumption on endogenous grid and implied asset income under tauW_at
	        do xp_ind = 1,nx 	
				P = P_tr(ti+1); R = R_tr(ti+1) 
				MB_aprime_t(xp_ind) = MB_a_at(agrid_t(ai),zi,xp_ind)
    		enddo 
	        EndoCons(na_t+sw)  = (beta*survP(age)*	&
	        	& sum(pr_x(xi,:,zi,age)*MB_aprime_t*Cons_t_pr(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power)) &
				& + (1.0_dp-survP(age))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
		    	& *   chi_bq*((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ) **euler_power
	    	EndoYgrid(na_t+sw) = agrid_t(ai) +  EndoCons(na_t+sw) - RetY_lambda_e(lambdai,ei)

	  !   	print*, ' '
			! print*, ' Threshold test - Retirement'
			! print*, ' Current State', age, ai, zi, lambdai, ei, xi
			! print*, ' ', (pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold))
			! print*, ' ', (any((pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold)).lt.1e-8))
			! print*, ' ', MBGRID_t(ai,zi,:)
			! print*, ' ', MB_a_bt(agrid(ai),zi,xi)
			! print*, ' ', MB_a_at(agrid(ai),zi,xi)
			! print*, ' ', EndoCons(ai)
			! print*, ' ', EndoCons(na_t+sw)
			! print*, ' '
	    else 
	    	! Consumption on endogenous grid and implied asset income
	    	EndoCons(ai)  = (beta*survP(age)* 	&
	    		& sum(pr_x(xi,:,zi,age)*MBGRID_t(ai,zi,:)*Cons_t_pr(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power))&
				& + (1.0_dp-survP(age))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
		    	& *   chi_bq*((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) )**euler_power
	        EndoYgrid(ai) = agrid_t(ai) +  EndoCons(ai) - RetY_lambda_e(lambdai,ei)
	        ! !$omp critical
	        ! print*,' Standard EGM - State:',age,ai,zi,lambdai,ei,xi,ti
	        ! print*,' 	EndoCons(ai)=',EndoCons(ai)
	        ! print*,' 	Expected value=',&
	        ! 	& sum(pr_x(xi,:,zi,age)*MBGRID_t(ai,zi,:)*Cons_t_pr(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power))
	        ! print*,' 	Cons_t+1=',Cons_t_pr(age+1,ai,zi,lambdai,ei,:)
	        ! print*,' 	Cons_t+1=',Cons_exp(age+1,ai,zi,lambdai,ei,:)
	        ! !$omp end critical
	    end if 


	 !    if (any(isnan(EndoCons))) then 
	 !    	print*,' '
		! 	print*,' '
		! 	print*,' '
		! 	print*,' Endo Consumption in retirement', age,ai,zi,lambdai,ei,xi
		! 	print*,' ',EndoCons(ai),(beta*survP(age)* 	&
	 !    				& sum(pr_x(xi,:,zi,age)*MBGRID_t(ai,zi,:)*Cons_t(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power)) ) **euler_power
		! 	print*,' MBGRID_t=',MBGRID_t(ai,zi,:)
		! 	print*,' cons(t+1)=',Cons_t(age+1,ai,zi,lambdai,ei,:)
		! 	print*,' Size(endoCons)=',size(EndoCons)
		! 	print*,' sw=',sw,'na_t=',na_t
		! 	print*,' Threshold_test=',(any((pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold)).lt.1e-8))
		! 	print*,' '
		! 	print*,' ',EndoCons
		! 	print*,' '
		! 	STOP
		! endif 
	ENDDO ! ai

	    !$omp critical
		if (any(isnan(EndoCons))) then 
			print*, "isnan - Consumption endogenous"
			print*, age,zi,ai,lambdai,ei,xi,ti
			print*, EndoCons
			STOP 
		end if 
		!$omp end critical

		
	
	! Sort endogenous grid for interpolation
	call Sort(na_t+nx*nz+1,EndoYgrid,EndoYgrid_sort,sort_ind)
	! print*,' Yendo'
	! print*, EndoYgrid
	! print*,' Yendo_sort'
	! print*, EndoYgrid_sort
	! print*,' indices'
	! print*, sort_ind

	EndoYgrid = EndoYgrid_sort
	EndoCons = EndoCons(sort_ind)

	! print*, ' '
	! print*, ' isnan(endocons)', any(isnan(EndoCons))

	! Find  decision rules on exogenous grids
		! decision rules are obtained taking care of extrapolations
		tempai=1           
		DO WHILE ( YGRID_t(tempai,zi,xi) .lt. EndoYgrid(1) )
            tempai = tempai + 1
        ENDDO
	                
		DO ai=tempai,na_t              
		    ! CONSUMPTION ON EXOGENOUS GRIDS 
		    Cons_t(age,ai,zi,lambdai, ei,xi)  = Linear_Int(EndoYgrid(1:na_t+sw), EndoCons(1:na_t+sw),na_t+sw, YGRID_t(ai,zi,xi))
		    Aprime_t(age,ai,zi,lambdai,ei,xi) = &
		    								& YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei)-Cons_t(age, ai, zi, lambdai, ei,xi)
		    
		    If (Aprime_t(age,ai,zi,lambdai,ei,xi).lt.amin) then
		    	Aprime_t(age,ai,zi,lambdai,ei,xi) = amin
	            Cons_t(age,ai,zi,lambdai,ei,xi)   = &
	            	& YGRID_t(ai,zi,xi) + RetY_lambda_e(lambdai,ei) - Aprime_t(age,ai,zi,lambdai,ei,xi)
				IF (Cons_t(age,ai,zi,lambdai,ei,xi) .le. 0.0_DP)  THEN
				    print*,'r1: Cons(age, ai, zi, lambdai,ei)=',Cons_t(age, ai, zi, lambdai, ei,ti)
				ENDIF                   
	        endif 

	        ! ! !$omp critical 
	        ! if (isnan(Cons_t(age,ai,zi,lambdai,ei,xi))) then
	        ! 	print*,' '
	        ! 	print*,' '
	        ! 	print*,' '
	        ! 	print*,' isnan Consumption', age,ai,zi,lambdai,ei,xi
	        ! 	print*,' sw=',sw 
	        ! 	print*,' Yendo'
	        ! 	print*, EndoYgrid(1:na_t+sw)
	        ! 	print*,' Cendo'
	        ! 	print*, EndoCons(1:na_t+sw)
	        ! 	print*,' YGRID'
	        ! 	print*, YGRID_t(ai,zi,xi)
	        ! 	print*,' Linear Interpolation',Linear_Int(EndoYgrid(1:na_t+sw), EndoCons(1:na_t+sw),na_t+sw, YGRID_t(ai,zi,xi))
	        ! 	print*,' Aprime=',Aprime_t(age,ai,zi,lambdai,ei,xi)
	        ! 	call Sort(na_t+nx*nz+1,EndoYgrid,EndoYgrid_sort,sort_ind)
	        ! 	print*,' Yendo'
	        ! 	print*, EndoYgrid_sort
	        ! 	print*,' indices'
	        ! 	print*, sort_ind
	        ! 	print*, ' '
	        ! 	print*, ' ',minval(EndoYgrid),maxval(EndoYgrid)
	        ! 	print*, ' ',minval(EndoYgrid_sort),maxval(EndoYgrid_sort)
	        ! 	print*, ' The end'
	        ! 	STOP
	        ! endif 
	        ! ! !$omp end critical
		ENDDO ! ai  

        ai=1           
        DO WHILE ( YGRID_t(ai,zi,xi) .lt. EndoYgrid(1) )
			! Solve for a' directly by solving the Euler equation for retirement FOC_R
			state_FOC  = (/age,ai,zi,lambdai,ei,xi,ti/)
			brentvalue = brent_p(min(amin,YGRID_t(ai,zi,xi)), (amin+YGRID_t(ai,zi,xi))/2.0_DP , &
			                & YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei) *0.95_DP, &
			                & FOC_R_Transition, brent_tol, Aprime_t(age,ai,zi,lambdai,ei,xi),state_FOC)
			
			Cons_t(age,ai,zi,lambdai,ei,xi) = &
					& YGRID_t(ai,zi,xi) + RetY_lambda_e(lambdai,ei) - Aprime_t(age,ai,zi,lambdai,ei,xi)
			
			IF (Cons_t(age, ai, zi, lambdai, ei,xi) .le. 0.0_DP)  THEN
			    print*,'r2: Cons(age,ai,zi,lambdai,ei,xi,ti)=',Cons_t(age,ai,zi,lambdai,ei,xi)
			    print*,'Aprime(age,ai,zi,lambdai,ei,xi,ti)=',Aprime_t(age,ai,zi,lambdai,ei,xi)
			    print*,'YGRID(ai,zi,xi)+RetY_lambda_e(lambdai,ei)=',YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei)
			    print*,'YGRID(ai,zi,xi)=',YGRID(ai,zi,xi),'EndoYgrid(1)=',EndoYgrid(1)
			    print*,'RetY_lambda_e(lambdai,ei)=',RetY_lambda_e(lambdai,ei)
			    print*,'lambdai=',lambdai
			ENDIF                   
			ai = ai + 1
		ENDDO  
	              
    ENDDO ! ei     
    ENDDO ! lambda
    ENDDO ! xi 
	ENDDO ! zi
    ENDDO !age


    ! print*, ' 			Retirement Period Ends'
	!------RETIREMENT PERIOD ENDS------------------------------------------------------------
	!========================================================================================
	
	!========================================================================================
	!------Working Period Starts-------------------------------------------------------------
	! print*, ' 			Working Period Starts'

	DO age=RetAge-1,1,-1
		! print*,' 	Age=',age
	!$omp parallel do private(lambdai,ei,ai,xi,xp_ind,EndoCons,EndoHours,EndoYgrid,sw,sort_ind,tempai, &
	!$omp& C_foc,state_FOC,par_FOC,MB_aprime_t)
    DO zi=1,nz
    !$omp parallel do private(lambdai,ei,ai,xp_ind,EndoCons,EndoHours,EndoYgrid,sw,sort_ind,tempai, &
	!$omp& C_foc,state_FOC,par_FOC,MB_aprime_t)
    DO xi=1,nx
    DO lambdai=1,nlambda
    DO ei=1,ne	
    	! Endogenous grid and consumption are initialized
	    EndoCons  = big_p 
	    EndoYgrid = big_p 
	    EndoHours = big_p 
	    sw 		  = 0                    
        DO ai=1,na_t
        state_FOC  = (/age,ai,zi,lambdai,ei,xi,ti/)
		if (abs(agrid_t(ai)-Y_a_threshold).lt.1e-8) then 
			sw 			  = sw+1	

	    	! Below threshold
    		do xp_ind = 1,nx 	
				P = P_tr(ti+1); R = R_tr(ti+1) 
				MB_aprime_t(xp_ind) = MB_a_bt(agrid_t(ai),zi,xp_ind)
    		enddo 
			call EGM_Working_Period_Transition( MB_aprime_t , H_min , state_FOC , & 
			      & EndoCons(ai), EndoHours(ai) , EndoYgrid(ai)  )
			
			! Above threshold
			do xp_ind = 1,nx 	
				P = P_tr(ti+1); R = R_tr(ti+1) 
				MB_aprime_t(xp_ind) = MB_a_at(agrid_t(ai),zi,xp_ind)
    		enddo 
			call EGM_Working_Period_Transition( MB_aprime_t , H_min , state_FOC , & 
			      & EndoCons(na_t+sw), EndoHours(na_t+sw) , EndoYgrid(na_t+sw)  )

	    	!print*, ' '
	    	!print*, State_FOC 
	    	!print*, MB_a_bt(agrid_t(ai),zgrid(zi)), MB_a_at(agrid_t(ai),zgrid(zi))
	  !   	print*, ' '
			! print*, ' Threshold test - Working Period'
			! print*, ' Current State', age, ai, zi, lambdai, ei, xi
			! print*, ' ', (pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold))
			! print*, ' ', (any((pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold)).lt.1e-8))
			! print*, ' ', MBGRID_t(ai,zi,:)
			! print*, ' ', MB_a_bt(agrid(ai),zi,xi)
			! print*, ' ', MB_a_at(agrid(ai),zi,xi)
			! print*, ' ', MB_aprime_t
			! print*, ' ', EndoCons(ai)
			! print*, ' ', EndoCons(na_t+sw)
			! print*, ' '
		else 
			! Usual EGM
			call EGM_Working_Period_Transition( MBGRID_t(ai,zi,:) , H_min , state_FOC , & 
			      & EndoCons(ai), EndoHours(ai) , EndoYgrid(ai)  )
	
		end if 

    	ENDDO ! ai

	    ! !$omp critical
		if (any(isnan(EndoCons))) then 
			print*, "isnan - Consumption endogenous"
			print*, age,zi,ai,lambdai,ei,xi,ti
			print*, EndoCons
			STOP 
		end if 
		! !$omp end critical

    ! Sort endogenous grid for interpolation
	call Sort(na_t+nx*nz+1,EndoYgrid,EndoYgrid,sort_ind)
	EndoHours = EndoHours(sort_ind)
	EndoCons  = EndoCons(sort_ind)
	! 	print*, ' '
	! 	do ai=1,na_t+1
	! 		print*, EndoHours(ai), EndoCons(ai), EndoYgrid(ai)
	! 	end do
	! 	print*, ' '

		if (any(isnan(EndoCons)).or.any(isnan(EndoHours)).or.any(isnan(EndoYgrid))) then 
			print*, "isnan - Consumption working 4"
			print*, age,lambdai,ai,zi,ei,xi,ti
			STOP 
		end if 

    	! Find  decision rules on exogenous grids
        tempai=1           
        DO WHILE ( YGRID_t(tempai,zi,xi) .lt. EndoYgrid(1) )
              tempai = tempai + 1
        ENDDO
	                
		! decision rules are obtained taking care of extrapolations
		DO ai=tempai,na_t 
			! Interpolate for value of consumption in exogenous grid
			Cons_t(age,ai,zi,lambdai,ei,xi) = Linear_Int(EndoYgrid(1:na_t+sw),EndoCons(1:na_t+sw),na_t+sw,YGRID_t(ai,zi,xi))

			if (Progressive_Tax_Switch.eqv..true.) then
			! Check FOC for implied consumption 
				if (NSU_Switch.eqv..true.) then 
					! Non-Separable Utility
					C_foc = (gamma/(1.0_dp-gamma))*(1.0_dp-H_min)*MB_h(H_min,age,lambdai,ei,wage_tr(ti))
				else
					! Separable Utility
					C_foc = (MB_h(H_min,age,lambdai,ei,wage_tr(ti))*(1.0_dp-H_min)**(gamma)/phi)**(1.0_dp/sigma)
				end if 
			! Hours
				if (Cons_t(age,ai,zi,lambdai,ei,xi).ge.C_foc) then
					Hours_t(age,ai,zi,lambdai,ei,xi) = 0.0_dp
				else
					! Auxiliary variables for solving FOC for hours 
					par_FOC(1:6) = (/age,ai,zi,lambdai,ei,xi/) 
					par_FOC(7)   = Cons_t(age,ai,zi,lambdai,ei,xi)
					! FOC (set wage to current wage)
					wage = wage_tr(ti)
					brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_H, brent_tol,Hours_t(age,ai,zi,lambdai,ei,xi) , par_FOC) 
				end if
			else 
				if (NSU_Switch.eqv..true.) then 
					Hours_t(age,ai,zi,lambdai,ei,xi) = max( 0.0_dp , &
						&  1.0_DP - (1.0_DP-gamma)*Cons_t(age,ai,zi,lambdai,ei,xi)/(gamma*psi*yh(age,lambdai,ei)) )
				else 
					Hours_t(age,ai,zi,lambdai,ei,xi) = max( 0.0_DP , &
						&  1.0_DP - phi*Cons_t(age,ai,zi,lambdai,ei,xi)/(psi*yh(age, lambdai,ei)) )
				end if 
			end if 

			if (Hours_t(age,ai,zi,lambdai,ei,xi).ge.1.0_dp) then 
				print*, ' '
				print*, 'Hours highers than 1', age,ai,zi,lambdai,ei,xi,ti
				print*, Cons_t(age,ai,zi,lambdai,ei,xi)
				STOP
			endif 


			! Savings 
				Aprime_t(age,ai,zi,lambdai,ei,xi) = & 
								& YGRID_t(ai,zi,xi)  + Y_h(Hours_t(age,ai,zi,lambdai,ei,xi),age,lambdai,ei,wage_tr(ti))  & 
            					& - Cons_t(age,ai,zi,lambdai,ei,xi) 
		                    
		    If (Aprime_t(age,ai,zi,lambdai,ei,xi)  .lt. amin) then

		    	! print*, ' Aprime was below minimum!!!!
            	Aprime_t(age,ai,zi,lambdai,ei,xi) = amin
		         
	           	! Compute hours
		        if (NSU_Switch.and.(Progressive_Tax_Switch.eqv..false.)) then 
		        	Hours_t(age,ai,zi,lambdai,ei,xi)  = max( 0.0_dp , &
		        	&  gamma - (1.0_DP-gamma)*(YGRID_t(ai,zi,xi)-Aprime_t(age,ai,zi,lambdai,ei,xi))/(psi*yh(age,lambdai,ei)))
                else               	        
                	!compute  hours using FOC_HA                              
		        	par_FOC(1:6) = (/age,ai,zi,lambdai,ei,xi/) 
		        	par_FOC(7)   = amin
					wage = wage_tr(ti)
		        	brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_HA, brent_tol, Hours_t(age,ai,zi,lambdai,ei,xi),par_FOC)
		        end if  

	            Cons_t(age,ai,zi,lambdai,ei,xi) = &
	            					& YGRID_t(ai,zi,xi)+  Y_h(Hours_t(age,ai,zi,lambdai,ei,xi),age,lambdai,ei,wage_tr(ti))  &
		                            & - Aprime_t(age,ai,zi,lambdai,ei,xi)

				IF (Cons_t(age,ai,zi,lambdai,ei,xi) .le. 0.0_DP)  THEN
					print*,'w1: Cons(age,ai,zi,lambdai,ei,xi,ti)=',Cons_t(age,ai,zi,lambdai,ei,xi)
					STOP
				ENDIF 
				if (Hours_t(age,ai,zi,lambdai,ei,xi).ge.1.0_dp) then 
					print*, ' '
					print*, 'w1 Hours highers than 1', age,ai,zi,lambdai,ei,xi,ti
					print*, Cons_t(age,ai,zi,lambdai,ei,xi)
					STOP
				endif                   
		     endif  
		    !$omp critical
		    IF (Cons_t(age,ai,zi,lambdai,ei,xi) .le. 0.0_DP)  THEN
				print*,'w1: Cons(age,ai,zi,lambdai,ei,xi)=',Cons_t(age,ai,zi,lambdai,ei,xi)
				STOP
			ENDIF 
			if (Hours_t(age,ai,zi,lambdai,ei,xi).ge.1.0_dp) then 
				print*, ' '
				print*, 'w1 Hours highers than 1', age,ai,zi,lambdai,ei,xi,ti
				print*, Cons_t(age,ai,zi,lambdai,ei,xi)
				STOP
			endif    
			!$omp end critical
		ENDDO ! ai   

		!$omp critical
		if (any(isnan(Cons_t(age,:,zi,lambdai,ei,xi)))) then 
			print*, "isnan - Consumption working 3"
			print*, age,lambdai,ai,zi,ei,xi,ti
			print*, 'Cons'
			print*, Cons_t(age,:,zi,lambdai,ei,xi)
			print*, 'Hours'
			print*, Hours_t(age,:,zi,lambdai,ei,xi) 
			STOP 
		end if         
		!$omp end critical        

		ai=1           
        DO WHILE ( YGRID_t(ai,zi,xi) .lt. EndoYgrid(1) )
        	! print*, ' Extrapolation between YGRID and EndoYgrid!!!!'
	        ! Solve for the Euler equation directly
	        state_FOC  = (/age,ai,zi,lambdai,ei,xi,ti/)
			brentvalue = brent_p( min(amin,YGRID_t(ai,zi,xi))   ,  (amin+YGRID_t(ai,zi,xi))/2.0_DP  ,  &
                             	& YGRID_t(ai,zi,xi)  + psi * ( yh(age, lambdai,ei)*0.95_DP )**(1.0_DP-tauPL)  ,  &
	                             & FOC_WH_Transition, brent_tol, Aprime_t(age,ai,zi,lambdai,ei,xi) , state_FOC )

			! Compute hours
	        if (NSU_Switch.and.(Progressive_Tax_Switch.eqv..false.)) then 
	        	Hours_t(age,ai,zi,lambdai,ei,xi)  = max( 0.0_dp , &
	        	&  gamma - (1.0_DP-gamma) *( YGRID_t(ai,zi,xi) - Aprime_t(age,ai,zi,lambdai,ei,xi))/(psi*yh(age,lambdai,ei)))
            else               	        
				!compute  hours using FOC_HA                              
				par_FOC(1:6) = (/age,ai,zi,lambdai,ei,xi/)
				par_FOC(7)   = Aprime_t(age,ai,zi,lambdai,ei,xi)
	            brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_HA, brent_tol, Hours_t(age,ai,zi,lambdai,ei,xi), par_FOC) 
 			end if 

            Cons_t(age,ai,zi,lambdai,ei,xi)=  &
            		& YGRID_t(ai,zi,xi) + Y_h(Hours_t(age,ai,zi,lambdai,ei,xi),age,lambdai,ei,wage_tr(ti))  &
					& - Aprime_t(age,ai,zi,lambdai,ei,xi)

           IF (Cons_t(age,ai,zi,lambdai,ei,xi) .le. 0.0_DP)  THEN
                print*,'w2:',age,zi,lambdai,ei,ai,xi,ti, 'Cons=',Cons_t(age,ai,zi,lambdai,ei,xi), &
                    & 'Aprime=',Aprime_t(age,ai,zi,lambdai,ei,xi), &
                    & 'Hours=', Hours_t(age,ai,zi,lambdai,ei,xi), &
                    & 'yh(age, lambdai,ei)=', yh(age, lambdai,ei),'YGRID(ai,zi)=',YGRID_t(ai,zi,xi)
                !pause
                print*, "there is a problem in line 2063 (Transition)"
                STOP
            ENDIF 

            if (Hours_t(age,ai,zi,lambdai,ei,xi).ge.1.0_dp) then 
				print*, ' '
				print*, 'w2 Hours highers than 1', age,ai,zi,lambdai,ei,xi,ti
				print*, Cons_t(age,ai,zi,lambdai,ei,xi)
				STOP
			endif 

            ai = ai + 1
        ENDDO  

	    if (any(isnan(Cons_t(age,:,zi,lambdai,ei,xi)))) then 
			print*, "isnan - Consumption working 2"
			print*, age,lambdai,ai,zi,ei,xi,ti
			print*, 'Cons'
			print*, Cons_t(age,:,zi,lambdai,ei,xi)
			print*, 'Hours'
			print*, Hours_t(age,:,zi,lambdai,ei,xi) 
			STOP 
		end if 

	                 
    ENDDO !ei         
    ENDDO !lambdai
    ENDDO ! xi 
	ENDDO !zi
    ENDDO !age


    ! print*, ' 			Clean-up and Interpolation'

    ! Update value for "next period's" policy functions
    Cons_t_pr   = Cons_t   ;
    Hours_t_pr  = Hours_t  ;

	! Interpolate to get values of policy functions on agrid (note that policy functions are defined on agrid_t)
	
	if (Y_a_threshold.eq.0.0_dp) then 
		! print*,' Assigning value of policy functions'
		Cons_tr(:,:,:,:,:,:,ti)   = Cons_t
		Hours_tr(:,:,:,:,:,:,ti)  = Hours_t
		Aprime_tr(:,:,:,:,:,:,ti) = Aprime_t
	else 
		! print*,' Interpolating policy functions'
		DO xi=1,nx
		DO age=1,MaxAge
	    DO lambdai=1,nlambda
	    DO zi=1,nz
	    DO ei=1,ne	                
	    DO ai=1,na
			Cons_tr(age,ai,zi,lambdai,ei,xi,ti)   = &
				& Linear_Int(YGRID_t(:,zi,xi) , Cons_t(age,:,zi,lambdai,ei,xi)   , na_t , YGRID(ai,zi,xi))
	    	Hours_tr(age,ai,zi,lambdai,ei,xi,ti)  = &
	    		& Linear_Int(YGRID_t(:,zi,xi) , Hours_t(age,:,zi,lambdai,ei,xi)  , na_t , YGRID(ai,zi,xi))
	    	Aprime_tr(age,ai,zi,lambdai,ei,xi,ti) = &
	    		& Linear_Int(YGRID_t(:,zi,xi) , Aprime_t(age,:,zi,lambdai,ei,xi) , na_t , YGRID(ai,zi,xi))
		ENDDO !ai         
		ENDDO !ei         
	    ENDDO !zi
	    ENDDO !lambdai
		ENDDO !age
		ENDDO !xi
	end if 

	enddo ! Time 

	! Adjust consumption by taxes
	Cons_tr = Cons_tr/(1.0_DP+tauC)

	! if (any(isnan(Cons_tr))) then 
	! 	print*, "isnan - Consumption"
	! 	STOP 
	! end if 
	! if (any(isnan(Hours_tr))) then 
	! 	print*, "isnan - Hours"
	! 	STOP 
	! end if 
	! if (any(isnan(Aprime_tr))) then 
	! 	print*, "isnan - Hours"
	! 	STOP 
	! end if 

	! 	print*, 'Policy Functions'
	! 	print*, ' YGRID ',' Cons ',' Hours ',' Aprime ',' A '
	! 	do ai=1,na
	! 		! print*, YGRID(ai,4), Cons(40,ai,4,3,3), Hours(40,ai,4,3,3), Aprime(40,ai,4,3,3), agrid(ai),&
	! 		! 		& Y_h(Hours(40,ai,4,3,3),40,3,3,wage),&
	! 		! 		& (1+tauC)*Cons(40,ai,4,3,3) + Aprime(40,ai,4,3,3) - YGRID(ai,4) -Y_h(Hours(40,ai,4,3,3),40,3,3,wage)
	! 	end do 
	! 	print*, ' '

	print*,' 	----------------------------'
	print*,' 	End of EGM Transition'
	print*,' 	----------------------------'
	print*,' '

END SUBROUTINE EGM_Transition


!========================================================================================
!========================================================================================
!========================================================================================
! For use with Agg_Debt to solve for equilibrium interest rate
SUBROUTINE EGM_Transition_aux(Ap_aux,time)
	use omp_lib
	IMPLICIT NONE
 	INTEGER , intent(in) :: time
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable, intent(out) :: Ap_aux
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: Ap_t_aux, Cons_aux, Hours_aux
	REAL(DP) :: brentvalue, C_foc, H_min, euler_power, chi_aux
	INTEGER  :: tempai, sw 
	REAL(DP), DIMENSION(na_t+nz*nx+1)  	:: EndoCons, EndoYgrid, EndoHours
	INTEGER , DIMENSION(na_t+nz*nx+1)   :: sort_ind 
	REAL(DP), DIMENSION(na_t,nz,nx) :: Wealth_mat
	REAL(DP), DIMENSION(7)       	:: state_FOC
	REAL(DP), DIMENSION(7)       	:: par_FOC
	REAL(DP), DIMENSION(nx)       	:: MB_aprime_t
	integer  :: age, ai, zi, lambdai, ei, xi, xp_ind
	real(dp), dimension(na_t+nz*nx+1) :: EndoYgrid_sort

	allocate( Ap_aux(MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Ap_t_aux(MaxAge,na_t,nz,nlambda,ne,nx) )
	allocate( Cons_aux(MaxAge,na_t,nz,nlambda,ne,nx) )
	allocate( Hours_aux(MaxAge,na_t,nz,nlambda,ne,nx) )

	! print*,' '
	! print*,' 	----------------------------------------'
	! print*,' 	Starting EGM Transition inside Agg_Debt'
	! print*,' 	----------------------------------------'

	!$ call omp_set_num_threads(nz)

	! Set a minimum value for labor to check in the FOC
		H_min = 0.000001_dp

	! Set the power used in the Euler equation for the retirement period
		if (NSU_Switch.eqv..true.) then 
			! Non-Separable Utility
				euler_power = (1.0_dp/((1.0_dp-sigma)*gamma-1.0_dp))
		else 
			! Separable Utility
				euler_power = (-1.0_dp/sigma)
		end if 

	! Set value of auxiliary chi parameter for last period's saving
		chi_aux = ((1.0_dp+tauC)*chi_bq)**(1.0_dp/(1-gamma*(1.0_dp-sigma)))/(1.0_dp+tauC)

	! Policy Functions for period t (time)
		! Note: these arrays should be adjusted to agrid_t 
		!       currently they work on agrid and we assume that
		!       agrid_t = agrid, which happens without threshold
		Cons_t_pr  = Cons_tr(:,:,:,:,:,:,time)*(1.0_DP+tauC);
		Hours_t_pr = Hours_tr(:,:,:,:,:,:,time) ;


	! Grids and auxiliary variables, must be solved per period 

		! Solve the model at current aggregate values
		wage = wage_tr(time-1)
		! Compute labor units 
			CALL ComputeLaborUnits(Ebar_tr(time-1), wage_tr(time-1)) 
		! Compute Capital demand and Profits by (a,z)
			K_mat  = K_Matrix(R_tr(time-1),P_tr(time-1))
			Pr_mat = Profit_Matrix(R_tr(time-1),P_tr(time-1))
		! Compute wealth given current R and P
			Wealth_mat = Wealth_Matrix_t(R_tr(time),P_tr(time))
		! Form YGRID for the capital income economy given interest rate "P"
			CALL FORM_Y_MB_GRID_Transition(YGRID,MBGRID,YGRID_t,MBGRID_t,time-1)
			! print*, 'YGRID_t(:,5,1)=',YGRID_t(:,5,1)
			! print*, 'RetY_lambda_e=',RetY_lambda_e
			! print*,"Period",ti,"|YGRID-YGRID_bench|=",maxval(abs(YGRID_t-YGRID_aux)),&
			! 	& "|RetY-RetY_bench|=",maxval(abs(RetY_lambda_e-RetY_lambda_e_aux))
			! print*,"Period",ti,"max(YGRID)=",maxval(YGRID_t),"max(YGRID_bench)=",maxval(YGRID_aux),&
			! 	& "max(RetY)=",maxval(RetY_lambda_e),"max(RetY_bench)=",minval(RetY_lambda_e_aux)

		! print*, 'R=',R,'P=',P, 'W=',wage, 'na=', na, 'na_t=', na_t
	!========================================================================================
	!------RETIREMENT PERIOD-----------------------------------------------------------------

	! SET HOURS TO ZERO SO THAT RETIRED HAS ZERO HOURS
	Hours_aux = 0.0_DP

	! Last period of life
	age=MaxAge
		! print*,' 	Age=',age
	!$omp parallel do private(lambdai,ei,ai,xi)
	DO zi=1,nz
	!$omp parallel do private(lambdai,ei,ai)
	DO xi=1,nx
    DO lambdai=1,nlambda
    DO ei=1,ne
    DO ai=1,na_t
        if ((YGRID_t(ai,zi,xi) + RetY_lambda_e(lambdai,ei)).le.(bq_0/chi_aux) ) then 
        Cons_aux(age,ai,zi,lambdai,ei,xi)   =  YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei) 
        else
        Cons_aux(age,ai,zi,lambdai,ei,xi)   = (YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei)+bq_0)/(1.0_dp+chi_aux)
        endif
        Ap_t_aux(age,ai,zi,lambdai,ei,xi) = YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei)-Cons_t(age,ai,zi,lambdai,ei,xi)
	ENDDO ! ai
    ENDDO ! ei
	ENDDO ! lambdai
	ENDDO ! xi
	ENDDO ! zi
	! print*,"Period",ti,"|Const_t-Const_bench|=",maxval(abs(Cons_t(MaxAge,:,:,:,:,:)-Cons_t_pr(MaxAge,:,:,:,:,:)))
	! print*,"Period",ti,"Test-Const_bench=",(YGRID_aux(50,5,2) + RetY_lambda_e_aux(3,3))/(1.0_DP+tauC)-Cons_bench(MaxAge,50,5,3,3,2)
	
	! Rest of retirement
	DO age=MaxAge-1,RetAge,-1
		! print*,' 	Age=',age
	!$omp parallel do private(lambdai,ei,ai,xi,xp_ind,EndoCons,EndoYgrid,sw,sort_ind,tempai, &
	!$omp& state_FOC,par_FOC,MB_aprime_t,EndoYgrid_sort)
    DO zi=1,nz
    !$omp parallel do private(lambdai,ei,ai,xp_ind,EndoCons,EndoYgrid,sw,sort_ind,tempai, &
	!$omp& state_FOC,par_FOC,MB_aprime_t,EndoYgrid_sort)
    DO xi=1,nx
    DO lambdai=1,nlambda
    DO ei=1,ne
    	! Endogenous grid and consumption are initialized
	    EndoCons  = big_p 
	    EndoYgrid = big_p 
	    sw 		  = 0
    DO ai=1,na_t 
		if (any((pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold)).lt.1e-8)) then 
			! print*, ' Threshold section - EGM Retirement'
			sw 			  = sw+1	
    		MB_aprime_t   = MBGRID_t(ai,zi,:)
    		! Consumption on endogenous grid and implied asset income under tauW_bt
    		do xp_ind = 1,nx 	
    			if (pr_x(xi,xp_ind,zi,age)/pr_x(xi,xp_ind,zi,age)*abs(Wealth_mat(ai,zi,xp_ind)-Y_a_threshold).lt.1e-8) then 
    				P = P_tr(time); R = R_tr(time) 
    				MB_aprime_t(xp_ind) = MB_a_bt(agrid_t(ai),zi,xp_ind)
    			endif
    		enddo 
    		EndoCons(ai)  =  (beta*survP(age)* 	&
    			& sum(pr_x(xi,:,zi,age)*MB_aprime_t*Cons_t_pr(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power)) &
    			& + (1.0_dp-survP(age))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
		    	& *   chi_bq*((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ) **euler_power
	        EndoYgrid(ai) = agrid_t(ai) +  EndoCons(ai) - RetY_lambda_e(lambdai,ei)
	        ! Consumption on endogenous grid and implied asset income under tauW_at
	        do xp_ind = 1,nx 	
    			if (pr_x(xi,xp_ind,zi,age)/pr_x(xi,xp_ind,zi,age)*abs(Wealth_mat(ai,zi,xp_ind)-Y_a_threshold).lt.1e-8) then 
    				P = P_tr(time); R = R_tr(time) 
    				MB_aprime_t(xp_ind) = MB_a_at(agrid_t(ai),zi,xp_ind)
    			endif
    		enddo 
	        EndoCons(na_t+sw)  = (beta*survP(age)*	&
	        	& sum(pr_x(xi,:,zi,age)*MB_aprime_t*Cons_t_pr(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power)) &
    			& + (1.0_dp-survP(age))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
		    	& *   chi_bq*((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ) **euler_power
	    	EndoYgrid(na_t+sw) = agrid_t(ai) +  EndoCons(na_t+sw) - RetY_lambda_e(lambdai,ei)

	  !   	print*, ' '
			! print*, ' Threshold test - Retirement'
			! print*, ' Current State', age, ai, zi, lambdai, ei, xi
			! print*, ' ', (pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold))
			! print*, ' ', (any((pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold)).lt.1e-8))
			! print*, ' ', MBGRID_t(ai,zi,:)
			! print*, ' ', MB_a_bt(agrid(ai),zi,xi)
			! print*, ' ', MB_a_at(agrid(ai),zi,xi)
			! print*, ' ', EndoCons(ai)
			! print*, ' ', EndoCons(na_t+sw)
			! print*, ' '
	    else 
	    	! Consumption on endogenous grid and implied asset income
	    	EndoCons(ai)  = (beta*survP(age)* 	&
	    		& sum(pr_x(xi,:,zi,age)*MBGRID_t(ai,zi,:)*Cons_t_pr(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power))&
    			& + (1.0_dp-survP(age))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
		    	& *   chi_bq*((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ) **euler_power
	        EndoYgrid(ai) = agrid_t(ai) +  EndoCons(ai) - RetY_lambda_e(lambdai,ei)
	        ! !$omp critical
	        ! print*,' Standard EGM - State:',age,ai,zi,lambdai,ei,xi,ti
	        ! print*,' 	EndoCons(ai)=',EndoCons(ai)
	        ! print*,' 	Expected value=',&
	        ! 	& sum(pr_x(xi,:,zi,age)*MBGRID_t(ai,zi,:)*Cons_t_pr(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power))
	        ! print*,' 	Cons_t+1=',Cons_t_pr(age+1,ai,zi,lambdai,ei,:)
	        ! print*,' 	Cons_t+1=',Cons_exp(age+1,ai,zi,lambdai,ei,:)
	        ! !$omp end critical
	    end if 


	 !    if (any(isnan(EndoCons))) then 
	 !    	print*,' '
		! 	print*,' '
		! 	print*,' '
		! 	print*,' Endo Consumption in retirement', age,ai,zi,lambdai,ei,xi
		! 	print*,' ',EndoCons(ai),(beta*survP(age)* 	&
	 !    				& sum(pr_x(xi,:,zi,age)*MBGRID_t(ai,zi,:)*Cons_t(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power)) ) **euler_power
		! 	print*,' MBGRID_t=',MBGRID_t(ai,zi,:)
		! 	print*,' cons(t+1)=',Cons_t(age+1,ai,zi,lambdai,ei,:)
		! 	print*,' Size(endoCons)=',size(EndoCons)
		! 	print*,' sw=',sw,'na_t=',na_t
		! 	print*,' Threshold_test=',(any((pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold)).lt.1e-8))
		! 	print*,' '
		! 	print*,' ',EndoCons
		! 	print*,' '
		! 	STOP
		! endif 
	ENDDO ! ai

	    !$omp critical
		if (any(isnan(EndoCons))) then 
			print*, "isnan - Consumption endogenous"
			print*, age,zi,ai,lambdai,ei,xi,time
			print*, EndoCons
			STOP 
		end if 
		!$omp end critical

		
	
	! Sort endogenous grid for interpolation
	call Sort(na_t+nx*nz+1,EndoYgrid,EndoYgrid_sort,sort_ind)
	! print*,' Yendo'
	! print*, EndoYgrid
	! print*,' Yendo_sort'
	! print*, EndoYgrid_sort
	! print*,' indices'
	! print*, sort_ind

	EndoYgrid = EndoYgrid_sort
	EndoCons = EndoCons(sort_ind)

	! print*, ' '
	! print*, ' isnan(endocons)', any(isnan(EndoCons))

	! Find  decision rules on exogenous grids
		! decision rules are obtained taking care of extrapolations
		tempai=1           
		DO WHILE ( YGRID_t(tempai,zi,xi) .lt. EndoYgrid(1) )
            tempai = tempai + 1
        ENDDO
	                
		DO ai=tempai,na_t              
		    ! CONSUMPTION ON EXOGENOUS GRIDS 
		    Cons_aux(age,ai,zi,lambdai, ei,xi)  = Linear_Int(EndoYgrid(1:na_t+sw), EndoCons(1:na_t+sw),na_t+sw, YGRID_t(ai,zi,xi))
		    Ap_t_aux(age,ai,zi,lambdai,ei,xi)   = &
		    								& YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei)-Cons_t(age, ai, zi, lambdai, ei,xi)
		    
		    If (Ap_t_aux(age,ai,zi,lambdai,ei,xi).lt.amin) then
		    	Ap_t_aux(age,ai,zi,lambdai,ei,xi) = amin
	            Cons_aux(age,ai,zi,lambdai,ei,xi) = &
	            	& YGRID_t(ai,zi,xi) + RetY_lambda_e(lambdai,ei) - Aprime_t(age,ai,zi,lambdai,ei,xi)
				IF (Cons_aux(age,ai,zi,lambdai,ei,xi) .le. 0.0_DP)  THEN
				    print*,'r1: Cons(age, ai, zi, lambdai,ei)=',Cons_aux(age, ai, zi, lambdai, ei,xi)
				ENDIF                   
	        endif 

	        ! ! !$omp critical 
	        ! if (isnan(Cons_t(age,ai,zi,lambdai,ei,xi))) then
	        ! 	print*,' '
	        ! 	print*,' '
	        ! 	print*,' '
	        ! 	print*,' isnan Consumption', age,ai,zi,lambdai,ei,xi
	        ! 	print*,' sw=',sw 
	        ! 	print*,' Yendo'
	        ! 	print*, EndoYgrid(1:na_t+sw)
	        ! 	print*,' Cendo'
	        ! 	print*, EndoCons(1:na_t+sw)
	        ! 	print*,' YGRID'
	        ! 	print*, YGRID_t(ai,zi,xi)
	        ! 	print*,' Linear Interpolation',Linear_Int(EndoYgrid(1:na_t+sw), EndoCons(1:na_t+sw),na_t+sw, YGRID_t(ai,zi,xi))
	        ! 	print*,' Aprime=',Aprime_t(age,ai,zi,lambdai,ei,xi)
	        ! 	call Sort(na_t+nx*nz+1,EndoYgrid,EndoYgrid_sort,sort_ind)
	        ! 	print*,' Yendo'
	        ! 	print*, EndoYgrid_sort
	        ! 	print*,' indices'
	        ! 	print*, sort_ind
	        ! 	print*, ' '
	        ! 	print*, ' ',minval(EndoYgrid),maxval(EndoYgrid)
	        ! 	print*, ' ',minval(EndoYgrid_sort),maxval(EndoYgrid_sort)
	        ! 	print*, ' The end'
	        ! 	STOP
	        ! endif 
	        ! ! !$omp end critical
		ENDDO ! ai  

        ai=1           
        DO WHILE ( YGRID_t(ai,zi,xi) .lt. EndoYgrid(1) )
			! Solve for a' directly by solving the Euler equation for retirement FOC_R
			state_FOC  = (/age,ai,zi,lambdai,ei,xi,time-1/)
			brentvalue = brent_p(min(amin,YGRID_t(ai,zi,xi)), (amin+YGRID_t(ai,zi,xi))/2.0_DP , &
			                & YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei) *0.95_DP, &
			                & FOC_R_Transition, brent_tol, Ap_t_aux(age,ai,zi,lambdai,ei,xi),state_FOC)
			
			Cons_aux(age,ai,zi,lambdai,ei,xi) = &
					& YGRID_t(ai,zi,xi) + RetY_lambda_e(lambdai,ei) - Ap_t_aux(age,ai,zi,lambdai,ei,xi)
			
			IF (Cons_aux(age, ai, zi, lambdai, ei,xi) .le. 0.0_DP)  THEN
			    print*,'r2: Cons(age,ai,zi,lambdai,ei,xi,ti)=',Cons_aux(age,ai,zi,lambdai,ei,xi)
			    print*,'Aprime(age,ai,zi,lambdai,ei,xi,ti)=',Ap_t_aux(age,ai,zi,lambdai,ei,xi)
			    print*,'YGRID(ai,zi,xi)+RetY_lambda_e(lambdai,ei)=',YGRID_t(ai,zi,xi)+RetY_lambda_e(lambdai,ei)
			    print*,'YGRID(ai,zi,xi)=',YGRID(ai,zi,xi),'EndoYgrid(1)=',EndoYgrid(1)
			    print*,'RetY_lambda_e(lambdai,ei)=',RetY_lambda_e(lambdai,ei)
			    print*,'lambdai=',lambdai
			ENDIF                   
			ai = ai + 1
		ENDDO  
	              
    ENDDO ! ei     
    ENDDO ! lambda
    ENDDO ! xi 
	ENDDO ! zi
    ENDDO !age

	!------RETIREMENT PERIOD ENDS------------------------------------------------------------
	!========================================================================================
	
	!========================================================================================
	!------Working Period Starts-------------------------------------------------------------

	DO age=RetAge-1,1,-1
		! print*,' 	Age=',age
	!$omp parallel do private(lambdai,ei,ai,xi,xp_ind,EndoCons,EndoHours,EndoYgrid,sw,sort_ind,tempai, &
	!$omp& C_foc,state_FOC,par_FOC,MB_aprime_t)
    DO zi=1,nz
    !$omp parallel do private(lambdai,ei,ai,xp_ind,EndoCons,EndoHours,EndoYgrid,sw,sort_ind,tempai, &
	!$omp& C_foc,state_FOC,par_FOC,MB_aprime_t)
    DO xi=1,nx
    DO lambdai=1,nlambda
    DO ei=1,ne	
    	! Endogenous grid and consumption are initialized
	    EndoCons  = big_p 
	    EndoYgrid = big_p 
	    EndoHours = big_p 
	    sw 		  = 0                    
        DO ai=1,na_t
        state_FOC  = (/age,ai,zi,lambdai,ei,xi,time-1/)
		if (any(pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold).lt.1e-8)) then 
			sw 			  = sw+1	
    		MB_aprime_t   = MBGRID_t(ai,zi,:)
	    	! Below threshold
    		do xp_ind = 1,nx 	
    			if (pr_x(xi,xp_ind,zi,age)/pr_x(xi,xp_ind,zi,age)*abs(Wealth_mat(ai,zi,xp_ind)-Y_a_threshold).lt.1e-8) then 
    				P = P_tr(time); R = R_tr(time) 
    				MB_aprime_t(xp_ind) = MB_a_bt(agrid_t(ai),zi,xp_ind)
    			endif
    		enddo 
			call EGM_Working_Period_Transition( MB_aprime_t, H_min, state_FOC, EndoCons(ai), EndoHours(ai), EndoYgrid(ai) )
			
			! Above threshold
			do xp_ind = 1,nx 	
    			if (pr_x(xi,xp_ind,zi,age)/pr_x(xi,xp_ind,zi,age)*abs(Wealth_mat(ai,zi,xp_ind)-Y_a_threshold).lt.1e-8) then 
    				P = P_tr(time); R = R_tr(time) 
    				MB_aprime_t(xp_ind) = MB_a_at(agrid_t(ai),zi,xp_ind)
    			endif
    		enddo 
			call EGM_Working_Period_Transition( MB_aprime_t, H_min, state_FOC, EndoCons(na_t+sw), EndoHours(na_t+sw), EndoYgrid(na_t+sw) )

	    	!print*, ' '
	    	!print*, State_FOC 
	    	!print*, MB_a_bt(agrid_t(ai),zgrid(zi)), MB_a_at(agrid_t(ai),zgrid(zi))
	  !   	print*, ' '
			! print*, ' Threshold test - Working Period'
			! print*, ' Current State', age, ai, zi, lambdai, ei, xi
			! print*, ' ', (pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold))
			! print*, ' ', (any((pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold)).lt.1e-8))
			! print*, ' ', MBGRID_t(ai,zi,:)
			! print*, ' ', MB_a_bt(agrid(ai),zi,xi)
			! print*, ' ', MB_a_at(agrid(ai),zi,xi)
			! print*, ' ', MB_aprime_t
			! print*, ' ', EndoCons(ai)
			! print*, ' ', EndoCons(na_t+sw)
			! print*, ' '
		else 
			! Usual EGM
			call EGM_Working_Period_Transition( MBGRID_t(ai,zi,:) , H_min , state_FOC , & 
			      & EndoCons(ai), EndoHours(ai) , EndoYgrid(ai)  )
	
		end if 

    	ENDDO ! ai

	    ! !$omp critical
		if (any(isnan(EndoCons))) then 
			print*, "isnan - Consumption endogenous"
			print*, age,zi,ai,lambdai,ei,xi,time
			print*, EndoCons
			STOP 
		end if 
		! !$omp end critical

    ! Sort endogenous grid for interpolation
	call Sort(na_t+nx*nz+1,EndoYgrid,EndoYgrid,sort_ind)
	EndoHours = EndoHours(sort_ind)
	EndoCons  = EndoCons(sort_ind)
	! 	print*, ' '
	! 	do ai=1,na_t+1
	! 		print*, EndoHours(ai), EndoCons(ai), EndoYgrid(ai)
	! 	end do
	! 	print*, ' '

		if (any(isnan(EndoCons)).or.any(isnan(EndoHours)).or.any(isnan(EndoYgrid))) then 
			print*, "isnan - Consumption working 4"
			print*, age,lambdai,ai,zi,ei,xi,time
			STOP 
		end if 

    	! Find  decision rules on exogenous grids
        tempai=1           
        DO WHILE ( YGRID_t(tempai,zi,xi) .lt. EndoYgrid(1) )
              tempai = tempai + 1
        ENDDO
	                
		! decision rules are obtained taking care of extrapolations
		DO ai=tempai,na_t 
			! Interpolate for value of consumption in exogenous grid
			Cons_aux(age,ai,zi,lambdai,ei,xi) = Linear_Int(EndoYgrid(1:na_t+sw),EndoCons(1:na_t+sw),na_t+sw,YGRID_t(ai,zi,xi))

			if (Progressive_Tax_Switch.eqv..true.) then
			! Check FOC for implied consumption 
				if (NSU_Switch.eqv..true.) then 
					! Non-Separable Utility
					C_foc = (gamma/(1.0_dp-gamma))*(1.0_dp-H_min)*MB_h(H_min,age,lambdai,ei,wage_tr(time-1))
				else
					! Separable Utility
					C_foc = (MB_h(H_min,age,lambdai,ei,wage_tr(time-1))*(1.0_dp-H_min)**(gamma)/phi)**(1.0_dp/sigma)
				end if 
			! Hours
				if (Cons_t(age,ai,zi,lambdai,ei,xi).ge.C_foc) then
					Hours_aux(age,ai,zi,lambdai,ei,xi) = 0.0_dp
				else
					! Auxiliary variables for solving FOC for hours 
					par_FOC(1:6) = (/age,ai,zi,lambdai,ei,xi/) 
					par_FOC(7)   = Cons_aux(age,ai,zi,lambdai,ei,xi)
					! FOC (set wage to current wage)
					wage = wage_tr(time-1)
					brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_H, brent_tol,Hours_aux(age,ai,zi,lambdai,ei,xi) , par_FOC) 
				end if
			else 
				if (NSU_Switch.eqv..true.) then 
					Hours_aux(age,ai,zi,lambdai,ei,xi) = max( 0.0_dp , &
						&  1.0_DP - (1.0_DP-gamma)*Cons_aux(age,ai,zi,lambdai,ei,xi)/(gamma*psi*yh(age,lambdai,ei)) )
				else 
					Hours_aux(age,ai,zi,lambdai,ei,xi) = max( 0.0_DP , &
						&  1.0_DP - phi*Cons_aux(age,ai,zi,lambdai,ei,xi)/(psi*yh(age, lambdai,ei)) )
				end if 
			end if 

			if (Hours_aux(age,ai,zi,lambdai,ei,xi).ge.1.0_dp) then 
				print*, ' '
				print*, 'Hours highers than 1', age,ai,zi,lambdai,ei,xi,time
				print*, Cons_aux(age,ai,zi,lambdai,ei,xi)
				STOP
			endif 


			! Savings 
				Ap_t_aux(age,ai,zi,lambdai,ei,xi) = & 
								& YGRID_t(ai,zi,xi)  + Y_h(Hours_aux(age,ai,zi,lambdai,ei,xi),age,lambdai,ei,wage_tr(time-1))  & 
            					& - Cons_aux(age,ai,zi,lambdai,ei,xi) 
		                    
		    If (Ap_t_aux(age,ai,zi,lambdai,ei,xi)  .lt. amin) then

		    	! print*, ' Aprime was below minimum!!!!
            	Ap_t_aux(age,ai,zi,lambdai,ei,xi) = amin
		         
	           	! Compute hours
		        if (NSU_Switch.and.(Progressive_Tax_Switch.eqv..false.)) then 
		        	Hours_aux(age,ai,zi,lambdai,ei,xi)  = max( 0.0_dp , &
		        	&  gamma - (1.0_DP-gamma)*(YGRID_t(ai,zi,xi)-Ap_t_aux(age,ai,zi,lambdai,ei,xi))/(psi*yh(age,lambdai,ei)))
                else               	        
                	!compute  hours using FOC_HA                              
		        	par_FOC(1:6) = (/age,ai,zi,lambdai,ei,xi/) 
		        	par_FOC(7)   = amin
					wage = wage_tr(time-1)
		        	brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_HA, brent_tol, Hours_aux(age,ai,zi,lambdai,ei,xi),par_FOC)
		        end if  

	            Cons_aux(age,ai,zi,lambdai,ei,xi) = &
	            					& YGRID_t(ai,zi,xi)+  Y_h(Hours_aux(age,ai,zi,lambdai,ei,xi),age,lambdai,ei,wage_tr(time-1))  &
		                            & - Ap_t_aux(age,ai,zi,lambdai,ei,xi)

				IF (Cons_aux(age,ai,zi,lambdai,ei,xi) .le. 0.0_DP)  THEN
					print*,'w1: Cons(age,ai,zi,lambdai,ei,xi,ti)=',Cons_aux(age,ai,zi,lambdai,ei,xi)
					STOP
				ENDIF 
				if (Hours_aux(age,ai,zi,lambdai,ei,xi).ge.1.0_dp) then 
					print*, ' '
					print*, 'w1 Hours highers than 1', age,ai,zi,lambdai,ei,xi,time
					print*, Cons_aux(age,ai,zi,lambdai,ei,xi)
					STOP
				endif                   
		     endif  
		    !$omp critical
		    IF (Cons_aux(age,ai,zi,lambdai,ei,xi) .le. 0.0_DP)  THEN
				print*,'w1: Cons(age,ai,zi,lambdai,ei,xi)=',Cons_aux(age,ai,zi,lambdai,ei,xi)
				STOP
			ENDIF 
			if (Hours_aux(age,ai,zi,lambdai,ei,xi).ge.1.0_dp) then 
				print*, ' '
				print*, 'w1 Hours highers than 1', age,ai,zi,lambdai,ei,xi,time
				print*, Cons_aux(age,ai,zi,lambdai,ei,xi)
				STOP
			endif    
			!$omp end critical
		ENDDO ! ai   

		!$omp critical
		if (any(isnan(Cons_aux(age,:,zi,lambdai,ei,xi)))) then 
			print*, "isnan - Consumption working 3"
			print*, age,lambdai,ai,zi,ei,xi,time
			print*, 'Cons'
			print*, Cons_aux(age,:,zi,lambdai,ei,xi)
			print*, 'Hours'
			print*, Hours_aux(age,:,zi,lambdai,ei,xi) 
			STOP 
		end if         
		!$omp end critical        

		ai=1           
        DO WHILE ( YGRID_t(ai,zi,xi) .lt. EndoYgrid(1) )
        	! print*, ' Extrapolation between YGRID and EndoYgrid!!!!'
	        ! Solve for the Euler equation directly
	        state_FOC  = (/age,ai,zi,lambdai,ei,xi,time-1/)
			brentvalue = brent_p( min(amin,YGRID_t(ai,zi,xi))   ,  (amin+YGRID_t(ai,zi,xi))/2.0_DP  ,  &
                             	& YGRID_t(ai,zi,xi)  + psi * ( yh(age, lambdai,ei)*0.95_DP )**(1.0_DP-tauPL)  ,  &
	                             & FOC_WH_Transition, brent_tol, Ap_t_aux(age,ai,zi,lambdai,ei,xi) , state_FOC )

			! Compute hours
	        if (NSU_Switch.and.(Progressive_Tax_Switch.eqv..false.)) then 
	        	Hours_aux(age,ai,zi,lambdai,ei,xi)  = max( 0.0_dp , &
	        	&  gamma - (1.0_DP-gamma) *( YGRID_t(ai,zi,xi) - Ap_t_aux(age,ai,zi,lambdai,ei,xi))/(psi*yh(age,lambdai,ei)))
            else               	        
				!compute  hours using FOC_HA                              
				par_FOC(1:6) = (/age,ai,zi,lambdai,ei,xi/)
				par_FOC(7)   = Ap_t_aux(age,ai,zi,lambdai,ei,xi)
				wage = wage_tr(time-1)
	            brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_HA, brent_tol, Hours_aux(age,ai,zi,lambdai,ei,xi), par_FOC) 
 			end if 

            Cons_aux(age,ai,zi,lambdai,ei,xi)=  &
            		& YGRID_t(ai,zi,xi) + Y_h(Hours_aux(age,ai,zi,lambdai,ei,xi),age,lambdai,ei,wage_tr(time-1))  &
					& - Ap_t_aux(age,ai,zi,lambdai,ei,xi)

           IF (Cons_aux(age,ai,zi,lambdai,ei,xi) .le. 0.0_DP)  THEN
                print*,'w2:',age,zi,lambdai,ei,ai,xi,time, 'Cons=',Cons_aux(age,ai,zi,lambdai,ei,xi), &
                    & 'Aprime=',Ap_t_aux(age,ai,zi,lambdai,ei,xi), &
                    & 'Hours=', Hours_aux(age,ai,zi,lambdai,ei,xi), &
                    & 'yh(age, lambdai,ei)=', yh(age, lambdai,ei),'YGRID(ai,zi)=',YGRID_t(ai,zi,xi)
                !pause
                print*, "there is a problem in EGM Transition with AGG_Debt"
                STOP
            ENDIF 

            if (Hours_aux(age,ai,zi,lambdai,ei,xi).ge.1.0_dp) then 
				print*, ' '
				print*, 'w2 Hours highers than 1', age,ai,zi,lambdai,ei,xi,time
				print*, Cons_aux(age,ai,zi,lambdai,ei,xi)
				STOP
			endif 

            ai = ai + 1
        ENDDO  

	    if (any(isnan(Cons_aux(age,:,zi,lambdai,ei,xi)))) then 
			print*, "isnan - Consumption working 2"
			print*, age,lambdai,ai,zi,ei,xi,time
			print*, 'Cons'
			print*, Cons_aux(age,:,zi,lambdai,ei,xi)
			print*, 'Hours'
			print*, Hours_aux(age,:,zi,lambdai,ei,xi) 
			STOP 
		end if 

	                 
    ENDDO !ei         
    ENDDO !lambdai
    ENDDO ! xi 
	ENDDO !zi
    ENDDO !age


	! Interpolate to get values of policy functions on agrid (note that policy functions are defined on agrid_t)
	
	if (Y_a_threshold.eq.0.0_dp) then 
		! print*,' Assigning value of policy functions'
		Ap_aux = Ap_t_aux
	else 
		! print*,' Interpolating policy functions'
		DO xi=1,nx
		DO age=1,MaxAge
	    DO lambdai=1,nlambda
	    DO zi=1,nz
	    DO ei=1,ne	                
	    DO ai=1,na
			Ap_aux(age,ai,zi,lambdai,ei,xi) = &
	    		& Linear_Int(YGRID_t(:,zi,xi) , Ap_t_aux(age,:,zi,lambdai,ei,xi) , na_t , YGRID(ai,zi,xi))
		ENDDO !ai         
		ENDDO !ei         
	    ENDDO !zi
	    ENDDO !lambdai
		ENDDO !age
		ENDDO !xi
	end if 

	
	! print*,' 	----------------------------------------'
	! print*,' 	End of EGM Transition inside Agg_Debt'
	! print*,' 	----------------------------------------'
	! print*,' '

END SUBROUTINE EGM_Transition_Aux

!========================================================================================
!========================================================================================
!========================================================================================
! EGM_Working_Period: Solves for the endogenous grid during the working period
!
! Usage: call EGM_Working_Period(MB_in,H_min,C_endo,H_endo,Y_endo)
!
! Input: MB_in   , real(dp), marginal benefit from assets in t+1
!        H_min   , real(dp), minimum value of hours
!
! Implicit inputs: ai     , integer , index of current level of assets
!				   zi     , integer , index of agent's permanent entreprenurial ability
!				   lambdai, integer , index of agent's permanent labor income component
!				   ei     , integer , index of agent's transitory labor income component
!				   age    , integer , agent's age
!
! Output: C_endo , real(dp), Endogenous value of consumption in t
! 		  H_endo , real(dp), Endogenous value of hours in t
! 		  Y_endo , real(dp), Endogenous value of Y grid in t
!
Subroutine EGM_Working_Period_Transition(MB_in,H_min,state_FOC,C_endo,H_endo,Y_endo)
	Implicit None 
	real(dp), intent(in)   	  :: MB_in(nx), H_min, state_FOC(7)
	real(dp), intent(out)  	  :: C_endo, H_endo, Y_endo
	real(dp)               	  :: C_euler, C_foc, brentvalue, E_MU_cp(nx), MB_a_vec(nx)
	REAL(DP), DIMENSION(7+nx) :: par_FOC
	integer                	  :: age, ai, zi, lambdai, ei, xi, ti, xp_ind

	! Allocate state and get indeces
	age 	= int(state_FOC(1))
	ai   	= int(state_FOC(2))
	zi      = int(state_FOC(3))
	lambdai = int(state_FOC(4))
	ei      = int(state_FOC(5))
	xi      = int(state_FOC(6))
	ti      = int(state_FOC(7))

	! Set current wage
		wage = wage_tr(ti)

	if (Progressive_Tax_Switch.eqv..true.) then !Progressive labor taxes 
		if (NSU_Switch.eqv..true.) then 
			! Non-Separable Utility
			do xp_ind=1,nx
				E_MU_cp(xp_ind) = SUM(pr_e(ei,:) * &
			        & Cons_t_pr(age+1,ai,zi,lambdai,:,xp_ind)**((1.0_dp-sigma)*gamma-1.0_dp)    * &
			        & (1.0_dp-Hours_t_pr(age+1,ai,zi,lambdai,:,xp_ind))**((1.0_dp-sigma)*(1.0_dp-gamma))  )
			enddo

			C_euler = ( (beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_MU_cp)) &
	        			& + (1.0_dp-survP(age))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
				    	& *   chi_bq*((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ) &
						& **(1.0_dp/((1.0_dp-sigma)*gamma-1.0_dp))
			C_foc   = (gamma/(1.0_dp-gamma))*(1.0_dp-H_min)*MB_h(H_min,age,lambdai,ei,wage)

			if (C_euler.ge.C_foc) then
				C_endo  = C_euler 
				H_endo = 0.0_dp
			else 
				if (Log_Switch.eqv..true.) then
		  			! Auxiliary consumption variable for FOC_H 
		  			par_FOC(1:6) = state_FOC(1:6)       
				    par_FOC(7) = C_euler
				    ! Solution for hours from Euler or FOC  equation according to sigma 
				    brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_H, brent_tol, H_endo , par_FOC ) 
				    C_endo = C_euler
				else 
				    ! Set Marginal benefit of assets to the below threshold level
				    par_FOC(1:7) = state_FOC
				    par_FOC(8:7+nx) = MB_in
				    ! Solution for hours from Euler or FOC  equation according to sigma 
				    brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_H_NSU_Transition, brent_tol, H_endo , par_FOC ) 
				    ! Implied consumption by hours from Labor FOC
				    C_endo = (gamma/(1.0_dp-gamma))*(1.0_dp-H_endo)*MB_h(H_endo,age,lambdai,ei,wage)
				end if 
			end if 

		else 
			! Separable Utility
			do xp_ind=1,nx
				E_MU_cp(xp_ind) = SUM(pr_e(ei,:) * Cons_t_pr(age+1,ai,zi,lambdai,:,xp_ind)**(-sigma) )
			enddo
			C_endo = 1.0_dp/( (beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_mu_cp)  &
	        			& + (1.0_dp-survP(age))*(1.0_dp+tauC)**(1.0_dp-sigma)&
				    	& *   chi_bq*((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**(-sigma) ) **(1.0_dp/sigma) )
			C_foc  = (MB_h(H_min,age,lambdai,ei,wage)*(1.0_dp-H_min)**(gamma)/phi)**(1.0_dp/sigma)

			if (C_endo.ge.C_foc) then
			  H_endo = 0.0_dp
			else 
			  ! Auxiliary consumption variable for FOC_H 
			  par_FOC(1:6) = state_FOC(1:6)              
			  par_FOC(7) = C_endo
			  ! Solution for hours from Euler or FOC  equation according to sigma 
			  brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_H, brent_tol, H_endo , par_FOC ) 
			end if  

			end if 

	else ! Linear labor taxes 
		if (NSU_Switch.eqv..true.) then 
			! Non-Separable Utility
			do xp_ind=1,nx
				E_MU_cp(xp_ind) = &
				& sum( pr_e(ei,:) * (Cons_t_pr(age+1,ai,zi,lambdai,:,xp_ind)**(gamma*(1.0_DP-sigma)-1.0_DP)) &
			    & *  ( (1.0_DP-Hours_t_pr(age+1, ai, zi, lambdai,:,xp_ind))**((1.0_DP-gamma)*(1.0_DP-sigma))))
			enddo
			  C_endo = ((gamma*psi*yh(age, lambdai,ei)/(1.0_DP-gamma))**((1.0_DP-gamma)*(1.0_DP-sigma)) &
			    & *  (beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_MU_cp)  &
	        		& + (1.0_dp-survP(age))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
					& *   chi_bq*((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) ) )**(-1.0_DP/sigma)

			  H_endo = 1.0_DP - (1.0_DP-gamma)*C_endo/(gamma*psi*yh(age, lambdai,ei))   

			If (H_endo .lt. 0.0_DP) then
			    H_endo = 0.0_DP 
			    C_endo = ( beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_MU_cp)  &
	        		& + (1.0_dp-survP(age))*(1.0_dp+tauC)**((1.0_dp-sigma)*gamma)&
					& *   chi_bq*((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**((1.0_dp-sigma)*gamma-1.0_dp) )&
			    	& **(1.0_DP/(gamma*(1.0_DP-sigma)-1.0_DP))
			endif 

			! print*,' '
			! print*,' 	Inside EGM'
			! print*,' 	Consumption_t+1=',Cons_t_pr(age+1,ai,zi,lambdai,:,:)
		else 
			! Separable Utility
			do xp_ind=1,nx
				E_MU_cp(xp_ind) = sum( pr_e(ei,:) * (Cons_t_pr(age+1,ai,zi,lambdai,:,xp_ind)**(-sigma)) )
			enddo
			C_endo  = 1.0_DP/( beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_MU_cp)  &
	        			& + (1.0_dp-survP(age))*(1.0_dp+tauC)**(1.0_dp-sigma)&
				    	& *   chi_bq*((1.0_dp-tau_bq)*agrid_t(ai)+bq_0)**(-sigma) )**(1.0_DP/sigma)

			H_endo = max(0.0_DP , 1.0_DP - (phi*C_endo**sigma/(psi*yh(age, lambdai,ei)))**(1.0_dp/gamma) )  

		end if 
	end if


	! Endogenous grid for asset income
	Y_endo = agrid_t(ai) + C_endo - Y_h(H_endo,age,lambdai,ei,wage)

	! 		!$omp critical
	! 		if ((zi.eq.1).and.(ai.eq.16)) then 
	! 		print*, "EGM Working Periods"
	! 		print*, C_endo,H_endo,Y_endo
	! 		print*, MB_in,state_FOC
	! 		print*, ' '
	! 		endif 
	! 		!$omp end critical
end Subroutine EGM_Working_Period_Transition


!========================================================================================
!========================================================================================
!========================================================================================

! Under construction!!!!!!
Subroutine Find_TauW_Threshold(DBN_in,Y_a_threshold_out)
	real(dp), intent(in)  :: DBN_in(MaxAge, na, nz, nlambda, ne,nx)
	real(dp), intent(out) :: Y_a_threshold_out
	real(dp), dimension(na,nz,nx) :: DBN_azx, Wealth
	real(dp)                   :: Mean_Wealth
	integer                    :: prctile, a_ind, z_ind

	if (bench_indx.eq.1) then 
		Y_a_threshold_out = 0
	else 
		! Compute distribution of agents by (a,z)
		DBN_azx = sum(sum(sum(DBN_in,5),4),1)
		! Compute mean before tax wealth
		Wealth = Wealth_Matrix(R,P)
		! Mean Wealth
		Mean_Wealth = sum( Wealth*DBN_azx )
		! Set threshold
		Y_a_threshold_out = Mean_Wealth
			!print*, "Current Threshold for wealth taxes", Y_a_threshold_out
	end if 
			! 	! Compute threshold for wealth tax
			! 	! Find distribution of assets
			! 		DO a_ind=1,na
			! 		    pr_az_dbn(a_ind)  = sum( DBN_in(:,a_ind,:,:,:) ) 
			! 		    cdf_az_dbn(a_ind) = sum( pr_a_dbn(1:a_ind)  )      
			! 		ENDDO
			! 		cdf_a_dbn = cdf_a_dbn + 1.0_DP - cdf_a_dbn(na)

			! 	! Find ai (grid node) that corresponds to percentile of wealth distribution
			! 		DO prctile=1,100
			! 		    a_ind=1
			! 		    DO while (cdf_a_dbn(a_ind) .lt. (REAL(prctile,8)/100.0_DP-0.000000000000001))
			! 		        a_ind=a_ind+1
			! 		    ENDDO
			! 		    prctile_ai_ind(prctile) = a_ind
			! 		ENDDO

			! 	! Set threshold for wealth tax as:
			! 		if (tauW_threshold_in.eq.0) then 
			! 			! Threshold is minimum
			! 			Y_a_threshold_out = 0.0_dp
			! 		else 
			! 			! Threshold at some percentile
			! 			Y_a_threshold_out = agrid( prctile_ai_ind(tauW_threshold_in) )
			! 		end if 

end Subroutine Find_TauW_Threshold



!========================================================================================
!========================================================================================
!========================================================================================

Function Labor_Market_Clearing(W_in)
	Implicit None 
	real(dp), intent(in) :: W_in
	real(dp)             :: Labor_Market_Clearing
	real(dp) 			 :: Ld_P, Ld_C 

	Ld_P  = ( (1.0_dp-alpha  )*Aprod/W_in )**(1.0_dp/alpha  ) * QBAR 
	Ld_C  = ( (1.0_dp-alpha_C)*A_C  /W_in )**(1.0_dp/alpha_C) * K_C

	Labor_Market_Clearing = ( ( Ld_P+Ld_C )/NBAR - 1.0_dp )**2.0_dp 


end Function Labor_Market_Clearing



!========================================================================================
!========================================================================================
!========================================================================================
! This function computes the capital demand matrix for agents of type (a,z)
! The function uses the interest rate R, price of intermediate good P

Function K_Matrix(R_in,P_in)
	Implicit None 
	real(dp), intent(in)          :: R_in, P_in
	real(dp), dimension(na,nz,nx) :: K_Matrix
	integer :: i,j,h

	! K_Matrix = min( theta*spread(agrid,2,nz) , (mu*P_in*spread(zgrid,1,na)**mu/(R_in+DepRate))**(1.0_dp/(1.0_dp-mu)) )

	do h=1,nx
	do j=1,nz
	do i=1,na 
		K_Matrix(i,j,h) = min( theta(j)*agrid(i) , (mu*P_in*xz_grid(h,j)**mu/(R_in+DepRate))**(1.0_dp/(1.0_dp-mu)) )
	enddo
	enddo
	enddo 

end Function K_Matrix

Function K_Matrix_t(R_in,P_in)
	Implicit None 
	real(dp), intent(in)            :: R_in, P_in
	real(dp), dimension(na_t,nz,nx) :: K_Matrix_t
	integer :: i,j,h

	! K_Matrix_t = min( theta*spread(agrid_t,2,nz) , (mu*P_in*spread(zgrid,1,na_t)**mu/(R_in+DepRate))**(1.0_dp/(1.0_dp-mu)) )

	do h=1,nx
	do j=1,nz
	do i=1,na_t
		K_Matrix_t(i,j,h) = min( theta(j)*agrid_t(i) , (mu*P_in*xz_grid(h,j)**mu/(R_in+DepRate))**(1.0_dp/(1.0_dp-mu)) )
	enddo
	enddo
	enddo 

end Function K_Matrix_t

!========================================================================================
!========================================================================================
!========================================================================================
! This function computes the profit matrix for agents of type (a,z)
! The function uses the interest rate R, price of intermediate good P

Function Profit_Matrix(R_in,P_in)
	Implicit None 
	real(dp), intent(in)          :: R_in, P_in
	real(dp), dimension(na,nz,nx) :: Profit_Matrix
	real(dp), dimension(na,nz,nx) :: K
	integer :: i,j,h

	K = K_matrix(R_in,P_in)

	! Profit_Matrix = P_in*(spread(zgrid,1,na)*K)**mu - (R_in+DepRate)*K

	do h=1,nx
	do j=1,nz 
	do i=1,na 
		Profit_Matrix(i,j,h) = P_in*(xz_grid(h,j)*K(i,j,h))**mu - (R_in+DepRate)*K(i,j,h)
	enddo 
	enddo
	enddo 

end Function Profit_Matrix

Function Profit_Matrix_t(R_in,P_in)
	Implicit None 
	real(dp), intent(in)            :: R_in, P_in
	real(dp), dimension(na_t,nz,nx) :: Profit_Matrix_t
	real(dp), dimension(na_t,nz,nx) :: K
	integer :: i,j,h

	K = K_matrix_t(R_in,P_in)

	! Profit_Matrix_t = P_in*(spread(zgrid,1,na_t)*K)**mu - (R_in+DepRate)*K
	do h=1,nx 
	do j=1,nz 
	do i=1,na_t 
		Profit_Matrix_t(i,j,h) = P_in*(xz_grid(h,j)*K(i,j,h))**mu - (R_in+DepRate)*K(i,j,h)
	enddo 
	enddo 
	enddo 

end Function Profit_Matrix_t

!========================================================================================
!========================================================================================
!========================================================================================
! This function computes the before wealth tax wealth matrix for agents of type (a,z) on agrid_t
! The function uses the interest rate R, price of intermediate good P

Function Wealth_Matrix(R_in,P_in)
	Implicit None 
	real(dp), intent(in)          :: R_in, P_in
	real(dp), dimension(na,nz,nx) :: Wealth_Matrix
	real(dp), dimension(na,nz,nx) :: Pr
	integer :: i,j,h

	Pr = Profit_Matrix(R_in,P_in)

	Wealth_Matrix = ( (1.0_dp+R_in*(1.0_DP-tauK))*spread(spread(agrid,2,nz),3,nx) +  Pr *(1.0_DP-tauK) )

end Function Wealth_Matrix

Function Wealth_Matrix_t(R_in,P_in)
	Implicit None 
	real(dp), intent(in)            :: R_in, P_in
	real(dp), dimension(na_t,nz,nx) :: Wealth_Matrix_t
	real(dp), dimension(na_t,nz,nx) :: Pr
	integer :: i,j,h

	Pr = Profit_Matrix_t(R_in,P_in)

	Wealth_Matrix_t = ( (1.0_dp+R_in*(1.0_DP-tauK))*spread(spread(agrid_t,2,nz),3,nx) +  Pr *(1.0_DP-tauK) )

end Function Wealth_Matrix_t


!========================================================================================
!========================================================================================
! This function computes the difference between capital demand and supply
! The function uses the interest rate R and implicitely the price of intermediate good P and distribution DBN1

Function Agg_Debt(R_in)
	Implicit None 
	real(dp), intent(in) :: R_in
	real(dp)             :: Agg_Debt
	real(dp), dimension(na,nz,nx) :: DBN_azx, K_mat
	real(dp)             :: Wealth , Kd

	DBN_azx  = sum(sum(sum(DBN1,5),4),1)

	Wealth   = sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid )

	K_mat = K_matrix(R_in,P)

	Agg_Debt = (sum(DBN_azx*(K_mat-spread(spread(agrid,2,nz),3,nx)))/Wealth)**2.0_DP 
	Kd       = (sum(DBN_azx*(K_mat)))

	
	Agg_Debt = 0.0_dp
	do xi=1,nx
	do zi=1,nz 
	do ai=1,na
		Agg_Debt = Agg_Debt + sum(DBN1(:,ai,zi,:,:,xi)*(K_mat(ai,zi,xi)-agrid(ai)))
	enddo 
	enddo 
	enddo 

	! Adjust with Government Debt
		Agg_Debt = Agg_Debt + Debt_Absorption*Debt_SS 

	! print*, mu, P, R_in, DepRate
	! print*, xz_grid(1,5:)
	! print*, (mu*P*xz_grid(1,5:)**mu/(R_in+DepRate))**(1.0_dp/(1.0_dp-mu))
	! print*, '------------',Wealth, Kd, Agg_Debt, R_in, P
	Agg_Debt = (Agg_Debt/Wealth)**2.0_dp


end Function Agg_Debt

!========================================================================================
!========================================================================================
! This function computes the difference between capital demand and supply
! The function uses the interest rate R and implicitely the price of intermediate good P and distribution DBN1

Function Agg_Debt_C(R_in)
	Implicit None 
	real(dp), intent(in) :: R_in
	real(dp)             :: Agg_Debt_C
	real(dp), dimension(na,nz,nx) :: DBN_azx, K_mat
	real(dp)             :: Wealth , Kd, K_Corp

	DBN_azx  = sum(sum(sum(DBN1,5),4),1)

	Wealth   = sum( sum(sum(sum(sum(sum(DBN1                ,6),5),4),3),1)*agrid )

	K_Corp   = sum( sum(sum(sum(sum(sum(DBN1(:,:,z_C:,:,:,:),6),5),4),3),1)*agrid ) 

	K_mat    = K_matrix(R_in,P)

	Kd       = (sum(DBN_azx*(K_mat)))

	! Debt as difference between private demand supply of assets
	Agg_Debt_C = Kd - (Wealth - K_Corp)

	! Adjust with Government Debt
		Agg_Debt_C = Agg_Debt_C + Debt_Absorption*Debt_SS 

	! print*, mu, P, R_in, DepRate
	! print*, xz_grid(1,5:)
	! print*, (mu*P*xz_grid(1,5:)**mu/(R_in+DepRate))**(1.0_dp/(1.0_dp-mu))
	! print*, '------------',Wealth, Kd, Agg_Debt, R_in, P
	Agg_Debt_C = (Agg_Debt_C/(Wealth-K_Corp))**2.0_dp

end Function Agg_Debt_C

!========================================================================================
!========================================================================================
! This function computes the difference between capital demand and supply
! The function is modified for the transition experiments
! Government debt is included
! The function uses the interest rate R and implicitely the price of intermediate good P and distribution DBN1
Function Agg_Debt_Tr(R_in)
	use omp_lib
	Implicit None 
	real(dp), intent(in) :: R_in
	real(dp)             :: Agg_Debt_Tr, Private_Demand, Public_Demand, Agg_Demand
	real(dp), dimension(na,nz,nx) :: K_mat
	real(dp)             :: Wealth
	INTEGER :: tklo, tkhi, age1, age2, z1, z2, a1, a2, lambda1, lambda2, e1, e2, x1, x2
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: Ap_aux, DBN_aux
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: PrAprimelo, PrAprimehi, PrBqlo, PrBqhi
	INTEGER , DIMENSION(:,:,:,:,:,:), allocatable :: Aplo, Aphi, Bqlo, Bqhi

	!$ call omp_set_num_threads(nz)

	allocate( Ap_aux(     MaxAge,na,nz,nlambda,ne,nx) )
	allocate( DBN_aux(    MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrAprimehi( MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrAprimelo( MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Aplo(       MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Aphi(       MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrBqlo(     MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrBqhi(     MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Bqlo(       MaxAge,na,nz,nlambda,ne,nx) )
	allocate( BQhi(       MaxAge,na,nz,nlambda,ne,nx) )

	! Assign new interest rate to R_tr
	R_tr(ti) = R_in 

	! Initialize new distribution
	DBN_aux = 0.0_dp

	! Compute auxiliary distribution of assets 
	if (ti.ge.2) then
		! Solve auxiliary DP problem 
		call EGM_Transition_aux(Ap_aux,ti)

		! Update t-1 distribution using Ap_aux


		! Discretize policy function for assets (a') for current period
		! For each age and state vector bracket optimal a' between two grid points
		! When at that age and state the optimal decision is approximated by selecting one the grid points
		! The grid points are selected with probability proportional to their distance to the optimal a'
		! print*,' Discretizing Policy Functions'
		DO age=1,MaxAge
			! print*,' 		Age ',age
		!$omp parallel do private(lambdai,ei,ai,xi,tklo,tkhi)
		DO zi=1,nz
		DO xi=1,nx
		DO ai=1,na
		DO lambdai=1,nlambda
		DO ei=1, ne
	        if ( Ap_aux(age,ai,zi,lambdai,ei,xi) .ge. amax) then
	            tklo =na-1
	        elseif (Ap_aux(age,ai,zi,lambdai,ei,xi) .lt. amin) then
	            tklo = 1
	        else
	            tklo = ((Ap_aux(age,ai,zi,lambdai,ei,xi) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
	        endif
	        tkhi = tklo + 1        
	        Aplo(age,ai,zi,lambdai,ei,xi)  	   = tklo
	        Aphi(age,ai,zi,lambdai,ei,xi)  	   = tkhi        
	        PrAprimelo(age,ai,zi,lambdai,ei,xi) = (agrid(tkhi) - Ap_aux(age,ai,zi,lambdai,ei,xi))/(agrid(tkhi)-agrid(tklo))
	        PrAprimehi(age,ai,zi,lambdai,ei,xi) = (Ap_aux(age,ai,zi,lambdai,ei,xi) - agrid(tklo))/(agrid(tkhi)-agrid(tklo))

        if ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Ap_aux(age,ai,zi,lambdai,ei,xi) .ge. amax) then
            tklo =na-1
        elseif ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Ap_aux(age,ai,zi,lambdai,ei,xi) .lt. amin) then
            tklo = 1
        else
            tklo = (((1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Ap_aux(age,ai,zi,lambdai,ei,xi)-amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
        endif
        tkhi = tklo + 1        
        Bqlo(age,ai,zi,lambdai,ei,xi)  	= tklo
        Bqhi(age,ai,zi,lambdai,ei,xi)  	= tkhi        
        PrBqlo(age,ai,zi,lambdai,ei,xi) = &
        	& ( agrid(tkhi) - (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Ap_aux(age,ai,zi,lambdai,ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
        PrBqhi(age,ai,zi,lambdai,ei,xi) = &
        	& ( (1.0_dp-tau_bq)*(1.0_dp-bq_fee)*Ap_aux(age,ai,zi,lambdai,ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		! print*, ' Discretizing Policy Functions Completed'

		! Probablities are adjusted to lie in [0,1]
		PrAprimelo = min(PrAprimelo, 1.0_DP); PrAprimelo = max(PrAprimelo, 0.0_DP)
		PrAprimehi = min(PrAprimehi, 1.0_DP); PrAprimehi = max(PrAprimehi, 0.0_DP)

		PrBqlo = min(PrBqlo, 1.0_DP); PrBqlo = max(PrBqlo, 0.0_DP)
		PrBqhi = min(PrBqhi, 1.0_DP); PrBqhi = max(PrBqhi, 0.0_DP)


		! Everyone in MaxAge dies. Those who die, switch to z2, lambda2 and start at ne/2+1 and x=1
	    age1=MaxAge
	    !$omp parallel do reduction(+:DBN_aux) private(x1,a1,lambda1,e1,z2,lambda2)
	    DO z1=1,nz
	    DO x1=1,nx
	    DO a1=1,na
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN_aux(1,Bqlo(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN_aux(1,Bqlo(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1) + DBN_tr(age1,a1,z1,lambda1,e1,x1,ti-1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN_aux(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN_aux(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1) + DBN_tr(age1,a1,z1,lambda1,e1,x1,ti-1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO    
	    ENDDO
	    ENDDO 
	    !$omp barrier   

		! retirees "e" stays the same for benefit retirement calculation purposes
		!$omp parallel do reduction(+:DBN_aux) private(x1,age1,a1,lambda1,e1,z2,lambda2,x2)
	    DO z1=1,nz
	    DO x1=1,nx
	    DO age1=RetAge-1, MaxAge-1
	    DO a1=1,na
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        ! Those who die, switch to z2, lambda2 and start at ne/2+1
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN_aux(1,Bqlo(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN_aux(1,Bqlo(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1) + DBN_tr(age1,a1,z1,lambda1,e1,x1,ti-1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN_aux(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN_aux(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1) + DBN_tr(age1,a1,z1,lambda1,e1,x1,ti-1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, and also e1 since they are retired
	        !e2=e1
	        DO x2=1,nx
		        DBN_aux(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) =  &
		        	& DBN_aux(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) + &
		        	& DBN_tr(age1, a1, z1, lambda1, e1,x1,ti-1) &
		            & * survP(age1) * PrAprimelo(age1,a1,z1,lambda1,e1,x1) * pr_x(x1,x2,z1,age1)     
		        DBN_aux(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) =  &
		          	& DBN_aux(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2) + &
		          	& DBN_tr(age1, a1, z1, lambda1, e1,x1,ti-1) &
		            & * survP(age1) * PrAprimehi(age1,a1,z1,lambda1,e1,x1) * pr_x(x1,x2,z1,age1)  
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    !$omp barrier
	    
	    ! Working age agents
	    !$omp parallel do reduction(+:DBN_aux) private(x1,age1,a1,lambda1,e1,z2,lambda2,x2,e2)
	    DO z1=1,nz
	    DO x1=1,nx
	    DO age1=1,RetAge-2
	    DO a1=1,na
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        ! Those who die, switch to z2, lambda2 and start at ne/2+1
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN_aux(1,Bqlo(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN_aux(1,Bqlo(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1) + DBN_tr(age1,a1,z1,lambda1,e1,x1,ti-1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqlo(age1,a1,z1,lambda1,e1,x1)
	            DBN_aux(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1)   =  &
	           		& DBN_aux(1,Bqhi(age1,a1,z1,lambda1,e1,x1),z2,lambda2,ne/2+1,1) + DBN_tr(age1,a1,z1,lambda1,e1,x1,ti-1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrBqhi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, but switch to e2 and x2
	        DO x2=1,nx
	        DO e2=1,ne
				DBN_aux(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) =  &
	          		& DBN_aux(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) + &
	          		& DBN_tr(age1, a1, z1, lambda1, e1,x1,ti-1) &
	                & * survP(age1) * pr_e(e1,e2) * pr_x(x1,x2,z1,age1) * PrAprimelo(age1,a1,z1,lambda1,e1,x1)
	            DBN_aux(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) =  &
	          		& DBN_aux(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2) + &
	          		& DBN_tr(age1, a1, z1, lambda1, e1,x1,ti-1) &
	                & * survP(age1) * pr_e(e1,e2) * pr_x(x1,x2,z1,age1) * PrAprimehi(age1,a1,z1,lambda1,e1,x1) 
	        ENDDO
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    !$omp barrier

	else 
		! Distribution is obtained from benchmark
		! In first period there is no updating with R 
		! Change is a shock and not anticipated
		DBN_aux = DBN_bench

	endif 


	! Total wealth - Supply of assets in the economy
	Wealth   = sum( sum(sum(sum(sum(sum(DBN_aux,6),5),4),3),1)*agrid )

	! Private demand is computed as total demand for capital in the economy 
	K_mat = K_matrix(R_in,P_tr(ti))
	Private_Demand = 0.0_dp
	do xi=1,nx
	do zi=1,nz 
	do ai=1,na
		Private_Demand = Private_Demand + sum(DBN_aux(:,ai,zi,:,:,xi)*K_mat(ai,zi,xi))
	enddo 
	enddo 
	enddo 

	! Public demand is computed as current debt plus governemnt deficit
		! Debt_tr(t) = (1+R)*Debt_tr(t-1) + Deficit(t)
		! Deficit(t) = GBAR_bench - GBAR_Tr(t) 
	if (ti.eq.1) then 
		Public_Demand = 0.0_dp
	else
		Public_Demand = Debt_Absorption*Debt_tr(ti-1)
	endif 


	! Aggregate debt is the sum of private and public demand for funds 
	Agg_Demand  = Private_Demand + Public_Demand
		! Agg_Demand = Private_Demand
	
	! Function outputs aggregate demand relative to total wealth (squared to have a min at 0)
	Agg_Debt_Tr = ((Agg_Demand-Wealth)/Wealth)**2.0_dp

	! print*, ' 	Agg_Debt_Error=',Agg_Debt_Tr,'	R_in=',R_in

end Function Agg_Debt_Tr


!========================================================================================
!========================================================================================
!========================================================================================
! THIS YGRID and Marginal Benefit of Investment GRID 
! NEEDS TO BE COMPUTED FOR EACH TIME THE INTEREST RATE "P" IS UPDATED. 

SUBROUTINE FORM_Y_MB_GRID(TYGRID,TMBGRID,TYGRID_t,TMBGRID_t)
	IMPLICIT NONE
	REAL(DP), DIMENSION(na,nz,nx),   INTENT(OUT) :: TYGRID, TMBGRID
	REAL(DP), DIMENSION(na_t,nz,nx), INTENT(OUT) :: TYGRID_t, TMBGRID_t
	!REAL(DP), INTENT(IN) :: P
	integer :: a_ind, z_ind, x_ind

	DO x_ind=1,nx 
	DO z_ind=1,nz
		DO a_ind=1,na
			TYGRID(a_ind,z_ind,x_ind)  = Y_a(agrid(a_ind),z_ind,x_ind)
			TMBGRID(a_ind,z_ind,x_ind) = MB_a(agrid(a_ind),z_ind,x_ind)
		ENDDO 
		if (Y_a_threshold.eq.0.0_dp) then
			TYGRID_t(:,z_ind,x_ind)  = TYGRID(:,z_ind,x_ind)
			TMBGRID_t(:,z_ind,x_ind) = TMBGRID(:,z_ind,x_ind) 
		else 
			DO a_ind=1,na_t
				TYGRID_t(a_ind,z_ind,x_ind)  = Y_a(agrid_t(a_ind),z_ind,x_ind)
				TMBGRID_t(a_ind,z_ind,x_ind) = MB_a(agrid_t(a_ind),z_ind,x_ind)
			ENDDO
		endif 
	ENDDO
	ENDDO

	!print *, "Grid for asset income"
	!do ai=1,na
	!	write(*,*) TMBGRID_t(ai,:)
 	!end do
	!pause

END SUBROUTINE FORM_Y_MB_GRID

SUBROUTINE FORM_Y_MB_GRID_Transition(TYGRID,TMBGRID,TYGRID_t,TMBGRID_t,time)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: time
	REAL(DP), DIMENSION(na,nz,nx),   INTENT(OUT) :: TYGRID, TMBGRID
	REAL(DP), DIMENSION(na_t,nz,nx), INTENT(OUT) :: TYGRID_t, TMBGRID_t
	!REAL(DP), INTENT(IN) :: P
	!integer :: ai, zi

	DO xi=1,nx 
	DO zi=1,nz
		DO ai=1,na
			P = P_tr(time); R = R_tr(time) 
			TYGRID(ai,zi,xi)  = Y_a(agrid(ai),zi,xi)
			P = P_tr(time+1); R = R_tr(time+1) 
			TMBGRID(ai,zi,xi) = MB_a(agrid(ai),zi,xi)
		ENDDO 
		if (Y_a_threshold.eq.0.0_dp) then
			TYGRID_t  = TYGRID
			TMBGRID_t = TMBGRID 
		else 
			DO ai=1,na_t
				P = P_tr(time); R = R_tr(time) 
				TYGRID_t(ai,zi,xi)  = Y_a(agrid_t(ai),zi,xi)
				P = P_tr(time+1); R = R_tr(time+1) 
				TMBGRID_t(ai,zi,xi) = MB_a(agrid_t(ai),zi,xi)
			ENDDO
		endif 
	ENDDO
	ENDDO
	! Return P and R to current time
	P = P_tr(time); R = R_tr(time) 

	!print *, "Grid for asset income"
	!do ai=1,na
	!	write(*,*) TMBGRID_t(ai,:)
 	!end do
	!pause

END SUBROUTINE FORM_Y_MB_GRID_Transition


!========================================================================================
!========================================================================================
!========================================================================================
! THIS SUBROUTINE COMPUTE EFFICIENCY UNITS OF LABOR DURING WORKING TIMES (NOT LABOR INCOME SINCE
! LABOR INCOME DEPENDS ON LABOR SUPPLY. IT ALSO fS RETIREMENT INCOME).
! THESE VARIABLES DEPEND ON EQUILIBRIUM WAGE AND EARNINGS AND AS A RESULT, WE NEED TO COMPUTE THESE
! FOR ALL EQUILIBRIUM WAGES AND EARNINGS

SUBROUTINE ComputeLaborUnits(Ebart,Waget)
	IMPLICIT NONE
	!integer:: lambdai  
	REAL(DP), INTENT(IN) :: Ebart, Waget 

	! This computes efficiency units times wage
	yh = Waget * eff_un

	! This part computes Retirement Income
	RetY_lambda_e = phi_lambda_e  * Ebart !* psi / psi_bench
	IF ((KeepSSatBench .eq. 1) .AND. (solving_bench .eq. 0)) THEN
    	RetY_lambda_e = phi_lambda_e  * Ebar_bench
	ENDIF

	if (solving_bench .eq. 1) then 
	OPEN  (UNIT=1,  FILE=trim(Result_Folder)//'Ret_Y'  , STATUS='replace')
	WRITE (UNIT=1,  FMT=*) RetY_lambda_e
	CLOSE (unit=1)
	endif 

END SUBROUTINE ComputeLaborUnits

!========================================================================================
!========================================================================================
!========================================================================================
! This subroutine computes after tax income given current taxes for each point in the state space

SUBROUTINE Compute_After_Tax_Income
	real(dp) :: capital_income

		DO xi=1,nx
		DO zi=1,nz			   
	    DO ai=1,na
	    	capital_income = YGRID(ai,zi,xi)
	    DO lambdai=1,nlambda
	    DO ei=1,ne
	    	DO age=RetAge,MaxAge
			Income_AT(age,ai,zi,lambdai,ei,xi)   = capital_income + RetY_lambda_e(lambdai,ei) 
			ENDDO !age        

			DO age=1,RetAge-1
			Income_AT(age,ai,zi,lambdai,ei,xi)   = capital_income + Y_h(hours(age,ai,zi,lambdai,ei,xi),age,lambdai,ei,wage)
			ENDDO !age        
		ENDDO !ei         
	    ENDDO !lambdai
		ENDDO !ai
		ENDDO !zi
		ENDDO !xi



END SUBROUTINE Compute_After_Tax_Income


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE  INITIALIZE()
	IMPLICIT NONE

	!integer::  lambdai
	REAL(DP) :: lambdaBAR, tempno 
	REAL(DP) :: m, Rh, start_timet, finish_timet
	INTEGER  :: ee0, ee1, ee2, zindx1, zindx2, lambdaindx1, lambdaindx2, diff_array, eindx1, eindx2, zi
	INTEGER, DIMENSION(RetAge) :: agevec
	! Entrepreneurial ability
		! grid (zgrid), invariant distribution (Glz)
		REAL(DP), DIMENSION(nz_aux)    :: zgrid_aux , Gz_aux
		! transition matrix (pr_z)
		REAL(DP), DIMENSION(nz_aux,nz_aux) :: pr_z_aux
	
	! Initiliaze grids for z, lamda and e	
		CALL tauchen(mtauchen_z,rho_z,sigma_z_eps,nz_aux,zgrid_aux,pr_z_aux,Gz_aux)
		CALL tauchen(mtauchen,rho_e,sigma_e_eps,ne,egrid,pr_e,Ge)
		CALL tauchen(mtauchen,rho_lambda,sigma_lambda_eps,nlambda,lambdagrid,pr_lambda,Glambda)

		! Tauchen gives grids for the log of the variables. Exponentiate to adjust
		zgrid_aux  = exp(zgrid_aux) + mu_z
		egrid      = exp(egrid)
		lambdagrid = exp(lambdagrid) 

		! Cut bottom elements of zgrid 
		CALL Markov_Cut(nz_aux,zgrid_aux,pr_z_aux,Gz_aux,nz_aux-nz,zgrid,pr_z,Gz)
		! print*,'G_z = ',Gz

	! Transitory investment productivity x
		if (nx.gt.1) then
			! print*, 'X probability '
			xgrid = (/x_hi , x_lo , x_0/)
			! Low z types stay in x=1 until retirement
				pr_x(1,1,1:4,:) = 1.00_dp - p2_x
				pr_x(1,2,1:4,:) = 0.00_dp 
				pr_x(1,3,1:4,:) = p2_x
			! High z types have 5% probability of going from x=1 to x=2
				pr_x(1,1,5:nz,:) = 1.00_dp - p1_x - p2_x
				pr_x(1,2,5:nz,:) = p1_x 
				pr_x(1,3,5:nz,:) = p2_x
			! x=2 goes to x=3 with probability 3%
				pr_x(2,1,:,:) = 0.00_dp 
				pr_x(2,2,:,:) = 1.00_dp - p2_x
				pr_x(2,3,:,:) = p2_x
			! x=3 is an absorbing state
				pr_x(3,1,:,:) = 0.00_dp 
				pr_x(3,2,:,:) = 0.00_dp 
				pr_x(3,3,:,:) = 1.00_dp
			! Gx is not used. So it is initialized to an arbitrary value
				Gx(1,:,:) = 0.50_dp ; Gx(2,:,:) = 0.50_dp ; Gx(3,:,:) = 0.00_dp ;
			! xz grid
				xz_grid(1,:)   = exp(log(zgrid)*xgrid(1))
				xz_grid(2,1:4) = xz_grid(1,1:4); xz_grid(2,5:) = exp(log(zgrid(5:))*xgrid(2));
				xz_grid(3,:)   = 0.0_dp
				! xz_grid = spread(zgrid,1,nx)*spread(xgrid,2,nz)
				! xz_grid(1,:)   = zgrid 	; xz_grid(2,1:3) = zgrid(1:3)	;	xz_grid(2,4:)  = zgrid(4)
				! xz_grid(1,:) = zgrid; xz_grid(2,:) = 0.00_dp*zgrid
				! print*, ' xgrid', xgrid
				! print*, ' zgrid', zgrid 
				! do xi=1,nx
				! print*, 'xzgrid', xz_grid(xi,:)
				! enddo
				! print*, 'xgrid error'
				! STOP 
		else 
			xgrid = 1.0_dp
			xz_grid(1,:) = zgrid 
			pr_x(1,1,:,:)  = 1.0_dp
			Gx(1,:,:)      = 1.0_dp
		endif 

		! Obtain CDF of invariant distribution and transition matrix
			! Entrepreneurial ability
			DO ee0 = 1,nz
			    cdf_Gz(ee0) = sum(Gz(1:ee0))
			    DO ee1 = 1,nz
			        cdf_pr_z(ee0,ee1) = sum(pr_z(ee0,1:ee1))
			    ENDDO
			ENDDO
			do age=1,MaxAge
			do zi=1,nz
				DO ee0 = 1,nx
				    cdf_Gx(ee0,zi,age) = sum(Gx(1:ee0,zi,age))
				    DO ee1 = 1,nx
				        cdf_pr_x(ee0,ee1,zi,age) = sum(pr_x(ee0,1:ee1,zi,age))
				    ENDDO
				ENDDO
			enddo 
			enddo

			! print*, ' '
			! print*, 'CDF of X'
			! do age=1,MaxAge 
			! 	print*, ' '
			! 	do zi=1,nz 
			! 		print*, 'Age=',age,'zi=',zi
			! 		do ee0=1,nx 
			! 			print*, cdf_pr_x(ee0,:,zi,age)
			! 		enddo 
			! 	enddo 
			! enddo 

			! Labor income permanent component
			DO ee0 = 1,nlambda
			    cdf_Glambda(ee0) = sum(Glambda(1:ee0))
			    DO ee1 = 1,nlambda
			        cdf_pr_lambda(ee0,ee1) = sum(pr_lambda(ee0,1:ee1))
			    ENDDO
			ENDDO
			! Labor income transitory component
			DO ee0 = 1,ne
			    cdf_Ge(ee0) = sum(Ge(1:ee0))
			    DO ee1 = 1,ne
			        cdf_pr_e(ee0,ee1) = sum(pr_e(ee0,1:ee1))
			    ENDDO
			ENDDO

			OPEN  (UNIT=1,  FILE=trim(Result_Folder)//'P_e'  , STATUS='replace')
			WRITE (UNIT=1,  FMT=*) cdf_pr_e
			CLOSE (unit=1)

	! Distribution of "e" by age: Ge_byage(age,e)
		! We assume that in the first period, everyone starts at median "e"
		Ge_byage(:,:)      = 0.0_DP
		Ge_byage(1,ne/2+1) = 1.0_DP  
		DO age=1, MaxAge-1
		      DO ee1=1,ne
		           DO ee0=1,ne
			            Ge_byage(age+1,ee0) = Ge_byage(age+1,ee0)  + Ge_byage(age,ee1)*pr_e(ee1,ee0)
		            END DO
		       END DO
		!       print*,Ge_byage(age,:)
		END DO
		!print*,Ge_byage(age,:)
	    
	    ! Compute CDF of distribution of "e" by age
		DO age=1, MaxAge
		      DO ee0 = 1, ne
		            cdf_Ge_byage(age,ee0) = sum(Ge_byage(age,1:ee0))
		      ENDDO
		ENDDO

	! Efficiency units of labor
		! NOTE THAT THIS IS NOT LABOR INCOME
		eff_un=0.0_DP
		Rh=1.0_DP
		DO age=1,RetAge
	        agevec(age) = age
	        ! Compute life cycle component
	        kappagrid(age) = exp(  (60.0_DP *(Rh-1.0_DP)-(Rh-1.0_DP)**2.0_DP)/1800.0_DP )
	        ! kappagrid(age) = 1.0_DP
	        DO lambdai=1,nlambda        
                DO ei=1,ne                
                    eff_un(age,lambdai,ei) = kappagrid(age) * lambdagrid(lambdai) * egrid(ei)        
                ENDDO
	        ENDDO
	        Rh=Rh+1.0_DP
		ENDDO
		!print*,'kappagrid='
		!print*,kappagrid-exp(  (60.0_DP *(agevec-1.0_DP)-(agevec-1.0_DP)**2.0_DP)/1800.0_DP )
		!PAUSE
		OPEN  (UNIT=1,  FILE=trim(Result_Folder)//'eff_un'  , STATUS='replace')
		WRITE (UNIT=1,  FMT=*) eff_un
		CLOSE (unit=1)

	! Initialize asset grid
		! Normal grid
		m=(amax-amin)**(1.0_DP/a_theta)/REAL(na-1,DP)
		DO ai=1,na
			agrid(ai)=REAL(ai-1,DP)*m
		END DO
		agrid=amin+agrid**a_theta
		!	print*,'agrid=',agrid
		!!pause

		OPEN(UNIT=12, FILE=trim(Result_Folder)//'agrid', STATUS='replace')
		WRITE(unit=12, FMT=*) agrid
		CLOSE (unit=12)

		OPEN(UNIT=12, FILE=trim(Result_Folder)//'zgrid', STATUS='replace')
		WRITE(unit=12, FMT=*) zgrid
		CLOSE (unit=12)

		OPEN(UNIT=12, FILE=trim(Result_Folder)//'lambdagrid', STATUS='replace')
		WRITE(unit=12, FMT=*) lambdagrid
		CLOSE (unit=12)

		OPEN(UNIT=12, FILE=trim(Result_Folder)//'egrid', STATUS='replace')
		WRITE(unit=12, FMT=*) egrid
		CLOSE (unit=12)

		! Fine grid
		m=(amax-amin)**(1.0_DP/a_theta)/REAL(fine_na-1,DP)
		DO ai=1,fine_na
			fine_agrid(ai)=REAL(ai-1,DP)*m
		END DO
		fine_agrid=amin+fine_agrid**a_theta

	!----------------------------------------------
	! life-cycle component - Survival probablity
	!----------------------------------------------
		
	! Population Numbers from Bell and Miller (2002)
		pop(1)=	197316.0_DP
		pop(2)=	197141.0_DP
		pop(3)=	196959.0_DP
		pop(4)=	196770.0_DP
		pop(5)=	196580.0_DP
		pop(6)=	196392.0_DP
		pop(7)=	196205.0_DP
		pop(8)=	196019.0_DP
		pop(9)=	195830.0_DP
		pop(10)=195634.0_DP
		pop(11)=195429.0_DP
		pop(12)=195211.0_DP
		pop(13)=194982.0_DP
		pop(14)=194739.0_DP
		pop(15)=194482.0_DP
		pop(16)=194211.0_DP
		pop(17)=193924.0_DP
		pop(18)=193619.0_DP
		pop(19)=193294.0_DP
		pop(20)=192945.0_DP
		pop(21)=192571.0_DP
		pop(22)=192169.0_DP
		pop(23)=191736.0_DP
		pop(24)=191271.0_DP
		pop(25)=190774.0_DP
		pop(26)=190243.0_DP
		pop(27)=189673.0_DP
		pop(28)=189060.0_DP
		pop(29)=188402.0_DP
		pop(30)=187699.0_DP
		pop(31)=186944.0_DP
		pop(32)=186133.0_DP
		pop(33)=185258.0_DP
		pop(34)=184313.0_DP
		pop(35)=183290.0_DP
		pop(36)=182181.0_DP
		pop(37)=180976.0_DP
		pop(38)=179665.0_DP
		pop(39)=178238.0_DP
		pop(40)=176689.0_DP
		pop(41)=175009.0_DP
		pop(42)=173187.0_DP
		pop(43)=171214.0_DP
		pop(44)=169064.0_DP
		pop(45)=166714.0_DP
		pop(46)=164147.0_DP
		pop(47)=161343.0_DP
		pop(48)=158304.0_DP
		pop(49)=155048.0_DP
		pop(50)=151604.0_DP
		pop(51)=147990.0_DP
		pop(52)=144189.0_DP
		pop(53)=140180.0_DP
		pop(54)=135960.0_DP
		pop(55)=131532.0_DP
		pop(56)=126888.0_DP
		pop(57)=122012.0_DP
		pop(58)=116888.0_DP
		pop(59)=111506.0_DP
		pop(60)=105861.0_DP
		pop(61)=99957.0_DP
		pop(62)=93806.0_DP
		pop(63)=87434.0_DP
		pop(64)=80882.0_DP
		pop(65)=74204.0_DP
		pop(66)=67462.0_DP
		pop(67)=60721.0_DP
		pop(68)=54053.0_DP
		pop(69)=47533.0_DP
		pop(70)=41241.0_DP
		pop(71)=35259.0_DP
		pop(72)=29663.0_DP
		pop(73)=24522.0_DP
		pop(74)=19890.0_DP
		pop(75)=15805.0_DP
		pop(76)=12284.0_DP
		pop(77)=9331.0_DP
		pop(78)=6924.0_DP
		pop(79)=5016.0_DP
		pop(80)=3550.0_DP

		if (MaxAge .gt. 80)	 then	
			pop(MaxAge)=2454.0_DP
		endif	
		! pop = 1.0_DP	

	
	! Survival probabilities: surv(i)=prob(alive in i+1|alive in i)
		FORALL (age=1:maxAge-1) survP(age)= pop(age+1)/pop(age)
		survP(maxAge)=0.0_DP
		OPEN  (UNIT=1,  FILE=trim(Result_Folder)//'survP'  , STATUS='replace')
		WRITE (UNIT=1,  FMT=*) survP
		CLOSE (unit=1)
			


	! Set the initial distribution
	DBN1=0.0_DP
		DO xi=1,nx
		DO age=1,MaxAge
		DO zi=1,nz
		DO lambdai=1,nlambda
		DO ei=1, ne
		    DBN1(age,1,zi,lambdai,ei,xi) = (pop(age)/sum(pop))*Gz(zi)*Glambda(lambdai)*Ge_byage(age,ei)
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		ENDDO  
		! print*,'Initial Distrbution', sum(DBN1)
		DBN1 = DBN1/sum(DBN1)

	CALL LIFETIME_Y_ESTIMATE

END SUBROUTINE INITIALIZE


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE  LIFETIME_Y_ESTIMATE()
	IMPLICIT NONE
	integer  :: currentzi, currentlambdai, currentei
	REAL(DP) :: tempnoage, tempnoz, tempnolambda, tempnoe, tempno, currenta, currentY
	REAL(DP) :: start_timet, finish_timet
	INTEGER  :: paneli, cohortsize
	REAL(DP) :: mean_panel_lifetime_eff_unit , lambdaBAR
	INTEGER , DIMENSION(:), allocatable :: panele, panellambda
	REAL(DP), DIMENSION(:), allocatable :: panel_lifetime_eff_unit
	REAL(DP), DIMENSION(nlambda,ne) :: lifetime_eff_unit_by_lambda_e, size_by_lambda_e 

	cohortsize=10000000
	allocate(panele(cohortsize),panellambda(cohortsize),panel_lifetime_eff_unit(cohortsize))

	newiseed=-1

	panel_lifetime_eff_unit = 0.0_DP
	age=1

	DO paneli=1, cohortsize   
		! LAMBDA  
	   tempnolambda = ran1(newiseed) 
	   lambdai=1
	   DO WHILE (tempnolambda .gt. cdf_Glambda(lambdai))
	       lambdai=lambdai+1
	   ENDDO
	   panellambda(paneli)=lambdai

	   ei=ne/2+1
	   panele(paneli)=ei
	   panel_lifetime_eff_unit(paneli) =  eff_un(age,lambdai,ei)   
	ENDDO ! paneli

	DO age=2, RetAge-1
        DO paneli=1, cohortsize  
            currentei = panele(paneli)   
            tempno = ran1(newiseed)   
            ei=1
            DO WHILE (tempno .gt. cdf_pr_e(currentei,ei))
                ei=ei+1
            ENDDO
            panele(paneli)=ei    
            panel_lifetime_eff_unit(paneli) =  panel_lifetime_eff_unit(paneli) + eff_un(age,panellambda(paneli),ei)   
        ENDDO ! paneli     
	ENDDO ! age
	panel_lifetime_eff_unit = panel_lifetime_eff_unit/REAL(RetAge-1,8)
	mean_panel_lifetime_eff_unit = sum(panel_lifetime_eff_unit)/cohortsize  


	lifetime_eff_unit_by_lambda_e=0.0_DP
	size_by_lambda_e = 0.0_DP
	DO paneli=1,cohortsize
		ei = panele(paneli)   
		lambdai = panellambda(paneli)   
		lifetime_eff_unit_by_lambda_e(lambdai, ei) = lifetime_eff_unit_by_lambda_e(lambdai, ei) + panel_lifetime_eff_unit(paneli)
		size_by_lambda_e(lambdai, ei) = size_by_lambda_e(lambdai, ei) +1.0_DP       
	ENDDO
	lifetime_eff_unit_by_lambda_e = lifetime_eff_unit_by_lambda_e/ size_by_lambda_e
	lifetime_eff_unit_by_lambda_e = lifetime_eff_unit_by_lambda_e/mean_panel_lifetime_eff_unit 


	! compute retirement income replacement rate for each lambda type
	!print*,'phi_lambda_e='
	DO lambdai=1,nlambda
	DO ei=1,ne
	    IF (lifetime_eff_unit_by_lambda_e(lambdai, ei) .le. 0.3_DP) THEN
		
			phi_lambda_e(lambdai,ei) = 0.9_DP*  lifetime_eff_unit_by_lambda_e(lambdai, ei)
	    
	    ELSE IF   (lifetime_eff_unit_by_lambda_e(lambdai, ei) .le. 2.0_DP) THEN
        
            phi_lambda_e(lambdai,ei)  = 0.27_DP +  0.32_DP*  (lifetime_eff_unit_by_lambda_e(lambdai, ei)-0.3_DP) 
		
		ELSE IF (lifetime_eff_unit_by_lambda_e(lambdai, ei) .le. 4.1_DP) THEN
            
            phi_lambda_e(lambdai,ei)  = 0.81_DP +  0.15_DP*  (lifetime_eff_unit_by_lambda_e(lambdai, ei)-2.0_DP)
		
		ELSE
            
            phi_lambda_e(lambdai,ei)  =1.1_DP	   
        
        ENDIF
	ENDDO 
	!print*,phi_lambda_e(lambdai,:)       
	ENDDO
	!PAUSE

END SUBROUTINE LIFETIME_Y_ESTIMATE

!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE  Firm_Value()
	use omp_lib

	IMPLICIT NONE
	integer  :: age, ai, zi, lambdai, ei, xi_p, tklo, tkhi, zp_ind
	real(dp) :: dV_low, dV_high, spline_coeff(na), V_spline_R, sp_coeff_W(na,ne), V_spline_W(nx), Prob_lo, Prob_hi, V_Pr_aux(nz)

	!$ call omp_set_num_threads(5)

	K_mat  = K_Matrix(R,P)
	Pr_mat = Profit_Matrix(R,P)

	! Final period of life
	age = MaxAge 
	do xi=1,nx
	do zi=1,nz
	do ai=1,na 
	do ei=1,ne 
	do lambdai=1,nlambda
		V_Pr(age,ai,zi,lambdai,ei,xi) = Pr_mat(ai,zi,xi)
	enddo 
	enddo 
	enddo
	enddo 
	enddo 

	! Retirement Periods
	do zi=1,nz
		! ! Prepare spline interpolation
		! 	! Derivative of value function in first grid point
		! 	if (K_mat(1,zi).lt.theta*agrid(1)) then 
		! 		dV_low  = 0.0_dp
		! 	else 
		! 		dV_low  = P*mu*((theta*zgrid(zi))**mu)*agrid(1)**(mu-1.0_DP)-(R+DepRate)*theta
		! 	endif
		! 	! Derivative of value function in last grid point
		! 	if (K_mat(na,zi).lt.theta*agrid(na)) then 
		! 		dV_high = 0.0_dp
		! 	else 
		! 		dV_high = P*mu*((theta*zgrid(zi))**mu)*agrid(na)**(mu-1.0_DP)-(R+DepRate)*theta
		! 	endif

		do age=MaxAge-1,RetAge,-1
		!$omp parallel do private(lambdai,ei,ai,xi,spline_coeff,V_spline_R,tklo,tkhi,Prob_lo,Prob_hi)
		do xi=1,nx
		do ei=1,ne 
		do lambdai=1,nlambda

			!CALL spline( agrid, V_Pr(age+1, :, zi, lambdai, ei) , na , dV_low , dV_high , spline_coeff)  
		
			do ai=1,na
				! call splint( agrid, V_Pr(age+1, :, zi, lambdai, ei), spline_coeff, na, Aprime(age,ai,zi,lambdai, ei), V_spline_R ) 
				! V_Pr(age,ai,zi,lambdai,ei) = Pr_mat(ai,zi) + survP(age)/(1.0_dp+R) * V_spline_R

				! Prepare linear interpolation
				if ( Aprime(age,ai,zi,lambdai, ei,xi) .ge. amax) then
					tklo =na-1
				elseif (Aprime(age,ai,zi,lambdai, ei,xi) .lt. amin) then
				    tklo = 1
				else
				    tklo = ((Aprime(age,ai,zi,lambdai, ei,xi) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
				endif            
				tkhi = tklo + 1        
				Prob_lo = ( agrid(tkhi) - Aprime(age,ai,zi,lambdai, ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
				Prob_hi = ( Aprime(age,ai,zi,lambdai, ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )        
				Prob_lo = min(Prob_lo, 1.0_DP)
				Prob_lo = max(Prob_lo, 0.0_DP)
				Prob_hi = min(Prob_hi, 1.0_DP)
				Prob_hi = max(Prob_hi, 0.0_DP)    

				V_Pr(age,ai,zi,lambdai,ei,xi) = Pr_mat(ai,zi,xi) + survP(age)/(1.0_dp+MeanReturn) &
			  		&  * sum( pr_x(xi,:,zi,age) * (Prob_lo*V_Pr(age+1,tklo,zi,lambdai,ei,:)+Prob_hi*V_Pr(age+1,tkhi,zi,lambdai,ei,:) ) ) 

			enddo
	enddo
	enddo
	enddo
	enddo
	enddo

	! Working Periods
	do zi=1,nz 
		! ! Prepare spline interpolation
		! 	! Derivative of value function in first grid point
		! 	if (K_mat(1,zi).lt.theta*agrid(1)) then 
		! 		dV_low  = 0.0_dp
		! 	else 
		! 		dV_low  = P*mu*((theta*zgrid(zi))**mu)*agrid(1)**(mu-1.0_DP)-(R+DepRate)*theta
		! 	endif
		! 	! Derivative of value function in last grid point
		! 	if (K_mat(na,zi).lt.theta*agrid(na)) then 
		! 		dV_high = 0.0_dp
		! 	else 
		! 		dV_high = P*mu*((theta*zgrid(zi))**mu)*agrid(na)**(mu-1.0_DP)-(R+DepRate)*theta
		! 	endif

		do age=RetAge-1,1,-1
		!$omp parallel do private(lambdai,ei,ai,xi,sp_coeff_W,V_spline_W,xi_p,tklo,tkhi,Prob_lo,Prob_hi)
		do xi=1,nx
		do lambdai=1,nlambda
			! do ei_p=1,ne
			! 	CALL spline( agrid, V_Pr(age+1, :, zi, lambdai, ei_p) , na , dV_low , dV_high , sp_coeff_W(:,ei_p))
			! enddo   
		
		do ei=1,ne 
		do ai=1,na
			if ( Aprime(age,ai,zi,lambdai,ei,xi) .ge. amax) then
				tklo =na-1
			elseif (Aprime(age,ai,zi,lambdai,ei,xi) .lt. amin) then
			    tklo = 1
			else
			    tklo = ((Aprime(age,ai,zi,lambdai,ei,xi) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
			endif            
			tkhi = tklo + 1        
			Prob_lo = ( agrid(tkhi) - Aprime(age,ai,zi,lambdai,ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
			Prob_hi = ( Aprime(age,ai,zi,lambdai,ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )        
			Prob_lo = min(Prob_lo, 1.0_DP)
			Prob_lo = max(Prob_lo, 0.0_DP)
			Prob_hi = min(Prob_hi, 1.0_DP)
			Prob_hi = max(Prob_hi, 0.0_DP) 

			do xi_p=1,nx
				! call splint( agrid, V_Pr(age+1, :, zi, lambdai, ei_p), sp_coeff_W, na, Aprime(age,ai,zi,lambdai, ei), V_spline_W(ei_p) )
				V_spline_W(xi_p) = sum( pr_e(ei,:)*( Prob_lo*V_Pr(age+1,tklo,zi,lambdai,:,xi_p) + Prob_hi*V_Pr(age+1,tkhi,zi,lambdai,:,xi_p)) )
			enddo

			V_Pr(age,ai,zi,lambdai,ei,xi) = Pr_mat(ai,zi,xi) + survP(age)/(1.0_dp+MeanReturn) * sum(pr_x(xi,:,zi,age)*V_spline_W)

		enddo
	enddo
	enddo
	enddo
	enddo
	enddo


	! Define Firm based Wealth measure
	do ai=1,na
		Firm_Wealth(:,ai,:,:,:,:) = V_Pr(:,ai,:,:,:,:) + (1.0_dp+R)*agrid(ai)
	enddo

	! Value for new borns
	!$omp parallel do private(lambdai,ai,zp_ind,V_Pr_aux)
	do zi=1,nz 
	do lambdai=1,nlambda
	do ai=1,na 
		do zp_ind=1,nz 
			V_Pr_aux(zp_ind) = sum(pr_lambda(lambdai,:)*Firm_Wealth(1,ai,zp_ind,:,ne/2+1,1))
		enddo 
		V_Pr_nb(ai,zi,lambdai) = sum(pr_z(zi,:) * V_Pr_aux )
	enddo 
	enddo 
	enddo 


END SUBROUTINE Firm_Value


!========================================================================================
!========================================================================================
!========================================================================================

    function omp_ran1(idum)
        ! returns a uniform random number between 0 and 1
        ! see numerical recipes
        ! press,flannery,teukolsky & vetterling
        ! cambridge university press 1986 pp 191-2-3
        ! http://geco.mines.edu/prototype/Show_me_some_local_HPC_tutorials/examples/openmp/ranmod.f90
        ! http://inside.mines.edu/~tkaiser/fortran//
        ! https://people.sc.fsu.edu/~jburkardt/f_src/random_openmp/random_openmp.html
        ! http://jblevins.org/log/openmp
        implicit none
        integer, intent(inout), optional :: idum
        real(DP) :: omp_ran1
        integer  :: m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3
        integer  :: iff
        integer  :: ix1,ix2,ix3,j
        real(DP) :: r(97),rm1,rm2
        
        parameter (m1=259200,ia1=7141,ic1=54773)
        parameter (m2=134456,ia2=8121,ic2=28411)
        parameter (m3=243000,ia3=4561,ic3=51349)
        data iff /0/
        
        !$OMP THREADPRIVATE(iff,ix1,ix2,ix3,j,r,rm1,rm2) 
        save iff,ix1,ix2,ix3,j,r,rm1,rm2
        if(present(idum))then
            if (idum<0.or.iff.eq.0)then
                rm1=1.0_dp/m1
                rm2=1.0_dp/m2
                iff=1
                ix1=mod(ic1-idum,m1)
                ix1=mod(ia1*ix1+ic1,m1)
                ix2=mod(ix1,m2)
                ix1=mod(ia1*ix1+ic1,m1)
                ix3=mod(ix1,m3)
                do  j=1,97
                    ix1=mod(ia1*ix1+ic1,m1)
                    ix2=mod(ia2*ix2+ic2,m2)
                    r(j)=(float(ix1)+float(ix2)*rm2)*rm1
                enddo 
                idum=1
            endif
        endif
        
        ix1=mod(ia1*ix1+ic1,m1)
        ix2=mod(ia2*ix2+ic2,m2)
        ix3=mod(ia3*ix3+ic3,m3)
        j=1+(97*ix3)/m3

        if(j>97.or.j<1)then
            write(*,*)' error in omp_ran1 j=',j
            stop
        endif

        omp_ran1=r(j)
        r(j)=(float(ix1)+float(ix2)*rm2)*rm1
        return
    end function omp_ran1



!========================================================================================
!========================================================================================
!========================================================================================

end module programfunctions



