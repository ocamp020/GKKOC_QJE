
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
		Y_a = ( a_in +  ( Pr + R*a_in ) *(1.0_DP-tauK) )

		! Compute after tax wealth according to threshold
		if (a_in.le.Y_a_threshold) then 
			Y_a = a_in* (1.0_dp-tauW_bt) + ( Pr + R*a_in ) *(1.0_DP-tauK)  
		else
			Y_a = Y_a_threshold*(1.0_dp-tauW_bt) + (a_in - Y_a_threshold) * (1.0_dp-tauW_at) + ( Pr + R*a_in ) *(1.0_DP-tauK)
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
		Y_a = ( a_in +  ( Pr + R*a_in ) *(1.0_DP-tauK) )
		if (a_in.le.Y_a_threshold) then 
			tauW = tauW_bt 
		else
			tauW = tauW_at 
		end if

		! After tax marginal benefit of assets
		if (K.lt.theta(z_in)*a_in) then 
			MB_a = (1.0_dp*(1.0_dp-tauW) + R*(1.0_dp-tauK))
		else 
			MB_a = (1.0_dp*(1.0_dp-tauW) + R*(1.0_dp-tauK)) &
         	& + (P*mu*((theta(z_in)*xz_grid(x_in,z_in))**mu)*a_in**(mu-1.0_DP)-(R+DepRate)*theta(z_in))*(1.0_dp-tauK)
		endif 

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
			MB_a_at = (1.0_dp*(1.0_dp-tauW_at) + R*(1.0_dp-tauK))
		else 
			MB_a_at = (1.0_dp*(1.0_dp-tauW_at) + R*(1.0_dp-tauK)) &
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
			MB_a_bt = (1.0_dp*(1.0_dp-tauW_bt) + R*(1.0_dp-tauK))
		else 
			MB_a_bt = (1.0_dp*(1.0_dp-tauW_bt) + R*(1.0_dp-tauK)) &
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
		           & - ( beta*survP(age_in) * sum( pr_x(x_in,:,z_in,age_in) * MB_aprime * cprime**(1.0_dp/euler_power) ) ) &
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
		           & - ( beta*survP(age_in) * sum( pr_x(x_in,:,z_in,age_in) * MB_aprime * cprime**(1.0_dp/euler_power) ) ) &
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
				FOC_WH   = ( (1.0_dp/ctemp)- (beta*survP(age_in) * sum( pr_x(x_in,:,z_in,age_in) * MB_aprime * E_MU_cp) ) ) **2.0_DP 
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
					         & - beta*survP(age_in)* sum(pr_x(x_in,:,z_in,age_in)*MB_aprime*E_MU_cp)  )**2.0_DP
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
					& - beta*survP(age_in)* sum( pr_x(x_in,:,z_in,age_in)*MB_aprime*E_MU_cp ) )**2.0_DP
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
			FOC_WH   = (ctemp - 1.0_dp/(beta*survP(age_in)*sum(pr_x(x_in,:,z_in,age_in)*MB_aprime*E_MU_cp))**(1.0_dp/sigma)) **2.0_DP 
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
					& ( (1.0_dp/ctemp)- (beta*survP(age_in) * sum( pr_x(x_in,:,z_in,age_in) * MB_aprime * E_MU_cp) ) ) **2.0_DP 
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
					         & - beta*survP(age_in)* sum(pr_x(x_in,:,z_in,age_in)*MB_aprime*E_MU_cp)  )**2.0_DP
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
					& - beta*survP(age_in)* sum( pr_x(x_in,:,z_in,age_in)*MB_aprime*E_MU_cp ) )**2.0_DP
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
				& (ctemp - 1.0_dp/(beta*survP(age_in)*sum(pr_x(x_in,:,z_in,age_in)*MB_aprime*E_MU_cp))**(1.0_dp/sigma)) **2.0_DP
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
			         	& - beta*survP(age_in)*sum(pr_x(x_in,:,z_in,age_in)*MB_aprime*E_MU_cp)  )**2.0_DP

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
			         	& - beta*survP(age_in)*sum(pr_x(x_in,:,z_in,age_in)*MB_aprime*E_MU_cp)  )**2.0_DP

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

				C_euler = ( (beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_MU_cp))) **(1.0_dp/((1.0_dp-sigma)*gamma-1.0_dp))
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
				C_endo = 1.0_dp/( (beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_mu_cp)  ) **(1.0_dp/sigma) )
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
				    & *  beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_MU_cp)  )**(-1.0_DP/sigma)

				  H_endo = 1.0_DP - (1.0_DP-gamma)*C_endo/(gamma*psi*yh(age,lambdai,ei))   

				If (H_endo .lt. 0.0_DP) then
				    H_endo = 0.0_DP 
				    C_endo  = ( beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_MU_cp) )**(1.0_DP/(gamma*(1.0_DP-sigma)-1.0_DP))
				endif 

				! print*,' '
				! print*,' 	Inside EGM'
				! print*,' 	Consumption_t+1=',Cons_t(age+1,ai,zi,lambdai,:,:)
			else 
				! Separable Utility
				do xp_ind=1,nx
					E_MU_cp(xp_ind) = sum( pr_e(ei,:) * (Cons_t(age+1,ai,zi,lambdai,:,xp_ind)**(-sigma)) )
				enddo
				C_endo  = 1.0_DP/( beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_MU_cp) )**(1.0_DP/sigma)

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

!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE COMPUTE_WELFARE_GAIN()
	IMPLICIT NONE
	real(DP), dimension(MaxAge):: CumDiscountF
	! REAL(DP), dimension(MaxAge, na, nz, nlambda, ne) ::  Cons_Eq_Welfare ! ValueFunction_Bench, ValueFunction_Exp,
	REAL(DP), dimension(nz) ::  temp_ce_by_z, temp_cons_by_z, temp_leisure_by_z, temp_dbn_by_z 
	REAL(DP) :: frac_pos_welfare 
	REAL(DP), dimension(MaxAge, nz) :: frac_pos_welfare_by_age_z, size_pos_welfare_by_age_z, size_by_age_z_bench, size_by_age_z_exp
	INTEGER, dimension(max_age_category+1) :: age_limit
	INTEGER :: age_group_counter
	REAL(DP), dimension(max_age_category,nz) :: CE_by_agegroup_z 
	REAL(DP), dimension(max_age_category,nz) :: size_pos_welfare_by_agegroup_z, frac_pos_welfare_by_agegroup_z  
	REAL(DP), dimension(max_age_category,nz) :: tot_wealth_by_agegroup_z_bench, size_by_agegroup_z_bench, size_by_agegroup_z_exp
	REAL(DP), dimension(max_age_category,nz) :: tot_wealth_by_agegroup_z_exp
	REAL(DP), dimension(max_age_category,nz) :: tot_cons_by_agegroup_by_z_bench, tot_hours_by_agegroup_by_z_bench
	REAL(DP), dimension(max_age_category,nz) :: tot_aprime_by_agegroup_by_z_bench, tot_return_by_agegroup_by_z_bench
	REAL(DP), dimension(max_age_category,nz) :: tot_at_return_by_agegroup_by_z_bench, tot_cap_tax_by_agegroup_by_z_bench
	real(DP), dimension(max_age_category,nz) :: tot_lab_tax_by_agegroup_by_z_bench, tot_cons_tax_by_agegroup_by_z_bench
	REAL(DP), dimension(max_age_category,nz) :: tot_cons_by_agegroup_by_z_exp, tot_hours_by_agegroup_by_z_exp
	REAL(DP), dimension(max_age_category,nz) :: tot_aprime_by_agegroup_by_z_exp, tot_return_by_agegroup_by_z_exp
	REAL(DP), dimension(max_age_category,nz) :: tot_at_return_by_agegroup_by_z_exp, tot_cap_tax_by_agegroup_by_z_exp
	REAL(DP), dimension(max_age_category,nz) :: tot_lab_tax_by_agegroup_by_z_exp, tot_cons_tax_by_agegroup_by_z_exp
	REAL(DP), dimension(draft_age_category,nz) :: CE_draft_group_z,  size_draft_group_z, frac_pos_welfare_draft_group_z
	REAL(DP), dimension(draft_age_category,nz) :: wealth_draft_group_z
	REAL(DP), dimension(draft_age_category,nz) :: K_Tax_draft_group_z, L_Tax_draft_group_z, C_Tax_draft_group_z
	REAL(DP), dimension(draft_age_category,nz) :: Tot_Income_draft_group_z, K_Inc_draft_group_z, L_Inc_draft_group_z
	REAL(DP), dimension(draft_age_category,nz) :: K_Inc_frac_draft_group_z, L_Inc_frac_draft_group_z
	REAL(DP), dimension(draft_age_category,nz) :: K_Tax_Inc_draft_group_z, L_Tax_Inc_draft_group_z, C_Tax_Inc_draft_group_z
	REAL(DP), dimension(draft_age_category,nz) :: Tax_Increase_tk_draft_group_z, Tax_Increase_tl_draft_group_z
	REAL(DP), dimension(draft_age_category,nz) :: Tax_Increase_draft_group_z, Tax_Rate_Increase_draft_group_z
	REAL(DP), dimension(draft_age_category,nz) :: Tax_Rate_Increase_tk_draft_group_z, Tax_Rate_Increase_tl_draft_group_z
	REAL(DP), dimension(draft_age_category,draft_z_category) :: CE_draft_group,  size_draft_group, frac_pos_welfare_draft_group
	REAL(DP), dimension(draft_age_category,draft_z_category) :: wealth_draft_group,  av_wealth_draft_group, frac_wealth_draft_group
	REAL(DP), dimension(draft_age_category,draft_z_category) :: K_Tax_draft_group, L_Tax_draft_group, C_Tax_draft_group
	REAL(DP), dimension(draft_age_category,draft_z_category) :: Tot_Income_draft_group, K_Inc_draft_group, L_Inc_draft_group
	REAL(DP), dimension(draft_age_category,draft_z_category) :: K_Inc_frac_draft_group, L_Inc_frac_draft_group
	REAL(DP), dimension(draft_age_category,draft_z_category) :: K_Tax_Inc_draft_group, L_Tax_Inc_draft_group, C_Tax_Inc_draft_group
	REAL(DP), dimension(draft_age_category,draft_z_category) :: Tax_Increase_tk_draft_group, Tax_Increase_tl_draft_group
	REAL(DP), dimension(draft_age_category,draft_z_category) :: Tax_Increase_draft_group, Tax_Rate_Increase_draft_group
	REAL(DP), dimension(draft_age_category,draft_z_category) :: Tax_Rate_Increase_tk_draft_group, Tax_Rate_Increase_tl_draft_group
	REAL(DP), dimension(draft_age_category,nz,nx) :: CE_draft_group_xz,  size_draft_group_xz, frac_pos_welfare_draft_group_xz
	REAL(DP), dimension(nz) :: DBN_Z, CDF_Z 
	REAL(DP), dimension(nz,nx) :: DBN_XZ
	INTEGER , dimension(draft_age_category+1) :: draft_age_limit
	INTEGER :: age2, z2
	REAL(DP):: K_Inc_aux, L_Inc_aux, cdf_xz, cdf_xz_low
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: Income_bench, K_Tax_bench, L_Tax_bench
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: Income_exp,   K_Tax_exp,   L_Tax_exp
	allocate( Income_bench( MaxAge,na,nz,nlambda,ne,nx) )
	allocate( K_Tax_bench(  MaxAge,na,nz,nlambda,ne,nx) )
	allocate( L_Tax_bench(  MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Income_exp(   MaxAge,na,nz,nlambda,ne,nx) )
	allocate( K_Tax_exp(    MaxAge,na,nz,nlambda,ne,nx) )
	allocate( L_Tax_exp(    MaxAge,na,nz,nlambda,ne,nx) )



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
		
		! CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction_Bench)
		! CALL COMPUTE_VALUE_FUNCTION_SPLINE  
		! ValueFunction_Bench = ValueFunction


	! Profit Matrix
		Pr_mat = Profit_Matrix(R,P)

	! Print policy functions and distribution 
		OPEN (UNIT=70, FILE=trim(Result_Folder)//'cons_by_age_z_bench', STATUS='replace')    
		OPEN (UNIT=80, FILE=trim(Result_Folder)//'leisure_by_age_z_bench', STATUS='replace')    
		OPEN (UNIT=90, FILE=trim(Result_Folder)//'dbn_by_age_z_bench', STATUS='replace')   
		DO age=1,MaxAge    
		    DO zi=1,nz
		        temp_cons_by_z(zi)       	 = sum(Cons(age,:,zi,:,:,:)*DBN_bench(age,:,zi,:,:,:))/sum(DBN_bench(age,:,zi,:,:,:))
		        temp_leisure_by_z(zi)    	 = sum((1.0_DP-HOURS(age,:,zi,:,:,:))*DBN_bench(age,:,zi,:,:,:))/sum(DBN_bench(age,:,zi,:,:,:))
		        size_by_age_z_bench(age, zi) = sum(DBN_bench(age,:,zi,:,:,:))
		    ENDDO ! zi
		    WRITE  (UNIT=70, FMT=*) temp_cons_by_z
		    WRITE  (UNIT=80, FMT=*) temp_leisure_by_z
		    WRITE  (UNIT=90, FMT=*)  size_by_age_z_bench(age, :) 
		ENDDO
		close (unit=70)
		close (unit=80)
		close (unit=90)

	! This prints Aprime for different "z" for median lambda and e
		OPEN (UNIT=50, FILE=trim(Result_Folder)//'aprime_age1_bench', STATUS='replace')     
		DO zi=1,nz 
		    WRITE  (UNIT=50, FMT=*) Aprime(1, :, zi, nlambda/2+1, ne/2+1, 1)
		ENDDO
		close (unit=50)

		OPEN (UNIT=50, FILE=trim(Result_Folder)//'aprime_age16_bench', STATUS='replace')   
		DO zi=1,nz 
		    WRITE  (UNIT=50, FMT=*) Aprime(16, :, zi, nlambda/2+1, ne/2+1, 1)
		ENDDO
		close (unit=50)  

		OPEN (UNIT=50, FILE=trim(Result_Folder)//'aprime_age31_bench', STATUS='replace')   
		DO zi=1,nz 
		    WRITE  (UNIT=50, FMT=*) Aprime(31, :, zi, nlambda/2+1, ne/2+1, 1)
		ENDDO
		close (unit=50)  

		OPEN (UNIT=50, FILE=trim(Result_Folder)//'aprime_age46_bench', STATUS='replace')   
		DO zi=1,nz 
		    WRITE  (UNIT=50, FMT=*) Aprime(46, :, zi, nlambda/2+1, ne/2+1, 1)
		ENDDO
		close (unit=50)  

	! Wealth by age group
		age_group_counter=1
		tot_wealth_by_agegroup_z_bench=0.0_DP
		size_by_agegroup_z_bench =0.0_DP
		DO age=1,MaxAge 

		    DO while (age .gt.   age_limit(age_group_counter+1) )
		        age_group_counter=age_group_counter+1
		    ENDDO    
		 
		 	DO xi=1,nx
		    DO ai=1,na
	        DO zi=1,nz
            DO lambdai=1,nlambda
            DO ei=1,ne
                size_by_agegroup_z_bench(age_group_counter,zi) = size_by_agegroup_z_bench(age_group_counter,zi) + &
                         & DBN_bench(age,ai,zi,lambdai,ei,xi)       

                tot_wealth_by_agegroup_z_bench(age_group_counter,zi) = tot_wealth_by_agegroup_z_bench(age_group_counter,zi) + &
                         & agrid(ai)*DBN_bench(age,ai,zi,lambdai,ei,xi)                           
            ENDDO
            ENDDO
	        ENDDO
		    ENDDO
		    ENDDO
		ENDDO

		OPEN (UNIT=60, FILE=trim(Result_Folder)//'mean_wealth_by_agegroup_z_bench.txt', STATUS='replace')  
		OPEN (UNIT=70, FILE=trim(Result_Folder)//'size_by_agegroup_z_bench.txt', STATUS='replace')  
		DO age_group_counter=1,max_age_category
		    WRITE  (UNIT=60, FMT=*)   tot_wealth_by_agegroup_z_bench(age_group_counter,:)/ size_by_agegroup_z_bench(age_group_counter,:)
		    WRITE  (UNIT=70, FMT=*)    size_by_agegroup_z_bench(age_group_counter,:)
		ENDDO
		close (UNIT=60)
		close (UNIT=70)
		
	age_group_counter=1

	tot_cons_by_agegroup_by_z_bench 	 = 0.0_DP
	tot_hours_by_agegroup_by_z_bench 	 = 0.0_DP
	tot_aprime_by_agegroup_by_z_bench 	 = 0.0_DP
	tot_return_by_agegroup_by_z_bench    = 0.0_DP
	tot_at_return_by_agegroup_by_z_bench = 0.0_DP
	tot_cap_tax_by_agegroup_by_z_bench 	 = 0.0_DP
	tot_lab_tax_by_agegroup_by_z_bench 	 = 0.0_DP
	tot_cons_tax_by_agegroup_by_z_bench  = 0.0_DP

	DO age=1,MaxAge 

	    DO while (age .gt.   age_limit(age_group_counter+1) )
	        age_group_counter=age_group_counter+1
	    ENDDO    
	 
	 	DO xi=1,nx 
	    DO ai=1,na
	    DO zi=1,nz
	    DO lambdai=1,nlambda
	    DO ei=1,ne

	        tot_cons_by_agegroup_by_z_bench(age_group_counter,zi) =  tot_cons_by_agegroup_by_z_bench(age_group_counter,zi)+ &
	                         & Cons(age,ai,zi,lambdai,ei,xi)*DBN_bench(age,ai,zi,lambdai,ei,xi)

	        tot_hours_by_agegroup_by_z_bench(age_group_counter,zi) = tot_hours_by_agegroup_by_z_bench(age_group_counter,zi)+&
	                         & HOURS(age,ai,zi,lambdai,ei,xi)*DBN_bench(age,ai,zi,lambdai,ei,xi)

	        tot_aprime_by_agegroup_by_z_bench(age_group_counter,zi)=tot_aprime_by_agegroup_by_z_bench(age_group_counter,zi)+&
	                         & Aprime(age,ai,zi,lambdai,ei,xi)*DBN_bench(age,ai,zi,lambdai,ei,xi)

	        tot_return_by_agegroup_by_z_bench(age_group_counter,zi) = tot_at_return_by_agegroup_by_z_bench(age_group_counter,zi)+&
	                         & DBN_bench(age,ai,zi,lambdai,ei,xi) * (R*agrid(ai) + Pr_mat(ai,zi,xi) )

	        tot_at_return_by_agegroup_by_z_bench(age_group_counter,zi) = tot_at_return_by_agegroup_by_z_bench(age_group_counter,zi)+&
	                         & DBN_bench(age,ai,zi,lambdai,ei,xi) * (Y_a(agrid(ai),zi,xi)-1.0_dp)
	       
	        tot_lab_tax_by_agegroup_by_z_bench(age_group_counter,zi) =tot_lab_tax_by_agegroup_by_z_bench(age_group_counter,zi) + &
	                & DBN_bench(age,ai,zi,lambdai,ei,xi) * ( wage_bench*eff_un(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei,xi) - &
	                & Y_h(Hours(age,ai,zi,lambdai,ei,xi),age,lambdai,ei,wage_bench) )
	      
	        tot_cons_tax_by_agegroup_by_z_bench(age_group_counter,zi)=tot_cons_tax_by_agegroup_by_z_bench(age_group_counter,zi)+&
	                         & DBN_bench(age,ai,zi,lambdai,ei,xi) * (  tauC *cons(age,ai,zi,lambdai,ei,xi) )
	         
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	ENDDO

	OPEN (UNIT=10, FILE=trim(Result_Folder)//'mean_cons_by_agegroup_by_z_bench.txt'     , STATUS='replace')  
	OPEN (UNIT=11, FILE=trim(Result_Folder)//'mean_hours_by_agegroup_by_z_bench.txt'    , STATUS='replace')  
	OPEN (UNIT=12, FILE=trim(Result_Folder)//'mean_aprime_by_agegroup_by_z_bench.txt'   , STATUS='replace')  
	OPEN (UNIT=13, FILE=trim(Result_Folder)//'mean_weighted_return_by_agegroup_by_z_bench.txt'   , STATUS='replace')  
	OPEN (UNIT=14, FILE=trim(Result_Folder)//'mean_weighted_at_return_by_agegroup_by_z_bench.txt', STATUS='replace')  
	OPEN (UNIT=15, FILE=trim(Result_Folder)//'tot_lab_tax_by_agegroup_by_z_bench.txt'  , STATUS='replace')  
	OPEN (UNIT=16, FILE=trim(Result_Folder)//'tot_cons_tax_by_agegroup_by_z_bench.txt' , STATUS='replace')  
	!OPEN (UNIT=17, FILE=trim(Result_Folder)//' tot_cap_tax_by_agegroup_by_z_bench.txt', STATUS='replace')  
	DO age_group_counter=1,max_age_category
		WRITE (UNIT=10, FMT=*) tot_cons_by_agegroup_by_z_bench(age_group_counter,:)/ size_by_agegroup_z_bench(age_group_counter,:)
		WRITE (UNIT=11, FMT=*) tot_hours_by_agegroup_by_z_bench(age_group_counter,:)/ size_by_agegroup_z_bench(age_group_counter,:)
		WRITE (UNIT=12, FMT=*) tot_aprime_by_agegroup_by_z_bench(age_group_counter,:)/ size_by_agegroup_z_bench(age_group_counter,:)
		WRITE (UNIT=13, FMT=*) tot_return_by_agegroup_by_z_bench(age_group_counter,:)/ tot_wealth_by_agegroup_z_bench(age_group_counter,:)
		WRITE (UNIT=14, FMT=*) &
			&	tot_at_return_by_agegroup_by_z_bench(age_group_counter,:)/ tot_wealth_by_agegroup_z_bench(age_group_counter,:)
		WRITE (UNIT=15, FMT=*) tot_lab_tax_by_agegroup_by_z_bench(age_group_counter,:) 
		WRITE (UNIT=16, FMT=*) tot_cons_tax_by_agegroup_by_z_bench(age_group_counter,:) 
		!WRITE (UNIT=17, FMT=*) tot_cap_tax_by_agegroup_by_z_bench(age_group_counter,:)  
	ENDDO
	close (UNIT=10)
	close (UNIT=11)
	close (UNIT=12)
	close (UNIT=13)
	close (UNIT=14)
	close (UNIT=15)
	close (UNIT=16)
	!close (UNIT=17)



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
		! CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction_Exp)
		!CALL COMPUTE_VALUE_FUNCTION_SPLINE 
		! ValueFunction_Exp = ValueFunction

	! Profit Matrix
		Pr_mat = Profit_Matrix(R,P)

	! Print policy functions and distribution 
		OPEN (UNIT=70, FILE=trim(Result_Folder)//'cons_by_age_z_exp', STATUS='replace')    
		OPEN (UNIT=80, FILE=trim(Result_Folder)//'leisure_by_age_z_exp', STATUS='replace')   
		OPEN (UNIT=90, FILE=trim(Result_Folder)//'dbn_by_age_z_exp', STATUS='replace')   
		DO age=1,MaxAge 
		    DO zi=1,nz
		          temp_cons_by_z(zi)         = sum(Cons(age,:,zi,:,:,:)*DBN1(age,:,zi,:,:,:))/sum(DBN1(age,:,zi,:,:,:))
		          temp_leisure_by_z(zi)      = sum((1.0_DP-HOURS(age,:,zi,:,:,:))*DBN1(age,:,zi,:,:,:))/sum(DBN1(age,:,zi,:,:,:))
		          size_by_age_z_exp(age, zi) = sum(DBN1(age,:,zi,:,:,:))
		    ENDDO ! zi
		    WRITE  (UNIT=70, FMT=*) temp_cons_by_z
		    WRITE  (UNIT=80, FMT=*) temp_leisure_by_z
		    WRITE  (UNIT=90, FMT=*) size_by_age_z_exp(age, :)  
		ENDDO
		close (unit=70)
		close (unit=80)
		close (unit=90)

	! This prints Aprime for different "z" for median lambda and e
		OPEN (UNIT=50, FILE=trim(Result_Folder)//'aprime_age1_exp', STATUS='replace')     
		DO zi=1,nz 
		    WRITE  (UNIT=50, FMT=*) Aprime(1, :, zi, nlambda/2+1, ne/2+1, 1)
		ENDDO
		close (unit=50)

		OPEN (UNIT=50, FILE=trim(Result_Folder)//'aprime_age16_exp', STATUS='replace')   
		DO zi=1,nz 
		    WRITE  (UNIT=50, FMT=*) Aprime(16, :, zi, nlambda/2+1, ne/2+1, 1)
		ENDDO
		close (unit=50)  

		OPEN (UNIT=50, FILE=trim(Result_Folder)//'aprime_age31_exp', STATUS='replace')   
		DO zi=1,nz 
		    WRITE  (UNIT=50, FMT=*) Aprime(31, :, zi, nlambda/2+1, ne/2+1, 1)
		ENDDO
		close (unit=50)  

		OPEN (UNIT=50, FILE=trim(Result_Folder)//'aprime_age46_exp', STATUS='replace')   
		DO zi=1,nz 
		    WRITE  (UNIT=50, FMT=*) Aprime(46, :, zi, nlambda/2+1, ne/2+1, 1)
		ENDDO
		close (unit=50)  

	! Consumption Equivalent Welfare
		OPEN (UNIT=50, FILE=trim(Result_Folder)//'CE_alternative.txt', STATUS='replace') 
		    WRITE  (UNIT=50, FMT=*) 'CE2_NB =', 100*(( sum(ValueFunction_exp(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:)) / &
		                                &  sum(ValueFunction_Bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))  ) &
		                                &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
		    WRITE  (UNIT=50, FMT=*) 'CE2_Pop =', 100*(( sum(ValueFunction_exp*DBN1) / sum(ValueFunction_Bench*DBN_bench)  ) &
		                                &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
		    WRITE  (UNIT=50, FMT=*) 'CE2_NB_bench =', 100*(( sum(ValueFunction_exp(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:)) / &
		                                &  sum(ValueFunction_Bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))  ) &
		                                &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
		    WRITE  (UNIT=50, FMT=*) 'CE2_Pop_bench =', 100*(( sum(ValueFunction_exp*DBN_bench)/sum(ValueFunction_Bench*DBN_bench)) &
		                                &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
            WRITE  (UNIT=50, FMT=*) 'CE2_NB_exp =', 100*(( sum(ValueFunction_exp(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:)) / &
		                                &  sum(ValueFunction_Bench(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))  ) &
		                                &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
		    WRITE  (UNIT=50, FMT=*) 'CE2_Pop_exp =', 100*(( sum(ValueFunction_exp*DBN1) / sum(ValueFunction_Bench*DBN1)  ) &
		                                &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
            WRITE  (UNIT=50, FMT=*) ' '

		    WRITE  (UNIT=50, FMT=*) 'Av Utility NB (bench) =',sum(ValueFunction_Bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))  
		    WRITE  (UNIT=50, FMT=*) 'Av Utility NB (exp  ) =',sum(ValueFunction_exp(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))    
		    WRITE  (UNIT=50, FMT=*) 'Av Utility  (bench)   =',sum(ValueFunction_Bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:))  
		    WRITE  (UNIT=50, FMT=*) 'Av Utility  (exp)     =',sum(ValueFunction_exp(:,:,:,:,:,:)*DBN1(:,:,:,:,:,:))
		    print*,' '
		    print*, 'Av Utility NB (bench) =',sum(ValueFunction_Bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))  
		    print*, 'Av Utility NB (exp  ) =',sum(ValueFunction_exp(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))    
		    print*, sum(ValueFunction_exp(1,:,:,:,:,:)), sum(DBN1(1,:,:,:,:,:))    
		close (unit=50)		


		OPEN (UNIT=50, FILE=trim(Result_Folder)//'CE_NEWBORN', STATUS='replace')  
		OPEN (UNIT=60, FILE=trim(Result_Folder)//'CE', STATUS='replace')  
		OPEN (UNIT=70, FILE=trim(Result_Folder)//'CE_by_age', STATUS='replace')  
		OPEN (UNIT=80, FILE=trim(Result_Folder)//'CE_by_age_z', STATUS='replace')  

		DO age=1,MaxAge
			if (Log_Switch.eqv..true.) then 
		    	Cons_Eq_Welfare(age,:,:,:,:,:)= & 
		    		& exp((ValueFunction_exp(age,:,:,:,:,:)-ValueFunction_Bench(age,:,:,:,:,:))/CumDiscountF(age))-1.0_DP
		    else 
		    	Cons_Eq_Welfare(age,:,:,:,:,:)=(ValueFunction_exp(age,:,:,:,:,:)/ValueFunction_Bench(age,:,:,:,:,:)) &
                                				&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP
		    end if 

		    WRITE  (UNIT=70, FMT=*) 100*sum(Cons_Eq_Welfare(age,:,:,:,:,:)*DBN_bench(age,:,:,:,:,:))/sum(DBN_bench(age,:,:,:,:,:))
		    DO zi=1,nz
		         temp_ce_by_z(zi) = 100*sum(Cons_Eq_Welfare(age,:,zi,:,:,:)*DBN_bench(age,:,zi,:,:,:))/sum(DBN_bench(age,:,zi,:,:,:))
		    ENDDO
		    WRITE  (UNIT=80, FMT=*) temp_ce_by_z
		    !print*,'age=',age, temp_ce_by_z, ', mean:  ', &
		    !    & 100*sum(Cons_Eq_Welfare(age,:,:,:,:)*DBN_bench(age,:,:,:,:))/sum(DBN_bench(age,:,:,:,:))
		ENDDO

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
		OPEN (UNIT=80, FILE=trim(Result_Folder)//'CE_by_AgeGroup_z', STATUS='replace') 
		DO zi=1,nz
		    DO age_group_counter=1,max_age_category
		         CE_by_agegroup_z(age_group_counter,zi)= &
		            & 100*sum(Cons_Eq_Welfare(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:,:)* &
		            &                         DBN_bench(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:,:))/&
		            &                sum( DBN_bench(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:,:))
		    ENDDO
		ENDDO
		DO age_group_counter=1,max_age_category
		    WRITE  (UNIT=80, FMT=*)  CE_by_agegroup_z(age_group_counter,:)
		ENDDO
		close (unit=80)

		OPEN (UNIT=80, FILE=trim(Result_Folder)//'CE_by_Asset_z_med_E_Lambda_age1', STATUS='replace') 
		OPEN (UNIT=90, FILE=trim(Result_Folder)//'MeanAsset_by_z_med_E_Lambda_age1', STATUS='replace') 
		DO zi=1,nz
		    WRITE  (UNIT=80, FMT=*)   Cons_Eq_Welfare(1,:,zi,nlambda/2+1, ne/2+1, 1)
		    WRITE  (UNIT=90, FMT=*)   sum(agrid*DBN_bench(1,:,zi,nlambda/2+1, ne/2+1, 1)) &
		                                                        & /sum(DBN_bench(1,:,zi,nlambda/2+1, ne/2+1, 1))
		ENDDO           
		close (unit=80)
		close (unit=90)

		OPEN (UNIT=80, FILE=trim(Result_Folder)//'CE_by_Asset_z_med_E_Lambda_age16', STATUS='replace') 
		OPEN (UNIT=90, FILE=trim(Result_Folder)//'MeanAsset_by_z_med_E_Lambda_age16', STATUS='replace') 
		DO zi=1,nz
		    WRITE  (UNIT=80, FMT=*)   Cons_Eq_Welfare(16,:,zi,nlambda/2+1, ne/2+1, 1)
		    WRITE  (UNIT=90, FMT=*)   sum(agrid*DBN_bench(16,:,zi,nlambda/2+1, ne/2+1, 1)) &
		                                                        & /sum(DBN_bench(16,:,zi,nlambda/2+1, ne/2+1, 1))
		ENDDO           
		close (unit=80)
		close (unit=90)

		OPEN (UNIT=80, FILE=trim(Result_Folder)//'CE_by_Asset_z_med_E_Lambda_age31', STATUS='replace') 
		OPEN (UNIT=90, FILE=trim(Result_Folder)//'MeanAsset_by_z_med_E_Lambda_age31', STATUS='replace') 
		DO zi=1,nz
			WRITE  (UNIT=80, FMT=*)   Cons_Eq_Welfare(31,:,zi,nlambda/2+1, ne/2+1, 1)
			WRITE  (UNIT=90, FMT=*)   sum(agrid*DBN_bench(31,:,zi,nlambda/2+1, ne/2+1, 1)) &
			                                            & /sum(DBN_bench(31,:,zi,nlambda/2+1, ne/2+1, 1))
		ENDDO           
		close (unit=8)
		close (unit=9)


		! FRACTION POSITIVE WELFARE BY AGE-Z GROUP
		frac_pos_welfare=0.0_DP
		size_pos_welfare_by_age_z=0.0_DP

		size_by_agegroup_z_exp = 0.0_DP
		size_pos_welfare_by_agegroup_z = 0.0_DP

		tot_wealth_by_agegroup_z_exp = 0.0_DP

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
                size_by_agegroup_z_exp(age_group_counter,zi) = size_by_agegroup_z_exp(age_group_counter,zi) + &
                         & DBN1(age,ai,zi,lambdai,ei,xi)       

                tot_wealth_by_agegroup_z_exp(age_group_counter,zi) = tot_wealth_by_agegroup_z_exp(age_group_counter,zi) + &
                         & agrid(ai)*DBN1(age,ai,zi,lambdai,ei,xi)                           
		                    
	        ENDDO
	        ENDDO
	        ENDDO
		    ENDDO
		    ENDDO
		ENDDO


	OPEN (UNIT=60, FILE=trim(Result_Folder)//'mean_wealth_by_agegroup_z_exp.txt', STATUS='replace')  
	OPEN (UNIT=70, FILE=trim(Result_Folder)//'size_by_agegroup_z_exp.txt', STATUS='replace')  
	DO age_group_counter=1,max_age_category
	    WRITE  (UNIT=60, FMT=*)   tot_wealth_by_agegroup_z_exp(age_group_counter,:)/ size_by_agegroup_z_exp(age_group_counter,:)
	    WRITE  (UNIT=70, FMT=*)   size_by_agegroup_z_exp(age_group_counter,:)
	ENDDO
	close (UNIT=60)
	close (UNIT=70)

	OPEN (UNIT=60, FILE=trim(Result_Folder)//'frac_pos_welfare_by_agegroup_z', STATUS='replace')  
	DO age_group_counter=1,max_age_category
	    WRITE  (UNIT=60, FMT=*)  size_pos_welfare_by_agegroup_z(age_group_counter,:)/ size_by_agegroup_z_bench(age_group_counter,:)
	ENDDO
	close (UNIT=60)

	OPEN (UNIT=60, FILE=trim(Result_Folder)//'frac_pos_welfare', STATUS='replace')  
	WRITE  (UNIT=60, FMT=*) frac_pos_welfare
	close (unit=60)

	OPEN (UNIT=60, FILE=trim(Result_Folder)//'frac_pos_welfare_by_age_z', STATUS='replace')  
	DO age=1, MaxAge
	    WRITE  (UNIT=60, FMT=*) size_pos_welfare_by_age_z(age,:)/size_by_age_z_bench(age,:)
	ENDDO
	close (UNIT=60)


	age_group_counter=1

	tot_cons_by_agegroup_by_z_exp 	 	= 0.0_DP
	tot_hours_by_agegroup_by_z_exp 	 	= 0.0_DP
	tot_aprime_by_agegroup_by_z_exp 	= 0.0_DP
	tot_return_by_agegroup_by_z_exp    	= 0.0_DP
	tot_at_return_by_agegroup_by_z_exp 	= 0.0_DP
	tot_cap_tax_by_agegroup_by_z_exp 	= 0.0_DP
	tot_lab_tax_by_agegroup_by_z_exp 	= 0.0_DP
	tot_cons_tax_by_agegroup_by_z_exp  	= 0.0_DP

	DO age=1,MaxAge 

	    DO while (age .gt.   age_limit(age_group_counter+1) )
	        age_group_counter=age_group_counter+1
	    ENDDO    
	 
	 	DO xi=1,nx
	    DO ai=1,na
	    DO zi=1,nz
	    DO lambdai=1,nlambda
	    DO ei=1,ne

	        tot_cons_by_agegroup_by_z_exp(age_group_counter,zi) =  tot_cons_by_agegroup_by_z_exp(age_group_counter,zi)+ &
	                         & Cons(age,ai,zi,lambdai,ei,xi)*DBN1(age,ai,zi,lambdai,ei,xi)

	        tot_hours_by_agegroup_by_z_exp(age_group_counter,zi) = tot_hours_by_agegroup_by_z_exp(age_group_counter,zi)+&
	                         & HOURS(age,ai,zi,lambdai,ei,xi)*DBN1(age,ai,zi,lambdai,ei,xi)

	        tot_aprime_by_agegroup_by_z_exp(age_group_counter,zi)=tot_aprime_by_agegroup_by_z_exp(age_group_counter,zi)+&
	                         & Aprime(age,ai,zi,lambdai,ei,xi)*DBN1(age,ai,zi,lambdai,ei,xi)

	        tot_return_by_agegroup_by_z_exp(age_group_counter,zi) = tot_at_return_by_agegroup_by_z_exp(age_group_counter,zi)+&
	                         & DBN1(age,ai,zi,lambdai,ei,xi) * (R*agrid(ai) + Pr_mat(ai,zi,xi) )

	        tot_at_return_by_agegroup_by_z_exp(age_group_counter,zi) = tot_at_return_by_agegroup_by_z_exp(age_group_counter,zi)+&
	                         & DBN1(age,ai,zi,lambdai,ei,xi) * (Y_a(agrid(ai),zi,xi)-1.0_dp)
	       
	        tot_lab_tax_by_agegroup_by_z_exp(age_group_counter,zi) = tot_lab_tax_by_agegroup_by_z_exp(age_group_counter,zi) + &
	                         & DBN1(age,ai,zi,lambdai,ei,xi) * ( wage_exp * eff_un(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei,xi) - &
	                          & Y_h(Hours(age,ai,zi,lambdai,ei,xi),age,lambdai,ei,wage_exp) )
	      
	        tot_cons_tax_by_agegroup_by_z_exp(age_group_counter,zi)=tot_cons_tax_by_agegroup_by_z_exp(age_group_counter,zi)+&
	                         & DBN1(age,ai,zi,lambdai,ei,xi) * (  tauC *cons(age,ai,zi,lambdai,ei,xi) )
	         
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	ENDDO

	OPEN (UNIT=10, FILE=trim(Result_Folder)//'mean_cons_by_agegroup_by_z_exp.txt'     , STATUS='replace')  
	OPEN (UNIT=11, FILE=trim(Result_Folder)//'mean_hours_by_agegroup_by_z_exp.txt'    , STATUS='replace')  
	OPEN (UNIT=12, FILE=trim(Result_Folder)//'mean_aprime_by_agegroup_by_z_exp.txt'   , STATUS='replace')  
	OPEN (UNIT=13, FILE=trim(Result_Folder)//'mean_weighted_return_by_agegroup_by_z_exp.txt'   , STATUS='replace')  
	OPEN (UNIT=14, FILE=trim(Result_Folder)//'mean_weighted_at_return_by_agegroup_by_z_exp.txt', STATUS='replace')  
	OPEN (UNIT=15, FILE=trim(Result_Folder)//'tot_lab_tax_by_agegroup_by_z_exp.txt'  , STATUS='replace')  
	OPEN (UNIT=16, FILE=trim(Result_Folder)//'tot_cons_tax_by_agegroup_by_z_exp.txt' , STATUS='replace')  
	!OPEN (UNIT=17, FILE=trim(Result_Folder)//' tot_cap_tax_by_agegroup_by_z_exp.txt', STATUS='replace')  
	DO age_group_counter=1,max_age_category
		WRITE (UNIT=10, FMT=*) tot_cons_by_agegroup_by_z_exp(age_group_counter,:)/ size_by_agegroup_z_exp(age_group_counter,:)
		WRITE (UNIT=11, FMT=*) tot_hours_by_agegroup_by_z_exp(age_group_counter,:)/ size_by_agegroup_z_exp(age_group_counter,:)
		WRITE (UNIT=12, FMT=*) tot_aprime_by_agegroup_by_z_exp(age_group_counter,:)/ size_by_agegroup_z_exp(age_group_counter,:)
		WRITE (UNIT=13, FMT=*) tot_return_by_agegroup_by_z_exp(age_group_counter,:)/ tot_wealth_by_agegroup_z_exp(age_group_counter,:)
		WRITE (UNIT=14, FMT=*) tot_at_return_by_agegroup_by_z_exp(age_group_counter,:)/ tot_wealth_by_agegroup_z_exp(age_group_counter,:)
		WRITE (UNIT=15, FMT=*) tot_lab_tax_by_agegroup_by_z_exp(age_group_counter,:) 
		WRITE (UNIT=16, FMT=*) tot_cons_tax_by_agegroup_by_z_exp(age_group_counter,:) 
		!WRITE (UNIT=17, FMT=*) tot_cap_tax_by_agegroup_by_z_exp(age_group_counter,:)  
	ENDDO
	close (UNIT=10)
	close (UNIT=11)
	close (UNIT=12)
	close (UNIT=13)
	close (UNIT=14)
	close (UNIT=15)
	close (UNIT=16)
	!close (UNIT=17)



	! Compute average welfare
		Welfare_Gain_Pop_bench = 100.0_DP*sum(Cons_Eq_Welfare*DBN_bench)
		Welfare_Gain_Pop_exp   = 100.0_DP*sum(Cons_Eq_Welfare*DBN1)
		Welfare_Gain_NB_bench  = 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
		Welfare_Gain_NB_exp    = 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))

	print*,'---------------------------'
	print*,''
	print*,'Average Welfare Gain Whole Population (bench dbn) (prct)=',100.0_DP*sum(Cons_Eq_Welfare*DBN_bench)
	print*,'Average Welfare Gain Whole Population (exp dbn)     (prct)=',100.0_DP*sum(Cons_Eq_Welfare*DBN1)
	print*,'Average Welfare Gain NEW BORN (bench dbn) (prct)          =',&
	    & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
	print*,'Average Welfare Gain NEW BORN (exp dbn)     (prct)          =',&
	    & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:))/sum(DBN1(1,:,:,:,:,:))
	print*,''
	print*,'---------------------------'

	! Compute Average Utility - CE2

	Av_Util_NB  =  100.0_dp * (( sum(ValueFunction_exp(1,:,:,:,:,:)*DBN1(1,:,:,:,:,:)) / &
				&               sum(ValueFunction_Bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:)) ) &
				& ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

	Av_Util_Pop =  100.0_dp * (( sum(ValueFunction_exp(:,:,:,:,:,:)*DBN1(:,:,:,:,:,:)) / &
				& 	            sum(ValueFunction_Bench(:,:,:,:,:,:)*DBN_bench(:,:,:,:,:,:)) ) &
				& ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Draft Tables
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	frac_pos_welfare_draft_group_z = 0.0_dp 
	wealth_draft_group_z           = 0.0_dp 

	do zi  = 1,nz
	do age = 1,draft_age_category
		size_draft_group_z(age,zi) = &
			& sum(DBN_bench(draft_age_limit(age)+1:draft_age_limit(age+1),:,zi,:,:,:))

		CE_draft_group_z(age,zi) =  100* &
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
        	wealth_draft_group_z(age,zi) = wealth_draft_group_z(age,zi) + agrid(ai)*DBN_bench(age2,ai,zi,lambdai,ei,xi)
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
	size_draft_group(:,1) = size_draft_group_z(:,1) + size_draft_group_z(:,2) + size_draft_group_z(:,3) & 
							&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*size_draft_group_z(:,4)
	size_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*size_draft_group_z(:,4)+((0.80_dp-CDF_Z(4))/DBN_Z(5))*size_draft_group_z(:,5)
	size_draft_group(:,3) = (0.10_dp/DBN_Z(5))*size_draft_group_z(:,5)
	size_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*size_draft_group_z(:,5) + &
						& ((0.99_dp-CDF_Z(5))/DBN_Z(6))*size_draft_group_z(:,6) 
	size_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*size_draft_group_z(:,6) + & 
						& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*size_draft_group_z(:,7) 
	size_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*size_draft_group_z(:,7) + size_draft_group_z(:,8) + size_draft_group_z(:,9) 
	size_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*size_draft_group_z(:,7) + &
							& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*size_draft_group_z(:,8)
	size_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*size_draft_group_z(:,8) + size_draft_group_z(:,9)

	! CE of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
	CE_draft_group(:,1)   = ( DBN_Z(1)*CE_draft_group_z(:,1) + DBN_Z(2)*CE_draft_group_z(:,2) + & 
							& DBN_Z(3)*CE_draft_group_z(:,3) + (0.40_dp-CDF_Z(3))*CE_draft_group_z(:,4) )/0.40_dp
	CE_draft_group(:,2)   = ( (CDF_Z(4)-0.40_dp)*CE_draft_group_z(:,4) + (0.80_dp-CDF_Z(4))*CE_draft_group_z(:,5) )/0.40_dp
	CE_draft_group(:,3)   = CE_draft_group_z(:,5)
	CE_draft_group(:,4)   = ( (CDF_Z(5)-0.90_dp)*CE_draft_group_z(:,5) + (0.99_dp-CDF_Z(5))*CE_draft_group_z(:,6) )/0.09_dp
	CE_draft_group(:,5)   = ( (CDF_Z(6)-0.99_dp)*CE_draft_group_z(:,6) + (0.999_dp-CDF_Z(6))*CE_draft_group_z(:,7) )/0.009_dp
	CE_draft_group(:,6)   = ( (CDF_Z(7)-0.999_dp)*CE_draft_group_z(:,7) + &
							&  DBN_Z(8)*CE_draft_group_z(:,8) + DBN_Z(9)*CE_draft_group_z(:,9) )/0.001_dp
	CE_draft_group(:,7)   = ( (CDF_Z(7)-0.999_dp)*CE_draft_group_z(:,7) + (0.9999_dp-CDF_Z(7))*CE_draft_group_z(:,8) )/0.0009_dp
	CE_draft_group(:,8)   = ( (CDF_Z(8)-0.9999_dp)*CE_draft_group_z(:,8) + DBN_Z(9)*CE_draft_group_z(:,9) )/0.0001_dp

	! Frac. pos. welfare by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
	frac_pos_welfare_draft_group(:,1) = frac_pos_welfare_draft_group_z(:,1) + frac_pos_welfare_draft_group_z(:,2) + &
							&  frac_pos_welfare_draft_group_z(:,3) + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*frac_pos_welfare_draft_group_z(:,4)
	frac_pos_welfare_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*frac_pos_welfare_draft_group_z(:,4) + & 
							&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*frac_pos_welfare_draft_group_z(:,5)
	frac_pos_welfare_draft_group(:,3) = (0.10_dp/DBN_Z(5))*frac_pos_welfare_draft_group_z(:,5)
	frac_pos_welfare_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*frac_pos_welfare_draft_group_z(:,5) + & 
							& ((0.99_dp-CDF_Z(5))/DBN_Z(6))*frac_pos_welfare_draft_group_z(:,6) 
	frac_pos_welfare_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*frac_pos_welfare_draft_group_z(:,6) + &
							& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*frac_pos_welfare_draft_group_z(:,7) 
	frac_pos_welfare_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*frac_pos_welfare_draft_group_z(:,7) + &
							& frac_pos_welfare_draft_group_z(:,8) + frac_pos_welfare_draft_group_z(:,9) 
	frac_pos_welfare_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*frac_pos_welfare_draft_group_z(:,7) + &
							& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*frac_pos_welfare_draft_group_z(:,8)
	frac_pos_welfare_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*frac_pos_welfare_draft_group_z(:,8) + &
							&  frac_pos_welfare_draft_group_z(:,9)

	! Size of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
	wealth_draft_group(:,1) = wealth_draft_group_z(:,1) + wealth_draft_group_z(:,2) + wealth_draft_group_z(:,3) & 
							&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*wealth_draft_group_z(:,4)
	wealth_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*wealth_draft_group_z(:,4) +& 
							&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*wealth_draft_group_z(:,5)
	wealth_draft_group(:,3) = (0.10_dp/DBN_Z(5))*wealth_draft_group_z(:,5)
	wealth_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*wealth_draft_group_z(:,5) + & 
							&  ((0.99_dp-CDF_Z(5))/DBN_Z(6))*wealth_draft_group_z(:,6)
	wealth_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*wealth_draft_group_z(:,6) + & 
							& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*wealth_draft_group_z(:,7) 
	wealth_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*wealth_draft_group_z(:,7) +&
							&  wealth_draft_group_z(:,8) + wealth_draft_group_z(:,9) 
	wealth_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*wealth_draft_group_z(:,7) + &
							& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*wealth_draft_group_z(:,8)
	wealth_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*wealth_draft_group_z(:,8) + wealth_draft_group_z(:,9)

	! Fix fractions
	frac_pos_welfare_draft_group = 100*frac_pos_welfare_draft_group/size_draft_group
    av_wealth_draft_group        = (EBAR_data/(EBAR_bench*0.727853584919652_dp))*wealth_draft_group/size_draft_group
    frac_wealth_draft_group      = 100*wealth_draft_group/sum(wealth_draft_group)

    OPEN (UNIT=80, FILE=trim(Result_Folder)//'draft_group_CE.txt', STATUS='replace') 
    OPEN (UNIT=81, FILE=trim(Result_Folder)//'draft_group_size.txt', STATUS='replace') 
    OPEN (UNIT=82, FILE=trim(Result_Folder)//'draft_group_fpos_welfare.txt', STATUS='replace') 
    OPEN (UNIT=83, FILE=trim(Result_Folder)//'draft_group_wealth.txt', STATUS='replace') 
    OPEN (UNIT=84, FILE=trim(Result_Folder)//'draft_group_av_wealth.txt', STATUS='replace') 
    OPEN (UNIT=85, FILE=trim(Result_Folder)//'draft_group_frac_wealth.txt', STATUS='replace') 
	do age = 1,draft_age_category
	    WRITE  (UNIT=80, FMT=*)  CE_draft_group(age,:)
	    WRITE  (UNIT=81, FMT=*)  size_draft_group(age,:)
	    WRITE  (UNIT=82, FMT=*)  frac_pos_welfare_draft_group(age,:)
	    WRITE  (UNIT=83, FMT=*)  wealth_draft_group(age,:)
	    WRITE  (UNIT=84, FMT=*)  av_wealth_draft_group(age,:)
	    WRITE  (UNIT=85, FMT=*)  frac_wealth_draft_group(age,:)
	ENDDO
	close(unit=80)
	close(unit=81)
	close(unit=82)
	close(unit=83)
	close(unit=84)
	close(unit=85)

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
	!! Draft Tables by xz producitivity groups
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	frac_pos_welfare_draft_group_xz = 0.0_dp 
	size_draft_group_xz = 0.0_dp
	DBN_XZ = sum(sum(sum(sum(DBN_bench,5),4),2),1) 

	do xi  = 1,nx
	do zi  = 1,nz
	do age = 1,draft_age_category
		size_draft_group_xz(age,zi,xi) = &
			& sum(DBN_bench(draft_age_limit(age)+1:draft_age_limit(age+1),:,zi,:,:,xi))

		CE_draft_group_xz(age,zi,xi) =  100* &
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

	! CE of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
    	cdf_xz = sum(DBN_XZ(1:4,3))
	CE_draft_group(:,1)   = ( DBN_XZ(1,3)*CE_draft_group_xz(:,1,3) + DBN_XZ(2,3)*CE_draft_group_xz(:,2,3) + &
							& DBN_XZ(3,3)*CE_draft_group_xz(:,3,3) + DBN_XZ(4,3)*CE_draft_group_xz(:,4,3) + &
							& (0.40_dp-cdf_xz)*CE_draft_group_xz(:,5,3) )/0.40_dp 
		cdf_xz = cdf_xz + DBN_XZ(5,3)
		cdf_xz_low = cdf_xz 
		cdf_xz = cdf_xz + sum(DBN_XZ(6:,3)) + sum(DBN_XZ(1:3,1))
	CE_draft_group(:,2)   = ( (cdf_xz_low-0.40_dp)*CE_draft_group_xz(:,5,3) + &
							& DBN_XZ(6,3)*CE_draft_group_xz(:,6,3) + DBN_XZ(7,3)*CE_draft_group_xz(:,7,3) + &
							& DBN_XZ(8,3)*CE_draft_group_xz(:,8,3) + DBN_XZ(9,3)*CE_draft_group_xz(:,9,3) + &
							& DBN_XZ(1,1)*CE_draft_group_xz(:,1,1) + DBN_XZ(2,1)*CE_draft_group_xz(:,2,1) + &
							& DBN_XZ(3,1)*CE_draft_group_xz(:,3,1) + &
							& (0.80_dp-cdf_xz)*CE_draft_group_xz(:,4,1) )/0.40_dp
		cdf_xz = cdf_xz + DBN_XZ(4,1)
	CE_draft_group(:,3)   = ( (cdf_xz-0.80_dp)*CE_draft_group_xz(:,4,1) + (0.90_dp-cdf_xz)*CE_draft_group_xz(:,5,2) )/0.10_dp
		cdf_xz = cdf_xz + DBN_XZ(5,2)
		cdf_xz_low = cdf_xz 
		cdf_xz = cdf_xz + sum(DBN_XZ(6:,2)) + DBN_XZ(5,1)
	CE_draft_group(:,4)   = ( (cdf_xz_low-0.90_dp)*CE_draft_group_xz(:,5,2) + & 
							& DBN_XZ(6,2)*CE_draft_group_xz(:,6,2) + DBN_XZ(7,2)*CE_draft_group_xz(:,7,2) + &
							& DBN_XZ(8,2)*CE_draft_group_xz(:,8,2) + DBN_XZ(9,2)*CE_draft_group_xz(:,9,2) + &
							& DBN_XZ(5,1)*CE_draft_group_xz(:,5,1) + &
							& (0.99_dp-cdf_xz)*CE_draft_group_xz(:,6,1) )/0.09_dp
		cdf_xz = cdf_xz + DBN_XZ(6,1)
	CE_draft_group(:,5)   = ( (cdf_xz-0.99_dp)*CE_draft_group_xz(:,6,1) + (0.999_dp-cdf_xz)*CE_draft_group_xz(:,7,1) )/0.009_dp
		cdf_xz = cdf_xz + DBN_XZ(7,1)
	CE_draft_group(:,6)   = ( (CDF_xz-0.999_dp)*CE_draft_group_xz(:,7,1) + &
							&  DBN_XZ(8,1)*CE_draft_group_xz(:,8,1) + DBN_XZ(9,1)*CE_draft_group_xz(:,9,1) )/0.001_dp
	CE_draft_group(:,7)   = 0
	CE_draft_group(:,8)   = 0

		! Adjustment for first age group
		    ! % 40%-80%    of Current Productivity 
		    !     % (z1,x1), (z2,x1), (z3,x1), (z4,x1)
	    	CE_draft_group(1,2)   = ( DBN_XZ(1,1)*CE_draft_group_xz(1,1,1) + DBN_XZ(2,1)*CE_draft_group_xz(1,2,1) + &
						& DBN_XZ(3,1)*CE_draft_group_xz(1,3,1) + DBN_XZ(4,1)*CE_draft_group_xz(1,4,1) )/sum(DBN_XZ(1:4,1))
		    ! % 80%-90%    of Current Productivity
		    !     % (z4,x1)
		    CE_draft_group(1,3)   = CE_draft_group_xz(1,4,1)
		    ! % 90%-99%    of Current Productivity 
		    !     % (z5,x1), (z6,x1)
		    CE_draft_group(1,4)   = ( DBN_XZ(5,1)*CE_draft_group_xz(1,5,1) + DBN_XZ(6,1)*CE_draft_group_xz(1,6,1) )/sum(DBN_XZ(5:6,1))

	! Size of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
    	cdf_xz = sum(DBN_XZ(1:4,3))
	size_draft_group(:,1)   = size_draft_group_xz(:,1,3) + size_draft_group_xz(:,2,3) + &
							& size_draft_group_xz(:,3,3) + size_draft_group_xz(:,4,3) + &
							& (0.40_dp-cdf_xz)*size_draft_group_xz(:,5,3)/DBN_XZ(5,3) 
		cdf_xz = cdf_xz + DBN_XZ(5,3)
		cdf_xz_low = cdf_xz 
		cdf_xz = cdf_xz + sum(DBN_XZ(6:,3)) + sum(DBN_XZ(1:3,1))
	size_draft_group(:,2)   = (cdf_xz_low-0.40_dp)*size_draft_group_xz(:,5,3)/DBN_XZ(5,3) + &
							& size_draft_group_xz(:,6,3) + size_draft_group_xz(:,7,3) + &
							& size_draft_group_xz(:,8,3) + size_draft_group_xz(:,9,3) + &
							& size_draft_group_xz(:,1,1) + size_draft_group_xz(:,2,1) + &
							& size_draft_group_xz(:,3,1) + &
							& (0.80_dp-cdf_xz)*size_draft_group_xz(:,4,1)/DBN_XZ(4,1)
		cdf_xz = cdf_xz + DBN_XZ(4,1)
	size_draft_group(:,3)   = (cdf_xz-0.80_dp)*size_draft_group_xz(:,4,1)/DBN_XZ(4,1) + &
							& (0.90_dp-cdf_xz)*size_draft_group_xz(:,5,2)/DBN_XZ(5,2)
		cdf_xz = cdf_xz + DBN_XZ(5,2)
		cdf_xz_low = cdf_xz 
		cdf_xz = cdf_xz + sum(DBN_XZ(6:,2)) + DBN_XZ(5,1)
	size_draft_group(:,4)   = (cdf_xz_low-0.90_dp)*size_draft_group_xz(:,5,2)/DBN_XZ(5,2) + & 
							& size_draft_group_xz(:,6,2) + size_draft_group_xz(:,7,2) + &
							& size_draft_group_xz(:,8,2) + size_draft_group_xz(:,9,2) + &
							& size_draft_group_xz(:,5,1) + &
							& (0.99_dp-cdf_xz)*size_draft_group_xz(:,6,1)/DBN_XZ(6,1)
		cdf_xz = cdf_xz + DBN_XZ(6,1)
	size_draft_group(:,5)   = (cdf_xz-0.99_dp)*size_draft_group_xz(:,6,1)/DBN_XZ(6,1) + &
							& (0.999_dp-cdf_xz)*size_draft_group_xz(:,7,1)/DBN_XZ(7,1)
		cdf_xz = cdf_xz + DBN_XZ(7,1)
	size_draft_group(:,6)   = (CDF_xz-0.999_dp)*size_draft_group_xz(:,7,1)/DBN_XZ(7,1) + &
							&  size_draft_group_xz(:,8,1) + size_draft_group_xz(:,9,1)
	size_draft_group(:,7)   = 0
	size_draft_group(:,8)   = 0


	! Frac. pos. welfare by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
    	cdf_xz = sum(DBN_XZ(1:4,3))
	frac_pos_welfare_draft_group(:,1)   = frac_pos_welfare_draft_group_xz(:,1,3) + frac_pos_welfare_draft_group_xz(:,2,3) + &
							& frac_pos_welfare_draft_group_xz(:,3,3) + frac_pos_welfare_draft_group_xz(:,4,3) + &
							& (0.40_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(:,5,3)/DBN_XZ(5,3) 
		cdf_xz = cdf_xz + DBN_XZ(5,3)
		cdf_xz_low = cdf_xz 
		cdf_xz = cdf_xz + sum(DBN_XZ(6:,3)) + sum(DBN_XZ(1:3,1))
	frac_pos_welfare_draft_group(:,2)   = (cdf_xz_low-0.40_dp)*frac_pos_welfare_draft_group_xz(:,5,3)/DBN_XZ(5,3) + &
							& frac_pos_welfare_draft_group_xz(:,6,3) + frac_pos_welfare_draft_group_xz(:,7,3) + &
							& frac_pos_welfare_draft_group_xz(:,8,3) + frac_pos_welfare_draft_group_xz(:,9,3) + &
							& frac_pos_welfare_draft_group_xz(:,1,1) + frac_pos_welfare_draft_group_xz(:,2,1) + &
							& frac_pos_welfare_draft_group_xz(:,3,1) + &
							& (0.80_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(:,4,1)/DBN_XZ(4,1)
		cdf_xz = cdf_xz + DBN_XZ(4,1)
	frac_pos_welfare_draft_group(:,3)   = (cdf_xz-0.80_dp)*frac_pos_welfare_draft_group_xz(:,4,1)/DBN_XZ(4,1) + &
										& (0.90_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(:,5,2)/DBN_XZ(5,2)
		cdf_xz = cdf_xz + DBN_XZ(5,2)
		cdf_xz_low = cdf_xz 
		cdf_xz = cdf_xz + sum(DBN_XZ(6:,2)) + DBN_XZ(5,1)
	frac_pos_welfare_draft_group(:,4)   = (cdf_xz_low-0.90_dp)*frac_pos_welfare_draft_group_xz(:,5,2)/DBN_XZ(5,2) + & 
							& frac_pos_welfare_draft_group_xz(:,6,2) + frac_pos_welfare_draft_group_xz(:,7,2) + &
							& frac_pos_welfare_draft_group_xz(:,8,2) + frac_pos_welfare_draft_group_xz(:,9,2) + &
							& frac_pos_welfare_draft_group_xz(:,5,1) + &
							& (0.99_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(:,6,1)/DBN_XZ(6,1)
		cdf_xz = cdf_xz + DBN_XZ(6,1)
	frac_pos_welfare_draft_group(:,5)   = (cdf_xz-0.99_dp)*frac_pos_welfare_draft_group_xz(:,6,1)/DBN_XZ(6,1) + &
										& (0.999_dp-cdf_xz)*frac_pos_welfare_draft_group_xz(:,7,1)/DBN_XZ(7,1)
		cdf_xz = cdf_xz + DBN_XZ(7,1)
	frac_pos_welfare_draft_group(:,6)   = (CDF_xz-0.999_dp)*frac_pos_welfare_draft_group_xz(:,7,1)/DBN_XZ(7,1) + &
							&  frac_pos_welfare_draft_group_xz(:,8,1) + frac_pos_welfare_draft_group_xz(:,9,1)
	frac_pos_welfare_draft_group(:,7)   = 0
	frac_pos_welfare_draft_group(:,8)   = 0


	! Fix fractions
	frac_pos_welfare_draft_group = 100*frac_pos_welfare_draft_group/size_draft_group

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
	frac_pos_welfare_draft_group = 100*frac_pos_welfare_draft_group/size_draft_group

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
	!! Draft Tables for Tax Burden
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!! Benchmark 
		Pr_mat = Profit_Matrix(R_bench,P_bench)
		CALL ComputeLaborUnits(EBAR_bench,wage_bench)

		K_Tax_draft_group_z			= 0.0_dp
		L_Tax_draft_group_z			= 0.0_dp
		C_Tax_draft_group_z			= 0.0_dp
		Tot_Income_draft_group_z	= 0.0_dp
		K_Inc_draft_group_z         = 0.0_dp
		L_Inc_draft_group_z         = 0.0_dp
		K_Inc_frac_draft_group_z    = 0.0_dp
		L_Inc_frac_draft_group_z    = 0.0_dp
		K_Tax_Inc_draft_group_z		= 0.0_dp
		L_Tax_Inc_draft_group_z		= 0.0_dp
		C_Tax_Inc_draft_group_z		= 0.0_dp
		wealth_draft_group_z   		= 0.0_dp 

		do zi  = 1,nz
		do age = 1,draft_age_category
	        do xi=1,nx
		    do ei=1,ne
		    do lambdai=1,nlambda
		    do ai=1,na
		    do age2=draft_age_limit(age)+1,draft_age_limit(age+1)

		    	if (age.lt.draft_age_category) then
		    	L_Inc_aux   = yh(age2,lambdai,ei)*Hours_bench(age2,ai,zi,lambdai,ei,xi)
		    	else
		    	L_Inc_aux   = RetY_lambda_e(lambdai,ei) 
		    	endif 
		    	K_Inc_aux   = R_bench*agrid(ai) + Pr_mat(ai,zi,xi)

		    	Income_bench(age2,ai,zi,lambdai,ei,xi) = L_Inc_aux + K_Inc_aux

		    	K_Tax_bench(age2,ai,zi,lambdai,ei,xi)  = tauK_bench*( K_Inc_aux )
		    	if (age.lt.draft_age_category) then
		    	L_Tax_bench(age2,ai,zi,lambdai,ei,xi)  = L_Inc_aux - psi_bench*(L_Inc_aux)**(1.0_DP-tauPL_bench)
		    	endif 

	        	K_Tax_draft_group_z(age,zi) = K_Tax_draft_group_z(age,zi) + & 
	        		& ( tauK_bench*( K_Inc_aux ) )*DBN_bench(age2,ai,zi,lambdai,ei,xi)

	        	if (age.lt.draft_age_category) then
        		L_Tax_draft_group_z(age,zi) = L_Tax_draft_group_z(age,zi) + & 
	        		& ( L_Inc_aux - psi_bench*(L_Inc_aux)**(1.0_DP-tauPL_bench) )* DBN_bench(age2,ai,zi,lambdai,ei,xi)
	        	endif 

                C_Tax_draft_group_z(age,zi) = C_Tax_draft_group_z(age,zi) + & 
	        		& ( tauC*Cons_bench(age2,ai,zi,lambdai,ei,xi) )*DBN_bench(age2,ai,zi,lambdai,ei,xi)

	        	Tot_Income_draft_group_z(age,zi) = Tot_Income_draft_group_z(age,zi) + & 
	        		& ( K_Inc_aux + L_Inc_aux )*DBN_bench(age2,ai,zi,lambdai,ei,xi)

        		K_Inc_draft_group_z(age,zi) = K_Inc_draft_group_z(age,zi) + ( K_Inc_aux )*DBN_bench(age2,ai,zi,lambdai,ei,xi)

        		L_Inc_draft_group_z(age,zi) = L_Inc_draft_group_z(age,zi) + ( L_Inc_aux )*DBN_bench(age2,ai,zi,lambdai,ei,xi)

        		K_Inc_frac_draft_group_z(age,zi) = K_Inc_frac_draft_group_z(age,zi) + & 
	        		& ( K_Inc_aux )/( K_Inc_aux + L_Inc_aux )*DBN_bench(age2,ai,zi,lambdai,ei,xi)

        		L_Inc_frac_draft_group_z(age,zi) = L_Inc_frac_draft_group_z(age,zi) + & 
	        		& ( L_Inc_aux )/( K_Inc_aux + L_Inc_aux )*DBN_bench(age2,ai,zi,lambdai,ei,xi)

        		K_Tax_Inc_draft_group_z(age,zi) = K_Tax_Inc_draft_group_z(age,zi) + & 
	        		& ( tauK_bench*( K_Inc_aux ) )*DBN_bench(age2,ai,zi,lambdai,ei,xi)/( K_Inc_aux + L_Inc_aux )

	        	if (age.lt.draft_age_category) then
        		L_Tax_Inc_draft_group_z(age,zi) = L_Tax_Inc_draft_group_z(age,zi) + & 
	        		& ( L_Inc_aux - psi_bench*(L_Inc_aux)**(1.0_DP-tauPL_bench) )*DBN_bench(age2,ai,zi,lambdai,ei,xi)/&
	        		& ( K_Inc_aux + L_Inc_aux )
	        	endif 

                C_Tax_Inc_draft_group_z(age,zi) = C_Tax_Inc_draft_group_z(age,zi) + & 
	        		& ( tauC*Cons_bench(age2,ai,zi,lambdai,ei,xi) )*DBN_bench(age2,ai,zi,lambdai,ei,xi)/( K_Inc_aux + L_Inc_aux )

        		wealth_draft_group_z(age,zi) = wealth_draft_group_z(age,zi) + agrid(ai)*DBN_bench(age2,ai,zi,lambdai,ei,xi)
	    	enddo 
	    	enddo 
	    	enddo 
	    	enddo 
	    	enddo  
		enddo
		enddo

		! Total Capital Tax adjusted by productivity group
		K_Tax_draft_group(:,1) = K_Tax_draft_group_z(:,1) + K_Tax_draft_group_z(:,2) + K_Tax_draft_group_z(:,3) & 
								&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*K_Tax_draft_group_z(:,4)
		K_Tax_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*K_Tax_draft_group_z(:,4) +& 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*K_Tax_draft_group_z(:,5)
		K_Tax_draft_group(:,3) = (0.10_dp/DBN_Z(5))*K_Tax_draft_group_z(:,5)
		K_Tax_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*K_Tax_draft_group_z(:,5) + & 
								&  ((0.99_dp-CDF_Z(5))/DBN_Z(6))*K_Tax_draft_group_z(:,6)
		K_Tax_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*K_Tax_draft_group_z(:,6) + & 
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*K_Tax_draft_group_z(:,7) 
		K_Tax_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*K_Tax_draft_group_z(:,7) +&
								&  K_Tax_draft_group_z(:,8) + K_Tax_draft_group_z(:,9) 
		K_Tax_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*K_Tax_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*K_Tax_draft_group_z(:,8)
		K_Tax_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*K_Tax_draft_group_z(:,8) + K_Tax_draft_group_z(:,9)

		! Total Labor Tax adjusted by productivity group
		L_Tax_draft_group(:,1) = L_Tax_draft_group_z(:,1) + L_Tax_draft_group_z(:,2) + L_Tax_draft_group_z(:,3) & 
								&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*L_Tax_draft_group_z(:,4)
		L_Tax_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*L_Tax_draft_group_z(:,4) +& 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*L_Tax_draft_group_z(:,5)
		L_Tax_draft_group(:,3) = (0.10_dp/DBN_Z(5))*L_Tax_draft_group_z(:,5)
		L_Tax_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*L_Tax_draft_group_z(:,5) + & 
								&  ((0.99_dp-CDF_Z(5))/DBN_Z(6))*L_Tax_draft_group_z(:,6)
		L_Tax_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*L_Tax_draft_group_z(:,6) + & 
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*L_Tax_draft_group_z(:,7) 
		L_Tax_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*L_Tax_draft_group_z(:,7) +&
								&  L_Tax_draft_group_z(:,8) + L_Tax_draft_group_z(:,9) 
		L_Tax_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*L_Tax_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*L_Tax_draft_group_z(:,8)
		L_Tax_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*L_Tax_draft_group_z(:,8) + L_Tax_draft_group_z(:,9)


		! Total Capital Tax adjusted by productivity group
		C_Tax_draft_group(:,1) = C_Tax_draft_group_z(:,1) + C_Tax_draft_group_z(:,2) + C_Tax_draft_group_z(:,3) & 
								&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*C_Tax_draft_group_z(:,4)
		C_Tax_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*C_Tax_draft_group_z(:,4) +& 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*C_Tax_draft_group_z(:,5)
		C_Tax_draft_group(:,3) = (0.10_dp/DBN_Z(5))*C_Tax_draft_group_z(:,5)
		C_Tax_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*C_Tax_draft_group_z(:,5) + & 
								&  ((0.99_dp-CDF_Z(5))/DBN_Z(6))*C_Tax_draft_group_z(:,6)
		C_Tax_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*C_Tax_draft_group_z(:,6) + & 
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*C_Tax_draft_group_z(:,7) 
		C_Tax_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*C_Tax_draft_group_z(:,7) +&
								&  C_Tax_draft_group_z(:,8) + C_Tax_draft_group_z(:,9) 
		C_Tax_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*C_Tax_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*C_Tax_draft_group_z(:,8)
		C_Tax_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*C_Tax_draft_group_z(:,8) + C_Tax_draft_group_z(:,9)

		! Total Pre-Tax Income adjusted by productivity group
		Tot_Income_draft_group(:,1) = Tot_Income_draft_group_z(:,1) + Tot_Income_draft_group_z(:,2) + Tot_Income_draft_group_z(:,3) & 
								&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*Tot_Income_draft_group_z(:,4)
		Tot_Income_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*Tot_Income_draft_group_z(:,4) +& 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*Tot_Income_draft_group_z(:,5)
		Tot_Income_draft_group(:,3) = (0.10_dp/DBN_Z(5))*Tot_Income_draft_group_z(:,5)
		Tot_Income_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*Tot_Income_draft_group_z(:,5) + & 
								&  ((0.99_dp-CDF_Z(5))/DBN_Z(6))*Tot_Income_draft_group_z(:,6)
		Tot_Income_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*Tot_Income_draft_group_z(:,6) + & 
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*Tot_Income_draft_group_z(:,7) 
		Tot_Income_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Tot_Income_draft_group_z(:,7) +&
								&  Tot_Income_draft_group_z(:,8) + Tot_Income_draft_group_z(:,9) 
		Tot_Income_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Tot_Income_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*Tot_Income_draft_group_z(:,8)
		Tot_Income_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*Tot_Income_draft_group_z(:,8) + Tot_Income_draft_group_z(:,9)

		! Total Pre-Tax Capital Income adjusted by productivity group
		K_Inc_draft_group(:,1) = K_Inc_draft_group_z(:,1) + K_Inc_draft_group_z(:,2) + K_Inc_draft_group_z(:,3) & 
								&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*K_Inc_draft_group_z(:,4)
		K_Inc_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*K_Inc_draft_group_z(:,4) +& 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*K_Inc_draft_group_z(:,5)
		K_Inc_draft_group(:,3) = (0.10_dp/DBN_Z(5))*K_Inc_draft_group_z(:,5)
		K_Inc_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*K_Inc_draft_group_z(:,5) + & 
								&  ((0.99_dp-CDF_Z(5))/DBN_Z(6))*K_Inc_draft_group_z(:,6)
		K_Inc_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*K_Inc_draft_group_z(:,6) + & 
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*K_Inc_draft_group_z(:,7) 
		K_Inc_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*K_Inc_draft_group_z(:,7) +&
								&  K_Inc_draft_group_z(:,8) + K_Inc_draft_group_z(:,9) 
		K_Inc_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*K_Inc_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*K_Inc_draft_group_z(:,8)
		K_Inc_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*K_Inc_draft_group_z(:,8) + K_Inc_draft_group_z(:,9)

		! Total Pre-Tax Capital Income adjusted by productivity group
		L_Inc_draft_group(:,1) = L_Inc_draft_group_z(:,1) + L_Inc_draft_group_z(:,2) + L_Inc_draft_group_z(:,3) & 
								&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*L_Inc_draft_group_z(:,4)
		L_Inc_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*L_Inc_draft_group_z(:,4) +& 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*L_Inc_draft_group_z(:,5)
		L_Inc_draft_group(:,3) = (0.10_dp/DBN_Z(5))*L_Inc_draft_group_z(:,5)
		L_Inc_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*L_Inc_draft_group_z(:,5) + & 
								&  ((0.99_dp-CDF_Z(5))/DBN_Z(6))*L_Inc_draft_group_z(:,6)
		L_Inc_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*L_Inc_draft_group_z(:,6) + & 
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*L_Inc_draft_group_z(:,7) 
		L_Inc_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*L_Inc_draft_group_z(:,7) +&
								&  L_Inc_draft_group_z(:,8) + L_Inc_draft_group_z(:,9) 
		L_Inc_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*L_Inc_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*L_Inc_draft_group_z(:,8)
		L_Inc_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*L_Inc_draft_group_z(:,8) + L_Inc_draft_group_z(:,9)


		! Get ratios to total income in group
		K_Tax_draft_group = K_Tax_draft_group/Tot_Income_draft_group
		L_Tax_draft_group = L_Tax_draft_group/Tot_Income_draft_group
		C_Tax_draft_group = C_Tax_draft_group/Tot_Income_draft_group
		K_Inc_draft_group = K_Inc_draft_group/Tot_Income_draft_group
		L_Inc_draft_group = L_Inc_draft_group/Tot_Income_draft_group

		! Divide by mass in group
		K_Tax_Inc_draft_group_z  = K_Tax_Inc_draft_group_z/size_draft_group_z
		L_Tax_Inc_draft_group_z  = L_Tax_Inc_draft_group_z/size_draft_group_z
		C_Tax_Inc_draft_group_z  = C_Tax_Inc_draft_group_z/size_draft_group_z
		K_Inc_frac_draft_group_z = K_Inc_frac_draft_group_z/size_draft_group_z 
		L_Inc_frac_draft_group_z = L_Inc_frac_draft_group_z/size_draft_group_z 

		! Total Capital Tax to income ratio adjusted by productivity group
		K_Tax_Inc_draft_group(:,1)   = ( DBN_Z(1)*K_Tax_Inc_draft_group_z(:,1) + DBN_Z(2)*K_Tax_Inc_draft_group_z(:,2) + & 
								& DBN_Z(3)*K_Tax_Inc_draft_group_z(:,3) + (0.40_dp-CDF_Z(3))*K_Tax_Inc_draft_group_z(:,4) )/0.40_dp
		K_Tax_Inc_draft_group(:,2)   = ( (CDF_Z(4)-0.40_dp)*K_Tax_Inc_draft_group_z(:,4) + &
								& (0.80_dp-CDF_Z(4))*K_Tax_Inc_draft_group_z(:,5) )/0.40_dp
		K_Tax_Inc_draft_group(:,3)   = K_Tax_Inc_draft_group_z(:,5)
		K_Tax_Inc_draft_group(:,4)   = ( (CDF_Z(5)-0.90_dp)*K_Tax_Inc_draft_group_z(:,5) + &
								& (0.99_dp-CDF_Z(5))*K_Tax_Inc_draft_group_z(:,6) )/0.09_dp
		K_Tax_Inc_draft_group(:,5)   = ( (CDF_Z(6)-0.99_dp)*K_Tax_Inc_draft_group_z(:,6) + &
								& (0.999_dp-CDF_Z(6))*K_Tax_Inc_draft_group_z(:,7) )/0.009_dp
		K_Tax_Inc_draft_group(:,6)   = ( (CDF_Z(7)-0.999_dp)*K_Tax_Inc_draft_group_z(:,7) + &
								&  DBN_Z(8)*K_Tax_Inc_draft_group_z(:,8) + DBN_Z(9)*K_Tax_Inc_draft_group_z(:,9) )/0.001_dp
		K_Tax_Inc_draft_group(:,7)   = ( (CDF_Z(7)-0.999_dp)*K_Tax_Inc_draft_group_z(:,7) + &
								& (0.9999_dp-CDF_Z(7))*K_Tax_Inc_draft_group_z(:,8) )/0.0009_dp
		K_Tax_Inc_draft_group(:,8)   = ( (CDF_Z(8)-0.9999_dp)*K_Tax_Inc_draft_group_z(:,8) + &
								& DBN_Z(9)*K_Tax_Inc_draft_group_z(:,9) )/0.0001_dp

		! Total Labor Tax to income ratio adjusted by productivity group
		L_Tax_Inc_draft_group(:,1)   = ( DBN_Z(1)*L_Tax_Inc_draft_group_z(:,1) + DBN_Z(2)*L_Tax_Inc_draft_group_z(:,2) + & 
								& DBN_Z(3)*L_Tax_Inc_draft_group_z(:,3) + (0.40_dp-CDF_Z(3))*L_Tax_Inc_draft_group_z(:,4) )/0.40_dp
		L_Tax_Inc_draft_group(:,2)   = ( (CDF_Z(4)-0.40_dp)*L_Tax_Inc_draft_group_z(:,4) + &
								& (0.80_dp-CDF_Z(4))*L_Tax_Inc_draft_group_z(:,5) )/0.40_dp
		L_Tax_Inc_draft_group(:,3)   = L_Tax_Inc_draft_group_z(:,5)
		L_Tax_Inc_draft_group(:,4)   = ( (CDF_Z(5)-0.90_dp)*L_Tax_Inc_draft_group_z(:,5) + &
								& (0.99_dp-CDF_Z(5))*L_Tax_Inc_draft_group_z(:,6) )/0.09_dp
		L_Tax_Inc_draft_group(:,5)   = ( (CDF_Z(6)-0.99_dp)*L_Tax_Inc_draft_group_z(:,6) + &
								& (0.999_dp-CDF_Z(6))*L_Tax_Inc_draft_group_z(:,7) )/0.009_dp
		L_Tax_Inc_draft_group(:,6)   = ( (CDF_Z(7)-0.999_dp)*L_Tax_Inc_draft_group_z(:,7) + &
								&  DBN_Z(8)*L_Tax_Inc_draft_group_z(:,8) + DBN_Z(9)*L_Tax_Inc_draft_group_z(:,9) )/0.001_dp
		L_Tax_Inc_draft_group(:,7)   = ( (CDF_Z(7)-0.999_dp)*L_Tax_Inc_draft_group_z(:,7) + &
								& (0.9999_dp-CDF_Z(7))*L_Tax_Inc_draft_group_z(:,8) )/0.0009_dp
		L_Tax_Inc_draft_group(:,8)   = ( (CDF_Z(8)-0.9999_dp)*L_Tax_Inc_draft_group_z(:,8) + &
								& DBN_Z(9)*L_Tax_Inc_draft_group_z(:,9) )/0.0001_dp

		! Total Consumption Tax to income ratio adjusted by productivity group
		C_Tax_Inc_draft_group(:,1)   = ( DBN_Z(1)*C_Tax_Inc_draft_group_z(:,1) + DBN_Z(2)*C_Tax_Inc_draft_group_z(:,2) + & 
								& DBN_Z(3)*C_Tax_Inc_draft_group_z(:,3) + (0.40_dp-CDF_Z(3))*C_Tax_Inc_draft_group_z(:,4) )/0.40_dp
		C_Tax_Inc_draft_group(:,2)   = ( (CDF_Z(4)-0.40_dp)*C_Tax_Inc_draft_group_z(:,4) + &
								& (0.80_dp-CDF_Z(4))*C_Tax_Inc_draft_group_z(:,5) )/0.40_dp
		C_Tax_Inc_draft_group(:,3)   = C_Tax_Inc_draft_group_z(:,5)
		C_Tax_Inc_draft_group(:,4)   = ( (CDF_Z(5)-0.90_dp)*C_Tax_Inc_draft_group_z(:,5) + &
								& (0.99_dp-CDF_Z(5))*C_Tax_Inc_draft_group_z(:,6) )/0.09_dp
		C_Tax_Inc_draft_group(:,5)   = ( (CDF_Z(6)-0.99_dp)*C_Tax_Inc_draft_group_z(:,6) + &
								& (0.999_dp-CDF_Z(6))*C_Tax_Inc_draft_group_z(:,7) )/0.009_dp
		C_Tax_Inc_draft_group(:,6)   = ( (CDF_Z(7)-0.999_dp)*C_Tax_Inc_draft_group_z(:,7) + &
								&  DBN_Z(8)*C_Tax_Inc_draft_group_z(:,8) + DBN_Z(9)*C_Tax_Inc_draft_group_z(:,9) )/0.001_dp
		C_Tax_Inc_draft_group(:,7)   = ( (CDF_Z(7)-0.999_dp)*C_Tax_Inc_draft_group_z(:,7) + &
								& (0.9999_dp-CDF_Z(7))*C_Tax_Inc_draft_group_z(:,8) )/0.0009_dp
		C_Tax_Inc_draft_group(:,8)   = ( (CDF_Z(8)-0.9999_dp)*C_Tax_Inc_draft_group_z(:,8) + &
								& DBN_Z(9)*C_Tax_Inc_draft_group_z(:,9) )/0.0001_dp


		! Capital Income Share Tax adjusted by productivity group
		K_Inc_frac_draft_group(:,1)   = ( DBN_Z(1)*K_Inc_frac_draft_group_z(:,1) + DBN_Z(2)*K_Inc_frac_draft_group_z(:,2) + & 
								& DBN_Z(3)*K_Inc_frac_draft_group_z(:,3) + (0.40_dp-CDF_Z(3))*K_Inc_frac_draft_group_z(:,4) )/0.40_dp
		K_Inc_frac_draft_group(:,2)   = ( (CDF_Z(4)-0.40_dp)*K_Inc_frac_draft_group_z(:,4) + &
								& (0.80_dp-CDF_Z(4))*K_Inc_frac_draft_group_z(:,5) )/0.40_dp
		K_Inc_frac_draft_group(:,3)   = K_Inc_frac_draft_group_z(:,5)
		K_Inc_frac_draft_group(:,4)   = ( (CDF_Z(5)-0.90_dp)*K_Inc_frac_draft_group_z(:,5) + &
								& (0.99_dp-CDF_Z(5))*K_Inc_frac_draft_group_z(:,6) )/0.09_dp
		K_Inc_frac_draft_group(:,5)   = ( (CDF_Z(6)-0.99_dp)*K_Inc_frac_draft_group_z(:,6) + &
								& (0.999_dp-CDF_Z(6))*K_Inc_frac_draft_group_z(:,7) )/0.009_dp
		K_Inc_frac_draft_group(:,6)   = ( (CDF_Z(7)-0.999_dp)*K_Inc_frac_draft_group_z(:,7) + &
								&  DBN_Z(8)*K_Inc_frac_draft_group_z(:,8) + DBN_Z(9)*K_Inc_frac_draft_group_z(:,9) )/0.001_dp
		K_Inc_frac_draft_group(:,7)   = ( (CDF_Z(7)-0.999_dp)*K_Inc_frac_draft_group_z(:,7) + &
								& (0.9999_dp-CDF_Z(7))*K_Inc_frac_draft_group_z(:,8) )/0.0009_dp
		K_Inc_frac_draft_group(:,8)   = ( (CDF_Z(8)-0.9999_dp)*K_Inc_frac_draft_group_z(:,8) + &
								& DBN_Z(9)*K_Inc_frac_draft_group_z(:,9) )/0.0001_dp

		! Labor Income Share Tax adjusted by productivity group
		L_Inc_frac_draft_group(:,1)   = ( DBN_Z(1)*L_Inc_frac_draft_group_z(:,1) + DBN_Z(2)*L_Inc_frac_draft_group_z(:,2) + & 
								& DBN_Z(3)*L_Inc_frac_draft_group_z(:,3) + (0.40_dp-CDF_Z(3))*L_Inc_frac_draft_group_z(:,4) )/0.40_dp
		L_Inc_frac_draft_group(:,2)   = ( (CDF_Z(4)-0.40_dp)*L_Inc_frac_draft_group_z(:,4) + &
								& (0.80_dp-CDF_Z(4))*L_Inc_frac_draft_group_z(:,5) )/0.40_dp
		L_Inc_frac_draft_group(:,3)   = L_Inc_frac_draft_group_z(:,5)
		L_Inc_frac_draft_group(:,4)   = ( (CDF_Z(5)-0.90_dp)*L_Inc_frac_draft_group_z(:,5) + &
								& (0.99_dp-CDF_Z(5))*L_Inc_frac_draft_group_z(:,6) )/0.09_dp
		L_Inc_frac_draft_group(:,5)   = ( (CDF_Z(6)-0.99_dp)*L_Inc_frac_draft_group_z(:,6) + &
								& (0.999_dp-CDF_Z(6))*L_Inc_frac_draft_group_z(:,7) )/0.009_dp
		L_Inc_frac_draft_group(:,6)   = ( (CDF_Z(7)-0.999_dp)*L_Inc_frac_draft_group_z(:,7) + &
								&  DBN_Z(8)*L_Inc_frac_draft_group_z(:,8) + DBN_Z(9)*L_Inc_frac_draft_group_z(:,9) )/0.001_dp
		L_Inc_frac_draft_group(:,7)   = ( (CDF_Z(7)-0.999_dp)*L_Inc_frac_draft_group_z(:,7) + &
								& (0.9999_dp-CDF_Z(7))*L_Inc_frac_draft_group_z(:,8) )/0.0009_dp
		L_Inc_frac_draft_group(:,8)   = ( (CDF_Z(8)-0.9999_dp)*L_Inc_frac_draft_group_z(:,8) + &
								& DBN_Z(9)*L_Inc_frac_draft_group_z(:,9) )/0.0001_dp

		! Size of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		wealth_draft_group(:,1) = wealth_draft_group_z(:,1) + wealth_draft_group_z(:,2) + wealth_draft_group_z(:,3) & 
								&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*wealth_draft_group_z(:,4)
		wealth_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*wealth_draft_group_z(:,4) +& 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*wealth_draft_group_z(:,5)
		wealth_draft_group(:,3) = (0.10_dp/DBN_Z(5))*wealth_draft_group_z(:,5)
		wealth_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*wealth_draft_group_z(:,5) + & 
								&  ((0.99_dp-CDF_Z(5))/DBN_Z(6))*wealth_draft_group_z(:,6)
		wealth_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*wealth_draft_group_z(:,6) + & 
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*wealth_draft_group_z(:,7) 
		wealth_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*wealth_draft_group_z(:,7) +&
								&  wealth_draft_group_z(:,8) + wealth_draft_group_z(:,9) 
		wealth_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*wealth_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*wealth_draft_group_z(:,8)
		wealth_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*wealth_draft_group_z(:,8) + wealth_draft_group_z(:,9)


		OPEN (UNIT=80, FILE=trim(Result_Folder)//'draft_group_Tax_K_bench.txt', STATUS='replace') 
	    OPEN (UNIT=81, FILE=trim(Result_Folder)//'draft_group_Tax_L_bench.txt', STATUS='replace') 
	    OPEN (UNIT=82, FILE=trim(Result_Folder)//'draft_group_Tax_C_bench.txt', STATUS='replace') 
	    OPEN (UNIT=83, FILE=trim(Result_Folder)//'draft_group_Inc_bench.txt', STATUS='replace') 
	    OPEN (UNIT=84, FILE=trim(Result_Folder)//'draft_group_Tax_K_Inc_bench.txt', STATUS='replace') 
	    OPEN (UNIT=85, FILE=trim(Result_Folder)//'draft_group_Tax_L_Inc_bench.txt', STATUS='replace') 
	    OPEN (UNIT=86, FILE=trim(Result_Folder)//'draft_group_Tax_C_Inc_bench.txt', STATUS='replace') 
	    OPEN (UNIT=87, FILE=trim(Result_Folder)//'draft_group_K_Inc_bench.txt', STATUS='replace') 
	    OPEN (UNIT=88, FILE=trim(Result_Folder)//'draft_group_L_Inc_bench.txt', STATUS='replace') 
	    OPEN (UNIT=89, FILE=trim(Result_Folder)//'draft_group_K_Inc_frac_bench.txt', STATUS='replace') 
	    OPEN (UNIT=90, FILE=trim(Result_Folder)//'draft_group_L_Inc_frac_bench.txt', STATUS='replace') 
	    OPEN (UNIT=91, FILE=trim(Result_Folder)//'draft_group_Wealth_bench.txt', STATUS='replace') 
		do age = 1,draft_age_category
		    WRITE  (UNIT=80, FMT=*)  K_Tax_draft_group(age,:)
		    WRITE  (UNIT=81, FMT=*)  L_Tax_draft_group(age,:)
		    WRITE  (UNIT=82, FMT=*)  C_Tax_draft_group(age,:)
		    WRITE  (UNIT=83, FMT=*)  Tot_Income_draft_group(age,:)
		    WRITE  (UNIT=84, FMT=*)  K_Tax_Inc_draft_group(age,:)
		    WRITE  (UNIT=85, FMT=*)  L_Tax_Inc_draft_group(age,:)
		    WRITE  (UNIT=86, FMT=*)  C_Tax_Inc_draft_group(age,:)
		    WRITE  (UNIT=87, FMT=*)  K_Inc_draft_group(age,:)
		    WRITE  (UNIT=88, FMT=*)  L_Inc_draft_group(age,:)
		    WRITE  (UNIT=89, FMT=*)  K_Inc_frac_draft_group(age,:)
		    WRITE  (UNIT=90, FMT=*)  L_Inc_frac_draft_group(age,:)
		    WRITE  (UNIT=91, FMT=*)  wealth_draft_group(age,:)
		ENDDO
		close(unit=80); close(unit=81); close(unit=82); close(unit=83); close(unit=84); close(unit=85)
		close(unit=86); close(unit=87); close(unit=88); close(unit=89); close(unit=90); close(unit=91); 

		
		!! Experiment 
		Pr_mat = Profit_Matrix(R_exp,P_exp)
		CALL ComputeLaborUnits(EBAR_exp,wage_exp)
		CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)

		K_Tax_draft_group_z			= 0.0_dp
		L_Tax_draft_group_z			= 0.0_dp
		C_Tax_draft_group_z			= 0.0_dp
		Tot_Income_draft_group_z	= 0.0_dp
		K_Inc_draft_group_z         = 0.0_dp
		L_Inc_draft_group_z         = 0.0_dp
		K_Inc_frac_draft_group_z    = 0.0_dp
		L_Inc_frac_draft_group_z    = 0.0_dp
		K_Tax_Inc_draft_group_z		= 0.0_dp
		L_Tax_Inc_draft_group_z		= 0.0_dp
		C_Tax_Inc_draft_group_z		= 0.0_dp
		wealth_draft_group_z 	    = 0.0_dp

		do zi  = 1,nz
		do age = 1,draft_age_category
	        do xi=1,nx
		    do ei=1,ne
		    do lambdai=1,nlambda
		    do ai=1,na
		    do age2=draft_age_limit(age)+1,draft_age_limit(age+1)

		    	if (age.lt.draft_age_category) then
		    	L_Inc_aux   = yh(age2,lambdai,ei)*Hours_exp(age2,ai,zi,lambdai,ei,xi)
		    	else
		    	L_Inc_aux   = RetY_lambda_e(lambdai,ei) 
		    	endif 
		    	K_Inc_aux   = R_exp*agrid(ai) + Pr_mat(ai,zi,xi)

		    	Income_exp(age2,ai,zi,lambdai,ei,xi) = L_Inc_aux + K_Inc_aux

		    	K_Tax_exp(age2,ai,zi,lambdai,ei,xi)  = ((1.0_dp+R_exp)*agrid(ai) + Pr_mat(ai,zi,xi))-YGRID(ai,zi,xi) 
		    	if (age.lt.draft_age_category) then
		    	L_Tax_exp(age2,ai,zi,lambdai,ei,xi)  =  L_Inc_aux - psi_exp*(L_Inc_aux)**(1.0_DP-tauPL_exp) 
		    	endif 


	        	K_Tax_draft_group_z(age,zi) = K_Tax_draft_group_z(age,zi) + & 
	        		& (((1.0_dp+R_exp)*agrid(ai) + Pr_mat(ai,zi,xi))-YGRID(ai,zi,xi) )*DBN_exp(age2,ai,zi,lambdai,ei,xi)

	        	if (age.lt.draft_age_category) then
        		L_Tax_draft_group_z(age,zi) = L_Tax_draft_group_z(age,zi) + & 
	        		& ( L_Inc_aux - psi_exp*(L_Inc_aux)**(1.0_DP-tauPL_exp) )* DBN_exp(age2,ai,zi,lambdai,ei,xi)
	        	endif 

                C_Tax_draft_group_z(age,zi) = C_Tax_draft_group_z(age,zi) + & 
	        		& ( tauC*Cons_exp(age2,ai,zi,lambdai,ei,xi) )*DBN_exp(age2,ai,zi,lambdai,ei,xi)

	        	Tot_Income_draft_group_z(age,zi) = Tot_Income_draft_group_z(age,zi) + & 
	        		& ( K_Inc_aux + L_Inc_aux )*DBN_exp(age2,ai,zi,lambdai,ei,xi)

        		K_Inc_draft_group_z(age,zi) = K_Inc_draft_group_z(age,zi) + ( K_Inc_aux )*DBN_exp(age2,ai,zi,lambdai,ei,xi)

        		L_Inc_draft_group_z(age,zi) = L_Inc_draft_group_z(age,zi) + ( L_Inc_aux )*DBN_exp(age2,ai,zi,lambdai,ei,xi)

        		K_Inc_frac_draft_group_z(age,zi) = K_Inc_frac_draft_group_z(age,zi) + & 
	        		& ( K_Inc_aux )/( K_Inc_aux + L_Inc_aux )*DBN_exp(age2,ai,zi,lambdai,ei,xi)

        		L_Inc_frac_draft_group_z(age,zi) = L_Inc_frac_draft_group_z(age,zi) + & 
	        		& ( L_Inc_aux )/( K_Inc_aux + L_Inc_aux )*DBN_exp(age2,ai,zi,lambdai,ei,xi)

        		K_Tax_Inc_draft_group_z(age,zi) = K_Tax_Inc_draft_group_z(age,zi) + & 
	        		& (((1.0_dp+R_exp)*agrid(ai) + Pr_mat(ai,zi,xi))-YGRID(ai,zi,xi))*DBN_exp(age2,ai,zi,lambdai,ei,xi)/&
	        		& ( K_Inc_aux + L_Inc_aux )
	        		
	        	if (age.lt.draft_age_category) then
        		L_Tax_Inc_draft_group_z(age,zi) = L_Tax_Inc_draft_group_z(age,zi) + & 
	        		& ( L_Inc_aux - psi_exp*(L_Inc_aux)**(1.0_DP-tauPL_exp) )*DBN_exp(age2,ai,zi,lambdai,ei,xi)/&
	        		& ( K_Inc_aux + L_Inc_aux )
	        	endif 

                C_Tax_Inc_draft_group_z(age,zi) = C_Tax_Inc_draft_group_z(age,zi) + & 
	        		& ( tauC*Cons_exp(age2,ai,zi,lambdai,ei,xi) )*DBN_exp(age2,ai,zi,lambdai,ei,xi)/( K_Inc_aux + L_Inc_aux )

        		wealth_draft_group_z(age,zi) = wealth_draft_group_z(age,zi) + agrid(ai)*DBN_exp(age2,ai,zi,lambdai,ei,xi)

	    	enddo 
	    	enddo 
	    	enddo 
	    	enddo 
	    	enddo  
		enddo
		enddo

		! Total Capital Tax adjusted by productivity group
		K_Tax_draft_group(:,1) = K_Tax_draft_group_z(:,1) + K_Tax_draft_group_z(:,2) + K_Tax_draft_group_z(:,3) & 
								&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*K_Tax_draft_group_z(:,4)
		K_Tax_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*K_Tax_draft_group_z(:,4) +& 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*K_Tax_draft_group_z(:,5)
		K_Tax_draft_group(:,3) = (0.10_dp/DBN_Z(5))*K_Tax_draft_group_z(:,5)
		K_Tax_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*K_Tax_draft_group_z(:,5) + & 
								&  ((0.99_dp-CDF_Z(5))/DBN_Z(6))*K_Tax_draft_group_z(:,6)
		K_Tax_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*K_Tax_draft_group_z(:,6) + & 
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*K_Tax_draft_group_z(:,7) 
		K_Tax_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*K_Tax_draft_group_z(:,7) +&
								&  K_Tax_draft_group_z(:,8) + K_Tax_draft_group_z(:,9) 
		K_Tax_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*K_Tax_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*K_Tax_draft_group_z(:,8)
		K_Tax_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*K_Tax_draft_group_z(:,8) + K_Tax_draft_group_z(:,9)

		! Total Labor Tax adjusted by productivity group
		L_Tax_draft_group(:,1) = L_Tax_draft_group_z(:,1) + L_Tax_draft_group_z(:,2) + L_Tax_draft_group_z(:,3) & 
								&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*L_Tax_draft_group_z(:,4)
		L_Tax_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*L_Tax_draft_group_z(:,4) +& 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*L_Tax_draft_group_z(:,5)
		L_Tax_draft_group(:,3) = (0.10_dp/DBN_Z(5))*L_Tax_draft_group_z(:,5)
		L_Tax_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*L_Tax_draft_group_z(:,5) + & 
								&  ((0.99_dp-CDF_Z(5))/DBN_Z(6))*L_Tax_draft_group_z(:,6)
		L_Tax_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*L_Tax_draft_group_z(:,6) + & 
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*L_Tax_draft_group_z(:,7) 
		L_Tax_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*L_Tax_draft_group_z(:,7) +&
								&  L_Tax_draft_group_z(:,8) + L_Tax_draft_group_z(:,9) 
		L_Tax_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*L_Tax_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*L_Tax_draft_group_z(:,8)
		L_Tax_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*L_Tax_draft_group_z(:,8) + L_Tax_draft_group_z(:,9)


		! Total Capital Tax adjusted by productivity group
		C_Tax_draft_group(:,1) = C_Tax_draft_group_z(:,1) + C_Tax_draft_group_z(:,2) + C_Tax_draft_group_z(:,3) & 
								&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*C_Tax_draft_group_z(:,4)
		C_Tax_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*C_Tax_draft_group_z(:,4) +& 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*C_Tax_draft_group_z(:,5)
		C_Tax_draft_group(:,3) = (0.10_dp/DBN_Z(5))*C_Tax_draft_group_z(:,5)
		C_Tax_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*C_Tax_draft_group_z(:,5) + & 
								&  ((0.99_dp-CDF_Z(5))/DBN_Z(6))*C_Tax_draft_group_z(:,6)
		C_Tax_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*C_Tax_draft_group_z(:,6) + & 
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*C_Tax_draft_group_z(:,7) 
		C_Tax_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*C_Tax_draft_group_z(:,7) +&
								&  C_Tax_draft_group_z(:,8) + C_Tax_draft_group_z(:,9) 
		C_Tax_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*C_Tax_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*C_Tax_draft_group_z(:,8)
		C_Tax_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*C_Tax_draft_group_z(:,8) + C_Tax_draft_group_z(:,9)

		! Total Pre-Tax Income adjusted by productivity group
		Tot_Income_draft_group(:,1) = Tot_Income_draft_group_z(:,1) + Tot_Income_draft_group_z(:,2) + Tot_Income_draft_group_z(:,3) & 
								&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*Tot_Income_draft_group_z(:,4)
		Tot_Income_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*Tot_Income_draft_group_z(:,4) +& 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*Tot_Income_draft_group_z(:,5)
		Tot_Income_draft_group(:,3) = (0.10_dp/DBN_Z(5))*Tot_Income_draft_group_z(:,5)
		Tot_Income_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*Tot_Income_draft_group_z(:,5) + & 
								&  ((0.99_dp-CDF_Z(5))/DBN_Z(6))*Tot_Income_draft_group_z(:,6)
		Tot_Income_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*Tot_Income_draft_group_z(:,6) + & 
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*Tot_Income_draft_group_z(:,7) 
		Tot_Income_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Tot_Income_draft_group_z(:,7) +&
								&  Tot_Income_draft_group_z(:,8) + Tot_Income_draft_group_z(:,9) 
		Tot_Income_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Tot_Income_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*Tot_Income_draft_group_z(:,8)
		Tot_Income_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*Tot_Income_draft_group_z(:,8) + Tot_Income_draft_group_z(:,9)

		! Total Pre-Tax Capital Income adjusted by productivity group
		K_Inc_draft_group(:,1) = K_Inc_draft_group_z(:,1) + K_Inc_draft_group_z(:,2) + K_Inc_draft_group_z(:,3) & 
								&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*K_Inc_draft_group_z(:,4)
		K_Inc_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*K_Inc_draft_group_z(:,4) +& 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*K_Inc_draft_group_z(:,5)
		K_Inc_draft_group(:,3) = (0.10_dp/DBN_Z(5))*K_Inc_draft_group_z(:,5)
		K_Inc_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*K_Inc_draft_group_z(:,5) + & 
								&  ((0.99_dp-CDF_Z(5))/DBN_Z(6))*K_Inc_draft_group_z(:,6)
		K_Inc_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*K_Inc_draft_group_z(:,6) + & 
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*K_Inc_draft_group_z(:,7) 
		K_Inc_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*K_Inc_draft_group_z(:,7) +&
								&  K_Inc_draft_group_z(:,8) + K_Inc_draft_group_z(:,9) 
		K_Inc_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*K_Inc_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*K_Inc_draft_group_z(:,8)
		K_Inc_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*K_Inc_draft_group_z(:,8) + K_Inc_draft_group_z(:,9)

		! Total Pre-Tax Capital Income adjusted by productivity group
		L_Inc_draft_group(:,1) = L_Inc_draft_group_z(:,1) + L_Inc_draft_group_z(:,2) + L_Inc_draft_group_z(:,3) & 
								&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*L_Inc_draft_group_z(:,4)
		L_Inc_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*L_Inc_draft_group_z(:,4) +& 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*L_Inc_draft_group_z(:,5)
		L_Inc_draft_group(:,3) = (0.10_dp/DBN_Z(5))*L_Inc_draft_group_z(:,5)
		L_Inc_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*L_Inc_draft_group_z(:,5) + & 
								&  ((0.99_dp-CDF_Z(5))/DBN_Z(6))*L_Inc_draft_group_z(:,6)
		L_Inc_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*L_Inc_draft_group_z(:,6) + & 
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*L_Inc_draft_group_z(:,7) 
		L_Inc_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*L_Inc_draft_group_z(:,7) +&
								&  L_Inc_draft_group_z(:,8) + L_Inc_draft_group_z(:,9) 
		L_Inc_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*L_Inc_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*L_Inc_draft_group_z(:,8)
		L_Inc_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*L_Inc_draft_group_z(:,8) + L_Inc_draft_group_z(:,9)


		! Get ratios to total income in group
		K_Tax_draft_group = K_Tax_draft_group/Tot_Income_draft_group
		L_Tax_draft_group = L_Tax_draft_group/Tot_Income_draft_group
		C_Tax_draft_group = C_Tax_draft_group/Tot_Income_draft_group
		K_Inc_draft_group = K_Inc_draft_group/Tot_Income_draft_group
		L_Inc_draft_group = L_Inc_draft_group/Tot_Income_draft_group


		! Divide by mass in group
		K_Tax_Inc_draft_group_z  = K_Tax_Inc_draft_group_z/size_draft_group_z
		L_Tax_Inc_draft_group_z  = L_Tax_Inc_draft_group_z/size_draft_group_z
		C_Tax_Inc_draft_group_z  = C_Tax_Inc_draft_group_z/size_draft_group_z
		K_Inc_frac_draft_group_z = K_Inc_frac_draft_group_z/size_draft_group_z 
		L_Inc_frac_draft_group_z = L_Inc_frac_draft_group_z/size_draft_group_z 

		! Total Capital Tax to income ratio adjusted by productivity group
		K_Tax_Inc_draft_group(:,1)   = ( DBN_Z(1)*K_Tax_Inc_draft_group_z(:,1) + DBN_Z(2)*K_Tax_Inc_draft_group_z(:,2) + & 
								& DBN_Z(3)*K_Tax_Inc_draft_group_z(:,3) + (0.40_dp-CDF_Z(3))*K_Tax_Inc_draft_group_z(:,4) )/0.40_dp
		K_Tax_Inc_draft_group(:,2)   = ( (CDF_Z(4)-0.40_dp)*K_Tax_Inc_draft_group_z(:,4) + &
								& (0.80_dp-CDF_Z(4))*K_Tax_Inc_draft_group_z(:,5) )/0.40_dp
		K_Tax_Inc_draft_group(:,3)   = K_Tax_Inc_draft_group_z(:,5)
		K_Tax_Inc_draft_group(:,4)   = ( (CDF_Z(5)-0.90_dp)*K_Tax_Inc_draft_group_z(:,5) + &
								& (0.99_dp-CDF_Z(5))*K_Tax_Inc_draft_group_z(:,6) )/0.09_dp
		K_Tax_Inc_draft_group(:,5)   = ( (CDF_Z(6)-0.99_dp)*K_Tax_Inc_draft_group_z(:,6) + &
								& (0.999_dp-CDF_Z(6))*K_Tax_Inc_draft_group_z(:,7) )/0.009_dp
		K_Tax_Inc_draft_group(:,6)   = ( (CDF_Z(7)-0.999_dp)*K_Tax_Inc_draft_group_z(:,7) + &
								&  DBN_Z(8)*K_Tax_Inc_draft_group_z(:,8) + DBN_Z(9)*K_Tax_Inc_draft_group_z(:,9) )/0.001_dp
		K_Tax_Inc_draft_group(:,7)   = ( (CDF_Z(7)-0.999_dp)*K_Tax_Inc_draft_group_z(:,7) + &
								& (0.9999_dp-CDF_Z(7))*K_Tax_Inc_draft_group_z(:,8) )/0.0009_dp
		K_Tax_Inc_draft_group(:,8)   = ( (CDF_Z(8)-0.9999_dp)*K_Tax_Inc_draft_group_z(:,8) + &
								& DBN_Z(9)*K_Tax_Inc_draft_group_z(:,9) )/0.0001_dp

		! Total Labor Tax to income ratio adjusted by productivity group
		L_Tax_Inc_draft_group(:,1)   = ( DBN_Z(1)*L_Tax_Inc_draft_group_z(:,1) + DBN_Z(2)*L_Tax_Inc_draft_group_z(:,2) + & 
								& DBN_Z(3)*L_Tax_Inc_draft_group_z(:,3) + (0.40_dp-CDF_Z(3))*L_Tax_Inc_draft_group_z(:,4) )/0.40_dp
		L_Tax_Inc_draft_group(:,2)   = ( (CDF_Z(4)-0.40_dp)*L_Tax_Inc_draft_group_z(:,4) + &
								& (0.80_dp-CDF_Z(4))*L_Tax_Inc_draft_group_z(:,5) )/0.40_dp
		L_Tax_Inc_draft_group(:,3)   = L_Tax_Inc_draft_group_z(:,5)
		L_Tax_Inc_draft_group(:,4)   = ( (CDF_Z(5)-0.90_dp)*L_Tax_Inc_draft_group_z(:,5) + &
								& (0.99_dp-CDF_Z(5))*L_Tax_Inc_draft_group_z(:,6) )/0.09_dp
		L_Tax_Inc_draft_group(:,5)   = ( (CDF_Z(6)-0.99_dp)*L_Tax_Inc_draft_group_z(:,6) + &
								& (0.999_dp-CDF_Z(6))*L_Tax_Inc_draft_group_z(:,7) )/0.009_dp
		L_Tax_Inc_draft_group(:,6)   = ( (CDF_Z(7)-0.999_dp)*L_Tax_Inc_draft_group_z(:,7) + &
								&  DBN_Z(8)*L_Tax_Inc_draft_group_z(:,8) + DBN_Z(9)*L_Tax_Inc_draft_group_z(:,9) )/0.001_dp
		L_Tax_Inc_draft_group(:,7)   = ( (CDF_Z(7)-0.999_dp)*L_Tax_Inc_draft_group_z(:,7) + &
								& (0.9999_dp-CDF_Z(7))*L_Tax_Inc_draft_group_z(:,8) )/0.0009_dp
		L_Tax_Inc_draft_group(:,8)   = ( (CDF_Z(8)-0.9999_dp)*L_Tax_Inc_draft_group_z(:,8) + &
								& DBN_Z(9)*L_Tax_Inc_draft_group_z(:,9) )/0.0001_dp


		! Total Capital Tax adjusted by productivity group
		C_Tax_Inc_draft_group(:,1)   = ( DBN_Z(1)*C_Tax_Inc_draft_group_z(:,1) + DBN_Z(2)*C_Tax_Inc_draft_group_z(:,2) + & 
								& DBN_Z(3)*C_Tax_Inc_draft_group_z(:,3) + (0.40_dp-CDF_Z(3))*C_Tax_Inc_draft_group_z(:,4) )/0.40_dp
		C_Tax_Inc_draft_group(:,2)   = ( (CDF_Z(4)-0.40_dp)*C_Tax_Inc_draft_group_z(:,4) + &
								& (0.80_dp-CDF_Z(4))*C_Tax_Inc_draft_group_z(:,5) )/0.40_dp
		C_Tax_Inc_draft_group(:,3)   = C_Tax_Inc_draft_group_z(:,5)
		C_Tax_Inc_draft_group(:,4)   = ( (CDF_Z(5)-0.90_dp)*C_Tax_Inc_draft_group_z(:,5) + &
								& (0.99_dp-CDF_Z(5))*C_Tax_Inc_draft_group_z(:,6) )/0.09_dp
		C_Tax_Inc_draft_group(:,5)   = ( (CDF_Z(6)-0.99_dp)*C_Tax_Inc_draft_group_z(:,6) + &
								& (0.999_dp-CDF_Z(6))*C_Tax_Inc_draft_group_z(:,7) )/0.009_dp
		C_Tax_Inc_draft_group(:,6)   = ( (CDF_Z(7)-0.999_dp)*C_Tax_Inc_draft_group_z(:,7) + &
								&  DBN_Z(8)*C_Tax_Inc_draft_group_z(:,8) + DBN_Z(9)*C_Tax_Inc_draft_group_z(:,9) )/0.001_dp
		C_Tax_Inc_draft_group(:,7)   = ( (CDF_Z(7)-0.999_dp)*C_Tax_Inc_draft_group_z(:,7) + &
								& (0.9999_dp-CDF_Z(7))*C_Tax_Inc_draft_group_z(:,8) )/0.0009_dp
		C_Tax_Inc_draft_group(:,8)   = ( (CDF_Z(8)-0.9999_dp)*C_Tax_Inc_draft_group_z(:,8) + &
								& DBN_Z(9)*C_Tax_Inc_draft_group_z(:,9) )/0.0001_dp

				! Capital Income Share Tax adjusted by productivity group
		K_Inc_frac_draft_group(:,1)   = ( DBN_Z(1)*K_Inc_frac_draft_group_z(:,1) + DBN_Z(2)*K_Inc_frac_draft_group_z(:,2) + & 
								& DBN_Z(3)*K_Inc_frac_draft_group_z(:,3) + (0.40_dp-CDF_Z(3))*K_Inc_frac_draft_group_z(:,4) )/0.40_dp
		K_Inc_frac_draft_group(:,2)   = ( (CDF_Z(4)-0.40_dp)*K_Inc_frac_draft_group_z(:,4) + &
								& (0.80_dp-CDF_Z(4))*K_Inc_frac_draft_group_z(:,5) )/0.40_dp
		K_Inc_frac_draft_group(:,3)   = K_Inc_frac_draft_group_z(:,5)
		K_Inc_frac_draft_group(:,4)   = ( (CDF_Z(5)-0.90_dp)*K_Inc_frac_draft_group_z(:,5) + &
								& (0.99_dp-CDF_Z(5))*K_Inc_frac_draft_group_z(:,6) )/0.09_dp
		K_Inc_frac_draft_group(:,5)   = ( (CDF_Z(6)-0.99_dp)*K_Inc_frac_draft_group_z(:,6) + &
								& (0.999_dp-CDF_Z(6))*K_Inc_frac_draft_group_z(:,7) )/0.009_dp
		K_Inc_frac_draft_group(:,6)   = ( (CDF_Z(7)-0.999_dp)*K_Inc_frac_draft_group_z(:,7) + &
								&  DBN_Z(8)*K_Inc_frac_draft_group_z(:,8) + DBN_Z(9)*K_Inc_frac_draft_group_z(:,9) )/0.001_dp
		K_Inc_frac_draft_group(:,7)   = ( (CDF_Z(7)-0.999_dp)*K_Inc_frac_draft_group_z(:,7) + &
								& (0.9999_dp-CDF_Z(7))*K_Inc_frac_draft_group_z(:,8) )/0.0009_dp
		K_Inc_frac_draft_group(:,8)   = ( (CDF_Z(8)-0.9999_dp)*K_Inc_frac_draft_group_z(:,8) + &
								& DBN_Z(9)*K_Inc_frac_draft_group_z(:,9) )/0.0001_dp

		! Labor Income Share Tax adjusted by productivity group
		L_Inc_frac_draft_group(:,1)   = ( DBN_Z(1)*L_Inc_frac_draft_group_z(:,1) + DBN_Z(2)*L_Inc_frac_draft_group_z(:,2) + & 
								& DBN_Z(3)*L_Inc_frac_draft_group_z(:,3) + (0.40_dp-CDF_Z(3))*L_Inc_frac_draft_group_z(:,4) )/0.40_dp
		L_Inc_frac_draft_group(:,2)   = ( (CDF_Z(4)-0.40_dp)*L_Inc_frac_draft_group_z(:,4) + &
								& (0.80_dp-CDF_Z(4))*L_Inc_frac_draft_group_z(:,5) )/0.40_dp
		L_Inc_frac_draft_group(:,3)   = L_Inc_frac_draft_group_z(:,5)
		L_Inc_frac_draft_group(:,4)   = ( (CDF_Z(5)-0.90_dp)*L_Inc_frac_draft_group_z(:,5) + &
								& (0.99_dp-CDF_Z(5))*L_Inc_frac_draft_group_z(:,6) )/0.09_dp
		L_Inc_frac_draft_group(:,5)   = ( (CDF_Z(6)-0.99_dp)*L_Inc_frac_draft_group_z(:,6) + &
								& (0.999_dp-CDF_Z(6))*L_Inc_frac_draft_group_z(:,7) )/0.009_dp
		L_Inc_frac_draft_group(:,6)   = ( (CDF_Z(7)-0.999_dp)*L_Inc_frac_draft_group_z(:,7) + &
								&  DBN_Z(8)*L_Inc_frac_draft_group_z(:,8) + DBN_Z(9)*L_Inc_frac_draft_group_z(:,9) )/0.001_dp
		L_Inc_frac_draft_group(:,7)   = ( (CDF_Z(7)-0.999_dp)*L_Inc_frac_draft_group_z(:,7) + &
								& (0.9999_dp-CDF_Z(7))*L_Inc_frac_draft_group_z(:,8) )/0.0009_dp
		L_Inc_frac_draft_group(:,8)   = ( (CDF_Z(8)-0.9999_dp)*L_Inc_frac_draft_group_z(:,8) + &
								& DBN_Z(9)*L_Inc_frac_draft_group_z(:,9) )/0.0001_dp
		L_Inc_frac_draft_group(:,6)   = ( (CDF_Z(7)-0.999_dp)*L_Inc_frac_draft_group_z(:,7) + &
								&  DBN_Z(8)*L_Inc_frac_draft_group_z(:,8) + DBN_Z(9)*L_Inc_frac_draft_group_z(:,9) )/0.001_dp
		L_Inc_frac_draft_group(:,7)   = ( (CDF_Z(7)-0.999_dp)*L_Inc_frac_draft_group_z(:,7) + &
								& (0.9999_dp-CDF_Z(7))*L_Inc_frac_draft_group_z(:,8) )/0.0009_dp
		L_Inc_frac_draft_group(:,8)   = ( (CDF_Z(8)-0.9999_dp)*L_Inc_frac_draft_group_z(:,8) + &
								& DBN_Z(9)*L_Inc_frac_draft_group_z(:,9) )/0.0001_dp

		! Size of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		wealth_draft_group(:,1) = wealth_draft_group_z(:,1) + wealth_draft_group_z(:,2) + wealth_draft_group_z(:,3) & 
								&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*wealth_draft_group_z(:,4)
		wealth_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*wealth_draft_group_z(:,4) +& 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*wealth_draft_group_z(:,5)
		wealth_draft_group(:,3) = (0.10_dp/DBN_Z(5))*wealth_draft_group_z(:,5)
		wealth_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*wealth_draft_group_z(:,5) + & 
								&  ((0.99_dp-CDF_Z(5))/DBN_Z(6))*wealth_draft_group_z(:,6)
		wealth_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*wealth_draft_group_z(:,6) + & 
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*wealth_draft_group_z(:,7) 
		wealth_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*wealth_draft_group_z(:,7) +&
								&  wealth_draft_group_z(:,8) + wealth_draft_group_z(:,9) 
		wealth_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*wealth_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*wealth_draft_group_z(:,8)
		wealth_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*wealth_draft_group_z(:,8) + wealth_draft_group_z(:,9)



		OPEN (UNIT=80, FILE=trim(Result_Folder)//'draft_group_Tax_K_exp.txt', STATUS='replace') 
	    OPEN (UNIT=81, FILE=trim(Result_Folder)//'draft_group_Tax_L_exp.txt', STATUS='replace') 
	    OPEN (UNIT=82, FILE=trim(Result_Folder)//'draft_group_Tax_C_exp.txt', STATUS='replace') 
	    OPEN (UNIT=83, FILE=trim(Result_Folder)//'draft_group_Inc_exp.txt', STATUS='replace') 
	    OPEN (UNIT=84, FILE=trim(Result_Folder)//'draft_group_Tax_K_Inc_exp.txt', STATUS='replace') 
	    OPEN (UNIT=85, FILE=trim(Result_Folder)//'draft_group_Tax_L_Inc_exp.txt', STATUS='replace') 
	    OPEN (UNIT=86, FILE=trim(Result_Folder)//'draft_group_Tax_C_Inc_exp.txt', STATUS='replace') 
	    OPEN (UNIT=87, FILE=trim(Result_Folder)//'draft_group_K_Inc_exp.txt', STATUS='replace') 
	    OPEN (UNIT=88, FILE=trim(Result_Folder)//'draft_group_L_Inc_exp.txt', STATUS='replace') 
	    OPEN (UNIT=89, FILE=trim(Result_Folder)//'draft_group_K_Inc_frac_exp.txt', STATUS='replace') 
	    OPEN (UNIT=90, FILE=trim(Result_Folder)//'draft_group_L_Inc_frac_exp.txt', STATUS='replace')
	    OPEN (UNIT=91, FILE=trim(Result_Folder)//'draft_group_Wealth_exp.txt', STATUS='replace') 
		do age = 1,draft_age_category
		    WRITE  (UNIT=80, FMT=*)  K_Tax_draft_group(age,:)
		    WRITE  (UNIT=81, FMT=*)  L_Tax_draft_group(age,:)
		    WRITE  (UNIT=82, FMT=*)  C_Tax_draft_group(age,:)
		    WRITE  (UNIT=83, FMT=*)  Tot_Income_draft_group(age,:)
		    WRITE  (UNIT=84, FMT=*)  K_Tax_Inc_draft_group(age,:)
		    WRITE  (UNIT=85, FMT=*)  L_Tax_Inc_draft_group(age,:)
		    WRITE  (UNIT=86, FMT=*)  C_Tax_Inc_draft_group(age,:)
		    WRITE  (UNIT=87, FMT=*)  K_Inc_draft_group(age,:)
		    WRITE  (UNIT=88, FMT=*)  L_Inc_draft_group(age,:)
		    WRITE  (UNIT=89, FMT=*)  K_Inc_frac_draft_group(age,:)
		    WRITE  (UNIT=90, FMT=*)  L_Inc_frac_draft_group(age,:)
		    WRITE  (UNIT=91, FMT=*)  wealth_draft_group(age,:)
		ENDDO
		close(unit=80); close(unit=81); close(unit=82); close(unit=83); close(unit=84); close(unit=85)
		close(unit=86); close(unit=87); close(unit=88); close(unit=89); close(unit=90); close(unit=91); 


		!! Fraction of tax increses
		Tax_Increase_tk_draft_group_z = 0.0_dp ; Tax_Increase_tl_draft_group_z = 0.0_dp ;
		Tax_Increase_draft_group_z = 0.0_dp ; Tax_Rate_Increase_draft_group_z = 0.0_dp ;
		Tax_Rate_Increase_tk_draft_group_z = 0.0_dp ; Tax_Rate_Increase_tl_draft_group_z = 0.0_dp ;

		do zi  = 1,nz
		do age = 1,draft_age_category
	        do xi=1,nx
		    do ei=1,ne
		    do lambdai=1,nlambda
		    do ai=1,na
		    do age2=draft_age_limit(age)+1,draft_age_limit(age+1)

		    	! Capital taxes
		    	If ( K_Tax_exp(age2,ai,zi,lambdai,ei,xi) .gt. K_Tax_bench(age2,ai,zi,lambdai,ei,xi)) then
		    	Tax_Increase_tk_draft_group_z(age,zi) = Tax_Increase_tk_draft_group_z(age,zi) + DBN_bench(age2,ai,zi,lambdai,ei,xi)
		    	endif 

		    	! Labor taxes
		    	If ( L_Tax_exp(age2,ai,zi,lambdai,ei,xi) .gt. L_Tax_bench(age2,ai,zi,lambdai,ei,xi)) then
		    	Tax_Increase_tl_draft_group_z(age,zi) = Tax_Increase_tl_draft_group_z(age,zi) + DBN_bench(age2,ai,zi,lambdai,ei,xi)
		    	endif 

		    	! Total taxes
		    	If ( (K_Tax_exp(age2,ai,zi,lambdai,ei,xi)+L_Tax_exp(age2,ai,zi,lambdai,ei,xi)) .gt. &
		    	  &  (K_Tax_bench(age2,ai,zi,lambdai,ei,xi)+L_Tax_bench(age2,ai,zi,lambdai,ei,xi)) ) then
		    	Tax_Increase_draft_group_z(age,zi) = Tax_Increase_draft_group_z(age,zi) + DBN_bench(age2,ai,zi,lambdai,ei,xi)
		    	endif 

				! Capital tax rate
		    	If ( K_Tax_exp(age2,ai,zi,lambdai,ei,xi)/Income_exp(age2,ai,zi,lambdai,ei,xi) .gt. & 
		    	  &  K_Tax_bench(age2,ai,zi,lambdai,ei,xi)/Income_bench(age2,ai,zi,lambdai,ei,xi)) then
		    	Tax_Rate_Increase_tk_draft_group_z(age,zi) = Tax_Rate_Increase_tk_draft_group_z(age,zi) + & 
		    		&  DBN_bench(age2,ai,zi,lambdai,ei,xi)
		    	endif 

		    	! Labor tax rate
		    	If ( L_Tax_exp(age2,ai,zi,lambdai,ei,xi)/Income_exp(age2,ai,zi,lambdai,ei,xi) .gt. & 
		    	  &  L_Tax_bench(age2,ai,zi,lambdai,ei,xi)/Income_bench(age2,ai,zi,lambdai,ei,xi)) then
		    	Tax_Rate_Increase_tl_draft_group_z(age,zi) = Tax_Rate_Increase_tl_draft_group_z(age,zi) + & 
		    		& DBN_bench(age2,ai,zi,lambdai,ei,xi)
		    	endif 

		    	! Total tax rate
		    	If ( (K_Tax_exp(age2,ai,zi,lambdai,ei,xi)+L_Tax_exp(age2,ai,zi,lambdai,ei,xi))/Income_exp(age2,ai,zi,lambdai,ei,xi) &
		    		& .gt. (K_Tax_bench(age2,ai,zi,lambdai,ei,xi)+L_Tax_bench(age2,ai,zi,lambdai,ei,xi))/&
		    		& Income_bench(age2,ai,zi,lambdai,ei,xi) ) then
		    	Tax_Rate_Increase_draft_group_z(age,zi) = Tax_Rate_Increase_draft_group_z(age,zi)+DBN_bench(age2,ai,zi,lambdai,ei,xi)
		    	endif 
	    	enddo 
	    	enddo 
	    	enddo 
	    	enddo 
	    	enddo  
		enddo
		enddo

		! Frac. capital tax increase by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		Tax_Increase_tk_draft_group(:,1) = Tax_Increase_tk_draft_group_z(:,1) + Tax_Increase_tk_draft_group_z(:,2) + &
						&  Tax_Increase_tk_draft_group_z(:,3) + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*Tax_Increase_tk_draft_group_z(:,4)
		Tax_Increase_tk_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*Tax_Increase_tk_draft_group_z(:,4) + & 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*Tax_Increase_tk_draft_group_z(:,5)
		Tax_Increase_tk_draft_group(:,3) = (0.10_dp/DBN_Z(5))*Tax_Increase_tk_draft_group_z(:,5)
		Tax_Increase_tk_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*Tax_Increase_tk_draft_group_z(:,5) + & 
								& ((0.99_dp-CDF_Z(5))/DBN_Z(6))*Tax_Increase_tk_draft_group_z(:,6) 
		Tax_Increase_tk_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*Tax_Increase_tk_draft_group_z(:,6) + &
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*Tax_Increase_tk_draft_group_z(:,7) 
		Tax_Increase_tk_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Tax_Increase_tk_draft_group_z(:,7) + &
								& Tax_Increase_tk_draft_group_z(:,8) + Tax_Increase_tk_draft_group_z(:,9) 
		Tax_Increase_tk_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Tax_Increase_tk_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*Tax_Increase_tk_draft_group_z(:,8)
		Tax_Increase_tk_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*Tax_Increase_tk_draft_group_z(:,8) + &
								&  Tax_Increase_tk_draft_group_z(:,9)

			! Fix fractions
			Tax_Increase_tk_draft_group = 100*Tax_Increase_tk_draft_group/size_draft_group

		! Frac. labor tax increase by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		Tax_Increase_tl_draft_group(:,1) = Tax_Increase_tl_draft_group_z(:,1) + Tax_Increase_tl_draft_group_z(:,2) + &
						&  Tax_Increase_tl_draft_group_z(:,3) + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*Tax_Increase_tl_draft_group_z(:,4)
		Tax_Increase_tl_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*Tax_Increase_tl_draft_group_z(:,4) + & 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*Tax_Increase_tl_draft_group_z(:,5)
		Tax_Increase_tl_draft_group(:,3) = (0.10_dp/DBN_Z(5))*Tax_Increase_tl_draft_group_z(:,5)
		Tax_Increase_tl_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*Tax_Increase_tl_draft_group_z(:,5) + & 
								& ((0.99_dp-CDF_Z(5))/DBN_Z(6))*Tax_Increase_tl_draft_group_z(:,6) 
		Tax_Increase_tl_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*Tax_Increase_tl_draft_group_z(:,6) + &
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*Tax_Increase_tl_draft_group_z(:,7) 
		Tax_Increase_tl_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Tax_Increase_tl_draft_group_z(:,7) + &
								& Tax_Increase_tl_draft_group_z(:,8) + Tax_Increase_tl_draft_group_z(:,9) 
		Tax_Increase_tl_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Tax_Increase_tl_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*Tax_Increase_tl_draft_group_z(:,8)
		Tax_Increase_tl_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*Tax_Increase_tl_draft_group_z(:,8) + &
								&  Tax_Increase_tl_draft_group_z(:,9)

			! Fix fractions
			Tax_Increase_tl_draft_group = 100*Tax_Increase_tl_draft_group/size_draft_group

		! Frac. total tax increase by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		Tax_Increase_draft_group(:,1) = Tax_Increase_draft_group_z(:,1) + Tax_Increase_draft_group_z(:,2) + &
						&  Tax_Increase_draft_group_z(:,3) + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*Tax_Increase_draft_group_z(:,4)
		Tax_Increase_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*Tax_Increase_draft_group_z(:,4) + & 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*Tax_Increase_draft_group_z(:,5)
		Tax_Increase_draft_group(:,3) = (0.10_dp/DBN_Z(5))*Tax_Increase_draft_group_z(:,5)
		Tax_Increase_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*Tax_Increase_draft_group_z(:,5) + & 
								& ((0.99_dp-CDF_Z(5))/DBN_Z(6))*Tax_Increase_draft_group_z(:,6) 
		Tax_Increase_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*Tax_Increase_draft_group_z(:,6) + &
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*Tax_Increase_draft_group_z(:,7) 
		Tax_Increase_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Tax_Increase_draft_group_z(:,7) + &
								& Tax_Increase_draft_group_z(:,8) + Tax_Increase_draft_group_z(:,9) 
		Tax_Increase_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Tax_Increase_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*Tax_Increase_draft_group_z(:,8)
		Tax_Increase_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*Tax_Increase_draft_group_z(:,8) + &
								&  Tax_Increase_draft_group_z(:,9)

			! Fix fractions
			Tax_Increase_draft_group = 100*Tax_Increase_draft_group/size_draft_group

		! Frac. capital tax rate increase by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		Tax_Rate_Increase_tk_draft_group(:,1) = Tax_Rate_Increase_tk_draft_group_z(:,1) + Tax_Rate_Increase_tk_draft_group_z(:,2) + &
						&  Tax_Rate_Increase_tk_draft_group_z(:,3) + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*Tax_Rate_Increase_tk_draft_group_z(:,4)
		Tax_Rate_Increase_tk_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*Tax_Rate_Increase_tk_draft_group_z(:,4) + & 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*Tax_Rate_Increase_tk_draft_group_z(:,5)
		Tax_Rate_Increase_tk_draft_group(:,3) = (0.10_dp/DBN_Z(5))*Tax_Rate_Increase_tk_draft_group_z(:,5)
		Tax_Rate_Increase_tk_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*Tax_Rate_Increase_tk_draft_group_z(:,5) + & 
								& ((0.99_dp-CDF_Z(5))/DBN_Z(6))*Tax_Rate_Increase_tk_draft_group_z(:,6) 
		Tax_Rate_Increase_tk_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*Tax_Rate_Increase_tk_draft_group_z(:,6) + &
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*Tax_Rate_Increase_tk_draft_group_z(:,7) 
		Tax_Rate_Increase_tk_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Tax_Rate_Increase_tk_draft_group_z(:,7) + &
								& Tax_Rate_Increase_tk_draft_group_z(:,8) + Tax_Rate_Increase_tk_draft_group_z(:,9) 
		Tax_Rate_Increase_tk_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Tax_Rate_Increase_tk_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*Tax_Rate_Increase_tk_draft_group_z(:,8)
		Tax_Rate_Increase_tk_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*Tax_Rate_Increase_tk_draft_group_z(:,8) + &
								&  Tax_Rate_Increase_tk_draft_group_z(:,9)

			! Fix fractions
			Tax_Rate_Increase_tk_draft_group = 100*Tax_Rate_Increase_tk_draft_group/size_draft_group

		! Frac. labor tax rate increase by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		Tax_Rate_Increase_tl_draft_group(:,1) = Tax_Rate_Increase_tl_draft_group_z(:,1) + Tax_Rate_Increase_tl_draft_group_z(:,2) + &
						&  Tax_Rate_Increase_tl_draft_group_z(:,3) + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*Tax_Rate_Increase_tl_draft_group_z(:,4)
		Tax_Rate_Increase_tl_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*Tax_Rate_Increase_tl_draft_group_z(:,4) + & 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*Tax_Rate_Increase_tl_draft_group_z(:,5)
		Tax_Rate_Increase_tl_draft_group(:,3) = (0.10_dp/DBN_Z(5))*Tax_Rate_Increase_tl_draft_group_z(:,5)
		Tax_Rate_Increase_tl_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*Tax_Rate_Increase_tl_draft_group_z(:,5) + & 
								& ((0.99_dp-CDF_Z(5))/DBN_Z(6))*Tax_Rate_Increase_tl_draft_group_z(:,6) 
		Tax_Rate_Increase_tl_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*Tax_Rate_Increase_tl_draft_group_z(:,6) + &
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*Tax_Rate_Increase_tl_draft_group_z(:,7) 
		Tax_Rate_Increase_tl_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Tax_Rate_Increase_tl_draft_group_z(:,7) + &
								& Tax_Rate_Increase_tl_draft_group_z(:,8) + Tax_Rate_Increase_tl_draft_group_z(:,9) 
		Tax_Rate_Increase_tl_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Tax_Rate_Increase_tl_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*Tax_Rate_Increase_tl_draft_group_z(:,8)
		Tax_Rate_Increase_tl_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*Tax_Rate_Increase_tl_draft_group_z(:,8) + &
								&  Tax_Rate_Increase_tl_draft_group_z(:,9)

			! Fix fractions
			Tax_Rate_Increase_tl_draft_group = 100*Tax_Rate_Increase_tl_draft_group/size_draft_group

		! Frac. total tax rate increase by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
		Tax_Rate_Increase_draft_group(:,1) = Tax_Rate_Increase_draft_group_z(:,1) + Tax_Rate_Increase_draft_group_z(:,2) + &
						&  Tax_Rate_Increase_draft_group_z(:,3) + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*Tax_Rate_Increase_draft_group_z(:,4)
		Tax_Rate_Increase_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*Tax_Rate_Increase_draft_group_z(:,4) + & 
								&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*Tax_Rate_Increase_draft_group_z(:,5)
		Tax_Rate_Increase_draft_group(:,3) = (0.10_dp/DBN_Z(5))*Tax_Rate_Increase_draft_group_z(:,5)
		Tax_Rate_Increase_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*Tax_Rate_Increase_draft_group_z(:,5) + & 
								& ((0.99_dp-CDF_Z(5))/DBN_Z(6))*Tax_Rate_Increase_draft_group_z(:,6) 
		Tax_Rate_Increase_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*Tax_Rate_Increase_draft_group_z(:,6) + &
								& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*Tax_Rate_Increase_draft_group_z(:,7) 
		Tax_Rate_Increase_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Tax_Rate_Increase_draft_group_z(:,7) + &
								& Tax_Rate_Increase_draft_group_z(:,8) + Tax_Rate_Increase_draft_group_z(:,9) 
		Tax_Rate_Increase_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*Tax_Rate_Increase_draft_group_z(:,7) + &
								& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*Tax_Rate_Increase_draft_group_z(:,8)
		Tax_Rate_Increase_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*Tax_Rate_Increase_draft_group_z(:,8) + &
								&  Tax_Rate_Increase_draft_group_z(:,9)

			! Fix fractions
			Tax_Rate_Increase_draft_group = 100*Tax_Rate_Increase_draft_group/size_draft_group

		OPEN (UNIT=80, FILE=trim(Result_Folder)//'draft_group_frac_Tax_K.txt', STATUS='replace') 
	    OPEN (UNIT=81, FILE=trim(Result_Folder)//'draft_group_frac_Tax_L.txt', STATUS='replace') 
	    OPEN (UNIT=82, FILE=trim(Result_Folder)//'draft_group_frac_Tax.txt', STATUS='replace') 
	    OPEN (UNIT=83, FILE=trim(Result_Folder)//'draft_group_frac_Tax_Rate_K.txt', STATUS='replace') 
	    OPEN (UNIT=84, FILE=trim(Result_Folder)//'draft_group_frac_Tax_Rate_L.txt', STATUS='replace') 
	    OPEN (UNIT=85, FILE=trim(Result_Folder)//'draft_group_frac_Tax_Rate.txt', STATUS='replace') 
		do age = 1,draft_age_category
		    WRITE  (UNIT=80, FMT=*)  Tax_Increase_tk_draft_group(age,:)
		    WRITE  (UNIT=81, FMT=*)  Tax_Increase_tl_draft_group(age,:)
		    WRITE  (UNIT=82, FMT=*)  Tax_Increase_draft_group(age,:)
		    WRITE  (UNIT=83, FMT=*)  Tax_Rate_Increase_tk_draft_group(age,:)
		    WRITE  (UNIT=84, FMT=*)  Tax_Rate_Increase_tl_draft_group(age,:)
		    WRITE  (UNIT=85, FMT=*)  Tax_Rate_Increase_draft_group(age,:)
		ENDDO
		close(unit=80); close(unit=81); close(unit=82); close(unit=83); close(unit=84); close(unit=85)



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
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: Value_mat

	allocate( Value_mat(  MaxAge,na,nz,nlambda,ne,nx) )

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Routine Set Up
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Age Brackets
		draft_age_limit = [0, 1, 15, 30, 45, MaxAge ] 

	! Discount factor
		CumDiscountF(MaxAge)=1.0_DP
		DO age=MaxAge-1,1,-1
		    CumDiscountF(age) = 1.0_DP + beta * survP(age) *CumDiscountF(age+1) 
		ENDDO

	! Select relevant values of value function
		DO age=1,MaxAge
			Value_mat(age,:,:,:,:,:) = ValueFunction_tr(age,:,:,:,:,:,age)
		ENDDO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Measuring Consumption Equivalent Welfare
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Consumption Equivalent Welfare: CE 1
		DO age=1,MaxAge
			if (Log_Switch.eqv..true.) then 
		    	CE1_tr(age,:,:,:,:,:)= & 
		    		& exp((Value_mat(age,:,:,:,:,:)-ValueFunction_Bench(age,:,:,:,:,:))/CumDiscountF(age))-1.0_DP
		    else 
		    	CE1_tr(age,:,:,:,:,:)=(Value_mat(age,:,:,:,:,:)/ValueFunction_Bench(age,:,:,:,:,:)) &
                                				&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP
		    end if 
		ENDDO

		! Aggregates 
		CE1_nb_tr  = 100.0_DP*sum(CE1_tr(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:))/sum(DBN_bench(1,:,:,:,:,:))
		CE1_pop_tr = 100.0_DP*sum(CE1_tr*DBN_bench)/sum(DBN_bench)

	! Consumption Equivalent Welfare: CE 2
		CE2_nb_tr  = 100.0_dp * (( sum(Value_mat(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:)) / &
				&               sum(ValueFunction_Bench(1,:,:,:,:,:)*DBN_bench(1,:,:,:,:,:)) ) &
				& ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)
		CE2_pop_tr = 100*(( sum(Value_mat*DBN_bench) / sum(ValueFunction_Bench*DBN_bench)  ) &
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


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Draft Tables
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	frac_pos_welfare_draft_group_z = 0.0_dp 

	do zi  = 1,nz
	do age = 1,draft_age_category
		size_draft_group_z(age,zi) = &
			& sum(DBN_bench(draft_age_limit(age)+1:draft_age_limit(age+1),:,zi,:,:,:))

		CE_draft_group_z(age,zi) =  100* &
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
	do zi=1,nz 
		CDF_Z(zi) = sum(DBN_Z(1:zi))
	enddo 
	! print*,' '
	! print*,'DBN_Z=',DBN_Z
	! print*,'CDF_Z=',CDF_Z 
	! print*,' '

	! Size of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
	size_draft_group(:,1) = size_draft_group_z(:,1) + size_draft_group_z(:,2) + size_draft_group_z(:,3) & 
							&  + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*size_draft_group_z(:,4)
	size_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*size_draft_group_z(:,4)+((0.80_dp-CDF_Z(4))/DBN_Z(5))*size_draft_group_z(:,5)
	size_draft_group(:,3) = (0.10_dp/DBN_Z(5))*size_draft_group_z(:,5)
	size_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*size_draft_group_z(:,5) + &
						& ((0.99_dp-CDF_Z(5))/DBN_Z(6))*size_draft_group_z(:,6) 
	size_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*size_draft_group_z(:,6) + & 
						& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*size_draft_group_z(:,7) 
	size_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*size_draft_group_z(:,7) + size_draft_group_z(:,8) + size_draft_group_z(:,9) 
	size_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*size_draft_group_z(:,7) + &
							& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*size_draft_group_z(:,8)
	size_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*size_draft_group_z(:,8) + size_draft_group_z(:,9)

	! CE of groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
	CE_draft_group(:,1)   = ( DBN_Z(1)*CE_draft_group_z(:,1) + DBN_Z(2)*CE_draft_group_z(:,2) + & 
							& DBN_Z(3)*CE_draft_group_z(:,3) + (0.40_dp-CDF_Z(3))*CE_draft_group_z(:,4) )/0.40_dp
	CE_draft_group(:,2)   = ( (CDF_Z(4)-0.40_dp)*CE_draft_group_z(:,4) + (0.80_dp-CDF_Z(4))*CE_draft_group_z(:,5) )/0.40_dp
	CE_draft_group(:,3)   = CE_draft_group_z(:,5)
	CE_draft_group(:,4)   = ( (CDF_Z(5)-0.90_dp)*CE_draft_group_z(:,5) + (0.99_dp-CDF_Z(5))*CE_draft_group_z(:,6) )/0.09_dp
	CE_draft_group(:,5)   = ( (CDF_Z(6)-0.99_dp)*CE_draft_group_z(:,6) + (0.999_dp-CDF_Z(6))*CE_draft_group_z(:,7) )/0.009_dp
	CE_draft_group(:,6)   = ( (CDF_Z(7)-0.999_dp)*CE_draft_group_z(:,7) + &
							&  DBN_Z(8)*CE_draft_group_z(:,8) + DBN_Z(9)*CE_draft_group_z(:,9) )/0.001_dp
	CE_draft_group(:,7)   = ( (CDF_Z(7)-0.999_dp)*CE_draft_group_z(:,7) + (0.9999_dp-CDF_Z(7))*CE_draft_group_z(:,8) )/0.0009_dp
	CE_draft_group(:,8)   = ( (CDF_Z(8)-0.9999_dp)*CE_draft_group_z(:,8) + DBN_Z(9)*CE_draft_group_z(:,9) )/0.0001_dp

	! Frac. pos. welfare by groups adjusting by z group: 0%-40% - 40%-80% - 80%-90% - 90%-99% - 99%-99.9% - 99.9%-100% - (99.9%-99.99% - 99.99%-100%)
	frac_pos_welfare_draft_group(:,1) = frac_pos_welfare_draft_group_z(:,1) + frac_pos_welfare_draft_group_z(:,2) + &
							&  frac_pos_welfare_draft_group_z(:,3) + ((0.40_dp-CDF_Z(3))/DBN_Z(4))*frac_pos_welfare_draft_group_z(:,4)
	frac_pos_welfare_draft_group(:,2) = ((CDF_Z(4)-0.40_dp)/DBN_Z(4))*frac_pos_welfare_draft_group_z(:,4) + & 
							&  ((0.80_dp-CDF_Z(4))/DBN_Z(5))*frac_pos_welfare_draft_group_z(:,5)
	frac_pos_welfare_draft_group(:,3) = (0.10_dp/DBN_Z(5))*frac_pos_welfare_draft_group_z(:,5)
	frac_pos_welfare_draft_group(:,4) = ((CDF_Z(5)-0.90_dp)/DBN_Z(5))*frac_pos_welfare_draft_group_z(:,5) + & 
							& ((0.99_dp-CDF_Z(5))/DBN_Z(6))*frac_pos_welfare_draft_group_z(:,6) 
	frac_pos_welfare_draft_group(:,5) = ((CDF_Z(6)-0.99_dp)/DBN_Z(6))*frac_pos_welfare_draft_group_z(:,6) + &
							& ((0.999_dp-CDF_Z(6))/DBN_Z(7))*frac_pos_welfare_draft_group_z(:,7) 
	frac_pos_welfare_draft_group(:,6) = ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*frac_pos_welfare_draft_group_z(:,7) + &
							& frac_pos_welfare_draft_group_z(:,8) + frac_pos_welfare_draft_group_z(:,9) 
	frac_pos_welfare_draft_group(:,7) =   ((CDF_Z(7)-0.999_dp)/DBN_Z(7))*frac_pos_welfare_draft_group_z(:,7) + &
							& ((0.9999_dp-CDF_Z(7))/DBN_Z(8))*frac_pos_welfare_draft_group_z(:,8)
	frac_pos_welfare_draft_group(:,8) = ((CDF_Z(8)-0.9999_dp)/DBN_Z(8))*frac_pos_welfare_draft_group_z(:,8) + &
							&  frac_pos_welfare_draft_group_z(:,9)

	! Fix fractions
	frac_pos_welfare_draft_group = 100*frac_pos_welfare_draft_group/size_draft_group

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


SUBROUTINE COMPUTE_VALUE_FUNCTION_LINEAR(Cons_mat,Hours_mat,Aprime_mat,Value_mat)
	IMPLICIT NONE
	REAL(DP), DIMENSION(MaxAge,na,nz,nlambda,ne,nx), INTENT(in)  :: Cons_mat, Hours_mat, Aprime_mat
	REAL(DP), DIMENSION(MaxAge,na,nz,nlambda,ne,nx), INTENT(out) :: Value_mat
	INTEGER  :: tklo, tkhi, xp_ind, age, xi, ai, zi, lambdai, ei
	REAL(DP) :: PrAprimelo, PrAprimehi, E_MU_cp(nx), aux

	! print*,'VALUE FUNCTION LINEAR'

	age=MaxAge
	DO xi=1,nx
	DO ai=1,na
    DO zi=1,nz
    DO lambdai=1,nlambda          
	DO ei=1,ne
      	Value_mat(age, ai, zi, lambdai, ei, xi) = &
      				& Utility(Cons_mat(age, ai, zi, lambdai, ei, xi),Hours_mat(age, ai, zi, lambdai, ei, xi))
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
			  & 				                     +  PrAprimehi*Value_mat(age+1, tkhi, zi, lambdai, ei, :)) ) 
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
			E_MU_cp(xp_ind) = sum( ( PrAprimelo * Value_mat(age+1, tklo, zi, lambdai,:,xp_ind)  &
		   & 					+    PrAprimehi * Value_mat(age+1, tkhi, zi, lambdai,:,xp_ind)) * pr_e(ei,:))
		enddo    

		Value_mat(age, ai, zi, lambdai, ei, xi) = Utility(Cons_mat(age,ai,zi,lambdai,ei,xi),Hours_mat(age,ai,zi,lambdai,ei,xi))  &
		   & + beta*survP(age)* sum( pr_x(xi,:,zi,age)*E_mu_cp )
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

END SUBROUTINE COMPUTE_VALUE_FUNCTION_LINEAR 


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE COMPUTE_VALUE_FUNCTION_TRANSITION
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: Cons_mat, Hours_mat, Aprime_mat, Value_mat
	INTEGER :: aux_T

	allocate( Cons_mat(   MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Hours_mat(  MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Aprime_mat( MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Value_mat(  MaxAge,na,nz,nlambda,ne,nx) )

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
			call COMPUTE_VALUE_FUNCTION_LINEAR(Cons_mat,Hours_mat,Aprime_mat,Value_mat)

		! Save value function
			ValueFunction_tr(:,:,:,:,:,:,age) = Value_mat
	enddo 

END SUBROUTINE COMPUTE_VALUE_FUNCTION_TRANSITION




!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE GOVNT_BUDGET(print_flag)
	IMPLICIT NONE
	LOGICAL, INTENT(IN) :: print_flag

	! Initialize variables
	GBAR 		 = 0.0_DP
	GBAR_K 		 = 0.0_DP
	GBAR_W		 = 0.0_DP
	GBAR_L 		 = 0.0_DP
	GBAR_C       = 0.0_DP
	SSC_Payments = 0.0_DP
	Tot_Lab_Inc  = 0.0_DP

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
	    GBAR = GBAR + DBN1(age,ai,zi,lambdai,ei,xi) * ( tauK*( R*agrid(ai) + Pr_mat(ai,zi,xi) )  	    &
	          & + ( agrid(ai) + ( R*agrid(ai) + Pr_mat(ai,zi,xi) ) *(1.0_DP-tauK)  ) - YGRID(ai,zi,xi)  &	
	          & + yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei,xi)  									&
	          & - psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei,xi))**(1.0_DP-tauPL)  			&
	          & + tauC * cons(age, ai, zi, lambdai,ei,xi)  )   

	    GBAR_L = GBAR_L  + DBN1(age,ai,zi,lambdai,ei,xi) * (  yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei,xi) &
	               &- psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei,xi))**(1.0_DP-tauPL) )

	    GBAR_K = GBAR_K + DBN1(age,ai,zi,lambdai,ei,xi) * tauK*( R*agrid(ai) + Pr_mat(ai,zi,xi) )

	    GBAR_W = GBAR_W + DBN1(age,ai,zi,lambdai,ei,xi) * &
	            & (( agrid(ai) + ( R*agrid(ai) + Pr_mat(ai,zi,xi) ) *(1.0_DP-tauK)  ) - YGRID(ai,zi,xi) )

      	GBAR_C = GBAR_C +  DBN1(age,ai,zi,lambdai,ei,xi) * tauC * cons(age, ai, zi, lambdai,ei,xi)

	    Tot_Lab_Inc = Tot_Lab_Inc + DBN1(age,ai,zi,lambdai,ei,xi) * yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei,xi) 
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
	print*, ' '
	print*, "Government Budget - Revenues and taxes"
	print*,'GBAR=',GBAR,'SSC_Payments=', SSC_Payments, 'GBAR_L=',GBAR_L,'Av. Labor Tax=', GBAR_L/Ebar 
	print*, 'GBAR_K=', GBAR_K, "GBAR_W=", GBAR_W, 'GBAR_C=', GBAR_C
	print*, 'Tau_K=', tauK, 'Tau_W_bt=', tauW_bt, 'Tau_W_at=', tauW_at, 'Tau_C=', tauC, "Threshold", Y_a_threshold
	print*, 'Tax Revenue over GDP', (GBAR_K+GBAR_W+GBAR_L+GBAR_C)/YBAR
	print*, 'Capital Tax / Total Tax', GBAR_K/(GBAR_K+GBAR_W+GBAR_L+GBAR_C)
	print*, 'Labor Tax / Total Tax', GBAR_L/(GBAR_K+GBAR_W+GBAR_L+GBAR_C)
	print*, 'Labor Tax / GDP', GBAR_L/YBAR
	print*, 'Average Labor Tax', GBAR_L/Tot_Lab_Inc
	print*, 'Total Labor Income', Tot_Lab_Inc , 'EBAR', EBAR
	print*, ' '

	if (solving_bench.eq.1) then 
	OPEN(UNIT=11, FILE=trim(Result_Folder)//'Govnt_Budget_Bench.txt', STATUS='replace')
	else
	OPEN(UNIT=11, FILE=trim(Result_Folder)//'Govnt_Budget_Exp.txt', STATUS='replace')
	endif 
	WRITE(UNIT=11, FMT=*) ' '
	WRITE(UNIT=11, FMT=*) "Government Budget - Revenues and taxes"
	WRITE(UNIT=11, FMT=*) 'GBAR=',GBAR,'SSC_Payments=', SSC_Payments, 'GBAR_L=',GBAR_L,'Av. Labor Tax=', GBAR_L/Ebar 
	WRITE(UNIT=11, FMT=*) 'GBAR_K=', GBAR_K, "GBAR_W=", GBAR_W, 'GBAR_C=', GBAR_C
	WRITE(UNIT=11, FMT=*) 'Tau_K=', tauK, 'Tau_W=', tauW_at, 'Tau_C=', tauC, "Threshold", Y_a_threshold
	WRITE(UNIT=11, FMT=*) 'Tax Revenue over GDP', (GBAR_K+GBAR_W+GBAR_L+GBAR_C)/YBAR
	WRITE(UNIT=11, FMT=*) 'Capital Tax / Total Tax', GBAR_K/(GBAR_K+GBAR_W+GBAR_L+GBAR_C)
	WRITE(UNIT=11, FMT=*) 'Capital Tax / GDP', GBAR_K/YBAR
	WRITE(UNIT=11, FMT=*) 'Labor Tax / Total Tax', GBAR_L/(GBAR_K+GBAR_W+GBAR_L+GBAR_C)
	WRITE(UNIT=11, FMT=*) 'Labor Tax / GDP', GBAR_L/YBAR
	WRITE(UNIT=11, FMT=*) 'Average Labor Tax', GBAR_L/Tot_Lab_Inc
	WRITE(UNIT=11, FMT=*) 'Total Labor Income', Tot_Lab_Inc , 'EBAR', EBAR
	Close(UNIT=11)
	ENDIF  
END  SUBROUTINE GOVNT_BUDGET

!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE FIND_DBN_EQ()
	use omp_lib
	IMPLICIT NONE
	INTEGER:: tklo, tkhi, age1, age2, z1, z2, a1, a2, lambda1, lambda2, e1, e2, DBN_iter, simutime, iter_indx, x1, x2
	REAL   :: DBN_dist, DBN_criteria
	REAL(DP)   ::BBAR, MeanWealth, brent_value
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable ::  PrAprimelo, PrAprimehi, DBN2
	INTEGER , DIMENSION(:,:,:,:,:,:), allocatable ::  Aplo, Aphi

	allocate( PrAprimelo(   MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrAprimehi(   MaxAge,na,nz,nlambda,ne,nx) )
	allocate( DBN2(         MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Aplo(         MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Aphi(         MaxAge,na,nz,nlambda,ne,nx) )

	!$ call omp_set_num_threads(nz)
	DBN_criteria = 1.0E-07_DP

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
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	! Probablities are adjusted to lie in [0,1]
		PrAprimelo = min(PrAprimelo, 1.0_DP)
		PrAprimelo = max(PrAprimelo, 0.0_DP)
		PrAprimehi = min(PrAprimehi, 1.0_DP)
		PrAprimehi = max(PrAprimehi, 0.0_DP)

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
	        	DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1, a1, z1, lambda1, e1, x1)
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2)*pr_lambda(lambda1,lambda2)    *  PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
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
	        	DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1, a1, z1, lambda1, e1, x1)
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2)*pr_lambda(lambda1,lambda2)    *  PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
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
	        	DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1, a1, z1, lambda1, e1, x1)
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2)*pr_lambda(lambda1,lambda2)    *  PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
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

	    	!!
	    	print*, 'DBN_diff=', DBN_dist, 'R=',R,'P=',P
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
				ENDDO
				ENDDO
				ENDDO
				ENDDO
				ENDDO
				ENDDO
	
	        ! Probablities are adjusted to lie in [0,1]
		        PrAprimelo = min(PrAprimelo, 1.0_DP)
		        PrAprimelo = max(PrAprimelo, 0.0_DP)
		        PrAprimehi = min(PrAprimehi, 1.0_DP)
		        PrAprimehi = max(PrAprimehi, 0.0_DP)  

		    ! Reset counter for next update of policy functions
	        iter_indx=0
	    ENDIF
	    
	    iter_indx = iter_indx + 1
	    simutime  = simutime +1 
	 
	ENDDO ! WHILE

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
	REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne, nx) :: PrAprimelo, PrAprimehi, DBN2
	INTEGER,  DIMENSION(MaxAge, na, nz, nlambda, ne, nx) :: Aplo, Aphi
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
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	! Probablities are adjusted to lie in [0,1]
		PrAprimelo = min(PrAprimelo, 1.0_DP)
		PrAprimelo = max(PrAprimelo, 0.0_DP)
		PrAprimehi = min(PrAprimehi, 1.0_DP)
		PrAprimehi = max(PrAprimehi, 0.0_DP)

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
	        	DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1, a1, z1, lambda1, e1, x1)
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2)*pr_lambda(lambda1,lambda2)    *  PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
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
	        	DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1, a1, z1, lambda1, e1, x1)
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2)*pr_lambda(lambda1,lambda2)    *  PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
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
	        	DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1, a1, z1, lambda1, e1, x1)
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2)*pr_lambda(lambda1,lambda2)    *  PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
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
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: PrAprimelo, PrAprimehi, DBN2
	INTEGER , DIMENSION(:,:,:,:,:,:), allocatable :: Aplo, Aphi
	real(dp), dimension(na,nz,nx) :: YGRID_old

	! Allocate
	allocate( PrAprimehi( MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrAprimelo( MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Aplo(       MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Aphi(       MaxAge,na,nz,nlambda,ne,nx) )
	allocate( DBN2(       MaxAge,na,nz,nlambda,ne,nx) )

	!$ call omp_set_num_threads(nz)
	DBN_criteria = 1.0E-07_DP
	DBN1 = DBN_bench

		! Form YGRID for the capital income economy given interest rate "P"
			YGRID_old = YGRID 
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
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	! Probablities are adjusted to lie in [0,1]
		PrAprimelo = min(PrAprimelo, 1.0_DP)
		PrAprimelo = max(PrAprimelo, 0.0_DP)
		PrAprimehi = min(PrAprimehi, 1.0_DP)
		PrAprimehi = max(PrAprimehi, 0.0_DP)

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
	        	DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1, a1, z1, lambda1, e1, x1)
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2)*pr_lambda(lambda1,lambda2)    *  PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
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
	        	DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1, a1, z1, lambda1, e1, x1)
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2)*pr_lambda(lambda1,lambda2)    *  PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
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
	        	DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1, a1, z1, lambda1, e1, x1)
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2)*pr_lambda(lambda1,lambda2)    *  PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
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
	real(dp)   :: BBAR, MeanWealth, brent_value
	REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne, nx) :: PrAprimelo, PrAprimehi, DBN2
	INTEGER,  DIMENSION(MaxAge, na, nz, nlambda, ne, nx) :: Aplo, Aphi
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
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	! Probablities are adjusted to lie in [0,1]
		PrAprimelo = min(PrAprimelo, 1.0_DP)
		PrAprimelo = max(PrAprimelo, 0.0_DP)
		PrAprimehi = min(PrAprimehi, 1.0_DP)
		PrAprimehi = max(PrAprimehi, 0.0_DP)

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
	        	DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1, a1, z1, lambda1, e1, x1)
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2)*pr_lambda(lambda1,lambda2)    *  PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
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
	        	DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1, a1, z1, lambda1, e1, x1)
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2)*pr_lambda(lambda1,lambda2)    *  PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
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
	        	DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1, a1, z1, lambda1, e1, x1)
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2)*pr_lambda(lambda1,lambda2)    *  PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
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

	    	!!
	    	print*,'Iter',simutime,'DBN_diff=',DBN_dist,'R=',R,'P=',P,'Assets',sum(sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid)
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
			        PrAprimelo(age,ai,zi,lambdai,ei,xi) = ( agrid(tkhi) - Aprime(age,ai,zi,lambdai,ei,xi) ) / ( agrid(tkhi) -agrid(tklo) )
			        PrAprimehi(age,ai,zi,lambdai,ei,xi) = ( Aprime(age,ai,zi,lambdai,ei,xi) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )
				ENDDO
				ENDDO
				ENDDO
				ENDDO
				ENDDO
				ENDDO
	
	        ! Probablities are adjusted to lie in [0,1]
		        PrAprimelo = min(PrAprimelo, 1.0_DP)
		        PrAprimelo = max(PrAprimelo, 0.0_DP)
		        PrAprimehi = min(PrAprimehi, 1.0_DP)
		        PrAprimehi = max(PrAprimehi, 0.0_DP)  

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
	real(dp)   ::BBAR, MeanWealth, brent_value
	REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne, nx) :: PrAprimelo, PrAprimehi, DBN2
	INTEGER,  DIMENSION(MaxAge, na, nz, nlambda, ne, nx) :: Aplo, Aphi

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
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	! Probablities are adjusted to lie in [0,1]
		PrAprimelo = min(PrAprimelo, 1.0_DP)
		PrAprimelo = max(PrAprimelo, 0.0_DP)
		PrAprimehi = min(PrAprimehi, 1.0_DP)
		PrAprimehi = max(PrAprimehi, 0.0_DP)

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
	        	DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1, a1, z1, lambda1, e1, x1)
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2)*pr_lambda(lambda1,lambda2)    *  PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
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
	        	DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1, a1, z1, lambda1, e1, x1)
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2)*pr_lambda(lambda1,lambda2)    *  PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
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
	        	DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1, a1, z1, lambda1, e1, x1)
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1) + DBN1(age1, a1, z1, lambda1, e1, x1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2)*pr_lambda(lambda1,lambda2)    *  PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
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

	    	!!
	    	print*, 'DBN_diff=',DBN_dist,'R=',R,'P=',P,'W=',wage,'Assets',sum(sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid)
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
				ENDDO
				ENDDO
				ENDDO
				ENDDO
				ENDDO
				ENDDO
	
	        ! Probablities are adjusted to lie in [0,1]
		        PrAprimelo = min(PrAprimelo, 1.0_DP)
		        PrAprimelo = max(PrAprimelo, 0.0_DP)
		        PrAprimehi = min(PrAprimehi, 1.0_DP)
		        PrAprimehi = max(PrAprimehi, 0.0_DP)  

		    ! Reset counter for next update of policy functions
	        iter_indx=0
	    ENDIF
	    
	    iter_indx = iter_indx + 1
	    simutime  = simutime +1 
	 
	ENDDO ! WHILE



END SUBROUTINE FIND_DBN_EQ_Prices



!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE FIND_DBN_Transition()
	use omp_lib
	IMPLICIT NONE
	INTEGER    :: tklo, tkhi, age1, age2, z1, z2, a1, a2, lambda1, lambda2, e1, e2, DBN_iter, simutime, iter_indx, x1, x2
	REAL       :: DBN_dist, DBN_criteria, Q_dist, N_dist, Price_criteria, Chg_criteria, Old_DBN_dist, Chg_dist
	REAL(DP)   :: BBAR, MeanWealth, brent_value, K_bench, C_bench, K_exp, C_exp
	REAL(DP), DIMENSION(:,:,:,:,:,:), allocatable :: PrAprimelo, PrAprimehi
	INTEGER , DIMENSION(:,:,:,:,:,:), allocatable :: Aplo, Aphi
	REAL(DP), DIMENSION(T+1) :: QBAR2_tr, NBAR2_tr


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Set Up
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Allocate
	allocate( PrAprimehi( MaxAge,na,nz,nlambda,ne,nx) )
	allocate( PrAprimelo( MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Aplo(       MaxAge,na,nz,nlambda,ne,nx) )
	allocate( Aphi(       MaxAge,na,nz,nlambda,ne,nx) )

	!$ call omp_set_num_threads(nz)
	DBN_criteria    = 1.0E-06_DP
	Price_criteria  = 1.0E-05_DP
	Chg_criteria    = 1.5E-07_DP

	! Set grids that depend on wealth tax threshold
		! Adjust agrid to include breaking points
		CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)

	! Wealth and consumption in benchmark and experiment
		K_bench = sum( sum(sum(sum(sum(sum(DBN_bench,6),5),4),3),1)*agrid )
		K_exp   = sum( sum(sum(sum(sum(sum(DBN_exp  ,6),5),4),3),1)*agrid )
		C_bench = sum( DBN_bench*Cons_bench )
		C_exp   = sum( DBN_exp  *Cons_exp   )


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Initial guess for transition path variables
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Guess NBAR, QBAR and R as a linear combination of starting and end values
			NBAR_tr(1)   = NBAR_bench ; NBAR_tr(T+1) = NBAR_exp   ;
			QBAR_tr(1)   = QBAR_bench ; QBAR_tr(T+1) = QBAR_exp   ;
			R_tr(1)      = R_bench    ; R_tr(T+1)    = R_exp      ;
			do ti=2,T
				NBAR_tr(ti) = NBAR_tr(ti-1) + (NBAR_tr(T+1)-NBAR_tr(1))/T
				QBAR_tr(ti) = QBAR_tr(ti-1) + (QBAR_tr(T+1)-QBAR_tr(1))/T
				R_tr(ti)    = R_tr(ti-1) + (R_tr(T+1)-R_tr(1))/T
			enddo 
		! Load Guess From Files
			! print*, 'Loading initial variables from file'
			! OPEN (UNIT=1,  FILE=trim(Result_Folder)//'QBAR_tr'   , STATUS='old', ACTION='read')
			! OPEN (UNIT=2,  FILE=trim(Result_Folder)//'NBAR_tr'	 , STATUS='old', ACTION='read')
			! OPEN (UNIT=3,  FILE=trim(Result_Folder)//'R_tr'		 , STATUS='old', ACTION='read')
			! READ (UNIT=1,  FMT=*), QBAR_tr
			! READ (UNIT=2,  FMT=*), NBAR_tr
			! READ (UNIT=3,  FMT=*), R_tr
			! CLOSE (unit=1); CLOSE (unit=2); CLOSE (unit=3);
			! print*, 'Reading completed'

		! Choose YBAR, EBAR, P and Wage to be consistent
		P_tr    = alpha* QBAR_tr**(alpha-mu) * NBAR_tr**(1.0_DP-alpha)
        YBAR_tr = QBAR_tr ** alpha * NBAR_tr **(1.0_DP-alpha)
        wage_tr = (1.0_DP-alpha)*QBAR_tr**alpha * NBAR_tr**(-alpha)
        Ebar_tr = wage_tr  * NBAR_tr  * sum(pop)/sum(pop(1:RetAge-1))
        	print*, "Test Prices"
	        print*, "P   ", P_bench, P_tr(1), P_tr(T), P_tr(T+1)
	        print*, "wage", wage_bench, wage_tr(1), wage_tr(T), wage_tr(T+1)
	        print*, "EBAR", EBAR_bench, EBAR_tr(1), EBAR_tr(T), EBAR_tr(T+1)
	        print*, "Tau_K", tauK,"Tau_W_at",tauW_at,"Tau_bt",tauW_bt

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
     		
     		WRITE(UNIT=77, FMT=*) 'Iteration, DBN_dist, Q_dist, N_dist, Q(T)/Q(SS), N(T)/N(SS)'
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

        	WRITE(UNIT=77, FMT=*) NBAR_bench, NBAR_tr, NBAR_exp
        	WRITE(UNIT=78, FMT=*) QBAR_bench, QBAR_tr, QBAR_exp
        	WRITE(UNIT=79, FMT=*) R_bench, R_tr, R_exp

    	CLOSE (unit=76); CLOSE (unit=77); CLOSE (unit=78); CLOSE (unit=79);
    	CLOSE (unit=80); CLOSE (unit=81); CLOSE (unit=82); CLOSE (unit=83); CLOSE (unit=84); CLOSE (unit=85);
    	CLOSE (unit=86); CLOSE (unit=87);

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
	N_dist       = 1.0_DP
	Chg_dist     = 1.0_DP 
	Old_DBN_dist = 0.0_DP 
	simutime     = 1
	iter_indx    = 1
	!print*, 'Computing Equilibrium Distribution'
	DO WHILE ((DBN_dist.ge.DBN_criteria).and.(max(Q_dist,N_dist).ge.Price_criteria)&
			& .and.(simutime.le.25).and.(Chg_dist.ge.Chg_criteria) )
		! print*, 'DBN_dist=', DBN_dist

		! Start Q_dist and N_dist
		Q_dist = 0.0 ; N_dist = 0.0 ;

		! Start Debt to zero
		Debt_tr = 0.0_dp


	    ! Solve for policy functions by backwards induction and EGM
	    	! Output is policy functions for all times and all ages
	    	CALL EGM_Transition

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


		! Discretize policy function for assets (a') for current period
			! For each age and state vector bracket optimal a' between two grid points
			! When at that age and state the optimal decision is approximated by selecting one the grid points
			! The grid points are selected with probability proportional to their distance to the optimal a'
		! print*,' Discretizing Policy Functions'
		DO age=1,MaxAge
			! print*,' 		Age ',age
		!$omp parallel do private(lambdai,ei,ai,xi,tklo,tkhi)
		DO zi=1,nz
			! print*,' 			z ',zi 
		DO xi=1,nx
			! print*,' 				x ',xi 
		DO ai=1,na
			! print*,' 					a ',ai 
		DO lambdai=1,nlambda
			! print*,' 						lambda ',lambdai 
		DO ei=1, ne
			! print*, '  							e', ei
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
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		ENDDO
		! print*, ' Discretizing Policy Functions Completed'

		! Probablities are adjusted to lie in [0,1]
			PrAprimelo = min(PrAprimelo, 1.0_DP)
			PrAprimelo = max(PrAprimelo, 0.0_DP)
			PrAprimehi = min(PrAprimehi, 1.0_DP)
			PrAprimehi = max(PrAprimehi, 0.0_DP)


		! Everyone in MaxAge dies. Those who die, switch to z2, lambda2 and start at ne/2+1 and x=1
	    age1=MaxAge
	    DO x1=1,nx
	    DO z1=1,nz
	    DO a1=1,na
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN_tr(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1 , ti+1)   =  &
	           		& DBN_tr(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1 , ti+1) + &
	           		& DBN_tr(age1, a1, z1, lambda1, e1, x1,ti) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1,a1,z1,lambda1,e1,x1)
	            DBN_tr(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1 , ti+1)   =  &
	           		& DBN_tr(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1 , ti+1) + &
	           		&  DBN_tr(age1, a1, z1, lambda1, e1, x1,ti) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
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
	        	DBN_tr(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1,ti+1)   =  &
	           		& DBN_tr(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1,ti+1) + &
	           		& DBN_tr(age1, a1, z1, lambda1, e1, x1,ti) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1,a1,z1,lambda1,e1,x1)
	            DBN_tr(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1,ti+1)   =  &
	           		& DBN_tr(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1,ti+1) + &
	           		& DBN_tr(age1, a1, z1, lambda1, e1, x1,ti) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, and also e1 since they are retired
	        !e2=e1
	        DO x2=1,nx
		        DBN_tr(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2,ti+1) =  &
		        	& DBN_tr(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2,ti+1) + &
		        	& DBN_tr(age1, a1, z1, lambda1, e1,x1,ti) &
		            & * survP(age1) * PrAprimelo(age1,a1,z1,lambda1,e1,x1) * pr_x(x1,x2,z1,age1)     
		        DBN_tr(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2,ti+1) =  &
		          	& DBN_tr(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e1,x2,ti+1) + &
		          	& DBN_tr(age1, a1, z1, lambda1, e1,x1,ti) &
		            & * survP(age1) * PrAprimehi(age1,a1,z1,lambda1,e1,x1) * pr_x(x1,x2,z1,age1)  
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
	        	DBN_tr(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1,ti+1)   =  &
	           		& DBN_tr(1,Aplo(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1,ti+1) + & 
	           		& DBN_tr(age1, a1, z1, lambda1, e1, x1,ti) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimelo(age1,a1,z1,lambda1,e1,x1)
	            DBN_tr(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2,lambda2,ne/2+1 , 1,ti+1)   =  &
	           		& DBN_tr(1,Aphi(age1, a1, z1, lambda1, e1,x1), z2, lambda2, ne/2+1, 1,ti+1) + & 
	           		& DBN_tr(age1, a1, z1, lambda1, e1, x1,ti) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2) * PrAprimehi(age1,a1,z1,lambda1,e1,x1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, but switch to e2 and x2
	        DO x2=1,nx
	        DO e2=1,ne
				DBN_tr(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2,ti+1) =  &
	          		& DBN_tr(age1+1, Aplo(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2,ti+1) + &
	          		& DBN_tr(age1, a1, z1, lambda1, e1,x1,ti) &
	                & * survP(age1) * pr_e(e1,e2) * pr_x(x1,x2,z1,age1) * PrAprimelo(age1,a1,z1,lambda1,e1,x1)
	            DBN_tr(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2,ti+1) =  &
	          		& DBN_tr(age1+1, Aphi(age1, a1, z1, lambda1, e1,x1), z1,lambda1,e2,x2,ti+1) + &
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

	    ! Set global variables to current period  (Time: ti)
	    	R = R_tr(ti) ; P = P_tr(ti) ; wage = wage_tr(ti) ;
	    	K_mat  = K_Matrix(R,P)
			Pr_mat = Profit_Matrix(R,P)
			CALL ComputeLaborUnits(Ebar_tr(ti), wage_tr(ti)) 
			CALL FORM_Y_MB_GRID(YGRID,MBGRID,YGRID_t,MBGRID_t)
			DBN1  = DBN_tr(:,:,:,:,:,:,ti)
			Cons  = Cons_tr(:,:,:,:,:,:,ti)
			Hours = Hours_tr(:,:,:,:,:,:,ti)


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
			    Debt_tr             = Debt_tr*(1+R_tr(ti)) + (GBAR_bench-GBAR_tr(ti)) 


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

	        	! Get Q_dist and N_dist before dampening 
	        	Q_dist = max(Q_dist,abs(QBAR2_tr(ti)/QBAR_tr(ti)-1))
	        	N_dist = max(N_dist,abs(NBAR2_tr(ti)/NBAR_tr(ti)-1))

            	! Dampened Update of QBAR and NBAR
	        	QBAR_tr(ti)  = 0.7*QBAR_tr(ti) + 0.3*QBAR2_tr(ti)
	        	NBAR_tr(ti)  = 0.7*NBAR_tr(ti) + 0.3*NBAR2_tr(ti)

        	! Update other prices and quantities             
	        P_tr(ti)     = alpha* QBAR_tr(ti)**(alpha-mu) * NBAR_tr(ti)**(1.0_DP-alpha)
	        YBAR_tr(ti)  = QBAR_tr(ti)**alpha * NBAR_tr(ti)**(1.0_DP-alpha)
	        wage_tr(ti)  = (1.0_DP-alpha)*QBAR_tr(ti)**alpha * NBAR_tr(ti)**(-alpha)
	        Ebar_tr(ti)  = wage_tr(ti)  * NBAR_tr(ti) * sum(pop)/sum(pop(1:RetAge-1))

	    	! Solve for new R (that clears market under new guess for prices)
	    	if (sum(theta)/nz .gt. 1.0_DP) then
	    		! Set price 
	    		P = min(P_tr(ti),1.0_dp)
	    		! Set DBN1 as the distribution for the current period (Time: ti)
	    		DBN1  = DBN_tr(:,:,:,:,:,:,ti)
	            brent_value = brent(-0.1_DP,0.01_DP,10.0_DP,Agg_Debt, brent_tol,R_tr(ti))
            else
                R_tr(ti) = 0.0_DP
	        endif

	        ! Compute total assets and consumption
	        K_tr(ti) = sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid )
	        C_tr(ti) = sum( DBN1*Cons )


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
			    ! Note that Debt_tr is not computed for T+1 since that would increase interest payments


	    ! Compute prices and aggregates for the current period (Time: T+1)
	    ! print*,' Updating Prices and Quantities fot T+1'
    		! Compute aggregates with current distribution (Time: T+1)
	        QBAR2_tr(ti) =0.0
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

	        	! Get Q_dist and N_dist before dampening 
	        	Q_dist = max(Q_dist,abs(QBAR2_tr(T+1)/QBAR_tr(T+1)-1))
	        	N_dist = max(N_dist,abs(NBAR2_tr(T+1)/NBAR_tr(T+1)-1))

            	! Dampened Update of QBAR and NBAR
	        	QBAR_tr(T+1)  = 0.7*QBAR_tr(T+1) + 0.3*QBAR2_tr(T+1)
	        	NBAR_tr(T+1)  = 0.7*NBAR_tr(T+1) + 0.3*NBAR2_tr(T+1)

        	! Update other prices and quantities             
	        P_tr(T+1)     = alpha* QBAR_tr(T+1)**(alpha-mu) * NBAR_tr(T+1)**(1.0_DP-alpha)
	        YBAR_tr(T+1)  = QBAR_tr(T+1)**alpha * NBAR_tr(T+1)**(1.0_DP-alpha)
	        wage_tr(T+1)  = (1.0_DP-alpha)*QBAR_tr(T+1)**alpha * NBAR_tr(T+1)**(-alpha)
	        Ebar_tr(T+1)  = wage_tr(T+1)  * NBAR_tr(T+1) * sum(pop)/sum(pop(1:RetAge-1))

	    	! Solve for new R (that clears market under new guess for prices)
	    	if (sum(theta)/nz .gt. 1.0_DP) then
	    		! Set price 
	    		P = min(P_tr(T+1),1.0_dp)
	    		! Set DBN1 as the distribution for the current period (Time: T+1)
	    		DBN1  = DBN_tr(:,:,:,:,:,:,T+1)
	            brent_value = brent(-0.1_DP,0.01_DP,10.0_DP,Agg_Debt, brent_tol,R_tr(T+1))
            else
                R_tr(T+1) = 0.0_DP
	        endif

	        ! Compute total assets and consumption
	        K_tr(T+1) = sum( sum(sum(sum(sum(sum(DBN1,6),5),4),3),1)*agrid )
	        C_tr(T+1) = sum( DBN1*Cons )
	    

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
	        	WRITE(UNIT=77, FMT=*) NBAR_bench, NBAR2_tr, NBAR_exp
	        	WRITE(UNIT=78, FMT=*) QBAR_bench, QBAR2_tr, QBAR_exp
	        	WRITE(UNIT=79, FMT=*) R_bench, R_tr, R_exp
	        	WRITE(UNIT=80, FMT=*) GBAR_bench, GBAR_tr, GBAR_exp
	        	WRITE(UNIT=81, FMT=*) GBAR_K_tr
	        	WRITE(UNIT=82, FMT=*) GBAR_W_tr
	        	WRITE(UNIT=83, FMT=*) GBAR_L_tr
	        	WRITE(UNIT=84, FMT=*) GBAR_C_tr
	        	WRITE(UNIT=85, FMT=*) SSC_Payments_tr
	        	WRITE(UNIT=86, FMT=*) K_tr
	        	WRITE(UNIT=87, FMT=*) C_tr
	    	CLOSE (unit=77); CLOSE (unit=78); CLOSE (unit=79);
	    	CLOSE (unit=80); CLOSE (unit=81); CLOSE (unit=82); CLOSE (unit=83); CLOSE (unit=84); CLOSE (unit=85);
	    	CLOSE (unit=86); CLOSE (unit=87);


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
			! CLOSE (unit=1); CLOSE (unit=2); CLOSE (unit=3); 
			CLOSE (unit=77); CLOSE (unit=78); CLOSE (unit=79);
	    	CLOSE (unit=80); CLOSE (unit=81); CLOSE (unit=82); CLOSE (unit=83); 
	    	CLOSE (unit=84); CLOSE (unit=85); CLOSE (unit=86); CLOSE (unit=87);
			print*,' 	--------------------------------------'
			print*,' 	Variable Printing Completed'
			print*,' 	--------------------------------------'
			print*,' '

		! Compare distance to tax reform distribution
		    DBN_dist = maxval(abs(DBN_tr(:,:,:,:,:,:,T+1)-DBN_exp))
		    Chg_dist = abs(DBN_dist-Old_DBN_dist)
		    Old_DBN_dist = DBN_dist

		    print*, 'Iteration=',simutime,' DBN_diff=', DBN_dist,' Q_dist=',Q_dist,' N_dist=',N_dist,&
		    	&' Q(T)/Q(SS)=',100*(QBAR2_tr(T+1)/QBAR_exp-1),' N(T)/N(SS)=',100*(NBAR2_tr(T+1)/NBAR_exp-1),&
		    	&' Chg_dist=',Chg_dist
	    	print*,'	GBAR_exp=',GBAR_tr(T+1),'Debt/GDP=',Debt_tr/YBAR_tr(T+1),'Deficit=',GBAR_tr(T+1)-GBAR_bench-R_tr(T+1)*Debt_tr

	    	OPEN (UNIT=76, FILE=trim(Result_Folder)//'Transition_Distance.txt', STATUS='old', POSITION='append')
	    	WRITE(UNIT=76, FMT=*) simutime,DBN_dist,Q_dist,N_dist,100*(QBAR2_tr(T+1)/QBAR_exp-1),100*(NBAR2_tr(T+1)/NBAR_exp-1)
	    	CLOSE(UNIT=76)


	    simutime  = simutime +1 
	 
	ENDDO ! WHILE

	print*,'---------------------------------------------------'
	print*,' 	Transition Completed'
	print*,'---------------------------------------------------'
	print*,' '

	! Print Summary File
	OPEN  (UNIT=78,  FILE=trim(Result_Folder)//'Transition_Summary.txt'   , STATUS='replace')
	WRITE (UNIT=78,  FMT=*) 'Period, Q, N, R, Wage, Y, K, C, GBAR, GBAR_K, GBAR_W, GBAR_L, GBAR_C, SSC'
	WRITE (UNIT=78,  FMT=*) 'SS_1,',QBAR_bench,',',NBAR_bench,',',R_bench,',',wage_bench,',' & 
								&  ,Y_bench,',',K_bench,',',C_bench,',',GBAR_bench
	do ti=1,T+1
	WRITE (UNIT=78,  FMT=*) ti,',',QBAR_tr(ti),',',NBAR_tr(ti),',',R_tr(ti),',',Wage_tr(ti),',' & 
								& ,YBAR_tr(ti),',',K_tr(ti),',',C_tr(ti),',',GBAR_tr(ti),',' &
								& ,GBAR_K_tr(ti),',',GBAR_W_tr(ti),',',GBAR_L_tr(ti),',',GBAR_C_tr(ti),',',SSC_Payments_tr(ti)
	enddo 
	WRITE (UNIT=78,  FMT=*) 'SS_2,',QBAR_exp,',',NBAR_exp,',',R_exp,',',wage_exp,',' & 
								&  ,Y_exp,',',K_exp,',',C_exp,',',GBAR_exp
	CLOSE (UNIT=78);



END SUBROUTINE FIND_DBN_Transition

!========================================================================================
!========================================================================================
!========================================================================================

SUBROUTINE EGM_Transition()
	use omp_lib
	IMPLICIT NONE
	REAL(DP) :: brentvalue, C_foc, H_min, euler_power
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
	DO xi=1,nx
    DO lambdai=1,nlambda
    DO ei=1,ne
    DO ai=1,na_t
        Cons_t(age,ai,zi,lambdai,ei,xi) =  YGRID_t(ai,zi,xi) + RetY_lambda_e(lambdai,ei) 
	ENDDO ! ai
    ENDDO ! ei
	ENDDO ! lambdai
	ENDDO ! xi
	ENDDO ! zi
	Aprime_t(age, :, :, :, :,:) = 0.0_DP
	! print*,"Period",ti,"|Const_t-Const_bench|=",maxval(abs(Cons_t(MaxAge,:,:,:,:,:)-Cons_t_pr(MaxAge,:,:,:,:,:)))
	! print*,"Period",ti,"Test-Const_bench=",(YGRID_aux(50,5,2) + RetY_lambda_e_aux(3,3))/(1.0_DP+tauC)-Cons_bench(MaxAge,50,5,3,3,2)
	
	! Rest of retirement
	DO age=MaxAge-1,RetAge,-1
		! print*,' 	Age=',age
	!$omp parallel do private(lambdai,ei,ai,xi,xp_ind,EndoCons,EndoYgrid,sw,sort_ind,tempai, &
	!$omp& state_FOC,par_FOC,MB_aprime_t,EndoYgrid_sort)
    DO zi=1,nz
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
    				P = P_tr(ti+1); R = R_tr(ti+1) 
    				MB_aprime_t(xp_ind) = MB_a_bt(agrid_t(ai),zi,xp_ind)
    			endif
    		enddo 
    		EndoCons(ai)  =  (beta*survP(age)* 	&
    			& sum(pr_x(xi,:,zi,age)*MB_aprime_t*Cons_t_pr(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power)) )**euler_power
	        EndoYgrid(ai) = agrid_t(ai) +  EndoCons(ai) - RetY_lambda_e(lambdai,ei)
	        ! Consumption on endogenous grid and implied asset income under tauW_at
	        do xp_ind = 1,nx 	
    			if (pr_x(xi,xp_ind,zi,age)/pr_x(xi,xp_ind,zi,age)*abs(Wealth_mat(ai,zi,xp_ind)-Y_a_threshold).lt.1e-8) then 
    				P = P_tr(ti+1); R = R_tr(ti+1) 
    				MB_aprime_t(xp_ind) = MB_a_at(agrid_t(ai),zi,xp_ind)
    			endif
    		enddo 
	        EndoCons(na_t+sw)  = (beta*survP(age)*	&
	        	& sum(pr_x(xi,:,zi,age)*MB_aprime_t*Cons_t_pr(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power)) )**euler_power
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
	    		& sum(pr_x(xi,:,zi,age)*MBGRID_t(ai,zi,:)*Cons_t_pr(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power)))**euler_power
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

	!------RETIREMENT PERIOD ENDS------------------------------------------------------------
	!========================================================================================
	
	!========================================================================================
	!------Working Period Starts-------------------------------------------------------------

	DO age=RetAge-1,1,-1
		! print*,' 	Age=',age
	!$omp parallel do private(lambdai,ei,ai,xi,xp_ind,EndoCons,EndoHours,EndoYgrid,sw,sort_ind,tempai, &
	!$omp& C_foc,state_FOC,par_FOC,MB_aprime_t)
    DO zi=1,nz
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
		if (any(pr_x(xi,:,zi,age)/pr_x(xi,:,zi,age)*abs(Wealth_mat(ai,zi,:)-Y_a_threshold).lt.1e-8)) then 
			sw 			  = sw+1	
    		MB_aprime_t   = MBGRID_t(ai,zi,:)
	    	! Below threshold
    		do xp_ind = 1,nx 	
    			if (pr_x(xi,xp_ind,zi,age)/pr_x(xi,xp_ind,zi,age)*abs(Wealth_mat(ai,zi,xp_ind)-Y_a_threshold).lt.1e-8) then 
    				P = P_tr(ti+1); R = R_tr(ti+1) 
    				MB_aprime_t(xp_ind) = MB_a_bt(agrid_t(ai),zi,xp_ind)
    			endif
    		enddo 
			call EGM_Working_Period_Transition( MB_aprime_t , H_min , state_FOC , & 
			      & EndoCons(ai), EndoHours(ai) , EndoYgrid(ai)  )
			
			! Above threshold
			do xp_ind = 1,nx 	
    			if (pr_x(xi,xp_ind,zi,age)/pr_x(xi,xp_ind,zi,age)*abs(Wealth_mat(ai,zi,xp_ind)-Y_a_threshold).lt.1e-8) then 
    				P = P_tr(ti+1); R = R_tr(ti+1) 
    				MB_aprime_t(xp_ind) = MB_a_at(agrid_t(ai),zi,xp_ind)
    			endif
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

			C_euler = ( (beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_MU_cp))) **(1.0_dp/((1.0_dp-sigma)*gamma-1.0_dp))
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
			C_endo = 1.0_dp/( (beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_mu_cp)  ) **(1.0_dp/sigma) )
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
			    & *  beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_MU_cp)  )**(-1.0_DP/sigma)

			  H_endo = 1.0_DP - (1.0_DP-gamma)*C_endo/(gamma*psi*yh(age, lambdai,ei))   

			If (H_endo .lt. 0.0_DP) then
			    H_endo = 0.0_DP 
			    C_endo = ( beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_MU_cp) )**(1.0_DP/(gamma*(1.0_DP-sigma)-1.0_DP))
			endif 

			! print*,' '
			! print*,' 	Inside EGM'
			! print*,' 	Consumption_t+1=',Cons_t_pr(age+1,ai,zi,lambdai,:,:)
		else 
			! Separable Utility
			do xp_ind=1,nx
				E_MU_cp(xp_ind) = sum( pr_e(ei,:) * (Cons_t_pr(age+1,ai,zi,lambdai,:,xp_ind)**(-sigma)) )
			enddo
			C_endo  = 1.0_DP/( beta*survP(age)*sum(pr_x(xi,:,zi,age)*MB_in*E_MU_cp) )**(1.0_DP/sigma)

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


SUBROUTINE COMPUTE_STATS()
	use omp_lib

	IMPLICIT NONE
	INTEGER  :: prctile, group, i
	REAL(DP), DIMENSION(nz)    :: cdf_Gz_DBN, Capital_by_z 
	REAL(dp), DIMENSION(na,nz,nx) :: DBN_azx, Wealth_mat
	REAL(DP) :: MeanATReturn, StdATReturn, VarATReturn, MeanATReturn_by_z(nz), Mean_Capital
	REAL(DP) :: Std_k_Return,    Var_K_Return,    Mean_K_Return_by_z(nz)
	REAL(DP) :: Std_AT_K_Return, Var_AT_K_Return, Mean_AT_K_Return_by_z(nz)
	REAL(DP) :: A_Age(max_age_category), A_AZ(max_age_category,nz), A_W(3)
	REAL(DP) :: Ap_Age(max_age_category), Ap_AZ(max_age_category,nz), Ap_W(3)
	REAL(DP) :: Y_Age(max_age_category), Y_AZ(max_age_category,nz), Y_W(3)
	REAL(DP) :: S_Age(max_age_category), S_AZ(max_age_category,nz), S_W(3)
	REAL(DP) :: S_Rate_A_Age(max_age_category), S_Rate_A_AZ(max_age_category,nz), S_Rate_A_W(3)
	REAL(DP) :: S_Rate_Y_Age(max_age_category), S_Rate_Y_AZ(max_age_category,nz), S_Rate_Y_W(3)
	REAL(DP) :: size_Age(max_age_category), size_AZ(max_age_category,nz), size_W(3)
	real(DP) :: leverage_age_z(MaxAge,nz), size_by_age_z(MaxAge,nz), constrained_firms_age_z(MaxAge,nz)
	real(DP) :: constrained_firms_age(MaxAge), size_by_age(MaxAge)
	real(DP) :: FW_top_x(6),  prctile_FW(6), prctile_bq(7), a, b, c, CCDF_c
	REAL(DP), DIMENSION(:), allocatable :: DBN_vec, Firm_Wealth_vec, CDF_Firm_Wealth, BQ_vec, DBN_bq_vec, CDF_bq
	real(DP), dimension(:,:,:,:,:,:), allocatable :: Firm_Output, Firm_Profit, DBN_bq
	real(DP), dimension(:,:,:,:,:,:), allocatable :: Labor_Income, Total_Income, K_L_Income, K_T_Income
	integer , dimension(:,:,:,:,:,:), allocatable :: constrained_firm_ind
	real(DP) :: Frisch_Aux, Frisch_Aux_2
	integer , dimension(max_age_category+1) :: age_limit
	character(100) :: rowname

	allocate(DBN_vec(size(DBN1)))
	allocate(Firm_Wealth_vec(size(DBN1)))
	allocate(CDF_Firm_Wealth(size(DBN1)))
	allocate(BQ_vec(size(DBN1)))
	allocate(DBN_bq_vec(size(DBN1)))
	allocate(CDF_bq(size(DBN1)))
	allocate(Firm_Output(MaxAge,na,nz,nlambda,ne,nx))
	allocate(Firm_Profit(MaxAge,na,nz,nlambda,ne,nx))
	allocate(DBN_bq(MaxAge,na,nz,nlambda,ne,nx))
	allocate(Labor_Income(MaxAge,na,nz,nlambda,ne,nx))
	allocate(Total_Income(MaxAge,na,nz,nlambda,ne,nx))
	allocate(K_L_Income(MaxAge,na,nz,nlambda,ne,nx))
	allocate(K_T_Income(MaxAge,na,nz,nlambda,ne,nx))
	allocate(constrained_firm_ind(MaxAge,na,nz,nlambda,ne,nx))


	!$ call omp_set_num_threads(20)
	!$ print *, "OMP Test Message"
	
	! Age Brackets
		age_limit = [0, 5, 15, 25, 35, 45, 55, MaxAge ]




	! Distribution of Assets
		DO ai=1,na
		     pr_a_dbn(ai)          = sum(DBN1(:,ai,:,:,:,:)) 
		     cdf_a_dbn(ai)         = sum( pr_a_dbn(1:ai) )      
		     tot_a_by_grid(ai)     = sum(DBN1(:,ai,:,:,:,:) * agrid(ai) )
		     cdf_tot_a_by_grid(ai) = sum(tot_a_by_grid(1:ai))   
		!     print*, pr_a_dbn(ai), cdf_a_dbn(ai)
		ENDDO
		cdf_a_dbn = cdf_a_dbn + 1.0_DP - cdf_a_dbn(na)

		! Percentage of the population above wealth tax threshold
		! Compute distribution of agents by (a,z,x)
		DBN_azx = sum(sum(sum(DBN1,5),4),1)
		! Compute mean before tax wealth
		Wealth_mat = Wealth_Matrix(R,P)
		! Compute share of agents above threshold
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
		print*,''
		prct1_wealth  = 1.0_DP-cdf_tot_a_by_prctile(99)/cdf_tot_a_by_prctile(100)
		prct10_wealth = 1.0_DP-cdf_tot_a_by_prctile(90)/cdf_tot_a_by_prctile(100)
		prct20_wealth = 1.0_DP-cdf_tot_a_by_prctile(80)/cdf_tot_a_by_prctile(100)
		prct40_wealth = 1.0_DP-cdf_tot_a_by_prctile(60)/cdf_tot_a_by_prctile(100)
	

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


    ! Debt to GDP Ratio
    External_Debt_GDP = 0.0_DP
	DO xi=1,nx
	DO zi=1,nz
	DO ai=1,na
	    External_Debt_GDP = External_Debt_GDP + sum(DBN1(:, ai, zi, :, :,xi))*abs(K_mat(ai,zi,xi)-agrid(ai))
	ENDDO
	ENDDO
	ENDDO
	External_Debt_GDP = 0.5_dp*External_Debt_GDP / YBAR

	! Savings Rate
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
	size_Age = 0.0_dp
	size_AZ  = 0.0_dp
	size_W   = 0.0_dp
	DO age=1,MaxAge 

	    DO while (age.gt.age_limit(group+1))
	        group = group+1
	    ENDDO    
	 
	 	DO xi=1,nx
	    DO ai=1,na
        DO zi=1,nz
        DO lambdai=1,nlambda
        DO ei=1,ne
        	size_Age(group)   = size_Age(group)   + DBN1(age,ai,zi,lambdai,ei,xi)
        	size_AZ(group,zi) = size_AZ(group,zi) + DBN1(age,ai,zi,lambdai,ei,xi)

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
        		size_W(1) = size_W(1) + DBN1(age,ai,zi,lambdai,ei,xi)
        		A_W(1)    = A_W(1)    + DBN1(age,ai,zi,lambdai,ei,xi)*agrid(ai)
        		Ap_W(1)   = Ap_W(1)   + DBN1(age,ai,zi,lambdai,ei,xi)*Aprime(age,ai,zi,lambdai,ei,xi)
        		if (age.lt.RetAge) then 
        		Y_W(1)    = Y_W(1)    + DBN1(age,ai,zi,lambdai,ei,xi)*&
        					& (YGRID(ai,zi,xi)+ Y_h(Hours(age,ai,zi,lambdai,ei,xi),age,lambdai,ei,Wage))
        		else 
        		Y_W(1)    = Y_W(1)    + DBN1(age,ai,zi,lambdai,ei,xi)*(YGRID(ai,zi,xi)+ RetY_lambda_e(lambdai,ei))
        		endif 
        	else if  ((ai.gt.prctile_ai_ind(90)).and.(ai.le.prctile_ai_ind(99))) then
        		size_W(2) = size_W(2) + DBN1(age,ai,zi,lambdai,ei,xi)
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
        		size_W(3) = size_W(3) + DBN1(age,ai,zi,lambdai,ei,xi)
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

	! Leverage Ratio and fraction of constrainted firms 
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
			Firm_Output(age,ai,zi,lambdai,ei,xi) = xz_grid(xi,zi)*K_mat(ai,zi,xi)
			Firm_Profit(age,ai,zi,lambdai,ei,xi) = Pr_mat(ai,zi,xi)
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
		




	! Distribution of firm wealth
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
			FW_top_x_share(i) = 100*sum(Firm_Wealth_vec*DBN_vec,Firm_Wealth_vec>=c)/Mean_Firm_Wealth
			WRITE(UNIT=11, FMT=*) 100_dp*prctile_FW(i),FW_top_x(i),FW_top_x_share(i), CCDF_c
		enddo 

			CLOSE(UNIT=11)

	! Distribution of bequest
		Bequest_Wealth=0.0_DP
		DO xi=1,nx
		DO zi=1,nz
		DO ai=1,na
		DO lambdai=1,nlambda
		DO ei=1, ne
		   Bequest_Wealth = Bequest_Wealth  +   DBN1(1, ai, zi, lambdai, ei, xi) * agrid(ai)
		ENDDO
		ENDDO
		ENDDO    
		ENDDO 
		ENDDO  


		! Distribution of bequest (matrix)	
		do ai=1,MaxAge
			DBN_bq(age,:,:,:,:,:) = DBN1(age,:,:,:,:,:)*(1.0_DP-survP(age))
		enddo 
		DBN_bq = DBN_bq/sum(DBN_bq)
		
		! Vectorization
		DBN_bq_vec        = reshape(DBN_bq,(/size(DBN1)/))
		BQ_vec            = reshape(Aprime,(/size(DBN1)/))

		! Mean Bequest
		Mean_Bequest      = sum(BQ_vec*DBN_bq_vec)

		if (solving_bench.eq.1) then
			OPEN(UNIT=11, FILE=trim(Result_Folder)//'Bequest_Stats_Bench.txt', STATUS='replace')
		else
			OPEN(UNIT=11, FILE=trim(Result_Folder)//'Bequest_Stats_Exp.txt', STATUS='replace')
		end if 
			WRITE(UNIT=11, FMT=*) ' '
			WRITE(UNIT=11, FMT=*) 'Bequest_Stats'
			WRITE(UNIT=11, FMT=*) 'Mean_Bequest/Wealth= '		, Bequest_Wealth/MeanWealth 
			WRITE(UNIT=11, FMT=*) 'Mean_Bequest/PV_Wealth= '	, Bequest_Wealth/Mean_Firm_Wealth 
			WRITE(UNIT=11, FMT=*) 'Bequests_Above_Threshold= '	, Threshold_Share_bq
			WRITE(UNIT=11, FMT=*) 'Bequest_Revenue/YBAR= '		, 0
			WRITE(UNIT=11, FMT=*) 'Top_x% ','x_percentile ','x_percentile/YBAR'

		prctile_bq = (/0.90_dp, 0.70_dp, 0.5_dp, 0.30_dp, 0.10_dp, 0.02_dp, 0.01_dp/)
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
			do while ((abs(CCDF_c-prctile_bq(i))>0.0001_dp).and.(b-a>1e-8))
				if (CCDF_c<prctile_bq(i)) then 
					b = c 
					c = (a+b)/2.0_dp
					CCDF_c = sum(DBN_bq_vec,BQ_vec>=c)
				else 
					a = c 
					c = (a+b)/2.0_dp
					CCDF_c = sum(DBN_bq_vec,BQ_vec>=c)
				endif
				!print*, 'a',a,'c',c,'b',b,'CCDF',CCDF_c,'Error', abs(CCDF_c-prctile_bq(i))
			enddo 
			BQ_top_x(i) = c 
			WRITE(UNIT=11, FMT=*) 100_dp*prctile_bq(i),BQ_top_x(i),BQ_top_x(i)/YBAR, CCDF_C
		enddo 

			CLOSE(UNIT=11)


	! Frisch Elasticity 
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
		Frisch_Elasticity = Frisch_Elasticity/Size_Frisch
		Hours_Frisch 	  = Hours_Frisch/Size_Frisch
		Frisch_Elasticity_2 = (1.0_dp-tauPL)/( sigma/(1.0_dp-(1.0_dp-sigma)*gamma) * Hours_Frisch/(1-Hours_Frisch) - tauPL )
		print*, 'Frisch_Elasticity',Frisch_Elasticity,Frisch_Elasticity_2,Hours_Frisch, Size_Frisch
		



	!print*, 'MeanReturn=',MeanReturn, 'StdReturn=', StdReturn
	!print*,'MeanReturn_by_z=',MeanReturn_by_z

	SSE_Moments = (1.0-Wealth_Output/3.0_DP)**2.0_DP  + (1.0_DP-prct1_wealth/0.34_DP)**2.0_DP  + (prct10_wealth-0.69_DP)**2.0_DP &
                   & + (1.0_DP-Std_Log_Earnings_25_60 / 0.8_DP)**2.0_DP + (1.0_DP-meanhours_25_60/0.4_DP)**2.0_DP &
                   & + (1.0_DP-MeanReturn/0.069_DP)**2.0_DP
	!print*,''
	!print*,"Current parameters"
	!print*,'beta',beta,'rho_z',rho_z,'sigma_z',sigma_z_eps,'sigma_lam',sigma_lambda_eps,'phi',phi
	print*,"Statistics"
	print*,'Debt/GDP',External_Debt_GDP,'W/GDP',Wealth_Output,'Top 1% A',prct1_wealth,'Top 10% A',prct10_wealth
	print*,'STD Labor Earnings',Std_Log_Earnings_25_60,'Mean Labor (hours 25-60)',meanhours_25_60,'MeanReturn',MeanReturn
	print*,'PV_Wealth_Top_1%', FW_top_x_share(4), 'PV_Top_10%', FW_top_x_share(3)
	print*,'Z','Constrained_firms_by_z: ','Capital_high_shock','Capital_low_shock'
	do zi=1,nz
		print*, zi, & 
			100.0_dp*sum(constrained_firm_ind(:,:,zi,:,:,:)*DBN1(:,:,zi,:,:,:))/sum(DBN1(:,:,zi,:,:,:)), &
			100.0_dp*sum(constrained_firm_ind(:,:,zi,:,:,1)*DBN1(:,:,zi,:,:,1))/sum(DBN1(:,:,zi,:,:,1)), &
			100.0_dp*sum(constrained_firm_ind(:,:,zi,:,:,2)*DBN1(:,:,zi,:,:,2))/sum(DBN1(:,:,zi,:,:,2)), &
			(EBAR_data/(EBAR*0.727853584919652_dp))*(mu*P*xz_grid(1,zi)**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) , & 
			(EBAR_data/(EBAR*0.727853584919652_dp))*(mu*P*xz_grid(2,zi)**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu))  
	enddo 
	print*, 'Total Constrained', 100.0_dp*sum(constrained_firm_ind*DBN1)
	print*,'Moments',SSE_Moments 
	print*,''

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


END SUBROUTINE COMPUTE_STATS


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
			FW_top_x_share(i) = 100*sum(Firm_Wealth_vec*DBN_vec,Firm_Wealth_vec>=c)/Mean_Firm_Wealth
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

SUBROUTINE EGM_RETIREMENT_WORKING_PERIOD()
	use omp_lib
	IMPLICIT NONE
	REAL(DP) :: brentvalue, C_foc, H_min, euler_power
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
	DO xi=1,nx
    DO lambdai=1,nlambda
    DO ei=1,ne
    DO ai=1,na_t
        Cons_t(age,ai,zi,lambdai,ei,xi) =  YGRID_t(ai,zi,xi) + RetY_lambda_e(lambdai,ei) 
	ENDDO ! ai
    ENDDO ! ei
	ENDDO ! lambdai
	ENDDO ! xi
	ENDDO ! zi
	Aprime_t(age, :, :, :, :,:) = 0.0_DP
	
	! Rest of retirement
	DO age=MaxAge-1,RetAge,-1
	!$omp parallel do private(lambdai,ei,ai,xi,xp_ind,EndoCons,EndoYgrid,sw,sort_ind,tempai, &
	!$omp& state_FOC,par_FOC,MB_aprime_t,EndoYgrid_sort)
    DO zi=1,nz
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
    					& sum(pr_x(xi,:,zi,age)*MB_aprime_t*Cons_t(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power)) ) **euler_power
	        EndoYgrid(ai) = agrid_t(ai) +  EndoCons(ai) - RetY_lambda_e(lambdai,ei)
	        
	        ! Consumption on endogenous grid and implied asset income under tauW_at
	        do xp_ind = 1,nx 	
				MB_aprime_t(xp_ind) = MB_a_at(agrid_t(ai),zi,xp_ind)
    		enddo 
	        EndoCons(na_t+sw)  = (beta*survP(age)*	&
	        			& sum(pr_x(xi,:,zi,age)*MB_aprime_t*Cons_t(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power)) ) **euler_power
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
	    				& sum(pr_x(xi,:,zi,age)*MBGRID_t(ai,zi,:)*Cons_t(age+1,ai,zi,lambdai,ei,:)**(1.0_dp/euler_power)) ) **euler_power
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
	
		end if 

    	ENDDO ! ai

	    ! !$omp critical
		if (any(isnan(EndoCons))) then 
			print*, "isnan - Consumption endogenous"
			print*, age,zi,ai,lambdai,ei,xi
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
	! print*, mu, P, R_in, DepRate
	! print*, xz_grid(1,5:)
	! print*, (mu*P*xz_grid(1,5:)**mu/(R_in+DepRate))**(1.0_dp/(1.0_dp-mu))
	! print*, '------------',Wealth, Kd, Agg_Debt, R_in, P
	Agg_Debt = (Agg_Debt/Wealth)**2.0_dp


end Function Agg_Debt


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
	!integer :: ai, zi

	DO xi=1,nx 
	DO zi=1,nz
		DO ai=1,na
			TYGRID(ai,zi,xi)  = Y_a(agrid(ai),zi,xi)
			TMBGRID(ai,zi,xi) = MB_a(agrid(ai),zi,xi)
		ENDDO 
		if (Y_a_threshold.eq.0.0_dp) then
			TYGRID_t  = TYGRID
			TMBGRID_t = TMBGRID 
		else 
			DO ai=1,na_t
				TYGRID_t(ai,zi,xi)  = Y_a(agrid_t(ai),zi,xi)
				TMBGRID_t(ai,zi,xi) = MB_a(agrid_t(ai),zi,xi)
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

	! Transitory investment productivity x
		if (nx.gt.1) then
			print*, 'X probability '
			xgrid = (/x_hi , x_lo , x_0/)
			! Low z types stay in x=1 until retirement
				pr_x(1,1,1:4,:) = 0.97_dp 
				pr_x(1,2,1:4,:) = 0.00_dp 
				pr_x(1,3,1:4,:) = 0.03_dp
			! High z types have 5% probability of going from x=1 to x=2
				pr_x(1,1,5:nz,:) = 0.92_dp 
				pr_x(1,2,5:nz,:) = 0.05_dp 
				pr_x(1,3,5:nz,:) = 0.03_dp
			! x=2 goes to x=3 with probability 3%
				pr_x(2,1,:,:) = 0.00_dp 
				pr_x(2,2,:,:) = 0.97_dp 
				pr_x(2,3,:,:) = 0.03_dp
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
				print*, ' xgrid', xgrid
				print*, ' zgrid', zgrid 
				do xi=1,nx
				print*, 'xzgrid', xz_grid(xi,:)
				enddo
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
		print*,'Initial Distrbution', sum(DBN1)
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
		CALL COMPUTE_VALUE_FUNCTION_LINEAR(Cons,Hours,Aprime,ValueFunction)
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
		deallocate( YGRID_t, MBGRID_t, Cons_t, Hours_t, Aprime_t )

	! Output variable
	Tax_Reform_Welfare = Av_Util_NB 

	print*, 'Tax Reform: tauK=',tauK,'tauW=',tauW_at,'Gain',Tax_Reform_Welfare


end Function Tax_Reform_Welfare

! !========================================================================================
! !========================================================================================
! !========================================================================================



SUBROUTINE WRITE_VARIABLES(bench_indx)
	IMPLICIT NONE
	integer, intent(in) :: bench_indx
	integer :: prctile, status, zi

	
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
		OPEN (UNIT=19, FILE=trim(Result_Folder)//'output.txt', STATUS='old', POSITION='append', iostat=status) 
			if (status.ne.0) then 
			OPEN (UNIT=19, FILE=trim(Result_Folder)//'output.txt', STATUS='replace') 
			end if 
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
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) 'GBAR'	, GBAR
			WRITE(UNIT=19, FMT=*) 'QBAR'	, QBAR
			WRITE(UNIT=19, FMT=*) 'NBAR'	, NBAR
			WRITE(UNIT=19, FMT=*) 'EBAR'	, EBAR
			WRITE(UNIT=19, FMT=*) 'YBAR'	, YBAR
			WRITE(UNIT=19, FMT=*) 'CBAR'    , MeanCons
			WRITE(UNIT=19, FMT=*) 'P'		, 100.0_dp*P
			WRITE(UNIT=19, FMT=*) 'wage'	, wage
			WRITE(UNIT=19, FMT=*) 'R'		, 100.0_dp*R
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) "Moments:"
			WRITE(UNIT=19, FMT=*) " "
			WRITE(UNIT=19, FMT=*) "Debt_Output"		  		, External_Debt_GDP
			WRITE(UNIT=19, FMT=*) "Wealth_Output"		  	, Wealth_Output
			WRITE(UNIT=19, FMT=*) "Mean_Assets"				, MeanWealth
			WRITE(UNIT=19, FMT=*) "Bequest_Wealth"	        , Bequest_Wealth/MeanWealth 
			WRITE(UNIT=19, FMT=*) 'Wealth_held_by_Top_1%' 	, prct1_wealth
			WRITE(UNIT=19, FMT=*) 'Wealth_held_by_Top_10%'	, prct10_wealth
			WRITE(UNIT=19, FMT=*) 'Wealth_held_by_Top_20%'	, prct20_wealth
			WRITE(UNIT=19, FMT=*) 'Wealth_held_by_Top_40%'	, prct40_wealth
			WRITE(UNIT=19, FMT=*) 'Mean_Labor_Earnings'	  	, mean_log_earnings_25_60
			WRITE(UNIT=19, FMT=*) 'STD_Labor_Earnings'	  	, Std_Log_Earnings_25_60
			WRITE(UNIT=19, FMT=*) 'Mean_Labor_25_60'	   	, meanhours_25_60
			WRITE(UNIT=19, FMT=*) 'Mean_Return'				, 100.0_dp*MeanReturn
			WRITE(UNIT=19, FMT=*) 'Std_Return'				, StdReturn
			do zi=1,nz
			WRITE(UNIT=19, FMT=*) 'Mean_Return_by_z'		, 100.0_dp*MeanReturn_by_z(zi)
			enddo 
			WRITE(UNIT=19, FMT=*) 'Moments'				  	, SSE_Moments 
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) 'Present_Value_Wealth'
			WRITE(UNIT=19, FMT=*) "Mean_PV_Wealth"		    , Mean_Firm_Wealth
			WRITE(UNIT=19, FMT=*) 'PV_Wealth_Top_0.01%' 	, FW_top_x_share(6)
			WRITE(UNIT=19, FMT=*) 'PV_Wealth_Top_0.1%' 		, FW_top_x_share(5)
			WRITE(UNIT=19, FMT=*) 'PV_Wealth_Top_1%' 		, FW_top_x_share(4)
			WRITE(UNIT=19, FMT=*) 'PV_Wealth_Top_10%'		, FW_top_x_share(3)
			WRITE(UNIT=19, FMT=*) 'PV_Wealth_Top_20%'		, FW_top_x_share(2)
			WRITE(UNIT=19, FMT=*) 'PV_Wealth_Top_40%'		, FW_top_x_share(1)
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) 'Bequest'
			WRITE(UNIT=19, FMT=*) 'Mean_Bequest/Wealth'		, Bequest_Wealth/MeanWealth 
			WRITE(UNIT=19, FMT=*) 'Mean_Bequest/PV_Wealth'	, Bequest_Wealth/Mean_Firm_Wealth 
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) 'Labor'
			WRITE(UNIT=19, FMT=*) 'Fraction_of_workers'		, Size_Frisch/sum(DBN1(1:RetAge-1,:,:,:,:,:))
			WRITE(UNIT=19, FMT=*) 'Av.Hours'				, Hours_Frisch
			WRITE(UNIT=19, FMT=*) 'Av.Frisch_Elasticity'   	, Frisch_Elasticity
			WRITE(UNIT=19, FMT=*) 'Frisch_Elasticity_Av'   	, Frisch_Elasticity_2
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) 'Taxes'
			WRITE(UNIT=19, FMT=*) 'Tax_Rev/GDP'				, (GBAR_K+GBAR_W+GBAR_L+GBAR_C)/YBAR
			WRITE(UNIT=19, FMT=*) 'Capital_Tax/Total_Tax'	, (GBAR_K+GBAR_W)/(GBAR_K+GBAR_W+GBAR_L+GBAR_C)
			WRITE(UNIT=19, FMT=*) 'Capital_Tax/_GDP'		, (GBAR_K+GBAR_W)/YBAR
			WRITE(UNIT=19, FMT=*) 'Labor_Tax/Total_Tax'		, GBAR_L/(GBAR_K+GBAR_W+GBAR_L+GBAR_C)
			WRITE(UNIT=19, FMT=*) 'Labor_Tax/GDP'			, GBAR_L/YBAR
			WRITE(UNIT=19, FMT=*) 'Average_Labor_Tax'		, GBAR_L/Tot_Lab_Inc
			WRITE(UNIT=19, FMT=*) 'Cons_Tax/Total_Tax'		, GBAR_C/(GBAR_K+GBAR_W+GBAR_L+GBAR_C)
			WRITE(UNIT=19, FMT=*) 'Cons_Tax/GDP'			, GBAR_C/YBAR
			WRITE(UNIT=19, FMT=*) ' '
		CLOSE(Unit=19)
	if (bench_indx.ne.1) then
		OPEN (UNIT=19, FILE=trim(Result_Folder)//'output.txt', STATUS='old', POSITION='append') 
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) 'Welfare and output gain'
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) "CE_Pop(bench)" , Welfare_Gain_Pop_bench
			WRITE(UNIT=19, FMT=*) "CE_Pop(exp)"   , Welfare_Gain_Pop_exp
			WRITE(UNIT=19, FMT=*) "CE_NB(bench)"  , Welfare_Gain_NB_bench
			WRITE(UNIT=19, FMT=*) "CE_NB(exp)"    , Welfare_Gain_NB_exp
			WRITE(UNIT=19, FMT=*) "Output_Gain(prct)"	  	, 100.0_DP*(Y_exp/Y_bench-1.0) 
			WRITE(UNIT=19, FMT=*) "Av_Util_Pop(exp)"		, Av_Util_Pop
			WRITE(UNIT=19, FMT=*) "Av_Util_NB(exp)"			, Av_Util_NB
		CLOSE(Unit=19)
	end if 

			
END SUBROUTINE WRITE_VARIABLES

!========================================================================================
!========================================================================================
!========================================================================================

SUBROUTINE Write_Benchmark_Results(Compute_bench)
	IMPLICIT NONE
	logical :: Compute_bench
	character(100) :: bench_folder, string_theta

	write(string_theta,'(f4.2)')  theta_folder

	if ((Progressive_Tax_Switch.eqv..false.).and.(NSU_Switch.eqv..true.)) then 
		bench_folder = './NSU_ZS_LT_Results/Theta_'//trim(string_theta)//'/Bench_Files/'
	else if ((Progressive_Tax_Switch.eqv..true.).and.(NSU_Switch.eqv..true.)) then 
		bench_folder = './NSU_ZS_PT_Results/Theta_'//trim(string_theta)//'/Bench_Files/'
	else if ((Progressive_Tax_Switch.eqv..false.).and.(NSU_Switch.eqv..false.)) then 
		bench_folder = './SU_ZS_LT_Results/Theta_'//trim(string_theta)//'/Bench_Files/'
	else if ((Progressive_Tax_Switch.eqv..true.).and.(NSU_Switch.eqv..false.)) then 
		bench_folder = './SU_ZS_PT_Results/Theta_'//trim(string_theta)//'/Bench_Files/'
	end if 

	bench_folder = trim(Result_Folder)//'Bench_Files/'
		call system( 'mkdir -p ' // trim(bench_folder) )
		print*, "Bench Files Folder:", bench_folder
	
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

		print*, "Writing of benchmark results completed"
	ELSE 
		OPEN (UNIT=1,  FILE=trim(bench_folder)//'cons'  , STATUS='old', ACTION='read')
		OPEN (UNIT=2,  FILE=trim(bench_folder)//'aprime', STATUS='old', ACTION='read')
		OPEN (UNIT=3,  FILE=trim(bench_folder)//'hours' , STATUS='old', ACTION='read')
		OPEN (UNIT=4,  FILE=trim(bench_folder)//'value' , STATUS='old', ACTION='read')
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

		READ (UNIT=1,  FMT=*), cons
		READ (UNIT=2,  FMT=*), aprime
		READ (UNIT=3,  FMT=*), hours
		READ (UNIT=4,  FMT=*), ValueFunction
		READ (UNIT=5,  FMT=*), DBN1 
		READ (UNIT=60, FMT=*), GBAR 
		READ (UNIT=7,  FMT=*), EBAR
		READ (UNIT=8,  FMT=*), NBAR
		READ (UNIT=9,  FMT=*), QBAR
		READ (UNIT=10, FMT=*), P
		READ (UNIT=11, FMT=*), R
		READ (UNIT=12, FMT=*), wage 
		READ (UNIT=13, FMT=*), YBAR

		READ (UNIT=14, FMT=*), tauK
		READ (UNIT=15, FMT=*), tauPL
		READ (UNIT=16, FMT=*), psi
		READ (UNIT=17, FMT=*), tauW_bt
		READ (UNIT=18, FMT=*), tauW_at

		READ (UNIT=19, FMT=*), V_Pr
		! READ (UNIT=20, FMT=*), V_Pr_nb
		READ (UNIT=21, FMT=*), SSC_Payments
		READ (UNIT=22, FMT=*), Income_AT

		CLOSE (unit=1); CLOSE (unit=2); CLOSE (unit=3); CLOSE (unit=4); CLOSE (unit=5)
		CLOSE (unit=60); CLOSE (unit=7); CLOSE (unit=8); CLOSE (unit=9); CLOSE (unit=10)
		CLOSE (unit=11); CLOSE (unit=12); CLOSE (unit=13); CLOSE (unit=14); CLOSE (unit=15)
		CLOSE (unit=16); CLOSE (unit=17); CLOSE (unit=18); CLOSE (unit=19); !CLOSE (unit=20)
		CLOSE (unit=21); CLOSE (unit=22);

		print*, "Reading of benchmark results completed"
	END IF 
END SUBROUTINE Write_Benchmark_Results


SUBROUTINE Write_Experimental_Results(compute_exp)
	IMPLICIT NONE
	logical, intent(in) :: compute_exp

	call system( 'mkdir -p ' // trim(Result_Folder) // 'Exp_Files/' )

	if (compute_exp) then 
		print*, "Writing experimental results in folder", trim(Result_Folder) // 'Exp_Files/'
		OPEN  (UNIT=1,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_cons'  , STATUS='replace')
		WRITE (UNIT=1,  FMT=*) cons
		OPEN  (UNIT=2,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_aprime', STATUS='replace')
		WRITE (UNIT=2,  FMT=*) aprime
		OPEN  (UNIT=3,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_hours' , STATUS='replace')
		WRITE (UNIT=3,  FMT=*) hours
		OPEN  (UNIT=4,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_value' , STATUS='replace')
		WRITE (UNIT=4,  FMT=*) ValueFunction

		OPEN  (UNIT=5,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_DBN'   , STATUS='replace')
		WRITE (UNIT=5,  FMT=*) DBN1 
		OPEN  (UNIT=60,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_GBAR' , STATUS='replace')
		WRITE (UNIT=60,  FMT=*) GBAR
		OPEN  (UNIT=7,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_EBAR'  , STATUS='replace')
		WRITE (UNIT=7,  FMT=*) EBAR
		OPEN  (UNIT=8,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_NBAR'  , STATUS='replace')
		WRITE (UNIT=8,  FMT=*) NBAR
		OPEN  (UNIT=9,  FILE=trim(Result_Folder)//'Exp_Files/Exp_results_QBAR'  , STATUS='replace')
		WRITE (UNIT=9,  FMT=*) QBAR
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

		print*, "Writing of experimental results completed"

	else 
		print*, "Reading experimental results from folder", trim(Result_Folder) // 'Exp_Files/'
		OPEN (UNIT=1 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_cons'  	, STATUS='old', ACTION='read')
		OPEN (UNIT=2 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_aprime'	, STATUS='old', ACTION='read')
		OPEN (UNIT=3 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_hours' 	, STATUS='old', ACTION='read')
		OPEN (UNIT=4 , FILE=trim(Result_Folder)//'Exp_Files/Exp_results_value' 	, STATUS='old', ACTION='read')
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

		READ (UNIT=1,  FMT=*), cons
		READ (UNIT=2,  FMT=*), aprime
		READ (UNIT=3,  FMT=*), hours
		READ (UNIT=4,  FMT=*), ValueFunction
		READ (UNIT=5,  FMT=*), DBN1 
		READ (UNIT=60, FMT=*), GBAR 
		READ (UNIT=7,  FMT=*), EBAR
		READ (UNIT=8,  FMT=*), NBAR
		READ (UNIT=9,  FMT=*), QBAR
		READ (UNIT=10, FMT=*), P
		READ (UNIT=11, FMT=*), R
		READ (UNIT=12, FMT=*), wage 
		READ (UNIT=13, FMT=*), YBAR
		READ (UNIT=14, FMT=*), tauK
		READ (UNIT=15, FMT=*), psi
		READ (UNIT=16, FMT=*), tauPL
		READ (UNIT=17, FMT=*), tauW_bt
		READ (UNIT=18, FMT=*), tauW_at
		READ (UNIT=19, FMT=*), V_Pr
		READ (UNIT=20, FMT=*), V_Pr_nb
		READ (UNIT=21, FMT=*), Income_AT
		print*, "Reading of experimental results completed"
	endif 

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
	CLOSE (unit=13)
	CLOSE (unit=14)
	CLOSE (unit=15)
	CLOSE (unit=16)
	CLOSE (unit=17)
	CLOSE (unit=18)
	CLOSE (unit=19)
	CLOSE (unit=20); CLOSE (unit=21)

END SUBROUTINE Write_Experimental_Results


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



