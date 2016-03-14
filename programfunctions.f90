
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
	real(dp), intent(in)    :: Y_a_threshold_in
	integer , intent(out)   :: na_t
	real(dp), dimension(:), allocatable, intent(out) :: agrid_t
	real(dp), dimension(nz) :: a_aux
	integer                 :: a_ind
	integer , dimension(:), allocatable :: agrid_t_ind
	real(dp), dimension(:), allocatable :: par
	!real(dp), dimension(2)  :: par
	real(dp)                :: max_wealth, K

	allocate( par(2) )
	a_ind = 0
	! If the threshold for wealth taxes is positive then agrid is adjusted
	if (Y_a_threshold_in.gt.0.0_dp) then 
		par(1) = Y_a_threshold_in
 		do zi=1,nz 
			! New points are added to agrid if there is an "a" st Y(a,z))=Y_threshold
			K = min( theta*agrid(na) , (mu*P*zgrid(zi)/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )
			max_wealth = (1.0_dp+R)*agrid(na) + P*(zgrid(zi)*K)**mu - (R+DepRate)*K 
			if (Y_a_threshold_in.lt.max_wealth) then
				a_ind		 = a_ind + 1 
				par(2)       = zgrid(zi)
				a_aux(a_ind) = zbrent_p(Y_a_res,0.0_dp,agrid(na),brent_tol,par) 
				!a_aux(a_ind) = zbrent(Y_a_res,0.0_dp,agrid(na),brent_tol)
			else 
				print*, 'Error in forming a grid with threshold'
				STOP
			end if 
 		end do 

 		na_t = na + a_ind
 		allocate( agrid_t(1:na_t) )
 		allocate( agrid_t_ind(1:na_t) )
 		agrid_t = [agrid,a_aux(1:a_ind)]
 		call Sort(na_t,agrid_t,agrid_t,agrid_t_ind)
 	else 
 		na_t    = na
 		allocate( agrid_t(1:na_t) )
 		agrid_t = agrid
	end if

	! Allocate variables that depend on na_t
		! Grids for Y and MB
		allocate( YGRID_t(na_t,nz)  )
		allocate( MBGRID_t(na_t,nz) )
		! Allocate size of policy function on adjusted grid:
		allocate( Cons_t(MaxAge,na_t,nz,nlambda,ne) )
		allocate( Hours_t(MaxAge,na_t,nz,nlambda,ne) )
		allocate( Aprime_t(MaxAge,na_t,nz,nlambda,ne) )

	!contains 

		
end Subroutine Asset_Grid_Threshold

	function Y_a_res(a_in,par)
	!function Y_a_res(a_in)
		IMPLICIT NONE
		real(dp), intent(in)  :: a_in
		real(dp), dimension(:), allocatable, intent(in) :: par
		real(dp) :: Y_a_res
		real(dp) :: Y_a_th, z_in, K

		Y_a_th = par(1)
		z_in   = par(2)

		K = min( theta*a_in , (mu*P*z_in/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )

		Y_a_res = (1.0_dp+R)*a_in + P*(z_in*K)**mu - (R+DepRate)*K  - Y_a_th

	end function Y_a_res


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
	FUNCTION Y_a(a_in,z_in)
		IMPLICIT NONE   
		real(DP), intent(in)  :: a_in, z_in
		real(DP)              :: Y_a, K, Pr

		! Capital demand 
		K   = min( theta*a_in , (mu*P*z_in**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )
		! Profits 
		Pr  = P*(z_in*K)**mu - (R+DepRate)*K
		! Before tax wealth
		Y_a = ( a_in +  ( Pr + R*a_in ) *(1.0_DP-tauK) )

		! Compute after tax wealth according to threshold
		if (Y_a.le.Y_a_threshold) then 
			Y_a = Y_a* (1.0_dp-tauW_bt)
		else
			Y_a = Y_a_threshold*(1.0_dp-tauW_bt) + (Y_a - Y_a_threshold) * (1.0_dp-tauW_at)
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
	FUNCTION MB_a(a_in,z_in)
		IMPLICIT NONE   
		real(DP), intent(in)  :: a_in, z_in
		real(DP) 			  :: MB_a
		real(DP) :: K, Pr, Y_a, tauW

		! Capital demand 
		K   = min( theta*a_in , (mu*P*z_in**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )
		! Profits 
		Pr  = P*(z_in*K)**mu - (R+DepRate)*K
		! Before tax wealth
		Y_a = ( a_in +  ( Pr + R*a_in ) *(1.0_DP-tauK) )
		if (Y_a.le.Y_a_threshold) then 
			tauW = tauW_bt 
		else
			tauW = tauW_at 
		end if

		! After tax marginal benefit of assets
		if (K.lt.theta*a_in) then 
			MB_a = (1.0_dp+R*(1.0_dp-tauK))*(1.0_dp-tauW) 
		else 
			MB_a = (1.0_dp+R*(1.0_dp-tauK))*(1.0_dp-tauW) &
         		& + (P*mu*((theta*z_in)**mu)*a_in**(mu-1.0_DP)-(R+DepRate)*theta)*(1.0_dp-tauK)*(1.0_dp-tauW)
		endif 

	END  FUNCTION MB_a

	FUNCTION MB_a_at(a_in,z_in)
		IMPLICIT NONE   
		real(DP), intent(in)  :: a_in, z_in
		real(DP)			  :: MB_a_at, K

		! Capital demand 
		K   = min( theta*a_in , (mu*P*z_in**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )
		! Compute asset marginal benefit - subject to taxes
		if (K.lt.theta*a_in) then 
			MB_a_at = (1.0_dp+R*(1.0_dp-tauK))*(1.0_dp-tauW_at) 
		else 
			MB_a_at = (1.0_dp+R*(1.0_dp-tauK))*(1.0_dp-tauW_at) &
         		& + (P*mu*((theta*z_in)**mu)*a_in**(mu-1.0_DP)-(R+DepRate)*theta)*(1.0_dp-tauK)*(1.0_dp-tauW_at)
		endif 

	END  FUNCTION MB_a_at

	FUNCTION MB_a_bt(a_in,z_in)
		IMPLICIT NONE   
		real(DP), intent(in) :: a_in, z_in
		real(DP)             :: MB_a_bt, K

		! Capital demand 
		K   = min( theta*a_in , (mu*P*z_in**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )
		! Compute asset marginal benefit - subject to taxes
		if (K.lt.theta*a_in) then 
			MB_a_bt = (1.0_dp+R*(1.0_dp-tauK))*(1.0_dp-tauW_bt) 
		else 
			MB_a_bt = (1.0_dp+R*(1.0_dp-tauK))*(1.0_dp-tauW_bt) &
         		& + (P*mu*((theta*z_in)**mu)*a_in**(mu-1.0_DP)-(R+DepRate)*theta)*(1.0_dp-tauK)*(1.0_dp-tauW_bt)
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
		real(DP)              :: MB_aprime, yprime, cprime, euler_power
		integer               :: age_in, a_in, z_in, l_in, e_in

		! Allocate state and get indeces
		age_in = int(state(1))
		a_in   = int(state(2))
		z_in   = int(state(3))
		l_in   = int(state(4))
		e_in   = int(state(5))

		! Set the power used in the Euler equation for the retirement period
		if (NSU_Switch.eqv..true.) then 
			! Non-Separable Utility
			euler_power = (1.0_dp/((1.0_dp-sigma)*gamma-1.0_dp))
		else 
			! Separable Utility
			euler_power = (-1.0_dp/sigma)
		end if 

		! Compute marginal benefit of next period at a' (given zi)
		MB_aprime = MB_a(aprimet,zgrid(z_in))
		 
		! Compute asset income of next period at a' (given zi)
		yprime   = Y_a(aprimet,zgrid(z_in))
		 
		! Compute consumption of next period given a' (given zi, lambdai and ei)
			! The value of c' comes from interpolating next period's consumption policy function
			! The policy function is defined over a grid of asset income Y, and is interpolated for y'
		cprime =   Linear_Int(Ygrid_t(:,z_in), Cons_t(age_in+1,:,z_in,l_in,e_in),na_t, yprime)    

		! Evaluate square residual of Euler equation at current state (given by (ai,zi,lambdai,ei)) and savings given by a'
		FOC_R	= ( (YGRID_t(a_in,z_in)+RetY_lambda_e(l_in,e_in)-aprimet)    &
		           & - ((beta * survP(age_in) *  MB_aprime)**euler_power)  * cprime  ) ** 2.0_DP

	END  FUNCTION FOC_R

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
		real(DP)             :: ctemp, ntemp, yprime, MB_aprime, FOC_WH, exp1overcprime, E_MU_cp
		REAL(DP)             :: brentvaluet, consin, H_min, c_foc, c_budget, cp, hp
		integer              :: ep_ind
		real(DP), dimension(ne) :: cprime, nprime, MU_cp
		REAL(DP), DIMENSION(6)  :: par_FOC
		integer               :: age_in, a_in, z_in, l_in, e_in

		! Allocate state and get indeces
		age_in = int(state(1))
		a_in   = int(state(2))
		z_in   = int(state(3))
		l_in   = int(state(4))
		e_in   = int(state(5))

		par_FOC(1:5) = state

		H_min = 0.000001 

		! Compute marginal benefit of next period at a' (given zi)
		MB_aprime = MB_a(aprimet,zgrid(z_in))
		! Compute asset income of next period at a' (given zi)
		yprime   = Y_a(aprimet,zgrid(z_in))

		! Set auxiliary variable for FOC_HA
		par_FOC(6) = aprimet
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
			ctemp   = YGRID_t(a_in,z_in)+  Y_h(ntemp,age_in,l_in,e_in,wage) - aprimet   

		if (NSU_Switch.eqv..true.) then
			if (Progressive_Tax_Switch) then 
				! Non-Separable Utility
				if (Log_Switch.eqv..true.) then
					! I have to evaluate the FOC in expectation over eindx prime given eindx
					! Compute c' for each value of e'
					DO ep_ind=1,ne
					    cprime(ep_ind) = Linear_Int(Ygrid_t(:,z_in),Cons_t(age_in+1,:,z_in,l_in,ep_ind), na_t, yprime  )
					ENDDO
					! Compute the expected value of 1/c' conditional on current ei
					E_MU_cp = SUM( pr_e(e_in,:) / cprime )

					! Evaluate the squared residual of the Euler equation for working period
					FOC_WH   = ( (1.0_dp/ctemp)  - beta * survP(age_in) * MB_aprime * E_MU_cp ) **2.0_DP 
				else
					! I have to evaluate the FOC in expectation over eindx prime given eindx
					! Compute consumption and labor for eachvalue of eindx prime
					DO ep_ind=1,ne
					    cp     = Linear_Int(Ygrid_t(:,z_in),Cons_t(age_in+1,:,z_in,l_in,ep_ind), na_t, yprime  )
						c_foc  = (gamma/(1.0_dp-gamma))*(1.0_dp-H_min)*MB_h(H_min,age_in,l_in,e_in,wage)
							if (cp.ge.c_foc) then
								hp = 0.0_dp
							else
								par_FOC(4) = ep_ind
								par_FOC(6) = cp 
								brentvaluet = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_H, brent_tol, hp, par_FOC ) 
							end if 
					    MU_cp(ep_ind) = cp*((1.0_dp-sigma)*gamma-1.0_dp) * (1.0_dp-hp)**((1.0_dp-sigma)*(1.0_dp-gamma))
					END DO
					! Compute the expected value of 1/c' conditional on current ei
					E_MU_cp = SUM( pr_e(e_in,:) * MU_cp )

					! Evaluate the squared residual of the Euler equation for working period
					FOC_WH = ( ctemp**((1.0_dp-sigma)*gamma-1.0_dp) * (1.0_dp-ntemp)**((1.0_dp-sigma)*(1.0_dp-gamma)) & 
						         & - beta*survP(age_in)*MB_aprime*E_MU_cp  )**2.0_DP
				end if 
			else ! Linear Taxes 
				ntemp = max(0.0_DP , gamma - (1.0_DP-gamma)*(YGRID_t(a_in,z_in) - aprimet)/(psi*yh(age_in,l_in,e_in)) )
				ctemp = YGRID_t(a_in,z_in) + psi*yh(age_in,l_in,e_in) * ntemp - aprimet

				DO ep_ind=1,ne
				      cprime(ep_ind) = Linear_Int(Ygrid(:,z_in), Cons_t(age_in+1,:,z_in,l_in,ep_ind), na,    yprime  )
				      nprime(ep_ind) = 1.0_DP-(1.0_DP-gamma)*cprime(ep_ind)/(gamma*psi*yh(age_in,l_in,e_in))
				ENDDO

				nprime = max(0.0_DP,nprime)

				FOC_WH = ((ctemp**(gamma*(1.0_DP-sigma)-1))*((1.0_DP-ntemp)**((1.0_DP-gamma)*(1.0_DP-sigma))) &
						& - beta*survP(age_in)*MB_aprime  &
						& * sum( pr_e(e_in,:) * (cprime**(gamma*(1.0_DP-sigma)-1.0_DP)) &
						& *((1.0_DP-nprime)**((1.0_DP-gamma)*(1.0_DP-sigma)))) )**2.0_DP
			end if 
		else 
			! Separable Utility
				! I have to evaluate the FOC in expectation over eindx prime given eindx
				! Compute c' for each value of e'
				DO ep_ind=1,ne
				    cprime(ep_ind) = Linear_Int(Ygrid_t(:,z_in),Cons_t(age+1,:,z_in,l_in,ep_ind), na_t, yprime  )
				ENDDO
				! Compute the expected value of 1/c' conditional on current ei
				E_MU_cp = SUM( pr_e(e_in,:) / cprime**sigma )

				! Evaluate the squared residual of the Euler equation for working period
				FOC_WH   = (ctemp - 1.0_dp/(beta*survP(age_in)*MB_aprime*E_MU_cp)**(1.0_dp/sigma)) **2.0_DP 
		end if 

	END  FUNCTION FOC_WH

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
		integer                 :: age_in, a_in, z_in, l_in, e_in

		! Allocate state and get indeces
		age_in = int(par(1))
		a_in   = int(par(2))
		z_in   = int(par(3))
		l_in   = int(par(4))
		e_in   = int(par(5))
		cons_in = par(6)

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
		real(DP)             	:: FOC_H_NSU, cons, E_MU_cp
		real(DP), dimension(ne) :: MU_cp
		real(DP)                :: MB_aprime
		integer                 :: age_in, a_in, z_in, l_in, e_in

		! Allocate state and get indeces
		age_in = int(par(1))
		a_in   = int(par(2))
		z_in   = int(par(3))
		l_in   = int(par(4))
		e_in   = int(par(5))
		MB_aprime = par(6)

		! Compute current consumption implied by hours and labor FOC
			cons  = (gamma/(1.0_dp-gamma))*(1.0_dp-hoursin)*MB_h(hoursin,age_in,l_in,e_in,wage)
		! Compute marginal utility of consumption for next period at a' (for all values of e)
	    	MU_cp = Cons_t(age_in+1,a_in,z_in,l_in,:)**((1.0_dp-sigma)*gamma-1.0_dp) &
	    	        * (1.0_dp-Hours_t(age_in+1,a_in,z_in,l_in,:))**((1.0_dp-sigma)*(1.0_dp-gamma))
	    ! Compute expected value of marginal uitility of consumption
	    	E_MU_cp = SUM( pr_e(e_in,:) * MU_cp )
		! Compute square residual of Euler FOC
			FOC_H_NSU = ( cons**((1.0_dp-sigma)*gamma-1.0_dp) * (1.0_dp-hoursin)**((1.0_dp-sigma)*(1.0_dp-gamma)) & 
			         	& - beta*survP(age_in)*MB_aprime*E_MU_cp  )**2.0_DP

	END  FUNCTION FOC_H_NSU

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
		integer              :: age_in, a_in, z_in, l_in, e_in

		! Allocate state and get indeces
		age_in = int(par(1))
		a_in   = int(par(2))
		z_in   = int(par(3))
		l_in   = int(par(4))
		e_in   = int(par(5))
		ap_in  = par(6)

		! Consumption given ain and hoursin
			cons   = YGRID_t(a_in,z_in)+  Y_h(hoursin,age_in,l_in,e_in,wage) - ap_in
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
		real(dp), intent(in)   :: MB_in , H_min, state_FOC(5)
		real(dp), intent(out)  :: C_endo, H_endo, Y_endo
		real(dp)               :: C_euler, C_foc, brentvalue
		REAL(DP), DIMENSION(6) :: par_FOC
		integer                :: age, ai, zi, lambdai, ei

		! Allocate state and get indeces
		age 	= int(state_FOC(1))
		ai   	= int(state_FOC(2))
		zi      = int(state_FOC(3))
		lambdai = int(state_FOC(4))
		ei      = int(state_FOC(5))

		par_FOC(1:5) = state_FOC

		if (Progressive_Tax_Switch.eqv..true.) then !Progressive labor taxes 
			if (NSU_Switch.eqv..true.) then 
				! Non-Separable Utility
				C_euler = ( (beta*survP(age)*MB_in) * SUM(pr_e(ei,:) * &
				        & Cons_t(age+1,ai,zi,lambdai,:)**((1.0_dp-sigma)*gamma-1.0_dp)    * &
				        & (1.0_dp-Hours_t(age+1,ai,zi,lambdai,:))**((1.0_dp-sigma)*(1.0_dp-gamma))  ) ) **(1.0_dp/((1.0_dp-sigma)*gamma-1.0_dp))
				C_foc   = (gamma/(1.0_dp-gamma))*(1.0_dp-H_min)*MB_h(H_min,age,lambdai,ei,wage)

				if (C_euler.ge.C_foc) then
					C_endo  = C_euler 
					H_endo = 0.0_dp
				else 
					if (Log_Switch.eqv..true.) then
			  			! Auxiliary consumption variable for FOC_H        
					    par_FOC(6) = C_euler
					    ! Solution for hours from Euler or FOC  equation according to sigma 
					    brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_H, brent_tol, H_endo , par_FOC ) 
					    C_endo = C_euler
					else 
					    ! Set Marginal benefit of assets to the below threshold level
					    par_FOC(6) = MB_in
					    ! Solution for hours from Euler or FOC  equation according to sigma 
					    brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_H_NSU, brent_tol, H_endo , par_FOC ) 
					    ! Implied consumption by hours from Labor FOC
					    C_endo = (gamma/(1.0_dp-gamma))*(1.0_dp-H_endo)*MB_h(H_endo,age,lambdai,ei,wage)
					end if 
				end if 

			else 
				! Separable Utility
				C_endo = 1.0_dp/( (beta*survP(age)*MB_in) * SUM(pr_e(ei,:) * &
				          & Cons_t(age+1,ai,zi,lambdai,:)**(-sigma) )  ) **(1.0_dp/sigma)
				C_foc  = (MB_h(H_min,age,lambdai,ei,wage)*(1.0_dp-H_min)**(gamma)/phi)**(1.0_dp/sigma)

				if (C_endo.ge.C_foc) then
				  H_endo = 0.0_dp
				else 
				  ! Auxiliary consumption variable for FOC_H        
				  par_FOC(6) = C_endo
				  ! Solution for hours from Euler or FOC  equation according to sigma 
				  brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_H, brent_tol, H_endo , par_FOC ) 
				end if  

				end if 

		else ! Linear labor taxes 
			if (NSU_Switch.eqv..true.) then 
				! Non-Separable Utility
				  C_endo = ((gamma*psi*yh(age, lambdai,ei)/(1.0_DP-gamma))**((1.0_DP-gamma)*(1.0_DP-sigma)) &
				    & *  beta*survP(age)*MB_in  &
				    & *  sum( pr_e(ei,:) * (Cons_t(age+1,ai,zi,lambdai,:)**(gamma*(1.0_DP-sigma)-1.0_DP)) &
				    & *  ( (1.0_DP-Hours_t(age+1, ai, zi, lambdai,:))**((1.0_DP-gamma)*(1.0_DP-sigma)))))**(-1.0_DP/sigma)

				  H_endo = 1.0_DP - (1.0_DP-gamma)*C_endo/(gamma*psi*yh(age, lambdai,ei))   

				  If (H_endo .lt. 0.0_DP) then
				    H_endo = 0.0_DP 
				    C_endo  = ( beta*survP(age)*MB_in  &
				    & * sum( pr_e(ei,:) * (Cons_t(age+1,ai,zi,lambdai,:)**(gamma*(1.0_DP-sigma)-1.0_DP)) &
				    & *((1.0_DP-Hours_t(age+1,ai,zi,lambdai,:))**((1.0_DP-gamma)*(1.0_DP-sigma)))) )**(1.0_DP/(gamma*(1.0_DP-sigma)-1.0_DP))
				  endif 
			else 
				! Separable Utility
				C_endo  = 1.0_DP/( beta*survP(age)*MB_in &
					& * sum( pr_e(ei,:) * (Cons_t(age+1,ai,zi,lambdai,:)**(-sigma)) ) )**(1.0_DP/sigma)

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


	! Age Brackets
		age_limit = [0, 5, 15, 25, 35, 45, 55, MaxAge ]

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
		CALL COMPUTE_VALUE_FUNCTION_LINEAR
		!CALL COMPUTE_VALUE_FUNCTION_SPLINE  
		ValueFunction_Bench = ValueFunction

	! Profit Matrix
		Pr_mat = Profit_Matrix(R,P)

	! Print policy functions and distribution 
		OPEN (UNIT=70, FILE=trim(Result_Folder)//'cons_by_age_z_bench', STATUS='replace')    
		OPEN (UNIT=80, FILE=trim(Result_Folder)//'leisure_by_age_z_bench', STATUS='replace')    
		OPEN (UNIT=90, FILE=trim(Result_Folder)//'dbn_by_age_z_bench', STATUS='replace')   
		DO age=1,MaxAge    
		    DO zi=1,nz
		        temp_cons_by_z(zi)       	 = sum(Cons(age,:,zi,:,:)*DBN_bench(age,:,zi,:,:))/sum(DBN_bench(age,:,zi,:,:))
		        temp_leisure_by_z(zi)    	 = sum((1.0_DP-HOURS(age,:,zi,:,:))*DBN_bench(age,:,zi,:,:))/sum(DBN_bench(age,:,zi,:,:))
		        size_by_age_z_bench(age, zi) = sum(DBN_bench(age,:,zi,:,:))
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
		    WRITE  (UNIT=50, FMT=*) Aprime(1, :, zi, nlambda/2+1, ne/2+1)
		ENDDO
		close (unit=50)

		OPEN (UNIT=50, FILE=trim(Result_Folder)//'aprime_age16_bench', STATUS='replace')   
		DO zi=1,nz 
		    WRITE  (UNIT=50, FMT=*) Aprime(16, :, zi, nlambda/2+1, ne/2+1)
		ENDDO
		close (unit=50)  

		OPEN (UNIT=50, FILE=trim(Result_Folder)//'aprime_age31_bench', STATUS='replace')   
		DO zi=1,nz 
		    WRITE  (UNIT=50, FMT=*) Aprime(31, :, zi, nlambda/2+1, ne/2+1)
		ENDDO
		close (unit=50)  

		OPEN (UNIT=50, FILE=trim(Result_Folder)//'aprime_age46_bench', STATUS='replace')   
		DO zi=1,nz 
		    WRITE  (UNIT=50, FMT=*) Aprime(46, :, zi, nlambda/2+1, ne/2+1)
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
		 
		    DO ai=1,na
	        DO zi=1,nz
            DO lambdai=1,nlambda
            DO ei=1,ne
                size_by_agegroup_z_bench(age_group_counter,zi) = size_by_agegroup_z_bench(age_group_counter,zi) + &
                         & DBN_bench(age,ai,zi,lambdai,ei)       

                tot_wealth_by_agegroup_z_bench(age_group_counter,zi) = tot_wealth_by_agegroup_z_bench(age_group_counter,zi) + &
                         & agrid(ai)*DBN_bench(age,ai,zi,lambdai,ei)                           
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
	 
	    DO ai=1,na
	    DO zi=1,nz
	    DO lambdai=1,nlambda
	    DO ei=1,ne

	        tot_cons_by_agegroup_by_z_bench(age_group_counter,zi) =  tot_cons_by_agegroup_by_z_bench(age_group_counter,zi)+ &
	                         & Cons(age,ai,zi,lambdai,ei)*DBN_bench(age,ai,zi,lambdai,ei)

	        tot_hours_by_agegroup_by_z_bench(age_group_counter,zi) = tot_hours_by_agegroup_by_z_bench(age_group_counter,zi)+&
	                         & HOURS(age,ai,zi,lambdai,ei)*DBN_bench(age,ai,zi,lambdai,ei)

	        tot_aprime_by_agegroup_by_z_bench(age_group_counter,zi)=tot_aprime_by_agegroup_by_z_bench(age_group_counter,zi)+&
	                         & Aprime(age,ai,zi,lambdai,ei)*DBN_bench(age,ai,zi,lambdai,ei)

	        tot_return_by_agegroup_by_z_bench(age_group_counter,zi) = tot_at_return_by_agegroup_by_z_bench(age_group_counter,zi)+&
	                         & DBN_bench(age,ai,zi,lambdai,ei) * (R*agrid(ai) + Pr_mat(ai,zi) )

	        tot_at_return_by_agegroup_by_z_bench(age_group_counter,zi) = tot_at_return_by_agegroup_by_z_bench(age_group_counter,zi)+&
	                         & DBN_bench(age,ai,zi,lambdai,ei) * (Y_a(agrid(ai),zgrid(zi))-1.0_dp)
	       
	        tot_lab_tax_by_agegroup_by_z_bench(age_group_counter,zi) =tot_lab_tax_by_agegroup_by_z_bench(age_group_counter,zi) + &
	                         & DBN_bench(age,ai,zi,lambdai,ei) * ( wage_bench * eff_un(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei) - &
	                          & Y_h(Hours(age,ai,zi,lambdai,ei),age,lambdai,ei,wage_bench) )
	      
	        tot_cons_tax_by_agegroup_by_z_bench(age_group_counter,zi)=tot_cons_tax_by_agegroup_by_z_bench(age_group_counter,zi)+&
	                         & DBN_bench(age,ai,zi,lambdai,ei) * (  tauC *cons(age,ai,zi,lambdai,ei) )
	         
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
		!CALL COMPUTE_VALUE_FUNCTION_LINEAR
		!CALL COMPUTE_VALUE_FUNCTION_SPLINE 
		!ValueFunction_Exp = ValueFunction

	! Profit Matrix
		Pr_mat = Profit_Matrix(R,P)

	! Print policy functions and distribution 
		OPEN (UNIT=70, FILE=trim(Result_Folder)//'cons_by_age_z_exp', STATUS='replace')    
		OPEN (UNIT=80, FILE=trim(Result_Folder)//'leisure_by_age_z_exp', STATUS='replace')   
		OPEN (UNIT=90, FILE=trim(Result_Folder)//'dbn_by_age_z_exp', STATUS='replace')   
		DO age=1,MaxAge 
		    DO zi=1,nz
		          temp_cons_by_z(zi)         = sum(Cons(age,:,zi,:,:)*DBN1(age,:,zi,:,:))/sum(DBN1(age,:,zi,:,:))
		          temp_leisure_by_z(zi)      = sum((1.0_DP-HOURS(age,:,zi,:,:))*DBN1(age,:,zi,:,:))/sum(DBN1(age,:,zi,:,:))
		          size_by_age_z_exp(age, zi) = sum(DBN1(age,:,zi,:,:))
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
		    WRITE  (UNIT=50, FMT=*) Aprime(1, :, zi, nlambda/2+1, ne/2+1)
		ENDDO
		close (unit=50)

		OPEN (UNIT=50, FILE=trim(Result_Folder)//'aprime_age16_exp', STATUS='replace')   
		DO zi=1,nz 
		    WRITE  (UNIT=50, FMT=*) Aprime(16, :, zi, nlambda/2+1, ne/2+1)
		ENDDO
		close (unit=50)  

		OPEN (UNIT=50, FILE=trim(Result_Folder)//'aprime_age31_exp', STATUS='replace')   
		DO zi=1,nz 
		    WRITE  (UNIT=50, FMT=*) Aprime(31, :, zi, nlambda/2+1, ne/2+1)
		ENDDO
		close (unit=50)  

		OPEN (UNIT=50, FILE=trim(Result_Folder)//'aprime_age46_exp', STATUS='replace')   
		DO zi=1,nz 
		    WRITE  (UNIT=50, FMT=*) Aprime(46, :, zi, nlambda/2+1, ne/2+1)
		ENDDO
		close (unit=50)  

	! Consumption Equivalent Welfare
		OPEN (UNIT=50, FILE=trim(Result_Folder)//'CE_alternative', STATUS='replace') 
		    WRITE  (UNIT=50, FMT=*) 'CE Newborn=', ( sum(ValueFunction_exp(1,:,:,:,:)*DBN1(1,:,:,:,:)) / &
		                                &  sum(ValueFunction_Bench(1,:,:,:,:)*DBN_bench(1,:,:,:,:))  ) &
		                                &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP
		    WRITE  (UNIT=50, FMT=*) 'CE =', ( sum(ValueFunction_exp(:,:,:,:,:)*DBN1(:,:,:,:,:)) / &
		                                &  sum(ValueFunction_Bench(:,:,:,:,:)*DBN_bench(:,:,:,:,:))  ) &
		                                &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP

		    WRITE  (UNIT=50, FMT=*) 'Av Utility NB (bench) =',sum(ValueFunction_Bench(1,:,:,:,:)*DBN_bench(1,:,:,:,:))  
		    WRITE  (UNIT=50, FMT=*) 'Av Utility NB (exp  ) =',sum(ValueFunction_exp(1,:,:,:,:)*DBN1(1,:,:,:,:))    
		    WRITE  (UNIT=50, FMT=*) 'Av Utility  (bench)   =',sum(ValueFunction_Bench(:,:,:,:,:)*DBN_bench(:,:,:,:,:))  
		    WRITE  (UNIT=50, FMT=*) 'Av Utility  (exp)     =',sum(ValueFunction_exp(:,:,:,:,:)*DBN1(:,:,:,:,:))
		    print*,' '
		    print*, 'Av Utility NB (bench) =',sum(ValueFunction_Bench(1,:,:,:,:)*DBN_bench(1,:,:,:,:))  
		    print*, 'Av Utility NB (exp  ) =',sum(ValueFunction_exp(1,:,:,:,:)*DBN1(1,:,:,:,:))    
		close (unit=50)		


		OPEN (UNIT=50, FILE=trim(Result_Folder)//'CE_NEWBORN', STATUS='replace')  
		OPEN (UNIT=60, FILE=trim(Result_Folder)//'CE', STATUS='replace')  
		OPEN (UNIT=70, FILE=trim(Result_Folder)//'CE_by_age', STATUS='replace')  
		OPEN (UNIT=80, FILE=trim(Result_Folder)//'CE_by_age_z', STATUS='replace')  

		DO age=1,MaxAge
			if (Log_Switch.eqv..true.) then 
		    	Cons_Eq_Welfare(age,:,:,:,:)=exp((ValueFunction_exp(age,:,:,:,:)-ValueFunction_Bench(age,:,:,:,:))/CumDiscountF(age))-1.0_DP
		    else 
		    	Cons_Eq_Welfare(age,:,:,:,:)=(ValueFunction_exp(age,:,:,:,:)/ValueFunction_Bench(age,:,:,:,:)) &
                                				&  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP
		    end if 

		    WRITE  (UNIT=70, FMT=*) 100*sum(Cons_Eq_Welfare(age,:,:,:,:)*DBN_bench(age,:,:,:,:))/sum(DBN_bench(age,:,:,:,:))
		    DO zi=1,nz
		         temp_ce_by_z(zi) = 100*sum(Cons_Eq_Welfare(age,:,zi,:,:)*DBN_bench(age,:,zi,:,:))/sum(DBN_bench(age,:,zi,:,:))
		    ENDDO
		    WRITE  (UNIT=80, FMT=*) temp_ce_by_z
		    !print*,'age=',age, temp_ce_by_z, ', mean:  ', &
		    !    & 100*sum(Cons_Eq_Welfare(age,:,:,:,:)*DBN_bench(age,:,:,:,:))/sum(DBN_bench(age,:,:,:,:))
		ENDDO

		CE_NEWBORN = 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN_bench(1,:,:,:,:))/sum(DBN_bench(1,:,:,:,:))

		WRITE  (UNIT=50, FMT=*) 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN_bench(1,:,:,:,:))/sum(DBN_bench(1,:,:,:,:))
		WRITE  (UNIT=50, FMT=*) 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))
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
		            & 100*sum(Cons_Eq_Welfare(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:)* &
		            &                         DBN_bench(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:))/&
		            &                sum( DBN_bench(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:))
		    ENDDO
		ENDDO
		DO age_group_counter=1,max_age_category
		    WRITE  (UNIT=80, FMT=*)  CE_by_agegroup_z(age_group_counter,:)
		ENDDO
		close (unit=80)

		OPEN (UNIT=80, FILE=trim(Result_Folder)//'CE_by_Asset_z_med_E_Lambda_age1', STATUS='replace') 
		OPEN (UNIT=90, FILE=trim(Result_Folder)//'MeanAsset_by_z_med_E_Lambda_age1', STATUS='replace') 
		DO zi=1,nz
		    WRITE  (UNIT=80, FMT=*)   Cons_Eq_Welfare(1,:,zi,nlambda/2+1, ne/2+1)
		    WRITE  (UNIT=90, FMT=*)   sum(agrid*DBN_bench(1,:,zi,nlambda/2+1, ne/2+1)) &
		                                                        & /sum(DBN_bench(1,:,zi,nlambda/2+1, ne/2+1))
		ENDDO           
		close (unit=80)
		close (unit=90)

		OPEN (UNIT=80, FILE=trim(Result_Folder)//'CE_by_Asset_z_med_E_Lambda_age16', STATUS='replace') 
		OPEN (UNIT=90, FILE=trim(Result_Folder)//'MeanAsset_by_z_med_E_Lambda_age16', STATUS='replace') 
		DO zi=1,nz
		    WRITE  (UNIT=80, FMT=*)   Cons_Eq_Welfare(16,:,zi,nlambda/2+1, ne/2+1)
		    WRITE  (UNIT=90, FMT=*)   sum(agrid*DBN_bench(16,:,zi,nlambda/2+1, ne/2+1)) &
		                                                        & /sum(DBN_bench(16,:,zi,nlambda/2+1, ne/2+1))
		ENDDO           
		close (unit=80)
		close (unit=90)

		OPEN (UNIT=80, FILE=trim(Result_Folder)//'CE_by_Asset_z_med_E_Lambda_age31', STATUS='replace') 
		OPEN (UNIT=90, FILE=trim(Result_Folder)//'MeanAsset_by_z_med_E_Lambda_age31', STATUS='replace') 
		DO zi=1,nz
			WRITE  (UNIT=80, FMT=*)   Cons_Eq_Welfare(31,:,zi,nlambda/2+1, ne/2+1)
			WRITE  (UNIT=90, FMT=*)   sum(agrid*DBN_bench(31,:,zi,nlambda/2+1, ne/2+1)) &
			                                            & /sum(DBN_bench(31,:,zi,nlambda/2+1, ne/2+1))
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

		    DO while (age .gt.   age_limit(age_group_counter+1) )
		        age_group_counter=age_group_counter+1
		    ENDDO    
		 
		    DO ai=1,na
		    DO zi=1,nz
            DO lambdai=1,nlambda
            DO ei=1,ne
		                    
                If ( Cons_Eq_Welfare(age,ai,zi,lambdai,ei) .ge. 0.0_DP) then
                    frac_pos_welfare = frac_pos_welfare +DBN_bench(age,ai,zi,lambdai,ei)
                    size_pos_welfare_by_age_z(age,zi) = size_pos_welfare_by_age_z(age,zi) + DBN_bench(age,ai,zi,lambdai,ei)
                    size_pos_welfare_by_agegroup_z(age_group_counter,zi)=size_pos_welfare_by_agegroup_z(age_group_counter,zi) &
                         &+  DBN_bench(age,ai,zi,lambdai,ei)       
                endif 
                size_by_agegroup_z_exp(age_group_counter,zi) = size_by_agegroup_z_exp(age_group_counter,zi) + &
                         & DBN1(age,ai,zi,lambdai,ei)       

                tot_wealth_by_agegroup_z_exp(age_group_counter,zi) = tot_wealth_by_agegroup_z_exp(age_group_counter,zi) + &
                         & agrid(ai)*DBN1(age,ai,zi,lambdai,ei)                           
		                    
	        ENDDO
	        ENDDO
	        ENDDO
		    ENDDO
		ENDDO


	OPEN (UNIT=60, FILE=trim(Result_Folder)//'mean_wealth_by_agegroup_z_exp', STATUS='replace')  
	OPEN (UNIT=70, FILE=trim(Result_Folder)//'size_by_agegroup_z_exp', STATUS='replace')  
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
	 
	    DO ai=1,na
	    DO zi=1,nz
	    DO lambdai=1,nlambda
	    DO ei=1,ne

	        tot_cons_by_agegroup_by_z_exp(age_group_counter,zi) =  tot_cons_by_agegroup_by_z_exp(age_group_counter,zi)+ &
	                         & Cons(age,ai,zi,lambdai,ei)*DBN1(age,ai,zi,lambdai,ei)

	        tot_hours_by_agegroup_by_z_exp(age_group_counter,zi) = tot_hours_by_agegroup_by_z_exp(age_group_counter,zi)+&
	                         & HOURS(age,ai,zi,lambdai,ei)*DBN1(age,ai,zi,lambdai,ei)

	        tot_aprime_by_agegroup_by_z_exp(age_group_counter,zi)=tot_aprime_by_agegroup_by_z_exp(age_group_counter,zi)+&
	                         & Aprime(age,ai,zi,lambdai,ei)*DBN1(age,ai,zi,lambdai,ei)

	        tot_return_by_agegroup_by_z_exp(age_group_counter,zi) = tot_at_return_by_agegroup_by_z_exp(age_group_counter,zi)+&
	                         & DBN1(age,ai,zi,lambdai,ei) * (R*agrid(ai) + Pr_mat(ai,zi) )

	        tot_at_return_by_agegroup_by_z_exp(age_group_counter,zi) = tot_at_return_by_agegroup_by_z_exp(age_group_counter,zi)+&
	                         & DBN1(age,ai,zi,lambdai,ei) * (Y_a(agrid(ai),zgrid(zi))-1.0_dp)
	       
	        tot_lab_tax_by_agegroup_by_z_exp(age_group_counter,zi) =tot_lab_tax_by_agegroup_by_z_exp(age_group_counter,zi) + &
	                         & DBN1(age,ai,zi,lambdai,ei) * ( wage_exp * eff_un(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei) - &
	                          & Y_h(Hours(age,ai,zi,lambdai,ei),age,lambdai,ei,wage_exp) )
	      
	        tot_cons_tax_by_agegroup_by_z_exp(age_group_counter,zi)=tot_cons_tax_by_agegroup_by_z_exp(age_group_counter,zi)+&
	                         & DBN1(age,ai,zi,lambdai,ei) * (  tauC *cons(age,ai,zi,lambdai,ei) )
	         
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
		Welfare_Gain_NB_bench  = 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN_bench(1,:,:,:,:))/sum(DBN_bench(1,:,:,:,:))
		Welfare_Gain_NB_exp    = 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))

	print*,'---------------------------'
	print*,''
	print*,'Average Welfare Gain Whole Population (bench dbn) (prct)=',100.0_DP*sum(Cons_Eq_Welfare*DBN_bench)
	print*,'Average Welfare Gain Whole Population (exp dbn)     (prct)=',100.0_DP*sum(Cons_Eq_Welfare*DBN1)
	print*,'Average Welfare Gain NEW BORN (bench dbn) (prct)          =',&
	    & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN_bench(1,:,:,:,:))/sum(DBN_bench(1,:,:,:,:))
	print*,'Average Welfare Gain NEW BORN (exp dbn)     (prct)          =',&
	    & 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))
	print*,''
	print*,'---------------------------'

	! Compute Average Utility

	Av_Util_NB  =  100.0_dp * (( sum(ValueFunction_exp(1,:,:,:,:)*DBN1(1,:,:,:,:)) / &
				&               sum(ValueFunction_Bench(1,:,:,:,:)*DBN_bench(1,:,:,:,:)) ) &
				& ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)

	Av_Util_Pop =  100.0_dp * (( sum(ValueFunction_exp(:,:,:,:,:)*DBN1(:,:,:,:,:)) / &
				& 	            sum(ValueFunction_Bench(:,:,:,:,:)*DBN_bench(:,:,:,:,:)) ) &
				& ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP)



END SUBROUTINE  COMPUTE_WELFARE_GAIN


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE COMPUTE_VALUE_FUNCTION_SPLINE()
	IMPLICIT NONE
	INTEGER :: tklo, tkhi
	REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) :: PrAprimelo, PrAprimehi
	REAL(DP), DIMENSION(na) :: ValueP1, ValueP2, ValueP, ExpValueP

	! Announce method of interpolation
	! Interpolation is used to get value function at optimal a' (when computing expectations)
	!print*,'VALUE FUNCTION SPLINE'

	! Compute value function at grid nodes
	! Final age
	age=MaxAge
	DO ai=1,na    
    DO zi=1,nz
    DO lambdai=1,nlambda          
	DO ei=1,ne
      	ValueFunction(age, ai, zi, lambdai, ei) = Utility(Cons(age, ai, zi, lambdai, ei),Hours(age, ai, zi, lambdai, ei))
		! print*,Cons(age, ai, zi, lambdai, ei),  ValueFunction(age, ai, zi, lambdai, ei) 
		! pause
	ENDDO ! ei          
    ENDDO ! lambdai
    ENDDO ! zi
	ENDDO ! ai

	! Retirement Period
	DO age=MaxAge-1,RetAge,-1
	DO zi=1,nz
	DO lambdai=1,nlambda          
	DO ei=1,ne            
  		if (NSU_Switch.eqv..true.) then
            CALL spline( agrid, ValueFunction(age+1, :, zi, lambdai, ei) , na , &
      			& gamma*MBGRID(1,zi) *Cons(age+1, 1, zi, lambdai, ei) **((1.0_DP-sigma)*gamma-1.0_DP)/(1_DP+tauC), &
        		& gamma*MBGRID(na,zi)*Cons(age+1, na, zi, lambdai, ei)**((1.0_DP-sigma)*gamma-1.0_DP)/(1_DP+tauC), ValueP2)  
        else 
        	CALL spline( agrid, ValueFunction(age+1, :, zi, lambdai, ei) , na , &
      			& gamma*MBGRID(1,zi) /Cons(age+1, 1, zi, lambdai, ei) **(sigma)/(1_DP+tauC), &
        		& gamma*MBGRID(na,zi)/Cons(age+1, na, zi, lambdai, ei)**(sigma)/(1_DP+tauC), ValueP2)  
      	end if 
	                  
        DO ai=1,na    
            call splint( agrid, ValueFunction(age+1, :, zi, lambdai, ei), &
                    & ValueP2, na, Aprime(age,ai,zi,lambdai, ei), ValueP(ai))  

            ValueFunction(age,ai,zi,lambdai,ei) = Utility(Cons(age,ai,zi,lambdai,ei),Hours(age,ai,zi,lambdai,ei)) &
              										& + beta*survP(age)* ValueP(ai)      
        ENDDO ! ai
	              
    ENDDO ! ei          
	ENDDO ! lambdai
    ENDDO ! zi
	ENDDO ! age
   

	! Working Period
	DO age=RetAge-1,1,-1
    DO zi=1,nz
    DO lambdai=1,nlambda          
	DO ei=1,ne
        DO ai=1,na    
              ExpValueP(ai) = sum(ValueFunction(age+1, ai, zi, lambdai, :) * pr_e(ei,:))
        ENDDO

        if (NSU_Switch.eqv..true.) then
        	CALL spline( agrid, ExpValueP , na , &
		            & (gamma*MBGRID(1,zi)/(1.0_DP+tauC)) * sum(pr_e(ei,:)* &
		            & Cons(age+1, 1, zi, lambdai, :)**((1.0_DP-sigma)*gamma-1.0_DP) * &
		            & (1.0_DP-Hours(age+1,1,zi,lambdai,:))**((1.0_DP-gamma)*(1.0_DP-sigma))),& 
		            & (gamma*MBGRID(na,zi)/(1.0_DP+tauC)) * sum(pr_e(ei,:)* &
		            & Cons(age+1, na, zi, lambdai, :)**((1.0_DP-sigma)*gamma-1.0_DP) * &
		            & (1.0_DP-Hours(age+1,na,zi,lambdai,:))**((1.0_DP-gamma)*(1.0_DP-sigma))),&
		            & ValueP2)   
        else 
	        CALL spline( agrid, ExpValueP , na , &
		            & (gamma*MBGRID(1,zi)/(1.0_DP+tauC))  * sum(pr_e(ei,:) / Cons(age+1,  1, zi, lambdai, :)**(sigma) ) ,&
		            & (gamma*MBGRID(na,zi)/(1.0_DP+tauC)) * sum(pr_e(ei,:) / Cons(age+1, na, zi, lambdai, :)**(sigma) ) ,&
		            & ValueP2)  
        end if 

        DO ai=1,na 
        	call splint( agrid, ExpValueP, ValueP2, na, Aprime(age,ai,zi,lambdai, ei), ValueP(ai))   
        
        	ValueFunction(age, ai, zi, lambdai, ei) = Utility(Cons(age,ai,zi,lambdai,ei),Hours(age,ai,zi,lambdai,ei)) &
        	                                         & + beta*survP(age)*ValueP(ai)
		ENDDO ! ai
	ENDDO ! ei          
    ENDDO ! lambdai
    ENDDO ! zi
	ENDDO ! age

END SUBROUTINE COMPUTE_VALUE_FUNCTION_SPLINE 


!========================================================================================
!========================================================================================
!========================================================================================


! SUBROUTINE COMPUTE_VALUE_FUNCTION_LINEAR()
! 	IMPLICIT NONE
! 	INTEGER :: tklo, tkhi
! 	REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) :: PrAprimelo, PrAprimehi

! 	print*,'VALUE FUNCTION LINEAR'

! 	age=MaxAge
! 	DO ai=1,na    
!     DO zi=1,nz
!     DO lambdai=1,nlambda          
! 	DO ei=1,ne
!       	! ValueFunction(age, ai, zi, lambdai, ei) = Utility(Cons(age, ai, zi, lambdai, ei),Hours(age, ai, zi, lambdai, ei))
!       	ValueFunction(age, ai, zi, lambdai, ei) = ((Cons(age,ai,zi,lambdai,ei)**gamma) &
!                    & * (1.0_DP-Hours(age,ai,zi,lambdai,ei))**(1.0_DP-gamma))**(1.0_DP-sigma)/(1.0_DP-sigma)  
! 		! print*,Cons(age, ai, zi, lambdai, ei),  ValueFunction(age, ai, zi, lambdai, ei) 
! 		! pause
! 	ENDDO ! ei          
!     ENDDO ! lambdai
!     ENDDO ! zi
! 	ENDDO ! ai

! 	! Retirement Period
! 	DO age=MaxAge-1,RetAge,-1
! 	DO ai=1,na    
!     DO zi=1,nz
!     DO lambdai=1,nlambda          
! 	DO ei=1,ne
! 		if ( Aprime(age,ai,zi,lambdai, ei) .ge. amax) then
! 		    tklo =na-1
! 		else if (Aprime(age,ai,zi,lambdai, ei) .lt. amin) then
! 		    tklo = 1
! 		else
! 		    tklo = ((Aprime(age,ai,zi,lambdai, ei) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
! 		end if            
! 		tkhi = tklo + 1        
! 		PrAprimelo(age,ai,zi,lambdai, ei) = ( agrid(tkhi) - Aprime(age,ai,zi,lambdai, ei) ) / ( agrid(tkhi) -agrid(tklo) )
! 		PrAprimehi(age,ai,zi,lambdai, ei) = ( Aprime(age,ai,zi,lambdai, ei) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )        
! 		PrAprimelo(age,ai,zi,lambdai, ei) = min (PrAprimelo(age,ai,zi,lambdai, ei), 1.0_DP)
! 		PrAprimelo(age,ai,zi,lambdai, ei) = max(PrAprimelo(age,ai,zi,lambdai, ei), 0.0_DP)
! 		PrAprimehi(age,ai,zi,lambdai, ei) = min (PrAprimehi(age,ai,zi,lambdai, ei), 1.0_DP)
! 		PrAprimehi(age,ai,zi,lambdai, ei) = max(PrAprimehi(age,ai,zi,lambdai, ei), 0.0_DP)    

! 		! ValueFunction(age, ai, zi, lambdai, ei) = Utility(Cons(age,ai,zi,lambdai,ei),Hours(age,ai,zi,lambdai,ei)) &
! 		! 	  & + beta*survP(age)* (PrAprimelo(age,ai,zi,lambdai, ei)*ValueFunction(age+1, tklo, zi, lambdai, ei)&
! 		! 	  & 				 +  PrAprimehi(age,ai,zi,lambdai, ei)*ValueFunction(age+1, tkhi, zi, lambdai, ei))

!         ValueFunction(age, ai, zi, lambdai, ei) = ((Cons(age,ai,zi,lambdai,ei)**gamma) &
!                       & * (1.0_DP-Hours(age,ai,zi,lambdai,ei))**(1.0_DP-gamma))**(1.0_DP-sigma)/(1.0_DP-sigma) &
!                       & + beta*survP(age)* (PrAprimelo(age,ai,zi,lambdai, ei)*ValueFunction(age+1, tklo, zi, lambdai, ei)&
!                       & +                   PrAprimehi(age,ai,zi,lambdai, ei)*ValueFunction(age+1, tkhi, zi, lambdai, ei))
! 	ENDDO ! ei          
!     ENDDO ! lambdai
!     ENDDO ! zi
! 	ENDDO ! ai
! 	ENDDO ! age
! 	!print*,ValueFunction


! 	! Working Period
! 	DO age=RetAge-1,1,-1
! 	DO ai=1,na    
!     DO zi=1,nz
!     DO lambdai=1,nlambda          
! 	DO ei=1,ne
! 		if ( Aprime(age,ai,zi,lambdai, ei) .ge. amax) then
! 		    tklo =na-1
! 	    elseif (Aprime(age,ai,zi,lambdai, ei) .lt. amin) then
! 	        tklo = 1
!         else
!             tklo = ((Aprime(age,ai,zi,lambdai, ei) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
! 		endif  

! 		tkhi = tklo + 1        
! 		PrAprimelo(age,ai,zi,lambdai, ei) = ( agrid(tkhi) - Aprime(age,ai,zi,lambdai, ei) ) / ( agrid(tkhi) -agrid(tklo) )
! 		PrAprimehi(age,ai,zi,lambdai, ei) = ( Aprime(age,ai,zi,lambdai, ei) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )        
! 		PrAprimelo(age,ai,zi,lambdai, ei) = min (PrAprimelo(age,ai,zi,lambdai, ei), 1.0_DP)
! 		PrAprimelo(age,ai,zi,lambdai, ei) = max(PrAprimelo(age,ai,zi,lambdai, ei), 0.0_DP)
! 		PrAprimehi(age,ai,zi,lambdai, ei) = min (PrAprimehi(age,ai,zi,lambdai, ei), 1.0_DP)
! 		PrAprimehi(age,ai,zi,lambdai, ei) = max(PrAprimehi(age,ai,zi,lambdai, ei), 0.0_DP)    

! 		! ValueFunction(age, ai, zi, lambdai, ei) = Utility(Cons(age,ai,zi,lambdai,ei),Hours(age,ai,zi,lambdai,ei))  &
! 		!    & + beta*survP(age)* sum( ( PrAprimelo(age,ai,zi,lambdai, ei) * ValueFunction(age+1, tklo, zi, lambdai,:)  &
! 		!    & 						 + PrAprimehi(age,ai,zi,lambdai, ei) * ValueFunction(age+1, tkhi, zi, lambdai,:)) * pr_e(ei,:))
! 		ValueFunction(age, ai, zi, lambdai, ei) = ((Cons(age,ai,zi,lambdai,ei)**gamma) &
!                        & * (1.0_DP-Hours(age,ai,zi,lambdai,ei))**(1.0_DP-gamma))**(1.0_DP-sigma)/(1.0_DP-sigma) &
!                        & + beta*survP(age)* sum( ( PrAprimelo(age,ai,zi,lambdai, ei) * ValueFunction(age+1, tklo, zi, lambdai,:)  &
!                        & + PrAprimehi(age,ai,zi,lambdai, ei) * ValueFunction(age+1, tkhi, zi, lambdai,:)) * pr_e(ei,:) )
! 		! if ( ValueFunction(age, ai, zi, lambdai, ei) .lt. (-100.0_DP) ) then
! 		!    print*,'ValueFunction(age, ai, zi, lambdai, ei)=',ValueFunction(age, ai, zi, lambdai, ei)
! 		! endif
! 	ENDDO ! ei          
!     ENDDO ! lambdai
!     ENDDO ! zi
! 	ENDDO ! ai
! 	ENDDO ! age

! END SUBROUTINE COMPUTE_VALUE_FUNCTION_LINEAR 

SUBROUTINE COMPUTE_VALUE_FUNCTION_LINEAR
IMPLICIT NONE
INTEGER :: tklo, tkhi
REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) :: PrAprimelo, PrAprimehi

age=MaxAge
DO ai=1,na    
DO zi=1,nz
DO lambdai=1,nlambda          
DO ei=1,ne
	ValueFunction(age, ai, zi, lambdai, ei) = ((Cons(age,ai,zi,lambdai,ei)**gamma) &
	& * (1.0_DP-Hours(age,ai,zi,lambdai,ei))**(1.0_DP-gamma))**(1.0_DP-sigma)/(1.0_DP-sigma)  
	!                  print*,Cons(age, ai, zi, lambdai, ei),  ValueFunction(age, ai, zi, lambdai, ei) 
	!                  pause
ENDDO ! ei          
ENDDO ! lambdai
ENDDO ! zi
ENDDO ! ai

! Retirement Period
DO age=MaxAge-1,RetAge,-1
DO ai=1,na    
DO zi=1,nz
DO lambdai=1,nlambda          
DO ei=1,ne
	if ( Aprime(age,ai,zi,lambdai, ei) .ge. amax) then
		tklo =na-1
	elseif (Aprime(age,ai,zi,lambdai, ei) .lt. amin) then
	    tklo = 1
    else
        tklo = ((Aprime(age,ai,zi,lambdai, ei) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
	endif            
	tkhi = tklo + 1        
	PrAprimelo(age,ai,zi,lambdai, ei) = ( agrid(tkhi) - Aprime(age,ai,zi,lambdai, ei) ) / ( agrid(tkhi) -agrid(tklo) )
	PrAprimehi(age,ai,zi,lambdai, ei) = ( Aprime(age,ai,zi,lambdai, ei) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )        
	PrAprimelo(age,ai,zi,lambdai, ei) = min (PrAprimelo(age,ai,zi,lambdai, ei), 1.0_DP)
	PrAprimelo(age,ai,zi,lambdai, ei) = max(PrAprimelo(age,ai,zi,lambdai, ei), 0.0_DP)
	PrAprimehi(age,ai,zi,lambdai, ei) = min (PrAprimehi(age,ai,zi,lambdai, ei), 1.0_DP)
	PrAprimehi(age,ai,zi,lambdai, ei) = max(PrAprimehi(age,ai,zi,lambdai, ei), 0.0_DP)    

	ValueFunction(age, ai, zi, lambdai, ei) = ((Cons(age,ai,zi,lambdai,ei)**gamma) &
	& * (1.0_DP-Hours(age,ai,zi,lambdai,ei))**(1.0_DP-gamma))**(1.0_DP-sigma)/(1.0_DP-sigma) &
	& + beta*survP(age)* (PrAprimelo(age,ai,zi,lambdai, ei)*ValueFunction(age+1, tklo, zi, lambdai, ei)&
	& +                   PrAprimehi(age,ai,zi,lambdai, ei)*ValueFunction(age+1, tkhi, zi, lambdai, ei))
ENDDO ! ei          
ENDDO ! lambdai
ENDDO ! zi
ENDDO ! ai
ENDDO ! age
!print*,ValueFunction


! Working Period
DO age=RetAge-1,1,-1
DO ai=1,na    
DO zi=1,nz
DO lambdai=1,nlambda          
DO ei=1,ne
	if ( Aprime(age,ai,zi,lambdai, ei) .ge. amax) then
		tklo =na-1
    elseif (Aprime(age,ai,zi,lambdai, ei) .lt. amin) then
		tklo = 1
	else
		tklo = ((Aprime(age,ai,zi,lambdai, ei) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
	endif            
	tkhi = tklo + 1        
	PrAprimelo(age,ai,zi,lambdai, ei) = ( agrid(tkhi) - Aprime(age,ai,zi,lambdai, ei) ) / ( agrid(tkhi) -agrid(tklo) )
	PrAprimehi(age,ai,zi,lambdai, ei) = ( Aprime(age,ai,zi,lambdai, ei) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )        
	PrAprimelo(age,ai,zi,lambdai, ei) = min (PrAprimelo(age,ai,zi,lambdai, ei), 1.0_DP)
	PrAprimelo(age,ai,zi,lambdai, ei) = max(PrAprimelo(age,ai,zi,lambdai, ei), 0.0_DP)
	PrAprimehi(age,ai,zi,lambdai, ei) = min (PrAprimehi(age,ai,zi,lambdai, ei), 1.0_DP)
	PrAprimehi(age,ai,zi,lambdai, ei) = max(PrAprimehi(age,ai,zi,lambdai, ei), 0.0_DP)    

	ValueFunction(age, ai, zi, lambdai, ei) = ((Cons(age,ai,zi,lambdai,ei)**gamma) &
	       & * (1.0_DP-Hours(age,ai,zi,lambdai,ei))**(1.0_DP-gamma))**(1.0_DP-sigma)/(1.0_DP-sigma) &
	       & + beta*survP(age)* sum( ( PrAprimelo(age,ai,zi,lambdai, ei) * ValueFunction(age+1, tklo, zi, lambdai,:)  &
	       & + PrAprimehi(age,ai,zi,lambdai, ei) * ValueFunction(age+1, tkhi, zi, lambdai,:)) * pr_e(ei,:) )
	!                  if ( ValueFunction(age, ai, zi, lambdai, ei) .lt. (-100.0_DP) ) then
	!                       print*,'ValueFunction(age, ai, zi, lambdai, ei)=',ValueFunction(age, ai, zi, lambdai, ei)
	!                  endif
ENDDO ! ei          
ENDDO ! lambdai
ENDDO ! zi
ENDDO ! ai
ENDDO ! age


END SUBROUTINE COMPUTE_VALUE_FUNCTION_LINEAR 



!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE GOVNT_BUDGET()
	IMPLICIT NONE
	real(DP) ::  GBAR_K,  GBAR_W,  GBAR_L, GBAR_C

	! Initialize variables
	GBAR 		 = 0.0_DP
	GBAR_K 		 = 0.0_DP
	GBAR_W		 = 0.0_DP
	GBAR_L 		 = 0.0_DP
	GBAR_C       = 0.0_DP
	SSC_Payments = 0.0_DP

	! Compute total expenditure = total revenue
		! For each state accumulate tax income weighted by equilibrium distribution
		! Taxes come from capital income, wealth, labor income and consumption
		! G_L is the tax revenue that comes exclusively from labor income
	DO age=1, MaxAge
	DO ai=1,na
	DO zi=1,nz
	DO lambdai=1,nlambda
	DO ei=1,ne
	    GBAR = GBAR + DBN1(age,ai,zi,lambdai,ei) * ( tauK*( R*agrid(ai) + Pr_mat(ai,zi) )  	     &
	          & + ( agrid(ai) + ( R*agrid(ai) + Pr_mat(ai,zi) ) *(1.0_DP-tauK)  ) - YGRID(ai,zi) &	
	          & + yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei)  												&
	          & - psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei))**(1.0_DP-tauPL)  					&
	          & + tauC * cons(age, ai, zi, lambdai,ei)  )   

	    GBAR_L = GBAR_L  + DBN1(age,ai,zi,lambdai,ei) * (  yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei) &
	               &- psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei))**(1.0_DP-tauPL) )

	    GBAR_K = GBAR_K + DBN1(age,ai,zi,lambdai,ei) * tauK*( R*agrid(ai) + Pr_mat(ai,zi) )

	    GBAR_W = GBAR_W + DBN1(age,ai,zi,lambdai,ei) * &
	            & (( agrid(ai) + ( R*agrid(ai) + Pr_mat(ai,zi) ) *(1.0_DP-tauK)  ) - YGRID(ai,zi) )

      	GBAR_C = GBAR_C +  DBN1(age,ai,zi,lambdai,ei) * tauC * cons(age, ai, zi, lambdai,ei)
	    
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	! Compute social security expenditure
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

	! Set government expenditures as total expenditures minus social security payments
	GBAR = GBAR -  SSC_Payments


	! Print Results
	print*, ' '
	print*, "Government Budget - Revenues and taxes"
	print*,'GBAR=',GBAR,'SSC_Payments=', SSC_Payments, 'GBAR_L=',GBAR_L,'Av. Labor Tax=', GBAR_L/Ebar 
	print*, 'GBAR_K=', GBAR_K, "GBAR_W=", GBAR_W, 'GBAR_C=', GBAR_C
	print*, 'Tau_K=', tauK, 'Tau_W=', tauW_at, 'Tau_C=', tauC, "Threshold", Y_a_threshold
	print*, ' '
END  SUBROUTINE GOVNT_BUDGET

!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE FIND_DBN_EQ()
	use omp_lib
	IMPLICIT NONE
	INTEGER:: tklo, tkhi, age1, age2, z1, z2, a1, a2, lambda1, lambda2, e1, e2, DBN_iter, simutime, iter_indx
	REAL   :: DBN_dist, DBN_criteria
	real(dp)   ::BBAR, MeanWealth, brent_value
	REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) :: PrAprimelo, PrAprimehi, DBN2
	INTEGER,  DIMENSION(MaxAge, na, nz, nlambda, ne) :: Aplo, Aphi

	!$ call omp_set_num_threads(7)
	DBN_criteria = 1.0E-08_DP

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
	!$omp parallel do private(lambdai,ei,ai,tklo,tkhi)
	DO zi=1,nz
	DO ai=1,na
	DO lambdai=1,nlambda
	DO ei=1, ne
	        if ( Aprime(age,ai,zi,lambdai, ei) .ge. amax) then
	            tklo =na-1
	            elseif (Aprime(age,ai,zi,lambdai, ei) .lt. amin) then
	                 tklo = 1
	                else
	                    tklo = ((Aprime(age,ai,zi,lambdai, ei) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
	        endif            
	        tkhi = tklo + 1        
	        Aplo(age,ai,zi,lambdai, ei)  = tklo
	        Aphi(age,ai,zi,lambdai, ei)  = tkhi        
	        PrAprimelo(age,ai,zi,lambdai, ei) = ( agrid(tkhi) - Aprime(age,ai,zi,lambdai, ei) ) / ( agrid(tkhi) -agrid(tklo) )
	        PrAprimehi(age,ai,zi,lambdai, ei) = ( Aprime(age,ai,zi,lambdai, ei) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )        
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO

	! Probablities are adjusted to lie in [0,1]
		PrAprimelo = min (PrAprimelo, 1.0_DP)
		PrAprimelo = max(PrAprimelo, 0.0_DP)
		PrAprimehi = min (PrAprimehi, 1.0_DP)
		PrAprimehi = max(PrAprimehi, 0.0_DP)


	! The following reads the distribution from the file
		!OPEN   (UNIT=2, FILE='dbn1')
		!DO age=1,MaxAge
		!DO ai=1,na
		!DO zi=1,nz
		!DO lambdai=1,nlambda
		!DO ei=1,ne
		!     READ(unit=2, FMT=*) DBN1(age,ai,zi,lambdai,ei)       
		!ENDDO
		!ENDDO
		!ENDDO
		!ENDDO
		!ENDDO
		!CLOSE (unit=2)


	! Compute distribution of assets by age and state
		! Distribution is obtained by iterating over an initial distribution using policy functions
	DBN_dist=1.0_DP
	simutime = 1
	iter_indx = 1
	!print*, 'Computing Equilibrium Distribution'
	DO WHILE ( ( DBN_dist .ge. DBN_criteria ) .and. ( simutime .le. MaxSimuTime ) )
		!print*, 'Eq. Distribution difference=', DBN_dist
		!    print*, 'sum DBN1=', sum(DBN1)
	    DBN2=0.0_DP

		! Everyone in MaxAge dies. Those who die, switch to z2, lambda2 and start at ne/2+1
	    age1=MaxAge
	    DO z1=1,nz
	    DO a1=1,na
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN2(1,Aplo(age1, a1, z1, lambda1, e1), z2,lambda2,ne/2+1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1) + DBN1(age1, a1, z1, lambda1, e1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2)  *  PrAprimelo(age1, a1, z1, lambda1, e1)
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1) + DBN1(age1, a1, z1, lambda1, e1) & 
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2)*pr_lambda(lambda1,lambda2)    *  PrAprimehi(age1,a1,z1,lambda1,e1)   
	        ENDDO
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO    
	    ENDDO    

		! retirees "e" stays the same for benefit retirement calculation purposes
	    DO age1=RetAge-1, MaxAge-1
	    DO a1=1,na
	    DO z1=1,nz
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        ! Those who die, switch to z2, lambda2 and start at ne/2+1
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	            DBN2(1,Aplo(age1, a1, z1, lambda1, e1), z2,lambda2,ne/2+1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1) + DBN1(age1, a1, z1, lambda1, e1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2)  *  PrAprimelo(age1, a1, z1, lambda1, e1)     
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1) + DBN1(age1, a1, z1, lambda1, e1) & 
	                & *(1.0_DP-survP(age1)) * pr_z(z1,z2) *pr_lambda(lambda1,lambda2)*PrAprimehi(age1,a1,z1,lambda1,e1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, and also e1 since they are retired
	        !e2=e1
	        DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1), z1,lambda1,e1) =  &
	        	& DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1), z1,lambda1,e1) + DBN1(age1, a1, z1, lambda1, e1) &
	            & * survP(age1) * PrAprimelo(age1, a1, z1, lambda1, e1)     
	        DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1), z1,lambda1,e1) =  &
	          	& DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1), z1,lambda1,e1) + DBN1(age1, a1, z1, lambda1, e1) &
	            & * survP(age1) * PrAprimehi(age1, a1, z1, lambda1, e1) 
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    
	    ! Working age agents
	    DO age1=1,RetAge-2
	    DO a1=1,na
	    DO z1=1,nz
	    DO lambda1=1,nlambda
	    DO e1=1, ne
	        ! Those who die, switch to z2, lambda2 and start at ne/2+1
	        DO z2=1,nz
	        DO lambda2=1,nlambda
	        	DBN2(1,Aplo(age1, a1, z1, lambda1, e1), z2,lambda2,ne/2+1)   =  &
	           		& DBN2(1,Aplo(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1) + DBN1(age1, a1, z1, lambda1, e1) &
	                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2)  *  PrAprimelo(age1, a1, z1, lambda1, e1)     
	            DBN2(1,Aphi(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1)   =  &
	           		& DBN2(1,Aphi(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1) + DBN1(age1, a1, z1, lambda1, e1) & 
	                & *(1.0_DP-survP(age1)) * pr_z(z1,z2) *pr_lambda(lambda1,lambda2)*PrAprimehi(age1,a1,z1,lambda1,e1)   
	        ENDDO
	        ENDDO
	        
	        ! Those who live stay at z1, lambda1, but switch to e2
	        DO e2=1, ne
				DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1), z1,lambda1,e2) =  &
	          		& DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1), z1,lambda1,e2) + DBN1(age1, a1, z1, lambda1, e1) &
	                & * survP(age1) * pr_e(e1,e2)   * PrAprimelo(age1, a1, z1, lambda1, e1)     
	            DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1), z1,lambda1,e2) =  &
	          		& DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1), z1,lambda1,e2) + DBN1(age1, a1, z1, lambda1, e1) &
	                & * survP(age1) * pr_e(e1,e2)  *  PrAprimehi(age1, a1, z1, lambda1, e1) 
	        ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    ENDDO
	    
	    
	    DBN_dist = maxval(abs(DBN2-DBN1))
	    DBN1 = DBN2

	    ! Instead of adjusting policy functions to be consistent with current DBN at each iteration
	    ! The code iterates the distribution and only updates policy functions every "update_period"

	    ! Update of policy function with current aggregates
	    IF (iter_indx .ge. update_period) THEN

	    	! Compute aggregates with current distribution
	        QBAR =0.0
	        NBAR =0.0
	        DO age1=1,MaxAge
	        DO z1=1,nz
	        DO a1=1,na
	        DO lambda1=1,nlambda
	        DO e1=1, ne
	             QBAR= QBAR+ DBN1(age1, a1, z1, lambda1, e1) * ( zgrid(z1) *K_mat(a1,z1) )**mu
	             NBAR= NBAR+ DBN1(age1, a1, z1, lambda1, e1) * eff_un(age1, lambda1, e1) * Hours(age1, a1, z1, lambda1,e1)
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
	    	!R = zbrent(Agg_Debt,0.0_dp,0.50_dp,brent_tol) 
	    	if (theta .gt. 1.0_DP) then
	           brent_value = brent(-0.1_DP,0.1_DP,0.5_DP,Agg_Debt, brent_tol,R)
            else
                R = 0.0_DP
	        endif

	    	!!
	    	print*, 'DBN_diff=', DBN_dist, "Agg_Debt", Agg_Debt(R)**0.5_dp, 'R=',R,'P=',P
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
	        !$omp parallel do private(lambdai,ei,ai,tklo,tkhi)
	        DO zi=1,nz
	        DO ai=1,na
	        DO lambdai=1,nlambda
	        DO ei=1, ne
                if ( Aprime(age,ai,zi,lambdai, ei) .ge. amax) then
                    tklo =na-1
                    elseif (Aprime(age,ai,zi,lambdai, ei) .lt. amin) then
                         tklo = 1
                        else
                            tklo = ((Aprime(age,ai,zi,lambdai, ei) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
                endif            
                tkhi = tklo + 1        
                Aplo(age,ai,zi,lambdai, ei)  = tklo
                Aphi(age,ai,zi,lambdai, ei)  = tkhi        
                PrAprimelo(age,ai,zi,lambdai, ei) = ( agrid(tkhi) - Aprime(age,ai,zi,lambdai, ei) ) / ( agrid(tkhi) -agrid(tklo) )
                PrAprimehi(age,ai,zi,lambdai, ei) = ( Aprime(age,ai,zi,lambdai, ei) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )
	        ENDDO
	        ENDDO
	        ENDDO
	        ENDDO
	        ENDDO
	
	        ! Probablities are adjusted to lie in [0,1]
		        PrAprimelo = min (PrAprimelo, 1.0_DP)
		        PrAprimelo = max(PrAprimelo, 0.0_DP)
		        PrAprimehi = min (PrAprimehi, 1.0_DP)
		        PrAprimehi = max(PrAprimehi, 0.0_DP)  

		    ! Reset counter for next update of policy functions
	        iter_indx=0
	    ENDIF
	    
	    iter_indx = iter_indx + 1
	    simutime  = simutime +1 
	 
	ENDDO ! WHILE

	! Write
		!OPEN   (UNIT=3, FILE='agrid', STATUS='replace')
		!WRITE(unit=3, FMT=*) agrid
		!CLOSE (unit=3)
		!
		!OPEN   (UNIT=3, FILE='aprime', STATUS='replace')
		!WRITE(unit=3, FMT=*) Aprime(1, :, nz/2+1, nlambda/2+1, ne/2+1)
		!CLOSE (unit=3)
		!
		!OPEN   (UNIT=3, FILE='aplo', STATUS='replace')
		!WRITE(unit=3, FMT=*) Aplo(1, :,nz/2+1, nlambda/2+1, ne/2+1)
		!CLOSE (unit=3)
		!
		!OPEN   (UNIT=3, FILE='aphi', STATUS='replace')
		!WRITE(unit=3, FMT=*) Aphi(1, :, nz/2+1, nlambda/2+1, ne/2+1)
		!CLOSE (unit=3)
		!
		!OPEN   (UNIT=3, FILE='Pr_aplo', STATUS='replace')
		!WRITE(unit=3, FMT=*) PrAprimelo(1, :, nz/2+1, nlambda/2+1, ne/2+1)
		!CLOSE (unit=3)
		!
		!OPEN   (UNIT=3, FILE='pr_aphi', STATUS='replace')
		!WRITE(unit=3, FMT=*) PrAprimehi(1, :, nz/2+1, nlambda/2+1, ne/2+1)
		!CLOSE (unit=3)

END SUBROUTINE FIND_DBN_EQ


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE COMPUTE_STATS()
	use omp_lib

	IMPLICIT NONE
	INTEGER  :: prctile, group
	REAL(DP), DIMENSION(nz)    :: cdf_Gz_DBN, Capital_by_z 
	REAL(dp), DIMENSION(na,nz) :: DBN_az, wealth
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
	integer  :: constrained_firm_ind(MaxAge,na,nz,nlambda,ne), i
	real(DP) :: Firm_Output(MaxAge,na,nz,nlambda,ne), Firm_Profit(MaxAge,na,nz,nlambda,ne)
	real(DP), dimension(size(DBN1)) :: DBN_vec, Firm_Wealth_vec, CDF_Firm_Wealth, BQ_vec, DBN_bq_vec, CDF_bq
	real(DP) :: FW_top_x(4),  prctile_FW(4), prctile_bq(7), a, b, c, CCDF_c
	real(DP) :: DBN_bq(MaxAge,na,nz,nlambda,ne)
	character(100) :: rowname
	INTEGER, dimension(max_age_category+1) :: age_limit

	!$ call omp_set_num_threads(20)

	
	! Age Brackets
		age_limit = [0, 5, 15, 25, 35, 45, 55, MaxAge ]
	

	! COMPUTE AVERAGE HOURS FOR AGES 25-60 (5-40 IN THE MODEL) INCLUDING NON-WORKERS
	! COMPUTE VARIANCE OF LOG EARNINGS FOR 25-60 FOR THOSE WHO WORK MORE THAN 260 HOURS
	! WHICH CORRESPOND TO 0.055 IN THE MODEL
	pop_25_60        	   = 0.0_DP
	tothours_25_60         = 0.0_DP
	pop_pos_earn_25_60     = 0.0_DP
	tot_log_earnings_25_60 = 0.0_DP 
	Var_Log_Earnings_25_60 = 0.0_DP
	DO age=5,40
	DO ai=1,na
	DO zi=1,nz
	DO lambdai=1,nlambda
	DO ei=1,ne
		tothours_25_60 = tothours_25_60 + DBN1(age, ai, zi, lambdai, ei)  * Hours(age, ai, zi, lambdai,ei)
		pop_25_60      = pop_25_60 +  DBN1(age, ai, zi, lambdai, ei)
		IF (Hours(age, ai, zi, lambdai, ei) .ge. 0.055) THEN
		tot_log_earnings_25_60 = tot_log_earnings_25_60 + DBN1(age, ai, zi, lambdai, ei)  &
		                 		& *  log( Y_h(Hours(age, ai, zi, lambdai, ei),age,lambdai,ei,wage) )
		pop_pos_earn_25_60     = pop_pos_earn_25_60 +  DBN1(age, ai, zi, lambdai, ei)
		ENDIF
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	meanhours_25_60         = tothours_25_60 / pop_25_60
	mean_log_earnings_25_60 = tot_log_earnings_25_60 / pop_pos_earn_25_60

	DO age=5,40
	DO ai=1,na
	DO zi=1,nz
	DO lambdai=1,nlambda
	DO ei=1,ne
		IF (Hours(age, ai, zi, lambdai, ei) .ge. 0.055) THEN
		    Var_Log_Earnings_25_60 =  Var_Log_Earnings_25_60 + DBN1(age, ai, zi, lambdai, ei)  &
		                 			& * ( log( Y_h(Hours(age, ai, zi, lambdai, ei),age,lambdai,ei,wage) ) &
		                 			& -   mean_log_earnings_25_60 ) ** 2.0_DP
		ENDIF
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
		do ai=1,na 
			YGRID(ai,zi) = Y_a(agrid(ai),zgrid(zi))
		enddo 
		enddo
		CALL ComputeLaborUnits(Ebar, wage) 

	MeanWealth 	 = 0.0_DP
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
	DO age=1,MaxAge
	DO zi=1,nz
	DO ai=1,na
	DO lambdai=1,nlambda
	DO ei=1, ne
	    MeanWealth   = MeanWealth   + DBN1(age, ai, zi, lambdai, ei)*agrid(ai)
	    Mean_Capital = Mean_Capital + DBN1(age, ai, zi, lambdai, ei)*K_mat(ai,zi)
	    MeanCons     = MeanCons     + DBN1(age, ai, zi, lambdai, ei)*cons(age, ai, zi, lambdai, ei)   

	    size_by_z(zi)    = size_by_z(zi)    + DBN1(age, ai, zi, lambdai, ei) 
	    Wealth_by_z(zi)  = Wealth_by_z(zi)  + DBN1(age, ai, zi, lambdai, ei) * agrid(ai)
	    Capital_by_z(zi) = Capital_by_z(zi) + DBN1(age, ai, zi, lambdai, ei) * K_mat(ai,zi)

	    MeanReturn           = MeanReturn          + DBN1(age, ai, zi, lambdai, ei) * &
	    							& (R*agrid(ai) + Pr_mat(ai,zi))
	    MeanReturn_by_z(zi)  = MeanReturn_by_z(zi) + DBN1(age, ai, zi, lambdai, ei) * &
	    							& (R*agrid(ai) + Pr_mat(ai,zi)) 
	    
	    MeanATReturn           = MeanATReturn          + DBN1(age, ai, zi, lambdai, ei) * (YGRID(ai,zi)-agrid(ai))
	    MeanATReturn_by_z(zi)  = MeanATReturn_by_z(zi) + DBN1(age, ai, zi, lambdai, ei) * (YGRID(ai,zi)-agrid(ai))
	    

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
	DO age=1,MaxAge
	DO zi=1,nz
	DO ai=1,na
	DO lambdai=1,nlambda
	DO ei=1, ne  

	    VarReturn    = VarReturn +  DBN1(age, ai, zi, lambdai, ei) * agrid(ai)/MeanWealth * &
	    				& ((R*agrid(ai) + Pr_mat(ai,zi))/agrid(ai)-MeanReturn)**2.0_dp

	    Var_K_Return = Var_K_Return +  DBN1(age, ai, zi, lambdai, ei) * K_mat(ai,zi)/MeanWealth * &
	    				& ((R*agrid(ai) + Pr_mat(ai,zi))/K_mat(ai,zi)-MeanATReturn)**2.0_dp

	    VarATReturn  = VarATReturn +  DBN1(age, ai, zi, lambdai, ei) * agrid(ai)/MeanWealth * &
	    				& ((YGRID(ai,zi)-agrid(ai))/agrid(ai)-MeanATReturn)**2.0_dp 

	    Var_AT_K_Return  = Var_AT_K_Return +  DBN1(age, ai, zi, lambdai, ei) * K_mat(ai,zi)/MeanWealth * &
	    				& ((YGRID(ai,zi)-agrid(ai))/K_mat(ai,zi)-MeanATReturn)**2.0_dp 
	    
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
	StdATReturn     = VarATReturn**0.5_DP
	StdReturn       = VarReturn**0.5_DP
	Std_AT_K_Return = Var_AT_K_Return**0.5_DP
	Std_K_Return    = Var_K_Return**0.5_DP


    ! Debt to GDP Ratio
    External_Debt_GDP = 0.0_DP
	DO zi=1,nz
	DO ai=1,na
	    External_Debt_GDP = External_Debt_GDP + sum(DBN1(:, ai, zi, :, :))*abs(K_mat(ai,zi)-agrid(ai))
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
	 
	    DO ai=1,na
        DO zi=1,nz
        DO lambdai=1,nlambda
        DO ei=1,ne
        	size_Age(group)   = size_Age(group)   + DBN1(age,ai,zi,lambdai,ei)
        	size_AZ(group,zi) = size_AZ(group,zi) + DBN1(age,ai,zi,lambdai,ei)

        	A_Age(group)      = A_Age(group)   + DBN1(age,ai,zi,lambdai,ei)*agrid(ai)
        	A_AZ(group,zi)    = A_AZ(group,zi) + DBN1(age,ai,zi,lambdai,ei)*agrid(ai)

        	Ap_Age(group)     = Ap_Age(group)   + DBN1(age,ai,zi,lambdai,ei)*Aprime(age,ai,zi,lambdai,ei)
        	Ap_AZ(group,zi)   = Ap_AZ(group,zi) + DBN1(age,ai,zi,lambdai,ei)*Aprime(age,ai,zi,lambdai,ei)

        	if (age.lt.RetAge) then
        	Y_Age(group)      = Y_Age(group)   + DBN1(age,ai,zi,lambdai,ei)*&
        						& ( YGRID(ai,zi)+ Y_h(Hours(age,ai,zi,lambdai,ei),age,lambdai,ei,Wage))
        	Y_AZ(group,zi)    = Y_AZ(group,zi) + DBN1(age,ai,zi,lambdai,ei)*&
        						& ( YGRID(ai,zi)+ Y_h(Hours(age,ai,zi,lambdai,ei),age,lambdai,ei,Wage))
        	else 
        	Y_Age(group)      = Y_Age(group)   + DBN1(age,ai,zi,lambdai,ei)*&
        						& ( YGRID(ai,zi)+ RetY_lambda_e(lambdai,ei) )
        	Y_AZ(group,zi)    = Y_AZ(group,zi) + DBN1(age,ai,zi,lambdai,ei)*&
        						& ( YGRID(ai,zi)+ RetY_lambda_e(lambdai,ei) )
        	endif

        	if (ai.le.prctile_ai_ind(90)) then 
        		size_W(1) = size_W(1) + DBN1(age,ai,zi,lambdai,ei)
        		A_W(1)    = A_W(1)    + DBN1(age,ai,zi,lambdai,ei)*agrid(ai)
        		Ap_W(1)   = Ap_W(1)   + DBN1(age,ai,zi,lambdai,ei)*Aprime(age,ai,zi,lambdai,ei)
        		if (age.lt.RetAge) then 
        		Y_W(1)    = Y_W(1)    + DBN1(age,ai,zi,lambdai,ei)*&
        					(YGRID(ai,zi)+ Y_h(Hours(age,ai,zi,lambdai,ei),age,lambdai,ei,Wage))
        		else 
        		Y_W(1)    = Y_W(1)    + DBN1(age,ai,zi,lambdai,ei)*&
        					(YGRID(ai,zi)+ RetY_lambda_e(lambdai,ei))
        		endif 
        	else if  ((ai.gt.prctile_ai_ind(90)).and.(ai.le.prctile_ai_ind(99))) then
        		size_W(2) = size_W(2) + DBN1(age,ai,zi,lambdai,ei)
        		A_W(2)    = A_W(2)    + DBN1(age,ai,zi,lambdai,ei)*agrid(ai)
        		Ap_W(2)   = Ap_W(2)   + DBN1(age,ai,zi,lambdai,ei)*Aprime(age,ai,zi,lambdai,ei)
        		if (age.lt.RetAge) then 
        		Y_W(2)    = Y_W(2)    + DBN1(age,ai,zi,lambdai,ei)*&
        					(YGRID(ai,zi)+ Y_h(Hours(age,ai,zi,lambdai,ei),age,lambdai,ei,Wage))
        		else 
        		Y_W(2)    = Y_W(2)    + DBN1(age,ai,zi,lambdai,ei)*&
        					(YGRID(ai,zi)+ RetY_lambda_e(lambdai,ei))
        		endif
        	else 
        		size_W(3) = size_W(3) + DBN1(age,ai,zi,lambdai,ei)
        		A_W(3)    = A_W(3)    + DBN1(age,ai,zi,lambdai,ei)*agrid(ai)
        		Ap_W(3)   = Ap_W(3)   + DBN1(age,ai,zi,lambdai,ei)*Aprime(age,ai,zi,lambdai,ei)
        		if (age.lt.RetAge) then 
        		Y_W(3)    = Y_W(3)    + DBN1(age,ai,zi,lambdai,ei)*&
        					(YGRID(ai,zi)+ Y_h(Hours(age,ai,zi,lambdai,ei),age,lambdai,ei,Wage))
        		else 
        		Y_W(3)    = Y_W(3)    + DBN1(age,ai,zi,lambdai,ei)*&
        					(YGRID(ai,zi)+ RetY_lambda_e(lambdai,ei))
        		endif
        	endif 
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

	! Leverage Ratio and fraction of contrainted firms 
	leverage_age_z = 0.0_dp 
	size_by_age_z  = 0.0_dp 
	constrained_firms_age_z = 0.0_dp
	constrained_firm_ind = 0
	do age = 1,MaxAge 
	do zi  = 1,nz 
        DO lambdai=1,nlambda
        DO ei=1,ne
        DO ai=1,na
        	size_by_age(age)       = size_by_age(age)       + DBN1(age,ai,zi,lambdai,ei)
			size_by_age_z(age,zi)  = size_by_age_z(age,zi)  + DBN1(age,ai,zi,lambdai,ei)
			leverage_age_z(age,zi) = leverage_age_z(age,zi) + DBN1(age,ai,zi,lambdai,ei)*K_mat(ai,zi)/agrid(ai)
			if (K_mat(ai,zi).ge.(theta*agrid(ai))) then 
				constrained_firms_age(age)      = constrained_firms_age(age)      + DBN1(age,ai,zi,lambdai,ei)
				constrained_firms_age_z(age,zi) = constrained_firms_age_z(age,zi) + DBN1(age,ai,zi,lambdai,ei)
				constrained_firm_ind(age,ai,zi,lambdai,ei) = 1
			endif 
			Firm_Output(age,ai,zi,lambdai,ei) = zgrid(zi)*K_mat(ai,zi)
			Firm_Profit(age,ai,zi,lambdai,ei) = Pr_mat(ai,zi)
		enddo
		enddo 
		enddo 
	enddo 
	enddo 
	leverage_age_z = leverage_age_z/size_by_age_z
	constrained_firms_age_z = constrained_firms_age_z/size_by_age_z 
	constrained_firms_age   = constrained_firms_age/size_by_age 


	! Percentage of the population above wealth tax threshold
		! Compute distribution of agents by (a,z)
		DBN_az = sum(sum(sum(DBN1,5),4),1)
		! Compute mean before tax wealth
		Wealth = Wealth_Matrix(R,P)
		! Compute share of agents above threshold
		Threshold_Share = 0.0_dp
		do ai=1,na
		do zi=1,nz 
			if (Wealth(ai,zi).gt.Y_a_threshold) then 
				Threshold_Share = Threshold_Share + DBN_az(ai,zi)
			end if 
		end do 
		end do 

	! Distribution of Assets
		DO ai=1,na
		     pr_a_dbn(ai)          = sum(DBN1(:,ai,:,:,:)) 
		     cdf_a_dbn(ai)         = sum( pr_a_dbn(1:ai) )      
		     tot_a_by_grid(ai)     = sum(DBN1(:,ai,:,:,:) * agrid(ai) )
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
			Firm_Wealth(:,ai,:,:,:) = V_Pr(:,ai,:,:,:) + (1.0_dp+R)*agrid(ai)
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
			WRITE(UNIT=11, FMT=*) 100_dp*prctile_FW(i),FW_top_x(i),FW_top_x_share(i), CCDF_c
		enddo 

			CLOSE(UNIT=11)

	! Distribution of bequest
		! Distribution of bequest (matrix)	
		do ai=1,MaxAge
			DBN_bq(age,:,:,:,:) = DBN1(age,:,:,:,:)*(1.0_DP-survP(age))
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
			WRITE(UNIT=11, FMT=*) 'Mean_Bequest/Wealth= '		, Mean_Bequest/MeanWealth 
			WRITE(UNIT=11, FMT=*) 'Mean_Bequest/PV_Wealth= '	, Mean_Bequest/Mean_Firm_Wealth 
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
	





	!print*, 'MeanReturn=',MeanReturn, 'StdReturn=', StdReturn
	!print*,'MeanReturn_by_z=',MeanReturn_by_z

	SSE_Moments = (1.0-Wealth_Output/3.0_DP)**2.0_DP  + (1.0_DP-prct1_wealth/0.34_DP)**2.0_DP  + (prct10_wealth-0.69_DP)**2.0_DP &
                   & + (1.0_DP-Std_Log_Earnings_25_60 / 0.8_DP)**2.0_DP + (1.0_DP-meanhours_25_60/0.4_DP)**2.0_DP &
                   & + (1.0_DP-MeanReturn/0.069_DP)**2.0_DP
	!print*,''
	!print*,"Current parameters"
	!print*,'beta',beta,'rho_z',rho_z,'sigma_z',sigma_z_eps,'sigma_lam',sigma_lambda_eps,'phi',phi
	print*,"Statistics"
	print*,'Debt/GDP',External_Debt_GDP,'W/GDP',Wealth_Output,'Top 1%',prct1_wealth,'Top 10%',prct10_wealth
	print*,'STD Labor Earnings',Std_Log_Earnings_25_60,'Mean Labor (hours 25-60)',meanhours_25_60,'MeanReturn',MeanReturn
	print*,'Moments',SSE_Moments 
	!print*,''

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

	! Save files of constrained index, output and profits
	if (solving_bench.eq.1) then
		OPEN(UNIT=1,  FILE=trim(Result_Folder)//'constrained_ind_bench'  , STATUS='replace')
		OPEN(UNIT=2,  FILE=trim(Result_Folder)//'firm_output_bench'      , STATUS='replace')
		OPEN(UNIT=3,  FILE=trim(Result_Folder)//'firm_profit_bench' 	    , STATUS='replace')
	else 
		OPEN(UNIT=1,  FILE=trim(Result_Folder)//'constrained_ind_exp'    , STATUS='replace')
		OPEN(UNIT=2,  FILE=trim(Result_Folder)//'firm_output_exp'        , STATUS='replace')
		OPEN(UNIT=3,  FILE=trim(Result_Folder)//'firm_profit_exp'  	    , STATUS='replace')
	endif
		WRITE(UNIT=1,FMT=*) constrained_firm_ind
		WRITE(UNIT=2,FMT=*) Firm_Output
		WRITE(UNIT=3,FMT=*) Firm_Profit
		CLOSE(UNIT=1)
		CLOSE(UNIT=2)
		CLOSE(UNIT=3)


END SUBROUTINE COMPUTE_STATS


!========================================================================================
!========================================================================================
!========================================================================================

SUBROUTINE EGM_RETIREMENT_WORKING_PERIOD()
	use omp_lib
	IMPLICIT NONE
	REAL(DP) :: brentvalue, C_foc, H_min, euler_power
	INTEGER  :: tempai, sw 
	REAL(DP), DIMENSION(na_t+1)  :: EndoCons, EndoYgrid, EndoHours
	INTEGER , DIMENSION(na_t+1)  :: sort_ind 
	REAL(DP), DIMENSION(na_t,nz) :: Wealth
	REAL(DP), DIMENSION(5)       :: state_FOC
	REAL(DP), DIMENSION(6)       :: par_FOC
	integer  :: age, ai, zi, lambdai, ei

	!$ call omp_set_num_threads(11)

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
		Wealth = Wealth_Matrix_t(R,P)

		! print*, 'R=',R,'P=',P, 'W=',wage, 'na=', na, 'na_t=', na_t
	!========================================================================================
	!------RETIREMENT PERIOD-----------------------------------------------------------------

	! SET HOURS TO ZERO SO THAT RETIRED HAS ZERO HOURS
	Hours_t = 0.0_DP

	! Last period of life
	age=MaxAge
	!$omp parallel do private(lambdai,ei,ai)
	DO zi=1,nz
    DO lambdai=1,nlambda
    DO ei=1,ne
    DO ai=1,na_t
        Cons_t(age,ai,zi,lambdai,ei) =  YGRID_t(ai, zi) + RetY_lambda_e(lambdai,ei) 
	ENDDO ! ai
    ENDDO ! ei
	ENDDO ! lambdai
	ENDDO ! zi
	Aprime_t(age, :, :, :, :) = 0.0_DP
	
	! Rest of retirement
	DO age=MaxAge-1,RetAge,-1
	!$omp parallel do private(lambdai,ei,ai,EndoCons,EndoYgrid,sw,sort_ind,tempai,state_FOC,par_FOC)
    DO zi=1,nz
    DO lambdai=1,nlambda
    DO ei=1,ne
    	! Endogenous grid and consumption are initialized
	    EndoCons  = big_p 
	    EndoYgrid = big_p 
	    sw 		  = 0
    DO ai=1,na_t
		if (abs(Wealth(ai,zi)-Y_a_threshold).lt.1e-8) then 
    		! Consumption on endogenous grid and implied asset income under tauW_bt
    		EndoCons(ai)  = Cons_t(age+1,ai,zi,lambdai,ei) * (beta*survP(age)*MB_a_bt(agrid(ai),zgrid(zi)) )**euler_power
	        EndoYgrid(ai) = agrid_t(ai) +  EndoCons(ai) - RetY_lambda_e(lambdai,ei)
	        ! Consumption on endogenous grid and implied asset income under tauW_at
	        EndoCons(na_t+1)  = Cons_t(age+1,ai,zi,lambdai,ei)*(beta*survP(age)*MB_a_at(agrid(ai),zgrid(zi)))**euler_power
	    	EndoYgrid(na_t+1) = agrid_t(ai) +  EndoCons(na_t+1) - RetY_lambda_e(lambdai,ei)
	    	! Set the flag!
	    	sw                = 1 
	    else 
	    	! Consumption on endogenous grid and implied asset income
	    	EndoCons(ai)  = Cons_t(age+1,ai,zi,lambdai,ei) * (beta*survP(age)*MBGRID_t(ai,zi))**euler_power
	        EndoYgrid(ai) = agrid_t(ai) +  EndoCons(ai) - RetY_lambda_e(lambdai,ei)
	    end if 
	ENDDO ! ai
	
	! Sort endogenous grid for interpolation
	call Sort(na_t+1,EndoYgrid,EndoYgrid,sort_ind)
	EndoCons = EndoCons(sort_ind)

	! Find  decision rules on exogenous grids
		! decision rules are obtained taking care of extrapolations
		tempai=1           
		DO WHILE ( YGRID_t(tempai,zi) .lt. EndoYgrid(1) )
            tempai = tempai + 1
        ENDDO
	                
		DO ai=tempai,na_t              
		    ! CONSUMPTION ON EXOGENOUS GRIDS 
		    Cons_t(age, ai, zi, lambdai, ei)  = Linear_Int(EndoYgrid(1:na_t+sw), EndoCons(1:na_t+sw),na_t+sw, YGRID_t(ai,zi))
		    Aprime_t(age, ai, zi, lambdai,ei) = YGRID_t(ai,zi)+ RetY_lambda_e(lambdai,ei) - Cons_t(age, ai, zi, lambdai, ei)
		    
		    If (Aprime_t(age, ai, zi, lambdai,ei)  .lt. amin) then
		    	Aprime_t(age, ai, zi, lambdai,ei) = amin
	            Cons_t(age, ai, zi, lambdai,ei)   = YGRID_t(ai,zi) + RetY_lambda_e(lambdai,ei) - Aprime_t(age, ai, zi, lambdai,ei) 
				IF (Cons_t(age, ai, zi, lambdai, ei) .le. 0.0_DP)  THEN
				    print*,'r1: Cons(age, ai, zi, lambdai,ei)=',Cons_t(age, ai, zi, lambdai, ei)
				ENDIF                   
	        endif                                                         
		ENDDO ! ai  

        ai=1           
        DO WHILE ( YGRID_t(ai,zi) .lt. EndoYgrid(1) )
			! Solve for a' directly by solving the Euler equation for retirement FOC_R
			state_FOC  = (/age,ai,zi,lambdai,ei/)
			brentvalue = brent_p(min(amin,YGRID_t(ai,zi)), (amin+YGRID_t(ai,zi))/2.0_DP , &
			                & YGRID_t(ai,zi)+RetY_lambda_e(lambdai,ei) *0.95_DP, &
			                & FOC_R, brent_tol, Aprime_t(age, ai, zi, lambdai,ei),state_FOC)
			
			Cons_t(age, ai, zi, lambdai, ei) =  YGRID_t(ai,zi)  + RetY_lambda_e(lambdai,ei)  - Aprime_t(age, ai, zi, lambdai,ei)
			
			IF (Cons_t(age, ai, zi, lambdai, ei) .le. 0.0_DP)  THEN
			    print*,'r2: Cons(age, ai, zi, lambdai,ei)=',Cons_t(age, ai, zi, lambdai,ei)
			    print*,'Aprime(age, ai, zi, lambdai,ei)=',Aprime_t(age, ai, zi, lambdai,ei)
			    print*,'YGRID(ai,zi)+RetY_lambda_e(lambdai,ei)=',YGRID_t(ai,zi)+RetY_lambda_e(lambdai,ei)
			    print*,'YGRID(ai,zi)=',YGRID(ai,zi),'EndoYgrid(1)=',EndoYgrid(1)
			    print*,'RetY_lambda_e(lambdai,ei)=',RetY_lambda_e(lambdai,ei)
			    print*,'lambdai=',lambdai
			ENDIF                   
			ai = ai + 1
		ENDDO  
	               
    ENDDO ! ei     
    ENDDO ! lambda
	ENDDO ! zi
    ENDDO !age


	!------RETIREMENT PERIOD ENDS------------------------------------------------------------
	!========================================================================================
	
	!========================================================================================
	!------Working Period Starts-------------------------------------------------------------

	DO age=RetAge-1,1,-1
	! !$omp parallel do private(lambdai,ei,ai,EndoCons,EndoHours,EndoYgrid,sw,sort_ind,tempai,C_foc,state_FOC,par_FOC)
    DO zi=1,nz
    DO lambdai=1,nlambda
    DO ei=1,ne	
    	! Endogenous grid and consumption are initialized
	    EndoCons  = big_p 
	    EndoYgrid = big_p 
	    EndoHours = big_p 
	    sw 		  = 0                    
        DO ai=1,na_t
        state_FOC  = (/age,ai,zi,lambdai,ei/)
		if (abs(Wealth(ai,zi)-Y_a_threshold).lt.1e-8) then 
	    	! Below threshold
			call EGM_Working_Period( MB_a_bt(agrid_t(ai),zgrid(zi)) , H_min , state_FOC , & 
			      & EndoCons(ai), EndoHours(ai) , EndoYgrid(ai)  )
			
			! Above threshold
			call EGM_Working_Period( MB_a_at(agrid_t(ai),zgrid(zi)) , H_min , state_FOC , & 
			      & EndoCons(na_t+1), EndoHours(na_t+1) , EndoYgrid(na_t+1)  )

			! Set the flag!
	    	sw = 1 

	    	!print*, ' '
	    	!print*, State_FOC 
	    	!print*, MB_a_bt(agrid_t(ai),zgrid(zi)), MB_a_at(agrid_t(ai),zgrid(zi))
		else 
			! Usual EGM
			call EGM_Working_Period( MBGRID_t(ai,zi) , H_min , state_FOC , & 
			      & EndoCons(ai), EndoHours(ai) , EndoYgrid(ai)  )
	
		end if 

    	ENDDO ! ai

	    ! !$omp critical
		if (any(isnan(EndoCons))) then 
			print*, "isnan - Consumption endogenous"
			print*, age,zi,ai,lambdai,ei
			print*, EndoCons
			STOP 
		end if 
		! !$omp end critical

    ! Sort endogenous grid for interpolation
	call Sort(na_t+1,EndoYgrid,EndoYgrid,sort_ind)
	EndoHours = EndoHours(sort_ind)
	EndoCons  = EndoCons(sort_ind)
! 	print*, ' '
! 	do ai=1,na_t+1
! 		print*, EndoHours(ai), EndoCons(ai), EndoYgrid(ai)
! 	end do
! 	print*, ' '

		if (any(isnan(EndoCons)).or.any(isnan(EndoHours)).or.any(isnan(EndoYgrid))) then 
			print*, "isnan - Consumption working 4"
			print*, age,lambdai,ai,zi,ei
			STOP 
		end if 

    	! Find  decision rules on exogenous grids
        tempai=1           
        DO WHILE ( YGRID_t(tempai,zi) .lt. EndoYgrid(1) )
              tempai = tempai + 1
        ENDDO
	                
		! decision rules are obtained taking care of extrapolations
		DO ai=tempai,na_t 
			! Interpolate for value of consumption in exogenous grid
				Cons_t(age,ai,zi,lambdai,ei) = Linear_Int(EndoYgrid(1:na_t+sw), EndoCons(1:na_t+sw),na_t+sw, YGRID_t(ai,zi)) 

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
				if (Cons_t(age,ai,zi,lambdai,ei).ge.C_foc) then
					Hours_t(age,ai,zi,lambdai,ei) = 0.0_dp
				else
					! Auxiliary variables for solving FOC for hours 
					par_FOC(1:5) = (/age,ai,zi,lambdai,ei/) 
					par_FOC(6)   = Cons_t(age, ai, zi, lambdai,ei)
					! FOC
					brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_H, brent_tol, Hours_t(age,ai,zi,lambdai,ei) , par_FOC ) 
				end if
			else 
				if (NSU_Switch.eqv..true.) then 
					Hours_t(age,ai,zi,lambdai,ei) = max( 0.0_dp , &
						&  1.0_DP - (1.0_DP-gamma)*Cons_t(age,ai,zi,lambdai,ei)/(gamma*psi*yh(age,lambdai,ei)) )
				else 
					Hours_t(age,ai,zi,lambdai,ei) = max( 0.0_DP , &
						&  1.0_DP - phi*Cons_t(age,ai,zi,lambdai,ei)/(psi*yh(age, lambdai,ei)) )
				end if 
			end if 

			if (Hours_t(age,ai,zi,lambdai,ei).ge.1.0_dp) then 
				print*, ' '
				print*, 'Hours highers than 1', age,ai,zi,lambdai,ei
				print*, Cons_t(age,ai,zi,lambdai,ei)
				STOP
			endif 


			! Savings 
				Aprime_t(age, ai, zi, lambdai,ei) = YGRID_t(ai,zi)  + Y_h(Hours_t(age, ai, zi, lambdai,ei),age,lambdai,ei,wage)  & 
		                    					& - Cons_t(age, ai, zi, lambdai,ei) 
		                    
		    If (Aprime_t(age, ai, zi, lambdai,ei)  .lt. amin) then

		    	print*, ' Aprime was below minimum!!!!'
            	Aprime_t(age, ai, zi, lambdai,ei) = amin
		         
	           	! Compute hours
		        if (NSU_Switch.and.(Progressive_Tax_Switch.eqv..false.)) then 
		        	Hours_t(age, ai, zi, lambdai,ei)  = max( 0.0_dp , &
		        			&  gamma - (1.0_DP-gamma) *( YGRID_t(ai,zi) - Aprime_t(age, ai, zi, lambdai,ei)) / (psi*yh(age, lambdai,ei)) )
                else               	        
                	!compute  hours using FOC_HA                              
		        	par_FOC(1:5) = (/age,ai,zi,lambdai,ei/) 
		        	par_FOC(6)   = amin
		        	brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_HA, brent_tol, Hours_t(age, ai, zi, lambdai,ei), par_FOC)
		        end if  

	            Cons_t(age,ai,zi,lambdai,ei) = YGRID_t(ai,zi)+  Y_h(Hours_t(age, ai, zi, lambdai,ei),age,lambdai,ei,wage)  &
		                            		  & - Aprime_t(age, ai, zi, lambdai,ei)      

				IF (Cons_t(age, ai, zi, lambdai,ei) .le. 0.0_DP)  THEN
					print*,'w1: Cons(age, ai, zi, lambdai,ei)=',Cons_t(age, ai, zi, lambdai,ei)
					STOP
				ENDIF 
				if (Hours_t(age,ai,zi,lambdai,ei).ge.1.0_dp) then 
					print*, ' '
					print*, 'w1 Hours highers than 1', age,ai,zi,lambdai,ei
					print*, Cons_t(age,ai,zi,lambdai,ei)
					STOP
				endif                   
		     endif  
		    !$omp critical
		    IF (Cons_t(age, ai, zi, lambdai,ei) .le. 0.0_DP)  THEN
				print*,'w1: Cons(age, ai, zi, lambdai,ei)=',Cons_t(age, ai, zi, lambdai,ei)
				STOP
			ENDIF 
			if (Hours_t(age,ai,zi,lambdai,ei).ge.1.0_dp) then 
				print*, ' '
				print*, 'w1 Hours highers than 1', age,ai,zi,lambdai,ei
				print*, Cons_t(age,ai,zi,lambdai,ei)
				STOP
			endif    
			!$omp end critical
		ENDDO ! ai   

		!$omp critical
		if (any(isnan(Cons_t(age,:,zi,lambdai,ei)))) then 
			print*, "isnan - Consumption working 3"
			print*, age,lambdai,ai,zi,ei
			print*, 'Cons'
			print*, Cons_t(age,:,zi,lambdai,ei)
			print*, 'Hours'
			print*, Hours_t(age,:,zi,lambdai,ei) 
			STOP 
		end if         
		!$omp end critical        

		ai=1           
        DO WHILE ( YGRID_t(ai,zi) .lt. EndoYgrid(1) )
        	! print*, ' Extrapolation between YGRID and EndoYgrid!!!!'
	        ! Solve for the Euler equation directly
	        state_FOC  = (/age,ai,zi,lambdai,ei/)
			brentvalue = brent_p( min(amin,YGRID_t(ai,zi))   ,  (amin+YGRID_t(ai,zi))/2.0_DP  ,  &
                             	& YGRID_t(ai,zi)  + psi * ( yh(age, lambdai,ei)*0.95_DP )**(1.0_DP-tauPL)  ,  &
	                             & FOC_WH, brent_tol, Aprime_t(age, ai, zi, lambdai,ei) , state_FOC )

			! Compute hours
	        if (NSU_Switch.and.(Progressive_Tax_Switch.eqv..false.)) then 
	        	Hours_t(age, ai, zi, lambdai,ei)  = max( 0.0_dp , &
	        			&  gamma - (1.0_DP-gamma) *( YGRID_t(ai,zi) - Aprime_t(age, ai, zi, lambdai,ei)) / (psi*yh(age, lambdai,ei)) )
            else               	        
				!compute  hours using FOC_HA                              
				par_FOC(1:5) = (/age,ai,zi,lambdai,ei/)
				par_FOC(6)   = Aprime_t(age, ai, zi, lambdai,ei)
	            brentvalue = brent_p(H_min, 0.4_DP, 0.99_DP, FOC_HA, brent_tol, Hours_t(age, ai, zi, lambdai,ei), par_FOC) 
 			end if 

            Cons_t(age, ai, zi, lambdai,ei)=  YGRID_t(ai,zi) + Y_h(Hours_t(age, ai, zi, lambdai,ei),age,lambdai,ei,wage)  &
                           						 & - Aprime_t(age, ai, zi, lambdai,ei)
           IF (Cons_t(age, ai, zi, lambdai,ei) .le. 0.0_DP)  THEN
                print*,'w2:',age,zi,lambdai,ei,ai, 'Cons(age, ai, zi, lambdai,ei)=',Cons_t(age, ai, zi, lambdai,ei), &
                    & 'Aprime(age, ai, zi, lambdai,ei)=',Aprime_t(age, ai, zi, lambdai,ei), &
                    & 'Hours(age, ai, zi, lambdai,ei)=', Hours_t(age, ai, zi, lambdai,ei), &
                    & 'yh(age, lambdai,ei)=', yh(age, lambdai,ei),'YGRID(ai,zi)=',YGRID_t(ai,zi)
                !pause
                print*, "there is a problem in line 2063"
                STOP
            ENDIF 

            if (Hours_t(age,ai,zi,lambdai,ei).ge.1.0_dp) then 
				print*, ' '
				print*, 'w2 Hours highers than 1', age,ai,zi,lambdai,ei
				print*, Cons_t(age,ai,zi,lambdai,ei)
				STOP
			endif 

            ai = ai + 1
        ENDDO  

	    if (any(isnan(Cons_t(age,:,zi,lambdai,ei)))) then 
			print*, "isnan - Consumption working 2"
			print*, age,lambdai,ai,zi,ei
			print*, 'Cons'
			print*, Cons_t(age,:,zi,lambdai,ei)
			print*, 'Hours'
			print*, Hours_t(age,:,zi,lambdai,ei) 
			STOP 
		end if 

	                 
    ENDDO !ei         
    ENDDO !lambdai
	ENDDO !zi
    ENDDO !age


	! Interpolate to get values of policy functions on agrid (note that policy functions are defined on agrid_t)
	
	if (Y_a_threshold.eq.0.0_dp) then 
		Cons   = Cons_t
		Hours  = Hours_t
		Aprime = Aprime_t
	else 
		DO age=1,MaxAge
	    DO lambdai=1,nlambda
	    DO zi=1,nz
	    DO ei=1,ne	                
	    DO ai=1,na
			Cons(age, ai, zi, lambdai,ei)   = Linear_Int(YGRID_t(:,zi) , Cons_t(age,:,zi,lambdai,ei)  , na_t , YGRID(ai,zi))
	    	Hours(age, ai, zi, lambdai,ei)  = Linear_Int(YGRID_t(:,zi) , Hours_t(age,:,zi,lambdai,ei) , na_t , YGRID(ai,zi))
	    	Aprime(age, ai, zi, lambdai,ei) = Linear_Int(YGRID_t(:,zi) , Aprime_t(age,:,zi,lambdai,ei) , na_t , YGRID(ai,zi))
		ENDDO !ai         
		ENDDO !ei         
	    ENDDO !zi
	    ENDDO !lambdai
		ENDDO !age
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
	real(dp), intent(in)  :: DBN_in(MaxAge, na, nz, nlambda, ne)
	real(dp), intent(out) :: Y_a_threshold_out
	real(dp), dimension(na,nz) :: DBN_az, Wealth
	real(dp)                   :: Mean_Wealth
	integer                    :: prctile, a_ind, z_ind

	if (bench_indx.eq.1) then 
		Y_a_threshold_out = 0
	else 
		! Compute distribution of agents by (a,z)
		DBN_az = sum(sum(sum(DBN_in,5),4),1)
		! Compute mean before tax wealth
		Wealth = Wealth_Matrix(R,P)
		! Mean Wealth
		Mean_Wealth = sum( Wealth*DBN_az )
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
	real(dp), intent(in)       :: R_in, P_in
	real(dp), dimension(na,nz) :: K_Matrix
	integer :: i,j

	! K_Matrix = min( theta*spread(agrid,2,nz) , (mu*P_in*spread(zgrid,1,na)**mu/(R_in+DepRate))**(1.0_dp/(1.0_dp-mu)) )

	do j=1,nz
	do i=1,na 
		K_Matrix(i,j) = min( theta*agrid(i) , (mu*P_in*zgrid(j)**mu/(R_in+DepRate))**(1.0_dp/(1.0_dp-mu)) )
	enddo
	enddo 

end Function K_Matrix

Function K_Matrix_t(R_in,P_in)
	Implicit None 
	real(dp), intent(in)         :: R_in, P_in
	real(dp), dimension(na_t,nz) :: K_Matrix_t
	integer :: i,j

	! K_Matrix_t = min( theta*spread(agrid_t,2,nz) , (mu*P_in*spread(zgrid,1,na_t)**mu/(R_in+DepRate))**(1.0_dp/(1.0_dp-mu)) )

	do j=1,nz
	do i=1,na_t
		K_Matrix_t(i,j) = min( theta*agrid_t(i) , (mu*P_in*zgrid(j)**mu/(R_in+DepRate))**(1.0_dp/(1.0_dp-mu)) )
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
	real(dp), intent(in)       :: R_in, P_in
	real(dp), dimension(na,nz) :: Profit_Matrix
	real(dp), dimension(na,nz) :: K
	integer :: i,j

	K = K_matrix(R_in,P_in)

	! Profit_Matrix = P_in*(spread(zgrid,1,na)*K)**mu - (R_in+DepRate)*K

	do j=1,nz 
	do i=1,na 
		Profit_Matrix(i,j) = P_in*(zgrid(j)*K(i,j))**mu - (R_in+DepRate)*K(i,j)
	enddo 
	enddo 

end Function Profit_Matrix

Function Profit_Matrix_t(R_in,P_in)
	Implicit None 
	real(dp), intent(in)         :: R_in, P_in
	real(dp), dimension(na_t,nz) :: Profit_Matrix_t
	real(dp), dimension(na_t,nz) :: K
	integer :: i,j

	K = K_matrix_t(R_in,P_in)

	! Profit_Matrix_t = P_in*(spread(zgrid,1,na_t)*K)**mu - (R_in+DepRate)*K
	do j=1,nz 
	do i=1,na_t 
		Profit_Matrix_t(i,j) = P_in*(zgrid(j)*K(i,j))**mu - (R_in+DepRate)*K(i,j)
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
	real(dp), intent(in)       :: R_in, P_in
	real(dp), dimension(na,nz) :: Wealth_Matrix
	real(dp), dimension(na,nz) :: K, Pr
	integer :: i,j

	K = K_matrix(R_in,P_in)

	Pr = Profit_Matrix(R_in,P_in)

	Wealth_Matrix = ( (1.0_dp+R_in*(1.0_DP-tauK))*spread(agrid,2,nz) +  Pr *(1.0_DP-tauK) )

end Function Wealth_Matrix

Function Wealth_Matrix_t(R_in,P_in)
	Implicit None 
	real(dp), intent(in)         :: R_in, P_in
	real(dp), dimension(na_t,nz) :: Wealth_Matrix_t
	real(dp), dimension(na_t,nz) :: K, Pr
	integer :: i,j

	K = K_matrix_t(R_in,P_in)

	Pr = Profit_Matrix_t(R_in,P_in)

	Wealth_Matrix_t = ( (1.0_dp+R_in*(1.0_DP-tauK))*spread(agrid_t,2,nz) +  Pr *(1.0_DP-tauK) )

end Function Wealth_Matrix_t

!========================================================================================
!========================================================================================
! This function computes the difference between capital demand and supply
! The function uses the interest rate R and implicitely the price of intermediate good P and distribution DBN1

Function Agg_Debt(R_in)
	Implicit None 
	real(dp), intent(in) :: R_in
	real(dp)             :: Agg_Debt
	real(dp), dimension(na,nz) :: DBN_az
	real(dp)             :: Wealth 

	DBN_az   = sum(sum(sum(DBN1,5),4),1)

	Wealth   = sum( sum(sum(sum(sum(DBN1,5),4),3),1)*agrid )

	Agg_Debt = (sum(DBN_az*(K_matrix(R_in,P)-spread(agrid,2,nz)))/Wealth)**2.0_DP 

end Function Agg_Debt


!========================================================================================
!========================================================================================
!========================================================================================
! THIS YGRID and Marginal Benefit of Investment GRID 
! NEEDS TO BE COMPUTED FOR EACH TIME THE INTEREST RATE "P" IS UPDATED. 

SUBROUTINE FORM_Y_MB_GRID(TYGRID,TMBGRID,TYGRID_t,TMBGRID_t)
	IMPLICIT NONE
	REAL(DP), DIMENSION(na,nz),   INTENT(OUT) :: TYGRID, TMBGRID
	REAL(DP), DIMENSION(na_t,nz), INTENT(OUT) :: TYGRID_t, TMBGRID_t
	!REAL(DP), INTENT(IN) :: P
	!integer :: ai, zi

	DO zi=1,nz
		DO ai=1,na
			TYGRID(ai,zi)  = Y_a(agrid(ai),zgrid(zi))
			TMBGRID(ai,zi) = MB_a(agrid(ai),zgrid(zi))
		ENDDO 
		if (Y_a_threshold.eq.0.0_dp) then
			TYGRID_t  = TYGRID
			TMBGRID_t = TMBGRID 
		else 
			DO ai=1,na_t
				TYGRID_t(ai,zi)  = Y_a(agrid_t(ai),zgrid(zi))
				TMBGRID_t(ai,zi) = MB_a(agrid_t(ai),zgrid(zi))
			ENDDO
		endif 
	ENDDO

	!print *, "Grid for asset income"
	!do ai=1,na
	!	write(*,*) TMBGRID_t(ai,:)
 	!end do
	!pause

END SUBROUTINE FORM_Y_MB_GRID


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
	RetY_lambda_e = phi_lambda_e  * Ebart 
	IF ((KeepSSatBench .eq. 1) .AND. (solving_bench .eq. 0)) THEN
    	RetY_lambda_e = phi_lambda_e  * Ebar_bench
	ENDIF

	OPEN  (UNIT=1,  FILE=trim(Result_Folder)//'Ret_Y'  , STATUS='replace')
	WRITE (UNIT=1,  FMT=*) RetY_lambda_e
	CLOSE (unit=1)

END SUBROUTINE ComputeLaborUnits


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE  INITIALIZE()
	IMPLICIT NONE

	!integer::  lambdai
	REAL(DP) :: lambdaBAR, tempno 
	REAL(DP) :: m, Rh, start_timet, finish_timet
	INTEGER  :: ee0, ee1, ee2, zindx1, zindx2, lambdaindx1, lambdaindx2, diff_array, eindx1, eindx2
	INTEGER, DIMENSION(RetAge) :: agevec
	
	! Initiliaze grids for z, lamda and e	
		CALL tauchen(mtauchen_z,rho_z,sigma_z_eps,nz,zgrid,pr_z,Gz)
		CALL tauchen(mtauchen,rho_E,sigma_e_eps,ne,egrid,pr_e,Ge)
		CALL tauchen(mtauchen,rho_lambda,sigma_lambda_eps,nlambda,lambdagrid,pr_lambda,Glambda)

		! Tauchen gives grids for the log of the variables. Exponentiate to adjust
		zgrid      = exp(zgrid) + mu_z
		egrid      = exp(egrid)
		lambdagrid = exp(lambdagrid) 

		! Obtain CDF of invariant distribution and transition matrix
			! Entrepreneurial ability
			DO ee0 = 1,nz
			    cdf_Gz(ee0) = sum(Gz(1:ee0))
			    DO ee1 = 1,nz
			        cdf_pr_z(ee0,ee1) = sum(pr_z(ee0,1:ee1))
			    ENDDO
			ENDDO
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

		OPEN   (UNIT=12, FILE=trim(Result_Folder)//'agrid', STATUS='replace')
		WRITE(unit=12, FMT=*) agrid
		CLOSE (unit=12)

		OPEN   (UNIT=12, FILE=trim(Result_Folder)//'zgrid', STATUS='replace')
		WRITE(unit=12, FMT=*) zgrid
		CLOSE (unit=12)

		OPEN   (UNIT=12, FILE=trim(Result_Folder)//'lambdagrid', STATUS='replace')
		WRITE(unit=12, FMT=*) lambdagrid
		CLOSE (unit=12)

		OPEN   (UNIT=12, FILE=trim(Result_Folder)//'egrid', STATUS='replace')
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
		DO age=1,MaxAge
		DO zi=1,nz
		DO lambdai=1,nlambda
		DO ei=1, ne
		      DBN1(age,1,zi,lambdai,ei) = (pop(age)/sum(pop))*Gz(zi)*Glambda(lambdai)*Ge_byage(age,ei)      
		ENDDO
		ENDDO
		ENDDO
		ENDDO  

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
	integer  :: age, ai, zi, lambdai, ei, ei_p, tklo, tkhi
	real(dp) :: dV_low, dV_high, spline_coeff(na), V_spline_R, sp_coeff_W(na,ne), V_spline_W(ne), Prob_lo, Prob_hi

	!$ call omp_set_num_threads(5)

	K_mat  = K_Matrix(R,P)
	Pr_mat = Profit_Matrix(R,P)

	! Final period of life
	age = MaxAge 
	do zi=1,nz
	do ai=1,na 
	do ei=1,ne 
	do lambdai=1,nlambda
		V_Pr(age,ai,zi,lambdai,ei) = Pr_mat(ai,zi)
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
		!$omp parallel do private(lambdai,ei,ai,spline_coeff,V_spline_R,tklo,tkhi,Prob_lo,Prob_hi)
		do ei=1,ne 
		do lambdai=1,nlambda

			!CALL spline( agrid, V_Pr(age+1, :, zi, lambdai, ei) , na , dV_low , dV_high , spline_coeff)  
		
			do ai=1,na
				! call splint( agrid, V_Pr(age+1, :, zi, lambdai, ei), spline_coeff, na, Aprime(age,ai,zi,lambdai, ei), V_spline_R ) 
				! V_Pr(age,ai,zi,lambdai,ei) = Pr_mat(ai,zi) + survP(age)/(1.0_dp+R) * V_spline_R

				! Prepare linear interpolation
				if ( Aprime(age,ai,zi,lambdai, ei) .ge. amax) then
					tklo =na-1
				elseif (Aprime(age,ai,zi,lambdai, ei) .lt. amin) then
				    tklo = 1
				else
				    tklo = ((Aprime(age,ai,zi,lambdai, ei) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
				endif            
				tkhi = tklo + 1        
				Prob_lo = ( agrid(tkhi) - Aprime(age,ai,zi,lambdai, ei) ) / ( agrid(tkhi) -agrid(tklo) )
				Prob_hi = ( Aprime(age,ai,zi,lambdai, ei) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )        
				Prob_lo = min(Prob_lo, 1.0_DP)
				Prob_lo = max(Prob_lo, 0.0_DP)
				Prob_hi = min(Prob_hi, 1.0_DP)
				Prob_hi = max(Prob_hi, 0.0_DP)    

				V_Pr(age,ai,zi,lambdai,ei) = Pr_mat(ai,zi) + survP(age)/(1.0_dp+R) &
			  		&  * ( Prob_lo* V_Pr(age+1, tklo, zi, lambdai, ei)  + Prob_hi* V_Pr(age+1, tkhi, zi, lambdai, ei) ) 

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
		!$omp parallel do private(lambdai,ei,ai,sp_coeff_W,V_spline_W,ei_p,tklo,tkhi,Prob_lo,Prob_hi)
		do lambdai=1,nlambda
			! do ei_p=1,ne
			! 	CALL spline( agrid, V_Pr(age+1, :, zi, lambdai, ei_p) , na , dV_low , dV_high , sp_coeff_W(:,ei_p))
			! enddo   
		
		do ei=1,ne 
		do ai=1,na
			do ei_p=1,ne
				! call splint( agrid, V_Pr(age+1, :, zi, lambdai, ei_p), sp_coeff_W, na, Aprime(age,ai,zi,lambdai, ei), V_spline_W(ei_p) )
				if ( Aprime(age,ai,zi,lambdai, ei_p) .ge. amax) then
					tklo =na-1
				elseif (Aprime(age,ai,zi,lambdai, ei_p) .lt. amin) then
				    tklo = 1
				else
				    tklo = ((Aprime(age,ai,zi,lambdai, ei_p) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
				endif            
				tkhi = tklo + 1        
				Prob_lo = ( agrid(tkhi) - Aprime(age,ai,zi,lambdai, ei_p) ) / ( agrid(tkhi) -agrid(tklo) )
				Prob_hi = ( Aprime(age,ai,zi,lambdai, ei_p) - agrid(tklo) ) / ( agrid(tkhi) -agrid(tklo) )        
				Prob_lo = min(Prob_lo, 1.0_DP)
				Prob_lo = max(Prob_lo, 0.0_DP)
				Prob_hi = min(Prob_hi, 1.0_DP)
				Prob_hi = max(Prob_hi, 0.0_DP) 

				V_spline_W(ei_p) = Prob_lo* V_Pr(age+1, tklo, zi, lambdai, ei_p)  + Prob_hi* V_Pr(age+1, tkhi, zi, lambdai, ei_p)   
			enddo

			V_Pr(age,ai,zi,lambdai,ei) = Pr_mat(ai,zi) + survP(age)/(1.0_dp+R) * sum(pr_e(ei,:)*V_spline_W)

		enddo
	enddo
	enddo
	enddo
	enddo


	! Define Firm based Wealth measure
	do ai=1,na
		Firm_Wealth(:,ai,:,:,:) = V_Pr(:,ai,:,:,:) + (1.0_dp+R)*agrid(ai)
	enddo

END SUBROUTINE Firm_Value


!========================================================================================
!========================================================================================
!========================================================================================

SUBROUTINE  SIMULATION(bench_indx)
	use parameters
	use global
	use omp_lib
	IMPLICIT NONE
	integer, intent(in) :: bench_indx
	integer  :: currentzi, currentlambdai, currentei
	REAL(DP) :: tempnoage, tempnoz, tempnolambda, tempnoe, tempno, currenta, currentY
	REAL(DP) :: start_timet, finish_timet
	INTEGER  :: agecounter, agesign, tage, tzi, tlambdai, tei, tklo, tkhi, paneli, simutime
	INTEGER , DIMENSION(MaxAge) :: requirednumberby_age, cdfrequirednumberby_age
	INTEGER , DIMENSION(totpop) :: panelage, panelz, panellambda, panele 
	REAL(DP), DIMENSION(totpop) :: panela, panelPV_a, panelK  

	!$ call omp_set_num_threads(20)

	print*,'SIMULATION STARTED'

    age=1
    requirednumberby_age(age)    = NINT(totpop*pop(age)/sum(pop))
    cdfrequirednumberby_age(age) = requirednumberby_age(age)
    DO age=2,MaxAge
        requirednumberby_age(age)    = NINT(totpop*pop(age)/sum(pop))
        cdfrequirednumberby_age(age) = requirednumberby_age(age) + cdfrequirednumberby_age(age-1)
    ENDDO
    ! If the total number of people are not equal to the total population, then I will add the remainder to the last age
    requirednumberby_age(MaxAge) =  requirednumberby_age(MaxAge)-cdfrequirednumberby_age(MaxAge) + totpop
    cdfrequirednumberby_age(MaxAge) = totpop

	!=====================================================================
	!                     GENERATE   INITIAL   PANEL
	!=====================================================================


	newiseed=-1

	DO paneli=1,totpop

	! AGE
	   tempnoage = ran1(newiseed)
	   age=1
	   DO WHILE (tempnoage*totpop .gt. cdfrequirednumberby_age(age))
	       age=age+1
	   ENDDO

	! Z   
	   tempnoz = ran1(newiseed)
	   zi=1
	   DO WHILE (tempnoz .gt. cdf_Gz(zi))
	       zi=zi+1
	   ENDDO
	 
	! LAMBDA  
	   tempnolambda = ran1(newiseed) 
	   lambdai=1
	   DO WHILE (tempnolambda .gt. cdf_Glambda(lambdai))
	       lambdai=lambdai+1
	   ENDDO

	! E   
	   tempnoe = ran1(newiseed)   
	   ei=1
	   DO WHILE (tempnoe .gt. cdf_Ge_byage(age,ei))
	       ei=ei+1
	   ENDDO

	   panelage(paneli)=age
	   panelz(paneli)=zi
	   panellambda(paneli)=lambdai
	   panele(paneli)=ei
	   
	ENDDO

	! SET INITIAL ASSET DISTRIBUTION
	panela            = 1.0_DP

	!=============================================================================
	!
	! SIMULATE FROM THE SECOND PERIOD SHOCKS AND UPDATE NEW DISTRIBUTIONS
	!
	!=============================================================================

	!call cpu_time(start_timet) 

	DO simutime=1, MaxSimuTime

	   DO paneli=1,totpop
	    
	       	currenta  		= panela(paneli)
	       	age       		= panelage(paneli)
	       	currentzi      	= panelz(paneli)
	       	currentlambdai 	= panellambda(paneli) 
	       	currentei 		= panele(paneli)
	        
			!  COMPUTE NEXT PERIOD'S ASSET
	       	if (age .lt. MaxAge) then
	            ! do linear interpolation here to find aprime. calling the function takes much more time
	            if (currenta .ge. amax) then
                	tklo = na-1
                elseif (currenta .lt. amin) then
                    tklo = 1
                else
                    tklo = ((currenta - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
	            endif 
	            tkhi = tklo + 1        

	            panela(paneli) = ((agrid(tkhi) - currenta)*Aprime(age,tklo,currentzi,currentlambdai,currentei) 		&
	                           	&  + (currenta - agrid(tklo))*Aprime(age,tkhi,currentzi,currentlambdai, currentei)) &
	                            &  / ( agrid(tkhi) - agrid(tklo) )    
	            
	            if (panela(paneli)  .ge. amax) then
	                panela(paneli) = min(panela(paneli), amax) 
	            endif      
	            if (panela(paneli)  .lt. amin) then
	                panela(paneli) = max(panela(paneli), amin) 
	            endif      

	       	endif !age .lt. MaxAge
			!  NEXT PERIOD'S ASSET IS COMPUTED

	             
			! DRAW NEXT PERIOD'S AGE DBN
	      	tempnoage = ran1(newiseed)  
	  
	      	IF (tempnoage .gt. survP(age)) THEN
				panelage(paneli)   = 1
			ELSE
				panelage(paneli)   = age+1
				panelz(paneli)     = currentzi
				panellambda(paneli)= currentlambdai   
	      	ENDIF
	    
			! DRAW Z and LAMBDA DISTRIBUTION FOR ONE-YEAR OLDS
	     	age = panelage(paneli)   
	     	IF (age .eq. 1) THEN    
				! Z      
		       	tempnoz = ran1(newiseed) 
		       	zi=1
				DO WHILE (tempnoz .gt. cdf_pr_z(currentzi,zi))
					zi=zi+1
				ENDDO
	       
				! LAMBDA  
				tempnolambda = ran1(newiseed) 
				lambdai=1
				DO WHILE (tempnolambda .gt. cdf_pr_lambda(currentlambdai,lambdai))
					lambdai=lambdai+1
				ENDDO

				! E       
				currentei = panele(paneli)
				ei = ne/2+1 ! ALL NEWBORN START FROM THE MEDIAN E  but if too many people died and started from median E, draw a new E for them

				panelz(paneli)		= zi    
				panellambda(paneli)	= lambdai
				panele(paneli) 		= ei 

       		else  IF (age .gt. 1) THEN
	            currentei = panele(paneli)   
	            tempno 	  = ran1(newiseed)   
	            ei 		  = 1
	            DO WHILE (tempno .gt. cdf_pr_e(currentei,ei))
	               ei = ei+1
	            ENDDO            
	            panele(paneli)=ei            
	     	ENDIF ! new age==1
	    ENDDO ! paneli
	ENDDO ! simutime

	!$omp parallel do private(currenta,age,currentzi,currentlambdai,currentei,tklo,tkhi)
	DO paneli=1,totpop
	    currenta = panela(paneli)
	    age=panelage(paneli)
	    currentzi = panelz(paneli)
	    currentlambdai = panellambda(paneli) 
	    currentei = panele(paneli)

	    if (currenta .ge. amax) then
	        tklo = na-1
	    elseif (currenta .lt. amin) then
            tklo = 1
        else
            tklo = ((currenta - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
	    endif    
	    tkhi = tklo + 1        

	    panelPV_a(paneli) = (   (agrid(tkhi) - currenta)     * Firm_Wealth(age,tklo,currentzi,currentlambdai, currentei)  &
	                       &  + (currenta - agrid(tklo))     * Firm_Wealth(age,tkhi,currentzi,currentlambdai, currentei)) &
	                       &  / ( agrid(tkhi) - agrid(tklo) ) 

	    panelK(paneli)    = min( theta*currenta , (mu*P*zgrid(currentzi)**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )             
	           
	ENDDO ! paneli


	print*, ' '
	print*, 'Writing simulation results'
	call system( 'mkdir -p ' // trim(Result_Folder) // 'Simul/' )

	if (bench_indx.eq.1) then
		OPEN   (UNIT=10, FILE=trim(Result_Folder)//'Simul/panela_bench'		   	, STATUS='replace')
		OPEN   (UNIT=11, FILE=trim(Result_Folder)//'Simul/panelage_bench'	 	, STATUS='replace')
		OPEN   (UNIT=12, FILE=trim(Result_Folder)//'Simul/panelz_bench'		    , STATUS='replace')
		OPEN   (UNIT=13, FILE=trim(Result_Folder)//'Simul/panellambda_bench'    , STATUS='replace')
		OPEN   (UNIT=14, FILE=trim(Result_Folder)//'Simul/panele_bench'        	, STATUS='replace')
		OPEN   (UNIT=26, FILE=trim(Result_Folder)//'Simul/panelPV_a_bench'      , STATUS='replace')
		OPEN   (UNIT=27, FILE=trim(Result_Folder)//'Simul/panelK_bench'         , STATUS='replace')
	else 
		OPEN   (UNIT=10, FILE=trim(Result_Folder)//'Simul/panela_exp'		 	, STATUS='replace')
		OPEN   (UNIT=11, FILE=trim(Result_Folder)//'Simul/panelage_exp'		    , STATUS='replace')
		OPEN   (UNIT=12, FILE=trim(Result_Folder)//'Simul/panelz_exp'		 	, STATUS='replace')
		OPEN   (UNIT=13, FILE=trim(Result_Folder)//'Simul/panellambda_exp'   	, STATUS='replace')
		OPEN   (UNIT=14, FILE=trim(Result_Folder)//'Simul/panele_exp'        	, STATUS='replace')
		OPEN   (UNIT=26, FILE=trim(Result_Folder)//'Simul/panelPV_a_exo'        , STATUS='replace')
		OPEN   (UNIT=27, FILE=trim(Result_Folder)//'Simul/panelK_exp' 	        , STATUS='replace')
	endif 


	WRITE  (UNIT=10, FMT=*) panela
	WRITE  (UNIT=11, FMT=*) panelage 
	WRITE  (UNIT=12, FMT=*) panelz 
	WRITE  (UNIT=13, FMT=*) panellambda 
	WRITE  (UNIT=14, FMT=*) panele 
	WRITE  (UNIT=26, FMT=*) panelPV_a
	WRITE  (UNIT=27, FMT=*) panelK

	close (unit=10)
	close (unit=11)
	close (unit=12)
	close (unit=13)
	close (unit=14)
	close (unit=26)
	close (unit=27)


END SUBROUTINE SIMULATION


!========================================================================================
!========================================================================================
!========================================================================================

SUBROUTINE  SIMULATION_TOP(bench_indx)
	use parameters
	use global
	use omp_lib
	IMPLICIT NONE
	integer, intent(in) :: bench_indx
	integer  :: currentzi, currentlambdai, currentei
	REAL(DP) :: tempnoage, tempnoz, tempnolambda, tempnoe, tempno, currenta, currentY
	REAL(DP) :: start_timet, finish_timet
	INTEGER  :: agecounter, agesign, tage, tzi, tlambdai, tei, tklo, tkhi, paneli, simutime
	INTEGER , DIMENSION(MaxAge) :: requirednumberby_age, cdfrequirednumberby_age
	INTEGER , DIMENSION(totpop) :: panelage, panelz, panellambda, panele 
	REAL(DP), DIMENSION(totpop) :: panela, panelPV_a, panelK 
	INTEGER , DIMENSION(150,80) :: panelage_top, panelz_top
	REAL(DP), DIMENSION(150,80) :: panela_top, panelK_top 
	INTEGER 				    :: top_ind(80)

	!$ call omp_set_num_threads(20)

	! top_ind = (/ 12989503 , 16359499 , 17210684 , 2112354  , 2286434  , 16526136 , 7438210  , 6869588  , 19078203 , 15606788 ,&
	! 		&    2225856  , 12531737 , 11062763 , 19793497 , 16188809 , 18407385 , 11304514 , 6412028  , 1148323 , 8595835 ,&
	! 		&    14145270 , 18478397 , 2737060  , 1032430  , 7474392  , 12080222 , 7661199  , 8814829  , 7723389 , 168007 ,&
	! 		&    17158184 , 9059254  , 4597720  , 11328153 , 9623802  , 877507   , 9864235  , 5559945  , 16345916 , 7342133 ,&
	! 		&    16218846 , 13141854 , 9039802  , 16373055 , 14633731 , 3949499  , 8253465  , 48763    , 18557196 , 19577804 ,&
	! 		&    19590912 , 6845111  , 10804332 , 10151506 , 3053807  , 18567572 , 6728883  , 13320399 , 2160040 , 7239959 ,&
	! 		&    366537   , 15071340 , 1194864  , 11098872 , 12735017 , 2876272  , 5195770  , 4720772  , 12421325 , 12027769 ,&
	! 		&    6114885  , 14643750 , 11624454 , 13317533 , 10460044 , 9292769  , 177466   , 9419824  , 17555631 , 16615490/)

	!!!!! mu=0.95
	top_in = (/ 10563514 , 15606788 , 13548896 , 16730950 , 16766492 , 11328153 , 18525999 , 12531737 , &
			&	19793497 , 7438210 , 5976746 , 16526136 , 6412028 , 18478397 , 14145270 , 9623802 , 2286434 , &
			& 	16359499 , 12080222 , 2737060 , 19078203 , 16188809 , 877507 , 2225856 , 16218846 , 11062763 , &
			&	8595835 , 7723389 , 9864235 , 7474392 , 4597720 , 18407385 , 8814829 , 14633731 , 17158184 , &
			&   16345916 , 19577804 , 18557196 , 1032430 , 16373055 , 7661199 , 9059254 , 13141854 , 168007 , &
			&   9039802 , 5559945 , 7342133 , 10151506 , 6728883 , 8253465 , 3053807 , 6845111 , 3949499 , &
			&	15071340 , 5195770 , 48763 , 10804332 , 19590912 , 2160040 , 7239959 , 13320399 , 11098872 , &
			&	12735017 , 18567572 , 12027769 , 1194864 , 366537 , 4720772 , 11624454 , 13317533 , 2876272 , &
			& 	6114885 , 12421325 , 14643750 , 9292769 , 10460044 , 9419824 , 177466 , 17555631 , 16615490 /)

	print*,'SIMULATION STARTED (for top agents)'

    age=1
    requirednumberby_age(age)    = NINT(totpop*pop(age)/sum(pop))
    cdfrequirednumberby_age(age) = requirednumberby_age(age)
    DO age=2,MaxAge
        requirednumberby_age(age)    = NINT(totpop*pop(age)/sum(pop))
        cdfrequirednumberby_age(age) = requirednumberby_age(age) + cdfrequirednumberby_age(age-1)
    ENDDO
    ! If the total number of people are not equal to the total population, then I will add the remainder to the last age
    requirednumberby_age(MaxAge) =  requirednumberby_age(MaxAge)-cdfrequirednumberby_age(MaxAge) + totpop
    cdfrequirednumberby_age(MaxAge) = totpop

	!=====================================================================
	!                     GENERATE   INITIAL   PANEL
	!=====================================================================


	newiseed=-1

	DO paneli=1,totpop

	! AGE
	   tempnoage = ran1(newiseed)
	   age=1
	   DO WHILE (tempnoage*totpop .gt. cdfrequirednumberby_age(age))
	       age=age+1
	   ENDDO

	! Z   
	   tempnoz = ran1(newiseed)
	   zi=1
	   DO WHILE (tempnoz .gt. cdf_Gz(zi))
	       zi=zi+1
	   ENDDO
	 
	! LAMBDA  
	   tempnolambda = ran1(newiseed) 
	   lambdai=1
	   DO WHILE (tempnolambda .gt. cdf_Glambda(lambdai))
	       lambdai=lambdai+1
	   ENDDO

	! E   
	   tempnoe = ran1(newiseed)   
	   ei=1
	   DO WHILE (tempnoe .gt. cdf_Ge_byage(age,ei))
	       ei=ei+1
	   ENDDO

	   panelage(paneli)=age
	   panelz(paneli)=zi
	   panellambda(paneli)=lambdai
	   panele(paneli)=ei
	   
	ENDDO

	! SET INITIAL ASSET DISTRIBUTION
	panela            = 1.0_DP

	!=============================================================================
	!
	! SIMULATE FROM THE SECOND PERIOD SHOCKS AND UPDATE NEW DISTRIBUTIONS
	!
	!=============================================================================

	!call cpu_time(start_timet) 

	DO simutime=1, MaxSimuTime
	    DO paneli=1,totpop
	    
	       	currenta  		= panela(paneli)
	       	age       		= panelage(paneli)
	       	currentzi      	= panelz(paneli)
	       	currentlambdai 	= panellambda(paneli) 
	       	currentei 		= panele(paneli)
	        
			!  COMPUTE NEXT PERIOD'S ASSET
	       	if (age .lt. MaxAge) then
	            ! do linear interpolation here to find aprime. calling the function takes much more time
	            if (currenta .ge. amax) then
                	tklo = na-1
                elseif (currenta .lt. amin) then
                    tklo = 1
                else
                    tklo = ((currenta - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
	            endif 
	            tkhi = tklo + 1        

	            panela(paneli) = ((agrid(tkhi) - currenta)*Aprime(age,tklo,currentzi,currentlambdai,currentei) 		&
	                           	&  + (currenta - agrid(tklo))*Aprime(age,tkhi,currentzi,currentlambdai, currentei)) &
	                            &  / ( agrid(tkhi) - agrid(tklo) )    
	            
	            if (panela(paneli)  .ge. amax) then
	                panela(paneli) = min(panela(paneli), amax) 
	            endif      
	            if (panela(paneli)  .lt. amin) then
	                panela(paneli) = max(panela(paneli), amin) 
	            endif      

	       	endif !age .lt. MaxAge
			!  NEXT PERIOD'S ASSET IS COMPUTED

	             
			! DRAW NEXT PERIOD'S AGE DBN
	      	tempnoage = ran1(newiseed)  
	  
	      	IF (tempnoage .gt. survP(age)) THEN
				panelage(paneli)   = 1
			ELSE
				panelage(paneli)   = age+1
				panelz(paneli)     = currentzi
				panellambda(paneli)= currentlambdai   
	      	ENDIF
	    
			! DRAW Z and LAMBDA DISTRIBUTION FOR ONE-YEAR OLDS
	     	age = panelage(paneli)   
	     	IF (age .eq. 1) THEN    
				! Z      
		       	tempnoz = ran1(newiseed) 
		       	zi=1
				DO WHILE (tempnoz .gt. cdf_pr_z(currentzi,zi))
					zi=zi+1
				ENDDO
	       
				! LAMBDA  
				tempnolambda = ran1(newiseed) 
				lambdai=1
				DO WHILE (tempnolambda .gt. cdf_pr_lambda(currentlambdai,lambdai))
					lambdai=lambdai+1
				ENDDO

				! E       
				currentei = panele(paneli)
				ei = ne/2+1 ! ALL NEWBORN START FROM THE MEDIAN E  but if too many people died and started from median E, draw a new E for them

				panelz(paneli)		= zi    
				panellambda(paneli)	= lambdai
				panele(paneli) 		= ei 

       		else  IF (age .gt. 1) THEN
	            currentei = panele(paneli)   
	            tempno 	  = ran1(newiseed)   
	            ei 		  = 1
	            DO WHILE (tempno .gt. cdf_pr_e(currentei,ei))
	               ei = ei+1
	            ENDDO            
	            panele(paneli)=ei            
	     	ENDIF ! new age==1

	    ENDDO ! paneli

	    ! Record data of top agents
     	if (simutime.ge.(MaxSimuTime-149)) then 
     		print*, "Selecting top agents - Period", simutime-MaxSimuTime+150
     		panelage_top(simutime-MaxSimuTime+150,:) = panelage(top_ind)
     		panelz_top(simutime-MaxSimuTime+150,:)   = panelz(top_ind)
     		panela_top(simutime-MaxSimuTime+150,:)   = panela(top_ind)
     		panelk_top(simutime-MaxSimuTime+150,:)   = min( theta*panela_top(simutime-MaxSimuTime+150,:) ,&
     					& (mu*P*zgrid(panelz_top(simutime-MaxSimuTime+150,:))**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )
     	endif
	    print*, "Simulation period", simutime
	ENDDO ! simutime

	print*, ' '
	print*, 'Writing simulation results for top agents'
	call system( 'mkdir -p ' // trim(Result_Folder) // 'Simul/' )

	if (bench_indx.eq.1) then
		OPEN   (UNIT=10, FILE=trim(Result_Folder)//'Simul/panela_top_bench'		   	, STATUS='replace')
		OPEN   (UNIT=11, FILE=trim(Result_Folder)//'Simul/panelage_top_bench'	 	, STATUS='replace')
		OPEN   (UNIT=12, FILE=trim(Result_Folder)//'Simul/panelz_top_bench'		    , STATUS='replace')
		OPEN   (UNIT=27, FILE=trim(Result_Folder)//'Simul/panelK_top_bench'         , STATUS='replace')
	else 
		OPEN   (UNIT=10, FILE=trim(Result_Folder)//'Simul/panela_top_exp'		 	, STATUS='replace')
		OPEN   (UNIT=11, FILE=trim(Result_Folder)//'Simul/panelage_top_exp'		    , STATUS='replace')
		OPEN   (UNIT=12, FILE=trim(Result_Folder)//'Simul/panelz_top_exp'		 	, STATUS='replace')
		OPEN   (UNIT=27, FILE=trim(Result_Folder)//'Simul/panelK_top_exp' 	        , STATUS='replace')
	endif 


	WRITE  (UNIT=10, FMT=*) panela_top
	WRITE  (UNIT=11, FMT=*) panelage_top  
	WRITE  (UNIT=12, FMT=*) panelz_top 
	WRITE  (UNIT=27, FMT=*) panelK_top

	close (unit=10)
	close (unit=11)
	close (unit=12)
	close (unit=27)


END SUBROUTINE SIMULATION_TOP


!========================================================================================
!========================================================================================
!========================================================================================



! SUBROUTINE  SIMULATION_STATS(bench_indx)
! 	use parameters
! 	use global
! 	use omp_lib

! 	IMPLICIT NONE
! 	integer, intent(in) :: bench_indx
! 	integer  :: currentzi, currentlambdai, currentei
! 	REAL(DP) :: tempnoage, tempnoz, tempnolambda, tempnoe, tempno, currenta, currentY, K
! 	REAL(DP) :: start_timet, finish_timet
! 	INTEGER  :: agecounter, agesign, tage, tzi, tlambdai, tei, tklo, tkhi, paneli, simutime
! 	INTEGER,  DIMENSION(MaxAge) :: requirednumberby_age, cdfrequirednumberby_age
! 	INTEGER,  DIMENSION(totpop) :: panelage , panelz , panellambda, panele,   newpanelage , newpanelz , newpanellambda, newpanele
! 	REAL(DP), DIMENSION(totpop) :: panela,  newpanela,  panel_return, panelcons, panelhours, panelaprime, panel_at_return
! 	REAL(DP), DIMENSION(totpop) :: panel_firm_wealth
! 	INTEGER,  DIMENSION(totpop) :: eligible
! 	REAL(DP), DIMENSION(totpop) :: panela_old_1, panela_old_2, panela_old_3, panela_new_1, panela_new_2, panela_new_3 
! 	Real(DP), allocatable       :: eligible_panela_old_1(:), eligible_panela_old_2(:), eligible_panela_old_3(:)
! 	Real(DP), allocatable       :: eligible_panela_new_1(:), eligible_panela_new_2(:), eligible_panela_new_3(:)
! 	INTEGER                     :: n_eligible

! 	!$ call omp_set_num_threads(20)

!     age=1
!     requirednumberby_age(age)    = NINT(totpop*pop(age)/sum(pop))
!     cdfrequirednumberby_age(age) = requirednumberby_age(age)
!     DO age=2,MaxAge
!         requirednumberby_age(age)    = NINT(totpop*pop(age)/sum(pop))
!         cdfrequirednumberby_age(age) = requirednumberby_age(age) + cdfrequirednumberby_age(age-1)
!     ENDDO
!     ! If the total number of people are not equal to the total population, then I will add the remainder to the last age
!     requirednumberby_age(MaxAge)    = requirednumberby_age(MaxAge)-cdfrequirednumberby_age(MaxAge) + totpop
!     cdfrequirednumberby_age(MaxAge) = totpop

! 	!=====================================================================
! 	!                     GENERATE   INITIAL   PANEL
! 	!=====================================================================

	
! 		newiseed=-1

! 		!numberby_age_z_lambda=0
! 		!numberby_age_e =0

! 		! !$omp parallel do private(age,zi,lambdai,ei,tempnoage,tempnoz,tempnolambda,tempnoe)
! 		DO paneli=1,totpop

! 			! AGE
! 		   	tempnoage = ran1(newiseed)
! 		   	age=1
! 		   	DO WHILE (tempnoage*totpop .gt. cdfrequirednumberby_age(age))
! 		    	age=age+1
! 		   	ENDDO

! 			! Z   
! 		   	tempnoz = ran1(newiseed)
! 		   	zi=1
! 		   	DO WHILE (tempnoz .gt. cdf_Gz(zi))
! 		    	zi=zi+1
! 		   	ENDDO
		 
! 			! LAMBDA  
! 		   	tempnolambda = ran1(newiseed) 
! 		   	lambdai=1
! 		   	DO WHILE (tempnolambda .gt. cdf_Glambda(lambdai))
! 		    	lambdai=lambdai+1
! 		   	ENDDO

! 			! E   
! 		   	tempnoe = ran1(newiseed)   
! 		   	ei=1
! 		   	DO WHILE (tempnoe .gt. cdf_Ge_byage(age,ei))
! 		    	ei=ei+1
! 		   	ENDDO

! 			! CORRECT THE NUMBER OF PEOPLE IF THERE ARE EXCESS
! 				!
! 				!   if (age .gt. 1) then
! 				!        if ( (cdfrequirednumberby_age(age)-tempnoage*totpop) .gt.  (tempnoage*totpop-cdfrequirednumberby_age(age-1)) ) then
! 				!            agesign=1
! 				!            else
! 				!                 agesign=-1
! 				!        endif      
! 				!    else
! 				!        agesign=1
! 				!    endif 
! 				!   agecounter=1
! 				!   tage=age        
! 				!111 IF (sum(numberby_age_z_lambda(age,:,:)) .ge. sum(requirednumberby_age_z_lambda(age,:,:))) then
! 				!            age = tage + agecounter * agesign
! 				!            age = max(age,1)
! 				!            age = min(age,MaxAge)
! 				!            if (sum(numberby_age_z_lambda(age,:,:)) .ge. sum(requirednumberby_age_z_lambda(age,:,:))) then
! 				!                age = age - agecounter * agesign
! 				!                age = max(age,1)
! 				!                age = min(age,MaxAge)
! 				!                if (sum(numberby_age_z_lambda(age,:,:)) .ge. sum(requirednumberby_age_z_lambda(age,:,:))) then
! 				!                    agecounter = agecounter +1
! 				!                    go to 111
! 				!                endif    
! 				!            endif
! 				!       ENDIF
! 				!   
! 				!   if (zi .gt. 1) then 
! 				!       if ( (cdf_Gz(zi) -tempnoz) .gt.  (tempnoz-cdf_Gz(zi-1)) )    then
! 				!           agesign=1
! 				!           else
! 				!                agesign=-1
! 				!       endif      
! 				!   else
! 				!       agesign=1
! 				!   endif
! 				!   agecounter=1  
! 				!   tzi=zi       
! 				!112 IF (sum(numberby_age_z_lambda(age,zi,:)) .ge. sum(requirednumberby_age_z_lambda(age,zi,:))) then
! 				!           zi = tzi + agecounter * agesign
! 				!           zi = max(zi,1)
! 				!           zi=min(zi,nz)
! 				!           IF (sum(numberby_age_z_lambda(age,zi,:)) .ge. sum(requirednumberby_age_z_lambda(age,zi,:))) then
! 				!               zi = zi - agecounter * agesign
! 				!               zi = max(zi,1)
! 				!               zi=min(zi,nz)
! 				!               IF (sum(numberby_age_z_lambda(age,zi,:)) .ge. sum(requirednumberby_age_z_lambda(age,zi,:))) then
! 				!                   agecounter = agecounter +1
! 				!                   go to 112
! 				!               ENDIF
! 				!           ENDIF               
! 				!       ENDIF    
! 				! 
! 				!   if (lambdai .gt. 1) then      
! 				!       if ( (cdf_Glambda(lambdai) -tempnolambda) .gt.  (tempnolambda-cdf_Glambda(lambdai-1)) )    then
! 				!           agesign=1
! 				!           else
! 				!                agesign=-1
! 				!       endif  
! 				!   else
! 				!       agesign=1
! 				!   endif 
! 				!    
! 				!   agecounter=1  
! 				!   tlambdai=lambdai  
! 				!113 IF  (numberby_age_z_lambda(age,zi,lambdai) .ge. requirednumberby_age_z_lambda(age,zi,lambdai)) then
! 				!           lambdai = tlambdai + agecounter * agesign
! 				!           lambdai = max(lambdai,1)
! 				!           lambdai=min(lambdai,nlambda)
! 				!           IF  (numberby_age_z_lambda(age,zi,lambdai) .ge. requirednumberby_age_z_lambda(age,zi,lambdai)) then
! 				!               lambdai = lambdai - agecounter * agesign
! 				!               lambdai = max(lambdai,1)
! 				!               lambdai=min(lambdai,nlambda)
! 				!               IF  (numberby_age_z_lambda(age,zi,lambdai) .ge. requirednumberby_age_z_lambda(age,zi,lambdai)) then
! 				!                   agecounter = agecounter +1
! 				!                   go to 113
! 				!               ENDIF
! 				!           ENDIF               
! 				!       ENDIF
! 				!
! 				!   if (ei .gt. 1) then
! 				!       if ( (Ge_byage(age,ei) -tempnoe) .gt.  (tempnolambda-Ge_byage(age,ei-1) ) )    then
! 				!           agesign=1
! 				!           else
! 				!                agesign=-1
! 				!       endif     
! 				!    else
! 				!        agesign=1
! 				!    endif 
! 				!      
! 				!   agecounter=1  
! 				!   tei=ei      
! 				!114  IF (numberby_age_e(age,ei) .ge. requirednumberby_age_e(age,ei)) THEN
! 				!           ei = tei + agecounter * agesign
! 				!           ei = max(ei,1)
! 				!           ei=min(ei,ne)
! 				!           IF (numberby_age_e(age,ei) .ge. requirednumberby_age_e(age,ei)) THEN
! 				!               ei = tei -  agecounter * agesign
! 				!               ei = max(ei,1)
! 				!               ei=min(ei,ne)
! 				!               IF (numberby_age_e(age,ei) .ge. requirednumberby_age_e(age,ei)) THEN
! 				!                   agecounter = agecounter +1
! 				!                   go to 114
! 				!              ENDIF
! 				!           ENDIF        
! 				!        ENDIF
! 				!   numberby_age_e(age,ei) = numberby_age_e(age,ei)+1    
! 				!   numberby_age_z_lambda(age,zi,lambdai) = numberby_age_z_lambda(age,zi,lambdai)+1
! 				!    
! 			! CORRECTION ENDED	

! 		   panelage(paneli)    = age
! 		   panelz(paneli)      = zi
! 		   panellambda(paneli) = lambdai
! 		   panele(paneli)      = ei
		   
! 		ENDDO
	
! 		!print '("INITIAL = ",f6.3," seconds.")',finish_timet-start_timet

! 		newpanelage    = panelage
! 		newpanelz      = panelz
! 		newpanele      = panele
! 		newpanellambda = panellambda

! 		! SET INITIAL ASSET DISTRIBUTION
! 		panela     = 1.0_DP
! 		newpanela  = 1.0_DP


! 		! Set default value for eligibility indicator to one
! 		eligible = 1 
	

! 	!=============================================================================
! 	!
! 	! SIMULATE FROM THE SECOND PERIOD SHOCKS AND UPDATE NEW DISTRIBUTIONS
! 	!
! 	!=============================================================================

! 		!call cpu_time(start_timet) 


! 		DO simutime=1, MaxSimuTime

! 			panelage = newpanelage
! 			panelz   = newpanelz
! 			panele   = newpanele
! 			panellambda = newpanellambda
! 			panela   = newpanela

! 			!print*,'simutime=',simutime

! 			!numberby_age=0
! 			!deathby_age=0
! 			!survivingby_age =0
! 			!
! 			!deathby_age_z_lambda = 0
! 			!survivingby_age_z_lambda=0
! 			!numberby_age_z_lambda=0
! 			!numberby_e_e=0

! 			newpanela = amin

! 			! !$omp parallel do &
! 			! !$omp private(age,zi,ei,lambdai,currenta,currentzi,currentlambdai,currentei,tklo,tkhi,tempnoage,tempnoz,tempnolambda,tempnoe)
! 			DO paneli=1,totpop
			    
! 				currenta  = panela(paneli)
! 				age       = panelage(paneli)
! 				currentzi = panelz(paneli)
! 				currentlambdai = panellambda(paneli) 
! 				currentei = panele(paneli)
			       
! 				! currentY=  Y_a(currenta,zgrid(currentzi)) 
			       
! 				!  COMPUTE NEXT PERIOD'S ASSET
! 			    if (age .lt. MaxAge) then
! 					! newpanela(paneli) = Linear_Int(agrid, Aprime(age,:,currentzi,currentlambdai, currentei),na,currenta)
! 					! newpanela(paneli) = Linear_Int(YGRID(:,currentzi), Aprime(age,:,currentzi,currentlambdai, currentei),na,currentY)
! 					! newpanela(paneli) = Linear_Int_Aprime(agrid, Aprime(age,:,currentzi,currentlambdai, currentei),na,currenta)

! 			        ! do linear interpolation here to find aprime. calling the function takes much more time
! 		            if (currenta .ge. amax) then
! 		                tklo = na-1
! 		            else if (currenta .lt. amin) then
! 		                tklo = 1
! 		            else
! 		                tklo = ((currenta - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
! 		            endif 
			            
! 			        tkhi = tklo + 1        

! 			        newpanela(paneli) = ((agrid(tkhi) - currenta)*Aprime(age,tklo,currentzi,currentlambdai, currentei) &
! 			                            &  + (currenta - agrid(tklo))*Aprime(age,tkhi,currentzi,currentlambdai, currentei)) &
! 			                            &  / ( agrid(tkhi) - agrid(tklo) )    
			            
! 		            if (newpanela(paneli)  .ge. amax) then
! 		                newpanela(paneli) = min(newpanela(paneli), amax) 
! 		            endif      
! 		            if (newpanela(paneli)  .lt. amin) then
! 		                newpanela(paneli) = max(newpanela(paneli), amin) 
! 		            endif      

! 				endif !age .lt. MaxAge
! 				!  NEXT PERIOD'S ASSET IS COMPUTED

			             
! 				! DRAW NEXT PERIOD'S AGE DBN
! 			    tempnoage = ran1(newiseed)  
			  
! 			    IF (tempnoage .gt. survP(age)) THEN
! 			        newpanelage(paneli) = 1
! 			    ELSE
! 			        newpanelage(paneli)    = age+1
! 			        newpanelz(paneli)      = currentzi
! 			        newpanellambda(paneli) = currentlambdai   
! 			    ENDIF
			  
! 				! CORRECT AGE-Z-LAMBDA DEATH DISTRIBUTION 
! 					!     
! 					!       IF (tempnoage .gt. survP(age)) THEN
! 					!           IF (deathby_age_z_lambda(age,currentzi,currentlambdai) .lt. &
! 					!                                    & requireddeathby_age_z_lambda(age,currentzi,currentlambdai)) THEN
! 					!                newpanelage(paneli)=1
! 					!                deathby_age_z_lambda(age,currentzi,currentlambdai) = deathby_age_z_lambda(age,currentzi,currentlambdai)+1
! 					!                ELSE
! 					!                       newpanelage(paneli) =age+1
! 					!                       survivingby_age_z_lambda(age,currentzi,currentlambdai) = &
! 					!                                    &  survivingby_age_z_lambda(age,currentzi,currentlambdai)+1
! 					!                       newpanelz(paneli)=currentzi
! 					!                       newpanellambda(paneli)=currentlambdai
! 					!                       numberby_age_z_lambda(age+1,currentzi,currentlambdai) = & 
! 					!                                    & numberby_age_z_lambda(age+1,currentzi,currentlambdai) +1
! 					!           ENDIF  
! 					!      ELSE
! 					!          IF (survivingby_age_z_lambda(age,currentzi,currentlambdai) .lt. &
! 					!                                    & requiredsurvivingby_age_z_lambda(age,currentzi,currentlambdai)) THEN
! 					!              newpanelage(paneli)=age+1
! 					!              survivingby_age_z_lambda(age,currentzi,currentlambdai) = &
! 					!                                    & survivingby_age_z_lambda(age,currentzi,currentlambdai)+1
! 					!              newpanelz(paneli)=currentzi
! 					!              newpanellambda(paneli)=currentlambdai
! 					!              numberby_age_z_lambda(age+1,currentzi,currentlambdai) = &
! 					!                                    & numberby_age_z_lambda(age+1,currentzi,currentlambdai) +1  
! 					!              ELSE
! 					!                  newpanelage(paneli)=1
! 					!                  deathby_age_z_lambda(age,currentzi,currentlambdai) = deathby_age_z_lambda(age,currentzi,currentlambdai) +1               
! 					!          ENDIF           
! 					!      ENDIF
! 					!      
! 				! CORRECT AGE-Z-LAMBDA DEATH DISTRIBUTION ENDED
			  

			 
! 				! DRAW Z and LAMBDA DISTRIBUTION FOR ONE-YEAR OLDS
! 			    age = newpanelage(paneli)   
			 
! 			   	IF (age .eq. 1) THEN    

! 					! Z      
! 			       	tempnoz = ran1(newiseed) 
! 			       	zi=1
! 			       	DO WHILE (tempnoz .gt. cdf_pr_z(currentzi,zi))
! 			            zi=zi+1
! 			       	ENDDO
			       
! 					! LAMBDA  
! 			       	tempnolambda = ran1(newiseed) 
! 			       	lambdai=1
! 			       	DO WHILE (tempnolambda .gt. cdf_pr_lambda(currentlambdai,lambdai))
! 			           lambdai=lambdai+1
! 			       	ENDDO

! 					! E       
! 			       	currentei = panele(paneli)
! 			       	ei=ne/2+1 ! ALL NEWBORN START FROM THE MEDIAN E  but if too many people died and started from median E, draw a new E for them

! 					! CORRECT AGE-Z-LAMBDA  DISTRIBUTIONS
! 						!
! 						!       if (zi  .gt. 1) then
! 						!           if ( (cdf_pr_z(currentzi,zi) -tempnoz) .lt.  (tempnoz-cdf_pr_z(currentzi,zi-1) ) )    then
! 						!               agesign=1
! 						!               else
! 						!                    agesign=-1
! 						!           endif
! 						!       else
! 						!            agesign=1          
! 						!       endif
! 						!       agecounter=1  
! 						!       tzi=zi       
! 						!115 IF (sum(numberby_age_z_lambda(age,zi,:)) .ge. sum(requirednumberby_age_z_lambda(age,zi,:))) then
! 						!           zi = tzi + agecounter * agesign
! 						!           zi = max(zi,1)
! 						!           zi=min(zi,nz)
! 						!           IF (sum(numberby_age_z_lambda(age,zi,:)) .ge. sum(requirednumberby_age_z_lambda(age,zi,:))) then
! 						!               zi = zi - agecounter * agesign
! 						!               zi = max(zi,1)
! 						!               zi=min(zi,nz)
! 						!               IF (sum(numberby_age_z_lambda(age,zi,:)) .ge. sum(requirednumberby_age_z_lambda(age,zi,:))) then
! 						!                   agecounter = agecounter +1
! 						!                   go to 115
! 						!               ENDIF
! 						!           ENDIF               
! 						!       ENDIF    
! 						!
! 						!       if (lambdai .gt. 1) then
! 						!            if ( (cdf_pr_lambda(currentlambdai,lambdai) -tempnolambda) .gt.  &
! 						!                            & (tempnolambda-cdf_pr_lambda(currentlambdai,lambdai-1)) )   then
! 						!               agesign=1
! 						!               else
! 						!                    agesign=-1
! 						!           endif      
! 						!       else
! 						!            agesign=1          
! 						!       endif       
! 						!       agecounter=1  
! 						!       tlambdai=lambdai  
! 						!116 IF  (numberby_age_z_lambda(age,zi,lambdai) .ge. requirednumberby_age_z_lambda(age,zi,lambdai)) then
! 						!               lambdai = tlambdai + agecounter * agesign
! 						!               lambdai = max(lambdai,1)
! 						!               lambdai=min(lambdai,nlambda)
! 						!               IF  (numberby_age_z_lambda(age,zi,lambdai) .ge. requirednumberby_age_z_lambda(age,zi,lambdai)) then
! 						!                   lambdai = lambdai - agecounter * agesign
! 						!                   lambdai = max(lambdai,1)
! 						!                   lambdai=min(lambdai,nlambda)
! 						!                   IF  (numberby_age_z_lambda(age,zi,lambdai) .ge. requirednumberby_age_z_lambda(age,zi,lambdai)) then
! 						!                       agecounter = agecounter +1
! 						!                       go to 116
! 						!                   ENDIF
! 						!               ENDIF               
! 						!         ENDIF
! 						!          
! 						!        if ( numberby_e_e(currentei,ei) .ge. requirednumberby_e_e(currentei,ei) ) then
! 						!            tempnoe = ran1(newiseed) 
! 						!            ei=1
! 						!            DO WHILE (tempnoe .gt. cdf_Ge(ei))
! 						!               ei=ei+1
! 						!            ENDDO    
! 						!
! 						!           if (ei .gt. 1) then 
! 						!               if ( (cdf_Ge(ei) -tempnoe) .gt.  (tempnolambda-cdf_Ge(ei-1)) )    then
! 						!                   agesign=1
! 						!                   else
! 						!                        agesign=-1
! 						!               endif      
! 						!           else
! 						!               agesign=1 
! 						!           endif        
! 						!           agecounter=1  
! 						!           tei=ei      
! 						! 117    IF (numberby_e_e(currentei,ei) .ge. requirednumberby_e_e(currentei,ei)) THEN
! 						!                   ei = tei + agecounter * agesign
! 						!                   ei = max(ei,1)
! 						!                   ei=min(ei,ne)
! 						!                   IF (numberby_e_e(currentei,ei) .ge. requirednumberby_e_e(currentei,ei)) THEN
! 						!                       ei = tei -  agecounter * agesign
! 						!                       ei = max(ei,1)
! 						!                       ei=min(ei,ne)
! 						!                       IF (numberby_e_e(currentei,ei) .ge. requirednumberby_e_e(currentei,ei)) THEN
! 						!                           agecounter = agecounter +1
! 						!                           go to 117
! 						!                      ENDIF
! 						!                   ENDIF        
! 						!             ENDIF                
! 						!        endif
! 						!
! 						!       numberby_age_z_lambda(age,zi,lambdai) = numberby_age_z_lambda(age,zi,lambdai) +1        
! 						!       numberby_e_e(currentei, ei)  = numberby_e_e(currentei,ei) +1
! 						!
! 					!  CORRECTING DISTRIBUTIONS ENDED

! 			        newpanelz(paneli)      = zi    
! 			        newpanellambda(paneli) = lambdai
! 			        newpanele(paneli)      = ei  
			        
! 			    ENDIF ! new age==1
			 
! 			    ! DRAW NEW E FOR THOSE WHO ARE NOT NEWBORN
! 			    IF (age .gt. 1) THEN             
! 			        currentei = panele(paneli)   
! 			        tempno = ran1(newiseed)   
! 			        ei=1
! 			        DO WHILE (tempno .gt. cdf_pr_e(currentei,ei))
! 						ei=ei+1
! 			        ENDDO
			            
! 					! CORRECT E DISTRIBUTION 
! 						!
! 						!           if (ei .gt. 1) then 
! 						!                if ( (cdf_pr_e(currentei,ei) -tempnoe) .gt.  (tempnolambda-cdf_pr_e(currentei,ei-1)) )    then
! 						!                   agesign=1
! 						!                   else
! 						!                        agesign=-1
! 						!               endif    
! 						!           else
! 						!               agesign=1 
! 						!           endif  
! 						!           agecounter=1  
! 						!           tei=ei      
! 						! 118    IF (numberby_e_e(currentei,ei) .ge. requirednumberby_e_e(currentei,ei)) THEN
! 						!                   ei = tei + agecounter * agesign
! 						!                   ei = max(ei,1)
! 						!                   ei=min(ei,ne)
! 						!                   IF (numberby_e_e(currentei,ei) .ge. requirednumberby_e_e(currentei,ei)) THEN
! 						!                       ei = tei -  agecounter * agesign
! 						!                       ei = max(ei,1)
! 						!                       ei=min(ei,ne)
! 						!                       IF (numberby_e_e(currentei,ei) .ge. requirednumberby_e_e(currentei,ei)) THEN
! 						!                           agecounter = agecounter +1
! 						!                           go to 118
! 						!                      ENDIF
! 						!                   ENDIF        
! 						!             ENDIF                
! 						!             numberby_e_e(currentei,ei) = numberby_e_e(currentei,ei) + 1
! 						!
! 					! CORRECT E DISTRIBUTION ENDED

! 			        newpanele(paneli)=ei            
! 			    ENDIF ! age .gt. 1        
! 			ENDDO ! paneli


! 			panelage     = newpanelage
! 			panela       = newpanela
! 			panelz       = newpanelz
! 			panellambda  = newpanellambda
! 			panele       = newpanele

! 			! Save data on assets for the last periods
! 			! Agents are eligible if:
! 				! 1) They don't die during the first two recording periods
! 				! 2) They they die in the third recording period
! 				! 3) They don't die again
! 			if (simutime.eq.(MaxSimuTime-15)) then 
! 		    	panela_old_1 = panela
! 		    	where(panelage==1) eligible = 0 
! 	        endif 
! 	        if (simutime.eq.(MaxSimuTime-14)) then 
! 		    	panela_old_2 = panela
! 		    	where(panelage==1) eligible = 0 
! 	        endif 
! 	        if (simutime.eq.(MaxSimuTime-13)) then 
! 		    	panela_old_3 = panela
! 		    	where(panelage==1) eligible = 0 
! 	        endif 
! 	        if (simutime.eq.(MaxSimuTime-12)) then
! 	        	where(panelage>1) eligible = 0
! 	        endif 
! 	        if ((simutime.gt.(MaxSimuTime-12)).and.(simutime.lt.(MaxSimuTime-2))) then
! 	        	where(panelage==1) eligible = 0 
! 	        endif
! 	        if (simutime.eq.(MaxSimuTime-2)) then 
! 		    	panela_new_1 = panela
! 		    	where(panelage==1) eligible = 0 
! 	        endif 
! 	        if (simutime.eq.(MaxSimuTime-1)) then 
! 		    	panela_new_2 = panela
! 		    	where(panelage==1) eligible = 0 
! 	        endif 
! 	        if (simutime.eq.(MaxSimuTime)) then 
! 		    	panela_new_3 = panela
! 		    	where(panelage==1) eligible = 0 
! 	        endif 
! 		ENDDO ! simutime

! 		n_eligible = sum(eligible)

! 		allocate( eligible_panela_old_1(n_eligible), eligible_panela_old_2(n_eligible), eligible_panela_old_3(n_eligible) )
! 		allocate( eligible_panela_new_1(n_eligible), eligible_panela_new_2(n_eligible), eligible_panela_new_3(n_eligible) )

! 		eligible_panela_old_1 = pack(panela_old_1 , (eligible.eq.1) )
! 		eligible_panela_old_2 = pack(panela_old_2 , (eligible.eq.1) )
! 		eligible_panela_old_3 = pack(panela_old_3 , (eligible.eq.1) )
! 		eligible_panela_new_1 = pack(panela_new_1 , (eligible.eq.1) )
! 		eligible_panela_new_2 = pack(panela_new_2 , (eligible.eq.1) )
! 		eligible_panela_new_3 = pack(panela_new_3 , (eligible.eq.1) )

! 		print*, ' '
! 		print*, 'n_eligible', sum(eligible)
! 		print*, ' '


! 	!=============================================================================
! 	!
! 	! Panel on hours, returns, etc for final period
! 	!
! 	!=============================================================================
! 		!$omp parallel do private(currenta,age,currentzi,currentlambdai,currentei,tklo,tkhi,K)
! 		DO paneli=1,totpop

! 			currenta  = panela(paneli)
! 			age       = panelage(paneli)
! 			currentzi = panelz(paneli)
! 			currentlambdai = panellambda(paneli) 
! 			currentei = panele(paneli)
		       

! 	        ! do linear interpolation here to find aprime. calling the function takes much more time
! 	        if (currenta .ge. amax) then
! 	            tklo = na-1
! 	        else if (currenta .lt. amin) then
! 	            tklo = 1
! 	        else
! 	            tklo = ((currenta - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
! 	        endif 
	        
! 	        tkhi = tklo + 1        

! 	        panelcons(paneli) = ((agrid(tkhi) - currenta)*cons(age,tklo,currentzi,currentlambdai, currentei) &
! 	                                &  + (currenta - agrid(tklo))*cons(age,tkhi,currentzi,currentlambdai, currentei)) &
! 	                                &  / ( agrid(tkhi) - agrid(tklo) )              

! 	        panelhours(paneli) = ((agrid(tkhi) - currenta)*hours(age,tklo,currentzi,currentlambdai, currentei) &
! 	                                &  + (currenta - agrid(tklo))*hours(age,tkhi,currentzi,currentlambdai, currentei)) &
! 	                                &  / ( agrid(tkhi) - agrid(tklo) )  

! 	        panelaprime(paneli) =((agrid(tkhi) - currenta)*Aprime(age,tklo,currentzi,currentlambdai, currentei) &
! 	                                &  + (currenta - agrid(tklo))*Aprime(age,tkhi,currentzi,currentlambdai, currentei)) &
! 	                                &  / ( agrid(tkhi) - agrid(tklo) )  
		
! 			K = min( theta*currenta , (mu*P*zgrid(currentzi)**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )

! 			panel_return(paneli)    = R*currenta + ( P*(zgrid(currentzi)*K)**mu - (R+DepRate)*K )

! 			if ( ( currenta + panel_return(paneli)*(1.0_DP-tauK) ).le.Y_a_threshold) then  
! 				panel_at_return(paneli) = ( currenta + panel_return(paneli)*(1.0_DP-tauK) )*(1-tauW_bt) - currenta
! 			else 
! 				panel_at_return(paneli) = ( currenta + panel_return(paneli)*(1.0_DP-tauK) )*(1-tauW_at) - currenta
! 			endif 

! 			panel_firm_wealth(paneli) = (1.0_dp+R)*currenta + &
! 									& ((agrid(tkhi) - currenta)*V_Pr(age,tklo,currentzi,currentlambdai, currentei) &
! 	                                &  + (currenta - agrid(tklo))*V_Pr(age,tkhi,currentzi,currentlambdai, currentei)) &
! 	                                &  / ( agrid(tkhi) - agrid(tklo) )              

! 		ENDDO ! paneli


! 	print*, ' '
! 	print*, 'Writing simulation results'
! 	call system( 'mkdir -p ' // trim(Result_Folder) // 'Simul/' )

! 	if (bench_indx.eq.1) then
! 		OPEN   (UNIT=10, FILE=trim(Result_Folder)//'Simul/panela_bench'		   	  , STATUS='replace')
! 		OPEN   (UNIT=11, FILE=trim(Result_Folder)//'Simul/panelage_bench'	  	  , STATUS='replace')
! 		OPEN   (UNIT=12, FILE=trim(Result_Folder)//'Simul/panelz_bench'		 	  , STATUS='replace')
! 		OPEN   (UNIT=13, FILE=trim(Result_Folder)//'Simul/panellambda_bench'   	  , STATUS='replace')
! 		OPEN   (UNIT=14, FILE=trim(Result_Folder)//'Simul/panele_bench'        	  , STATUS='replace')
! 		OPEN   (UNIT=15, FILE=trim(Result_Folder)//'Simul/panel_return_bench'  	  , STATUS='replace')
! 		OPEN   (UNIT=16, FILE=trim(Result_Folder)//'Simul/panel_cons_bench'		  , STATUS='replace') 
! 		OPEN   (UNIT=17, FILE=trim(Result_Folder)//'Simul/panel_hours_bench'	  , STATUS='replace') 
! 		OPEN   (UNIT=18, FILE=trim(Result_Folder)//'Simul/panel_aprime_bench' 	  , STATUS='replace') 
! 		OPEN   (UNIT=19, FILE=trim(Result_Folder)//'Simul/panel_at_return_bench'  , STATUS='replace')

! 		OPEN   (UNIT=20, FILE=trim(Result_Folder)//'Simul/panela_old_1'           , STATUS='replace')
! 		OPEN   (UNIT=21, FILE=trim(Result_Folder)//'Simul/panela_old_2'           , STATUS='replace')
! 		OPEN   (UNIT=22, FILE=trim(Result_Folder)//'Simul/panela_old_3'           , STATUS='replace')
! 		OPEN   (UNIT=23, FILE=trim(Result_Folder)//'Simul/panela_new_1'           , STATUS='replace')
! 		OPEN   (UNIT=24, FILE=trim(Result_Folder)//'Simul/panela_new_2'           , STATUS='replace')
! 		OPEN   (UNIT=25, FILE=trim(Result_Folder)//'Simul/panela_new_3'           , STATUS='replace')

! 		OPEN   (UNIT=26, FILE=trim(Result_Folder)//'Simul/panel_firm_wealth_bench', STATUS='replace')
! 	else 
! 		OPEN   (UNIT=10, FILE=trim(Result_Folder)//'Simul/panela_exp'		 	, STATUS='replace')
! 		OPEN   (UNIT=11, FILE=trim(Result_Folder)//'Simul/panelage_exp'		    , STATUS='replace')
! 		OPEN   (UNIT=12, FILE=trim(Result_Folder)//'Simul/panelz_exp'		 	, STATUS='replace')
! 		OPEN   (UNIT=13, FILE=trim(Result_Folder)//'Simul/panellambda_exp'   	, STATUS='replace')
! 		OPEN   (UNIT=14, FILE=trim(Result_Folder)//'Simul/panele_exp'        	, STATUS='replace')
! 		OPEN   (UNIT=15, FILE=trim(Result_Folder)//'Simul/panel_return_exp'  	, STATUS='replace')
! 		OPEN   (UNIT=16, FILE=trim(Result_Folder)//'Simul/panel_cons_exp'		, STATUS='replace') 
! 		OPEN   (UNIT=17, FILE=trim(Result_Folder)//'Simul/panel_hours_exp'	 	, STATUS='replace') 
! 		OPEN   (UNIT=18, FILE=trim(Result_Folder)//'Simul/panel_aprime_exp' 	, STATUS='replace') 
! 		OPEN   (UNIT=19, FILE=trim(Result_Folder)//'Simul/panel_at_return_exp'	, STATUS='replace') 

! 		OPEN   (UNIT=26, FILE=trim(Result_Folder)//'Simul/panel_firm_wealth_exp', STATUS='replace')
! 	endif 


! 	WRITE  (UNIT=10, FMT=*) panela
! 	WRITE  (UNIT=11, FMT=*) panelage 
! 	WRITE  (UNIT=12, FMT=*) panelz 
! 	WRITE  (UNIT=13, FMT=*) panellambda 
! 	WRITE  (UNIT=14, FMT=*) panele 
! 	WRITE  (UNIT=15, FMT=*) panel_return 
! 	WRITE  (UNIT=16, FMT=*) panelcons
! 	WRITE  (UNIT=17, FMT=*) panelhours
! 	WRITE  (UNIT=18, FMT=*) panelaprime
! 	WRITE  (UNIT=19, FMT=*) panel_at_return 

! 	WRITE  (UNIT=26, FMT=*) panel_firm_wealth

! 	close (unit=10)
! 	close (unit=11)
! 	close (unit=12)
! 	close (unit=13)
! 	close (unit=14)
! 	close (unit=15)
! 	close (unit=16)
! 	close (unit=17)
! 	close (unit=18)
! 	close (unit=19)

! 	close (unit=26)

! 	if (bench_indx.eq.1) then
! 		WRITE (UNIT=20, FMT=*) eligible_panela_old_1
! 		WRITE (UNIT=21, FMT=*) eligible_panela_old_2
! 		WRITE (UNIT=22, FMT=*) eligible_panela_old_3
! 		WRITE (UNIT=23, FMT=*) eligible_panela_new_1
! 		WRITE (UNIT=24, FMT=*) eligible_panela_new_2
! 		WRITE (UNIT=25, FMT=*) eligible_panela_new_3

! 		close (unit=20)
! 		close (unit=21)
! 		close (unit=22)
! 		close (unit=23)
! 		close (unit=24)
! 		close (unit=25)
! 	endif

! 	print*, 'Averages from simulation'
! 	print*, sum(panela)/totpop, sum(panelage)/totpop, sum(panel_return)/totpop, sum(panelhours)/totpop
! 	print*, 'Mean Firm Wealth:', 47000/EBAR*sum(panel_firm_wealth)/totpop, 47000/EBAR*sum(Firm_Wealth*DBN1) 
! 	print*, 'Number of eligible agents for dynamics', sum(eligible)
! 	print*, ' '

! END SUBROUTINE SIMULATION_STATS


!========================================================================================
!========================================================================================
!========================================================================================


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
			WRITE(UNIT=19, FMT=*) 'Theta',theta
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
			WRITE(UNIT=19, FMT=*) "Bequest_Wealth"	        , Bequest_Wealth
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
			WRITE(UNIT=19, FMT=*) 'PV_Wealth_Top_1%' 		, FW_top_x_share(4)
			WRITE(UNIT=19, FMT=*) 'PV_Wealth_Top_10%'		, FW_top_x_share(3)
			WRITE(UNIT=19, FMT=*) 'PV_Wealth_Top_20%'		, FW_top_x_share(2)
			WRITE(UNIT=19, FMT=*) 'PV_Wealth_Top_40%'		, FW_top_x_share(1)
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) 'Bequest'
			WRITE(UNIT=19, FMT=*) 'Mean_Bequest/Wealth='	, Mean_Bequest/MeanWealth 
			WRITE(UNIT=19, FMT=*) 'Mean_Bequest/PV_Wealth='	, Mean_Bequest/Mean_Firm_Wealth 
			WRITE(UNIT=19, FMT=*) 'BQ_Above_Threshold='		, Threshold_Share_bq
			WRITE(UNIT=19, FMT=*) 'Bequest_Revenue/YBAR='	, 0
			WRITE(UNIT=19, FMT=*) 'BQ_p10/YBAR'				, BQ_top_x(1)/YBAR
			WRITE(UNIT=19, FMT=*) 'BQ_p30/YBAR'				, BQ_top_x(2)/YBAR
			WRITE(UNIT=19, FMT=*) 'BQ_p50/YBAR'				, BQ_top_x(3)/YBAR
			WRITE(UNIT=19, FMT=*) 'BQ_p70/YBAR'				, BQ_top_x(4)/YBAR
			WRITE(UNIT=19, FMT=*) 'BQ_p90/YBAR'				, BQ_top_x(5)/YBAR
			WRITE(UNIT=19, FMT=*) 'BQ_p98/YBAR'				, BQ_top_x(6)/YBAR
			WRITE(UNIT=19, FMT=*) 'BQ_p99/YBAR'				, BQ_top_x(7)/YBAR
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

	write(string_theta,'(f4.2)')  theta

	if ((Progressive_Tax_Switch.eqv..false.).and.(NSU_Switch.eqv..true.)) then 
		bench_folder = './NSU_F_LT_Results/Theta_'//trim(string_theta)//'/Bench_Files/'
	else if ((Progressive_Tax_Switch.eqv..true.).and.(NSU_Switch.eqv..true.)) then 
		bench_folder = './NSU_F_PT_Results/Theta_'//trim(string_theta)//'/Bench_Files/'
	else if ((Progressive_Tax_Switch.eqv..false.).and.(NSU_Switch.eqv..false.)) then 
		bench_folder = './SU_F_LT_Results/Theta_'//trim(string_theta)//'/Bench_Files/'
	else if ((Progressive_Tax_Switch.eqv..true.).and.(NSU_Switch.eqv..false.)) then 
		bench_folder = './SU_F_PT_Results/Theta_'//trim(string_theta)//'/Bench_Files/'
	end if 

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

END SUBROUTINE Write_Experimental_Results


!========================================================================================
!========================================================================================
!========================================================================================

end module programfunctions



