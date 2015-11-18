

! This program solves for the wealth taxation model of 
! Guvenen, Kuruscu and Kambourov

!========================================================================================
!========================================================================================
!========================================================================================


MODULE parameters

    use nrtype
    use nrutil

    ! Huge number
    real(dp), parameter         ::  big_p   = HUGE(1.0_dp)

    ! Switch for computing social security with benchmark earnings
    	! If KeepSSatBench=1 then E_bar is kept at E_bar_bench for experiments
    INTEGER(I4B),  PARAMETER :: KeepSSatBench=1

	! Labor efficiency shocks
		! log(y)=  lambda + kappa + e 
		! lambda: inidividual fixed effect (fixed within generation)
		! kappa: life-cycle component
		! e: indiveidual transitory component (varies across time)
    ! Permanent labor earnings componenet (lambda)
    REAL(DP), PARAMETER      :: rho_lambda=0.5_DP 
    REAL(DP)             	 :: sigma_lambda_eps
    INTEGER(I4B),  PARAMETER :: nlambda=1         ! Number of grid points

    ! Transitory labor earnings componenet (e)
    REAL(DP), PARAMETER  	 :: rho_e=0.9_DP, sigma_e_eps=0.20_DP
    INTEGER(I4B),  PARAMETER :: ne=5              ! Number of grid points

    ! Entrepreneurial ability (z)
    REAL(DP)         	     :: rho_z, sigma_z_eps 
    INTEGER(I4B),  PARAMETER :: nz=3              ! Number of grid points

 

    ! Utility: Discount factor (beta) and Disutility from labor (phi) 
	REAL(DP)                 :: beta, phi
    
	! Production 
		! Final good producer
		REAL(DP), PARAMETER  :: alpha=0.33_DP, Aprod=1.0_DP
		! Intermediate good (or home production)
		REAL(DP), PARAMETER  :: mu=0.9_DP
		! Depreciation rate
		REAL(DP), PARAMETER  :: DepRate=0.0_DP
	

	! Life cycle: retirement age, maximum age
	INTEGER(I4B), PARAMETER  :: MaxAge=40

	! Asset grid: nodes (na,fine_na), min (amin), max (amax), curvature (a_theta) 
	INTEGER(I4B), PARAMETER  :: na=201
	REAL(DP)    , PARAMETER  :: a_theta=4.0_DP , amax=100000.0_DP, amin=0.0001_DP

	

	! Control for updates on stationary distribution
		! Every "update_period" iterations policy functions are updated
	INTEGER(I4B), PARAMETER  :: update_period=5
		! The distribution is iterated until convergence or until "MaxSimuTime" iterations
	INTEGER(I4B),  PARAMETER :: MaxSimuTime=500 

	! Parameters for external functions and procedures
		! Number of std std away from mean in tauchen
		REAL(DP), PARAMETER  :: mtauchen=3.0_DP
		! Tolerance for brent algorithm
		REAL(DP), PARAMETER  :: brent_tol=0.00000001_DP
		

	! Taxes
		! Wealth tax: minimum wealth tax to consider and increments for balancing budget
		REAL(DP), PARAMETER  :: tauWmin_bt=0.02_DP, tauWinc_bt=0.005_DP ! Minimum tax below threshold and increments
		REAL(DP), PARAMETER  :: tauWmin_at=0.02_DP, tauWinc_at=0.005_DP ! Minimum tax above threshold and increments
		INTEGER , PARAMETER  :: Y_a_threshold_pct = 0.0_dp ! Percentile of the distribution of wealth
		! Consumption tax
		REAL(DP), PARAMETER  :: tauC=0.075_DP
		! Labor income tax: This is a progresive tax.
			! 1-psi controls the level of tax, and tauPL controls progressivity
		REAL(DP), PARAMETER  :: tauPL=0.185_DP, psi_PL=0.77_DP  

END MODULE parameters


!========================================================================================
!========================================================================================
!========================================================================================


MODULE global
    USE parameters

    ! "params" determines the economy to be simulated. 
    	! It contains the values of: beta, rho_z, sigma_z_eps, sigma_lambda_eps and phi. In order.
    real(DP) , dimension(5)     :: params

    ! Population size (pop) and survival probability by age (survP)
    REAL(DP), DIMENSION(MaxAge) :: pop, survP          
	
    ! Labor efficiency shocks: 
    	! Lyfe cycle component (by age)
		REAL(DP), DIMENSION(MaxAge) :: kappagrid

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

	! Retirement income 
	REAL(DP), DIMENSION(nlambda,ne) :: phi_lambda_e   ! phi_lambda_e is the income replacement ratio used to compute SS payments
	REAL(DP), DIMENSION(nlambda,ne) :: RetY_lambda_e  ! Retirement income = phi_lambda_e*Ebar, the first 45 years are assigned zero. 

	! Labor efficiency units
		! eff_un(age,lambda,e) = kappa(age)*lambda*e
		! yh(age,lambda,e)     = Wage * eff_un
		REAL(DP), DIMENSION(MaxAge  , nlambda, ne) :: eff_un,  yh

	! Policy function and value function (defined on the exogenous grid)
    REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) :: Cons, Hours, Aprime
    REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) :: ValueFunction
		! Analytical solution for mu=1 for all the lifecycle not just retirement period
    	REAL(DP), DIMENSION(MaxAge, na, nz, nlambda) :: AnRetCons,   AnRetValue, AnRetHours 
	! Policy function and value function (defined on the adjusted grid for breakpoints)
	REAL(DP), DIMENSION(:,:,:,:,:), allocatable :: Cons_t, Hours_t, Aprime_t
 
 	! Aggregate variables
	 	! Benchmark values of Q, N, E, Wage, R, G, Y
	    REAL(DP) :: QBAR_bench, NBAR_bench, Ebar_bench, wage_bench, rr_bench, GBAR_bench, Y_bench
	    ! Experiment values of Q, N, E, Wage, R, G, Y
	    REAL(DP) :: QBAR_exp,   NBAR_exp,   Ebar_exp,   wage_exp,   rr_exp,   GBAR_exp, GBAR_exp_old, Y_exp
	    ! Values for aggregate variables (used when solving a given economy)
	    REAL(DP) :: rr, Ebar , wage, NBAR, QBAR, YBAR, GBAR

	! Asset, resources and marginal benefit of wealth grids
		! Note that these grids will change as we change tax rates and for each intermediate good price
		! For each tax rate we weill have different Y grids
		REAL(DP), DIMENSION(na)      :: agrid
	    REAL(DP), DIMENSION(na,nz)   :: YGRID, MBGRID
	    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: YGRID_t, MBGRID_t
	    REAL(DP), DIMENSION(:),   ALLOCATABLE :: agrid_t
	    INTEGER                      :: na_t
	
	! Values for taxes in benchmark and experiment
    REAL(DP) :: tauk_bench, tauL_bench, tauw_bt_bench, tauw_at_bench, Y_a_threshold_bench 
    REAL(DP) :: tauk_exp,   tauL_exp,   tauw_bt_exp,   tauw_at_exp,   Y_a_threshold_exp

    ! Values for taxes (when solving an economy)
    REAL(DP) :: tauK, tauL
    REAL(DP) :: tauW_bt, tauW_at ! Wealth taxes below threshold and above threshold
    REAL(DP) :: Y_a_threshold = 0.0_dp ! Value of the threshold for change in tauW

    ! Auxiliary variables to find wealth tax that balances the budget in experiment economy
    REAL(DP) :: tauWindx, tauW_low_bt, tauW_up_bt, tauW_low_at, tauW_up_at

	! Counters for the age, and index of lamnbda, z, a and e
    INTEGER :: age, lambdai, zi, ai, ei    

    ! Distribution of population by age, a, z, lambda, e
    REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) ::DBN1, DBN_bench, DBN_exp

    ! Stats and distribution in equilibrium
	    ! Distribution of assets
	    REAL(DP), DIMENSION(na)  :: pr_a_dbn, cdf_a_dbn
	    ! Mass of assets by a_grid (defined as a_grid(i)sum(DBN(:,i,:,:,:)))
	    REAL(DP), DIMENSION(na)  :: tot_a_by_grid, cdf_tot_a_by_grid
	    ! Percentiles of the asset distribution and mass of total assets by them
	    INTEGER , DIMENSION(100) ::  prctile_ai
	    REAL(DP), DIMENSION(100) ::  cdf_tot_a_by_prctile
    	! Other stats
	    REAL(DP) :: pop_25_60 , tothours_25_60, pop_pos_earn_25_60, tot_log_earnings_25_60, mean_log_earnings_25_60 
	    REAL(DP) :: meanhours_25_60, Var_Log_Earnings_25_60, Std_Log_Earnings_25_60, MeanWealth, Wealth_Output
	    REAL(DP) :: prct1_wealth, prct10_wealth, SSE_Moments, Min_SSE_Moments
	
	! Auxiliary variables for evaluating FOC: Consumption and assets
    	REAL(DP) :: consin, ain

    ! Auxiliary variable for labor earnings (it is set to psi=psi_PL)
    	REAL(DP) :: psi

    ! Objective moments (not currently in use)	
    	REAL(DP), DIMENSION(5)  ::  Min_Moments
    
    ! Switch to solve benchmark or experiment
    INTEGER :: solving_bench
    
    ! Seed for random number generator
    INTEGER :: newiseed
   
END MODULE global


!========================================================================================
!========================================================================================
!========================================================================================


Module programfunctions
	use parameters
	use global
	    
	contains

!========================================================================================
!========================================================================================

Subroutine Asset_Grid_Threshold(Y_a_threshold_in,agrid_t,na_t)
	real(dp), intent(in)    :: Y_a_threshold_in
	integer , intent(out)   :: na_t
	real(dp), dimension(:), allocatable, intent(out) :: agrid_t
	real(dp), dimension(nz) :: a_aux
	integer                 :: a_ind
	integer , dimension(:), allocatable :: agrid_t_ind
	real(dp), dimension(:), allocatable :: p
	real(dp)                :: max_wealth

	allocate( p(2) )
	a_ind = 0
	! If the threshold for wealth taxes is positive then agrid is adjusted
	if (Y_a_threshold_in.gt.0.0_dp) then 
		p(1) = Y_a_threshold_in
 		do zi=1,nz 
			! New points are added to agrid if there is an "a" st Y(a,z))=Y_threshold
			max_wealth = (1.0_dp-DepRate)*agrid(na)+rr*(agrid(na)*zgrid(zi))**mu
			if (Y_a_threshold_in.lt.max_wealth) then
				a_ind		 = a_ind + 1 
				p(2)         = zgrid(zi)
				a_aux(a_ind) = zbrent_p(Y_a_res,0.0_dp,agrid(na),brent_tol,p) 
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

	contains 

		function Y_a_res(a_in,p)
			real(dp), intent(in) :: a_in
			real(dp), dimension(:), allocatable, intent(in) :: p
			real(dp) :: Y_a_res, Y_a_th, z_in

			Y_a_th = p(1)
			z_in   = p(2)

			Y_a_res = ( a_in + ( rr * (z_in * a_in )**mu - DepRate*a_in ) ) - Y_a_th

		end function Y_a_res

end Subroutine Asset_Grid_Threshold


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
		real(DP), intent(in) :: a_in, z_in
		real(DP)             :: Y_a, wealth 

		! Before tax wealth
		Y_a = ( a_in + ( rr * (z_in * a_in )**mu - DepRate*a_in ) *(1.0_DP-tauK) )

		! Compute after tax wealth according to threshold
		if (Y_a.le.Y_a_threshold) then 
			Y_a = Y_a * (1-tauW_bt)
		else
			Y_a = Y_a * (1-tauW_at)
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
		real(DP), intent(in) :: a_in, z_in
		real(DP)             :: MB_a, Y_a

		! Before tax wealth
		Y_a = ( a_in + ( rr * (z_in * a_in )**mu - DepRate*a_in ) *(1.0_DP-tauK) )

		! After tax marginal benefit of assets
		if (Y_a.le.Y_a_threshold) then 
			! Compute asset marginal benefit - tax free
			MB_a = ( 1.0_DP + ( rr *mu* (z_in**mu) * (a_in**(mu-1.0_DP)) - DepRate ) *(1.0_DP-tauK) )*(1.0_DP-tauW_bt)
		else
			! Compute asset marginal benefit - subject to taxes
			MB_a = ( 1.0_DP + ( rr *mu* (z_in**mu) * (a_in**(mu-1.0_DP)) - DepRate ) *(1.0_DP-tauK) )*(1.0_DP-tauW_at)
		end if
	END  FUNCTION MB_a

	FUNCTION MB_a_at(a_in,z_in)
		IMPLICIT NONE   
		real(DP), intent(in) :: a_in, z_in
		real(DP)             :: MB_a_at

		! Compute asset marginal benefit - subject to taxes
		MB_a_at = ( 1.0_DP + ( rr *mu* (z_in**mu) * (a_in**(mu-1.0_DP)) - DepRate ) *(1.0_DP-tauK) )*(1.0_DP-tauW_at)

	END  FUNCTION MB_a_at

	FUNCTION MB_a_bt(a_in,z_in)
		IMPLICIT NONE   
		real(DP), intent(in) :: a_in, z_in
		real(DP)             :: MB_a_bt

		! Compute asset marginal benefit - subject to taxes
		MB_a_bt = ( 1.0_DP + ( rr *mu* (z_in**mu) * (a_in**(mu-1.0_DP)) - DepRate ) *(1.0_DP-tauK) )*(1.0_DP-tauW_bt)

	END  FUNCTION MB_a_bt


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
	FUNCTION FOC_WH(aprimet)
		IMPLICIT NONE   
		real(DP), intent(in) :: aprimet
		real(DP)             :: ntemp, MBaprime, FOC_WH,  yprime, exp1overcprime
		REAL(DP)             :: brentvaluet
		integer              :: epindx
		real(DP), dimension(ne):: cprime
		

		! Set auxiliary variable for FOC_HA
		ain=aprimet
		! Solve for hours choice by solving the FOC for labor
		brentvaluet = brent(0.000001_DP, 0.4_DP, 0.99_DP, FOC_HA, brent_tol, ntemp)           
		
		! Compute marginal benefit of next period at a' (given zi)
		MBaprime = MB_a(aprimet,zgrid(zi))
		 
		! Compute asset income of next period at a' (given zi)
		yprime   = Y_a(aprimet,zgrid(zi))

		! I have to evaluate the FOC in expectation over eindx prime given eindx
		! Compute c' for each value of e'
		DO epindx=1,ne
		    cprime(epindx) = Linear_Int(Ygrid_t(:,zi),Cons_t(age+1,:,zi,lambdai,epindx), na_t, yprime  )
		ENDDO
		! Compute the expected value of 1/c' conditional on current ei
		exp1overcprime = SUM( pr_e(ei,:) / cprime )

		! Evaluate the squared residual of the Euler equation for working period
		FOC_WH   = (1.0_DP / (YGRID_t(ai,zi)  + psi* (yh(age, lambdai,ei) * ntemp)**(1.0_DP-tauPL ) - aprimet )  &
		           & - beta *  survP(age) *MBaprime* exp1overcprime) **2.0_DP 

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
	FUNCTION FOC_H(hoursin)
		IMPLICIT NONE   
		real(DP), intent(in) :: hoursin
		real(DP)             :: FOC_H 

		FOC_H  = ( (psi*(1.0_DP-tauPL)*yh(age, lambdai,ei)**(1.0_DP-tauPL) )* (1.0_DP-hoursin)*(hoursin**(-tauPL)) - phi*consin )**2.0_DP

	END  FUNCTION FOC_H

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
	FUNCTION FOC_HA(hoursin)
		IMPLICIT NONE   
		real(DP), intent(in) :: hoursin
		real(DP)   			 :: FOC_HA 

		FOC_HA  = ( (psi*(1.0_DP-tauPL)*yh(age, lambdai,ei)**(1.0_DP-tauPL) )* (1.0_DP-hoursin)*(hoursin**(-tauPL)) &
		        &  - phi*  ( YGRID_t(ai,zi) + psi*(yh(age, lambdai,ei)*hoursin)**(1.0_DP-tauPL) - ain  )   )**2.0_DP

	END  FUNCTION FOC_HA

!========================================================================================
!========================================================================================
! Given arrays xa(1:n) and ya(1:n) of length n; this subroutine returns  linear interpolated value "y" at point "x"
	FUNCTION Linear_Int(xa,ya,n,x)  
		USE parameters
		IMPLICIT NONE      
		INTEGER    :: n  
		REAL(DP)   :: x, xa(n),ya(n)  
		INTEGER    :: k,khi,klo  
		REAL(DP)   :: a,b,h, Linear_Int
		      
		   
		klo=1  
		khi=n  

		1     if (khi-klo.gt.1) then  
		        k=(khi+klo)/2  
		        
		        
		        if(xa(k).gt.x)then  
		          khi=k  
		        else  
		          klo=k  
		        endif  
		      goto 1  
		      endif
		 
		      h=xa(khi)-xa(klo)  
		      if (h.eq.0.) then 
		      	print*,'bad xa input in linear int'  
		      end if 
		      a=(xa(khi)-x)/h  
		      b=(x-xa(klo))/h  
		      Linear_Int=a*ya(klo)+b*ya(khi)     
		     return  

	END  Function Linear_Int

!========================================================================================
!========================================================================================

	FUNCTION brent(ax,bx,cx,func,tol,xmin)
		USE nrtype; USE nrutil, ONLY : nrerror
		IMPLICIT NONE
		
		REAL(DP), INTENT(IN) :: ax,bx,cx,tol
		REAL(DP), INTENT(OUT) :: xmin
		REAL(DP) :: brent
		INTERFACE
			FUNCTION func(x)
				USE nrtype
				IMPLICIT NONE
				REAL(DP), INTENT(IN) :: x
				REAL(DP) :: func
			END FUNCTION func
		END INTERFACE
		INTEGER(I4B), PARAMETER :: ITMAX=100
		REAL(DP), PARAMETER :: CGOLD=0.3819660_DP,ZEPS=1.0e-3_DP*epsilon(ax)
		INTEGER(I4B) :: iter
		REAL(DP) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
		
		a=min(ax,cx)
		b=max(ax,cx)
		v=bx
		w=v
		x=v
		e=0.0
		fx=func(x)
		fv=fx
		fw=fx
		do iter=1,ITMAX

			xm=0.5_DP*(a+b)
			tol1=tol*abs(x)+ZEPS
			tol2=2.0_DP*tol1
			if (abs(x-xm) <= (tol2-0.5_DP*(b-a))) then
				xmin=x
				brent=fx
				RETURN
			end if
			if (abs(e) > tol1) then
				r=(x-w)*(fx-fv)
				q=(x-v)*(fx-fw)
				p=(x-v)*q-(x-w)*r
				q=2.0_DP*(q-r)
				if (q > 0.0) p=-p
				q=abs(q)
				etemp=e
				e=d
				if (abs(p) >= abs(0.5_DP*q*etemp) .or. &
					p <= q*(a-x) .or. p >= q*(b-x)) then
					e=merge(a-x,b-x, x >= xm )
					d=CGOLD*e
				else
					d=p/q
					u=x+d
					if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
				end if
			else
				e=merge(a-x,b-x, x >= xm )
				d=CGOLD*e
			end if
			u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
			fu=func(u)
			if (fu <= fx) then
				if (u >= x) then
					a=x
				else
					b=x
				end if
				call shft(v,w,x,u)
				call shft(fv,fw,fx,fu)
			else
				if (u < x) then
					a=u
				else
					b=u
				end if
				if (fu <= fw .or. w == x) then
					v=w
					fv=fw
					w=u
					fw=fu
				else if (fu <= fv .or. v == x .or. v == w) then
					v=u
					fv=fu
				end if
			end if
		end do
		call nrerror('brent: exceed maximum iterations')
		
		CONTAINS
		SUBROUTINE shft(a,b,c,d)
			REAL(DP), INTENT(OUT) :: a
			REAL(DP), INTENT(INOUT) :: b,c
			REAL(DP), INTENT(IN) :: d
			a=b
			b=c
			c=d
		END SUBROUTINE shft
	END FUNCTION brent

!========================================================================================
!========================================================================================

	FUNCTION ran1(idum)  
	      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV  
	      REAL(DP):: ran1,AM,EPS,RNMX  
	      PARAMETER (IA=16807.0_DP,IM=2147483647.0_DP,AM=1.0_DP/IM,IQ=127773,IR=2836,  &
	     & NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.0_DP-EPS)  
	      INTEGER j,k,iv(NTAB),iy  
	      SAVE iv,iy  
	      DATA iv /NTAB*0/, iy /0/  
	      
	      if (idum.le.0.or.iy.eq.0) then  
	        idum=max(-idum,1)  
	        do 11 j=NTAB+8,1,-1  
	          k=idum/IQ  
	          idum=IA*(idum-k*IQ)-IR*k  
	          if (idum.lt.0) idum=idum+IM  
	          if (j.le.NTAB) iv(j)=idum  
	11      continue  
	        iy=iv(1)  
	      endif  
	      k=idum/IQ  
	      idum=IA*(idum-k*IQ)-IR*k  
	      if (idum.lt.0) idum=idum+IM  
	      j=1+iy/NDIV  
	      iy=iv(j)  
	      iv(j)=idum        
	      ran1=min(AM*iy,RNMX)  
	      
	      return  
	END  FUNCTION ran1

!========================================================================================
!========================================================================================
! Sort: Sorts in ascending the elements of a one dimensional array of type real(8) 
!       It also gives the original indeces of the sorted elements
!
! Usage: Call Sort(n,A,A_sort,Sort_ind)
!
! Input: n   , integer(4), the number of elements in A
!        A   , real(8), the array whose elements are to be sorted
!
! Output: A_sort, real(8), array with sorted elements (in ascending order)
!		  Sort_ind, integer(4), array with original indeces of the elements in A_sort
!

	Subroutine Sort(n,A,A_sort,Sort_ind)
		integer , intent(in) :: n    !Number of elements in A
		real(dp), intent(in) , dimension(n) :: A
		real(dp), intent(out), dimension(n) :: A_sort
		integer , intent(out), dimension(n) :: Sort_ind
		integer :: i,j

		A_sort = A
		do i=1,n
			Sort_ind(i)=i
		end do

		do i=1,(n-1)
			do j=i+1,n
				if (A_sort(i) .ge. A_sort(j)) then
					A_sort((/i,j/))   = A_sort((/j,i/))
					Sort_ind((/i,j/)) = Sort_ind((/j,i/))
				end if
			end do
		end do

		return
	End Subroutine Sort	

!========================================================================================
!========================================================================================

	FUNCTION zbrent_p(func,x1,x2,tol,par)
		USE nrtype; USE nrutil, ONLY : nrerror
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x1,x2,tol
		REAL(DP), dimension(:), allocatable, INTENT(IN) :: par
		REAL(DP) :: zbrent_p
		INTERFACE
			FUNCTION func(x,par)
			USE nrtype
			IMPLICIT NONE
			REAL(DP), INTENT(IN) :: x
			REAL(DP), dimension(:), allocatable, INTENT(IN) :: par
			REAL(DP) :: func
			END FUNCTION func
		END INTERFACE
		INTEGER(I4B), PARAMETER :: ITMAX=100
		REAL(DP), PARAMETER :: EPS=epsilon(x1)
		INTEGER(I4B) :: iter
		REAL(DP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
		a=x1
		b=x2
		fa=func(a,par)
		fb=func(b,par)
		if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) &
			call nrerror('root must be bracketed for zbrent_p')
		c=b
		fc=fb
		do iter=1,ITMAX
			if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
				c=a
				fc=fa
				d=b-a
				e=d
			end if
			if (abs(fc) < abs(fb)) then
				a=b
				b=c
				c=a
				fa=fb
				fb=fc
				fc=fa
			end if
			tol1=2.0_DP*EPS*abs(b)+0.5_DP*tol
			xm=0.5_DP*(c-b)
			if (abs(xm) <= tol1 .or. fb == 0.0) then
				zbrent_p=b
				RETURN
			end if
			if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
				s=fb/fa
				if (a == c) then
					p=2.0_DP*xm*s
					q=1.0_DP-s
				else
					q=fa/fc
					r=fb/fc
					p=s*(2.0_DP*xm*q*(q-r)-(b-a)*(r-1.0_DP))
					q=(q-1.0_DP)*(r-1.0_DP)*(s-1.0_DP)
				end if
				if (p > 0.0) q=-q
				p=abs(p)
				if (2.0_DP*p  <  min(3.0_DP*xm*q-abs(tol1*q),abs(e*q))) then
					e=d
					d=p/q
				else
					d=xm
					e=d
				end if
			else
				d=xm
				e=d
			end if
			a=b
			fa=fb
			b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
			fb=func(b,par)
		end do
		call nrerror('zbrent_p: exceeded maximum iterations')
		zbrent_p=b
		
	END FUNCTION zbrent_p

end module programfunctions


!========================================================================================
!========================================================================================
!========================================================================================


PROGRAM main
	USE parameters
	USE GLOBAL
	use  programfunctions

	IMPLICIT NONE
	! Variables to measure running time
		REAL(DP) :: start_time, finish_time
	! Values of parameters for different parametrizations (Low and High)
		REAL(DP) :: betaL, betaH, phiL, phiH, sigmalambL,sigmalambH, sigmazL, sigmazH, rhozL, rhozH
	! Number of different values to be considered for each parameter 
		INTEGER  ::nbeta, nphi, nsigmalambda, nrhoz, nsigmaz
	! Counters for runs of the model
		INTEGER  :: parindx1,  parindx2, parindx3, parindx4, parindx5
	! Compute benchmark or load results
		INTEGER  :: read_write_bench
	
	! Unused values of parameters
		! the following solves equilibrium of capital tax economy
		params=[  0.9500,    0.7000,    0.2996,    0.5673,    1.2280]
		params=[  0.9510,    0.5000,    0.4060,    0.5680,    1.2250]
		params=[  0.9503,    0.6500,    0.3288,    0.5667,    1.2280]
		params=[  0.9510,    0.5250,    0.3942,    0.5680,    1.2247]
		params=[  0.9511,    0.4200,    0.4400,    0.5680,    1.2250]
		params=[  0.9522,    0.3400,    0.4740,    0.5690,    1.2240]
		params=[  0.9506,    0.6000,    0.3564,    0.5667,    1.2280]

		! NEW PARAMETERS

		params=[ 0.947  ,  0.4 , 0.490 , 0.340 , 1.01 ]
		params=[ 0.9455 ,  0.6 , 0.381 , 0.335 , 1.00 ]
		params=[ 0.9455 ,  0.8 , 0.255 , 0.34  , 1.00 ]
		params=[ 0.948  ,  0.2 , 0.56  , 0.34  , 1.02 ]


	!print*,'---------------------------       PSI    NOT  ADJUSTED   ---------------------------'
	!print*,'------------------------- RETIREMENT BENEFITS ADJUSTED - DBN ADJUSTED ------------------------'
	!print*,'------------------------- CONS TAX SIMPLE ------------------------'
	print*,'na=',na,'update_period=',update_period

	! Set parameters to be used in all simulations economy
		OPEN   (UNIT=3, FILE='params', STATUS='replace')
		WRITE(unit=3, FMT=*) params
		CLOSE (unit=3)
		
		beta             = params(1)
		rho_z            = params(2)
		sigma_z_eps      = params(3)
		sigma_lambda_eps = params(4)
		phi              = params(5)
		!print*,  beta, rho_z, sigma_z_eps, sigma_lambda_eps,  phi

	! Start timining of  the process
		call cpu_time(start_time) 

		
	! Set initia lvalues of R, Wage, Ebar to find equilibrium
		! ------- DO NOT REMOVE THE LINES BELOW

		rr    =  4.906133597851297E-002 
		wage  =  1.97429920063330 
		Ebar  =  1.82928004963637  
		Ebar_bench = Ebar
		 

		! ------- DO NOT REMOVE THE LINES ABOVE

	!--------------------------------------------------------------------------------------------------------------------------

		!call cpu_time(start_time) 
		! 
		!
		!
		!betaL =0.947
		!betaH=0.948
		!rhozL = 0.4
		!rhozH = 0.4
		!sigmazL = 0.48
		!sigmazH = 0.49
		!sigmalambL = 0.34
		!sigmalambH =0.35
		!phiL    = 0.99
		!phiH   = 1.01
		!
		!
		!!betaL =0.9455
		!!betaH=0.946
		!!rhozL = 0.6
		!!rhozH = 0.6
		!!sigmazL = 0.381
		!!sigmazH = 0.385
		!!sigmalambL = 0.335
		!!sigmalambH =0.34
		!!phiL    = 0.99
		!!phiH   = 1.00
		!
		!
		!betaL =0.9455
		!betaH=0.946
		!rhozL = 0.8
		!rhozH = 0.8
		!sigmazL = 0.253
		!sigmazH = 0.255
		!sigmalambL = 0.335
		!sigmalambH =0.34
		!phiL    = 1.00
		!phiH   = 1.02
		!
		!betaL =0.948
		!betaH=0.949
		!rhozL = 0.2
		!rhozH = 0.2
		!sigmazL = 0.56
		!sigmazH = 0.57
		!sigmalambL = 0.34
		!sigmalambH =0.35
		!phiL    = 1.02
		!phiH   = 1.04
		!
		!
		!nbeta =2
		!nrhoz=1
		!nsigmaz=2
		!nsigmalambda=2
		!nphi=2
		!
		!Min_SSE_Moments=1000.0_DP
		!
		!DO parindx3=1,nsigmalambda
		!DO parindx2=1,nphi
		!DO parindx1=1,nbeta
		!DO parindx4=1,nrhoz
		!DO parindx5=1,nsigmaz
		!
		!    beta = betaL + real(parindx1-1,8) *(betaH-betaL)/max(real(nbeta-1,8),1.0_DP)
		!    phi   = phiL   + real(parindx2-1,8) *(phiH-phiL)/max(real(nphi-1,8),1.0_DP)
		!    sigma_lambda_eps = sigmalambL + real(parindx3-1,8)*(sigmalambH -sigmalambL) / max(real(nsigmalambda-1,8),1.0_DP)
		!    rho_z= rhozL   +  real(parindx4-1,8)*(rhozH-rhozL) / max(real(nrhoz-1,8),1.0_DP)
		!    sigma_z_eps = sigmazL +  real(parindx5-1,8)*(sigmazH-sigmazL) / max(real(nsigmaz-1,8),1.0_DP)
		!
		!    CALL  INITIALIZE
		!    CALL FIND_DBN_EQ
		!    CALL COMPUTE_STATS
		!    IF (SSE_Moments .lt. Min_SSE_Moments ) THEN
		!        Min_SSE_Moments =SSE_Moments
		!        params= [ beta, rho_z, sigma_z_eps, sigma_lambda_eps, phi ]
		!        Min_Moments = [  Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60  ]
		!    ENDIF
		!    !CALL WRITE_TO_FILE
		!
		!ENDDO
		!ENDDO
		!ENDDO
		!ENDDO
		!ENDDO
		!print*, params
		!print*, Min_Moments
	!

	PRINT*,''
	Print*,'--------------- SOLVING BENCHMARK WITH BEST PARAMETERS -----------------'
	PRINT*,''
	print*,'CAPITAL TAX ECONOMY'

	! Benchmark economy
		solving_bench=1
	! Set taxes for benchmark economy
		tauK = 0.25_DP
		tauL = 0.30_DP
		tauW_bt = 0.00_DP
		tauW_at = 0.00_DP
		Y_a_threshold = 0.00_DP 

	! Solve for the model and compute stats
	read_write_bench = 0
	print*,"	Initializing program"
		CALL INITIALIZE
	if (read_write_bench.eq.0) then
		print*,"	Computing equilibrium distribution"
		CALL FIND_DBN_EQ
		print*,"	Computing satitics"
		CALL COMPUTE_STATS
		print*,"	Computing government spending"
		CALL GOVNT_BUDGET
		print*,"	Writing variables"
		CALL WRITE_VARIABLES(1)
		print*,"	Saving results in text files to be read later"
		CALL Write_Benchmark_Results(read_write_bench)
	else
		print*,"	Reading benchmark results from files"
		CALL Write_Benchmark_Results(read_write_bench)
	end if 

	! Aggregate variables in benchmark economy
		GBAR_bench  = GBAR
		QBAR_bench  = QBAR 
		NBAR_bench  = NBAR 
		Ebar_bench  = EBAR
		rr_bench    = rr
		wage_bench  = wage
		Y_bench     = YBAR
		tauK_bench  = tauK
		tauL_bench  = tauL
		DBN_bench   = DBN1
		tauw_bt_bench = tauW_bt
		tauw_at_bench = tauW_at
		Y_a_threshold_bench = Y_a_threshold

		write(*,*) "Benchmark variables"
		write(*,*) "GBAR=",GBAR,"EBAR=",EBAR,"NBAR=",NBAR,"QBAR=",QBAR,"rr=",rr,"wage=",wage

	PRINT*,''
	Print*,'--------------- SOLVING EXPERIMENT WITH BEST PARAMETERS -----------------'
	PRINT*,''
	print*,'Wealth Tax Economy'
	
	! Experiment economy
		solving_bench=0
	! Set capital taxes to zero
		tauK = 0.0_DP
	! Set Y_a_threshold
		write(*,*) "Y_a threshold is set to a proportion of the mean wealth under current distribution"
		Y_a_threshold = 0.0_dp ! big_p   !8.1812138704441200  
		!call Find_TauW_Threshold(DBN1,Y_a_threshold)
		!Y_a_threshold = 10.0_dp

	! Find wealth taxes that balances budget
	print*, "	Computing Wealth Tax to balance the budget"
		! Set initial value for G in experimental economy and for wealth taxes
		GBAR_exp = 0.0_DP
		tauW_bt  = tauWmin_bt
		tauW_at  = tauWmin_at
		tauWindx = 0.0_DP
		! Solve for the model increasing wealth taxes until revenue is enough to finance G_benchamark
		DO WHILE (GBAR_exp .lt. GBAR_bench)
			! Set old G and new value of tauW
			GBAR_exp_old = GBAR_exp
			tauW_bt = tauWmin_bt + tauWindx * tauWinc_bt
			tauW_at = tauWmin_at + tauWindx * tauWinc_at
			! Solve the model
			CALL FIND_DBN_EQ
			CALL COMPUTE_STATS
			CALL GOVNT_BUDGET
			! Get new G
			GBAR_exp = GBAR 
			! Iteratioins  
			tauWindx = tauWindx + 1.0_DP   
			write(*,*) "Bracketing GBAR: tauW_bt=", tauW_bt*100, "And tauW_at=", tauW_at*100
			print*, "Current Threshold for wealth taxes", Y_a_threshold
			print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
		ENDDO

		! Set tauW as weighted average of point in  the grid to balance budget more precisely
			tauW_up_bt  = tauW_bt
			tauW_low_bt = tauW_bt  -  tauWinc_bt
			tauW_bt     = tauW_low_bt + tauWinc_bt * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
			tauW_up_at  = tauW_at
			tauW_low_at = tauW_at  -  tauWinc_at  
			tauW_at     = tauW_low_at + tauWinc_at * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
			print*,''
			print*,'GBAR bracketed by taxes:'
			print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
			print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
			print*,''

		! Solve (again) experimental economy
			CALL FIND_DBN_EQ
			CALL COMPUTE_STATS
			CALL GOVNT_BUDGET

		! Find tauW that exactly balances the budget (up to precisioin 0.1) using bisection
			GBAR_exp = GBAR
			print*,"Gbar at midpoint of bracket and GBAR at benchmark"
			print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
			print*,''
			print*,'Bisection for TauW:'
			DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.1 ) ! as long as the difference is greater than 0.1% continue
			    if (GBAR_exp .gt. GBAR_bench ) then
			        tauW_up_bt  = tauW_bt 
			        tauW_up_at  = tauW_at 
			    else
			        tauW_low_bt = tauW_bt
			        tauW_low_at = tauW_at
			    endif
			    tauW_bt = (tauW_low_bt + tauW_up_bt)/2.0_DP
			    tauW_at = (tauW_low_at + tauW_up_at)/2.0_DP
			    CALL FIND_DBN_EQ
			    CALL COMPUTE_STATS
			    CALL GOVNT_BUDGET
			    GBAR_exp = GBAR
			    print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
				print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
				print*, "Current Threshold for wealth taxes", Y_a_threshold
				print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
			ENDDO

	! AGgregate variable in experimental economy
		GBAR_exp  = GBAR
		QBAR_exp  = QBAR 
		NBAR_exp  = NBAR  
		Y_exp 	  = YBAR
		Ebar_exp  = EBAR
		rr_exp    = rr
		wage_exp  = wage
		tauK_exp  = tauK
		tauL_exp  = tauL
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold

	CALL WRITE_VARIABLES(0)

	! Compute welfare gain between economies
	CALL COMPUTE_WELFARE_GAIN
	print*,'---------------------------'
	print*,''
	print*,'Output Gain Prct=', 100.0_DP*(Y_exp/Y_bench-1.0) 
	print*,''
	print*,'---------------------------'

	call cpu_time(finish_time)
	print*,'Total time =',finish_time-start_time

	print*, " "
	print*, "Asset Function"
	do ai=1,na 
		print*,agrid(ai), Aprime(10,ai,:,1,3)
	end do 

END PROGRAM main


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE COMPUTE_WELFARE_GAIN
	use GLOBAL 
	use programfunctions
	IMPLICIT NONE
	real(DP), DIMENSION(MaxAge):: CumDiscountF
	REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) ::  ValueFunction_Bench, ValueFunction_Exp, Cons_Eq_Welfare
	REAL(DP), dimension(nz) ::  temp_ce_by_z

	! Discount factor
		CumDiscountF(MaxAge)=1.0_DP
		DO age=MaxAge-1,1,-1
		    CumDiscountF(age)   = 1.0_DP + beta * survP(age) *CumDiscountF(age+1) 
		ENDDO
		!print*,CumDiscountF
		!PAUSE

	! Solve for the benchmark economy 
		solving_bench = 1
		tauK  = tauK_bench
		tauL  = tauL_bench
		rr    = rr_bench
		wage  = wage_bench
		Ebar  = Ebar_bench
		tauW_bt  = tauW_bt_bench
		tauW_at  = tauW_at_bench
		Y_a_threshold = Y_a_threshold_bench
		print*,'BENCH: rr=',rr,'wage=',wage,'Ebar=',Ebar
		CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
		CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
		CALL ComputeLaborUnits(Ebar, wage) 
		CALL EGM_RETIREMENT_WORKING_PERIOD 

		! Compute the value function using interpolation and save it
			!CALL COMPUTE_VALUE_FUNCTION_LINEAR
			CALL COMPUTE_VALUE_FUNCTION_SPLINE  
			ValueFunction_Bench = ValueFunction
		
		!CALL COMPUTE_STATS
		
		! Print policy and value function benchmark economy
			!OPEN (UNIT=7, FILE='c_bench', STATUS='replace')    
			!OPEN (UNIT=8, FILE='n_bench', STATUS='replace')    
			!OPEN (UNIT=9, FILE='ap_bench', STATUS='replace')    
			!OPEN (UNIT=10, FILE='v_bench', STATUS='replace')    
			!DO age=1,MaxAge 
			!DO ai=1,na    
			!    DO zi=1,nz
			!        DO lambdai=1,nlambda          
			!             DO ei=1,ne
			!                  WRITE  (UNIT=7, FMT=*) cons(age, ai, zi, lambdai, ei)
			!                  WRITE  (UNIT=8, FMT=*) HOURS(age, ai, zi, lambdai, ei)
			!                  WRITE  (UNIT=9, FMT=*) Aprime(age, ai, zi, lambdai, ei)
			!                  WRITE  (UNIT=10, FMT=*) ValueFunction_Bench(age, ai, zi, lambdai, ei)
			!               ENDDO ! ei          
			!        ENDDO ! lambdai
			!    ENDDO ! zi
			!ENDDO ! ai
			!ENDDO
			!close (unit=7)
			!close (unit=8)
			!close (unit=9)
			!close (unit=10)

	! Solve the experimental economy  
		solving_bench = 0  
		tauK  = tauK_exp
		tauL  = tauL_exp
		rr    = rr_exp
		wage  = wage_exp
		Ebar  = Ebar_exp
		tauW_bt  = tauW_bt_exp
		tauW_at  = tauW_at_exp
		Y_a_threshold = Y_a_threshold_exp
		print*,' EXP: rr=',rr,'wage=',wage,'Ebar=',Ebar
		CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
		CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
		CALL ComputeLaborUnits(Ebar, wage) 
		CALL EGM_RETIREMENT_WORKING_PERIOD 

		! Compute the value function using interpolation and save it
			!CALL COMPUTE_VALUE_FUNCTION_LINEAR
			CALL COMPUTE_VALUE_FUNCTION_SPLINE 
			ValueFunction_Exp = ValueFunction

		!CALL COMPUTE_STATS

		! Print policy and value function benchmark economy
			!OPEN (UNIT=7, FILE='c_exp', STATUS='replace')    
			!OPEN (UNIT=8, FILE='n_exp', STATUS='replace')    
			!OPEN (UNIT=9, FILE='ap_exp', STATUS='replace')    
			!OPEN (UNIT=10, FILE='v_exp', STATUS='replace')    
			!DO age=1,MaxAge 
			!DO ai=1,na    
			!    DO zi=1,nz
			!        DO lambdai=1,nlambda          
			!             DO ei=1,ne
			!                  WRITE  (UNIT=7, FMT=*) cons(age, ai, zi, lambdai, ei)
			!                  WRITE  (UNIT=8, FMT=*) HOURS(age, ai, zi, lambdai, ei)
			!                  WRITE  (UNIT=9, FMT=*) Aprime(age, ai, zi, lambdai, ei)
			!                  WRITE  (UNIT=10, FMT=*) ValueFunction_exp(age, ai, zi, lambdai, ei)
			!               ENDDO ! ei          
			!        ENDDO ! lambdai
			!    ENDDO ! zi
			!ENDDO ! ai
			!ENDDO
			!close (unit=7)
			!close (unit=8)
			!close (unit=9)
			!close (unit=10)

	! Compute consumption equivalent - Average consumption equivalents are stored
	OPEN (UNIT=10, FILE='CE', STATUS='replace')  
	OPEN (UNIT=11, FILE='CE_by_age', STATUS='replace')  
	OPEN (UNIT=12, FILE='CE_by_age_z', STATUS='replace')  

	WRITE  (UNIT=10, FMT=*) 100.0_DP*sum(Cons_Eq_Welfare*DBN_bench)
	DO age=MaxAge,1,-1
	    Cons_Eq_Welfare(age,:,:,:,:)=exp((ValueFunction_exp(age,:,:,:,:)-ValueFunction_Bench(age,:,:,:,:))/CumDiscountF(age))-1.0_DP
	    WRITE  (UNIT=11, FMT=*) 100*sum(Cons_Eq_Welfare(age,:,:,:,:)*DBN_bench(age,:,:,:,:))/sum(DBN_bench(age,:,:,:,:))
	    DO zi=1,nz
	         temp_ce_by_z(zi) = 100*sum(Cons_Eq_Welfare(age,:,zi,:,:)*DBN_bench(age,:,zi,:,:))/sum(DBN_bench(age,:,zi,:,:))
	    ENDDO
	    WRITE  (UNIT=12, FMT=*) temp_ce_by_z
	ENDDO

	close (unit=10)
	close (unit=11)
	close (unit=12)

	print*,'---------------------------'
	print*,''
	print*,'Average Welfare Gain Whole Population (prct)=',100.0_DP*sum(Cons_Eq_Welfare*DBN_bench)
	print*,''
	print*,'---------------------------'

	! 
		!OPEN (UNIT=7, FILE='value_bench', STATUS='replace')    
		!WRITE  (UNIT=7, FMT=*) ValueFunction_Bench
		!close (unit=7)
		!
		!OPEN (UNIT=7, FILE='value_exp', STATUS='replace')    
		!WRITE  (UNIT=7, FMT=*) ValueFunction_Exp
		!close (unit=7)
		!
		!OPEN (UNIT=7, FILE=' cons_eq_welfare', STATUS='replace')    
		!WRITE  (UNIT=7, FMT=*)  Cons_Eq_Welfare
		!close (unit=7)
		!
		!OPEN (UNIT=1, FILE=' vbench', STATUS='replace')    
		!OPEN (UNIT=2, FILE=' vexp', STATUS='replace')    
		!OPEN (UNIT=3, FILE=' vdbn', STATUS='replace')    
		!OPEN (UNIT=4, FILE=' v_CE', STATUS='replace')    
		!DO age=1,MaxAge 
		!DO ai=1,na    
		!    DO zi=1,nz
		!        DO lambdai=1,nlambda          
		!             DO ei=1,ne
		!                    WRITE  (UNIT=1, FMT=*) ValueFunction_Bench(age, ai, zi, lambdai, ei)
		!                    WRITE  (UNIT=2, FMT=*) ValueFunction_exp(age, ai, zi, lambdai, ei)
		!                    WRITE  (UNIT=3, FMT=*) DBN1(age, ai, zi, lambdai, ei)                  
		!                    WRITE  (UNIT=4, FMT=*) Cons_Eq_Welfare(age, ai, zi, lambdai, ei)                 
		!             ENDDO ! ei          
		!        ENDDO ! lambdai
		!    ENDDO ! zi
		!ENDDO ! ai
		!ENDDO
		!close (unit=1)
		!close (unit=2)
		!close (unit=3)
		!close (unit=4)
		!
		!

END SUBROUTINE  COMPUTE_WELFARE_GAIN


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE COMPUTE_VALUE_FUNCTION_SPLINE
	use global 
	use  programfunctions
	IMPLICIT NONE
	INTEGER :: tklo, tkhi
	REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) :: PrAprimelo, PrAprimehi
	REAL(DP), DIMENSION(na) :: ValueP1, ValueP2, ValueP, ExpValueP

	! Announce method of interpolation
	! Interpolation is used to get value function at optimal a' (when computing expectations)
	print*,'VALUE FUNCTION SPLINE'

	! Compute value function at grid nodes
	! Final age
	age=MaxAge
	DO ai=1,na    
	    DO zi=1,nz
	        DO lambdai=1,nlambda          
	              DO ei=1,ne
	                  ValueFunction(age, ai, zi, lambdai, ei) = log(Cons(age, ai, zi, lambdai, ei)) &
	                       & + phi * log(1.0_DP-Hours(age, ai, zi, lambdai, ei))
	                  ! print*,Cons(age, ai, zi, lambdai, ei),  ValueFunction(age, ai, zi, lambdai, ei) 
	                  ! pause
	              ENDDO ! ei          
	        ENDDO ! lambdai
	    ENDDO ! zi
	ENDDO ! ai

	! Working Period
	DO age=MaxAge-1,1,-1
	    DO zi=1,nz
	        DO lambdai=1,nlambda          
	              DO ei=1,ne
	                    DO ai=1,na    
	                          ExpValueP(ai) = sum(ValueFunction(age+1, ai, zi, lambdai, :) * pr_e(ei,:))
	                    ENDDO

	                    CALL spline( agrid, ExpValueP , na , &
	                    & sum(pr_e(ei,:)/Cons(age+1, 1, zi, lambdai,:)) , sum(pr_e(ei,:)/Cons(age+1, na, zi, lambdai,:)) , ValueP2)  
	                    
	                    DO ai=1,na                        
	                         call splint( agrid, ExpValueP, ValueP2, na, Aprime(age,ai,zi,lambdai, ei), ValueP(ai))   
	                         ValueFunction(age, ai, zi, lambdai, ei) = log(Cons(age, ai, zi, lambdai, ei)) &
	                               & + phi * log(1.0_DP-Hours(age, ai, zi, lambdai, ei)) + beta*survP(age)* ValueP(ai)                        
	                    ENDDO ! ai
	               ENDDO ! ei          
	        ENDDO ! lambdai
	    ENDDO ! zi
	ENDDO ! age

END SUBROUTINE COMPUTE_VALUE_FUNCTION_SPLINE 


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE COMPUTE_VALUE_FUNCTION_LINEAR
	use global
	use  programfunctions
	IMPLICIT NONE
	INTEGER :: tklo, tkhi
	REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) :: PrAprimelo, PrAprimehi

	print*,'VALUE FUNCTION LINEAR'

	age=MaxAge
	DO ai=1,na    
	    DO zi=1,nz
	        DO lambdai=1,nlambda          
	              DO ei=1,ne
	                  ValueFunction(age, ai, zi, lambdai, ei) = log(Cons(age, ai, zi, lambdai, ei)) &
	                  		& + phi * log(1.0_DP-Hours(age, ai, zi, lambdai, ei))
	!                  print*,Cons(age, ai, zi, lambdai, ei),  ValueFunction(age, ai, zi, lambdai, ei) 
	!                  pause
	              ENDDO ! ei          
	        ENDDO ! lambdai
	    ENDDO ! zi
	ENDDO ! ai

	! Working Period
	DO age=MaxAge-1,1,-1
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

		ValueFunction(age, ai, zi, lambdai, ei) = log(Cons(age, ai, zi, lambdai, ei)) &
		   & + phi * log(1.0_DP-Hours(age, ai, zi, lambdai, ei))  &
		   & + beta*survP(age)* sum( ( PrAprimelo(age,ai,zi,lambdai, ei) * ValueFunction(age+1, tklo, zi, lambdai,:)  &
		   & + PrAprimehi(age,ai,zi,lambdai, ei) * ValueFunction(age+1, tkhi, zi, lambdai,:)) * pr_e(ei,:) )
		if ( ValueFunction(age, ai, zi, lambdai, ei) .lt. (-100.0_DP) ) then
		   print*,'ValueFunction(age, ai, zi, lambdai, ei)=',ValueFunction(age, ai, zi, lambdai, ei)
		endif
	ENDDO ! ei          
    ENDDO ! lambdai
    ENDDO ! zi
	ENDDO ! ai
	ENDDO ! age

END SUBROUTINE COMPUTE_VALUE_FUNCTION_LINEAR 


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE GOVNT_BUDGET
	USE PARAMETERS
	USE GLOBAL
	use programfunctions
	IMPLICIT NONE
	real(DP) ::  GBAR_K,  GBAR_W,  GBAR_L, GBAR_C, SSC_Payments

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
	    GBAR = GBAR + DBN1(age,ai,zi,lambdai,ei) * ( tauK*( rr*(agrid(ai)*zgrid(zi))**mu-DepRate*agrid(ai) )  	&
	          & + (agrid(ai) + ( rr * (zgrid(zi) * agrid(ai) )**mu - DepRate*agrid(ai) ) *(1.0_DP-tauK)  )		&
	          & - Y_a(agrid(ai),zgrid(zi)) 																		&	
	          & + yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei)  												&
	          & - psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei))**(1.0_DP-tauPL)  					&
	          & + tauC * cons(age, ai, zi, lambdai,ei)  )   

	    GBAR_L = GBAR_L  + DBN1(age,ai,zi,lambdai,ei) * (  yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei) &
	               &- psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei))**(1.0_DP-tauPL) )

	    GBAR_K = GBAR_K + DBN1(age,ai,zi,lambdai,ei) * ( tauK*( rr*(agrid(ai)*zgrid(zi))**mu-DepRate*agrid(ai) ) )

	    GBAR_W = GBAR_W + DBN1(age,ai,zi,lambdai,ei) * &
	            & ((agrid(ai) + ( rr * (zgrid(zi) * agrid(ai) )**mu - DepRate*agrid(ai) ) *(1.0_DP-tauK)  )		&
	          	& - Y_a(agrid(ai),zgrid(zi)) )

      	GBAR_C = GBAR_C +  DBN1(age,ai,zi,lambdai,ei) * tauC * cons(age, ai, zi, lambdai,ei)
	    
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO


	! Print Results
	print*, ' '
	print*, "Government Budget - Revenues and taxes"
	print*,'GBAR=',GBAR,'GBAR_L=',GBAR_L,'GBAR_K=', GBAR_K, "GBAR_W=", GBAR_W, 'GBAR_C=', GBAR_C
	print*, 'Tau_K=', tauK, 'Tau_W=', tauW_at, 'Tau_C=', tauC, "Threshold", Y_a_threshold
	print*, ' '
END  SUBROUTINE GOVNT_BUDGET

!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE FIND_DBN_EQ
	USE PARAMETERS
	USE GLOBAL
	use  programfunctions
	IMPLICIT NONE
	INTEGER:: tklo, tkhi, age1, age2, z1, z2, a1, a2, lambda1, lambda2, e1, e2, DBN_iter, simutime, iter_indx
	REAL   :: DBN_dist, DBN_criteria
	REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) :: PrAprimelo, PrAprimehi, DBN2
	INTEGER,  DIMENSION(MaxAge, na, nz, nlambda, ne) :: Aplo, Aphi

	DBN_criteria = 1.0E-07_DP

	! Solve the model at current aggregate values
		! Find the threshold for wealth taxes (a_bar)
			!call Find_TauW_Threshold(DBN1,Y_a_threshold)
		! Adjust grid to include breaking points
			CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
		! Compute labor units 
			CALL ComputeLaborUnits(Ebar, wage) 
		! Form YGRID for the capital income economy given interest rate "rr"
			CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
		! Solve for policy and value functions 
			CALL EGM_RETIREMENT_WORKING_PERIOD 

	! Discretize policy function for assets (a')
		! For each age and state vector bracket optimal a' between two grid points
		! When at that age and state the optimal decision is approximated by selecting one the grid points
		! The grid points are selected with probability proportional to their distance to the optimal a'
	DO age=1,MaxAge
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
	DO WHILE ( ( DBN_dist .ge. DBN_criteria ) .and. ( simutime .le. MaxSimuTime ) )
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

		
	    ! Working age agents
	    DO age1=1,MaxAge-1
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
	             QBAR= QBAR+ DBN1(age1, a1, z1, lambda1, e1) * ( zgrid(z1) *agrid(a1) )**mu
	             NBAR= NBAR+ DBN1(age1, a1, z1, lambda1, e1) * eff_un(age1, lambda1, e1) * Hours(age1, a1, z1, lambda1,e1)
	        ENDDO
	        ENDDO
	        ENDDO    
	        ENDDO    
	        ENDDO    
	    
	        QBAR = ( QBAR)**(1.0_DP/mu)                
	        rr   = alpha* QBAR **(alpha-mu) * NBAR **(1.0_DP-alpha)
	        YBAR = QBAR ** alpha * NBAR **(1.0_DP-alpha)
	        wage = (1.0_DP-alpha)*QBAR **alpha * NBAR  **(-alpha)
	        Ebar = wage  * NBAR  * sum(pop)/sum(pop(1:MaxAge))
	    	! print*,'DBN_dist=',DBN_dist, 'QBAR=', QBAR ,  'NBAR=', NBAR 


	    	! Solve the model at current aggregate values
				! Find the threshold for wealth taxes (a_bar)
					!call Find_TauW_Threshold(DBN1,Y_a_threshold)
				! Adjust grid to include breaking points
					CALL Asset_Grid_Threshold(Y_a_threshold,agrid_t,na_t)
				! Compute labor units 
					CALL ComputeLaborUnits(Ebar, wage) 
				! Form YGRID for the capital income economy given interest rate "rr"
					CALL FORM_Y_MB_GRID(YGRID, MBGRID,YGRID_t,MBGRID_t)
				! Solve for policy and value functions 
					CALL EGM_RETIREMENT_WORKING_PERIOD 
	        
	     	! Discretize policy function for assets (a')
				! For each age and state vector bracket optimal a' between two grid points
				! When at that age and state the optimal decision is approximated by selecting one the grid points
				! The grid points are selected with probability proportional to their distance to the optimal a'
	        DO age=1,MaxAge
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


SUBROUTINE COMPUTE_STATS
	use global
	use programfunctions
	IMPLICIT NONE
	INTEGER :: prctile
	REAL(DP):: MeanReturn, StdReturn, VarReturn 
	REAL(DP), DIMENSION(nz) :: cdf_Gz_DBN 
	REAL(DP), DIMENSION(nz):: MeanReturn_by_z, size_by_z

	DO zi=1,nz
	    cdf_Gz_DBN(zi) = sum(DBN1(:,:,zi,:,:))
	ENDDO
	!print*,'cdf_Gz_DBN ='
	!print*,cdf_Gz_DBN

	DO ai=1,na
	     pr_a_dbn(ai)          = sum(DBN1(:,ai,:,:,:)) 
	     cdf_a_dbn(ai)         = sum( pr_a_dbn(1:ai) )      
	     tot_a_by_grid(ai)     = sum(DBN1(:,ai,:,:,:) * agrid(ai) )
	     cdf_tot_a_by_grid(ai) = sum(tot_a_by_grid(1:ai))   
	!     print*, pr_a_dbn(ai), cdf_a_dbn(ai)
	ENDDO
	cdf_a_dbn = cdf_a_dbn + 1.0_DP - cdf_a_dbn(na)

	!DO ai=1,na
	!     print*, pr_a_dbn(ai), cdf_a_dbn(ai)
	!ENDDO

	!print*,''

	! FIND THE ai THAT CORRESPONDS TO EACH PRCTILE OF WEALTH DBN & WEALTH HELD BY PEOPLE LOWER THAN THAT PRCTILE
	DO prctile=1,100
	    ai=1
	    DO while (cdf_a_dbn(ai) .lt. (REAL(prctile,8)/100.0_DP-0.000000000000001))
	        ai=ai+1
	    ENDDO
	    prctile_ai(prctile) = ai
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
		                 		& *  log( wage * yh(age, lambdai, ei) * Hours(age, ai, zi, lambdai, ei) )
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
		                 			& * ( log( wage * yh(age, lambdai, ei) * Hours(age, ai, zi, lambdai, ei) ) &
		                 			& -   mean_log_earnings_25_60 ) ** 2.0_DP
		ENDIF
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	Var_Log_Earnings_25_60 = Var_Log_Earnings_25_60 / pop_pos_earn_25_60
	Std_Log_Earnings_25_60 = Var_Log_Earnings_25_60 ** 0.5_DP

	MeanWealth = 0.0_DP
	MeanReturn = 0.0_DP
	DO age=1,MaxAge
	DO zi=1,nz
	DO ai=1,na
	DO lambdai=1,nlambda
	DO ei=1, ne
	     MeanWealth = MeanWealth  +   DBN1(age, ai, zi, lambdai, ei) * agrid(ai)         
	     MeanReturn = MeanReturn  +   DBN1(age, ai, zi, lambdai, ei) * (MBGRID(ai,zi)-1.0_DP)
	ENDDO
	ENDDO
	ENDDO    
	ENDDO    
	ENDDO    
	Wealth_Output = MeanWealth/YBAR 

	VarReturn = 0.0_DP
	DO age=1,MaxAge
	DO zi=1,nz
	DO ai=1,na
	DO lambdai=1,nlambda
	DO ei=1, ne      
	     VarReturn = VarReturn  +   DBN1(age, ai, zi, lambdai, ei) * (MBGRID(ai,zi)-1.0_DP-MeanReturn)**2.0_DP 
	ENDDO
	ENDDO
	ENDDO    
	ENDDO    
	ENDDO  
	StdReturn=VarReturn**0.5_DP

	MeanReturn_by_z=0.0_DP
	size_by_z = 0.0_DP
	DO age=1,MaxAge
	DO zi=1,nz
	DO ai=1,na
	DO lambdai=1,nlambda
	DO ei=1, ne
	     MeanReturn_by_z(zi) = MeanReturn_by_z(zi)  +   DBN1(age, ai, zi, lambdai, ei) * (MBGRID(ai,zi)-1.0_DP)
	     size_by_z(zi)       = size_by_z(zi) + DBN1(age, ai, zi, lambdai, ei) 
	ENDDO
	ENDDO
	ENDDO    
	ENDDO    
	ENDDO    
	MeanReturn_by_z = MeanReturn_by_z / size_by_z

	!print*, 'MeanReturn=',MeanReturn, 'StdReturn=', StdReturn
	!print*,'MeanReturn_by_z=',MeanReturn_by_z

	SSE_Moments = (Wealth_Output-3.0_DP)**2.0_DP + (prct1_wealth-0.34_DP)**2.0_DP  + (prct10_wealth-0.71_DP)**2.0_DP &
	                   & + (Std_Log_Earnings_25_60 -0.8_DP)**2.0_DP + (meanhours_25_60-0.4_DP)**2.0_DP
	!print*,''
	print*,"Current parameters"
	print*,'beta',beta,'rho_z',rho_z,'sigma_z',sigma_z_eps,'sigma_lam',sigma_lambda_eps,'phi',phi
	print*,"Statistics"
	print*,'W/GDP',Wealth_Output,'Top 1%',prct1_wealth,'Top 10%',prct10_wealth
	print*,'STD Labor Earnings',Std_Log_Earnings_25_60,'Mean Labor Earnings',meanhours_25_60,'Moments',SSE_Moments 
	!print*,''

END SUBROUTINE COMPUTE_STATS


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE WRITE_VARIABLES(bench_indx)
	USE GLOBAL
	use  programfunctions
	IMPLICIT NONE
	integer :: bench_indx,  prctile
		!IF (bench_indx .gt. 0) then 
		!    OPEN (UNIT=2, FILE='cons_bench', STATUS='replace')
		!    OPEN (UNIT=3, FILE='aprime_bench', STATUS='replace')
		!    OPEN (UNIT=4, FILE='agrid_bench', STATUS='replace')
		!    OPEN (UNIT=6, FILE='hours_bench', STATUS='replace')
		!    OPEN (UNIT=7, FILE='Value_Bench', STATUS='replace')    
		!    OPEN   (UNIT=8, FILE='DBN_bench', STATUS='replace')
		!    OPEN   (UNIT=9, FILE='QBAR_bench', STATUS='replace')
		!    OPEN   (UNIT=10, FILE='NBAR_bench', STATUS='replace')
		!ELSE
		!    OPEN (UNIT=2, FILE='cons_exp', STATUS='replace')
		!    OPEN (UNIT=3, FILE='aprime_exp', STATUS='replace')
		!    OPEN (UNIT=4, FILE='agrid_exp', STATUS='replace')
		!    OPEN (UNIT=6, FILE='hours_exp', STATUS='replace')
		!    OPEN (UNIT=7, FILE='Value_exp', STATUS='replace')    
		!    OPEN   (UNIT=8, FILE='DBN_exp', STATUS='replace')
		!    OPEN   (UNIT=9, FILE='QBAR_exp', STATUS='replace')
		!    OPEN   (UNIT=10, FILE='NBAR_exp', STATUS='replace')
		!    OPEN   (UNIT=11, FILE='Cons_Eq_Welfare', STATUS='replace')
		!    WRITE(unit=11, FMT=*) Cons_Eq_Welfare
		!    CLOSE (unit=11)
		!ENDIF
		!
		!WRITE  (UNIT=2, FMT=*) cons
		!WRITE  (UNIT=3, FMT=*) Aprime
		!WRITE  (UNIT=4, FMT=*) agrid
		!WRITE  (UNIT=6, FMT=*) Hours
		!WRITE  (UNIT=7, FMT=*) ValueFunction
		!WRITE  (UNIT=8, FMT=*) DBN1
		!WRITE  (UNIT=9, FMT=*) QBAR 
		!WRITE  (UNIT=10, FMT=*) NBAR
		!
		!close (unit=2)
		!close (unit=3)
		!close (unit=4)
		!close (unit=6)
		!close (unit=7)
		!close (unit=8)
		!close (unit=9)
		!close (unit=10)
		!
		!OPEN   (UNIT=2, FILE='pr_a_dbn', STATUS='replace')
		!DO ai=1,na
		!     WRITE(unit=2, FMT=*) pr_a_dbn(ai), cdf_a_dbn(ai), tot_a_by_grid(ai), cdf_tot_a_by_grid(ai)
		!ENDDO
		!CLOSE (unit=2)
		!
		!OPEN   (UNIT=2, FILE='prctile_a_dbn', STATUS='replace')
		!DO prctile=1,100
		!    WRITE(unit=2, FMT=*) prctile, cdf_a_dbn(prctile_ai(prctile)),   prctile_ai(prctile),  cdf_tot_a_by_prctile(prctile)
		!ENDDO
		!CLOSE (unit=2)

	! If benchmark economy then write as follows
	IF (bench_indx .gt. 0) then 
		OPEN   (UNIT=2, FILE='output_bench.txt', STATUS='replace')

		WRITE  (UNIT=2, FMT=*)  'Params=[', params,']'
		WRITE  (UNIT=2, FMT=*)  'GBAR_bench=', GBAR_bench
		WRITE  (UNIT=2, FMT=*)  'QBAR_bench=',QBAR_bench
		WRITE  (UNIT=2, FMT=*)  'NBAR_bench=',NBAR_bench
		WRITE  (UNIT=2, FMT=*)  'EBAR_bench=',EBAR_bench
		WRITE  (UNIT=2, FMT=*)  'rr_bench=',rr_bench
		WRITE  (UNIT=2, FMT=*)  'wage_bench=',wage_bench
		WRITE  (UNIT=2, FMT=*)  'Y_bench=',Y_bench
		WRITE  (UNIT=2, FMT=*)  'MOMENTS'
		WRITE  (UNIT=2, FMT=*)  Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60
	! If experimental economy then write as follows
	else

		OPEN   (UNIT=2, FILE='output_exp.txt', STATUS='replace')

		WRITE  (UNIT=2, FMT=*)  'GBAR_exp=', GBAR_exp
		WRITE  (UNIT=2, FMT=*)  'QBAR_exp=',QBAR_exp
		WRITE  (UNIT=2, FMT=*)  'NBAR_exp=',NBAR_exp
		WRITE  (UNIT=2, FMT=*)  'EBAR_exp=',EBAR_exp
		WRITE  (UNIT=2, FMT=*)  'rr_exp=',rr_exp
		WRITE  (UNIT=2, FMT=*)  'wage_exp=',wage_exp
		WRITE  (UNIT=2, FMT=*)  'Y_exp=',Y_exp
		WRITE  (UNIT=2, FMT=*)  ''
		WRITE  (UNIT=2, FMT=*)  'Y_exp_Y_bench=',Y_exp/Y_bench
		WRITE  (UNIT=2, FMT=*)  'tauW_bt =',tauW_bt_exp
		WRITE  (UNIT=2, FMT=*)  'tauW_at =',tauW_at_exp
		WRITE  (UNIT=2, FMT=*)  'Y_a_threshold =',Y_a_threshold_exp
		WRITE  (UNIT=2, FMT=*)  ''
		WRITE  (UNIT=2, FMT=*)  'MOMENTS'
		WRITE  (UNIT=2, FMT=*)  Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60
	ENDIF

	close (unit=2)

END SUBROUTINE WRITE_VARIABLES


!========================================================================================
!========================================================================================
!========================================================================================

SUBROUTINE EGM_RETIREMENT_WORKING_PERIOD
	use global
	use programfunctions
	IMPLICIT NONE
	!REAL(DP), DIMENSION(fine_na, nz) :: FineYGRID
	REAL(DP), DIMENSION(MaxAge, na, nz, nlambda) :: CorrectRetValueP1, CorrectRetValueP2
	REAL(DP) :: AnRetCCtemp, cprime, AntRetVVtemp, disty, slopeV, slopeC, linearC, linearV
	REAL(DP) :: ap_temp, brentvalue
	REAL(DP) :: tempvar1, tempvar2, tempvar3
	INTEGER  :: na1, na2, tempai
	REAL(DP), DIMENSION(na_t+1) :: EndoCons, EndoYgrid, EndoHours
	INTEGER , DIMENSION(na_t+1) :: sort_ind 
	INTEGER  :: sw 
	REAL(DP) :: Wealth

	! These lines are note being used!!!!!!!!
		! Set auxiliary value ofr psi 
			psi = psi_PL * Ebar ** tauPL
		! If solving experimental economy use Ebar_bench
			if (solving_bench .eq. 0) then
	            psi = psi_PL * Ebar_bench ** tauPL
			endif          
	! These lines are note being used!!!!!!!!      
	
	! Set psi to psi_PL
		psi = psi_PL 

	! SET HOURS TO ZERO SO THAT RETIRED HAS ZERO HOURS
	Hours_t = 0.0_DP

	! Last period of life
	age=MaxAge
	DO lambdai=1,nlambda
	DO zi=1,nz
    DO ei=1,ne
    DO ai=1,na_t
    	! Find hours under no savings 
    	ain = 0.0_dp         
		brentvalue = brent(0.000001_DP, 0.4_DP, 0.99_DP, FOC_HA, brent_tol, Hours_t(age, ai, zi, lambdai,ei))           
        ! Consumption given hours worked
        Cons_t(age,ai,zi,lambdai,ei) = YGRID_t(ai,zi)+psi*(yh(age,lambdai,ei)* Hours_t(age, ai, zi, lambdai,ei))**(1.0_DP-tauPL)
	ENDDO ! ai
    ENDDO ! ei
    ENDDO ! zi
	ENDDO ! lambdai
	Aprime_t(age, :, :, :, :) = 0.0_DP

	
	DO age=MaxAge-1,1,-1
    DO lambdai=1,nlambda
    DO zi=1,nz
    DO ei=1,ne	
    	! Endogenous grid and consumption are initialized
	    EndoCons  = big_p 
	    EndoYgrid = big_p 
	    EndoHours = big_p 
	    sw 		  = 0                    
    DO ai=1,na_t
    	Wealth = agrid_t(ai)+(rr*(zgrid(zi)*agrid_t(ai))**mu-DepRate*agrid_t(ai))*(1.0_DP-tauK)
		if (abs(Wealth-Y_a_threshold).lt.1e-8) then 
	    	! Consumption on endogenous grid and implied asset income under tauW_bt
			EndoCons(ai) = 1.0_DP/( beta*survP(age)*MB_a_bt(agrid_t(ai),zgrid(zi))*SUM(pr_e(ei,:)/Cons_t(age+1,ai,zi,lambdai,:)))    
			    
			! Auxiliary consumption variable for FOC_H        
			consin =  EndoCons(ai)     
			! Solution of Labor FOC for hours
			brentvalue = brent(0.000001_DP, 0.4_DP, 0.99_DP, FOC_H, brent_tol, EndoHours(ai) )           
			
			EndoYgrid(ai) = agrid_t(ai) + EndoCons(ai) - psi*(yh(age, lambdai,ei)* EndoHours(ai))**(1.0_DP-tauPL)

	    	! Consumption, hours and asset income in endogenous grid with above threshold tax
	    	EndoCons(na_t+1) = 1.0_DP/( beta*survP(age)*MB_a_at(agrid_t(ai),zgrid(zi))*SUM(pr_e(ei,:)/Cons_t(age+1,ai,zi,lambdai,:)))   
	    	! Auxiliary consumption variable for FOC_H        
			consin =  EndoCons(na_t+1)     
			! Solution of Labor FOC for hours
			brentvalue = brent(0.000001_DP, 0.4_DP, 0.99_DP, FOC_H, brent_tol, EndoHours(na_t+1) )           
			
			EndoYgrid(na_t+1) = agrid_t(ai) + EndoCons(na_t+1) - psi*(yh(age, lambdai,ei)* EndoHours(na_t+1))**(1.0_DP-tauPL)
	    	sw                = 1 
	    else 
	    	! Consumption, hours and asset income in endogenous grid
			EndoCons(ai) =   1.0_DP /( beta*survP(age)*MBGRID_t(ai,zi)*SUM(pr_e(ei,:) /Cons_t(age+1,ai,zi,lambdai,:)))    
			    
			! Auxiliary consumption variable for FOC_H        
			consin =  EndoCons(ai)     
			! Solution of Labor FOC for hours
			brentvalue = brent(0.000001_DP, 0.4_DP, 0.99_DP, FOC_H, brent_tol, EndoHours(ai) )           
			
			EndoYgrid(ai) = agrid_t(ai) + EndoCons(ai) - psi*(yh(age, lambdai,ei)* EndoHours(ai))**(1.0_DP-tauPL)
	    end if 

    ENDDO ! ai  

    ! Sort endogenous grid for interpolation
	call Sort(na_t+1,EndoYgrid,EndoYgrid,sort_ind)
	EndoCons = EndoCons(sort_ind)

! 	print*, "Y_grid"
! 	do ai=1,na
! 		print*, YGRID(ai,zi),YGRID_t(ai,zi),EndoYgrid(ai),EndoCons(ai)
! 	end do

    	! Find  decision rules on exogenous grids
        tempai=1           
        DO WHILE ( YGRID_t(tempai,zi) .lt. EndoYgrid(1) )
              tempai = tempai + 1
        ENDDO
       
		! decision rules are obtained taking care of extrapolations
		DO ai=tempai,na_t               
		    Cons_t(age, ai, zi, lambdai,ei)= Linear_Int(EndoYgrid(1:na_t+sw), EndoCons(1:na_t+sw),na_t+sw, YGRID_t(ai,zi))  
		     
		    consin = Cons_t(age, ai, zi, lambdai,ei)

		    brentvalue = brent(0.000001_DP, 0.4_DP, 0.99_DP, FOC_H, brent_tol, Hours_t(age, ai, zi, lambdai,ei))           
		    
		    Aprime_t(age, ai, zi, lambdai,ei) = YGRID_t(ai,zi)  &
		                    & + psi*(yh(age, lambdai,ei)*Hours_t(age, ai, zi, lambdai,ei))**(1.0_DP-tauPL)  & 
		                    & - Cons_t(age, ai, zi, lambdai,ei)   
		                    
		    If (Aprime_t(age, ai, zi, lambdai,ei)  .lt. amin) then
            	Aprime_t(age, ai, zi, lambdai,ei) = amin
		         
	           	!compute  hours using FOC_HA                              
		        ain = amin
		           
		        brentvalue = brent(0.000001_DP, 0.4_DP, 0.99_DP, FOC_HA, brent_tol, Hours_t(age, ai, zi, lambdai,ei))           
		        
	            Cons_t(age,ai,zi,lambdai,ei) = YGRID_t(ai,zi)+psi*(yh(age, lambdai,ei)&
		                            		& * Hours_t(age, ai, zi, lambdai,ei))**(1.0_DP-tauPL) &
		                            		& - Aprime_t(age, ai, zi, lambdai,ei)      

				IF (Cons_t(age, ai, zi, lambdai,ei) .le. 0.0_DP)  THEN
					print*,'w1: Cons(age, ai, zi, lambdai,ei)=',Cons_t(age, ai, zi, lambdai,ei)
				ENDIF                   
		     endif      
		ENDDO ! ai                  

		ai=1           
        DO WHILE ( YGRID_t(ai,zi) .lt. EndoYgrid(1) )
	        ! Solve for the Euler equation directly         
			brentvalue = brent( min(amin,YGRID_t(ai,zi))   ,  (amin+YGRID_t(ai,zi))/2.0_DP  ,  &
                             	& YGRID_t(ai,zi)  + psi * ( yh(age, lambdai,ei)*0.95_DP )**(1.0_DP-tauPL)  ,  &
	                             & FOC_WH, brent_tol, Aprime_t(age, ai, zi, lambdai,ei) )
			
			!compute  hours using FOC_HA                              
			ain = Aprime_t(age, ai, zi, lambdai,ei)
	                       
            brentvalue = brent(0.000001_DP, 0.4_DP, 0.99_DP, FOC_HA, brent_tol, Hours_t(age, ai, zi, lambdai,ei))           
	 		
            Cons_t(age, ai, zi, lambdai,ei)=  YGRID_t(ai,zi)  &
                            & + psi*(yh(age, lambdai,ei)*Hours_t(age, ai, zi, lambdai,ei))**(1.0_DP-tauPL)  &
                            & - Aprime_t(age, ai, zi, lambdai,ei)
           IF (Cons_t(age, ai, zi, lambdai,ei) .le. 0.0_DP)  THEN
                print*,'w2:',age,zi,lambdai,ei,ai, 'Cons(age, ai, zi, lambdai,ei)=',Cons_t(age, ai, zi, lambdai,ei), &
                    & 'Aprime(age, ai, zi, lambdai,ei)=',Aprime_t(age, ai, zi, lambdai,ei), &
                    & 'Hours(age, ai, zi, lambdai,ei)=', Hours_t(age, ai, zi, lambdai,ei), &
                    & 'yh(age, lambdai,ei)=', yh(age, lambdai,ei),'YGRID(ai,zi)=',YGRID_t(ai,zi)
                !pause
                print*, "there is a problem in line 2063"
            ENDIF                   
            ai = ai + 1
        ENDDO  
	                 
    ENDDO !ei         
    ENDDO !zi
    ENDDO !lambdai
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
	cons = cons/(1.0_DP+tauC)

END SUBROUTINE EGM_RETIREMENT_WORKING_PERIOD

!========================================================================================
!========================================================================================
!========================================================================================

Subroutine Find_TauW_Threshold(DBN_in,Y_a_threshold_out)
	use NRTYPE
	use parameters
	use global
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
		Wealth = spread(agrid,2,nz)+(rr*(spread(zgrid,1,na)*spread(agrid,2,nz))**mu-DepRate*spread(agrid,2,nz))*(1.0_DP-tauK)
		! Mean Wealth
		Mean_Wealth = sum( Wealth*DBN_az )
		! Set threshold
		Y_a_threshold_out = 0.75_dp*Mean_Wealth
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
			! 		    prctile_ai(prctile) = a_ind
			! 		ENDDO

			! 	! Set threshold for wealth tax as:
			! 		if (tauW_threshold_in.eq.0) then 
			! 			! Threshold is minimum
			! 			Y_a_threshold_out = 0.0_dp
			! 		else 
			! 			! Threshold at some percentile
			! 			Y_a_threshold_out = agrid( prctile_ai(tauW_threshold_in) )
			! 		end if 

end Subroutine Find_TauW_Threshold


!========================================================================================
!========================================================================================
!========================================================================================
! THIS YGRID and Marginal Benefit of Investment GRID 
! NEEDS TO BE COMPUTED FOR EACH TIME THE INTEREST RATE "rr" IS UPDATED. 

SUBROUTINE FORM_Y_MB_GRID(TYGRID, TMBGRID,TYGRID_t,TMBGRID_t)
	USE global
	USE programfunctions
	IMPLICIT NONE
	REAL(DP), DIMENSION(na,nz),   INTENT(OUT) :: TYGRID, TMBGRID
	REAL(DP), DIMENSION(na_t,nz), INTENT(OUT) :: TYGRID_t, TMBGRID_t
	!REAL(DP), INTENT(IN) :: rr
	!integer :: ai, zi

	DO zi=1,nz
		DO ai=1,na
		! Asset income grid (by current asset -ai- and entrepreneurial ability -zi-)
		TYGRID(ai,zi)  = Y_a(agrid(ai),zgrid(zi))
		! Asset marginal benefit grid (by current asset -ai- and entrepreneurial ability -zi-)
	    TMBGRID(ai,zi) = MB_a(agrid(ai),zgrid(zi))
		ENDDO
		DO ai=1,na_t
		! Asset income grid (by current asset -ai- and entrepreneurial ability -zi-)
		TYGRID_t(ai,zi)  = Y_a(agrid_t(ai),zgrid(zi))
		! Asset marginal benefit grid (by current asset -ai- and entrepreneurial ability -zi-)
	    TMBGRID_t(ai,zi) = MB_a(agrid_t(ai),zgrid(zi))
		ENDDO
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
	use global
	IMPLICIT NONE
	!integer:: lambdai  
	REAL(DP), INTENT(IN) :: Ebart, Waget 

	! This computes efficiency units times wage
	yh = Waget * eff_un

END SUBROUTINE ComputeLaborUnits


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE  INITIALIZE
	use parameters
	use global
	use programfunctions
	IMPLICIT NONE

	!integer::  lambdai
	REAL(DP) :: lambdaBAR, tempno 
	REAL(DP) :: m, Rh, start_timet, finish_timet
	INTEGER  :: ee0, ee1, ee2, zindx1, zindx2, lambdaindx1, lambdaindx2, diff_array, eindx1, eindx2
	INTEGER, DIMENSION(MaxAge) :: agevec
	!INTEGER :: ageij, ageii, eij, eii, zij, zii, lambdaij, lambdaii
	!integer, dimension(ne,ne,nz,nz,nlambda,nlambda,MaxAge,MaxAge) :: globaltransition
	!
	!DO eij=1,ne
	!DO eii=1,ne
	!DO zij=1,nz
	!DO zii=1,nz
	!DO lambdaij=1,nlambda
	!DO lambdaii=1,nlambda
	!DO ageij=1,MaxAge
	!DO ageii=1,MaxAge    
	!    globaltransition(ei,eii,zi,zii,lambdai,lambdaii,agei,ageii) = 1
	!ENDDO    
	!ENDDO   
	!ENDDO   
	!ENDDO   
	!ENDDO   
	!ENDDO   
	!ENDDO   
	!ENDDO   
	
	! Initiliaze grids for z, lamda and e	
		CALL tauchen(mtauchen,rho_z,sigma_z_eps,nz,zgrid,pr_z,Gz)
		CALL tauchen(mtauchen,rho_E,sigma_e_eps,ne,egrid,pr_e,Ge)
		CALL tauchen(mtauchen,rho_lambda,sigma_lambda_eps,nlambda,lambdagrid,pr_lambda,Glambda)

		! Tauchen gives grids for the log of the variables. Exponentiate to adjust
		zgrid      = exp(zgrid)
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
		DO age=1,MaxAge
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

	! Initialize asset grid
		! Normal grid
		m=(amax-amin)**(1.0_DP/a_theta)/REAL(na-1,DP)
		DO ai=1,na
			agrid(ai)=REAL(ai-1,DP)*m
		END DO
		agrid=amin+agrid**a_theta
		!	print*,'agrid=',agrid
		!!pause

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
! 		pop(41)=175009.0_DP
! 		pop(42)=173187.0_DP
! 		pop(43)=171214.0_DP
! 		pop(44)=169064.0_DP
! 		pop(45)=166714.0_DP
! 		pop(46)=164147.0_DP
! 		pop(47)=161343.0_DP
! 		pop(48)=158304.0_DP
! 		pop(49)=155048.0_DP
! 		pop(50)=151604.0_DP
! 		pop(51)=147990.0_DP
! 		pop(52)=144189.0_DP
! 		pop(53)=140180.0_DP
! 		pop(54)=135960.0_DP
! 		pop(55)=131532.0_DP
! 		pop(56)=126888.0_DP
! 		pop(57)=122012.0_DP
! 		pop(58)=116888.0_DP
! 		pop(59)=111506.0_DP
! 		pop(60)=105861.0_DP
! 		pop(61)=99957.0_DP
! 		pop(62)=93806.0_DP
! 		pop(63)=87434.0_DP
! 		pop(64)=80882.0_DP
! 		pop(65)=74204.0_DP
! 		pop(66)=67462.0_DP
! 		pop(67)=60721.0_DP
! 		pop(68)=54053.0_DP
! 		pop(69)=47533.0_DP
! 		pop(70)=41241.0_DP
! 		pop(71)=35259.0_DP
! 		pop(72)=29663.0_DP
! 		pop(73)=24522.0_DP
! 		pop(74)=19890.0_DP
! 		pop(75)=15805.0_DP
! 		pop(76)=12284.0_DP
! 		pop(77)=9331.0_DP
! 		pop(78)=6924.0_DP
! 		pop(79)=5016.0_DP
! 		pop(80)=3550.0_DP

		if (MaxAge .gt. 80)	 then	
			pop(MaxAge)=2454.0_DP
		endif	
		! pop = 1.0_DP	

	
	! Survival probabilities: surv(i)=prob(alive in i+1|alive in i)
		FORALL (age=1:maxAge-1) survP(age)= pop(age+1)/pop(age)
		survP(maxAge)=0.0_DP
			


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

END SUBROUTINE INITIALIZE


!========================================================================================
!========================================================================================
!========================================================================================


SUBROUTINE Write_Benchmark_Results(read_write)
	use global
	IMPLICIT NONE
	integer :: read_write
	
	IF (read_write .eq. 0) then 
		OPEN  (UNIT=1,  FILE='Simple_Bench/simple_results_cons'  , STATUS='replace')
		WRITE (UNIT=1,  FMT=*) cons
		CLOSE (unit=1)
		OPEN  (UNIT=2,  FILE='Simple_Bench/simple_results_aprime', STATUS='replace')
		WRITE (UNIT=2,  FMT=*) aprime
		CLOSE (unit=2)
		OPEN  (UNIT=3,  FILE='Simple_Bench/simple_results_hours' , STATUS='replace')
		WRITE (UNIT=3,  FMT=*) hours
		CLOSE (unit=3)
		OPEN  (UNIT=4,  FILE='Simple_Bench/simple_results_value' , STATUS='replace')
		WRITE (UNIT=4,  FMT=*) ValueFunction
		CLOSE (unit=4)

		OPEN  (UNIT=5,  FILE='Simple_Bench/simple_results_DBN'   , STATUS='replace')
		WRITE (UNIT=5,  FMT=*) DBN1 
		CLOSE (UNIT=5)
		OPEN  (UNIT=60,  FILE='Simple_Bench/simple_results_GBAR'  , STATUS='replace')
		WRITE (UNIT=60,  FMT=*) GBAR
		CLOSE (UNIT=60)
		OPEN  (UNIT=7,  FILE='Simple_Bench/simple_results_EBAR'  , STATUS='replace')
		WRITE (UNIT=7,  FMT=*) EBAR
		CLOSE (UNIT=7)
		OPEN  (UNIT=8,  FILE='Simple_Bench/simple_results_NBAR'  , STATUS='replace')
		WRITE (UNIT=8,  FMT=*) NBAR
		CLOSE (UNIT=8)
		OPEN  (UNIT=9,  FILE='Simple_Bench/simple_results_QBAR'  , STATUS='replace')
		WRITE (UNIT=9,  FMT=*) QBAR
		CLOSE (UNIT=9)
		OPEN  (UNIT=10, FILE='Simple_Bench/simple_results_rr'    , STATUS='replace')
		WRITE (UNIT=10, FMT=*) rr
		CLOSE (UNIT=10)
		OPEN  (UNIT=11, FILE='Simple_Bench/simple_results_wage'  , STATUS='replace')
		WRITE (UNIT=11, FMT=*) wage 
		CLOSE (UNIT=11)
		OPEN  (UNIT=12, FILE='Simple_Bench/simple_results_YBAR'  , STATUS='replace')
		WRITE (UNIT=12, FMT=*) YBAR
		CLOSE (UNIT=12)

		print*, "Writing of benchmark results completed"
	ELSE 
		OPEN (UNIT=1,  FILE='Simple_Bench/simple_results_cons'  , STATUS='old', ACTION='read')
		OPEN (UNIT=2,  FILE='Simple_Bench/simple_results_aprime', STATUS='old', ACTION='read')
		OPEN (UNIT=3,  FILE='Simple_Bench/simple_results_hours' , STATUS='old', ACTION='read')
		OPEN (UNIT=4,  FILE='Simple_Bench/simple_results_value' , STATUS='old', ACTION='read')
		OPEN (UNIT=5,  FILE='Simple_Bench/simple_results_DBN'   , STATUS='old', ACTION='read')
		OPEN (UNIT=60, FILE='Simple_Bench/simple_results_GBAR'  , STATUS='old', ACTION='read')
		OPEN (UNIT=7,  FILE='Simple_Bench/simple_results_EBAR'  , STATUS='old', ACTION='read')
		OPEN (UNIT=8,  FILE='Simple_Bench/simple_results_NBAR'  , STATUS='old', ACTION='read')
		OPEN (UNIT=9,  FILE='Simple_Bench/simple_results_QBAR'  , STATUS='old', ACTION='read')
		OPEN (UNIT=10, FILE='Simple_Bench/simple_results_rr'    , STATUS='old', ACTION='read')
		OPEN (UNIT=11, FILE='Simple_Bench/simple_results_wage'  , STATUS='old', ACTION='read')
		OPEN (UNIT=12, FILE='Simple_Bench/simple_results_YBAR'  , STATUS='old', ACTION='read')

		READ (UNIT=1,  FMT=*), cons
		READ (UNIT=2,  FMT=*), aprime
		READ (UNIT=3,  FMT=*), hours
		READ (UNIT=4,  FMT=*), ValueFunction
		READ (UNIT=5,  FMT=*), DBN1 
		READ (UNIT=60, FMT=*), GBAR 
		READ (UNIT=7,  FMT=*), EBAR
		READ (UNIT=8,  FMT=*), NBAR
		READ (UNIT=9,  FMT=*), QBAR
		READ (UNIT=10, FMT=*), rr
		READ (UNIT=11, FMT=*), wage 
		READ (UNIT=12, FMT=*), YBAR

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

		print*, "Reading of benchmark results completed"
	END IF 
END SUBROUTINE Write_Benchmark_Results

!========================================================================================
!========================================================================================
!========================================================================================

SUBROUTINE tauchen(mt,rhot,sigmat,nt,gridt,prt,Gt)
	USE parameters
	IMPLICIT NONE
	
	REAL(DP), INTENT(IN) :: rhot, sigmat
	INTEGER(I4B),  INTENT(IN) :: nt
	REAL(DP), INTENT(OUT), DIMENSION(nt)    :: gridt, Gt
	REAL(DP), INTENT(OUT), DIMENSION(nt,nt) :: prt
	REAL(DP), DIMENSION(nt)    :: Gt_new
	REAL(DP) :: a, stept, cdf_normal, mut
	REAL(DP), INTENT(IN) ::  mt
	INTEGER(I4B)  :: i, j, k, zi, zz



	mut=0.0_DP
	gridt=0.0_DP
	prt=0.0_DP
	Gt = 0.0_DP
	a=(1.0_DP-rhot)*mut;
        
    if (nt .gt. 1) then  
		gridt(nt)=mt*sqrt(sigmat**2.0_DP/(1.0_DP-rhot**2))
		gridt(1)=-gridt(nt)

		stept=(gridt(nt)-gridt(1))/REAL(nt-1,DP)
		DO i=2,nt-1
			gridt(i)=gridt(1)+stept*REAL(i-1,DP)
		END DO
		gridt=gridt+a/(1.0_DP-rhot)
		
		DO j=1,nt
			DO k=1,nt
				IF (k==1) THEN
					prt(j,k)=cdf_normal((gridt(1)-a-rhot*gridt(j)+stept/2.0_DP)/sigmat)
				ELSE IF (k==nt) THEN
					prt(j,k)=1.0_DP-cdf_normal((gridt(nt)-a-rhot*gridt(j)-stept/2.0_DP)/sigmat)
				ELSE
	                prt(j,k)=cdf_normal((gridt(k)-a-rhot*gridt(j)+stept/2.0_DP)/sigmat)- &
						& cdf_normal((gridt(k)-a-rhot*gridt(j)-stept/2.0_DP)/sigmat)
				END IF
			END DO
		END DO
		
		Gt(1)=cdf_normal((gridt(1)+stept/2.0_DP)/sigmat)
		DO zi=2,nt-1
			Gt(zi)=cdf_normal((gridt(zi)+stept/2.0_DP)/sigmat)- &
						& cdf_normal((gridt(zi)-stept/2.0_DP)/sigmat)
		END DO
		Gt(nt)=1.0_DP-cdf_normal((gridt(nt)-stept/2.0_DP)/sigmat)
	 	! print*, 'Gt', Gt, 'sum', sum(Gt)

		DO i=1,1000
			Gt_new=0.0_DP
			DO zi=1,nt
				DO zz=1,nt
					Gt_new(zz)=Gt_new(zz)+Gt(zi)*prt(zi,zz)
				END DO
			END DO
			Gt=Gt_new
		END DO
		
	 	! print*, 'Gt', Gt, 'sum', sum(Gt)
	 	! pause

  	ELSE
  		prt = 1.0_DP
		Gt  = 1.0_DP
	endif       
	
END SUBROUTINE tauchen 

!====================================================================

REAL(DP) FUNCTION cdf_normal(x)
	USE parameters
	IMPLICIT NONE

	real(DP), parameter :: a1 = 0.398942280444D+00
	real(DP), parameter :: a2 = 0.399903438504D+00
	real(DP), parameter :: a3 = 5.75885480458D+00
	real(DP), parameter :: a4 = 29.8213557808D+00
	real(DP), parameter :: a5 = 2.62433121679D+00
	real(DP), parameter :: a6 = 48.6959930692D+00
	real(DP), parameter :: a7 = 5.92885724438D+00
	real(DP), parameter :: b0 = 0.398942280385D+00
	real(DP), parameter :: b1 = 3.8052D-08
	real(DP), parameter :: b2 = 1.00000615302D+00
	real(DP), parameter :: b3 = 3.98064794D-04
	real(DP), parameter :: b4 = 1.98615381364D+00
	real(DP), parameter :: b5 = 0.151679116635D+00
	real(DP), parameter :: b6 = 5.29330324926D+00
	real(DP), parameter :: b7 = 4.8385912808D+00
	real(DP), parameter :: b8 = 15.1508972451D+00
	real(DP), parameter :: b9 = 0.742380924027D+00
	real(DP), parameter :: b10 = 30.789933034D+00
	real(DP), parameter :: b11 = 3.99019417011D+00
	real(DP) q
	real(DP) x
	real(DP) y

	!  |X| <= 1.28.
	!
	  if ( abs ( x ) <= 1.28D+00 ) then
	
	    y = 0.5D+00 * x * x
	
	    q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
	      + a6 / ( y + a7 ) ) ) )
	!
	!  1.28 < |X| <= 12.7
	!
	  else if ( abs ( x ) <= 12.7D+00 ) then
	
	    y = 0.5D+00 * x * x
	
	    q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
	      + b2 / ( abs ( x ) + b3 &
	      + b4 / ( abs ( x ) - b5 &
	      + b6 / ( abs ( x ) + b7 &
	      - b8 / ( abs ( x ) + b9 &
	      + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
	!
	!  12.7 < |X|
	!
	  else
	
	    q = 0.0D+00
	
	  end if
	!
	!  Take account of negative X.
	!
	if ( x < 0.0D+00 ) then
	    cdf_normal = q
	  else
	    cdf_normal = 1.0D+00 - q
	  end if
	
	  return
	end


!=======================================================================


SUBROUTINE spline(x,y,n,yp1,ypn,y2)  
	! Given arrays x(1:n) and y(1:n) of length n; and first derivatives yp1 and ypn at points 1 and n, 
	! this subroutine returns an array y2(1:n) (2nd derivative) at points x(1:n) parameter NMAX is the largest anticipated value of "n"
	 
	USE parameters
	IMPLICIT NONE      

	INTEGER, PARAMETER :: NMAX=1000
	INTEGER  :: n 
	INTEGER  :: i,k  
	REAL(DP) ::  yp1,ypn,x(n),y(n),y2(n), y1(n)  
	REAL(DP) :: p,qn,sig,un,u(NMAX)  
	  
	   
	if (yp1.gt..99e30) then  
		y2(1)=0.0_DP  
		u(1)=0.0_DP  
	else  
		y2(1)=-0.5_DP  
		u(1)=(3.0_DP/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)  
	endif  

	do i=2,n-1  
	    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))  
	    p=sig*y2(i-1)+2.0_DP  
	    y2(i)=(sig-1.0_DP)/p  
	    u(i)=(6.0_DP*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p  
	enddo
	 
	if (ypn .gt. .99e30) then  
		qn=0.0_DP  
		un=0.0_DP  
	else  
		qn=0.50_DP  
		un=(3.0_DP/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))  
	endif  
	      
	y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0_DP)  
	      
	do k=n-1,1,-1  
        y2(k)=y2(k)*y2(k+1)+u(k)   
	enddo
	 
	! If you want y1(1:n) (1st derivative) at points x(1:n), then do the following
	!do k=1,n-1 
	!        y1(k)= (y(k+1)-y(k))/(x(k+1)-x(k)) - (x(k+1)-x(k)) * y2(k)/3.0_DP - (x(k+1)-x(k))*y2(k+1)/6.0_DP
	!enddo
	!
	!y1(n) = (y(n)-y(n-1))/(x(n)-x(n-1)) + (x(n)-x(n-1)) * y2(n)/3.0_DP + (x(n)-x(n-1))*y2(n-1)/6.0_DP

	!print*,'yp1',yp1,y1(1)
	!print*,'ypn',ypn,y1(n)     

	return  
  
END SUBROUTINE spline
      
!=======================================================================
       
 
SUBROUTINE splint(xa,ya,y2a,n,x,y)  
	! Given arrays xa(1:n) and ya(1:n) of length n; and second derivatives y2a(1:n), which is the output of spline
	! given the value of "x" this subroutine returns  a cubic-spline interpolated value "y" at point "x"

	USE parameters
	IMPLICIT NONE      
	INTEGER  :: n  
	INTEGER  :: k,khi,klo  
	REAL(DP) :: x, y, yprime, xa(n),y2a(n),ya(n)  
	REAL(DP) :: a,b,h  
      
	!      y2a=0.0_DP
	      
	klo=1  
	khi=n  

1   if (khi-klo.gt.1) then  
    	k=(khi+klo)/2  
        if(xa(k).gt.x)then  
        	khi=k  
        else  
        	klo=k  
        endif  
      	
      	goto 1  
    
    endif
  
	h = xa(khi)-xa(klo)  
	
	if (h.eq.0.) then 
		print*,'bad xa input in splint'  
	end if
	
	a = (xa(khi)-x)/h  
	b = (x-xa(klo))/h  
	y = a*ya(klo)+b*ya(khi)     +((a**3.0_DP-a)*y2a(klo)+(b**3.0_DP-b)*y2a(khi))*(h**2.0_DP)/6.0_DP  
	! it also returns "yprime" (first derivative) at point "x"
	!      yprime = ( ya(khi) - ya(klo) ) / h  & 
	!                 & + h * ( ( 1.0_DP- 3.0_DP * a**2.0_DP  ) *  y2a(klo)  +  ( 3.0_DP * b**2.0_DP -1.0_DP ) *  y2a(khi)  )/6.0_DP
	return  

END  
 

!===========================================================================

FUNCTION zbrent(func,x1,x2,tol)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x1,x2,tol
	REAL(DP) :: zbrent
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(DP), PARAMETER :: EPS=epsilon(x1)
	INTEGER(I4B) :: iter
	REAL(DP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
	a=x1
	b=x2
	fa=func(a)
	fb=func(b)
	if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) &
		call nrerror('root must be bracketed for zbrent')
	c=b
	fc=fb
	do iter=1,ITMAX
		if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
			c=a
			fc=fa
			d=b-a
			e=d
		end if
		if (abs(fc) < abs(fb)) then
			a=b
			b=c
			c=a
			fa=fb
			fb=fc
			fc=fa
		end if
		tol1=2.0_DP*EPS*abs(b)+0.5_DP*tol
		xm=0.5_DP*(c-b)
		if (abs(xm) <= tol1 .or. fb == 0.0) then
			zbrent=b
			RETURN
		end if
		if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
			s=fb/fa
			if (a == c) then
				p=2.0_DP*xm*s
				q=1.0_DP-s
			else
				q=fa/fc
				r=fb/fc
				p=s*(2.0_DP*xm*q*(q-r)-(b-a)*(r-1.0_DP))
				q=(q-1.0_DP)*(r-1.0_DP)*(s-1.0_DP)
			end if
			if (p > 0.0) q=-q
			p=abs(p)
			if (2.0_DP*p  <  min(3.0_DP*xm*q-abs(tol1*q),abs(e*q))) then
				e=d
				d=p/q
			else
				d=xm
				e=d
			end if
		else
			d=xm
			e=d
		end if
		a=b
		fa=fb
		b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
		fb=func(b)
	end do
	call nrerror('zbrent: exceeded maximum iterations')
	zbrent=b
	
END FUNCTION zbrent
