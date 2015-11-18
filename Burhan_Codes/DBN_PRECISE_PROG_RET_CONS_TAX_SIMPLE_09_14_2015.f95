MODULE parameters

    use nrtype
    use nrutil

!	REAL(DP), PARAMETER :: beta=0.9515_DP
!	REAL(DP), PARAMETER :: phi=1.25_DP ! parameter for leisure utility ! when phi=0, do not value leisure.
!    REAL(DP), PARAMETER  :: rho_z= 0.50651_DP, sigma_z_eps=0.624_DP	! parameters for z process
!	REAL(DP), PARAMETER  :: rho_lambda=0.5_DP, sigma_lambda_eps=0.34_DP

    INTEGER(I4B),  PARAMETER ::  KeepSSatBench=1, SolveInelasticLabor=0, update_period=5
    
    REAL(DP), PARAMETER  :: rho_lambda=0.5_DP 
	REAL(DP) :: beta, phi, sigma_lambda_eps,  rho_z, sigma_z_eps    
	REAL(DP), PARAMETER  :: rho_e=0.9_DP, sigma_e_eps=0.20_DP    
	
	! production 
	REAL(DP), PARAMETER :: alpha=0.33_DP, Aprod=1.0_DP, DepRate=0.085_DP
	REAL(DP), PARAMETER :: mu=0.9_DP, hp=0.002_DP ! hp=home production (min cons)
		
	INTEGER(I4B), PARAMETER       :: MaxAge=81, RetAge=45 

	! asset grid
	REAL(DP), PARAMETER        ::  a_theta=4.0_DP , amax=100000.0_DP, amin=0.0001_DP, azero=0.0001_DP
	INTEGER(I4B),  PARAMETER :: na=201, fine_na=801, MaxSimuTime=500, max_age_category=7 
		
	! labor efficiency shocks
	! log(y)=  lambda + kappa + e 
	! lambda: inidividual fixed effect (fixed within generation)
	! kappa: life-cycle component
	! grid size
	INTEGER(I4B),  PARAMETER  :: ne=5       ! labor idiosyncratic shocks
	INTEGER(I4B),  PARAMETER  :: nlambda=5  ! labor efficiency fixed effect
	INTEGER(I4B),  PARAMETER  :: nz=7       ! inv ability fixed component
	REAL(DP),      PARAMETER    :: tauchen_m=3.0_DP ! how many std away from mean in tauchen
	REAL(DP), PARAMETER  :: mtauchen=3.0_DP,  brent_tol=0.00000001_DP
	
	REAL(DP), PARAMETER  :: tauWmin=0.01_DP, tauWinc=0.005_DP, tauC=0.075_DP

! progressive labor income tax parameters
	REAL(DP), PARAMETER  :: tauPL=0.185_DP, psi_PL=0.77_DP  ! 1-psi controls the level of tax, and tauPL controls progressivity
!	REAL(DP), PARAMETER  :: tauPL=0.0_DP, psi_PL=0.7_DP  ! this is flat tax that corresponds to our earlier case
	
END MODULE parameters


MODULE global
    USE parameters

    real(DP) , dimension(5) :: params

    REAL(DP), DIMENSION(MaxAge) :: pop, survP          ! population size by age and survival prob by age
	REAL(DP), DIMENSION(RetAge)                     :: kappagrid

	REAL(DP), DIMENSION(ne)                            :: egrid, Ge	, cdf_Ge
 	REAL(DP), DIMENSION(ne,ne)                       :: pr_e, cdf_pr_e
    REAL(DP), DIMENSION(MaxAge,ne)              :: Ge_byage, cdf_Ge_byage
                
	REAL(DP), DIMENSION(nlambda)                  :: lambdagrid,  Glambda, cdf_Glambda
	REAL(DP), DIMENSION(nlambda,nlambda)    :: pr_lambda, cdf_pr_lambda
	
    REAL(DP), DIMENSION(nlambda, ne)                  :: phi_lambda_e, RetY_lambda_e  ! Retirement income, the first 45 years are assigned zero. 

	! investment abilities
	REAL(DP), DIMENSION(nz)                  :: zgrid , Gz, cdf_Gz
	REAL(DP), DIMENSION(nz,nz)             :: pr_z, cdf_pr_z

                                                                                   ! I keep this format because I do not want to mess the ages
	REAL(DP), DIMENSION(MaxAge  , nlambda, ne) :: eff_un,  yh       ! Labor efficiency units, Labor efficiency units x wage

    REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) :: Cons, Hours, Aprime
    REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) ::ValueFunction
  
! Analytical solution for mu=1 for all the lifecycle not just retirement period
    REAL(DP), DIMENSION(MaxAge, na, nz, nlambda) :: AnRetCons,   AnRetValue, AnRetHours 
 
    REAL(DP) :: QBAR_bench,  QBAR_exp, NBAR_exp, NBAR_bench, Ebar_bench, Ebar_exp, wage_bench, wage_exp, rr_bench, rr_exp


! ASSET AND RESOURCE GRIDS
! NOTE THAT THESE GRIDS WILL CHANGE AS WE CHANGE TAX RATES AND FOR EACH INTERMEDIATE GOOD PRICE
! FOR EACH TAX RATE WE WILL HAVE DIFFERENT Y GRIDS
	REAL(DP), DIMENSION(na)                            :: agrid
	REAL(DP), DIMENSION(na,nz)                       :: YK, YW, MBK, MBW ! Y grids for CAPITAL INCOME tax and WEALTH taxes. 
    REAL(DP), DIMENSION(na,nz)                       :: YGRID, MBGRID
	
	REAL(DP), DIMENSION(fine_na)                     :: fine_agrid
    REAL(DP) :: tauk_bench, tauw_bench, tauL_bench
    REAL(DP) :: tauk_exp,     tauw_exp,    tauL_exp
    REAL(DP) :: tauK, tauW, tauL , rr, Ebar , wage, tauWindx, tauW_low, tauW_up

    INTEGER :: age, lambdai, zi, ai, ei    

    REAL(DP) :: NBAR, QBAR
    
    REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) ::DBN1, DBN_bench
    REAL(DP), DIMENSION(na) :: pr_a_dbn, cdf_a_dbn, tot_a_by_grid, cdf_tot_a_by_grid
    REAL(DP), DIMENSION(100) ::  cdf_tot_a_by_prctile
    INTEGER, DIMENSION(100) ::  prctile_ai
    REAL(DP) :: pop_25_60 , tothours_25_60, pop_pos_earn_25_60, tot_log_earnings_25_60, mean_log_earnings_25_60 
    REAL(DP) :: meanhours_25_60, Var_Log_Earnings_25_60, Std_Log_Earnings_25_60, MeanWealth, Wealth_Output
    REAL(DP) :: prct1_wealth, prct10_wealth, SSE_Moments, Min_SSE_Moments
    REAL(DP) :: YBAR, GBAR, GBAR_exp,GBAR_exp_old, GBAR_bench, Y_bench, Y_exp
    REAL(DP) :: consin, ain, psi
    REAL(DP) , DIMENSION(5) ::  Min_Moments
    INTEGER:: solving_bench, newiseed
   
END MODULE global

!====================================================================
!====================================================================

Module programfunctions
use parameters
use global
    
contains

FUNCTION FOC_R(aprimet)
IMPLICIT NONE   
real(DP):: aprimet,  MBaprime, FOC_R, yprime, cprime


MBaprime = ( 1.0_DP + ( rr *mu* (zgrid(zi)**mu) * (aprimet**(mu-1.0_DP)) - DepRate ) *(1.0_DP-tauK) )*(1.0_DP-tauW)

!print*,'mu=', mu,'tauK=',tauK
!print*,'(aprimet**(mu-1.0_DP))', (aprimet**(mu-1.0_DP))

!print*,'rr',rr,'zi',zi,'zgrid(zi)',zgrid(zi)
 

yprime =  ( aprimet+ ( rr * (zgrid(zi)* aprimet )**mu - DepRate* aprimet ) *(1.0_DP-tauK) )*(1.0_DP-tauW)
 

!cprime =   Linear_Int(End_Ret_Ygrid(age+1,:,zi,lambdai), &
!                                    & End_RetCons(age+1,:,zi,lambdai),na, yprime)    
 
! print*,Ygrid(:,zi)
! print*,RetCons(age+1,:,zi,lambdai)
 
cprime =   Linear_Int(Ygrid(:,zi), &
                                    & Cons(age+1,:,zi,lambdai,ei),na, yprime)    


FOC_R   = (1.0_DP / (YGRID(ai,zi)  + RetY_lambda_e(lambdai,ei) - aprimet )  &
           & - beta *  survP(age) *  MBaprime /cprime ) **2.0_DP
 

!print*,'aprimet=',aprimet,'MBaprime=',MBaprime,'cprime',cprime,'FOC_R=',FOC_R
!print*
!print*,'yprime=',yprime
!pause
END  FUNCTION


!====================================================================

FUNCTION FOC_WH(aprimet)
IMPLICIT NONE   
real(DP):: aprimet, ntemp, MBaprime, FOC_WH,  yprime, exp1overcprime
integer:: epindx
real(DP), dimension(ne):: cprime
REAL(DP):: brentvaluet

ain=aprimet

IF  (SolveInelasticLabor .eq. 1) THEN
      ntemp = 0.4_DP
      ELSE
         brentvaluet = brent(0.000001_DP, 0.4_DP, 0.99_DP, FOC_HA, brent_tol, ntemp)           
ENDIF

MBaprime = ( 1.0_DP + ( rr *mu* (zgrid(zi)**mu) * (aprimet**(mu-1.0_DP)) - DepRate ) *(1.0_DP-tauK) )*(1.0_DP-tauW)

yprime =  ( aprimet+ ( rr * (zgrid(zi)* aprimet )**mu - DepRate* aprimet ) *(1.0_DP-tauK) )*(1.0_DP-tauW)

! I have to evaluate the FOC in expectation over eindx prime given eindx

DO epindx=1,ne

      cprime(epindx)  =   Linear_Int(Ygrid(:,zi),&
                                    & Cons(age+1,:,zi,lambdai,epindx), na,    yprime  )
ENDDO

exp1overcprime = SUM( pr_e(ei,:) / cprime )

FOC_WH   = (1.0_DP / (YGRID(ai,zi)  + psi* (yh(age, lambdai,ei) * ntemp)**(1.0_DP-tauPL ) - aprimet )  &
           & - beta *  survP(age) *MBaprime* exp1overcprime) **2.0_DP 

END  FUNCTION


!====================================================================

FUNCTION FOC_H(hoursin)
IMPLICIT NONE   
real(DP):: hoursin, FOC_H 

FOC_H  = ( (psi*(1.0_DP-tauPL)*yh(age, lambdai,ei)**(1.0_DP-tauPL) )* (1.0_DP-hoursin)*(hoursin**(-tauPL)) - phi*consin )**2.0_DP

END  FUNCTION

!====================================================================

FUNCTION FOC_HA(hoursin)
IMPLICIT NONE   
real(DP):: hoursin, FOC_HA 

FOC_HA  = ( (psi*(1.0_DP-tauPL)*yh(age, lambdai,ei)**(1.0_DP-tauPL) )* (1.0_DP-hoursin)*(hoursin**(-tauPL)) &
        &  - phi*  ( YGRID(ai,zi) + psi*(yh(age, lambdai,ei)*hoursin)**(1.0_DP-tauPL) - ain  )   )**2.0_DP

END  FUNCTION

!================================================================================



FUNCTION Linear_Int(xa,ya,n,x)  
! Given arrays xa(1:n) and ya(1:n) of length n; this subroutine returns  linear interpolated value "y" at point "x"

USE parameters
IMPLICIT NONE      
      INTEGER:: n  
      REAL(DP)   :: x, yprime, xa(n),ya(n)  
      INTEGER k,khi,klo  
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
      if (h.eq.0.) pause 'bad xa input in linear int'  
      a=(xa(khi)-x)/h  
      b=(x-xa(klo))/h  
      Linear_Int=a*ya(klo)+b*ya(khi)     
     return  

END  Function Linear_Int



!================================================================================


FUNCTION Linear_Int_Aprime(xa,ya,n,x)  
! Given arrays xa(1:n) and ya(1:n) of length n; this subroutine returns  linear interpolated value "y" at point "x"

USE parameters
IMPLICIT NONE      
      INTEGER:: n  
      REAL(DP)   :: x, yprime, xa(n),ya(n)  
      INTEGER k,khi,klo  
      REAL(DP)   :: a,b,h, Linear_Int_Aprime
      

!if (x .gt. amax) then
!    klo = na-1
!    elseif (x .lt. amin) then
!        klo = 1
!        else
!            klo = ((x - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
!endif      
klo=1  
khi = klo + 1   
 
      h=xa(khi)-xa(klo)  
      if (h.eq.0.) pause 'bad xa input in linear int'  
      a=(xa(khi)-x)/h  
      b=(x-xa(klo))/h  
      Linear_Int_Aprime=a*ya(klo)+b*ya(khi)     
     return  

END  Function Linear_Int_Aprime



!================================================================================


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
!BL
	SUBROUTINE shft(a,b,c,d)
	REAL(DP), INTENT(OUT) :: a
	REAL(DP), INTENT(INOUT) :: b,c
	REAL(DP), INTENT(IN) :: d
	a=b
	b=c
	c=d
	END SUBROUTINE shft
	END FUNCTION brent

!====================================================================

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


end module programfunctions

!====================================================================


PROGRAM main
USE parameters
USE GLOBAL
use  programfunctions
IMPLICIT NONE
real :: start_time, finish_time, betaL, betaH, phiL, phiH, sigmalambL,sigmalambH, sigmazL, sigmazH, rhozL, rhozH
integer :: parindx1,  parindx2, parindx3, parindx4, parindx5, nbeta, nphi, nsigmalambda, nrhoz, nsigmaz

PRINT*,' **** PSI CORRECTED BY EBAR ----------'


! PARAMETERS FOR ELASTIC LABOR
params=[0.9455,  0.8,  0.255,  0.34,    1.00]
params=[0.9455,  0.6,  0.381,  0.335,  1.00]
!params=[0.947  ,  0.4,  0.490,  0.340,  1.01]
!params=[0.948,    0.2,  0.55,    0.34,    1.02]

! PARAMETERS FOR ELASTIC LABOR: NEW CALIBRATION WITH HIGHER DISPERSION
params=[0.9455,  0.8,  0.26058,  0.34,  1.00]
params=[0.9473463893,  0.6,   0.39882,  0.335,  1.00]
params=[0.9489479959,  0.40,  0.514595031738281,  0.34,  1.01]
!params=[0.949949,  0.20,  0.59388,  0.34,  1.02]


! PARAMETERS FOR INELASTIC LABOR
!params=[0.946,   0.80,   0.2325,   0.4484,  0.00]
!params=[0.9485, 0.40,   0.45,   0.4484,  0.0 ]
!params=[ 0.95,    0.20,   0.53,   0.4484,  0.00]
!params=[0.9471, 0.60,   0.35,    0.4484,   0.00]

PRINT*,''
If (DepRate .gt. 0.0_DP) then
    print*,'Positive Depretiation: Dep Rate = ', DepRate
    ELSE
        print*,'Zero Depretiation'
ENDIF
PRINT*,''

IF (SolveInelasticLabor .eq. 1) THEN
    params(5) = 0.0_DP
    print*,'Solving inelastic labor'
    else
    print*,'Solving elastic labor'
ENDIF

if (KeepSSatBench .eq. 1) then
    print*,''
    print*,'SS kept at benchmark level'
    print*,''
else
    print*,''
    print*,'SS adjusted with tax change'
    print*,''
endif    

print*,'na=',na,'update_period=',update_period

beta=params(1)
rho_z=params(2)
sigma_z_eps =params(3)
sigma_lambda_eps = params(4)
phi   =  params(5)
call cpu_time(start_time) 
tauK = 0.25_DP
tauL = 0.30_DP
tauW= 0.00_DP

!print*,  beta, rho_z, sigma_z_eps, sigma_lambda_eps,  phi

! ------- DO NOT REMOVE THE LINES BELOW

rr=  4.906133597851297E-002 
wage=  1.97429920063330 
Ebar=  1.82928004963637
Ebar_bench = Ebar

! ------- DO NOT REMOVE THE LINES ABOVES

!---------------------calibration begins-----------------------------------------------------------------------------------------------------

!
!!call cpu_time(start_time) 

! INELASTIC CASE PARAMETERS

!sigmalambL = 0.4484
!sigmalambH = 0.4484
!phiL  = 0.0_DP
!phiH = 0.0_DP
!
!
!betaL =0.9506
!betaH=0.9508
!rhozL = 0.2
!rhozH = 0.2
!sigmazL = 0.55
!sigmazH = 0.56
!
!betaL =0.949
!betaH=0.9495
!rhozL = 0.4
!rhozH = 0.4
!sigmazL = 0.47
!sigmazH = 0.475
!
!betaL =0.948
!betaH=0.949
!rhozL = 0.6
!rhozH = 0.6
!sigmazL = 0.37
!sigmazH = 0.375

!betaL =0.946
!betaH=0.948
!rhozL = 0.8
!rhozH = 0.8
!sigmazL = 0.245
!sigmazH = 0.255


params=[0.9455,  0.8,  0.26058,  0.34,  1.00]
params=[1.01606,  0.6,   0.46,  0.335,  1.1]
params=[1.01606,  0.40,  0.56,  0.34,  1.1]
!params=[0.949949,  0.20,  0.59388,  0.34,  1.02]

rhozL =PARAMS(2)
rhozH =PARAMS(2)
sigmalambL = PARAMS(4)
sigmalambH = PARAMS(4)

phiL  = PARAMS(5)  
phiH = PARAMS(5)  

sigmazL = PARAMS(3)  *1.00
sigmazH = PARAMS(3) *1.05

betaL   = PARAMS(1)  *1.00
betaH  = PARAMS(1)  *1.006


nbeta =3
nrhoz=1
nsigmaz=3 
nsigmalambda=1
nphi=1

Min_SSE_Moments=1000.0_DP

DO parindx3=1,nsigmalambda
DO parindx2=1,nphi
DO parindx1=1,nbeta
DO parindx4=1,nrhoz
DO parindx5=1,nsigmaz

    beta = betaL + real(parindx1-1,8) *(betaH-betaL)/max(real(nbeta-1,8),1.0_DP)
    phi   = phiL   + real(parindx2-1,8) *(phiH-phiL)/max(real(nphi-1,8),1.0_DP)
    sigma_lambda_eps = sigmalambL + real(parindx3-1,8)*(sigmalambH -sigmalambL) / max(real(nsigmalambda-1,8),1.0_DP)
    rho_z= rhozL   +  real(parindx4-1,8)*(rhozH-rhozL) / max(real(nrhoz-1,8),1.0_DP)
    sigma_z_eps = sigmazL +  real(parindx5-1,8)*(sigmazH-sigmazL) / max(real(nsigmaz-1,8),1.0_DP)

    print*,'parameters=',beta, rho_z,sigma_z_eps,sigma_lambda_eps, phi 

    solving_bench=1
    CALL  INITIALIZE
    CALL FIND_DBN_EQ
    print*,'parameters=',beta, rho_z,sigma_z_eps,sigma_lambda_eps, phi    
    CALL COMPUTE_STATS
 
   IF (SSE_Moments .lt. Min_SSE_Moments ) THEN
        Min_SSE_Moments =SSE_Moments
        params= [ beta, rho_z, sigma_z_eps, sigma_lambda_eps, phi ]
        Min_Moments = [  Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60  ]
    ENDIF
    !CALL WRITE_TO_FILE
ENDDO
ENDDO
ENDDO
ENDDO
ENDDO
print*, params
print*, Min_Moments
!!
!----------------------- calibration ends --------------------------------------

Print*,'--------------- SOLVING EXPERIMENT WITH BEST PARAMETERS -----------------'
PRINT*,''
print*,'CAPITAL TAX ECONOMY'
PRINT*,''

beta=params(1)
rho_z=params(2)
sigma_z_eps =params(3)
sigma_lambda_eps = params(4)
phi   =  params(5)

OPEN   (UNIT=3, FILE='params', STATUS='replace')
WRITE(unit=3, FMT=*) params
CLOSE (unit=3)

print*,'parameters=',beta, rho_z,sigma_z_eps,sigma_lambda_eps, phi

rr=  4.906133597851297E-002 
wage=  1.97429920063330 
Ebar=  1.82928004963637
Ebar_bench = Ebar

solving_bench=1
CALL  INITIALIZE
CALL FIND_DBN_EQ
CALL COMPUTE_STATS
CALL GOVNT_BUDGET

GBAR_bench=GBAR
QBAR_bench = QBAR 
NBAR_bench = NBAR 
Ebar_bench  = EBAR
rr_bench      = rr
wage_bench = wage
Y_bench       = YBAR
tauK_bench = tauK
tauw_bench = tauW
tauL_bench  = tauL
DBN_bench = DBN1
CALL WRITE_VARIABLES(1)


 !the following solves equilibrium of WEALTH tax economy
print*,'WEALTH TAX ECONOMY'
solving_bench=0 

tauK = 0.0_DP
GBAR_exp=0.0_DP
tauW=tauWmin
tauWindx=0.0_DP

DO WHILE (GBAR_exp .lt. GBAR_bench)
      GBAR_exp_old=GBAR_exp
      tauW = tauWmin +    tauWindx * tauWinc       
      CALL FIND_DBN_EQ
      CALL COMPUTE_STATS
      CALL GOVNT_BUDGET
      GBAR_exp = GBAR   
      tauWindx = tauWindx + 1.0_DP   
      print*,'tauW=', tauW, 'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
ENDDO

tauW_up  = tauW 
tauW_low = tauW  -  tauWinc  
tauW = tauW_low + tauWinc * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
print*,''
print*,'tauW_low =',tauW_low , 'tauW_up=', tauW_up, 'tauW=', tauW
print*,''
CALL FIND_DBN_EQ
CALL GOVNT_BUDGET

GBAR_exp = GBAR
print*,'tauW_low =', tauW_low, 'tauW_up=', tauW_up, 'tauW=', tauW, 'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.01 ) ! as long as the difference is greater than 0.05% continue
    if (GBAR_exp .gt. GBAR_bench ) then
        tauW_up  = tauW 
    else
        tauW_low = tauW 
    endif
    tauW = (tauW_low + tauW_up)/2.0_DP
    CALL FIND_DBN_EQ
    CALL GOVNT_BUDGET
    GBAR_exp = GBAR
    print*,'tauW_low =', tauW_low, 'tauW_up=', tauW_up, 'tauW=', tauW, 'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
ENDDO

CALL COMPUTE_STATS

GBAR_exp=GBAR
QBAR_exp = QBAR 
NBAR_exp = NBAR  
Y_exp = YBAR
Ebar_exp  = EBAR
rr_exp      = rr
wage_exp = wage
tauW_exp = tauW
tauK_exp = tauK
tauL_exp = tauL

CALL WRITE_VARIABLES(0)

CALL COMPUTE_WELFARE_GAIN
print*,'---------------------------'
print*,''
print*,'Output Gain Prct=', 100.0_DP*(Y_exp/Y_bench-1.0) 
print*,''
print*,'---------------------------'

call cpu_time(finish_time)
print*,'Total time =',finish_time-start_time

END PROGRAM

!====================================================================

SUBROUTINE COMPUTE_WELFARE_GAIN
USE GLOBAL 
use  programfunctions
IMPLICIT NONE
real(DP), DIMENSION(MaxAge):: CumDiscountF
REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) ::  ValueFunction_Bench, ValueFunction_Exp, Cons_Eq_Welfare
REAL(DP), dimension(nz) ::  temp_ce_by_z, temp_cons_by_z, temp_leisure_by_z, temp_dbn_by_z 
REAL(DP) :: frac_pos_welfare 
REAL(DP), dimension(MaxAge, nz) :: frac_pos_welfare_by_age_z, size_pos_welfare_by_age_z, size_by_age_z_bench, size_by_age_z_exp
INTEGER, dimension(max_age_category+1) :: age_limit
INTEGER :: age_group_counter
REAL(DP), dimension(max_age_category,nz) :: CE_by_agegroup_z 
REAL(DP), dimension(max_age_category,nz) :: size_pos_welfare_by_agegroup_z, frac_pos_welfare_by_agegroup_z  
REAL(DP), dimension(max_age_category,nz) ::  tot_wealth_by_agegroup_z_bench, size_by_agegroup_z_bench, size_by_agegroup_z_exp
REAL(DP), dimension(max_age_category,nz) ::  tot_wealth_by_agegroup_z_exp

age_limit = [0, 5, 15, 25, 35, 45, 55, MaxAge ]

CumDiscountF(MaxAge)=1.0_DP
DO age=MaxAge-1,1,-1
    CumDiscountF(age)   = 1.0_DP + beta * survP(age) *CumDiscountF(age+1) 
ENDDO
!print*,CumDiscountF
!PAUSE

solving_bench = 1
tauK=tauK_bench
tauL=tauL_bench
tauW=tauW_bench
rr = rr_bench
wage = wage_bench
Ebar = Ebar_bench
print*,'BENCH: rr=',rr,'wage=',wage,'Ebar=',Ebar
CALL FORM_Y_MB_GRID(YGRID,MBGRID) 
CALL ComputeLaborUnits(Ebar, wage) 
CALL EGM_RETIREMENT_WORKING_PERIOD 
!CALL COMPUTE_VALUE_FUNCTION_LINEAR
CALL COMPUTE_VALUE_FUNCTION_SPLINE  !-----------------.................................................
!CALL COMPUTE_STATS

ValueFunction_Bench = ValueFunction

OPEN (UNIT=7, FILE='cons_by_age_z_bench', STATUS='replace')    
OPEN (UNIT=8, FILE='leisure_by_age_z_bench', STATUS='replace')    
OPEN (UNIT=9, FILE='dbn_by_age_z_bench', STATUS='replace')   
DO age=1,MaxAge    
    DO zi=1,nz
          temp_cons_by_z(zi)       = sum(Cons(age,:,zi,:,:)*DBN_bench(age,:,zi,:,:))/sum(DBN_bench(age,:,zi,:,:))
          temp_leisure_by_z(zi)    = sum((1.0_DP-HOURS(age,:,zi,:,:))*DBN_bench(age,:,zi,:,:))/sum(DBN_bench(age,:,zi,:,:))
          size_by_age_z_bench(age, zi)        = sum(DBN_bench(age,:,zi,:,:))
    ENDDO ! zi
    WRITE  (UNIT=7, FMT=*) temp_cons_by_z
    WRITE  (UNIT=8, FMT=*) temp_leisure_by_z
    WRITE  (UNIT=9, FMT=*)  size_by_age_z_bench(age, :) 
ENDDO
close (unit=7)
close (unit=8)
close (unit=9)
  
! This prints Aprime for different "z" for median lambda and e
OPEN (UNIT=5, FILE='aprime_age1_bench', STATUS='replace')     
DO zi=1,nz 
    WRITE  (UNIT=5, FMT=*) Aprime(1, :, zi, nlambda/2+1, ne/2+1)
ENDDO
close (unit=5)

OPEN (UNIT=5, FILE='aprime_age16_bench', STATUS='replace')   
DO zi=1,nz 
    WRITE  (UNIT=5, FMT=*) Aprime(16, :, zi, nlambda/2+1, ne/2+1)
ENDDO
close (unit=5)  

OPEN (UNIT=5, FILE='aprime_age31_bench', STATUS='replace')   
DO zi=1,nz 
    WRITE  (UNIT=5, FMT=*) Aprime(31, :, zi, nlambda/2+1, ne/2+1)
ENDDO
close (unit=5)  

OPEN (UNIT=5, FILE='aprime_age46_bench', STATUS='replace')   
DO zi=1,nz 
    WRITE  (UNIT=5, FMT=*) Aprime(46, :, zi, nlambda/2+1, ne/2+1)
ENDDO
close (unit=5)  



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

OPEN (UNIT=6, FILE='mean_wealth_by_agegroup_z_bench', STATUS='replace')  
OPEN (UNIT=7, FILE='size_by_agegroup_z_bench', STATUS='replace')  
DO age_group_counter=1,max_age_category
    WRITE  (UNIT=6, FMT=*)   tot_wealth_by_agegroup_z_bench(age_group_counter,:)/ size_by_agegroup_z_bench(age_group_counter,:)
    WRITE  (UNIT=7, FMT=*)    size_by_agegroup_z_bench(age_group_counter,:)
ENDDO
close (UNIT=6)
close (UNIT=7)


!=========================================== SOLVING EXP NEXT =================================

solving_bench = 0  
tauK=tauK_exp
tauL=tauL_exp
tauW=tauW_exp
rr = rr_exp
wage = wage_exp
Ebar = Ebar_exp
print*,' EXP: rr=',rr,'wage=',wage,'Ebar=',Ebar
CALL FORM_Y_MB_GRID(YGRID,MBGRID) 
CALL ComputeLaborUnits(Ebar, wage) 
CALL EGM_RETIREMENT_WORKING_PERIOD 
!CALL COMPUTE_VALUE_FUNCTION_LINEAR
CALL COMPUTE_VALUE_FUNCTION_SPLINE  !-----------------.................................................
!CALL COMPUTE_STATS

ValueFunction_Exp = ValueFunction

OPEN (UNIT=7, FILE='cons_by_age_z_exp', STATUS='replace')    
OPEN (UNIT=8, FILE='leisure_by_age_z_exp', STATUS='replace')   
OPEN (UNIT=9, FILE='dbn_by_age_z_exp', STATUS='replace')   
DO age=1,MaxAge 
    DO zi=1,nz
          temp_cons_by_z(zi)       = sum(Cons(age,:,zi,:,:)*DBN1(age,:,zi,:,:))/sum(DBN1(age,:,zi,:,:))
          temp_leisure_by_z(zi)    = sum((1.0_DP-HOURS(age,:,zi,:,:))*DBN1(age,:,zi,:,:))/sum(DBN1(age,:,zi,:,:))
          size_by_age_z_exp(age, zi)         = sum(DBN1(age,:,zi,:,:))
    ENDDO ! zi
    WRITE  (UNIT=7, FMT=*) temp_cons_by_z
    WRITE  (UNIT=8, FMT=*) temp_leisure_by_z
    WRITE  (UNIT=9, FMT=*) size_by_age_z_exp(age, :)  
ENDDO
close (unit=7)
close (unit=8)
close (unit=9)

! This prints Aprime for different "z" for median lambda and e
OPEN (UNIT=5, FILE='aprime_age1_exp', STATUS='replace')     
DO zi=1,nz 
    WRITE  (UNIT=5, FMT=*) Aprime(1, :, zi, nlambda/2+1, ne/2+1)
ENDDO
close (unit=5)

OPEN (UNIT=5, FILE='aprime_age16_exp', STATUS='replace')   
DO zi=1,nz 
    WRITE  (UNIT=5, FMT=*) Aprime(16, :, zi, nlambda/2+1, ne/2+1)
ENDDO
close (unit=5)  

OPEN (UNIT=5, FILE='aprime_age31_exp', STATUS='replace')   
DO zi=1,nz 
    WRITE  (UNIT=5, FMT=*) Aprime(31, :, zi, nlambda/2+1, ne/2+1)
ENDDO
close (unit=5)  

OPEN (UNIT=5, FILE='aprime_age46_exp', STATUS='replace')   
DO zi=1,nz 
    WRITE  (UNIT=5, FMT=*) Aprime(46, :, zi, nlambda/2+1, ne/2+1)
ENDDO
close (unit=5)  


OPEN (UNIT=5, FILE='CE_NEWBORN', STATUS='replace')  
OPEN (UNIT=6, FILE='CE', STATUS='replace')  
OPEN (UNIT=7, FILE='CE_by_age', STATUS='replace')  
OPEN (UNIT=8, FILE='CE_by_age_z', STATUS='replace')  

DO age=1,MaxAge
    Cons_Eq_Welfare(age,:,:,:,:)=exp((ValueFunction_exp(age,:,:,:,:)-ValueFunction_Bench(age,:,:,:,:))/CumDiscountF(age))-1.0_DP
    WRITE  (UNIT=7, FMT=*) 100*sum(Cons_Eq_Welfare(age,:,:,:,:)*DBN_bench(age,:,:,:,:))/sum(DBN_bench(age,:,:,:,:))
    DO zi=1,nz
         temp_ce_by_z(zi) = 100*sum(Cons_Eq_Welfare(age,:,zi,:,:)*DBN_bench(age,:,zi,:,:))/sum(DBN_bench(age,:,zi,:,:))
    ENDDO
    WRITE  (UNIT=8, FMT=*) temp_ce_by_z
    print*,'age=',age, temp_ce_by_z, ', mean:  ', &
        & 100*sum(Cons_Eq_Welfare(age,:,:,:,:)*DBN_bench(age,:,:,:,:))/sum(DBN_bench(age,:,:,:,:))
ENDDO

WRITE  (UNIT=5, FMT=*) 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN_bench(1,:,:,:,:))/sum(DBN_bench(1,:,:,:,:))
WRITE  (UNIT=5, FMT=*) 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN1(1,:,:,:,:))/sum(DBN1(1,:,:,:,:))
WRITE  (UNIT=6, FMT=*) 100.0_DP*sum(Cons_Eq_Welfare*DBN_bench)
WRITE  (UNIT=6, FMT=*) 100.0_DP*sum(Cons_Eq_Welfare*DBN1)

close (unit=5)
close (unit=6)
close (unit=7)
close (unit=8)

! CE by AGE-Z GROUP
OPEN (UNIT=8, FILE='CE_by_AgeGroup_z', STATUS='replace') 
DO zi=1,nz
    DO age_group_counter=1,max_age_category
         CE_by_agegroup_z(age_group_counter,zi)= &
            & 100*sum(Cons_Eq_Welfare(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:)* &
            &                         DBN_bench(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:))/&
            &                sum( DBN_bench(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:))
    ENDDO
ENDDO
DO age_group_counter=1,max_age_category
    WRITE  (UNIT=8, FMT=*)  CE_by_agegroup_z(age_group_counter,:)
ENDDO
close (unit=8)

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


OPEN (UNIT=6, FILE='mean_wealth_by_agegroup_z_exp', STATUS='replace')  
OPEN (UNIT=7, FILE='size_by_agegroup_z_exp', STATUS='replace')  
DO age_group_counter=1,max_age_category
    WRITE  (UNIT=6, FMT=*)   tot_wealth_by_agegroup_z_exp(age_group_counter,:)/ size_by_agegroup_z_exp(age_group_counter,:)
    WRITE  (UNIT=7, FMT=*)   size_by_agegroup_z_exp(age_group_counter,:)
ENDDO
close (UNIT=6)
close (UNIT=7)

OPEN (UNIT=6, FILE='frac_pos_welfare_by_agegroup_z', STATUS='replace')  
DO age_group_counter=1,max_age_category
    WRITE  (UNIT=6, FMT=*)  size_pos_welfare_by_agegroup_z(age_group_counter,:)/ size_by_agegroup_z_bench(age_group_counter,:)
ENDDO
close (UNIT=6)

OPEN (UNIT=6, FILE='frac_pos_welfare', STATUS='replace')  
WRITE  (UNIT=6, FMT=*) frac_pos_welfare
close (unit=6)

OPEN (UNIT=6, FILE='frac_pos_welfare_by_age_z', STATUS='replace')  
DO age=1, MaxAge
    WRITE  (UNIT=6, FMT=*) size_pos_welfare_by_age_z(age,:)/size_by_age_z_bench(age,:)
ENDDO
close (UNIT=6)


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

!====================================================================

SUBROUTINE COMPUTE_VALUE_FUNCTION_SPLINE
USE GLOBAL 
use  programfunctions
IMPLICIT NONE
INTEGER :: tklo, tkhi
REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) :: PrAprimelo, PrAprimehi
REAL(DP), DIMENSION(na) :: ValueP1, ValueP2, ValueP, ExpValueP

print*,'VALUE FUNCTION SPLINE'


age=MaxAge
DO ai=1,na    
    DO zi=1,nz
        DO lambdai=1,nlambda          
              DO ei=1,ne
                  ValueFunction(age, ai, zi, lambdai, ei) = log(Cons(age, ai, zi, lambdai, ei)) 
!                  print*,Cons(age, ai, zi, lambdai, ei),  ValueFunction(age, ai, zi, lambdai, ei) 
!                  pause
              ENDDO ! ei          
        ENDDO ! lambdai
    ENDDO ! zi
ENDDO ! ai

! Retirement Period
DO age=MaxAge-1,RetAge,-1
    DO zi=1,nz
        DO lambdai=1,nlambda          
              DO ei=1,ne            
          
                    CALL spline( agrid, ValueFunction(age+1, :, zi, lambdai, ei) , na , &
                    & 1.0_DP/Cons(age+1, 1, zi, lambdai,ei) , 1.0_DP/Cons(age+1, na, zi, lambdai,ei) , ValueP2)  
                  
                    DO ai=1,na    
                         call splint( agrid, ValueFunction(age+1, :, zi, lambdai, ei), &
                                & ValueP2, na, Aprime(age,ai,zi,lambdai, ei), ValueP(ai))   
    
                         ValueFunction(age, ai, zi, lambdai, ei) = log(Cons(age, ai, zi, lambdai, ei)) &
                              & + beta*survP(age)* ValueP(ai)
                    ENDDO ! ai
              
            ENDDO ! ei          
        ENDDO ! lambdai
    ENDDO ! zi
ENDDO ! age
!print*,ValueFunction


! Working Period
DO age=RetAge-1,1,-1
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

!====================================================================

SUBROUTINE COMPUTE_VALUE_FUNCTION_LINEAR
USE GLOBAL 
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
                  ValueFunction(age, ai, zi, lambdai, ei) = log(Cons(age, ai, zi, lambdai, ei)) 
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
             
                  ValueFunction(age, ai, zi, lambdai, ei) = log(Cons(age, ai, zi, lambdai, ei)) &
                      & + beta*survP(age)* (PrAprimelo(age,ai,zi,lambdai, ei)*ValueFunction(age+1, tklo, zi, lambdai, ei)&
                      & +                             PrAprimehi(age,ai,zi,lambdai, ei)*ValueFunction(age+1, tkhi, zi, lambdai, ei))
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

!====================================================================

SUBROUTINE GOVNT_BUDGET
USE PARAMETERS
USE GLOBAL
use  programfunctions
IMPLICIT NONE

real(DP) ::  GBAR_K,  GBAR_W,  GBAR_L, SSC_Payments

GBAR=0.0_DP
GBAR_K =0.0_DP
GBAR_W=0.0_DP
GBAR_L =0.0_DP

DO age=1, MaxAge
DO ai=1,na
DO zi=1,nz
DO lambdai=1,nlambda
DO ei=1,ne
    GBAR = GBAR + DBN1(age,ai,zi,lambdai,ei) * ( tauK*( rr*(agrid(ai)*zgrid(zi))**mu-DepRate*agrid(ai) )  &
          & + tauW * ( rr*(agrid(ai)*zgrid(zi))**mu-DepRate*agrid(ai) + agrid(ai) )  &
          & - tauK * tauW * ( rr*(agrid(ai)*zgrid(zi))**mu-DepRate*agrid(ai) ) &
          & + yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei)  &
          & -  psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei))**(1.0_DP-tauPL)  &
          & + tauC * cons(age, ai, zi, lambdai,ei)  )         

    GBAR_L = GBAR_L  + DBN1(age,ai,zi,lambdai,ei) * (  yh(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei) &
               &- psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei))**(1.0_DP-tauPL) )
    
ENDDO
ENDDO
ENDDO
ENDDO
ENDDO

SSC_Payments=0.0_DP
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

print*,'GBAR=',GBAR,'SSC_Payments=', SSC_Payments, 'GBAR_L=',GBAR_L,'Av. Labor Tax=',GBAR_L/Ebar 

END  SUBROUTINE GOVNT_BUDGET

!====================================================================

SUBROUTINE FIND_DBN_EQ
USE PARAMETERS
USE GLOBAL
use  programfunctions
IMPLICIT NONE
INTEGER:: tklo, tkhi, age1, age2, z1, z2, a1, a2, lambda1, lambda2, e1, e2, DBN_iter, simutime, iter_indx
REAL :: DBN_dist, DBN_criteria
REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) :: PrAprimelo, PrAprimehi, DBN2
INTEGER,  DIMENSION(MaxAge, na, nz, nlambda, ne) :: Aplo, Aphi

DBN_criteria = 1.0E-07_DP

! Form YGRID for the capital income economy given interest rate "rr"
CALL FORM_Y_MB_GRID(YGRID,MBGRID) 
CALL ComputeLaborUnits(Ebar, wage) 
CALL EGM_RETIREMENT_WORKING_PERIOD 
!
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
           &DBN2(1,Aplo(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1) + DBN1(age1, a1, z1, lambda1, e1) &
                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2)  *  PrAprimelo(age1, a1, z1, lambda1, e1)     
             DBN2(1,Aphi(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1)   =  &
           &DBN2(1,Aphi(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1) + DBN1(age1, a1, z1, lambda1, e1) & 
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
           &DBN2(1,Aplo(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1) + DBN1(age1, a1, z1, lambda1, e1) &
                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2)  *  PrAprimelo(age1, a1, z1, lambda1, e1)     
             DBN2(1,Aphi(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1)   =  &
           &DBN2(1,Aphi(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1) + DBN1(age1, a1, z1, lambda1, e1) & 
                & *(1.0_DP-survP(age1)) * pr_z(z1,z2) *pr_lambda(lambda1,lambda2)*PrAprimehi(age1,a1,z1,lambda1,e1)   
        ENDDO
        ENDDO
        
        ! Those who live stay at z1, lambda1, and also e1 since they are retired
        !e2=e1
            DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1), z1,lambda1,e1) =  &
          &DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1), z1,lambda1,e1) + DBN1(age1, a1, z1, lambda1, e1) &
                & * survP(age1) * PrAprimelo(age1, a1, z1, lambda1, e1)     
            DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1), z1,lambda1,e1) =  &
          &DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1), z1,lambda1,e1) + DBN1(age1, a1, z1, lambda1, e1) &
                & * survP(age1) * PrAprimehi(age1, a1, z1, lambda1, e1) 
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    
    DO age1=1,RetAge-2
    DO a1=1,na
    DO z1=1,nz
    DO lambda1=1,nlambda
    DO e1=1, ne
        ! Those who die, switch to z2, lambda2 and start at ne/2+1
        DO z2=1,nz
        DO lambda2=1,nlambda
             DBN2(1,Aplo(age1, a1, z1, lambda1, e1), z2,lambda2,ne/2+1)   =  &
           &DBN2(1,Aplo(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1) + DBN1(age1, a1, z1, lambda1, e1) &
                & * (1.0_DP-survP(age1)) * pr_z(z1,z2) * pr_lambda(lambda1,lambda2)  *  PrAprimelo(age1, a1, z1, lambda1, e1)     
             DBN2(1,Aphi(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1)   =  &
           &DBN2(1,Aphi(age1, a1, z1, lambda1, e1), z2, lambda2, ne/2+1) + DBN1(age1, a1, z1, lambda1, e1) & 
                & *(1.0_DP-survP(age1)) * pr_z(z1,z2) *pr_lambda(lambda1,lambda2)*PrAprimehi(age1,a1,z1,lambda1,e1)   
        ENDDO
        ENDDO
        ! Those who live stay at z1, lambda1, but switch to e2
        DO e2=1, ne
            DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1), z1,lambda1,e2) =  &
          &DBN2(age1+1, Aplo(age1, a1, z1, lambda1, e1), z1,lambda1,e2) + DBN1(age1, a1, z1, lambda1, e1) &
                & * survP(age1) * pr_e(e1,e2)   * PrAprimelo(age1, a1, z1, lambda1, e1)     
            DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1), z1,lambda1,e2) =  &
          &DBN2(age1+1, Aphi(age1, a1, z1, lambda1, e1), z1,lambda1,e2) + DBN1(age1, a1, z1, lambda1, e1) &
                & * survP(age1) * pr_e(e1,e2)  *  PrAprimehi(age1, a1, z1, lambda1, e1) 
        ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    
    
    DBN_dist = maxval(abs(DBN2-DBN1))
    DBN1 = DBN2
!    print*,'DBN_dist=',DBN_dist

!---------- update decision rules if iter_indx is greater than update period

    IF (iter_indx .ge. update_period) THEN
        QBAR =0.0
        NBAR =0.0
        DO age1=1,MaxAge
        DO z1=1,nz
        DO a1=1,na
        DO lambda1=1,nlambda
        DO e1=1, ne
             QBAR= QBAR+  DBN1(age1, a1, z1, lambda1, e1) * ( zgrid(z1) *agrid(a1) )**mu
             NBAR= NBAR+  DBN1(age1, a1, z1, lambda1, e1)  &
                             & * eff_un(age1, lambda1, e1) * Hours(age1, a1, z1, lambda1,e1)
        ENDDO
        ENDDO
        ENDDO    
        ENDDO    
        ENDDO    
    
        QBAR= ( QBAR)**(1.0_DP/mu)                
        rr = alpha* QBAR **(alpha-mu) * NBAR **(1.0_DP-alpha)
        YBAR  = QBAR ** alpha * NBAR **(1.0_DP-alpha)
        wage = (1.0_DP-alpha)*QBAR **alpha * NBAR  **(-alpha)
        Ebar  = wage  * NBAR  * sum(pop)/sum(pop(1:RetAge-1))
    !    print*,'DBN_dist=',DBN_dist, 'QBAR=', QBAR ,  'NBAR=', NBAR 
     
        CALL FORM_Y_MB_GRID(YGRID,MBGRID) 
        CALL ComputeLaborUnits(Ebar, wage) 
        CALL EGM_RETIREMENT_WORKING_PERIOD
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
        
        PrAprimelo = min (PrAprimelo, 1.0_DP)
        PrAprimelo = max(PrAprimelo, 0.0_DP)
        PrAprimehi = min (PrAprimehi, 1.0_DP)
        PrAprimehi = max(PrAprimehi, 0.0_DP)   
        iter_indx=0
    ENDIF
    
    iter_indx  = iter_indx + 1
    simutime = simutime +1
   
!call COMPUTE_STATS   
 
ENDDO ! WHILE
!
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

!===============================================================================

SUBROUTINE COMPUTE_STATS
USE GLOBAL
use  programfunctions
IMPLICIT NONE
INTEGER :: prctile
REAL(DP), DIMENSION(nz) :: cdf_Gz_DBN 
REAL(DP):: MeanReturn, StdReturn, VarReturn 
REAL(DP), DIMENSION(nz):: MeanReturn_by_z, size_by_z

DO zi=1,nz
    cdf_Gz_DBN(zi) = sum(DBN1(:,:,zi,:,:))
ENDDO
!print*,'cdf_Gz_DBN ='
!print*,cdf_Gz_DBN

DO ai=1,na
     pr_a_dbn(ai)   = sum(DBN1(:,ai,:,:,:)) 
     cdf_a_dbn(ai) = sum( pr_a_dbn(1:ai) )      
     tot_a_by_grid(ai) = sum(DBN1(:,ai,:,:,:) * agrid(ai) )
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
!    print*,prctile, REAL(prctile,8)/100.0_DP,  ai
    IF (ai .gt. 1) THEN
        cdf_tot_a_by_prctile(prctile)  =   cdf_tot_a_by_grid(ai-1) + (REAL(prctile,8)/100.0_DP - cdf_a_dbn(ai-1))*agrid(ai) 
        else
             cdf_tot_a_by_prctile(prctile)  = (REAL(prctile,8)/100.0_DP )*agrid(ai)     
    ENDIF
ENDDO
print*,''
prct1_wealth   =  1.0_DP-cdf_tot_a_by_prctile(99)/cdf_tot_a_by_prctile(100)
prct10_wealth =  1.0_DP-cdf_tot_a_by_prctile(90)/cdf_tot_a_by_prctile(100)

! COMPUTE AVERAGE HOURS FOR AGES 25-60 (5-40 IN THE MODEL) INCLUDING NON-WORKERS
! COMPUTE VARIANCE OF LOG EARNINGS FOR 25-60 FOR THOSE WHO WORK MORE THAN 260 HOURS
! WHICH CORRESPOND TO 0.055 IN THE MODEL
pop_25_60        = 0.0_DP
tothours_25_60 = 0.0_DP
pop_pos_earn_25_60  = 0.0_DP
tot_log_earnings_25_60 = 0.0_DP 
Var_Log_Earnings_25_60=0.0_DP
DO age=5,40
DO ai=1,na
DO zi=1,nz
DO lambdai=1,nlambda
DO ei=1,ne
      tothours_25_60 = tothours_25_60 + DBN1(age, ai, zi, lambdai, ei)  * Hours(age, ai, zi, lambdai,ei)
      pop_25_60  = pop_25_60 +  DBN1(age, ai, zi, lambdai, ei)
      IF (Hours(age, ai, zi, lambdai, ei) .ge. 0.055) THEN
            tot_log_earnings_25_60 =  tot_log_earnings_25_60 + DBN1(age, ai, zi, lambdai, ei)  &
                         & *  log( wage * yh(age, lambdai, ei) * Hours(age, ai, zi, lambdai, ei) )
             pop_pos_earn_25_60  = pop_pos_earn_25_60 +  DBN1(age, ai, zi, lambdai, ei)
      ENDIF
ENDDO
ENDDO
ENDDO
ENDDO
ENDDO
meanhours_25_60  = tothours_25_60 / pop_25_60
mean_log_earnings_25_60 = tot_log_earnings_25_60 / pop_pos_earn_25_60

!print*,'pop_pos_earn_25_60=',pop_pos_earn_25_60
!print*,'mean_log_earnings_25_60 =',mean_log_earnings_25_60 
!print*,' tot_log_earnings_25_60=', tot_log_earnings_25_60
!print*,'sum DBN=',sum(DBN1)
!print*,'yh=',yh
!print*,'wage=',wage

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
!print*,'------------------'
!!print*,'Var_Log_Earnings_25_60=',Var_Log_Earnings_25_60
!print*,'std of log earnings=',Std_Log_Earnings_25_60
!print*,'------------------'

MeanWealth =0.0_DP
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
     size_by_z(zi) = size_by_z(zi) + DBN1(age, ai, zi, lambdai, ei) 
ENDDO
ENDDO
ENDDO    
ENDDO    
ENDDO    
MeanReturn_by_z = MeanReturn_by_z / size_by_z

!print*, 'MeanReturn=',MeanReturn, 'StdReturn=', StdReturn
!print*,'MeanReturn_by_z=',MeanReturn_by_z

SSE_Moments = (Wealth_Output-3.0_DP)**2.0_DP + (prct1_wealth-0.35_DP)**2.0_DP  + (prct10_wealth-0.75_DP)**2.0_DP &
                   & + (Std_Log_Earnings_25_60 -0.8_DP)**2.0_DP + (meanhours_25_60-0.4_DP)**2.0_DP

print*,'MOMENTS=', Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60 
print*,'--------------------------------------------------'

END SUBROUTINE COMPUTE_STATS

!====================================================================

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
WRITE  (UNIT=2, FMT=*)  'tauW =',tauW
WRITE  (UNIT=2, FMT=*)  ''
WRITE  (UNIT=2, FMT=*)  'MOMENTS'
WRITE  (UNIT=2, FMT=*)  Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60

ENDIF

close (unit=2)

END SUBROUTINE WRITE_VARIABLES

!====================================================================

SUBROUTINE EGM_RETIREMENT_WORKING_PERIOD !(YGRID, MBGRID)
USE GLOBAL 
use  programfunctions
IMPLICIT NONE
!REAL(DP), DIMENSION(na,nz), INTENT(IN)  :: YGRID, MBGRID
!REAL(DP), DIMENSION(fine_na, nz) :: FineYGRID
REAL(DP), DIMENSION(MaxAge, na, nz, nlambda) :: CorrectRetValueP1, CorrectRetValueP2
REAL(DP) :: AnRetCCtemp, cprime, AntRetVVtemp, disty, slopeV, slopeC, linearC, linearV
REAL(DP) ::  ap_temp, brentvalue
INTEGER :: na1, na2, tempai
real(DP) :: tempvar1, tempvar2, tempvar3
REAL(DP), DIMENSION(na) :: EndoCons, EndoYgrid, EndoHours

! SET HOURS TO ZERO SO THAT RETIRED HAS ZERO HOURS
Hours(RetAge:MaxAge, :, :, :,: ) = 0.0_DP
Hours(1:RetAge-1       , :, :, :,: ) =0.4_DP

psi=psi_PL

age=MaxAge

DO lambdai=1,nlambda
    DO zi=1,nz
        DO ei=1,ne
                DO ai=1,na
                        Cons(age,ai,zi,lambdai,ei) =  YGRID(ai, zi) + RetY_lambda_e(lambdai,ei) 
                  ENDDO ! ai
        ENDDO ! ei
    ENDDO ! zi
ENDDO ! lambdai
Aprime(age, :, :, :, :) = 0.0_DP


!------RETIREMENT PERIOD---------------------------

DO age=MaxAge-1,RetAge,-1
!    print*,'age=',age
    DO lambdai=1,nlambda
        DO zi=1,nz
        DO ei=1,ne
              DO ai=1,na
                    EndoCons(ai) =   Cons(age+1, ai, zi, lambdai,ei)/( beta*survP(age)*MBGRID(ai,zi))    
                    EndoYgrid(ai) =  agrid(ai) +  EndoCons(ai)   - RetY_lambda_e(lambdai,ei)
              ENDDO ! ai

              tempai=1           
              DO WHILE ( YGRID(tempai,zi) .lt. EndoYgrid(1) )
                      tempai = tempai + 1
              ENDDO
                
              ! Find  decision rules on exogenous grids
              DO ai=tempai,na               
                    ! CONSUMPTION ON EXOGENOUS GRIDS
                    Cons(age, ai, zi, lambdai, ei) = Linear_Int(EndoYgrid, &
                                        & EndoCons,na, YGRID(ai,zi))                                                                                 
                    Aprime(age, ai, zi, lambdai,ei) = YGRID(ai,zi)+ RetY_lambda_e(lambdai,ei)  &
                                        & - Cons(age, ai, zi, lambdai, ei)                                         
                       If   (Aprime(age, ai, zi, lambdai,ei)  .lt. amin) then
                              Aprime(age, ai, zi, lambdai,ei)=amin
                              Cons(age, ai, zi, lambdai,ei) =  YGRID(ai,zi)  + RetY_lambda_e(lambdai,ei)  &
                                                                       & - Aprime(age, ai, zi, lambdai,ei) 
                              IF (Cons(age, ai, zi, lambdai, ei) .le. 0.0_DP)  THEN
                                    print*,'r1: Cons(age, ai, zi, lambdai,ei)=',Cons(age, ai, zi, lambdai, ei)
                              ENDIF                   
                        endif                                                         
              ENDDO ! ai  

               ai=1           
               DO WHILE ( YGRID(ai,zi) .lt. EndoYgrid(1) )
!                      ap_temp    = Aprime(age, tempai, zi, lambdai,1) 
                      brentvalue = brent(min(amin,YGRID(ai,zi)), (amin+YGRID(ai,zi))/2.0_DP , &
                                        & YGRID(ai,zi)+RetY_lambda_e(lambdai,ei) *0.95_DP, &
                                        & FOC_R, brent_tol, Aprime(age, ai, zi, lambdai,ei))

                      Cons(age, ai, zi, lambdai, ei) =  YGRID(ai,zi)  + RetY_lambda_e(lambdai,ei)  - Aprime(age, ai, zi, lambdai,ei)
                      IF (Cons(age, ai, zi, lambdai, ei) .le. 0.0_DP)  THEN
                            print*,'r2: Cons(age, ai, zi, lambdai,ei)=',Cons(age, ai, zi, lambdai,ei)
                            print*,'Aprime(age, ai, zi, lambdai,ei)=',Aprime(age, ai, zi, lambdai,ei)
                            print*,'YGRID(ai,zi)+RetY_lambda_e(lambdai,ei)=',YGRID(ai,zi)+RetY_lambda_e(lambdai,ei)
                            print*,'YGRID(ai,zi)=',YGRID(ai,zi),'EndoYgrid(1)=',EndoYgrid(1)
                            print*,'RetY_lambda_e(lambdai,ei)=',RetY_lambda_e(lambdai,ei)
                            print*,'lambdai=',lambdai
                       ENDIF                   
                      ai = ai + 1
               ENDDO  
               
        ENDDO ! ei     
        ENDDO ! zi
    ENDDO ! lambdai

    
ENDDO !age


!------RETIREMENT PERIOD ENDS---------------------------

!------WORKING PERIOD STARTS---------------------------

DO age=RetAge-1,1,-1

    DO lambdai=1,nlambda
        DO zi=1,nz
            DO ei=1,ne
                
                DO ai=1,na
                      EndoCons(ai) =   1.0_DP &
                                    & /( beta*survP(age)*MBGRID(ai,zi)*SUM(pr_e(ei,:) /Cons(age+1,ai,zi,lambdai,:)))    
                                    
                      consin =  EndoCons(ai)     
                      
                      EndoHours(ai) =0.4_DP                        
                      IF  (SolveInelasticLabor .ne. 1) THEN
                              brentvalue = brent(0.000001_DP, 0.4_DP, 0.99_DP,&
                                 & FOC_H, brent_tol, EndoHours(ai) )          
                      ENDIF
                      
                      EndoYgrid(ai) =  agrid(ai) +  EndoCons(ai)   &
                                          &                 -   psi*(yh(age, lambdai,ei)* EndoHours(ai))**(1.0_DP-tauPL)
                 ENDDO ! ai  

                tempai=1           
                DO WHILE ( YGRID(tempai,zi) .lt. EndoYgrid(1) )
                      tempai = tempai + 1
                ENDDO
                
                 ! Find  decision rules on exogenous grids
                 DO ai=tempai,na               
                         Cons(age, ai, zi, lambdai,ei)= Linear_Int(EndoYgrid,&
                                    & EndoCons,na, YGRID(ai,zi))   
                         
                         consin = Cons(age, ai, zi, lambdai,ei)
                                           
                         IF  (SolveInelasticLabor .ne. 1) THEN
                                 brentvalue = brent(0.000001_DP, 0.4_DP, 0.99_DP,&
                                     & FOC_H, brent_tol, Hours(age, ai, zi, lambdai,ei))           
                          ENDIF
 
                         Aprime(age, ai, zi, lambdai,ei) = YGRID(ai,zi)  &
                                        & + psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei))**(1.0_DP-tauPL)  & 
                                        & - Cons(age, ai, zi, lambdai,ei)   
                                        
                         If   (Aprime(age, ai, zi, lambdai,ei)  .lt. amin) then
                                   Aprime(age, ai, zi, lambdai,ei) = amin
                                 
                                   !compute  hours using FOC_HA                              
                                   ain = amin
         
                                  IF  (SolveInelasticLabor .ne. 1) THEN
                                        brentvalue = brent(0.000001_DP, 0.4_DP, 0.99_DP,&
                                             & FOC_HA, brent_tol, Hours(age, ai, zi, lambdai,ei))           
                                 ENDIF
                                 
                                 Cons(age,ai,zi,lambdai,ei)=YGRID(ai,zi)+psi*(yh(age, lambdai,ei)&
                                                    & * Hours(age, ai, zi, lambdai,ei))**(1.0_DP-tauPL) &
                                                    & - Aprime(age, ai, zi, lambdai,ei)      
      
                                  IF (Cons(age, ai, zi, lambdai,ei) .le. 0.0_DP)  THEN
                                        print*,'w1: Cons(age, ai, zi, lambdai,ei)=',Cons(age, ai, zi, lambdai,ei)
                                  ENDIF                   
                         endif      
                  ENDDO ! ai                  
                
                 ai=1           
                 DO WHILE ( YGRID(ai,zi) .lt. EndoYgrid(1) )
                 
                       brentvalue = brent( min(amin,YGRID(ai,zi))   ,  (amin+YGRID(ai,zi))/2.0_DP  ,  &
                             & YGRID(ai,zi)  + psi * ( yh(age, lambdai,ei)*0.95_DP )**(1.0_DP-tauPL)  ,  &
                             & FOC_WH, brent_tol, Aprime(age, ai, zi, lambdai,ei) )

                       !compute  hours using FOC_HA                              
                       ain = Aprime(age, ai, zi, lambdai,ei)
 
                       IF  (SolveInelasticLabor .ne. 1) THEN
                           brentvalue = brent(0.000001_DP, 0.4_DP, 0.99_DP,&
                                 & FOC_HA, brent_tol, Hours(age, ai, zi, lambdai,ei))           
                      ENDIF
   
                       Cons(age, ai, zi, lambdai,ei)=  YGRID(ai,zi)  &
                                        & + psi*(yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei))**(1.0_DP-tauPL)  &
                                        & - Aprime(age, ai, zi, lambdai,ei)
                       IF (Cons(age, ai, zi, lambdai,ei) .le. 0.0_DP)  THEN
                            print*,'w2:',age,zi,lambdai,ei,ai, 'Cons(age, ai, zi, lambdai,ei)=',Cons(age, ai, zi, lambdai,ei), &
                                & 'Aprime(age, ai, zi, lambdai,ei)=',Aprime(age, ai, zi, lambdai,ei), &
                                & 'Hours(age, ai, zi, lambdai,ei)=', Hours(age, ai, zi, lambdai,ei), &
                                & 'yh(age, lambdai,ei)=', yh(age, lambdai,ei),'YGRID(ai,zi)=',YGRID(ai,zi)
                            pause
                        ENDIF                   
                        ai = ai + 1
              ENDDO  
                 
            ENDDO  ! ei         
         ENDDO ! zi
    ENDDO ! lambdai
   
!    print*,'age=',age 
!    print*,'max deviation EXO Cons=',    MAXLOC(ABS(Cons(age, :, :, :,1)-AnRetCons(age, :, :, :) )), &
!                & MAXVAL(ABS(Cons(age, :, :, :,1)-AnRetCons(age, :, :, :) ))
!          
!    print*,'max % deviation EXO Cons=',    MAXLOC(ABS(100*(Cons(age, :, :, :,1)/AnRetCons(age, :, :, :)-1.0_DP) )), &
!                & MAXVAL(ABS(100*(Cons(age, :, :, :,1)/AnRetCons(age, :, :, :)-1.0_DP) ))
!    print*,Cons(age, :, 1, nlambda,ne)          
!    print*,AnRetCons(age, :, 1, nlambda)  
!!    print*,'' 
!    print*,'max deviation EXO ret Hours=',    MAXLOC(ABS(Hours(age, :, :, :,1)-AnRetHours(age, :, :, :) )), &
!                & MAXVAL(ABS(Hours(age, :, :, :,1)-AnRetHours(age, :, :, :) ))
!    print*,Hours(age, :, 1, nlambda,ne)
!    print*,AnRetHours(age, :, 1, nlambda)
!!    print*,''
!    pause

    ! Compute expected value function and expected value function derivative
!   DO lambdai=1,nlambda
!        DO zi=1,nz
!            DO ai=1,na
!                DO ei=1,ne
!                    ExpLValueP1(age,ai,zi,lambdai,ei) = SUM(pr_e(ei,:)/Cons(age, ai, zi, lambdai,:))
!                ENDDO
!            ENDDO
!        ENDDO
!    ENDDO    
ENDDO !age

cons = cons/(1.0_DP+tauC)

!print*,'min(Aprime) =', MINVAL(Aprime) 

END SUBROUTINE




!====================================================================
! THIS YGRID and Marginal Benefit of Investment GRID 
! NEEDS TO BE COMPUTED FOR EACH TIME THE INTEREST RATE "rr" IS UPDATED. 

SUBROUTINE FORM_Y_MB_GRID(TYGRID, TMBGRID)
USE GLOBAL
IMPLICIT NONE
!REAL(DP), INTENT(IN) :: rr
REAL(DP), DIMENSION(na,nz), INTENT(OUT)  :: TYGRID, TMBGRID
!integer :: ai, zi


!print*, 'TMBGRID(ai,zi)*beta*surP(1)='
DO ai=1,na
DO zi=1,nz
      TYGRID(ai,zi)   =  ( agrid(ai)+ ( rr * (zgrid(zi)*agrid(ai))**mu - DepRate*agrid(ai) ) *(1.0_DP-tauK) )*(1.0_DP-tauW)
      TMBGRID(ai,zi) =  ( 1.0_DP + ( rr *mu* (zgrid(zi)**mu) * (agrid(ai)**(mu-1.0_DP)) - DepRate ) *(1.0_DP-tauK) )*(1.0_DP-tauW)
ENDDO
!print*, TMBGRID(ai,:)*beta*survP(1)
!print*, 'TMBGRID(ai,nz)*beta*surP(age)=',TMBGRID(ai,nz)*beta*survP(1),TMBGRID(ai,nz)*beta*survP(11),TMBGRID(ai,nz)*beta*survP(21), &
!    & TMBGRID(ai,nz)*beta*survP(31), TMBGRID(ai,nz)*beta*survP(41),TMBGRID(ai,nz)*beta*survP(51),TMBGRID(ai,nz)*beta*survP(61)
ENDDO

!pause

END SUBROUTINE FORM_Y_MB_GRID

!====================================================================

! THIS SUBROUTINE COMPUTE EFFICIENCY UNITS OF LABOR DURING WORKING TIMES (NOT LABOR INCOME SINCE
! LABOR INCOME DEPENDS ON LABOR SUPPLY. IT ALSO fS RETIREMENT INCOME).
! THESE VARIABLES DEPEND ON EQUILIBRIUM WAGE AND EARNINGS AND AS A RESULT, WE NEED TO COMPUTE THESE
! FOR ALL EQUILIBRIUM WAGES AND EARNINGS

SUBROUTINE ComputeLaborUnits(Ebart,Waget)
USE GLOBAL
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

END SUBROUTINE 



!====================================================================

SUBROUTINE  INITIALIZE
use parameters
USE GLOBAL
use  programfunctions
IMPLICIT NONE
!integer::  lambdai
REAL(DP):: lambdaBAR, tempno 
REAL(DP) :: m, Rh, start_timet, finish_timet
INTEGER :: ee0, ee1, ee2, zindx1, zindx2, lambdaindx1, lambdaindx2, diff_array, eindx1, eindx2
INTEGER, DIMENSION(RetAge) :: agevec
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


	m=(amax-amin)**(1.0_DP/a_theta)/REAL(na-1,DP)
	DO ai=1,na
		agrid(ai)=REAL(ai-1,DP)*m
	END DO
	agrid=amin+agrid**a_theta

OPEN   (UNIT=2, FILE='agrid', STATUS='replace')
WRITE(unit=2, FMT=*) agrid
CLOSE (unit=2)

	m=(amax-amin)**(1.0_DP/a_theta)/REAL(fine_na-1,DP)
	DO ai=1,fine_na
		fine_agrid(ai)=REAL(ai-1,DP)*m
	END DO
	fine_agrid=amin+fine_agrid**a_theta
	
CALL tauchen(mtauchen,rho_z,sigma_z_eps,nz,zgrid,pr_z,Gz)
CALL tauchen(mtauchen,rho_E,sigma_e_eps,ne,egrid,pr_e,Ge)
CALL tauchen(mtauchen,rho_lambda,sigma_lambda_eps,nlambda,lambdagrid,pr_lambda,Glambda)

zgrid=exp(zgrid)
egrid=exp(egrid)
lambdagrid=exp(lambdagrid) 

!print*,'pr_z='
DO ee0 = 1,nz
    cdf_Gz(ee0) = sum(Gz(1:ee0))
    DO ee1 = 1,nz
        cdf_pr_z(ee0,ee1) = sum(pr_z(ee0,1:ee1))
    ENDDO
!    print*,pr_z(ee0,:) 
ENDDO
!print*,''
!print*,'Gz=',Gz
!print*,''
!PAUSE

!print*,'pr_lambda='
DO ee0 = 1,nlambda
    cdf_Glambda(ee0) = sum(Glambda(1:ee0))
    DO ee1 = 1,nlambda
        cdf_pr_lambda(ee0,ee1) = sum(pr_lambda(ee0,1:ee1))
    ENDDO
!    print*,pr_lambda(ee0,:) 
ENDDO
!print*,''
!print*,'Glambda=',Glambda
!print*,''
!PAUSE


DO ee0 = 1,ne
    cdf_Ge(ee0) = sum(Ge(1:ee0))
    DO ee1 = 1,ne
        cdf_pr_e(ee0,ee1) = sum(pr_e(ee0,1:ee1))
    ENDDO
!    print*,cdf_pr_e(ee0,:) ,pr_e(ee0,:) 
ENDDO

!   First construct dbn of "e" by age since we assume that in the first period, everyone starts at median "e"
!print*,'Ge_byage='
Ge_byage(:,:) =0.0_DP
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
    
DO age=1, MaxAge
      DO ee0 = 1, ne
            cdf_Ge_byage(age,ee0) = sum(Ge_byage(age,1:ee0))
      ENDDO
ENDDO

! this part computes efficiency units of labor, NOTE THAT THIS IS NOT LABOR INCOME
eff_un=0.0_DP
Rh=1.0_DP
DO age=1,RetAge
        agevec(age) = age
        kappagrid(age) = exp(  (60.0_DP *(Rh-1.0_DP)-(Rh-1.0_DP)**2.0_DP)/1800.0_DP )
!        kappagrid(age) = 1.0_DP
        DO lambdai=1,nlambda        
                DO ei=1,ne                
 
                     eff_un(age,lambdai,ei) =  kappagrid(age) *  &
                                                           & lambdagrid(lambdai) * egrid(ei)        
                ENDDO
        ENDDO
        Rh=Rh+1.0_DP
ENDDO

!print*,'kappagrid='
!print*,kappagrid-exp(  (60.0_DP *(agevec-1.0_DP)-(agevec-1.0_DP)**2.0_DP)/1800.0_DP )
!PAUSE

	!----------------------------------------------
	! life-cycle component
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
	

!pop=1.0_DP
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


CALL LIFETIME_Y_ESTIMATE

END SUBROUTINE INITIALIZE


!============================================================================


SUBROUTINE  LIFETIME_Y_ESTIMATE
use parameters
USE GLOBAL
use  programfunctions
IMPLICIT NONE
integer::  currentzi, currentlambdai, currentei
REAL(DP)::  tempnoage, tempnoz, tempnolambda, tempnoe, tempno, currenta, currentY
REAL(DP) ::  start_timet, finish_timet
INTEGER :: paneli
INTEGER , DIMENSION(10000000) :: panele, panellambda
REAL(DP) , DIMENSION(10000000) :: panel_lifetime_eff_unit
REAL(DP) :: mean_panel_lifetime_eff_unit , cohortsize, lambdaBAR
REAL(DP) , DIMENSION(nlambda,ne) :: lifetime_eff_unit_by_lambda_e, size_by_lambda_e 

cohortsize=10000000
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
lifetime_eff_unit_by_lambda_e =lifetime_eff_unit_by_lambda_e/ size_by_lambda_e
lifetime_eff_unit_by_lambda_e =lifetime_eff_unit_by_lambda_e/mean_panel_lifetime_eff_unit 


!------------------OLD RETIREMENT SCHEME  '
!lambdaBAR= sum(lambdagrid*Glambda)
!DO lambdai=1,nlambda
!DO ei=1,ne
!                lifetime_eff_unit_by_lambda_e(lambdai, ei) = lambdagrid(lambdai)/ lambdaBAR
!ENDDO
!ENDDO
!print*,''
!print*,'  OLD RETIREMENT SCHEME  '
!print*,''

!print*,'LIFETIME AVERAGE Y RELATIVE TO THE MEAN LIFETIME AVERAGE Y'
!DO lambdai=1,nlambda
!                print*, lifetime_eff_unit_by_lambda_e(lambdai, :)
!ENDDO
!PRINT*,''
!DO lambdai=1,nlambda
!                print*, size_by_lambda_e(lambdai, :)
!ENDDO
!
!PRINT*,''
!print*,'egrid=',egrid
!
!PRINT*,''
!print*,'lambdagrid=',lambdagrid


! compute retirement income replacement rate for each lambda type
!print*,'phi_lambda_e='
DO lambdai=1,nlambda
DO ei=1,ne
      IF (lifetime_eff_unit_by_lambda_e(lambdai, ei) .le. 0.3_DP) THEN
          phi_lambda_e(lambdai,ei) = 0.9_DP*  lifetime_eff_unit_by_lambda_e(lambdai, ei)
          ELSEIF   (lifetime_eff_unit_by_lambda_e(lambdai, ei) .le. 2.0_DP) THEN
	           phi_lambda_e(lambdai,ei)  = 0.27_DP +  0.32_DP*  (lifetime_eff_unit_by_lambda_e(lambdai, ei)-0.3_DP) 
	           ELSEIF (lifetime_eff_unit_by_lambda_e(lambdai, ei) .le. 4.1_DP) THEN
	                  phi_lambda_e(lambdai,ei)  = 0.81_DP +  0.15_DP*  (lifetime_eff_unit_by_lambda_e(lambdai, ei)-2.0_DP)
                          ELSE
                                phi_lambda_e(lambdai,ei)  =1.1_DP	   
        ENDIF
ENDDO 
!print*,phi_lambda_e(lambdai,:)       
ENDDO
!PAUSE

END SUBROUTINE


!====================================================================


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
	a=(1.0_DP-rhot)*mut;
        
        if (nt .gt. 1) then  
	        gridt(nt)=mt*sqrt(sigmat**2.0_DP/(1.0_DP-rhot**2))
	        gridt(1)=-gridt(nt)
	        ELSE
	        print*,'only one grid. For this case transition probabilities might not be right. check it.'
	 endif       

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
! 	print*, 'Gt', Gt, 'sum', sum(Gt)

	DO i=1,1000
		Gt_new=0.0_DP
		DO zi=1,nt
			DO zz=1,nt
				Gt_new(zz)=Gt_new(zz)+Gt(zi)*prt(zi,zz)
			END DO
		END DO
		Gt=Gt_new
	END DO
	
! 	print*, 'Gt', Gt, 'sum', sum(Gt)
! 	!pause

  if (nt .eq. 1) then
      Gt=1.0_DP
  endif
  
	
END SUBROUTINE

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
!
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

INTEGER :: n,NMAX  
      REAL(DP)   ::  yp1,ypn,x(n),y(n),y2(n), y1(n)  
      PARAMETER (NMAX=1000)  
      INTEGER :: i,k  
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
      
  
END  
      
!=======================================================================
       
 
SUBROUTINE splint(xa,ya,y2a,n,x,y)  
! Given arrays xa(1:n) and ya(1:n) of length n; and second derivatives y2a(1:n), which is the output of spline
! given the value of "x" this subroutine returns  a cubic-spline interpolated value "y" at point "x"

USE parameters
IMPLICIT NONE      
      INTEGER:: n  
      REAL(DP)   :: x, y, yprime, xa(n),y2a(n),ya(n)  
      INTEGER k,khi,klo  
      REAL(DP)   :: a,b,h  
      
      
      
!      y2a=0.0_DP
      
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
      if (h.eq.0.) pause 'bad xa input in splint'  
      a=(xa(khi)-x)/h  
      b=(x-xa(klo))/h  
      y=a*ya(klo)+b*ya(khi)     +((a**3.0_DP-a)*y2a(klo)+(b**3.0_DP-b)*y2a(khi))*(h**2.0_DP)/6.0_DP  
! it also returns "yprime" (first derivative) at point "x"
!      yprime = ( ya(khi) - ya(klo) ) / h  & 
!                 & + h * ( ( 1.0_DP- 3.0_DP * a**2.0_DP  ) *  y2a(klo)  +  ( 3.0_DP * b**2.0_DP -1.0_DP ) *  y2a(khi)  )/6.0_DP
      return  

END  
 
!=======================================================================




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
