MODULE parameters

    use nrtype
    use nrutil

    INTEGER(I4B),  PARAMETER ::  KeepSSatBench=1 , update_period=5 
    
! REAL(DP), PARAMETER :: beta=0.9515_DP
! REAL(DP), PARAMETER :: phi=1.25_DP ! parameter for leisure utility ! when phi=0, do not value leisure.
!    REAL(DP), PARAMETER  :: rho_z= 0.50651_DP, sigma_z_eps=0.624_DP  ! parameters for z process
! REAL(DP), PARAMETER  :: rho_lambda=0.5_DP, sigma_lambda_eps=0.34_DP 
        REAL(DP), PARAMETER  :: rho_lambda=0.5_DP, sigma=4.0_DP
  REAL(DP) :: beta, phi, sigma_lambda_eps,  rho_z, sigma_z_eps   , mu_z, gamma 
  REAL(DP), PARAMETER  :: rho_e=0.9_DP, sigma_e_eps=0.20_DP    
  
  ! production  
  REAL(DP), PARAMETER :: alpha=0.33_DP, Aprod=1.0_DP, DepRate=0.00_DP
  REAL(DP), PARAMETER :: mu=0.9_DP, hp=0.002_DP ! hp=home production (min cons)
    
  INTEGER(I4B), PARAMETER       :: MaxAge=81, RetAge=45 , totpop=4200000

  ! asset grid
  REAL(DP), PARAMETER        ::  a_theta=4.0_DP , amax=100000.0_DP, amin=0.0001_DP, azero=0.0001_DP
  INTEGER(I4B),  PARAMETER :: na=201, fine_na=801, MaxSimuTime=2000 ,  max_age_category=7 



   ! Financial constraint parameter
  REAL(DP), PARAMETER :: vartheta=5.0_DP 


    
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
  
  REAL(DP), PARAMETER  :: tauWmin=0.005_DP, tauWinc=0.005_DP
  
END MODULE parameters


MODULE global
    USE parameters
   real(DP) :: tauC
    real(DP) , dimension(6) :: params

    REAL(DP), DIMENSION(MaxAge) :: pop, survP          ! population size by age and survival prob by age
  REAL(DP), DIMENSION(RetAge)                     :: kappagrid

  REAL(DP), DIMENSION(ne)                            :: egrid, Ge , cdf_Ge
  REAL(DP), DIMENSION(ne,ne)                       :: pr_e, cdf_pr_e
    REAL(DP), DIMENSION(MaxAge,ne)              :: Ge_byage, cdf_Ge_byage
                
  REAL(DP), DIMENSION(nlambda)                  :: lambdagrid,  Glambda, cdf_Glambda
  REAL(DP), DIMENSION(nlambda,nlambda)    :: pr_lambda  , cdf_pr_lambda
    REAL(DP), DIMENSION(nlambda,ne)                  :: phi_lambda_e, RetY_lambda_e

  ! investment abilities
  REAL(DP), DIMENSION(nz)                  :: zgrid , Gz, cdf_Gz
  REAL(DP), DIMENSION(nz,nz)             :: pr_z, cdf_pr_z   
                                                                                   ! I keep this format because I do not want to mess the ages
  REAL(DP), DIMENSION(MaxAge  , nlambda, ne) :: eff_un,  yh       ! Labor efficiency units, Labor efficiency units x wage

!    REAL(DP), DIMENSION(MaxAge, na, nz, nlambda) :: RetCons !, RetConsP1
    REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) :: Cons, Hours, Aprime
  REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) ::ValueFunction !, ValueFunction_Bench, ValueFunction_Exp, Cons_Eq_Welfare

! Analytical solution for mu=1 for all the lifecycle not just retirement period
    REAL(DP), DIMENSION(MaxAge, na, nz, nlambda) :: AnRetCons,   AnRetValue, AnRetHours 
 
 
! ASSET AND RESOURCE GRIDS
! NOTE THAT THESE GRIDS WILL CHANGE AS WE CHANGE TAX RATES AND FOR EACH INTERMEDIATE GOOD PRICE
! FOR EACH TAX RATE WE WILL HAVE DIFFERENT Y GRIDS
  REAL(DP), DIMENSION(na)                            :: agrid
  REAL(DP), DIMENSION(na,nz)                       :: YK, YW, MBK, MBW ! Y grids for CAPITAL INCOME tax and WEALTH taxes. 
    REAL(DP), DIMENSION(na,nz)                       :: YGRID, MBGRID, MBGRID_bench, MBGRID_exp
  
  REAL(DP), DIMENSION(fine_na)                     :: fine_agrid
    REAL(DP) :: tauk_bench, tauw_bench, tauL_bench
    REAL(DP) :: tauk_exp,     tauw_exp,    tauL_exp



    REAL(DP) :: tauK, tauW, tauL , tauWindx, tauW_low, tauW_up,  Ebar , wage
    REAL(DP) :: QBAR_bench,  QBAR_exp, NBAR_exp, NBAR_bench, Ebar_bench, Ebar_exp, wage_bench, wage_exp, BBAR 
    REAL(DP) :: rr, rr_bench, rr_exp, int_rate, int_rate_bench, int_rate_exp, External_Debt_GDP, BBAR_POS 



    INTEGER :: age, lambdai, zi, ai, ei    

    real(DP) :: NBAR, QBAR, CE_NEWBORN
    
    REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) ::DBN1, DBN_bench
    REAL(DP), DIMENSION(na) :: pr_a_dbn, cdf_a_dbn, tot_a_by_grid, cdf_tot_a_by_grid
    REAL(DP), DIMENSION(100) ::  cdf_tot_a_by_prctile
    INTEGER, DIMENSION(100) ::  prctile_ai
    REAL(DP) :: pop_25_60 , tothours_25_60, pop_pos_earn_25_60, tot_log_earnings_25_60, mean_log_earnings_25_60 
    REAL(DP) :: meanhours_25_60, Var_Log_Earnings_25_60, Std_Log_Earnings_25_60, MeanWealth, Wealth_Output
    REAL(DP) :: prct1_wealth, prct10_wealth, prct20_wealth, prct40_wealth, Bequest_Wealth
    REAL(DP) :: YBAR, GBAR, GBAR_exp,GBAR_exp_old, GBAR_bench, Y_bench, Y_exp

    REAL(DP) :: SSE_Moments, Min_SSE_Moments, MeanReturn, StdReturn,  MeanCons
    REAL(DP) , DIMENSION(6) ::  Min_Moments
    INTEGER:: solving_bench, newiseed
    
    real(DP) :: SSC_Payments, GBAR_K, GBAR_L, GBAR_C
 
    REAL(DP) , DIMENSION(nz) ::  kstar, pistar
    REAL(DP) , DIMENSION(na, nz) ::  kdemand, bdemand
END MODULE global

!====================================================================
!====================================================================

Module programfunctions
use parameters
use global
    
contains

Function Agg_Debt(R_in)
    use NRTYPE
    use parameters
    use global
    Implicit None 
    real(dp), intent(in) :: R_in
    real(dp)             :: Agg_Debt
    integer              :: z1, a1

    Agg_Debt =0.0_DP          
            DO z1=1,nz
            DO a1=1,na
              Agg_Debt = Agg_Debt + sum(DBN1(:, a1, z1, :, :))* &
                & ( min((mu*rr*zgrid(z1)**mu/(R_in+DepRate))**(1.0_DP/(1.0_DP-mu)),vartheta*agrid(a1)) &
                & - agrid(a1)) 
            ENDDO
            ENDDO

    Agg_Debt = Agg_Debt/MeanWealth 

  end Function Agg_Debt



FUNCTION FOC_R(aprimet)
IMPLICIT NONE   
real(DP), intent(in) :: aprimet
real(DP) :: MBaprime, FOC_R, yprime, cprime

if (kstar(zi) .le. (vartheta*aprimet) ) then
   MBaprime = (1.0_DP+int_rate*(1.0_DP-tauK))*(1.0_DP-tauW) 
   yprime   = (1.0_DP+int_rate*(1.0_DP-tauK))*(1.0_DP-tauW)*aprimet + pistar(zi)*(1.0_DP-tauK)*(1.0_DP-tauW) 
else
   MBaprime = (1.0_DP+int_rate*(1.0_DP-tauK))*(1.0_DP-tauW) &
      & + (rr*mu*((vartheta*zgrid(zi))**mu)*aprimet**(mu-1.0_DP)-(int_rate+DepRate)*vartheta)*(1.0_DP-tauK)*(1.0_DP-tauW)
   yprime   = (1.0_DP+int_rate*(1.0_DP-tauK))*(1.0_DP-tauW)*aprimet &
      & + (rr*(vartheta*zgrid(zi)*aprimet)**mu-(int_rate+DepRate)*vartheta*aprimet)*(1.0_DP-tauK)*(1.0_DP-tauW)
endif

cprime =   Linear_Int(Ygrid(:,zi), &
                                    & Cons(age+1,:,zi,lambdai,ei),na, yprime)    

FOC_R   =  ( YGRID(ai,zi) + RetY_lambda_e(lambdai,ei) - aprimet -  cprime * & 
        &  ( beta*survP(age)* MBaprime  ) ** ( 1 / (gamma*(1.0_DP-sigma)-1) ) )**2.0_DP 

END  FUNCTION


!====================================================================

FUNCTION FOC_W(aprimet)
IMPLICIT NONE   
real(DP), intent(in) :: aprimet
real(DP) :: ntemp, ctemp, MBaprime, FOC_W,  yprime 
integer:: epindx
real(DP), dimension(ne):: cprime, nprime

ntemp = max(0.0_DP, gamma - (1.0_DP-gamma)*(YGRID(ai,zi) - aprimet)/yh(age, lambdai,ei) )
ctemp = YGRID(ai,zi) + yh(age, lambdai,ei) * ntemp - aprimet

if (kstar(zi) .le. (vartheta*aprimet) ) then
   MBaprime = (1.0_DP+int_rate*(1.0_DP-tauK))*(1.0_DP-tauW) 
   yprime   = (1.0_DP+int_rate*(1.0_DP-tauK))*(1.0_DP-tauW)*aprimet + pistar(zi)*(1.0_DP-tauK)*(1.0_DP-tauW) 
else
   MBaprime = (1.0_DP+int_rate*(1.0_DP-tauK))*(1.0_DP-tauW) &
      & + (rr*mu*((vartheta*zgrid(zi))**mu)*aprimet**(mu-1.0_DP)-(int_rate+DepRate)*vartheta)*(1.0_DP-tauK)*(1.0_DP-tauW)
   yprime   = (1.0_DP+int_rate*(1.0_DP-tauK))*(1.0_DP-tauW)*aprimet &
      & + (rr*(vartheta*zgrid(zi)*aprimet)**mu-(int_rate+DepRate)*vartheta*aprimet)*(1.0_DP-tauK)*(1.0_DP-tauW)
endif

DO epindx=1,ne
      cprime(epindx) = Linear_Int(Ygrid(:,zi),&
                                    & Cons(age+1,:,zi,lambdai,epindx), na,    yprime  )
      nprime(epindx) = 1.0_DP-(1.0_DP-gamma)*cprime(epindx)/(gamma*yh(age, lambdai,ei))                            
ENDDO

nprime = max(0.0_DP,nprime)

FOC_W = ((ctemp**(gamma*(1.0_DP-sigma)-1))*((1.0_DP-ntemp)**((1.0_DP-gamma)*(1.0_DP-sigma))) &
& - beta*survP(age)*MBaprime  &
& * sum( pr_e(ei,:) * (cprime**(gamma*(1.0_DP-sigma)-1.0_DP)) &
& *((1.0_DP-nprime)**((1.0_DP-gamma)*(1.0_DP-sigma)))) )**2.0_DP

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
      !if (h.eq.0.) pause 'bad xa input in linear int'  
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
      !if (h.eq.0.) pause 'bad xa input in linear int'  
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



end module programfunctions

!====================================================================




PROGRAM main
USE parameters
USE GLOBAL
use  programfunctions
IMPLICIT NONE
real(DP):: start_time, finish_time, betaL, betaH, sigmalambL,sigmalambH, sigmazL, sigmazH, rhozL, rhozH, gammaL, gammaH, muzL, muzH
integer :: parindx1,  parindx2, parindx3, parindx4, parindx5, parindx6, nbeta, ngamma, nsigmalambda, nrhoz, nsigmaz, nmuz



call cpu_time(start_time) 
tauK = 0.25_DP
tauL = 0.224_DP
tauW= 0.00_DP
tauC = 0.075_DP

! ------- DO NOT REMOVE THE LINES BELOW

rr=  4.906133597851297E-002 
wage=  1.97429920063330 
Ebar=  1.82928004963637
Ebar_bench = Ebar
int_rate =0.000000001_DP


! ------- DO NOT REMOVE THE LINES ABOVES

! the following solves equilibrium of capital tax economy
print*,'---- CRRA + FLAT LABOR INCOME TAX SOLVING  + FINANCIAL MARKETS ------ '

if (KeepSSatBench .eq. 1) then
    print*,''
    print*,'SS kept at benchmark level'
    print*,''
else
    print*,''
    print*,'SS adjusted with tax change'
    print*,''
endif    
if (DepRate .gt. 0.0 ) then
        print*,'Pos Dep Rate=',DepRate
endif
print*,'VARTHETA=',vartheta
print*,'tauC=',tauC
print*,'tauK=',tauK
print*,'tauL=',tauL


Params=[ 0.944,  0.00, 0.50,  0.60,  0.35,  0.4494]
Params=[  0.9438,  0.00,  0.50,  0.63,  0.35,  0.4494 ] ! tauL=0.3, tauC=0.075 calibration
Params=[  0.9436,  0.00, 0.50,  0.70444445,   0.34,  0.4494 ] ! tauL=0.224, tauC=0.075 calibration
!Params=[  0.944,  0.00,  0.70,  0.523684194,  0.34,  0.4494]
!params=[0.9321,  0.00,  0.50, 0.00,  0.2565,  0.4866] ! equal z calibration with zero depreciation


!Params=[  0.969,  7.7777777777777,  0.5,  0.704444468021393,  0.34,  0.4494 ] ! calibration that matches W/Y and MeanReturn
!Params=[  0.968500003218651,  6.11111111111111,  0.7,  0.523684203624725,  0.340000003576279,  0.449400007724762 ]
!Params=[  0.969,  0.0,  0.5,  0.704444468021393,  0.34,  0.4494 ] 
!Params=[  0.968500003218651,  0.0,  0.7,  0.523684203624725,  0.340000003576279,  0.449400007724762 ]

!----------------------- calibration starts--------------------------------------
!
!!rho_z=0.3 -- sigma_z=0.77
!!rho_z=0.4 -- sigma_z=0.74
!!rho_z=0.6 -- sigma_z=0.65
!!rho_z=0.7 -- sigma_z=0.57
!!rho_z=0.8 -- sigma_z=0.49
!
!betaL=0.9322
!betaH=0.9322
!
!gammaL=0.46
!gammaH=0.5
!
!sigmalambL= sigma_lambda_eps
!sigmalambH= sigma_lambda_eps
!
!muzL = mu_z
!muzH = mu_z
!
!rhozL = rho_z
!rhozH= rho_z
! 
!sigmazL = 0.0
!sigmazH = 0.0
!
!
!nbeta =1
!nrhoz=1
!nsigmaz=1
!nsigmalambda=1
!ngamma=4
!nmuz=1
!
!Min_SSE_Moments=1000.0_DP
!
!DO parindx2=1,ngamma
!DO parindx1=1,nbeta
!DO parindx4=1,nrhoz
!DO parindx5=1,nsigmaz
!DO parindx6=1,nmuz
!DO parindx3=1,nsigmalambda
!
!    beta  = betaL  + real(parindx1-1,8) *(betaH-betaL)/max(real(nbeta-1,8),1.0_DP)
!    gamma = gammaL + real(parindx2-1,8) *(gammaH-gammaL)/max(real(ngamma-1,8),1.0_DP)
!    sigma_lambda_eps = sigmalambL + real(parindx3-1,8)*(sigmalambH -sigmalambL) / max(real(nsigmalambda-1,8),1.0_DP)
!    rho_z= rhozL   +  real(parindx4-1,8)*(rhozH-rhozL) / max(real(nrhoz-1,8),1.0_DP)
!    sigma_z_eps = sigmazL +  real(parindx5-1,8)*(sigmazH-sigmazL) / max(real(nsigmaz-1,8),1.0_DP)
!    mu_z = muzL +  real(parindx6-1,8)*(muzH-muzL) / max(real(nmuz-1,8),1.0_DP)
!        
!    print*,'parameters=',beta, mu_z, rho_z, sigma_z_eps, sigma_lambda_eps , gamma 
!
!    solving_bench=1
!    CALL  INITIALIZE
!    CALL FIND_DBN_EQ
!    print*,'parameters=',beta, mu_z, rho_z,sigma_z_eps,sigma_lambda_eps, gamma
!    CALL COMPUTE_STATS
! 
!    IF (SSE_Moments .lt. Min_SSE_Moments ) THEN
!        Min_SSE_Moments =SSE_Moments
!        params= [ beta, mu_z, rho_z, sigma_z_eps, sigma_lambda_eps, gamma ]
!        Min_Moments = [  Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, MeanReturn  ]
!    ENDIF
!    !CALL WRITE_TO_FILE
!    
!ENDDO
!ENDDO
!ENDDO
!ENDDO
!ENDDO
!ENDDO
!print*, params
!print*, Min_Moments

!----------------------- calibration ends --------------------------------------


beta=params(1)
mu_z=params(2)
rho_z=params(3)
sigma_z_eps =params(4)
sigma_lambda_eps = params(5)
gamma= params(6)

OPEN   (UNIT=3, FILE='params', STATUS='replace')
WRITE(unit=3, FMT=*) params
CLOSE (unit=3)

print*,''
print*,''
print*,'---------------------------------------------------'
print*,''
print*,'CAPITAL TAX ECONOMY nnnnnnnnn'
print*,'parameters=', params 

solving_bench=1
CALL  INITIALIZE
CALL FIND_DBN_EQ
CALL COMPUTE_STATS
CALL GOVNT_BUDGET

MBGRID_bench=MBGRID

GBAR_bench=GBAR
QBAR_bench = QBAR 
NBAR_bench = NBAR 
Ebar_bench  = EBAR
rr_bench      = rr
int_rate_bench = int_rate
wage_bench = wage
Y_bench       = YBAR
tauK_bench = tauK
tauw_bench = tauW
tauL_bench  = tauL
DBN_bench = DBN1
CALL WRITE_VARIABLES(1)

OPEN   (UNIT=3, FILE='panela_bench', STATUS='replace')
OPEN   (UNIT=4, FILE='panelage_bench', STATUS='replace')
OPEN   (UNIT=5, FILE='panelz_bench', STATUS='replace')
OPEN   (UNIT=6, FILE='panellambda_bench', STATUS='replace')
OPEN   (UNIT=7, FILE='panele_bench', STATUS='replace')
OPEN   (UNIT=8, FILE='panel_return_bench', STATUS='replace')
OPEN   (UNIT=9, FILE='panel_cons_bench', STATUS='replace') 
OPEN   (UNIT=10, FILE='panel_hours_bench', STATUS='replace') 
OPEN   (UNIT=11, FILE='panel_aprime_bench', STATUS='replace') 
OPEN   (UNIT=12, FILE='panel_at_return_bench', STATUS='replace') 

call SIMULATION

close (unit=3)
close (unit=4)
close (unit=5)
close (unit=6)
close (unit=7)
close (unit=8)
close (unit=9)
close (unit=10)
close (unit=11)
close (unit=12)


 !the following solves equilibrium of WEALTH tax economy
print*,'WEALTH TAX ECONOMY'
!-------DO NOT ERASE THE LINE BELOW--------
solving_bench=0
!-------DO NOT ERASE THE LINE ABOVE--------
 
tauK = 0.0_DP 
 
! policy experiment
!tauK = 0.0_DP
!tauW=  1.719379006132605E-002  
!tauL = 0.224_DP
 
!! Optimal capital income tax rates
!tauK=  1.055728085076737E-002 
!tauW=  0.000000000000000 
!tauL=  0.295893735802321 
!
!! Optimal wealth tax rates
!tauK=  0.000000000000000 
!tauW=  1.661803400000000E-002 
!tauL=  0.224762614881983


! =========solving for the tauW that generates the same G

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
      print*,'tauW_low =', tauW_low, 'tauW_up=', tauW_up, 'tauW=', tauW, 'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
ENDDO

tauW_up  = tauW 
tauW_low = tauW  -  tauWinc  
print*, 'tauW=', tauW, 'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench,'tauW_low =', tauW_low, 'tauW_up=', tauW_up
tauW = tauW_low + tauWinc * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)

CALL FIND_DBN_EQ
CALL COMPUTE_STATS
CALL GOVNT_BUDGET
GBAR_exp = GBAR

DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.001 ) ! as long as the difference is greater than 0.1% continue
    if (GBAR_exp .gt. GBAR_bench ) then
        tauW_up  = tauW 
    else
        tauW_low = tauW 
    endif
    tauW = (tauW_low + tauW_up)/2.0_DP
    CALL FIND_DBN_EQ
    CALL COMPUTE_STATS
    CALL GOVNT_BUDGET
    GBAR_exp = GBAR
    print*, 'tauW=', tauW, 'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench,'tauW_low =', tauW_low, 'tauW_up=', tauW_up
ENDDO


! =========solving for the tauW that generates the same G ENDED


CALL FIND_DBN_EQ
CALL COMPUTE_STATS
CALL GOVNT_BUDGET
MBGRID_exp=MBGRID

GBAR_exp=GBAR
QBAR_exp = QBAR 
NBAR_exp = NBAR  
Y_exp = YBAR
Ebar_exp  = EBAR
rr_exp    = rr
int_rate_exp = int_rate
wage_exp = wage
tauW_exp = tauW
tauK_exp = tauK
tauL_exp = tauL
CALL WRITE_VARIABLES(0)

OPEN   (UNIT=3, FILE='panela_exp', STATUS='replace')
OPEN   (UNIT=4, FILE='panelage_exp', STATUS='replace')
OPEN   (UNIT=5, FILE='panelz_exp', STATUS='replace')
OPEN   (UNIT=6, FILE='panellambda_exp', STATUS='replace')
OPEN   (UNIT=7, FILE='panele_exp', STATUS='replace')
OPEN   (UNIT=8, FILE='panel_return_exp', STATUS='replace')
OPEN   (UNIT=9, FILE='panel_cons_exp', STATUS='replace') 
OPEN   (UNIT=10, FILE='panel_hours_exp', STATUS='replace') 
OPEN   (UNIT=11, FILE='panel_aprime_exp', STATUS='replace') 
OPEN   (UNIT=12, FILE='panel_at_return_bench', STATUS='replace') 

call SIMULATION

close (unit=3)
close (unit=4)
close (unit=5)
close (unit=6)
close (unit=7)
close (unit=8)
close (unit=9)
close (unit=10)
close (unit=11)
close (unit=12)

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
REAL(DP), dimension(max_age_category,nz) :: tot_cons_by_agegroup_by_z_bench,  tot_hours_by_agegroup_by_z_bench , &
       &tot_aprime_by_agegroup_by_z_bench,  tot_at_return_by_agegroup_by_z_bench ,tot_cap_tax_by_agegroup_by_z_bench ,&
       & tot_lab_tax_by_agegroup_by_z_bench , tot_cons_tax_by_agegroup_by_z_bench
REAL(DP), dimension(max_age_category,nz) :: tot_cons_by_agegroup_by_z_exp , tot_hours_by_agegroup_by_z_exp  , &
       & tot_aprime_by_agegroup_by_z_exp, tot_cons_tax_by_agegroup_by_z_exp, &
       & tot_at_return_by_agegroup_by_z_exp  , tot_cap_tax_by_agegroup_by_z_exp ,tot_lab_tax_by_agegroup_by_z_exp 


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

OPEN (UNIT=6, FILE='mean_wealth_by_agegroup_z_bench.txt', STATUS='replace')  
OPEN (UNIT=7, FILE='size_by_agegroup_z_bench.txt', STATUS='replace')  
DO age_group_counter=1,max_age_category
    WRITE  (UNIT=6, FMT=*)   tot_wealth_by_agegroup_z_bench(age_group_counter,:)/ size_by_agegroup_z_bench(age_group_counter,:)
    WRITE  (UNIT=7, FMT=*)    size_by_agegroup_z_bench(age_group_counter,:)
ENDDO
close (UNIT=6)
close (UNIT=7)


age_group_counter=1
tot_cons_by_agegroup_by_z_bench=0.0_DP
tot_hours_by_agegroup_by_z_bench =0.0_DP
tot_aprime_by_agegroup_by_z_bench =0.0_DP
tot_at_return_by_agegroup_by_z_bench =0.0_DP
tot_cap_tax_by_agegroup_by_z_bench=0.0_DP
tot_lab_tax_by_agegroup_by_z_bench=0.0_DP
tot_cons_tax_by_agegroup_by_z_bench=0.0_DP

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

         tot_at_return_by_agegroup_by_z_bench(age_group_counter,zi) = tot_at_return_by_agegroup_by_z_bench(age_group_counter,zi)+&
                         & MBGRID_bench(ai,zi)*agrid(ai)*DBN_bench(age,ai,zi,lambdai,ei)
       
         tot_lab_tax_by_agegroup_by_z_bench(age_group_counter,zi) =tot_lab_tax_by_agegroup_by_z_bench(age_group_counter,zi) + &
                         & DBN_bench(age,ai,zi,lambdai,ei) * ( &
                          & tauL_bench * wage_bench * eff_un(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei) )
      
         tot_cons_tax_by_agegroup_by_z_bench(age_group_counter,zi)=tot_cons_tax_by_agegroup_by_z_bench(age_group_counter,zi)+&
                         & DBN_bench(age,ai,zi,lambdai,ei) * ( &
                          & tauC *cons(age,ai,zi,lambdai,ei) )
      
         
             ENDDO
         ENDDO
     ENDDO
   ENDDO
ENDDO

OPEN (UNIT=6, FILE='mean_cons_by_agegroup_by_z_bench.txt', STATUS='replace')  
OPEN (UNIT=7, FILE='mean_hours_by_agegroup_by_z_bench.txt', STATUS='replace')  
OPEN (UNIT=8, FILE='mean_aprime_by_agegroup_by_z_bench.txt', STATUS='replace')  
OPEN (UNIT=9, FILE='mean_weighted_at_return_by_agegroup_by_z_bench.txt', STATUS='replace')  
!OPEN (UNIT=10, FILE=' tot_cap_tax_by_agegroup_by_z_bench.txt', STATUS='replace')  
OPEN (UNIT=11, FILE=' tot_lab_tax_by_agegroup_by_z_bench.txt', STATUS='replace')  
OPEN (UNIT=12, FILE=' tot_cons_tax_by_agegroup_by_z_bench.txt', STATUS='replace')  
DO age_group_counter=1,max_age_category
WRITE  (UNIT=6, FMT=*)  tot_cons_by_agegroup_by_z_bench(age_group_counter,:)/ size_by_agegroup_z_bench(age_group_counter,:)
WRITE  (UNIT=7, FMT=*)  tot_hours_by_agegroup_by_z_bench(age_group_counter,:)/ size_by_agegroup_z_bench(age_group_counter,:)
WRITE  (UNIT=8, FMT=*)  tot_aprime_by_agegroup_by_z_bench(age_group_counter,:)/ size_by_agegroup_z_bench(age_group_counter,:)
WRITE (UNIT=9, FMT=*) tot_at_return_by_agegroup_by_z_bench(age_group_counter,:)/ tot_wealth_by_agegroup_z_bench(age_group_counter,:)
!WRITE  (UNIT=10, FMT=*) tot_cap_tax_by_agegroup_by_z_bench(age_group_counter,:)  
WRITE  (UNIT=11, FMT=*) tot_lab_tax_by_agegroup_by_z_bench(age_group_counter,:) 
WRITE  (UNIT=12, FMT=*) tot_cons_tax_by_agegroup_by_z_bench(age_group_counter,:) 
ENDDO
close (UNIT=6)
close (UNIT=7)
close (UNIT=8)
close (UNIT=9)
!close (UNIT=10)
close (UNIT=11)
close (UNIT=12)
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


OPEN (UNIT=5, FILE='CE_alternative', STATUS='replace') 
    WRITE  (UNIT=5, FMT=*) 'CE Newborn=', ( sum(ValueFunction_exp(1,:,:,:,:)*DBN1(1,:,:,:,:)) / &
                                &  sum(ValueFunction_Bench(1,:,:,:,:)*DBN_bench(1,:,:,:,:))  ) &
                                &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP
    WRITE  (UNIT=5, FMT=*) 'CE =', ( sum(ValueFunction_exp(:,:,:,:,:)*DBN1(:,:,:,:,:)) / &
                                &  sum(ValueFunction_Bench(:,:,:,:,:)*DBN_bench(:,:,:,:,:))  ) &
                                &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP

    WRITE  (UNIT=5, FMT=*) 'Av Utility NB (bench) =',sum(ValueFunction_Bench(1,:,:,:,:)*DBN_bench(1,:,:,:,:))  
    WRITE  (UNIT=5, FMT=*) 'Av Utility NB (exp  ) =',sum(ValueFunction_exp(1,:,:,:,:)*DBN1(1,:,:,:,:))    
    WRITE  (UNIT=5, FMT=*) 'Av Utility  (bench)   =',sum(ValueFunction_Bench(:,:,:,:,:)*DBN_bench(:,:,:,:,:))  
    WRITE  (UNIT=5, FMT=*) 'Av Utility  (exp)     =',sum(ValueFunction_exp(:,:,:,:,:)*DBN1(:,:,:,:,:))
close (unit=5)

OPEN (UNIT=5, FILE='CE_NEWBORN', STATUS='replace')  
OPEN (UNIT=6, FILE='CE', STATUS='replace')  
OPEN (UNIT=7, FILE='CE_by_age', STATUS='replace')  
OPEN (UNIT=8, FILE='CE_by_age_z', STATUS='replace')  

DO age=1,MaxAge
  Cons_Eq_Welfare(age,:,:,:,:)=(ValueFunction_exp(age,:,:,:,:)/ValueFunction_Bench(age,:,:,:,:)) &
                                &  ** ( 1.0_DP / ( gamma* (1.0_DP-sigma)) )-1.0_DP

    WRITE  (UNIT=7, FMT=*) 100*sum(Cons_Eq_Welfare(age,:,:,:,:)*DBN_bench(age,:,:,:,:))/sum(DBN_bench(age,:,:,:,:))

    DO zi=1,nz
         temp_ce_by_z(zi) = 100*sum(Cons_Eq_Welfare(age,:,zi,:,:)*DBN_bench(age,:,zi,:,:))/sum(DBN_bench(age,:,zi,:,:))
    ENDDO
    WRITE  (UNIT=8, FMT=*) temp_ce_by_z
!    print*,'age=',age, temp_ce_by_z, ', mean:  ', &
!        & 100*sum(Cons_Eq_Welfare(age,:,:,:,:)*DBN_bench(age,:,:,:,:))/sum(DBN_bench(age,:,:,:,:))
ENDDO

CE_NEWBORN = 100.0_DP*sum(Cons_Eq_Welfare(1,:,:,:,:)*DBN_bench(1,:,:,:,:))/sum(DBN_bench(1,:,:,:,:))

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

OPEN (UNIT=8, FILE='CE_by_Asset_z_med_E_Lambda_age1', STATUS='replace') 
OPEN (UNIT=9, FILE='MeanAsset_by_z_med_E_Lambda_age1', STATUS='replace') 
DO zi=1,nz
           WRITE  (UNIT=8, FMT=*)   Cons_Eq_Welfare(1,:,zi,nlambda/2+1, ne/2+1)
           WRITE  (UNIT=9, FMT=*)   sum(agrid*DBN_bench(1,:,zi,nlambda/2+1, ne/2+1)) &
                                                        & /sum(DBN_bench(1,:,zi,nlambda/2+1, ne/2+1))
ENDDO           
close (unit=8)
close (unit=9)

OPEN (UNIT=8, FILE='CE_by_Asset_z_med_E_Lambda_age16', STATUS='replace') 
OPEN (UNIT=9, FILE='MeanAsset_by_z_med_E_Lambda_age16', STATUS='replace') 
DO zi=1,nz
           WRITE  (UNIT=8, FMT=*)   Cons_Eq_Welfare(16,:,zi,nlambda/2+1, ne/2+1)
           WRITE  (UNIT=9, FMT=*)   sum(agrid*DBN_bench(16,:,zi,nlambda/2+1, ne/2+1)) &
                                                        & /sum(DBN_bench(16,:,zi,nlambda/2+1, ne/2+1))
ENDDO           
close (unit=8)
close (unit=9)

OPEN (UNIT=8, FILE='CE_by_Asset_z_med_E_Lambda_age31', STATUS='replace') 
OPEN (UNIT=9, FILE='MeanAsset_by_z_med_E_Lambda_age31', STATUS='replace') 
DO zi=1,nz
           WRITE  (UNIT=8, FMT=*)   Cons_Eq_Welfare(31,:,zi,nlambda/2+1, ne/2+1)
           WRITE  (UNIT=9, FMT=*)   sum(agrid*DBN_bench(31,:,zi,nlambda/2+1, ne/2+1)) &
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

age_group_counter=1
tot_cons_by_agegroup_by_z_exp=0.0_DP
tot_hours_by_agegroup_by_z_exp =0.0_DP
tot_aprime_by_agegroup_by_z_exp =0.0_DP
tot_at_return_by_agegroup_by_z_exp =0.0_DP
tot_cap_tax_by_agegroup_by_z_exp=0.0_DP
tot_lab_tax_by_agegroup_by_z_exp=0.0_DP
tot_cons_tax_by_agegroup_by_z_exp=0.0_DP

DO age=1,MaxAge 

    DO while (age .gt.   age_limit(age_group_counter+1) )
        age_group_counter=age_group_counter+1
    ENDDO    
 
    DO ai=1,na
        DO zi=1,nz
            DO lambdai=1,nlambda
                DO ei=1,ne

         tot_cons_by_agegroup_by_z_exp(age_group_counter,zi) =      tot_cons_by_agegroup_by_z_exp(age_group_counter,zi)+ &
                    & Cons(age,ai,zi,lambdai,ei)*DBN1(age,ai,zi,lambdai,ei)

         tot_hours_by_agegroup_by_z_exp(age_group_counter,zi) =    tot_hours_by_agegroup_by_z_exp(age_group_counter,zi) + &
                    & HOURS(age,ai,zi,lambdai,ei)*DBN1(age,ai,zi,lambdai,ei)

        tot_aprime_by_agegroup_by_z_exp(age_group_counter,zi) =    tot_aprime_by_agegroup_by_z_exp(age_group_counter,zi) + &
                    & Aprime(age,ai,zi,lambdai,ei)*DBN1(age,ai,zi,lambdai,ei)

        tot_at_return_by_agegroup_by_z_exp(age_group_counter,zi) =    tot_at_return_by_agegroup_by_z_exp(age_group_counter,zi) + &
                    & MBGRID_exp(ai,zi)*agrid(ai)*DBN1(age,ai,zi,lambdai,ei)
  
        tot_lab_tax_by_agegroup_by_z_exp(age_group_counter,zi) =     tot_lab_tax_by_agegroup_by_z_exp(age_group_counter,zi) + &
                    & DBN1(age,ai,zi,lambdai,ei) * ( &
                     & tauL_exp * wage_exp * eff_un(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei) )
 
        tot_cons_tax_by_agegroup_by_z_exp(age_group_counter,zi) =     tot_cons_tax_by_agegroup_by_z_exp(age_group_counter,zi) + &
                    & DBN1(age,ai,zi,lambdai,ei) * ( &
                     & tauC *cons(age,ai,zi,lambdai,ei) )
         
            
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO

OPEN (UNIT=6, FILE='mean_cons_by_agegroup_by_z_exp.txt', STATUS='replace')  
OPEN (UNIT=7, FILE='mean_hours_by_agegroup_by_z_exp.txt', STATUS='replace')  
OPEN (UNIT=8, FILE='mean_aprime_by_agegroup_by_z_exp.txt', STATUS='replace')  
OPEN (UNIT=9, FILE='mean_at_weighted_return_by_agegroup_by_z_exp.txt', STATUS='replace')  
!OPEN (UNIT=10, FILE=' tot_cap_tax_by_agegroup_by_z_exp.txt', STATUS='replace')  
OPEN (UNIT=11, FILE=' tot_lab_tax_by_agegroup_by_z_exp.txt', STATUS='replace')  
OPEN (UNIT=12, FILE=' tot_cons_tax_by_agegroup_by_z_exp.txt', STATUS='replace')  
DO age_group_counter=1,max_age_category
    WRITE  (UNIT=6, FMT=*)  tot_cons_by_agegroup_by_z_exp(age_group_counter,:)/ size_by_agegroup_z_exp(age_group_counter,:)
    WRITE  (UNIT=7, FMT=*)  tot_hours_by_agegroup_by_z_exp(age_group_counter,:)/ size_by_agegroup_z_exp(age_group_counter,:)
    WRITE  (UNIT=8, FMT=*)  tot_aprime_by_agegroup_by_z_exp(age_group_counter,:)/ size_by_agegroup_z_exp(age_group_counter,:)
  WRITE(UNIT=9, FMT=*) tot_at_return_by_agegroup_by_z_exp(age_group_counter,:)/tot_wealth_by_agegroup_z_exp(age_group_counter,:)-1
!    WRITE  (UNIT=10, FMT=*) tot_cap_tax_by_agegroup_by_z_exp(age_group_counter,:)  
    WRITE  (UNIT=11, FMT=*) tot_lab_tax_by_agegroup_by_z_exp(age_group_counter,:) 
    WRITE  (UNIT=12, FMT=*) tot_cons_tax_by_agegroup_by_z_exp(age_group_counter,:) 
ENDDO
close (UNIT=6)
close (UNIT=7)
close (UNIT=8)
close (UNIT=9)
!close (UNIT=10)
close (UNIT=11)
close (UNIT=12)



print*,'---------------------------'
print*,''
print*,'tauK=', tauK
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
    DO zi=1,nz
        DO lambdai=1,nlambda          
              DO ei=1,ne            
          
                    CALL spline( agrid, ValueFunction(age+1, :, zi, lambdai, ei) , na , &
                    & gamma*MBGRID(1,zi) *Cons(age+1, 1, zi, lambdai, ei) **((1.0_DP-sigma)*gamma-1.0_DP)/(1_DP+tauC), &
                    & gamma*MBGRID(na,zi)*Cons(age+1, na, zi, lambdai, ei)**((1.0_DP-sigma)*gamma-1.0_DP)/(1_DP+tauC), ValueP2)  
                  
                    DO ai=1,na    
                         call splint( agrid, ValueFunction(age+1, :, zi, lambdai, ei), &
                                & ValueP2, na, Aprime(age,ai,zi,lambdai, ei), ValueP(ai))   
    
                         ValueFunction(age, ai, zi, lambdai, ei) = ((Cons(age,ai,zi,lambdai,ei)**gamma) &
                           & * (1.0_DP-Hours(age,ai,zi,lambdai,ei))**(1.0_DP-gamma))**(1.0_DP-sigma)/(1.0_DP-sigma) &
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
                    & (gamma*MBGRID(1,zi)/(1.0_DP+tauC)) * sum(pr_e(ei,:)* &
                    & Cons(age+1, 1, zi, lambdai, :)**((1.0_DP-sigma)*gamma-1.0_DP) * &
                    & (1.0_DP-Hours(age+1,1,zi,lambdai,:))**((1.0_DP-gamma)*(1.0_DP-sigma))),&                    
                    & (gamma*MBGRID(na,zi)/(1.0_DP+tauC)) * sum(pr_e(ei,:)* &
                    & Cons(age+1, na, zi, lambdai, :)**((1.0_DP-sigma)*gamma-1.0_DP) * &
                    & (1.0_DP-Hours(age+1,na,zi,lambdai,:))**((1.0_DP-gamma)*(1.0_DP-sigma))),&
                    & ValueP2)  
                    
                    DO ai=1,na                        
                         call splint( agrid, ExpValueP, ValueP2, na, Aprime(age,ai,zi,lambdai, ei), ValueP(ai))   
                         ValueFunction(age, ai, zi, lambdai, ei) = ((Cons(age,ai,zi,lambdai,ei)**gamma) &
                           & * (1.0_DP-Hours(age,ai,zi,lambdai,ei))**(1.0_DP-gamma))**(1.0_DP-sigma)/(1.0_DP-sigma) &
                           & + beta*survP(age)* ValueP(ai)
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
real(DP) ::  GBAR_NL, BT_EARNINGS , new_tauL

GBAR=0.0_DP
GBAR_K=0.0_DP
GBAR_L=0.0_DP
GBAR_C=0.0_DP

DO age=1, MaxAge
DO ai=1,na
DO zi=1,nz
DO lambdai=1,nlambda
DO ei=1,ne

   if (kstar(zi) .le. (vartheta*agrid(ai)) ) then
      GBAR = GBAR + DBN1(age,ai,zi,lambdai,ei) * ( tauK * int_rate * agrid(ai)    &
          & + tauW * (1.0_DP+int_rate) * agrid(ai) - tauK*tauW*int_rate*agrid(ai)   &
          & + tauL * wage * eff_un(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei) &
          & + tauC * cons(age, ai, zi, lambdai,ei)    )

      GBAR_K = GBAR_K + DBN1(age,ai,zi,lambdai,ei) * ( tauK * int_rate * agrid(ai)    &
          & + tauW * (1.0_DP+int_rate) * agrid(ai) - tauK*tauW*int_rate*agrid(ai) )             
   else 
      GBAR = GBAR + DBN1(age,ai,zi,lambdai,ei) * ( tauK * int_rate * agrid(ai)    &
          & + tauW * (1.0_DP+int_rate) * agrid(ai) - tauK*tauW*int_rate*agrid(ai)   &
          & + (rr*(vartheta*zgrid(zi)*agrid(ai))**mu-(int_rate+DepRate)*vartheta*agrid(ai))*(tauK+tauW-tauK*tauW) &
          & + tauL * wage * eff_un(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei) &
          & + tauC * cons(age, ai, zi, lambdai,ei)    )

      GBAR_K = GBAR_K + DBN1(age,ai,zi,lambdai,ei) * ( tauK * int_rate * agrid(ai)    &
          & + tauW * (1.0_DP+int_rate) * agrid(ai) - tauK*tauW*int_rate*agrid(ai)   &
          & + (rr*(vartheta*zgrid(zi)*agrid(ai))**mu-(int_rate+DepRate)*vartheta*agrid(ai))*(tauK+tauW-tauK*tauW)  )  
   endif
          
   GBAR_L = GBAR_L + DBN1(age,ai,zi,lambdai,ei) * ( &
          & tauL * wage * eff_un(age,lambdai,ei)*Hours(age,ai,zi,lambdai,ei) )
          
   GBAR_C = GBAR_C + DBN1(age,ai,zi,lambdai,ei) * (tauC * cons(age, ai, zi, lambdai,ei)    )
                 
ENDDO
ENDDO
ENDDO
ENDDO
ENDDO


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


print*, 'GBAR_bench',GBAR_bench, 'GBAR=',GBAR

END  SUBROUTINE GOVNT_BUDGET

!====================================================================

SUBROUTINE FIND_DBN_EQ
USE PARAMETERS
USE GLOBAL
use  programfunctions
IMPLICIT NONE
INTEGER:: tklo, tkhi, age1, age2, z1, z2, a1, a2, lambda1, lambda2, e1, e2, DBN_iter, simutime,  iter_indx
REAL :: DBN_dist, DBN_criteria , DIST_INT_RATE, new_int_rate, old_int_rate
REAL(DP), DIMENSION(MaxAge, na, nz, nlambda, ne) :: PrAprimelo, PrAprimehi, DBN2
INTEGER,  DIMENSION(MaxAge, na, nz, nlambda, ne) :: Aplo, Aphi

DBN_criteria = 1.0E-8_DP

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
simutime =1
DO WHILE ( ( DBN_dist .ge. DBN_criteria ) .and. ( simutime .le. MaxSimuTime ) )
!    print*, 'sum DBN1=', sum(DBN1)
    DBN2=0.0_DP
    
! Everyone in MaxAge dies. Those who die, switch to z2, lambda2 and start at ne/2+1
    age1=MaxAge
    DO z1=1,nz
    DO a1=1,na
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
             QBAR= QBAR+  DBN1(age1, a1, z1, lambda1, e1) * ( zgrid(z1) * kdemand(a1,z1) )**mu
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

!        print*,'DBN_dist=',DBN_dist, 'QBAR=', QBAR ,  'NBAR=', NBAR 

        MeanWealth=0.0_DP
        DO a1=1,na
               MeanWealth = MeanWealth + sum(DBN1(:, a1, :, :, :))*agrid(a1)
        ENDDO
        
!         old_int_rate = int_rate
!         DIST_INT_RATE=1.0_DP
!         DO WHILE (DIST_INT_RATE .gt. 0.000001)
        
!           BBAR =0.0_DP          
!           DO z1=1,nz
!           DO a1=1,na
!             BBAR = BBAR + sum(DBN1(:, a1, z1, :, :))* &
!               & ( min((mu*rr*zgrid(z1)**mu/(int_rate+DepRate))**(1.0_DP/(1.0_DP-mu)),vartheta*agrid(a1)) &
!               & - agrid(a1)) 
!           ENDDO
!           ENDDO
          
!           new_int_rate = int_rate * ( 1.0_DP + 0.1*BBAR / MeanWealth )
!           DIST_INT_RATE = abs(BBAR/MeanWealth)
!           int_rate  = new_int_rate
! !          print*,'BBAR/MeanWealth=',BBAR/MeanWealth, 'int_rate=',int_rate,'new_int_rate=',new_int_rate !,'dist=',DIST_INT_RATE
          
!         ENDDO  
!         int_rate = 0.1_DP * min(int_rate,0.5_DP) + 0.9_DP * old_int_rate
        
!         print*,'DBN_dist=',DBN_dist, 'BBAR/MeanWealth=',BBAR/MeanWealth, 'int_rate=',int_rate
!           BBAR =0.0_DP          
!           DO z1=1,nz
!           DO a1=1,na
!             BBAR = BBAR + sum(DBN1(:, a1, z1, :, :))* &
!               & ( min((mu*rr*zgrid(z1)**mu/(int_rate+DepRate))**(1.0_DP/(1.0_DP-mu)),vartheta*agrid(a1)) &
!               & - agrid(a1)) 
!           ENDDO
!           ENDDO
!         print*,'New BBAR/MeanWealth=',BBAR/MeanWealth

        int_rate = zbrent(Agg_Debt,0.0_dp,0.50_dp,brent_tol) 
          BBAR =0.0_DP          
          DO z1=1,nz
          DO a1=1,na
            BBAR = BBAR + sum(DBN1(:, a1, z1, :, :))* &
              & ( min((mu*rr*zgrid(z1)**mu/(int_rate+DepRate))**(1.0_DP/(1.0_DP-mu)),vartheta*agrid(a1)) &
              & - agrid(a1)) 
          ENDDO
          ENDDO
        print*,'DBN_dist=',DBN_dist, 'BBAR/MeanWealth=',BBAR/MeanWealth, Agg_Debt(int_rate), 'int_rate=',int_rate

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
 
ENDDO ! WHILE

BBAR_POS = 0.0_DP
DO z1=1,nz
DO a1=1,na
    BBAR_POS = BBAR_POS + sum(DBN1(:, a1, z1, :, :))* &
      & abs( min((mu*rr*zgrid(z1)**mu/(int_rate+DepRate))**(1.0_DP/(1.0_DP-mu)),vartheta*agrid(a1)) &
      & - agrid(a1)) 
ENDDO
ENDDO
External_Debt_GDP = 0.5*BBAR_POS / YBAR


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
REAL(DP):: MeanATReturn, StdATReturn, VarATReturn , VarReturn 
REAL(DP), DIMENSION(nz):: MeanATReturn_by_z, MeanReturn_by_z, size_by_z

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
prct20_wealth =  1.0_DP-cdf_tot_a_by_prctile(80)/cdf_tot_a_by_prctile(100)
prct40_wealth =  1.0_DP-cdf_tot_a_by_prctile(60)/cdf_tot_a_by_prctile(100)

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
                         & *  log(  yh(age, lambdai, ei) * Hours(age, ai, zi, lambdai, ei) )
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
                         & * ( log(  yh(age, lambdai, ei) * Hours(age, ai, zi, lambdai, ei) ) &
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
MeanATReturn = 0.0_DP
MeanReturn = 0.0_DP
MeanCons  = 0.0_DP

DO age=1,MaxAge
DO zi=1,nz
DO ai=1,na
DO lambdai=1,nlambda
DO ei=1, ne
   MeanWealth = MeanWealth  +   DBN1(age, ai, zi, lambdai, ei) * agrid(ai)         
   MeanCons  =  MeanCons  + DBN1(age, ai, zi, lambdai, ei) * cons(age, ai, zi, lambdai, ei)

   if (kstar(zi) .le. (vartheta*agrid(ai)) ) then
      MeanReturn = MeanReturn+ DBN1(age, ai, zi, lambdai, ei) * int_rate * agrid(ai)
      MeanATReturn = MeanATReturn  +   DBN1(age, ai, zi, lambdai, ei) &
            & *  ((1.0_DP+int_rate*(1.0_DP-tauK))*(1.0_DP-tauW)-1.0_DP)* agrid(ai) 
   else
      MeanReturn = MeanReturn+ DBN1(age, ai, zi, lambdai, ei) * agrid(ai) * (int_rate &
         & + (rr*mu*((vartheta*zgrid(zi))**mu)*agrid(ai)**(mu-1.0_DP)-(int_rate+DepRate)*vartheta)) 
      MeanATReturn = MeanATReturn  +   DBN1(age, ai, zi, lambdai, ei) &
         & * ((1.0_DP+int_rate*(1.0_DP-tauK))*(1.0_DP-tauW) &
         & + (rr*mu*((vartheta*zgrid(zi))**mu)*agrid(ai)**(mu-1.0_DP) &
         & -(int_rate+DepRate)*vartheta)*(1.0_DP-tauK)*(1.0_DP-tauW) -1.0_DP)*agrid(ai) 
   endif

ENDDO
ENDDO
ENDDO    
ENDDO    
ENDDO    
Wealth_Output = MeanWealth/YBAR 
MeanReturn = MeanReturn/MeanWealth

Bequest_Wealth=0.0_DP
DO zi=1,nz
DO ai=1,na
DO lambdai=1,nlambda
DO ei=1, ne
     Bequest_Wealth = Bequest_Wealth  +   DBN1(1, ai, zi, lambdai, ei) * agrid(ai)         
ENDDO
ENDDO
ENDDO    
ENDDO  
Bequest_Wealth =Bequest_Wealth/MeanWealth

SSE_Moments = (1.0-Wealth_Output/3.0_DP)**2.0_DP  + (1.0_DP-prct1_wealth/0.34_DP)**2.0_DP  + (prct10_wealth-0.71_DP)**2.0_DP &
                   & + (1.0_DP-Std_Log_Earnings_25_60 / 0.8_DP)**2.0_DP + (1.0_DP-meanhours_25_60/0.4_DP)**2.0_DP &
                   & + (1.0_DP-MeanReturn/0.069_DP)**2.0_DP
                   
print*,'Moments=',Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, MeanReturn
print*,''
print*,'--------------------------------------------------------------------------------------'


END SUBROUTINE COMPUTE_STATS

!====================================================================

SUBROUTINE WRITE_VARIABLES(bench_indx)
USE GLOBAL
use  programfunctions
IMPLICIT NONE
integer :: bench_indx,  prctile


IF (bench_indx .gt. 0) then 
OPEN   (UNIT=2, FILE='output_bench.txt', STATUS='replace')

WRITE  (UNIT=2, FMT=*)  'Params=[', params,']'
WRITE  (UNIT=2, FMT=*)  'tauK=',tauK,'tauW=',tauW,'tauL=',tauL,'tauC=',tauC
WRITE  (UNIT=2, FMT=*)  'MeanWealth_bench=',MeanWealth
WRITE  (UNIT=2, FMT=*)  'QBAR_bench=',QBAR_bench
WRITE  (UNIT=2, FMT=*)  'NBAR_bench=',NBAR_bench
WRITE  (UNIT=2, FMT=*)  'EBAR_bench=',EBAR_bench
WRITE  (UNIT=2, FMT=*)  'Pr_bench=',rr_bench
WRITE  (UNIT=2, FMT=*)  'wage_bench=',wage_bench
WRITE  (UNIT=2, FMT=*)  'Y_bench=',Y_bench
WRITE  (UNIT=2, FMT=*)  'Cons_bench=',meancons
WRITE  (UNIT=2, FMT=*)  'External_Debt_GDP_bench',External_Debt_GDP
WRITE  (UNIT=2, FMT=*)  'Int_rate_bench=',int_rate
WRITE  (UNIT=2, FMT=*)  ''
WRITE  (UNIT=2, FMT=*)  'GBAR_bench=', GBAR_bench
WRITE  (UNIT=2, FMT=*)  'SSC_Payments_bench=', SSC_Payments
WRITE  (UNIT=2, FMT=*)  'GBAR_K_bench=', GBAR_K
WRITE  (UNIT=2, FMT=*)  'GBAR_L_bench=', GBAR_L
WRITE  (UNIT=2, FMT=*)  'GBAR_C_bench=', GBAR_C
WRITE  (UNIT=2, FMT=*)  ''
WRITE  (UNIT=2, FMT=*)  'MOMENTS'
WRITE  (UNIT=2, FMT=*)  'Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, MeanReturn'
WRITE  (UNIT=2, FMT=*)  Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, MeanReturn
WRITE  (UNIT=2, FMT=*)  ''
WRITE  (UNIT=2, FMT=*)  'prct20_wealth=',prct20_wealth,'prct40_wealth=',prct40_wealth
WRITE  (UNIT=2, FMT=*)  ''
!WRITE  (UNIT=2, FMT=*)  'StdReturn=',StdReturn
WRITE  (UNIT=2, FMT=*)  'Bequest_Wealth=',Bequest_Wealth


else

OPEN   (UNIT=2, FILE='output_exp.txt', STATUS='replace')

WRITE  (UNIT=2, FMT=*)  'tauK=',tauK,'tauW=',tauW,'tauL=',tauL,'tauC=',tauC
WRITE  (UNIT=2, FMT=*)  'MeanWealth_exp=',MeanWealth
WRITE  (UNIT=2, FMT=*)  'QBAR_exp=',QBAR_exp
WRITE  (UNIT=2, FMT=*)  'NBAR_exp=',NBAR_exp
WRITE  (UNIT=2, FMT=*)  'EBAR_exp=',EBAR_exp
WRITE  (UNIT=2, FMT=*)  'Pr_exp=',rr_exp
WRITE  (UNIT=2, FMT=*)  'wage_exp=',wage_exp
WRITE  (UNIT=2, FMT=*)  'Y_exp=',Y_exp
WRITE  (UNIT=2, FMT=*)  'Cons_exp=',meancons
WRITE  (UNIT=2, FMT=*)  'External_Debt_GDP_exp',External_Debt_GDP
WRITE  (UNIT=2, FMT=*)  'Int_rate_exp=',int_rate
WRITE  (UNIT=2, FMT=*)  ''
WRITE  (UNIT=2, FMT=*)  'Y_exp/Y_bench-1=',Y_exp/Y_bench-1.0_DP
WRITE  (UNIT=2, FMT=*)  ''
WRITE  (UNIT=2, FMT=*)  'GBAR_exp=', GBAR_exp
WRITE  (UNIT=2, FMT=*)  'SSC_Payments_exp=', SSC_Payments
WRITE  (UNIT=2, FMT=*)  'GBAR_K_exp=', GBAR_K
WRITE  (UNIT=2, FMT=*)  'GBAR_L_exp=', GBAR_L
WRITE  (UNIT=2, FMT=*)  'GBAR_C_exp=', GBAR_C
WRITE  (UNIT=2, FMT=*)  ''
WRITE  (UNIT=2, FMT=*)  'MOMENTS'
WRITE  (UNIT=2, FMT=*)  'Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, MeanReturn'
WRITE  (UNIT=2, FMT=*)  Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60, MeanReturn
WRITE  (UNIT=2, FMT=*)  ''
WRITE  (UNIT=2, FMT=*)  'prct20_wealth=',prct20_wealth,'prct40_wealth=',prct40_wealth
WRITE  (UNIT=2, FMT=*)  ''
!WRITE  (UNIT=2, FMT=*)  'StdReturn=',StdReturn
WRITE  (UNIT=2, FMT=*)  'Bequest_Wealth=',Bequest_Wealth


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
!Hours(1:RetAge-1       , :, :, :,: ) =0.4_DP

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
                    EndoCons(ai) =    Cons(age+1, ai, zi, lambdai,ei) * & 
                                 &  ( beta*survP(age)*MBGRID(ai,zi) ) ** ( 1.0_DP / (gamma*(1.0_DP-sigma)-1.0_DP))     
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
                   EndoCons(ai) = ((gamma*yh(age, lambdai,ei)/(1.0_DP-gamma))**((1.0_DP-gamma)*(1.0_DP-sigma)) &
                        & *  beta*survP(age)*MBGRID(ai,zi)  &
                        & *  sum( pr_e(ei,:) * (Cons(age+1,ai,zi,lambdai,:)**(gamma*(1.0_DP-sigma)-1.0_DP)) &
                        & *  ( (1.0_DP-Hours(age+1, ai, zi, lambdai,:))**((1.0_DP-gamma)*(1.0_DP-sigma)))))**(-1.0_DP/sigma)                    
                     
                   EndoHours(ai) = 1.0_DP - (1.0_DP-gamma)*EndoCons(ai)/(gamma*yh(age, lambdai,ei))   
                      
                   If (EndoHours(ai) .lt. 0.0_DP) then
                    EndoHours(ai) = 0.0_DP 
                    EndoCons(ai)  = ( beta*survP(age)*MBGRID(ai,zi)  &
                    & * sum( pr_e(ei,:) * (Cons(age+1,ai,zi,lambdai,:)**(gamma*(1.0_DP-sigma)-1.0_DP)) &
                    & *((1.0_DP-Hours(age+1,ai,zi,lambdai,:))**((1.0_DP-gamma)*(1.0_DP-sigma)))) )&
                    & **(1.0_DP/(gamma*(1.0_DP-sigma)-1.0_DP))                                                 
                   endif                 
                   EndoYgrid(ai) =  agrid(ai) +  EndoCons(ai) - yh(age, lambdai,ei)* EndoHours(ai)                   
                ENDDO ! ai  

                tempai=1           
                DO WHILE ( YGRID(tempai,zi) .lt. EndoYgrid(1) )
                      tempai = tempai + 1
                ENDDO
                
                ! Find  decision rules on exogenous grids
                DO ai=tempai,na               
                         Cons(age, ai, zi, lambdai,ei)= Linear_Int(EndoYgrid,&
                                    & EndoCons,na, YGRID(ai,zi))                          
                         Hours(age, ai, zi, lambdai,ei) = max(0.0_DP, &                          
                                    & 1.0_DP - (1.0_DP-gamma)*Cons(age,ai,zi,lambdai,ei)/(gamma*yh(age, lambdai,ei)) )                                
                         Aprime(age, ai, zi, lambdai,ei) = YGRID(ai,zi)  &
                                        & + yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei)  & 
                                        & - Cons(age, ai, zi, lambdai,ei)   
                         If   (Aprime(age, ai, zi, lambdai,ei)  .lt. amin) then
                               Aprime(age, ai, zi, lambdai,ei) = amin
                               
                               Hours(age, ai, zi, lambdai,ei)  = max( gamma - &
                                & (1.0_DP-gamma) *( YGRID(ai,zi) - Aprime(age, ai, zi, lambdai,ei)) / yh(age, lambdai,ei) ,0.0_DP)       
                              
                                Cons(age, ai, zi, lambdai,ei)= YGRID(ai,zi)  + yh(age, lambdai,ei)*Hours(age, ai, zi, lambdai,ei)&
                                                & -Aprime(age, ai, zi, lambdai,ei)        
                                IF (Cons(age, ai, zi, lambdai,ei) .le. 0.0_DP)  THEN
                                    print*,'w1: Cons(age, ai, zi, lambdai,ei)=',Cons(age, ai, zi, lambdai,ei)
                                ENDIF                   
                         endif      
                 ENDDO ! ai                  

                 ai=1           
                 DO WHILE ( YGRID(ai,zi) .lt. EndoYgrid(1) )
!                       ap_temp    = Aprime(age, tempai, zi, lambdai,ei) 
                       brentvalue = brent(min(amin,YGRID(ai,zi)), (amin+YGRID(ai,zi))/2.0_DP , YGRID(ai,zi)+yh(age, lambdai,ei),&
                             & FOC_W, brent_tol, Aprime(age, ai, zi, lambdai,ei) )
                        
                       Hours(age, ai, zi, lambdai,ei)  =  max(0.0_DP, &
                                & gamma - (1.0_DP-gamma)*(YGRID(ai,zi) - Aprime(age, ai, zi, lambdai,ei))/yh(age, lambdai,ei) )            
 
                       Cons(age, ai, zi, lambdai,ei)=  YGRID(ai,zi)  + yh(age, lambdai,ei) * Hours(age, ai, zi, lambdai,ei) &
                                        & -Aprime(age, ai, zi, lambdai,ei)
                                        
                       IF (Cons(age, ai, zi, lambdai,ei) .le. 0.0_DP)  THEN
                            print*,'w2:',age,zi,lambdai,ei,ai, 'Cons(age, ai, zi, lambdai,ei)=',Cons(age, ai, zi, lambdai,ei), &
                                & 'Aprime(age, ai, zi, lambdai,ei)=',Aprime(age, ai, zi, lambdai,ei), &
                                & 'Hours(age, ai, zi, lambdai,ei)=', Hours(age, ai, zi, lambdai,ei), &
                                & 'yh(age, lambdai,ei)=', yh(age, lambdai,ei),'YGRID(ai,zi)=',YGRID(ai,zi)
                            !pause
                        ENDIF                   
                        ai = ai + 1
              ENDDO  
                 
            ENDDO  ! ei         
         ENDDO ! zi
    ENDDO ! lambdai   
ENDDO !age


cons = cons/(1.0_DP+tauC)

END SUBROUTINE




!====================================================================
! THIS YGRID and Marginal Benefit of Investment GRID 
! NEEDS TO BE COMPUTED FOR EACH TIME THE INTEREST RATE "rr" IS UPDATED. 

SUBROUTINE FORM_Y_MB_GRID(TYGRID, TMBGRID)
USE GLOBAL
IMPLICIT NONE
REAL(DP), DIMENSION(na,nz), INTENT(OUT)  :: TYGRID, TMBGRID

DO zi=1,nz
        kstar(zi)  = (mu*rr*zgrid(zi)**mu/(int_rate+DepRate))**(1.0_DP/(1.0_DP-mu))
        
        pistar(zi) = (1.0_DP-mu)*(mu*rr*zgrid(zi)**mu)**(1.0_DP/(1.0_DP-mu)) &
                &  / ( mu* (int_rate+DepRate)**(mu/(1.0_DP-mu))  ) 
                
!       print*,'kstar=',kstar(zi),'pistar=',pistar(zi), rr*(zgrid(zi)*kstar(zi))**mu-(int_rate+DepRate)*kstar(zi)       
ENDDO
!pause

DO ai=1,na
DO zi=1,nz   
   
   kdemand(ai,zi) = min( kstar(zi), vartheta*agrid(ai) )
   bdemand(ai,zi) = kdemand(ai,zi) - agrid(ai)
   
   if (kstar(zi) .le. (vartheta*agrid(ai)) ) then
      TYGRID(ai,zi) = (1.0_DP+int_rate*(1.0_DP-tauK))*(1.0_DP-tauW)*agrid(ai) +pistar(zi)*(1.0_DP-tauK)*(1.0_DP-tauW) 
      TMBGRID(ai,zi)= (1.0_DP+int_rate*(1.0_DP-tauK))*(1.0_DP-tauW) 
   else
      TYGRID(ai,zi) = (1.0_DP+int_rate*(1.0_DP-tauK))*(1.0_DP-tauW)*agrid(ai) &
         & + (rr*(vartheta*zgrid(zi)*agrid(ai))**mu-(int_rate+DepRate)*vartheta*agrid(ai))*(1.0_DP-tauK)*(1.0_DP-tauW)
      TMBGRID(ai,zi)= (1.0_DP+int_rate*(1.0_DP-tauK))*(1.0_DP-tauW) &
         & + (rr*mu*((vartheta*zgrid(zi))**mu)*agrid(ai)**(mu-1.0_DP)-(int_rate+DepRate)*vartheta)*(1.0_DP-tauK)*(1.0_DP-tauW)
   endif
ENDDO
ENDDO

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
yh = (1.0_DP-tauL) * Waget * eff_un

! This part computes Retirement Income
RetY_lambda_e = phi_lambda_e  * Ebart 
IF ((KeepSSatBench .eq. 1) .AND. (solving_bench .eq. 0)) THEN
    RetY_lambda_e = phi_lambda_e  * Ebar_bench
ENDIF

!print*,'yh=',yh
!print*, 'RetY =',RetY_lambda_e
!pause

END SUBROUTINE 

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
REAL(DP) :: mean_panel_lifetime_eff_unit , lambdaBAR
INTEGER  :: cohortsize
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
! print*,'agrid=',agrid
!pause

OPEN   (UNIT=2, FILE='agrid', STATUS='replace')
WRITE  (UNIT=2, FMT=*)  agrid
CLOSE (UNIT=2)

!pause


  m=(amax-amin)**(1.0_DP/a_theta)/REAL(fine_na-1,DP)
  DO ai=1,fine_na
    fine_agrid(ai)=REAL(ai-1,DP)*m
  END DO
  fine_agrid=amin+fine_agrid**a_theta
  
CALL tauchen(mtauchen, 0.0_DP, rho_E,sigma_e_eps,ne,egrid,pr_e,Ge)
CALL tauchen(mtauchen, 0.0_DP, rho_lambda,sigma_lambda_eps,nlambda,lambdagrid,pr_lambda,Glambda)
CALL tauchen(mtauchen, 0.0_DP, rho_z,sigma_z_eps,nz,zgrid,pr_z,Gz)

zgrid=exp(zgrid)
zgrid=zgrid+mu_z

egrid=exp(egrid)
lambdagrid=exp(lambdagrid) 

OPEN   (UNIT=2, FILE='zgrid', STATUS='replace')
WRITE  (UNIT=2, FMT=*)  zgrid
CLOSE (UNIT=2)

OPEN   (UNIT=2, FILE='egrid', STATUS='replace')
WRITE  (UNIT=2, FMT=*)  egrid
CLOSE (UNIT=2)

OPEN   (UNIT=2, FILE='lambdagrid', STATUS='replace')
WRITE  (UNIT=2, FMT=*) lambdagrid
CLOSE (UNIT=2)


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


! Daphne's retirement income calculator
!     
!phi_lambda(1) =  0.251822117832046_DP
!phi_lambda(2) =  0.431345243929363_DP
!phi_lambda(3) =  0.560744163108887_DP
!phi_lambda(4) =  0.793921042376030_DP 
!phi_lambda(5) =  1.25254963402369_DP  
! 

!print*,'retirement income corrected'
!print*,''
!print*,'phi lambda=',phi_lambda 
!print*,'lambdabar=',lambdaBAR
!print*,'lambdagrid=',lambdagrid
!pause
     
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

print*,'kappagrid='
print*,kappagrid 
OPEN   (UNIT=2, FILE='kappagrid', STATUS='replace')
WRITE  (UNIT=2, FMT=*)  kappagrid
CLOSE (UNIT=2)
!PAUSE

  !----------------------------------------------
  ! life-cycle component
  !----------------------------------------------
    
  ! Population Numbers from Bell and Miller (2002)
  
  pop(1)= 197316.0_DP
  pop(2)= 197141.0_DP
  pop(3)= 196959.0_DP
  pop(4)= 196770.0_DP
  pop(5)= 196580.0_DP
  pop(6)= 196392.0_DP
  pop(7)= 196205.0_DP
  pop(8)= 196019.0_DP
  pop(9)= 195830.0_DP
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

if (MaxAge .gt. 80)  then 
  pop(MaxAge)=2454.0_DP
endif 
  

!pop=1.0_DP
  ! Survival probabilities: surv(i)=prob(alive in i+1|alive in i)
  FORALL (age=1:maxAge-1) survP(age)= pop(age+1)/pop(age)
  survP(maxAge)=0.0_DP
    
 
! Set the initial distribution

print*,'agrid=',agrid(1:20)
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


!====================================================================

SUBROUTINE  SIMULATION
use parameters
USE GLOBAL
use  programfunctions
IMPLICIT NONE
integer::  currentzi, currentlambdai, currentei
REAL(DP)::  tempnoage, tempnoz, tempnolambda, tempnoe, tempno, currenta, currentY
REAL(DP) ::  start_timet, finish_timet
INTEGER :: agecounter, agesign, tage, tzi, tlambdai, tei, tklo, tkhi, paneli, simutime
INTEGER, DIMENSION(MaxAge) :: requirednumberby_age, cdfrequirednumberby_age
INTEGER, DIMENSION(totpop)::    panelage , panelz , panellambda, panele,   newpanelage , newpanelz , newpanellambda, newpanele
REAL(DP), DIMENSION(totpop)::   panela,  newpanela,  panel_return, panelcons, panelhours, panelaprime, panel_at_return


    age=1
    requirednumberby_age(age) = NINT(totpop*pop(age)/sum(pop))
    cdfrequirednumberby_age(age) = requirednumberby_age(age)
    DO age=2,MaxAge
        requirednumberby_age(age) = NINT(totpop*pop(age)/sum(pop))
        cdfrequirednumberby_age(age) = requirednumberby_age(age) + cdfrequirednumberby_age(age-1)
    ENDDO
    ! If the total number of people are not equal to the total population, then I will add the remainder to the last age
    requirednumberby_age(MaxAge) =  requirednumberby_age(MaxAge)-cdfrequirednumberby_age(MaxAge) + totpop
    cdfrequirednumberby_age(MaxAge) = totpop

!=====================================================================
!                     GENERATE   INITIAL   PANEL
!=====================================================================


newiseed=-1

!numberby_age_z_lambda=0
!numberby_age_e =0

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

! CORRECT THE NUMBER OF PEOPLE IF THERE ARE EXCESS
!
!   if (age .gt. 1) then
!        if ( (cdfrequirednumberby_age(age)-tempnoage*totpop) .gt.  (tempnoage*totpop-cdfrequirednumberby_age(age-1)) ) then
!            agesign=1
!            else
!                 agesign=-1
!        endif      
!    else
!        agesign=1
!    endif 
!   agecounter=1
!   tage=age        
!111 IF (sum(numberby_age_z_lambda(age,:,:)) .ge. sum(requirednumberby_age_z_lambda(age,:,:))) then
!            age = tage + agecounter * agesign
!            age = max(age,1)
!            age = min(age,MaxAge)
!            if (sum(numberby_age_z_lambda(age,:,:)) .ge. sum(requirednumberby_age_z_lambda(age,:,:))) then
!                age = age - agecounter * agesign
!                age = max(age,1)
!                age = min(age,MaxAge)
!                if (sum(numberby_age_z_lambda(age,:,:)) .ge. sum(requirednumberby_age_z_lambda(age,:,:))) then
!                    agecounter = agecounter +1
!                    go to 111
!                endif    
!            endif
!       ENDIF
!   
!   if (zi .gt. 1) then 
!       if ( (cdf_Gz(zi) -tempnoz) .gt.  (tempnoz-cdf_Gz(zi-1)) )    then
!           agesign=1
!           else
!                agesign=-1
!       endif      
!   else
!       agesign=1
!   endif
!   agecounter=1  
!   tzi=zi       
!112 IF (sum(numberby_age_z_lambda(age,zi,:)) .ge. sum(requirednumberby_age_z_lambda(age,zi,:))) then
!           zi = tzi + agecounter * agesign
!           zi = max(zi,1)
!           zi=min(zi,nz)
!           IF (sum(numberby_age_z_lambda(age,zi,:)) .ge. sum(requirednumberby_age_z_lambda(age,zi,:))) then
!               zi = zi - agecounter * agesign
!               zi = max(zi,1)
!               zi=min(zi,nz)
!               IF (sum(numberby_age_z_lambda(age,zi,:)) .ge. sum(requirednumberby_age_z_lambda(age,zi,:))) then
!                   agecounter = agecounter +1
!                   go to 112
!               ENDIF
!           ENDIF               
!       ENDIF    
! 
!   if (lambdai .gt. 1) then      
!       if ( (cdf_Glambda(lambdai) -tempnolambda) .gt.  (tempnolambda-cdf_Glambda(lambdai-1)) )    then
!           agesign=1
!           else
!                agesign=-1
!       endif  
!   else
!       agesign=1
!   endif 
!    
!   agecounter=1  
!   tlambdai=lambdai  
!113 IF  (numberby_age_z_lambda(age,zi,lambdai) .ge. requirednumberby_age_z_lambda(age,zi,lambdai)) then
!           lambdai = tlambdai + agecounter * agesign
!           lambdai = max(lambdai,1)
!           lambdai=min(lambdai,nlambda)
!           IF  (numberby_age_z_lambda(age,zi,lambdai) .ge. requirednumberby_age_z_lambda(age,zi,lambdai)) then
!               lambdai = lambdai - agecounter * agesign
!               lambdai = max(lambdai,1)
!               lambdai=min(lambdai,nlambda)
!               IF  (numberby_age_z_lambda(age,zi,lambdai) .ge. requirednumberby_age_z_lambda(age,zi,lambdai)) then
!                   agecounter = agecounter +1
!                   go to 113
!               ENDIF
!           ENDIF               
!       ENDIF
!
!   if (ei .gt. 1) then
!       if ( (Ge_byage(age,ei) -tempnoe) .gt.  (tempnolambda-Ge_byage(age,ei-1) ) )    then
!           agesign=1
!           else
!                agesign=-1
!       endif     
!    else
!        agesign=1
!    endif 
!      
!   agecounter=1  
!   tei=ei      
!114  IF (numberby_age_e(age,ei) .ge. requirednumberby_age_e(age,ei)) THEN
!           ei = tei + agecounter * agesign
!           ei = max(ei,1)
!           ei=min(ei,ne)
!           IF (numberby_age_e(age,ei) .ge. requirednumberby_age_e(age,ei)) THEN
!               ei = tei -  agecounter * agesign
!               ei = max(ei,1)
!               ei=min(ei,ne)
!               IF (numberby_age_e(age,ei) .ge. requirednumberby_age_e(age,ei)) THEN
!                   agecounter = agecounter +1
!                   go to 114
!              ENDIF
!           ENDIF        
!        ENDIF
!   numberby_age_e(age,ei) = numberby_age_e(age,ei)+1    
!   numberby_age_z_lambda(age,zi,lambdai) = numberby_age_z_lambda(age,zi,lambdai)+1
!    
! CORRECTION ENDED

   panelage(paneli)=age
   panelz(paneli)=zi
   panellambda(paneli)=lambdai
   panele(paneli)=ei
   
ENDDO

!print '("INITIAL = ",f6.3," seconds.")',finish_timet-start_timet

newpanelage = panelage
newpanelz = panelz
newpanele = panele
newpanellambda = panellambda

! SET INITIAL ASSET DISTRIBUTION
panela            = 1.0_DP
newpanela      = 1.0_DP
!=============================================================================
!
! SIMULATE FROM THE SECOND PERIOD SHOCKS AND UPDATE NEW DISTRIBUTIONS
!
!=============================================================================

!call cpu_time(start_timet) 

DO simutime=1, MaxSimuTime

panelage = newpanelage
panelz = newpanelz
panele = newpanele
panellambda = newpanellambda
panela = newpanela

!print*,'simutime=',simutime

!numberby_age=0
!deathby_age=0
!survivingby_age =0
!
!deathby_age_z_lambda = 0
!survivingby_age_z_lambda=0
!numberby_age_z_lambda=0
!numberby_e_e=0

newpanela=amin

DO paneli=1,totpop
    
       currenta = panela(paneli)
       age=panelage(paneli)
       currentzi = panelz(paneli)
       currentlambdai = panellambda(paneli) 
       currentei = panele(paneli)
       
!       currentY=  (currenta+ ( rr * (zgrid(currentzi)*currenta)**mu - DepRate*currenta ) *(1.0_DP-tauK) )*(1.0_DP-tauW)
       
!  COMPUTE NEXT PERIOD'S ASSET
       if (age .lt. MaxAge) then
!            newpanela(paneli) = Linear_Int(agrid, Aprime(age,:,currentzi,currentlambdai, currentei),na,currenta)
!            newpanela(paneli) = Linear_Int(YGRID(:,currentzi), Aprime(age,:,currentzi,currentlambdai, currentei),na,currentY)
!            newpanela(paneli) = Linear_Int_Aprime(agrid, Aprime(age,:,currentzi,currentlambdai, currentei),na,currenta)

            ! do linear interpolation here to find aprime. calling the function takes much more time
            if (currenta .ge. amax) then
                tklo = na-1
                elseif (currenta .lt. amin) then
                    tklo = 1
                    else
                        tklo = ((currenta - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
            endif 
            
            tkhi = tklo + 1        

            newpanela(paneli) = ((agrid(tkhi) - currenta)*Aprime(age,tklo,currentzi,currentlambdai, currentei) &
                                    &  + (currenta - agrid(tklo))*Aprime(age,tkhi,currentzi,currentlambdai, currentei)) &
                                    &  / ( agrid(tkhi) - agrid(tklo) )    
            
            if (newpanela(paneli)  .ge. amax) then
                newpanela(paneli) = min(newpanela(paneli), amax) 
            endif      
            if (newpanela(paneli)  .lt. amin) then
                newpanela(paneli) = max(newpanela(paneli), amin) 
            endif      

       endif !age .lt. MaxAge
!  NEXT PERIOD'S ASSET IS COMPUTED

             
! DRAW NEXT PERIOD'S AGE DBN
       tempnoage = ran1(newiseed)  
  
     IF (tempnoage .gt. survP(age)) THEN
           newpanelage(paneli)=1
           ELSE
           newpanelage(paneli)=age+1
           newpanelz(paneli)=currentzi
           newpanellambda(paneli)=currentlambdai   
      ENDIF
  
! CORRECT AGE-Z-LAMBDA DEATH DISTRIBUTION 
!     
!       IF (tempnoage .gt. survP(age)) THEN
!           IF (deathby_age_z_lambda(age,currentzi,currentlambdai) .lt. &
!                                    & requireddeathby_age_z_lambda(age,currentzi,currentlambdai)) THEN
!                newpanelage(paneli)=1
!                deathby_age_z_lambda(age,currentzi,currentlambdai) = deathby_age_z_lambda(age,currentzi,currentlambdai)+1
!                ELSE
!                       newpanelage(paneli) =age+1
!                       survivingby_age_z_lambda(age,currentzi,currentlambdai) = &
!                                    &  survivingby_age_z_lambda(age,currentzi,currentlambdai)+1
!                       newpanelz(paneli)=currentzi
!                       newpanellambda(paneli)=currentlambdai
!                       numberby_age_z_lambda(age+1,currentzi,currentlambdai) = & 
!                                    & numberby_age_z_lambda(age+1,currentzi,currentlambdai) +1
!           ENDIF  
!      ELSE
!          IF (survivingby_age_z_lambda(age,currentzi,currentlambdai) .lt. &
!                                    & requiredsurvivingby_age_z_lambda(age,currentzi,currentlambdai)) THEN
!              newpanelage(paneli)=age+1
!              survivingby_age_z_lambda(age,currentzi,currentlambdai) = &
!                                    & survivingby_age_z_lambda(age,currentzi,currentlambdai)+1
!              newpanelz(paneli)=currentzi
!              newpanellambda(paneli)=currentlambdai
!              numberby_age_z_lambda(age+1,currentzi,currentlambdai) = &
!                                    & numberby_age_z_lambda(age+1,currentzi,currentlambdai) +1  
!              ELSE
!                  newpanelage(paneli)=1
!                  deathby_age_z_lambda(age,currentzi,currentlambdai) = deathby_age_z_lambda(age,currentzi,currentlambdai) +1               
!          ENDIF           
!      ENDIF
!      
! CORRECT AGE-Z-LAMBDA DEATH DISTRIBUTION ENDED
  

 
! DRAW Z and LAMBDA DISTRIBUTION FOR ONE-YEAR OLDS

    age = newpanelage(paneli)   
 
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
       ei=ne/2+1 ! ALL NEWBORN START FROM THE MEDIAN E  but if too many people died and started from median E, draw a new E for them

! CORRECT AGE-Z-LAMBDA  DISTRIBUTIONS
!
!       if (zi  .gt. 1) then
!           if ( (cdf_pr_z(currentzi,zi) -tempnoz) .lt.  (tempnoz-cdf_pr_z(currentzi,zi-1) ) )    then
!               agesign=1
!               else
!                    agesign=-1
!           endif
!       else
!            agesign=1          
!       endif
!       agecounter=1  
!       tzi=zi       
!115 IF (sum(numberby_age_z_lambda(age,zi,:)) .ge. sum(requirednumberby_age_z_lambda(age,zi,:))) then
!           zi = tzi + agecounter * agesign
!           zi = max(zi,1)
!           zi=min(zi,nz)
!           IF (sum(numberby_age_z_lambda(age,zi,:)) .ge. sum(requirednumberby_age_z_lambda(age,zi,:))) then
!               zi = zi - agecounter * agesign
!               zi = max(zi,1)
!               zi=min(zi,nz)
!               IF (sum(numberby_age_z_lambda(age,zi,:)) .ge. sum(requirednumberby_age_z_lambda(age,zi,:))) then
!                   agecounter = agecounter +1
!                   go to 115
!               ENDIF
!           ENDIF               
!       ENDIF    
!
!       if (lambdai .gt. 1) then
!            if ( (cdf_pr_lambda(currentlambdai,lambdai) -tempnolambda) .gt.  &
!                            & (tempnolambda-cdf_pr_lambda(currentlambdai,lambdai-1)) )   then
!               agesign=1
!               else
!                    agesign=-1
!           endif      
!       else
!            agesign=1          
!       endif       
!       agecounter=1  
!       tlambdai=lambdai  
!116 IF  (numberby_age_z_lambda(age,zi,lambdai) .ge. requirednumberby_age_z_lambda(age,zi,lambdai)) then
!               lambdai = tlambdai + agecounter * agesign
!               lambdai = max(lambdai,1)
!               lambdai=min(lambdai,nlambda)
!               IF  (numberby_age_z_lambda(age,zi,lambdai) .ge. requirednumberby_age_z_lambda(age,zi,lambdai)) then
!                   lambdai = lambdai - agecounter * agesign
!                   lambdai = max(lambdai,1)
!                   lambdai=min(lambdai,nlambda)
!                   IF  (numberby_age_z_lambda(age,zi,lambdai) .ge. requirednumberby_age_z_lambda(age,zi,lambdai)) then
!                       agecounter = agecounter +1
!                       go to 116
!                   ENDIF
!               ENDIF               
!         ENDIF
!          
!        if ( numberby_e_e(currentei,ei) .ge. requirednumberby_e_e(currentei,ei) ) then
!            tempnoe = ran1(newiseed) 
!            ei=1
!            DO WHILE (tempnoe .gt. cdf_Ge(ei))
!               ei=ei+1
!            ENDDO    
!
!           if (ei .gt. 1) then 
!               if ( (cdf_Ge(ei) -tempnoe) .gt.  (tempnolambda-cdf_Ge(ei-1)) )    then
!                   agesign=1
!                   else
!                        agesign=-1
!               endif      
!           else
!               agesign=1 
!           endif        
!           agecounter=1  
!           tei=ei      
! 117    IF (numberby_e_e(currentei,ei) .ge. requirednumberby_e_e(currentei,ei)) THEN
!                   ei = tei + agecounter * agesign
!                   ei = max(ei,1)
!                   ei=min(ei,ne)
!                   IF (numberby_e_e(currentei,ei) .ge. requirednumberby_e_e(currentei,ei)) THEN
!                       ei = tei -  agecounter * agesign
!                       ei = max(ei,1)
!                       ei=min(ei,ne)
!                       IF (numberby_e_e(currentei,ei) .ge. requirednumberby_e_e(currentei,ei)) THEN
!                           agecounter = agecounter +1
!                           go to 117
!                      ENDIF
!                   ENDIF        
!             ENDIF                
!        endif
!
!       numberby_age_z_lambda(age,zi,lambdai) = numberby_age_z_lambda(age,zi,lambdai) +1        
!       numberby_e_e(currentei, ei)  = numberby_e_e(currentei,ei) +1
!
!  CORRECTING DISTRIBUTIONS ENDED

        newpanelz(paneli)=zi    
        newpanellambda(paneli)=lambdai
        newpanele(paneli) = ei  
        
    ENDIF ! new age==1
 
ENDDO ! paneli


!E for surviving
!print*,'E'
DO paneli=1,totpop
!print*,'paneli',paneli
        age = newpanelage(paneli)  
        ! DRAW NEW E FOR THOSE WHO ARE NOT NEWBORN
        IF (age .gt. 1) THEN             
            currentei = panele(paneli)   
            tempno = ran1(newiseed)   
            ei=1
            DO WHILE (tempno .gt. cdf_pr_e(currentei,ei))
               ei=ei+1
            ENDDO
            
! CORRECT E DISTRIBUTION 
!
!           if (ei .gt. 1) then 
!                if ( (cdf_pr_e(currentei,ei) -tempnoe) .gt.  (tempnolambda-cdf_pr_e(currentei,ei-1)) )    then
!                   agesign=1
!                   else
!                        agesign=-1
!               endif    
!           else
!               agesign=1 
!           endif  
!           agecounter=1  
!           tei=ei      
! 118    IF (numberby_e_e(currentei,ei) .ge. requirednumberby_e_e(currentei,ei)) THEN
!                   ei = tei + agecounter * agesign
!                   ei = max(ei,1)
!                   ei=min(ei,ne)
!                   IF (numberby_e_e(currentei,ei) .ge. requirednumberby_e_e(currentei,ei)) THEN
!                       ei = tei -  agecounter * agesign
!                       ei = max(ei,1)
!                       ei=min(ei,ne)
!                       IF (numberby_e_e(currentei,ei) .ge. requirednumberby_e_e(currentei,ei)) THEN
!                           agecounter = agecounter +1
!                           go to 118
!                      ENDIF
!                   ENDIF        
!             ENDIF                
!             numberby_e_e(currentei,ei) = numberby_e_e(currentei,ei) + 1
!
! CORRECT E DISTRIBUTION ENDED

             newpanele(paneli)=ei            
       ENDIF ! age .gt. 1        
ENDDO ! paneli


panelage     = newpanelage
panela       = newpanela
panelz       = newpanelz
panellambda  = newpanellambda
panele       = newpanele

ENDDO ! simutime

DO paneli=1,totpop

panel_return(paneli)=(rr*mu*(zgrid(panelz(paneli))**mu)*(panela(paneli)**(mu-1.0_DP))-DepRate)   

   if (kstar(zi) .le. (vartheta*panela(paneli)) ) then
      panel_return(paneli)= int_rate 
      panel_at_return(paneli)= ((1.0_DP+int_rate*(1.0_DP-tauK))*(1.0_DP-tauW)-1.0_DP)
   else
      panel_return(paneli)= (int_rate &
         & + (rr*mu*((vartheta*zgrid(zi))**mu)*agrid(ai)**(mu-1.0_DP)-(int_rate+DepRate)*vartheta)) 
      panel_at_return(paneli)= ((1.0_DP+int_rate*(1.0_DP-tauK))*(1.0_DP-tauW) &
         & + (rr*mu*((vartheta*zgrid(zi))**mu)*agrid(ai)**(mu-1.0_DP) &
         & -(int_rate+DepRate)*vartheta)*(1.0_DP-tauK)*(1.0_DP-tauW) -1.0_DP)
   endif


       currenta = panela(paneli)
       age=panelage(paneli)
       currentzi = panelz(paneli)
       currentlambdai = panellambda(paneli) 
       currentei = panele(paneli)
       

            ! do linear interpolation here to find aprime. calling the function takes much more time
            if (currenta .ge. amax) then
                tklo = na-1
                elseif (currenta .lt. amin) then
                    tklo = 1
                    else
                        tklo = ((currenta - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
            endif 
            
            tkhi = tklo + 1        

            panelcons(paneli) = ((agrid(tkhi) - currenta)*cons(age,tklo,currentzi,currentlambdai, currentei) &
                                    &  + (currenta - agrid(tklo))*cons(age,tkhi,currentzi,currentlambdai, currentei)) &
                                    &  / ( agrid(tkhi) - agrid(tklo) )              

            panelhours(paneli) = ((agrid(tkhi) - currenta)*hours(age,tklo,currentzi,currentlambdai, currentei) &
                                    &  + (currenta - agrid(tklo))*hours(age,tkhi,currentzi,currentlambdai, currentei)) &
                                    &  / ( agrid(tkhi) - agrid(tklo) )  

            panelaprime(paneli) =((agrid(tkhi) - currenta)*Aprime(age,tklo,currentzi,currentlambdai, currentei) &
                                    &  + (currenta - agrid(tklo))*Aprime(age,tkhi,currentzi,currentlambdai, currentei)) &
                                    &  / ( agrid(tkhi) - agrid(tklo) )             
ENDDO ! paneli


WRITE  (UNIT=3, FMT=*) panela
WRITE  (UNIT=4, FMT=*) panelage 
WRITE  (UNIT=5, FMT=*) panelz 
WRITE  (UNIT=6, FMT=*) panellambda 
WRITE  (UNIT=7, FMT=*) panele 
WRITE  (UNIT=8, FMT=*) panel_return 
WRITE  (UNIT=9, FMT=*) panelcons
WRITE  (UNIT=10, FMT=*) panelhours
WRITE  (UNIT=11, FMT=*) panelaprime
WRITE  (UNIT=12, FMT=*) panel_at_return 

END SUBROUTINE SIMULATION

!====================================================================

!============================================================================

SUBROUTINE tauchen(mt,mut,rhot,sigmat,nt,gridt,prt,Gt)
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



! mut=0.0_DP
  gridt=0.0_DP
  prt=0.0_DP
  a=(1.0_DP-rhot)*mut;
        
        if (nt .gt. 1) then  
          gridt(nt)=mt*sqrt(sigmat**2.0_DP/(1.0_DP-rhot**2))
          gridt(1)=-gridt(nt)
          ELSE
!         print*,'only one grid. For this case transition probabilities might not be right. check it.'
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
!   print*, 'Gt', Gt, 'sum', sum(Gt)

  DO i=1,1000
    Gt_new=0.0_DP
    DO zi=1,nt
      DO zz=1,nt
        Gt_new(zz)=Gt_new(zz)+Gt(zi)*prt(zi,zz)
      END DO
    END DO
    Gt=Gt_new
  END DO
  
!   print*, 'Gt', Gt, 'sum', sum(Gt)
!   !pause

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
      !if (h.eq.0.) pause 'bad xa input in splint'  
      a=(xa(khi)-x)/h  
      b=(x-xa(klo))/h  
      y=a*ya(klo)+b*ya(khi)     +((a**3.0_DP-a)*y2a(klo)+(b**3.0_DP-b)*y2a(khi))*(h**2.0_DP)/6.0_DP  
! it also returns "yprime" (first derivative) at point "x"
!      yprime = ( ya(khi) - ya(klo) ) / h  & 
!                 & + h * ( ( 1.0_DP- 3.0_DP * a**2.0_DP  ) *  y2a(klo)  +  ( 3.0_DP * b**2.0_DP -1.0_DP ) *  y2a(khi)  )/6.0_DP
      return  

END  
 
!=======================================================================





