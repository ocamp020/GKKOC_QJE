module OBJECTIVE
    use genericParams
    Use nrtype
    Use Globals
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    PRIVATE
    PUBLIC objFun, dfovec, initial0
contains
    FUNCTION objFun(theta)
        use genericParams
        !use Parameters
        !use Globals
        implicit none
        integer :: i
        REAL(DP) :: objFun
        REAL(DP),DIMENSION(p_nx),INTENT(IN) :: theta

        REAL(DP) :: objFun0(1,1)
        REAL(DP),DIMENSION(nmoments,1) :: Fnout
        REAL(DP),DIMENSION(nmoments,nmoments) :: W
        REAL(DP) :: PAR(npar),pen_ub(npar), pen_lb(npar),pencons(npar),penalty

        CALL dfovec(npar, nmoments, theta, Fnout)

        PAR = theta
        PAR(2) = 1.-exp(-exp(theta(2))) 
        pen_ub=(/.5,0.9,1,1,.5,.5,1/)
        pen_Lb = (/0.001,0.1,-1,-1,0.001,0.001,0/)
        pencons = (/1000000D0,1000D0,10000D0,10000D0,1000000D0,1000000D0,1000D0/)
        penalty = 0D0

        do i = 1,npar
            penalty = penalty + pencons(i)*(min(0.D0,PAR(i)-pen_Lb(i)))**2D0
            penalty = penalty + pencons(i)*(min(0.D0,pen_Ub(i)-PAR(i)))**2d0 
        enddo

        objFun0=matmul(matmul(transpose(Fnout),W),Fnout)
        objFun = objFun0(1,1) + penalty;

 
    END FUNCTION objFun

    SUBROUTINE dfovec(n, mv, x, Fnout)
    !use nrtype
        use genericParams
        !use Parameters
        USE Globals
        IMPLICIT NONE

        INTEGER, INTENT(IN)     :: n, mv
        REAL(DP), DIMENSION(n), INTENT(IN)  :: x
        REAL(DP), DIMENSION(nmoments)  :: C_x, Cx
        REAL(DP), DIMENSION(mv,1),INTENT(OUT) :: Fnout
        REAL(DP), DIMENSION(mv,nsim) ::C_THETA_M
        REAL(DP), DIMENSION(mv) ::C_THETA_BAR, C_x_abs
        INTEGER :: i,i_m
        !Variables for sleep call
        !integer(c_int) :: mytime, dur
        real(dp) :: gama(mv)

        common Cx


        call Estimation

        C_x=Cx
!        allocate(C_THETA_M(size(C_x),nsim))
        !OPEN(2, FILE = datadir // '/Cx.dat')
        !do i_m=1,nmoments
         !   read(2,*) C_x(i_m)
        !enddo
        !CLOSE(2)
        
        
        !seed missing;
        do i_m = 1,nsim
            call f_moments(x,mv,C_THETA_M(:,i_m)) 
        enddo
        
        C_THETA_BAR=(1./real(nsim))*sum(C_THETA_M,2)
        C_x_abs= abs(C_x)

        ! Distance
        gama = 0.D0
        FORALL (i=1:mv) Fnout(i,1)=(C_THETA_BAR(i)-C_x(i))/(C_x_abs(i)+gama(i))


        !mytime=1;
        !dur=myFortSleep(mytime);

    END SUBROUTINE dfovec

    SUBROUTINE f_moments(theta,nm,C_M)
        Use Globals
        Use nrtype
        use genericParams
        use prctile
        USE random
        implicit none
        

        integer, intent(in) :: nm
        REAL(DP), DIMENSION(nm),intent(out) :: C_M
        REAL(DP),DIMENSION(npar),INTENT(IN) :: theta
        
        integer :: ii,jj,i_t,lag,isize, iseed(2)
        ! parameters
        real(dp) :: sig_eps, p1_eta, mu2_eta, mu3_eta, sig1_eta, sig2_eta, phi
        ! other constants
        real(dp) :: p2_eta, p3_eta, sig3_eta, mu_eps, mu_bar(T_sim)
        ! X series
        real(dp) :: yrgdp(T_sim+1),rgdp(T_sim+1), GDPgr(T_sim)
        real(dp) :: BC_X_1(T_sim), BC_X(T_data,nhhsim), mu_bar_t(T_sim,nhhsim)
        integer :: yrsample(T_sim)
        ! shocks
        real(dp), dimension(T_sim,nhhsim) :: eps, uniprob, eta1, eta2, eta3, eta
        real(dp) :: INCCHANGE_SRsim(T_sim-L_sr,nhhsim), &
                    INCCHANGE_MRsim(T_sim-L_mr,nhhsim), &
                    INCCHANGE_LRsim(T_sim-L_lr,nhhsim)
        real(dp) :: INCCHANGE_SR(T_data,nhhsim), &
                    INCCHANGE_MR(T_data+L_sr-L_mr,nhhsim), &
                    INCCHANGE_LR(T_data+L_sr-L_lr,nhhsim)
        real(dp) ::p10_SR(T_data),p50_SR(T_data),p90_SR(T_data),mea_SR(T_data),std_SR(T_data),ske_SR(T_data),kur_SR(T_data)
        real(dp) ::p10_MR(T_data+L_sr-L_mr),p50_MR(T_data+L_sr-L_mr),p90_MR(T_data+L_sr-L_mr),mea_MR(T_data+L_sr-L_mr),std_MR(T_data+L_sr-L_mr),ske_MR(T_data+L_sr-L_mr),kur_MR(T_data+L_sr-L_mr)
        real(dp) ::p10_LR(T_data+L_sr-L_lr),p50_LR(T_data+L_sr-L_lr),p90_LR(T_data+L_sr-L_lr),mea_LR(T_data+L_sr-L_lr),std_LR(T_data+L_sr-L_lr),ske_LR(T_data+L_sr-L_lr),kur_LR(T_data+L_sr-L_lr)
!integer, parameter :: leng
integer,parameter:: np=3
!real(8) :: xx(5),pout(np)
real(8) :: pp(np), pout(np)


        ! Read in parameters
        

        sig_eps = theta(1)
        p1_eta = 1.-exp(-exp(theta(2)))
        mu2_eta = theta(3)
        mu3_eta = theta(4)
        sig1_eta = theta(5)
        sig2_eta = theta(6)
        phi = theta(7)

        ! Done reading in parameters

        ! 1. X series
        ! -- GDP growth
        OPEN(1, FILE = datadir // '/Aggr_data/rgdp_y.dat')
        do ii=1,T_sim+1
            READ(1,*,iostat=iostatus) yrgdp(ii),rgdp(ii)
        enddo
        CLOSE(1)
        forall (ii=1:T_sim) GDPgr(ii) = (rgdp(ii+1)/rgdp(ii) - 1.0D0)
        
        BC_X_1 = -phi * GDPgr
        !NORMALIZE X HERE

        !forall (ii=1:nhhsim) BC_X(:,ii)=BC_X_1

        
        ! Temporary shock
        mu_eps = - log(exp(sig_eps**2/2))

        ! Persistent shock
        p2_eta = (1./2.) - p1_eta/2.
        p3_eta = p2_eta

        forall (ii=1:T_sim) mu_bar(ii) = - log( p1_eta*exp((sig1_eta**2.)/2.) + & 
                                                p2_eta*exp(mu2_eta - BC_X_1(ii) + (sig2_eta**2./2.))+ & 
                                                p3_eta*exp(mu3_eta - BC_X_1(ii) + (sig3_eta**2./2.)))
        !draw random numbers
        isize = 2
        iseed(1) = 6555
        iseed(2) = 5444
        CALL RANDOM_SEED(size = isize)
        CALL RANDOM_SEED(put = iseed)   
        call random_number(uniprob)

        
        do i_t = 1,T_sim
           do ii = 1,nhhsim
                eps(i_t,ii)  = mu_eps + sig_eps*random_normal()
                eta1(i_t,ii) = mu_bar(i_t) + sig1_eta*random_normal()
                eta2(i_t,ii) = mu_bar(i_t) + mu2_eta - BC_X_1(i_t) + sig2_eta*random_normal()
                eta3(i_t,ii) = mu_bar(i_t) + mu3_eta - BC_X_1(i_t) + sig3_eta*random_normal()

                if (uniprob(i_t,ii)<=p1_eta) then
                    eta(i_t,ii) = eta1(i_t,ii)
                elseif (uniprob(i_t,ii)<=p1_eta+p2_eta) then
                    eta(i_t,ii) = eta2(i_t,ii)
                else
                    eta(i_t,ii) = eta3(i_t,ii)
                endif
            enddo
        enddo
!        write(*,'(f7.4)'),mu_bar

        ! -- write eta
        !OPEN(2, FILE = datadir // '/eta.dat')
        !do i_t=1,T_sim
        !    write(2,'(100000f17.5)') eta(i_t,:)
        !enddo
        !CLOSE(2)

        INCCHANGE_SRsim=0.D0
        INCCHANGE_MRsim=0.D0
        INCCHANGE_LRsim=0.D0
        do i_t = 1,T_sim-L_sr
            INCCHANGE_SRsim(i_t,:) = eps(i_t+L_sr,:) - eps(i_t,:)
            do lag = 1,L_sr
                INCCHANGE_SRsim(i_t,:) = INCCHANGE_SRsim(i_t,:) + eta(i_t+lag,:)
            enddo
        enddo
        do i_t = 1,T_sim-L_mr
            INCCHANGE_MRsim(i_t,:) = eps(i_t+L_mr,:) - eps(i_t,:)
            do lag = 1,L_mr
                INCCHANGE_MRsim(i_t,:) = INCCHANGE_MRsim(i_t,:) + eta(i_t+lag,:)
            enddo
        enddo
        do i_t = 1,T_sim-L_lr
            INCCHANGE_LRsim(i_t,:) = eps(i_t+L_lr,:) - eps(i_t,:)
            do lag = 1,L_lr
                INCCHANGE_LRsim(i_t,:) = INCCHANGE_LRsim(i_t,:) + eta(i_t+lag,:)
            enddo
        enddo
        !if (country==1 .or. country==3) then
        !    INCCHANGE_SR = INCCHANGE_SRsim
        !    INCCHANGE_MR = INCCHANGE_MRsim
        !    INCCHANGE_LR = INCCHANGE_LRsim
        !elseif (country == 2) then
            i_t=19
            jj=i_t
            INCCHANGE_SR(1:i_t-1,:)=INCCHANGE_SRsim(1:i_t-1,:)
            do 
                if (i_t>size(INCCHANGE_SRsim,1)) exit
                INCCHANGE_SR(jj,:) = INCCHANGE_SRsim(i_t,:)
                jj=jj+1
                i_t=i_t+2
            enddo
            i_t=17
            jj=i_t
            INCCHANGE_MR(1:i_t-1,:)=INCCHANGE_MRsim(1:i_t-1,:)
            do 
                if (i_t>size(INCCHANGE_MRsim,1)) exit
                INCCHANGE_MR(jj,:) = INCCHANGE_MRsim(i_t,:)
                jj=jj+1
                i_t=i_t+2
            enddo
            i_t=15
            jj=i_t
            INCCHANGE_LR(1:i_t-1,:)=INCCHANGE_LRsim(1:i_t-1,:)
            do 
                if (i_t>size(INCCHANGE_LRsim,1)) exit
                INCCHANGE_LR(jj,:) = INCCHANGE_LRsim(i_t,:)
                jj=jj+1
                i_t=i_t+2
            enddo
        !endif ! country

        pp=(/.1,.5,.9/)
        do i_t = 1,T_data
            call percentile(INCCHANGE_SR(i_t,:),nhhsim,np,pp,pout) 
            p10_SR(i_t)=pout(1)
            p50_SR(i_t)=pout(2)
            p90_SR(i_t)=pout(3)
        enddo
        do i_t = 1,T_data+L_sr-L_mr
            call percentile(INCCHANGE_MR(i_t,:),nhhsim,np,pp,pout) 
            p10_MR(i_t)=pout(1)
            p50_MR(i_t)=pout(2)
            p90_MR(i_t)=pout(3)
        enddo
        do i_t = 1,T_data+L_sr-L_lr
            call percentile(INCCHANGE_LR(i_t,:),nhhsim,np,pp,pout) 
            p10_LR(i_t)=pout(1)
            p50_LR(i_t)=pout(2)
            p90_LR(i_t)=pout(3)
        enddo
    C_M = (/p10_SR, p10_MR, p10_LR, &
            p50_SR, p50_MR, p50_LR, &
            p90_SR, p90_MR, p90_LR/)
!print*,C_M


        ! -- write eta
        !OPEN(2, FILE = datadir // '/c_m.dat')
        !do i_t=1,nmoments
        !    write(2,*) C_M(i_t)
        !    enddo
        !CLOSE(2)


    END SUBROUTINE f_moments
    SUBROUTINE initial0
        USE global
        IMPLICIT NONE
    END SUBROUTINE initial0
    
SUBROUTINE Estimation
use genericParams
        Use nrtype

IMPLICIT NONE


!Declare Local Variables
!----------------------------------
INTEGER                                 :: info,ip,ia,it,ii,jj
REAL(dp), DIMENSION(:), ALLOCATABLE      :: Xguess, Xsol
REAL(dp), DIMENSION(nmoments)            :: Fsol
REAL(dp), dimension(:,:), allocatable :: ones, X, Y, YY, XX, Yres
integer, dimension(:,:), allocatable :: isnanY
REAL(dp) :: bet(2,1), meanY
integer :: Nnan
REAL(dp) :: Cx0(nn), Cx(nm)
real(dp):: C_THETA(nmoments)
common Cx
! Country-specific simulation and estimation options


    allocate(table(T_data,22))
    allocate(data_d(T_data,22))
    allocate(X(T_data,2))
    allocate(Y(T_data,1),isnanY(T_data,1))

    OPEN(1, FILE = datadir // '/moment_est_pre.txt')
    read(1,*)
    do ii=1,T_data
        READ(1,*,iostat=iostatus) table(ii,:)
    enddo
    CLOSE(1)

    ! Detrend
    data_d = table
    where(data_d == -1.0)
        data_d =  sqrt(-1.0)
    endwhere 
    allocate(ones(T_data,1))
    ones = 1.0

    X(:,1) = ones(:,1)
    X(:,2) = data_d(:,1)
    do ii = 2,size(data_d,2)
        Y(:,1) = data_d(:,ii)
        where (isnan(Y) == .TRUE.)
            isnanY = 1
        elsewhere
            isnanY = 0
        endwhere
        Nnan = sum(isnanY(:,1))
        allocate(YY(T_data-Nnan,1),Yres(T_data-Nnan,1),XX(T_data-Nnan,2))
        ia = 1
        it = 1
        do 
            if ((ia>size(Y,1))) exit
                if (isnanY(ia,1)==0) then
                    YY(it,1) = Y(ia,1)
                    XX(it,1) = X(ia,1)
                    XX(it,2) = X(ia,2)
                    ia = ia + 1
                    it = it + 1
                else
                    ia = ia + 1
                endif
        enddo
        bet = matmul(inver(matmul(transpose(XX),XX)),(matmul(transpose(XX),YY)))
        meanY = sum(YY(:,1))/size(YY(:,1))
        Yres = YY - matmul(XX,bet)
        ia = 1
        it = 1
        do 
            if ((ia>size(Y,1))) exit
                if (isnanY(ia,1)==0) then
                    Y(ia,1) = Yres(it,1)+meanY
                    ia = ia + 1
                    it = it + 1
                else
                    ia = ia + 1
                endif
        enddo
        deallocate(YY,XX,Yres)
        data_d(:,ii) = Y(:,1)
    enddo
    ! Done detrending

       !allocate(Cx0(nn),Cx(nm))

    Cx0=reshape(data_d(:,2:10),(/nn/))
    ia=1
    it=1
    do
        if (ia>nn) exit
            if (isnan(Cx0(ia))==.FALSE.) then
                Cx(it) = Cx0(ia)
                ia = ia+1
                it = it+1
            else
                ia = ia+1
            endif
    enddo

    ! Initial Value
!print*,theta0
!theta0=(/0.1580D0,0.8948D0,0.355D0,-.298D0,0.0143D0, 0.1041D0,1.D0/)
!call f_moments(theta0,nmoments,C_THETA)
! call dfovec(npar,nmoments,theta0,Cx,Fout)
!Fval= objFun(THETA0,Cx)
    !call f_moments(THETA0)
    deallocate(table,data_d,X,Y,isnanY)

        ! -- write eta
        !OPEN(2, FILE = datadir // '/Cx.dat')
        !do ia=1,nmoments
         !   write(2,*) Cx(ia)
        !enddo
        !CLOSE(2)
END SUBROUTINE Estimation
function inver(matrix)
	        !This function computes the inverse of matrix 'matrix(lda,n)' calling LAPACK routines
	        ! dgetrf and dgetri
	        !
        implicit none

        REAL(dp):: matrix(:,:)
        integer:: info
        integer:: lwork
        integer:: lda,n
        integer, dimension(size(matrix,2)):: ipiv
        REAL(dp), dimension(size(matrix,2))::work
        REAL(dp), dimension(size(matrix,2),size(matrix,2)):: inver, M

        lda=size(matrix,1)
        n=size(matrix,2)
        M=matrix

        !First compute optimal lwork
        call DGETRF( n, n, M, lda, ipiv, info )
        call DGETRI( n, M, lda, ipiv, work, -1, info)
        !Now, use it to get the inverse
        lwork=min(1000,int(work(1)))
        M=matrix
        call DGETRF( n, n, M, lda, ipiv, info )
        call DGETRI( n, M, lda, ipiv, work, lwork, info)

        if (info/=0) then
            write(*,*) 'Careful with inverse for',matrix,'!, info=', info
            stop
        else
            inver=M
        end if

    end function inver

end MODULE objective
