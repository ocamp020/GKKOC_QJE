MODULE Toolbox
	USE nrtype
	use nrutil
	IMPLICIT NONE

	Contains
		! Linear_Int
		! Brent
		! ran1
		! Sort
		! ZBrent
		! ZBrent_p
		! cdf_normal
		! Tauchen
		! Spline
		! Splint
		
! =============================================================================
!
! Linear_IP: Conducts Linear Interpolation to compute f(x) given a grid
!              Given a function y=f(x) tabulated in xvec and yvec and a point x0,
!              return the function value y0=f(x0) obtained from linear interpolations
!              between x1<x<x2.
!
! Usage: y0 = Linear_IP(n,xvec,yvec,x0)
!
! Input: n_grid, integer ,                    Size of the grid
!        x_grid, real(dp), dimension(n_grid), Values of the grid
!        f_grid, real(dp), dimension(n_grid), Values of the objective function at the grid
!        x     , real(dp),                  , Point for interpolation
!
! Output: y0, real(dp), Value of function at x.
!
! Remarks: Taken from Heer & Maussner (2nd Edition) - Original source in file Funcion.for

    Function Linear_IP(n_grid,x_grid,f_grid,x)
        implicit none
        integer  :: n_grid, j
        real(dp) :: x, x_grid(n_grid), f_grid(n_grid), Linear_IP

    
        ! if ((x.lt.x_grid(1)) .or. (x.gt.x_grid(n_grid))) then
        !     Print *, "Linear Interpolation Error: Input off of grid!"
        !     Linear_IP = x
        !     Return
        ! end if
    
        if (x.lt.x_grid(1)) then
            Linear_IP = f_grid(1)
        elseif (x.gt.x_grid(n_grid)) then
            Linear_IP = f_grid(n_grid)
        else
            j = count(x_grid.le.x) ! this determines the lower bracket for x
            Linear_IP=f_grid(j)+((f_grid(j+1)-f_grid(j))/(x_grid(j+1)-x_grid(j)))*(x-x_grid(j))
        end if

        Return

    End Function Linear_IP

!========================================================================================
!========================================================================================
! Given arrays xa(1:n) and ya(1:n) of length n; this subroutine returns  linear interpolated value "y" at point "x"
	FUNCTION Linear_Int(xa,ya,n,x)  
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
		      	print*,'bad xa input in linear int', khi,klo, xa(khi), xa(klo),x
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

	FUNCTION brent_p(ax,bx,cx,func,tol,xmin,par)
		USE nrtype; USE nrutil, ONLY : nrerror
		IMPLICIT NONE
		
		REAL(DP), INTENT(IN) :: ax,bx,cx,tol
		REAL(DP), INTENT(OUT) :: xmin
		REAL(DP), dimension(:), INTENT(IN) :: par
		REAL(DP) :: brent_p
		INTERFACE
			FUNCTION func(x,par)
				USE nrtype
				IMPLICIT NONE
				REAL(DP), INTENT(IN) :: x
				REAL(DP), dimension(:), INTENT(IN) :: par
				REAL(DP) :: func
			END FUNCTION func
		END INTERFACE
		INTEGER(I4B), PARAMETER :: ITMAX=1000
		REAL(DP), PARAMETER :: CGOLD=0.3819660_DP,ZEPS=1.0e-3_DP*epsilon(ax)
		INTEGER(I4B) :: iter
		REAL(DP) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
		
		a=min(ax,cx)
		b=max(ax,cx)
		v=bx
		w=v
		x=v
		e=0.0
		fx=func(x,par)
		fv=fx
		fw=fx
		do iter=1,ITMAX

			xm=0.5_DP*(a+b)
			tol1=tol*abs(x)+ZEPS
			tol2=2.0_DP*tol1
			if (abs(x-xm) <= (tol2-0.5_DP*(b-a))) then
				xmin=x
				brent_p=fx
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
			fu=func(u,par)
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
		call nrerror('brent_p: exceed maximum iterations')
		
		CONTAINS
		SUBROUTINE shft(a,b,c,d)
			REAL(DP), INTENT(OUT) :: a
			REAL(DP), INTENT(INOUT) :: b,c
			REAL(DP), INTENT(IN) :: d
			a=b
			b=c
			c=d
		END SUBROUTINE shft
	END FUNCTION brent_p

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
		if (fa<0.0_dp) then
			zbrent = a 
			print*, 'Lower bound is not low enough - This message is for Agg_Debt', a,fa
			Return 
		endif 
		fb=func(b)
		if (fb>0.0_dp) then 
			zbrent = b
			print*, 'Higher bound is not high enough - This message is for Agg_Debt' , b, fb 
			Return 
		endif 
		! if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) &
		! 	call nrerror('root must be bracketed for zbrent')
	
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

!========================================================================================
!========================================================================================

	FUNCTION cdf_normal(x)
		use nrtype
		use nrutil
		IMPLICIT NONE

		real(DP) :: cdf_normal

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

	END FUNCTION cdf_normal

!========================================================================================
!========================================================================================


	SUBROUTINE tauchen(mt,rhot,sigmat,nt,gridt,prt,Gt)
		use nrtype
    	use nrutil
		IMPLICIT NONE
		
		REAL(DP), INTENT(IN) :: rhot, sigmat
		INTEGER(I4B),  INTENT(IN) :: nt
		REAL(DP), INTENT(OUT), DIMENSION(nt)    :: gridt, Gt
		REAL(DP), INTENT(OUT), DIMENSION(nt,nt) :: prt
		REAL(DP), DIMENSION(nt)    :: Gt_new
		REAL(DP) :: a, stept, mut
		REAL(DP), INTENT(IN) ::  mt
		INTEGER(I4B)  :: i, j, k, zi, zz


		mut=0.0_DP
		gridt=0.0_DP
		prt=0.0_DP
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
		elseif (nt.eq.1) then
			prt = 1.0_DP
			Gt  = 1.0_DP 
		else
			print*, "Tauchen: Fatal error, nt is not greater than or equal to 1"
			call EXIT
		endif 
		
	END SUBROUTINE tauchen 

! =============================================================================
! MarkovAR_95: Approximates a contiuous AR(1) process with a discrete Markov process - Rowenhorst 95
!
! Usage: Call MarkovAR_95(n_z,rho,sigma,z_grid,P)
!
! Input: n_z    , integer , Grid dimension for Markov process
!        rho    , real(dp), Persistence of x where x'=rho*x+sigma*e and e~N(0,1)
!        sigma  , real(dp), Standard Deviation of shocks to x where x'=rho*x+sigma*e and e~N(0,1)
!
! Output: z_grid , real(dp), dimension(n_z),     Grid of values of the Markov process. Size n_z.
!         P      , real(dp), dimension(n_z,n_z), Transition Matrix of z. Prob sum across rows.
!
! Remarks: This routine generates a grid for z (z_grid) and a transition probability matrix (P)
!          The transition probability matrix is organized (z,z') so that it sums to 1 by rows
! 		   P(i,j) = p_ij is the transition prob of going from z=z_i to z'=z_j
!          Note that this creates matrixes in the opposite order than MarkovAR that follows Tauchen
!          The method follows Rouwenhorst (1995) as shown in Kopecky and Suen (2010)
	
	Recursive Subroutine MarkovAR_95(n_z,rho,sigma,z_grid,P_z)
	    integer, intent(in) :: n_z
	    real(dp), intent(in) :: rho, sigma
	    real(dp), intent(out), dimension(n_z)     :: z_grid
	    real(dp), intent(out), dimension(n_z,n_z) :: P_z
	    integer :: i, j
	    real(dp) :: step, p, q, psi, P_2(2,2)
	    real(dp) :: z_grid_aux(n_z-1), P_aux(n_z-1,n_z-1), P_a(n_z+1,n_z+1), P_b(n_z,n_z), P_half(n_z,n_z)

	    ! Parameters p, q and psi
	        p = (1.0_dp+rho)/2.0_dp
	        q = (1.0_dp+rho)/2.0_dp
	        psi = sqrt(real(n_z)-1.0_dp)*sigma/sqrt(1.0_dp-rho**2.0_dp)

	    ! Step of grid
	        step = 2.0_dp*psi/(real(n_z)-1.0_dp)
	    
	    ! Compute z_grid
	        z_grid(1) = -psi
	        do i=2,n_z
	            z_grid(i) = z_grid(i-1) + step
	        end do
	    
	    ! Compute transition matrix for n_z=2
	        P_2 = reshape((/p,1.0_dp-q,1.0_dp-p,q/),(/2,2/))

	    ! Compute transition matrix for arbitrary n_z 
	        if (n_z>2) then
	            Call MarkovAR_95(n_z-1,rho,sigma,z_grid_aux,P_aux)
		        ! To create P_n you take P which is n_z -1 by n_z - 1 and create
		        ! 4 matrixes which have P in one corner and 0's in the row and column
		        ! touching the other corner.  For example,
		        ! [1 2; 3 4] => [1 2 0; 3 4 0; 0 0 0], [ 0 1 2; 0 3 4; 0 0 0] ...
		        ! plus
		        P_a = 0.0_dp
		        P_a(2:n_z, 2:n_z) = P_aux
		        P_b = ( p*P_a(2:n_z+1,2:n_z+1) + (1.0_dp-p)*P_a(2:n_z+1,1:n_z) + &
		              & (1.0_dp-q)*P_a(1:n_z,2:n_z+1) + q*P_a(1:n_z,1:n_z) )
		        P_half(1,:)       =  1.0_dp
		        P_half(2:n_z-1,:) = 0.50_dp
		        P_half(n_z,:)     = 1.0_dp
		        P_z = P_b*P_half
	        else
	            P_z = P_2
	        end if

	    return
	end subroutine  MarkovAR_95

! =============================================================================
! Markov_Cut: Cuts down the first n_cut nodes of a discrete markov process. 
!			  This subroutine takes in an "n" state markov process with transition matrix P_in
! 			  and invariant distribution G_in, and returns a process with "n-n_cut" states and
!			  transition matrix P and invariant distribution G
!			  The transition matrix is assumed to sum trhough columns so that p(i,j) = Pr(z'=j|z=i)
!			  The first n_cut+1 states of the original grid are merged into the first state of the 
!			  new grid
!
! Usage: Call Markov_Cut(n,grid_in,P_in,G_in,n_cut,grid,P,G)
!
! Input:  n      , integer , 			    Grid dimension for original Markov process
! 		  n_cut  , integer , 			    Number of states to eliminate. New process has n-n_cut states
!		  grid_in, real(dp), dimension(n)  , Grid of states of original markov process
!		  P_in   , real(dp), dimension(n,n), Transition matrix of original markov process
!		  G_in   , real(dp), dimension(n)  , Invariant distribution of original markov process
!
! Output: grid   , real(dp), dimension(n-n_cut)        , Grid of states of new markov process
!		  P      , real(dp), dimension(n-n_cut,n-n_cut), Transition matrix of new markov process
!		  G      , real(dp), dimension(n-n_cut)        , Invariant distribution of new markov process
!
	
	Subroutine Markov_Cut(n,grid_in,P_in,G_in,n_cut,grid,P,G)
		use nrtype
		IMPLICIT NONE
		integer , intent(in)  :: n, n_cut
		real(dp), intent(in)  :: grid_in(n), P_in(n,n), G_in(n)
		real(dp), intent(out) :: grid(n-n_cut), P(n-n_cut,n-n_cut), G(n-n_cut)
		real(dp)              :: G_cut(n_cut+1), P_aux(n,n-n_cut)


		if (n_cut.ge.n) then
			print*, "Error in Markov_Cut: n_cut>=n"
			STOP
		elseif (n_cut.eq.0) then
			print*, "Warning in Markov_Cut: Nothing to cut (n_cut=0)"
			print*, "Output is unchanged."
			grid = grid_in
			P 	 = P_in
			G 	 = G_in
		else
			! Distribution of nodes to be eliminated
			G_cut = G_in(1:n_cut+1)/sum(G_in(1:n_cut+1))

			! Grid 
			grid(1)  = sum(grid_in(1:n_cut+1)*G_cut)
			grid(2:) = grid_in(n_cut+2:)

			! Transition Matrix
				! Sum transition probabilities into first n_cut+1 states
				P_aux(:,1)  = sum(P_in(:,1:n_cut+1),2)
				P_aux(:,2:) = P_in(:,n_cut+2:)
				! Average across the states to be merged with G_cut
				P(1,:)      = sum(P_aux(1:n_cut+1,:)*spread(G_cut,2,n-n_cut),1)
				P(2:,:)     = P_aux(n_cut+2:,:)

			! Satationary Distribution
			G(1)  = sum(G_in(1:n_cut+1))
			G(2:) = G_in(n_cut+2:)
		endif

	end Subroutine Markov_Cut


!========================================================================================
!========================================================================================


	SUBROUTINE spline(x,y,n,yp1,ypn,y2)  
		! Given arrays x(1:n) and y(1:n) of length n; and first derivatives yp1 and ypn at points 1 and n, 
		! this subroutine returns an array y2(1:n) (2nd derivative) at points x(1:n) parameter NMAX is the largest anticipated value of "n"
		 
		use nrtype
    	use nrutil
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
      

!========================================================================================
!========================================================================================
       
 
	SUBROUTINE splint(xa,ya,y2a,n,x,y)  
		! Given arrays xa(1:n) and ya(1:n) of length n; and second derivatives y2a(1:n), which is the output of spline
		! given the value of "x" this subroutine returns  a cubic-spline interpolated value "y" at point "x"

		use nrtype
    	use nrutil
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

	END Subroutine splint
 

!========================================================================================
!========================================================================================

end MODULE Toolbox