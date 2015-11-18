module OBJECTIVE
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    PRIVATE
    PUBLIC objFun, dfovec, initial0

   interface

      function myFortSleep (seconds)  bind ( C, name="sleep" )
          import
          integer (c_int) :: myFortSleep
          integer (c_int), intent (in), VALUE :: seconds
      end function myFortSleep

   end interface

contains

    FUNCTION objFun(theta)
        use genericParams
        implicit none
        REAL(DP),DIMENSION(p_nx),INTENT(IN) :: theta
        REAL(DP) :: objFun
        REAL(DP),DIMENSION(p_nmom) :: v_err

        CALL dfovec(p_nx, p_nmom, theta, v_err)
        objFun=v_err(1)
        print *,v_err(1)


    END FUNCTION objFun

    SUBROUTINE dfovec(n, mv, x, v_err)
        USE nrtype
        USE global
        use genericParams
        IMPLICIT NONE

        INTEGER, INTENT(IN)     :: n, mv
        REAL(DP), DIMENSION(n), INTENT(IN)  :: x
        REAL(DP), DIMENSION(mv),INTENT(OUT) :: v_err
        INTEGER :: i
        REAL(DP), DIMENSION(p_nx) :: values

        !Variables for sleep call
        integer(c_int) :: mytime, dur

        !values(1)=0.7
        !values(2)=0.2
        !values(3)=12
        !values(4)=0.87

        !FORALL (i=1:mv) v_err(i)=(x(i)-values(i))/values(i)

        v_err(1)=(SIN(3.141592653589793D0*x(1)*2.0D0)**2+1.0D0)*(x(2)-1.0D0)**2+ &
               (SIN(3.141592653589793D0*x(2)*3.0D0)**2+1.0D0)*(x(1)-1.0D0)**2+ &
                SIN(3.141592653589793D0*x(1)*3.0D0)**2 
        
        mytime=1;
        !dur=myFortSleep(mytime);

    END SUBROUTINE dfovec

    SUBROUTINE initial0
        USE global
        IMPLICIT NONE
    END SUBROUTINE initial0
end MODULE objective
