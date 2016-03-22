module OBJECTIVE
    use, intrinsic :: iso_c_binding, only: c_int
    use omp_lib
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
        use omp_lib
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
        use GKK_Calibration
        use omp_lib
        IMPLICIT NONE

        INTEGER, INTENT(IN)     :: n, mv
        REAL(DP), DIMENSION(n), INTENT(IN)  :: x
        REAL(DP), DIMENSION(mv),INTENT(OUT) :: v_err
        INTEGER :: i
        REAL(DP), DIMENSION(p_nx) :: values

        !Variables for sleep call
        integer(c_int) :: mytime, dur


        v_err(1) = Moments_Objective(x)
        
        mytime=1;
        !dur=myFortSleep(mytime);

    END SUBROUTINE dfovec

    SUBROUTINE initial0
        USE global
        IMPLICIT NONE
    END SUBROUTINE initial0
end MODULE objective
