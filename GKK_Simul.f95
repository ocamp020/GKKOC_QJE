Program GKK_Simul
  use parameters
  use global
  use programfunctions
  use Toolbox
  Implicit None
  character(100) :: Bench_Folder

  ! Set Parameters 
    Params =[ 0.9436, 0.00, 0.50, 0.70444445, 0.34, 0.4494 ] ! tauL=0.224, tauC=0.075 calibration
    beta   = params(1)
    mu_z   = params(2) ! this is just shifting the z grids. it is zero now.
    rho_z  = params(3) 
    sigma_z_eps      = params(4)
    sigma_lambda_eps = params(5)
    gamma  = params(6)
    sigma  = 4.0_dp

  ! Taxes
  ! Wealth tax: minimum wealth tax to consider and increments for balancing budget
    TauW_bt          = 0.00_dp
    Threshold_Factor = 0.00_dp 
  ! Consumption tax
    tauC=0.075_DP
  ! Set Labor Tax Regime
    !tauPL=0.185_DP
    !psi=0.77_DP  
    tauPL=0.0_DP
    psi=0.776_DP    

  ! Resutls Folder
    write(Result_Folder,'(f4.2)') Threshold_Factor

    if ((TauPL.eq.0.0_dp).and.(sigma.ne.1.0_dp)) then 
      Result_Folder = './NSU_LT_Results/Factor_'//trim(Result_Folder)//'/'
    else if ((TauPL.ne.0.0_dp).and.(sigma.ne.1.0_dp)) then 
      Result_Folder = './NSU_PT_Results/Factor_'//trim(Result_Folder)//'/'
    else if ((TauPL.eq.0.0_dp).and.(sigma.eq.1.0_dp)) then 
      Result_Folder = './SU_LT_Results/Factor_'//trim(Result_Folder)//'/'
    else if ((TauPL.ne.0.0_dp).and.(sigma.eq.1.0_dp)) then 
      Result_Folder = './SU_PT_Results/Factor_'//trim(Result_Folder)//'/'
    end if 
    print*, "Results are stored in directory: ", Result_Folder

  ! Bench_Folder
    if ((TauPL.eq.0.0_dp).and.(sigma.ne.1.0_dp)) then 
      Bench_Folder = './NSU_LT_Results/Bench_Files/'
    else if ((TauPL.ne.0.0_dp).and.(sigma.ne.1.0_dp)) then 
      Bench_Folder = './NSU_PT_Results/Bench_Files/'
    else if ((TauPL.eq.0.0_dp).and.(sigma.eq.1.0_dp)) then 
      Bench_Folder = './SU_LT_Results/Bench_Files/'
    else if ((TauPL.ne.0.0_dp).and.(sigma.eq.1.0_dp)) then 
      Bench_Folder = './SU_PT_Results/Bench_Files/'
    end if 

  ! Initialize program and load functions
    print*,"  Initializing program"
    CALL INITIALIZE

!====================================================================================================
  PRINT*,''
  Print*,'Loading benchmark'
  PRINT*,''

  ! Set taxes for benchmark economy
    tauK = 0.25_DP
    tauW_at = 0.00_DP
    Y_a_threshold = 0.00_DP 

  print*,"  Reading benchmark results from files"
    CALL Write_Benchmark_Results(1)

  print*,"  Simulation benchmark"
  CALL SIMULATION(1)
  print*,"  End of Simulation benchmark"
 

!====================================================================================================
  PRINT*,''
  Print*,'Loading experiment'
  PRINT*,''

  tauK = 0.0_DP
  Y_a_threshold = Threshold_Factor*Ebar_bench 
  tauW_at = 0.017072675596579098_dp

  OPEN  (UNIT=1,  FILE=trim(Result_Folder)//'Exp_results_cons'  , STATUS='old', ACTION='read')
  OPEN  (UNIT=2,  FILE=trim(Result_Folder)//'Exp_results_aprime', STATUS='old', ACTION='read')
  OPEN  (UNIT=3,  FILE=trim(Result_Folder)//'Exp_results_hours' , STATUS='old', ACTION='read')
  OPEN  (UNIT=4,  FILE=trim(Result_Folder)//'Exp_results_value' , STATUS='old', ACTION='read')
  OPEN  (UNIT=5,  FILE=trim(Result_Folder)//'Exp_results_DBN'   , STATUS='old', ACTION='read')
  OPEN  (UNIT=60, FILE=trim(Result_Folder)//'Exp_results_GBAR'  , STATUS='old', ACTION='read')
  OPEN  (UNIT=7,  FILE=trim(Result_Folder)//'Exp_results_EBAR'  , STATUS='old', ACTION='read')
  OPEN  (UNIT=8,  FILE=trim(Result_Folder)//'Exp_results_NBAR'  , STATUS='old', ACTION='read')
  OPEN  (UNIT=9,  FILE=trim(Result_Folder)//'Exp_results_QBAR'  , STATUS='old', ACTION='read')
  OPEN  (UNIT=10, FILE=trim(Result_Folder)//'Exp_results_rr'    , STATUS='old', ACTION='read')
  OPEN  (UNIT=11, FILE=trim(Result_Folder)//'Exp_results_wage'  , STATUS='old', ACTION='read')
  OPEN  (UNIT=12, FILE=trim(Result_Folder)//'Exp_results_YBAR'  , STATUS='old', ACTION='read')

  READ (UNIT=1,  FMT=*), Cons
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

  print*,"  Simulation experiment"
  CALL SIMULATION(0)
  print*,"  End Simulation experiment"



Contains 
  !====================================================================

  SUBROUTINE SIMULATION(ben_switch)
    use parameters
    USE GLOBAL
    use programfunctions
    use Toolbox
    IMPLICIT NONE
    integer  :: currentzi, currentlambdai, currentei
    REAL(DP) :: tempnoage, tempnoz, tempnolambda, tempnoe, tempno, currenta, currentY
    REAL(DP) :: start_timet, finish_timet
    INTEGER  :: agecounter, agesign, tage, tzi, tlambdai, tei, tklo, tkhi
    INTEGER, PARAMETER  :: totpop=4000000, cohortsize=50000,  Max_Simu_Time=1100, LastTperiod=100
    INTEGER, intent(in) :: ben_switch
    INTEGER :: i


    REAL(DP):: lambdaBAR
    REAL(DP) :: m, Rh
    INTEGER :: ee0, ee1, ee2, zindx1, zindx2, lambdaindx1, lambdaindx2, diff_array, eindx1, eindx2
    INTEGER, DIMENSION(RetAge) :: agevec

    ! PANEL VARIABLES
      integer,DIMENSION(MaxAge)    :: requirednumberby_age, cdfrequirednumberby_age, numberby_age
      integer,DIMENSION(MaxAge)    :: requireddeathby_age, deathby_age, requiredsurvivingby_age, survivingby_age
      integer,DIMENSION(MaxAge,nz) :: requirednumberby_age_z, Numberby_age_z
      !integer,DIMENSION(MaxAge,nlambda)   ::  requirednumberby_age_lambda, numberby_age_lambda

      integer,DIMENSION(MaxAge,nz,nlambda)   :: requirednumberby_age_z_lambda, Numberby_age_z_lambda
      integer,DIMENSION(MaxAge,nz,nlambda)   :: requireddeathby_age_z_lambda, deathby_age_z_lambda
      integer,DIMENSION(MaxAge,nz,nlambda)   :: requiredsurvivingby_age_z_lambda, survivingby_age_z_lambda

      !integer, DIMENSION(MaxAge,nz,nlambda,ne) ::  requirednumberby_age_z_lambda_e, numberby_age_z_lambda_e
      !integer, DIMENSION(MaxAge,nz,nlambda,ne) ::  requireddeathby_age_z_lambda_e, deathby_age_z_lambda_e
      !integer, DIMENSION(MaxAge,nz,nlambda,ne) ::  requiredsurvivingby_age_z_lambd_e, survivingby_age_z_lambda_e


      integer,DIMENSION(MaxAge,ne)   ::  requirednumberby_age_e, numberby_age_e
      integer, DIMENSION(ne)      :: requirednumberby_e   , numberby_e
      integer, DIMENSION(MaxAge, ne)      :: requireddeathby_age_e, requiredsurvivingby_age_e, deathby_age_e, survivingby_age_e
      integer, DIMENSION(ne,ne) :: requirednumberby_e_e, numberby_e_e


      integer :: paneli, oldiseed, simutime 
      real(DP), DIMENSION(Max_Simu_Time) :: NBAR_seq, QBAR_seq
      integer, DIMENSION(totpop) :: panelage, newpanelage, panelz, newpanelz, newpanele, panele, newpanellambda, panellambda
      real(DP), DIMENSION(totpop) :: panela, newpanela, newpaneln, panel_c, panel_h, panel_r, panel_r_at, panel_ap, panel_Yh !,   newpanelc, , 

    !=====================================================================
    !                     INITIALIZE VARIABLES
    !=====================================================================

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
      
      DO age=1,MaxAge-1
          requireddeathby_age(age) = - requirednumberby_age(age+1)+requirednumberby_age(age)
          requiredsurvivingby_age(age) = requirednumberby_age(age+1)
      ENDDO
      requireddeathby_age(MaxAge) = requirednumberby_age(MaxAge)
      requiredsurvivingby_age(MaxAge) = 0
      
      !  Construct requirednumberby_age_x;  x is either z, lambda, e 
      !PRINT*, 'DISTRIBUTIONLAR NASIL GORUNUYOR?'   
      DO age=1,MaxAge
               requirednumberby_age_z(age,:)           = NINT( REAL(requirednumberby_age(age),8)*Gz)
               requirednumberby_age_e(age,:)           = NINT(requirednumberby_age(age)*Ge_byage(age,:))
      ENDDO

      ! Make sure that  requirednumberby_age_z,lambda,e are consistent with requirednumberby_age for each age
      ! Here I am adding the excess to the median shock for each age to make sure that for each age we have the right number of people
      ! If adding to the median shock makes negative death for that bin, I add to the other shocks


         DO age=1,MaxAge
         
       !      print*,'age=',age       
             requirednumberby_age_z(age,nz/2+1) = requirednumberby_age_z(age,nz/2+1) &
                          & +  requirednumberby_age(age)-sum(requirednumberby_age_z(age,:))
              requirednumberby_age_e(age,ne/2+1) = requirednumberby_age_e(age,ne/2+1) &
                          & +  requirednumberby_age(age)-sum(requirednumberby_age_e(age,:))          
              if (abs(requirednumberby_age(age) -SUM(requirednumberby_age_z(age,:))) .gt. 0)  then        
                  PRINT*,'required numbers not compatible'
                  PRINT*,requirednumberby_age(age) -SUM(requirednumberby_age_z(age,:))
                  stop
              endif    
          ENDDO

      !    print*,'required number by age z and lambda'
         
          DO age=1,MaxAge
                   DO ee1=1,nz
                       requirednumberby_age_z_lambda(age,ee1,:) = NINT(requirednumberby_age_z(age,ee1)*Glambda)
                        DO ee0=1,nlambda/2+1
                               lambdaindx1 = nlambda/2+1 + ee0 - 1
                               lambdaindx2 = nlambda/2+1 - ee0 + 1 
                               
                               requirednumberby_age_z_lambda(age,ee1,lambdaindx1) =  & 
                                  & requirednumberby_age_z_lambda(age,ee1,lambdaindx1) &
                                  & + requirednumberby_age_z(age,ee1) - sum(requirednumberby_age_z_lambda(age,ee1,:)) 

                               requirednumberby_age_z_lambda(age,ee1,lambdaindx2) =  & 
                                  & requirednumberby_age_z_lambda(age,ee1,lambdaindx2) &
                                  & + requirednumberby_age_z(age,ee1) - sum(requirednumberby_age_z_lambda(age,ee1,:)) 
                                                              
                                  if (age .gt. 1) then
                                       requirednumberby_age_z_lambda(age,ee1,lambdaindx1) = &
                                           & min (requirednumberby_age_z_lambda(age,ee1,lambdaindx1), &
                                           & requirednumberby_age_z_lambda(age-1,ee1,lambdaindx1) )
                                      requirednumberby_age_z_lambda(age,ee1,lambdaindx2) = & 
                                           & min (requirednumberby_age_z_lambda(age,ee1,lambdaindx2), &
                                           & requirednumberby_age_z_lambda(age-1,ee1,lambdaindx2) )
                                  ENDIF     
                        ENDDO  !ee0   lambda   
                        
                        if (abs(requirednumberby_age_z(age,ee1) - sum(requirednumberby_age_z_lambda(age,ee1,:)) ) .gt. 0)  then        
                            PRINT*,'required numbers not compatible'
                            print*, requirednumberby_age_z(age,ee1) - sum(requirednumberby_age_z_lambda(age,ee1,:))   
                            stop
                        endif    
                                                                                 
      !                 DO ee0=1,nlambda
      !                     requirednumberby_age_z_lambda_e(age,ee1,ee0,:) = &
      !                            & NINT(requirednumberby_age_z_lambda(age,ee1,ee0)*Ge_byage(age))
      !                     DO ee2=1,ne
      !                         eindx1 = ne/2 + 1 + ee2 - 1
      !                         eindx2 = ne/2 + 1 - ee2 + 1
      !                        
      !                         requirednumberby_age_z_lambda_e(age,ee1,ee0,eindx1) = &
      !                            & requirednumberby_age_z_lambda_e(age,ee1,ee0,eindx1) + &
      !                            & requirednumberby_age_z_lambda(age,ee1,ee0) - &
      !                            & sum(requirednumberby_age_z_lambda_e(age,ee1,ee0,:))
      !
      !                         requirednumberby_age_z_lambda_e(age,ee1,ee0,eindx2) = &
      !                            & requirednumberby_age_z_lambda_e(age,ee1,ee0,eindx2) + &
      !                            & requirednumberby_age_z_lambda(age,ee1,ee0) - &
      !                            & sum(requirednumberby_age_z_lambda_e(age,ee1,ee0,:))
      !                           
                               ! I do not need to make sure that the number of people in a bin at age should be less than the ones in age+1 since they can switch
                              
      !                     ENDDO ! ee2  e
      !                 ENDDO ! ee0  lambda                 
                    ENDDO  ! ee1  z

      ! The following prints the required number of people 
      !              print*,'age=',age
      !              DO ee1=1,nz
      !                  DO ee0=1,nlambda
      !                        print*,NINT(requirednumberby_age_z_lambda(age,ee1,ee0) * Ge_byage(age,:))
      !                  ENDDO
      !              ENDDO
      !              print*,''
      !              print*,requirednumberby_age_e(age,:)
      !              stop    
         ENDDO ! age

      ! Compute required death and surviving by age, z, and lambda
          DO age=1,MaxAge-1
          
              requireddeathby_age_z_lambda(age,:,:) = - requirednumberby_age_z_lambda(age+1,:,:) & 
                                                                           & + requirednumberby_age_z_lambda(age,:,:)
              requiredsurvivingby_age_z_lambda(age,:,:) = requirednumberby_age_z_lambda(age+1,:,:)

      !        print*,'age=',age
      !        do lambdai=1,nlambda
      !            print*,requirednumberby_age_z_lambda(age,:,lambdai) 
      !        ENDDO    
      !        print*,''
      !        stop
          ENDDO
          requireddeathby_age_z_lambda(MaxAge,:,:) = requirednumberby_age_z_lambda(MaxAge,:,:)
          requiredsurvivingby_age_z_lambda(MaxAge,:,:) = 0
      !        print*,'age=',age
      !        do lambdai=1,nlambda
      !            print*,requirednumberby_age_z_lambda(age,:,lambdai) 
      !        ENDDO    

      !    requireddeathby_age_z_lambda_e(MaxAge,:,:,:) = requirednumberby_age_z_lambda_e(MaxAge,:,:,:)
      !    requiredsurvivingby_age_z_lambd_e(MaxAge,:,:,:) = 0


      IF  (minval(requireddeathby_age_z_lambda) .lt. 0) THEN
          print*,'-----------------------------------------'
          print*,'-----------------------------------------'
          print*,'-----------------------------------------'
          print*,'minval(requireddeathby_age_z_lambda)=',minval(requireddeathby_age_z_lambda)
          print*,'I need to make a fix for age_z dbn similar to what I did for age_z_lambda dbn'
          print*,'-----------------------------------------'
          print*,'-----------------------------------------'
          print*,'-----------------------------------------'
          stop
      ENDIF


      ! E -------------------------------------------   
      ! Compute required death and surviving by age and e

          DO age=1,MaxAge-1
              DO ee1=1,ne
                  requiredsurvivingby_age_e(age,ee1) =  NINT(requirednumberby_age_e(age,ee1) *  survP(age))
              ENDDO
              requiredsurvivingby_age_e(age,ne/2+1) =  requiredsurvivingby_age_e(age,ne/2+1)  &
                      & +  sum( requirednumberby_age_e(age+1,:) - requiredsurvivingby_age_e(age,:) )
              requireddeathby_age_e(age,:) =  requirednumberby_age_e(age,:) -requiredsurvivingby_age_e(age,:)
       
      !       print*,age,'num=',requirednumberby_age_e(age,:),'suv=',requiredsurvivingby_age_e(age,:),'death=',&
      !                & requireddeathby_age_e(age,:)         
          ENDDO !age
          requireddeathby_age_e(MaxAge,:) =  requirednumberby_age_e(MaxAge,:)
          requiredsurvivingby_age_e(MaxAge,:) = 0
      !    

      !
      !print*,'required death by age='
      !DO age=1,MaxAge
      !    print*,requireddeathby_age_e(age,:)
      !ENDDO
      !stop
      !
      !    print*,''
      !    print*,'pr_e(ee1,:)='
      !    DO ee1=1,ne 
      !        print*,pr_e(ee1,:)
      !    ENDDO    
      !    print*,''
      !    
      !    PRINT*,'requirednumberby_e_e='
      !    
          requirednumberby_e_e(:,ne/2+1) = sum(requireddeathby_age_e,1)
          DO ee1=1,ne     
               requirednumberby_e(ee1) = sum(requirednumberby_age_e(:,ee1))        
               requirednumberby_e_e(ee1,:) = NINT( pr_e(ee1,:)*sum(requiredsurvivingby_age_e(:,ee1)) )
      !         PRINT*,requirednumberby_e_e(ee1,:)
         ENDDO

      !print*,'requirednumberby_e=',requirednumberby_e
      !stop

          requirednumberby_e_e(:,ne/2+1) =  requirednumberby_e_e(:,ne/2+1) + sum(requireddeathby_age_e,1)
      !    print*,''
      !    print*,'death correction'
      !    print*,''
      !    DO ee1=1,ne     
      !          PRINT*,requirednumberby_e_e(ee1,:)
      !    ENDDO

      !    print*,''
      !    print*,'horizontal correction'
      !    print*,''
          DO ee1=1,ne     
               requirednumberby_e_e(ee1,ne/2+1) = requirednumberby_e_e(ee1,ne/2+1 ) &
                                                  &  + requirednumberby_e(ee1)  - sum( requirednumberby_e_e(ee1,:) )                            
               PRINT*,requirednumberby_e_e(ee1,:)
         ENDDO   
         print*,''
      !   print*,'requirednumberby_e=',requirednumberby_e   
      !   print*,''
      !   
      !   print*,'requirednumberby_e_e corrected'
      !  print*,''
      ! correction below produces negative number of people
      !   requirednumberby_e_e(ne/2+1,:) = requirednumberby_e_e(ne/2+1,:)  +  &
      !            & sum(requirednumberby_e_e,2) -  sum(requirednumberby_e_e,1)

      ! ALL CORRECTION
          DO ee1=1,ne-1     
                requirednumberby_e_e(ee1+1,ee1) =  requirednumberby_e_e(ee1+1,ee1) &
                    & + sum(requirednumberby_e_e(ee1,:)) -  sum(requirednumberby_e_e(:,ee1))

      !          DO ee0 =1,ne     
      !                 PRINT*,requirednumberby_e_e(ee0,:)
      !           ENDDO 
      !           stop


                requirednumberby_e_e(ee1+1,ee1+1) =     requirednumberby_e_e(ee1+1,ee1+1) &
                    & + requirednumberby_e(ee1+1) - sum(requirednumberby_e_e(ee1+1,:) )
          
      !           DO ee0 =1,ne     
      !                 PRINT*,requirednumberby_e_e(ee0,:)
      !           ENDDO 
      !           stop
      !    
          ENDDO

     
    NBAR_seq = 0.0_DP
    QBAR_seq = 0.0_DP
    !=====================================================================
    !                     GENERATE   INITIAL   PANEL
    !=====================================================================


    newiseed=-1

    numberby_age_z_lambda=0
    numberby_age_e =0

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
    newpaneln      = 0.0_DP
    !=============================================================================
    !
    ! SIMULATE FROM THE SECOND PERIOD SHOCKS AND UPDATE NEW DISTRIBUTIONS
    !
    !=============================================================================

    !call cpu_time(start_timet) 

    DO simutime=1, Max_Simu_Time

      panelage = newpanelage
      panelz = newpanelz
      panele = newpanele
      panellambda = newpanellambda
      panela = newpanela

      !print*,'simutime=',simutime

      numberby_age=0
      deathby_age=0
      survivingby_age =0

      deathby_age_z_lambda = 0
      survivingby_age_z_lambda=0
      numberby_age_z_lambda=0
      numberby_e_e=0

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
                else if (currenta .lt. amin) then
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

    QBAR_seq(simutime) = 0.0_DP 
    DO paneli=1,totpop 
        QBAR_seq(simutime) = QBAR_seq(simutime) +  ( newpanela(paneli)*zgrid(newpanelz(paneli)) )**mu
    ENDDO

    QBAR_seq(simutime) = ( QBAR_seq(simutime) / totpop )**(1.0_DP/mu)  !**(1.0_DP/mu)

    !print*,'simutime=',simutime, 'QBAR_seq(simutime) =',QBAR_seq(simutime)


    if (simutime .gt. (Max_Simu_Time-LastTperiod)) then
        NBAR_seq(simutime) = 0.0_DP
        DO paneli=1,totpop    
               currenta              = newpanela(paneli)
               age                   = newpanelage(paneli)
               currentzi             = newpanelz(paneli)
               currentlambdai        = newpanellambda(paneli) 
               currentei             = newpanele(paneli)
         
              if (age .lt. RetAge) then
        
                    ! do linear interpolation here to find aprime. calling the function takes much more time
                    if (currenta .ge. amax) then
                        tklo = na-1
                        elseif (currenta .lt. amin) then
                            tklo = 1
                            else
                                tklo = ((currenta - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
                    endif 
                    
                    tkhi = tklo + 1        
        
                    newpaneln(paneli) = ((agrid(tkhi) - currenta)*Hours(age,tklo,currentzi,currentlambdai, currentei) &
                                            &  + (currenta - agrid(tklo))*Hours(age,tkhi,currentzi,currentlambdai, currentei)) &
                                            &  / ( agrid(tkhi) - agrid(tklo) )    
                    
                    if (currenta .ge. amax) then
                        newpaneln(paneli) = Hours(age,na,currentzi,currentlambdai, currentei)
                    endif      
                   NBAR_seq(simutime) =  NBAR_seq(simutime) + kappagrid(age) *  lambdagrid(currentlambdai) * &
                                         & egrid(currentei) *newpaneln(paneli)
              endif !age .lt. MaxAge     
                
        ENDDO !paneli
        NBAR_seq(simutime) =  NBAR_seq(simutime) / totpop    
    endif 


    panelage    = newpanelage
    panela      = newpanela
    panelz      = newpanelz
    panellambda = newpanellambda
    panele      = newpanele


    ENDDO ! simutime

    ! Policy functions 
    print*, 'Computing endogenous variables'
    do i=1,totpop
      panel_c(i)  = Linear_Int(agrid, Cons(panelage(i),:,panelz(i),panellambda(i),panele(i)),na, panela(i))
      panel_r(i)  = 1 + ( mu*  rr * (zgrid(panelz(i)))**mu * panela(i)**(mu-1)  - DepRate) ;
      panel_r_at(i) = (1 + ( mu * rr * (zgrid(panelz(i)))**mu * panela(i)**(mu-1)  - DepRate)*(1.0_dp-tauK) )*(1.0_dp-tauW_at) ;
      
      if (panelage(i).lt.RetAge) then 
        panel_h(i)  = max(0.0_DP, &
              & 1.0_DP - (1.0_DP-gamma)*(1.0_dp+tauC)*panel_c(i)/(gamma*psi*wage*eff_un(panelage(i),panellambda(i),panele(i))) )                                
        panel_ap(i) = Y_a(panela(i),zgrid(panelz(i)))+Y_h(panel_h(i),panelage(i),panellambda(i),panele(i),wage)-(1.0_dp+tauC)*panel_c(i)
      else 
        panel_h(i)  = 0.0_DP
        panel_ap(i) = Y_a(panela(i),zgrid(panelz(i))) + RetY_lambda_e(panellambda(i),panele(i))  - (1.0_dp+tauC)*panel_c(i)
      end if 

      if (panel_ap(i).lt.agrid(1)) then
        panel_ap(i) = agrid(1)
        if (panelage(i).lt.RetAge) then                   
          panel_h(i)  = max( 0.0_dp , &
                        & gamma - (1-gamma)*(Y_a(panela(i),zgrid(panelz(i)))-panel_ap(i)) &
                                    & /(psi*wage*eff_un(panelage(i),panellambda(i),panele(i))) ) 
          panel_c(i)  = (Y_a(panela(i),zgrid(panelz(i))) + Y_h(panel_h(i),panelage(i),panellambda(i),panele(i),wage)  - panel_ap(i))/(1.0_dp+tauC)
        else 
          panel_h(i) = 0.0_DP
          panel_c(i) = (Y_a(panela(i),zgrid(panelz(i))) + RetY_lambda_e(panellambda(i),panele(i))  - panel_ap(i))/(1.0_dp+tauC)
        end if 
      end if

      if (panelage(i).lt.RetAge) then  
        panel_Yh(i) =  Y_h(panel_h(i),panelage(i),panellambda(i),panele(i),wage)
      else 
        panel_Yh(i) =  RetY_lambda_e(panellambda(i),panele(i))
      end if 
    end do



    ! Save results
    print*, 'Saving Results'
    if (ben_switch.eq.1) then 
      OPEN  (UNIT=1,  FILE=trim(Result_Folder)//'Simul/Sim_age_ben' , STATUS='replace') 
      OPEN  (UNIT=2,  FILE=trim(Result_Folder)//'Simul/Sim_A_ben'   , STATUS='replace')  
      OPEN  (UNIT=3,  FILE=trim(Result_Folder)//'Simul/Sim_Z_ben'   , STATUS='replace')  
      OPEN  (UNIT=4,  FILE=trim(Result_Folder)//'Simul/Sim_C_ben'   , STATUS='replace')  
      OPEN  (UNIT=5,  FILE=trim(Result_Folder)//'Simul/Sim_H_ben'   , STATUS='replace') 
      OPEN  (UNIT=6,  FILE=trim(Result_Folder)//'Simul/Sim_R_ben'   , STATUS='replace')  
      OPEN  (UNIT=7,  FILE=trim(Result_Folder)//'Simul/Sim_R_at_ben', STATUS='replace') 
      OPEN  (UNIT=8,  FILE=trim(Result_Folder)//'Simul/Sim_Ap_ben'  , STATUS='replace') 
      OPEN  (UNIT=9,  FILE=trim(Result_Folder)//'Simul/Sim_Yh_ben'  , STATUS='replace') 
    else
      OPEN  (UNIT=1,  FILE=trim(Result_Folder)//'Simul/Sim_age_exp' , STATUS='replace') 
      OPEN  (UNIT=2,  FILE=trim(Result_Folder)//'Simul/Sim_A_exp'   , STATUS='replace')  
      OPEN  (UNIT=3,  FILE=trim(Result_Folder)//'Simul/Sim_Z_exp'   , STATUS='replace')  
      OPEN  (UNIT=4,  FILE=trim(Result_Folder)//'Simul/Sim_C_exp'   , STATUS='replace')  
      OPEN  (UNIT=5,  FILE=trim(Result_Folder)//'Simul/Sim_H_exp'   , STATUS='replace')
      OPEN  (UNIT=6,  FILE=trim(Result_Folder)//'Simul/Sim_R_exp'   , STATUS='replace')  
      OPEN  (UNIT=7,  FILE=trim(Result_Folder)//'Simul/Sim_R_at_exp', STATUS='replace') 
      OPEN  (UNIT=8,  FILE=trim(Result_Folder)//'Simul/Sim_Ap_exp'  , STATUS='replace') 
      OPEN  (UNIT=9,  FILE=trim(Result_Folder)//'Simul/Sim_Yh_exp'  , STATUS='replace')    
    end if 

    WRITE (UNIT=1,  FMT=*) panelage
    WRITE (UNIT=2,  FMT=*) panela 
    WRITE (UNIT=3,  FMT=*) panelz
    WRITE (UNIT=4,  FMT=*) panel_c 
    WRITE (UNIT=5,  FMT=*) panel_h
    WRITE (UNIT=6,  FMT=*) panel_r
    WRITE (UNIT=7,  FMT=*) panel_r_at
    WRITE (UNIT=8,  FMT=*) panel_ap
    WRITE (UNIT=9,  FMT=*) panel_Yh

    CLOSE (unit=1)
    CLOSE (unit=2)
    CLOSE (unit=3)
    CLOSE (unit=4)
    CLOSE (unit=5)
    CLOSE (unit=6)
    CLOSE (unit=7)
    CLOSE (unit=8)
    CLOSE (unit=9)

    !call cpu_time(finish_timet)
    !!print '("SIMULATION TIME = ",f6.3," seconds.")',finish_timet-start_timet
    !print*,'SIMULATION TIME =', finish_timet-start_timet
    !
    !stop

  END SUBROUTINE SIMULATION

  !====================================================================

end Program GKK_Simul
