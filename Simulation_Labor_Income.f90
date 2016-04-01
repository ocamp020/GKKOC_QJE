Program Simulation_Labor_Income
	USE parameters
	USE GLOBAL
	USE omp_lib
	
	Implicit None
	real(dp), dimension(totpop) :: Y_L, panel_a
	integer , dimension(totpop) :: panel_age, panel_z, panel_lambda, panel_e, panel_x
	real(dp) :: h_i
	integer  :: tklo, tkhi, i
	character(100) :: Results_Folder, Bench_Folder, Simul_Folder

	!$ call omp_set_num_threads(20)

	print*, ' '
	print*, 'Simulation of labor income'
	print*, ' '

	! Set Folders
	print*, 'Setting up folders'
	Results_Folder = './NSU_ZS_LT_Results/Theta_1.50/Factor_0.00/Exp_Shock_mu90/'
	Bench_Folder  = './NSU_ZS_LT_Results/Theta_1.50/Factor_0.00/Exp_Shock_mu90/Bench_Files/'
	Simul_Folder  = './NSU_ZS_LT_Results/Theta_1.50/Factor_0.00/Exp_Shock_mu90/Simul/'
	print*, 'Result_Folder', Results_Folder
	print*, 'Bench_Folder' , Bench_Folder
	print*, 'Simul_Folder' , Simul_Folder

	! Load Variables
	print*, ' '
	print*, 'Loading Variables'
		! Wage
		OPEN (UNIT=1,  FILE=trim(Bench_Folder)//'wage'  , STATUS='old', ACTION='read')
		READ (UNIT=1,  FMT=*), wage
		CLOSE(UNIT=1)
		! psi
		OPEN (UNIT=1,  FILE=trim(Bench_Folder)//'psi'   , STATUS='old', ACTION='read')
		READ (UNIT=1,  FMT=*), psi
		CLOSE(UNIT=1)
		! tauPL
		OPEN (UNIT=1,  FILE=trim(Bench_Folder)//'tauPL' , STATUS='old', ACTION='read')
		READ (UNIT=1,  FMT=*), tauPL
		CLOSE(UNIT=1)
		! Hours
		OPEN (UNIT=1,  FILE=trim(Bench_Folder)//'hours'  , STATUS='old', ACTION='read')
		READ (UNIT=1,  FMT=*), Hours
		CLOSE(UNIT=1)
		! agrid
		OPEN (UNIT=1,  FILE=trim(Results_Folder)//'agrid'  , STATUS='old', ACTION='read')
		READ (UNIT=1,  FMT=*), agrid
		CLOSE(UNIT=1)
		! eff_un
		OPEN (UNIT=1,  FILE=trim(Results_Folder)//'eff_un'  , STATUS='old', ACTION='read')
		READ (UNIT=1,  FMT=*), eff_un
		CLOSE(UNIT=1)
		! RetY
		OPEN (UNIT=1,  FILE=trim(Results_Folder)//'Ret_Y'  , STATUS='old', ACTION='read')
		READ (UNIT=1,  FMT=*), RetY_lambda_e
		CLOSE(UNIT=1)
		! panel_age
		OPEN (UNIT=1,  FILE=trim(Simul_Folder)//'panelage_bench'  , STATUS='old', ACTION='read')
		READ (UNIT=1,  FMT=*), panel_age
		CLOSE(UNIT=1)
		! panel_a
		OPEN (UNIT=1,  FILE=trim(Simul_Folder)//'panela_bench'  , STATUS='old', ACTION='read')
		READ (UNIT=1,  FMT=*), panel_a
		CLOSE(UNIT=1)
		! panel_z
		OPEN (UNIT=1,  FILE=trim(Simul_Folder)//'panelz_bench'  , STATUS='old', ACTION='read')
		READ (UNIT=1,  FMT=*), panel_z
		CLOSE(UNIT=1)
		! panel_lambda
		OPEN (UNIT=1,  FILE=trim(Simul_Folder)//'panellambda_bench'  , STATUS='old', ACTION='read')
		READ (UNIT=1,  FMT=*), panel_lambda
		CLOSE(UNIT=1)
		! panel_e
		OPEN (UNIT=1,  FILE=trim(Simul_Folder)//'panele_bench'  , STATUS='old', ACTION='read')
		READ (UNIT=1,  FMT=*), panel_e
		CLOSE(UNIT=1)
		! panel_x
		OPEN (UNIT=1,  FILE=trim(Simul_Folder)//'panelx_bench'  , STATUS='old', ACTION='read')
		READ (UNIT=1,  FMT=*), panel_x
		CLOSE(UNIT=1)

	! Simulation
	print*, ' '
	print*, 'Simulating'
	!$omp parallel do private(tklo,tkhi,h_i)
	do i=1,totpop
		if (panel_age(i).lt.RetAge) then 
			if (panel_a(i) .ge. amax) then
	            tklo = na-1
	        else if (panel_a(i) .lt. amin) then
	            tklo = 1
	        else
	            tklo = ((panel_a(i) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
	        endif 
	        
	        tkhi = tklo + 1        

	        h_i  = ((agrid(tkhi) - panel_a(i))*hours(panel_age(i),tklo,panel_z(i),panel_lambda(i),panel_e(i),panel_x(i)) &
	           &  + (panel_a(i) - agrid(tklo))*hours(panel_age(i),tkhi,panel_z(i),panel_lambda(i),panel_e(i),panel_x(i))) &
	                                &  / ( agrid(tkhi) - agrid(tklo) )  

			Y_L(i) = psi*( Wage*eff_un(panel_age(i),panel_lambda(i),panel_e(i))*h_i)**(1.0_dp-tauPL)
		else 
			Y_L(i) = RetY_lambda_e(panel_lambda(i),panel_e(i))
		endif 
 		!$omp critical
		print*, i, Y_L(i)
		!$omp end critical
	enddo

	! Write Results 
	print*, ' '
	print*, 'Writing Results'
	OPEN  (UNIT=1,  FILE=trim(Simul_Folder)//'Y_L_bench'  , STATUS='replace')
	WRITE (UNIT=1,  FMT=*) Y_L
	CLOSE (unit=1)


end Program Simulation_Labor_Income