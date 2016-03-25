! Calibration Loop

Program Calibration_Loop
	use nrtype
	Implicit None
	real(DP) :: beta, beta_L, beta_H
	real(DP) :: sigmaz, sigmaz_L, sigmaz_H
	real(DP) :: x_hi, x_hi_L, x_hi_H
	real(DP) :: a_x, a_x_L, a_x_H
	real(DP) :: b_x, b_x_L, b_x_H
	real(DP) :: par(5)
	integer  :: n_beta, n_sigmaz, n_x_hi, n_a, n_b, i, j, k, l, m, ind
	character(100) :: Result_Folder, log_file

	Result_Folder = './Calibration_Loop/'
	call system( 'mkdir -p ' // trim(Result_Folder) )

	OPEN(UNIT=3, FILE=trim(Result_Folder)//'Calibration_Loop_Restuls.txt', STATUS='replace')
	WRITE(unit=3, FMT=*) ' '
	WRITE(unit=3, FMT=*) 'This file contains the results of the benchmark model with different calibrations'
	WRITE(unit=3, FMT=*) ' '
	WRITE(unit=3, FMT=*) 'beta ','sigmaz ','x_hi ','a_x ','b_x ','W/GDP ','STD_Earnings ','Mean_Labor ','MeanReturn ', &
			& 'PV_Top_1% ','PV_Top_10% '
	CLOSE(unit=3)

	call system( 'make GKK_Calibration_Loop.a' )

	beta_L   = 0.95_dp
	beta_H   = 0.95_dp
	n_beta   = 1

	sigmaz_L = 0.3350_dp
	sigmaz_H = 0.35_dp
	n_sigmaz = 1

	x_hi_L   = 4.0_dp 
	x_hi_H   = 5.0_dp
	n_x_hi   = 1

	a_x_L    = 0.05_dp
	a_x_H    = 0.20_dp
	n_a      = 5

	b_x_L    = 0.00_dp
	b_x_H    = 1.00_dp
	n_b      = 5

	ind =1
	do i=1,n_beta
		beta   = beta_L   + real(i-1,8)*(beta_H-beta_L)/real(max((n_beta-1),1),8)
	do j=1,n_sigmaz
		sigmaz = sigmaz_L + real(j-1,8)*(sigmaz_H-sigmaz_L)/real(max((n_sigmaz-1),1),8)
	do k=1,n_x_hi
		x_hi   = x_hi_L   + real(k-1,8)*(x_hi_H-x_hi_L)/real(max((n_x_hi-1),1),8)
	do l=1,n_a
		a_x    = a_x_L    + real(l-1,8)*(a_x_H-a_x_L)/real(max((n_a-1),1),8)
	do m=1,n_b
		b_x    = b_x_L    + real(l-1,8)*(b_x_H-b_x_L)/real(max((n_b-1),1),8)

		! Set Log File 
			write(log_file,'(I3.1)') ind 
			log_file = 'log_ind_'//trim(log_file)//'.txt'
			ind = ind + 1

		! Write parameters into file
			par = (/beta,sigmaz,x_hi,a_x,b_x/)
			OPEN(UNIT=3, FILE=trim(Result_Folder)//'Loop_Par', STATUS='replace')
			WRITE(unit=3, FMT=*) par
			CLOSE(unit=3)

		! Call Main Program
			call system( 'nohup Compiled_Files/GKK_Calibration_Loop.a >' // trim(log_file) //' &' )
			call system( 'sleep 15' )

	enddo 
	enddo
	enddo 
	enddo
	enddo 


end Program Calibration_Loop
