Program Test_Rowenhurts 
	Implicit None 
	integer, parameter          ::  DP      = KIND(1.0D0)
    integer, parameter          ::  QP      = selected_real_kind(32)
    integer, parameter          ::  LG      = KIND(.true.)
    real(dp), parameter         ::  eps_p   = EPSILON(1.0_dp)
    real(dp), parameter         ::  small_p = 1E-7_dp
    real(dp), parameter         ::  big_p   = HUGE(1.0_dp)
	
	integer, parameter  :: n_z=5
	real(dp) :: z_grid(n_z), P_z(n_z,n_z), rho, sigma
	integer :: i

	rho = 0.99_dp
	sigma = 0.04_dp
	call MarkovAR_95(n_z,rho,sigma,z_grid,P_z)

	print*,"z_grid"
	print*,z_grid
	print*," "
	print*,"P_z"
	do i=1,n_z
		print*,P_z(i,:)
	end do



Contains 

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
	        p = (1+rho)/2
	        q = (1+rho)/2
	        psi = sqrt(real(n_z)-1)*sigma/sqrt(1-rho**2)

	    ! Step of grid
	        step = 2*psi/(n_z-1)
	    
	    ! Compute z_grid
	        z_grid(1) = -psi
	        do i=2,n_z
	            z_grid(i) = z_grid(i-1) + step
	        end do
	    
	    ! Compute transition matrix for n_z=2
	        P_2 = reshape((/p,1-q,1-p,q/),(/2,2/))

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
		        P_b = ( p*P_a(2:n_z+1,2:n_z+1) + (1-p)*P_a(2:n_z+1,1:n_z) + &
		              & (1-q)*P_a(1:n_z,2:n_z+1) + q*P_a(1:n_z,1:n_z) )
		        P_half(1,:)       =  1.0_dp
		        P_half(2:n_z-1,:) = 0.50_dp
		        P_half(n_z,:)     = 1.0_dp
		        P_z = P_b*P_half
	        else
	            P_z = P_2
	        end if

	    return
	end subroutine  MarkovAR_95

end Program Test_Rowenhurts