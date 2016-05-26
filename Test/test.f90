module my_subs

contains

Subroutine select(x,y) !result(y)
    implicit none
    integer, dimension (:), intent (in) :: x
    integer, dimension (:), allocatable :: y
    integer :: i, j

    j = 0
    do i=1, size (x)
        if (x(i)/=0) j = j+1
    enddo

    allocate ( y (1:j) )

    j = 0
    do i=1, size (x)
        if (x(i)/=0) then
            j = j+1
            y(j) = x(i)
        endif
    enddo

    return

end Subroutine select



	Subroutine Sort(n,A,A_sort,Sort_ind)
		integer , intent(in) :: n    !Number of elements in A
		real, intent(in) , dimension(n) :: A
		real, intent(out), dimension(n) :: A_sort
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

end module my_subs

program test

use my_subs
!use, intrinsic :: iso_fortran_env
use, intrinsic :: ieee_arithmetic
implicit none
integer, dimension (6) :: array = [ 5, 0, 3, 0, 6, 1 ]
integer, dimension (:), allocatable :: answer
real, dimension(4) :: aa,bb,dd
integer, dimension(4) :: bb_ind
real, dimension(2,3,4) :: matrix, matrix_in
real, dimension(2,3,4,5,6) :: dbn
real, dimension(3,4)       :: dbn_1, dbn_2, a_mat, z_mat
real, dimension(3,3)       :: zz_mat
real, dimension(3)         :: agrid
real, dimension(4)         :: zgrid
integer :: i,j,k
character(100) :: mydir 
!real :: aaa = 3.25

real, allocatable  :: zz(:)
real, dimension(15) :: yy
integer, dimension(15) :: xx
real, dimension(2,2) :: aaa 
real, dimension(4)   :: bbb
real(8), dimension(2) :: ccc
logical :: log_1, log_2(2)

type mytype
   real(8), dimension(2,3,4,5,6,7) :: a
end type mytype

type(mytype), dimension(2,2,2) :: kk 
real(8), dimension(2,2,2,2,3,4,5,6,7) :: kkk

! !answer = select (array)
! call select(array,answer)

! write (*, *) size (array), size (answer)
! write (*, *) array
! write (*, *) answer

! !aa = 0
! !aa = 0/aa 
! aa = HUGE(1.0)
! !aa = IEEE_VALUE(aa,IEEE_QUIET_NAN)
! aa(4) = 1.0
! !aa(4) = 2.98
! aa(1) = 5.1
! write(*,*) aa
! call sort(4,aa,bb,bb_ind)

! dd = aa(bb_ind)


! write(*,*) aa, bb_ind
! write(*,*) dd

! matrix(:,:,1) = 1
! matrix(:,:,2) = 2
! matrix(:,:,3) = 3
! matrix(:,:,4) = 4

! OPEN (UNIT=2, FILE='test_matrix', STATUS='replace')
! WRITE (UNIT=2, FMT=*) matrix
! CLOSE (unit=2)
! OPEN (UNIT=2, FILE='test_matrix', STATUS='old', ACTION='read')
! READ(Unit=2,FMT=*), matrix_in

! print*,"print matrix"
! do i=1,4
! 	print*,matrix_in(1,:,i)
! 	print*,matrix_in(2,:,i)
! end do 

! dbn(:,1,:,:,:)= 1.0
! dbn(:,2,:,:,:)= 2.0
! dbn(:,3,:,:,:)= 3.0

! dbn_1 = sum(sum(sum(dbn,5),4),1) 
! do i=1,3
! 	do j=1,4
! 		dbn_2(i,j) = sum(dbn(:,i,j,:,:))
! 	end do 
! end do 

! do i=1,3 
! 	print*,"DBN_1", dbn_1(i,:), "DBN_2", dbn_2(i,:)
! end do 

! agrid=[1,2,3]
! a_mat= spread(agrid,2,4)
! zgrid=[3,4,5,6]
! z_mat= spread(zgrid,1,3)

! do i=1,3
! 	print*, "aa",a_mat(i,:),"zz",z_mat(i,:)
! end do 

! print*, shape( sum(sum(sum(dbn,5),4),1)  )
! print*, shape( spread(agrid,2,4) )
! print*, shape( spread(zgrid,1,3) )


! print*, " "
! print*, " "
! print*, "Make new directory"
! print*, aaa
! write(mydir,'(f4.2)') aaa
! print*, mydir
! print*, trim(mydir)
! mydir = 'Sub_Test_'//trim(mydir)
! call execute_command_line( 'mkdir -p ' // trim(mydir) )


xx = 0
xx(1)  = 1
xx(4)  = 1
xx(8)  = 1
xx(11) = 1
xx(14) = 1

yy = 3.1416
where(xx.eq.1) yy=0.11
yy(4) = 0.99
yy(14) = 3.1111

zz = pack(yy,(xx.eq.1))

print *, zz

aaa(1,1)=22.1
aaa(1,2)=55.1
aaa(2,1)=99.0
aaa(2,2)=34.22

bbb=reshape(aaa,(/size(aaa)/))
print*, bbb
print*, size(aaa)
print*, bbb( (/4,1/)  )
print*, min(bbb,(/30.0,1.0,1.0,10000.0/))


z_mat(:,1) = 1.0
z_mat(:,2) = 2.0
z_mat(:,3) = 3.0
z_mat(:,4) = 4.0

zz_mat(:,1) = sum(z_mat(:,1:2),2)
zz_mat(:,2:3) = z_mat(:,3:4)

print*,zz_mat(1,:)
print*,zz_mat(2,:)
print*,zz_mat(3,:)

	
print*, ' '
print*, (1e-9.lt.1e-8)
	ccc(1)=0.0
	ccc(2)=0.3
print*, any(ccc/ccc*(/1e-9,1e-7/).lt.1e-8)
log_1 = any(ccc/ccc*(/1e-9,1e-7/).lt.1e-8)
log_2 =    (ccc/ccc*(/1e-9,1e-7/).lt.1e-8)
if (log_1) then 
	print*, 'Test 1'
	endif 
if (any(log_2)) then 
	print*, 'Test 2'
	endif 

print*, ' '
print*, 'Test Count'
print*, xx
print*, count(xx.eq.0), count(xx.eq.1), count(xx.gt.0), count(xx.gt.2)

print*, ' '
print*, z_mat(1,:)
print*, z_mat(1,(/2,1,2/))


do i=1,2; do j=1,2; do k=1,2;
kk(i,j,k)%a = 0.0
enddo; enddo; enddo

kk(1,1,1)%a(:,:,:,:,:,:)=10
kk(1,1,2)%a(:,:,:,:,:,:)=11
kk(1,2,1)%a(:,:,:,:,:,:)=12
kk(1,2,2)%a(:,:,:,:,:,:)=13

print*, ' '
print*,'test kk'
print*, kk



end program test