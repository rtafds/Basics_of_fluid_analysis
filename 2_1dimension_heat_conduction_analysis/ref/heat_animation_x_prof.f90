module constants
  implicit none
  integer, parameter :: m=61, nmax=5000
  integer, parameter :: SP = kind(1.0)
  integer, parameter :: DP = selected_real_kind(2*precision(1.0_SP))
end module constants
module print
  use constants
  implicit none
contains
  subroutine print__profile(f)
    real(DP), dimension(0:m+1), intent(in) :: f
    integer :: i
    integer :: counter = 0  ! automatic save attribute
    character(len=*), parameter :: base = "./data/temp.j=middle."
    character(len=4) :: serial_num
    write(serial_num,'(i4.4)') counter
    open(10,file=base//serial_num)
    do i = 0 , m+1
       write(10,*) i, f(i)
    end do
    close(10)
    counter = counter + 1
  end subroutine print__profile
end module print
program heat1_animation_x_prof
  use constants
  use print
  implicit none
  integer :: i,j,n
  real(DP), dimension(:,:), allocatable :: u, un
  real(DP) :: h, heat=1.0_DP
  real(DP) :: random
  allocate(u(0:m+1,0:m+1))
  allocate(un(m,m))
  h = 1.0_DP/(m+1)
  u(:,:) = 0.0_DP
  call print__profile(u(:,m/2+1))
  do n=1, nmax
    do j=1, m
      do i=1, m
        un(i,j)=(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1))*0.25_DP+heat*h*h
      end do
    end do
    u(1:m,1:m)=un(1:m,1:m)
    if (mod(n,100)==0) call print__profile(u(:,m/2+1))
  end do
  deallocate(u,un)
end program heat1_animation_x_prof