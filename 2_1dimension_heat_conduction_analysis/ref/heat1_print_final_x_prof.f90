program heat1_print_final_x_prof
  implicit none
  integer, parameter :: m=31, nmax=5000
  integer :: i,j,n
  integer, parameter :: SP = kind(1.0)
  integer, parameter :: DP = selected_real_kind(2*precision(1.0_SP))
  real(DP), dimension(:,:), allocatable :: u, un
  real(DP) :: h, heat=1.0_DP, heat_hsq
  allocate(u(0:m+1,0:m+1))
  allocate(un(m,m))
  h = 1.0_DP/(m+1)
  heat_hsq = heat*h*h
  u(:,:) = 0.0_DP
  open(10,file='./data/temp.final_profile_x')
  do n = 1 , nmax
    do j = 1 , m
      do i = 1 , m
        un(i,j)=(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1))*0.25_DP+heat_hsq
      end do
    end do
    u(1:m,1:m)=un(1:m,1:m)
  end do
  do i = 0 , m+1
     write(10,*) i, u(i,m/2+1)
  end do
  deallocate(u,un)
  close(10)
end program heat1_print_final_x_prof