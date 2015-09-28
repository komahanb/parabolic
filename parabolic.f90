program heat_equation_one_dim
  ! solve a u_xx = u_t for 0 <= x <= xf, 0 <= t <= T; u = u(x,t)
  ! u(x, 0) = sin πx,
  ! u(0, t) = 0,
  ! u(1, t) = 0

  implicit none
  real*8, parameter  :: pi = dble(22)/dble(7)
  integer, parameter :: M = 20, N = 100
  real*8  :: dx, dt
  real*8  :: xf, T
  real*8  :: u(M+1, N+1) ! u(spatial, time)
  real*8  :: x(M+1)
  real*8  :: r, a, r1, norm
  integer :: i,j,k

  print*, "Parabolic ODE solver using Explicit Euler"

  a  = 1.0d0
  xf = 1.0d0
  T  = 1.0d-1

  dx = xf/dble(M) ! spatial stepping
  dt = T/dble(N)  ! time stepping

  ! generate the spatial grid
  do i = 1, M+1
     x(i) = dble(i-1)*dx
  end do

  ! initial conditions
  do i = 1, M+1
     u(i,1) = sin(pi*x(i))
  end do

  ! boundary conditions
  ! left end
  do j = 1, N+1
     u(1, j)   = 0
  end do
  ! right end
  do j = 1, N+1
     u(M+1, j) = 0
  end do
 
  ! Explicit Euler

  r = A*dt/dx**2
  if (r.gt.0.5) stop"Explicit Euler will be unstable... change parameters"
  r1= 1.0d0 - 2.0d0*r

  do k = 2, N + 1 ! ignore the initial condition (marching until final time)
     do j = 2, M  ! ignore left and right ends
        u(j, k) = r*(u(j+1, k-1) + u(j-1, k-1)) + r1*u(j, k-1)
     end do
     write(*,'(F15.3,E15.3)'), k*dt, norm2(u(:,k-1) - u(:,k))/sqrt(dble(M+1))
  end do 

  print*, "Writing solution to heat_soln.dat"

  ! write the solution
  open (unit = 2, file = "heat_soln.dat")
  write(2,*) 'variables= "x","t","u"'
  write(2,*) 'zone i=', M+1,' j=',N+1,' datapacking=point'
  do i = 1, N+1
     do j =1, M+1
        write (2,*) dble(j-1)*dx, dble(i-1)*dt, u(j,i)
     end do
  end do
  close(2)

end program heat_equation_one_dim

