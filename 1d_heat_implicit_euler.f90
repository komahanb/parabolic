program heat_equation_one_dim_implicit_euler
  ! solve a u_xx = u_t for 0 <= x <= xf, 0 <= t <= T; u = u(x,t)
  ! u(x, 0) = sin πx,
  ! u(0, t) = 0.15,
  ! u(1, t) = 0.15

  implicit none
  real*8, parameter  :: pi = dble(22)/dble(7)
  integer, parameter :: M = 10, N = 100

  real*8  :: dx, dt
  real*8  :: xf, T

  ! solution
  real*8  :: u(M+1, N+1) =0.0d0! u(spatial, time)

  ! vandermonte system and its inverse
  real*8  :: V(M+1, M+1) =0.0d0, rhs(M+1)
  real*8  :: inverse_V(M+1, M+1) =0.0d0! u(spatial, time)

  ! spatial grid
  real*8  :: x(M+1)
  real*8  :: r, a, r1,d
  integer :: i, j, k, errorflag

  print*, "Parabolic ODE solver using Implicit Euler"

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
     u(1, j)   = 0.15d0
  end do
  ! right end
  do j = 1, N+1
     u(M+1, j) = 0.15d0
  end do
 
  ! Implicit Euler setup vandermonte system
  r = A*dt/dx**2
  d=1.0d0+2.0d0*r
  do i = 1, M +1
     do j = 1, M +1
        if (j.eq.i) then
           V(j,i) = d
        else if (j.eq.(i-1)) then
           V(j,i)= -r
        else if (j.eq.(i+1)) then
           V(j,i)= -r
        end if
     end do
  end do
  
  ! change the system to fit the boundary conditions
  V(1,1)     = 1.0d0
  V(M+1,M+1) = 1.0d0
  V(1,2)     = 0.0d0
  V(M+1,M) = 0.0d0

  ! find the inverse of vandermonte matrix
  CALL FINDInv(V,Inverse_V,m+1,errorflag)
  if (errorflag .eq. -1) stop" Computing Inverse failed"

  do k = 2, N+1
     U(: , k) = matmul(Inverse_V,U(:,k-1))
     write(*,'(F15.3,E15.3)'), (k-1)*dt, norm2(u(:,k-1) - u(:,k))/sqrt(dble(M+1))
  end do

  print*, "Writing solution to heat_soln.dat"
  
  ! write the solution
  open (unit = 2, file = "results/heat_soln.dat")
  write(2,*) 'variables= "x","t","u"'
  write(2,*) 'zone i=', M+1,' j=',N+1,' datapacking=point'
  do i = 1, N+1
     do j =1, M+1
        write (2,*) dble(j-1)*dx, dble(i-1)*dt, u(j,i)
     end do
  end do
  close(2)

end program heat_equation_one_dim_implicit_euler
