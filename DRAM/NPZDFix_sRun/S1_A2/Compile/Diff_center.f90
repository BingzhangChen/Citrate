!subroutine diff_center(N,dt,cnpar,posconc,h,Bcup,Bcdw, &
!                       Yup,Ydw,nuY,Lsour,Qsour,Taur,Yobs,Yin,Yout)
subroutine diff_center(N,dt,cnpar,posconc,h,nuY,Yin,Yout)

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!   integer,parameter                    :: Dirichlet      = 0
!   integer,parameter                    :: Neumann        = 1

!  number of vertical layers
   integer,  intent(in)               :: N

!  time step (s)
   real,    intent(in)                :: dt

!  "implicitness" parameter
   real,    intent(in)                :: cnpar

!  1: non-negative concentration, 0: else
   integer, intent(in)                :: posconc

!  layer thickness (m)
   real,    intent(in)                :: h(N)

!  type of upper BC
!   integer,  intent(in)               :: Bcup

!  type of lower BC
!   integer,  intent(in)               :: Bcdw

!  value of upper BC
!   real(8), intent(in)                :: Yup

!  value of lower BC
!   real(8), intent(in)                :: Ydw

!  diffusivity of Y
   real,    intent(in)                :: nuY(0:N)

!  linear source term
!  (treated implicitly)
!   real(8), intent(in)                :: Lsour(N)

!  constant source term
!  (treated explicitly)
!   real(8), intent(in)                :: Qsour(N)

!  relaxation time (s)
!   real(8), intent(in)                :: Taur(N)

!  observed value of Y
!   real(8), intent(in)                :: Yobs(N)
!
! !INPUT/OUTPUT PARAMETERS:
   real(8), intent(in)                :: Yin(N)
   real(8), intent(out)               :: Yout(N)

   real(8), dimension(N)              :: au,bu,cu,du
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
   real(8)                   :: a,c,l
!
!-----------------------------------------------------------------------
!BOC
!
!   Debug:
!   write(6,*) 'Yup    = ',Yup
!   write(6,*) 'Ydw    = ',Ydw

!  Initialize au, bu, cu, du
   au(:) = 0d0; bu(:) = 0d0; cu(:) = 0d0; du(:) = 0d0

!  set up matrix
   do i=2,N-1
      c     = 2.0*dt*nuY(i)  /(h(i)+h(i+1))/h(i)
      a     = 2.0*dt*nuY(i-1)/(h(i)+h(i-1))/h(i)
!      l     =     dt*Lsour(i)

      cu(i) =-cnpar*c
      au(i) =-cnpar*a
     ! bu(i) = 1d0 + cnpar*(a + c) - l
      bu(i) = 1d0 + cnpar*(a + c)

     ! du(i) = (1d0 - (1d0-cnpar)*(a + c))*Yin(i)                  &
     !       + (1d0 - cnpar)*( a*Yin(i-1) + c*Yin(i+1) ) + dt*Qsour(i)

      du(i) = (1d0 - (1d0-cnpar)*(a + c))*Yin(i)                  &
            + (1d0 - cnpar)*( a*Yin(i-1) + c*Yin(i+1) )

   enddo

!   set up upper boundary condition
!   select case(Bcup)
!   case(Neumann)
      a     = 2.0*dt*nuY(N-1)/(h(N)+h(N-1))/h(N)
!      l     = dt*Lsour(N)

      au(N) =-cnpar*a
!      if (posconc .eq. 1 .and. Yup.lt.0d0) then ! Patankar (1980) trick
!         bu(N) =  1d0 - au(N) - l  - dt*Yup/Yin(N)/h(N)
!         du(N) = Yin(N) + dt*Qsour(N)   &
!               + (1d0 - cnpar)*a*(Yin(N-1)-Yin(N))
!      else
!       bu(N) =  1d0 - au(N) - l
       bu(N) = 1d0 - au(N)
!       du(N) = Yin(N) + dt*(Qsour(N)+Yup/h(N))   &
!             + (1d0 - cnpar)*a*(Yin(N-1)-Yin(N))
       du(N) = Yin(N) + (1d0 - cnpar)*a*(Yin(N-1)-Yin(N))
!      end if
!   case(Dirichlet)
!      au(N) = 0d0
!      bu(N) = 1d0
!      du(N) = Yup
!   case default
!      write(6,*) 'Fatal error: invalid boundary condition type for upper boundary'
!      stop  'diff_center.F90'
!   end select

!   set up lower boundary condition
!   select case(Bcdw)
!   case(Neumann)
      c     = 2.0*dt*nuY(1)/(h(1)+h(2))/h(1)
!      l     = dt*Lsour(1)

      cu(1) =-cnpar*c
     ! if (posconc.eq.1 .and. Ydw.lt.0d0) then ! Patankar (1980) trick
     !    bu(1) = 1d0 - cu(1) - l - dt*Ydw/Yin(1)/h(1)
     !    du(1) = Yin(1) + dt*(Qsour(1))   &
     !          + (1d0 - cnpar)*c*(Yin(2)-Yin(1))
     ! else
     !    bu(1) = 1d0 - cu(1) - l
     !    du(1) = Yin(1) + dt*(Qsour(1)+Ydw/h(1))   &
     !          + (1d0 - cnpar)*c*(Yin(2)-Yin(1))

         bu(1) = 1d0 - cu(1)
         du(1) = Yin(1) + (1d0 - cnpar)*c*(Yin(2)-Yin(1))

     ! endif
!   case(Dirichlet)
!      cu(1) = 0d0
!      bu(1) = 1d0
!      du(1) = Ydw
!   case default
!      write(6,*) 'Fatal error: invalid boundary condition type for lower boundary'
!      stop  'diff_center.F90'
!   end select


!  relaxation to observed value
!   if (minval(Taur).lt.1.E10) then
!      do i=1,N
!         bu(i)=bu(i)+dt/Taur(i)
!         du(i)=du(i)+dt/Taur(i)*Yobs(i)
!      end do
!   end if


!  solve linear system
   call tridiagonal(N,au,bu,cu,du,1,N,Yout)

   return
end subroutine diff_center
!-----------------------------------------------------------------------
