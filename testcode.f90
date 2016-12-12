module general
implicit none
integer :: step, i, j, p
integer, parameter :: a = 10, npart = 100, nstep = 1000
real(8), dimension(npart, 2) :: pos
real(8), parameter :: R = 0.1, dt=0.01, m = 0.00001, U0 = 1, Pi = atan(1.)*4
end module general


program deuxD
use general
implicit none
real(8) :: v1, v2, r0, r1, r2, a0, a1, a2, x, y, t, theta
real(8) :: LennardJones, gauss
real(8), dimension(npart - 1, npart - 1) :: dx, dy


!Verlet
call PosPart      !crée positions initiale


open(20, file = 'newpos.res')

do step = 1, nstep
  t = step*dt
!call newpos.res if it doesn't exist call use PosPart

  do i = 1, npart - 1; do p = i + 1, npart                 !pot ressenti par la part i des p autres part
    dx(i, p) = pos(i, 1) - pos(p, 1)                    !distance en x de la ième part avec les p autrees part
    dy(i, p) = pos(i, 2) - pos(p, 2)
  end do

  if (i == 1 ) then
    r1 = sqrt(minval(dx)**2 + minval(dy)**2)
    LennardJones(r)      !!distance en y de la ième part avec les p autrees part
    a0 = f / m
    v0 = gauss()
    v1 = v0
    theta = random_number()*2*Pi
  endif

  do j = 1, npart

    r2 = r1 + dt * v1 +((dt ** 2) / 2) * a0               !1=t et 2=t+1
    a1 = 24 / m * (2 * (R ** 12 / r2 ** 13) - ( R ** 6 / r2 ** 7))
    v2 = v1 + (dt / 2 ) * (a0 + a1)

    if (t > 1) then

      x = pos(j,1) + cos(theta)*v1*dt
      y = pos(j,2) + sin(theta)*v1*dt
      write(20,*) x, y    !nouvelles positions apres dt pour toutes les particules
    endif
  enddo

r1 = r2                                    !nouvelles valeurs
v1 = v2
a0 = a1
enddo

end program deuxD


!Position aleatoire des particules
function PosPart(pos)
use general
implicit none
intent(out) :: pos
real(8) xi, yi, xp, yp



createpoint: do i = 1, npart            !tirage de position aleatoire
call random_number(xi)               !0<x<1, jaimerais faire un vrai nombre aléatoire entre 0 et a
call random_number(yi)               !0<y<1
pos(i, 1) = xi
pos(i, 2) = yi
  if (i > 1) then
    do p = 1, i - 1
    xp = pos(p, 1)
    yp = pos(p, 2)
      if ((xi - xp) ** 2 + (yi - yp) ** 2 <= 4 * R ** 2) then
        exit createpoint
      endif
    enddo
  endif

enddo createpoint
close(10)

end function PosPart

!Verlet


! subroutine Pot_LennardJones
! use general
! implicit none
! call PosPart
!     do i = 1, npart
!     rx = x - rx
!     ry = y - ry
!     r = sqrt(rx ** 2 + ry ** 2)
!     U(i) = 4 * U0 * [(R / r) ** 12 - (R / r) ** 6]
! end do
! end

subroutine LennardJones(r,U,f)
use general
implicit none
real(8), intent(in) :: r
real(8), intent(out) :: U, f
real(8), dimension(npart - 1, npart - 1) :: dx, dy
real(8), dimension(1:npart - 1) :: Upar, d

f = 24 * (2 * (R ** 12 / r ** 13) - (R ** 6 / r ** 7))

do i = 1, npart - 1; do p = i + 1, npart                 !pot ressenti par la part i des p autres part
  dx(i, p) = pos(i, 1) - pos(p, 1)                    !distance en x de la ième part avec les p autrees part
  dy(i, p) = pos(i, 2) - pos(p, 2)
  d(i) = sqrt(minval(dx)**2 + minval(dy)**2)                   !!distance en y de la ième part avec les p autrees part
end do

Upar(i) = 4 * U0 * ((R / d(i) ** 12 - (R / d(i)) ** 6)
Ueff = sum(U(i))
U = 4 * U0 * ((R / r ) ** 12 - (R / r) ** 6)

end

real function gauss()
!    *******************************************************************
!    ** RANDOM VARIATE FROM THE STANDARD NORMAL DISTRIBUTION.         **
!    **                                                               **
!    ** THE DISTRIBUTION IS GAUSSIAN WITH ZERO MEAN AND UNIT VARIANCE.**
!    **                                                               **
!    ** REFERENCE:                                                    **
!    **                                                               **
!    ** KNUTH D, THE ART OF COMPUTER PROGRAMMING, (2ND EDITION        **
!    **    ADDISON-WESLEY), 1978                                      **
!    **                                                               **
!    *******************************************************************
implicit none
real, parameter :: a1 = 3.949846138, a3 = 0.252408784, a5 = 0.076542912, &
                   a7 = 0.008355968, a9 = 0.029899776
real :: r, r2, dummy
real, dimension(12) :: s

call random_number(s) ; r = (sum(s)-6.0)/4.0 ; r2 = r*r
gauss =  (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1 )*r
end function gauss

