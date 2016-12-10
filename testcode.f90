module general
implicit none
integer :: i, p
integer, parameter :: a = 10, npart = 100, T = 500
real(8), dimension(npart, 2) :: pos
real(8), parameter :: R = 0.1, m = 0.00001, v0 = 1, U0 = 1


end module general

program deuxD
use general
implicit none


call PosPart
end program deuxD

!Position al√©atoire des particules
subroutine PosPart
use general
implicit none
real(8) xi, yi, xp, yp

open(10, file = 'positions.res')

createpoint: do i = 1, npart            !tirage de position al√©atoire
call random_number(xi)               !0<x<1, jaimerais faire un vrai nombre alÈatoire entre 0 et a
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
write(10, *) pos(i, 1), pos(i, 2)

enddo createpoint
close(10)

end

subroutine Verlet
use general
implicit none
real(8) :: v1, v2, r0, r1, r2, a0, a1, a2, xi, yi, dt = 0.01

open(20, file = 'verlet.res')


do i = 1, npart, 1
read(10, *) xi, yi

r0 = sqrt(xi ** 2 + yi ** 2)

a0 = 24 / m * (2 * (R ** 12 / r0 ** 13) - (R ** 6 / r0 ** 7))
r1 = r0
v1 = v0

    do p = 0, T, 1

    r2 = r1 + dt * v1 +((dt ** 2) / 2) * a0               !1=t et 2=t+1
    a1 = 24 / m * (2 * (R ** 12 / r2 ** 13) - ( R ** 6 / r2 ** 7))
    v2 = v1 + (dt / 2 ) * (a0 + a1)

    r1 = r2                                    !nouvelles valeurs
    v1 = v2
    a0 = a1

    enddo

write(20, *) r1, v1, a0
enddo
close(20)
end

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

subroutine LennardJones
use general
implicit none
real(8), dimension(npart - 1, npart - 1) :: U, dx, dy
real(8), dimension(1:npart - 1) :: Ueff
call PosPart


do i = 1, npart - 1; do p = i + 1, npart                 !pot ressenti par la part i des p autres part
  dx(i, p) = pos(i, 1) - pos(p, 1)                    !distance en x de la iËme part avec les p autrees part
  dy(i, p) = pos(i, 2) - pos(p, 2)                    !!distance en y de la iËme part avec les p autrees part  
  U(i, p) = 4 * U0 * ((R / sqrt(dx(i, p) ** 2 + dy(i, p) ** 2)) ** 12 - (R / sqrt(dx(i, p) ** 2 + dy(i, p) ** 2)) ** 6)
end do
Ueff(i) = sum(U(i,:))      !somme de tous les elements de la ligne i
end do

end
