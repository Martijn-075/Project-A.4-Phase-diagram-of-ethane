        subroutine MCROT(En, Vir, Attempt, Nacc, Iseed)

        IMPLICIT NONE
        INCLUDE 'parameter.inc'
        INCLUDE 'conf.inc'
        INCLUDE 'system.inc'
        INCLUDE "potential.inc"

        DOUBLE PRECISION En(*), Vir(*), RANF, eno, enn,
     &  viro, virn, rand, rn, ro,
     &  xn(chainlength), yn(chainlength), zn(chainlength)
        INTEGER Attempt, Nacc, Iseed, ic, ie, o, ido, jb

        ! return if there is no chain
        if (chainlength .lt. 2) RETURN

        Attempt = Attempt + 1
        jb = 1

c     ---select a particle at random
        o = INT(NPART*RANF(Iseed)) + 1
        ido = ID(o)

c     ---calculate energy old configuration
      CALL ENERI(X(o,:), Y(o,:), Z(o,:), o, jb, eno, viro, ido,
     &           chainlength)

        ! selct rendom wich part of the chain will be rotated
        if (RANF(Iseed) .lt. 0.5D0) then
                ic = 1
                ie = 2
        else
                ic = 2
                ie = 1
        end if

        xn(1) = x(o,ic)
        yn(1) = y(o,ic)
        zn(1) = z(o,ic)

        ! random on a sphere
        rand = RANF(iseed)
        xn(2) = xn(1) + optbondlength * 
     &  sin(rand * PI) * 
     &  cos(rand * 2.D0 * PI)
        yn(2) = yn(1) + optbondlength * 
     &  sin(rand * PI) * 
     &  sin(rand * 2.D0 * PI)
        zn(2) = zn(1) + optbondlength * 
     &  cos(rand * PI)

c     ---calculate energy new configuration:
        CALL ENERI(xn(:), yn(:), zn(:), o, jb, enn, virn, ido, 
     &  chainlength)
c     ---acceptance test
        IF (RANF(Iseed).LT.EXP(-BETA*(enn-eno))) THEN
c        --accepted
                Nacc = Nacc + 1
                En(ido) = En(ido) + (enn-eno)
                Vir(ido) = Vir(ido) + (virn-viro)
c        ---put particle in simulation box 
!                 ro = sqrt((x(o,1) - 
!      &  x(o,2))**2 + (y(o,1) - y(o,2))**2 + (z(o,1) - z(o,2))**2)

!                 rn = sqrt((xn(1) - 
!      &          xn(2))**2 + (yn(1) - yn(2))**2 + (zn(1) - zn(2))**2)
!                 write(6,*) o, 'ro', ro, 'rn', rn, 'del', ro - rn
!                 write(6,*) ido, box(1), box(2)
!                 write(6,*) "xyz before", xn, yn, zn
                IF (xn(ie).LT.0.) xn(ie) = xn(ie) + BOX(ido)
                IF (xn(ie).GT.BOX(ido)) xn(ie) = xn(ie) - BOX(ido)
                IF (yn(ie).LT.0.) yn(ie) = yn(ie) + BOX(ido)
                IF (yn(ie).GT.BOX(ido)) yn(ie) = yn(ie) - BOX(ido)
                IF (zn(ie).LT.0.) zn(ie) = zn(ie) + BOX(ido)
                IF (zn(ie).GT.BOX(ido)) zn(ie) = zn(ie) - BOX(ido)
!                 write(6,*) "xyz after", xn, yn, zn
!                 write(6,*) "after"
!                 ro = sqrt((x(o,1) - 
!      &    x(o,2))**2 + (y(o,1) - y(o,2))**2 + (z(o,1) - z(o,2))**2)

!                 rn = sqrt((xn(1) - 
!      &          xn(2))**2 + (yn(1) - yn(2))**2 + (zn(1) - zn(2))**2)
!                 write(6,*) o, 'ro', ro, 'rn', rn, 'del', ro - rn

                X(o,:) = xn(:)
                Y(o,:) = yn(:)
                Z(o,:) = zn(:)
        END IF





        END