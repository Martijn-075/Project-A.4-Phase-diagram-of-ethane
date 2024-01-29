! ==mcmove.spg  processed by SPAG 4.52O  at 18:10 on 19 Jul 1996
      SUBROUTINE MCMOVE(En, Vir, Attempt, Nacc, Dr, Iseed)
c     attempts to displace a randomly selected particle
      IMPLICIT NONE
      INCLUDE 'parameter.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'system.inc'
      INCLUDE "potential.inc"
      DOUBLE PRECISION enn, eno, En, RANF, xn(chainlength),
     &                 yn(chainlength), zn(chainlength), 
     &                 viro, virn, Dr, Vir, ro, rn
      DIMENSION En(*), Vir(*)
      INTEGER o, Attempt, Nacc, jb, Iseed, ido, i
 
      Attempt = Attempt + 1
      jb = 1
c     ---select a particle at random
      o = INT(NPART*RANF(Iseed)) + 1
      ido = ID(o)
c     ---calculate energy old configuration
      CALL ENERI(X(o,:), Y(o,:), Z(o,:), o, jb, eno, viro, ido,
     &           chainlength)
c     ---give particle a random displacement

      xn(:) = X(o,:) + (RANF(Iseed)-0.5D0)*Dr
      yn(:) = Y(o,:) + (RANF(Iseed)-0.5D0)*Dr
      zn(:) = Z(o,:) + (RANF(Iseed)-0.5D0)*Dr

c     ---calculate energy new configuration:
      CALL ENERI(xn, yn, zn, o, jb, enn, virn, ido, chainlength)
c     ---acceptance test
      IF (RANF(Iseed).LT.EXP(-BETA*(enn-eno))) THEN
c        --accepted
         Nacc = Nacc + 1
         En(ido) = En(ido) + (enn-eno)
         Vir(ido) = Vir(ido) + (virn-viro)
c        ---put particle in simulation box
            do i = 1, chainlength
                  IF (xn(i).LT.0.) xn(i) = xn(i) + BOX(ido)
                  IF (xn(i).GT.BOX(ido)) xn(i) = xn(i) - BOX(ido)
                  IF (yn(i).LT.0.) yn(i) = yn(i) + BOX(ido)
                  IF (yn(i).GT.BOX(ido)) yn(i) = yn(i) - BOX(ido)
                  IF (zn(i).LT.0.) zn(i) = zn(i) + BOX(ido)
                  IF (zn(i).GT.BOX(ido)) zn(i) = zn(i) - BOX(ido)
            end do

         X(o,:) = xn(:)
         Y(o,:) = yn(:)
         Z(o,:) = zn(:)
      END IF
      RETURN
      END
