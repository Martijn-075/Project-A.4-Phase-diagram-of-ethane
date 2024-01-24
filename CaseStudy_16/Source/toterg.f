**==toterg.spg  processed by SPAG 4.52O  at 18:10 on 19 Jul 1996
      SUBROUTINE TOTERG(Ener, Vir, Ib)
c     ---calculates total energy
      IMPLICIT NONE
      INCLUDE 'parameter.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'potential.inc'
      INCLUDE 'system.inc'
 
      DOUBLE PRECISION xi, yi, zi, Ener, eni, CORU, viri, Vir, rho
      INTEGER i, jb, Ib
 
      Ener = 0
      Vir = 0
      DO i = 1, NPART - 1
         IF (ID(i).EQ.Ib) THEN
            xi = X(i)
            yi = Y(i)
            zi = Z(i)
            jb = i + 1
            CALL ENERI(xi, yi, zi, i, jb, eni, viri, Ib)
            Ener = Ener + eni
            Vir = Vir + viri
         END IF
      END DO
c     ---add tail corrections
      IF (TAILCO) THEN
         rho = NPBOX(Ib)/(BOX(Ib)**3)
         Ener = Ener + NPBOX(Ib)*CORU(RC(Ib), rho)
      END IF
      RETURN
      END
 
 
 
