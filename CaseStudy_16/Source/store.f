**==store.spg  processed by SPAG 4.52O  at 18:10 on 19 Jul 1996
      SUBROUTINE STORE(Iout, Dr, Vmax)
c     writes configuration to disk
      IMPLICIT NONE
      INCLUDE 'parameter.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'system.inc'
      INTEGER Iout, i
      DOUBLE PRECISION Dr, Vmax
 

! svaes to the lj.res file on iout 8 (gibbs program)
      write (iout, *) "box 1, hbox 1, box 2, hbox 2"
      WRITE (Iout, *) BOX(1), HBOX(1), BOX(2), HBOX(2)
      write (iout, *) "N particals, particales box 1, particles box 2"
      WRITE (Iout, *) NPART, NPBOX(1), NPBOX(2)
      write (iout, *) "Dr, Vmax"
      WRITE (Iout, *) Dr, Vmax
      write (iout, *) "X, Y, Z, box id"
      DO i = 1, NPART
         WRITE (Iout, *) X(i), Y(i), Z(i), ID(i)
      END DO
      REWIND (Iout)
      RETURN
      END
