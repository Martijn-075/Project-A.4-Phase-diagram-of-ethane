**==eneri.spg  processed by SPAG 4.52O  at 18:10 on 19 Jul 1996
      SUBROUTINE ENERI(Xi, Yi, Zi, I, Jb, En, Vir, Ib, cl)
      IMPLICIT NONE
      INCLUDE 'parameter.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'system.inc'
      INCLUDE 'potential.inc'
 
      DOUBLE PRECISION Xi(cl), Yi(cl), Zi(cl)
      DOUBLE PRECISION En, dx, dy, dz, r2, Vir, virij, enij
      INTEGER I, j, Jb, Ib, k, cl
c
      En = 0.D0
      Vir = 0.D0
      DO j = Jb, NPART
         ! check if particle is in the same box
         IF (ID(j).EQ.Ib) THEN
            ! check if it is not the same particle or chain
            IF (j.NE.I) THEN
               do k = 1,cl
                  dx = Xi(k) - X(j,k)
                  dy = Yi(k) - Y(j,k)
                  dz = Zi(k) - Z(j,k)

                  ! check dx
                  IF (dx.GT.HBOX(Ib)) THEN
                     dx = dx - BOX(Ib)
                  ELSE
                     IF (dx.LT.-HBOX(Ib)) dx = dx + BOX(Ib)
                  END IF
                  ! check dy
                  IF (dy.GT.HBOX(Ib)) THEN
                     dy = dy - BOX(Ib)
                  ELSE
                     IF (dy.LT.-HBOX(Ib)) dy = dy + BOX(Ib)
                  END IF
                  ! check dz
                  IF (dz.GT.HBOX(Ib)) THEN
                     dz = dz - BOX(Ib)
                  ELSE
                     IF (dz.LT.-HBOX(Ib)) dz = dz + BOX(Ib)
                  END IF
                  ! calcualte r2
                  r2 = dx*dx + dy*dy + dz*dz

                  ! calcualte energy

                  CALL ENER(enij, virij, r2, Ib)
                  En = En + enij
                  Vir = Vir + virij
               end do

            END IF
         END IF
      END DO
      RETURN
      END
