**==mcswap.spg  processed by SPAG 4.52O  at 18:10 on 19 Jul 1996
      SUBROUTINE MCSWAP(En, Vir, Attempt, Acc, Iseed)
c     ---exchange a particle bewteen the two boxes
 
      IMPLICIT NONE
 
      INCLUDE 'parameter.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'potential.inc'
      INCLUDE 'system.inc'
      INCLUDE 'chem.inc'
 
      DOUBLE PRECISION En, Vir, RANF, enn, virn, eno, viro,
     &                 xn(chainlength), yn(chainlength), sum,
     &                 zn(chainlength), arg, CORU, vola, vold, rhoan,
     &                 rhoao, rhodn, rhodo, dele, dtaila, dtaild,
     &                 rosenbluthn, xk(ntrialor), yk(ntrialor), 
     &                 zk(ntrialor), enk(ntrialor), prob, treshold,
     &                 rosenblutho, enkfirst, rbn(ntrialor), 
     &                 rbo(ntrialor), rand, ro, rn
      INTEGER Attempt, Iseed, o, iadd, idel, jb, idi, Acc, k, cn, co
      DIMENSION En(*), Vir(*)
 
 
      Attempt = Attempt + 1
c     ===select a box at random
      IF (RANF(Iseed).LT.0.5D0) THEN
         iadd = 1
         idel = 2
      ELSE
         iadd = 2
         idel = 1
      END IF

      !  check if there is a particle to remove
      IF (NPBOX(idel).EQ.0) RETURN

      vola = BOX(iadd)**3
      vold = BOX(idel)**3

c     ---add first particle to box iadd
      xn(1) = BOX(iadd)*RANF(Iseed)
      yn(1) = BOX(iadd)*RANF(Iseed)
      zn(1) = BOX(iadd)*RANF(Iseed)
c     ---calculate energy of this particle
      jb = 1
      o = NPART + 1
      CALL ENERI(xn(1), yn(1), zn(1), o, jb,
     &           enkfirst, virn, iadd, 1)

!  check if there is a chain to be created
      if (chainlength .ge. 2) then
            rosenbluthn = 0.D0
            xk = 0.D0
            yk = 0.D0
            zk = 0.D0
            enk = 0.D0

            !  chain with one particle has an internal energy of zero thus the probaility of creating
            ! of a trial insertion is equeal
            do k = 1,ntrialor
                  ! random on a sphere
<<<<<<< HEAD
<<<<<<< HEAD
                  rand = RANF(iseed)
                  xk(k) = xn(1) + optbondlength * 
     &            sin(rand * PI) * 
     &            cos(rand * 2.D0 * PI)
                  yk(k) = yn(1) + optbondlength * 
     &             sin(rand * PI) * 
     &            sin(rand * 2.D0 * PI)
                  zk(k) = zn(1) + optbondlength * 
     &            cos(rand * PI)
=======
=======
>>>>>>> ecf5840b54e9e69ae5125152c9fbc3a39cd297e7
                  xk(k) = xn(1) + optbondlength * 
     &            sin(RANF(iseed) * PI) * 
     &            cos(RANF(iseed) * 2.D0 * PI)
                  yk(k) = yn(1) + optbondlength * 
     &             sin(RANF(iseed) * PI) * 
     &            sin(RANF(iseed) * 2.D0 * PI)
                  zk(k) = zn(1) + optbondlength * 
     &            cos(RANF(iseed) * PI)
<<<<<<< HEAD
>>>>>>> ecf5840b54e9e69ae5125152c9fbc3a39cd297e7
=======
>>>>>>> ecf5840b54e9e69ae5125152c9fbc3a39cd297e7

                  CALL ENERI(xk(k), yk(k), zk(k), o, jb,
     &                 enk(k) ,virn, iadd, 1)
                  enk(k) = enk(k) + enkfirst
                  ! savegaurd against NAN with high energy
                  if (enk(k) .lt. 200.D0) then
                        rbn(k) = exp(-1.D0 * BETA * enk(k))
                  else 
                        rbn(k) = 0.D0
                  end if
                  rosenbluthn = rosenbluthn + rbn(k)
            end do


            ! save gaurd against NAN
            if (rosenbluthn .le. 1.d-10) RETURN

            ! get the second particle with the rosenbloth probability
            
            treshold = RANF(iseed)
            k = 0
            sum = 0.D0
            prob = 0.D0
            do while (sum .lt. treshold)
                  k = k + 1
                  ! not convered protection
                  if (k .gt. ntrialor) RETURN
                  sum = sum + rbn(k)
            end do

      if (rosenbluthn .le. 1.d-10) RETURN

            xn(2) = xk(k)
            yn(2) = yk(k)
            zn(2) = zk(k)
            CALL ENERI(xn(:), yn(:), zn(:), o, jb, enn, virn,
     &                       iadd, chainlength)


      else
            ! one particle no chain (old calculation)
            CALL ENERI(xn, yn, zn, o, jb, enn, virn,
     &                  iadd, 1)
      end if

c     ---calculate contibution to the chemical potential:
      arg = -BETA*enn
      IF (TAILCO) THEN
         rhoan = (NPBOX(iadd)+1)/vola
         arg = -BETA*(enn+2*CORU(RC(iadd),rhoan))
      END IF
      CHP(iadd) = CHP(iadd) + vola*EXP(arg)/DBLE(NPBOX(iadd)+1)
      IF (NPBOX(iadd).EQ.NPART) CHP(iadd) = CHP(iadd) + vola*EXP(arg)
     &    /DBLE(NPBOX(iadd)+1)
      ICHP(iadd) = ICHP(iadd) + 1

      idi = 0
      ! get particle random particle from the other box
      DO WHILE (idi.NE.idel)
         o = INT(NPART*RANF(Iseed)) + 1
         idi = ID(o)
      END DO

! rosenbluth of old configuration
      if (chainlength .ge. 2) then
            rosenblutho = 0.D0
            xk = 0.D0
            yk = 0.D0
            zk = 0.D0
            enk = 0.D0
      
            ! select wich carbon in the chain is sued for the rosenbluth
      
            if (RANF(iseed) .gt. 0.5D0) then
                  cn = 1
                  co = 2
            else
                  cn = 2
                  co = 1
            end if

            CALL ENERI(x(o,:), y(o,:), z(o,:), o, jb,
     &                 enk(1), virn, idi, chainlength)
     
            CALL ENERI(x(o,co), y(o,co), z(o,co), o, jb,
     &                 enkfirst, virn, idi, 1)



            rbo(1) = exp(-1.D0 * BETA * enk(1))
            xk(1) = x(o,cn)
            yk(1) = y(o,cn)
            zk(1) = z(o,cn)
            
            do k = 2,ntrialor
                  ! random on a sphere
                  rand = RANF(Iseed)
                  xk(k) = x(o,co) + optbondlength * 
     &             sin(rand * PI) * 
     &             cos(rand * 2.D0 * PI)
                  yk(k) = y(o,co) + optbondlength * 
     &            sin(rand * PI) * 
     &            sin(rand * 2.D0 * PI)
                  zk(k) = z(o,co) + optbondlength * 
     &            cos(rand * PI)

                  CALL ENERI(xk(k), yk(k), zk(k), o, jb,
     &            enk(k),virn, idi, 1)
                  enk(k) = enk(k) + enkfirst

<<<<<<< HEAD
<<<<<<< HEAD
                  if (enk(k) .lt. 200.D0) then
                        rbo(k) = exp(-1.D0 * BETA * enk(k))
                  else
                        rbo(k) = 0.D0
                  end if
                  
=======
                  rbo(k) = exp(-1.D0 * BETA * enk(k))
>>>>>>> ecf5840b54e9e69ae5125152c9fbc3a39cd297e7
=======
                  rbo(k) = exp(-1.D0 * BETA * enk(k))
>>>>>>> ecf5840b54e9e69ae5125152c9fbc3a39cd297e7
                  rosenblutho = rosenblutho + rbo(k)
            end do
      else
            IF (NPBOX(idel).EQ.0) RETURN

            CALL ENERI(X(o,1), Y(o,1), Z(o,1), o, jb, eno, viro, 
     &                  idel, 1)
      end if

 
! not used due to rosenbluth
! c     ---delete particle from box b:
!       IF (NPBOX(idel).EQ.0) THEN
!          RETURN
!       END IF
!       idi = 0
!       ! get particle to delete from the right box 
!       DO WHILE (idi.NE.idel)
!          o = INT(NPART*RANF(Iseed)) + 1
!          idi = ID(o)
!       END DO
!       ! calculate energy of partcile to be removed
!       CALL ENERI(X(o,:), Y(o,:), Z(o,:), o, jb, eno, viro, idel,
!      &           chainlength)
 
c     ---acceptence test:
! first check if the rosenbluth (chain) ore a normal acceptnace rule needs to be used
      if (chainlength.ge. 2) then
            dele = (vola * NPBOX(idel)) / (vold * (NPBOX(iadd) + 1)) *
     &            (rosenbluthn / rosenblutho) 
            IF (RANF(Iseed).LT.dele) THEN
c        ---accepted:

                  CALL ENERI(x(o,:), y(o,:), z(o,:), o, jb,
     &                 eno, viro, idel, chainlength)

                  Acc = Acc + 1
                  NPBOX(iadd) = NPBOX(iadd) + 1
                  X(o,:) = xn(:)
                  Y(o,:) = yn(:)
                  Z(o,:) = zn(:)
                  ID(o) = iadd

                  En(iadd) = En(iadd) + enn

                  IF (TAILCO) En(iadd) = En(iadd) + dtaila
                  Vir(iadd) = Vir(iadd) + virn
                  NPBOX(idel) = NPBOX(idel) - 1
                  En(idel) = En(idel) - eno
                  IF (TAILCO) En(idel) = En(idel) + dtaild
                  Vir(idel) = Vir(idel) - viro
            END IF
      else
            dele = enn - eno + LOG(vold*(NPBOX(iadd)+1)/
     &       (vola*NPBOX(idel)))/BETA
            IF (TAILCO) THEN
c        ---tail corrections:
            rhoao = NPBOX(iadd)/vola
            dtaila = (NPBOX(iadd)+1)*CORU(RC(iadd), rhoan) - NPBOX(iadd)
     &            *CORU(RC(iadd), rhoao)
            rhodn = (NPBOX(idel)-1)/vold
            rhodo = NPBOX(idel)/vold
            dtaild = (NPBOX(idel)-1)*CORU(RC(idel), rhodn) - NPBOX(idel)
     &            *CORU(RC(idel), rhodo)
            dele = dele + dtaila + dtaild
            END IF

            IF (RANF(Iseed).LT.EXP(-BETA*dele)) THEN
c        ---accepted:
                  Acc = Acc + 1
                  NPBOX(iadd) = NPBOX(iadd) + 1
                  X(o,1) = xn(1)
                  Y(o,1) = yn(1)
                  Z(o,1) = zn(1)
                  ID(o) = iadd
                  En(iadd) = En(iadd) + enn
                  IF (TAILCO) En(iadd) = En(iadd) + dtaila
                  Vir(iadd) = Vir(iadd) + virn
                  NPBOX(idel) = NPBOX(idel) - 1
                  En(idel) = En(idel) - eno
                  IF (TAILCO) En(idel) = En(idel) + dtaild
                  Vir(idel) = Vir(idel) - viro
            END IF
      end if

      RETURN
      END
