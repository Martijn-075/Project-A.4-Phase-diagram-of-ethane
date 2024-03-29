c potential.inc 
      double precision eps4,sig2,mass,ecut,rc2,pi
     $             ,eps48,rc, optbondlength
      logical tailco,shift
      integer chainlength, ntrialor
      character bondtype*5
      common/pot1/pi,eps4,eps48,sig2,mass,rc(2),rc2(2),ecut(2)
     $        ,tailco,shift
      common/pot2/optbondlength, ntrialor, chainlength, bondtype
c
c     pi        : pi = 3.14 ..
c     eps4      : 4 * epsilon 
c     eps48     : 48 * epsilon
c     (epsilon) : energy parameter Lennard-Jones potential
c     sig2      : sigma*sigma
c     (sigma)   : size parameter Lennard-Jones potenital
c     mass      : mass of the molecules
c     rc        : cut-off radius of the potenial
c     rc2       : rc * rc
c     ecut      : energy at cut-off radius
c     tailco    : .true. apply tail corrections
c     shift     : .true. shift the potential
c     chainlength: number of crabon atoms in alkane
c     optbondlength: optimal bond length of c-c in alkane
c     ntrialor    : number of trail postion of grwoing of chain
