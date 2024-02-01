c conf.inc
      double precision x,y,z, rchain
      integer npart,id,npbox
      common/conf1/x(npmax,2),y(npmax,2),z(npmax,2),id(npmax),npart
      common/conf2/npbox(2), rchain(npmax, 1)
c
c     x(i,cn),y(i,cn),z(i,cn)    : position particle i and carnon number cn
c     id(i)             : box 1 or 2
c     npart             : actual number of particles
c     npbox(i)          : number of particles in box i
c     rchain(i)         : bond lenght in chain
