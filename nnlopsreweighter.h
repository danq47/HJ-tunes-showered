c -*- Fortran -*-

c maxrwgtfiles.h
      integer maxRwgtFiles
      parameter (maxRwgtFiles=50)

c mxpart.h
      integer    bigmxpart
      parameter (bigmxpart=2048)

c jetcuts.h
      double precision ptjetmin,etajetmin,etajetmax,ptbjetmin,etabjetmax

c jetlabel.h
      integer jets
      character*2 jetlabel(bigmxpart)
c      common/parts_int/jets
c      common/parts_char/jetlabel
      
c npart.h
c----Number of partons in final state
      integer npart
c      common/npart/npart


      common/cnnlopsreweighter/ptjetmin,etajetmin,etajetmax,ptbjetmin,
     1 etabjetmax,jets,jetlabel,npart
