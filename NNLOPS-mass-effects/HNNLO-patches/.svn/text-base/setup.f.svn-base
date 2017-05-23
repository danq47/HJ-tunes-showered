      subroutine setup
      implicit none
      include 'constants.f'
      include 'virtonly.f'
      include 'realonly.f'
      include 'noglue.f'
      include 'lc.f'
      include 'cutoff.f'
      include 'maxwt.f'
      include 'masses.f'
      include 'process.f'
      include 'scale.f'
      include 'facscale.f'
      include 'zerowidth.f'
      include 'vzerowidth.f'
      include 'removebr.f'
      include 'clustering.f'
      include 'gridinfo.f'
      include 'limits.f'
      include 'alfacut.f'
      include 'pdfiset.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'nlooprun.f'
      include 'rescoeff.f'
      include 'jetcuts.f'
CC
      include 'vegas_common.f'
      include 'nwz.f'
      include 'ewcharge.f'
      include 'lhapdf.f'
      character *2 plabel(mxpart)
      common/plabel/plabel
      integer order,notag,nqcdjets,nqcdstart,isub,higgsdec,ndec
      common/nnlo/order
      common/notag/notag
      common/nqcdjets/nqcdjets,nqcdstart
      common/isub/isub
      double precision BrnRat,gamgambr,wwbr,zzbr,br0
      common/BrnRat/BrnRat
      common/higgsdec/higgsdec,ndec
      common/rescale/resfac

      double precision qtcut
      common/qtcut/qtcut      

      double precision beta1,H2ggdelta,H2ggD0
      common/Hstcoeff/beta1,H2ggdelta,H2ggD0

      logical isol
      common/isol/isol

      logical int
      common/int/int

      logical dorebin
      common/dorebin/dorebin

      character *50 prefix
      character *36 pdfstring
      integer nset
      common/prefix/nset,prefix
      common/pdfstring/pdfstring
ch
      character*4 part
      character*30 runstring
      integer j
      logical makecuts
      integer nmin,nmax,n2,n3,approxim
      integer ih1,ih2,itmx1,itmx2,ncall1,ncall2,idum,rseed
      double precision rtsmin,sroot,LT,c1ggdeltaf
      double precision Mwmin,Mwmax
      double precision Rcut,resfac
      double precision ran2,randummy
      double precision cmass,bmass,mtop,mbot,sigmatb,sigmat
      double precision mass2,width2,mass3,width3
      double precision amz,alphas
      double precision brwen,brzee,brtau,brtop
      double precision abisq,abisqnnlo
      
      character *3 str1
      character *8 str2
      character *38 str3
      character *50 string

ch      commom/flags/approxim
      common/couple/amz
      
      common/breit/n2,n3,mass2,width2,mass3,width3
      
      common/nmin/nmin
      common/nmax/nmax
      common/rtsmin/rtsmin
 

      common/part/part
      common/runstring/runstring
      common/energy/sroot
      common/density/ih1,ih2
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/ranno/idum
      common/flag/approxim
      common/tbmass/mtop,mbot


CH    Declaration of new common block for rescaling the Born term
      common/born/sigmatb,sigmat
CH

      
      common/Rcut/Rcut
      common/makecuts/makecuts

      common/qmass/cmass,bmass

      common/rseed/rseed
      save /ranno/

      logical lhapdfs
      common/lhapdfs/lhapdfs

      lhapdfs=.false.

      isol=.false.

      virtonly=.false.
      realonly=.false.

      noglue=.false.
      ggonly=.false.
      gqonly=.false.

      nmin=1
      nmax=2

      clustering=.true.
      colourchoice=0
      rtsmin=40d0
      cutoff=0.01d0


CC
      aii=1d0
      aif=1d0
      afi=1d0
      aff=1d0

CC Finite width effects for the Higgs boson            

      zerowidth=.true.

CC Finite width effects for W and Z decay

      vzerowidth=.false.

C     Limits on invariant mass of Higgs decay products
C     (irrelevant if zerowidth=.true.)

      
      Mwmin=40d0
      Mwmax=14000d0

      inclusive=.true.
      algorithm='ktal'


CC   Parameters used to define jets

      ptjetmin=0d0
      etajetmin=0d0
      etajetmax=20d0

      Rcut=0.7d0

CC

C -      removebr=.false.
      removebr=.true.   ! Now we remove the Higgs branching ratios from results
      makecuts=.true.     


CC    Adjust the grid at each iteration

      dorebin=.true.


CC    Read a previously saved grid

      readin=.false.


      writeout=.false.
      ingridfile=''
      outgridfile=''

CC    Read inputfile


      read(*,*) sroot
      read(*,*) ih1,ih2
      read(*,*) hmass
      read(*,*) scale,facscale  ! mur,muf
      read(*,*) order        
      read(*,*) higgsdec !decay mode
      read(*,*) part
      read(*,*) itmx1,ncall1
      read(*,*) itmx2,ncall2
      read(*,*) rseed          ! seed
      read(*,*) iset,nset ! pdf set, eigenvector set
CH
      read(*,*) mtop
      read(*,*) mbot
CH
      read(*,*) approxim      
      read(*,*) PDFname,PDFmember
      read(*,*) runstring
      
      

CC    Cut on qt

      qtcut=0.8d0


      wsqmin=Mwmin**2
      wsqmax=Mwmax**2

C     Check if the limits are compatible with sroot
C     The minum x cannot go below 10^(-5):
C     At the LHC the minimum invariant mass cannot go below 44.27 GeV

      if(wsqmin.lt.(sroot**2*1d-5)) wsqmin=(sroot**2*1d-5)
      if(wsqmax.gt.(sroot**2)) wsqmax=sroot**2


      do j=1,mxpart
      plabel(j)=''
      enddo


      plabel(1)='pp'
      plabel(2)='pp'

c--- the default behaviour is to remove no branching ratio

      BrnRat=1d0

      call coupling
 
      call cstring(pdfstring)

      if(lhapdfs.eqv.(.false.)) then
      write(6,*)'CCCCCCCCCC      Parton Distributions    CCCCCCCCCCCC'
      write(6,*)'C                                                  C'
      write(6,*)'C       ', pdfstring,'       C'
      write(6,*)'C                                                  C'
      endif
      write(6,*)'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'

      notag=0
      nqcdjets=0
      isub=0


      
      mb=0
    
CC    Compute Higgs BR ratios

C     br0 is tree level H->bbar BR which is needed here

      call sethparams(br0,gamgambr,wwbr,zzbr)

      case='ggfus1'
    
     
      if((ih1.eq.1).and.(ih2.eq.-1)) then
       str1='ppb'
      elseif((ih1.eq.1).and.(ih2.eq.1)) then
       str1='pp'
      else
       write(*,*)'Initial state not allowed'
      endif

      str2=' -> H -> ' 

CC    Here ndim is the number of dimensions for Born

      if(higgsdec.eq.1) then
C
C     H->2gamma
C
       ndec=2
       ndim=4    
       plabel(3)='tl'
       plabel(4)='ta'
       plabel(5)='pp'
       plabel(6)='pp'
       n2=0
       n3=1
       mass3=hmass
       width3=hwidth

C - If removebr is true then all MEs will have this Brnrat=br0 divided
C - out of them. br0 is the Higgs->bbbar branching fraction. When higgsdec
C - equals 1, i.e. when Higgs->gamma gamma production is selected (default)
C - then the MEs that are calculated are the H->bbbar ones. So by dividing
C - out Brnrat=br0 you go back to the total production cross section. You
C - can see that H->bbbar is what is actually calculated when H->gamma gamma
C - is selected if you look in e.g. Matrix/gg_h.f . By the same token, if
C - removebr is false then dividing out Brnrat=br0/gamgambr divides out the
C - H->bbbar branching fraction and multiplies in the gamma gamma one.
       if (removebr) then
        Brnrat=br0
       else
        Brnrat=br0/gamgambr
       endif

      str3=' gamma(p3)+gamma(p4)'

      elseif(higgsdec.eq.2) then
C
C      H->2W->4l
C
       ndec=4
       ndim=10
       plabel(3)='nl'
       plabel(4)='ea'
       plabel(5)='el'
       plabel(6)='na'
       plabel(7)='pp'
       plabel(8)='pp'
       n2=1
       n3=1
       mass2=wmass
       width2=wwidth
       mass3=wmass
       width3=wwidth
       
       if(zerowidth.and.(hmass.lt.2*wmass).and.vzerowidth)then
        write(*,*)'MH < 2MW !'
        stop
       endif

       if (removebr) then
        call branch(brwen,brzee,brtau,brtop)
        BrnRat=brwen**2*wwbr
       endif

      str3=' WW -> nu(p3)+e+(p4)+e-(p5)+nubar(p6)'

      elseif(higgsdec.eq.31) then
C
C      H->2Z->mu-mu+e-e+
C
       ndec=4
       ndim=10

       l1=le
       r1=re
       l2=le
       r2=re

       plabel(3)='el'
       plabel(4)='ea'
       plabel(5)='ml'
       plabel(6)='ma'
       plabel(7)='pp'
       plabel(8)='pp'
       n2=1
       n3=1
       mass2=zmass
       width2=zwidth
       mass3=zmass
       width3=zwidth

       int=.false.

       str3=' ZZ -> mu-(p3)+mu+(p4)+e-(p5)+e+(p6)'
       
       if(zerowidth.and.(hmass.lt.2*zmass).and.vzerowidth) then
        write(*,*)'MH < 2MZ !'
        stop
       endif

       if (removebr) then
        call branch(brwen,brzee,brtau,brtop)
        BrnRat=2d0*brzee**2*zzbr ! factor of 2 for identical particles      
       endif


C
C      H->2Z->e-e+e-e+
C

       elseif(higgsdec.eq.32) then

       ndec=4
       ndim=10

       l1=le
       r1=re
       l2=le
       r2=re

       plabel(3)='el'
       plabel(4)='ea'
       plabel(5)='el'
       plabel(6)='ea'
       plabel(7)='pp'
       plabel(8)='pp'
       n2=1
       n3=1
       mass2=zmass
       width2=zwidth
       mass3=zmass
       width3=zwidth


       str3=' ZZ -> e-(p3)+e+(p4)+e-(p5)+e+(p6)'

       if((hmass.lt.2*zmass).and.vzerowidth.eqv..true.) then
       write(*,*)'MH < 2MZ !'
       stop
       endif
        
       int=.true.   

      else
       write(*,*)'Wrong decay channel'
       stop
      endif

CCCCCCCCCCCCC

      call strcat(str1,str2,string)
      call strcat(string,str3,string)

     
      call cstring(string)      
     

      write(6,*)'C                                                  C'

      if(order.eq.0) then
      write(6,*)'C         Computing LO cross section for           C'
      elseif(order.eq.1) then
      write(6,*)'C        Computing NLO cross section for           C'
      elseif(order.eq.2) then
      write(6,*)'C        Computing NNLO cross section for          C'
      else
      write(*,*)'Order can be 0,1 or 2 !'
       stop
      endif
  
      write(6,*)'C                                                  C' 
      write(*,*)'C',string,'C'
      write(6,*)'C                                                  C' 
      write(6,96)sroot
      write(6,*)'C                                                  C'
      write(6,*)'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC' 


 96   format(' C           at  Sqrt(s)=',f8.2,' GeV               C')

      nqcdjets=0
    

      call ckmfill(nwz)

    

CCCCCCCCCCCCCCCCCC


c--- set-up the random number generator with a negative seed
      idum=-abs(rseed)
      randummy=ran2()

c--- initialize masses for alpha_s routine
      cmass=dsqrt(mcsq)
      bmass=mbot


c--- check that we have a valid value of 'part'
      if ( (part .ne. 'lord') .and. (part .ne. 'real') .and.
     .     (part .ne. 'virt') .and. (part .ne. 'tota') ) then
          write(6,*) 'part=',part,' is not a valid option'
          stop     
      endif      



        as=alphas(scale,amz,nlooprun)
        ason2pi=as/twopi
        ason4pi=as/fourpi
        gsq=fourpi*as
        musq=scale**2

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC    Resummation coefficients

      beta0=(33-2*nf)/12d0

      beta1=(153d0-19*nf)/24d0

      Kappa=67/6d0-(pi**2)/2d0-5d0/9d0*nf

      A1g=3d0
      A2g=0.5d0*A1g*Kappa
      B1g=-2*beta0

C     B2g from qt code (checked !)

      B2g=9*(23/24d0+(11*pi**2)/18d0-3*Z3/2d0)+
     /    2*nf/3d0-3*nf*(1/12d0+pi**2/9d0)-11/2d0

      

C     Delta term in c1gg coefficient

      C1ggdelta=(11+3*pi**2)/4d0
      if(approxim.eq.0)then
        C1ggdeltaex=(11+3*pi**2)/4d0
      else
       C1ggdeltaex=c1ggdeltaf(0,hmass,mtop,mbot)
      endif

C     Delta term in P2gg splitting function (as/pi normalization)

      Delta2gg=9*(8d0/3+3*Z3)-2d0/3*nf-2*nf

      Delta2gg=Delta2gg/4d0

CC    Coefficients of D0 and D1 in P*P (as/pi normalization)

      D0gggg=6*beta0
      D1gggg=18d0

CC    Coefficients of delta(1-z) in P*P

      Deltagggg=beta0**2-3d0/2*pi**2

C     H2gg contribution: coefficient of delta(1-z)

      LT=2*dlog(hmass/mtop)

      H2ggdelta=11399d0/144+19d0/8*LT+133*Pi**2/8+13d0/16*Pi**4
     #          -55d0/2*Z3+nf*(-1189d0/144+2d0*LT/3-5*Pi**2/12)

C     H2gg contribution: coefficient of D0(z)

      H2ggD0=-101d0/3+14d0*nf/9d0+63d0/2*Z3


CH    
      sigmatb=abisq(mtop**2/hmass**2,mbot**2/hmass**2)
      sigmat=abisqnnlo(mtop**2/hmass**2)
      if (approxim.ne.0) then
       resfac=sigmat
      elseif(approxim.eq.0) then
       resfac=1d0
      endif
CH
      
      return

      end
      
