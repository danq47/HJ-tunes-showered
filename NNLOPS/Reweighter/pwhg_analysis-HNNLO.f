      subroutine init_hist
      implicit none
      include 'LesHouches.h'
      include 'pwhg_bookhist-multi.h'
      real*8   dyh,dpth,dHT,dptj1
      call inihists

C - Set bin widths
      dyh=0.1d0
      dpth=3d0
      dHT=5d0
      dptj1=3d0

      call bookupeqbins('yh',             dyh, -5d0,  5d0)
      call bookupeqbins('yh-pth-1',       dyh, -5d0,  5d0)
      call bookupeqbins('yh-pth-2',       dyh, -5d0,  5d0)
      call bookupeqbins('yh-pth-3',       dyh, -5d0,  5d0)
      call bookupeqbins('yh-pth-4',       dyh, -5d0,  5d0)
      call bookupeqbins('yh-HT-1',        dyh, -5d0,  5d0)
      call bookupeqbins('yh-HT-2',        dyh, -5d0,  5d0)
      call bookupeqbins('yh-HT-3',        dyh, -5d0,  5d0)
      call bookupeqbins('yh-HT-4',        dyh, -5d0,  5d0)
      call bookupeqbins('yh-ptj1-1',      dyh, -5d0,  5d0)
      call bookupeqbins('yh-ptj1-2',      dyh, -5d0,  5d0)
      call bookupeqbins('yh-ptj1-3',      dyh, -5d0,  5d0)
      call bookupeqbins('yh-ptj1-4',      dyh, -5d0,  5d0)
      call bookupeqbins('pth',           dpth,  0d0,300d0)
      call bookupeqbins('pth-hfact-3',   dpth,  0d0,300d0)
      call bookupeqbins('HT',             dHT,100d0,600d0)
      call bookupeqbins('HT-hfact-2',     dHT,100d0,600d0)
      call bookupeqbins('ptj1',         dptj1,  0d0,300d0)
      call bookupeqbins('ptj1-hfact-3', dptj1,  0d0,300d0)
      call bookupeqbins('ptj2',         dptj1,  0d0,300d0)

      write(*,*) 'booked ',jhist,' histograms !'
      end


      subroutine analysis(dsig0)
      implicit  none
      include  'hepevt.h'
      include  'pwhg_weights.h'
      real * 8  dsig0,dsig(1:weights_max)
C - Maximum number of particles to send to jet finder: 
      include  'mxpart.f'
      character*6 WHCPRG
      common/cWHCPRG/WHCPRG
      data    WHCPRG/'NLO   '/
C - Higgs kinematics
      real*8  hmass
      common/hmass/hmass
      real*8  ph(4),yh,etah,pth,mh
C - Jet radius: 
      real*8  Rmin !:
      data    Rmin/0.4d0/
      save    Rmin
C - HNNLO jet finding:
      real*8  q(bigmxpart,4),p(bigmxpart,4)
      real*8  HT,ptj1,ptj2,pjm,ptj
      real*8  damping_vec(20) !:
C - Auxiliary counters:
      integer ixx,nu,i,njets
      integer nBosons
C - Kinematic functions from Need/ptyrap.f
      real*8  pt
C - Initialisation indicator (for new hisogramming setup).
      logical ini
      data ini/.true./
      save ini

      if(ini) then
         hmass = 125.5d0  
         write(*,*) 'Higgs mass set to', hmass 
         if(hmass.le.0d0) then
            write(*,*) 'pwhg_analysis-HNNLO.f:'
            write(*,*) 'hmass = ',hmass
            write(*,*) 'Have not got acceptable Higgs mass'
            stop
         endif
         if(WHCPRG.eq.'NLO') then
            weights_num=0
         endif
         if(weights_num.eq.0) then
            call setupmulti(1)
         else
            call setupmulti(weights_num)
         endif
         ini=.false.
      endif

      dsig=0
      if(weights_num.eq.0) then
         dsig(1)=dsig0
      else
         dsig(1:weights_num)=weights_val(1:weights_num)
      endif
      if(sum(abs(dsig)).eq.0) return

C     find jets 
      call find_jets(Rmin,p,njets) ! p are output jet momenta 

C     find the Higgs
      do ixx=1,nhep
         if(idhep(ixx).eq.25) then
            ph=phep(1:4,ixx)
            cycle
         endif
      enddo
     
C     get Higgs kinematics
      call getyetaptmass(ph,yh,etah,pth,mh)
      HT=sqrt(pth*pth+mh*mh)
      call get_ptj1j2(p,ptj1,ptj2)
      HT=HT+ptj1+ptj2

C     compute damping and fill histograms 
      call get_hdamp(yH,pth,HT,ptj1,ptj2,damping_vec)

C - yh
      call filld('yh',            yh,damping_vec(1)*dsig)
      
C - yh with pth damping function
      call filld('yh-pth-1',      yh,damping_vec(2)*dsig)
      call filld('yh-pth-2',      yh,damping_vec(3)*dsig)
      call filld('yh-pth-3',      yh,damping_vec(4)*dsig)
      call filld('yh-pth-4',      yh,damping_vec(5)*dsig)

C - yh with HT damping function
      call filld('yh-HT-1',       yh,damping_vec(6)*dsig)
      call filld('yh-HT-2',       yh,damping_vec(7)*dsig)
      call filld('yh-HT-3',       yh,damping_vec(8)*dsig)
      call filld('yh-HT-4',       yh,damping_vec(9)*dsig)

C - yh with ptj1 damping function
      call filld('yh-ptj1-1',     yh,damping_vec(10)*dsig)
      call filld('yh-ptj1-2',     yh,damping_vec(11)*dsig)
      call filld('yh-ptj1-3',     yh,damping_vec(12)*dsig)
      call filld('yh-ptj1-4',     yh,damping_vec(13)*dsig)

C - pth
      call filld('pth',           pth,damping_vec(14)*dsig)
C - pth with canonical damping factor
      call filld('pth-hfact-3',   pth,damping_vec(15)*dsig)

C - HT
      call filld('HT',            HT,damping_vec(16)*dsig)
C - HT with canonical damping factor
      call filld('HT-hfact-2',    HT,damping_vec(17)*dsig)

C - ptj1
      call filld('ptj1',        ptj1,damping_vec(18)*dsig)
C - ptj1 with canonical damping factor
      call filld('ptj1-hfact-3',ptj1,damping_vec(19)*dsig)

C - ptj2
      call filld('ptj2',        ptj2,damping_vec(20)*dsig)

      end

c$$$      subroutine pwhgfinalopshist
c$$$      implicit none
c$$$      include     'pwhg_bookhist-multi.h'
c$$$      end

