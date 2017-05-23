      subroutine init_hist
      implicit none
      real*8  dyh,dpth,dHT,dptj1!,dy,dpt,dptzm,dptzmzm

      call inihists

cccccccccccccccccccccccccccccccccccccccccccccccccc
c     HNNLO-like plots
c     Plots as those produced with HNNLO, using
c     the HNNLO kt-algo.
cccccccccccccccccccccccccccccccccccccccccccccccccc
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

ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Cross-check plots for jets, using FastJet
ccccccccccccccccccccccccccccccccccccccccccccccccc
      call bookupeqbins('ptj1_check',   dptj1,  0d0,300d0)
      call bookupeqbins('ptj2_check',   dptj1,  0d0,300d0)

      end

     
      subroutine analysis(dsig0)
      implicit none
      include 'hepevt.h'
      include 'pwhg_weights.h'
      include 'pwhg_bookhist-multi.h'
c     allow multiweights 
      real * 8 dsig0,dsig(1:weights_max)
c     tell to the analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
c     Initialisation indicator (for new hisogramming setup)
      logical ini
      data ini/.true./
      save ini

ccccccccccccccccccccccccc
c     HNNLO-like plots
ccccccccccccccccccccccccc
C - Maximum number of particles to send to jet finder: 
      include  'mxpart.f'
C - Higgs kinematics
      real*8  hmass
      common/hmass/hmass
      real*8  ph(4),yh,etah,pth,mh
C - Jet radius: 
      real*8  Rmin !:
      common/Rmin/Rmin
C - HNNLO jet finding:
      real*8  p(bigmxpart,4)
      real*8  HT,ptj1,ptj2
      real*8  damping_vec(20) !:
C - Auxiliary counters:
      integer ixx,njets

ccccccccccccccccccccccccccc
c     Cross-check plots for jets, using FastJet
ccccccccccccccccccccccccccc
      integer maxjet
      parameter (maxjet=2048)
      real *8 palg,rr,ptmin
      real *8  ktj(maxjet),etaj(maxjet),rapj(maxjet),
     1    phij(maxjet),pj(4,maxjet),ptrel(4)
      integer mjets
 
      if (ini) then
         if(hmass.le.0d0) then
            write(*,*) 'pwhg_analysis-release.f:'
            write(*,*) 'hmass = ',hmass
            write(*,*) 'Have not got acceptable Higgs mass'
            stop
         endif
         write(*,*) '*****************************'
         if(WHCPRG.ne.'LHE   ') then
            write(*,*) 'pwhg_analysis-release has to be'
            write(*,*) 'used only at the LHE file level.'
            write(*,*) 'Quitting ...'
            stop
         else
            write(*,*) '       LHE ANALYSIS'
         endif
         write(*,'(a,f10.5)') 
     C        ' pwhg_analysis: Higgs mass (hmass) value: ', hmass
         write(*,'(a,f6.4)') 
     C        ' pwhg_analysis: Using kt-algo with R = ',Rmin
         write(*,*) '*****************************'
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

cccccccccccccccccccccc
c     HNNLO-like plots
cccccccccccccccccccccc
c     Settings in jet_finder.f: ptjetmin,etajetmin,etajetmax,ptbjetmin,etabjetmax

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

cccccccccccccccccccccccccccccc
c     find jets with FastJet
cccccccccccccccccccccccccccccc
c     palg: -1=antikt, +1=kt
      palg=1
      rr=Rmin
      ptmin=0d0 !:
      mjets=0
      ktj(:)=0
      etaj(:)=0
      rapj(:)=0
      call buildjetspalg(1,palg,rr,ptmin,mjets,ktj,etaj,rapj,phij,
     $     ptrel,pj)

      if(mjets.lt.1) then
         write(*,*) '>>> EVENT WITH NO JETS: ',ktj(:5)
         write(*,*) 'It cannot happen from HJ-MiNLO'
         stop
      endif
c      print*, ktj(1:2)

ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Cross-check plots for jets, using FastJet
ccccccccccccccccccccccccccccccccccccccccccccccccc
      call filld('ptj1_check',ktj(1),dsig)
      call filld('ptj2_check',ktj(2),dsig)

      end
      

      subroutine buildjetspalg(iflag,palg,rr,ptmin,mjets,kt,eta,rap,phi,
     $     ptrel,pjet)
c     arrays to reconstruct jets, radius parameter rr
      implicit none
      integer iflag,mjets
      real * 8  rr,ptmin,kt(*),eta(*),rap(*),
     1     phi(*),ptrel(3),pjet(4,*)
      include   'hepevt.h'
      include  'LesHouches.h'
      integer   maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real * 8  ptrack(4,maxtrack),pj(4,maxjet)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      integer   ntracks,njets
      integer   j,k,mu,i
      real * 8 r,palg,tmp
      logical islept
      external islept
      real * 8 vec(3),pjetin(0:3),pjetout(0:3),beta,
     $     ptrackin(0:3),ptrackout(0:3)
      real * 8 get_ptrel
      external get_ptrel
C - Initialize arrays and counters for output jets
      do j=1,maxtrack
         do mu=1,4
            ptrack(mu,j)=0d0
         enddo
         jetvec(j)=0
      enddo      
      ntracks=0
      do j=1,maxjet
         do mu=1,4
            pjet(mu,j)=0d0
            pj(mu,j)=0d0
         enddo
      enddo
      if(iflag.eq.1) then
C     - Extract final state particles to feed to jet finder
         do j=1,nhep
c all but the Higgs
            if (isthep(j).eq.1.and..not.idhep(j).eq.25) then
               if(ntracks.eq.maxtrack) then
                  write(*,*) 'analyze: need to increase maxtrack!'
                  write(*,*) 'ntracks: ',ntracks
                  stop
               endif
               ntracks=ntracks+1
               do mu=1,4
                  ptrack(mu,ntracks)=phep(mu,j)
               enddo
               itrackhep(ntracks)=j
            endif
         enddo
      else
         do j=1,nup
            if (istup(j).eq.1.and..not.islept(idup(j))) then
               if(ntracks.eq.maxtrack) then
                  write(*,*) 'analyze: need to increase maxtrack!'
                  write(*,*) 'ntracks: ',ntracks
                  stop
               endif
               ntracks=ntracks+1
               do mu=1,4
                  ptrack(mu,ntracks)=pup(mu,j)
               enddo
               itrackhep(ntracks)=j
            endif
         enddo
      endif
      if (ntracks.eq.0) then
         mjets=0
         return
      endif
C --------------------------------------------------------------------- C
C     R = 0.7   radius parameter
c palg=1 is standard kt, -1 is antikt
c      palg=1 !:
      r=rr
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,
     $                        jetvec)
      mjets=njets
      if(njets.eq.0) return
c check consistency
      do k=1,ntracks
         if(jetvec(k).gt.0) then
            do mu=1,4
               pj(mu,jetvec(k))=pj(mu,jetvec(k))+ptrack(mu,k)
            enddo
         endif
c         print*, 'FJ: track ',k, 'is in jet ',jetvec(k)
      enddo
      tmp=0
      do j=1,mjets
         do mu=1,4
            tmp=tmp+abs(pj(mu,j)-pjet(mu,j))
         enddo
      enddo
      if(tmp.gt.1d-4) then
         write(*,*) ' bug!'
      endif
C --------------------------------------------------------------------- C
C - Computing arrays of useful kinematics quantities for hardest jets - C
C --------------------------------------------------------------------- C
      do j=1,mjets
         call getyetaptmass(pjet(:,j),rap(j),eta(j),kt(j),tmp)
         phi(j)=atan2(pjet(2,j),pjet(1,j))
      enddo

c     loop over the hardest 3 jets
      do j=1,min(njets,3)
         do mu=1,3
            pjetin(mu) = pjet(mu,j)
         enddo
         pjetin(0) = pjet(4,j)         
         vec(1)=0d0
         vec(2)=0d0
         vec(3)=1d0
         beta = -pjet(3,j)/pjet(4,j)
         call mboost(1,vec,beta,pjetin,pjetout)         
c     write(*,*) pjetout
         ptrel(j) = 0
         do i=1,ntracks
            if (jetvec(i).eq.j) then
               do mu=1,3
                  ptrackin(mu) = ptrack(mu,i)
               enddo
               ptrackin(0) = ptrack(4,i)
               call mboost(1,vec,beta,ptrackin,ptrackout) 
               ptrel(j) = ptrel(j) + get_ptrel(ptrackout,pjetout)
            endif
         enddo
      enddo
      end
