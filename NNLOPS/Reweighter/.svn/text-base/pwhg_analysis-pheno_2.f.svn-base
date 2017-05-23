c     This analysis should replace pwhg_analysis-HNNLO.f
c     for all the runs involving POWHEG.
c     Besides results with the kt-algo, it also 
c     outputs jets spectra obtained with the antikt algorithm.
c     Plots 1-20 are the same, and have the same name,
c     as those produced when running HNNLO.
c     It also includes a check that for results at NLO or LH
c     level, results are exactly the same as what one obtains
c     using the original pwhg_analysis-HNNLO.f.
c     This can also be seen by inspecting the plots with 'check' 
c     in their name.
c     When pwhg_analysis-HNNLO.f is used on showered results,
c     there are small differences between the kt jet-algo 
c     implemented in HNNLO files and the implementation if FastJet.
c     In particular, in FJ
c     dR=sqrt(dy^2+dphi^2)
c     whereas in the MCFM/HNNLO implementation
c     dR=sqrt(deta^2+dphi^2)
c     To reproduce exactly the FJ results using the MCFM
c     implementation, changes should be made in the file r.f.
      subroutine init_hist
      implicit none
      real*8  dyh,dpth,dHT,dptj1,dy,dpt,dptzm,dptzmzm
      integer ialg,iRcut,ipTcut
      real*8  pTcut(2)
      real*8  Rcut(4)
      integer alg(2)
      character*4 algsuffix(2)
      character*4 Rsuffix(4)
      character*6 pTsuffix(2)
      common/pTcut/pTcut,Rcut,alg,algsuffix,Rsuffix,pTsuffix

      call inihists

cccccccccccccccccccccccccccccccccccccccccccccccccc
c     HNNLO-like plots
c     Plots as those produced with HNNLO, using
c     the HNNLO kt-algo.
c     This should produce results identical to plots
c     produced with POWHEG until July 24th.
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
      call bookupeqbins('ptj1_check',3d0,0d0,300d0)
      call bookupeqbins('ptj2_check',3d0,0d0,300d0)

ccccccccccccccccccccccccccccccc
c     Other plots
ccccccccccccccccccccccccccccccc

      dy=0.1d0
      dpt=5d0
      dptzm=3d0
      dptzmzm=1d0

      call bookupeqbins('total',1d0,0d0,1d0)

      call bookupeqbins('yH',dy,-8d0,8d0)
      call bookupeqbins('ptHzm',dptzm,0d0,300d0)
      call bookupeqbins('ptHzmzm',dptzmzm,0d0,300d0)

      algsuffix(1) = '--kt'
      alg(1) =  1
      algsuffix(2) = '-akt'
      alg(2) = -1

      pTsuffix(1) = '-ptj20'
      pTcut(1) = 20d0 
      pTsuffix(2) = '-ptj30'
      pTcut(2) = 30d0 

      Rsuffix(1) = '-r04'
      Rcut(1)    = 0.4d0
      Rsuffix(2) = '-r05'
      Rcut(2)    = 0.5d0
      Rsuffix(3) = '-r07'
      Rcut(3)    = 0.7d0
      Rsuffix(4) = '-r10'
      Rcut(4)    = 1.0d0

      do ialg=1,2
         do iRcut=1,4

            call bookupeqbins('ptJ1zm'//algsuffix(ialg)
     $                               //Rsuffix(iRcut),dptzm,0d0,300d0)
            call bookupeqbins('ptJ1zmzm'//algsuffix(ialg)
     $                               //Rsuffix(iRcut),dptzmzm,0d0,300d0)
            call bookupeqbins('ptJ2zm'//algsuffix(ialg)
     $                               //Rsuffix(iRcut),dptzm,0d0,300d0)
            call bookupeqbins('ptJ2zmzm'//algsuffix(ialg)
     $                               //Rsuffix(iRcut),dptzmzm,0d0,300d0)

C     -- loop over two sets of jet cuts (20 and 30 GeV) 
            do ipTcut=1,2
               call bookupeqbins('yJ1'//algsuffix(ialg)
     $                   //pTsuffix(ipTcut)//Rsuffix(iRcut),dy,-8d0,8d0)
               call bookupeqbins('yJ2'//algsuffix(ialg)
     $                   //pTsuffix(ipTcut)//Rsuffix(iRcut),dy,-8d0,8d0)
C     - Inclusive jet multiplicties
               call bookupeqbins('Njet-inc'//algsuffix(ialg)
     $              //pTsuffix(ipTcut)//Rsuffix(iRcut),1d0,-0.5d0,5.5d0)
C     - Exclusive jet multiplicties
               call bookupeqbins('Njet-exc'//algsuffix(ialg)
     $              //pTsuffix(ipTcut)//Rsuffix(iRcut),1d0,-0.5d0,5.5d0)
            enddo
         enddo
      enddo

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
      data    Rmin/0.7d0/
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

ccccccccccccccccccccccccccc
c     Other plots
ccccccccccccccccccccccccccc
      real*8  pTcut(2)
      real*8  Rcut(4)
      integer alg(2)
      character*4 algsuffix(2)
      character*4 Rsuffix(4)
      character*6 pTsuffix(2)
      common/pTcut/pTcut,Rcut,alg,algsuffix,Rsuffix,pTsuffix
c     arrays to reconstruct jets
      integer maxjet
      parameter (maxjet=2048)
      real *8 ptmin,palg
      real *8  ktj(maxjet),etaj(maxjet),rapj(maxjet),
     1    phij(maxjet),pj(4,maxjet),rr,ptrel(4)
      integer j1,found,mjets,ihep,icut
      integer ialg,iRcut,ipTcut
      integer jxx

      integer  indexhist
      external indexhist

      logical debug
      parameter (debug=.false.)
      real *8 Etot,p1(4),p2(4),dy,deta,dphi,dr
 
      if (ini) then
         hmass = 125.5d0  
         write(*,*) 'Higgs mass set to', hmass 
         if(hmass.le.0d0) then
            write(*,*) 'pwhg_analysis-HNNLO.f:'
            write(*,*) 'hmass = ',hmass
            write(*,*) 'Have not got acceptable Higgs mass'
            stop
         endif
         write(*,*) '*****************************'
         if(whcprg.eq.'NLO') then
            write(*,*) '       NLO ANALYSIS'
            weights_num=0
         elseif(WHCPRG.eq.'LHE   ') then
            write(*,*) '       LHE ANALYSIS'
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS            '
         elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS            '
         endif
c$$$         write(*,*) 'kt-algo run with R= ',Rmin
c$$$         write(*,*) 'antikt-algo run with R= ',Rmin_akt
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
c     Be careful on settings in jet_finder.f,
c     in particular
c     ptjetmin
c     and
c     etajetmin,etajetmax,ptbjetmin,etabjetmax
c     Also check Rmin above!

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
      do ialg=1,2

         palg=alg(ialg)

         do iRcut=1,4
ccccccc      rr=Rmin
            rr=Rcut(iRcut)
            ptmin=0d0  !:
            mjets=0
            ktj(:)=0
            etaj(:)=0
            rapj(:)=0
            call buildjetspalg(1,palg,rr,ptmin,mjets,ktj,etaj,rapj,phij,
     $                         ptrel,pj)

c - Start 'if(alg.eq.1.and.Rcut(iRcut).eq.0.4d0) ...' block 
            if((palg.eq.1).and.(Rcut(iRcut).eq.0.4d0)) then

ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Cross-check plots for jets, using FastJet
ccccccccccccccccccccccccccccccccccccccccccccccccc
               if(debug) then
                  if((mjets.ge.1).and.((dabs(ktj(1)/ptj1-1).gt.0.00001)
     $                 .or.(dabs(ktj(2)/ptj2-1).gt.0.00001))) then
c     if(mjets.ne.njets) then
                     Etot=0d0
                     print*, '-----------'
                     print*, 'FJ vs HNNLO: njets ',mjets,njets
                     print*, 'FJ vs HNNLO: ptj1 ',ktj(1),ptj1
                     print*, 'FJ vs HNNLO: ptj2 ',ktj(2),ptj2
                     print*, 'FJ jets: '
                     do ixx=1,mjets
                        print*, ktj(ixx),' | ',pj(:,ixx)
                        if(ktj(ixx).gt.0) Etot=Etot+pj(4,ixx)
                     enddo
                     print*, 'FJ jets: total E ',Etot
                     stop
                  endif
               endif
c$$$      call getdydetadphidr(p1,p2,dy,deta,dphi,dr)
               call filld('ptj1_check',ktj(1),dsig)
               call filld('ptj2_check',ktj(2),dsig)

ccccccccccccccccccccccc
c     Other plots
ccccccccccccccccccccccc
c     Loop over final state particles to find Higgs
               found=0
               do ihep=1,nhep
                  if (((isthep(ihep).eq.1).or.(isthep(ihep).eq.2)
     $               .or.(isthep(ihep).eq.155).or.(isthep(ihep).eq.195))
     $               .and.(idhep(ihep).eq.25)) then
                     j1=ihep
                     found=found+1
                     pH(1:4) = phep(1:4,j1)
                  endif
               enddo
               if(found.lt.1) then
                  write(*,*) 'ERROR: Higgs not found'
                  call exit(1)
               elseif(found.gt.1) then
                  write(*,*) 'ERROR: more Higgs-like particles found'
                  call exit(1)
               endif

c     total sigma
               call filld('total',0.5d0,dsig)

c     Higgs histograms 
               call getyetaptmass(pH,yH,etaH,ptH,mH)
               call filld('yH',yH,dsig)
               call filld('ptHzm',ptH,dsig)
               call filld('ptHzmzm',ptH,dsig)

c - End of 'if(alg.eq.1.and.Rcut(irad).eq.0.4d0) ...' block 
            endif

            if(mjets.lt.1) then
               print*, '>>> EVENT WITH NO JETS, mjets....',ktj(:5)
c               stop
            endif
            if(mjets.ge.1) then
               call filld('ptJ1zm'
     $              //algsuffix(ialg)//Rsuffix(iRcut),ktj(1),dsig)
               call filld('ptJ1zmzm'
     $              //algsuffix(ialg)//Rsuffix(iRcut),ktj(1),dsig)
            endif
            if(mjets.ge.2) then
               call filld('ptJ2zm'
     $              //algsuffix(ialg)//Rsuffix(iRcut),ktj(2),dsig)
               call filld('ptJ2zmzm'
     $              //algsuffix(ialg)//Rsuffix(iRcut),ktj(2),dsig)
            endif


            do ipTcut = 1,2
C - Count number of jets above current pT threshold (assumes ktj array
C - is ordered from highest to lowest pT):
               mjets=0
               do ixx=1,10 ! mjets stops at 10 jets i.e. mjets=10 really means
                           ! there were 10 or more jets above the pT threshold.
                  if(ktj(ixx).gt.pTcut(ipTcut)) mjets=mjets+1
               enddo

C - Plot jet rapidities:
               if(mjets.ge.1) call filld('yJ1'//algsuffix(ialg)
     $                  //pTsuffix(ipTcut)//Rsuffix(iRcut),rapj(1),dsig)
               if(mjets.ge.2) call filld('yJ2'//algsuffix(ialg)
     $                  //pTsuffix(ipTcut)//Rsuffix(iRcut),rapj(2),dsig)

C - Plot inclusive jet multiplicities:
               do ixx=0,mjets
                  call filld('Njet-inc'
     $                        //algsuffix(ialg)
     $                        //pTsuffix(ipTcut)//Rsuffix(iRcut),
     $                        ixx*1d0,dsig)
               enddo

C - Plot exclusive jet multiplicities:
               call filld('Njet-exc'
     $                     //algsuffix(ialg)
     $                     //pTsuffix(ipTcut)//Rsuffix(iRcut),
     $                     mjets*1d0,dsig)

            enddo

         enddo

      enddo
 
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
