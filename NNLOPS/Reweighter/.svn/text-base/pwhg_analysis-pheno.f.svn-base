      subroutine init_hist
      implicit none
      real * 8 dy,dpt,dptzoom
      integer icut
      real *8 jetcut(2)
      character * 6 suffix(2) 
      common/jcut/jetcut,suffix

      call inihists

      dy=0.1d0
      dpt=5d0
      dptzoom=1d0
      suffix(1) = '-ptj20'
      jetcut(1) = 20d0 
      suffix(2) = '-ptj30'
      jetcut(2) = 30d0 


      call bookupeqbins('total',1d0,0d0,1d0)
      call bookupeqbins('yH',dy,-8d0,8d0)

      call bookupeqbins('ptH',dpt,0d0,500d0)
      call bookupeqbins('ptHzoom',dptzoom,0d0,50d0)

c      call bookupeqbins('ptj1_check',3d0,0d0,300d0)
      call bookupeqbins('ptj1',dpt,0d0,500d0)
      call bookupeqbins('ptj1zoom',dptzoom,0d0,50d0)

c      call bookupeqbins('ptj2_check',3d0,0d0,300d0)
      call bookupeqbins('ptj2',dpt,0d0,500d0)
      call bookupeqbins('ptj2zoom',dptzoom,0d0,50d0)

C     -- loop over two sets of jet cuts (20 and 30 GeV) 
      do icut=1,2
C     - Inclusive jet multiplicties
         call bookupeqbins('Njet-inc'//suffix(icut),1d0,-0.5d0,5.5d0)
C     - Exclusive jet multiplicties
         call bookupeqbins('Njet-exc'//suffix(icut),1d0,-0.5d0,5.5d0)
      enddo
         
      end

      
     
      subroutine analysis(dsig0)
      implicit none
      include 'hepevt.h'
      include 'pwhg_weights.h'
      logical ini
      data ini/.true./
      save ini
C     allow multiweights 
      real * 8 dsig0,dsig(1:weights_max)
c     tell to the analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      real *8 jetcut(2)
      character * 6 suffix(2) 
      common/jcut/jetcut,suffix
C     Higgs variables 
      real*8 pH(4),mH,ptH,yH,etaH
c     arrays to reconstruct jets
      integer maxjet
      parameter (maxjet=2048)
      real *8 ptmin
      real *8  ktj(maxjet),etaj(maxjet),rapj(maxjet),
     1    phij(maxjet),pj(4,maxjet),rr,ptrel(4)
      integer j1,found,mjets,ihep,icut,i
      
 
      if (ini) then
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

c     Loop over final state particles to find Higgs
      found=0
      do ihep=1,nhep
         if (((isthep(ihep).eq.1).or.(isthep(ihep).eq.2)
     #.or.(isthep(ihep).eq.155).or.(isthep(ihep).eq.195))
     #.and.(idhep(ihep).eq.25)) then
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
      call filld('ptH',ptH,dsig)
      call filld('ptHzoom',ptH,dsig)

C     find jets 
      rr=0.4d0       
      ptmin=0d0
      call buildjets(1,rr,ptmin,mjets,ktj,etaj,rapj,phij,ptrel,pj)

      if(mjets.ge.1) then
c         call filld('ptj1_check',ktj(1),dsig)
         call filld('ptj1',ktj(1),dsig)
         call filld('ptj1zoom',ktj(1),dsig)
      endif
      if(mjets.ge.2) then
c         call filld('ptj2_check',ktj(2),dsig)
         call filld('ptj2',ktj(2),dsig)
         call filld('ptj2zoom',ktj(2),dsig)
      endif
      
c     following plots have not been tested carefully...

C     -- Plot inclusive and exclusive jet multiplicties at two jet-pT thresholds
      do icut = 1,2


c     no-jets OR hardest-jet softer than threshold: 
c     ->fill 0-th bin of Njet-inc and 0-th bin of Njet-exc
         if((mjets.eq.0) .or. 
     $        (mjets.gt.0 .and. ktj(1).lt.jetcut(icut))) then
            call filld('Njet-inc'//suffix(icut),0d0,dsig)
            call filld('Njet-exc'//suffix(icut),0d0,dsig)
         endif

c     at least 1 jet: loop over jets
         do i=1,5
c     ->fill inclusive i-th bin if i-th jet is harder than threshold
            if(mjets.ge.i.and.ktj(i).ge.jetcut(icut)) 
     $           call filld('Njet-inc'//suffix(icut),i*1d0,dsig)
c     ->fill exclusive i-th bin if 
c     1) I have exactly i jets and i-th jet is harder than threshold
c     OR
c     2) I have more than i jets, i-th jet is harder than threshold
c        but (i+1)-jet (and all the others) is too soft
            if( (mjets.eq.i.and.ktj(i).ge.jetcut(icut)) .or. 
     $          (mjets.gt.i.and.ktj(i).ge.jetcut(icut) .and. 
     $           ktj(i+1).lt.jetcut(icut)))
     $           call filld('Njet-exc'//suffix(icut),i*1d0,dsig)
         enddo

      enddo
 
      end
      


      subroutine buildjets(iflag,rr,ptmin,mjets,kt,eta,rap,phi,
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
      palg=-1
      r=rr
c      ptmin=20d0 
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
