c A trial analysis to test out the HJ generator

c Initialise histograms
      subroutine init_hist
      implicit none
      include 'LesHouches.h'
      include 'pwhg_math.h'
      integer 			 ixx,jxx,kxx,len_pt1,len_pt2,len_y,t1,t2,ls1
C     integer maxbins1,maxbins2
C     parameter (maxbins1=30,maxbins2=8)
      real*8           pt_bins1(300),pt_bins2(300),y_bins1(300)
      character * 60   s1(100)
      character * 60   tmp,tmp2
      common/bins/pt_bins1,y_bins1,len_pt1,len_pt2,len_y
      integer          lenocc
      external         lenocc
      logical          negative ! check whether then bin edges (which will become strings for naming the plots) is negative
      call inihists

c     total plots
      call bookupeqbins('sigmatot',1d0,-1d0,2d0)
      call bookupeqbins('hpt-5GeV',5d0,0d0,300d0)
      call bookupeqbins('hpt-2GeV',2d0,0d0,400d0)

      call bookupeqbins('hpt-5GeV-y-lt-0.4',5d0,0d0,300d0)
      call bookupeqbins('hpt-5GeV-0.4-y-1.2',5d0,0d0,300d0)
      call bookupeqbins('hpt-5GeV-1.2-y-2.4',5d0,0d0,300d0)
      call bookupeqbins('hpt-5GeV-2.4-y-4.0',5d0,0d0,300d0)

      call bookupeqbins('hpt-2GeV-y-lt-0.4',2d0,0d0,400d0)
      call bookupeqbins('hpt-2GeV-0.4-y-1.2',2d0,0d0,400d0)
      call bookupeqbins('hpt-2GeV-1.2-y-2.4',2d0,0d0,400d0)
      call bookupeqbins('hpt-2GeV-2.4-y-4.0',2d0,0d0,400d0)


      end
      
c     Analysis subroutine      
      subroutine analysis(dsig0)
      implicit none
      include 'hepevt.h'
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'pwhg_weights.h'
      include 'pwhg_lhrwgt.h'
      real * 8     dsig0,dsig(7)
      integer      nweights
      logical      iniwgts
      data         iniwgts/.true./
      save         iniwgts
      character*6  WHCPRG
      common/cWHCPRG/WHCPRG
      data         WHCPRG/'NLO   '/
      integer      ihep,ixx,i_higgs,i_j1,i_j2,id,kxx
      integer      mjets,njets,maxjet
      parameter   (maxjet=2048)
      logical      IsForClustering(maxjet)
      real * 8     ktj(maxjet),etaj(maxjet),rapj(maxjet),
     1     phij(maxjet),pj(4,maxjet),jetRadius,ptrel(4),ptmin
      real * 8     ph(4)
      real * 8     ptj1,ptj2,pt_test,pt_higgs,y_higgs
      real * 8     y,eta,pt,m
      real * 8     powheginput
      real * 8 	 pt_bins1(300),y_bins1(300),pt_bins2(300)
      integer      len_pt1,len_pt2,len_y,lenocc,t1,t2,ls1
      character*60 s1(100),tmp,tmp2
      common/bins/ pt_bins1,y_bins1,len_pt1,len_pt2,len_y
      external     powheginput,lenocc
      logical      negative
      
      if (iniwgts) then
         write(*,*) '*********************'
         if(whcprg.eq.'NLO    ') then
            write(*,*) ' NLO ANALYSIS      '
            weights_num=0
         elseif(WHCPRG.eq.'LHE    ') then
            write(*,*) ' LHE ANALYSIS      '
         elseif(WHCPRG.eq.'HERWIG ') then
            write(*,*) ' HERWIG ANALYSIS   '
         elseif(WHCPRG.eq.'PYTHIA ') then
            write(*,*) ' PYTHIA ANALYSIS   '
         elseif(WHCPRG.eq.'PYTHIA8') then
            write(*,*) ' PYTHIA8 ANALYSIS   '
         endif
         write(*,*) '*********************'
         if(weights_num.eq.0) then
            call setupmulti(1)
         else
            call setupmulti(weights_num)
         endif
         iniwgts=.false.
      endif
      
      dsig=0
      if(weights_num.eq.0) then
         dsig(1)=dsig0
      else
         dsig(1:weights_num)=weights_val(1:weights_num)
      endif
      if(sum(abs(dsig)).eq.0) return
      
c     Initialise everything to zero
      i_higgs=0 		! index of the higgs
      do ixx=1,100
         s1(ixx) = ''
      enddo
      do ixx=1,4
         ph(:)  =0.0
      enddo

c     Find the Higgs - we will choose the final higgs
      do ihep=1,nhep
         id=abs(idhep(ihep))
         if(idhep(ihep).eq.25) i_higgs = ihep
      enddo

      ph=phep(1:4,i_higgs)
      
c     Call Fastjet to build jets
      jetRadius= 0.4d0       
      ptmin = 1d0
      call buildjets(1,jetRadius,ptmin,mjets,ktj,
     1        etaj,rapj,phij,ptrel,pj) 	
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                        c
c                        Make plots                      c
c                                                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c    Total plots
      call getyetaptmass(ph,y,eta,pt,m)
      pt_higgs=pt
      y_higgs=y
      call filld('sigmatot',0.5d0,dsig)
      call filld('hpt-5GeV',pt_higgs,dsig)
      call filld('hpt-2GeV',pt_higgs,dsig)

      y_higgs=abs(y_higgs)

      if(y_higgs.lt.0.4) then
         call filld('hpt-5GeV-y-lt-0.4',pt_higgs,dsig)
         call filld('hpt-2GeV-y-lt-0.4',pt_higgs,dsig)
      elseif(y_higgs.ge.0.4.and.y_higgs.lt.1.2) then
         call filld('hpt-5GeV-0.4-y-1.2',pt_higgs,dsig)
         call filld('hpt-2GeV-0.4-y-1.2',pt_higgs,dsig)
      elseif(y_higgs.ge.1.2.and.y_higgs.lt.2.4) then
         call filld('hpt-5GeV-1.2-y-2.4',pt_higgs,dsig)
         call filld('hpt-2GeV-1.2-y-2.4',pt_higgs,dsig)
      elseif(y_higgs.ge.2.4.and.y_higgs.lt.4.0) then
         call filld('hpt-5GeV-2.4-y-4.0',pt_higgs,dsig)
         call filld('hpt-2GeV-2.4-y-4.0',pt_higgs,dsig)
      endif

      end
      
      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(4),y,eta,pt,mass,pv
      real *8 tiny
      parameter (tiny=1.d-5)
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      pt=sqrt(p(1)**2+p(2)**2)
      pv=sqrt(pt**2+p(3)**2)
      if(pt.lt.tiny)then
         eta=sign(1.d0,p(3))*1.d8
      else
         eta=0.5d0*log((pv+p(3))/(pv-p(3)))
      endif
      mass=sqrt(abs(p(4)**2-pv**2))
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
      integer   j,k,mu,jb,i
      real * 8 r,palg,pp,tmp
      logical islept
      external islept
      real * 8 vec(3),pjetin(0:3),pjetout(0:3),beta,
     $     ptrackin(0:3),ptrackout(0:3)
      real * 8 get_ptrel
      external get_ptrel
C     - Initialize arrays and counters for output jets
      do j=1,maxtrack
         do mu=1,4
            ptrack(mu,j)=0d0
         enddo
         jetvec(j)=0
         ptrel(j) = 0d0
      enddo      
      ntracks=0
      do j=1,maxjet
         do mu=1,4
            pjet(mu,j)=0d0
            pj(mu,j)=0d0
         enddo
         kt(j)  = 0d0
         eta(j) = 0d0
         rap(j) = 0d0
         phi(j) = 0d0
      enddo
      if(iflag.eq.1) then
C     - Extract final state particles to feed to jet finder
         do j=1,nhep
c     all but the Higgs
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
C     --------------------------------------------------------------------- C
C     R = 0.7   radius parameter
c     palg=1 is standard kt, -1 is antikt
      palg=-1
      r=rr
c     ptmin=20d0 
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,
     $     jetvec)
      mjets=njets
      if(njets.eq.0) return
c     check consistency
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
C     --------------------------------------------------------------------- C
C     - Computing arrays of useful kinematics quantities for hardest jets - C
C     --------------------------------------------------------------------- C
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
      
      function islept(j)
      implicit none
      logical islept
      integer j
      if(abs(j).ge.11.and.abs(j).le.15) then
         islept = .true.
      else
         islept = .false.
      endif
      end
      
      function get_ptrel(pin,pjet)
      implicit none
      real * 8 get_ptrel,pin(0:3),pjet(0:3)
      real * 8 pin2,pjet2,cth2,scalprod
      pin2  = pin(1)**2 + pin(2)**2 + pin(3)**2
      pjet2 = pjet(1)**2 + pjet(2)**2 + pjet(3)**2
      scalprod = pin(1)*pjet(1) + pin(2)*pjet(2) + pin(3)*pjet(3)
      cth2 = scalprod**2/pin2/pjet2
      get_ptrel = sqrt(pin2*abs(1d0 - cth2))
      end

      
