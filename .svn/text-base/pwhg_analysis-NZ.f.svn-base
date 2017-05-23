c This analysis was proposed as benchmark for the Higgs workshop at
c the end of 2012, by Bruce Mellado and others.
c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'

      call inihists

C - Total cross-section 
      call bookupeqbins('sigmatot', 1d0,0d0,1d0)
      call bookupeqbins('H-y',0.1d0,-5d0,5d0)
      call bookupeqbins('H-pt',3d0,0d0,300d0)
      call bookupeqbins('Yj1',0.1d0,-5d0,5d0)
      call bookupeqbins('ptj1',3d0,0d0,300d0)

C - Inclusive boson pt in pseudorapidity bins:
      call bookupeqbins('H-pt--inf-yH--004',3d0,0d0,300d0)
      call bookupeqbins('H-pt--004-yH--003',3d0,0d0,300d0)
      call bookupeqbins('H-pt--003-yH--002',3d0,0d0,300d0)
      call bookupeqbins('H-pt--002-yH--001',3d0,0d0,300d0)
      call bookupeqbins('H-pt--001-yH-+000',3d0,0d0,300d0)
      call bookupeqbins('H-pt-+000-yH-+001',3d0,0d0,300d0)
      call bookupeqbins('H-pt-+001-yH-+002',3d0,0d0,300d0)
      call bookupeqbins('H-pt-+002-yH-+003',3d0,0d0,300d0)
      call bookupeqbins('H-pt-+003-yH-+004',3d0,0d0,300d0)
      call bookupeqbins('H-pt-+004-yH-+inf',3d0,0d0,300d0)

C - Inclusive boson pt in pseudorapidity bins (cumulative):
      call bookupeqbins('H-ptint',3d0,0d0,300d0)
      call bookupeqbins('H-ptint--inf-yH--004',3d0,0d0,300d0)
      call bookupeqbins('H-ptint--004-yH--003',3d0,0d0,300d0)
      call bookupeqbins('H-ptint--003-yH--002',3d0,0d0,300d0)
      call bookupeqbins('H-ptint--002-yH--001',3d0,0d0,300d0)
      call bookupeqbins('H-ptint--001-yH-+000',3d0,0d0,300d0)
      call bookupeqbins('H-ptint-+000-yH-+001',3d0,0d0,300d0)
      call bookupeqbins('H-ptint-+001-yH-+002',3d0,0d0,300d0)
      call bookupeqbins('H-ptint-+002-yH-+003',3d0,0d0,300d0)
      call bookupeqbins('H-ptint-+003-yH-+004',3d0,0d0,300d0)
      call bookupeqbins('H-ptint-+004-yH-+inf',3d0,0d0,300d0)

      end
    
      subroutine analysis(dsig0)
      implicit  none
      real * 8     dsig0,dsig(20)
      include     'hepevt.h'
      include     'nlegborn.h'
      include     'pwhg_kn.h'
      include     'pwhg_math.h' 
      include     'pwhg_weights.h'
      include     'pwhg_lhrwgt.h'
      integer      maxjet
      parameter   (maxjet=2048)
      character*6 WHCPRG
      common/cWHCPRG/WHCPRG
      data    WHCPRG/'NLO   '/
      integer      mjets,njets
      real * 8       ktj(maxjet),etaj(maxjet),rapj(maxjet),
     1             phij(maxjet),pj(4,maxjet),jetRadius,ptrel(4)

      real * 8       ph(4)
      real * 8       y,eta,pt,m
      integer      ihep,ixx,nBosons
      real * 8       powheginput
      external     powheginput
      real * 8     low_pt_bin_edge,ptmin

      logical ini
      data    ini/.true./
      save    ini

      if(ini) then
         if(lhrwgt_nids.gt.0) then
            call setupmulti(lhrwgt_nids)
         elseif(weights_num.gt.0) then
            call setupmulti(weights_num+1)
         else
            call setupmulti(1)
         endif
         ini=.false.
      endif

      dsig=0
      if(lhrwgt_nids.gt.0) then
         dsig(1:lhrwgt_nids)=lhrwgt_weights(1:lhrwgt_nids)
      elseif(weights_num.gt.0) then
         dsig(1)=dsig0
         dsig(2:weights_num+1)=weights_val(1:weights_num)
      else
         dsig(1)=dsig0
      endif

      if(sum(abs(dsig)).eq.0) return



C - Find the boson:
      nBosons = 0
C - For HERWIG or PYTHIA, the last H boson is kept, which is likely the
C - Higgs before decaying, if it decays, or the final Higgs, if undecayed.
      if (WHCPRG == 'NLO') then 
         ph = phep(1:4,3)
         nBosons = 1
      else
         do ihep=1,nhep
            if(idhep(ihep).eq.25) then
               ph = phep(1:4,ihep)
               nBosons = 1
            endif
         enddo
         if(nBosons.eq.0) then
            write(*,*) 'pwhg_analysis-NZ.f: Fatal error.'
            write(*,*) ' Higgs boson not found'
            call exit(-1)
         endif
      endif


C - Boson kinematics fully inclusive w.r.t radiation:
      call getyetaptmass(ph,y,eta,pt,m)
      call filld('sigmatot',0.5d0,dsig)
      call filld("H-pt",pt,dsig)
      call filld("H-y",y,dsig)

      jetRadius= 0.4d0       
      ptmin = 1d0

      call buildjets(1,jetRadius,ptmin,mjets,ktj,etaj,rapj,phij,ptrel,pj)

      if(mjets.ge.1) then
         if(ktj(1).gt.25) then
            call filld("Yj1",rapj(1),dsig)
         endif
         call filld("ptj1",ktj(1),dsig)
      endif

C - Inclusive boson pt in pseudorapidity bins:
      if(y.lt.-4d0) then
         call filld("H-pt--inf-yH--004",
     $        pt,dsig)
      elseif(y.gt.-4d0.and.y.lt.-3d0) then
         call filld("H-pt--004-yH--003",
     $        pt,dsig)
      elseif(y.gt.-3d0.and.y.lt.-2d0) then
         call filld("H-pt--003-yH--002",
     $        pt,dsig)
      elseif(y.gt.-2d0.and.y.lt.-1d0) then
         call filld("H-pt--002-yH--001",
     $        pt,dsig)
      elseif(y.gt.-1d0.and.y.lt. 0d0) then
         call filld("H-pt--001-yH-+000",
     $        pt,dsig)
      elseif(y.gt. 0d0.and.y.lt. 1d0) then
         call filld("H-pt-+000-yH-+001",
     $        pt,dsig)
      elseif(y.gt. 1d0.and.y.lt. 2d0) then
         call filld("H-pt-+001-yH-+002",
     $        pt,dsig)
      elseif(y.gt. 2d0.and.y.lt. 3d0) then
         call filld("H-pt-+002-yH-+003",
     $        pt,dsig)
      elseif(y.gt. 3d0.and.y.lt. 4d0) then
         call filld("H-pt-+003-yH-+004",
     $        pt,dsig)
      elseif(y.gt. 4d0) then
         call filld("H-pt-+004-yH-+inf",
     $        pt,dsig)
      endif

C - Inclusive boson pt in rapidity bins (cumulative):
      do ixx=1,100
         low_pt_bin_edge = (ixx-1d0)*3d0
         if(pt.gt.low_pt_bin_edge) then
            call filld("H-ptint",
     $             low_pt_bin_edge+1d-3,dsig)
           if(y.lt.-4d0) then
             call filld("H-ptint--inf-yH--004",
     $             low_pt_bin_edge+1d-3,dsig)
           elseif(y.gt.-4d0.and.y.lt.-3d0) then
             call filld("H-ptint--004-yH--003",
     $             low_pt_bin_edge+1d-3,dsig)
           elseif(y.gt.-3d0.and.y.lt.-2d0) then
             call filld("H-ptint--003-yH--002",
     $             low_pt_bin_edge+1d-3,dsig)
           elseif(y.gt.-2d0.and.y.lt.-1d0) then
             call filld("H-ptint--002-yH--001",
     $             low_pt_bin_edge+1d-3,dsig)
           elseif(y.gt.-1d0.and.y.lt. 0d0) then
             call filld("H-ptint--001-yH-+000",
     $             low_pt_bin_edge+1d-3,dsig)
           elseif(y.gt. 0d0.and.y.lt. 1d0) then
             call filld("H-ptint-+000-yH-+001",
     $             low_pt_bin_edge+1d-3,dsig)
           elseif(y.gt. 1d0.and.y.lt. 2d0) then
             call filld("H-ptint-+001-yH-+002",
     $             low_pt_bin_edge+1d-3,dsig)
           elseif(y.gt. 2d0.and.y.lt. 3d0) then
             call filld("H-ptint-+002-yH-+003",
     $             low_pt_bin_edge+1d-3,dsig)
           elseif(y.gt. 3d0.and.y.lt. 4d0) then
             call filld("H-ptint-+003-yH-+004",
     $             low_pt_bin_edge+1d-3,dsig)
           elseif(y.gt. 4d0) then
             call filld("H-ptint-+004-yH-+inf",
     $             low_pt_bin_edge+1d-3,dsig)
           endif
         endif
      enddo



      end

 






c      subroutine yetaptmassplot(p,dsig,prefix)
c      implicit none
c      real * 8 p(4),dsig(*)
c      character *(*) prefix
c      real * 8 y,eta,pt,m
c      call getyetaptmass(p,y,eta,pt,m)
c      call filld(prefix//'-y',y,dsig)
c      call filld(prefix//'-eta',eta,dsig)
c      call filld(prefix//'-pt',pt,dsig)
c      call filld(prefix//'-m',m,dsig)
c      end

      subroutine deltaplot(p1,p2,dsig,prefix,postfix)
      implicit none
      real * 8 p1(4),p2(4),dsig(*)
      character *(*) prefix,postfix
      real * 8 dy,deta,delphi,dr
      call getdydetadphidr(p1,p2,dy,deta,delphi,dr)
      call filld(prefix//'-dy'//postfix,dy,dsig)
      call filld(prefix//'-deta'//postfix,deta,dsig)
      call filld(prefix//'-delphi'//postfix,delphi,dsig)
      call filld(prefix//'-dr'//postfix,dr,dsig)
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

      subroutine getdydetadphidr(p1,p2,dy,deta,dphi,dr)
      implicit none
      include 'pwhg_math.h' 
      real * 8 p1(*),p2(*),dy,deta,dphi,dr
      real * 8 y1,eta1,pt1,mass1,phi1
      real * 8 y2,eta2,pt2,mass2,phi2
      call getyetaptmass(p1,y1,eta1,pt1,mass1)
      call getyetaptmass(p2,y2,eta2,pt2,mass2)
      dy=y1-y2
      deta=eta1-eta2
      phi1=atan2(p1(1),p1(2))
      phi2=atan2(p2(1),p2(2))
      dphi=abs(phi1-phi2)
      dphi=min(dphi,2d0*pi-dphi)
      dr=sqrt(deta**2+dphi**2)
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

      subroutine sortbypt(n,iarr)
      implicit none
      integer n,iarr(n)
      include 'hepevt.h'
      integer j,k
      real * 8 tmp,pt(nmxhep)
      logical touched
      do j=1,n
         pt(j)=sqrt(phep(1,iarr(j))**2+phep(2,iarr(j))**2)
      enddo
c bubble sort
      touched=.true.
      do while(touched)
         touched=.false.
         do j=1,n-1
            if(pt(j).lt.pt(j+1)) then
               k=iarr(j)
               iarr(j)=iarr(j+1)
               iarr(j+1)=k
               tmp=pt(j)
               pt(j)=pt(j+1)
               pt(j+1)=tmp
               touched=.true.
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

      subroutine pwhgfinalopshist
      end


