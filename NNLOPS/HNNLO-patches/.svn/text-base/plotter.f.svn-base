      subroutine bookplot(n,tag,titlex,var,wt,xmin,xmax,dx,llplot) 
      implicit none
      include 'nplot.f'
      integer n,higgsdec
      character*(*) titlex
      character*3 llplot
      character*4 tag,mypart
      double precision var,wt,xmin,xmax,dx
      common/mypart/mypart

      if (tag.eq.'book') then
          call mbook(n,titlex,dx,xmin,xmax)
          call mbook(20+n,titlex,dx,xmin,xmax)
          call mbook(40+n,titlex,dx,xmin,xmax)
          call mbook(60+n,titlex,dx,xmin,xmax)
          call mbook(80+n,titlex,dx,xmin,xmax)
          call mbook(100+n,titlex,dx,xmin,xmax)
          call mbook(120+n,titlex,dx,xmin,xmax)
      elseif (tag .eq. 'plot') then
          call mfill(n,var,wt)
        linlog(n)=llplot
        titlearray(n)=titlex
      endif

      return
      end

    
      subroutine plotter(p,wt,switch)
      implicit none
      include 'clustering.f'
      include 'constants.f'
      include 'cutoff.f'
      include 'jetlabel.f'
      include 'npart.f'
      include 'mxdim.f'
      include 'process.f'
      include 'removebr.f'
      include 'masses.f'

      integer n,switch,nplotmax
      character tag*4
  
      double precision dot
      double precision wt
      double precision p(mxpart,4)
      double precision pt,yraptwo,pttwo
      integer eventpart

      logical first,jetmerge
      common/nplotmax/nplotmax
      common/jetmerge/jetmerge
      integer order,higgsdec,ndec
      common/higgsdec/higgsdec,ndec
      common/nnlo/order
      data first/.true./
      save first

      double precision yh,pth
      double precision HT,ptj1,ptj2,pjm
      double precision HT_exp,pth_exp,ptj1_exp
      double precision HT_min,pth_min,ptj1_min
      double precision HT_hfact(4)
      double precision pth_hfact(4)
      double precision ptj1_hfact(4)
      double precision damping

      wt=wt/1000d0  !  Convert to Powheg-Box units.

      if (first) then
         tag='book'
c--- ensure we initialize all possible histograms
         eventpart=npart+3
         
         yh=0d0
         pth=0d0
         HT=0d0
         ptj1=0d0
         ptj2=0d0
         
         jetmerge=.true.
CC      Here set jets to the maximum number of jets 
CC      to book the necessary histograms: 0 at LO, 1 at NLO and 2 at NNLO
         jets=order
CC
         goto 99
      else
         tag='plot'
      endif

C - Damping factors.
C - These are generically defined as
C - damping=(X_hfact-Xmin)^X_exp/((X_hfact-Xmin)^X_exp+(X-Xmin)^X_exp)
C - where X refers to some kinematic variable e.g. HT, pth, ptj1,
C - which typically increases with increasing hardness of radiation
C - in the event. Xmin is the minimum value kinematically accessible
C - to X. So damping --->1 for soft events and damping ---> 0 for
C - hard events. 
      pth_exp=2d0
      pth_min=0d0
      pth_hfact(1)=0.5*hmass
      pth_hfact(2)=0.7*hmass
      pth_hfact(3)=1.0*hmass
      pth_hfact(4)=1.4*hmass

      HT_exp=2d0
      HT_min=hmass
      HT_hfact(1)=1.5*hmass
      HT_hfact(2)=2.0*hmass
      HT_hfact(3)=2.5*hmass
      HT_hfact(4)=3.0*hmass

      ptj1_exp=2d0
      ptj1_min=0d0
      ptj1_hfact(1)=0.5*hmass
      ptj1_hfact(2)=0.7*hmass
      ptj1_hfact(3)=1.0*hmass
      ptj1_hfact(4)=1.4*hmass

C     eventpart is the total number of four momenta in the event
C     2 for initial state + ndec for the H decay + number of jets 

      eventpart=2+ndec+jets
    

      yh=yraptwo(3,4,p)
      pth=pttwo(3,4,p)
      HT=sqrt(pth*pth+2*dot(p,3,4))

 99   continue

C - Order the jets
      if(jets.eq.2) then
         ptj1=pt(5,p)
         ptj2=pt(6,p)
         HT=HT+ptj1+ptj2
         if(ptj1.lt.ptj2) then
            pjm=ptj2
            ptj2=ptj1
            ptj1=pjm
         endif
      elseif(jets.eq.1) then
         ptj1=pt(5,p)
         ptj2=0d0
         HT=HT+ptj1
      elseif(jets.eq.0) then
         ptj1=0d0
         ptj2=0d0
      endif
      
      n=1                  
C - yh
      call bookplot(n,tag,'yh',yh,wt,-5d0,5d0,0.1d0,'lin')
      n=n+1

C - yh with pth damping function
      damping=((pth_hfact(1)-pth_min)**pth_exp)
     $       /((pth_hfact(1)-pth_min)**pth_exp
     $        +(pth-pth_min)**pth_exp)
      call bookplot(n,tag,'yh-pth-1',yh,
     $              damping*wt,-5d0,5d0,0.1d0,'lin')
      n=n+1
      damping=((pth_hfact(2)-pth_min)**pth_exp)
     $       /((pth_hfact(2)-pth_min)**pth_exp
     $        +(pth-pth_min)**pth_exp)
      call bookplot(n,tag,'yh-pth-2',yh,
     $              damping*wt,-5d0,5d0,0.1d0,'lin')
      n=n+1
      damping=((pth_hfact(3)-pth_min)**pth_exp)
     $       /((pth_hfact(3)-pth_min)**pth_exp
     $        +(pth-pth_min)**pth_exp)
      call bookplot(n,tag,'yh-pth-3',yh,
     $              damping*wt,-5d0,5d0,0.1d0,'lin')
      n=n+1
      damping=((pth_hfact(4)-pth_min)**pth_exp)
     $       /((pth_hfact(4)-pth_min)**pth_exp
     $        +(pth-pth_min)**pth_exp)
      call bookplot(n,tag,'yh-pth-4',yh,
     $              damping*wt,-5d0,5d0,0.1d0,'lin')
      n=n+1

C - yh with HT damping function
      damping=((HT_hfact(1)-HT_min)**HT_exp)
     $       /((HT_hfact(1)-HT_min)**HT_exp
     $        +(HT-HT_min)**HT_exp)
      call bookplot(n,tag,'yh-HT-1',yh,
     $              damping*wt,-5d0,5d0,0.1d0,'lin')
      n=n+1
      damping=((HT_hfact(2)-HT_min)**HT_exp)
     $       /((HT_hfact(2)-HT_min)**HT_exp
     $        +(HT-HT_min)**HT_exp)
      call bookplot(n,tag,'yh-HT-2',yh,
     $              damping*wt,-5d0,5d0,0.1d0,'lin')
      n=n+1
      damping=((HT_hfact(3)-HT_min)**HT_exp)
     $       /((HT_hfact(3)-HT_min)**HT_exp
     $        +(HT-HT_min)**HT_exp)
      call bookplot(n,tag,'yh-HT-3',yh,
     $              damping*wt,-5d0,5d0,0.1d0,'lin')
      n=n+1
      damping=((HT_hfact(4)-HT_min)**HT_exp)
     $       /((HT_hfact(4)-HT_min)**HT_exp
     $        +(HT-HT_min)**HT_exp)
      call bookplot(n,tag,'yh-HT-4',yh,
     $              damping*wt,-5d0,5d0,0.1d0,'lin')
      n=n+1

C - yh with ptj1 damping function
      damping=((ptj1_hfact(1)-ptj1_min)**ptj1_exp)
     $       /((ptj1_hfact(1)-ptj1_min)**ptj1_exp
     $        +(ptj1-ptj1_min)**ptj1_exp)
      call bookplot(n,tag,'yh-ptj1-1',yh,
     $              damping*wt,-5d0,5d0,0.1d0,'lin')
      n=n+1
      damping=((ptj1_hfact(2)-ptj1_min)**ptj1_exp)
     $       /((ptj1_hfact(2)-ptj1_min)**ptj1_exp
     $        +(ptj1-ptj1_min)**ptj1_exp)
      call bookplot(n,tag,'yh-ptj1-2',yh,
     $              damping*wt,-5d0,5d0,0.1d0,'lin')
      n=n+1
      damping=((ptj1_hfact(3)-ptj1_min)**ptj1_exp)
     $       /((ptj1_hfact(3)-ptj1_min)**ptj1_exp
     $        +(ptj1-ptj1_min)**ptj1_exp)
      call bookplot(n,tag,'yh-ptj1-3',yh,
     $              damping*wt,-5d0,5d0,0.1d0,'lin')
      n=n+1
      damping=((ptj1_hfact(4)-ptj1_min)**ptj1_exp)
     $       /((ptj1_hfact(4)-ptj1_min)**ptj1_exp
     $        +(ptj1-ptj1_min)**ptj1_exp)
      call bookplot(n,tag,'yh-ptj1-4',yh,
     $              damping*wt,-5d0,5d0,0.1d0,'lin')
      n=n+1

C - pth
      call bookplot(n,tag,'pth',pth,wt,0d0,300d0,3d0,'lin')
      n=n+1
C - pth with canonical damping factor
      damping=((pth_hfact(3)-pth_min)**pth_exp)
     $       /((pth_hfact(3)-pth_min)**pth_exp
     $        +(pth-pth_min)**pth_exp)
      call bookplot(n,tag,'pth-hfact-3',pth,
     $              damping*wt,0d0,300d0,3d0,'lin')
      n=n+1
C - HT
      call bookplot(n,tag,'HT',HT,wt,100d0,600d0,5d0,'lin')
      n=n+1
C - HT with canonical damping factor
      damping=((HT_hfact(2)-HT_min)**HT_exp)
     $       /((HT_hfact(2)-HT_min)**HT_exp
     $        +(HT-HT_min)**HT_exp)
      call bookplot(n,tag,'HT-hfact-2',HT,
     $              damping*wt,100d0,600d0,5d0,'lin')
      n=n+1
C - ptj1
      call bookplot(n,tag,'ptj1',ptj1,wt,0d0,300d0,3d0,'lin')
      n=n+1
C - ptj1 with canonical damping factor
      damping=((ptj1_hfact(3)-ptj1_min)**ptj1_exp)
     $       /((ptj1_hfact(3)-ptj1_min)**ptj1_exp
     $        +(ptj1-ptj1_min)**ptj1_exp)
      call bookplot(n,tag,'ptj1-hfact-3',ptj1,
     $              damping*wt,0d0,300d0,3d0,'lin')
      n=n+1
C - ptj2
      call bookplot(n,tag,'ptj2',ptj2,wt,0d0,300d0,3d0,'lin')
      n=n+1


      n=n-1

      if (n.gt.20) then
        write(6,*) 'WARNING - TOO MANY HISTOGRAMS!'
        write(6,*) n,' > 20, which is the built-in maximum'
        stop
      endif

c--- set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif
      

      return 
      end
      

      
