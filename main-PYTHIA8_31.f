      program main_pythia8
      implicit none
      include 'LesHouches.h'
      include 'hepevt.h'
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      integer j,k,l,m,iret
      integer maxev
      common/mcmaxev/maxev
      real * 8 powheginput,scalupfac
      external powheginput

      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix

      logical analysis_jetveto
      common/canalysis_jetveto/analysis_jetveto

      character *20 eventfile

      WHCPRG='PYTHIA'

      call init_hist

      call getmaxev(maxev)

c      eventfile=pwgprefix(1:lprefix)//'events.lhe'
c      eventfile="w-events.lhe"

c      maxev=49000
c      maxev=1000
c      maxev=50000
c      maxev=949000
c      maxev=100

      call lhefreadhdr(97)

c      call pythia_init(trim(eventfile(1:)))
      call pythia_init

      nevhep=0

      do l=1,maxev

         call lhefreadev(97)

         nevhep=nevhep+1

         do m=1,1 !ER:
c Insist to shower this event;
            call pythia_next(iret)

            if(iret.eq.-1) then
               nevhep=nevhep-1
               print*, 'EOF'
               print*, nevhep
               print*, l
               goto 123
            endif
           
            if(iret.ne.1) then
               write(*,*) ' return code ',iret
               if(m.eq.1) then
                  write(*,*) ' Pythia could not shower this event'
                  call printleshouches
               endif
c               write(*,*) ' retry: ',l
            else
c               write(*,*) ' done: ',l
               exit
            endif
         enddo

         if(iret.eq.1) then
            call pythia_to_hepevt(nmxhep,nhep,isthep,idhep,jmohep,
     1           jdahep,phep,vhep)
            if(nevhep.lt.2) then
               do j=1,nhep
                  write(*,100)j,isthep(j),idhep(j),jmohep(1,j),
     1           jmohep(2,j),jdahep(1,j),jdahep(2,j), (phep(k,j),k=1,5)
               enddo
            endif

            call pyanal

            if(nevhep.gt.0.and.mod(nevhep,20000).eq.0) then
               write(*,*)'# of events processed=',nevhep
               write(*,*)'# of events generated=',nevhep
               call pyaend
            endif 
            if(nevhep.gt.0.and.mod(nevhep,20000).eq.0) then
               write(*,*)'# of events processed=',nevhep
               write(*,*)'# of events generated=',nevhep
            endif 
         endif


      enddo

 123  continue

      write(*,*) 'At the end NEVHEP is ',nevhep

      call pythia_stat
!:      write(*,*) 'At the end: #warnings= ',mstu(27),' #errors= ',mstu(23)
c---user's terminal calculations


      call pyaend

 100  format(i4,2x,i5,2x,i5,2x,i4,1x,i4,2x,i4,1x,i4,2x,5(d10.4,1x))
      end


      subroutine pyanal
      implicit none
      include 'LesHouches.h'
      include 'hepevt.h'
c     check parameters
      logical verbose
      parameter (verbose=.false.)
      real * 8 powheginput
      external powheginput
c     !ER:
c      nevhep=nevhep+1
c     print*, xwgtup,xsecup(1),nup
      if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
      call analysis(xwgtup)
      call pwhgaccumup 
      end

      subroutine getmaxev(maxev)
      integer maxev
C--- Opens input file and counts number of events, setting MAXEV;
      call opencount(maxev)
      end

      subroutine pyaend
      character * 6 vetoname
      character * 20 pwgprefix
      character * 100 filename
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      include 'pwhg_rnd.h'
      logical analysis_jetveto
      common/canalysis_jetveto/analysis_jetveto
      
      if(analysis_jetveto) then
         vetoname=trim(adjustl('_veto'))
      else
         vetoname=trim(adjustl(''))
      endif
      if(rnd_cwhichseed.ne.'none') then
         filename=pwgprefix(1:lprefix)//'POWHEG+PYTHIA8-output-'
     1        //rnd_cwhichseed//vetoname
      else
         filename=pwgprefix(1:lprefix)//'POWHEG+PYTHIA8-output'//vetoname
      endif

      call pwhgsetout
      call pwhgtopout(filename)
      end

      subroutine printleshouches
c useful for debugging
      call lhefwritev(6)
      end

c...lhefeader(nlf)
c...writes event information to a les houches events file on unit nlf. 
      subroutine lhefwritev(nlf)
      implicit none
      integer nlf
      include 'LesHouches.h'
      include 'pwhg_flg.h'
      integer i,j
      write(nlf,'(a)')'<event>'
      write(nlf,210) nup,idprup,xwgtup,scalup,aqedup,aqcdup
      do 200 i=1,nup
         write(nlf,220) idup(i),istup(i),mothup(1,i),
     & mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     & vtimup(i),spinup(i)
 200  continue
      write(nlf,'(a)')'</event>'      
 210  format(1p,2(1x,i8),4(1x,e12.5))
 220  format(1p,i8,5(1x,i5),5(1x,e16.9),1x,e12.5,1x,e10.3)
      end
