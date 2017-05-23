C ************************************************************** C
C ************************************************************** C
C ************************************************************** C
C Usage: e.g. 
C ./minnlo <lhe file> <string> <hnnlo output> [<index (1:13)>]
C input arguments are: 
C 1) file to be rescaled
C 2) list of HNNLO output files, each one preceded by an identifying string
C 3) index of distribution used for rescaling (not mandatory)
C
C >>> We recommend to use 10 as value for 'index'. This corresponds
C >>> to use beta=0.5 in eq. 2.9 of 1309.0017. If no value is
C >>> entered, the default is used.
C
C >>> Variables defined locally in this file: Rmin, prescription.
C >>> We recommend to leave them as they are (Rmin=0.7, prescription=2)
C >>> 
C >>> More details can be found in the following.
C ************************************************************** C
C ************************************************************** C
C ************************************************************** C
C
C
C Rmin is used as jet radius to generate the plots which we 
C extract the denumerator weights from. In other words, it 
C is used to analyze the LH files from HJ-MiNLO.
C >>>
C >>> Therefore, it MUST match the jet radius used in HNNLO.
C >>> Look at all the occurences of !:
C >>>
C 
C prescription is the switch controlling the W(y,pt) reweighing
C function we use.
C prescription=1 corresponds to eq. 2.10 in 1309.0017
C prescription=2 corresponds to eq. 2.12 in 1309.0017
C i.e.
C 1: W = h*(NNLO_A/MINLO_A)       + (1-h)
C 2: W = h*(NNLO-MINLO_B)/MINLO_A + (1-h)
C >>>
C >>> The suggested value, which was used in 1309.0017, is 2.
C >>>
      program minnlo
      implicit none
      include 'pwhg_weights.h'
      include 'LesHouches.h'
      include 'pwhg_bookhist-multi.h'
      include 'maxrwgtfiles.h' 
      include 'mxpart.f'

      character*6 WHCPRG
      common/cWHCPRG/WHCPRG

      real*8  hmass
      common/hmass/hmass

      integer ios,nlf,nlo,nev,iev

      character*(1000) lhFile
      integer nFiles
C     maxgwgtfiles: max number of "input" files for reweighting
      character*(1000) rwgtFiles(maxRwgtFiles)
      character*(20) rwgtstring(maxRwgtFiles)
      integer dataIndex
      character*(2) cdataIndex

C     nhist: max number of distributions 
C     maxbins: max number of bins per distribution 
      integer  nDataSets(maxRwgtFiles)
      integer  nLines(maxRwgtFiles,nhist)
      character * 100 dataSetTitles(maxRwgtFiles,nhist)
      real * 8 Xlo(maxRwgtFiles,nhist,maxbins)
      real * 8 Xhi(maxRwgtFiles,nhist,maxbins)
      real * 8 num_Y(maxRwgtFiles,nhist,maxbins)
      real * 8 num_Yerr(maxRwgtFiles,nhist,maxbins)

      real * 8 den_Y(maxmulti,nhist,maxbins)
      real * 8 den_Yerr(maxmulti,nhist,maxbins)
     
      real*8  Rmin !:
      data    Rmin/0.7d0/
      common/Rmin/Rmin
      character*(3) cRmin

      integer prescription
      data prescription/2/
      save prescription
      character*(1) cprescription

      real * 8 resc(maxRwgtFiles,maxmulti)
      real * 8 yH,pth,HT,ptj1,ptj2,hdamp
      real * 8 etah,ph(4),mh
      real * 8 damping_vec(20)  !:
      real * 8 p(bigmxpart,4)
      integer  ibin,iline

      integer i,j,njets

      character*(100) tmpLine,tmpLine2

      character*(1000) lhFile2
      integer lprefix

      logical addnewweightline
      real * 8 renscfact,facscfact
      integer lhans1,lhans2
      character*(3) PDFp
      common/creweightinfo/renscfact,facscfact,lhans1,lhans2,
     $     addnewweightline,PDFp


C     Which level are we reweighting at ('LHE   ' / 'PYTHIA' / 'HERWIG')?
      WHCPRG='LHE   '

C     What jet radius do the jets have in the RWGT denominator / numerator files?
c$$$      Rmin=0.7d0 !:
      write(cRmin,'(f3.1)') Rmin
      cRmin=trim(adjustl(cRmin))

C     Which reweighting prescription are we using 
c     ( 1: h*(NNLO_A/MINLO_A)       + (1-h) )
c     ( 2: h*(NNLO-MINLO_B)/MINLO_A + (1-h) )
c$$$      prescription=2 !:
      write(cprescription,'(i1)') prescription
      cprescription=trim(adjustl(cprescription))

      write(*,*) ''
      write(*,*) '*************************************'
      write(*,*) '*                                   *'
      write(*,*) '* minnlo NNLO reweighting:          *'
      write(*,*) '* ------------------------          *'
      write(*,*) '*                                   *'
      write(*,*) '* Reweighting will be carried out   *'
      write(*,*) '* at the ',WHCPRG,' level.              *'
      write(*,*) '*                                   *'
      write(*,*) '* Reweighting functions profiled    *'
      write(*,*) '* w.r.t jet pTs in input data for   *'
      write(*,*) '* weight numerator and denominators *'
      write(*,*) '* are assumed to be clustered with  *'
      write(*,'(a,f9.6,a)') 
     $     ' * radius Rmin = ',Rmin,'           *'
      write(*,'(a,i2,a)') 
     $     ' * Reweighting prescription used:',prescription,'  *'
      write(*,*) '*                                   *'
      write(*,*) '*************************************'
      write(*,*) ''

C     Get the inputs from the command line and check for consistency
      call readCommandLine(lhFile,nFiles,rwgtFiles,rwgtstring,dataIndex)
      write(cdataIndex,'(i2)') dataIndex
      cdataindex=trim(adjustl(cdataIndex))

      if(dataIndex.ne.10) then
         write(*,*) 'Please select 10 as data index to be used, as'
         write(*,*) 'this is the recommended value in 1309.0017.'
         write(*,*) 'Other values should be used only by developers'
         write(*,*) 'Quitting ...'
         stop
      endif

C     Load in the input data
      call loadInputData(nFiles,rwgtFiles,
     $     nDataSets,nLines,dataSetTitles,
     $     Xlo,Xhi,num_Y,num_Yerr)

C     Open the event file and count how many events
      call opencount(lhFile,nev)
C     Check if the events have reweight info or not. 
C     If they do, count how many weights.
C     If they don't, remember this in addnewweightline.
      nlf = 10
      open(unit=nlf,file=trim(lhFile),status='old',iostat=ios)
      call lhefreadhdr(nlf)
c     read first event to get the number of weights & rewind 
      call lhefreadev(nlf)   
      rewind(nlf)
      if(weights_num.ne.0) then
         addnewweightline=.false.
         write(*,'(/,a,i3,a,/)') 'Found ',weights_num,' weights per '//
     $        'event in LH input file'
      else
         write(*,*) 'In ',trim(lhFile),' the "#new weight" line '//
     $        'is missing. Will use the standard weight (xwgtup).'
         addnewweightline=.true.
c     Extract from header the values associated to the keywords
c     'renscfact','facscfact','lhans1','lhans2','PDF package'.
c     If any of them is missing, use -1 or '?' as value.
c     If this happens, some information is missing in the
c     header copied from the input file. In general, this
c     should not happen.       
         renscfact=-1
         facscfact=-1
         lhans1=-1
         lhans2=-1
         PDFp=' ? '
         i=0
         do while(i.eq.0)
            read(nlf,fmt='(a)',iostat=ios) tmpLine
            do j=1,100
               if(tmpline(j:j).eq.'#'.or.tmpline(j:j).eq.'!') then
                  tmpline(j:)=' '
               endif
            enddo
            if(index(tmpLine,'</eve').ne.0) i=1
            if(index(tmpLine,'renscfact').ne.0) then
               tmpLine2 = tmpLine(index(trim(tmpLine),' ',back=.true.)
     $              :len(tmpLine))
               read(tmpLine2,*,iostat=ios) renscfact
            endif
            if(index(tmpLine,'facscfact').ne.0) then
               tmpLine2 = tmpLine(index(trim(tmpLine),' ',back=.true.)
     $              :len(tmpLine))
               read(tmpLine2,*,iostat=ios) facscfact
            endif
            if(index(tmpLine,'lhans1').ne.0) then
               tmpLine2 = tmpLine(index(trim(tmpLine),' ',back=.true.)
     $              :len(tmpLine))
               read(tmpLine2,*,iostat=ios) lhans1
            endif
            if(index(tmpLine,'lhans2').ne.0) then
               tmpLine2 = tmpLine(index(trim(tmpLine),' ',back=.true.)
     $              :len(tmpLine))
               read(tmpLine2,*,iostat=ios) lhans2
            endif
            if(index(tmpLine,'PDF package').ne.0) then
               tmpLine2 = tmpLine(index(trim(tmpLine),' ',back=.true.)
     $              :len(tmpLine))
               read(tmpLine2,*,iostat=ios) PDFp
            endif
         enddo
         rewind(nlf)
      endif

c     Find the value of hmass in the header of the event
c     file, and read its value. This is needed, because we 
c     will use this (fixed) value in the reweighting function,
c     and the same must have been used in HNNLO.
c     If this is not found, abort the program.
c     
C     Read in all lines until you find the first occurrence of hmass
      ios=0
      hmass=-1d0
      do while(ios.eq.0)
         read(nlf,fmt='(a)',iostat=ios) tmpLine
         if(index(tmpLine,'hmass').ne.0) then
            tmpLine2 = tmpLine(index(trim(tmpLine),' ',back=.true.)
     $                        :len(tmpLine))
            read(tmpLine2,*,iostat=ios) hmass
         endif
      enddo
      if(hmass.lt.0) then
         write(*,*) 'hmass token not found in LH input file'
         write(*,*) 'Quitting ...'
         stop
      else
         write(*,*) ' In LH input file header, found '//
     $        'hmass = ',hmass
      endif
      rewind(nlf)

cccccccccccccccccccccccccccc
C     Compute denominators 
cccccccccccccccccccccccccccc
      write(*,*) 'Starting step 1: computing denominators'
      call init_hist 
      do j=1,nev
         if(mod(j,10000) == 0) write(*,*) '  step 1: at event ',
     $        j,'/',nev
         call lhefreadev(nlf)   
         if(WHCPRG == 'LHE   ') then
            call lhuptohepevt(j)
            if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
            call analysis(xwgtup)
         else
c     here shower and then call analysis...
            write(*,*) 'denominators after showering'//
     $           'not implemented. Program stops.'
            stop
         endif
         call pwhgaccumup
      enddo
c     finish denominator plots
      call lheanend
      close(nlf) 

C     Now we have the denominators in yhistarr2 and errhistarr2
C     just swap indices to keep things simple 
      do i=1,nhist
         do j=1,maxbins
            den_Y(:,i,j) = yhistarr2(:,j,i) 
            den_Yerr(:,i,j) = errhistarr2(:,j,i) 
         enddo
      enddo

C     now open event file again to read events 
C     and do the rescaling 
      nlf = 10
      open(unit=nlf,file=lhFile,status='old',iostat=ios)
C     open also output file (will contain rescaled weights)
      nlo = 11

C      lhFile2=trim(lhfile)//'-fn'//trim(cdataIndex)//
C     $     '-p'//cprescription//
C     $     '-R'//cRmin//
C     $     '.nnlo'
C
C for general distribution use a simplified name 
      lhFile2=trim(lhfile)//'.nnlo'


      open(unit=nlo,file=trim(lhFile2),status='unknown',iostat=ios)
C     read and copy header
      call lhefreadwritehdr(nlf,nlo)

c     If we have a single NLO weight and no 'new weight' line,
c     then weights_num is 0, whereas the index to be used 
c     for the denominator plots is 1 [do j=1,weights_num]
      if(addnewweightline) weights_num=1

ccccccccccccccccc
C     main loop
ccccccccccccccccc
      write(*,*) ''
      write(*,*) 'Starting step 2: computing NNLO rescaling'
      do iev=1,nev
         if(mod(iev,10000) == 0) write(*,*) '  step 2: at event ',
     $        iev,'/',nev
C     read and copy event 
         call lhefreadwriteev(nlf,nlo) 

C     fill hepevt common block, neeed by find_jets
         call lhuptohepevt(iev)

C     find jets 
         call find_jets(Rmin,p,njets) ! p are output jet momenta

C     get Higgs and jet kinematics
         ph(1:4) = pup(1:4,3)
         call getyetaptmass(ph,yh,etah,pth,mh)
         HT=sqrt(pth*pth+mh*mh)
         call get_ptj1j2(p,ptj1,ptj2)
c$$$         if(nup.eq.4.and.njets.eq.1) print*,njets,pth,ptj1,ptj2
         HT = HT+ptj1+ptj2 ! NB HT includes only two hardest jets

C     Compute damping factor associated to the weight 
c     profile fn (dataindex) depending on pTH / pTJ1 / HT ...
         call get_hdamp(yH,pth,HT,ptj1,ptj2,damping_vec)

C     Find bin associated to current value of Higgs rapidity
         ibin = 0 
         do iline=1,nlines(1,dataindex) 
            if ( (Xlo(1,dataindex,iline) < yH) .and. 
     $           (yH <= Xhi(1,dataindex,iline))) then 
               ibin = iline 
               exit 
            endif
         enddo
         if (ibin == 0) then 
            write(*,*) 'yH out of range' 
         endif

C     Compute the rescaling matrix.
C     dataindex set the distribution to be used for the rescaling.

C     num_y(:,1,:)                        ==> (dsig_NNLO)
C     num_y(:,dataindex,:)                ==> (dsig_NNLO)_A

C     den_y(:,1,:)                        ==> (dsig_MINLO)
C     den_y(:,dataindex,:)                ==> (dsig_MINLO)_A
C     (den_y(:,1,:)-den_y(:,dataindex,:)) ==> (dsig_MINLO)_B

         resc = 0d0
         hdamp=damping_vec(dataindex)
         do i=1,nfiles          ! nfiles: number of HHNLO input files 
            do j=1,weights_num  ! number of weights in the event file 
               if(prescription.eq.1) then
                  resc(i,j) = 
     $                 num_y(i,dataindex,ibin)/den_y(j,dataindex,ibin)
     $                 * hdamp 
     $                 + 1d0-hdamp
               elseif(prescription.eq.2) then
c$$$                  print*, den_y(j,1,ibin)-den_y(j,dataindex,ibin)
                  resc(i,j) = 
     $                 (  num_y(i,1,ibin)-
     $                    (den_y(j,1,ibin)-den_y(j,dataindex,ibin)) 
     $                 )
     $                 /den_y(j,dataindex,ibin)
     $                 * hdamp 
     $                 + 1d0-hdamp
               else
                  write(*,*) 'reweighting prescription '//
     $                 'not implemented yet.'
                  write(*,*) 'Quitting ...'
                  stop
               endif
c$$$               if ((iev == 1) .and. (i == 1)) print*, 'ev 1 resc',
c$$$     $              j,ibin,den_y(j,dataindex,ibin) 
            enddo
         enddo
         
C     write event extra info with rescaled weights 
         call lhefreadwriteextrannlo(resc,nfiles,rwgtstring,nlf,nlo)
      enddo

C     close all units 
      close(nlf) 
      close(nlo) 

C     If ever needed, here add some sanity-check counters

C     write how the reweighted file is called
      write(*,*) 'LH file with nnlo weights: ',trim(lhFile2)
      end





C ************************************************************** C
C ************************************************************** C
C ************************************************************** C
C - Loops over all numerator files, then within that loop,
C - it further loops over all data sets (histograms),
C - loading them into the returned 3-D arrays
C - Xlo,Xhi,num_Y,num_Yerr,den_Y,den_Yerr
C - i.e. ALL data in all numerator files is loaded (in case down the line
C - we might want to compute reweightings for all different possible
C - weights that the data files allow).
      subroutine loadInputData(nFiles,rwgtFiles,
     $                         nDataSets,nLines,dataSetTitles,
     $                         Xlo,Xhi,num_Y,num_Yerr)
      implicit none
      include 'pwhg_bookhist-multi.h'
      include 'maxrwgtfiles.h'

      integer nFiles
      character*(1000) rwgtFiles(maxRwgtFiles)

      integer  nDataSets(maxRwgtFiles)
      integer  nLines(maxRwgtFiles,nhist)
      character * 100 dataSetTitles(maxRwgtFiles,nhist)
      real * 8 Xlo(maxRwgtFiles,nhist,maxbins)
      real * 8 Xhi(maxRwgtFiles,nhist,maxbins)
      real * 8 num_Y(maxRwgtFiles,nhist,maxbins)
      real * 8 num_Yerr(maxRwgtFiles,nhist,maxbins)


      integer jxx,ixx

      integer ios
      integer nDataSetsNum
      character*(100) tmpLine
      integer nLinesNum

      character * 100 tmp_title
      real * 8 tmp_Xlo(maxbins),tmp_Xhi(maxbins)
      real * 8 tmp_Y(maxbins),tmp_Yerr(maxbins)

      ndatasets = 0 

C - Read in all data in 'numerator' reweighting files
      do ixx=1,nFiles

         nDataSets(ixx)=-1
C - Count number of data sets in the numerator file
         ios=0
         nDataSetsNum=0
         open(unit=11,file=rwgtFiles(ixx),status='old')
         do while(ios.eq.0)
            read(11,fmt='(a)',iostat=ios) tmpLine
            if(index(tmpLine,'#').ne.0) nDataSetsNum=nDataSetsNum+1
         enddo
         close(11)

         nDataSets(ixx) = nDataSetsNum

C - Load in each data set for the current numerator file
         do jxx=1,nDataSetsNum
            nLinesNum=-1
            call loadDataSet(rwgtFiles(ixx),jxx-1,
     $                       nLinesNum,tmp_title,
     $                       tmp_Xlo,tmp_Xhi,tmp_Y,tmp_Yerr)
            nLines(ixx,jxx)        = nLinesNum
            dataSetTitles(ixx,jxx) = tmp_title
            Xlo(ixx,jxx,:)         = tmp_Xlo
            Xhi(ixx,jxx,:)         = tmp_Xhi
            num_Y(ixx,jxx,:)       = tmp_Y
            num_Yerr(ixx,jxx,:)    = tmp_Yerr
         enddo

      enddo

      end


C ************************************************************** C
C ************************************************************** C
C ************************************************************** C
C - Does a well-defined thing: looks through file rwgtFile, locates
C - histogram with index dataIndex, reads in the data and fills in
C - the 1-D output arrays Xlo,Xhi,Y,Yerr, with the relevant values.
      subroutine loadDataSet(rwgtFile,dataIndex,
     $                       nLines,title,Xlo,Xhi,Y,Yerr)
      implicit none
      include 'pwhg_bookhist-multi.h'
      include 'maxrwgtfiles.h'


      character*(1000) rwgtFile
      integer  dataIndex
      integer  nLines
      character * 100 title
      real * 8 Xlo(maxbins),Xhi(maxbins)
      real * 8 Y(maxbins),Yerr(maxbins)

      integer ios

      character*(100) tmpLine,tmpLine2
      integer tmpIndex

      open(unit=11,file=rwgtFile,status='old')

C - Read in all lines until you get to an error or the end of the file
      ios=0
      nLines=0
      do while(ios.eq.0)
         read(11,fmt='(a)',iostat=ios,end=111) tmpLine
C - If you get to a line with a # extract the integer data index in it
         if(index(tmpLine,'#').ne.0) then
            tmpLine2 = tmpLine(index(trim(tmpLine),' ',back=.true.)
     $                        :len(tmpLine))
            read(tmpLine2,*,iostat=ios) tmpIndex
         endif
C - If we already read in some lines and the current data set index is
C - no longer equal to that requested we must have read in all we asked
C - for so we exit.
         if(tmpIndex.ne.dataIndex.and.nLines.gt.0) goto 111
C - If the index extracted is not equal to the data index requested go
C - back to the start of the loop, equally, if the line read is blank,
C - or if the line read contains #. The latter two if's mean that we
C - don't get past here inside the while loop unless we read actual data,
C - moreover, the data that we asked for.
         if(tmpIndex.ne.dataIndex)   goto 101
         if(trim(tmpLine).eq.'')     goto 101
         if(index(tmpLine,'#').ne.0) then
            title=trim(tmpLine)
            goto 101
         endif
C - Increment the line counter for every line of requested data read
         nLines=nLines+1
C - Fill the next element of the returned arrays from the tmpLine string
         read(tmpLine,'(4(1x,e14.8))')
     $        Xlo(nLines),Xhi(nLines),Y(nLines),Yerr(nLines)
 101     continue
      enddo

 111  return
      close(11)
      
      end



C ************************************************************** C
C ************************************************************** C
C ************************************************************** C
C - Does a well defined thing. Reads the command line looking
C - for the LH event file to be reweighted, returned to the
C - program as lhfile, the "data" file names to be used to
C - construct the weights, returned as 1-D array rwgtFiles,
C - the number of data files nFiles, and also the data (histogram)
C - index to be used for the reweighting.
      subroutine readCommandLine(lhFile,nFiles,rwgtFiles,rwgtstring,
     $     dataIndex)           ! all output
      implicit none
      include 'pwhg_bookhist-multi.h'
      include 'maxrwgtfiles.h'

      character*(1000) lhFile
      integer nFiles
      character*(1000) rwgtFiles(maxRwgtFiles)
      character*(20) rwgtstring(maxRwgtFiles)
      integer       dataIndex, nstrings

      integer ixx

      character*(3) cdataIndex

C - Get the event file name containing the events to be reweighted.
      call getarg(1,lhFile)
      lhFile=trim(lhFile)
      if(lhfile.eq.'') then
         write(*,'(/,a,/,a,/,a,/,a,/,/,a,/)')
     $        'minnlo:',
     $        'No numerator file given.',
     $        'Program should be run with the following inputs:',
     $        'minnlo <event file> '//
     $        ' <HNNLO-identifier1> <HNNLO-outputfile1> '//
     $        '[<HNNLO-identifier2> <HNNLO-outputfile2> ...]'//
     $        'Quitting ... '
         stop
      endif

C - Since the reweighting data index is not mandatory,
C - set it to the default value, that could be changed
C - on the command line (but can't be changed now)
      dataIndex=10

C - Get the names of the file(s) assumed to contain the numerator
C - information for constructing the NNLO weights.
      nfiles = 0; nstrings = 0  
      do ixx=1,2*maxRwgtFiles
         if (mod(ixx,2).eq.0) then  
            nfiles = nfiles+1
            call getarg(ixx+1,rwgtFiles(nfiles))
         else
            nstrings = nstrings+1
            call getarg(ixx+1,rwgtstring(nstrings))
         endif
         if(trim(rwgtFiles(nfiles)).eq.'') then
            nfiles = nfiles-1
            goto 101
         endif
         if(trim(rwgtstring(nstrings)).eq.'') then
            nstrings = nstrings-1
            goto 101
         endif
      enddo

 101  continue
      if (nfiles == 0) then 
         write(*,'(/,a,/,a,/,a,/,a,/,/,a,/)')
     $        'minnlo:',
     $        'No numerator file given.',
     $        'Program should be run with the following inputs:',
     $        'minnlo <event file> '//
     $        ' <HNNLO-identifier1> <HNNLO-outputfile1> '//
     $        '[<HNNLO-identifier2> <HNNLO-outputfile2> ...]'//
     $        'Quitting ... '
         stop 
      endif


      if(nstrings.eq.nfiles+1) then
C - If we enter here, then it means we might have a 
C - non-default dataindex in the command line.
C - This has been read as a character in 
C     rwgtstring(nstrings)
C - Get the index of the data array in the file(s) with which to
C - construct the weight(s).

         cdataIndex=trim(rwgtstring(nstrings))
         nstrings=nstrings-1
         read(cdataIndex,*) dataIndex
         if (dataIndex > 13 .or. dataindex < 1) 
     $        stop 'Dataindex out of allowed range. Check analysis' 
      endif
      write(6,'(/,a,I3)')
     $     'The data index to be used to reweight with is ',dataIndex

      if(nFiles.lt.1) then
         write(*,'(/,a,/,a,/,a,/,a,/,a,/,/,a,/,a/)')
     $        'minnlo:',
     $        'No reweighting string or data files and / or',
     $        'no reweighting data '//
     $        'index were supplied on the command line.',
     $        'Program should be run with the following inputs:',
     $        'minnlo <event file> '//
     $        '<identifier numerator> <wgt numerator> ...'//
     $        '[<reweighting data index>]',
     $        'Quitting ... '
         stop
      endif
      if(nFiles.ne.nstrings) then
         write(*,'(/,a,/,a,/,a,/,/,a,/)')
     $        'minnlo:',
     $        'Mismatch between number of numerator files',
     $        'and number of identifier string',
     $        'Quitting ... '
         print*, nfiles,nstrings
         stop
      endif


C - Print out information on what was found on the command line.
      write(*,*) ''
      write(6,'(/,a,a,a)')
     $     'Found LH event file ',trim(lhFile),' on the command line.'
      write(6,'(/,a,i3,a,/,a)')
     $     'Found a further ',nFiles,
     $     ' files assumed to contain numerator ',
     $     ' data for event reweighting.'
      write(*,*), ''
      do ixx=1,nFiles
         write(*,*) ixx,' ',trim(rwgtFiles(ixx))
      enddo
      write(*,*), ''

      end


      subroutine lhuptohepevt(n)
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      integer ihep,mu,n
      nhep=nup
      nevhep=n
      do ihep=1,nhep
         isthep(ihep)=istup(ihep)
         idhep(ihep)=idup(ihep)
         do mu=1,2
            jmohep(mu,ihep)=mothup(mu,ihep)
         enddo
         do mu=1,5
            phep(mu,ihep)=pup(mu,ihep)
         enddo
      enddo
      end


      subroutine lheanend
      implicit none
      character * 20 pwgprefix
      character * 1000 filename
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
cccccccccccc
      include 'maxrwgtfiles.h'
      character*(1000) lhFile
      integer nFiles
C     maxgwgtfiles: max number of "input" files for reweighting
      character*(1000) rwgtFiles(maxRwgtFiles)
      character*(20) rwgtstring(maxRwgtFiles)
      integer dataIndex
cccccccccccc
      include 'pwhg_rnd.h'
      if(rnd_cwhichseed.ne.'none') then
         filename=pwgprefix(1:lprefix)//'LHEF_analysis-'
     1        //rnd_cwhichseed
      else
         filename=pwgprefix(1:lprefix)//'LHEF_analysis'
      endif
      call pwhgsetout
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Uncommenting the following will produce files 
C     containing LH-level plots from the HJ-MiNLO sample 
C     (i.e. the 'denominator' plots).
c$$$C     We don't really need to write results in a file, since 
c$$$C     eveything is in common/histnew/ in pwhg_boookhist-multi.h.
c$$$C     However, it's a good sanity check.
c$$$      call readCommandLine(lhFile,nFiles,rwgtFiles,rwgtstring,dataIndex)
c$$$      filename=trim(lhFILE)//'LHEF_analysis.top'
c$$$      call pwhgtopout(filename)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end



