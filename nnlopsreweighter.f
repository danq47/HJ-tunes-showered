C ************************************************************** C
C ************************************************************** C
C ************************************************************** C
C Usage: e.g. 
C ./minnlo-new
c It requires an input file nnlopsreweighter.input, containing <lhe file> <string> <hnnlo output> [<index (1:13)>]
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
      include 'pwhg_lhrwgt.h'
      include 'LesHouches.h'
      include 'pwhg_bookhist-multi.h'
      include 'nnlopsreweighter.h'
c      include 'maxrwgtfiles.h' 
c      include 'mxpart.h'

      character*6 WHCPRG
      common/cWHCPRG/WHCPRG

      real*8  hmass
      common/hmass/hmass

      integer ios,nlf,nlo,nev,iev,inputfile

      character*(100) lhFile
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

      real * 8 yH,pth,HT,ptj1,ptj2,hdamp,rwgtnnlo
      real * 8 etah,ph(4),mh
      real * 8 damping_vec(20)  !:
      real * 8 p(bigmxpart,4)
      integer  ibin,iline,n_minlo_wgts

      integer i,j,k,njets

      character*(100) tmpLine,tmpLine2

      character*(1000) lhFile2
      integer lprefix

      logical filled(maxids)

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
      call readInputFile(inputfile,lhFile,nFiles,rwgtFiles,rwgtstring,dataIndex)
      write(cdataIndex,'(i2)') dataIndex
      cdataindex=trim(adjustl(cdataIndex))

C     Open the event file and count how many events
      call opencount(lhFile,nev)

      call newunit(nlf)
      open(unit=nlf,file=trim(lhFile),status='old',iostat=ios)
      call lhrwgt_readheader(nlf)
      call lhrwgt_clearheader
      call lhrwgt_headerloader
c store number of weights from MiNLO run
      n_minlo_wgts = lhrwgt_nids
c load header information for new nnlo weights
      call lhrwgt_readheader(inputfile)
C     read and copy header
      call lhrwgt_headerloader
c Now check for consistency of nllo.input
      filled = .false.
      do i=1,nfiles
         do j=1,n_minlo_wgts
            call get_weight_index(n_minlo_wgts+1,
     1           trim(lhrwgt_id_arr(j)),trim(rwgtstring(i)),k)
            if(k.gt.0) filled(k) = .true.
         enddo
      enddo
      if(count(filled).ne. lhrwgt_nids - n_minlo_wgts) then
         write(*,*) ' could not fill the following required nnlo elements: '
         do k=n_minlo_wgts+1,lhrwgt_nids
            if(.not.filled(k)) then
               write(*,'(a)') '<weight id='//trim(lhrwgt_id_arr(k))//'/>'
     1              // trim(lhrwgt_descr_arr(k)) //'</weight>'
            endif
         enddo
         write(*,*) ' exiting ...'
         call exit(-1)
      endif

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


c     Find the value of hmass in the header of the event
c     file, and read its value. This is needed, because we 
c     will use this (fixed) value in the reweighting function,
c     and the same must have been used in HNNLO.
c     If this is not found, abort the program.
c     
C     Read in all lines until you find the first occurrence of hmass
      rewind(nlf)
      ios=0
      hmass=-1d0
      do while(ios.eq.0)
         read(nlf,fmt='(a)',iostat=ios) tmpLine
         tmpline = adjustl(tmpline)
         if(tmpline.eq.'<event>') exit
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
      call lhefreadhdr(nlf)
      do j=1,nev
         if(mod(j,10000) == 0) then
            write(*,*) '  step 1: at event ', j,'/',nev
            call pwhgsetout
            i = index(lhFile,'.lhe')
            call pwhgtopout(trim(adjustl(lhfile(:i-1)))//'-denominators')
         endif
         call lhefreadev(nlf)   
         if(WHCPRG == 'LHE   ') then
            call lhuptohepevt(j)
            if(abs(idwtup).eq.3) then
c We should also rescale the weights
c               xwgtup=xwgtup*xsecup(1)
               write(*,*) ' idwtup=',idwtup,' not handled properly!'
               write(*,*) ' exiting ...'
            endif
            call analysis(xwgtup)
         else
c     here shower and then call analysis...
            write(*,*) 'denominators after showering'//
     $           'not implemented. Program stops.'
            stop
         endif
         call pwhgaccumup
      enddo
      close(nlf) 
c     finish denominator plots
      call pwhgsetout
      i = index(lhFile,'.lhe')
      call pwhgtopout(trim(adjustl(lhfile(:i-1)))//'-denominators')

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
      call newunit(nlf)
      open(unit=nlf,file=lhFile,status='old',iostat=ios)



C     open also output file (will contain rescaled weights)
      call newunit(nlo)
C for general distribution use a simplified name 
      lhFile2=trim(lhfile)//'.nnlo'
      open(unit=nlo,file=trim(lhFile2),status='unknown',iostat=ios)
c load header information for new nnlo weights
      call lhrwgt_readheader(inputfile)
C     read and copy header
      call lhefreadwritehdr(nlf,nlo)
      call lhrwgt_headerloader

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

         hdamp=damping_vec(dataindex)
         do i=1,nfiles          ! nfiles: number of HHNLO input files 
            do j=1,n_minlo_wgts  ! number of weights in the event file
               call get_weight_index(n_minlo_wgts+1,
     1              trim(lhrwgt_id_arr(j)),trim(rwgtstring(i)),k)
               if(k.gt.0) then
                  filled(k) = .true.
                  if(prescription.eq.1) then
                     rwgtnnlo =
     $                    num_y(i,dataindex,ibin)/den_y(j,dataindex,ibin)
     $                    * hdamp 
     $                    + 1d0-hdamp
                     lhrwgt_weights(k) = rwgtnnlo * lhrwgt_weights(j)
                  elseif(prescription.eq.2) then
c$$$  print*, den_y(j,1,ibin)-den_y(j,dataindex,ibin)
                     rwgtnnlo =
     $                    (  num_y(i,1,ibin)-
     $                    (den_y(j,1,ibin)-den_y(j,dataindex,ibin)) 
     $                    )
     $                    /den_y(j,dataindex,ibin)
     $                    * hdamp 
     $                    + 1d0-hdamp
                     lhrwgt_weights(k) = rwgtnnlo * lhrwgt_weights(j)
                  else
                     write(*,*) 'reweighting prescription '//
     $                    'not implemented yet.'
                     write(*,*) 'Quitting ...'
                     stop
                  endif
               endif
c$$$               if ((iev == 1) .and. (i == 1)) print*, 'ev 1 resc',
c$$$     $              j,ibin,den_y(j,dataindex,ibin) 
            enddo
         enddo
C     write event extra info with rescaled weights;
c     nnlo weights start with the n_minlo_wgts+1 entry
         call lhefreadwriteextrannlo(n_minlo_wgts+1,nlf,nlo)
      enddo

C     close all units 
      close(nlf) 
      close(nlo) 

C     If ever needed, here add some sanity-check counters

C     write how the reweighted file is called
      write(*,*) 'LH file with nnlo weights: ',trim(lhFile2)
      contains
      end



      subroutine get_weight_index(istart,id1,id2,ind)
      implicit none
      character * (*) id1,id2
      integer istart,ind
      include 'pwhg_lhrwgt.h'
      integer j
      do j=istart,lhrwgt_nids
         if(index(lhrwgt_descr_arr(j),id1) .gt. 0 .and.
     1      index(lhrwgt_descr_arr(j),id2) .gt. 0 ) then
            ind = j
            return
         endif
      enddo
      ind = -1
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
      include 'nnlopsreweighter.h'
c      include 'maxrwgtfiles.h'

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

      integer ios,iun
      integer nDataSetsNum
      character*(100) tmpLine
      integer nLinesNum

      character * 100 tmp_title
      real * 8 tmp_Xlo(maxbins),tmp_Xhi(maxbins)
      real * 8 tmp_Y(maxbins),tmp_Yerr(maxbins)

      ndatasets = 0 

      call newunit(iun)

C - Read in all data in 'numerator' reweighting files
      do ixx=1,nFiles

         nDataSets(ixx)=-1
C - Count number of data sets in the numerator file
         ios=0
         nDataSetsNum=0
         open(unit=iun,file=rwgtFiles(ixx),status='old')
         do while(ios.eq.0)
            read(iun,fmt='(a)',iostat=ios) tmpLine
            if(index(tmpLine,'#').ne.0) nDataSetsNum=nDataSetsNum+1
         enddo
         close(iun)

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
      include 'nnlopsreweighter.h'
c      include 'maxrwgtfiles.h'


      character*(1000) rwgtFile
      integer  dataIndex
      integer  nLines
      character * 100 title
      real * 8 Xlo(maxbins),Xhi(maxbins)
      real * 8 Y(maxbins),Yerr(maxbins)

      integer ios,iun

      character*(100) tmpLine,tmpLine2
      integer tmpIndex

      call newunit(iun)

      open(unit=iun,file=rwgtFile,status='old')

C - Read in all lines until you get to an error or the end of the file
      ios=0
      nLines=0
      do while(ios.eq.0)
         read(iun,fmt='(a)',iostat=ios,end=111) tmpLine
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

 111  continue
      close(iun)      
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
      include 'nnlopsreweighter.h'
c      include 'maxrwgtfiles.h'

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


C ************************************************************** C
C ************************************************************** C
C ************************************************************** C
C - Does a well defined thing. Reads the command line looking
C - for the LH event file to be reweighted, returned to the
C - program as lhfile, the "data" file names to be used to
C - construct the weights, returned as 1-D array rwgtFiles,
C - the number of data files nFiles, and also the data (histogram)
C - index to be used for the reweighting.
      subroutine readInputFile(inputfile,lhFile,nFiles,rwgtFiles,rwgtstring,
     $     dataIndex)           ! all output
      implicit none
      include 'nnlopsreweighter.h'
c      include 'maxrwgtfiles.h'
      character*(*) lhFile
      integer inputfile,nFiles
      character*(*) rwgtFiles(maxRwgtFiles)
      character*(*) rwgtstring(maxRwgtFiles)
      integer dataindex
      character * 200 string, token
      call newunit(inputfile)
      dataindex = 10
      open(unit=inputfile,file='nnlopsreweighter.input',status='old')
      lhfile=' '
      nFiles=0
 1    continue

      read(inputfile,'(a)',end=777) string
      string = adjustl(string)
      if(string(1:6).eq.'lhfile') then
         if(lhfile.ne.' ') then
            write(*,*) 'readInputFile: <lhfile> has already appeared'
            write(*,*) ' exiting ...'
            call exit(-1)
         endif
         lhfile = adjustl(string(7:))
      endif
      if(string.eq.'nnlofiles') then
         if(nFiles.ne.0) then
            write(*,*) 'nnlofiles appears more than once'
            write(*,*) 'exiting ...'
            call exit(-1)
         endif
         do nFiles=1,maxRwgtFiles
            read(inputfile,'(a)',end=666) string
            if(string.eq.' ') exit
            string = adjustl(string)
            call next_token_from_string(string,token)
            select case(token(1:1))
            case('"',"'")
               continue
            case default
               write(*,*) 'readInputFile: after nnlofiles first line'
               write(*,*) 'must be a string; exiting ...'
               call exit(-1)
            end select
            rwgtstring(nFiles) = token
            rwgtFiles(nFiles) = adjustl(string)
         enddo
         if(nFiles.eq.maxRwgtFiles+1) then
            read(inputfile,'(a)',end=666) string
            if(string.ne.' ') then
               write(*,*) 'readInputFile: after nnlofiles'
               write(*,*) 'too many filenames; increase maxRwgtFiles'
               write(*,*) 'exiting ...'
               call exit(-1)
            endif
         endif
         nFiles = nFiles - 1
      endif
      goto 1
 666  continue
      write(*,*) ' readInputFile: some error while reading'
      write(*,*) ' exiting ...'
      call exit(-1)
 777  continue
      rewind(inputfile)
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
      include 'nnlopsreweighter.h'
c      include 'maxrwgtfiles.h'
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



C     a number of handy functions used by generic analyses 

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
      subroutine opencount(filename,maxev)
      implicit none
      include 'pwhg_rnd.h'
      integer maxev
      character * (*) filename
      integer ios
      character * 7 string
      real * 8 powheginput
      external powheginput
      maxev=0
      open(unit=98,file=trim(filename),status='old',iostat=ios)
      if(ios.ne.0) then
         write(*,*)' file not found:',trim(filename)
         stop 
      endif
      
      write(*,*) ' Opened event file ',trim(filename)
      write(*,*) ' Counting events in ', trim(filename)
      write(*,*) ' This may take some time...'
 1    continue
      read(unit=98,fmt='(a)',end=2) string
      if(string.eq.'</event') then
         maxev=maxev+1
         goto 1
      endif
      goto 1
 2    continue
      write(*,*) ' Found ',maxev,' events in file ',trim(filename)
      if (maxev.eq.0) then
         write(*,*) ' NO EVENTS!! Program exits'
         call exit(3)
      endif
      close(98)
      end
c...lhefheader(nlf)
c...reads initialization information from a les houches events file on unit nlf. 
      subroutine lhefreadwritehdr(nlf,nlo)
      implicit none
      integer nlf,nlo
      character * 200 string
      integer ipr
      include 'LesHouches.h'
 1    read(nlf,fmt='(a)',err=998,end=998) string
      if(string.eq.'</initrwgt>') then
c Write LH infor for new NLL weights
         call lhrwgt_writeheader(nlo)
      endif
      write(nlo,fmt='(a)') trim(string)

      if(string(1:6).eq.'<init>') then
         read(nlf,*) idbmup(1),idbmup(2),ebmup(1),ebmup(2),
     &        pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2),idwtup,nprup
         write(nlo,110) idbmup(1),idbmup(2),ebmup(1),ebmup(2),
     &        pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2),idwtup,nprup
         do ipr=1,nprup
            read(nlf,*) xsecup(ipr),xerrup(ipr),xmaxup(ipr),
     &           lprup(ipr)
            write(nlo,120) xsecup(ipr),xerrup(ipr),xmaxup(ipr),
     &           lprup(ipr)
         enddo
         goto 999
      else
         goto 1
      endif
 998  write(*,*) 'lhef_readwritehdr: could not find <init> data'
      call exit(1)
 110  format(1p,2(1x,i8),2(1x,e12.5),6(1x,i6))
 120  format(1p,3(1x,e12.5),1x,i6)
 999  write(nlo,'(a)') '</init>'
      end


Cc...reads event information from a les houches events file on unit nlf. 
      subroutine lhefreadwriteev(nlf,nlo)
      implicit none
      integer nlf,nlo
      character * 200 string
      include 'LesHouches.h'
      integer i,j
 1    continue
      string=' '
      read(nlf,fmt='(a)',err=777,end=666) string
      if(string.eq.'</LesHouchesEvents>') then
         write(nlo,fmt='(a)') trim(string)
         goto 998
      endif
      if(string(1:6).eq.'<event') then
c on error try next event. The error may be cause by merging
c truncated event files. On EOF return with no event found
         write(nlo,'(a)') trim(string) 
         read(nlf,*,end=998,err=1)nup,idprup,xwgtup,scalup,aqedup,aqcdup
         write(nlo,210)nup,idprup,xwgtup,scalup,aqedup,aqcdup
         do i=1,nup
            read(nlf,*,end=998,err=1) idup(i),istup(i),mothup(1,i),
     &           mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     &           vtimup(i),spinup(i)
         enddo
         do i=1,nup
            write(nlo,220) idup(i),istup(i),mothup(1,i),
     &           mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     &           vtimup(i),spinup(i)
         enddo
         call lhefreadwriteextra(nlf,nlo)
         goto 999
      else
         goto 1
      endif
c no event found:
 777  continue
      print *,"Error in reading in lhefreadev"
      print *,string
      nup=0
      return
 666  continue
      print *,"reached EOF"
      print *,string
      nup=0
      return
 998  continue
      print *,"read </LesHouchesEvents>"
      nup=0   
 210  format(1p,2(1x,i6),4(1x,e12.5))   
 220  format(1p,i8,5(1x,i5),5(1x,e16.9),1x,e12.5,1x,e10.3)

 999  end


c     Read extra info (in prticular the multiple weights infos)
c     and change the NLO weights to the NNLO ones.
      subroutine lhefreadwriteextra(nlf,nlo)
      implicit none
      include 'LesHouches.h'
      character * 200 string
      integer istart,nlf,nlo
      integer k,iret
 1    continue
      string=' '
      read(unit=nlf,fmt='(a)',end=998) string
      string=adjustl(string)
      if(string.eq.'<event>') then
         backspace nlf
         return
      else
         if(string.eq.'<rwgt>') then
            write(nlo,'(a)') trim(string)
            call lhrwgt_readwriteweights(nlf,nlo,iret)
            if(iret.ne.0) then
               write(*,*) ' lhefreadwriteextra: Worry:'
     1              //'iret=',iret
            endif
            return
         endif
      endif
      write(nlo,'(a)') trim(string)
      goto 1
      return
 998  continue
      write(*,*) 'lhefreadwriteextrannlo: last event in file'
      end




      subroutine lhrwgt_readwriteweights(nlf,nlo,iret)
      implicit none
      integer nlf,nlo,iret
c assumes that the '<rwgt>' line has already been red
      include 'pwhg_lhrwgt.h'
      integer nweights,iid,iendid
      real * 8 value
      character * 100 string
      nweights = 0
 1    read(unit=nlf,fmt='(a)',end=998) string
      call lhrwgt_id_value_ind(string,iid,iendid,value)
      if(iid.lt.0) goto 996
      write(nlo,'(a)') trim(string)
      nweights = nweights+1
      if(nweights.gt.lhrwgt_nids) then
         write(*,*) ' lhrwgt_loadweights: more weights in event'
         write(*,*) ' than in the header'
         write(*,*) ' exiting ...'
         call exit(-1)
      endif
      if(string(iid:iendid).eq.lhrwgt_id_arr(nweights)) then
         lhrwgt_weights(nweights) = value
      else
         write(*,*) ' lhrwgt_loadweights: id in event'
         write(*,*) ' does not match id in declaration'
         write(*,*) ' exiting ...'
         call exit(-1)
      endif
      goto 1
 996  continue
c did we find </rwgt>?
      if(adjustl(string).eq.'</rwgt>') then
         backspace nlf
         iret = 0
         return
      endif
 998  continue
c end of file
      iret = -1
      end




c     Read extra info (in prticular the multiple weights infos)
c     and change the NLO weights to the NNLO ones.
      subroutine lhefreadwriteextrannlo(istart,nlf,nlo)
      implicit none
      include 'LesHouches.h'
      include 'pwhg_lhrwgt.h'
      character * 200 string
      integer istart,nlf,nlo
      integer k
 1    continue
      string=' '
      read(unit=nlf,fmt='(a)',end=998) string
      string=adjustl(string)
      if(string.eq.'<event>') then
         backspace nlf
         return
      else
         if(string.eq.'</rwgt>') then
            do k=istart,lhrwgt_nids
               write(nlo,'(a,e16.10,a)') '<wgt id='//
     1      trim(lhrwgt_id_arr(k))// '> ',lhrwgt_weights(k),' </wgt>'
            enddo
         endif
      endif
      write(nlo,'(a)') trim(string)
      goto 1
      return
 998  continue
      write(*,*) 'lhefreadwriteextrannlo: last event in file'
      end

c     Read and copy all event file, and add after each event
c     record a reweight line, starting with '#new weight,...',
c     containing the same weight as xwgtup, and infos
c     read from the header, where a copy of the powheg
c     input file was pasted.
      subroutine lhefcopyaddweightline(nlf,nlo,nev)
      implicit none
      include 'LesHouches.h'
      integer ios,nlf,nev,nlo,i,j,iev
      real * 8 renscfact,facscfact
      integer lhans1,lhans2
      character*(3) PDFp
      character*(200) tmpLine,tmpLine2

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
     $           :len(tmpLine))
            read(tmpLine2,*,iostat=ios) renscfact
         endif
         if(index(tmpLine,'facscfact').ne.0) then
            tmpLine2 = tmpLine(index(trim(tmpLine),' ',back=.true.)
     $           :len(tmpLine))
            read(tmpLine2,*,iostat=ios) facscfact
         endif
         if(index(tmpLine,'lhans1').ne.0) then
            tmpLine2 = tmpLine(index(trim(tmpLine),' ',back=.true.)
     $           :len(tmpLine))
            read(tmpLine2,*,iostat=ios) lhans1
         endif
         if(index(tmpLine,'lhans2').ne.0) then
            tmpLine2 = tmpLine(index(trim(tmpLine),' ',back=.true.)
     $           :len(tmpLine))
            read(tmpLine2,*,iostat=ios) lhans2
         endif
         if(index(tmpLine,'PDF package').ne.0) then
            tmpLine2 = tmpLine(index(trim(tmpLine),' ',back=.true.)
     $           :len(tmpLine))
            read(tmpLine2,*,iostat=ios) PDFp
         endif
      enddo
      rewind(nlf)

c     Value for renscfact,facscfact,lhans1,lhans2,PDFp now
c     have been read
c      print*, renscfact,facscfact,lhans1,lhans2,PDFp

c     Rewind, go through all the file, copy each line,
c     and add the extra lines

C     read and copy header
      call lhefreadwritehdr(nlf,nlo)
      
      iev=0
 123  continue
      iev=iev+1
      if(iev.eq.nev+1) then
         return
      endif
      call lhefreadwriteev(nlf,nlo)
 456  continue
      read(nlf,fmt='(a)',iostat=ios,end=123) tmpLine
      if(index(tmpLine,'<event>').ne.0) then
c     new event, so go back one line and use lhefreadwriteev
         backspace nlf
         goto 123
      elseif(index(tmpLine,'</event>').ne.0) then
c     last line of event, so add the new reweight line
         write(nlo,*) '#new weight,renfact,facfact,pdf1,pdf2',
     1        xwgtup,renscfact,facscfact,lhans1,lhans2,PDFp
         write(nlo,'(a)') '</event>'
         goto 456
      else
c     generic line, copy from old file leaving unchanged
         write(nlo,fmt='(a)') trim(tmpLine)
         goto 456
      endif
      end
C     computes damping function 
      subroutine get_hdamp(yH,pth,HT,ptj1,ptj2,damping_vec)
      implicit  none
      real*8  yh,pth,HT,ptj1,ptj2
      real*8  damping_vec(20)!:
C - Higgs kinematics
      real*8  hmass
      common/hmass/hmass
C - Damping factor related quantities:
      real*8  HT_exp,pth_exp,ptj1_exp
      real*8  HT_min,pth_min,ptj1_min
      real*8  HT_hfact(4)
      real*8  pth_hfact(4)
      real*8  ptj1_hfact(4)

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



C - yh
      damping_vec(1)=1d0

C - yh with pth damping function
      damping_vec(2)=((pth_hfact(1)-pth_min)**pth_exp)
     $       /((pth_hfact(1)-pth_min)**pth_exp
     $        +(pth-pth_min)**pth_exp)

      damping_vec(3)=((pth_hfact(2)-pth_min)**pth_exp)
     $       /((pth_hfact(2)-pth_min)**pth_exp
     $        +(pth-pth_min)**pth_exp)

      damping_vec(4)=((pth_hfact(3)-pth_min)**pth_exp)
     $       /((pth_hfact(3)-pth_min)**pth_exp
     $        +(pth-pth_min)**pth_exp)

      damping_vec(5)=((pth_hfact(4)-pth_min)**pth_exp)
     $       /((pth_hfact(4)-pth_min)**pth_exp
     $        +(pth-pth_min)**pth_exp)

C - yh with HT damping function
      damping_vec(6)=((HT_hfact(1)-HT_min)**HT_exp)
     $       /((HT_hfact(1)-HT_min)**HT_exp
     $        +(HT-HT_min)**HT_exp)

      damping_vec(7)=((HT_hfact(2)-HT_min)**HT_exp)
     $       /((HT_hfact(2)-HT_min)**HT_exp
     $        +(HT-HT_min)**HT_exp)

      damping_vec(8)=((HT_hfact(3)-HT_min)**HT_exp)
     $       /((HT_hfact(3)-HT_min)**HT_exp
     $        +(HT-HT_min)**HT_exp)

      damping_vec(9)=((HT_hfact(4)-HT_min)**HT_exp)
     $       /((HT_hfact(4)-HT_min)**HT_exp
     $        +(HT-HT_min)**HT_exp)

C - yh with ptj1 damping function
      damping_vec(10)=((ptj1_hfact(1)-ptj1_min)**ptj1_exp)
     $       /((ptj1_hfact(1)-ptj1_min)**ptj1_exp
     $        +(ptj1-ptj1_min)**ptj1_exp)

      damping_vec(11)=((ptj1_hfact(2)-ptj1_min)**ptj1_exp)
     $       /((ptj1_hfact(2)-ptj1_min)**ptj1_exp
     $        +(ptj1-ptj1_min)**ptj1_exp)

      damping_vec(12)=((ptj1_hfact(3)-ptj1_min)**ptj1_exp)
     $       /((ptj1_hfact(3)-ptj1_min)**ptj1_exp
     $        +(ptj1-ptj1_min)**ptj1_exp)

      damping_vec(13)=((ptj1_hfact(4)-ptj1_min)**ptj1_exp)
     $       /((ptj1_hfact(4)-ptj1_min)**ptj1_exp
     $        +(ptj1-ptj1_min)**ptj1_exp)

C - pth
      damping_vec(14)=1d0

C - pth with canonical damping factor
      damping_vec(15)=((pth_hfact(3)-pth_min)**pth_exp)
     $       /((pth_hfact(3)-pth_min)**pth_exp
     $        +(pth-pth_min)**pth_exp)

C - HT
      damping_vec(16)=1d0

C - HT with canonical damping factor
      damping_vec(17)=((HT_hfact(2)-HT_min)**HT_exp)
     $       /((HT_hfact(2)-HT_min)**HT_exp
     $        +(HT-HT_min)**HT_exp)

C - ptj1
      damping_vec(18)=1d0

C - ptj1 with canonical damping factor
      damping_vec(19)=((ptj1_hfact(3)-ptj1_min)**ptj1_exp)
     $       /((ptj1_hfact(3)-ptj1_min)**ptj1_exp
     $        +(ptj1-ptj1_min)**ptj1_exp)

C - ptj2
      damping_vec(20)=1d0

      end

      subroutine genclust_kt(q,Rmin,qfinal,isub)
c--- clusters momenta using plabel to determine which 
c--- particles should be clustered. Forms 'jets' jets according to
c--- the standard kT algorithm with cone size Rmin.
c--- Furthermore, the clustered jets are only observed if
c--- pT(jet) > ptjetmin and eta(jet) < etajetmax
c--- 
c--- qfinal is the final vector q1,.... q(4+jets)
c--- where non-jet four vectors are set equal to the incoming q 
      implicit none
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
c      include 'npart.h'
c      include 'jetcuts.h'
c      include 'jetlabel.h'
      double precision q(bigmxpart,4),qjet(bigmxpart,4),
     $                                qfinal(bigmxpart,4)
      double precision pt,Rmin,dijmin,dkmin,aetarap
      integer i,nu,iter,nmin1,nmin2,maxjet,nk,
     . ajet,jetindex(bigmxpart),isub
      character*2 plabel(bigmxpart)
      logical jetmerge,failed
      common/plabel/plabel
      common/jetmerge/jetmerge

      jets=0
      maxjet=0
      jetmerge=.false.
      
      do i=1,bigmxpart
        do nu=1,4
        qfinal(i,nu)=0d0
        enddo
      enddo

c--- pick out jets: note that we search to npart+2-isub, to get the
c--- number of particles right. Note that isub=0 for all calls except
c--- the dipole contributions, where isub=1.   

      do i=3,npart+2-isub
      if ( (plabel(i) .eq. 'pp') .or. (plabel(i) .eq. 'pj')
     . .or.(plabel(i) .eq. 'bq') .or. (plabel(i) .eq. 'ba')
     . .or.(plabel(i) .eq. 'qj') ) then
        maxjet=maxjet+1
        jetindex(maxjet)=i
        jetlabel(maxjet)=plabel(i)
        do nu=1,4
          qjet(maxjet,nu)=q(i,nu)
        enddo
      endif
      enddo

c--- for no partons, just switch q into qfinal
      if (maxjet .eq. 0) then
        do i=1,bigmxpart
          do nu=1,4
            qfinal(i,nu)=q(i,nu)
          enddo
        enddo
        jets=0
        return
      endif

c--- skip clustering if we only have one parton  
      if (maxjet .eq. 1) goto 2

c--- for W+bbj, skip if b and b-bar are too close together
c      if ( ((nproc.eq.292) .or. (nproc.eq.297))
c     .     .and. (isub .eq.0) .and. (R(q,5,6) .lt. Rmin) ) then
c        jets=-1
c        return
c      endif

      iter=0
c--- loops through all the iterations of the algorithm      
    1 iter=iter+1

c      write(*,*) 'iter ',iter
c      write(*,*) 'jets ',jets
c      write(*,*) 'maxjet ',maxjet

c--- step1: find (i,j) pair with lowest measure of all non-jets so far
      call findmind(q,qjet,iter,maxjet,dijmin,nmin1,nmin2)
 
c--- step2: find jet K with lowest Et
      call findminet(q,qjet,iter,maxjet,dkmin,nk)
      dkmin=dkmin*Rmin

c      write(*,*) 'Comparing pair (',nmin1,',',nmin2,') value of'
c      write(*,*) 'dijmin = ',dijmin,' with ',nk,' value of dk = ',dkmin
      
c--- step3: compare the two ...      
      if (dijmin .lt. dkmin) then
c---  ... if we should combine, go ahead
c        write(*,*) 'Clustered ',nmin1,nmin2
        jetmerge=.true.
        call combine(qjet,nmin1,nmin2)
c--- combined object goes into nmin1, now shuffle nmin2 off the end 
        call swapjet(qjet,jetindex,nmin2,maxjet)        
        maxjet=maxjet-1
        iter=iter-1
c        do i=1,maxjet
c          do j=1,4
c            write(*,*) 'qjet(',i,',',nu,') = ',qjet(i,nu)
c          enddo
c        enddo
      else
c---  ... we've finished a jet
        jets=jets+1
c        write(*,*) 'Now swapping ',jets,' and ',nk
        call swapjet(qjet,jetindex,jets,nk)
      endif

c--- in the next iteration we search for jets in pjet from iter+1...maxjet
c--- so if this condition isn't true then there's one jet left at maxjet

      if (iter .lt. maxjet-1) goto 1
      
 2    continue      
      jets=jets+1

c--- restore incoming partons
      do i=1,2
        do nu=1,4
          qfinal(i,nu)=q(i,nu)
        enddo
      enddo
c--- set all other momenta to zero and restore leptons
      do i=3,npart+2
        do nu=1,4
          qfinal(i,nu)=0d0
          if ((plabel(i) .ne. 'pp') .and. (plabel(i) .ne. 'pj')
     .   .and.(plabel(i) .ne. 'bq') .and. (plabel(i) .ne. 'ba')
     .   .and.(plabel(i) .ne. 'qj')) then
            qfinal(i,nu)=q(i,nu)
          endif
        enddo
      enddo
      
c----remove jets that are below the pT threhold or which lie outside
c----the observable rapidity region
     
c      write(*,*) 'AFTER CLUSTERING: Obtained ',jets,' jets'
     
c--- restore jets
      ajet=0
      do i=1,jets

c        write(*,*) 'Jet ',i,'(',jetlabel(i),')',jetindex(i)
c        write(*,*) 'pt: ',pt(i,qjet),' vs min. ',ptjetmin
c        write(*,*) 'aeta: ',aetarap(i,qjet),' vs min. ',etajetmin
c        write(*,*) 'aeta: ',aetarap(i,qjet),' vs max. ',etajetmax
        if ((pt(i,qjet) .gt. ptjetmin) .and.
     .      (aetarap(i,qjet) .gt. etajetmin) .and.
     .      (aetarap(i,qjet) .lt. etajetmax)) then 
        ajet=ajet+1
        do nu=1,4
          qfinal(jetindex(ajet),nu)=qjet(i,nu)
        enddo
        jetlabel(ajet)=jetlabel(i)
        endif
      enddo
      
c--- if no jets are removed by eta and pt cuts, then jets=ajet
      if (ajet .lt. jets) then
        do i=ajet+1,jets
          do nu=1,4
            qfinal(jetindex(i),nu)=0d0
          enddo
        enddo
        jets=ajet
      endif


      return
      end
      subroutine swapjet(pjet,jetindex,i,j)
c--- swaps jets i..j in pjet
      implicit none
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
c      include 'jetlabel.h'
      integer i,j,k,jetindex(bigmxpart)
      double precision pjet(bigmxpart,4),tmp
      character*2 chartmp
 
      do k=1,4
        tmp=pjet(i,k)
        pjet(i,k)=pjet(j,k)
        pjet(j,k)=tmp
      enddo
 
      chartmp=jetlabel(i)
      jetlabel(i)=jetlabel(j)
      jetlabel(j)=chartmp
      
      return
      end
      subroutine findmind(p,pjet,pjetmin,pjetmax,dijmin,nmin1,nmin2)
c--- this finds the minimum dij for pjet indices pjetmin through pjetmax
c--- returns dijmin and indices of minimum in (nmin1,nmin2)
      implicit none
      include 'nnlopsreweighter.h'
      double precision p(bigmxpart,4),pjet(bigmxpart,4),dijmin,dij,d
      integer pjetmin,pjetmax,nmin1,nmin2,i,j
      logical dijerror

      dijmin=1d9
      dijerror=.true.

      do i=pjetmin,pjetmax
        do j=i+1,pjetmax
          if (i .ne. j) then
            d=dij(p,pjet,i,j)
            if (d .lt. dijmin) then
              dijmin=d
              nmin1=i
              nmin2=j
              dijerror=.false.
            endif
          endif
        enddo
      enddo
      
      if (dijerror) then
        write(*,*) 'Error in dij minimum-finding routine'
CC        stop
      endif
      
      return
      end
      
      subroutine findminet(p,pjet,pjetmin,pjetmax,dkmin,nk)
c--- this finds the minimum dkmin for pjet indices pjetmin through pjetmax
c--- returns dijmin and indices of minimum in (nmin1,nmin2)
C--- calculate the beam proto-jet separation see NPB406(1993)187, Eqn. 7
C--- in  practice this is just the minimum ptsq of protojets       
      implicit none
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
      double precision p(bigmxpart,4),pjet(bigmxpart,4),dkmin,dk,pt
      integer pjetmin,pjetmax,nk,i
      logical dkerror
      
      dkmin=1d9
      dkerror=.true.
      
      do i=pjetmin,pjetmax
        dk=pt(i,pjet)
        if (dk .lt. dkmin) then
          dkmin=dk
          nk=i
          dkerror=.false.
        endif
      enddo
      
      if (dkerror) then
        write(*,*) 'Error in dk minimum-finding routine'
CC        stop
      endif
      
      return
      end
      
      double precision function dij(p,pjet,i,j)
C---calculate the proto-jet separation see NPB406(1993)187, Eqn. 7
      implicit none
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
      integer i,j
      double precision p(bigmxpart,4),pjet(bigmxpart,4),pti,ptj,pt,r
c      double precision etarap,yi,yj,phii,phij
      
      pti=pt(i,pjet)
      ptj=pt(j,pjet)

c--- old method - bad because (phii-phij) can be > pi       
c      yi=etarap(i,pjet)
c      yj=etarap(j,pjet)

c      phii=atan2(pjet(i,1),pjet(i,2))
c      phij=atan2(pjet(j,1),pjet(j,2))
      
c      dij=dsqrt((yi-yj)**2+(phii-phij)**2)

c--- new method - r() calculates true value of 0 < (phi-phij) < pi
      dij=r(pjet,i,j)
            
      dij=dij*min(pti,ptj)
      
      return
      end
      
      subroutine combine(pjet,i,j)
      implicit none
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
c      include 'jetlabel.h'
      integer i,j
      double precision pjet(bigmxpart,4)
      
c--Run II prescription
      pjet(i,1)=pjet(i,1)+pjet(j,1)
      pjet(i,2)=pjet(i,2)+pjet(j,2)
      pjet(i,3)=pjet(i,3)+pjet(j,3)
      pjet(i,4)=pjet(i,4)+pjet(j,4)


      if (((jetlabel(i) .eq. 'bq') .and. (jetlabel(j) .eq. 'pp'))
     ..or.((jetlabel(j) .eq. 'bq') .and. (jetlabel(i) .eq. 'pp'))) then
        jetlabel(i)='bq'
        return
      endif
      if (((jetlabel(i) .eq. 'ba') .and. (jetlabel(j) .eq. 'pp'))
     ..or.((jetlabel(j) .eq. 'ba') .and. (jetlabel(i) .eq. 'pp'))) then
        jetlabel(i)='ba'
        return
      endif
      if (((jetlabel(i) .eq. 'bq') .and. (jetlabel(j) .eq. 'ba'))
     ..or.((jetlabel(j) .eq. 'bq') .and. (jetlabel(i) .eq. 'ba'))) then
        jetlabel(i)='bq'
        return
      endif
      if (((jetlabel(i) .eq. 'bq') .and. (jetlabel(j) .eq. 'qj'))
     ..or.((jetlabel(j) .eq. 'bq') .and. (jetlabel(i) .eq. 'qj'))) then
        jetlabel(i)='bq'
        return
      endif
      if (((jetlabel(i) .eq. 'ba') .and. (jetlabel(j) .eq. 'qj'))
     ..or.((jetlabel(j) .eq. 'ba') .and. (jetlabel(i) .eq. 'qj'))) then
        jetlabel(i)='ba'
        return
      endif
      if (((jetlabel(i) .eq. 'qj') .and. (jetlabel(j) .eq. 'pp'))
     ..or.((jetlabel(j) .eq. 'pp') .and. (jetlabel(i) .eq. 'qj'))) then
        jetlabel(i)='qj'
        return
      endif

      return
      end
      
c      subroutine combine_snowmass(p,pjet,i,j)
c      implicit none
c      include 'mxpart.h'
c      include 'jetlabel.h'
c      integer i,j
c      double precision p(bigmxpart,4),pjet(bigmxpart,4),ptjetij,yjet,phijet,
c     . ejet,pt,etarap,pti,ptj,yi,yj,phii,phij
      
C----Snowmass style prescripton 
c      pti=pt(i,pjet)
c      ptj=pt(j,pjet)
      
c      yi=etarap(i,pjet)
c      yj=etarap(j,pjet)
      
c      phii=atan2(pjet(i,1),pjet(i,2))
c      phij=atan2(pjet(j,1),pjet(j,2))
c

c      ptjetij=pti+ptj
c      yjet=(pti*yi+ptj*yj)/ptjetij
c      phijet=(pti*phii+ptj*phij)/ptjetij
c      ejet=exp(yjet)
      
c      pjet(i,1)=ptjetij*dsin(phijet)
c      pjet(i,2)=ptjetij*dcos(phijet)
c      pjet(i,3)=ptjetij*(ejet-1d0/ejet)/2d0
c      pjet(i,4)=ptjetij*(ejet+1d0/ejet)/2d0


c      if (((jetlabel(i) .eq. 'bq') .and. (jetlabel(j) .eq. 'pp'))
c     ..or.((jetlabel(j) .eq. 'bq') .and. (jetlabel(i) .eq. 'pp'))) then
c        jetlabel(i)='bq'
c        return
c      endif
c      if (((jetlabel(i) .eq. 'ba') .and. (jetlabel(j) .eq. 'pp'))
c     ..or.((jetlabel(j) .eq. 'ba') .and. (jetlabel(i) .eq. 'pp'))) then
c        jetlabel(i)='ba'
c        return
c      endif
c      if (((jetlabel(i) .eq. 'bq') .and. (jetlabel(j) .eq. 'ba'))
c     ..or.((jetlabel(j) .eq. 'bq') .and. (jetlabel(i) .eq. 'ba'))) then
c        jetlabel(i)='pp'
c        return
c      endif
c
c      return
c      end
      
c      subroutine shuffle(pjet,nmin,nmax)
c--- shuffles jets nmin..nmax-1 in pjet down by 1 index
c      implicit none
c      include 'mxpart.h'
c      integer i,j,nmin,nmax
c      double precision pjet(bigmxpart,4)
      
c      if (nmin .eq. nmax) return
      
c      do i=nmin,nmax-1
c        do j=1,4
c          pjet(i,j)=pjet(i+1,j)
c        enddo
c      enddo
      
c      return
c      end
      
      subroutine swap(pjet,i,j)
c--- swaps jets i..j in pjet
      implicit none
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
c      include 'jetlabel.h'
      integer i,j,k
      double precision pjet(bigmxpart,4),tmp
      character*2 chartmp
      
      do k=1,4
        tmp=pjet(i,k)
        pjet(i,k)=pjet(j,k)
        pjet(j,k)=tmp
      enddo
      
      chartmp=jetlabel(i)
      jetlabel(i)=jetlabel(j)
      jetlabel(j)=chartmp
      
      return
      end
      
c      double precision function ptjet(j,p,pjet)
c      implicit none
c      include 'mxpart.h'
c      integer j
c      double precision p(bigmxpart,4),pjet(bigmxpart,4)
c--- This is the formula for pt
c      ptjet=dsqrt(pjet(j,1)**2+pjet(j,2)**2)
c--- This is the formula for Et
c      ptjet=dsqrt(pjet(j,1)**2+pjet(j,2)**2)
c     . *pjet(j,4)/dsqrt(pjet(j,1)**2+pjet(j,2)**2+pjet(j,3)**2)
c      return
c      end

      double precision function dotjet(p,i,pjet,j)
C---Dot the ith vector p with the jth vector pjet
      implicit none
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
      integer i,j
      double precision p(bigmxpart,4),pjet(bigmxpart,4)
       
      dotjet=p(i,4)*pjet(j,4)-p(i,1)*pjet(j,1)
     .      -p(i,2)*pjet(j,2)-p(i,3)*pjet(j,3)

      return
      end
      
      double precision function bclustmass(pjet)
      implicit none
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
c      include 'jetlabel.h'
      integer i,nbq,nba
      double precision pjet(bigmxpart,4)
      
c--- note: this function ASSUMES that there is at most one b-quark
c--- and one anti-b-quark, returning zero if there are less than this

      bclustmass=0d0
      nbq=0
      nba=0
      
      do i=1,jets
        if (jetlabel(i) .eq. 'bq') nbq=i+4
        if (jetlabel(i) .eq. 'ba') nba=i+4
      enddo

      if ((nbq .eq. 0) .or. (nba .eq. 0)) return
      
      bclustmass=(pjet(nbq,4)+pjet(nba,4))**2
      do i=1,3
        bclustmass=bclustmass-(pjet(nbq,i)+pjet(nba,i))**2
      enddo

      return
      end
       
      double precision function etarap(j,p)
      implicit none
C---returns the value of the pseudorapidity
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
      integer j
      double precision p(bigmxpart,4)
      etarap=dsqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)
      etarap=(etarap+p(j,3))/(etarap-p(j,3))
      if (etarap .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      etarap=100d0
      else
      etarap=0.5d0*dlog(etarap)
      endif
      return
      end

      double precision function aetarap(j,p)
      implicit none
C---returns the absolute value of the pseudorapidity
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
      integer j
      double precision p(bigmxpart,4)
      aetarap=dsqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)
      aetarap=(aetarap+p(j,3))/(aetarap-p(j,3))
      if (aetarap .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      aetarap=100d0
      else
      aetarap=0.5d0*abs(dlog(aetarap))
      endif
      return
      end
 
      double precision function yrap(j,p)
      implicit none
C---returns the value of the rapidity
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
      integer j
      double precision p(bigmxpart,4)
      yrap=(p(j,4)+p(j,3))/(p(j,4)-p(j,3))
      if (yrap .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrap=100d0
      else
      yrap=0.5d0*dlog(yrap)
      endif
      return
      end

      double precision function ayrap(j,p)
      implicit none
C---returns the absolute value of the rapidity
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
      integer j
      double precision p(bigmxpart,4)
      ayrap=(p(j,4)+p(j,3))/(p(j,4)-p(j,3))
      if (ayrap .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      ayrap=100d0
      else
      ayrap=0.5d0*dabs(dlog(ayrap))
      endif
      return
      end
 
      double precision function pt(j,p)
      implicit none
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
      integer j
      double precision p(bigmxpart,4)
c--- This is the formula for pt
      pt=dsqrt(p(j,1)**2+p(j,2)**2)
c--- This is the formula for Et 
c      pt=dsqrt(p(j,1)**2+p(j,2)**2)
c     . *p(j,4)/dsqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)
      return
      end

      double precision function pttwo(j,k,p)
      implicit none
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
      integer j,k
      double precision p(bigmxpart,4)
      pttwo=dsqrt((p(j,1)+p(k,1))**2+(p(j,2)+p(k,2))**2)
      return
      end

c--- this is the rapidity of pair j,k
      double precision function yraptwo(j,k,p)
      implicit none
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
      integer j,k
      double precision p(bigmxpart,4)
      yraptwo=(p(j,4)+p(k,4)+p(j,3)+p(k,3))
     .       /(p(j,4)+p(k,4)-p(j,3)-p(k,3))
      if (yraptwo .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yraptwo=100d0
      else 
      yraptwo=0.5d0*dlog(yraptwo)
      endif
            
      return
      end

c--- this is the pseudo-rapidity
      double precision function etaraptwo(j,k,p)
      implicit none
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
      integer j,k
      double precision p(bigmxpart,4)
      
      etaraptwo=dsqrt((p(j,1)+p(k,1))**2+(p(j,2)+p(k,2))**2
     .               +(p(j,3)+p(k,3))**2)
      if (abs(etaraptwo)-abs(p(j,3)+p(k,3)) .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      etaraptwo=100d0
      else 
      etaraptwo=(etaraptwo+p(j,3)+p(k,3))
     .         /(etaraptwo-p(j,3)-p(k,3))
      etaraptwo=0.5d0*dlog(etaraptwo)
      endif
      
      return
      end


CC    NEW

c--- this is the rapidity of pair j,k,l,m
      double precision function yrapfour(j,k,l,m,p)
      implicit none
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
      integer j,k,l,m
      double precision p(bigmxpart,4)
      yrapfour=(p(j,4)+p(k,4)+p(l,4)+p(m,4)+p(j,3)+p(k,3)+p(l,3)+p(m,3))
     .        /(p(j,4)+p(k,4)+p(l,4)+p(m,4)-p(j,3)-p(k,3)-p(l,3)-p(m,3))
      if (yrapfour .lt. 1d-13) then
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrapfour=100d0
      else 
      yrapfour=0.5d0*dlog(yrapfour)
      endif
            
      return
      end
      double precision function r(p,i,j)
c----calculate the jets separation between p(i) and p(j)
      implicit none
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
      double precision p(bigmxpart,4),r1,r2,dely,delphi,ei,ej
      integer i,j


c     !: original code
      ei=dsqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)
      ej=dsqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)

c     !: to do the same as fastjet kt does
c$$$      ei=p(i,4)
c$$$      ej=p(j,4)
      
      r1= (ei+p(i,3))*(ej-p(j,3))/
     .   ((ej+p(j,3))*(ei-p(i,3)))
      dely=0.5d0*dlog(r1)

      r2= (p(i,1)*p(j,1)+p(i,2)*p(j,2))
     .   /dsqrt((p(i,1)**2+p(i,2)**2)*(p(j,1)**2+p(j,2)**2))
      if (r2 .lt. -0.999999999D0) r2=-1D0
      if (r2 .gt.  0.999999999D0) r2=1D0
      delphi=dacos(r2)
      
      r=dsqrt(dely**2+delphi**2)
      
      return
      end
      


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
      include 'pwhg_lhrwgt.h'
      include 'pwhg_bookhist-multi.h'
c     allow multiweights 
      real * 8 dsig0,dsig(1:maxmulti)
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
      include 'nnlopsreweighter.h'
c      include  'mxpart.h'
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
         if(lhrwgt_nids.gt.0) then
            call setupmulti(lhrwgt_nids)
         elseif(weights_num.gt.0) then
            call setupmulti(weights_num)
         else
            call setupmulti(1)
         endif
         ini=.false.
      endif

      dsig=0

      if(lhrwgt_nids.gt.0) then
         dsig(1:lhrwgt_nids)=lhrwgt_weights(1:lhrwgt_nids)
      elseif(weights_num.gt.0) then
         dsig(:weights_num)=weights_val(:weights_num)
      else
         dsig(1)=dsig0
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


      subroutine find_jets(Rmin,p,mjets) ! p are output jet momenta 
      implicit  none
      include  'hepevt.h'
C - Maximum number of particles to send to jet finder: 
      include 'nnlopsreweighter.h'
c      include 'mxpart.h'
c      include 'npart.h'
c      include 'jetcuts.h'
c      include 'jetlabel.h'
      real*8  ph(4)
      character*6 WHCPRG
      common/cWHCPRG/WHCPRG
C - Jet radius:
      real*8  Rmin
C - HNNLO jet finding:
      real*8  q(bigmxpart,4),p(bigmxpart,4)
      integer isub
      character*2 plabel(bigmxpart)
      logical     jetmerge
      common/plabel/plabel
      common/jetmerge/jetmerge
      integer mjets
C - Auxiliary counters:
      integer ixx,nu
      integer nBosons
C - First call indicator:
      logical first
      data    first/.true./
      save    first

      if(WHCPRG.ne.'NLO   '.and.WHCPRG.ne.'LHE   ') then
         write(*,*) '******************************************'
         write(*,*) '******************************************'
         write(*,*) 'find_jets called but NOT at NLO/LHE level'
         write(*,*) 'this has not been tested carefully'
         write(*,*) '******************************************'
         write(*,*) '******************************************'
         stop
      endif

C - Do this because HNNLO does it in its plotter.f:
      if (first) then
        jetmerge=.true.
        first=.false.
      endif

C - Initialise inputs & outputs for the HNNLO genclust_kt subroutine
      jets=0
      do ixx=1,bigmxpart
        do nu=1,4
           q(ixx,nu)=0d0
           p(ixx,nu)=0d0
C - We always have isub=0 as in 'dipole' contributions we have npart 
C - set to the no. of final-state particles in the Born/virtual (of HJ 
C - see below) wheras HNNLO has it defined the same way but using the
C - real (again, in the case of the subtraction terms) then it uses
C - isub inside genclust_kt to effectively redefine npart->npart-isub. 
           isub=0
C - Label everything as blanks by default
           plabel(ixx)='  '
C - From include 'jetlabel.f'
           jetlabel(ixx)='  '
        enddo
      enddo

C - HNNLO's genclust_kt subroutine assumes npart has been set in the
C - associated include 'npart.f'. HNNLO code says npart=2 for 2-loop
C - counterterm, 3 for Virtual, 4 for real. Hence we set in all cases
C - accordingly to the relevant number of final state particles be it
C - in real or virtual kinematics or showered
      npart=nhep-2
C - It also assumes that the following jet variables cuts have been
C - set in include 'jetcuts.f' (we set them as in HNNLO Need/setup.f):
      ptjetmin=0d0    !:
      etajetmin=0d0
      etajetmax=1d6   !:
      ptbjetmin=0d0
      etabjetmax=1d6  !:

C - Get the Higgs and set up the input vector of input four-vectors
C - and particle labels needed for the HNNLO genclust_kt subroutine
      nBosons = 0
C - For HERWIG or PYTHIA, the last H boson is kept, which is likely
C - the Higgs before decaying, if it decays, or the final Higgs, if
C - undecayed.
      do ixx=1,nhep
         do nu=1,4
            q(ixx,nu)=phep(nu,ixx)
         enddo
C - If the particle is a Higgs label it as such and store the momentum
C - as ph, then exit to the next iteration of the ixx loop.
         if(idhep(ixx).eq.25) then
            ph=phep(1:4,ixx)
            nBosons=1
            plabel(ixx)='H0'
            cycle
         endif
C - We get to here for all non-Higgs particles
         if(WHCPRG.eq.'NLO   '.or.WHCPRG.eq.'LHE   ') then
            plabel(ixx)='pp'    !  N.B. This means EVERYTHING but Higgs goes
                                !       into the genclust_kt jet finder.
                                !       This is what HNNLO does.
         endif
      enddo

      if(nBosons.eq.0) then
         write(*,*) 'find_jets: Fatal error.'
         write(*,*) 'Higgs boson not found'
         call exit(-1)
      endif

C - Get HNNLO style jets:
      call genclust_kt(q,Rmin,p,isub)
      mjets = jets ! jets is stored in a common
C - N.B. in HNNLO plotter gets called as "call plotter(pjet,val,0/1)"
C - yet plotter is defined "subroutine plotter(p,wt,switch)" (where
C - switch is not used). So when HNNLO's plotter routine is doing
C - analysis with p(..,..) it is actually doing it with pjet(..,..),
C - as one would hope. What is pjet though - where did it come from?
C - The pjet fed to plotter in HNNLO, in realint.f and virtint.f,
C - results from a call getptildejet(..,pjet). subroutine getptildejet(nd,pjet)
C - simply copies array ptildejet(nd,i,j) to pjet(i,j) where the 
C - former lives in Inc/ptilde.f and the latter is the output. Now
C - the ptildejet array which it copies is set in function includedipole
C - to another LOCALLY DEFINED pjet(i,j), where THIS local pjet(i,j)
C - array results from a call to genclust2(ptrans,rcut,pjet,isub) two
C - lines previous to the setting. Hence we call genclust_kt as above.

      end

C     Order the jets where there are jets to order
      subroutine get_ptj1j2(p,ptj1,ptj2) 
      implicit  none
      include  'hepevt.h'
C     Maximum number of particles to send to jet finder: 
      include 'nnlopsreweighter.h'
c      include  'mxpart.h'
c      include  'jetlabel.h'
      character*2 plabel(bigmxpart)
      common/plabel/plabel
      character*6 WHCPRG
      common/cWHCPRG/WHCPRG
      real *8 p(bigmxpart,4)
      real *8 ptj1,ptj2,pt,pjm

      if(WHCPRG.ne.'NLO   '.and.WHCPRG.ne.'LHE   ') then
         write(*,*) '******************************************'
         write(*,*) '******************************************'
         write(*,*) 'get_ptj1j2 called but NOT at NLO/LHE level'
         write(*,*) 'this has not been tested carefully'
         write(*,*) '******************************************'
         write(*,*) '******************************************'
         stop
      endif
      if(jets.eq.2) then
         ptj1=pt(4,p)  !  This is pt(5,p) in HNNLO as there H decays to YY. 
         ptj2=pt(5,p)  !  This is pt(6,p) in HNNLO as there H decays to YY. 
         if(ptj1.lt.ptj2) then
            pjm=ptj2
            ptj2=ptj1
            ptj1=pjm
         endif
      elseif(jets.eq.1) then
         ptj1=pt(4,p)  !  This is pt(5,p) in HNNLO as there H decays to YY. 
         ptj2=0d0
      elseif(jets.eq.0) then
         ptj1=0d0
         ptj2=0d0
      endif
      end

