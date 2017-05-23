c...lhefheader(nlf)
c...reads initialization information from a les houches events file on unit nlf. 
      subroutine lhefreadwritehdr(nlf,nlo)
      implicit none
      integer nlf,nlo
      character * 200 string
      integer ipr
      include 'LesHouches.h'
 1    read(nlf,fmt='(a)',err=998,end=998) string
      write(nlo,fmt='(a)') trim(string)

      if(string(1:5).eq.'<init') then
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
C         call lhefreadwriteextra(nlf,nlo)
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
      subroutine lhefreadwriteextrannlo(resc,nfiles,rwgtstring,nlf,nlo)
      implicit none
      include 'pwhg_weights.h'
      include 'pwhg_bookhist-multi.h'
C     number of HNNLO scales 
      include 'maxrwgtfiles.h'
      real * 8 resc(maxRwgtFiles,maxmulti)
      character*(20) rwgtstring(maxRwgtFiles)
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      character * 200 string
      integer nlf,nlo,ifile,nfiles
      integer iweights_num
      real *8 tmp 
cccc
      logical addnewweightline
      real * 8 renscfact,facscfact
      integer lhans1,lhans2
      character*(3) PDFp
      common/creweightinfo/renscfact,facscfact,lhans1,lhans2,
     $     addnewweightline,PDFp
cccc
      iweights_num = 0
 1    continue
      string=' '
      read(unit=nlf,fmt='(a)',end=998) string
      string=adjustl(string)
      if(string.eq.'<event>') then
         backspace nlf
         return
      else
         if(string.eq.'</event>') then
            if(addnewweightline) then
               iweights_num=1
               do ifile=1,nfiles
                  tmp = xwgtup * resc(ifile,iweights_num)
                  write(nlo,*) '#new weight,renfact,facfact,pdf1,pdf2',
     1                 tmp,renscfact,facscfact,lhans1,lhans2,PDFp,
     2                 "  ",trim(rwgtstring(ifile))
               enddo
            endif
            write(nlo,'(a)') trim(string)
            return
         elseif(string(1:11).ne."#new weight") then
            write(nlo,'(a)') trim(string)
         elseif(string(1:11).eq.'#new weight') then
            if(iweights_num.eq.weights_max) then
               write(*,*) ' too many weights!'
               write(*,*) ' increase weights_max'
               call exit(-1)
            endif
            iweights_num = iweights_num + 1
            read(string(38:),*) weights_val(iweights_num),
     1           weights_renfac(iweights_num),
     2           weights_facfac(iweights_num),
     3           weights_npdf1(iweights_num),
     4           weights_npdf2(iweights_num),
     5           weights_whichpdf(iweights_num)
            
C     perform rescaling 
            do ifile=1,nfiles 
               tmp = 
     $              weights_val(iweights_num)*resc(ifile,iweights_num)
               
               write(nlo,*) string(1:37),
     $              tmp,
     1              weights_renfac(iweights_num),
     2              weights_facfac(iweights_num),
     3              weights_npdf1(iweights_num),
     4              weights_npdf2(iweights_num),
     5              weights_whichpdf(iweights_num),
     6              "  ",trim(rwgtstring(ifile))
               
            enddo
         endif

      endif
      goto 1
      return
 998  continue
      if (iweights_num /= weights_num) then 
         write(*,*) 'lhefreadwriteextra: weights_num do not match', 
     $        iweights_num,weights_num 
         stop 
      endif
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
