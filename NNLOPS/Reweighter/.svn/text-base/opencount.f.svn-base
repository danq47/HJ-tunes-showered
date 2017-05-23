      subroutine opencount(filename,maxev)
      implicit none
      include 'pwhg_rnd.h'
      integer maxev
      character * 1000 filename
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
