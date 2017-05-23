      subroutine writeinfo(unitno,xsec,xsec_err)
************************************************************************
*   Routine to write out run information to a desired unit             *
************************************************************************
      implicit none
      include 'maxwt.f'
      include 'masses.f'
      include 'facscale.f'
      include 'scale.f'
      include 'zerowidth.f'
      include 'flags.f'
      include 'clustering.f'
      include 'gridinfo.f'
      include 'limits.f'
      include 'pdfiset.f'
      integer unitno
      double precision xsec,xsec_err
      
      character*4 part
      character*30 runstring
      logical creatent,dswhisto,makecuts
      integer nproc,ih1,ih2,itmx1,itmx2,ncall1,ncall2,rseed
      integer order
      double precision sqrts
      double precision Rcut
      double precision leptpt,leptrap,misspt,Rjlmin,Rllmin,delyjjmin,
     . leptpt2,leptrap2,gammpt,gammrap,gammcone,gammcut
      integer lbjscheme
      logical jetsopphem
 
      common/outputflags/creatent,dswhisto      

      common/nnlo/order
      common/part/part
      common/runstring/runstring
      common/energy/sqrts
      common/density/ih1,ih2
      common/iterat/itmx1,ncall1,itmx2,ncall2
         
    
      
      common/Rcut/Rcut
      common/makecuts/makecuts
      common/leptcuts/leptpt,leptrap,misspt,Rjlmin,Rllmin,delyjjmin,
     . leptpt2,leptrap2,gammpt,gammrap,gammcone,gammcut,
     . lbjscheme,jetsopphem

      common/rseed/rseed

C - Commented out so that gnuplot scripts work
C - when comparing Powheg-Box and HNNLO output
C - and to have a single file format for both.
c$$$      write(unitno,'(a,e14.8,a,e14.8)')
c$$$     $     '! Cross-section = ',xsec,'+/-',xsec_err
c$$$      write(unitno,'(a)')
c$$$     $     '! Run corresponds to this input file:'
c$$$      write(unitno,'(a,e14.8)')
c$$$     $     '! sqrts = ',sqrts
c$$$      write(unitno,'(a,i3,a,i3)')
c$$$     $     '! ih1   = ',ih1,'   ih2   = ',ih2
c$$$      write(unitno,'(a,e14.8)')
c$$$     $     '! Mh        = ',hmass
c$$$      write(unitno,'(a,e14.8)')
c$$$     $     '! muf       = ',facscale
c$$$      write(unitno,'(a,e14.8)')
c$$$     $     '! mur       = ',scale
c$$$      write(unitno,'(a,i14)')
c$$$     $     '! order     = ',order
c$$$      write(unitno,'(a,a)')
c$$$     $     '! part      = ',part
c$$$      write(unitno,'(a,i14)')
c$$$     $     '! itmx1     = ',itmx1
c$$$      write(unitno,'(a,i14)')
c$$$     $     '! ncall1    = ',ncall1
c$$$      write(unitno,'(a,i14)')
c$$$     $     '! itmx2     = ',itmx2
c$$$      write(unitno,'(a,i14)')
c$$$     $     '! ncall2    = ',ncall2
c$$$      write(unitno,'(a,i14)')
c$$$     $     '! rnd seed  = ',rseed
c$$$      write(unitno,'(a,i14)')
c$$$     $     '! iset      = ',iset
c$$$      write(unitno,'(a,a)')
c$$$     $     '! runstring = ',runstring
c$$$      write(unitno,'(a)')
c$$$     $     '!'

      return
      
      end
      
      
