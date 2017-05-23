      double precision function r(p,i,j)
c----calculate the jets separation between p(i) and p(j)
      implicit none
      include 'mxpart.f'
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
      
