      subroutine bmass_in_minlo(result,alphas)
      implicit none
      real * 8 result,alphas
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      include 'deltasud.h'
      include 'PhysPars.h'
      include 'pwhg_pdf.h'
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1         pdf2(-pdf_nparton:pdf_nparton)
      double precision save_result,oyh,opth 
      double precision deltasud,x1b,x2b,tau0,tau1,eta
      save save_result,oyh,opth
      data opth/-1d0/
      logical ini
      data ini/.true./
      save ini
c Setup gaussian weights and points
      if(ini) then
         if(ph_bottommass.le.0d0) then
            write(*,*) 'bmass_in_minlo:'
            write(*,*) "cannot set option bmass_in_minlo to 1"
            write(*,*) "unless also bottommass is set > 0 in "
            write(*,*) "the input file. Exiting."
            call exit(-1)
         endif
         call dgset(0d0,1d0,ngauss,xgauss,wgauss)
         ini = .false.
      endif
      yh = 0.5d0*log((kn_pborn(0,3)+kn_pborn(3,3))/
     1               (kn_pborn(0,3)-kn_pborn(3,3)))
      pth = sqrt(kn_pborn(1,3)**2+kn_pborn(2,3)**2)
      if(yh.eq.oyh.and.pth.eq.opth) then
         result = save_result
         return
      endif
      mb = ph_bottommass
      q2 = kn_sbeams
      q  = sqrt(q2)
      pth2 = pth**2
      mh2 = mh**2
      mth2 = mh**2+pth**2
      mth = sqrt(mth2)
      eta=exp(-2*abs(yh))
      tau0 = (pth+mth)**2/q2
      tau1 = ((eta*q2-mh2)*sqrt(eta*q2*mth2)+eta*q2*pth2)
     1     /(q2*(eta*q2-mth2))
      if(tau0.lt.eta) then
         if(tau1.lt.eta) then
            logtaumin = log(tau0)
            logtaumax = 0
         else
            write(*,*) 'bmass_in_minlo: tau0<eta and tau1>eta !!!'
            write(*,*) "can't be, exiting ..."
            call exit(-1)
         endif
      else
         if(tau1.lt.tau0) then
            write(*,*) 'bmass_in_minlo: tau1<tau0 !!!'
            write(*,*) "can't be, exiting ..."
            call exit(-1)
         endif
         if(tau1.gt.1) then
            write(*,*) 'bmass_in_minlo: tau1>1 !!!'
            write(*,*) "can't be, exiting ..."
            call exit(-1)
         endif
         logtaumin = log(tau1)
         logtaumax = 0
      endif
      alphads = alphas
C      alphads = st_alpha ! fixed coupling 
      call dointdeltasud(result)

c now result is the integral of the non-factorizable term
c divided by the Born term. The Born pdf were not divided out;
c do it now
      x1b = mh/q*exp(yh)
      x2b = mh2/q2/x1b
      call pdfcall(1,x1b,pdf1) 
      call pdfcall(2,x2b,pdf2)
      if(pdf1(0).eq.0.or.pdf2(0).eq.0) then
         result = 0
      else
         result = result / (pdf1(0) * pdf2(0))
      endif
c Now compute the exponential
      result = exp(-result)
c      call plotdeltasud
      save_result = result
      oyh = yh
      opth = pth
      end

      function deltasud0(ycm,tau) result(res)
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      double precision, intent(in) :: ycm,tau
      double precision res
      include 'deltasud.h'
      include 'PhysPars.h'
      double precision kt2,xjac,x1,x2,taub,gg,gq,qg,att,auu,kl,alphas
      include 'pwhg_pdf.h'
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1         pdf2(-pdf_nparton:pdf_nparton)
      integer j

      real * 8 born
      integer  bflav(4)
      real * 8 ppp(0:3,4) ! Partonic COM frame momenta
      real * 8 R_t,R_tb,R_tb_int,R_tb_int_fact,R_tb_int_non_fact
      real * 8 B_t,B_tb,B_t_inf,B_tb_inf
      real * 8 powheginput
      external powheginput
      real * 8 savemfer
      real * 8 approx_result
      integer  saveafer
      logical ini
      data    ini/.true./
      save    ini
      integer approx
      save    approx
      double precision pwhg_alphas
C - Some initialisation ...
      if(ini) then
         approx = powheginput("#bmass_sud_approx")
         if(approx.lt.0d0) approx = powheginput("#bmass_sudakov_apprx")
         if(approx.lt.0d0) approx=2   
         ini=.false.
      endif


C - Compute partonic COM momenta now
      kl = eh * tanh(yh-ycm)
      kt2 = kmom**2-kl**2

      ppp(0,1)=sqrt(shat)/2d0    ! Incoming parton going right
      ppp(1,1)=0d0               ! Incoming parton going right
      ppp(2,1)=0d0               ! Incoming parton going right
      ppp(3,1)=sqrt(shat)/2d0    ! Incoming parton going right
      ppp(0,2)=ppp(0,1)          ! Incoming parton going left
      ppp(1,2)=0d0               ! Incoming parton going left
      ppp(2,2)=0d0               ! Incoming parton going left
      ppp(3,2)=-ppp(3,1)         ! Incoming parton going left
      ppp(0,3)=eh                ! Higgs
      ppp(1,3)=sqrt(kt2)         ! Higgs
      ppp(2,3)=0d0               ! Higgs
      ppp(3,3)=kl                ! Higgs
      ppp(0,4)=kmom              ! Jet
      ppp(1,4)=-ppp(1,3)         ! Jet
      ppp(2,4)=-ppp(2,3)         ! Jet
      ppp(3,4)=-ppp(3,3)         ! Jet

c Compute alpha at kt2
      if(sqrt(kt2).gt.st_lambda5MSB) then
c the -1 argument means do matching at the flavour thresholds
         alphas = pwhg_alphas(kt2,st_lambda5MSB,-1)
         if(alphas.lt.0.or.alphas.gt.alphads) then
            alphas = alphads
         endif
      else
         alphas = alphads
      endif

      if(approx /= 0) then
C - Save current number of fermions and fermion(1)'s (assumed to be top) mass
      saveafer=afer
      savemfer = mfer(1)
C - Compute relevant 1-loop 2->1 matrix element(s), t & b ...
      afer=2
      call m2_gg_h(mh2,B_tb)
C - Compute relevant 1-loop 2->1 matrix element(s), t only ...
      afer=1
      call m2_gg_h(mh2,B_t)
C - Compute relevant 1-loop 2->1 matrix element(s), t & b,  m_top -> infty...
      mfer(1)  = 200000d0       ! i.e. m_top -> infinity
      afer=2
      call m2_gg_h(mh2,B_tb_inf)
C - Compute relevant 1-loop 2->1 matrix element(s), t only, m_top -> infty...
      afer=1
      call m2_gg_h(mh2,B_t_inf)
C - Restore current number of fermions and fermion(1)'s (assumed to be top) mass
      afer=saveafer
      mfer(1) = savemfer
      endif

C - <<<--- End of KH mod.

c
      xjac = z/(4*pi)*(mh2+kt2)/(mh2*(1+z))
c
      x1 = sqrt(tau)*exp(ycm)
      x2 = sqrt(tau)*exp(-ycm)
c
C      call checkkin(x1,x2,kt2)

      call pdfcall(1,x1,pdf1)
      call pdfcall(2,x2,pdf2)
c
      taub = mh**2/(4*mb**2)

      gg = pdf1(0)*pdf2(0)

      gq = 0d0
      do j = 1,5
         gq = gq + pdf1(0)*(pdf2(j)+pdf2(-j))
      enddo

      qg = 0d0
      do j = 1,5
         qg = qg + (pdf1(j)+pdf1(-j))*pdf2(0)
      enddo

      if(approx.eq.3) then ! If neglect finite m_top in R_t do same in B_t(b)
C     ask KH to clarify logic behind approx = 3 
         call gghgApprox(ppp,approx,R_t,R_tb_int)
         gg = gg * ( (R_t+R_tb_int)/B_tb_inf-R_t/B_t_inf )
         gg = gg * alphas/st_alpha
      else if(approx.eq.2) then
         call gghgApprox(ppp,approx,R_t,R_tb_int)
         gg = gg * ( (R_t+R_tb_int)/B_tb-R_t/B_t )
         gg = gg * alphas/st_alpha
      else if(approx.eq.1) then
         bflav=(/0,0,25,0/)
         afer = 2
         call setborn_tm(ppp,bflav,R_tb)
         R_tb = R_tb*st_alpha/(2*pi)
         afer = 1
         call setborn_tm(ppp,bflav,R_t)
         R_t  = R_t *st_alpha/(2*pi)
         gg = gg * ( (R_tb)/B_tb-R_t/B_t )
         afer = saveafer
         gg = gg * alphas/st_alpha
      else if(approx.eq.0) then
         gg = gg
     1     *(8*pi*alphas*6/kt2)*(1-z)/z/taub * 3d0/4/z*(1-z+5d0/4*z**2)
     2     *log(mb**2/kt2)**2
         gg = gg * kt2/(kt2+mb**2)
      endif

      if(approx.eq.3) then ! If neglect finite m_top in R_t do same in B_t(b)
         call qghqApprox(ppp,approx,
     $        R_t,R_tb_int)
         qg = qg * ( (R_t+R_tb_int)/B_tb_inf-R_t/B_t_inf )
         qg = qg * alphas/st_alpha
      else if(approx.eq.2) then
         call qghqApprox(ppp,approx,
     $        R_t,R_tb_int)
         qg = qg * ( (R_t+R_tb_int)/B_tb-R_t/B_t )
         qg = qg * alphas/st_alpha
      else if(approx.eq.1) then
         bflav=(/1,0,25,1/)
         afer = 2
         call setborn_tm(ppp,bflav,R_tb)
         R_tb = R_tb*st_alpha/(2*pi)
         afer = 1
         call setborn_tm(ppp,bflav,R_t)
         R_t  = R_t *st_alpha/(2*pi)
         qg = qg * ( R_tb/B_tb-R_t/B_t )
         afer = saveafer
         qg = qg * alphas/st_alpha
      else if(approx.eq.0) then
c Absolute value of t
         att = sqrt(shat)*(kmom-kl)
         qg = qg *  8d0/2*pi*alphas/taub*(1+(1-z)**2)/z*1/(z*att)*3d0/2
     1        *log(mb**2/att)**2
         qg = qg * kt2/(kt2+mb**2)
      endif
      
      if(approx.eq.3) then ! If neglect finite m_top in R_t do same in B_t(b)
         call gqhqApprox(ppp,approx,
     $        R_t,R_tb_int)
         gq = gq * ( (R_t+R_tb_int)/B_tb_inf-R_t/B_t_inf )
         gq = gq * alphas/st_alpha
      else if(approx.eq.2) then
         call gqhqApprox(ppp,approx,
     $        R_t,R_tb_int)
         gq = gq * ( (R_t+R_tb_int)/B_tb-R_t/B_t )
         gq = gq * alphas/st_alpha
      else if(approx.eq.1) then
         bflav=(/0,1,25,1/)
         afer = 2
         call setborn_tm(ppp,bflav,R_tb)
         R_tb = R_tb*st_alpha/(2*pi)
         afer = 1
         call setborn_tm(ppp,bflav,R_t)
         R_t  = R_t *st_alpha/(2*pi)
         gq = gq * ( R_tb/B_tb-R_t/B_t )
         afer = saveafer
         gq = gq * alphas/st_alpha
      else if(approx.eq.0) then
c Absolute value of u
         auu = sqrt(shat)*(kmom+kl)
         gq = gq *  8d0/2*pi*alphas/taub*(1+(1-z)**2)/z*1/(z*auu)*3d0/2
     1        *log(mb**2/auu)**2
         gq = gq * kt2/(kt2+mb**2)
      endif

      res = (gg + qg + gq) * xjac 
c flux factor
     1 / (2*shat) 
c It has to be divided out by the Born integrated at fixed rapidity, i.e.
c \int  dtau 2pi delta(q2*tau-mh2) * mborn/(2 mh2) = 2pi/(2mh2 * q2) * mborn
c with mborn already divided out in our formula
     2 * (2*mh2*q2)/(2*pi)
c

      end

      function deltasud1(xx) result(res)
c this constructs the kinematics starting from x(1) and (x2),
c and yh, pth, pth2 (in deltasud.h), computes ycm,tau 
c and a jacobian to go from xx(1:2) to ycm,tau
      double precision, intent(in) :: xx(2)
      double precision res
      include 'deltasud.h'
      include 'PhysPars.h'
      double precision xjac,tau,logtau,ycm,zzz
      double precision deltasud0,deltasud0_NZ
c integration limit
      tau = exp(xx(1)**2*(logtaumax-logtaumin)+logtaumin)
      xjac = tau*(logtaumax-logtaumin)*2*xx(1)
c range in ycm
      shat = tau*q**2
      z = mh**2/shat
      eh = mh*(1+z)/(2*sqrt(z))
      kmom = mh*(1-z)/(2*sqrt(z))
      aslim = atanh(sqrt(kmom**2-pth2)/eh)
      ycmmin = max(log(tau)/2,
     1     yh-aslim)
      ycmmax = min(-log(tau)/2,
     1     yh+aslim)
      if(ycmmin.gt.ycmmax) then
         write(*,*) 'deltasud1:  Worry!'
         write(*,*) 'ycmmmin>ycmmax!!'
         call exit(-1)
      endif
      zzz = 6*(xx(2)**2/2-xx(2)**3/3)
      xjac = xjac * 6*xx(2)*(1-xx(2))
      ycm = zzz*(ycmmax-ycmmin)+ycmmin
      xjac = xjac*(ycmmax-ycmmin)
      res = deltasud0(ycm,tau)*xjac
      end


      
      
      subroutine dointdeltasud(res)
      implicit none
      real * 8 res
      include 'deltasud.h'
      real * 8 xx(2),deltasud1,random
      integer j,k,ncalls
      parameter (ncalls=100000)
      res = 0
c      do j=1,ncalls
c         xx(1)=random()
c         xx(2)=random()
c         res = res + deltasud1(xx)/ncalls
c      enddo
c      return

      do j = 1,ngauss
         xx(1) = xgauss(j)
         do k = 1,ngauss
            xx(2) = xgauss(k)
            res = res + deltasud1(xx)*wgauss(j)*wgauss(k)
         enddo
      enddo
      end

      subroutine plotdeltasud
      implicit none
      real * 8 res
      include 'deltasud.h'
      integer ncalls
      parameter (ncalls=100000)
      character * 20 str
      real * 8 xx(2),deltasud1,random
      integer j,k
      logical ini
      data ini/.true./
      save ini
      write(str,'(f4.1,"-",i3)') yh,int(pth)
      do j=1,15
         do k=1,2
            if(str(j:j).eq.' ') str(j:)=str(j+1:)
         enddo
      enddo
      write(*,*) '!'//str//'!'

      call inihists
      call bookupeqbins('x1',0.01d0,0d0,1d0)
      call bookupeqbins('x2',0.01d0,0d0,1d0)

      res = 0
      do k = 1,ncalls
         xx(1) = random()
         xx(2) = random()
         res = deltasud1(xx)
         call filld('x1',xx(1),res)
         call filld('x2',xx(2),res)
         call pwhgaccumup
      enddo
      call pwhgsetout
      call pwhgtopout("deltasudhists-"//trim(str))
      end


      subroutine checkkin(x1,x2,kt2)
      implicit none
      real * 8 x1,x2,kt2,ycm,p,p0,ppar,ycmh,shat_lcl
      include 'deltasud.h'
      include 'PhysPars.h'
      if(x1.lt.0.or.x1.gt.1) then
         write(*,*) ' checkkin: x1=',x1
         call exit(-1)
      endif
      if(x2.lt.0.or.x2.gt.1) then
         write(*,*) ' checkkin: x2=',x2
         call exit(-1)
      endif
      if(kt2.lt.0) then
         write(*,*) ' checkkin: kt2=',kt2
         call exit(-1)
      endif
      ycm = 0.5d0*log(x1/x2)
      shat_lcl = q2*x1*x2
      ycmh = yh-ycm
      p0 = (shat_lcl + mh2)/(2*sqrt(shat_lcl))
      p = (shat_lcl - mh2)/(2*sqrt(shat_lcl))
      ppar = p0*tanh(ycmh)
      if(abs(p**2-kt2-ppar**2).gt.1d-6) then
         write(*,*) ' checkkin: kinematics do not check'
         call exit(-1)
      endif
      end
      



C ******************************************************** C
C ******************************************************** C
C ******************************************************** C
C -                                                      - C
C - Below here is all that is needed for the approximate - C
C - H+1 jet matrix elements                              - C
C -                                                      - C
C ******************************************************** C
C ******************************************************** C
C ******************************************************** C

      subroutine gghgApprox(ppp,approx,
     $           R_t,R_tb_int)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      real * 8 ppp(0:3,nlegborn) ! These are only used if the exact finite
      integer  approx
      real * 8 R_t,R_t_factor
      real * 8 s,t,u,pt,z
      real * 8 mb_2
      real * 8 taub,F02
      real * 8 R_tb_int,R_tb_int_fact,R_tb_int_non_fact
      real * 8 B_t,B_tb,B_t_inf,B_tb_inf
      integer  saveafer
      real * 8 savemfer
      real * 8 dotp
      external dotp

C - Mandelstams
      s = 2 * dotp(ppp(:,1),ppp(:,2)) ! BMZ conventions
      t =-2 * dotp(ppp(:,2),ppp(:,4)) ! BMZ conventions
      u =-2 * dotp(ppp(:,1),ppp(:,4)) ! BMZ conventions

C - b quark mass stuff
      mb_2 = ph_bottommass**2
      taub = mh2/4/mb_2

C - Kinematics
      pt = sqrt(u*t/s)
      z  = mh2/s

C - Temporarily store mfer(1) and afer and set mfer(1) (mt) -> infty
      savemfer = mfer(1)
      saveafer=afer
      mfer(1)  = 200000d0
C - Get Born matrix element with mtop = infty top & bottom contributions
      call m2_gg_h(mh2,B_tb_inf)
C - Get Born matrix element with mtop = infty contribution alone
      afer=1
      call m2_gg_h(mh2,B_t_inf)      
C - Restore mfer(1) (mt) and afer
      mfer(1)  = savemfer
      afer=saveafer

C - Temporarily store mfer(1) and afer and set mfer(1) (mt) -> infty
      saveafer=afer
C - Get Born matrix element with mtop = infty top & bottom contributions
      call m2_gg_h(mh2,B_tb)
C - Get Born matrix element with mtop = infty contribution alone
      afer=1
      call m2_gg_h(mh2,B_t)      
C - Restore mfer(1) (mt) and afer
      afer=saveafer

C - |F_0|^2
      F02  = B_tb/B_t  ! Different to notes |F_0|^2 by b-loop^2 contribs in
                       ! here (numerically 1.43%).

C - gg->Hg ME in m_top->infty limit and no b-loop
      R_t = B_t_inf * (8*pi*st_alpha*6/mh2**2)
     $              * (mh2**4+s**4+t**4+u**4)/(2*s*t*u)

C - Get finite top mass correction factor for 2->2 ME?
      if(approx.eq.2) then
C - Rescaling factor |M^t-only,finite|^2/|M^t-only,inf|^2
C - (next to leading 1/mt^2 part only) ...
C - WARNING! If you use this and finite m_top in computing B_t
C -          above, you will over count 1/mt^2 corrections in
C -          the collinear limit with the code as it is now!

         R_t_factor = B_t/B_t_inf
      elseif(approx.eq.3) then
C - Don't do finite top-mass corrections (OK for gg, not for gq qg though) ...
         R_t_factor = 1d0
      else
         write(*,*) 'gghggApprox called with wrong approx,', approx 
         call pwhg_exit(-1) 
      endif
      
C - |M^(t,b interference)_gg->Hg|^2 factorisable bit ...
      R_tb_int_fact = R_t * (F02-1) 
C - |M^(t,b interference)_gg->Hg|^2 add factorisation breaking bit ...
      R_tb_int_non_fact = R_t * 1/(1-(-pt**2)/mh2)
     $                  * ( (-pt**2)/mh2 * (F02-1)
     $                    + (3d0/4)/taub
     $                    * ( 1/z*(1-z+5d0/4*z**2)
     $                      / (z*(1-z)+z/(1-z)+(1-z)/z)
     $                      * log(mb_2/(pt**2+mb_2))**2
     $                      - 4*pt**2/(pt**2+mh2)*log(pt**2/mh2)
     $                      )
     $                    )



C - Do finite top mass corrections on real ME
      R_t = R_t * R_t_factor
      R_tb_int_fact = R_tb_int_fact * R_t_factor
      R_tb_int_non_fact = R_tb_int_non_fact * R_t_factor

C - Add factorisable and non-factorisable to get total interference 
      R_tb_int = R_tb_int_fact + R_tb_int_non_fact

      return
      end


      subroutine qghqApprox(ppp,approx,
     $           R_t,R_tb_int)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      real * 8 ppp(0:3,nlegborn)
      integer  bflav(nlegborn)
      integer  approx
      real * 8 mh_2
      real * 8 s,t,u
      real * 8 R_t,R_tb_int,R_tb_int_fact,R_tb_int_non_fact
      real * 8 dotp
      external dotp
C - Mandelstams
      s = 2 * dotp(ppp(:,1),ppp(:,2)) ! BMZ conventions
      t =-2 * dotp(ppp(:,2),ppp(:,4)) ! BMZ conventions
      u =-2 * dotp(ppp(:,1),ppp(:,4)) ! BMZ conventions
      call qqhgApproxCore(ppp,approx,u,t,s,
     $     R_t,R_tb_int)
      R_t      = 3d0/8*R_t
      R_tb_int = 3d0/8*R_tb_int
      end


      subroutine gqhqApprox(ppp,approx,
     $           R_t,R_tb_int)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      real * 8 ppp(0:3,nlegborn)
      integer  bflav(nlegborn)
      integer  approx
      real * 8 mh_2
      real * 8 s,t,u
      real * 8 R_t,R_tb_int,R_tb_int_fact,R_tb_int_non_fact
      real * 8 dotp
      external dotp
C - Mandelstams
      s = 2 * dotp(ppp(:,1),ppp(:,2)) ! BMZ conventions
      t =-2 * dotp(ppp(:,2),ppp(:,4)) ! BMZ conventions
      u =-2 * dotp(ppp(:,1),ppp(:,4)) ! BMZ conventions
      call qqhgApproxCore(ppp,approx,t,s,u,
     $     R_t,R_tb_int)
      R_t      = 3d0/8*R_t
      R_tb_int = 3d0/8*R_tb_int
      end

      subroutine qqhgApprox(ppp,approx,
     $           R_t,R_tb_int)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      real * 8 ppp(0:3,nlegborn)
      integer  bflav(nlegborn)
      integer  approx
      real * 8 mh_2
      real * 8 s,t,u
      real * 8 R_t,R_tb_int,R_tb_int_fact,R_tb_int_non_fact
      real * 8 dotp
      external dotp
C - Mandelstams
      s = 2 * dotp(ppp(:,1),ppp(:,2)) ! BMZ conventions
      t =-2 * dotp(ppp(:,2),ppp(:,4)) ! BMZ conventions
      u =-2 * dotp(ppp(:,1),ppp(:,4)) ! BMZ conventions
      call qqhgApproxCore(ppp,approx,s,t,u,
     $     R_t,R_tb_int)
      R_t      = 3d0/8*R_t
      R_tb_int = 3d0/8*R_tb_int
      end

      subroutine qqhgApproxCore(ppp,approx,s,t,u,R_t,R_tb_int)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      real * 8 ppp(0:3,nlegborn)
      integer  approx
      real * 8 mb_2,mh_2,mt_2
      real * 8 s,t,u
      real * 8 taub,F02
      real * 8 R_t,R_t_factor
      real * 8 R_tb_int,R_tb_int_fact,R_tb_int_non_fact
      integer  saveafer
      real * 8 B_t,B_tb,B_t_inf,B_tb_inf
      real * 8 auxiliary
      real * 8 savemfer

C - Higgs mass stuff
      mh_2 = s+t+u

C - top quark mass stuff
      mt_2 = ph_topmass**2

C - b quark mass stuff
      mb_2 = ph_bottommass**2
      taub = mh_2/4/mb_2

C - Temporarily store mfer(1) and afer and set mfer(1) (mt) -> infty
      savemfer = mfer(1)
      saveafer=afer
      mfer(1)  = 200000d0
C - Get Born matrix element with mtop = infty top & bottom contributions
      call m2_gg_h(mh_2,B_tb_inf)
C - Get Born matrix element with mtop = infty contribution alone
      afer=1
      call m2_gg_h(mh_2,B_t_inf)      
C - Restore mfer(1) (mt) and afer
      mfer(1)  = savemfer
      afer=saveafer

C - Temporarily store mfer(1) and afer and set mfer(1) (mt) -> infty
      saveafer=afer
C - Get Born matrix element with mtop = infty top & bottom contributions
      call m2_gg_h(mh_2,B_tb)
C - Get Born matrix element with mtop = infty contribution alone
      afer=1
      call m2_gg_h(mh_2,B_t)      
C - Restore mfer(1) (mt) and afer
      afer=saveafer

C - |F_0|^2
      F02  = B_tb/B_t  ! Different to notes |F_0|^2 by b-loop^2 contribs in
                       ! here (numerically 1.43%).

C - qq->Hg ME in m_top->infty limit and no b-loop
      R_t = B_t_inf * (8*pi*st_alpha*4d0/3/mh_2**2)
     $              * (8d0/3)
     $              * (t**2+u**2)/sqrt(s*s)

C - Get finite top mass correction factor for 2->2 ME?

      if(approx.eq.2) then
C - Same as last bit of 'if' block but only 6 lines of code ...
         if(s.lt.0d0) then
            R_t_factor =
     $      s - mh_2 + 2*s*(sqrt(4*mt_2/mh_2-1)*asin(sqrt(mh_2/mt_2)/2)
     $                     -sqrt(1d0-4*mt_2/s) *asinh(sqrt( -s/mt_2)/2))
     $               + (4*mt_2-u-t)
     $               * (asin (sqrt(mh_2/mt_2)/2)**2
     $                 +asinh(sqrt(  -s/mt_2)/2)**2)
            R_t_factor = R_t_factor**2 * 36*mt_2**2/(mh_2-s)**4
         else
            R_t_factor = 1d0
         endif
      elseif(approx.eq.3) then
         R_t_factor = 1d0
      else
         write(*,*) 'qqhgApproxCore called with wrong approx,', approx 
         call pwhg_exit(-1) 
      endif

C - |M^(t,b interference)_qq->Hg|^2 factorisable bit ...
      R_tb_int_fact = R_t * (F02-1) 
C - |M^(t,b interference)_qq->Hg|^2 add factorisation breaking bit ...
      if(s.ge.0d0) then
         auxiliary = log( mb_2/s)**2 - 4*s/(s-mh_2)*log( s/mh_2) - pi**2
      else
         auxiliary = log(mb_2/(-s+mb_2))**2 - 4*s/(s-mh_2)*log(-s/mh_2)
      endif
      auxiliary = auxiliary * 3d0/4/taub
      R_tb_int_non_fact = R_t * 1/(1-s/mh_2)
     $                  * ( (s/mh_2) * (F02-1) + auxiliary )

C - Do finite top mass corrections on real ME
      R_t = R_t * R_t_factor
      R_tb_int_fact = R_tb_int_fact * R_t_factor
      R_tb_int_non_fact = R_tb_int_non_fact * R_t_factor

C - Add factorisable and non-factorisable to get total interference 
      R_tb_int = R_tb_int_fact + R_tb_int_non_fact

      return
      end


      subroutine M2_gg_h(s,amp2)
c     Born matrix element times normalizations and averages.
c     IMPORTANT: the flux factor 1/2s is intentionally missing 
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
      include 'Flags.h'
      real * 8 s,v2,tiny,tmp,xnorm
      parameter (v2=0.70710678118654757d0)
      parameter (tiny=1.d-8)
      real * 8 amp2
      integer nlegs
      parameter (nlegs=nlegborn)
      real * 8 pphy(0:3,nlegs)
      integer i
      real * 8 m12,y12,bornqcd
      complex * 16 ampl, reduced,x12,aux
      common /bornampl/ampl,bornqcd
      external reduced
      real * 8 getdeltaew,ewcorr
      external getdeltaew
      common /ew/ewcorr

c      s=(pphy(0,1)+pphy(0,2))**2-(pphy(1,1)+pphy(1,2))**2
c     $     -(pphy(2,1)+pphy(2,2))**2 -(pphy(3,1)+pphy(3,2))**2
      xnorm=1/(2.d0*pi)*1/(pi)*s/(256*v2)

c     1/(2pi) comes from the 2*pi*delta(s-m^2) of phase space
      ampl = dcmplx(0d0)

      do i=1,afer
         m12=mfer(i)
         y12=m12**2/mh2
         x12 = reduced(1d0/y12)
         aux = lambdafer(i)*trfer(i)*
     &        (-4d0)*y12*(2d0-(1d0-4d0*y12)*0.5d0*log(x12)**2)
         ampl = ampl+aux
c         write(*,*) 'fermion ', i , ' ampl ', aux
      end do
c     If enabled, here we add the EW corrections

      tmp = ampl * dconjg(ampl)
      if (flg_ew.eq.1) then
         ewcorr = getdeltaew(ampl)
      else
         ewcorr = 0d0
      endif

      bornqcd = xnorm*st_alpha*st_alpha*ph_GF*2d0*s*tmp
      amp2 = bornqcd*(1d0+ewcorr)
c     the multiplication for 2s is needed to remove the flux factor
       end


