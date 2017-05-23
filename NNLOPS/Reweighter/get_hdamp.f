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

