      subroutine find_jets(Rmin,p,mjets) ! p are output jet momenta 
      implicit  none
      include  'hepevt.h'
C - Maximum number of particles to send to jet finder: 
      include 'mxpart.f'
      include 'npart.f'
      include 'jetcuts.f'
      include 'jetlabel.f'
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
      include  'mxpart.f'
      include  'jetlabel.f'
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

