      subroutine setborn(p,bflav,born,bornjk,bmunu)
c Wrapper subroutine to call the MadGraph Borns
c and set the event-by-event couplings constant
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer equivto(maxprocborn)
      common/cequivtoborn/equivto
      include 'pwhg_flg.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
      integer iborn
      real * 8 p(0:3,nlegborn),bornjk(nlegborn,nlegborn)
      integer bflav(nlegborn)
      real * 8 bmunu(0:3,0:3,nlegborn),born,resc,born_tm,born_tm_nob
      integer saveafer
      include 'resc_array.h'
      call set_ebe_couplings
      call sborn_proc(p,bflav,born,bornjk,bmunu)
      if(.not.flg_in_smartsig.and.quarkmasseffects) then
         call setborn_tm(p,bflav,born_tm)
         born_tm = born_tm * st_alpha/(2*pi)
         saveafer = afer
         if(rescalenlotoponly) then
c afer is the number of fermions in the loop; exclude bottom:
            afer = 1
            call setborn_tm(p,bflav,born_tm_nob)
            afer = saveafer

            born_tm_nob = born_tm_nob * st_alpha/(2*pi)
c rescale factor with top only in this case
            resc = born_tm_nob/born
            iborn = flst_cur_iborn
            do while(equivto(iborn).ne.-1)
               iborn = equivto(iborn)
            enddo
c store also the full rescaling factor. The Born term, at
c the end, must be rescaled to this
            resc_array_tm_all(iborn) = born_tm/born
         else
C     full rescale factor
            resc   = born_tm/born
         endif
c rescalebornonly is on, the born is not rescaled here (it is
c also used for computing subtraction terms, and in the
c denominator of the sudakov exponent. It is rescaled in the
c btildeborn routine, using the pwhg_born_rescaling.h hook,
c that calls the rescaleborn routine
         if(.not.rescalebornonly) then
            born   = born * resc
            bornjk = bornjk * resc
            bmunu  = bmunu  * resc
         endif
         iborn = flst_cur_iborn
         do while(equivto(iborn).ne.-1)
            iborn = equivto(iborn)
         enddo
         resc_array_tm(iborn) = resc
      endif
      end

      subroutine rescaleborn(res)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer equivto(maxprocborn)
      common/cequivtoborn/equivto
      include 'resc_array.h'
      real * 8 res(*)
      integer j,iborn
      if(quarkmasseffects) then
         if(rescalebornonly) then
            do j=1,flst_nborn
               iborn = j
               do while(equivto(iborn).ne.-1)
                  iborn = equivto(iborn)
               enddo
               res(j) = res(j) * resc_array_tm(iborn)
            enddo
         elseif(rescalenlotoponly) then
            do j=1,flst_nborn
               iborn = j
               do while(equivto(iborn).ne.-1)
                  iborn = equivto(iborn)
               enddo
               res(j) = res(j) * resc_array_tm_all(iborn)/resc_array_tm(iborn)
            enddo
         endif
      endif
      end

      subroutine borncolour_lh
c Wrapper subroutine to call the MadGraph code to associate
c a (leading) color structure to an event.
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      integer equivto(maxprocborn)
      common/cequivtoborn/equivto
      include 'LesHouches.h'
      integer bflav0(nlegborn),bflav(nlegborn),color(2,nlegborn)
      integer i,j,iborn
      logical samecol,conjcol
c We should reach the madgraph flavour configuration that
c was actually computed, in case smartsig is on
      iborn = rad_ubornidx
      bflav0 = flst_born(:,iborn)
      do while(equivto(iborn).ne.-1)
         iborn=equivto(iborn)
      enddo
      bflav = flst_born(:,iborn)
      call born_color(bflav,color)
c Now we have the colour configuration associated with the
c amplitude that was computed instead of rad_ubornidx.
c However, that amplitude may differ from the original one
c by charge conjugation. Check if this is the case
      call matchcolour(nlegborn,bflav0,color)
      icolup(:,1:nlegborn)=color(:,1:nlegborn)
      end


      subroutine finalize_lh
      implicit none
      include 'pwhg_flg.h'
c Specify here if resonances need be written in the event file.
      if(flg_ckkwscalup) call change_scalup
      end



      subroutine change_scalup
      implicit none
      include 'LesHouches.h'
c      include 'hepevt.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      real * 8 ptmin,ptmin2,pcm(0:3,maxnup),beta,vec(3),ptkj2
      integer k,mu,j
      logical ini
      save ini
      data ini/.true./
      real * 8 dotp
      external dotp
      integer npart
      npart=0
      do k=3,nup
c     only light partons
         if (idup(k).eq.21.or.abs(idup(k)).le.5) then
            npart=npart+1
            do mu=1,3
               pcm(mu,npart)=pup(mu,k)
            enddo
            pcm(0,npart)=pup(4,k)
         endif
      enddo
c     compute min pt of light partons with respect to the incoming beam
      ptmin2=1d30
      do k=1,npart
         ptmin2=min(ptmin2,pcm(1,k)**2+pcm(2,k)**2)
      enddo

c     compute pt's of the final state partons with respect to each other
      beta=-(pup(3,1)+pup(3,2))/(pup(4,1)+pup(4,2))
      vec(1)=0
      vec(2)=0
      vec(3)=1
c     go in the CM frame   
      call mboost(npart,vec,beta,pcm,pcm)
      do k=1,npart-1
         do j=k+1,npart
            ptkj2 = 2*dotp(pcm(0,k),pcm(0,j))*
     $           pcm(0,k)*pcm(0,j)/(pcm(0,k)+pcm(0,j))**2
            ptmin2=min(ptmin2,ptkj2)
         enddo
      enddo
      ptmin=sqrt(ptmin2)
      if(scalup.gt.ptmin) then
         if (ini) then
            write(*,*) '*************************************'
            write(*,*) 'scalup set to the min pt in the event'
            write(*,*) '*************************************'
            ini=.false.
         endif         
         scalup = ptmin
      endif
      end
