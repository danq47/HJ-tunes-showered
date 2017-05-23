      subroutine setreal(p,rflav,amp2)
c Wrapper subroutine to call the MadGraph real emission matrix
c elements and set the event-by-event couplings constant
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer equivto(maxprocborn)
      common/cequivtoborn/equivto
      include 'pwhg_flg.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'resc_array.h'
      real * 8 p(0:3,nlegreal)
      integer rflav(nlegreal),iborn
      real * 8 amp2
      call set_ebe_couplings
      call sreal_proc(p,rflav,amp2)
c Cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
      amp2 = amp2/(st_alpha/(2d0*pi))     

c rescale to finite top mass; find uborn
c not to be done while looking for identical amplitude
c not to be done while computing regular remnants, that
c do not have an underlying Born (flst_cur_alr=-1 in this case)
      if(.not.flg_in_smartsig.and.flst_cur_alr.gt.0
     1     .and.quarkmasseffects.and..not.rescalebornonly) then
         iborn = flst_alr2born(flst_cur_alr)
         do while(equivto(iborn).ne.-1)
            iborn = equivto(iborn)
         enddo
         amp2 = amp2 * resc_array_tm(iborn)
      endif

      end



      subroutine realcolour_lh
c Wrapper subroutine to call the MadGraph code to associate
c a (leading) color structure to an event.
      implicit none
      include 'nlegborn.h'
      include 'LesHouches.h'
      integer rflav(nlegreal),color(2,nlegreal)
      integer i,j
      do i=1,nlegreal
         rflav(i)=idup(i)
         if (rflav(i).eq.21) rflav(i)=0
      enddo
      call real_color(rflav,color)
      do i=1,2
         do j=1,nlegreal
            icolup(i,j)=color(i,j)
         enddo
      enddo
      end

