c dummy subroutine to take in p_T and y for the Higgs, and the 4 scale variation factors, and give an output weight.
c Later, we will use something like this to link to SC's interpolating function,
c but for now I'm just trying to see where it fits into the code
      function tune_reweight(y_h,pt_h,renscfact_pwg,facscfact_pwg,
     1         renscfact_mrt,facscfact_mrt)
      implicit none
      real * 8 y_h, pt_h, tune_reweight
      real * 8 renscfact_pwg,facscfact_pwg,renscfact_mrt,facscfact_mrt ! KmuR and KmuF for POWHEG and MRT respectively

# when we get the interpolating function from Stefano we can put it
# here to get the reweighting. For now this is just a placeholder where
# we reweight by 1

      tune_reweight = 1.0
      if(renscfact_mrt.gt.1.9) then
         tune_reweight = tune_reweight + 0.1
      elseif(renscfact_mrt.lt.0.6) then
         tune_reweight = tune_reweight + 0.2
      endif

      if(facscfact_mrt.gt.1.9) then
         tune_reweight = tune_reweight + 0.3
      elseif(facscfact_mrt.lt.0.6) then
         tune_reweight = tune_reweight + 0.4
      endif

      end
