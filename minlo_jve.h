c -*- Fortran -*-
c contains minlo options to use ptj rather than ptH 

      logical minlo_jve
      double precision minlo_deltaR 

      common/cminlo_jve/minlo_deltaR,minlo_jve
      save /cminlo_jve/ 
