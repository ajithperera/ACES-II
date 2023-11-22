      subroutine initmax(maxocc,iuhf)
      implicit none
 
      integer iuhf
 
      double precision maxocc

      maxocc = 2.0D0

      if (iuhf.eq.1) then
         maxocc = 1.0D0
      end if

      return
      end
 
