
      subroutine cmpute_grad_analtcaly

#include "aces2.h"

      integer icore, is, ius
      logical do_derint
      character*79 szExport

#ifdef _TRACE_ACES2
      write(*,*) '@ACES2: Entered ', __FILE__
#endif

      ndrgeo=0
      szGExport = 'xa2proc grad_extrp'
      call a2getrec(20,'JOBARC','NDROPGEO',1,NDRGEO)

      call scf

      if (plain_scf) then
         call runit(der_integral_package)
         call c_gtod(is,ius)
         print '(a,f10.1,a)', 'ACES2: Total elapsed time is ',
     &                        is+1.d-6*ius-dTimeStart,' seconds'
         
      else
         if (ndrgeo.ne.0)  then
            call runit('cp JOBARC JOBARC_AM')
            call runit('cp JAINDX JAINDX_AM')
         end if
C
         call runit('xvtran')
         if (iflags2(152).eq.1) call runit('cp HF2 HF2N')
         call runit('xintprc')
         if (iflags2(152).eq.1) call runit('mv HF2N HF2')
         if (iflags(87).eq.0) then
            call runvcc
            call runit('xlambda')
         else if (iflags(87).eq.1) then
            call runit('xvee')
         else
            call runvcc
            call runit('xlambda')
            call runit('xvee')
         end if
         if (ndrgeo.ne.0) call runchg
         call runit('xdens')
         call runit('xanti')
         call runit('xbcktrn')
         call runit(der_integral_package)
c
c If the gradients need to modified externally (ie. extrapolation) 
c do it here (note that no extrapolation for HF-SCF calculations),
c 01/2006, Ajith Perera. 
c
         If (bGExport .or. bCExport) then
            call runit(szGExport) 
         endif
c
         call c_gtod(is,ius)
         print '(a,f10.1,a)', 'ACES2: Total elapsed time is ',
     &                        is+1.d-6*ius-dTimeStart,' seconds'
         
      return
      end

