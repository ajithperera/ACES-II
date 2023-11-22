      subroutine Post_scf_anlytcl_grad(Do_derint)

#include "aces2.h"

      integer icore, is, ius
      logical Do_derint
CSSS Dropmo_analy_grads
      character*79 szGExtrap

#ifdef _TRACE_ACES2
      write(*,*) '@ACES2: Entered ', __FILE__
#endif

      ndrgeo=0
      szGExtrap = 'xa2proc grad_extrp'
      call a2getrec(20,'JOBARC','NDROPGEO',1,NDRGEO)
c
c   o use NUMDROPA since ndrgeo is 0 for single-point calcs.
      Call A2getrec(-1, 'JOBARC', 'NUMDROPA', 1, ndropa)
#ifdef _DEBUG_LVLM
      Print*, "Ndropa in ACES II driver", Ndropa,
     &         Analytical_gradient
#endif
CSSS      dropmo_analy_grads = (Analytical_gradient .AND.
CSSS     &                      Ndropa .ne. 0)
c
c We have not yet run lambda for non excitation energy calcualtions
c
      If ((Iflags(87).eq.0).and.(.not. BlCCSDT)) then
         call runit('xlambda')
      Endif
c
c This rather tedious process for DROPMO gradients should be handled
c differently. The DROPMO procedure need to be reimplemented.
c
      if (ndrgeo.ne.0) call runchg
c
      call runit('xdens')

C The HF-DFT like runs with correlated densities; 06/2018.

      If (Hf_Dft) Then
         Call Runit('xintgrt')
         Return 
      End if
c
c Return for First-order properties.
c
      If (.NOT. Do_derint) Return
      call runit('xanti')
      call runit('xbcktrn')
      call runit(der_integral_package)
c
c If the gradients need to modified externally (ie. extrapolation)
c do it here (note that no extrapolation for HF-SCF calculations),
c 01/2006, Ajith Perera.
c
      If (bGExtrap .or. bCExtrap) then
         call runit(szGExtrap)
      endif
c
      call c_gtod(is,ius)
      print '(a,f10.1,a)', 'ACES2: Total elapsed time is ',
     &                      is+1.d-6*ius-dTimeStart,' seconds'

      return
      end

