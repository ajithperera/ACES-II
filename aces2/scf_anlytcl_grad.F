      Subroutine Scf_anlytcl_grad(Do_derint)

#include "aces2.h"

      integer icore, is, ius
      logical Do_derint

#ifdef _TRACE_ACES2
      write(*,*) '@ACES2: Entered ', __FILE__
#endif

CSSS      IF (Do_derint) call Runit(der_integral_package)
C
      If (Hf_Dft .and. Do_derint) Then
C
         Call Runit('xvtran')
         Call Runit('xintprc')
         Call Runit(der_integral_package)
         Call Runit('xvksdint')
C       
      Else If (Do_derint .and. ks_scf) Then 
C
         Call Runit(der_integral_package)
         Call Runit('xvksdint')
C
      Else If (Do_derint) Then
C            
         Call Runit(der_integral_package)
C
      End If
     
      call c_gtod(is,ius)
      print '(a,f10.1,a)', 'ACES2: Total elapsed time is ',
     &                      is+1.d-6*ius-dTimeStart,' seconds'
         
         
      Return
      End


