
c This routine is the common error handler for xcc.

#define _CALL_ACES_FIN

      subroutine xcc_errex
      implicit none

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _CALL_ACES_FIN
      call aces_fin
#endif /* _CALL_ACES_FIN */

      call aces_exit(1)

c ----------------------------------------------------------------------

      end

