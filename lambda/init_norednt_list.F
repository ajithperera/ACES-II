      SUBROUTINE INIT_NOREDNT_LIST(IRREPX, SYTYPL, SYTYPR, LIST)
      IMPLICIT NONE

      integer irrepx, sytypl, sytypr, list

c#include "lists.com"
#include "syminf.com"
#include "sympop.com"

CSSS      DO IRREPR = 1, NIRREP
CSSS         IRREPL = DIRPRD(IRREPR, IRREPX)
CSSS         NCOLS  = IRPDPD(IRREPR, SYTYPR)
CSSS         NROWS  = IRPDPD(IRREPR, SYTYPL)
CSSS         ILEFT  = IRREPR
CSSS         IRIGHT = LIST
CSSS         MOIOSZ(ILEFT, IRIGHT) = NROWS
CSSS         MOIODS(ILEFT, IRIGHT) = NCOLS
CSSS      END DO

      ISYTYP(1, LIST) = SYTYPL
      ISYTYP(2, LIST) = SYTYPR

      RETURN
      END
