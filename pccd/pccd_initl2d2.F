      SUBROUTINE PCCD_INITL2D2(IUHF)
      IMPLICIT INTEGER(A-Z)
      COMMON/FLAGS/IFLAGS(100)
      COMMON /SYM/ POP1(8),POP2(8),VRT1(8),VRT2(8),NT(2),
     &             NF1(2),NF2(2)
#include "syminf.com" /* for nirrep */
C
C CREATE LISTS FOR DOUBLES IN GRADIENT CALCULATIONS
C
      IMODE=0
      CALL INIPCK(1,13,14,146,IMODE,0,1)

      DO 20 ISPIN=1,IUHF+1
        CALL UPDMOI(1,NT(ISPIN),ISPIN,190,0,0)
        CALL UPDMOI(1,NT(ISPIN),ISPIN+2,190,0,0)
        call aces_list_memset(ispin,  190,0)
        call aces_list_memset(ispin+2,190,0)
20     CONTINUE

      CALL INIPCK(1,13,14,3,IMODE,0,1)
      do iGrp = 1, nirrep
         call aces_list_memset(iGrp,3,0)
      end do
      CALL INIPCK(1,13,14,116,IMODE,0,1)
      do iGrp = 1, nirrep
         call aces_list_memset(iGrp,116,0)
      enddo

      RETURN
      END

