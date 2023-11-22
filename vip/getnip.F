      SUBROUTINE GETNIP(NIP,SIRREP,ISPIN,IUHF)
C
      IMPLICIT INTEGER (A-Z)
      logical lefthand,  singonly, DROPCORE
      DIMENSION NS(8)
C
#include "syminf.com"
#include "sym.com"
#include "sympop.com"

      COMMON/IPCALC/LEFTHAND,SINGONLY,DROPCORE
C
C FIRST FIND THE DIMENSION OF THE MATRICES
C
      CALL IZERO(NS, 8)
      NS(SIRREP) = 1
      IF (IUHF.EQ.0) THEN
         CALL GETLEN_HHP(LENAB,POP(1,1),POP(1,1),VRT(1,1),NS(1))
         LENAA = 0
      ELSE
         MSPIN = 3 - ISPIN
         CALL GETLEN_HHP2(LENAA,IRPDPD(1,ISPIN+2),VRT(1,ISPIN),
     +                     NS(1))
         CALL GETLEN_HHP2(LENAB,IRPDPD(1,14),VRT(1,MSPIN),NS(1))
      ENDIF
      NIP = POP(SIRREP,ISPIN) + LENAA + LENAB
C
      RETURN
      END
