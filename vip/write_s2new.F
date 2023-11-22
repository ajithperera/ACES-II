       SUBROUTINE WRITE_S2NEW(ICORE,MAXCOR,IUHF,ISPIN)

#include "machsp.com"
#include "syminf.com"
#include  "sym.com"
      COMMON /SINFO/ NS(8), SIRREP
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
C
C S(m=MIXSPIN,aj=ISPIN)
C
      DO MIXSPIN = 1, IUHF+1 
         LISTS2EX = LS2OUT(ISPIN, MIXSPIN + 1 - IUHF)

         IF ((IUHF.NE.0) .AND. (MIXSPIN.EQ.ISPIN)) THEN
            CALL ASSYMALL(ICORE(I020),LENS,NUMSZS,IRPDPD(1,18+ISPIN),
     &                    VRT(1,ISPIN),ICORE(I030),MAXCOR-I030+1)
            CALL GETEXP2_HHP(ICORE(I010),LENS,NUMSZS,IRPDPD(1,20+ISPIN),
     &                       LISTS2EX,VRT(1,ISPIN),
     &                       IRPDPD(1,ISYTYP(1,LISTS2EX)))
       ELSE
          CALL GETALLS2_HHP(ICORE(I010),LENS,POP(1,SOUTSPIN),NS(1),
     &                      1,LISTS2EX)
       ENDIF
C
       CALL SAXPY (LENS,ONE,ICORE(I020),1,ICORE(I010),1)
C
       IF ((IUHF.NE.0).AND. (SOUTSPIN.EQ.ISPIN)) THEN
          CALL PUTSQZ_HHP(ICORE(I010),LENS,NUMSZS,IRPDPD(1,20+ISPIN),
     &                    LISTS2EX, VRT(1, ISPIN),
     &                    IRPDPD(1, ISYTYP(1, LISTS2EX)),ICORE(I020),
     &                     MAXCOR-I020+1)
       ELSE
          CALL PUTALLS2_HHP(ICORE(I010),LENS,VRT(1,ISPIN),NS(1),1,
     $                      LISTS2EX)
       ENDIF

      ENDDO
C
      RETURN
      END

