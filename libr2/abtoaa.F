      SUBROUTINE ABTOAA(ICORE, MAXCOR, IUHF, LISTT2)
C 
C Generate AA list from a AB list.
C 
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR)     
C
#include "machsp.com"
#include "sympop.com"
#include "syminf.com"
#include "sym.com"
C
      IF (IUHF .EQ. 0) THEN
   
         DO IRREP=1, NIRREP
            LISTAB   = LISTT2 + 2 
            LISTAA   = LISTT2
            NUMABCOL = IRPDPD(IRREP,ISYTYP(2,LISTAB))
            NUMABROW = IRPDPD(IRREP,ISYTYP(1,LISTAB))
C
            NUMAACOL = IRPDPD(IRREP,ISYTYP(2,LISTAA))
            NUMAAROW = IRPDPD(IRREP,ISYTYP(1,LISTAA))
C
                I000 = 1
                I010 = I000 + IINTFP*NUMABCOL*NUMABROW
                I020 = I010 + IINTFP*NUMAACOL*NUMAAROW 
                IF (I020 .GE. MAXCOR) CALL INSMEM("@ABTOAA",
     &                                I020, MAXCOR)
C
            CALL GETLST(ICORE(I010), 1, NUMABCOL, 1, IRREP, LISTAB)
            CALL ASSYM(IRREP, POP, NUMABROW, NUMABROW, ICORE(I000), 
     &                 ICORE(I010))
            CALL SQSYM(IRREP, VRT, NUMAAROW, NUMABROW, NUMAACOL, 
     &                 ICORE(I010), ICORE(I000))
            CALL PUTLST(ICORE(I010), 1, NUMAACOL, 1, IRREP, LISTAA)
          ENDDO
      ENDIF
C
      RETURN
      END
