       SUBROUTINE READ_S2OLD(ICORE,MAXCOR,IUHF,ISPIN)

       IMPLICIT INTEGER (A-Z)
       DIMENSION ICORE(MAXCOR), NUMSZS(8)

#include "machsp.com"
#include "syminf.com"
#include  "sym.com"
C
      COMMON /SINFO/ NS(8),SIRREP
      COMMON/SLISTS/LS1IN,LS1OUT,LS2IN(2,2),LS2OUT(2,2)
C
C S(m=MIXSPIN,aj=ISPIN)
C
      DO 5 XIRREP = 1, NIRREP
         MIRREP = DIRPRD(XIRREP, SIRREP)
         NUMSZS(XIRREP) = VRT(MIRREP,ISPIN) * NS(SIRREP)
 5    CONTINUE     
C
      IOFF = 1
      DO MIXSPIN = 1, 1+IUHF
         MSPIN = MIXSPIN
         LISTS2IN = LS2IN(ISPIN, MIXSPIN + 1 - IUHF)
         CALL GETLEN_HHP(LENS, POP(1,MIXSPIN), POP(1,ISPIN),
     $                   VRT(1,MIXSPIN), NS)
C
C AAA and BBB blocks
C
         IF (IOFF .GT. MAXCOR) CALL INSMEM("@-READS2",IOFF,MAXCOR)
C
         IF ((IUHF.NE.0) .AND. (MIXSPIN.EQ.ISPIN)) THEN
            CALL GETEXP2_HHP(ICORE(IOFF), LENS, NUMSZS,
     +                       IRPDPD(1,20+ISPIN),
     +                       LISTS2IN,VRT(1,ISPIN),
     +                       IRPDPD(1,ISYTYP(1,LISTS2IN)))
         ELSE
C
C ABB blocks
C
            IF (IOFF .GT. MAXCOR) CALL INSMEM("@-READS2",IOFF,MAXCOR)
C
            CALL GETALLS2_HHP(ICORE(IOFF),LENS,VRT(1,MIXSPIN),
     +                        NS(1),1,LISTS2IN)
         ENDIF
         IOFF = IOFF + LENS * IINTFP 

      ENDDO
C
      RETURN
      END

