#include "flags.h"
      SUBROUTINE SCRNC_FORM_RCC_HBAR_INV(W,MAXCOR,IUHF,ILOC)
C
C THIS ROUTINE DRIVES THE drCCD Hbar diagonalizations 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSIZ
      LOGICAL CIS
      DOUBLE PRECISION W(MAXCOR)
      DIMENSION ILOC(8)

      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
#include "flags.com"

      DATA TWO, TOL /2.0D0,1.0D-08/
#ifdef _DEBUG_LVL0
      Write(6,"(1x,a)") "-----Entering scrnc_form_rcc_hbar_inv------"
#endif
C
      CIS = .FALSE.

      IF (CIS) THEN
         LISTDIAG = 42
      ELSE
         LISTDIAG = 56
      ENDIF 

C LOOP OVER SYMMETRY BLOCKS, READ THE MATRIX IN AND DIAGONALIZE IT
C
      DO IRREP=1,NIRREP

C RHF block

         IF (IUHF.EQ.0)THEN
            NUMDIS=IRPDPD(IRREP,ISYTYP(2,56))
            MATDIM = NUMDIS
CSSS            I000 = 1
            I000=ILOC(IRREP)
            I010=I000+NUMDIS*NUMDIS
            CALL GETLST(W(I000),1,NUMDIS,1,IRREP,LISTDIAG)
         ELSE

C UHF AA and BB blocks 

            NUMAA=IRPDPD(IRREP,ISYTYP(1,19))
            NUMBB=IRPDPD(IRREP,ISYTYP(1,20))
            MATDIM=NUMAA+NUMBB

CSSS            I000 = 1
            I000=ILOC(IRREP)
            I010=I000+MATDIM*MATDIM

            IOFFAAAA=ILOC(IRREP)
CSSS            IOFFAAAA=1
            IOFFBBBB=IOFFAAAA+MATDIM*NUMAA+NUMAA
            IOFFBBAA=IOFFAAAA+NUMAA
            IOFFAABB=IOFFAAAA+MATDIM*NUMAA

            IF (CIS) THEN
               LISTW=23
            ELSE
               LISTW=54
            ENDIF 
         
            NUMDIS=IRPDPD(IRREP,ISYTYP(1,LISTW))        
            IOFF=IOFFAAAA
            DO 11 IDIS=1,NUMDIS
               CALL GETLST(W(IOFF),IDIS,1,1,IRREP,LISTW)
               IOFF=IOFF+MATDIM
11          CONTINUE

            IF (CIS) THEN
               LISTW=24
            ELSE
               LISTW=55
            ENDIF 

            NUMDIS=IRPDPD(IRREP,ISYTYP(1,LISTW))        
            IOFF=IOFFBBBB
            DO 12 IDIS=1,NUMDIS
               CALL GETLST(W(IOFF),IDIS,1,1,IRREP,LISTW)
               IOFF=IOFF+MATDIM
12          CONTINUE

            IF (CIS) THEN
               LISTW1=17
               LISTW2=17
            ELSE 
               LISTW1=56
               LISTW2=57
            ENDIF 

            DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))        
            NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))        
            IOFF1=IOFFBBAA
            IOFF2=IOFFAABB

            DO 13 IDIS=1,NUMDIS
               CALL GETLST(W(IOFF1),IDIS,1,1,IRREP,LISTW1)
               CALL GETLST(W(IOFF2),IDIS,1,1,IRREP,LISTW2)
               IOFF1=IOFF1+MATDIM
               IOFF2=IOFF2+MATDIM 
13          CONTINUE

            NUMDIS=MATDIM
         ENDIF
C       
#ifdef _DEBUG_LVL0
        Write(6,"(a)") " The drCCD Hbar matrix"
        Iprint=Iloc(Irrep)
        call output(W(Iprint),1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
#endif 
C Do a full non-symmetric diagonalization

CSSS         I000 = 1
         I000 = ILOC(IRREP)
         I010 = I000 + MATDIM * MATDIM
         I020 = I010 + MATDIM 
         IEND = I020
         IF (IEND .GT. MAXCOR) CALL INSMEM("SCRNC_FORM_RCC_HBAR_INV",
     +                                     IEND,MAXCOR)

         CALL MINV(W(I000),MATDIM,MATDIM,W(I010),DET,TOL,0)
#ifdef _DEBUG_LVL0
        Iprint=Iloc(Irrep)
        Write(6,"(a)") " The drCCD Hbar inverse matrix"
        call output(W(Iprint),1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
#endif 

      ENDDO 
      RETURN 
      END
