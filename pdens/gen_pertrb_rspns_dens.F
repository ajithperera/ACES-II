#include "flags.h"
      SUBROUTINE GEN_PERTRB_RSPNS_DENS(WORK, MAXCOR, IRREPX, IUHF)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      DIMENSION WORK(MAXCOR),LENVV(2),LENOO(2),LENVO(2)
      LOGICAL MBPT2, CC
C
#include "jodaflags.com"
#include "machsp.com"
C
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /REFTYPE/ MBPT2,CC
C
      DATA IONE, ONE /1, 1.0D0/    

C
      MBPT2 = .FALSE. 
      IF (iflags2(h_IFLAGS_calc).eq.2)  MBPT2 = .TRUE.
      IF (.NOT. MBPT2) CC = .TRUE.
C
      DO ISPIN = 1, 1 + IUHF
         LENVV(ISPIN) = IRPDPD(IRREPX, 18+ISPIN)
         LENOO(ISPIN) = IRPDPD(IRREPX, 20+ISPIN)
         LENVO(ISPIN) = IRPDPD(IRREPX,  8+ISPIN)
      ENDDO
C
      IDOO = IONE
      IDVV = IDOO + LENOO(1) + IUHF*LENOO(2)
      IDVO = IDVV + LENVV(1) + IUHF*LENVV(2)
      IDOV = IDVO + LENVO(1) + IUHF*LENVO(2)
      ITOP = IDOV + LENVO(1) + IUHF*LENVO(2)
C
      MXCOR = MAXCOR - ITOP + 1
C
         

CSSS      Print*, "Enter Resort"
      CALL RESORT(WORK,MAXCOR,IUHF,IRREPX,444,434)  
CSSS      Print*, "Out of resort"
C
          LSTR1 = 490
          LSTL1 = 190
       LSTR1OFF =   0 
       LSTL1OFF =   0
         LISTL2 = 144 
         LISTR2 = 444
        LSTL2RS = 134 
        LSTR2RS = 434
         LSTGRL = 400
         LSTGTL =   0
       LSTGRLOF =   0
       LSTGTLOF =   2
         LSTTMP = 490
       LSTTMPOF =   2
          LSTT1 =  90
       LSTT1OFF =   0
             R0 =   0.0D0
C
#ifdef _DEBUG_LVLM2
      print*, "Entering PTDENS"
#endif
      CALL PTDENS(1,IRREPX,IRREPX,WORK(IDOO),WORK(IDVV),WORK(IDVO),
     &            WORK(IDOV),WORK(ITOP),MXCOR,IUHF,1.0D0,R0,ONE,
     &            0.0D0,LSTGRL,LSTGTL,LSTGRLOF,LSTGTLOF,LSTTMP,
     &            LSTTMPOF,LSTR1,LSTL1,LSTR1OFF,LSTL1OFF,LISTR2,
     &            LISTL2,LSTR2RS,LSTL2RS,LSTT1, LSTT1OFF)
C
      CALL SAVE_PDENS(WORK(IDOO), WORK(IDOO + NFMI(1)), WORK(IDVV),
     &                WORK(IDVV + NFEA(1)), WORK(IDVO), 
     &                WORK(IDVO + NT(1)), WORK(IDOV), 
     &                WORK(IDOV + NT(1)), IUHF)
C 
      RETURN
      END



