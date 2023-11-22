      SUBROUTINE MKFOCK_BLCK(FOCK, TMP1, TMP2, NOCC, NVRT, NBAS, 
     &                       ISPIN)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION FOCK(NBAS,NBAS), TMP1(NOCC,NOCC), TMP2(NVRT,NVRT)
      INTEGER POP, VRT, DIRPRD
      CHARACTER*2 SP(2)
C
#include "machsp.com"
C
      COMMON /SYM/ POP(8,2),VRT(8,2),NFMI(2),NFEA(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      DATA SP /"A", "B"/
C      
      IOFF = 0
      JOFF = 0
      ITAR = 0
      JTAR = 0
      IBGN = 0
      CALL ZERO(TMP1, NOCC*NOCC)
      CALL ZERO(TMP2, NVRT*NVRT)
      DO IRREP = 1, NIRREP
         DO JPOP = 1, POP(IRREP, ISPIN)
               JTAR = JTAR + 1
            DO IPOP = 1, POP(IRREP, ISPIN)
               ITAR = ITAR + 1
               TMP1(ITAR, JTAR) = FOCK(IPOP+IOFF, JPOP+JOFF)
            ENDDO
            ITAR = IBGN
         ENDDO
         IOFF = POP(IRREP, ISPIN) + IOFF
         JOFF = POP(IRREP, ISPIN) + JOFF
         ITAR = JTAR
         IBGN = JTAR
      ENDDO
#ifdef _DEBUG_LVL0
         Write(6,*) "The symmetry packed occ-occ block"
         call output(tmp1, 1, nocc, 1, nocc, nocc, nocc, 1)
#endif
      CALL PUTREC(20, "JOBARC", "FOCKIJ"//SP(ISPIN), NOCC*NOCC*IINTFP,
     &            TMP1)
C
      IOFF = NOCCO(ISPIN)
      JOFF = NOCCO(ISPIN) 
      ITAR = 0
      JTAR = 0
      IBGN = 0
      DO IRREP = 1, NIRREP
         DO JPOP = 1, VRT(IRREP, ISPIN)
            JTAR = JTAR + 1
            DO IPOP = 1, VRT(IRREP, ISPIN)
               ITAR = ITAR + 1
               TMP2(ITAR, JTAR) = FOCK(IPOP+IOFF, JPOP+JOFF)
            ENDDO
            ITAR = IBGN
         ENDDO
         IOFF = VRT(IRREP, ISPIN) + IOFF
         JOFF = VRT(IRREP, ISPIN) + JOFF
         ITAR = JTAR
         IBGN = JTAR
      ENDDO
#ifdef _DEBUG_LVL0
         Write(6,*) "The symmetry packed vrt-vrt block"
         call output(tmp2, 1, nvrt, 1, nvrt, nvrt, nvrt, 1)
#endif
      CALL PUTREC(20, "JOBARC", "FOCKAB"//SP(ISPIN), NVRT*NVRT*IINTFP,
     &            TMP2)
     
      RETURN
      END
 
