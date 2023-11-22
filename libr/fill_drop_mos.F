
      SUBROUTINE FILL_DROP_MOS(NIRREP, IUHF, NOCC_EXPND, NVRT_EXPND,
     &                         NBAS_FULL, NOCCO_EXPND, NVRTO_EXPND)
C
C Recreate the new population vectors that correspond to the full space
C for dropmo calculations. Ajith Perera, 07/2005.
C
      IMPLICIT INTEGER (A-Z)
      INTEGER POP, VRT
      DIMENSION NDRPOP(8), NDRVRT(8), NOCC_EXPND(8,2),
     &          NVRT_EXPND(8,2), NOCCO_EXPND(2), NVRTO_EXPND(2)

      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF(4)

      CALL GETREC(20, 'JOBARC', 'NDROPPOP', NIRREP, NDRPOP)
      CALL GETREC(20, 'JOBARC', 'NDROPVRT', NIRREP, NDRVRT)

           NBAS_FULL = 0
      NOCCO_EXPND(1) = 0
      NVRTO_EXPND(1) = 0
      NOCCO_EXPND(2) = 0
      NVRTO_EXPND(2) = 0
      DO IRREP = 1, NIRREP
         NOCC_EXPND(IRREP,1) = POP(IRREP,1) + NDRPOP(IRREP)
         NVRT_EXPND(IRREP,1) = VRT(IRREP,1) + NDRVRT(IRREP)
     &
         NBAS_FULL = NBAS_FULL + NOCC_EXPND(IRREP,1)
     &                         + NVRT_EXPND(IRREP,1)
         NOCCO_EXPND(1) = NOCCO_EXPND(1) + NOCC_EXPND(IRREP,1)
         NVRTO_EXPND(1) = NVRTO_EXPND(1) + NVRT_EXPND(IRREP,1)
         IF (IUHF.GT.0) THEN
            NOCC_EXPND(IRREP,2) = POP(IRREP,2) + NDRPOP(IRREP)
            NVRT_EXPND(IRREP,2) = VRT(IRREP,2) + NDRVRT(IRREP)
         ELSE
            NOCC_EXPND(IRREP,2) = NOCC_EXPND(IRREP,1)
            NVRT_EXPND(IRREP,2) = NVRT_EXPND(IRREP,1)
         END IF
         NOCCO_EXPND(2) = NOCCO_EXPND(2) + NOCC_EXPND(IRREP,2)
         NVRTO_EXPND(2) = NVRTO_EXPND(2) + NVRT_EXPND(IRREP,2)
      END DO

#ifdef _DEBUG_LVLM
      Write(6,*) "@-FILL_DROP_MOS: expanded occ. vectors"
      do k=1,(iuhf+1)
      Write(6,"(8I4)")(NOCC_EXPND(i,k),i=1,8)
      Write(6,"(8I4)")(NVRT_EXPND(i,k),i=1,8)
      enddo
#endif 

C At this point we assume that number of both alpha and beta MOs
C are the same. This is assumed throughout the code, but Ken (Andrew's)
C frozen natural orbitals can break this assumption.
C This is incorrect comment. I am not sure why I thought of this.
C Perhaps what I meant is that total remain the same (unlike FNO),
C but for UHF clearly the occ-virtual distribution is different 
C (see above; Ajith Perera, 11,2015). 
C 
CSSS      NOCCO_EXPND(2) = NOCCO_EXPND(1)
CSSS      NVRTO_EXPND(2) = NVRTO_EXPND(1)

      RETURN
      END

