      SUBROUTINE BLTAO_QUNTITES(WORK, MAXCOR, IUHF, NONHF)
C
C This routine creates the appropriate AO basis density and
C I(mu,nu) intermediates for analytical gradient calculations.
C This is done in vdint for non-MOLCAS runs. Ajith Perera 07/2000
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      LOGICAL NONHF
C
      DIMENSION WORK(MAXCOR)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
C
      DATA IONE /1/
C
C Get all the basis set information. This is needed here because of the drop
C core gradients. The distinction between spherical or Cartesian is not
C really needed for MOLCAS runs. The Number of basis functions is the number
C of basis functions!
C
      CALL SETBAS
C
      I000 = IONE
      I010 = I000 + NBASIS*NBASIS
      I020 = I010 + NBASIS*NBASIS
      I030 = I020 + NBASIS*NBASIS
      I040 = I030 + NBASTT*2
      I050 = I040 + NBASTT
      I060 = I050 + NBASIS*NBASIS
C
      IF (I060 .GE. MAXCOR) CALL INSMEM(bltao_quntities, I060, MAXCOR)
C
      CALL TRANS_MO2AOS(WORK(I000), WORK(I010), WORK(I020), WORK(I030),
     &                  WORK(I040), WORK(I050), IUHF, NONHF)
C
      RETURN
      END
