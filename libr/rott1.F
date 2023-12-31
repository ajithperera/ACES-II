      SUBROUTINE ROTT1(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE ROTATES THE MOLECULAR ORBITALS BY THE T1-LIKE
C  AMPLITUDES ON LISTS 3,90 (AND 4,90 FOR UHF) CALCULATIONS.  THIS
C  ROUTINE GENERATES A NEW SET OF ORTHONORMALIZED MOLECULAR ORBITALS.
C
C Extended to work with dropmo. Ajith Perera, 07/2005.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      CHARACTER*1 ISP(2)
      DIMENSION ICORE(MAXCOR)
      DIMENSION NOCC_EXPND(8,2),NVRT_EXPND(8,2),
     &          NOCCO_EXPND(2),NVRTO_EXPND(2),
     &          NBF4IRREP(8)
      COMMON /SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF(4)
      COMMON /FLAGS/ IFLAGS(100)
C
      DATA ONE  /1.0/
      DATA ZILCH/0.0/
      DATA ONEM /-1.0/
      DATA ISP /'A','B'/
C
      CALL GETREC(20, "JOBARC", 'NUMDROPA', 1, NDROPMO)
      NBAS=NOCCO(1)+NVRTO(1)
      NBAS2=NBAS*NBAS
      IF (NDROPMO.GT.0) THEN
         NBAS_DROP = NOCCO(1)+NVRTO(1)
         CALL FILL_DROP_MOS(NIRREP, IUHF, NOCC_EXPND, NVRT_EXPND,
     &                      NBAS_FULL, NOCCO_EXPND, NVRTO_EXPND)
         NBAS  = NBAS_FULL
         NBAS2 = NBAS_FULL*NBAS_FULL
      END IF
C
C NOCC_EXPND, POP = The # of occ per irrep in the full and reduced spaces.
C NVRT_EXPND, VRT = The # of vrt per irrep in the full and reduced spaces.
C NOCCO_EXPND, NOCCO = The # of occ in the full and reduced spaces.
C NVRTO_EXPND, NVRTO = The # of vrt in the full and reduced spaces.
C NBAS = The number of basis functions in the full space.
C
      DO 10 ISPIN=1,1+IUHF
       NOCC=NOCCO(ISPIN)
       NVRT=NVRTO(ISPIN)
       IF (NDROPMO.GT.0) THEN
          NOCC=NOCCO_EXPND(ISPIN)
          NVRT=NVRTO_EXPND(ISPIN)
       END IF
       I000=1
       I010=I000+NBAS2*IINTFP
       I020=I010+NBAS2*IINTFP
       I030=I020+NOCC*NVRT*IINTFP
       I040=I030+NOCC*NVRT*IINTFP
       IOFF=IINTFP*(NBAS*NOCC)
C
C GET MO COEFFICIENTS AND T1 VECTOR
C (Always read the MO vectors in the full space.)
C
       IF (NDROPMO.GT.0) THEN
          CALL GETREC(20,'JOBARC','SCFEVC'//ISP(ISPIN)//'0',
     &                NBAS2*IINTFP,ICORE(I000))
       ELSE
          CALL GETREC(20,'JOBARC','SCFEVEC'//ISP(ISPIN),
     &                NBAS2*IINTFP,ICORE(I000))
       END IF
#ifdef _DEBUG_LVLM
        call checksum("@-ROTT1-INITIAL_VECS", ICORE(I000), NBAS2)
#endif
       CALL SCOPY (NBAS2,ICORE(I000),1,ICORE(I010),1)
       CALL GETLST(ICORE(I020),1,1,1,2+ISPIN,90)
C
       IF (NDROPMO.GT.0) THEN
          CALL EXPAND_T1(ICORE(I030),ICORE(I020),NIRREP,
     &                   NOCC_EXPND(1,ISPIN),NVRT_EXPND(1,ISPIN),
     &                   NOCCO(ISPIN),NVRTO(ISPIN),
     &                   NOCCO_EXPND(ISPIN),NVRTO_EXPND(ISPIN),
     &                   SYMP_OVFUL)
          CALL DCOPY(SYMP_OVFUL, ICORE(I030), 1, ICORE(I020), 1)
       END IF
C
C EXPAND T1 VECTOR TO FULL SIZE 
C
       CALL EXPT1(ICORE(I020),ICORE(I030),NOCC,NVRT,
     &            NOCC_EXPND,NVRT_EXPND,NOCCO_EXPND,NVRTO_EXPND,
     &            NDROPMO,ISPIN)
C
C DO ORBITAL ROTATIONS - FIRST GENERATE NEW VIRTUAL ORBITALS
C
       CALL XGEMM ('N','T',NBAS,NVRT,NOCC,ONEM,ICORE(I000),
     &             NBAS,ICORE(I020),NVRT,ONE,
     &             ICORE(I010+IOFF),NBAS)
C
C NOW GENERATE NEW OCCUPIED ORBITALS
C
       CALL XGEMM('N','N',NBAS,NOCC,NVRT,ONE,ICORE(IOFF+I000),
     &            NBAS,ICORE(I020),NVRT,ONE,ICORE(I010),NBAS)
#ifdef _DEBUG_LVLM
      write(*,*) "Then checksum of rotated vectors"
      call checksum("@-ROTT1", ICORE(I010), NBAS2)
#endif        

       CALL PUTREC(20,'JOBARC','SCFEVEC'//ISP(ISPIN),NBAS2*IINTFP,
     &             ICORE(I010))
C
C THESE ORBITALS ARE NOT ORTHONORMAL, AND NOW NEED TO BE ORTHONORMALIZED.
C  USE LOWDIN ORTHONORMALIZATION.
C
C
C GET AO OVERLAP MATRIX AND TRANSFORM TO MO BASIS
C
       CALL GETREC(20,'JOBARC','AOOVRLAP',NBAS2*IINTFP,ICORE(I000))
#ifdef _DEBUG_LVLM
      write(*,*) "Then checksum of the AO. overlap"
      call checksum("@-ROTT1", ICORE(I000), NBAS2)
       call output(ICORE(I000), 1, nbas, 1, nbas, nbas, nbas, 1)
#endif        
       CALL AO2MO2(ICORE(I000),ICORE(I000),ICORE(I010),ICORE(I020),
     &             NBAS,NBAS,ISPIN)
#ifdef _DEBUG_LVLM
      write(*,*) "Then checksum of the trans. overlap"
      call checksum("@-ROTT1", ICORE(I000), NBAS2)
       call output(ICORE(I000), 1, nbas, 1, nbas, nbas, nbas, 1)
#endif        
C
C NOW FORM S**(-1/2) IN THE MO BASIS
C
       CALL EIG(ICORE(I000),ICORE(I010),NBAS,NBAS,0)
       CALL INVSQT(ICORE(I000),NBAS+1,NBAS)
#ifdef _DEBUG_LVLM
       write(*,*) "Then checksum of the eig vals and vectors"
       call checksum("@-ROTT1", ICORE(I010), NBAS2)
       call checksum("@-ROTT1", ICORE(I000), NBAS2)
#endif 
       CALL XGEMM('N','N',NBAS,NBAS,NBAS,ONE,ICORE(I010),NBAS,
     &            ICORE(I000),NBAS,ZILCH,ICORE(I020),NBAS)
       CALL XGEMM('N','T',NBAS,NBAS,NBAS,ONE,ICORE(I020),NBAS,
     &            ICORE(I010),NBAS,ZILCH,ICORE(I000),NBAS)
#ifdef _DEBUG_LVLM
CSSSS       call dcopy (Nbas2, ICORE(I020), 1, ICORE(I000), 1)
       write(*,*) "Then canonical orth. transformation"
       call checksum("@-ROTT1", ICORE(I000), NBAS2)
       call output(ICORE(I000), 1, nbas, 1, nbas, nbas, nbas, 1)
#endif 
C
C NOW TRANSFORM MO VECTOR TO ORTHONORMAL REPRESENTATION
C
       CALL GETREC(20,'JOBARC','SCFEVEC'//ISP(ISPIN),NBAS2*IINTFP,
     &             ICORE(I010))
       CALL XGEMM('N','N',NBAS,NBAS,NBAS,ONE,ICORE(I010),NBAS,
     &            ICORE(I000),NBAS,ZILCH,ICORE(I020),NBAS)
C
C REORDER MO COEFFICIENTS SO THAT SCF CAN READ THEM
C
       CALL PUTREC(20,'JOBARC','SCFEVC'//ISP(ISPIN)//'0',
     &             NBAS2*IINTFP,ICORE(I020))
#ifdef _DEBUG_LVL0
       Write(*,*) "The final vectors written to JOBARC"
       call output(icore(i020), 1, nbas, 1, nbas, nbas, nbas, 1)
C       call checksum("@-ROTT1", ICORE(I020), NBAS2)
C       CALL GETREC(20,'JOBARC','AOOVRLAP',NBAS2*IINTFP,ICORE(I000))
C       call xgemm("T", "N",nbas,nbas,nbas,one,icore(i020),nbas,
C     &             icore(i000),nbas,zilch,icore(i010),nbas)
C
C       call xgemm("N", "N",nbas,nbas,nbas,one,icore(i010),nbas,
C     &             icore(i020),nbas,zilch,icore(i000),nbas)
C       call output(icore(i000), 1, nbas, 1, nbas, nbas, nbas, 1)
#endif
       IF(IUHF.EQ.0)THEN
        CALL PUTREC(20,'JOBARC','SCFEVCB0',
     &             NBAS2*IINTFP,ICORE(I020))
       ENDIF
C
10    CONTINUE
C
C Restore the occupation to the full mo basis. The records OCCUPY(A0|B0)
C are saved in vtran/dropvc.F and contain the original full occupation.
C The same goes for MODROPA0 and REORDEA0. See dropvc for further info.
C
      IF (NDROPMO.GT.0) THEN
cSSS         CALL GETREC(20,'JOBARC','OCCUPYA0',NIRREP,ICORE)
cSSS         CALL PUTREC(20,'JOBARC','OCCUPYA0',NIRREP,ICORE)
cSSS         IF (IUHF.GT.0) THEN
cSSS            CALL GETREC(20,'JOBARC','OCCUPYB0',NIRREP,ICORE)
cSSS            CALL PUTREC(20,'JOBARC','OCCUPYB0',NIRREP,ICORE)
cSSS         END IF
      ELSE
         CALL GETREC(20,'JOBARC','OCCUPYA ',NIRREP,ICORE)
         CALL PUTREC(20,'JOBARC','OCCUPYA0',NIRREP,ICORE)
         IF (IUHF.GT.0) THEN
            CALL GETREC(20,'JOBARC','OCCUPYB ',NIRREP,ICORE)
            CALL PUTREC(20,'JOBARC','OCCUPYB0',NIRREP,ICORE)
         END IF
      END IF
      IF (NDROPMO.GT.0) THEN
         CALL GETREC(20,'JOBARC','MODROPA0',NDROPMO,ICORE(I000))
         CALL PUTREC(20,'JOBARC','MODROPA ',NDROPMO,ICORE(I000))
         CALL GETREC(20,'JOBARC','REORDEA0',NBAS,ICORE(I000))
         CALL PUTREC(20,'JOBARC','REORDERA',NBAS,ICORE(I000))
         IF (IUHF.GT.0) THEN
            CALL GETREC(20,'JOBARC','MODROPB0',NDROPMO,ICORE(I000))
            CALL PUTREC(20,'JOBARC','MODROPB ',NDROPMO,ICORE(I000))
            CALL GETREC(20,'JOBARC','REORDEB0',NBAS,ICORE(I000))
            CALL PUTREC(20,'JOBARC','REORDERB',NBAS,ICORE(I000))
         END IF
      END IF
      CALL PUTREC(20,'JOBARC','IFLAGS  ',100,IFLAGS)
C
      RETURN
      END
