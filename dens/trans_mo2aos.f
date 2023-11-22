      SUBROUTINE TRANS_MO2AOS(CMO, FT, DV, DSO, FSO, SCR, IUHF, NONHF)
C
C This subroutine calculates the folded total AO one-electron relaxed density
C and I(mu,nu) intermediate. See below for the citations. Ajith Perera 07/2000
C
C    CMO ..... SCF eigen vectors
C    FT  ..... The I(p, q) intermediates.
C    DV  ..... Hold the total density matrix.
C    DSO ..... The folded symmetry adapted AO basis total density matrix. 
C    FSO ..... The folded symmetry adapted AO basis total I intermediate.
C    SCR ..... Scratch array.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 LABELC,LABELD,LABELI,LABELNHF, LABELFOCK, LABELDENS,
     &            LABELSDEN
      LOGICAL NONHF
      INTEGER R,S,RS,U,V,UV
      INTEGER POP,VRT,DIRPRD
      DIMENSION CMO(1), FT(1), DV(1), DSO(1), FSO(1),SCR(1)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
      COMMON/FLAGS/IFLAGS(100)
C
      DATA AZERO,HALF,ONE,TWO /0.0D0,0.5D0,1.0D0,2.0D0/
C
C Initilize the target indices and target arrays.
C The records  EFFCFOCK, TDENSITY, DDENSITY contains the I intermediate
C matrix, total relax density and spin relax density (UHF/ROHF only).

      LABELFOCK='EFFCFOCK'
      LABELDENS='TDENSITY'
      LABELSDEN='DDENSITY'

      CALL ZERO(DSO, 2*NBASTT)
      CALL ZERO(FSO, NBASTT)
C
      DO 10 ISPIN=1, 1 + IUHF
C
C Set the labels for correlated methods. We need SCF eiegnvectors, the relax
C density and the I(pq) intermediate (see JFS, JG and RJB 95, 2623, 1991)and
C the density for the non-HF cases (ROHF/QRHF/KS)
C
         IF(ISPIN .EQ. 1) THEN
            LABELC='SCFEVCA0'
            LABELI='IINTERMA'
            LABELD='RELDENSA'
            LABELNHF='NHFDENSA'
         ELSE
            LABELC='SCFEVCB0'
            LABELI='IINTERMB'
            LABELD='RELDENSB'
            LABELNHF='NHFDENSB'
         ENDIF
C
C Get the SCF eigen vectors from the JOBARC file, They are in symmetry
C Blocked form. Symmetry pack them for future use.
C
         CALL GETREC(20, 'JOBARC', LABELC, IINTFP*NBASIS*NBASIS, SCR)
         CALL SYMC(SCR, CMO, NBASIS, NBAS, .FALSE., .FALSE., ISPIN)
C
C Read the Relaxed Density from the JOBARC file (The density code put it
C in there). For non-HF we need additional part - read and add that to the
C relax density to get the total non-HF density. For a discussion of
C non-HF gradients see JG, JFS, RJB JCP, 
C
         CALL GETREC(20, 'JOBARC', LABELD, IINTFP*NBASIS*NBASIS, SCR)
C
         IF (NONHF) THEN
            CALL GETREC(-20, 'JOBARC', LABELNHF, IINTFP*NBASIS*NBASIS,
     &                  DV)
            CALL SAXPY(NBASIS*NBASIS, ONE, DV, 1, SCR, 1)
         ENDIF
C
C Symmetry pack the total one-particle density.
C
         CALL SYMD(SCR, DV, NBASIS, NBAS, ISPIN)
C
C Read in the I(p,q) intermediates for correlated gradients and symmetry pack 
C them.
C
         CALL GETREC(20, 'JOBARC', LABELI, IINTFP*NBASIS*NBASIS, SCR)
C
         CALL SYMI(SCR, FT, NBASIS, NBAS, ISPIN)
C
cSSS         Write(6,*) "SCF Molecular Orbitals"  
cSSS         CALL OUTPUT(CMO, 1, NBASIS, 1,  NBASIS, NBASIS, NBASIS, 1)
cSSS         Write(6,*) "Total Density in MO basis"  
cSSS         CALL EXPND2(DV, SCR, NBASIS)
cSSS         CALL OUTPUT(SCR, 1, NBASIS, 1,  NBASIS, NBASIS, NBASIS, 1)
cSSS         Write(6,*) "Intermediate I(p,q) in MO basis"  
cSSS         CALL EXPND2(FT, SCR, NBASIS)
cSSS         CALL OUTPUT(SCR, 1, NBASIS, 1,  NBASIS, NBASIS, NBASIS, 1)
C
         ISEND=0
         ICEND=0
         IFEND=0
         IOFFO=0     
         IOFFF=0
C
C Do the transformation to symm. adapted AO basis.
C 
         DO 110 ISYM =1, NIRREP
C
            NOCCI=POP(ISYM, ISPIN)
            NBASI=NBAS(ISYM)
C
            IF(NBASI .EQ. 0) GOTO 120
            RS=0
C
C Loop over the target indices R and S (these are final the indices of symm. adap.
C AO basis D and I quantities.)
C
            DO 100 R=1, NBASI
C
               DO 200 S=1, R
                  RS = RS + 1
C
                  DTRS=AZERO
                  FTRS=AZERO
C
                  IF (NBASI .GT. 0) THEN
C
                     UV = ISEND
                     IDVEND = ICEND 
                     ICENDU = IDVEND
C
C Loop over all pairs of Mos and transform Relaxed density (DV) and I(p,q) (FT) 
C intermediates.
C
                     DO 400 U=1, NBASI
                        ICENDV = IDVEND
C
                        DO 410 V=1, U
C
                           UV=UV+1
                           DUV=DV(UV)
                           TEMP=CMO(ICENDU+R)*CMO(ICENDV+S)
                           IF(U.NE.V) TEMP=TEMP+CMO(ICENDU+S)*
     &                                     CMO(ICENDV+R)
                           DTRS=DTRS+DUV*TEMP
                           ICENDV=ICENDV+NBASI
C
 410                    CONTINUE
C
                        ICENDU=ICENDU+NBASI
C
 400                 CONTINUE
C
                     UV=IFEND
                     IDVEND=ICEND
                     ICENDU=IDVEND 
C
                     DO 500 U=1,NBASI
                        ICENDV = IDVEND
C
                        DO 510 V=1,NBASI
                           UV=UV+1
                           FUV = FT(UV)
                           TEMP= CMO(ICENDU+R)*CMO(ICENDV+S)
                           IF(R.NE.S) TEMP=TEMP+CMO(ICENDU+S)*
     &                                     CMO(ICENDV+R)
                           FTRS=FTRS-FUV*TEMP
                           ICENDV = ICENDV + NBASI
C
 510                    CONTINUE
C
                        ICENDU=ICENDU+NBASI
C
 500                 CONTINUE  
                  ENDIF
C
                  IF (R .NE. S) THEN
cSSS                     DTRS = DTRS + DTRS
                     FTRS = HALF*FTRS 
                  ENDIF
C
                  DSO(ISEND+RS) = DSO(ISEND+RS) + DTRS
C
C Relaxed density difference for UHF/ROHF calcualtions.
C
                  IF (IUHF .NE. 0) THEN
                     IF (ISPIN .EQ. 1) THEN
                        DSO(ISEND+RS+NBASTT) = DTRS
                     ELSE
                        DSO(ISEND+RS+NBASTT) = DSO(ISEND+RS+NBASTT) 
     &                                         - DTRS
                     ENDIF
                  ENDIF
C
                  FSO(ISEND+RS) = FSO(ISEND+RS) + FTRS
C
 200           CONTINUE
C     
 100        CONTINUE
C
cSSS            WRITE (6,'(1X,A,I5)') ' Symmetry', ISYM
cSSS            CALL PRITRI(DSO(ISEND+1+NBASTT),NBASI,
cSSS     *           'Total density matrix (SO basis)')
cSSS            CALL PRITRI(FSO(ISEND+1),NBASI,
cSSS     *           'Total Fock matrix (SO basis)')
cSSS
 120        CONTINUE
C
            IFEND = IFEND + NBASI*NBASI
            ISEND = ISEND + (NBASI*(NBASI + 1))/2
            ICEND = ICEND + NBASI*NBASI
            IOFFO = IOFFO + NBASI
            IOFFF = IOFFF + NOCCI*NOCCI
C
 110     CONTINUE
C     
 10   CONTINUE
C
C Expand the the AO basis relaxed density and I(mu.nu) intermediates to the full square
C matrix form before dump into the JOBARC file. The reason is merely to stay 
C consitent with the SCF records TDENSITY, DDENSITY and EFFCFOK written for SCF
C gradient calculations in the vscf module. For correlated calculations those records 
C are not written by the vscf module. With these records gradient codes does not
C need to have any logic about whether it is SCF/Correlated
C or RHF/UHF/ROHF or something else. It only has to look for these three records. 
C We need to make sure that the right things are on these records, then gradients should 
C work (Two particle density is separtely calculated)
C
      CALL EXPND_STR_TO_SQ(DSO, SCR, NBASIS, NBAS, NIRREP)
CSSS      Write(6,*) "Total Density in SAO basis"  
CSSS      CALL OUTPUT(SCR, 1, NBASIS, 1,  NBASIS, NBASIS, NBASIS, 1)
      CALL PUTREC(20, 'JOBARC', LABELDENS, NBASIS*NBASIS*IINTFP, SCR)
C
      CALL EXPND_STR_TO_SQ(DSO(NBASTT+1), SCR, NBASIS, NBAS, NIRREP)
cSSS      Write(6,*) "Total Density difference in SAO basis"  
cSSS      CALL OUTPUT(SCR, 1, NBASIS, 1, NBASIS, NBASIS, NBASIS, 1)
      CALL PUTREC(20, 'JOBARC', LABELSDEN, NBASIS*NBASIS*IINTFP, SCR)

C
      CALL EXPND_STR_TO_SQ(FSO, SCR, NBASIS, NBAS, NIRREP)
cSSS      Write(6,*) "Intermediate I(p,q) in SAO basis"  
cSSS      CALL OUTPUT(SCR, 1, NBASIS, 1,  NBASIS, NBASIS, NBASIS, 1)
      CALL PUTREC(20, 'JOBARC', LABELFOCK, NBASIS*NBASIS*IINTFP, SCR)
C
      RETURN
      END

