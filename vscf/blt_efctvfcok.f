      SUBROUTINE BLT_EFCTVFOCK(EVEC, FOCK, FNFOCK, FNDENS, SCR, NBAS,
     &                         NOCC, NIRREP, NBFIRR, IUHF, POSTSCF)
C
C This routine built the occ-occ block of the Fock matrix and
C transform it to the SA-AO basis. 
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      DIMENSION EVEC(NBAS, NBAS), FOCK(NBAS, NBAS), SCR(NBAS, NBAS),
     &          FNFOCK(NBAS, NBAS), FNDENS(2*NBAS*NBAS),
     &          NBFIRR(8),NOCC(16), IOCC(2, 8), IVRT(2, 8)
C
      LOGICAL POSTSCF, ROHF, SEMICA
      CHARACTER*8 LABELC, LABELF, LABELD, LABELFOCK, LABELDENS, 
     &            LABELPSCF, LABELSDEN, LABELDSCF
C
      COMMON /MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/IFLAGS(100)
C
      DATA ONEM / -1.0D0/, ONE /1.0D0/, TWO /2.0D0/, IONE /1/
C
C Create the population vectors and define target labels. The records
C EFFCFOCK, TDENSITY, DDENSITY, HFDENSTY and HDDENSTY contains the effective 
C Fock matrix, total density, spin density (UHF/ROHF only) and the total
C HF density and spin density (for correlated calculations). 
C
      DO 10 ISPIN = 1, (IUHF + 1)
         DO 20 IRREP = 1, NIRREP
C
            IOCC(ISPIN, IRREP) = NOCC((ISPIN -1)*8 + IRREP)
            IVRT(ISPIN, IRREP) = NBFIRR(IRREP) - IOCC(ISPIN, IRREP)
C
 20      CONTINUE
 10   CONTINUE
C
      ROHF   = (IFLAGS(11) .EQ. 2)
      SEMICA = (IFLAGS(39) .EQ. 1) 
C
      LABELFOCK='EFFCFOCK'
      LABELDENS='TDENSITY'
      LABELSDEN='DDENSITY'
      LABELPSCF='HFDENSTY'
      LABELDSCF='HDDENSTY'
C 
C Post SCF gradient calculations we need SCF density in addition to the
C relaxed density, I(p,q) and Gamma elements. The relaxed density 
C I(p,q) and Gamma elements are created in the density module. Calculate
C the total SCF density here (Alpha and Beta densities are needed for
C UHF/ROHF)
C
      IF (POSTSCF) THEN
C
         DO 30 ISPIN = 1, (IUHF + 1) 
C
            IF (ISPIN .EQ. 1) THEN
               LABELD='SCFDENSA'
            ELSE
               LABELD='SCFDENSB'
            ENDIF
C
            CALL GETREC(20, 'JOBARC', LABELD, NBAS*NBAS*IINTFP, EVEC)
            IF (ISPIN .EQ. 1) CALL ZERO(FNDENS, 2*NBAS*NBAS)
            CALL SAXPY(NBAS*NBAS, ONE, EVEC, 1, FNDENS, 1)
            IF (IUHF .NE. 0) THEN
               IF (ISPIN .EQ. 1) THEN
                  CALL SAXPY(NBAS*NBAS, ONE, EVEC, 1, 
     &                       FNDENS(NBAS*NBAS+1), 1)
               ELSE
                  CALL SAXPY(NBAS*NBAS, ONEM, EVEC, 1,
     &                       FNDENS(NBAS*NBAS+1), 1)
               ENDIF
            ENDIF
C
 30      CONTINUE
C
         CALL PUTREC(20, 'JOBARC', LABELPSCF, NBAS*NBAS*IINTFP, FNDENS)
         CALL PUTREC(20, 'JOBARC', LABELDSCF, NBAS*NBAS*IINTFP, FNDENS
     &               (NBAS*NBAS+1))
C
         RETURN
C
      ENDIF
C
C The only SCF logic begin here. 
C
      DO 40 ISPIN = 1, (IUHF + 1)
C
         IF(ISPIN.EQ.1) THEN
            LABELC='SCFEVCA0'
            LABELF='FOCKA   '
            LABELD='SCFDENSA'
         ELSE
            LABELC='SCFEVCB0'
            LABELF='FOCKB   '
            LABELD='SCFDENSB'
         ENDIF
C
C Get the eigenvectors from the JOBARC file.
C     
         CALL GETREC(20, 'JOBARC', LABELC, NBAS*NBAS*IINTFP, EVEC)
C
C Now, get the Fock matrix. For ROHF standard orbitals (not semi canonical) we have 
C additional work to do after the transformation.
C
         CALL GETREC(20, 'JOBARC', LABELF, NBAS*NBAS*IINTFP, FOCK)
C
C Transform the Fock matrix to the MO basis         
C
         CALL XGEMM('T', 'N', NBAS, NBAS, NBAS, 1.0D0, EVEC, NBAS, 
     &               FOCK, NBAS, 0.0D0, SCR, NBAS)
         CALL XGEMM('N', 'N', NBAS, NBAS, NBAS, 1.0D0, SCR, NBAS, 
     &               EVEC, NBAS, 0.0D0, FOCK , NBAS)
C
C Now, zero out the diagonal virtual orbital eigenvalues for each irrep for RHF and UHF.
C
         IOFF = IONE
         IBEGN = IONE
         CALL ZERO(SCR, NBAS*NBAS)
C
         DO 5 IRREP = 1, NIRREP
C
            IF (.NOT. (ROHF .AND. (.NOT. SEMICA))) THEN 
C
               IOFF1 = IOFF + IOCC(ISPIN, IRREP)
C
               CALL SAXPY (IVRT(ISPIN, IRREP), ONEM, FOCK(IOFF1, IOFF1),
     &                     NBAS+1, FOCK(IOFF1, IOFF1), NBAS+1)
C
               IOFF = IOFF + NBFIRR(IRREP)
C
C For ROHF we need to do something clever. Unlike the RHF/UHF the 
C virtual-occupied block is not digonal. So, simply removing the
C virtual orbital eigenvalues do not work.  
C
            ELSE
C
               NELMNTS = IOCC(ISPIN, IRREP)
C     
               DO 99 IELMNTS = 1, NELMNTS
C
                  CALL SCOPY(NELMNTS, FOCK(IBEGN, (IBEGN+IELMNTS-1)),
     &                       IONE, SCR(IBEGN, (IBEGN+IELMNTS-1)), IONE)
C
 99            CONTINUE
C
               IBEGN = IBEGN + NBFIRR(IRREP)
C
            ENDIF
C
 5       CONTINUE
C
         IF (ROHF .AND. (.NOT. SEMICA)) THEN
            CALL ZERO(FOCK, NBAS*NBAS)
            CALL SCOPY(NBAS*NBAS, SCR, IONE, FOCK, IONE) 
         ENDIF
C
C Transform back to the AO basis.          
C
         CALL XGEMM('N', 'N', NBAS, NBAS, NBAS, 1.0D0, EVEC, NBAS, 
     &               FOCK, NBAS, 0.0D0, SCR, NBAS)
         CALL XGEMM('N', 'T', NBAS, NBAS, NBAS, 1.0D0, SCR, NBAS, 
     &               EVEC, NBAS, 0.0D0, FOCK, NBAS)
C
         IF (ISPIN .EQ. 1) CALL ZERO(FNFOCK, NBAS*NBAS)
         CALL SAXPY(NBAS*NBAS, ONE, FOCK, 1, FNFOCK, 1)
C
C For UHF gradient calculations we need Alpha and Beta density difference 
C in additions to the total effective density (for RHF one only needs the
C total density. Note that EVEC array is used  to store the density.    
C         
         CALL GETREC(20, 'JOBARC', LABELD, NBAS*NBAS*IINTFP, EVEC)

         IF (ISPIN .EQ. 1) CALL ZERO(FNDENS, 2*NBAS*NBAS)
         CALL SAXPY(NBAS*NBAS, ONE, EVEC, 1, FNDENS, 1)
         IF (IUHF .NE. 0) THEN
            IF (ISPIN .EQ. 1) THEN
               CALL SAXPY(NBAS*NBAS,ONE,EVEC,1,FNDENS(NBAS*NBAS+1),1)
            ELSE
               CALL SAXPY(NBAS*NBAS,ONEM,EVEC,1,FNDENS(NBAS*NBAS+1),1)
            ENDIF
         ENDIF
C     
 40   CONTINUE
C
      IF (IUHF .EQ. 0) CALL SSCAL(NBAS*NBAS, TWO, FNFOCK, 1)
C
C Write the total density (diffrence for UHF) and effective Fock matrix to the 
C JOBARC file. Note that the RHF density is not scaled by two.
C This is because the DENSA is already scaled (note the inconsistent 
C nomenclature of JOBARC records). The density code write the same records 
C for correlated calculations. 
C
      CALL PUTREC(20, 'JOBARC', LABELFOCK, NBAS*NBAS*IINTFP, FNFOCK)
      CALL PUTREC(20, 'JOBARC', LABELDENS, NBAS*NBAS*IINTFP, FNDENS)
      CALL PUTREC(20, 'JOBARC', LABELSDEN, NBAS*NBAS*IINTFP, FNDENS
     &            (NBAS*NBAS+1))
C
C Also write ''HFDENSTY'' and "HDDENSTY" density records. For SCF calculations
C these records are identical to "TDENSITY" and "DDENSITY" For correlated
C calculations "HFDENSTY" and HDDENSTY records remains as the HF density and
C spin density while the "TDENSITY" and DDENSITY record contains the 
C correlated "relaxed" density matrix and it's difference. This is simply
C to avoid additional logic in the ALASKA code.
C
      CALL PUTREC(20, 'JOBARC', LABELPSCF, NBAS*NBAS*IINTFP, FNDENS)
      CALL PUTREC(20, 'JOBARC', LABELDSCF, NBAS*NBAS*IINTFP, FNDENS
     &            (NBAS*NBAS+1))
C
cSSS      Write(6,*) "The final Fock"
cSSS      CALL OUTPUT(FNFOCK, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
CSSS      Write(6,*) "The final Dens"
CSSS      CALL OUTPUT(FNDENS, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
C
      RETURN
      END
