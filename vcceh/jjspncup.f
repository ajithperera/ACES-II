C
      SUBROUTINE JJSPNCUP(JJFC, SCR1, SCR2, NTPERT, NIRREP, IEOMPROP, 
     &                    ICALL)
C
C Calculated Fermi-contact NMR J-J coupling constants are transformed
C back to the non-symmetry adapted basis, and multiply with appropriate
C factors to report the coupling constant in MHZ.
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      DOUBLE PRECISION JJFC
      LOGICAL JFC, JPSO, JSD, NUCLEI
      LOGICAL ISOTOPES_PRESENT
      CHARACTER*80 FNAME 
C
      DIMENSION JJFC(NTPERT, NTPERT), SCR1(NTPERT, NTPERT), 
     &          SCR2(NTPERT, NTPERT)
C
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FILES/ LUOUT, MOINTS
      COMMON /NMR/JFC, JPSO, JSD, NUCLEI
C
      DATA AZERO, ONE /0.0D+00, 1.0D+00/
C
      CALL GFNAME("iSOTOPES",FNAME,LENTH)
      INQUIRE(FILE=FNAME(1:LENTH),EXIST=iSOTOPES_PRESENT)
      IUNITIS=1
      IF (iSOTOPES_PRESENT) THEN
         OPEN(UNIT=IUNITIS,FILE=FNAME(1:LENTH),FORM="FORMATTED")
      ENDIF

      IF (NIRREP .NE. 1) THEN
C
C Get the transformation matrix from the 'JOBARC' file. The
C transformation matrix is calculated in symadpt.f based on
C 'VPROP' NMR perturbation integrals.
C
      IF (JFC) THEN
        CALL GETREC (20, 'JOBARC', 'SMTRNFC ', NTPERT*NTPERT*IINTFP,
     &     SCR2)
      ELSE IF (JPSO) THEN
        CALL GETREC (20, 'JOBARC', 'SMTRNPSO', NTPERT*NTPERT*IINTFP,
     &     SCR2)
      ELSE IF (JSD) THEN
        CALL GETREC (20, 'JOBARC', 'SMTRNSD ', NTPERT*NTPERT*IINTFP,
     &     SCR2)
      ELSE IF (NUCLEI) THEN
        CALL GETREC (20, 'JOBARC', 'SMTRNNUC', NTPERT*NTPERT*IINTFP,
     &     SCR2)
      ENDIF
C      
C Now do the actual transformation
C
         CALL XGEMM('N', 'N', NTPERT, NTPERT, NTPERT, ONE, SCR2, 
     &               NTPERT, JJFC, NTPERT, AZERO, SCR1, NTPERT)
C
         CALL XGEMM('N', 'T', NTPERT, NTPERT, NTPERT, ONE, SCR1, 
     &               NTPERT, SCR2, NTPERT, AZERO, JJFC, NTPERT)
      ENDIF
C
      IF (IFLAGS(1) .GT. 20) THEN
            CALL HEADER('Nuclear Spin-Spin Coupling Tensor', -1, LUOUT)
            CALL TAB (6, JJFC, NTPERT, NTPERT, NTPERT, NTPERT)
      ENDIF
C Now we have transformed JJ coupling constants. Multiply by the
C appropriate factors to get the coupling constant in MHZ.
C
      IF (IFLAGS(18) .EQ. 8) THEN
         IFERMI = 0
         ISDIP  = 0
         IPSO   = 1
         IDSO   = 0
      ELSE IF (IFLAGS(18) .EQ. 9) THEN
         IFERMI = 1
         ISDIP  = 0
         IPSO   = 0
         IDSO   = 0
      ELSE IF (IFLAGS(18) .EQ. 10) THEN
         IFERMI = 0
         ISDIP  = 1
         IPSO   = 0
         IDSO   = 0
      ENDIF

      IF (ISOTOPES_PRESENT) THEN
         CALL FACTOR_IS(JJFC, NTPERT, IFERMI, ISDIP, IPSO, IDSO, 
     &                  IEOMPROP, ICALL, IUNITIS)
      ELSE
         CALL FACTOR(JJFC, NTPERT, IFERMI, ISDIP, IPSO, IDSO, 
     &               IEOMPROP, ICALL)
      ENDIF 
C  
      WRITE(6, *)
      CLOSE(IUNITIS)

      RETURN
      END
