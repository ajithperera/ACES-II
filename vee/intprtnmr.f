      SUBROUTINE INTPRTNMR (TRNDEN, SCR, TMNMR, MXCOR, NAO, NMO,
     &                      ISIDE, ISPIN, NTPERT, MAXCENT, NATOM,
     &                      FACT)
C
C     This routine calculates the transition moments over FC,
C     SD and PSO operators. This is implemented as a tool 
C     to interpret NMR spin-spin coupling constants.
C     S. Ajith Perera 04/02/96
C
C     TRNDEN : Left and right transition density matrix is 
C              stored depending on whether ISIDE = 1 or 2.
C     SCR    : Scratch space. Eventually, It will be used to
C              store NMR integrals.
C     MXCOR  : Available core space. 
C       NAO  : Number of atomic orbatils in the basis set. 
C       NMO  : Number of MOs. Differs from NAO only in the
C              case of dropcore calcualtions.
C     ISIDE  : The left/right hand side flag. ISIDE = 1 (LHS)
C            : and ISIDE =2 (RHS). 
C     IUHF   : The UHF/RHF flag. IUHF = 0 (RHF) and IUHF = 1 (UHF).
C     NTPERT : Total number of perturbations.
C     NATOM  : Total number of atoms.
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*8 LABELSD(6), LABELPSO(3)
      LOGICAL JFC, JPSO, JSD
C
      DIMENSION TRNDEN(NAO*NAO), SCR(MXCOR), TMNMR(6*MAXCENT)
C
      COMMON /FLAGS/ IFLAGS(100) 
      COMMON /FILES/ LUOUT, MOINTS
C

      DATA LABELSD/'  SDXX  ','  SDXY  ','  SDXZ  ','  SDYY  ',
     &             '  SDYZ  ','  SDZZ  '/
      DATA LABELPSO/'   OPX  ', '   OPY  ', '   OPZ  '/
C
      IONE =  1
      IRWND = 0
C
      IF (IFLAGS(18) .EQ. 8) THEN
         JPSO = .TRUE.
         JFC  = .FALSE.
         JSD  = .FALSE.
      ENDIF
      IF (IFLAGS(18) .EQ. 9) THEN
         JFC  = .TRUE.
         JPSO = .FALSE.
         JSD  = .FALSE.
      ENDIF
      IF (IFLAGS(18) .EQ. 10) THEN
         JSD  = .TRUE.
         JPSO = .FALSE.
         JFC  = .FALSE.
      ENDIF
C
      IF (ISPIN .EQ. 1 .AND. (JPSO .OR. JFC .OR. JSD)) THEN
         SIGN = 1.0D+00
      ELSE IF (ISPIN .EQ. 2 .AND. (JFC .OR. JSD)) THEN
         SIGN = - 1.0D+0
      ELSE IF (ISPIN .EQ. 2 .AND. JPSO) THEN
         SIGN = 0.0D+00
      ENDIF
C
      NSIZE  = NAO*NAO
C
      I000  = IONE
      I010  = I000 + NSIZE
C
C Open the vprop integral file 
C
      OPEN (UNIT=30, FILE='VPOUT', FORM='UNFORMATTED', STATUS='OLD')
C     
      IF (JFC) THEN 
C     
         DO 10 IPERT = 1, NATOM
C
            CALL SEEKLB ('   DEN  ', IERR, IRWND, 30)
            IF (IERR .NE. 0) CALL ERREX
            CALL LOADINT (SCR(I000), NATOM, NSIZE, NAO, IUHF)
C
            IF (IFLAGS(1) .GT. 40) THEN
               CALL HEADER ('PROPERTY INTEGRALS (FC)', 1, 6)
               CALL TAB (LUOUT, SCR(I000), NAO, NAO, NAO, NAO)
               CALL HEADER ('TRANSITION DENSITY', 1, 6)
               CALL TAB (LUOUT, TRNDEN, NAO, NAO, NAO, NAO)
            ENDIF
C     
            TMNMR(IPERT) = TMNMR(IPERT) + SIGN*SDOT(NAO*NAO, 
     &                     TRNDEN, 1, SCR(I000), 1)*FACT
C     
            IRWND = IONE
C     
  10      CONTINUE
C     
      ELSE IF (JSD) THEN
C     
         ICONT = IONE
C     
         DO 20  IATOMS = 1, NATOM
C     
            DO 30 IPERT = 1, 6
C     
               CALL SEEKLB (LABELSD(IPERT), IERR, IRWND, 30)
               IF (IERR .NE. 0) CALL ERREX
               CALL LOADINT (SCR(I000), NATOM, NSIZE, NAO, IUHF)
C     
               IF (IFLAGS(1) .GT. 40) THEN
                  CALL HEADER ('PROPERTY INTEGRALS (SD)', 1, 6)
                  CALL TAB (LUOUT, SCR(I000), NAO, NAO, NAO, NAO)
                  CALL HEADER ('TRANSITION DENSITY', 1, 6)
                  CALL TAB (LUOUT, TRNDEN, NAO, NAO, NAO, NAO)
               ENDIF
C     
               TMNMR(ICONT) = TMNMR(ICONT) + SIGN*SDOT(NAO*NAO, 
     &                        TRNDEN, 1, SCR(I000), 1)*FACT
C     
               ICONT = ICONT + IONE 
               IRWND = IONE
C     
 30         CONTINUE 
C     
 20      CONTINUE

      ELSE IF (JPSO) THEN
C     
         ICONT = IONE
C     
         DO 40  IATOMS = 1, NATOM
C     
            DO 50 IPERT = 1, 3
C     
               CALL SEEKLB (LABELPSO(IPERT), IERR, IRWND, 30)
               IF (IERR .NE. 0) CALL ERREX
               CALL LOADINT (SCR(I000), NATOM, NSIZE, NAO, IUHF)
C
               IF (IFLAGS(1) .GT. 40) THEN
                  CALL HEADER ('PROPERTY INTEGRALS (PSO)', 1, 6)
                  CALL TAB (LUOUT, SCR(I000), NAO, NAO, NAO, NAO)
                  CALL HEADER ('TRANSITION DENSITY', 1, 6)
                  CALL TAB (LUOUT, TRNDEN, NAO, NAO, NAO, NAO)
               ENDIF
C     
               TMNMR(ICONT) = TMNMR(ICONT) + SIGN*SDOT(NAO*NAO, 
     &                        TRNDEN, 1, SCR(I000), 1)*FACT
C     
               ICONT = ICONT + IONE
               IRWND = IONE
C     
 50         CONTINUE
C     
 40      CONTINUE
C     
      ENDIF
C     
      RETURN
      END

