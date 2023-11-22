      SUBROUTINE COMP_2NDORD_PROPS(PERTRB_DEN, SCR, NTPERT, MEMLEFT,  
     &                             NAO, NATOM, JPERT, IUHF, SPN_SPN,
     &                             POLAR, SPNC_TENSOR, POL_TENSOR)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*8 LABELSD(6), LABELPSO(3), LABELDPL(3)
      LOGICAL JFC, JPSO, JSD, POLAR, SPN_SPN
C
#include "mxatms.par"
C
      DIMENSION PERTRB_DEN(NAO*NAO), SCR(MEMLEFT), POL_TENSOR(3, 3),
     &          SPNC_TENSOR(NTPERT, NTPERT)
C
      COMMON /FLAGS/ IFLAGS(100) 
      COMMON /FILES/ LUOUT, MOINTS
      COMMON /SPN_SPN_TYPE/ JFC, JPSO, JSD 
C

      DATA LABELSD/'  SDXX  ','  SDXY  ','  SDXZ  ','  SDYY  ',
     &             '  SDYZ  ','  SDZZ  '/
      DATA LABELPSO/'   OPX  ', '   OPY  ', '   OPZ  '/
      DATA LABELDPL/'     X  ', '     Y  ', '     Z  '/
      DATA IONE, ONE, JZERO /1, 1.0D0, 0/
C
       IONE =  1
      IRWND = 0
      LUOUT = 6
C
      IF (IFLAGS(18) .EQ. 8) THEN
           JPSO = .TRUE.
           JFC  = .FALSE.
           JSD  = .FALSE.
         IFERMI = 0
         ISDIP  = 0
         IPSO   = 1
         IDSO   = 0
         SIGN   = ONE
      ENDIF
      IF (IFLAGS(18) .EQ. 9) THEN
           JFC  = .TRUE.
           JPSO = .FALSE.
           JSD  = .FALSE.
         IFERMI = 1
         ISDIP  = 0
         IPSO   = 0
         IDSO   = 0
         SIGN   = ONE
      ENDIF
      IF (IFLAGS(18) .EQ. 10) THEN
           JSD  = .TRUE.
           JPSO = .FALSE.
           JFC  = .FALSE.
         IFERMI = 0
         ISDIP  = 1
         IPSO   = 0
         IDSO   = 0
         SIGN   = ONE
      ENDIF
C
      SPN_SPN = JFC .OR. JPSO .OR. JSD
      IF (IFLAGS(18) .EQ. 11) THEN
         POLAR = .TRUE.
          FACT = 2.0D0 - DFLOAT(IUHF)
      ENDIF
C
      NSIZE  = NAO*NAO
C
      I000  = IONE + NSIZE
      I010  = I000 + NSIZE
C
C Open the vprop integral file 
C
      OPEN (UNIT=30, FILE='VPOUT', FORM='UNFORMATTED', STATUS='OLD')
C     
#ifdef _DEBUG_LVLM2
      Write(6,*)
      print*, "Perturbation of interest: FC,SD,PSO ", JFC, JSD, 
     &         JPSO,polar, natom, nsize
      Write(6,*) 
      Print*, "The perturbed density"
      call output(pertrb_den, 1, NAO, 1, NAO, NAO, NAO, 1)
#endif

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
               CALL HEADER ('PERTURBED DENSITY', 1, 6)
               CALL TAB (LUOUT, PERTRB_DEN, NAO, NAO, NAO, NAO)
            ENDIF
C     
            SPNC_TENSOR(IPERT, JPERT) = SPNC_TENSOR(IPERT, JPERT) +
     &                                  SIGN*SDOT(NAO*NAO, 
     &                                  PERTRB_DEN, 1, SCR(I000), 1)
#ifdef _DEBUG_LVL0
         write(6,*)
         Write(6,"(a)") "Building The NMR Tensor"
         Write(6,"(2I5,1x,F10.5)") ipert, jpert,
     &                             SPNC_TENSOR(ipert,jpert)
#endif

C
            IRWND = IONE
C     
  10      CONTINUE
C
#ifdef _DEBUG_LVLM2
          Write(6,*)
          Write(6,*) " The spn-spin tensor - intermediate"
          Write(6, '(4F10.5)'), ((SPNC_TENSOR(I, J), I=1, 
     &                            NTPERT), J=1, NTPERT)
#endif
C     
      ELSE IF (JSD) THEN
C     
         ICONT = JZERO
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
                  CALL HEADER ('PERTURBED DENSITY', 1, 6)
                  CALL TAB (LUOUT, PERTRB_DEN, NAO, NAO, NAO, NAO)
               ENDIF
C     
               SPNC_TENSOR(IPERT+ICONT, JPERT) = SIGN*SDOT(NAO*NAO, 
     &                                           PERTRB_DEN, 1, 
     &                                           SCR(I000), 1)
C
               IRWND = IONE
C     
 30         CONTINUE 
            ICONT = ICONT + IONE 
C     
 20      CONTINUE

      ELSE IF (JPSO) THEN
C     
         ICONT = JZERO
C     
         DO 40  IATOMS = 1, NATOM
C     
            DO 50 IPERT = 1, 3
C     
               CALL SEEKLB (LABELPSO(IPERT), IERR, IRWND, 30)
               IF (IERR .NE. 0) CALL ERREX
               CALL LOADINT (SCR(I000), NATOM, NSIZE, NAO, IUHF)
               CALL ANTI_SYM(SCR(I000), NAO)
C
               IF (IFLAGS(1) .GT. 40) THEN
                  CALL HEADER ('PROPERTY INTEGRALS (PSO)', 1, 6)
                  CALL TAB (LUOUT, SCR(I000), NAO, NAO, NAO, NAO)
                  CALL HEADER ('PERTURBED DENSITY', 1, 6)
                  CALL TAB (LUOUT, PERTRB_DEN, NAO, NAO, NAO, NAO)
               ENDIF
C     
               SPNC_TENSOR(IPERT+ICONT, JPERT) = SIGN*SDOT(NAO*NAO, 
     &                                           PERTRB_DEN, 1, 
     &                                           SCR(I000), 1)
               IRWND = IONE
C     
 50         CONTINUE

            ICONT = ICONT + 3
C     
 40      CONTINUE
      
C                 
      ELSE IF (POLAR) THEN

        DO IPERT = 1, 3
C
          CALL SEEKLB (LABELDPL(IPERT), IERR, IRWND, 30)
          IF (IERR .NE. 0) CALL ERREX
          CALL LOADINT (SCR(I000), NATOM, NSIZE, NAO, IUHF)
C
          IF (IFLAGS(1) .GT. 40) THEN
             CALL HEADER ('PROPERTY INTEGRALS (Dipole)', 1, 6)
             CALL TAB (6, SCR(I000), NAO, NAO, NAO, NAO)
             CALL HEADER ('PERTURBED DENSITY', 1, 6)
             CALL TAB (6, PERTRB_DEN, NAO, NAO, NAO, NAO)
          ENDIF
C
          POL_TENSOR(IPERT, JPERT) =  POL_TENSOR(IPERT, JPERT) + 
     &                                SDOT(NAO*NAO, PERTRB_DEN, 1, 
     &                                SCR(I000), 1)*FACT
          IRWND = IONE 
C
         END DO
C
      ENDIF
C  
      RETURN
      END

