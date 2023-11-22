C
      SUBROUTINE FACTOR(XMATX, N, IFERMI, ISDIP, IPSO, IDSO)
C
C Calculates the factors for NMR spin-spin coupling constant.
C The SI system of units has been used. Coded by Ajith 10/93.                             
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      CHARACTER*6 ANAME(45), LABEL(45)
      CHARACTER*51 STRINGFC
      CHARACTER*48 STRINGPSO
      CHARACTER*49 STRINGSD
C
      DIMENSION XMATX(N,N),SPINJJ(100, 100) 
      DIMENSION GTAB(45), ICHARGE(100), ITYPE(100)
C
      COMMON/FLAGS/IFLAGS(100)
C
      DATA IONE/1/
      DATA THREE/3.D0/
C  
C  Nuclei g-factors (per second per tesla). Extracted from CRC handbook
C  of chemistry and physics, 71'st edition 1990-1991
C  Update whenever possible.

C                      1H                3HE              7LI        
      DATA GTAB/ 2.675202874D+08, -2.037995287D+08,  1.039754808D+08,
C                      9BE               11B              13C
     &          -3.759981028D+07,  8.584481142D+07,  6.728214013D+07,
C                      15N               17O              19F       
     &          -2.712607916D+07, -3.628049629D+07,  2.518130433D+08,
C                      21NE              23Na             25Mg
     &          -2.113073577D+07,  7.080361014D+07,  1.638829367D+07,
C                      27AL              29SI             31P      
     &           6.976227164D+07, -5.319083216D+07,  1.083932031D+08,
C                      33S               35CL             39AR  
     &           2.055664899D+07,  2.624164069D+07,  1.778911832D+07,
C                      39K               43CA             45SC   
     &           1.249899943D+07, -1.802585043D+07,  6.508751029D+07,
C                      47TI              51V              53CR   
     &          -1.510531509D+07,  7.045456939D+07, -1.515167627D+07,
C                      55MN              57FE             59CO   
     &           6.645166199D+07,  8.678352950D+06,  6.335662908D+07,
C                      61NI              63CU             67ZN   
     &          -2.394752862D+07,  7.098816083D+07,  1.677240178D+07,
C                      69GA              73GE             75AS     
     &           6.438807865D+07, -9.360222020D+06,  4.596110641D+07,  
C                      77SE              79BR             83KR 
     &           5.125209194D+07,  6.725563891D+07, -1.033120987D+07,
C                                  SPECIAL NUCLEI
C                      14N               37CL             49TI  
     &           1.933759265D+07,  2.184339522D+07, -1.510939290D+07,
C                     2H                6LI              10B    
     &           4.106604279D+07,  3.937108191D+07,  2.587126212D+08,
C                      65CU              71GA             81BR 
     &           7.604574401D+07,  8.181119726D+07,  7.249719806D+07/
C
      DATA ANAME/ '[1] H ', '[3]HE ', '[7]LI ', '[9]BE ', '[11]B ',
     &            '[13]C ', '[15]N ', '[17]O ', '[19]F ', '[21]NE',
     &            '[23]NA', '[25]MG', '[27]AL', '[29]SI', '[31]P ',
     &            '[33]S ', '[35]CL', '[39]AR', '[39]K ', '[43]CA',
     &            '[45]SC', '[47]TI', '[51]V ', '[53]CR', '[55]MN',
     &            '[57]FE', '[59]CO', '[61]NI', '[63]CU', '[67]ZN',
     &            '[69]GA', '[73]GE', '[75]AS', '[77]SE', '[79]BR',
     &            '[83]KR', '[14]N ', '[37]CL', '[49]TI', '[2] H ',
     &            '[6]LI ', '[10]B ', '[65]CU', '[71]GA', '[81]BR'/ 
C
C Conversion of calculated coupling constant to usual units (HZ)
C Based on SI system of units. Conversion factors to Hz from
C atomic units is described as follows. Formulas are based on 
C Jan Geertsen et al. Chem. Phys. 90, 301-311, 1984 and 
C H. Fukui et al. JCP 97, 2299-2304. Physical constants are from 
C CRC handbook of chemistry and physics, 71'st edition 1990-1991.
C  
C  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  & Fermi-contact = -(2)(16)(mu_o)^2(mu_B)^2(hbar)^2/(9)(8)(PI)^2)(h)   &
C  &                                                                     &
C  & Spin-dipole   = -2(mu_O)^2(mu_B)^2(hbar)^2/(4)(8)(PI)^2)(h)         & 
C  &                                                                     &
C  & Para. magnetic spin-orbit = -2(mu_o)^2(mu_B)^2(hbar)^2/(2)(PI)^2)(h)&
C  &                                                                     &
C  & Diamgnetic spin-orbit = (2)(mu_o)^2(mu_B)(e)(hbar)/(8)(PI)^2        &
C  & &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      FACTFC  = -1.058311758D-14
      FACTPSO = -6.031653750D-16
      FACTDSO =  1.206313128D-15
      FACTSD  = -1.507913440D-16
C 
      WRITE(6,*)
C
C Get the atomic informations from JOBARC file
C
      CALL GETREC(20, 'JOBARC', 'NATOMS  ', IONE, NCENTR)
      CALL GETREC(20, 'JOBARC', 'NREALATM', 1, NATOMS) 
      CALL GETREC(20, 'JOBARC', 'VMOLORDR', NATOMS, ITYPE)
C
      K = 1
      DO 2 I = 1, NCENTR
         IF (ITYPE(I) .ne. 0) THEN 
            ICHARGE(K) = ITYPE(I)
            K = K + 1
         ENDIF
 2    CONTINUE
C
C Set the output strings depending on the method.
C
      IF (IFERMI .NE. 0) THEN
        STRINGFC='Transition moments over Fermi-contact operator (Hz)'
      ELSE IF (IPSO .NE. 0) THEN
        STRINGPSO='Transition moments over spin-orbit operator (Hz)'
      ELSE IF (ISDIP .NE. 0) THEN
        STRINGSD='Transition moments over spin-dipole operator (Hz)'
      ENDIF
C
C Take care of the g values for special nuclei
C
C      DO 5 I = 1, 100
C         IF (INPGVAL(I) .NE. 0) ICHARGE(I) = INPGVAL(I)
C 5    CONTINUE
C      
      IF(IFERMI .NE. 0) THEN
         NCOMPO = 1
         IF (IFLAGS(1) .GE. 20) THEN
         WRITE(6,1000) 'Conversion factor for Fermi-contact term = ',
     &                  FACTFC
         ENDIF
         NNUCL = N/NCOMPO
         DO 20 I = 1, NNUCL
            LABEL(I) = ANAME(ICHARGE(I))
            GFCI = GTAB(ICHARGE(I))
            IF (IFLAGS(1) .GE. 20) THEN
            WRITE(6, 1200)'Nuclear g-factor for nuclei', LABEL(I), 
     &                    ' = ', GFCI
            ENDIF
            DO 30 J = 1, I
               GFCJ = GTAB(ICHARGE(J))
               SPINJJ(I, J) = XMATX(I, J)*GFCI*GFCJ*FACTFC
 30         CONTINUE
 20      CONTINUE
         WRITE(6,*)
         CALL HEADER(STRINGFC, -1, 6)
         CALL PRNTLO (6, SPINJJ, LABEL, NNUCL, NNUCL, 100, 100)
      END IF
C
      IF(IPSO .NE. 0) THEN 
         NCOMPO = 3
         IF (IFLAGS(1) .GE. 20) THEN
         WRITE(6, 1000) 'Conversion factor for spin-orbit term = ', 
     &                   FACTPSO
         ENDIF
         NNUCL = N/NCOMPO
         DO 40 I = 1, NNUCL
            LABEL(I) = ANAME(ICHARGE(I))
            GFCI = GTAB(ICHARGE(I))
            IF (IFLAGS(1) .GE. 20) THEN
            WRITE(6, 1300)'Nuclear g-factor for nuclei', LABEL(I), 
     &                    ' = ', GFCI
            ENDIF
            DO 50 J = 1, I
               SPINJJ(I, J) = 0.0D+00
               GFCJ = GTAB(ICHARGE(J))
               DO 60 KK = 1, NCOMPO
                  II = (I - 1)*NCOMPO + KK
                  JJ = (J - 1)*NCOMPO + KK
                  SPINJJ(I, J) = SPINJJ(I, J) + XMATX(II, JJ)*GFCI*GFCJ
     &                           *FACTPSO/THREE 
 60            CONTINUE
 50         CONTINUE
 40      CONTINUE
               WRITE(6, *)
               CALL HEADER(STRINGPSO, -1, 6)
               CALL PRNTLO (6, SPINJJ, LABEL, NNUCL, NNUCL, 100, 100)
      END IF
C
      IF(ISDIP .NE. 0) THEN
         NCOMPO = 6
         IF (IFLAGS(1) .GE. 20) THEN
         WRITE(6,1000) 'Conversion factor for spin-dipole term = ', 
     &               FACTSD
         ENDIF
         NNUCL = N/NCOMPO
         DO 70 I = 1, NNUCL
            LABEL(I) = ANAME(ICHARGE(I))
            GFCI = GTAB(ICHARGE(I))
            IF (IFLAGS(1) .GE. 20) THEN
            WRITE(6, 1400)'Nuclear g-factor for nuclei', LABEL(I), 
     &                    ' = ', GFCI
            ENDIF
            DO 80 J = 1, I
               SPINJJ(I, J) = 0.0D+00
               GFCJ = GTAB(ICHARGE(J))
               DO 90 KK = 1, NCOMPO
                  II = (I - 1)*NCOMPO + KK
                  JJ = (J - 1)*NCOMPO + KK
                  SPINJJ(I, J) = SPINJJ(I, J) + XMATX(II, JJ)*GFCI*GFCJ
     &                           *FACTSD/THREE 
                  IF (KK .EQ. 2 .OR. KK .EQ. 3 .OR. KK .EQ. 5) THEN
                     SPINJJ(I, J) = SPINJJ(I, J) + XMATX(II, JJ)*GFCI*
     &                              GFCJ*FACTSD/THREE 
                  ENDIF
 90            CONTINUE
 80         CONTINUE
 70      CONTINUE 
         WRITE(6, *) 
         CALL HEADER(STRINGSD, -1, 6)
         CALL PRNTLO (6, SPINJJ, LABEL, NNUCL, NNUCL, 100, 100)
      END IF
C
      IF(IDSO .NE. 0) THEN 
         NCOMPO = 3
         IF (IFLAGS(1) .GE. 20) THEN
         WRITE(6, 1000)'Conversion factor for diam. spin-orbit term = '
     &                   , FACTDSO
         ENDIF
         NNUCL = N/NCOMPO
         DO 100 I = 1, NNUCL
            LABEL(I) = ANAME(ICHARGE(I))
            GFCI = GTAB(ICHARGE(I))
            IF (IFLAGS(1) .GE. 20) THEN
            WRITE(6, 1600)'Nuclear g-factor for nuclei', LABEL(I), 
     &                    ' = ', GFCI
            ENDIF
            DO 120 J = 1, I
               SPINJJ(I, J) = 0.0D+00
               GFCJ = GTAB(ICHARGE(J))
               DO 130 KK = 1, NCOMPO
                  II = (I - 1)*NCOMPO + KK
                  JJ = (J - 1)*NCOMPO + KK
                  SPINJJ(I, J) = SPINJJ(I, J) + XMATX(II, JJ)*GFCI*GFCJ
     &                           *FACTDSO/THREE 
 130           CONTINUE
 120        CONTINUE
 100     CONTINUE
               WRITE(6, *)
               CALL HEADER('DSO contribution to J (in Hz)', -1, 6)
               CALL PRNTLO (6, SPINJJ, LABEL, NNUCL, NNUCL, 100, 100)
      END IF
C
 1000 FORMAT (2X, A, D16.9)
 1200 FORMAT (2X, A, 1X, A6, 6X, A3, D16.9)
 1300 FORMAT (2X, A, 1X, A6, 3X, A3, D16.9)
 1400 FORMAT (2X, A, 1X, A6, 4X, A3, D16.9)
 1600 FORMAT (2X, A, 1X, A6, 9X, A3, D16.9)
 1500 FORMAT (2X, A)
C
      RETURN
      END
