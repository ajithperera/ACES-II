C
      SUBROUTINE FACTOR(XMATX, N, IFERMI, ISDIP, IPSO, IDSO, IEOMPROP,
     &                  ICALL)
C
C Calculates the factors for NMR spin-spin coupling constant.
C The SI system of units has been used. Here there is a very 
C reasonable hard coded limit on the number of atoms (100) in 
C the calcualtion. Coded by Ajith 10/93.               
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      DOUBLE PRECISION NINE, JJKEEP
C
      DIMENSION CENTER(100, 3)
      DIMENSION ZMTORDCENTER(100, 3)
      CHARACTER*6 ANAME(100), LABEL(100)
      CHARACTER*6 ZMTORDLABEL(100)
      CHARACTER*6 TMP1(100)
      CHARACTER*6 TMP2(100)
      CHARACTER*6 TMP3(100)

      CHARACTER*50 STRINGFC
      CHARACTER*40 STRINGPSO
      CHARACTER*48 STRINGSD
C
      DIMENSION XMATX(N,N), VMLORDSPINJJ(100,100),
     &          ZMTORDSPINJJ(100,100), JJKEEP(100,100)
      DIMENSION GTAB(45), ICHARGE(100), ITYPE(100), MAPCZ(100)
      DIMENSION COSYM(4,4)
      DIMENSION IMAP(100)
C
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FILES/LUOUT,MOINTS
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA IONE/1/
      DATA ONE,TWO,THREE,FOUR,EIGHT,NINE/1.D0,2.D0,3.D0,4.D0,8.D0,9.D0/
      DATA COSYM / 1.D0  , 1.D0  ,  1.D0 ,  1.D0
     &           , 1.D0  ,-1.D0  ,  1.D0 , -1.D0
     &           , 1.D0  , 1.D0  , -1.D0 , -1.D0
     &           , 1.D0  ,-1.D0  , -1.D0 ,  1.D0 /
C  
C  Nuclei g-factors (per second per tesla). Extracted from CRC handbook
C  of chemistry and physics, 71'st edition 1990-1991
C  Update whenever possible.

C                      1H                3He              7Li        
      DATA GTAB/ 2.675202874D+08, -2.037995287D+08,  1.039754808D+08,
C                      9BE               11B              13C
     &          -3.759981028D+07,  8.584481142D+07,  6.728214013D+07,
C                      15N               17O              19F       
     &          -2.712607916D+07, -3.628049629D+07,  2.518130433D+08,
C                      21Ne              23Na             25Mg
     &          -2.113073577D+07,  7.080361014D+07,  1.638829367D+07,
C                      27Al              29Si             31P      
     &           6.976227164D+07, -5.319083216D+07,  1.083932031D+08,
C                      33S               35Cl             39Ar  
     &           2.055664899D+07,  2.624164069D+07,  1.778911832D+07,
C                      39K               43Ca             45Sc   
     &           1.249899943D+07, -1.802585043D+07,  6.508751029D+07,
C                      47Ti              51V              53Cr   
     &          -1.510531509D+07,  7.045456939D+07, -1.515167627D+07,
C                      55Mn              57Fe             59Co   
     &           6.645166199D+07,  8.678352950D+06,  6.335662908D+07,
C                      61Ni              63Cu             67Zn   
     &          -2.394752862D+07,  7.098816083D+07,  1.677240178D+07,
C                      69Ga              73Ge             75As     
     &           6.438807865D+07, -9.360222020D+06,  4.596110641D+07,  
C                      77Se              79Br             83Kr 
     &           5.125209194D+07,  6.725563891D+07, -1.033120987D+07,
C                                  SPECIAL NUCLEI
C                      14N               37Cl             49Ti  
     &           1.933759265D+07,  2.184339522D+07, -1.510939290D+07,
C                      2H                6Li              10B    
     &           4.106604279D+07,  3.937108191D+07,  2.587126212D+08,
C                      65Cu              71Ga             81Br 
     &           7.604574401D+07,  8.181119726D+07,  7.249719806D+07/
C
      DATA ANAME/ '[1] H ', '[3]HE ', '[7]LI ', '[9]Be ', '[11]B ',
     &            '[13]C ', '[15]N ', '[17]O ', '[19]F ', '[21]NE',
     &            '[23]Na', '[25]Mg', '[27]Al', '[29]Si', '[31]P ',
     &            '[33]S ', '[35]Cl', '[39]Ar', '[39]K ', '[43]Ca',
     &            '[45]Sc', '[47]Ti', '[51]V ', '[53]Cr', '[55]Mn',
     &            '[57]Fe', '[59]Co', '[61]Ni', '[63]Cu', '[67]Zn',
     &            '[69]Ga', '[73]Ge', '[75]As', '[77]Se', '[79]Br',
     &            '[83]Kr', '[14]N ', '[37]Cl', '[49]Ti', '[2] H ',
     &            '[6]LI ', '[10]B ', '[65]Cu', '[71]Ga', '[81]Br'/ 
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
      CALL GETREC(20, 'JOBARC', 'VMORDCOR', 100*3*IINTFP, CENTER)

      IF (NCENTER .GT. 100) THEN
         Write(6,"(A,I3,A)") "The number of centers", NCENTER, 
     &                       " exceeded the currently allowed", 
     &                       " maximum of 100 atoms"
         Call ERREX
      ENDIF 
C   
      K = 1
      DO 2 I = 1, NCENTR
         IF (ITYPE(I) .ne. 0) THEN 
            ICHARGE(K) = ITYPE(I)
            K = K + 1
         ENDIF
 2    CONTINUE
C
C Write the ordering of atoms printed in the spin-spin coupling tensor.
C
      WRITE(LUOUT, 5000) 
      WRITE(LUOUT, *)
      WRITE(LUOUT, 5005)
      WRITE (LUOUT, 5010)
      WRITE (LUOUT, 5015)
      WRITE (LUOUT, 5005)
C
      DO 3 I = 1, NATOMS
         LABEL(I) = ANAME(ICHARGE(I))
         WRITE(LUOUT, 5020) LABEL(I), (CENTER(I, K), K = 1, 3)
 3    CONTINUE
C
      WRITE(LUOUT, 5005)
            

 5000 FORMAT (1X, 'The NMR spin-spin coupling tensors correspond to',
     &      /, 1X,'the following orientation and atomic ordering.')

 5005 FORMAT (T4, 56('-'))      
 5010 FORMAT (T23, 'Cartesian coordinates')
 5015 FORMAT (T5, '       ', T22, 'X', T34, 'Y', T46, 'Z')
 5020 FORMAT (T5, A, T17, F10.7, T29, F10.7, T41, F10.7)
C
C Set the output strings depending on the method.
C
      IF (IFERMI .NE. 0) THEN
C
         IF (IEOMPROP .EQ. 0)
     &   STRINGFC='Fermi-contact contribution to J (CI-like, in Hz)'
C
         IF (IEOMPROP .EQ. 1)
     &   STRINGFC='Fermi-contact contribution to J (Linear, in Hz)'
C
         IF (IEOMPROP .EQ. 2)
     &   STRINGFC='Fermi-contact contribution to J (Quadratic, in Hz)'

      ELSE IF (IPSO .NE. 0) THEN
C
         IF (IEOMPROP .EQ. 0)
     &   STRINGPSO='PSO contribution to J (CI-LIKE, in Hz)'
C
         IF (IEOMPROP .EQ. 1)
     &   STRINGPSO='PSO contribution to J (Linear, in Hz)'
C
         IF (IEOMPROP .EQ. 2)
     &   STRINGPSO='PSO contribution to J (Quadratic, in Hz)'

      ELSE IF (ISDIP .NE. 0) THEN
C
         IF (IEOMPROP .EQ. 0)
     &   STRINGSD='Spin-dipole contribution to J (CI-Like, in Hz)'
C
         IF (IEOMPROP .EQ. 1)
     &   STRINGSD='Spin-dipole contribution to J (Linear, in Hz)'
C 
        IF (IEOMPROP .EQ. 2)
     &   STRINGSD='Spin-dipole contribution to J (Quadratic, in Hz)'

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
         CALL PUTREC(20,"JOBARC","JCC_FC  ",N*N*IINTFP,XMATX)
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
               VMLORDSPINJJ(I, J) = XMATX(I, J)*GFCI*GFCJ*FACTFC
 30         CONTINUE
 20      CONTINUE
C
         WRITE(6,*)
         CALL HEADER(STRINGFC, -1, 6)
         CALL PRNTLO (6, VMLORDSPINJJ, LABEL, NNUCL, NNUCL, 100, 100)
C
      END IF
C
      IF(IPSO .NE. 0) THEN 
         NCOMPO = 3
         IF (IFLAGS(1) .GE. 20) THEN
         WRITE(6, 1000) 'Conversion factor for spin-orbit term = ', 
     &                   FACTPSO
         ENDIF
         CALL PUTREC(20,"JOBARC","JCC_PSO ",N*N*IINTFP,XMATX)
         NNUCL = N/NCOMPO
         DO 40 I = 1, NNUCL
            LABEL(I) = ANAME(ICHARGE(I))
            GFCI = GTAB(ICHARGE(I))
            IF (IFLAGS(1) .GE. 20) THEN
            WRITE(6, 1300)'Nuclear g-factor for nuclei', LABEL(I), 
     &                    ' = ', GFCI
            ENDIF
            DO 50 J = 1, I
               VMLORDSPINJJ(I, J) = 0.0D+00
               GFCJ = GTAB(ICHARGE(J))
               DO 60 KK = 1, NCOMPO
                  II = (I - 1)*NCOMPO + KK
                  JJ = (J - 1)*NCOMPO + KK
                  VMLORDSPINJJ(I, J) = VMLORDSPINJJ(I, J) + 
     &                            XMATX(II, JJ)*GFCI*GFCJ*FACTPSO/THREE 
 60            CONTINUE
 50         CONTINUE
 40      CONTINUE
         WRITE(6, *)
         CALL HEADER(STRINGPSO, -1, 6)
         CALL PRNTLO (6, VMLORDSPINJJ, LABEL, NNUCL, NNUCL, 100, 100)
C
      END IF
C
      IF(ISDIP .NE. 0) THEN
         NCOMPO = 6
         IF (IFLAGS(1) .GE. 20) THEN
         WRITE(6,1000) 'Conversion factor for spin-dipole term = ', 
     &               FACTSD
         ENDIF
         CALL PUTREC(20,"JOBARC","JCC_SDIP",N*N*IINTFP,XMATX)
         NNUCL = N/NCOMPO
         DO 70 I = 1, NNUCL
            LABEL(I) = ANAME(ICHARGE(I))
            GFCI = GTAB(ICHARGE(I))
            IF (IFLAGS(1) .GE. 20) THEN
            WRITE(6, 1400)'Nuclear g-factor for nuclei', LABEL(I), 
     &                    ' = ', GFCI
            ENDIF
            DO 80 J = 1, I
               VMLORDSPINJJ(I, J) = 0.0D+00
               GFCJ = GTAB(ICHARGE(J))
               DO 90 KK = 1, NCOMPO
                  II = (I - 1)*NCOMPO + KK
                  JJ = (J - 1)*NCOMPO + KK
                  VMLORDSPINJJ(I, J) = VMLORDSPINJJ(I, J) + 
     &                            XMATX(II, JJ)*GFCI*GFCJ*FACTSD/THREE 
C
                  IF (KK .EQ. 2 .OR. KK .EQ. 3 .OR. KK .EQ. 5) THEN
                     VMLORDSPINJJ(I, J) = VMLORDSPINJJ(I, J) +
     &                            XMATX(II, JJ)*GFCI*GFCJ*FACTSD/THREE 
                  ENDIF
 90            CONTINUE
 80         CONTINUE
 70      CONTINUE 
         WRITE(6, *) 
         CALL HEADER(STRINGSD, -1, 6)
         CALL PRNTLO (6, VMLORDSPINJJ, LABEL, NNUCL, NNUCL, 100, 100)
C
      END IF
C
      IF(IDSO .NE. 0) THEN 
         NCOMPO = 3
         IF (IFLAGS(1) .GE. 20) THEN
         WRITE(6, 1000)'Conversion factor for diam. spin-orbit term = '
     &                   , FACTDSO
         ENDIF
         CALL PUTREC(20,"JOBARC","JCC_DSO ",N*N*IINTFP,XMATX)
         NNUCL = N/NCOMPO
         DO 100 I = 1, NNUCL
            LABEL(I) = ANAME(ICHARGE(I))
            GFCI = GTAB(ICHARGE(I))
            IF (IFLAGS(1) .GE. 20) THEN
            WRITE(6, 1600)'Nuclear g-factor for nuclei', LABEL(I), 
     &                    ' = ', GFCI
            ENDIF
            DO 120 J = 1, I
               VMLORDSPINJJ(I, J) = 0.0D+00
               GFCJ = GTAB(ICHARGE(J))
               DO 130 KK = 1, NCOMPO
                  II = (I - 1)*NCOMPO + KK
                  JJ = (J - 1)*NCOMPO + KK
                  VMLORDSPINJJ(I, J) = VMLORDSPINJJ(I, J) +
     &                           XMATX(II, JJ)*GFCI*GFCJ*FACTDSO/THREE 
 130           CONTINUE
 120        CONTINUE
 100     CONTINUE
         WRITE(6, *)
         CALL HEADER('DSO contribution to J (in Hz)', -1, 6)
         CALL PRNTLO (6, VMLORDSPINJJ, LABEL, NNUCL, NNUCL, 100, 100)
C
      END IF

      NNUCL2 = 100*100
      IF (ICALL .EQ. 0) THEN
         CALL DCOPY(NNUCL2, VMLORDSPINJJ, 1, JJKEEP, 1)
      ELSE IF (ICALL .GT. 0) THEN
         CALL DAXPY(NNUCL2, ONE, VMLORDSPINJJ, 1, JJKEEP, 1) 
      ENDIF

C Change the ordering to ZMAT. This helps in managing the 
C vibrational corrections codes. 

      IF (ICALL .EQ. 3) THEN
         IF (NCENTR .NE. NATOMS) THEN
            CALL GETREC(20,'JOBARC','MAP2ZMAT',NCENTR,IMAP)
            CALL GETCREC(20,'JOBARC','ZSYM',5*NCENTR,ZMTORDLABEL)
            CALL FIXDUMMY(ZMTORDLABEL,TMP1,TMP2,TMP3,IMAP,ICHARGE,
     &                    NCENTR,NATOMS)
         ELSE
            CALL GETREC(20,'JOBARC','MAP2ZMAT',NCENTR,IMAP)
         ENDIF 
        
      CALL REORD_2ZMAT(JJKEEP,CENTER,IMAP,LABEL,ZMTORDSPINJJ,
     &                 ZMTORDLABEL,ZMTORDCENTER,NCENTR,100)

      Write(6,*)
      Write(6,"(2x,2a)") 'The total NMR spin-spin coupling tensors',
     &                   ' correspond to the'
      Write(6,"(2x,a)")  'following orientation and atomic ordering.'
      WRITE(LUOUT, *)
      WRITE(LUOUT, 5005)
      WRITE (LUOUT, 5010)
      WRITE (LUOUT, 5015)
      WRITE (LUOUT, 5005)
C
      DO I = 1, NATOMS
         WRITE(LUOUT, 5020) ZMTORDLABEL(I), 
     &                      (ZMTORDCENTER(I, K), K = 1, 3)
      ENDDO 
C
      WRITE(LUOUT, 5005)

C Note that IFLAGS(18) get reset from 13 (JSC_ALL) to 10 in vcceh. 
C Therefore JSC_ALL get set in setmet.f

      WRITE(6,*)
      CALL HEADER('The total spin-spin coupling tensor (in Hz)', 
     &             -1, 6)
      CALL PRNTLO_EXPANDED(6, ZMTORDSPINJJ, ZMTORDLABEL, NNUCL, 
     &                     NNUCL, 100, 100)
      CALL PUTREC(20, "JOBARC", "JJTENSOR",NNUCL2*IINTFP, ZMTORDSPINJJ) 
      ENDIF 
C
 1000 FORMAT (2X, A, D16.9)
 1200 FORMAT (2X, A, 1X, A6, 6X, A3, D16.9)
 1300 FORMAT (2X, A, 1X, A6, 3X, A3, D16.9)
 1400 FORMAT (2X, A, 1X, A6, 4X, A3, D16.9)
 1600 FORMAT (2X, A, 1X, A6, 9X, A3, D16.9)
 1500 FORMAT (2X, A)
C
      WRITE(6, *)
      WRITE(6, 1500)'@-PRFACT-I successful completion.'
C
      RETURN
      END
