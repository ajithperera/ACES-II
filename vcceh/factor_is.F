      SUBROUTINE FACTOR_IS(XMATX, N, IFERMI, ISDIP, IPSO, IDSO, 
     &                     IEOMPROP,ICALL,IUNITIS)
C
C Calculates the factors for NMR spin-spin coupling constant.
C The SI system of units has been used. Here there is a very 
C reasonable hard coded limit on the number of atoms (100) in 
C the calcualtion. Coded by Ajith 10/93.               
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      DOUBLE PRECISION NINE, JJKEEP
      INTEGER ATMCHRG,ZMTORDATMCHRG,ZMTORDISOTOPES
C
      DIMENSION CENTER(100, 3)
      DIMENSION ZMTORDCENTER(100, 3)
      DIMENSION ATMCHRG(100)
      DIMENSION ZMTORDATMCHRG(100)
      DIMENSION ATMMASS(100)
      DIMENSION ZMTORDATMMASS(100)
      DIMENSION ISOTOPES(100)
      DIMENSION ZMTORDISOTOPES(100)

      CHARACTER*6 ANAME_R(100) 
      CHARACTER*6 LABEL(100)
      CHARACTER*6 ANAME_S(100)
      CHARACTER*6 ZMTORDLABEL(100)
      CHARACTER*5 VMLORDLABEL(100)
      CHARACTER*5 TMP1(100)
      CHARACTER*5 TMP2(100)
      CHARACTER*5 TMP3(100)
      CHARACTER*6 NAME(100)

      CHARACTER*50 STRINGFC
      CHARACTER*40 STRINGPSO
      CHARACTER*48 STRINGSD
C
      DIMENSION XMATX(N,N), VMLORDSPINJJ(100,100),
     &          ZMTORDSPINJJ(100,100), JJKEEP(100,100)
      DIMENSION GTAB_S(100),GTAB_R(100)
      DIMENSION ITYPE(100), MAPCZ(100)
      DIMENSION IMAP(100)
      DIMENSION ICHARGE(100)
C
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FILES/LUOUT,MOINTS
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA IONE,IZERO/1,0/
      DATA ONE,TWO,THREE,FOUR,EIGHT,NINE/1.D0,2.D0,3.D0,4.D0,8.D0,9.D0/
C  
C  Nuclei g-factors (per second per tesla). Extracted from CRC handbook
C  of chemistry and physics, 71'st edition 1990-1991
C  Update whenever possible.

C                      1H                3He              7Li        
      DATA GTAB_R/ 2.675202874D+08, -2.037995287D+08,  1.039754808D+08,
C                      9BE               11B              13C
     &          -3.759981028D+07,  8.584481142D+07,  6.728214013D+07,
C                      15N               17O              19F       
     &          -2.712607916D+07, -3.628049629D+07,  2.518130433D+08,
C                      21Ne              23Na             25Mg
     &          -2.113073577D+07,  7.080361014D+07,  1.638829367D+07,
C                      27Al              29Si             31P      
     &           6.976227164D+07, -5.319083216D+07,  1.083932031D+08,
C                      33S               35Cl             39Ar  
     &           2.055664899D+07,  2.624164069D+07, -2.175763532D+07,
C                      39K               43Ca             45Sc
     &           1.249899943D+07, -1.802585043D+07,  6.508751029D+07,
C                      47Ti              51V              53Cr   
     &          -1.510531509D+07,  7.045456939D+07, -1.515167627D+07,
C                      55Mb              57Fe             59Co   
     &           6.645166199D+07,  8.678352950D+06,  6.335662908D+07,
C                      61Ni              63Cu             67Zn   
     &          -2.394752862D+07,  7.098816083D+07,  1.677240178D+07,
C                      69Ga              73Ge             75As     
     &           6.438807865D+07, -9.360222020D+06,  4.596110641D+07,  
C                      77Se              79Br             83Kr 
     &           5.125209194D+07,  6.725563891D+07, -1.033120987D+07,
C 
     &            0.0           ,  0.0            ,  0.0            ,
C 
     &            0.0           ,  0.0            ,  0.0            ,
C 
     &            0.0           ,  0.0            ,  0.0            ,
C 
     &            0.0           ,  0.0            ,  0.0            ,
C 
     &            0.0           ,  0.0            /
C                             ----SPECIAL NUCLEI-----
C                                
      DATA GTAB_S/ 4.1066287D+07,  2.853498382D+07, -4.470122142D+07,
C                     10B             14N                24Na
     &           2.874679354D+07,  1.933779354D+07,  2.023528505D+07,
C                     26Al            32P                35S
     &           2.685904819D+07, -1.208848745D+07,  3.192944387D+07,
C                     36Cl            37Cl               37Ar
     &           3.078325666D+07,  2.184357114D+07,  3.671886045D+07, 
C                     40K             41K                42K
     &          -1.554285415D+07,  6.860679605D+06, -2.735954222D+07, 
C                     43K             41Ca               45Ca
     &           5.204499351D+06, -2.182331875D+07, -1.815873086D+07,
C                     47Ca            43Sc               44Sc
     &          -1.923977061D+07,  6.198873460D+07,  5.986770726D+07,
C                     46Sc            47Sc               48Sc
     &           3.639956601D+07,  7.307281297D+07,  2.970000000D+07,
C                     45Ti            49Ti               48V
     &           1.299984505D+06, -1.510951458D+07,  2.406681832D+07,
C                     49V             50V                51Cr
     &           6.116768176D+07,  2.670650527D+07, -1.278090024D+07,
C                     52Mn            53Mn               59Fe
     &           2.449971646D+07,  6.889917852D+07, -1.072829314D+07,
C                     55Co            56Co               57Co
     &           6.598447644D+07,  4.609813463D+07,  4.609813459D+07,
C                     60Co            57Ni               63Ni
     &           3.600000000D+07, -2.547969621D+07,  4.789416581D+07,
C                     61Cu            65Cu               67Cu
     &           6.733919713D+07,  7.604635647D+07,  8.027062189D+07,  
C                     65Zn            69Zn               71Zn
     &           1.473224540D+07,  1.230000000D+07,  1.120000000D+07,
C                     67Ga            71Ga               72Ga
     &           5.765499680D+07,  8.181185615D+07, -2.111748289D+06,
C                     69Ge            71Ge               71As
     &           1.408088474D+07,  5.239621739D+07,  5.343392432D+07,
C                     72As            76As               77As
     &          -5.644278992D+07, -2.162421586D+07,  1.134862982D+07,
C                     73Se            75Se               79Se
     &           9.153107244D+06,  8.237796519D+07, -1.395772832D+07,
C                     76Br            77Br               81Br
     &           2.625606064D+07,  3.106734888D+07,  7.249778194D+07,
C                     82Br            79Kr               81Kr
     &           1.560000000D+07,  5.134254575D+07,  -1.242511502D+07/
C
      DATA ANAME_R/ '[1] H ', '[3]He ', '[7]Li ', '[9]Be ', '[11]B ',
     &              '[13]C ', '[15]N ', '[17]O ', '[19]F ', '[21]Ne',
     &              '[23]Na', '[25]Mg', '[27]Al', '[29]Si', '[31]P ',
     &              '[33]S ', '[35]Cl', '[39]Ar', '[39]K ', '[43]Ca',
     &              '[45]Sc', '[47]Ti', '[51]V ', '[53]Cr', '[55]Mn',
     &              '[57]Fe', '[59]Co', '[61]Ni', '[63]Cu', '[67]Zn',
     &              '[69]Ga', '[73]Ge', '[75]As', '[77]Se', '[79]Br',
     &              '[83]Kr', '[14]N ', '[37]Cl', '[49]Ti', '[2] H ',
     &              '[6]Li ', '[10]B ', '[65]Cu', '[71]Ga', '[81]Br'/ 

      DATA ANAME_S/ '2 [H] ', '[3] H ', '[7]Be ', '[10]BA', '[14]N ',
     &              '[24]Na', '[26]Al', '[32]P ', '[35]S ', '[36]Cl',
     &              '[41]K ', '[43]K ', '[41]Ca', '[45]Ca', '[47]Ca',
     &              '[43]Sc', '[44]Sc', '[45]Ti', '[49]Ti', '[48]V ',
     &              '[49]V ', '[51]Cr', '[52]Mn', '[53]Mn', '[59]Fe',
     &              '[59]Fe', '[55]Co', '[56]Co', '[57]Co', '[60]Co',
     &              '[63]Ni', '[61]Cu', '[65]Cu', '[67]Cu', '[65]Zn',
     &              '[69]Zn', '[71]Zn', '[67]Ga', '[71]Ga', '[72]Ga',
     &              '[69]Ge', '[71]Ge', '[71]As', '[72]As', '[76]As',
     &              '[77]As', '[73]Se', '[75]Se', '[79]Se', '[76]Br',
     &              '[77]Br', '[81]Br', '[82]Br', '[79]Kr', '[81]Kr'/


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
      CALL GETREC(20, 'JOBARC', 'ATOMMASS', NCENTR*IINTFP,ATMMASS)
      CALL GETREC(20, 'JOBARC', 'ATOMCHRG', NCENTR,ATMCHRG)

      IF (NCENTER .GT. 100) THEN
         Write(6,"(A,I3,A)") "The number of centers", NCENTER, 
     &                       " exceeded the currently allowed", 
     &                       " maximum of 100 atoms"
         Call ERREX
      ENDIF 

      K = 1
      DO 2 I = 1, NCENTR
         IF (ITYPE(I) .ne. 0) THEN
            ICHARGE(K) = ITYPE(I)
            K = K + 1
         ENDIF
 2    CONTINUE

C Change the ordering to ZMAT. This helps in managing the vibrational 
C corrections codes. 

      IF (NCENTR .NE. NATOMS) THEN
         CALL GETREC(20,'JOBARC','MAP2ZMAT',NCENTR,IMAP)
         NCOUNT = 0
         DO IATMVML = 1, NCENTR
            IATMZMT = IMAP(IATMVML)
            IF (IATMZMT .NE. IZERO) THEN
               NCOUNT = NCOUNT + 1
               ZMTORDATMCHRG(NCOUNT) = ATMCHRG(IATMZMT)
               ZMTORDATMMASS(NCOUNT) = ATMMASS(IATMZMT)
            ENDIF
         ENDDO

         DO IATOM = 1, NATOMS
            ATMCHRG(IATOM) =  ZMTORDATMCHRG(IATOM)
            ATMMASS(IATOM) =  ZMTORDATMMASS(IATOM)
         ENDDO 
        
         CALL FIXDUMMY(VMLORDLABEL,TMP1,TMP2,TMP3,IMAP,ICHARGE,
     &                 NCENTR,NATOMS)
         DO IATMVML = 1, NATOMS
            IATMZMT = IMAP(IATMVML)
            ZMTORDLABEL(IATMZMT)   = VMLORDLABEL(IATMVML)
            ZMTORDATMCHRG(IATMZMT) = ATMCHRG(IATMVML)
            ZMTORDATMMASS(IATMZMT) = ATMMASS(IATMVML)
         ENDDO
      ELSE
         CALL GETREC(20,'JOBARC','MAP2ZMAT',NCENTR,IMAP)
         CALL GETCREC(20,'JOBARC','ZSYM',5*NATOMS,LABEL)
         DO IATMZMT = 1, NATOMS
            ZMTORDLABEL(IATMZMT)   = LABEL(IATMZMT)
            ZMTORDATMCHRG(IATMZMT) = ATMCHRG(IATMZMT)
            ZMTORDATMMASS(IATMZMT) = ATMMASS(IATMZMT)
         ENDDO
      ENDIF 
#ifdef _DEBUG_LVL0
      write(6,"(a)") " Input ordered atom specific arrays"
      write(6,"(6(4x,6a))")    (ZMTORDLABEL(i),i=1,Natoms)
      write(6,"(6(1x,i4))")    (ZMTORDATMCHRG(i),i=1,Natoms)
      write(6,"(6(1x,F10.5))") (ZMTORDATMMASS(i),i=1,Natoms)  
#endif 

      CALL GET_ISODATA(ZMTORDATMMASS,ZMTORDLABEL,ZMTORDATMCHRG,
     &                 ZMTORDISOTOPES,NAME,IUNITIS,NATOMS)

C Reorder back to internal order since the NMR tensors are processed in 
C that order.

      DO IATMVML = 1, NCENTR 
         IATMZMT  = IMAP(IATMVML)
         ATMCHRG(IATMVML)  = ZMTORDATMCHRG(IATMZMT) 
         ATMMASS(IATMVML)  = ZMTORDATMMASS(IATMZMT)  
         ISOTOPES(IATMVML) = ZMTORDISOTOPES(IATMZMT)
         LABEL(IATMVML)    = ZMTORDLABEL(IATMZMT)
      ENDDO
#ifdef _DEBUG_LVL0
      write(6,"(a)") " Internally ordered atom specific arrays"
      write(6,"(6(4x,i4))")    (ISOTOPES(i),i=1,Natoms)
      write(6,"(6(1x,i4))")    (ATMCHRG(i),i=1,Natoms)
      write(6,"(6(1x,F10.5))") (ATMMASS(i),i=1,Natoms)   
      write(6,"(6(1x,5a))")    (LABEL(i),i=1,Natoms)  
      write(6,*)
#endif 

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
         IF (ISOTOPES(I) .EQ. 1) THEN
            LABEL(I) = LABEL(I)
         ELSE
            LABEL(I) = ANAME_R(ICHARGE(I))
         ENDIF 
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
            IF (ISOTOPES(I) .EQ. 1) THEN
               LABEL(I) = LABEL(I)
               GFCI = GTAB_S(ISOTOPES(I))
            ELSE
               LABEL(I) = ANAME_R(ICHARGE(I))
               GFCI = GTAB_R(ICHARGE(I))
            ENDIF 
            IF (IFLAGS(1) .GE. 20) THEN
            WRITE(6, 1200)'Nuclear g-factor for nuclei', LABEL(I), 
     &                    ' = ', GFCI
            ENDIF
            DO 30 J = 1, I
               IF (ISOTOPES(J) .EQ. 1) THEN
                  GFCJ = GTAB_S(ISOTOPES(J))
               ELSE
                  GFCJ = GTAB_R(ICHARGE(J))
               ENDIF 
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
            IF (ISOTOPES(I) .EQ. 1) THEN
               LABEL(I) = LABEL(I)
               GFCI = GTAB_S(ISOTOPES(I))
            ELSE
               LABEL(I) = ANAME_R(ICHARGE(I))
               GFCI = GTAB_R(ICHARGE(I))
            ENDIF 
            IF (IFLAGS(1) .GE. 20) THEN
            WRITE(6, 1300)'Nuclear g-factor for nuclei', LABEL(I), 
     &                    ' = ', GFCI
            ENDIF
            DO 50 J = 1, I
               VMLORDSPINJJ(I, J) = 0.0D+00
               IF (ISOTOPES(J) .EQ. 1) THEN
                  GFCJ = GTAB_S(ISOTOPES(J))
               ELSE
                  GFCJ = GTAB_R(ICHARGE(J))
               ENDIF
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
            IF (ISOTOPES(I) .EQ. 1) THEN
               LABEL(I) = LABEL(I)
               GFCI = GTAB_S(ISOTOPES(I))
            ELSE
               LABEL(I) = ANAME_R(ICHARGE(I))
               GFCI = GTAB_R(ICHARGE(I))
            ENDIF 
            IF (IFLAGS(1) .GE. 20) THEN
            WRITE(6, 1400)'Nuclear g-factor for nuclei', LABEL(I), 
     &                    ' = ', GFCI
            ENDIF
            DO 80 J = 1, I
               VMLORDSPINJJ(I, J) = 0.0D+00
               IF (ISOTOPES(J) .EQ. 1) THEN
                  GFCJ = GTAB_S(ISOTOPES(J))
               ELSE
                  GFCJ = GTAB_R(ICHARGE(J))
               ENDIF
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
            IF (ISOTOPES(I) .EQ. 1) THEN
               LABEL(I) = LABEL(I)
               GFCI = GTAB_S(ISOTOPES(I))
            ELSE
               LABEL(I) = ANAME_R(ICHARGE(I))
               GFCI = GTAB_R(ICHARGE(I))
            ENDIF 
            IF (IFLAGS(1) .GE. 20) THEN
            WRITE(6, 1600)'Nuclear g-factor for nuclei', LABEL(I), 
     &                    ' = ', GFCI
            ENDIF
            DO 120 J = 1, I
               VMLORDSPINJJ(I, J) = 0.0D+00
               IF (ISOTOPES(J) .EQ. 1) THEN
                  GFCJ = GTAB_S(ISOTOPES(J))
               ELSE
                  GFCJ = GTAB_R(ICHARGE(J))
               ENDIF
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
      CALL REORD_2ZMAT(JJKEEP,CENTER,IMAP,LABEL,ZMTORDSPINJJ,
     &                 NAME,ZMTORDCENTER,NCENTR,100)

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
