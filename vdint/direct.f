      SUBROUTINE DIRECT(BIGVEC,ICENTA,ICENTB,ICENTC,ICENTD,
     *                  NCENTA,NCENTB,NCENTC,NCENTD,
     *                  ISOA,ISOB,ISOC,ISOD,
     *                  SIGNAX,SIGNAY,SIGNAZ,SIGNBX,SIGNBY,SIGNBZ,
     *                  SIGNCX,SIGNCY,SIGNCZ,SIGNDX,SIGNDY,SIGNDZ,
     *                  MCCINT,MAXDER,DERECC,DERECD,EXPECT,
     *                  IATOM,MULATM,IPRINT)
C
C     This subroutine determines type of two-electron integral
C     and sets directives for integral calculation.
C
C     tuh
C
      IMPLICIT DOUBLE PRECISION (F,S), LOGICAL(D)
      DOUBLE PRECISION ONE, TWO
      PARAMETER (ONE = 1.00 D00, TWO = 2.00 D00)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      LOGICAL ONECEN, DERECC, DERECD, EXPECT, BIGVEC
      CHARACTER*4 ATYPE
      LOGICAL CROSS1, CROSS2
      COMMON /CROSSD/ CROSS1, CROSS2
      LOGICAL         DC101, DC1H1, DC1E1, DC2H1, DC2E1,
     *                DC102, DC1H2, DC1E2, DC2H2, DC2E2,
     *                DPATH1, DPATH2
      COMMON /SUBDIR/ DC101, DC1H1, DC1E1, DC2H1, DC2E1,
     *                DC102, DC1H2, DC1E2, DC2H2, DC2E2,
     *                DPATH1, DPATH2, NTOTAL
      COMMON /CRSDIR/ DHCHX, DHCHY,  DHCHZ,
     *                DHCEX, DHCEX1, DHCEX2,
     *                DHCEY, DHCEY1, DHCEY2,
     *                DHCEZ, DHCEZ1, DHCEZ2
      DOUBLE PRECISION
     *        SIGN1X, SIGN1Y, SIGN1Z, SIGN2X, SIGN2Y, SIGN2Z,
     *        SIGN3X, SIGN3Y, SIGN3Z, SIGN4X, SIGN4Y, SIGN4Z
      LOGICAL TWOCEN, THRCEN, FOUCEN, DERONE, DERTWO
      COMMON /EXPCOM/ SIGN1X, SIGN1Y, SIGN1Z, SIGN2X, SIGN2Y, SIGN2Z,
     *                SIGN3X, SIGN3Y, SIGN3Z, SIGN4X, SIGN4Y, SIGN4Z,
     *                NCENT1, NCENT2, NCENT3, NCENT4,
     *                ISO1,   ISO2,   ISO3,   ISO4,
     *                DERONE, DERTWO, TWOCEN, THRCEN, FOUCEN,
     *                NINTYP, NCCINT
      LOGICAL TRANS
      COMMON /DIRPRT/ SIGNDR(4,3), TRANS, NATOMS, IATOMS(4), ISOPDR(4)
      LOGICAL         DZER
      COMMON /DERZER/ FZERO, DZER, IZERO
      COMMON /DEODIR/ DEOVEC(54)
      COMMON /DETDIR/ DETVEC(36)
      COMMON /DHODIR/ DHOVEC(18)
      COMMON /DEOADR/ IEOVEC(54)
      COMMON /DETADR/ IETVEC(36)
      COMMON /DHOADR/ IHOVEC(18)
      COMMON /DEOFAC/ FEOVEC(54)
      COMMON /DETFAC/ FETVEC(36)
      COMMON /DHOFAC/ FHOVEC(18)
C
C     **********************************************
C     *************** Initialization ***************
C     **********************************************
C
      NCCINT = MCCINT
C
C     Number of centers
C
      ONECEN = .FALSE.
      TWOCEN = .FALSE.
      THRCEN = .FALSE.
      FOUCEN = .FALSE.
C
C     Expansion coefficient directives
C
      DERECC = .FALSE.
      DERECD = .FALSE.
C
C     Subroutine directives
C
      DZER  = .FALSE.
      DC101 = .FALSE.
      DC1H1 = .FALSE.
      DC1E1 = .FALSE.
      DC2H1 = .FALSE.
      DC2E1 = .FALSE.
      DC102 = .FALSE.
      DC1H2 = .FALSE.
      DC1E2 = .FALSE.
      DC2H2 = .FALSE.
      DC2E2 = .FALSE.
C
C     Directives for cross differentiation
C
      CROSS1 = .FALSE.
      CROSS2 = .FALSE.
      DHCHX  = .FALSE.
      DHCHY  = .FALSE.
      DHCHZ  = .FALSE.
      DHCEX  = .FALSE.
      DHCEY  = .FALSE.
      DHCEZ  = .FALSE.
      DHCEX1 = .FALSE.
      DHCEY1 = .FALSE.
      DHCEZ1 = .FALSE.
      DHCEX2 = .FALSE.
      DHCEY2 = .FALSE.
      DHCEZ2 = .FALSE.
C
C     Calculation paths
C
      DPATH1 = .FALSE.
      DPATH2 = .FALSE.
C
C     Derivative types DEOVEC, DETVEC and DHOVEC
C
      DO 100 I = 1, 54
         DEOVEC(I) = .FALSE.
  100 CONTINUE
      DO 110 I = 1, 36
         DETVEC(I) = .FALSE.
  110 CONTINUE
      DO 120 I = 1, 18
         DHOVEC(I) = .FALSE.
  120 CONTINUE
C
      IF (MAXDER .EQ. 0) THEN
         DPATH1 = .TRUE.
         DZER   = .TRUE.
         FZERO  = ONE
         IZERO  = 0
         DC101  = .TRUE.
         DC2H1  = .TRUE.
         NINTYP = 1
         NTOTAL = NCCINT
         RETURN
      ELSE IF (MAXDER .EQ. 1) THEN
         DERONE = .TRUE.
         DERTWO = .FALSE.
      ELSE
         DERONE = .TRUE.
         DERTWO = .TRUE.
      END IF
C
C     If only derivatives with respect to IATOM are needed,
C     check to see if this atom contributes to the integral.
C
      IF (.NOT. EXPECT) THEN
         IF (ICENTA .NE. IATOM .AND. ICENTB .NE. IATOM .AND.
     *       ICENTC .NE. IATOM .AND. ICENTD .NE. IATOM) RETURN
      END IF
C
C     Special case: PERTURB and differentiation on symmetry
C     related centers. All atoms are differentiated separately.
C
      IF (.NOT. EXPECT) THEN
         IF (MULATM .EQ. 1) THEN
            NATOMS = 1
            IATOMS(1) = IATOM
            TRANS = .FALSE.
            IF (IATOM .EQ. NCENTA) THEN
               ISOPDR(1)   = ISOA
               SIGNDR(1,1) = SIGNAX
               SIGNDR(1,2) = SIGNAY
               SIGNDR(1,3) = SIGNAZ
            ELSE IF (IATOM .EQ. NCENTB) THEN
               ISOPDR(1)   = ISOB
               SIGNDR(1,1) = SIGNBX
               SIGNDR(1,2) = SIGNBY
               SIGNDR(1,3) = SIGNBZ
            ELSE IF (IATOM .EQ. NCENTC) THEN
               ISOPDR(1)   = ISOC
               SIGNDR(1,1) = SIGNCX
               SIGNDR(1,2) = SIGNCY
               SIGNDR(1,3) = SIGNCZ
            ELSE
               ISOPDR(1)   = ISOD
               SIGNDR(1,1) = SIGNDX
               SIGNDR(1,2) = SIGNDY
               SIGNDR(1,3) = SIGNDZ
            END IF
         ELSE
            DOA = IATOM .EQ. NCENTA
            DOB = IATOM .EQ. NCENTB
            DOC = IATOM .EQ. NCENTC
            DOD = IATOM .EQ. NCENTD
            TRANS = .FALSE.
C           DOD = IATOM .EQ. NCENTD .AND. .NOT.(DOA.AND.DOB.AND.DOC)
C           TRANS = IATOM .EQ. NCENTD .AND. (DOA.AND.DOB.AND.DOC)
            INT = 0
            IF (DOA) THEN
               CALL DIREA0(1,2,3,ONE,ONE,ONE)
               INT           = INT +  1
               IATOMS(INT)   = IATOM
               ISOPDR(INT)   = ISOA
               SIGNDR(INT,1) = SIGNAX
               SIGNDR(INT,2) = SIGNAY
               SIGNDR(INT,3) = SIGNAZ
            END IF
            IF (DOB) THEN
               CALL DIREB0(3*INT+1,3*INT+2,3*INT+3,ONE,ONE,ONE)
               INT           = INT +  1
               IATOMS(INT)   = IATOM
               ISOPDR(INT)   = ISOB
               SIGNDR(INT,1) = SIGNBX
               SIGNDR(INT,2) = SIGNBY
               SIGNDR(INT,3) = SIGNBZ
            END IF
            IF (DOC) THEN
               CALL DIREC0(3*INT+1,3*INT+2,3*INT+3,ONE,ONE,ONE)
               INT           = INT +  1
               IATOMS(INT)   = IATOM
               ISOPDR(INT)   = ISOC
               SIGNDR(INT,1) = SIGNCX
               SIGNDR(INT,2) = SIGNCY
               SIGNDR(INT,3) = SIGNCZ
            END IF
            IF (DOD) THEN
               CALL DIRED0(3*INT+1,3*INT+2,3*INT+3,ONE,ONE,ONE)
               INT           = INT +  1
               IATOMS(INT)   = IATOM
               ISOPDR(INT)   = ISOD
               SIGNDR(INT,1) = SIGNDX
               SIGNDR(INT,2) = SIGNDY
               SIGNDR(INT,3) = SIGNDZ
            END IF
            NATOMS = INT
            NINTYP = 3*INT
            IF (DOA .OR. DOB) THEN
               DPATH2 = .TRUE.
               DC102  = .TRUE.
               DC2E2  = .TRUE.
            END IF
            IF (DOC .OR. DOD) THEN
               DPATH1 = .TRUE.
               DC101  = .TRUE.
               DC2E1  = .TRUE.
            END IF
            GO TO 200
         END IF
      END IF
C
C
C     ***********************************************************
C     ********* Determine type of multicenter integral **********
C     ***********************************************************
C
      IF(ICENTA.EQ.ICENTB) THEN
         IF(ICENTA.EQ.ICENTC) THEN
            IF(ICENTA.EQ.ICENTD) THEN
               ATYPE = '1111'
               GO TO 111
            ELSE
               ATYPE = '1112'
               GO TO 212
            END IF
         ELSE
            IF(ICENTA.EQ.ICENTD) THEN
               ATYPE = '1121'
               GO TO 212
            ELSE IF(ICENTC.EQ.ICENTD) THEN
               ATYPE = '1122'
               GO TO 211
            ELSE
               ATYPE = '1123'
               GO TO 312
            END IF
         END IF
      ELSE
         IF(ICENTA.EQ.ICENTC) THEN
            IF(ICENTA.EQ.ICENTD) THEN
               ATYPE = '1211'
               GO TO 221
            ELSE IF(ICENTB.EQ.ICENTD) THEN
               ATYPE = '1212'
               GO TO 222
            ELSE
               ATYPE = '1213'
               GO TO 322
            END IF
         ELSE IF(ICENTB.EQ.ICENTC) THEN
            IF(ICENTA.EQ.ICENTD) THEN
               ATYPE = '1221'
               GO TO 222
            ELSE IF(ICENTB.EQ.ICENTD) THEN
               ATYPE = '1222'
               GO TO 221
            ELSE
               ATYPE = '1223'
               GO TO 322
            END IF
         ELSE
            IF(ICENTA.EQ.ICENTD) THEN
               ATYPE = '1231'
               GO TO 322
            ELSE IF(ICENTB.EQ.ICENTD) THEN
               ATYPE = '1232'
               GO TO 322
            ELSE IF(ICENTC.EQ.ICENTD) THEN
               ATYPE = '1233'
               GO TO 321
            ELSE
               ATYPE = '1234'
               GO TO 422
            END IF
         END IF
      END IF
C
C     ************************************
C     **************** 111 ***************
C     ************************************
C
  111 CONTINUE
         ONECEN = .TRUE.
      RETURN
C
C     ************************************
C     *************** 212 ****************
C     ************************************
C
  212 CONTINUE
         TWOCEN = .TRUE.
         DPATH1 = .TRUE.
         DC101  = .TRUE.
         DC2E1  = .TRUE.
         NINTYP = 0
C
C        ***** 1112 *****
C
         IF (ATYPE .EQ. '1112') THEN
            DERECD = .TRUE.
            IF (EXPECT) THEN
               NCENT1 = NCENTD
               NCENT2 = NCENTA
               ISO1   = ISOD
               ISO2   = ISOA
               SIGN1X = SIGNDX
               SIGN2X = SIGNAX
               SIGN1Y = SIGNDY
               SIGN2Y = SIGNAY
               SIGN1Z = SIGNDZ
               SIGN2Z = SIGNAZ
               IF (DERONE) THEN
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
                  NINTYP = NINTYP + 3
               END IF
               IF (DERTWO) THEN
                  CALL DIREDD(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
            ELSE
               IF (IATOM .EQ. ICENTA) THEN
                  CALL DIRED0(1,2,3,-ONE,-ONE,-ONE)
               ELSE
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
C
C        ***** 1121 *****
C
         ELSE
            DERECC = .TRUE.
            IF (EXPECT) THEN
               NCENT1 = NCENTC
               NCENT2 = NCENTA
               ISO1   = ISOC
               ISO2   = ISOA
               SIGN1X = SIGNCX
               SIGN2X = SIGNAX
               SIGN1Y = SIGNCY
               SIGN2Y = SIGNAY
               SIGN1Z = SIGNCZ
               SIGN2Z = SIGNAZ
               IF (DERONE) THEN
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
                  NINTYP = NINTYP + 3
               END IF
               IF (DERTWO) THEN
                  CALL DIRECC(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
            ELSE
               IF (IATOM .EQ. ICENTA) THEN
                  CALL DIREC0(1,2,3,-ONE,-ONE,-ONE)
               ELSE
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
         END IF
      GO TO 200
C
C     ***********************************
C     *************** 221 ***************
C     ***********************************
C
  221 CONTINUE
         TWOCEN = .TRUE.
         DPATH2 = .TRUE.
         DC102 = .TRUE.
         DC2E2 = .TRUE.
         NINTYP = 0
C
C        ***** 1211 *****
C
         IF (ATYPE .EQ. '1211') THEN
            IF (EXPECT) THEN
               NCENT1 = NCENTB
               NCENT2 = NCENTA
               ISO1   = ISOB
               ISO2   = ISOA
               SIGN1X = SIGNBX
               SIGN2X = SIGNAX
               SIGN1Y = SIGNBY
               SIGN2Y = SIGNAY
               SIGN1Z = SIGNBZ
               SIGN2Z = SIGNAZ
               IF (DERONE) THEN
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
                  NINTYP = NINTYP + 3
               END IF
               IF (DERTWO) THEN
                  CALL DIREBB(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
            ELSE
               IF (IATOM .EQ. ICENTA) THEN
                  CALL DIREB0(1,2,3,-ONE,-ONE,-ONE)
               ELSE
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
C
C        ***** 1222 *****
C
         ELSE
            IF (EXPECT) THEN
               NCENT1 = NCENTA
               NCENT2 = NCENTB
               ISO1   = ISOA
               ISO2   = ISOB
               SIGN1X = SIGNAX
               SIGN2X = SIGNBX
               SIGN1Y = SIGNAY
               SIGN2Y = SIGNBY
               SIGN1Z = SIGNAZ
               SIGN2Z = SIGNBZ
               IF (DERONE) THEN
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  NINTYP = NINTYP + 3
               END IF
               IF (DERTWO) THEN
                  CALL DIREAA(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
            ELSE
               IF (IATOM .EQ. ICENTA) THEN
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
               ELSE
                  CALL DIREA0(1,2,3,-ONE,-ONE,-ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
         END IF
      GO TO 200
C
C     ***********************************
C     *************** 211 ***************
C     ***********************************
C
  211 CONTINUE
         TWOCEN = .TRUE.
         CALL DP211(DPATH1,DPATH2)
         NINTYP = 0
         IF (EXPECT) THEN
            IF (DPATH1) THEN
               NCENT1 = NCENTC
               NCENT2 = NCENTA
               ISO1   = ISOC
               ISO2   = ISOA
               SIGN1X = SIGNCX
               SIGN2X = SIGNAX
               SIGN1Y = SIGNCY
               SIGN2Y = SIGNAY
               SIGN1Z = SIGNCZ
               SIGN2Z = SIGNAZ
               DC101  = .TRUE.
               DC2H1  = .TRUE.
               IF (DERONE) THEN
                  CALL DIREQ0(1,2,3,ONE,ONE,ONE)
                  NINTYP = NINTYP + 3
               END IF
               IF (DERTWO) THEN
                  CALL DIREQQ(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
            ELSE
               NCENT1 = NCENTA
               NCENT2 = NCENTC
               ISO1   = ISOA
               ISO2   = ISOC
               SIGN1X = SIGNAX
               SIGN2X = SIGNCX
               SIGN1Y = SIGNAY
               SIGN2Y = SIGNCY
               SIGN1Z = SIGNAZ
               SIGN2Z = SIGNCZ
               DC102 = .TRUE.
               DC2H2 = .TRUE.
               IF (DERONE) THEN
                  CALL DIREP0(1,2,3,ONE,ONE,ONE)
                  NINTYP = NINTYP + 3
               END IF
               IF (DERTWO) THEN
                  CALL DIREPP(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
            END IF
         ELSE
            IF (DPATH1) THEN
               DC101 = .TRUE.
               DC2H1 = .TRUE.
               IF (IATOM .EQ. ICENTA) THEN
                  CALL DIREQ0(1,2,3,-ONE,-ONE,-ONE)
               ELSE
                  CALL DIREQ0(1,2,3,ONE,ONE,ONE)
               END IF
            ELSE
               DC102 = .TRUE.
               DC2H2 = .TRUE.
               IF (IATOM .EQ. ICENTA) THEN
                  CALL DIREP0(1,2,3,ONE,ONE,ONE)
               ELSE
                  CALL DIREP0(1,2,3,-ONE,-ONE,-ONE)
               END IF
            END IF
            NINTYP = NINTYP + 3
         END IF
      GO TO 200
C
C     ***********************************
C     *************** 222 ***************
C     ***********************************
C
  222 CONTINUE
         TWOCEN = .TRUE.
         DPATH1 = .TRUE.
         DPATH2 = .TRUE.
         DC101  = .TRUE.
         DC2E1  = .TRUE.
         DC102  = .TRUE.
         DC2E2  = .TRUE.
         IF (DERTWO) THEN
            CALL DCROSS(CROSS1,CROSS2)
            DC1E1 = CROSS1
            DC1E2 = CROSS2
         END IF
         NINTYP = 0
C
C        ***** 1212 *****
C
         IF (ATYPE .EQ. '1212') THEN
            IF (EXPECT) THEN
               NCENT1 = NCENTA
               NCENT2 = NCENTB
               ISO1   = ISOA
               ISO2   = ISOB
               SIGN1X = SIGNAX
               SIGN2X = SIGNBX
               SIGN1Y = SIGNAY
               SIGN2Y = SIGNBY
               SIGN1Z = SIGNAZ
               SIGN2Z = SIGNBZ
               DERECC = .TRUE.
               IF (DERONE) THEN
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
                  NINTYP = NINTYP + 3
               END IF
               IF (DERTWO) THEN
                  DHCEX = .TRUE.
                  DHCEY = .TRUE.
                  DHCEZ = .TRUE.
                  DHCEX1 = .TRUE.
                  DHCEY1 = .TRUE.
                  DHCEZ1 = .TRUE.
                  CALL DIREAA(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIRECC(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREAC(4,5,6,5,7,8,6,8,9,
     *                        TWO,ONE,ONE,ONE,TWO,ONE,ONE,ONE,TWO)
                  NINTYP = NINTYP + 6
               END IF
            ELSE
               IF (IATOM .EQ. ICENTA) THEN
                  DERECC = .TRUE.
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
               ELSE
                  DERECD = .TRUE.
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
C
C        ***** 1221 *****
C
         ELSE
            IF (EXPECT) THEN
               NCENT1 = NCENTA
               NCENT2 = NCENTB
               ISO1   = ISOA
               ISO2   = ISOB
               SIGN1X = SIGNAX
               SIGN2X = SIGNBX
               SIGN1Y = SIGNAY
               SIGN2Y = SIGNBY
               SIGN1Z = SIGNAZ
               SIGN2Z = SIGNBZ
               DERECD = .TRUE.
               IF (DERONE) THEN
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
                  NINTYP = NINTYP + 3
               END IF
               IF (DERTWO) THEN
                  DHCEX = .TRUE.
                  DHCEY = .TRUE.
                  DHCEZ = .TRUE.
                  IF (CROSS1) THEN
                     DHCEX1 = .TRUE.
                     DHCEY1 = .TRUE.
                     DHCEZ1 = .TRUE.
                  ELSE
                     DHCEX2 = .TRUE.
                     DHCEY2 = .TRUE.
                     DHCEZ2 = .TRUE.
                  END IF
                  CALL DIREAA(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREDD(4,5,6,7,8,9,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREAD(4,5,6,5,7,8,6,8,9,
     *                        TWO,ONE,ONE,ONE,TWO,ONE,ONE,ONE,TWO)
                  NINTYP = NINTYP + 6
               END IF
            ELSE
               IF (IATOM .EQ. ICENTA) THEN
                  DERECD = .TRUE.
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
               ELSE
                  DERECC = .TRUE.
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
         END IF
      GO TO 200
C
C     ***********************************
C     *************** 312 ***************
C     ***********************************
C
  312 CONTINUE
         THRCEN = .TRUE.
         NINTYP = 0
         IF (EXPECT) THEN
            NCENT1 = NCENTC
            NCENT2 = NCENTD
            NCENT3 = NCENTA
            ISO1   = ISOC
            ISO2   = ISOD
            ISO3   = ISOA
            SIGN1X = SIGNCX
            SIGN2X = SIGNDX
            SIGN3X = SIGNAX
            SIGN1Y = SIGNCY
            SIGN2Y = SIGNDY
            SIGN3Y = SIGNAY
            SIGN1Z = SIGNCZ
            SIGN2Z = SIGNDZ
            SIGN3Z = SIGNAZ
            DPATH1 = .TRUE.
            DC101  = .TRUE.
            DC2E1  = .TRUE.
            DERECC = .TRUE.
            DERECD = .TRUE.
            IF (DERONE) THEN
               CALL DIREC0(1,2,3,ONE,ONE,ONE)
               CALL DIRED0(4,5,6,ONE,ONE,ONE)
               NINTYP = NINTYP + 6
            END IF
            IF (DERTWO) THEN
               CALL DIRECC(7,8,9,10,11,12,ONE,ONE,ONE,ONE,ONE,ONE)
               CALL DIREDD(13,14,15,16,17,18,ONE,ONE,ONE,ONE,ONE,ONE)
               CALL DIRECD(19,20,21,22,23,24,25,26,27,
     *                     ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
               NINTYP = NINTYP + 21
            END IF
         ELSE
            IF (IATOM .EQ. ICENTA) THEN
               DPATH2 = .TRUE.
               DC102  = .TRUE.
               DC2H2  = .TRUE.
               CALL DIREP0(1,2,3,ONE,ONE,ONE)
            ELSE IF (IATOM .EQ. ICENTC) THEN
               DPATH1 = .TRUE.
               DC101  = .TRUE.
               DC2E1  = .TRUE.
               DERECC = .TRUE.
               CALL DIREC0(1,2,3,ONE,ONE,ONE)
            ELSE
               DPATH1 = .TRUE.
               DC101  = .TRUE.
               DC2E1  = .TRUE.
               DERECD = .TRUE.
               CALL DIRED0(1,2,3,ONE,ONE,ONE)
            END IF
            NINTYP = NINTYP + 3
         END IF
      GO TO 200
C
C     ***********************************
C     *************** 322 ***************
C     ***********************************
C
  322 CONTINUE
         THRCEN = .TRUE.
         NINTYP = 0
         IF (EXPECT) THEN
            DPATH1 = .TRUE.
            DPATH2 = .TRUE.
            DC101  = .TRUE.
            DC2E1  = .TRUE.
            DC102  = .TRUE.
            DC2E2  = .TRUE.
            IF (DERTWO) THEN
               CALL DCROSS(CROSS1,CROSS2)
               DC1E1 = CROSS1
               DC1E2 = CROSS2
            END IF
         END IF
C
C        ***** 1213 *****
C
         IF (ATYPE .EQ. '1213') THEN
            IF (EXPECT) THEN
               NCENT1 = NCENTB
               NCENT2 = NCENTD
               NCENT3 = NCENTA
               ISO1   = ISOB
               ISO2   = ISOD
               ISO3   = ISOA
               SIGN1X = SIGNBX
               SIGN2X = SIGNDX
               SIGN3X = SIGNAX
               SIGN1Y = SIGNBY
               SIGN2Y = SIGNDY
               SIGN3Y = SIGNAY
               SIGN1Z = SIGNBZ
               SIGN2Z = SIGNDZ
               SIGN3Z = SIGNAZ
               DERECD = .TRUE.
               IF (DERONE) THEN
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
                  CALL DIRED0(4,5,6,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
               IF (DERTWO) THEN
                  DHCEX = .TRUE.
                  DHCEY = .TRUE.
                  DHCEZ = .TRUE.
                  DHCEX2 = .TRUE.
                  DHCEY2 = .TRUE.
                  DHCEZ2 = .TRUE.
                  CALL DIREBB(7,8,9,10,11,12,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREDD(13,14,15,16,17,18,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREBD(19,20,21,22,23,24,25,26,27,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 21
               END IF
            ELSE
               IF (IATOM .EQ. ICENTA) THEN
                  DPATH1 = .TRUE.
                  DPATH2 = .TRUE.
                  DC101  = .TRUE.
                  DC2E1  = .TRUE.
                  DC102  = .TRUE.
                  DC2E2  = .TRUE.
                  DERECC = .TRUE.
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
               ELSE IF (IATOM .EQ. ICENTB) THEN
                  DPATH2 = .TRUE.
                  DC102  = .TRUE.
                  DC2E2  = .TRUE.
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
               ELSE
                  DPATH1 = .TRUE.
                  DC101  = .TRUE.
                  DC2E1  = .TRUE.
                  DERECD = .TRUE.
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
C
C        ***** 1223 *****
C
         ELSE IF (ATYPE .EQ. '1223') THEN
            IF (EXPECT) THEN
               NCENT1 = NCENTA
               NCENT2 = NCENTD
               NCENT3 = NCENTB
               ISO1   = ISOA
               ISO2   = ISOD
               ISO3   = ISOB
               SIGN1X = SIGNAX
               SIGN2X = SIGNDX
               SIGN3X = SIGNBX
               SIGN1Y = SIGNAY
               SIGN2Y = SIGNDY
               SIGN3Y = SIGNBY
               SIGN1Z = SIGNAZ
               SIGN2Z = SIGNDZ
               SIGN3Z = SIGNBZ
               DERECD = .TRUE.
               IF (DERONE) THEN
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIRED0(4,5,6,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
               IF (DERTWO) THEN
                  DHCEX = .TRUE.
                  DHCEY = .TRUE.
                  DHCEZ = .TRUE.
                  IF (CROSS1) THEN
                     DHCEX1 = .TRUE.
                     DHCEY1 = .TRUE.
                     DHCEZ1 = .TRUE.
                  ELSE
                     DHCEX2 = .TRUE.
                     DHCEY2 = .TRUE.
                     DHCEZ2 = .TRUE.
                  END IF
                  CALL DIREAA(7,8,9,10,11,12,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREDD(13,14,15,16,17,18,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREAD(19,20,21,22,23,24,25,26,27,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 21
               END IF
            ELSE
               IF (IATOM .EQ. ICENTA) THEN
                  DPATH2 = .TRUE.
                  DC102 = .TRUE.
                  DC2E2 = .TRUE.
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
               ELSE IF (IATOM .EQ. ICENTB) THEN
                  DPATH1 = .TRUE.
                  DPATH2 = .TRUE.
                  DC101  = .TRUE.
                  DC2E1  = .TRUE.
                  DC102  = .TRUE.
                  DC2E2  = .TRUE.
                  DERECC = .TRUE.
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
               ELSE
                  DPATH1 = .TRUE.
                  DC101  = .TRUE.
                  DC2E1  = .TRUE.
                  DERECD = .TRUE.
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
C
C        ***** 1231 *****
C
         ELSE IF (ATYPE .EQ. '1231') THEN
            IF (EXPECT) THEN
               NCENT1 = NCENTB
               NCENT2 = NCENTC
               NCENT3 = NCENTA
               ISO1   = ISOB
               ISO2   = ISOC
               ISO3   = ISOA
               SIGN1X = SIGNBX
               SIGN2X = SIGNCX
               SIGN3X = SIGNAX
               SIGN1Y = SIGNBY
               SIGN2Y = SIGNCY
               SIGN3Y = SIGNAY
               SIGN1Z = SIGNBZ
               SIGN2Z = SIGNCZ
               SIGN3Z = SIGNAZ
               DERECC = .TRUE.
               IF (DERONE) THEN
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
                  CALL DIREC0(4,5,6,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
               IF (DERTWO) THEN
                  DHCEX = .TRUE.
                  DHCEY = .TRUE.
                  DHCEZ = .TRUE.
                  IF (CROSS1) THEN
                     DHCEX2 = .TRUE.
                     DHCEY2 = .TRUE.
                     DHCEZ2 = .TRUE.
                  ELSE
                     DHCEX1 = .TRUE.
                     DHCEY1 = .TRUE.
                     DHCEZ1 = .TRUE.
                  END IF
                  CALL DIREBB(7,8,9,10,11,12,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIRECC(13,14,15,16,17,18,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREBC(19,20,21,22,23,24,25,26,27,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 21
               END IF
            ELSE
               IF (IATOM .EQ. ICENTA) THEN
                  DPATH1 = .TRUE.
                  DPATH2 = .TRUE.
                  DC101  = .TRUE.
                  DC102  = .TRUE.
                  DC2E1  = .TRUE.
                  DC2E2  = .TRUE.
                  DERECD = .TRUE.
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
               ELSE IF (IATOM .EQ. ICENTB) THEN
                  DPATH2 = .TRUE.
                  DC102  = .TRUE.
                  DC2E2  = .TRUE.
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
               ELSE
                  DPATH1 = .TRUE.
                  DC101  = .TRUE.
                  DC2E1  = .TRUE.
                  DERECC = .TRUE.
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
C
C        ***** 1232 *****
C
         ELSE
            IF (EXPECT) THEN
               NCENT1 = NCENTA
               NCENT2 = NCENTC
               NCENT3 = NCENTB
               ISO1   = ISOA
               ISO2   = ISOC
               ISO3   = ISOB
               SIGN1X = SIGNAX
               SIGN2X = SIGNCX
               SIGN3X = SIGNBX
               SIGN1Y = SIGNAY
               SIGN2Y = SIGNCY
               SIGN3Y = SIGNBY
               SIGN1Z = SIGNAZ
               SIGN2Z = SIGNCZ
               SIGN3Z = SIGNBZ
               DERECC = .TRUE.
               IF (DERONE) THEN
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIREC0(4,5,6,ONE,ONE,ONE)
                  NINTYP = NINTYP + 6
               END IF
               IF (DERTWO) THEN
                  DHCEX = .TRUE.
                  DHCEY = .TRUE.
                  DHCEZ = .TRUE.
                  DHCEX1 = .TRUE.
                  DHCEY1 = .TRUE.
                  DHCEZ1 = .TRUE.
                  CALL DIREAA(7,8,9,10,11,12,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIRECC(13,14,15,16,17,18,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREAC(19,20,21,22,23,24,25,26,27,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 21
               END IF
            ELSE
               IF (IATOM .EQ. ICENTA) THEN
                  DPATH2 = .TRUE.
                  DC102 = .TRUE.
                  DC2E2 = .TRUE.
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
               ELSE IF (IATOM .EQ. ICENTB) THEN
                  DPATH1 = .TRUE.
                  DPATH2 = .TRUE.
                  DC101  = .TRUE.
                  DC102  = .TRUE.
                  DC2E1  = .TRUE.
                  DC2E2  = .TRUE.
                  DERECD = .TRUE.
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
                  CALL DIRED0(1,2,3,ONE,ONE,ONE)
               ELSE
                  DPATH1 = .TRUE.
                  DC101  = .TRUE.
                  DC2E1  = .TRUE.
                  DERECC = .TRUE.
                  CALL DIREC0(1,2,3,ONE,ONE,ONE)
               END IF
               NINTYP = NINTYP + 3
            END IF
         END IF
      GO TO 200
C
C     ***********************************
C     *************** 321 ***************
C     ***********************************
C
  321 CONTINUE
         THRCEN = .TRUE.
         NINTYP = 0
         IF (EXPECT) THEN
            NCENT1 = NCENTA
            NCENT2 = NCENTB
            NCENT3 = NCENTC
            ISO1   = ISOA
            ISO2   = ISOB
            ISO3   = ISOC
            SIGN1X = SIGNAX
            SIGN2X = SIGNBX
            SIGN3X = SIGNCX
            SIGN1Y = SIGNAY
            SIGN2Y = SIGNBY
            SIGN3Y = SIGNCY
            SIGN1Z = SIGNAZ
            SIGN2Z = SIGNBZ
            SIGN3Z = SIGNCZ
            DPATH2 = .TRUE.
            DC102 = .TRUE.
            DC2E2 = .TRUE.
            IF (DERONE) THEN
               CALL DIREA0(1,2,3,ONE,ONE,ONE)
               CALL DIREB0(4,5,6,ONE,ONE,ONE)
               NINTYP = NINTYP + 6
            END IF
            IF (DERTWO) THEN
               CALL DIREAA(7,8,9,10,11,12,ONE,ONE,ONE,ONE,ONE,ONE)
               CALL DIREBB(13,14,15,16,17,18,ONE,ONE,ONE,ONE,ONE,ONE)
               CALL DIREAB(19,20,21,22,23,24,25,26,27,
     *                     ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
               NINTYP = NINTYP + 21
            END IF
         ELSE
            IF (IATOM .EQ. ICENTA) THEN
               DPATH2 = .TRUE.
               DC102 = .TRUE.
               DC2E2 = .TRUE.
               CALL DIREA0(1,2,3,ONE,ONE,ONE)
            ELSE IF (IATOM .EQ. ICENTB) THEN
               DPATH2 = .TRUE.
               DC102 = .TRUE.
               DC2E2 = .TRUE.
               CALL DIREB0(1,2,3,ONE,ONE,ONE)
            ELSE
               DPATH1 = .TRUE.
               DC101  = .TRUE.
               DC2H1  = .TRUE.
               CALL DIREQ0(1,2,3,ONE,ONE,ONE)
            END IF
            NINTYP = NINTYP + 3
         END IF
      GO TO 200
C
C     ***********************************
C     *************** 422 ***************
C     ***********************************
C
  422 CONTINUE
         FOUCEN = .TRUE.
         NINTYP = 0
         IF (EXPECT) THEN
            DPATH1 = .TRUE.
            DPATH2 = .TRUE.
            DC101  = .TRUE.
            DC2E1  = .TRUE.
            DC102  = .TRUE.
            DC2E2  = .TRUE.
            CALL DP422(CROSS1,CROSS2)
            IF (DERTWO) THEN
               DC1E1 = CROSS1
               DC1E2 = CROSS2
            END IF
            IF (CROSS1) THEN
               NCENT1 = NCENTB
               NCENT2 = NCENTC
               NCENT3 = NCENTD
               NCENT4 = NCENTA
               ISO1   = ISOB
               ISO2   = ISOC
               ISO3   = ISOD
               ISO4   = ISOA
               SIGN1X = SIGNBX
               SIGN2X = SIGNCX
               SIGN3X = SIGNDX
               SIGN4X = SIGNAX
               SIGN1Y = SIGNBY
               SIGN2Y = SIGNCY
               SIGN3Y = SIGNDY
               SIGN4Y = SIGNAY
               SIGN1Z = SIGNBZ
               SIGN2Z = SIGNCZ
               SIGN3Z = SIGNDZ
               SIGN4Z = SIGNAZ
               DERECC = .TRUE.
               DERECD = .TRUE.
               IF (DERONE) THEN
                  CALL DIREB0(1,2,3,ONE,ONE,ONE)
                  CALL DIREC0(4,5,6,ONE,ONE,ONE)
                  CALL DIRED0(7,8,9,ONE,ONE,ONE)
                  NINTYP = NINTYP + 9
               END IF
               IF (DERTWO) THEN
                  DHCEX = .TRUE.
                  DHCEY = .TRUE.
                  DHCEZ = .TRUE.
                  DHCEX2 = .TRUE.
                  DHCEY2 = .TRUE.
                  DHCEZ2 = .TRUE.
                  CALL DIREBB(10,11,12,13,14,15,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIRECC(16,17,18,19,20,21,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREDD(22,23,24,25,26,27,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREBC(28,29,30,31,32,33,34,35,36,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREBD(37,38,39,40,41,42,43,44,45,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIRECD(46,47,48,49,50,51,52,53,54,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 45
               END IF
            ELSE
               NCENT1 = NCENTA
               NCENT2 = NCENTB
               NCENT3 = NCENTC
               NCENT4 = NCENTD
               ISO1   = ISOA
               ISO2   = ISOB
               ISO3   = ISOC
               ISO4   = ISOD
               SIGN1X = SIGNAX
               SIGN2X = SIGNBX
               SIGN3X = SIGNCX
               SIGN4X = SIGNDX
               SIGN1Y = SIGNAY
               SIGN2Y = SIGNBY
               SIGN3Y = SIGNCY
               SIGN4Y = SIGNDY
               SIGN1Z = SIGNAZ
               SIGN2Z = SIGNBZ
               SIGN3Z = SIGNCZ
               SIGN4Z = SIGNDZ
               DERECC = .TRUE.
               IF (DERONE) THEN
                  CALL DIREA0(1,2,3,ONE,ONE,ONE)
                  CALL DIREB0(4,5,6,ONE,ONE,ONE)
                  CALL DIREC0(7,8,9,ONE,ONE,ONE)
                  NINTYP = NINTYP + 9
               END IF
               IF (DERTWO) THEN
                  DHCEX = .TRUE.
                  DHCEY = .TRUE.
                  DHCEZ = .TRUE.
                  DHCEX1 = .TRUE.
                  DHCEY1 = .TRUE.
                  DHCEZ1 = .TRUE.
                  CALL DIREAA(10,11,12,13,14,15,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREBB(16,17,18,19,20,21,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIRECC(22,23,24,25,26,27,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREAB(28,29,30,31,32,33,34,35,36,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREAC(37,38,39,40,41,42,43,44,45,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  CALL DIREBC(46,47,48,49,50,51,52,53,54,
     *                        ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE,ONE)
                  NINTYP = NINTYP + 45
               END IF
            END IF
         ELSE
            IF (IATOM .EQ. ICENTA) THEN
               DPATH2 = .TRUE.
               DC102 = .TRUE.
               DC2E2 = .TRUE.
               CALL DIREA0(1,2,3,ONE,ONE,ONE)
            ELSE IF (IATOM .EQ. ICENTB) THEN
               DPATH2 = .TRUE.
               DC102  = .TRUE.
               DC2E2  = .TRUE.
               CALL DIREB0(1,2,3,ONE,ONE,ONE)
            ELSE IF (IATOM .EQ. ICENTC) THEN
               DPATH1 = .TRUE.
               DC101  = .TRUE.
               DC2E1  = .TRUE.
               DERECC = .TRUE.
               CALL DIREC0(1,2,3,ONE,ONE,ONE)
            ELSE
               DPATH1 = .TRUE.
               DC101  = .TRUE.
               DC2E1  = .TRUE.
               DERECD = .TRUE.
               CALL DIRED0(1,2,3,ONE,ONE,ONE)
            END IF
            NINTYP = NINTYP + 3
         END IF
  200 CONTINUE
      NTOTAL = NINTYP*NCCINT
      IPRINT =0
      IF (IPRINT .GE. 10) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 1005) EXPECT, IATOM
         WRITE (LUPRI, 1010) ONECEN, TWOCEN, THRCEN, FOUCEN
         WRITE (LUPRI, 1020) ATYPE
         WRITE (LUPRI, 1030) DPATH1, DPATH2
         WRITE (LUPRI, 1035) NCCINT
         WRITE (LUPRI, 1040) NTOTAL
         WRITE (LUPRI, 1045) CROSS1, CROSS2
         WRITE (LUPRI, 1050) DC101, DC102
         WRITE (LUPRI, 1060) DC1H1, DC1H2
         WRITE (LUPRI, 1070) DC1E1, DC1E2
         WRITE (LUPRI, 1080) DC2H1, DC2H2
         WRITE (LUPRI, 1090) DC2E1, DC2E2
         WRITE (LUPRI, 1095) DERECC, DERECD
         WRITE (LUPRI, 1100) (DEOVEC(I),I=1,54)
         WRITE (LUPRI, 1110) (DETVEC(I),I=1,36)
         WRITE (LUPRI, 1120) (DHOVEC(I),I=1,18)
         WRITE (LUPRI, 1130) (IEOVEC(I),I=1,54)
         WRITE (LUPRI, 1140) (IETVEC(I),I=1,36)
         WRITE (LUPRI, 1150) (IHOVEC(I),I=1,18)
         WRITE (LUPRI, 1160) (IDNINT(FEOVEC(I)),I=1,54)
         WRITE (LUPRI, 1170) (IDNINT(FETVEC(I)),I=1,36)
         WRITE (LUPRI, 1180) (IDNINT(FHOVEC(I)),I=1,18)
      END IF
 1000 FORMAT (/,1X,'<<<<<<<<<< SUBROUTINE DIRECT >>>>>>>>>>',/)
 1005 FORMAT (1X,'EXPECT/IATOM     ',L5,I5)
 1010 FORMAT (1X,'ONECEN,... FOUCEN',4L5)
 1020 FORMAT (1X,'INTEGRAL TYPE    ',A4)
 1030 FORMAT (1X,'DPATH1/2         ',2L5)
 1035 FORMAT (1X,'NCCINT           ',I5)
 1040 FORMAT (1X,'NTOTAL           ',I5)
 1045 FORMAT (1X,'CROSS1/2         ',2L5)
 1050 FORMAT (1X,'DC101/2          ',2L5)
 1060 FORMAT (1X,'DC1H1/2          ',2L5)
 1070 FORMAT (1X,'DC1E1/2          ',2L5)
 1080 FORMAT (1X,'DC2H1/2          ',2L5)
 1090 FORMAT (1X,'DC2E1/2          ',2L5)
 1095 FORMAT (1X,'DERECC/D         ',2L5)
 1100 FORMAT (1X,'COMMON /DEODIR/',
     *        1X,4L5,/,2(17X,4L5,/),2(17X,3L5,/),4(17X,4L5,/),
     *        2(17X,3L5,/),2(17X,4L5,/),2(17X,3L5,/))
 1110 FORMAT (1X,'COMMON /DETDIR/',1X,3L5,/,11(17X,3L5,/))
 1120 FORMAT (1X,'COMMON /DHODIR/',1X,3L5,/,17X,3L5,/,2(17X,6L5,/))
 1130 FORMAT (1X,'COMMON /DEOADR/',
     *        1X,4I5,/,2(17X,4I5,/),2(17X,3I5,/),4(17X,4I5,/),
     *        2(17X,3I5,/),2(17X,4I5,/),2(17X,3I5,/))
 1140 FORMAT (1X,'COMMON /DETADR/',1X,3I5,/,11(17X,3I5,/))
 1150 FORMAT (1X,'COMMON /DHOADR/',1X,3I5,/,17X,3I5,/,2(17X,6I5,/))
 1160 FORMAT (1X,'COMMON /DEOFAC/',
     *        1X,4I5,/,2(17X,4I5,/),2(17X,3I5,/),4(17X,4I5,/),
     *        2(17X,3I5,/),2(17X,4I5,/),2(17X,3I5,/))
 1170 FORMAT (1X,'COMMON /DETFAC/',1X,3I5,/,11(17X,3I5,/))
 1180 FORMAT (1X,'COMMON /DHOFAC/',1X,3I5,/,17X,3I5,/,2(17X,6I5,/))
      IPRINT = 0
      RETURN
      END
