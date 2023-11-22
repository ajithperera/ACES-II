         SUBROUTINE  NLO__PRINT_NBO_RESULTS
     +
     +                    ( UNITID,
     +                      NBAS,NATOM,
     +                      MXNCEN,
     +                      NHB,NHA,
     +                      NCA,NRA,
     +                      NCB,NLB,NBB,NEB,NAB,NYB,NRB,
     +                      ZATOM,
     +                      ATNCB,ATCIDX,ATCOFF,
     +                      ATNRB,ATRIDX,ATROFF,
     +                      CSHELL,RSHELL,
     +                      BDNCEN,BDCEN,BDBAS,
     +                      NBOBD,
     +                      LOCAL,
     +                      W,B )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__PRINT_NBO_RESULTS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine prints details of the set of obtained
C                natural bond orbitals (NBOs) to a file specified
C                by its unit identification number. When entering this
C                routine, the NBO expansion coefficients in terms of
C                the atomic NHOs have been ordered such that all
C                Lone-pair, Bond, Empty-pair, Antibond and Rydberg NBOs
C                are clustered together and each set is in the order
C                mentioned. Also the overall weight vector has been
C                ordered to match the order of the NBO expansion
C                coefficients in the NBO part, and additionally the
C                remaining Core and Rydberg NAOs not used for NHO
C                construction were placed in front and at the end,
C                respectively. A picture of what is present at this
C                stage helps:
C
C                Order of NBO expansion coefficients (B matrix):
C
C                        |  LP  |   B  |  EP  |  AB  |RY(nbo)|
C                          (NLB)  (NBB)  (NEB)  (NAB)  (NYB)
C
C                Order of overall weights (W array):
C
C                 |   C  |  LP  |   B  |  EP  |  AB  |RY(nbo)|RY(nao)|
C                   (NCB)  (NLB)  (NBB)  (NEB)  (NAB)  (NYB)   (NRB)
C
C                where the variables below in parenthesis indicate the
C                number of elements present in each set.
C
C                The following is printed:
C
C                   1) The occupation numbers (weights) of the
C                      NBOs and their characterization as Core,
C                      Lone-pair (1-center), Bond and Antibond
C                      (x-center, x=2,3,...), and Rydberg type
C                      NBOs.
C
C                   2) A population chart indicating the NBO/NAO
C                      populations broken down in Lewis, Delocalized
C                      and Rest populations:
C
C                       Lewis: Core + Lone-pairs + 1- and 2-cen Bonds
C                       Delocalized: >2-cen Bonds
C                       Rest: Antibonds + Rydbergs (NBO and NAO)
C
C                  Input:
C
C                    UNITID       =  printout unit identification #
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atoms
C                    MXNCEN       =  maximum # of atomic valence
C                                    centers per bond.
C                    NHB          =  total # of Hybrid NAOs
C                    NHA          =  total # of hybrid atoms
C                    NCA          =  total # of NAO Core atoms
C                    NRA          =  total # of NAO Rydberg atoms
C                    NxB          =  total # of NAO Core, NBO Lone-
C                                    pairs, NBO Bonds, NBO Empty-pairs,
C                                    NBO Antibonds, NBO Rydbergs and NAO
C                                    Rydbergs found(x=C,L,B,E,A,Y,R).
C                    ZATOM (I)    =  atomic number for I-th atom.
C                    ATNCB (A)    =  # of Core NHOs on core atom A.
C                    ATCIDX (A)   =  atomic index for core atom A.
C                    ATCOFF (A)   =  index offset for Core NHOs for
C                                    core atom A. This index is equal
C                                    to the total number of Core NHOs
C                                    on all core atoms preceeding core
C                                    atom A.
C                    ATNRB (A)    =  # of Rydberg NHOs on core atom A.
C                    ATRIDX (A)   =  atomic index for Rydberg atom A.
C                    ATROFF (A)   =  index offset for Rydberg NHOs for
C                                    Rydberg atom A. This index is equal
C                                    to the total number of Rydberg NHOs
C                                    on all Rydberg atoms preceeding
C                                    Rydberg atom A.
C                    CSHELL (I)   =  l-shell type for the I-th Core NHO.
C                    RSHELL (I)   =  l-shell type for the I-th Rydberg
C                                    NHO.
C                    BDNCEN (I)   =  # of atomic centers for I-th NHO
C                                    bond.
C                    BDCEN (I,J)  =  I-th atomic center index for J-th
C                                    NHO bond.
C                    BDBAS (I,J)  =  I-th global basis (NHO) index for
C                                    J-th NHO bond.
C                    NBOBD (I)    =  contains the NHO bond index number
C                                    number for the I-th NBO in the NHB
C                                    part in NLB/NBB/NEB/NAB/NYB order.
C                    LOCAL (A,I)  =  NBO atom locality map for each
C                                    atom A and with NBOs in NCB/NLB/
C                                    NBB/NEB/NAB/NYB/NRB order.
C                    W            =  NBO weight vector with elements
C                                    in NCB/NLB/NBB/NEB/NAB/NYB/NRB
C                                    order.
C                    B (I,J)      =  MXNCEN x NVB matrix containing
C                                    the J-th NBO expansion coefficients
C                                    in terms of the I-th atomic NHOs
C                                    forming the J-th NBO. The NBO
C                                    column index is in NLB/NBB/NEB/
C                                    NAB/NYB order.
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         CHARACTER*2   ATCHAR
         CHARACTER*1   BAR,DASH
         CHARACTER*3   LCHAR
         CHARACTER*14  NBOCHAR

         CHARACTER*2   ATSYMB  (1:104)
         CHARACTER*3   LSYMB   (0:7)
         CHARACTER*14  NBOSYMB (1:7)

         LOGICAL     FIRST

         INTEGER     ATOM,ATOLD
         INTEGER     BDIDX
         INTEGER     I,N
         INTEGER     LTYPE
         INTEGER     MXNCEN
         INTEGER     NATOM
         INTEGER     NBAS
         INTEGER     NBEG,NEND
         INTEGER     NBO,NHO
         INTEGER     NCA,NRA
         INTEGER     NCB,NRB
         INTEGER     NCBA,NRBA
         INTEGER     NCEN
         INTEGER     NHA,NHB
         INTEGER     NLB,NBB,NEB,NAB,NYB
         INTEGER     OFF
         INTEGER     RYD
         INTEGER     UNITID
         INTEGER     ZVAL

         INTEGER     ATNCB  (1:NCA  )
         INTEGER     ATCIDX (1:NCA  )
         INTEGER     ATCOFF (1:NCA  )
         INTEGER     ATNRB  (1:NRA  )
         INTEGER     ATRIDX (1:NRA  )
         INTEGER     ATROFF (1:NRA  )
         INTEGER     BDNCEN (1:NHB  )
         INTEGER     CSHELL (1:NCB  )
         INTEGER     RSHELL (1:NRB  )
         INTEGER     NBOBD  (1:NHB  )
         INTEGER     ZATOM  (1:NATOM)

         INTEGER     BDBAS   (1:NHA,1:NHB)
         INTEGER     BDCEN   (1:NHA,1:NHB)

         DOUBLE PRECISION  COEFF
         DOUBLE PRECISION  LEWIS,DELOC,REST
         DOUBLE PRECISION  LOCVAL
         DOUBLE PRECISION  WEIGHT
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  W (1:NBAS)

         DOUBLE PRECISION  B     (1:MXNCEN,1:NHB )
         DOUBLE PRECISION  LOCAL (1:NATOM ,1:NBAS)

         DATA ATSYMB /' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     +                'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca',
     +                'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     +                'Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr',
     +                'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     +                'Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd',
     +                'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     +                'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     +                'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     +                'Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     +                'Md','No','Lr','xx'/
         DATA  BAR     /'|'/
         DATA  DASH    /'-'/
         DATA  LSYMB   /' s ',' p ',' d ',' f ',' g ',' h ',' i ','> i'/
         DATA  NBOSYMB /'   Core NAO   ',
     +                  ' Lone-pair NBO',
     +                  '   Bond NBO   ',
     +                  'Empty-pair NBO',
     +                  ' Anti-bond NBO',
     +                  '  Rydberg NBO ',
     +                  '  Rydberg NAO '/

         PARAMETER  (ZERO = 0.D0)

CCC hughes

C         integer alp, ppp

C         integer zimtx(natom,nbas)

C         do alp = 1, natom
C            do ppp = 1, nbas
C               zimtx(alp,ppp) = 0 
C            end do
C         end do

CCC hughes

C
C
C------------------------------------------------------------------------
C
C
C             ...print out title.
C
C
         WRITE (UNITID,8000) 'Natural bond orbital analysis'
 8000    FORMAT (//,30X,A29,//)
C
C
C             ...print out header line.
C
C
         WRITE (UNITID,9000) 'NBO #',BAR,'Atom',BAR,'Atom #',BAR,
     +                       'NHO #',BAR,'NHO coeffs',BAR,
     +                       'Orb-type',BAR,'Occupancy',BAR,'Locality'
         WRITE (UNITID,9010) (DASH,I=1,90)
C
C
C             ...print out Core NAO information.
C
C
         LEWIS = ZERO
         DELOC = ZERO
         REST  = ZERO

         IF (NCA.GT.0) THEN
             NBOCHAR = NBOSYMB (1)

             DO 1000 N = 1,NCA
                OFF = ATCOFF (N)
                NCBA = ATNCB (N)
                ATOM = ATCIDX (N)

                ZVAL = ZATOM (ATOM)
                IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                    ATCHAR = ATSYMB (104)
                ELSE
                    ATCHAR = ATSYMB (ZVAL)
                END IF

                FIRST = .TRUE.
                DO 100 I = 1,NCBA
                   NBO = OFF + I
                   WEIGHT = W (NBO)
                   LOCVAL = LOCAL (ATOM,NBO)
                   LEWIS = LEWIS + WEIGHT

                   LTYPE = CSHELL (NBO)
                   IF (LTYPE.GT.6) THEN
                       LCHAR = LSYMB (7)
                   ELSE
                       LCHAR = LSYMB (LTYPE)
                   END IF

                   IF (FIRST) THEN
                       WRITE (UNITID,9020) NBO,BAR,ATCHAR,BAR,ATOM,BAR,
     +                                     BAR,LCHAR,BAR,NBOCHAR,BAR,
     +                                     WEIGHT,BAR,LOCVAL

CCC hughes

C                       zimtx(atom,nbo) = 2
C                       print *, 'atom = ', atom, 'nbo = ',nbo,
C     & 'zi = ', zimtx(atom,nbo)

CCC hughes

                       FIRST = .FALSE.
                   ELSE
                       WRITE (UNITID,9030) NBO,BAR,BAR,BAR,BAR,LCHAR,
     +                                     BAR,NBOCHAR,BAR,WEIGHT,BAR,
     +                                     LOCVAL

CCC hughes

C                       print *, 'not sure'

CCC hughes

                   END IF
  100           CONTINUE

                WRITE (UNITID,9010) (DASH,I=1,90)
 1000        CONTINUE
         END IF
C
C
C             ...print out Lone-pair NBO information.
C
C
         NBO = NCB
         NEND = 0

         IF (NLB.GT.0) THEN
             NBEG = 1
             NEND = NLB
             NBOCHAR = NBOSYMB (2)

             DO 2000 N = NBEG,NEND
                NBO = NBO + 1
                BDIDX = NBOBD (N)
                WEIGHT = W (NBO)
                LEWIS = LEWIS + WEIGHT
                NHO = BDBAS (1,BDIDX)
                ATOM = BDCEN (1,BDIDX)
                COEFF = B (1,N)
                LOCVAL = LOCAL (ATOM,NBO)

                ZVAL = ZATOM (ATOM)
                IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                    ATCHAR = ATSYMB (104)
                ELSE
                    ATCHAR = ATSYMB (ZVAL)
                END IF

                WRITE (UNITID,9040) NBO,BAR,ATCHAR,BAR,ATOM,BAR,NHO,
     +                              BAR,COEFF,BAR,NBOCHAR,BAR,WEIGHT,
     +                              BAR,LOCVAL

CCC hughes

C                       zimtx(atom,nbo) = 2
C                       print *, 'atom = ', atom, 'nbo = ',nbo,
C     & 'zi = ', zimtx(atom,nbo)

CCC hughes


                WRITE (UNITID,9010) (DASH,I=1,90)
 2000        CONTINUE
         END IF
C
C
C             ...print out Bond NBO information.
C
C
         IF (NBB.GT.0) THEN
             NBEG = NEND + 1
             NEND = NEND + NBB
             NBOCHAR = NBOSYMB (3)

             DO 3000 N = NBEG,NEND
                NBO = NBO + 1
                BDIDX = NBOBD (N)
                NCEN = BDNCEN (BDIDX)
                WEIGHT = W (NBO)

                IF (NCEN.EQ.2) THEN
                    LEWIS = LEWIS + WEIGHT
                ELSE
                    DELOC = DELOC + WEIGHT
                END IF

                LOCVAL = ZERO
                DO 30 I = 1,NCEN
                   ATOM = BDCEN (I,BDIDX)
                   LOCVAL = LOCVAL + LOCAL (ATOM,NBO)
   30           CONTINUE

                FIRST = .TRUE.
                DO 300 I = 1,NCEN
                   NHO = BDBAS (I,BDIDX)
                   ATOM = BDCEN (I,BDIDX)
                   COEFF = B (I,N)

                   ZVAL = ZATOM (ATOM)
                   IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                       ATCHAR = ATSYMB (104)
                   ELSE
                       ATCHAR = ATSYMB (ZVAL)
                   END IF

                   IF (FIRST) THEN
                       WRITE (UNITID,9040) NBO,BAR,ATCHAR,BAR,ATOM,BAR,
     +                                     NHO,BAR,COEFF,BAR,NBOCHAR,
     +                                     BAR,WEIGHT,BAR,LOCVAL

CCC hughes

C                       zimtx(atom,nbo) = 1
C                       print *, 'atom = ', atom, 'nbo = ',nbo,
C     & 'zi = ', zimtx(atom,nbo)

CCC hughes

                       FIRST = .FALSE.
                   ELSE
                       WRITE (UNITID,9050) BAR,ATCHAR,BAR,ATOM,BAR,
     +                                     NHO,BAR,COEFF,BAR,BAR,BAR

CCC hughes

C                       zimtx(atom,nbo) = 1
C                       print *, 'atom = ', atom, 'nbo = ',nbo,
C     & 'zi = ', zimtx(atom,nbo)

CCC hughes

                   END IF
  300           CONTINUE

                WRITE (UNITID,9010) (DASH,I=1,90)
 3000        CONTINUE
         END IF
C
C
C             ...print out Empty-pair NBO information.
C
C
         IF (NEB.GT.0) THEN
             NBEG = NEND + 1
             NEND = NEND + NEB
             NBOCHAR = NBOSYMB (4)

             DO 4000 N = NBEG,NEND
                NBO = NBO + 1
                BDIDX = NBOBD (N)
                WEIGHT = W (NBO)
                REST = REST + WEIGHT
                NHO = BDBAS (1,BDIDX)
                ATOM = BDCEN (1,BDIDX)
                COEFF = B (1,N)
                LOCVAL = LOCAL (ATOM,NBO)

                ZVAL = ZATOM (ATOM)
                IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                    ATCHAR = ATSYMB (104)
                ELSE
                    ATCHAR = ATSYMB (ZVAL)
                END IF

                WRITE (UNITID,9040) NBO,BAR,ATCHAR,BAR,ATOM,BAR,NHO,
     +                              BAR,COEFF,BAR,NBOCHAR,BAR,WEIGHT,
     +                              BAR,LOCVAL

                WRITE (UNITID,9010) (DASH,I=1,90)
 4000        CONTINUE
         END IF
C
C
C             ...print out Antibond NBO information.
C
C
         IF (NAB.GT.0) THEN
             NBEG = NEND + 1
             NEND = NEND + NAB
             NBOCHAR = NBOSYMB (5)

             DO 5000 N = NBEG,NEND
                NBO = NBO + 1
                BDIDX = NBOBD (N)
                NCEN = BDNCEN (BDIDX)
                WEIGHT = W (NBO)
                REST = REST + WEIGHT

                LOCVAL = ZERO
                DO 50 I = 1,NCEN
                   ATOM = BDCEN (I,BDIDX)
                   LOCVAL = LOCVAL + LOCAL (ATOM,NBO)
   50           CONTINUE

                FIRST = .TRUE.
                DO 500 I = 1,NCEN
                   NHO = BDBAS (I,BDIDX)
                   ATOM = BDCEN (I,BDIDX)
                   COEFF = B (I,N)

                   ZVAL = ZATOM (ATOM)
                   IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                       ATCHAR = ATSYMB (104)
                   ELSE
                       ATCHAR = ATSYMB (ZVAL)
                   END IF

                   IF (FIRST) THEN
                       WRITE (UNITID,9040) NBO,BAR,ATCHAR,BAR,ATOM,BAR,
     +                                     NHO,BAR,COEFF,BAR,NBOCHAR,
     +                                     BAR,WEIGHT,BAR,LOCVAL
                       FIRST = .FALSE.
                   ELSE
                       WRITE (UNITID,9050) BAR,ATCHAR,BAR,ATOM,BAR,
     +                                     NHO,BAR,COEFF,BAR,BAR,BAR
                   END IF
  500           CONTINUE

                WRITE (UNITID,9010) (DASH,I=1,90)
 5000        CONTINUE
         END IF
C
C
C             ...print out Rydberg NBO information.
C
C
         IF (NYB.GT.0) THEN
             NBEG = NEND + 1
             NEND = NEND + NYB
             NBOCHAR = NBOSYMB (6)

             BDIDX = NBOBD (NBEG)
             ATOLD = BDCEN (1,BDIDX)

             DO 6000 N = NBEG,NEND
                NBO = NBO + 1
                BDIDX = NBOBD (N)
                WEIGHT = W (NBO)
                REST = REST + WEIGHT
                NHO = BDBAS (1,BDIDX)
                ATOM = BDCEN (1,BDIDX)
                COEFF = B (1,N)
                LOCVAL = LOCAL (ATOM,NBO)

                ZVAL = ZATOM (ATOM)
                IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                    ATCHAR = ATSYMB (104)
                ELSE
                    ATCHAR = ATSYMB (ZVAL)
                END IF

                IF (ATOM.NE.ATOLD) THEN
                    ATOLD = ATOM
                    WRITE (UNITID,9010) (DASH,I=1,90)
                END IF

                WRITE (UNITID,9040) NBO,BAR,ATCHAR,BAR,ATOM,BAR,NHO,
     +                              BAR,COEFF,BAR,NBOCHAR,BAR,WEIGHT,
     +                              BAR,LOCVAL

 6000        CONTINUE
             WRITE (UNITID,9010) (DASH,I=1,90)
         END IF
C
C
C             ...print out Rydberg NAO information.
C
C
         IF (NRA.GT.0) THEN
             NBOCHAR = NBOSYMB (7)

             DO 7000 N = 1,NRA
                OFF = ATROFF (N)
                NRBA = ATNRB (N)
                ATOM = ATRIDX (N)

                ZVAL = ZATOM (ATOM)
                IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                    ATCHAR = ATSYMB (104)
                ELSE
                    ATCHAR = ATSYMB (ZVAL)
                END IF

                FIRST = .TRUE.
                DO 700 I = 1,NRBA
                   RYD = OFF + I
                   NBO = NCB + NHB + RYD
                   WEIGHT = W (NBO)
                   LOCVAL = LOCAL (ATOM,NBO)
                   REST = REST + WEIGHT

                   LTYPE = RSHELL (RYD)
                   IF (LTYPE.GT.6) THEN
                       LCHAR = LSYMB (7)
                   ELSE
                       LCHAR = LSYMB (LTYPE)
                   END IF

                   IF (FIRST) THEN
                       WRITE (UNITID,9020) NBO,BAR,ATCHAR,BAR,ATOM,BAR,
     +                                     BAR,LCHAR,BAR,NBOCHAR,BAR,
     +                                     WEIGHT,BAR,LOCVAL
                       FIRST = .FALSE.
                   ELSE
                       WRITE (UNITID,9030) NBO,BAR,BAR,BAR,BAR,LCHAR,
     +                                     BAR,NBOCHAR,BAR,WEIGHT,BAR,
     +                                     LOCVAL
                   END IF
  700           CONTINUE

                WRITE (UNITID,9010) (DASH,I=1,90)
 7000        CONTINUE
         END IF

CCC hughes

C         call putrec(20,'JOBARC','ZIMATRIX',natom*nbas,zimtx)

CCC hughes

C
C
C             ...formats for printing the NBO information.
C
C
 9000    FORMAT (1X,A5,1X,     A1, 1X,A4,1X,    A1, 1X,A6,1X,  A1,
     +           1X,A5,1X,     A1, 1X,A10,1X,   A1, 4X,A8,4X,  A1,
     +           3X,A9,3X,     A1, 2X,A8,2X)
 9010    FORMAT (1X,90A1)
 9020    FORMAT (2X,I3,2X,     A1, 2X,A2,2X,    A1, 2X,I3,3X,  A1,
     +           7X,           A1, 5X,A3,4X,    A1, 1X,A14,1X, A1,
     +           2X,F12.10,1X, A1, 2X,F8.4,2X)
 9030    FORMAT (2X,I3,2X,     A1, 2X,2X,2X,    A1, 2X,3X,3X,  A1,
     +           7X,           A1, 5X,A3,4X,    A1, 1X,A14,1X, A1,
     +           2X,F12.10,1X, A1, 2X,F8.4,2X)
 9040    FORMAT (2X,I3,2X,     A1, 2X,A2,2X,    A1, 2X,I3,3X,  A1,
     +           2X,I3,2X,     A1, 1X,F10.6,1X, A1, 1X,A14,1X, A1,
     +           2X,F12.10,1X, A1, 2X,F8.4,2X)
 9050    FORMAT (2X,3X,2X,     A1, 2X,A2,2X,    A1, 2X,I3,3X,  A1,
     +           2X,I3,2X,     A1, 1X,F10.6,1X, A1, 1X,14X,1X, A1,
     +           2X,12X,1X,    A1)
C
C
C             ...print out population chart.
C
C
         WRITE (UNITID,9060) 'Population Chart'
 9060    FORMAT (//,20X,A16,//)

         WRITE (UNITID,9070) 'NBO-type',BAR,'Occupancy'
         WRITE (UNITID,9080) (DASH,I=1,28)
         WRITE (UNITID,9090) 'Lewis',BAR,LEWIS
         WRITE (UNITID,9100) 'Delocalized',BAR,DELOC
         WRITE (UNITID,9110) 'Rest',BAR,REST
         WRITE (UNITID,9080) (DASH,I=1,28)
         WRITE (UNITID,9120) 'Total',BAR,LEWIS+DELOC+REST

 9070    FORMAT (3X,A8,2X,  A1, 4X,A9)
 9080    FORMAT (2X,28A1)
 9090    FORMAT (7X,A5,1X,  A1, 1X,F12.8)
 9100    FORMAT (1X,A11,1X, A1, 1X,F12.8)
 9110    FORMAT (8X,A4,1X,  A1, 1X,F12.8)
 9120    FORMAT (7X,A5,1X,  A1, 1X,F12.8)
C
C
C             ...print out the NBO atomic localization map.
C
C
         WRITE (UNITID,9130) 'NBO Atomic Localization Map'
 9130    FORMAT (//,20X,A23)

         CALL    MAT__PRINT_A_FLOAT_3_NOZEROS
     +
     +                ( UNITID,
     +                  ' Rows = Atoms ; Columns = NBOs ',
     +                  NATOM,NBAS,
     +                  NATOM,NBAS,
     +                  LOCAL )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
