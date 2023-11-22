         SUBROUTINE  NLO__PRINT_NHO_RESULTS
     +
     +                    ( UNITID,
     +                      NATOM,NBOND,
     +                      MXNHBA,
     +                      NHB,NHA,
     +                      ZATOM,
     +                      ATNHB,ATHIDX,ATHOFF,
     +                      HSHELL,
     +                      BDNCEN,BDCEN,BDBAS,BDOCC,
     +                      H )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__PRINT_NHO_RESULTS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine prints details of the final set of
C                natural hybrid orbitals (NHOs) to be used to establish
C                the natural bond orbitals (NBOs) to a file specified
C                by its unit identification number.
C
C                For each atomic center the following is printed:
C
C                   1) Expansion of NHOs in terms of the valence NAOs.
C
C                   2) NHO hybridization pattern.
C
C                  Input:
C
C                    UNITID       =  printout unit identification #
C                    NATOM        =  total # of atoms
C                    NBOND        =  total # of bonds
C                    MXNHBA       =  maximum # of Hybrid NAOs per atom.
C                    NHB          =  # of Hybrid NAOs.
C                    NHA          =  # of atoms on which Hybrid NAOs
C                                    have been found.
C                    ZATOM (I)    =  atomic number for I-th atom.
C                    ATNHB (A)    =  # of Hybrid NAOs on hybrid atom A.
C                    ATHIDX (A)   =  atomic index for hybrid atom A.
C                    ATHOFF (A)   =  index offset for Hybrid NAOs for
C                                    hybrid atom A. This index is equal
C                                    to the total number of Hybrid
C                                    NAOs on all hybrid atoms preceeding
C                                    hybrid atom A.
C                    HSHELL (I)   =  l-shell type for the I-th Hybrid
C                                    NAO.
C                    BDNCEN (I)   =  # of atomic centers for I-th bond.
C                    BDCEN (I,J)  =  I-th atomic center index for
C                                    J-th bond.
C                    BDBAS (I,J)  =  I-th global basis (NHO) index for
C                                    J-th bond.
C                    BDOCC (I)    =  # of occupied levels for I-th bond.
C                    H (I,J)      =  MXNHBA x NHB matrix containing the
C                                    atomic NHOs. I is the local atomic
C                                    index labeling the atomic Hybrid
C                                    NAOs from which the NHOs are
C                                    constructed. J is the global NHO
C                                    index running over all NHB NHOs,
C                                    with all NHOs belonging to a
C                                    specific atomic center being
C                                    grouped together.
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

         CHARACTER*2   ATSYMB  (1:104)
         CHARACTER*1   HYBRID  (1:50)

         INTEGER     ATIDX
         INTEGER     ATOM,BOND
         INTEGER     I,J
         INTEGER     LENGTH
         INTEGER     LTYPE
         INTEGER     MXNHBA
         INTEGER     NAO,NHO
         INTEGER     NATOM,NBOND
         INTEGER     NOCC,NCEN
         INTEGER     NHB,NHA,NHBA
         INTEGER     OFF
         INTEGER     UNITID
         INTEGER     ZVAL

         INTEGER     ATNHB  (1:NHA  )
         INTEGER     ATHIDX (1:NHA  )
         INTEGER     ATHOFF (1:NHA  )
         INTEGER     BDNCEN (1:NHB  )
         INTEGER     BDOCC  (1:NHB  )
         INTEGER     HSHELL (1:NHB  )
         INTEGER     ZATOM  (1:NATOM)

         INTEGER     BDBAS   (1:NHA,1:NHB)
         INTEGER     BDCEN   (1:NHA,1:NHB)

         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  W (0:7)

         DOUBLE PRECISION  H (1:MXNHBA,1:NHB)

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

         PARAMETER  (ZERO = 0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...print out title.
C
C
         WRITE (UNITID,8000) 'Natural hybrid orbital analysis'
 8000    FORMAT (//,20X,A31,//)
C
C
C             ...print out header line.
C
C
         WRITE (UNITID,9000) 'NHO #',BAR,'Atom',BAR,'Atom #',BAR,
     +                       'NAO #',BAR,'NAO coeffs',BAR,
     +                       'Hybridization'
         WRITE (UNITID,9010) (DASH,I=1,62)
C
C
C             ...print out atomic NHO information.
C
C
         DO 1000 ATOM = 1,NHA
            OFF = ATHOFF (ATOM)
            NHBA = ATNHB (ATOM)
            ATIDX = ATHIDX (ATOM)

            ZVAL = ZATOM (ATIDX)

            IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                ATCHAR = ATSYMB (104)
            ELSE
                ATCHAR = ATSYMB (ZVAL)
            END IF

            DO 100 J = 1,NHBA
               NHO = OFF + J

               CALL  MAT__W_EQ_ZERO_FLOAT (8,8,W)

               DO 110 I = 1,NHBA
                  NAO = OFF + I
                  LTYPE = MIN0 (HSHELL (NAO),7)
                  W (LTYPE) = W (LTYPE) + H (I,NHO) ** 2
  110          CONTINUE

               CALL  NLO__FORM_NHO_HYBRID_PATTERN
     +
     +                    ( W,
     +
     +                          LENGTH,
     +                          HYBRID )
     +
     +
               WRITE (UNITID,9020) NHO,BAR,ATCHAR,BAR,ATIDX,BAR,1,BAR,
     +                             H (1,NHO),BAR,(HYBRID(I),I=1,LENGTH)

               DO 120 I = 2,NHBA
                  WRITE (UNITID,9030) BAR,BAR,BAR,I,BAR,H (I,NHO),BAR
  120          CONTINUE

               WRITE (UNITID,9010) (DASH,I=1,62)

  100       CONTINUE

 1000    CONTINUE
C
C
C             ...formats for printing the atomic NHO information.
C
C
 9000    FORMAT (1X,A5,1X, A1, 1X,A4,1X,  A1, 1X,A6,1X, A1,
     +           1X,A5,1X, A1, 1X,A10,1X, A1,
     +           2X,A13)
 9010    FORMAT (1X,62A1)
 9020    FORMAT (2X,I3,2X, A1, 2X,A2,2X,    A1, 2X,I3,3X, A1,
     +           1X,I4,2X, A1, 1X,F10.6,1X, A1, 2X,50(A1))
 9030    FORMAT (2X,3X,2X, A1, 2X,2X,2X,    A1, 2X,3X,3X, A1,
     +           1X,I4,2X, A1, 1X,F10.6,1X, A1)
C
C
C             ...print out title.
C
C
         WRITE (UNITID,8000) '    Bond formation analysis    '
C
C
C             ...print out header line and bond formation
C                information # 1.
C
C
         WRITE (UNITID,9040) 'BOND #',BAR,'Occ #',BAR,
     +                       '# of Centers',BAR,'Atomic Centers  '
         WRITE (UNITID,9050) (DASH,I=1,62)

         DO 200 BOND = 1,NBOND
            NOCC = BDOCC (BOND)
            NCEN = BDNCEN (BOND)
            WRITE (UNITID,9060) BOND,BAR,NOCC,BAR,NCEN,BAR,
     +                          (BDCEN (I,BOND),I=1,NCEN)
  200    CONTINUE
         WRITE (UNITID,9050) (DASH,I=1,62)
C
C
C             ...print out header line and bond formation
C                information # 2.
C
C
         WRITE (UNITID,9040) 'BOND #',BAR,'Occ #',BAR,
     +                       '# of Centers',BAR,'Bond NHO indices'
         WRITE (UNITID,9050) (DASH,I=1,62)

         DO 210 BOND = 1,NBOND
            NOCC = BDOCC (BOND)
            NCEN = BDNCEN (BOND)
            WRITE (UNITID,9060) BOND,BAR,NOCC,BAR,NCEN,BAR,
     +                          (BDBAS (I,BOND),I=1,NCEN)
  210    CONTINUE
         WRITE (UNITID,9050) (DASH,I=1,62)
C
C
C             ...formats for printing the bond formation information.
C
C
 9040    FORMAT (1X,A6,1X, A1, 1X,A5,1X,  A1, 1X,A12,1X, A1, 2X,A16)
 9050    FORMAT (1X,62A1)
 9060    FORMAT (2X,I4,2X, A1, 2X,I3,2X,  A1, 5X,I4,5X, A1, 2X,80(I4))
C
C
C             ...ready!
C
C
         RETURN
         END
