      SUBROUTINE INTEXP(PSO,PAO,ISOR,ISOS,ISOT,
     *                  NORBA,NORBB,NORBC,NORBD,NOABCD,
     *                  ICGTOA,ICGTOB,ICGTOC,ICGTOD,
     *                  DIAGAB,DIAGCD,SHABAB,SHAEQB,SHCEQD,
     *                  TCONAB,TCONCD,
     *                  KHKTA,KHKTB,KHKTC,KHKTD,
     *                  NHKTA,NHKTB,NHKTC,NHKTD,
     *                  THRESH,SYMFAC,IPRINT,BIGVEC,NOPMAT,NODV,
     *                  WORK1,LWORK1,ENEREE,SCF_GRADEE,GRADEE,
     *                  SCF_HESSEE,HESSEE,NCOORD,
     *                  DHFA,DHFB,DRELA,DRELB,DNHF,DREF,NBASIS,
     *                  NIR,PCHF)
C
C     Calculates expectation values of differentiated integrals
C
CEND
C                                      TUH 210488
C
C  ADAPTED TO THE ACES II ENVIRONMENT AND EXTENDED FOR UHF
C  OCT/90/JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
      LOGICAL TWOCEN, THRCEN, FOUCEN, DERONE, DERTWO
      LOGICAL NODER, DIAGAB, DIAGCD, SHABAB, SHAEQB, SHCEQD,
     *        BIGVEC, DCMPAB, DCMPCD, TCONAB, TCONCD, NOPMAT, NODV
      LOGICAL SCF,NONHF
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
#include "mxcent.par"
      DIMENSION DERIV(54), SCF_DERIV(54),
     *          ICGTOA(MXAOVC), ICGTOB(MXAOVC),
     *          ICGTOC(MXAOVC), ICGTOD(MXAOVC)
      DIMENSION PSO(1), PAO(1)
      DIMENSION GRADEE(NCOORD),HESSEE(NCOORD,NCOORD),
     &          SCF_GRADEE(NCOORD),SCF_HESSEE(NCOORD,NCOORD)
C
      COMMON /INTADR/ IWKAO, IWKSO, IWKHHS, IWK1HH, IWK1HC, IWKLST
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
      DIMENSION WORK1(LWORK1)
C
C  NEW COMMON BLOCK FOR HOLDING DENSITY MATRICES IN THE ACES II
C  ENVIRONMENT
C 
C    DHFA  .... ALPHA HF DENSITY
C    DHFB  .... BETA HF DENSITY
C    DRELA .... ALPHA RELAXED DENSITY
C    DRELB .... BETA RELAXED DENSITY
C    DNHF ..... TOTAL NON-DIRECT RELAXED DENSITY IN NON-HF METHODS
C    DREF ..... DENSITY OF THE REFERENCE FUNCTIONS IN NON HF METHODS
C
C  FOR DETAILS WITH RESPECT TO THE NON-HF METHODS SEE
C
C    J. GAUSS, J.F. STANTON, R.J. BARTLETT, J. CHEM. PHYS., SUBMITTED
C 
C    AND J. GAUSS, UNPUBLISHED NOTES
C
      DIMENSION DHFA(NBASIS,NBASIS,0:NIR),DHFB(NBASIS,NBASIS,0:NIR),
     &      DRELA(NBASIS,NBASIS,0:NIR),DRELB(NBASIS,NBASIS,0:NIR),
     &      DNHF(NBASIS,NBASIS,0:NIR),DREF(NBASIS,NBASIS,0:NIR)
C
C  NEW COMMON BLOCK IN THE ACES II ENVIRONMENT FOR KEEPINT TRACK ABOUT
C  THE METHOD 
C
      COMMON/METHOD/IUHF,SCF,NONHF
C
      COMMON /EXPCOM/ SIGN1X, SIGN1Y, SIGN1Z, SIGN2X, SIGN2Y, SIGN2Z,
     *                SIGN3X, SIGN3Y, SIGN3Z, SIGN4X, SIGN4Y, SIGN4Z,
     *                NCENT1, NCENT2, NCENT3, NCENT4,
     *                ISO1,   ISO2,   ISO3,   ISO4,
     *                DERONE, DERTWO, TWOCEN, THRCEN, FOUCEN,
     *                NINTYP, NCCINT
C
      COMMON /FLAGS2/ IFLAGS2(500)
C
      DATA FOURTH,HALF,ONE,TWO,FOUR,ZERADD,AZERO /0.25D0,0.5D0,1.0D0,
     &     2.0D0,4.0D0,1.D-16,0.D0/
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTXOR(I,J) = XOR(I,J)
C
      ITRI(I,J) = MAX(I,J)*(MAX(I,J) - 1)/2 + MIN(I,J)
      INDEX(I,J,K,L) = NORBA*(NORBB*(NORBC*(L-1)+K-1)+J-1)+I
C
c      IF (IPRINT .GT. 9) THEN
c       CALL HEADER('Output from INTEXP',-1)
c       WRITE (LUPRI, '(A,2L5)') ' NODV, NOPMAT', NODV, NOPMAT
c       WRITE (LUPRI, '(A,3I5)') ' ISOR/S/T  ', ISOR,ISOS,ISOT
c       WRITE (LUPRI, '(A,3F12.6)') 'SIGN1',SIGN1X,SIGN1Y,SIGN1Z
c       WRITE (LUPRI, '(A,3F12.6)') 'SIGN2',SIGN2X,SIGN2Y,SIGN2Z
c       WRITE (LUPRI, '(A,3F12.6)') 'SIGN3',SIGN3X,SIGN3Y,SIGN3Z
c       WRITE (LUPRI, '(A,3F12.6)') 'SIGN4',SIGN4X,SIGN4Y,SIGN4Z
c       WRITE (LUPRI, '(A,4I5)') ' NORB ', NORBA,NORBB,NORBC,NORBD
c       WRITE (LUPRI, '(A,4I5)') ' NCENT ',
c     *                              NCENT1,NCENT2,NCENT3,NCENT4
c       WRITE (LUPRI, '(A,2L5)') ' DIAGAB/CD ', DIAGAB,DIAGCD
c       WRITE (LUPRI, '(A,3L5)') ' SHAEQB, SHCEQD, SHABAB ',
c       WRITE (LUPRI, '(A,4I5)') ' KHKTA', KHKTA,KHKTB,KHKTC,KHKTD
c       WRITE (LUPRI, '(A,4I5)') ' NHKTA', NHKTA,NHKTB,NHKTC,NHKTD
c       WRITE (LUPRI, '(A,F12.6)') ' THRESH ', THRESH
c       WRITE (LUPRI, '(A,F12.6)') ' SYMFAC ', SYMFAC
c       WRITE (LUPRI, '(A,L5)')    ' BIGVEC ', BIGVEC
c       WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals A '
c       WRITE (LUPRI, '(20I5)') (ICGTOA(I),I=1, NORBA)
c       WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals B '
c       WRITE (LUPRI, '(20I5)') (ICGTOB(I),I=1, NORBB)
c       WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals C '
c       WRITE (LUPRI, '(20I5)') (ICGTOC(I),I=1, NORBC)
c       WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals D '
c       WRITE (LUPRI, '(20I5)') (ICGTOD(I),I=1, NORBD)
c      ENDIF
C
      ISOTS=IBTXOR(ISOT,ISOS)
      ISORT=IBTXOR(ISOR,ISOT)
      ISORTS=IBTXOR(ISOR,ISOTS)
C
      SFAC=SYMFAC
      IF(.NOT.SHABAB) SFAC=TWO*SFAC
      IF(.NOT.SHAEQB) SFAC=TWO*SFAC
      IF(.NOT.SHCEQD) SFAC=TWO*SFAC
C
      THRSH=MAX(THRESH,ZERADD)

      THRSH = 1.0D-18
      NODER=.TRUE.
      CALL ZERO(DERIV,NINTYP)
      CALL ZERO(SCF_DERIV,NINTYP)
CSS      Write(6,*)
CSS      Write(6, "(6(1x,F12.6))") (WORK1(IWKAO+i),i=1,2250)
C
C FIRST SET UP LOOP FOR RHF
C 
      IF(IUHF.EQ.0) THEN
C
C   SCF CASE
C
      IF(SCF) THEN
C
C     ***** Loop over shell components *****
C
      IAOFF=IWKAO+1
      DO 100 ICOMPA=1,KHKTA
       KHKTBB=KHKTB
       IF(DIAGAB) KHKTBB=ICOMPA
       DO 200 ICOMPB = 1,KHKTBB
        DCMPAB=DIAGAB.AND.ICOMPA.EQ.ICOMPB
        FACAB=PT(IBTAND(ITYPE(NHKTB,ICOMPB),ISOR))
        IF(DIAGAB.AND.ICOMPA.NE.ICOMPB) FACAB=TWO*FACAB
        DO 300 ICOMPC=1,KHKTC
         PTC=PT(IBTAND(ITYPE(NHKTC,ICOMPC),ISOT))
         KHKTDD=KHKTD
         IF(DIAGCD) KHKTDD=ICOMPC
         DO 400 ICOMPD=1,KHKTDD
          DCMPCD=DIAGCD.AND.ICOMPC.EQ.ICOMPD
          FACCD=PTC*PT(IBTAND(ITYPE(NHKTD,ICOMPD),ISOTS))
          IF(DIAGCD.AND.ICOMPC.NE.ICOMPD) FACCD=TWO*FACCD
          FCABAB=FACAB*FACCD*SFAC
          IAOADR=IAOFF
C
C  Loop over contracted functions *****
C
          DO 500 IORBA=1,NORBA
           INDXA=ICGTOA(IORBA)+IORBA+(ICOMPA-1)*NORBA
           NORBBB=NORBB
           IF (TCONAB) NORBBB=IORBA
           DO 600 IORBB=1,NORBBB
            INDXB=ICGTOB(IORBB)+IORBB+(ICOMPB-1)*NORBB
            DHFAB=DHFA(INDXA,INDXB,ISOR)
            FAB=ONE
            IF(TCONAB.AND.IORBA.NE.IORBB) FAB=TWO
            DO 700 IORBC=1,NORBC
             INDXC=ICGTOC(IORBC)+IORBC+(ICOMPC-1)*NORBC
             DHFAC=DHFA(INDXA,INDXC,ISOT)
             DHFBC=DHFA(INDXB,INDXC,ISORT)
             NORBDD=NORBD
             IF(TCONCD) NORBDD=IORBC
             DO 800 IORBD=1,NORBDD
              INDXD=ICGTOD(IORBD)+IORBD+(ICOMPD-1)*NORBD
              DHFAD=DHFA(INDXA,INDXD,ISOTS)
              DHFBD=DHFA(INDXB,INDXD,ISORTS)
              DHFCD=DHFA(INDXC,INDXD,ISOS)
              FCD=ONE
              IF(TCONCD.AND.IORBC.NE.IORBD) FCD=TWO
C
C              IF (IFLAGS2(153).EQ.0) THEN
C                 PVAL=HALF*FAB*FCD*FCABAB
C     &                *( DHFAB*DHFCD
C     &                  -FOURTH*(DHFAC*DHFBD+DHFAD*DHFBC))
C              ELSE
C                 PVAL=HALF*FAB*FCD*FCABAB*(DHFAB*DHFCD)
C              ENDIF
C
              PVAL=HALF*FAB*FCD*FCABAB
     &             *( DHFAB*DHFCD
     &               -PCHF*FOURTH*(DHFAC*DHFBD+DHFAD*DHFBC))
C
              IF (ABS(PVAL) .GT. THRSH) THEN
               NODER = .FALSE.
               INT   = IAOADR
CSSS               write(*,*)'a,b,c,d',indxa,indxb,indxc,indxd
CSSS               write(*,*) TCONAB,TCONCD,NORBA,NORBB,NORBC,NORBD
#ifdef _DEBUG
               if (int+NCCINT*NINTYP.gt.lwork1) then
                  print *, '@INTEXP: Out-of-bounds read.'
                  call errex
               end if
#endif
               DO 900 I = 1, NINTYP 
CSSS                write(6, "(1x,2F20.10)") pval,work1(int)
                DERIV(I) = DERIV(I)+PVAL*WORK1(INT)
                INT = INT + NCCINT
900            CONTINUE
              ENDIF
              IAOADR=IAOADR+1
800          CONTINUE
700         CONTINUE
600        CONTINUE
500       CONTINUE
          IAOFF=IAOFF+NOABCD
400      CONTINUE
300     CONTINUE
200    CONTINUE
100   CONTINUE
C
      ELSE
C   
C  CORRELATED METHODS
C
C     ***** Loop over shell components *****
C
      IAOFF=IWKAO+1
      DO 1100 ICOMPA=1,KHKTA
       KHKTBB=KHKTB
       IF(DIAGAB) KHKTBB=ICOMPA
       DO 1200 ICOMPB = 1,KHKTBB
        DCMPAB=DIAGAB.AND.ICOMPA.EQ.ICOMPB
        FACAB=PT(IBTAND(ITYPE(NHKTB,ICOMPB),ISOR))
        IF(DIAGAB.AND.ICOMPA.NE.ICOMPB) FACAB=TWO*FACAB
        DO 1300 ICOMPC=1,KHKTC
         PTC=PT(IBTAND(ITYPE(NHKTC,ICOMPC),ISOT))
         KHKTDD=KHKTD
         IF(DIAGCD) KHKTDD=ICOMPC
         DO 1400 ICOMPD=1,KHKTDD
          DCMPCD=DIAGCD.AND.ICOMPC.EQ.ICOMPD
          FACCD=PTC*PT(IBTAND(ITYPE(NHKTD,ICOMPD),ISOTS))
          IF(DIAGCD.AND.ICOMPC.NE.ICOMPD) FACCD=TWO*FACCD
          FCABAB=FACAB*FACCD*SFAC
C
C  Transform P-matrix block from SO basis to AO basis
C
          CALL PBLOCK(PSO,PAO,ICOMPA,ICOMPB,ICOMPC,ICOMPD,
     *                ISOR,ISOS,ISOT)
C          IF(IPRINT.GT.10) THEN
C           Write "Print the Shell INFo"
C           WRITE (LUPRI,'(A,I10)') ' IAOFF ', IAOFF
C           WRITE (LUPRI,'(A,4I5)') ' ICOMP ',
C     *             ICOMPA,ICOMPB,ICOMPC,ICOMPD
C           WRITE (LUPRI,'(A,2L5)') ' DCMPAB/CD ',DCMPAB,DCMPCD
C           WRITE (LUPRI,'(A,F12.6)') ' FCABAB ', FCABAB
C          ENDIF
          IAOADR=IAOFF
C
C  Loop over contracted functions *****
C
          DO 1500 IORBA=1,NORBA
           INDXA=ICGTOA(IORBA)+IORBA+(ICOMPA-1)*NORBA
           NORBBB=NORBB
           IF (TCONAB) NORBBB=IORBA
           DO 1600 IORBB=1,NORBBB
            INDXB=ICGTOB(IORBB)+IORBB+(ICOMPB-1)*NORBB
            DHFAB=DHFA(INDXA,INDXB,ISOR)
            DRELAB=DRELA(INDXA,INDXB,ISOR)
            FAB=ONE
            IF(TCONAB.AND.IORBA.NE.IORBB) FAB=TWO
            DO 1700 IORBC=1,NORBC
             INDXC=ICGTOC(IORBC)+IORBC+(ICOMPC-1)*NORBC
             DHFAC=DHFA(INDXA,INDXC,ISOT)
             DHFBC=DHFA(INDXB,INDXC,ISORT)
             DRELAC=DRELA(INDXA,INDXC,ISOT)
             DRELBC=DRELA(INDXB,INDXC,ISORT)
             NORBDD=NORBD
             IF(TCONCD) NORBDD=IORBC
             DO 1800 IORBD=1,NORBDD
              INDXD=ICGTOD(IORBD)+IORBD+(ICOMPD-1)*NORBD
              DHFAD=DHFA(INDXA,INDXD,ISOTS)
              DHFBD=DHFA(INDXB,INDXD,ISORTS)
              DHFCD=DHFA(INDXC,INDXD,ISOS)
              DRELAD=DRELA(INDXA,INDXD,ISOTS)
              DRELBD=DRELA(INDXB,INDXD,ISORTS)
              DRELCD=DRELA(INDXC,INDXD,ISOS)
              FCD=ONE
              IF(TCONCD.AND.IORBC.NE.IORBD) FCD=TWO
              PAOVAL = PAO(INDEX(IORBA,IORBB,IORBC,IORBD))
CSSS              PAOVAL = 0.0D0
C
C The HF-SCF gradients are apparently needed for a extrapolation
C scheme in rate constant work. Luckily we can compute it along
C with the total (correlated) gradeint. Ajith Perera, 01/2006.
C
              SCF_PVAL=HALF*FAB*FCD*FCABAB
     &                 *( DHFAB*DHFCD
     &                 -PCHF*FOURTH*(DHFAC*DHFBD+DHFAD*DHFBC))
C
              PVAL=HALF*FAB*FCD*FCABAB*
     *             ( DHFAB*(DHFCD+DRELCD)+DRELAB*DHFCD
     *               -FOURTH*( DHFAC*(DHFBD+DRELBD)
     *                        +DHFAD*(DHFBC+DRELBC)
     *                        +DRELAC*DHFBD
     *                        +DRELAD*DHFBC)
     *               +FOUR*PAOVAL)
CSSS                write(6, "(1x,4(1x,F20.10))")DRELAD, DRELBD, DRELCD,PAOVAL
C
C Let's use the only PVAL to check against the threshold and ignore
C SCF_PVAL
C
              IF (ABS(PVAL) .GT. THRSH) THEN
               NODER = .FALSE.
               INT   = IAOADR
#ifdef _DEBUG
               if (int+NCCINT*NINTYP.gt.lwork1) then
                  print *, '@INTEXP: Out-of-bounds read.'
                  call errex
               end if
#endif
               DO 1900 I = 1, NINTYP
                SCF_DERIV(I) = SCF_DERIV(I) + SCF_PVAL*WORK1(INT)
                DERIV(I) = DERIV(I)+PVAL*WORK1(INT)
CSSS                write(6, "(1x,2F20.10,1x,i4)") pval,work1(int),nccint
                INT = INT + NCCINT
1900           CONTINUE
              ENDIF
              IAOADR=IAOADR+1
1800         CONTINUE
1700        CONTINUE
1600       CONTINUE
1500      CONTINUE
          IAOFF=IAOFF+NOABCD
1400     CONTINUE
1300    CONTINUE
1200   CONTINUE
1100  CONTINUE
      ENDIF
C
CSSS      Write(6,*) DERIV(1), DERIV(2), DERIV(3) 
      ELSE
C
C  LOOP FOR UHF
C
C  FIRST SCF
C
      IF(SCF) THEN
C
C     ***** Loop over shell components *****
C
      IAOFF=IWKAO+1
      DO 2100 ICOMPA=1,KHKTA
       KHKTBB=KHKTB
       IF(DIAGAB) KHKTBB=ICOMPA
       DO 2200 ICOMPB = 1,KHKTBB
        DCMPAB=DIAGAB.AND.ICOMPA.EQ.ICOMPB
        FACAB=PT(IBTAND(ITYPE(NHKTB,ICOMPB),ISOR))
        IF(DIAGAB.AND.ICOMPA.NE.ICOMPB) FACAB=TWO*FACAB
        DO 2300 ICOMPC=1,KHKTC
         PTC=PT(IBTAND(ITYPE(NHKTC,ICOMPC),ISOT))
         KHKTDD=KHKTD
         IF(DIAGCD) KHKTDD=ICOMPC
         DO 2400 ICOMPD=1,KHKTDD
          DCMPCD=DIAGCD.AND.ICOMPC.EQ.ICOMPD
          FACCD=PTC*PT(IBTAND(ITYPE(NHKTD,ICOMPD),ISOTS))
          IF(DIAGCD.AND.ICOMPC.NE.ICOMPD) FACCD=TWO*FACCD
          FCABAB=FACAB*FACCD*SFAC
          IAOADR=IAOFF
C
C  Loop over contracted functions *****
C
          DO 2500 IORBA=1,NORBA
           INDXA=ICGTOA(IORBA)+IORBA+(ICOMPA-1)*NORBA
           NORBBB=NORBB
           IF (TCONAB) NORBBB=IORBA
           DO 2600 IORBB=1,NORBBB
            INDXB=ICGTOB(IORBB)+IORBB+(ICOMPB-1)*NORBB
            DHFAB=DHFA(INDXA,INDXB,ISOR)+DHFB(INDXA,INDXB,ISOR)
            FAB=ONE
            IF(TCONAB.AND.IORBA.NE.IORBB) FAB=TWO
            DO 2700 IORBC=1,NORBC
             INDXC=ICGTOC(IORBC)+IORBC+(ICOMPC-1)*NORBC
             DHFAAC=DHFA(INDXA,INDXC,ISOT)
             DHFBAC=DHFB(INDXA,INDXC,ISOT)
             DHFABC=DHFA(INDXB,INDXC,ISORT)
             DHFBBC=DHFB(INDXB,INDXC,ISORT)
             NORBDD=NORBD
             IF(TCONCD) NORBDD=IORBC
             DO 2800 IORBD=1,NORBDD
              INDXD=ICGTOD(IORBD)+IORBD+(ICOMPD-1)*NORBD
              DHFAAD=DHFA(INDXA,INDXD,ISOTS)
              DHFBAD=DHFB(INDXA,INDXD,ISOTS)
              DHFABD=DHFA(INDXB,INDXD,ISORTS)
              DHFBBD=DHFB(INDXB,INDXD,ISORTS)
              DHFCD=DHFA(INDXC,INDXD,ISOS)+DHFB(INDXC,INDXD,ISOS)
              FCD=ONE
              IF(TCONCD.AND.IORBC.NE.IORBD) FCD=TWO
C
C              IF (IFLAGS2(153).EQ.0) THEN
C                 PVAL=HALF*FAB*FCD*FCABAB
C     *                *( DHFAB*DHFCD
C     *                  -HALF*( DHFAAC*DHFABD
C     *                         +DHFAAD*DHFABC
C     *                         +DHFBAC*DHFBBD
C     *                         +DHFBAD*DHFBBC))
C              ELSE
C                 PVAL=HALF*FAB*FCD*FCABAB*(DHFAB*DHFCD)
C              ENDIF
C
              PVAL=HALF*FAB*FCD*FCABAB
     *             *( DHFAB*DHFCD
     *               -PCHF*HALF*( DHFAAC*DHFABD
     *                           +DHFAAD*DHFABC
     *                           +DHFBAC*DHFBBD
     *                           +DHFBAD*DHFBBC))
C
              IF (ABS(PVAL) .GT. THRSH) THEN
               NODER = .FALSE.
               INT   = IAOADR
#ifdef _DEBUG
               if (int+NCCINT*NINTYP.gt.lwork1) then
                  print *, '@INTEXP: Out-of-bounds read.'
                  call errex
               end if
#endif
               DO 2900 I = 1, NINTYP
                DERIV(I) = DERIV(I)+PVAL*WORK1(INT)
                INT = INT + NCCINT
2900           CONTINUE
              ENDIF
              IAOADR=IAOADR+1
2800         CONTINUE
2700        CONTINUE
2600       CONTINUE
2500      CONTINUE
          IAOFF=IAOFF+NOABCD
2400     CONTINUE
2300    CONTINUE
2200   CONTINUE
2100  CONTINUE
C
      ELSE
C
C CORRELATED METHODS HERE
C     
      IF(.NOT.NONHF) THEN
C
C   HF REFERENCE STATES
C
C     ***** Loop over shell components *****
C
      IAOFF=IWKAO+1
      DO 3100 ICOMPA=1,KHKTA
       KHKTBB=KHKTB
       IF(DIAGAB) KHKTBB=ICOMPA
       DO 3200 ICOMPB = 1,KHKTBB
        DCMPAB=DIAGAB.AND.ICOMPA.EQ.ICOMPB
        FACAB=PT(IBTAND(ITYPE(NHKTB,ICOMPB),ISOR))
        IF(DIAGAB.AND.ICOMPA.NE.ICOMPB) FACAB=TWO*FACAB
        DO 3300 ICOMPC=1,KHKTC
         PTC=PT(IBTAND(ITYPE(NHKTC,ICOMPC),ISOT))
         KHKTDD=KHKTD
         IF(DIAGCD) KHKTDD=ICOMPC
         DO 3400 ICOMPD=1,KHKTDD
          DCMPCD=DIAGCD.AND.ICOMPC.EQ.ICOMPD
          FACCD=PTC*PT(IBTAND(ITYPE(NHKTD,ICOMPD),ISOTS))
          IF(DIAGCD.AND.ICOMPC.NE.ICOMPD) FACCD=TWO*FACCD
          FCABAB=FACAB*FACCD*SFAC
C
C  Transform P-matrix block from SO basis to AO basis
C
          CALL PBLOCK(PSO,PAO,ICOMPA,ICOMPB,ICOMPC,ICOMPD,
     *                ISOR,ISOS,ISOT)
c          IF(IPRINT.GT.10) THEN
c           WRITE (LUPRI,'(A,I10)') ' IAOFF ', IAOFF
c           WRITE (LUPRI,'(A,4I5)') ' ICOMP ',
c     *             ICOMPA,ICOMPB,ICOMPC,ICOMPD
c           WRITE (LUPRI,'(A,2L5)') ' DCMPAB/CD ',DCMPAB,DCMPCD
c           WRITE (LUPRI,'(A,F12.6)') ' FCABAB ', FCABAB
c          ENDIF
          IAOADR=IAOFF
C
C  Loop over contracted functions *****
C
          DO 3500 IORBA=1,NORBA
           INDXA=ICGTOA(IORBA)+IORBA+(ICOMPA-1)*NORBA
           NORBBB=NORBB
           IF (TCONAB) NORBBB=IORBA
           DO 3600 IORBB=1,NORBBB
            INDXB=ICGTOB(IORBB)+IORBB+(ICOMPB-1)*NORBB
            DHFAB=DHFA(INDXA,INDXB,ISOR)+DHFB(INDXA,INDXB,ISOR)
            DRELAB=DRELA(INDXA,INDXB,ISOR)+DRELB(INDXA,INDXB,ISOR)
            FAB=ONE
            IF(TCONAB.AND.IORBA.NE.IORBB) FAB=TWO
            DO 3700 IORBC=1,NORBC
             INDXC=ICGTOC(IORBC)+IORBC+(ICOMPC-1)*NORBC
             DHFAAC=DHFA(INDXA,INDXC,ISOT)
             DHFBAC=DHFB(INDXA,INDXC,ISOT)
             DHFABC=DHFA(INDXB,INDXC,ISORT)
             DHFBBC=DHFB(INDXB,INDXC,ISORT)
             DRELAAC=DRELA(INDXA,INDXC,ISOT)
             DRELBAC=DRELB(INDXA,INDXC,ISOT)
             DRELABC=DRELA(INDXB,INDXC,ISORT)
             DRELBBC=DRELB(INDXB,INDXC,ISORT)
             NORBDD=NORBD
             IF(TCONCD) NORBDD=IORBC
             DO 3800 IORBD=1,NORBDD
              INDXD=ICGTOD(IORBD)+IORBD+(ICOMPD-1)*NORBD
              DHFAAD=DHFA(INDXA,INDXD,ISOTS)
              DHFBAD=DHFB(INDXA,INDXD,ISOTS)
              DHFABD=DHFA(INDXB,INDXD,ISORTS)
              DHFBBD=DHFB(INDXB,INDXD,ISORTS)
              DHFCD=DHFA(INDXC,INDXD,ISOS)+DHFB(INDXC,INDXD,ISOS)
              DRELAAD=DRELA(INDXA,INDXD,ISOTS)
              DRELBAD=DRELB(INDXA,INDXD,ISOTS)
              DRELABD=DRELA(INDXB,INDXD,ISORTS)
              DRELBBD=DRELB(INDXB,INDXD,ISORTS)
              DRELCD=DRELA(INDXC,INDXD,ISOS)+DRELB(INDXC,INDXD,ISOS)
              FCD=ONE
              IF(TCONCD.AND.IORBC.NE.IORBD) FCD=TWO
              PAOVAL = PAO(INDEX(IORBA,IORBB,IORBC,IORBD))
C
C The HF-SCF gradients are apparently needed for an extrapolation
C scheme in rate constant work. Luckily we can compute it along
C with the total (correlated) gradient. Ajith Perera, 01/2006.
C
              SCF_PVAL=HALF*FAB*FCD*FCABAB
     &                           *(DHFAB*DHFCD
     &                 -PCHF*HALF*(DHFAAC*DHFABD
     &                             +DHFAAD*DHFABC
     &                             +DHFBAC*DHFBBD
     &                             +DHFBAD*DHFBBC))
C
              PVAL=HALF*FAB*FCD*FCABAB
     *             *( DHFAB*(DHFCD+DRELCD)
     *               +DRELAB*DHFCD
     *               -HALF*( DHFAAC*(DHFABD+DRELABD)+DRELAAC*DHFABD
     *                      +DHFAAD*(DHFABC+DRELABC)+DRELAAD*DHFABC
     *                      +DHFBAC*(DHFBBD+DRELBBD)+DRELBAC*DHFBBD
     *                      +DHFBAD*(DHFBBC+DRELBBC)+DRELBAD*DHFBBC)
     *               +FOUR*PAOVAL)
              IF (ABS(PVAL) .GT. THRSH) THEN
               NODER = .FALSE.
               INT   = IAOADR
#ifdef _DEBUG
               if (int+NCCINT*NINTYP.gt.lwork1) then
                  print *, '@INTEXP: Out-of-bounds read.'
                  call errex
               end if
#endif
               DO 3900 I = 1, NINTYP
            SCF_DERIV(I) = SCF_DERIV(I) + SCF_PVAL*WORK1(INT)
                DERIV(I) = DERIV(I)+PVAL*WORK1(INT)
CSSS                write(6, "(1x,2F20.10)") pval,work1(int)
                INT = INT + NCCINT
3900           CONTINUE
              ENDIF
              IAOADR=IAOADR+1
3800         CONTINUE
3700        CONTINUE
3600       CONTINUE
3500      CONTINUE
          IAOFF=IAOFF+NOABCD
3400     CONTINUE
3300    CONTINUE
3200   CONTINUE
3100  CONTINUE
C 
      ELSE
C
C  NON HF METHODS
C
C     ***** Loop over shell components *****
C
      IAOFF=IWKAO+1
      DO 4100 ICOMPA=1,KHKTA
       KHKTBB=KHKTB
       IF(DIAGAB) KHKTBB=ICOMPA
       DO 4200 ICOMPB = 1,KHKTBB
        DCMPAB=DIAGAB.AND.ICOMPA.EQ.ICOMPB
        FACAB=PT(IBTAND(ITYPE(NHKTB,ICOMPB),ISOR))
        IF(DIAGAB.AND.ICOMPA.NE.ICOMPB) FACAB=TWO*FACAB
        DO 4300 ICOMPC=1,KHKTC
         PTC=PT(IBTAND(ITYPE(NHKTC,ICOMPC),ISOT))
         KHKTDD=KHKTD
         IF(DIAGCD) KHKTDD=ICOMPC
         DO 4400 ICOMPD=1,KHKTDD
          DCMPCD=DIAGCD.AND.ICOMPC.EQ.ICOMPD
          FACCD=PTC*PT(IBTAND(ITYPE(NHKTD,ICOMPD),ISOTS))
          IF(DIAGCD.AND.ICOMPC.NE.ICOMPD) FACCD=TWO*FACCD
          FCABAB=FACAB*FACCD*SFAC
C
C  Transform P-matrix block from SO basis to AO basis
C
          CALL PBLOCK(PSO,PAO,ICOMPA,ICOMPB,ICOMPC,ICOMPD,
     *                ISOR,ISOS,ISOT)
c          IF(IPRINT.GT.10) THEN
c           WRITE (LUPRI,'(A,I10)') ' IAOFF ', IAOFF
c           WRITE (LUPRI,'(A,4I5)') ' ICOMP ',
c     *             ICOMPA,ICOMPB,ICOMPC,ICOMPD
c           WRITE (LUPRI,'(A,2L5)') ' DCMPAB/CD ',DCMPAB,DCMPCD
c           WRITE (LUPRI,'(A,F12.6)') ' FCABAB ', FCABAB
c          ENDIF
          IAOADR=IAOFF
C
C  Loop over contracted functions *****
C
          DO 4500 IORBA=1,NORBA
           INDXA=ICGTOA(IORBA)+IORBA+(ICOMPA-1)*NORBA
           NORBBB=NORBB
           IF (TCONAB) NORBBB=IORBA
           DO 4600 IORBB=1,NORBBB
            INDXB=ICGTOB(IORBB)+IORBB+(ICOMPB-1)*NORBB
            DHFAB=DHFA(INDXA,INDXB,ISOR)+DHFB(INDXA,INDXB,ISOR)
            DRELAB=DRELA(INDXA,INDXB,ISOR)+DRELB(INDXA,INDXB,ISOR)
            DNHFAB=DNHF(INDXA,INDXB,ISOR)
            DREFAB=DREF(INDXA,INDXB,ISOR)
            FAB=ONE
            IF(TCONAB.AND.IORBA.NE.IORBB) FAB=TWO
            DO 4700 IORBC=1,NORBC
             INDXC=ICGTOC(IORBC)+IORBC+(ICOMPC-1)*NORBC
             DHFAAC=DHFA(INDXA,INDXC,ISOT)
             DHFBAC=DHFB(INDXA,INDXC,ISOT)
             DHFABC=DHFA(INDXB,INDXC,ISORT)
             DHFBBC=DHFB(INDXB,INDXC,ISORT)
             DRELAAC=DRELA(INDXA,INDXC,ISOT)
             DRELBAC=DRELB(INDXA,INDXC,ISOT)
             DRELABC=DRELA(INDXB,INDXC,ISORT)
             DRELBBC=DRELB(INDXB,INDXC,ISORT)
             DNHFAC=DNHF(INDXA,INDXC,ISOT)
             DNHFBC=DNHF(INDXB,INDXC,ISORT)
             DREFAC=DREF(INDXA,INDXC,ISOT)
             DREFBC=DREF(INDXB,INDXC,ISORT)
             NORBDD=NORBD
             IF(TCONCD) NORBDD=IORBC
             DO 4800 IORBD=1,NORBDD
              INDXD=ICGTOD(IORBD)+IORBD+(ICOMPD-1)*NORBD
              DHFAAD=DHFA(INDXA,INDXD,ISOTS)
              DHFBAD=DHFB(INDXA,INDXD,ISOTS)
              DHFABD=DHFA(INDXB,INDXD,ISORTS)
              DHFBBD=DHFB(INDXB,INDXD,ISORTS)
              DHFCD=DHFA(INDXC,INDXD,ISOS)+DHFB(INDXC,INDXD,ISOS)
              DRELAAD=DRELA(INDXA,INDXD,ISOTS)
              DRELBAD=DRELB(INDXA,INDXD,ISOTS)
              DRELABD=DRELA(INDXB,INDXD,ISORTS)
              DRELBBD=DRELB(INDXB,INDXD,ISORTS)
              DRELCD=DRELA(INDXC,INDXD,ISOS)+DRELB(INDXC,INDXD,ISOS)
              DNHFAD=DNHF(INDXA,INDXD,ISOTS)
              DNHFBD=DNHF(INDXB,INDXD,ISORTS)
              DNHFCD=DNHF(INDXC,INDXD,ISOS)
              DREFAD=DREF(INDXA,INDXD,ISOTS)
              DREFBD=DREF(INDXB,INDXD,ISORTS)
              DREFCD=DREF(INDXC,INDXD,ISOS)
              FCD=ONE
              IF(TCONCD.AND.IORBC.NE.IORBD) FCD=TWO
              PAOVAL = PAO(INDEX(IORBA,IORBB,IORBC,IORBD))
C
C The HF-SCF gradients are apparently needed for a extrapolation
C scheme in rate constant work. Luckily we can compute it along
C with the total (correlated) gradeint. Ajith Perera, 01/2006
C
              SCF_PVAL=HALF*FAB*FCD*FCABAB
     &                 *( DHFAB*DHFCD
     &                 -PCHF*FOURTH*(DHFAC*DHFBD+DHFAD*DHFBC))
C
              PVAL=HALF*FAB*FCD*FCABAB
     *             *( DHFAB*(DHFCD+DRELCD)
     *               +DRELAB*DHFCD
     *               +DREFAB*DNHFCD
     *               +DREFCD*DNHFAB 
     *               -HALF*( DHFAAC*(DHFABD+DRELABD)+DRELAAC*DHFABD
     *                      +DHFAAD*(DHFABC+DRELABC)+DRELAAD*DHFABC
     *                      +DHFBAC*(DHFBBD+DRELBBD)+DRELBAC*DHFBBD
     *                      +DHFBAD*(DHFBBC+DRELBBC)+DRELBAD*DHFBBC)
     *               -FOURTH*( DREFAC*DNHFBD+DREFBD*DNHFAC
     *                        +DREFAD*DNHFBC+DREFBC*DNHFAD)
     *               +FOUR*PAOVAL)
              IF (ABS(PVAL) .GT. THRSH) THEN
               NODER = .FALSE.
               INT   = IAOADR
#ifdef _DEBUG
               if (int+NCCINT*NINTYP.gt.lwork1) then
                  print *, '@INTEXP: Out-of-bounds read.'
                  call errex
               end if
#endif
               DO 4900 I = 1, NINTYP
            SCF_DERIV(I) = SCF_DERIV(I) + SCF_PVAL*WORK1(INT)
                DERIV(I) = DERIV(I)+PVAL*WORK1(INT)
CSSS                write(6, "(1x,2F20.10)") pval,work1(int)
                INT = INT + NCCINT
4900           CONTINUE
              ENDIF
              IAOADR=IAOADR+1
4800         CONTINUE
4700        CONTINUE
4600       CONTINUE
4500      CONTINUE
          IAOFF=IAOFF+NOABCD
4400     CONTINUE
4300    CONTINUE
4200   CONTINUE
4100  CONTINUE
C
      ENDIF
      ENDIF  
      ENDIF
C
      IF(.NOT.NODER) THEN
C
C See above comments about the SCF gradients and rate constant work.
C The GRADEE and HESSEE (if analytic Hessian is computed)
C are the gradient and Hessian.
C
      IF (.NOT.SCF) THEN
         CALL DEROUT(NINTYP,NCENT1,NCENT2,NCENT3,NCENT4,
     *             ISO1,ISO2,ISO3,ISO4,
     *             SIGN1X,SIGN1Y,SIGN1Z,SIGN2X,SIGN2Y,SIGN2Z,
     *             SIGN3X,SIGN3Y,SIGN3Z,SIGN4X,SIGN4Y,SIGN4Z,
     *             DERONE,DERTWO,TWOCEN,THRCEN,FOUCEN,
     *             IPRINT,SCF_DERIV,BIGVEC,ENEREE,SCF_GRADEE,
     &             SCF_HESSEE,NCOORD)
      END IF
C
       CALL DEROUT(NINTYP,NCENT1,NCENT2,NCENT3,NCENT4,
     *             ISO1,ISO2,ISO3,ISO4,
     *             SIGN1X,SIGN1Y,SIGN1Z,SIGN2X,SIGN2Y,SIGN2Z,
     *             SIGN3X,SIGN3Y,SIGN3Z,SIGN4X,SIGN4Y,SIGN4Z,
     *             DERONE,DERTWO,TWOCEN,THRCEN,FOUCEN,
     *             IPRINT,DERIV,BIGVEC,ENEREE,GRADEE,HESSEE,
     &             NCOORD)
      ENDIF
      RETURN
      END
