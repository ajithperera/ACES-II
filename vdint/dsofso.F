      SUBROUTINE DSOFSO(CMO,FT,DV,DSO,DSOSCF,FSO,SCR,DENSCF,DOSCF)
C
C     THIS SUBROUTINE CALCULATES THE FOLDED TOTAL ONE-ELECTRON DENSITY 
C     AND FOCK MATRICES IN AO BASIS.  THE REQUIRED INPUT IS 
C 
C     SCF : MO COEFFICIENTS AND ORBITAL ENERGIES
C
C     CORREL. METHODS : RELAXED DENSITY MATRIX, INTERMEDIATES I(P,Q)
C                       AND MO COEFFICIENTS
C
C     CMO ..... ARRAY ALLOCATED TO HOLD THE EIGEN VECTORS 
C               SIZE : FULL NBASIS*NBASIS
C     FT  ..... ARRAY ALLOCATED TO HOLD THE INTERMEDIATES
C               SIZE : FULL NBASIS*NBASIS
C     DV  ..... ARRAY ALLOCATED TO HOLD THE TOTAL DENSITY
C               MATRIX , SIZE FULL NBASIS*NBASIS
C     DSO ..... RETURNS THE FOLDED TOTAL ONE-ELECTRON
C               DENSITY MATRIX IN SYMMETRY ADAPTED AO BASIS
C     DSOSCF .. RETURNS THE FOLDED TOTAL ONE-ELECTRON 
C               SCF DENSITY MATRIX
C     FSO ..... RETURNS FOLDED INTERMEDIATE MATRIX I IN
C               SYMMETRY ADAPTED AO BASIS
C     SCR ..... SCRATCH ARRAY OF SIZE NBASI*NBASIS
C     DENSCF .. LOGICAL FLAG TURNING ON ADDITIONAL
C               CONTRUCTION OF SCF DENSITY MATRIX  
C
C The DOSCF logical variable controls formation of the SCF density
C during a post-SCF calculation. The DOSCF flag is set to true during
C a correlated calculation and FALSE for standard SCF calculations (it
C is redundant for SCF). The required JOBARC records (SCF-ordered
C eigenvalues and eigenvectors) are saved in vtran.
C Ajith Perera, 01/2006.
CEND
C
#include "flags.h"
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SCF,NONHF,ROHF,SEMI,IDFGH,DENSCF,DOSCF,HFDFT_GRAD
      CHARACTER*8 LABELC,LABELF,LABELD,LABELI,LABELNHF
      INTEGER R,S,RS,U,V,UV
      INTEGER POP,VRT,DIRPRD
      DIMENSION CMO(1),FT(1),DV(1),DSO(1),DSOSCF(1),
     &          FSO(1),SCR(1)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
      COMMON/BASSPH/NBAS5(8),NBASIS5,NBASSQ5,NBASTT5
      COMMON/METHOD/IUHF,SCF,NONHF
      COMMON/OPENSH/ROHF,SEMI
      COMMON/IPRI/IPRINT
      COMMON/DFGH/IDFGH  
C
#include "flags.com"
#include "flags2.com"
C
      DATA AZERO,HALF,ONE,TWO /0.0D0,0.5D0,1.0D0,2.0D0/
C
C  ZERO TARGET ARRAYS
C
      hfdft_grad=.false.
      call getrec(-1,'JOBARC','HFDFTGRA',1, ihfdft_grad)
      if (ihfdft_grad .eq. 1) hfdft_grad = .true.
C
      IF(DENSCF .OR. DOSCF) CALL ZERO(DSOSCF,NBASTT)
      CALL ZERO(DSO,NBASTT)
      CALL ZERO(FSO,NBASTT)
C
C  SET OCCUPATION NUMBER FOR SCF
C
      IF(IUHF.EQ.0) THEN
C
C  RHF :
C
       OCCUP=TWO
C
      ELSE
C
C  UHF :
C 
       OCCUP=ONE
C
      ENDIF
C
C  LOOP OVER SPIN CASES
C
      DO 10 ISPIN=1,1+IUHF
C
C
C   READ MO COEFFICIENTS IN A SYMMETRY PACKED FORM
C
C   GET FIRST THE LABELS FOR THE VARIOUS METHODS
C
       IF(SCF.OR.DOSCF) THEN
C
        IF(ISPIN.EQ.1) THEN
         LABELC='SCFEVCA0'
         LABELF='SCFEVLA0'
         IF (DOSCF) THEN
            LABELC='EVSCFORA'
            LABELF='EESCFORA'
         END IF
         IF(ROHF.AND.(.NOT.SEMI)) LABELF='FOCKA   '
         if(hfdft_grad) labelf='HFKSVDAA'
        ELSE
         LABELC='SCFEVCB0'
         LABELF='SCFEVLB0'
         IF (DOSCF) THEN
            LABELC='EVSCFORB'
            LABELF='EESCFORB'
         END IF
         IF(ROHF.AND.(.NOT.SEMI)) LABELF='FOCKB '
         if(hfdft_grad) labelf='HFKSVDBB'
        ENDIF
       ELSE
C
C   FOR CORRELATION METHODS READ ALSO TOTAL DENSITY MATRIX
C
        IF(ISPIN.EQ.1) THEN
         LABELC='SCFEVCA0'
         LABELI='IINTERMA'
         LABELD='RELDENSA'
         LABELNHF='NHFDENSA'
        ELSE
         LABELC='SCFEVCB0'
         LABELI='IINTERMB'
         LABELD='RELDENSB'
         LABELNHF='NHFDENSB'
        ENDIF
C
       ENDIF
C
C  GET THE EIGEN VECTORS FROM JOBARC
C
       CALL GETREC(20,'JOBARC',LABELC,IINTFP*NBASIS5*NBASIS5,SCR)
C
C  SYMMETRY PACK THE EIGEN VECTORS
C 
CSSS        CALL SYMC(SCR,CMO,NBASIS5,NBAS5,SCF.OR.DOSCF,.FALSE.,ISPIN)

      CALL SYMCEXT(SCR,CMO,NBASIS5,NBAS5,NBASIS5,NBAS5,SCF,
     &             .FALSE.,ISPIN)
C
C  DEAL WITH THE FOCK-MATRIX PART (INTERMEDIATES IN CC)
C
       IF(SCF.OR.DOSCF) THEN
C
C FOR ROHF DO SOME FANCY STUFF
C (BUT ONLY FOR STANDARD ORBITALS)
C
        IF(ROHF.AND.(.NOT.SEMI)) THEN
         CALL GETREC(20,'JOBARC',LABELF,IINTFP*NBASIS5*NBASIS5,SCR)
CSSS         CALL SYMC(SCR,FT,NBASIS5,NBAS5,SCF,.TRUE.,ISPIN)

      CALL SYMCEXT(SCR,FT,NBASIS5,NBAS5,NBASIS5,NBAS5,SCF,
     &             .TRUE.,ISPIN)

         CALL FORMFIJ(FT,CMO,SCR,ISPIN) 
        ELSE IF (hfdft_grad) then
         call getrec(20,'JOBARC',LABELF,iintfp*nbasis5*nbasis5,scr)
CSSS         call symc(scr,ft,nbasis5,nbas5,scf,.true.,ispin)
      CALL SYMCEXT(SCR,FT,NBASIS5,NBAS5,NBASIS5,NBAS5,SCF,
     &             .TRUE.,ISPIN)
         call formfij(ft,cmo,scr,ispin)
        else
C
C   READ EIGENVALUES OF THE FOCK MATRIX
C
         CALL GETREC(20,'JOBARC',LABELF,IINTFP*NBASIS5,FT)
C
         IF(IFLAGS(3).EQ.2) CALL SYMF(SCR,FT,NBASIS5,ISPIN)
C
        ENDIF
C
       ELSE
C
C   READ IN TOTAL DENSITY MATRIX IN MO BASIS AND INTERMEDIATES
C
        CALL GETREC(20,'JOBARC',LABELD,IINTFP*NBASIS5*NBASIS5,SCR)
C
        IF(NONHF) THEN
C
C   FOR NON-HF CASES, ADD HERE THE NON-DIRECT PART OF THE DENSITY
C   (FOR A DETAILED DISCUSSION OF GRADIENTS IN NON-HF CCSD METHODS
C   SEE J.GAUSS, J.F. STANTON, R.. BARTLETT, J.CHEM. PHYS, TO BE PUBLISHED)
C
         CALL GETREC(20,'JOBARC',LABELNHF,IINTFP*NBASIS5*NBASIS5,DV)
         CALL SAXPY(NBASIS5*NBASIS5,ONE,DV,1,SCR,1)
        ENDIF
C
C   SYMMETRY PACK THE DENSITY MATRIX
C
        CALL SYMD(SCR,DV,NBASIS5,NBAS5,ISPIN)
C
C   READ INTERMEDIATES I(P,Q) FOR CORELATED METHODS
C
        CALL GETREC(20,'JOBARC',LABELI,IINTFP*NBASIS5*NBASIS5,SCR)
C
C   SYMMETRY PACK THE INTERMEDIATES
C
        CALL SYMI(SCR,FT,NBASIS5,NBAS5,ISPIN)
C
       ENDIF
C
C  IN CASE OF SPHERICAL GAUSSIANS, EXPAND TO FULL SET OF CARTESIAN
C  GAUSSIAN
C
       IF(IDFGH) THEN
C
CSSS        CALL EXPSPH(CMO,SCR)

       CALL GETREC(20,'JOBARC',LABELC,IINTFP*NBASIS5*NBASIS5,CMO)
       CALL GETREC(20,'JOBARC','CMP2CART',NBASIS5*NBASIS*IINTFP,
     &            SCR)
       I010=1+NBASIS*NBASIS5
       CALL XGEMM('N','N',NBASIS,NBASIS5,NBASIS5,1.D0,SCR,NBASIS,
     &            CMO,NBASIS5,0.D0,SCR(I010),NBASIS)
       CALL SYMCEXT(SCR(I010),CMO,NBASIS,NBAS,NBASIS5,NBAS5,SCF,
     &             .FALSE.,ISPIN)
       ENDIF 
C
C  PRINT SECTION
       IF (IPRINT.GT.5) THEN
        WRITE (6, '(//A/)') ' <<<<< SUBROUTINE DSOFSO >>>>>'
        WRITE (6, '(A,8I5)') ' NOCC ', (POP(I,ISPIN),I = 1,NIRREP)
        WRITE (6, '(A,8I5)') ' NVRT ', (VRT(I,ISPIN),I = 1,NIRREP)
        WRITE (6, '(A,8I5)') ' NBAS ', (NBAS(I),I = 1,NIRREP)
        WRITE (6, '(A,8I5)') ' NBAS5 ', (NBAS5(I),I = 1,NIRREP)
        IF (IPRINT .GE. 10) THEN
         CALL HEADER('Molecular orbitals',0)
         IEND = 0
         DO 1000 ISYM = 1,NIRREP
          IF (NBAS(ISYM) .EQ. 0) GOTO 1000
          WRITE (6, '(//,A,I5,/)')
     *           ' Symmetry ', ISYM
          IENDI = 0
          DO 1200 I = 1, NBAS5(ISYM)
           WRITE (6, '(/,A,I5,/)')
     *            ' Molecular orbital ', I
           WRITE (6, '(6F12.6)')
     *           (CMO(IEND+IENDI+J), J = 1, NBAS(ISYM))
           IENDI = IENDI + NBAS(ISYM)
1200      CONTINUE
          IEND = IEND + NBAS5(ISYM)*NBAS(ISYM)
1000     CONTINUE
         IF(.NOT.SCF) THEN
          CALL PRITRI(FT,NBASTT5,'Intermediates I(p,q) ')
          CALL PRITRI(DV,NBASTT5,'Total density matrix (MO basis)')
         ENDIF
        ENDIF
       ENDIF
       IPRINT=0
C
C  CONSTRUCT SO MATRICES
C
       ISEND=0
       ISEND5=0
       ICEND=0
       IFEND=0
       IOFFO=0     
       IOFFF=0
C
C  LOOP OVER ALL IRREPS
C
       DO 110 ISYM=1,NIRREP
        NOCCI=POP(ISYM,ISPIN)
        NBASI=NBAS(ISYM)
        NBASI5=NBAS5(ISYM)
        IF(NBASI.EQ.0) GOTO 120
        RS=0
C
C  LOOP OVER TARGET INDICES R AND S (ACTUALLY THESE ARE INDICES IN
C  THE SAO BASIS   
C
        DO 100 R=1,NBASI
         DO 200 S=1,R
          RS=RS+1
C
C     CONSTRUCT FIRST ONE-ELECTRON DENSITY MATRIX
C
          DTRS=AZERO
          DTRSSCF=AZERO
          FTRS=AZERO
C
C  SCF : CONTRUCT DENSITY MATRIX FROM EIGEN VECTORS ONLY
C
          IF(SCF.OR.DOSCF)THEN
           IF(NOCCI.GT.0)THEN
            ICENDI=ICEND
C
C  LOOP OVER OCCUPIED ORBITALS 
C
            IF(ROHF.AND.(.NOT.SEMI)) THEN
             DO 310 I=1,NOCCI
              DTRS=DTRS+OCCUP*CMO(ICENDI+R)*CMO(ICENDI+S)
              ICENDI=ICENDI+NBASI
310          CONTINUE
             ICENDI=ICEND
             DO 320 I=1,NOCCI
              ICENDJ=ICEND
              DO 330 J=1,NOCCI
               FTRS=FTRS+OCCUP*CMO(ICENDI+R)*CMO(ICENDJ+S)
     &                        *FT(IOFFF+J+(I-1)*NOCCI)
               ICENDJ=ICENDJ+NBASI
330           CONTINUE
              ICENDI=ICENDI+NBASI
320          CONTINUE
C
C For HF-DFT grads, P. Verma and A. Perera 08/2008
C
             else if (hfdft_grad) then
               do 3101 i=1,nocci
                  dtrs=dtrs+occup*cmo(icendi+r)*cmo(icendi+s)
                  icendi=icendi+nbasi
3101           continue
               icendi=icend
               do 3201 i=1,nocci
                  icendj=icend
                  do 3301 j=1,nocci
                     ftrs=ftrs+occup*cmo(icendi+r)*cmo(icendj+S)
     &                    *ft(iofff+j+(i-1)*nocci)
                     icendj=icendj+nbasi
3301              continue
C
                  icendi=icendi+nbasi
3201           continue
C
            ELSE
             DO 300 I=1,NOCCI
              DTRS=DTRS+OCCUP*CMO(ICENDI+R)*CMO(ICENDI+S)
              FTRS=FTRS+OCCUP*CMO(ICENDI+R)*CMO(ICENDI+S)
     &                       *FT(IOFFO+I)
              ICENDI=ICENDI+NBASI
  300        CONTINUE
            ENDIF
            IF(R.NE.S)THEN
             FTRS=FTRS+FTRS
             DTRS=DTRS+DTRS
            ENDIF
           ENDIF
          ELSE
C
C  CORRELATION METHODS : TRANSFORM THE TOTAL MO DENSITY TO AO BASIS
C
           IF(NBASI.GT.0)THEN
            UV=ISEND5
            IDVEND=ICEND 
            ICENDU=IDVEND
C
C  LOOP OVER ALL MOS U AND V
C
            DO 400 U=1,NBASI5
             ICENDV=IDVEND
             DO 410 V=1,U
              UV=UV+1
              DUV=DV(UV)
              TEMP=CMO(ICENDU+R)*CMO(ICENDV+S)
              IF(U.NE.V)TEMP=TEMP+CMO(ICENDU+S)*CMO(ICENDV+R)
              DTRS=DTRS+DUV*TEMP
              ICENDV=ICENDV+NBASI
  410        CONTINUE
             ICENDU=ICENDU+NBASI
  400       CONTINUE
            UV=IFEND
            IDVEND=ICEND
            ICENDU=IDVEND 
            DO 500 U=1,NBASI5
             ICENDV = IDVEND
             DO 510 V=1,NBASI5
              UV=UV+1
              FUV = FT(UV)
              TEMP= CMO(ICENDU+R)*CMO(ICENDV+S)
              IF(R.NE.S) TEMP=TEMP+CMO(ICENDU+S)*CMO(ICENDV+R)
              FTRS=FTRS-FUV*TEMP
              ICENDV = ICENDV + NBASI
510          CONTINUE
             ICENDU=ICENDU+NBASI
500         CONTINUE  
           ENDIF
C 
C  FACTOR TWO FOR THE OFF-DIAGONAL ELEMENTS SINCE WE WANT TO
C  CALCULATE THE FOLDED DENSITY MATRIX
C
           IF (R .NE. S) THEN
            DTRS = DTRS + DTRS
           ENDIF
           IF(DENSCF)THEN
            IF(NOCCI.GT.0)THEN
             ICENDI=ICEND
C
C  LOOP OVER OCCUPIED ORBITALS 
C
             IF(ROHF.AND.(.NOT.SEMI)) THEN
              DO 610 I=1,NOCCI
               DTRSSCF=DTRSSCF+OCCUP*CMO(ICENDI+R)*CMO(ICENDI+S)
               ICENDI=ICENDI+NBASI
610           CONTINUE
             ELSE
              DO 600 I=1,NOCCI
               DTRSSCF=DTRSSCF+OCCUP*CMO(ICENDI+R)*CMO(ICENDI+S)
               ICENDI=ICENDI+NBASI
600           CONTINUE
             ENDIF
             IF(R.NE.S)THEN
              DTRSSCF=DTRSSCF+DTRSSCF
             ENDIF
            ENDIF
           ENDIF
          ENDIF
          DSO(ISEND+RS) = DSO(ISEND+RS) + DTRS
          FSO(ISEND+RS) = FSO(ISEND+RS) + FTRS
          IF(DENSCF)DSOSCF(ISEND+RS)=DSOSCF(ISEND+RS)+DTRSSCF
C
  200    CONTINUE
  100   CONTINUE
C
C  Print Section
        IF (IPRINT .GE. 10) THEN
         WRITE (6,'(1X,A,I5)') ' Symmetry', ISYM
         CALL PRITRI(DSO(ISEND+1),NBASI,
     *         'Total density matrix (SO basis)')
         CALL PRITRI(FSO(ISEND+1),NBASI,
     *         'Total Fock matrix (SO basis)')
        ENDIF
C
 120   CONTINUE
C
       IFEND = IFEND + NBASI5*NBASI5
       ISEND = ISEND + (NBASI*(NBASI + 1))/2
       ISEND5 = ISEND5 + (NBASI5*(NBASI5 + 1))/2
       ICEND = ICEND + NBASI*NBASI5
       IOFFO = IOFFO + NBASI5
       IOFFF = IOFFF + NOCCI*NOCCI
 110   CONTINUE
C
 10   CONTINUE
      IPRINT = 0
C
      RETURN
      END
