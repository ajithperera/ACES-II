      SUBROUTINE ONEINT(NTAP,H1,NSTO,POTNUC,BUF,IBUF,NSIZ,LEN2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 OVRLAP,ONEHAM,APOTL,SYMTRN,ONEEL,TWOEL,STAR
      DIMENSION TITLE(24),NSO(49),BUF(LEN2),IBUF(LEN2)
      DIMENSION H1(NSIZ)
C     EQUIVALENCE (LAB(1),IPACK)
C      DATA OVLAP/8HOVERLAP /, ONEHAM/8HONEHAMIL/
C      DATA FLAGM/8HTWOELSUP/, FLAGU/8HBASTWOEL/
      DATA OVRLAP/'OVERLAP '/,ONEHAM/'ONEHAMIL'/
      DATA APOTL/'LOCALPOT'/
      DATA SYMTRN/'SYMTRAN '/
      DATA ONEEL/'ONEHAMIL'/
      DATA TWOEL/'TWOELSUP'/
      DATA STAR /'********'/
      DATA ZERO/0.D00/
C      IUNPK1(IH2)=AND(ISHFT(IH2,-48),65535)
C      IUNPK2(IH2)=AND(ISHFT(IH2,-32),65535)
C      IUNPK3(IH2)=AND(ISHFT(IH2,-16),65535)
C      IUNPK4(IH2)=AND(IH2,65535)
C
C     65535 = 2**16 -1 ;  This puts 1's in the first 16 bits.
C     OPEN (UNIT=NTAP,FILE='INT',
      OPEN (UNIT=NTAP,FILE='IIII',
     1       STATUS='UNKNOWN',FORM='UNFORMATTED')
C GET NUMBER OF BASIS FUNCTIONS OFF OF MOLECULE.
CVAX     7
C     READ(NTAP) ILBM,(ILBM(I),I=1,22),NSYMHF,(NSO(I),I=1,NSYMHF),
C    *             REPNRG
C     WRITE(6,11)REPNRG,NSYMHF,(NSO(I),I=1,NSYMHF)
C  11 FORMAT(/,10X,' The nuclear repulsion energy is ',F20.10,
C    *       /,10X,' The number of symmetry types is ',I3,
C    *       /,10X,' The number of each sym. type are',
C    *       /,10X,8I6)
CFPS     7
      READ(NTAP) TITLE,NSYMHF,(NSO(I),I = 1,NSYMHF),POTNUC
      WRITE(6,711)TITLE,NSYMHF,(NSO(I),I = 1,NSYMHF)
  711 FORMAT(/,10X,' The AO integral title is ',12A6,
     *       /, 10X,'                         ',12A6,
     *       /,10X,' The number of symmetry types is ',I3,
     *       /,10X,' The number of each sym. type are',
     *       /,10X,8I6)
      WRITE(6,*) ' nuclear repulsion = ',POTNUC
      NSTO=0
      IVECT=0
      DO 12 I=1,NSYMHF
      IVECT=IVECT+((NSO(I)+1)*NSO(I))/2
   12 NSTO=NSTO+NSO(I)
      WRITE(6,13)NSTO,IVECT
   13 FORMAT(/,10X,' The number of basis functions is ',I3,
     *       /,10X,' The number of onel elements is ',I5)
      NIJ=(NSTO*(NSTO+1))/2
C     CALL SEARCH (OVRLAP,NTAP)
C 150 READ (NTAP) BUF,IBUF,NUM
C     IF (NUM.EQ.-1) GO TO 250
C     GO TO 150
C 250 CONTINUE
      DO 251 I=1,IVECT
  251 H1(I)=ZERO
      CALL SEARCH(ONEEL,NTAP)
  350 READ (NTAP) BUF,IBUF,NUM
      IF (NUM.EQ.-1) GO TO 450
      DO 360 I=1,NUM
      II=IBUF(I)
      IF (II.NE.0) H1(II)=BUF(I)
  360 CONTINUE
      GO TO 350
  450 CONTINUE
C     write(6,1000) (H1(II),II=1,IVECT)
 1000 FORMAT(1h ,6F10.4)
      CLOSE(NTAP)
      RETURN
      END
