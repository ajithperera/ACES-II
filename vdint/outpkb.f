      SUBROUTINE OUTPKB (AMATRX,NDIM,NBLK,NCTL,LUPRI)
C.......................................................................
C
C OUTPKB is OUTPAK modified for blocking as specified by
C        NDIM(NBLK).  16-Dec-1983 Hans Jorgen Aa. Jensen.
C
C 16-Jun-1986 hjaaj ( removed Hollerith )
C
C OUTPAK PRINTS A REAL OMSYMMETRIC MATRIX STORED IN ROW-PACKED LOWER
C TRIANGULAR FORM (SEE DIAGRAM BELOW) IN FORMATTED FORM WITH NUMBERED
C ROWS AND COLUMNS.  THE INPUT IS AS FOLLOWS:
C
C        AMATRX(')...........PACKED MATRIX, blocked
C
C        NDIM(').............dimension of each block
C
C        NBLK................number of blocks
C
C        NCTL................CARRIAGE CONTROL FLAG: 1 FOR SINGLE SPACE,
C                                                   2 FOR DOUBLE SPACE,
C                                                   3 FOR TRIPLE SPACE.
C
C THE MATRIX ELEMENTS in a block ARE ARRANGED IN STORAGE AS
C FOLLOWS:
C
C        1
C        2    3
C        4    5    6
C        7    8    9   10
C       11   12   13   14   15
C       16   17   18   19   20   21
C       22   23   24   25   26   27   28
C       AND SO ON.
C
C OUTPAK IS SET UP TO HANDLE 6 COLUMNS/PAGE WITH A 6F20.14 FORMAT
C FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED, CHANGE
C FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER OF
C COLUMNS.
C
C AUTHOR:  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
C          FLORIDA, GAINESVILLE
C..........OUTPAK VERSION = 09/05/73/03
C.......................................................................
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NDIM(NBLK),BEGIN
      DIMENSION AMATRX(2)
      CHARACTER*1 ASA(3),BLANK,CTL
      CHARACTER   PFMT*20, COLUMN*8
      PARAMETER (ZERO = 0.0D0, KCOLP=4, KCOLN=6)
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA COLUMN/'Column  '/, ASA/' ', '0', '-'/, BLANK/' '/
C
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
C
      MATLN = 0
      DO 200 IBLK = 1,NBLK
         MATLN = MATLN + NDIM(IBLK)*(NDIM(IBLK)+1)/2
  200 CONTINUE
C
      AMAX = ZERO
      DO 205 I=1,MATLN
         AMAX = MAX( AMAX, ABS(AMATRX(I)) )
  205 CONTINUE
      IF (AMAX .EQ. ZERO) THEN
         WRITE (LUPRI,3000) NBLK
         GO TO 800
      END IF
      IF (FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX) THEN
C        use F output format
         PFMT = '(A1,I7,2X,8F15.8)'
      ELSE
C        use 1PD output format
         PFMT = '(A1,I7,2X,1P8D15.6)'
      END IF
C
      IOFF = 0
      DO 100 IBLK = 1,NBLK
         IDIM = NDIM(IBLK)
      IF (IDIM.EQ.0) GO TO 100
         IIDIM = IDIM*(IDIM+1)/2
C
         DO 5 I=1,IIDIM
            IF (AMATRX(IOFF+I).NE.ZERO) GO TO 15
    5    CONTINUE
         WRITE (LUPRI,1100) IBLK
         GO TO 90
   15    CONTINUE
         WRITE (LUPRI,1200) IBLK
C
C LAST IS THE LAST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED
C
         LAST = MIN(IDIM,KCOL)
C
C BEGIN IS THE FIRST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED.
C
         BEGIN = 1
C.....BEGIN NON STANDARD DO LOOP.
 1050       NCOL = 1
            WRITE (LUPRI,1000) (COLUMN,I,I = BEGIN,LAST)
            KTOTAL = IOFF + BEGIN*(BEGIN+1)/2 - 1
            DO 40 K = BEGIN,IDIM
               DO 10 I = 1,NCOL
                  IF (AMATRX(KTOTAL+I) .NE. ZERO) GO TO 20
   10          CONTINUE
               GO TO 30
   20          WRITE (LUPRI,PFMT) CTL,K,(AMATRX(KTOTAL+J),J=1,NCOL)
   30          IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
               KTOTAL = KTOTAL + K
   40       CONTINUE
            LAST = MIN(LAST+KCOL,IDIM)
            BEGIN=BEGIN+NCOL
         IF (BEGIN.LE.IDIM) GO TO 1050
C
   90    IOFF = IOFF + IIDIM
  100 CONTINUE
C
  800 CONTINUE
      RETURN
 3000 FORMAT (/5X,'All',I3,' blocks zero matrices.')
 1100 FORMAT (/5X,'*** Block',I3,' zero matrix. ***')
 1200 FORMAT (/5X,'*** Block',I3,' ***')
 1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))
C2000 FORMAT (A1,'Row',I4,2X,1P8D15.6)
C2000 FORMAT (A1,I7,2X,1P8D15.6)
      END
