

      SUBROUTINE SETDEP(TMAT,AMAT,DONE,NPREF,ITBDEP)
C
C     Determines a set of Cartesian coordinates whose gradient and
C     Hessian elements can be determined by translational and rotational
C     symmetry.
C
C     Dec 1988 tuh
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mxcent.par"
#include "baslims.par"
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      PARAMETER (D0 = 0.0D0)
      LOGICAL OK, OKALL, ONCE, DONE(NCOOR)
      DIMENSION TMAT(6,NCOOR), AMAT(6,NCOOR), NPREF(NCOOR),
     *          ITBDEP(NCOOR)
C
      LOGICAL         COMPAR, SKIP, NOROT,
     *                HESTRO, GDTRO, RDTRO, TROGRD, TROHES, TRODIP
      COMMON /CBITRO/ THRESH, IPRINT, COMPAR, SKIP, NOROT,
     *                HESTRO, GDTRO, RDTRO, TROGRD, TROHES, TRODIP
      CHARACTER NAMEX*6
      LOGICAL DCORD, DCORGD, NOORBT, DOPERT
      COMMON /NUCLEIi/ NOORBT(MXCENT),
     &                NUCIND, NUCDEP, NUCPRE(MXCENT), NUCNUM(MXCENT,8),
     &                NUCDEG(MXCENT), ISTBNU(MXCENT), NDCORD,
     &                NDCOOR(MXCOOR), NTRACO, NROTCO, ITRACO(3),
     &                IROTCO(3),
     &                NATOMS, NFLOAT,
     &                IPTGDV(3*MXCENT),
     &                NGDVEC(8), IGDVEC(8)
      COMMON /NUCLEI/ CHARGE(MXCENT), CORD(MXCENT,3),
     &                DCORD(MXCENT,3),DCORGD(MXCENT,3),
     &                DOPERT(0:3*MXCENT)
      COMMON /NUCLEC/ NAMEX(MXCOOR)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
C
C     NCOOR  - total number of Cartesian coordinates
C     NCOORS - number of Cartesian coordinates of each symmetry
C     DEPEND - true for dependent symmetry coordinates
C     NCDEP  - total number of dependent/trarot coordinates
C     NCDEPS - number of dependent/trarot coordinates of each symmetry
C     NCIND  - total number of independent coordinates
C     NCINDS - number of independent coordinates of each symmetry
C     SYMTRO - symmetry of trarot coordinates Tx, Ty, Tz, Rx, Ry, Rz
C     IPTTRO - identifies symmetry-ordered trarot coordinate as
C              Tx, Ty, Tz, Rx, Ry, or Rz
C
      LOGICAL DEPEND
      COMMON /TRKOORi/ NCOOR, NCOORS(8),
     *                NCDEP, NCDEPS(8), NCIND, NCINDS(8),
     *                IPTTRO(6,8), NTRREP(0:7)
      COMMON /TRKOOR/ DEPEND(MXCOOR), SYMTRO(6)
C
      IF (IPRINT .GE. 5) THEN
         CALL TITLER('Output from SETDEP','*',103)
      END IF
      NSET = 0
      DO 100 I = 1, NCOOR
         DONE(I) = .FALSE.
  100 CONTINUE
C
C     (A) Eliminate if possible one atom completely
C     =============================================
C
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('Coordinates for atom completely eliminated',-1)
      END IF
C
C     Select atom of multiplicity one and highest charge
C
      MATOM  = 0
      CHAMAX = D0
      DO 200 IATOM = 1, NUCIND
         IF (MULT(ISTBNU(IATOM)) .EQ. 1) THEN
            CHARGI = CHARGE(IATOM)
            IF (CHARGI .GT. CHAMAX) THEN
               MATOM  = IATOM
               CHAMAX = CHARGI
            END IF
         END IF
  200 CONTINUE
C
C     Set coordinates of this atom - if any - as dependent
C
      IF (MATOM .GT. 0) THEN
         DO 210 ICOOR = 1, 3
            DO 220 IREP = 0, MAXLOP
               ISCOOR = IPTCNT(3*(MATOM - 1) + ICOOR,IREP)
               IF (ISCOOR .GT. 0) THEN
                  CALL SCOPY(NCDEP,TMAT(1,ISCOOR),1,AMAT(1,NSET+1),1)
                  IF (IPRINT .GE. 10) WRITE (LUPRI,'(A,I5)')
     *                 ' Calling LINDEP for coordinate ',ISCOOR
                  CALL LINDEP(AMAT,THRESH,OK,6,NCOOR,NSET,1,IPRINT)
                  IF (OK) THEN
                     NSET = NSET + 1
                     ITBDEP(NSET) = ISCOOR
                     DONE(ISCOOR) = .TRUE.
                     IF (NSET .EQ. NCDEP) GO TO 1000
                  END IF
               END IF
  220       CONTINUE
  210    CONTINUE
      END IF
      IF (IPRINT .GE. 5) THEN
         WRITE (LUPRI,'(/A,I5)')
     *      ' Number of dependent coordinates found so far: ',NSET
         IF (NSET .GT. 0) WRITE (LUPRI,'(A,6I5)') ' Coordinates found:',
     *                                              (ITBDEP(I),I=1,NSET)
      END IF
C
C     (B) Eliminate if possible one coordinate for all symmetries
C     ===========================================================
C
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('Coordinates eliminated for all symmetries',-1)
      END IF
C
C     Identify Cartesian directions which may be eliminated for
C     all symmetries
C
      NDIR = 0
      IDIR = 0
      DO 300 IATOM = 1, NUCIND
         DO 310 ICOOR = 1, 3
            IDIR = IDIR + 1
            OKALL = .TRUE.
            ONCE  = .FALSE.
            DO 320 IREP = 0, MAXLOP
               ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP)
               IF (ISCOOR .GT. 0) THEN
                  IF (.NOT.DONE(ISCOOR)) THEN
                     ONCE = .TRUE.
                     CALL SCOPY(NCDEP,TMAT(1,ISCOOR),1,AMAT(1,NSET+1),1)
                     IF (IPRINT .GE. 10) WRITE (LUPRI,'(A,I5)')
     *                    ' Calling LINDEP for coordinate ',ISCOOR
                     CALL LINDEP(AMAT,THRESH,OK,6,NCOOR,NSET,1,IPRINT)
                     OKALL = OKALL .AND. OK
                  END IF
               END IF
  320       CONTINUE
            IF (ONCE .AND. OKALL) THEN
               NDIR = NDIR + 1
               NPREF(NDIR) = IDIR
            END IF
  310    CONTINUE
  300 CONTINUE
      IF (IPRINT .GE. 5) THEN
         IF (NDIR .GT. 0) THEN
            CALL HEADER('Unsorted list of directions',1)
            WRITE (LUPRI,'(12I5)') (NPREF(I),I=1,NDIR)
         ELSE
            WRITE (LUPRI,'(A)') ' No coordinates found.'
         END IF
      END IF
C
C     Sort these coordinates by charge of their atoms
C
      DO 330 I = 1, NDIR - 1
         IDIR = NPREF(I)
         CHARGI = CHARGE((IDIR + 2)/3)
         DO 340 J = I + 1, NDIR
            JDIR = NPREF(J)
            CHARGJ = CHARGE((JDIR + 2)/3)
            IF (CHARGJ .GT. CHARGI) THEN
               NPREF(I) = JDIR
               NPREF(J) = IDIR
            END IF
  340    CONTINUE
  330 CONTINUE
      IF (IPRINT .GE. 5) THEN
         IF (NDIR .GT. 0) THEN
            CALL HEADER('Sorted list of directions',1)
            WRITE (LUPRI,'(12I5)') (NPREF(I),I=1,NDIR)
         END IF
      END IF
C
C     Set these coordinates - if any - as dependent
C
      DO 350 IDIR = 1, NDIR
         DO 360 IREP = 0, MAXLOP
            ISCOOR = IPTCNT(NPREF(IDIR),IREP)
            IF (ISCOOR .GT. 0) THEN
               IF (.NOT.DONE(ISCOOR)) THEN
                  CALL SCOPY(NCDEP,TMAT(1,ISCOOR),1,AMAT(1,NSET+1),1)
                  IF (IPRINT .GE. 10) WRITE (LUPRI,'(A,I5)')
     *                 ' Calling LINDEP for coordinate ',ISCOOR
                  CALL LINDEP(AMAT,THRESH,OK,6,NCOOR,NSET,1,IPRINT)
                  IF (OK) THEN
                     NSET = NSET + 1
                     ITBDEP(NSET) = ISCOOR
                     DONE(ISCOOR) = .TRUE.
                     IF (NSET .EQ. NCDEP) GO TO 1000
                  END IF
               END IF
            END IF
  360    CONTINUE
  350 CONTINUE
      IF (IPRINT .GE. 5) THEN
         WRITE (LUPRI,'(/A,I5)')
     *      ' Number of dependent coordinates found so far: ',NSET
         IF (NSET .GT. 0) WRITE (LUPRI,'(A,6I5)') ' Coordinates found:',
     *                                              (ITBDEP(I),I=1,NSET)
      END IF
C
C     (D) Set remaining coordinates until all have been set
C     =====================================================
C
      IF (IPRINT .GE. 5) CALL HEADER('Remaining coordinates',-1)
C
C     Sort coordinates by charge of their atoms
C
      NDIR = 0
      DO 400 I = 1, NCOOR
         IF (.NOT.DONE(I)) THEN
            NDIR = NDIR + 1
            NPREF(NDIR) = I
         END IF
  400 CONTINUE
      DO 410 I = 1, NDIR
         ISCOOR = NPREF(I)
         CHARGI = CHARGE((IPTCOR(ISCOOR) + 2)/3)
         DO 420 J = I + 1, NDIR
            JSCOOR = NPREF(J)
            CHARGJ = CHARGE((IPTCOR(JSCOOR) + 2)/3)
            IF (CHARGJ .GT. CHARGI) THEN
               NPREF(I) = JSCOOR
               NPREF(J) = ISCOOR
            END IF
  420    CONTINUE
  410 CONTINUE
      IF (IPRINT .GE. 5) THEN
         IF (NDIR .GT. 0) THEN
            CALL HEADER('Sorted list of directions',1)
            WRITE (LUPRI,'(12I5)') (NPREF(I),I=1,NDIR)
         END IF
      END IF
      DO 430 IDIR = 1, NDIR
         ISCOOR = NPREF(IDIR)
         CALL SCOPY(NCDEP,TMAT(1,ISCOOR),1,AMAT(1,NSET+1),1)
         IF (IPRINT .GE. 10) WRITE (LUPRI,'(A,I5)')
     *       ' Calling LINDEP for coordinate ',ISCOOR
         CALL LINDEP(AMAT,THRESH,OK,6,NCOOR,NSET,1,IPRINT)
         IF (OK) THEN
            NSET = NSET + 1
            ITBDEP(NSET) = ISCOOR
            DONE(ISCOOR) = .TRUE.
            IF (NSET .EQ. NCDEP) GO TO 1000
         END IF
  430 CONTINUE
 1000 CONTINUE
C
C     Set up DEPEND(NCOOR)
C
      ICOOR = 0
      DO 500 ISCOOR = 1, NCOOR
         DEPEND(ISCOOR) = .FALSE.
  500 CONTINUE
      DO 510 IDEP = 1, NCDEP
         DEPEND(ITBDEP(IDEP)) = .TRUE.
  510 CONTINUE
      IF (IPRINT .GE. 5) THEN
         WRITE (LUPRI,'(/A,I5)')
     *      ' Total number of dependent coordinates found: ',NSET
         IF (NSET .GT. 0) WRITE (LUPRI,'(A,6I5)') ' Coordinates found:',
     *                                              (ITBDEP(I),I=1,NSET)
         CALL HEADER('DEPEND in SETDEP',1)
         WRITE (LUPRI,'(2X,35L2)') (DEPEND(I),I=1,NCOOR)
      END IF
      IF (NSET .NE. NCDEP) THEN
         WRITE (LUPRI,'(/,2A,/,A,/)')
     *         ' ERROR in SETDEP - less than 6 (5)',
     *         ' independent Cartesian directions',' found.'
         STOP
      END IF
      RETURN
      END
