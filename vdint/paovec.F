      SUBROUTINE PAOVEC(MAXVEC,IATOM,IPRINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "baslims.par"
#include "mxcent.par"
      LOGICAL ORBOUT(MXSHEL), ALONE(MXSHEL), ALLSHR
      DIMENSION ISHELL(MXSHEL), IPRVOR(MXSHEL), IPRVPR(MXSHEL),
     *          JORBSH(MXSHEL,MXAOVC), KORBSH(MXSHEL,MXAOVC)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      LOGICAL SHARE
      COMMON /SHELLSi/ KMAX,
     &                NHKT(MXSHEL),   KHKT(MXSHEL), MHKT(MXSHEL),
     &                ISTBAO(MXSHEL), NUCO(MXSHEL), JSTRT(MXSHEL),
     &                NSTRT(MXSHEL),  MST(MXSHEL),  NCENT(MXSHEL),
     &                NRCO(MXSHEL), NUMCF(MXSHEL),
     &                NBCH(MXSHEL),   KSTRT(MXSHEL)
      COMMON /SHELLS/ CENT(MXSHEL,3), SHARE(MXSHEL)
      LOGICAL BIGVEC, SEGMEN
      COMMON /BLOCKS/ CENTSH(MXSHEL,3),
     *                MAXSHL, BIGVEC, SEGMEN,
     *                NHKTSH(MXSHEL), KHKTSH(MXSHEL), MHKTSH(MXSHEL),
     *                ISTBSH(MXSHEL), NUCOSH(MXSHEL), NORBSH(MXSHEL),
     *                NSTRSH(MXSHEL), NCNTSH(MXSHEL), NSETSH(MXSHEL),
     *                JSTRSH(MXSHEL,MXAOVC),
     *                NPRIMS(MXSHEL,MXAOVC),
     *                NCONTS(MXSHEL,MXAOVC),
     *                IORBSH(MXSHEL,MXAOVC),
     *                IORBSB(MXCORB), NRCSH(MXSHEL)
      COMMON /PINCOM/ IPIND(MXCORB), IBLOCK(MXCORB), INDGEN(MXCORB)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
      MXAOSH = 0
C
C     First assign each shell to a vector
C
      DO 100 I = 1, KMAX
         ORBOUT(I) = .FALSE.
 100  CONTINUE
C
C     Determine ALONE()
C
      ALLSHR = .TRUE.
      DO 110 I = 1, KMAX
         ALLSHR = ALLSHR .AND. .NOT.SHARE(I)
 110  CONTINUE
      DO 130 I = 1, KMAX
         ALONE(I) = .NOT.ALLSHR .AND. .NOT.SHARE(I)
 130  CONTINUE
C
      MAXSHL = 0
      NAOOUT = 0
 1000 CONTINUE
         DO 200 IA = 1, KMAX
            IF (.NOT.ORBOUT(IA)) THEN
               NHKTA  = NHKT(IA)
               KHKTA  = KHKT(IA)
               MHKTA  = MHKT(IA)
               NCENTA = NCENT(IA)
               NBCHA  = NBCH(IA)
               NUCA   = NUCO(IA)
               MAXSHL = MAXSHL + 1
               NUCSH  = 0
               NHKTSH(MAXSHL) = NHKTA
               KHKTSH(MAXSHL) = KHKTA
               MHKTSH(MAXSHL) = MHKTA
               ISTBSH(MAXSHL) = ISTBAO(IA)
               IF (BIGVEC) THEN
                  IF (NCENTA .EQ. IATOM) THEN
                     NCNTSH(MAXSHL) = NCENTA
                  ELSE
                     NCNTSH(MAXSHL) = 0
                  END IF
                  IF (MXAOSH .EQ. 0) THEN
                      IF (NHKTA .EQ. 1) MAXVC = 10
                      IF (NHKTA .EQ. 2) MAXVC = 5
                      IF (NHKTA .GE. 3) MAXVC = 3
                  ELSE
                      MAXVC = MXAOSH
                  END IF
               ELSE
                  NCNTSH(MAXSHL) = NCENTA
                  CENTSH(MAXSHL,1) = CENT(IA,1)
                  CENTSH(MAXSHL,2) = CENT(IA,2)
                  CENTSH(MAXSHL,3) = CENT(IA,3)
                  IF (MXAOSH .EQ. 0) THEN
                      MAXVC = NUCA
                  ELSE
                      MAXVC = MXAOSH
                  END IF
               END IF
               IORB = 1
               NSET = 0
               LBCH = 0
               DO 210 IB = 1, KMAX
                  IF (ORBOUT(IB))                              GO TO 210
                  IF (NHKT(IB) .NE. NHKTA)                     GO TO 210
                  IF (KHKT(IB) .NE. KHKTA)                     GO TO 210
                  IF (MHKT(IB) .NE. MHKTA)                     GO TO 210
                  IF (NCENT(IB) .NE. NCENTA) THEN
                     IF (.NOT.BIGVEC)                          GO TO 210
                     IF (BIGVEC .AND. ALONE(IA))               GO TO 210
                     IF (BIGVEC .AND. ALONE(IB))               GO TO 210
                  END IF
                  IF (NCENTA.EQ.IATOM.AND.NCENT(IB).NE.IATOM)  GO TO 210
                  IF (NCENTA.NE.IATOM.AND.NCENT(IB).EQ.IATOM)  GO TO 210
C
C                 Correct orbital type
C
                  NUCB = NUCO(IB)
                  IF (SEGMEN) THEN
                     IF (NUCSH + NUCB .LE. MAXVC) THEN
                        ORBOUT(IB) = .TRUE.
                        ISHELL(IB) = MAXSHL
                        JORBSH(MAXSHL,IORB) = IB
                        IORB       = IORB   + 1
                        NAOOUT     = NAOOUT + 1
                        NUCSH      = NUCSH  + NUCB
                        NSET       = 1
                     END IF
                  ELSE
                     NBCHB = NBCH(IB)
                     IF (NBCHB .EQ. LBCH .OR.
     *                  NUCSH + NUCB .LE. MAXVC) THEN
                        ORBOUT(IB) = .TRUE.
                        ISHELL(IB) = MAXSHL
                        DO 11 IR=1,NRCO(IB)
                        JORBSH(MAXSHL,IORB) = IB
                        IORB       = IORB+1
11                      CONTINUE
                        NAOOUT     = NAOOUT + 1
                        IF (NBCHB .NE. LBCH) THEN
                           NSET  = NSET + 1
                           LBCH  = NBCHB
                           NUCSH = NUCSH  + NUCB
                           KORBSH(MAXSHL,NSET) = IB
                        END IF
                     END IF
                  END IF
  210          CONTINUE
               NUCOSH(MAXSHL) = NUCSH
               NSETSH(MAXSHL) = NSET
               IF (NAOOUT .LT. KMAX) THEN
                   GO TO 1000
               ELSE
                   GO TO 1100
               END IF
            END IF
  200    CONTINUE
      IF (NAOOUT .LT. KMAX) GO TO 1000
 1100 CONTINUE
C
C     All orbitals have now been assigned to a shell
C
C
C     Determine NORBSH()
C
      DO 300 I = 1, MAXSHL
         NORBSH(I) = 0
         DO 310 J = 1, KMAX
            IF (ISHELL(J) .EQ. I) NORBSH(I) = NORBSH(I) + NRCO(J)
  310    CONTINUE
  300 CONTINUE
C
C     Determine NPRIMS(MAXSHL,NSETSH) and NCONTS(MAXSHL,NSETSH)
C
      DO 400 I = 1, MAXSHL
         DO 410 K = 1, NSETSH(I)
            NPRIMS(I,K) = NUCO (KORBSH(I,K))
            NCONTS(I,K) = NRCO (KORBSH(I,K))
  410    CONTINUE
  400 CONTINUE
C
C     Determine and NSTRSH()
C
      NSTRSH(1) = 0
      DO 500 IA = 1,MAXSHL - 1
         NSTRSH(IA + 1)   = NSTRSH(IA)   + KHKTSH(IA)*NORBSH(IA)
  500 CONTINUE
C
C     Determine IBLOCK
C
      IJ = 0
      DO 700 I = 1, KMAX
         DO 710 J = 1, KHKT(I)
            IJ = IJ + 1
            IBLOCK(IJ) = ISHELL(I)
  710    CONTINUE
  700 CONTINUE
C
C     Determine IORBSH(MAXSHL,NORBSH)
C
      IORB = 0
      IPRVOR(1) = 0
      DO 800 I = 1, KMAX - 1
         IORB = IORB + KHKT(I)*MULT(ISTBAO(I))*NRCO(I)
         IPRVOR(I + 1) = IORB
  800 CONTINUE
      DO 810 I = 1, MAXSHL
         DO 820 J = 1, NORBSH(I)
            IORBSH(I,J) = IPRVOR(JORBSH(I,J)) + 1
  820    CONTINUE
  810 CONTINUE
C
C     Determine JSTRSH(MAXSHL,NSETSH)
C
      INUC = 0
      IPRVPR(1) = 0
      DO 850 I = 1, KMAX - 1
         IF (NUMCF(I + 1) .EQ. 1) INUC = INUC + NUCO(I)
         IPRVPR(I + 1) = INUC
  850 CONTINUE
      DO 815 I = 1, MAXSHL
         DO 825 J = 1, NSETSH(I)
            JSTRSH(I,J) = IPRVPR(KORBSH(I,J))
  825    CONTINUE
  815 CONTINUE
C
C     IORBSB()
C
      IORBB = 1
      DO 900 I = 1, MAXSHL
         DO 910 J = 1, NORBSH(I)
            IORBS = IORBSH(I,J)
            DO 920 ICMP = 1, KHKTSH(I)
               IORBSB(IORBS) = IORBB
               IORBB = IORBB + 1
               IORBS = IORBS + 1
  920       CONTINUE
  910    CONTINUE
  900 CONTINUE
      NTOTAL = IORBB
C
C     Print Section
C
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('Output from PAOVEC',-1)
         WRITE (LUPRI,'(1X,A,L5)') ' BIGVEC ', BIGVEC
         WRITE (LUPRI,'(1X,A,I5)') ' MAXSHL ', MAXSHL
         WRITE (LUPRI,'(1X,A,(5I5))') ' NHKTSH ',(NHKTSH(I),I=1,MAXSHL)
         WRITE (LUPRI,'(1X,A,(5I5))') ' KHKTSH ',(KHKTSH(I),I=1,MAXSHL)
         WRITE (LUPRI,'(1X,A,(5I5))') ' MHKTSH ',(MHKTSH(I),I=1,MAXSHL)
         WRITE (LUPRI,'(1X,A,(5I5))') ' NUCOSH ',(NUCOSH(I),I=1,MAXSHL)
         WRITE (LUPRI,'(1X,A,(5I5))') ' NORBSH ',(NORBSH(I),I=1,MAXSHL)
         WRITE (LUPRI,'(1X,A,(5I5))') ' NSTRSH ',(NSTRSH(I),I=1,MAXSHL)
         WRITE (LUPRI,'(1X,A,(5I5))') ' IORBSB ',(IORBSB(I),I=1,NTOTAL)
         WRITE (LUPRI,'(1X,A,(5I5))') ' IBLOCK ',(IBLOCK(I),I=1,NTOTAL)
         IF (BIGVEC) THEN
            WRITE (LUPRI,'(1X,A,(10I5))')
     *            ' NCNTSH ',(NCNTSH(I),I=1,MAXSHL)
         END IF
         DO 2000 I = 1, MAXSHL
            WRITE (LUPRI,'(1X,A,I5)') ' JSTRSH for ISHELL =',I
            WRITE (LUPRI,'(1X,15I5)') (JSTRSH(I,J),J=1,NSETSH(I))
 2000    CONTINUE
         DO 2001 I = 1, MAXSHL
            WRITE (LUPRI,'(1X,A,I5)') ' NPRIMS for ISHELL =',I
            WRITE (LUPRI,'(1X,15I5)') (NPRIMS(I,J),J=1,NSETSH(I))
 2001    CONTINUE
         DO 2002 I = 1, MAXSHL
            WRITE (LUPRI,'(1X,A,I5)') ' CENTSH for ISHELL =',I
            WRITE (LUPRI,'(1X,3F20.10)')
     *             CENTSH(I,1),CENTSH(I,2),CENTSH(I,3)
 2002    CONTINUE
         DO 2003 I = 1, MAXSHL
            WRITE (LUPRI,'(1X,A,I5)') ' IORBSH for ISHELL =',I
            WRITE (LUPRI,'(1X,15I5)') (IORBSH(I,J),J=1,NORBSH(I))
 2003    CONTINUE
      END IF
      RETURN
      END
