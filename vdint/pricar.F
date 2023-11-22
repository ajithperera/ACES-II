      SUBROUTINE PRICAR
C
C     This subroutine prints information about atomic coordinates
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
#include "baslims.par"
      CHARACTER*1 CHRXYZ(-3:3)
      CHARACTER*1 CHRSGN(-1:1)
C
      CHARACTER*4 NAME
      INTEGER NBASE(8), NSIGN(8)
C
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
      COMMON /TRANUC/ TRCTOS(MXCOOR,MXCOOR), TRSTOC(MXCOOR,MXCOOR)
C
      DATA CHRSGN /'-',' ','+'/
      DATA CHRXYZ /'z','y','x',' ','X','Y','Z'/
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
      OPEN(UNIT=81,FILE='GRD',STATUS='UNKNOWN',FORM='FORMATTED')
      REWIND(81)
C
cmn      WRITE (LUPRI,'(/)')
      NCOOR = 3*NUCDEP
C
C  Cartesian Coordinates
C
      CALL HEADER('Cartesian Coordinates',1)
      WRITE (LUPRI,'(1X,A,I3,//)') ' Total number of coordinates:',NCOOR
      ENERGY=0.D0
      WRITE(81,3000) NUCDEP,ENERGY
3000  FORMAT(I5,F20.10)
C 
      ICRX = 1
      ICRY = 2
      ICRZ = 3
      DO 100 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         NAME   = NAMEX(3*ICENT)(1:4)
         IF (MULT(MULCNT) .EQ. 1) THEN
            WRITE (LUPRI,'(1X,I3,3X,A,5X,A,3X,F15.10)')
     *         ICRX, NAME,CHRXYZ(-1), CORD(ICENT,1)
            WRITE (LUPRI,'(1X,I3,12X,A,3X,F15.10)')
     *         ICRY, CHRXYZ(-2), CORD(ICENT,2)
             WRITE (LUPRI,'(1X,I3,12X,A,3X,F15.10,/)')
     *         ICRZ, CHRXYZ(-3), CORD(ICENT,3)
            WRITE(81,4000)CHARGE(ICENT),(CORD(ICENT,I),I=1,3)
4000        FORMAT(4F20.10)
            ICRX = ICRX + 3
            ICRY = ICRY + 3
            ICRZ = ICRZ + 3
         ELSE
           JATOM = 0
            DO 200 ISYMOP = 0, MAXLOT
               IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
                  JATOM = JATOM + 1
                  WRITE (LUPRI,'(1X,I3,3X,A,I2,3X,A,3X,F15.10)')
     *               ICRX, NAME,JATOM,CHRXYZ(-1),
     *               PT(IBTAND(ISYTYP(1),ISYMOP))*CORD(ICENT,1)
                  WRITE (LUPRI,'(1X,I3,12X,A,3X,F15.10)')
     *               ICRY,CHRXYZ(-2),
     *               PT(IBTAND(ISYTYP(2),ISYMOP))*CORD(ICENT,2)
                  WRITE (LUPRI,'(1X,I3,12X,A,3X,F15.10,/)')
     *               ICRZ,CHRXYZ(-3),
     *               PT(IBTAND(ISYTYP(3),ISYMOP))*CORD(ICENT,3)
                  WRITE(81,4000) CHARGE(ICENT),PT(IBTAND(ISYTYP(1),
     *                           ISYMOP))*CORD(ICENT,1),PT(IBTAND
     *                           (ISYTYP(2),ISYMOP))*CORD(ICENT,2), 
     *                           PT(IBTAND(ISYTYP(3),ISYMOP))*CORD 
     *                           (ICENT,3) 
                  ICRX = ICRX + 3
                  ICRY = ICRY + 3
                  ICRZ = ICRZ + 3
               END IF
  200       CONTINUE
         END IF
  100 CONTINUE
C
C  Symmetry Coordinates
C
      IF (MAXLOP .GT. 0) THEN
         CALL HEADER('Symmetry Coordinates',1)
         WRITE (LUPRI,'(1X,A,8I3)')
     *          ' Number of coordinates in each symmetry: ',
     *           (NCRREP(I),I=0,MAXLOP)
         DO 300 ISYM = 0, MAXLOP
         IF (NCRREP(ISYM) .GT. 0) THEN
cmn            WRITE (LUPRI,'(//1X,A,I2/)') ' Symmetry', ISYM + 1
            DO 400 IATOM = 1, NUCIND
               DO 500 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,ISYM)
                  IF (ISCOOR .GT. 0) THEN
                     NB = 0
                     DO 600 I = 1, NCOOR
                        NSGN = NINT(TRCTOS(ISCOOR,I))
                        IF (NSGN .NE. 0) THEN
                           NB = NB + 1
                           NBASE(NB) = I
                           NSIGN(NB) = NSGN
                        END IF
  600                CONTINUE
                     WRITE (LUPRI,
     *                     '(1X,I3,3X,A,2X,A,3X,I2,7(2X,A,1X,I2))')
     *                     ISCOOR, NAMEX(3*IATOM)(1:4), CHRXYZ(-ICOOR),
     *                     NBASE(1), (CHRSGN(NSIGN(I)),NBASE(I),I=2,NB)
                  END IF
  500          CONTINUE
  400       CONTINUE
         END IF
  300    CONTINUE
      END IF
      RETURN
      END
