      SUBROUTINE CABCI0(W,W1,BUF,IBUF,ICHAIN,
     &                  DISSIZ,NDBCK,NBKINT,NBUCK,
     &                  IRECL,ISPIN,NUMIRW,ISYM,
     &                  IPW,IPDIS,NWDIS,NMO,IRREPA)
C
C THIS ROUTINE CHAINS IN OFF OF THE PPPH SORT FILE AND
C  WRITES THEM TO MOINTS.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSIZ,DIRPRD,A,B,C,POP,VRT
      INTEGER AND
      CHARACTER*80 FNAME
      PARAMETER (LUSRT=15)
      DIMENSION W(DISSIZ*NDBCK),BUF(NBKINT),IBUF(NBKINT)
      DIMENSION ICHAIN(NBUCK),NUMIRW(2*NIRREP),ISYM(1),
     &          IPW(8),IPDIS(8),IPDSZ(8),W1(DISSIZ)
      DIMENSION IRREPA(NMO,2)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,
     &              NF2AA,NF2BB
      COMMON /SHIFT/ ISHIFT,NDRGEO
C
C STATEMENT FUNCTIONS
C
      IUPKI(IX)=AND(IX,IALONE)
      IUPKJ(IX)=AND(ISHFT(IX,-IBITWD),IALONE)
      IUPKK(IX)=AND(ISHFT(IX,-2*IBITWD),IALONE)
      IUPKL(IX)=AND(ISHFT(IX,-3*IBITWD),IALONE)
      INDX(I,J,N)=I+(J-1)*N
      IEXTI(IX,ISIZE)=1+INT((FLOAT(IX)-0.1)/FLOAT(ISIZE))
      IEXTC(IX,ISIZE)=IX-(IEXTI(IX,ISIZE)-1)*ISIZE
      CALL GFNAME('SRTFIL  ',FNAME,ILENGTH)
      OPEN(UNIT=LUSRT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=IRECL)
      NOCC=NOCCO(ISPIN)
      NVRT=NVRTO(ISPIN)
C
C     UPDATE POINTERS FOR WRITING
C
      IDIS=0
      IRREPW=1
      NDUMP=0
      NLEN=(NVRT*(NVRT-1))/2
      DO 10 IBK=1,NBUCK
       ISTART=(IBK-1)*NDBCK+1
       DO 1 IRREP=1,NIRREP
       IRREPST=IRREP
       IF(IPDIS(IRREP).GE.ISTART) THEN
        GO TO 2
       ENDIF
1      CONTINUE
       IRREPST=IRREPST+1
2      CONTINUE
       IRREPST=IRREPST-1
       IOFFST1=IPW(IRREPST)
       NNST=NUMIRW(IRREPST+NIRREP)
       IOFFST=IOFFST1+NNST*(ISTART-1-IPDIS(IRREPST))
       CALL ZERO(W,DISSIZ*NDBCK)
5      READ(LUSRT,REC=ICHAIN(IBK))BUF,IBUF,NUT,ICHAN
       ICHAIN(IBK)=ICHAN
       DO 20 INTGRL=1,NUT
        A=IUPKI(IBUF(INTGRL))
        B=IUPKJ(IBUF(INTGRL))
        C=IUPKK(IBUF(INTGRL))
        I=IUPKL(IBUF(INTGRL))
        IND1=INDX(A,B,NVRT)
        IND2=INDX(C,I,NVRT)
        IRREPAB=DIRPRD(IRREPA(A+NOCC,ISPIN),IRREPA(B+NOCC,ISPIN))
        IOFF=IPW(IRREPAB)
        NN=NUMIRW(IRREPAB+NIRREP)
        IADR=IOFF+ISYM(IND1+NWDIS)
     &       +NN*(ISYM(IND2)-1-IPDIS(IRREPAB))-IOFFST
        W(IADR)=BUF(INTGRL)
20     CONTINUE
       IF(ICHAIN(IBK).NE.0)GOTO 5

C IF WE REACH THIS POINT, THE ENTIRE BUCKET IS IN CORE.  ANTISYMMETRIZE
C  THE DISTRIBUTIONS AND DUMP THEM ONE AT A TIME TO DISK.  IF THE REFERE
C  FUNCTION IS RHF, THEN USE THE W ARRAY BEFORE ANTISYMMETRIZATION TO
C  WRITE THE AB INTEGRAL LISTS.
         IOFFW = 1
         IF (IFLAGS(11).GE.1) THEN
            DO ID = 1, NDBCK
28             CONTINUE
               IDIS = IDIS + 1
               IF (IDIS.GT.NUMIRW(IRREPW)) THEN
C                 SWITCH TO NEXT IRREP OR STOP
                  IRREPW = IRREPW + 1
                  IF (IRREPW.GT.NIRREP) GO TO 35
                  IDIS = 0
                  GO TO 28
               END IF
               CALL ASSYM(IRREPW,VRT(1,ISPIN),1,1,W1,W(IOFFW))
               CALL PUTLST(W1,IDIS,1,2,IRREPW,ISPIN+26+ISHIFT) 
               IOFFW = IOFFW + NUMIRW(IRREPW+NIRREP)
            END DO
         ELSE
            DO ID = 1, NDBCK
29             CONTINUE
               IDIS = IDIS + 1
               IF (IDIS.GT.NUMIRW(IRREPW)) THEN
                  IRREPW = IRREPW + 1
                  IF (IRREPW.GT.NIRREP) GO TO 35
                  IDIS = 0
                  GO TO 29
               END IF
               CALL PUTLST(W(IOFFW),IDIS,1,2,IRREPW,30+ISHIFT)
               IOFFW = IOFFW + NUMIRW(IRREPW+NIRREP)
            END DO
         END IF
35       CONTINUE

C     FINISHED PROCESSING THIS BUCKET.  GO BACK AND PICK UP ANOTHER ONE.
10    CONTINUE

      CLOSE(UNIT=LUSRT,STATUS='DELETE')

      WRITE(6,100)NDUMP
100   FORMAT(T3,'@CABCI0-I, Wrote ',I8,' logical records ',
     &          'of ordered PPPH integrals.')

      RETURN
      END

