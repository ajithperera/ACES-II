      SUBROUTINE SXEVAL1(IRREPX,FD,E,SD,NUM,IANTI)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD
      DIMENSION FD(1),E(1),SD(1),NUM(8)
C
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
      DATA HALFM /-0.5D0/
C
      IOFFD=0
      IOFFER=0
C
      DO 1000 IRREPR=1,NIRREP
C
       IRREPL=DIRPRD(IRREPX,IRREPR)
C
       IOFFEL=0
       DO 1 IRREP=1,IRREPL-1
        IOFFEL=IOFFEL+NUM(IRREP)
1      CONTINUE
C
       IOFFD2=0
       DO 110 IRREP1=1,IRREPL-1
        IRREP2=DIRPRD(IRREPX,IRREP1)
        IOFFD2=IOFFD2+NUM(IRREP1)*NUM(IRREP2)
110    CONTINUE
C
       NR=NUM(IRREPR)
       NL=NUM(IRREPL)
C
       IF(IANTI.EQ.0) THEN
C
        DO 100 I=1,NR
C
         CALL SAXPY(NL,HALFM*E(IOFFER+I),SD(IOFFD2+I),NR,
     &              FD(IOFFD+(I-1)*NL+1),1)
C
100     CONTINUE
C
       ELSE

        DO 150 I=1,NR
C
         CALL SAXPY(NL,HALFM*E(IOFFER+I),SD(IOFFD+(I-1)*NL+1),
     &              1,FD(IOFFD+(I-1)*NL+1),1)
C
150     CONTINUE
C
       ENDIF
C
       DO 200 I=1,NL
C
        CALL SAXPY(NR,HALFM*E(IOFFEL+I),SD(IOFFD+I),NL,
     &             FD(IOFFD+I),NL)
200    CONTINUE
C
       IOFFD=IOFFD+NR*NL
       IOFFER=IOFFER+NR
C
1000  CONTINUE
C
      RETURN
      END
