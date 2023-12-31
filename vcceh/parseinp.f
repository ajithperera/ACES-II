      SUBROUTINE PARSEINP(NBAS,SCR,ISCR,IMAPMO,IRREPORB,ILOC)
C
C READS EXCITATION INPUT GIVEN AS DOMINANT SINGLE EXCITATION
C IN GENERAL, THESE ARE IN THE FORM
C
C  2*
C    A 3 15 
C    B 11 103 
C
C WHICH MEANS THAT TWO ROOTS WILL BE SEARCHED FOR - ONE PRINCIPALLY
C DESCRIBED AS AN EXCITATION FROM ALPHA ORBITAL 3 TO ALPHA ORBITAL 15 
C WITH THE BEING BETA ORBITAL 11 TO BETA ORBITAL 103.  THE ORBITAL 
C NUMBERING USED HERE CORRESPONDS TO THE NUMBERING SCHEME OF ORBITALS 
C IN WHICH THE EIGENVALUES ARE SORTED FROM SMALLEST TO LARGEST.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP,VRT,DIRPRD
      CHARACTER*80 STRING
      CHARACTER SPIN,SPLABEL(2)
      DIMENSION SCR(NBAS,*),ISCR(NBAS*NBAS,*),IMAPMO(NBAS,*)
      DIMENSION IRREPORB(NBAS,*)
      DIMENSION IOFFO(8,2),IOFFV(8,2),ILOC(*)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/FLAGS/IFLAGS(100)
      COMMON/CALCINFO/NROOT(8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON/GUESS2/IMAP(100,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      DATA SPLABEL /'A','B'/
C
      INDX(I,J,N)=I+(J-1)*N
C
      IUHF=MIN(IFLAGS(11),1)
      CALL IZERO(NROOT,8)
      CALL IZERO(IMAP,800)
C
      WRITE(6,1000)
1000  FORMAT(T3,'@PARSEINP-I, Input particle-hole excitations used ',
     &          'as initial guesses.')
  
C
C READ IN EIGENVALUES, SORT THEM AND HOLD MAPPING VECTOR RELATING
C SYMMETRY ORDER AND UNPACKED ORDER
C
      DO 100 ISPIN=1,1+IUHF
       DO 101 I=1,NBAS
        IMAPMO(I,ISPIN)=I
101    CONTINUE
100   CONTINUE
      CALL GETREC(20,'JOBARC','SCFEVALA',IINTFP*NBAS,SCR)
      CALL GETREC(20,'JOBARC','IRREPALP',NBAS,IRREPORB(1,1))
      CALL PIKSR2(NBAS,SCR,IMAPMO(1,1))
      IF(IUHF.NE.0)THEN
       CALL GETREC(20,'JOBARC','SCFEVALB',IINTFP*NBAS,SCR(1,2))
       CALL GETREC(20,'JOBARC','IRREPBET',NBAS,IRREPORB(1,2))
       CALL PIKSR2(NBAS,SCR(1,2),IMAPMO(1,2))
      ENDIF
C
C CONSTRUCT SYMMETRY VECTOR FOR THIS IRREP SO THAT WE CAN ULTIMATELY
C GET THE OFFSET INTO THE AI OR ai VECTOR FROM THE INDIVIDUAL VALUES
C OF A AND I 
C
C
      DO 1 ISPIN=1,1+IUHF
       IOFFO(1,ISPIN)=0
       IOFFV(1,ISPIN)=0
       DO 2 IRREP=1,NIRREP-1
        IOFFO(IRREP+1,ISPIN)=IOFFO(IRREP,ISPIN)+POP(IRREP,ISPIN)
        IOFFV(IRREP+1,ISPIN)=IOFFV(IRREP,ISPIN)+VRT(IRREP,ISPIN)
2      CONTINUE
1     CONTINUE
C
       DO 10 IRREPX=1,NIRREP
        ITHRU=0
        DO 15 ISPIN=1,1+IUHF
         DO 11 IRREPI=1,NIRREP
          IRREPA=DIRPRD(IRREPI,IRREPX)
          NUMI=POP(IRREPI,ISPIN)
          NUMA=VRT(IRREPA,ISPIN)
          DO 12 INDI=1,NUMI
           DO 13 INDA=1,NUMA
            ITHRU=ITHRU+1
            INDABSI=INDI+IOFFO(IRREPI,ISPIN)
            INDABSA=INDA+IOFFV(IRREPA,ISPIN)
            ISCR(INDX(INDABSA,INDABSI,NVRTO(ISPIN)),ISPIN)=ITHRU
13         CONTINUE
12        CONTINUE
11       CONTINUE
15      CONTINUE
10     CONTINUE
C
C NOW PARSE THE INPUT
C
      WRITE(6,2000)
      WRITE(6,2001)
      WRITE(6,2000)
c500   READ(30,'(A1,2I5)',END=999)SPIN,INDEXI,INDEXA
500   READ(30,*,END=999)ISPIN,INDEXI,INDEXA
c      IF(SPIN.EQ.'A')ISPIN=1
c      IF(SPIN.EQ.'B')ISPIN=2
      IF(ISPIN.EQ.0.AND.INDEXI.EQ.0.AND.INDEXA.EQ.0)GOTO 999
      IPACKI=IMAPMO(INDEXI,ISPIN)
      IPACKA=IMAPMO(INDEXA,ISPIN)-NOCCO(ISPIN)
      ISYMI =IRREPORB(IPACKI,ISPIN)
      ISYMA =IRREPORB(IPACKA+NOCCO(ISPIN),ISPIN)
      ISYMAI=DIRPRD(ISYMA,ISYMI)
      INDXAI=INDX(IPACKA,IPACKI,NVRTO(ISPIN))
      NROOT(ISYMAI)=NROOT(ISYMAI)+1
      IMAP(NROOT(ISYMAI),ISYMAI)=ISCR(INDXAI,ISPIN)
c     &                          +(ISPIN-1)*IRPDPD(ISYMAI,9)
      WRITE(6,2002)IPACKI,SCR(INDEXI,ISPIN),IPACKA+NOCCO(ISPIN),
     &             SCR(INDEXA,ISPIN),ISYMAI,SPLABEL(ISPIN)
      GOTO 500
C
999   CONTINUE
      WRITE(6,2000)
C
      RETURN
2000  FORMAT(71('-'))
2001  FORMAT(T8,'Hole orbital',T32,'Particle orbital',T53,'Transition ',
     &       /,
     &       T5,'Offset',T14,'Eigenvalue',T30,'Offset',T39,'Eigenvalue',
     &       T54,'Symmetry',T66,'Spin') 
2002  FORMAT(T7,I3,T13,F12.6,T32,I3,T38,F12.6,T57,I1,T68,A1)
      END
