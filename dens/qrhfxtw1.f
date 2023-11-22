      SUBROUTINE QRHFXTW1(XOVA,XOVB,Z1IA,Z1IB,
     &                    IRREPZ,ICORE,MAXCOR)
C
C  THIS ROUTINE FORMS CONTRIBUTIONS TO THE MODIFIED X
C  INTERMEDIATE REQUIRED FOR QRHF-CCSD GRADIENTS. 
C
C  THIS HANDLES CASES WHERE THERE ARE ONLY QRHF_P TYPE
C  OPEN SHELL ORBITALS.
C
C   XTW(A,I) = X(A,I) - SUM Z(1,M) [<1M||AI> + <AM||1I>]
C                       1,M
C
C  HERE A,I ARE THE VIRTUAL AND OCCUPIED ORBITALS IN THE RHF
C       REFERENCE,
C       1,M ARE THE OPEN-SHELL AND DOUBLY OCCUPIED ORBITALS
C       IN THE QRHF REFERENCE
C
C  THE SPIN SYMMETRY OF THE ORBITALS IS USED, AS ONLY ONE
C  SPIN CASE OF THE Z*A CONTRACTION IS EVALUATED.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYW,DISFUL,POP,VRT,POPRHF,VRTRHF
      INTEGER POPDOC,VRTDOC
      DIMENSION XOVA(1),XOVB(1),ICORE(MAXCOR)
      DIMENSION Z1IA(1),Z1IB(1)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2) 
      COMMON/QRHFINF/POPRHF(8),VRTRHF(8),NOSH1(8),NOSH2(8),
     &               POPDOC(8),VRTDOC(8),NAI,N1I,NA2,
c&line mod
     &               NUMISCF,NUMASCF,ISPINP,ISPINM,IQRHF
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/INFO/NOCCO(2),NVRTO(2)
C
      INDXF(I,J,N)=I+(J-1)*N
C
      DATA ZERO,ONE,ONEM,TWO /0.D0,1.D0,-1.D0,2.D0/
C
C SET MAXIMUM CORE SIZE
C
      MXCOR=MAXCOR
C
      NDOCC=POPDOC(IRREPZ)
      NDVRT=VRTDOC(IRREPZ)
C
C TERM I:  SUM Z(1,M) * [<MI||1A> + <1I||MA>]  Z - AA   W - AAAA
C          1,M
C
c&start modif
      IF(ISPINP.EQ.2) THEN
C
C THE OPEN SHELL ORBITAL IS OCCUPIED IN THE ALPHA SPIN CASE,
C HENCE THE <1A||MI> INTEGRALS ARE OF TYPE <IJ||KA>
C
        LISTW=7
      ELSE
C
C THE OPEN SHELL ORBITAL IS OCCUPIED IN THE BETA SPIN CASE,
C HENCE THE <1A||MI> INTEGRALS ARE OF TYPE <ij||ka>
C
         LISTW=8
      ENDIF
C
C LOOP OVER IRREPS OF MI AND 1A
C
      DO 1000 IRREP=1,NIRREP
C
C GET THE SIZE OF THE EXPANDED INTEGRAL LIST
C
       NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
       DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
c&line mod
       DISFUL=IRPDPD(IRREP,23-ISPINP)
C
C ALOOCATE MEMORY
C
       I000=1
       I010=I000+IINTFP*DISFUL*NUMSYW
       I020=I010+IINTFP*NAI
       IF(I020.GE.MAXCOR)CALL INSMEM('QRHFXTW1',I020,MAXCOR)
C
C GET THE INTEGRALS 
C
       CALL GETLST(ICORE(I000),1,NUMSYW,2,IRREP,LISTW)
C
C EXPAND THE OCC.-OCC. BLOCK
C
c&line mod
       CALL SYMEXP2(IRREP,POP(1,3-ISPINP),DISFUL,DISSYW,NUMSYW,
     &              ICORE(I000),ICORE(I000))
C
C IRREPX IS THE IRREP WHICH BELONGS TO I AND A IN THE FOLLOWING
C
       IRREPX=DIRPRD(IRREP,IRREPZ)
       NUMI=POPRHF(IRREPX)
       NUMA=VRTRHF(IRREPX)
C
C  GET OFFSETS FOR 1I AND MA BLOCKS
C
       IOFFR=0
       IOFFL=0
       IOFFX=1
       DO 210 IRREPJ=1,IRREPX-1
        IRREPI=DIRPRD(IRREPJ,IRREP)
c&two lines mod
        IOFFL=IOFFL+POP(IRREPI,3-ISPINP)*POP(IRREPJ,3-ISPINP)
        IOFFR=IOFFR+POP(IRREPI,3-ISPINP)*VRT(IRREPJ,3-ISPINP)
        IOFFX=IOFFX+POPRHF(IRREPJ)*VRTRHF(IRREPJ)
210    CONTINUE
C
C   LOOP NOW OVER M (ONLY OVER DOUB. OCCS.) AND 1 (ALL TYPE I
C   OPEN-SHELL ORBITALS)
C
       IPOS=1
       DO 220 M=1,NDOCC
        DO 220 IOS1=1,NOSH1(IRREPZ)
C
         CALL IZERO(ICORE(I010),IINTFP*NAI)
C
C      OFFSETS FOR 1A AND MA (RIGHT HAND SIDE) :
C
         IOFF1AR=IOFFR+NDOCC+IOS1
         IOFFMAR=IOFFR+M
C
C COPY ONE A,I BLOCK TO SCRATCH
C
         IF(MIN(NUMI,NUMA).NE.0)THEN
          DO 230 I=1,NUMI
C
C      OFFSETS FOR 1I AND MI (LEFT HAND SIDE) :
C
c&two lines mod
           IOFF1IL=IOFFL+INDXF(NDOCC+IOS1,I,POP(IRREPZ,3-ISPINP))
           IOFFMIL=IOFFL+INDXF(M,I,POP(IRREPZ,3-ISPINP))
C
C  NOW COPY OVER BOTH REQUIRED PIECES 
C
           IOFF1IMA=I000+IINTFP*(INDXF(IOFF1IL,IOFFMAR,DISFUL)-1)
           IOFFMI1A=I000+IINTFP*(INDXF(IOFFMIL,IOFF1AR,DISFUL)-1)
           IOFFAI  =I010+IINTFP*(INDXF(1,I,NUMA)-1)
c&four lines mod
           CALL SAXPY(NUMA,ONE,ICORE(IOFF1IMA),
     &                DISFUL*POP(IRREPZ,3-ISPINP),ICORE(IOFFAI),1)
           CALL SAXPY(NUMA,ONE,ICORE(IOFFMI1A),
     &                DISFUL*POP(IRREPZ,3-ISPINP),ICORE(IOFFAI),1)
230       CONTINUE
C
C  MULTIPLY THE SCRATCH WITH Z1A OR Z1B
C
          LENGTH=NUMI*NUMA
          CALL SAXPY(LENGTH,Z1IA(IPOS),ICORE(I010),1,
     &               XOVA(IOFFX),1)
          CALL SAXPY(LENGTH,Z1IB(IPOS),ICORE(I010),1,
     &               XOVB(IOFFX),1)
         ENDIF
         IPOS=IPOS+1
C
220     CONTINUE
C
1000   CONTINUE
C
C TERM II:
C
C ALPHA-BETA CONTRIBUTION FROM Z(1,m) :
C
C         Z(1,m) * <A1|Im> 
C
C         FACTOR TWO SINCE <1A||mI> = <1I||mA>
C
c&lines modif
      IF(ISPINP.EQ.2) THEN
C
C THE OPEN SHELL ORBITAL IS VIRTUAL FOR SPIN CASE BETA, SO
C THE <A1|Im> INTEGRALS ARE OF TYPE <Ab|Ij>.  LIST 18
C ORDERS THEM (AI,1m) WHICH IS CONVENIENT FOR OUR PURPOSES.
C
         LISTW=18
      ELSE
C
C THE OPEN SHELL ORBITAL IS VIRTUAL FOR SPIN CASE ALPHA, SO
C THE <A1|Im> INTEGRALS ARE OF TYPE <Ab|Ij>.  LIST 17
C ORDERS THEM (ai,1M) WHICH IS CONVENIENT FOR OUR PURPOSES.
C
         LISTW=17
      ENDIF
C
C  WE NEED TO CONSIDER ONLY IRREP 1 OF THE INTEGRAL LIST
C
C CALCULATE OFFSET TO WHERE THE 1m PAIRS WITH 
C  IRREP(1)=IRREP(m)=IRREPZ BEGIN.
C
c&end modif
       IOFFDIS=0
       DO 3000 IRREP=1,IRREPZ-1
c&line mod
        IOFFDIS=IOFFDIS+VRT(IRREP,ISPINP)*POP(IRREP,ISPINP)
3000   CONTINUE
C
C  LOOP OVER m AND 1. 
C
       IPOS=1
       DO 3100 M=1,NDOCC
        DO 3200 IOS1=1,NOSH1(IRREPZ)
C
C  CALCULATE THE DISTRIBUTION NUMBER WE WANT HERE AND READ
C  IN THE CORRESPONDING W(AI,1m) ARRAY
C
c&line mod
         IDIS=IOFFDIS+INDXF(IOS1,M,VRT(IRREPZ,ISPINP))
         CALL GETLST(ICORE,IDIS,1,1,1,LISTW)
C
C  NOW SCALE ALL ELEMENTS OF W(AI) WITH Z(1,m) AND ADD THIS
C  TO EXISTING X PIECES
C
c&two lines mod
         CALL SAXPY(NT(IQRHF),TWO*Z1IB(IPOS),ICORE,1,XOVA,1)
         CALL SAXPY(NT(IQRHF),TWO*Z1IA(IPOS),ICORE,1,XOVB,1)
         IPOS=IPOS+1
C
3200    CONTINUE
3100   CONTINUE
C
      RETURN
      END