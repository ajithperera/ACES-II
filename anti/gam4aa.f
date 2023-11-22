      SUBROUTINE GAM4AA(GAMMA,BUF,BUF1,BUF2,GOFF,DOIRR,LENGAM)
C
C THIS ROUTINE DRIVES THE FORMATION OF THE ABCD SYMMETRY TYPE 
C  GAMMA LISTS FOR SPIN CASES AAAA AND BBBB AND THEN WRITES THEM 
C  TO THE BACKTRANSFORMATION MOINTS FILE.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION GAMMA(LENGAM),BUF(*),X,X1,X2,HALF
      DOUBLE PRECISION BUF1(*),BUF2(*)
      DOUBLE PRECISION FABIJ,FAIBJ,FABCD,FIJKL,FABCI,FIJKA
      INTEGER GOFF(8,8),IOFFSL(8),IOFFSR(8),DOIRR(42,4),IDID(8)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      LOGICAL ROHF
      LOGICAL TDA,EOM
      COMMON /EXCITE/ TDA,EOM
      COMMON /METH/ MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /MOPOPS/ MOPOP(8)
      COMMON /FACTORS/ FABIJ,FABCD,FIJKL,FAIBJ,FABCI,FIJKA
      COMMON /REF/ ROHF
C
      DATA HALF /0.5D0/
C
      INDXF(I,J,N)=I+(J-1)*N
      INDXT(I,J)  =I+(J*(J-1))/2
      NNP1O2(I)   =(I*(I+1))/2
C
      DO 100 ISPIN=1,2
       IF(TDA)THEN
        CALL G4AA4(GAMMA,BUF,GOFF,DOIRR,LENGAM,ISPIN,FAIBJ)
       ELSE
        CALL ZERO(GAMMA,LENGAM)
        CALL G4AA3(GAMMA,BUF,GOFF,DOIRR,LENGAM,ISPIN,FABIJ)
        IF(MBPT2)GOTO 2000
        CALL G4AA1(GAMMA,BUF,GOFF,DOIRR,LENGAM,ISPIN,FIJKL)
        CALL G4AA6(GAMMA,BUF,BUF1,BUF2,GOFF,DOIRR,LENGAM,ISPIN,FABCD)
        CALL G4AA4(GAMMA,BUF,GOFF,DOIRR,LENGAM,ISPIN,FAIBJ)
        IF((MBPT3.AND.(.NOT.ROHF)).OR.CCD)GOTO 2000
        CALL G4AA2(GAMMA,BUF,GOFF,DOIRR,LENGAM,ISPIN,FIJKA)
        CALL G4AA5(GAMMA,BUF,GOFF,DOIRR,LENGAM,ISPIN,FABCI)
       ENDIF
C      
2000   CONTINUE
       LIST1=9
       LIST2=1+(ISPIN-1)*50
       IOFF=1
       CALL SSCAL(LENGAM,HALF,GAMMA,1)
       DO 5000 IRREPDO=2,NIRREP 
        DO 5100 IRREPD=1,NIRREP
         IRREPC=DIRPRD(IRREPD,IRREPDO)
         IF(IRREPC.LT.IRREPD)GOTO 5100
         IBOT=MAX(IRREPC,IRREPD)+1
         CALL IZERO(IDID,NIRREP)
         DO 5200 IRREPTMP=IBOT,NIRREP
          IRREPA=DIRPRD(IRREPTMP,IRREPDO)
          IRREPB=MIN(IRREPTMP,IRREPA)
          IRREPA=MAX(IRREPTMP,IRREPA)
          IF(MAX(IDID(IRREPA),IDID(IRREPB)).NE.0)GOTO 5200
          IDID(IRREPA)=1
          IDID(IRREPB)=1
          NUMDIS=MOPOP(IRREPC)*MOPOP(IRREPD)
          DISSIZ=MOPOP(IRREPA)*MOPOP(IRREPB)
          CALL PUTLST(GAMMA(IOFF),1,NUMDIS,1,LIST1,LIST2)
          IOFF=IOFF+DISSIZ*NUMDIS
          LIST2=LIST2+1
5200     CONTINUE
5100    CONTINUE
5000   CONTINUE 
100   CONTINUE
      RETURN
      END
