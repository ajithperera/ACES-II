      SUBROUTINE GAM4AB(GAMMA,BUF,BUF1,BUF2,SCR1,SCR2,GOFF,
     &                  DOIRR,LENGAM,IUHF,ALASKA)
C
C THIS ROUTINE DRIVES THE FORMATION OF ABCD SYMMETRY TYPE AO
C  GAMMA LISTS FOR SPIN CASE AB AND WRITES THEM TO THE BACKTRANSFORMATION
C  MOINTS FILE.  THIS IS THE INCORE VERSION!
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION GAMMA(LENGAM),BUF(1),X,X1,X2 
      DOUBLE PRECISION FABIJ,FAIBJ,FABCD,FIJKL,FABCI,FIJKA
      DOUBLE PRECISION BUF1(*),BUF2(*),SCR1(*),SCR2(*)
      INTEGER GOFF(8,8,8,8),IDID(8),IOFFSL(8),DOIRR(42,4)
      INTEGER IOFFSR(8)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      LOGICAL ROHF,MOEQAO,ALASKA
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
      COMMON /SIZES/ MOEQAO
C
      INDXF(I,J,N)=I+(J-1)*N
      INDXT(I,J)  =I+(J*(J-1))/2
      NNP1O2(I)   =(I*(I+1))/2
C
      CALL ZERO(GAMMA,LENGAM)
C
C FIRST PROCESS THE G(Ab,Ij) LISTS
C
      IF(TDA)THEN
       CALL G4AB41(GAMMA,BUF,GOFF,DOIRR,LENGAM,FAIBJ)
       CALL G4AB42(GAMMA,BUF,GOFF,DOIRR,LENGAM,FAIBJ)
      ELSE
       CALL G4AB3(GAMMA,BUF,GOFF,LENGAM,FABIJ)
       IF(MBPT2)GOTO 2000
       CALL G4AB1(GAMMA,BUF,GOFF,DOIRR,LENGAM,FIJKL)
       CALL G4AB6(GAMMA,BUF,BUF1,BUF2,SCR1,SCR2,GOFF,
     &            DOIRR,LENGAM,FABCD,IUHF)
       CALL G4AB41(GAMMA,BUF,GOFF,DOIRR,LENGAM,FAIBJ)
       CALL G4AB42(GAMMA,BUF,GOFF,DOIRR,LENGAM,FAIBJ)
       IF(CCD.OR.(MBPT3.AND.(.NOT.ROHF)))GOTO 2000
       CALL G4AB2(GAMMA,BUF,GOFF,DOIRR,LENGAM,FIJKA)
       CALL G4AB5(GAMMA,BUF,GOFF,DOIRR,LENGAM,FABCI)
C
C WRITE THE ABCD GAMMAS TO DISK, IRREP BY IRREP
C
2000   CONTINUE
      ENDIF
      LIST1=4
      LIST2=1
      IOFF=1
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
         IF(MOEQAO .OR. ALASKA)THEN
          CALL PUTLST(GAMMA(IOFF),1,NUMDIS,1,LIST1,LIST2)
          IOFF=IOFF+DISSIZ*NUMDIS
         ELSE
C
C CODE BELOW REQUIRED FOR CASES WHERE NMO<NAO, SINCE NAO GAMMAS
C EVENTUALLY OVERWRITE LIST (4,1), WHICH HAS THE PHYSICAL DIMENSIONS
C OF THE CORRESPONDING AO GAMMA
C
          DO 5201 IDIS=1,NUMDIS
           CALL PUTLST(GAMMA(IOFF),IDIS,1,1,LIST1,LIST2)
           IOFF=IOFF+DISSIZ
5201      CONTINUE
         ENDIF        
         LIST2=LIST2+1
5200    CONTINUE
5100   CONTINUE
5000  CONTINUE 
      RETURN
      END