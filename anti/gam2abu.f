      SUBROUTINE GAM2ABU(GAMMA,BUF,BUF1,BUF2,GOFF,LENGAM,ALASKA)
C
C THIS ROUTINE FORMS THE AABB SYMMETRY TYPE AO GAMMA LISTS
C  FOR SPIN CASE AB AND WRITES THEM TO THE BACKTRANSFORMATION
C  MOINTS FILE FOR UHF CALCULATIONS.
C
C NOTE THAT AABB IN MULLIKEN CORRESPONDS TO ABAB IN DIRAC.  CONSEQUENTLY,
C  ALL ELEMENTS OF GAMMA WHICH ARE PROCESSED HERE COME FROM DPD LISTS
C  FOR IRREPS 2-h, SINCE A.NE.B.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION GAMMA(LENGAM),BUF(*),X,X1,X2,HALF
      DOUBLE PRECISION FABIJ,FAIBJ,FABCD,FIJKL,FABCI,FIJKA
      DOUBLE PRECISION BUF1(*),BUF2(*)
      INTEGER GOFF(8,8),IOFFSR(8),IOFFSL(8)
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
      DATA HALF /0.5D0/
C
      INDXF(I,J,N)=I+(J-1)*N
      INDXT(I,J)  =I+(J*(J-1))/2
      NNP1O2(I)   =(I*(I+1))/2
C
C SPIN CASE AABB
C
      CALL ZERO(GAMMA,LENGAM)
      IF(TDA)THEN
       CALL G2AB4A(GAMMA,BUF,GOFF,LENGAM,FAIBJ)
      ELSE
       CALL G2AB3A(GAMMA,BUF,GOFF,LENGAM,FABIJ)
       IF(MBPT2)GOTO 2000
       CALL G2AB1A(GAMMA,BUF,GOFF,LENGAM,FIJKL)
       CALL G2AB6A(GAMMA,BUF,BUF1,BUF2,GOFF,LENGAM,FABCD)
       CALL G2AB4A(GAMMA,BUF,GOFF,LENGAM,FAIBJ)
       IF(CCD.OR.(MBPT3.AND.(.NOT.ROHF)))GOTO 2000
       CALL G2AB2A(GAMMA,BUF,GOFF,LENGAM,FIJKA)
       CALL G2AB5A(GAMMA,BUF,GOFF,LENGAM,FABCI)
2000   CONTINUE
      ENDIF
      LIST1=2
      LIST2=1
      CALL SSCAL(LENGAM,HALF,GAMMA,1)
      IOFF=1
      DO 5000 IRREPR=2,NIRREP
       DO 5100 IRREPL=1,IRREPR-1
        NUMDIS=NNP1O2(MOPOP(IRREPR))
        DISSIZ=NNP1O2(MOPOP(IRREPL))
        IF(MOEQAO .OR. ALASKA)THEN
         CALL PUTLST(GAMMA(IOFF),1,NUMDIS,1,LIST1,LIST2)
         IOFF=IOFF+DISSIZ*NUMDIS
        ELSE
         DO 5101 IDIS=1,NUMDIS
          CALL PUTLST(GAMMA(IOFF),IDIS,1,1,LIST1,LIST2)
          IOFF=IOFF+DISSIZ
5101     CONTINUE
        ENDIF
        LIST2=LIST2+1
5100   CONTINUE
5000  CONTINUE
C
C SPIN CASE BBAA
C
      CALL ZERO(GAMMA,LENGAM)
      IF(TDA)THEN
       CALL G2AB4B(GAMMA,BUF,GOFF,LENGAM,FAIBJ)
      ELSE
       CALL G2AB3B(GAMMA,BUF,GOFF,LENGAM,FABIJ)
       IF(MBPT2)GOTO 3000
       CALL G2AB1B(GAMMA,BUF,GOFF,LENGAM,FIJKL)
       CALL G2AB6B(GAMMA,BUF,BUF1,BUF2,GOFF,LENGAM,FABCD)
       CALL G2AB4B(GAMMA,BUF,GOFF,LENGAM,FAIBJ)
       IF(CCD.OR.(MBPT3.AND.(.NOT.ROHF)))GOTO 3000
       CALL G2AB2B(GAMMA,BUF,GOFF,LENGAM,FIJKA)
       CALL G2AB5B(GAMMA,BUF,GOFF,LENGAM,FABCI)
3000   CONTINUE
      ENDIF
      CALL SSCAL(LENGAM,HALF,GAMMA,1)
      LIST1=2
      LIST2=51
      IOFF=1
      DO 6000 IRREPR=2,NIRREP
       DO 6100 IRREPL=1,IRREPR-1
        NUMDIS=NNP1O2(MOPOP(IRREPR))
        DISSIZ=NNP1O2(MOPOP(IRREPL))
        CALL PUTLST(GAMMA(IOFF),1,NUMDIS,1,LIST1,LIST2)
        IOFF=IOFF+DISSIZ*NUMDIS
        LIST2=LIST2+1
6100   CONTINUE
6000  CONTINUE
      RETURN
      END
