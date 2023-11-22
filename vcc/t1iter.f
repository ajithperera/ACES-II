      SUBROUTINE T1ITER(ICORE,MAXCOR,IUHF,ET1)
C
C SOLVES FOR THE T1[1] VECTOR ITERATIVELY FOR NONCANONICAL CASES.
C  EACH IRREP IS HANDLED INDIVIDUALLY
C
C D1 T1(a,i) = f(ai) + SUM t(e,i)*f(ea) - SUM t(a,m) * F(mi)
C                       e                  m
C
CEND
      IMPLICIT INTEGER (A-Z)
      CHARACTER*2 SPCASE(2)
      DOUBLE PRECISION ONE,ONEM,ZILCH,X,Z,SDOT,ET1,FNDLRGAB
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      DATA SPCASE /'AA','BB'/
      DATA ONE  /1.0/   
      DATA ONEM /-1.0/
      DATA ZILCH /0.0/
      ET1=0.0
      DO 5 ISPIN=1,1+IUHF
       WRITE(6,1000)SPCASE(ISPIN)
1000   FORMAT(T3,'@T1ITER-I, Beginning iterative solution for ',A2,
     &           ' amplitudes.')
       I000=1
       I010=I000+NT(ISPIN)*IINTFP
       I020=I010+NT(ISPIN)*IINTFP
       I030=I020+NT(ISPIN)*IINTFP
       I040=I030+NFEA(ISPIN)*IINTFP
       I050=I040+NFMI(ISPIN)*IINTFP
C
C LOAD INITIAL GUESS AMPLITUDES (NO F*T1 CONTRIBUTION), FEA AND FMI
C  ELEMENTS
C
       CALL GETLST(ICORE(I000),1,1,1,ISPIN,90)
       CALL GETLST(ICORE(I020),1,1,1,ISPIN+2,93)
       CALL GETLST(ICORE(I030),1,1,1,ISPIN+2,92)
       CALL GETLST(ICORE(I040),1,1,1,ISPIN+2,91) 
       CALL GETLST(ICORE(I050),1,1,1,9,63+ISPIN)
       CALL INVERS(ICORE(I050),NT(ISPIN))
       IOFFT1 =I000
       IOFFTAR=I010
       IOFFFAI=I020
       IOFFEA =I030
       IOFFMI =I040
       IOFFDEN=I050
       DO 10 IRREP=1,NIRREP
        ITER=0
C
C CALCULATE T(EI)*F(EA) CONTRIBUTION
C
        NVRT=VRT(IRREP,ISPIN)
        NPOP=POP(IRREP,ISPIN) 
        LENGTH=NVRT*NPOP
        IF(LENGTH.EQ.0)GOTO 2
1       CALL IZERO (ICORE(IOFFTAR),LENGTH*IINTFP)
        CALL XGEMM('N','N',NVRT,NPOP,NVRT,ONE,ICORE(IOFFEA),NVRT,
     &             ICORE(IOFFT1),NVRT,ONE,ICORE(IOFFTAR),NVRT)
C
C CALCULATE T(AM) * F(MI) CONTRIBUTION
C
        CALL XGEMM('N','N',NVRT,NPOP,NPOP,ONEM,ICORE(IOFFT1),NVRT,
     &             ICORE(IOFFMI),NPOP,ONE,ICORE(IOFFTAR),NVRT)
C
C ADD IN F(AI) PIECE AND DENOMINATOR WEIGHT
C
        CALL SAXPY(LENGTH,ONE,ICORE(IOFFFAI),1,ICORE(IOFFTAR),1)
        CALL VECPRD(ICORE(IOFFTAR),ICORE(IOFFDEN),ICORE(IOFFTAR),LENGTH)
C
C COMPARE WITH PREVIOUS ITERATE
C
        CALL SAXPY(LENGTH,ONEM,ICORE(IOFFTAR),1,ICORE(IOFFT1),1)
        X=FNDLRGAB(ICORE(IOFFT1),LENGTH)
        IF(X.GT.1.E-10)THEN
         ITER=ITER+1
c YAU : old
c        CALL ICOPY(LENGTH*IINTFP,ICORE(IOFFTAR),1,ICORE(IOFFT1),1)
c YAU : new
         CALL DCOPY(LENGTH,ICORE(IOFFTAR),1,ICORE(IOFFT1),1)
c YAU : end
         GOTO 1
        ELSE
         ITER=ITER+1
         WRITE(6,1001)IRREP,ITER
1001     FORMAT(T3,'Irrep ',I2,' amplitudes converged in ',I4,
     &         ' iterations.')
c YAU : old
c        CALL ICOPY(LENGTH*IINTFP,ICORE(IOFFTAR),1,ICORE(IOFFT1),1)
c YAU : new
         CALL DCOPY(LENGTH,ICORE(IOFFTAR),1,ICORE(IOFFT1),1)
c YAU : end
        ENDIF
2       IOFFDEN=IOFFDEN+LENGTH*IINTFP
        IOFFTAR=IOFFTAR+LENGTH*IINTFP
        IOFFT1 =IOFFT1 +LENGTH*IINTFP
        IOFFEA =IOFFEA +NVRT*NVRT*IINTFP
        IOFFMI =IOFFMI +NPOP*NPOP*IINTFP
        IOFFFAI=IOFFFAI+LENGTH*IINTFP
10     CONTINUE 
       CALL PUTLST(ICORE(I000),1,1,1,ISPIN,90)
       CALL GETLST(ICORE(I010),1,1,1,ISPIN+2,93)
       Z=SDOT(NT(ISPIN),ICORE(I000),1,ICORE(I010),1) 
       ET1=ET1+Z
       WRITE(6,1002)SPCASE(ISPIN),Z
1002   FORMAT(T3,'Spin case ',A2,' energy contribution ',F15.10,'.')
5     CONTINUE
      RETURN
      END
