      SUBROUTINE T2ITER(ICORE,MAXCOR,IUHF,ET2)
C
C SOLVES FOR THE MBPT[2] ENERGY AN T2 VECTOR ITERATIVELY
C  FOR NONCANONICAL CASES.  EACH IRREP IS HANDLED INDIVIDUALLY
C
C D2 T2(ab,ij)=<ab||ij>+P(ab) SUM T2(ae,ij)*F(be)-P(ij) SUM T2(ab,mj)*F(mj)
C                              e                         m
C
CEND
      IMPLICIT INTEGER (A-Z)
      CHARACTER*2 SPCASE(3)
      DOUBLE PRECISION ONE,ONEM,ZILCH,X,HALF,FNDLRGAB,Z,SDOT,ET2
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      DATA SPCASE /'AA','BB','AB'/
      DATA ONE  /1.0/   
      DATA ONEM /-1.0/
      DATA HALF /0.5/
      DATA ZILCH /0.0/
      ET2=0.0
C
C SPIN CASES AA AND BB
C
      DO 5 ISPIN=1,1+IUHF
       WRITE(6,1000)SPCASE(ISPIN)
1000   FORMAT(T3,'@T2ITER-I, Beginning iterative solution for ',A2,
     &           ' amplitudes.')
       Z=0.0
       DO 10 IRREP=1,NIRREP
        ITER=0
        NT2DSZ=IRPDPD(IRREP,ISYTYP(1,13+ISPIN))
        NT2DIS=IRPDPD(IRREP,ISYTYP(2,13+ISPIN))
        IF(MIN(NT2DIS,NT2DSZ).EQ.0)GOTO 10
        NSPACE=MAX(IRPDPD(IRREP,18+ISPIN)*NT2DIS,IRPDPD(IRREP,20+ISPIN)
     &             *NT2DSZ)
        DSZFUL=MAX(IRPDPD(IRREP,18+ISPIN),IRPDPD(IRREP,20+ISPIN)) 
        I000=1
        I010=I000+NSPACE*IINTFP
        I020=I010+NSPACE*IINTFP
        I030=I020+NFEA(ISPIN)*IINTFP
        I040=I030+NFMI(ISPIN)*IINTFP
        I050=I040+NT2DSZ*NT2DIS*IINTFP
        CALL IZERO(ICORE,I050)
C
C READ IN VV AND OO FOCK MATRICES
C
        CALL GETLST(ICORE(I020),1,1,1,2+ISPIN,92)
        CALL GETLST(ICORE(I030),1,1,1,2+ISPIN,91)
C
C READ IN INITIAL T2 AMPLITUDES (EXCLUDING F CONTRIBUTION)
C
        CALL GETLST(ICORE(I010),1,NT2DIS,1,IRREP,43+ISPIN)
C
C NOW TRANSPOSE AND EXPAND TO FORM T2(I<J,AE) BEFORE CALCULATING
C  INCREMENT.
C
1       CALL TRANSP(ICORE(I010),ICORE(I000),NT2DIS,NT2DSZ)
        CALL SYMEXP(IRREP,VRT(1,ISPIN),NT2DIS,ICORE(I000))
C
C PERFORM MULTIPLICATION
C                                 
C                    Z(I<J,AB) = T2(I<JA,E) * T(E,B)
C
        IOFFT  =I000
        IOFFTAR=I010
        IOFFF  =I020
        DO 20 IRREPE=1,NIRREP
         IRREPA=DIRPRD(IRREPE,IRREP)
         NUMA=VRT(IRREPA,ISPIN)
         NUME=VRT(IRREPE,ISPIN)
         NUMB=NUME
         IF(MAX(NT2DIS*NUMA,NUMB,NUME).NE.0)THEN
          CALL XGEMM('N','N',NT2DIS*NUMA,NUMB,NUME,ONE,ICORE(IOFFT),
     &               NT2DIS*NUMA,ICORE(IOFFF),NUME,ZILCH,
     &              ICORE(IOFFTAR),NT2DIS*NUMA)
         ENDIF
         IOFFF=IOFFF+NUME*NUMB*IINTFP
         IOFFT=IOFFT+NT2DIS*NUMA*NUME*IINTFP
         IOFFTAR=IOFFTAR+NT2DIS*NUMA*NUME*IINTFP
20      CONTINUE
C
C NOW ANTISYMMETRIZE - Z(I<J,AB) = Z(I<J,AB) - Z(I<J,BA) AND THEN
C  TRANSPOSE TO Z(A<B,I<J). DO THIS FOR BOTH THE T2 VECTOR AND
C  THE INCREMENT (Z)
C
        CALL ASSYM2(IRREP,VRT(1,ISPIN),NT2DIS,ICORE(I010))
        CALL SSCAL(NT2DIS*NT2DSZ,HALF,ICORE(I010),1)
        CALL TRANSP(ICORE(I010),ICORE(I040),NT2DSZ,NT2DIS)
c YAU : old
c       CALL ICOPY(NT2DIS*NT2DSZ*IINTFP,ICORE(I040),1,ICORE(I010),1)
c YAU : new
        CALL DCOPY(NT2DIS*NT2DSZ,ICORE(I040),1,ICORE(I010),1)
c YAU : end
        CALL ASSYM2(IRREP,VRT(1,ISPIN),NT2DIS,ICORE(I000))
        CALL SSCAL(NT2DIS*NT2DSZ,HALF,ICORE(I000),1)
        CALL TRANSP(ICORE(I000),ICORE(I040),NT2DSZ,NT2DIS)
c YAU : old
c       CALL ICOPY(NT2DIS*NT2DSZ*IINTFP,ICORE(I040),1,ICORE(I000),1)
c YAU : new
        CALL DCOPY(NT2DIS*NT2DSZ,ICORE(I040),1,ICORE(I000),1)
c YAU : end
C
C NOW DO THE F(MI) CONTRIBUTION
C
C
C FIRST EXPAND T2 AND Z VECTORS TO Z(A<B,IJ), T(A<B,IJ).
C
        CALL SYMEXP(IRREP,POP(1,ISPIN),NT2DSZ,ICORE(I000))
        CALL SYMEXP(IRREP,POP(1,ISPIN),NT2DSZ,ICORE(I010))
C
C PERFORM MULTIPLICATION
C                                 
C                    Z(A<B,IJ) = Z(A<B,IJ) + T2(A<BI,M) * F(MJ)
C
        IOFFT  =I000
        IOFFTAR=I010
        IOFFF  =I030
        DO 120 IRREPM=1,NIRREP
         IRREPI=DIRPRD(IRREPM,IRREP)
         NUMM=POP(IRREPM,ISPIN)
         NUMI=POP(IRREPI,ISPIN)
         NUMJ=NUMM
         IF(MAX(NT2DSZ*NUMI,NUMM,NUMJ).NE.0)THEN
          CALL XGEMM('N','N',NT2DSZ*NUMI,NUMJ,NUMM,ONEM,ICORE(IOFFT),
     &               NT2DSZ*NUMI,ICORE(IOFFF),NUMM,ONE,
     &               ICORE(IOFFTAR),NT2DSZ*NUMI)
         ENDIF
         IOFFF=IOFFF+NUMM*NUMJ*IINTFP
         IOFFT=IOFFT+NT2DSZ*NUMI*NUMM*IINTFP
         IOFFTAR=IOFFTAR+NT2DSZ*NUMI*NUMM*IINTFP
120     CONTINUE
C
C NOW ANTISYMMETRIZE - Z(A<B,I<J) = Z(A<B,IJ) - Z(A<B,JI), ADD IN
C  LEADING TERM AND DENOMINATOR WEIGHT
C
        CALL ASSYM2(IRREP,POP(1,ISPIN),NT2DSZ,ICORE(I000))
        CALL SSCAL(NT2DIS*NT2DSZ,HALF,ICORE(I000),1)
        CALL ASSYM2(IRREP,POP(1,ISPIN),NT2DSZ,ICORE(I010))
        CALL GETLST(ICORE(I040),1,NT2DIS,1,IRREP,13+ISPIN)
        CALL SAXPY(NT2DSZ*NT2DIS,ONE,ICORE(I040),1,ICORE(I010),1)
        CALL GETLST(ICORE(I040),1,NT2DIS,1,IRREP,47+ISPIN)
        CALL VECPRD(ICORE(I010),ICORE(I040),ICORE(I010),NT2DIS*NT2DSZ)
C
C COMPARE WITH PREVIOUS ITERATE
C
        CALL SAXPY(NT2DSZ*NT2DIS,ONEM,ICORE(I010),1,ICORE(I000),1)
        X=FNDLRGAB(ICORE(I000),NT2DSZ*NT2DIS)
        IF(ABS(X).Gt.1.E-08)THEN
         ITER=ITER+1
c YAU : old
c        CALL ICOPY(NT2DIS*NT2DSZ*IINTFP,ICORE(I010),1,ICORE(I000),1)
c YAU : new
         CALL DCOPY(NT2DIS*NT2DSZ,ICORE(I010),1,ICORE(I000),1)
c YAU : end
         GOTO 1
        ELSE
         ITER=ITER+1
         WRITE(6,1001)IRREP,ITER
         CALL PUTLST(ICORE(I010),1,NT2DIS,1,IRREP,43+ISPIN)
         CALL GETLST(ICORE(I000),1,NT2DIS,1,IRREP,13+ISPIN)
         Z=Z+SDOT(NT2DIS*NT2DSZ,ICORE(I000),1,ICORE(I010),1)
        ENDIF
1001    FORMAT(T3,'Irrep ',I2,' amplitudes converged in ',I4,
     &         ' iterations.')
10     CONTINUE
       WRITE(6,1002)SPCASE(ISPIN),Z
1002   FORMAT(T3,'Spin case ',A2,' T2-W energy contribution ',
     &        F15.10)
       ET2=ET2+Z
5     CONTINUE
C
C
C SPIN CASE AB
C
C     Z(Ab,Ij) = T(Ab,Ij) + T(Ae,Ij)*F(be) - T(Eb,Ij)*F(EA)
C
C                         - T(Ab,Im)*F(mj) + T(Ab,Mj)*F(MI)
C
      ISPIN=3
      WRITE(6,1000)SPCASE(ISPIN)
      Z=0.0
      DO 110 IRREP=1,NIRREP
       ITER=0
       NT2DSZ=IRPDPD(IRREP,ISYTYP(1,46))
       NT2DIS=IRPDPD(IRREP,ISYTYP(2,46))
       IF(MIN(NT2DSZ,NT2DIS).EQ.0)GOTO 110
       I000=1
       I010=I000+NT2DSZ*NT2DIS*IINTFP
       I020=I010+NT2DSZ*NT2DIS*IINTFP
       I030=I020+NT2DSZ*NT2DIS*IINTFP
C
C READ IN INITIAL GUESS T2 AMPLITUDES (NO F CONTRIBUTIONS) AND
C  TRANSPOSE FOR FIRST CONTRACTION
C
       CALL GETLST(ICORE(I000),1,NT2DIS,1,IRREP,46)
C
C DO FEA CONTRACTION FIRST
C
101    CONTINUE
       CALL TRANSP(ICORE(I000),ICORE(I010),NT2DIS,NT2DSZ)
c YAU : old
c      CALL ICOPY(IINTFP*NT2DIS*NT2DSZ,ICORE(I010),1,ICORE(I000),1)
c YAU : new
       CALL DCOPY(NT2DIS*NT2DSZ,ICORE(I010),1,ICORE(I000),1)
c YAU : end
       CALL IZERO(ICORE(I010),NT2DIS*NT2DSZ*IINTFP)
       DO 115 ISPIN=2,1,-1
C
C READ IN VV AND OO FOCK MATRICES
C
        I040=I030+NFEA(ISPIN)*IINTFP
        I050=I040+NFMI(ISPIN)*IINTFP
        CALL GETLST(ICORE(I030),1,1,1,2+ISPIN,92)
        CALL GETLST(ICORE(I040),1,1,1,2+ISPIN,91)
C
C TRANSPOSE KET INDICES IF ISPIN=1
C
C    T(Ij,Ab) => [for ISPIN=1] T(Ij,bA)
C
        IF(ISPIN.EQ.1)THEN
         I060=I050+NT2DIS*IINTFP
         I070=I060+NT2DIS*IINTFP
         I080=I070+NT2DIS*IINTFP
         CALL SYMTR1(IRREP,VRT(1,1),VRT(1,2),NT2DIS,ICORE(I010),
     &               ICORE(I050),ICORE(I060),ICORE(I070))
         CALL SYMTR1(IRREP,VRT(1,1),VRT(1,2),NT2DIS,ICORE(I000),
     &               ICORE(I050),ICORE(I060),ICORE(I070))
        ENDIF
C
C PERFORM MULTIPLICATION
C                                 
C                    Z(Ij,bA) = T2(Ijb,E) * F(EA) [ISPIN=1]
C
C                    Z(Ij,Ba) = T2(Ij,Be) * F(ea) [ISPIN=2]
C
        IOFFT  =I000
        IOFFTAR=I010
        IOFFF  =I030
        DO 220 IRREPE=1,NIRREP
         IRREPB=DIRPRD(IRREPE,IRREP)
         NUMB=VRT(IRREPB,3-ISPIN)
         NUME=VRT(IRREPE,ISPIN)
         NUMA=NUME
         IF(MAX(NT2DIS*NUMB,NUME,NUMA).NE.0)THEN
          CALL XGEMM('N','N',NT2DIS*NUMB,NUMA,NUME,ONE,ICORE(IOFFT),
     &               NT2DIS*NUMB,ICORE(IOFFF),NUME,ONE,
     &               ICORE(IOFFTAR),NT2DIS*NUMB)
         ENDIF
         IOFFF=IOFFF+NUME*NUMA*IINTFP
         IOFFT=IOFFT+NT2DIS*NUMB*NUME*IINTFP
         IOFFTAR=IOFFTAR+NT2DIS*NUMB*NUME*IINTFP
220     CONTINUE
115    CONTINUE
C
C ISPIN=1 RAN LAST SO WE HAVE TO TRANSPOSE KET INDICES BACK TO
C    T(Ij,Ab) AND Z(Ij,Ab).  THEN TRANSPOSE TO T(Ab,Ij) AND
C    Z(Ab,Ij) FOR FMI PART.
C
       I060=I050+NT2DIS*IINTFP
       I070=I060+NT2DIS*IINTFP
       I080=I070+NT2DIS*IINTFP
       CALL SYMTR1(IRREP,VRT(1,2),VRT(1,1),NT2DIS,ICORE(I000),
     &             ICORE(I050),ICORE(I060),ICORE(I070))
       CALL SYMTR1(IRREP,VRT(1,2),VRT(1,1),NT2DIS,ICORE(I010),
     &              ICORE(I050),ICORE(I060),ICORE(I070))
       CALL TRANSP(ICORE(I000),ICORE(I020),NT2DSZ,NT2DIS)
c YAU : old
c      CALL ICOPY(NT2DIS*NT2DSZ*IINTFP,ICORE(I020),1,ICORE(I000),1)
c YAU : new
       CALL DCOPY(NT2DIS*NT2DSZ,ICORE(I020),1,ICORE(I000),1)
c YAU : end
       CALL TRANSP(ICORE(I010),ICORE(I020),NT2DSZ,NT2DIS)
c YAU : old
c      CALL ICOPY(NT2DIS*NT2DSZ*IINTFP,ICORE(I020),1,ICORE(I010),1)
c YAU : new
       CALL DCOPY(NT2DIS*NT2DSZ,ICORE(I020),1,ICORE(I010),1)
c YAU : end
C
C NOW DO SECOND CONTRACTION - FMI PART
C
C         Z(Ab,Ij) = Z(Ab,Ij) + T(Ab,Mj) * F(MI)  [ISPIN=1]
C
C         Z(Ab,Ij) = Z(Ab,Ij) + T(Ab,Im) * F(mj)  [ISPIN=2]
C
C
C
C FOR ISPIN=1, WE HAVE TO DO A SYMTR1 ON THE T AMPLITUDES AND
C  INCREMENTS, GIVING T(Ab,jI) AND Z(Ab,jI).
C
       DO 116 ISPIN=2,1,-1
C
C READ IN VV AND OO FOCK MATRICES
C
        I040=I030+NFEA(ISPIN)*IINTFP
        I050=I040+NFMI(ISPIN)*IINTFP
        CALL GETLST(ICORE(I030),1,1,1,2+ISPIN,92)
        CALL GETLST(ICORE(I040),1,1,1,2+ISPIN,91)
        IF(ISPIN.EQ.1)THEN
         I060=I050+NT2DSZ*IINTFP
         I070=I060+NT2DSZ*IINTFP
         I080=I070+NT2DSZ*IINTFP
         CALL SYMTR1(IRREP,POP(1,1),POP(1,2),NT2DSZ,ICORE(I000),
     &               ICORE(I050),ICORE(I060),ICORE(I070))
         CALL SYMTR1(IRREP,POP(1,1),POP(1,2),NT2DSZ,ICORE(I010),
     &               ICORE(I050),ICORE(I060),ICORE(I070))
        ENDIF
C
C NOW PERFORM THE MATRIX MULTIPLICATION
C
C         Z(Ab,jI) = Z(Ab,jI) + T(Ab,jM) * F(MI)  [ISPIN=1]
C
C         Z(Ab,Ij) = Z(Ab,Ij) + T(Ab,Im) * F(mj)  [ISPIN=2]
C
        IOFFT=I000
        IOFFF=I040
        IOFFTAR=I010
        DO 320 IRREPM=1,NIRREP
         IRREPJ=DIRPRD(IRREPM,IRREP)
         IRREPI=IRREPM
         NUMM=POP(IRREPM,ISPIN)
         NUMI=POP(IRREPI,ISPIN)
         NUMJ=POP(IRREPJ,3-ISPIN)
         IF(MAX(NT2DIS*NUMJ,NUMM,NUMI).NE.0)THEN
          CALL XGEMM('N','N',NT2DSZ*NUMJ,NUMI,NUMM,ONEM,ICORE(IOFFT),
     &               NT2DSZ*NUMJ,ICORE(IOFFF),NUMM,ONE,ICORE(IOFFTAR),
     &               NT2DSZ*NUMJ)
         ENDIF
         IOFFT=IOFFT+NT2DSZ*NUMJ*NUMM*IINTFP
         IOFFF=IOFFF+NUMM*NUMI*IINTFP
         IOFFTAR=IOFFTAR+NT2DSZ*NUMJ*NUMM*IINTFP
320     CONTINUE
C
C NOW TRANSPOSE KET INDICES OF TARGET AND T2 VECTOR IF ISPIN IS 1
C
C        Z(Ab,jI) -> Z(Ab,Ij)
C
        IF(ISPIN.EQ.1)THEN
         I060=I050+NT2DSZ*IINTFP
         I070=I060+NT2DSZ*IINTFP
         I080=I070+NT2DSZ*IINTFP
         CALL SYMTR1(IRREP,POP(1,2),POP(1,1),NT2DSZ,ICORE(I000),
     &               ICORE(I050),ICORE(I060),ICORE(I070))
         CALL SYMTR1(IRREP,POP(1,2),POP(1,1),NT2DSZ,ICORE(I010),
     &               ICORE(I050),ICORE(I060),ICORE(I070))
        ENDIF
116    CONTINUE
C
C NOW ADD IN LEADING TERM AND DENOMINATOR WEIGHT
C
       CALL GETLST(ICORE(I020),1,NT2DIS,1,IRREP,16)
       CALL SAXPY(NT2DSZ*NT2DIS,ONE,ICORE(I020),1,ICORE(I010),1)
       CALL GETLST(ICORE(I020),1,NT2DIS,1,IRREP,50)
       CALL VECPRD(ICORE(I010),ICORE(I020),ICORE(I010),NT2DIS*NT2DSZ)
C
C COMPARE WITH PREVIOUS ITERATE
C
       CALL SAXPY(NT2DSZ*NT2DIS,ONEM,ICORE(I010),1,ICORE(I000),1)
       X=FNDLRGAB(ICORE(I000),NT2DSZ*NT2DIS)
       IF(ABS(X).Gt.1.E-12)THEN
        ITER=ITER+1
c YAU : old
c       CALL ICOPY(NT2DIS*NT2DSZ*IINTFP,ICORE(I010),1,ICORE(I000),1)
c YAU : new
        CALL DCOPY(NT2DIS*NT2DSZ,ICORE(I010),1,ICORE(I000),1)
c YAU : end
        GOTO 101
       ELSE
        ITER=ITER+1
        CALL PUTLST(ICORE(I010),1,NT2DIS,1,IRREP,46)
        CALL GETLST(ICORE(I000),1,NT2DIS,1,IRREP,16)
        Z=Z+SDOT(NT2DIS*NT2DSZ,ICORE(I000),1,ICORE(I010),1)
        WRITE(6,1001)IRREP,ITER
       ENDIF
110   CONTINUE
      WRITE(6,1002)SPCASE(3),Z
      ET2=ET2+Z
      RETURN
      END
