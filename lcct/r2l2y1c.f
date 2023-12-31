      SUBROUTINE R2L2Y1C(Y1,F1,ICORE,MAXCOR,IUHF)
C
C Y1(ai) = -1/2 F(if)*R(mn,ef)*L(mn,ea)
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE
      DIMENSION ICORE(MAXCOR),Y1(*),F1(*)
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA ONE/1.0D0/
C
      IOFFF=1
      IOFFY=1
      DO 10 ISPIN=1,1+IUHF
       I0G=1
       I000=I0G+IINTFP*NFEA(1)
C
C FIRST PICK UP G(fa) FROM LIST 192
C
       CALL GETLST(ICORE(I0G),1,1,1,ISPIN,192)
C
C FORM PRODUCT Y1(ai) = G(fa) * F(fi)
C
       IOFFG=I0G
       DO 11 IRREPI=1,NIRREP
        IRREPA=IRREPI
        IRREPF=IRREPA
        NUMA=VRT(IRREPA,ISPIN)
        NUMF=VRT(IRREPF,ISPIN)
        NUMI=POP(IRREPI,ISPIN)
        CALL XGEMM('T','N',NUMA,NUMI,NUMF,ONE,ICORE(IOFFG),NUMF,
     &             F1(IOFFF),NUMF,ONE,Y1(IOFFY),NUMA)
        IOFFF=IOFFF+IINTFP*NUMF*NUMI
        IOFFG=IOFFG+IINTFP*NUMF*NUMA
        IOFFY=IOFFY+IINTFP*NUMA*NUMI
11     CONTINUE
C
10    CONTINUE
C
      RETURN
      END
