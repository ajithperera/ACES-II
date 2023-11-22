      SUBROUTINE T2FINL1(L1,ICORE,MAXCOR,IUHF,ISPIN,BREDUNDANT)
C
C THIS ROUTINE COMPUTES THE CONTRIBUTION OF THE Fme INTERMEDIATE
C  IN THE T1 EQUATION.
C
C         L[2](i,a) = - T[1](im,ae) * f(me) 
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,L1
      LOGICAL BREDUNDANT
C
      DIMENSION ICORE(MAXCOR),L1(1)
C
      COMMON/INFO/NOCA,NOCB,NVRTA,NVRTB
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C
      DATA ONE,ONEM,ZILCH /1.0D0,-1.0D0,0.0D0/
C
C L[2](A;I) = - f(E,M) * T[1](E,M;A,I) - f(e,m) * T[1](e,m;A,I)
C
C L[2](a;i) = - f(e,m) * T[1](e,m;a,i) - f(E,M) * T[1](E,M;a,i)
C

      LISTTA=33+ISPIN
      LISTTB=35+ISPIN
      IF(IUHF.EQ.0)LISTTB=37
      LISTF =93
C
C WE ONLY NEED TO CONSIDER ONE DPD IRREP -- THE TOTALLY SYMMETRIC ONE.
C
      NLENT1=IRPDPD(1,8+ISPIN)
      NLENFA=IRPDPD(1,8+ISPIN)
      NLENFB=IRPDPD(1,11-ISPIN)
C
C I010 HOLDS FIRST F IN EQUATION ABOVE.
C I020 HOLDS SECOND F IN EQUATION ABOVE.
C I030 HOLDS INDIVIDUAL A1 IRREP OF BLOCKED T2 LIST.
C
      I010=1
      IF(IUHF.EQ.0)THEN
       I020=I010
      ELSE
       I020=I010+NLENFA*IINTFP
      ENDIF
      I030=I020+NLENFB*IINTFP
      CALL GETLST(ICORE(I010),1,1,1,2+ISPIN,LISTF)
      IF(IUHF.NE.0)CALL GETLST(ICORE(I020),1,1,1,5-ISPIN,LISTF)
      NTDSZA=IRPDPD(1,ISYTYP(1,LISTTA))
      NTDSZB=IRPDPD(1,ISYTYP(1,LISTTB))
      NTDISA=IRPDPD(1,ISYTYP(2,LISTTA))
      NTDISB=IRPDPD(1,ISYTYP(2,LISTTB))
      I040=I030+NTDSZA*NTDISA*IINTFP
      IF(I040.GT.MAXCOR)CALL INSMEM('T2FINL1',I040,MAXCOR)
      IF(BREDUNDANT) THEN
         CALL GETLST(ICORE(I030),1,NTDISA,1,1,LISTTA)
      ELSE
         CALL GETLST_NR(ICORE(I030),ICORE(I040),MAXCOR-I040,
     &                  LISTTA, 1)
      ENDIF
c      call chksums("t2finl1 listta 1: ",icore(i030),ntdisa*ntdsza)
      CALL XGEMM('N','N',1,NLENT1,NLENFA,ONE,ICORE(I010),
     &            1,ICORE(I030),NTDSZA,ZILCH,L1,
     &            1)
      I040=I030+NTDSZB*NTDISB*IINTFP
      IF(I040.GT.MAXCOR)CALL INSMEM('T2FINL1',I040,MAXCOR)
      IF(BREDUNDANT) THEN
         CALL GETLST(ICORE(I030),1,NTDISB,1,1,LISTTB)
      ELSE
         CALL GETLST_NR(ICORE(I030),ICORE(I040),MAXCOR-I040,
     &                  LISTTB, 1)
      ENDIF
c      call chksums("t2finl1 listtb 1: ",icore(i030),ntdisb*ntdszb)
      CALL XGEMM('N','N',1,NLENT1,NLENFB,ONEM,ICORE(I020),
     &            1,ICORE(I030),NTDSZB,ONE,L1,
     &            1)
20    CONTINUE
      RETURN
      END