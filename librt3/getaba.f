      SUBROUTINE GETABA(T2ABA,IRPK,K,IRREPX,LIST)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION T2ABA(1)
      DIMENSION IADT2(8),LENT2(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
C
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
      INDEX(I) = I*(I-1)/2
C
      ISPIN1 = 1
      ISPIN2 = 2
C
C     SET ADDRESSES FOR T2.
C
      DO   30 IRPL=1,NIRREP
      IRPKL = DIRPRD(IRPL,IRPK)
      IRPBC = DIRPRD(IRPKL,IRREPX)
      LENT2(IRPL) = IRPDPD(IRPBC,13) * POP(IRPL,ISPIN1)
      IF(IRPL.EQ.1)THEN
      IADT2(IRPL) = 1
      ELSE
      IADT2(IRPL) = IADT2(IRPL-1) + LENT2(IRPL-1)
      ENDIF
   30 CONTINUE
C
C     READ T2(BC,L).
C             AB A
C
      DO   70 IRPL=1,NIRREP
      IF(POP(IRPL,ISPIN1).EQ.0) GOTO 70
      IRPKL = DIRPRD(IRPK,IRPL)
C
      KL = IOFFOO(IRPK,IRPKL,5) + (K-1)*POP(IRPL,ISPIN1) + 1
      CALL GETLST(T2ABA(IADT2(IRPL)),
     1            KL,POP(IRPL,ISPIN1),1,IRPKL,LIST)
   70 CONTINUE
C
      RETURN
      END
