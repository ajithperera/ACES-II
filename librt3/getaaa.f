      SUBROUTINE GETAAA(SCR,T2AAA,IRPJ,J,ISPIN,IUHF,IRREPX,LISTU,LISTR)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Arguments :
C-----------------------------------------------------------------------
      DOUBLE PRECISION SCR,T2AAA
      INTEGER IRPJ,J,ISPIN,IUHF,IRREPX,LISTU,LISTR
C-----------------------------------------------------------------------
C     Common blocks :
C-----------------------------------------------------------------------
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD,
     &        NSTART,NIRREP,IRREPA,IRREPB,DIRPRD,
     &        IRPDPD,ISYTYP,ID,
     &        POP,VRT,NTAA,NTBB,NF1AA,NF2AA,NF1BB,NF2BB,
     &        IOFFVV,IOFFOO,IOFFVO
C-----------------------------------------------------------------------
C     Local variables :
C-----------------------------------------------------------------------
      INTEGER IRPL,IRPJL,IRPAB,LENT2,IADT2,L,JL,JLAB,LJAB,LENAB,I000,
     &        I010
C-----------------------------------------------------------------------
C     In-line function :
C-----------------------------------------------------------------------
      INTEGER INDEX,I
C-----------------------------------------------------------------------

      DIMENSION SCR(1),T2AAA(1)
      DIMENSION IADT2(8),LENT2(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     &                NF1BB,NF2BB
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
      INDEX(I) = I*(I-1)/2
C
C     SET ADDRESSES FOR T2.
C
      DO  130 IRPL=1,NIRREP
      IRPJL = DIRPRD(IRPL,IRPJ)
      IRPAB = DIRPRD(IRREPX,IRPJL)
      LENT2(IRPL) = IRPDPD(IRPAB,ISPIN) * POP(IRPL,ISPIN)
      IF(IRPL.EQ.1)THEN
      IADT2(IRPL) = 1
      ELSE
      IADT2(IRPL) = IADT2(IRPL-1) + LENT2(IRPL-1)
      ENDIF
  130 CONTINUE
C
      IF(IUHF.EQ.0) GOTO 200
C
C     READ T2(A<B,L).
C
      DO  170 IRPL=1,NIRREP
      IF(POP(IRPL,ISPIN).EQ.0) GOTO 170
      IRPJL = DIRPRD(IRPJ,IRPL)
      IRPAB = DIRPRD(IRREPX,IRPJL)
C
      IF(IRPJ.LT.IRPL)THEN
      DO  140    L=1,POP(IRPL,ISPIN)
      JL = IOFFOO(IRPL,IRPJL,ISPIN) + (L-1)*POP(IRPJ,ISPIN) + J
      CALL GETLST(T2AAA(IADT2(IRPL) + (L-1)*IRPDPD(IRPAB,ISPIN)),
     &            JL,1,1,IRPJL,LISTU)
  140 CONTINUE
      ENDIF
C
      IF(IRPJ.GT.IRPL)THEN
      JL = IOFFOO(IRPJ,IRPJL,ISPIN) + (J-1)*POP(IRPL,ISPIN) + 1
      CALL GETLST(T2AAA(IADT2(IRPL)),
     &            JL,POP(IRPL,ISPIN),1,IRPJL,LISTU)
      CALL VMINUS(T2AAA(IADT2(IRPL)),LENT2(IRPL))
      ENDIF
C
      IF(IRPJ.EQ.IRPL)THEN
      DO  160    L=1,POP(IRPL,ISPIN)
C
      IF(J.LT.L)THEN
      JL = IOFFOO(IRPL,IRPJL,ISPIN) + INDEX(L-1) + J
      CALL GETLST(T2AAA(IADT2(IRPL) + (L-1)*IRPDPD(IRPAB,ISPIN)),
     &            JL,1,1,IRPJL,LISTU)
      ENDIF
C
      IF(J.GT.L)THEN
      JL = IOFFOO(IRPJ,IRPJL,ISPIN) + INDEX(J-1) + L
      CALL GETLST(T2AAA(IADT2(IRPL) + (L-1)*IRPDPD(IRPAB,ISPIN)),
     &            JL,1,1,IRPJL,LISTU)
      CALL VMINUS(T2AAA(IADT2(IRPL) + (L-1)*IRPDPD(IRPAB,ISPIN)),
     &            IRPDPD(IRPAB,ISPIN))
      ENDIF
C
      IF(J.EQ.L)THEN
      CALL ZERO(T2AAA(IADT2(IRPL) + (L-1)*IRPDPD(IRPAB,ISPIN)),
     &          IRPDPD(IRPAB,ISPIN))
      ENDIF
  160 CONTINUE
C
      ENDIF
  170 CONTINUE
      RETURN
C
  200 CONTINUE
C
C     RHF case.
C
      DO  270 IRPL=1,NIRREP
      IF(POP(IRPL,ISPIN).EQ.0) GOTO 270
      IRPJL = DIRPRD(IRPJ,IRPL)
      IRPAB = DIRPRD(IRREPX,IRPJL)
      LENAB = IRPDPD(IRPAB,13)
      I000  = 1
      I010  = I000 + LENAB
C
      DO  260    L=1,POP(IRPL,ISPIN)
C
      JLAB = IOFFOO(IRPL,IRPJL,5) + (L-1)*POP(IRPJ,ISPIN) + J
      LJAB = IOFFOO(IRPJ,IRPJL,5) + (J-1)*POP(IRPL,ISPIN) + L
      CALL GETLST(SCR(I000),JLAB,1,1,IRPJL,LISTR)
      CALL GETLST(SCR(I010),LJAB,1,1,IRPJL,LISTR)
      CALL   VADD(SCR(I000),SCR(I000),SCR(I010),LENAB,-1.0D+00)
      CALL  SQSYM(IRPAB,VRT(1,1),IRPDPD(IRPAB,1),LENAB,1,
     &            T2AAA(IADT2(IRPL) + (L-1)*IRPDPD(IRPAB,1)),SCR(I000))
C
      IF(IRPJ.EQ.IRPL.AND.J.EQ.L)THEN
      CALL ZERO(T2AAA(IADT2(IRPL) + (L-1)*IRPDPD(IRPAB,ISPIN)),
     &          IRPDPD(IRPAB,ISPIN))
      ENDIF
  260 CONTINUE
  270 CONTINUE
      RETURN
      END
