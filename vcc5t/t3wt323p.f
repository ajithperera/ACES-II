      SUBROUTINE T3WT323P(D3T3,D3T3EXP,CORE,MAXCOR,IUHF,ISPIN1,ISPIN2,
     1                    LEN,LENEXP,IRPIJK,IJKVAL,IADT3,IADW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER VRT,POP,DIRPRD,DISSIZ
      INTEGER E,F,EF,FE,ELTF,ELOW,EHIGH,FLOW,FHIGH
      DIMENSION D3T3(LEN),D3T3EXP(LENEXP),CORE(1)
      DIMENSION IADT3(8),IADW(8),IADW2(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NT1(2),NFMI(2),NFEA(2)
      COMMON /FLAGS/  IFLAGS(100)
      EQUIVALENCE(ICLLVL,IFLAGS( 2))
      EQUIVALENCE(IDRLVL,IFLAGS( 3))
      COMMON /FILES/  LUOUT,MOINTS
C
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
      COMMON /T3IOOF/ IJKPOS(8,8,8,2),IJKLEN(36,8,4),IJKOFF(36,8,4),
     1                NCOMB(4)
C
      INDEX(I) = I*(I-1)/2
C
      I000  = 1
      I010  = I000 + LEN
      I020  = I010 + LENEXP
      I030  = I020 + LENEXP
      NEED  = I030 * IINTFP
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
 1010 FORMAT(' @T3WT323P-F, Insufficient memory. Need ',I10,' Got ',I10)
      CALL INSMEM('T3WT323P',NEED,MAXCOR)
      ENDIF
C
      CALL GETLIST(CORE(I000),IJKVAL,1,1,IRPIJK,ISPIN1 + 1 + 4)
c     write(6,*) ' in t3wt323p. 1st '
c      call sumblk(core(i000),len)
C
C     W(AB,EF) * T(EFc,IJk)
C     W(ab,ef) * T(efC,ijK)
C
      IOFF1 = I000
      IOFF2 = 1
      DO   10 IRPC=1,NIRREP
      IRPEF = DIRPRD(IRPC,IRPIJK)
      IRPAB = IRPEF
C
      IF(VRT(IRPC,ISPIN2).EQ.0.OR.IRPDPD(IRPAB,ISPIN1).EQ.0) GOTO 10
C
      IF(IUHF.EQ.0)THEN
      I040 = I030 + IRPDPD(IRPAB,ISPIN1) * IRPDPD(IRPEF,ISPIN1)
      I050 = I040 + IRPDPD(IRPAB,    13)
      I060 = I050 + IRPDPD(IRPAB,    13)
      NEED = I060 * IINTFP
      ELSE
      I040 = I030 + IRPDPD(IRPAB,ISPIN1) * IRPDPD(IRPEF,ISPIN1)
      I050 = I040 + IRPDPD(IRPAB,    13)
      I060 = I050 + IRPDPD(IRPAB,    13)
      NEED = I040 * IINTFP
      ENDIF
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
      CALL INSMEM('T3WT323P',NEED,MAXCOR)
      ENDIF
C
      IF(IUHF.EQ.0)THEN
      DO   5 IRPF=1,NIRREP
      IRPE = DIRPRD(IRPEF,IRPF)
      IF(IRPE.GT.IRPF) GOTO 5
C
      IF(IRPE.EQ.IRPF)THEN
      FHIGH = VRT(IRPF,1)
      FLOW  = 2
      ELSE
      FHIGH = VRT(IRPF,1)
      FLOW  = 1
      ENDIF
C
      IF(FLOW.GT.FHIGH) GOTO 5
      DO   3 F=FLOW,FHIGH
C
      IF(IRPE.EQ.IRPF)THEN
      EHIGH = F-1
      ELOW  = 1
      ELSE
      EHIGH = VRT(IRPE,1)
      ELOW  = 1
      ENDIF
C
      IF(ELOW.GT.EHIGH) GOTO 3
      DO   2 E=ELOW,EHIGH
C
      EF = IOFFVV(IRPF,IRPEF,5) + (F-1)*VRT(IRPE,1) + E
      FE = IOFFVV(IRPE,IRPEF,5) + (E-1)*VRT(IRPF,1) + F
C
      IF(IRPEF.EQ.1)THEN
      ELTF = IOFFVV(IRPF,IRPEF,1) + INDEX(F-1)        + E
      ELSE
      ELTF = IOFFVV(IRPF,IRPEF,1) + (F-1)*VRT(IRPE,1) + E
      ENDIF
C
      CALL GETLST(CORE(I040),EF,1,2,IRPEF,233)
      CALL GETLST(CORE(I050),FE,1,2,IRPEF,233)
      CALL   VADD(CORE(I040),CORE(I040),CORE(I050),IRPDPD(IRPAB,13),
     1            -1.0D+00)
      CALL  SQSYM(IRPAB,VRT(1,1),IRPDPD(IRPAB, 1),IRPDPD(IRPAB,13),
     1            1,CORE(I050),CORE(I040))
c YAU : old
c     CALL ICOPY(IINTFP*IRPDPD(IRPAB,1),
c    1           CORE(I050),1,CORE(I030+(ELTF-1)*IRPDPD(IRPAB,1)),1)
c YAU : new
      CALL DCOPY(IRPDPD(IRPAB,1),
     1           CORE(I050),1,CORE(I030+(ELTF-1)*IRPDPD(IRPAB,1)),1)
c YAU : end
    2 CONTINUE
    3 CONTINUE
    5 CONTINUE
      ELSE
      CALL GETLST(CORE(I030),1,IRPDPD(IRPEF,ISPIN1),2,IRPEF,230+ISPIN1)
      ENDIF
C
      DISSIZ = IRPDPD(IRPAB,ISPIN1)
      NSUM   = IRPDPD(IRPEF,ISPIN1)
      NDIS   = VRT(IRPC,ISPIN2)
c     write(6,*) ' before xgemm : d3t3, W, t3 '
c      call sumblk(d3t3(ioff2),dissiz*ndis)
c      call sumblk(core(i030),dissiz*nsum)
c      call sumblk(core(ioff1),nsum*ndis)
      CALL XGEMM('N','N',DISSIZ,NDIS,NSUM,1.0D+00,
     1           CORE(I030),DISSIZ,CORE(IOFF1),NSUM,1.0D+00,
     1           D3T3(IOFF2),DISSIZ)
c     write(6,*) ' after xgemm '
c      call sumblk(d3t3(ioff2),dissiz*ndis)
C
      IOFF1 = IOFF1 + NSUM   * NDIS
      IOFF2 = IOFF2 + DISSIZ * NDIS
   10 CONTINUE
C
C    T(Ef,A) * W(Bc,Ef) - T(Ef,B) * W(Ac,Ef)
C    T(eF,a) * W(Cb,Fe) - T(eF,b) * W(Ca,Fe)
C  [ T(eF,a) * W(bC,eF) - T(eF,b) * W(aC,eF) ]
C
C     Expand T3.
C
      CALL ZERO(CORE(I010),LENEXP)
      CALL ZERO(CORE(I020),LENEXP)
      CALL SYMTRW2(CORE(I000),CORE(I010),CORE(I020),IADT3,IADW2,
     1             IRPIJK,ISPIN1,ISPIN2)
C
C     At CORE(I010) we have T3(E,f,A)/T3(e,F,a), labelled by IRPA. Use
C     transpose so we form D3T3EXP(A,B,c)/D3T3EXP(a,b,C).
C
      IOFF1 = I010
      DO 20 IRPA=1,NIRREP
C
      IRPEF = DIRPRD(IRPA,IRPIJK)
      IRPBC = IRPEF
C
      DISSIZ = VRT(IRPA,ISPIN1)
      NDIS   = IRPDPD(IRPBC,13)
      NSUM   = IRPDPD(IRPEF,13)
C
      IF(DISSIZ.EQ.0.OR.NDIS.EQ.0.OR.NSUM.EQ.0) GOTO 20
C
      I040 = I030 + NDIS * NSUM
      I050 = I040 + MAX(NDIS,NSUM)
      I060 = I050 + MAX(NDIS,NSUM)
      I070 = I060 + MAX(NDIS,NSUM)
      NEED = I070 * IINTFP
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
      CALL INSMEM('T3WT323P',NEED,MAXCOR)
      ENDIF
C
      CALL GETLST(CORE(I030),1,NSUM,2,IRPEF,233)
C
C     If this is case 3 triples, we must convert the intermediate from
C     ABAB to BABA (see TRPS3).
C
      CALL XGEMM('T','T',DISSIZ,NDIS,NSUM,1.0D+00,
     1           CORE(IOFF1),NSUM,CORE(I030),NDIS,1.0D+00,
     1           D3T3EXP(IADW(IRPBC)),DISSIZ)
C
      IOFF1 = IOFF1 + NSUM * DISSIZ
   20 CONTINUE
      RETURN
      END
