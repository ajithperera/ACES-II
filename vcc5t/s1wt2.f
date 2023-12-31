      SUBROUTINE S1WT2(S1A,S1B,CORE,MAXCOR,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION SDOT
      INTEGER A,E,AM,EF
      INTEGER DISSIZ,POP,VRT,DIRPRD
      DIMENSION S1A(1),S1B(1),CORE(1)
      DIMENSION IOFFOV(8,8,4),IOFFT1(8,2)
C
      COMMON /FILES/  LUOUT,MOINTS
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
C     Computation of T1 like quantities, defined by
C
C     S1(A,I) =   Sum   <EF||AM> T(EF,IM) + Sum   <MN||IE> T(AE,MN)
C                 E<F                       M<N
C                  M                         E
C
      CALL  ZERO(S1A,NTAA)
      CALL  ZERO(S1B,NTBB)
      CALL IZERO(IOFFOV,256)
      CALL MKOFOV(IOFFOV)
      CALL MKOFOO
      CALL MKOFVV
      CALL MKOFVO
C
      CALL IZERO(IOFFT1,16)
      DO  10 ISPIN=1,2
      DO   5 IRREP=1,NIRREP
      IF(IRREP.GT.1)THEN
      IOFFT1(IRREP,ISPIN) = IOFFT1(IRREP-1,ISPIN) + 
     1                         VRT(IRREP-1,ISPIN) * POP(IRREP-1,ISPIN)
      ENDIF
    5 CONTINUE
   10 CONTINUE
C
      DO   100 ISPIN=1,IUHF+1
C
      LISTV = 26 + ISPIN
      LISTT = 43 + ISPIN
C
      DO    90 IRPEF=1,NIRREP
      IF(IRPDPD(IRPEF,  ISPIN).EQ.0.OR.IRPDPD(IRPEF,8+ISPIN).EQ.0.OR.
     1   IRPDPD(IRPEF,2+ISPIN).EQ.0)                          GOTO 90
C
      DISSIZ = IRPDPD(IRPEF,   ISPIN)
      NDISV  = IRPDPD(IRPEF, 8+ISPIN)
      NDIST  = IRPDPD(IRPEF, 2+ISPIN)
      NDISTE = IRPDPD(IRPEF,20+ISPIN)
C
      I000 = 1
      I010 = I000 + DISSIZ
      I020 = I010 + DISSIZ * NDISTE
      I030 = I020 + MAX(IRPDPD(IRPEF,18+ISPIN),IRPDPD(IRPEF,13))
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
 1010 FORMAT(' @S1WT2-F, Insufficient memory. Need ',I10,' Got ',I10)
      CALL ERREX
      ENDIF
C
      CALL GETLST(CORE(I010),1,NDIST,1,IRPEF,LISTT)
C
C     Expand the i<M part of T to i,m in place.
C
      CALL SYMEXP(IRPEF,POP(1,ISPIN),DISSIZ,CORE(I010))
C
      DO   50 IRPM=1,NIRREP
      IRPI = DIRPRD(IRPM,IRPEF)
      IRPA = IRPI
C
      IF(VRT(IRPA,ISPIN).EQ.0.OR.POP(IRPI,ISPIN).EQ.0.OR.
     1                           POP(IRPM,ISPIN).EQ.0) GOTO 50
C
      DO   40  M=1,POP(IRPM,ISPIN)
      DO   30  A=1,VRT(IRPA,ISPIN)
C
      IF(IUHF.EQ.0)THEN
      LISTV = 30
      AM = IOFFVO(IRPM,IRPEF,ISPIN) + (M-1)*VRT(IRPA,ISPIN) + A
      CALL GETLST(CORE(I020),AM,1,2,IRPEF,LISTV)
      CALL ASSYM3(IRPEF,VRT(1,1),IRPDPD(IRPEF,1),IRPDPD(IRPEF,13),
     1            1,CORE(I000),CORE(I020),IOFFVV,1)
      ELSE
      AM = IOFFVO(IRPM,IRPEF,ISPIN) + (M-1)*VRT(IRPA,ISPIN) + A
      CALL GETLST(CORE(I000),AM,1,2,IRPEF,LISTV)
      ENDIF
C
      DO   20  I=1,POP(IRPI,ISPIN)
C
      IM = IOFFOO(IRPM,IRPEF,2+ISPIN) + (M-1)*POP(IRPI,ISPIN) + I
C
      IF(ISPIN.EQ.1)THEN
       S1A(IOFFT1(IRPI,ISPIN) + (I-1)*VRT(IRPA,ISPIN) + A) =
     1 S1A(IOFFT1(IRPI,ISPIN) + (I-1)*VRT(IRPA,ISPIN) + A) +
     1 SDOT(DISSIZ,CORE(I000)                  ,1,
     1             CORE(I010 + (IM-1) * DISSIZ),1)
      ELSE
       S1B(IOFFT1(IRPI,ISPIN) + (I-1)*VRT(IRPA,ISPIN) + A) =
     1 S1B(IOFFT1(IRPI,ISPIN) + (I-1)*VRT(IRPA,ISPIN) + A) +
     1 SDOT(DISSIZ,CORE(I000)                  ,1,
     1             CORE(I010 + (IM-1) * DISSIZ),1)
      ENDIF
C
   20 CONTINUE
   30 CONTINUE
   40 CONTINUE
   50 CONTINUE
   90 CONTINUE
  100 CONTINUE
C
C
C     IPASS = 1  :  <Ef||Am> * T(Ef,Im)
C     IPASS = 2  :  <Fe||Ma> * T(Fe,Mi)
C
      DO   200 IPASS=1,IUHF+1
C
      IF(IPASS.EQ.1)THEN
      ISPIN1 = 1
      ISPIN2 = 2
      ELSE
      ISPIN1 = 2
      ISPIN2 = 1
      ENDIF
C
      LISTV = 31 - ISPIN1
      LISTT = 46
C
      DO   190 IRPEF=1,NIRREP
      IF(IRPDPD(IRPEF,13).EQ.0.OR.IRPDPD(IRPEF,10+IPASS).EQ.0.OR.
     1   IRPDPD(IRPEF,14).EQ.0)                          GOTO 190
C
      DISSIZ = IRPDPD(IRPEF,13)
      NDIST  = IRPDPD(IRPEF,14)
C
      I000 = 1
      I010 = I000 + DISSIZ
      I020 = I010 + DISSIZ * NDIST
C
      NEED = IINTFP * I020
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
      CALL ERREX
      ENDIF
C
      CALL GETLST(CORE(I010),1,NDIST,1,IRPEF,LISTT)
C
      DO  150 IRPM=1,NIRREP
      IRPI = DIRPRD(IRPM,IRPEF)
      IRPA = IRPI
C
      IF(VRT(IRPA,ISPIN1).EQ.0.OR.POP(IRPI,ISPIN1).EQ.0.OR.
     1                            POP(IRPM,ISPIN2).EQ.0) GOTO 150
C
      DO  140  M=1,POP(IRPM,ISPIN2)
      DO  130  A=1,VRT(IRPA,ISPIN1)
C
      IF(IPASS.EQ.1)THEN
      AM = IOFFVO(IRPM,IRPEF,4) + (M-1)*VRT(IRPA,ISPIN1) + A
      CALL GETLST(CORE(I000),AM,1,2,IRPEF,LISTV)
      ELSE
      MA = IOFFOV(IRPA,IRPEF,3) + (A-1)*POP(IRPM,ISPIN2) + M
      CALL GETLST(CORE(I000),MA,1,2,IRPEF,LISTV)
      ENDIF
C
      DO  120  I=1,POP(IRPI,ISPIN1)
C
      IF(IPASS.EQ.1)THEN
      IM = IOFFOO(IRPM,IRPEF,5) + (M-1)*POP(IRPI,ISPIN1) + I
C
       S1A(IOFFT1(IRPI,IPASS) + (I-1)*VRT(IRPA,ISPIN1) + A) =
     1 S1A(IOFFT1(IRPI,IPASS) + (I-1)*VRT(IRPA,ISPIN1) + A) +
     1 SDOT(DISSIZ,CORE(I000)                  ,1,
     1             CORE(I010 + (IM-1) * DISSIZ),1)
      ELSE
      MI = IOFFOO(IRPI,IRPEF,5) + (I-1)*POP(IRPM,ISPIN2) + M
       S1B(IOFFT1(IRPI,IPASS) + (I-1)*VRT(IRPA,ISPIN1) + A) =
     1 S1B(IOFFT1(IRPI,IPASS) + (I-1)*VRT(IRPA,ISPIN1) + A) +
     1 SDOT(DISSIZ,CORE(I000)                  ,1,
     1             CORE(I010 + (MI-1) * DISSIZ),1)
      ENDIF
C
  120 CONTINUE
  130 CONTINUE
  140 CONTINUE
  150 CONTINUE
  190 CONTINUE
  200 CONTINUE
C
C
      DO   300 ISPIN=1,IUHF+1
C
      LISTV =  6 + ISPIN
      LISTT = 43 + ISPIN
C
      DO   290 IRPMN=1,NIRREP
      IF(IRPDPD(IRPMN,  ISPIN).EQ.0.OR.IRPDPD(IRPMN,8+ISPIN).EQ.0.OR.
     1   IRPDPD(IRPMN,2+ISPIN).EQ.0)                         GOTO 290
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPMN, 2+ISPIN) * IRPDPD(IRPMN, 8+ISPIN)
      I020 = I010 +    MAX(
     1              IRPDPD(IRPMN, 2+ISPIN) , IRPDPD(IRPMN,18+ISPIN))
C
      NEED = IINTFP * I020
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
      CALL ERREX
      ENDIF
C
C     Put V(IE,MN) at CORE(I000).
C
      CALL GETTRN(CORE(I000),CORE(I010),IRPDPD(IRPMN, 2+ISPIN),
     1                                  IRPDPD(IRPMN, 8+ISPIN),
     1            2,IRPMN,LISTV)
C
      DO   280 MN=1,IRPDPD(IRPMN, 2+ISPIN)
C
      CALL GETLST(CORE(I010),MN,1,1,IRPMN,LISTT)
C
C     Expand A<E to A,E.
C
      CALL SYMEXP2(IRPMN,VRT(1,ISPIN),IRPDPD(IRPMN,18+ISPIN),
     1                                IRPDPD(IRPMN,   ISPIN),1,
     1             CORE(I010),CORE(I010))
C
      DO   270 IRPE=1,NIRREP
      IRPI = DIRPRD(IRPE,IRPMN)
      IRPA = IRPI
C
      IF(VRT(IRPA,ISPIN).EQ.0.OR.POP(IRPI,ISPIN).EQ.0.OR.
     1                           VRT(IRPE,ISPIN).EQ.0) GOTO 270
C
      IOFFV = I000 + (MN-1)*IRPDPD(IRPMN, 8+ISPIN) + 
     1               IOFFOV(IRPE,IRPMN,  ISPIN)
      IOFFT = I010 + IOFFVV(IRPE,IRPMN,2+ISPIN)
C
      IF(ISPIN.EQ.1)THEN
      CALL XGEMM('N','T',
     1           VRT(IRPA,ISPIN),POP(IRPI,ISPIN),VRT(IRPE,ISPIN),
     1           -1.0D+00,
     1           CORE(IOFFT),VRT(IRPA,ISPIN),
     1           CORE(IOFFV),POP(IRPI,ISPIN),1.0D+00,
     1            S1A(IOFFT1(IRPI,ISPIN)+1),VRT(IRPA,ISPIN))
      ELSE
      CALL XGEMM('N','T',
     1           VRT(IRPA,ISPIN),POP(IRPI,ISPIN),VRT(IRPE,ISPIN),
     1           -1.0D+00,
     1           CORE(IOFFT),VRT(IRPA,ISPIN),
     1           CORE(IOFFV),POP(IRPI,ISPIN),1.0D+00,
     1            S1B(IOFFT1(IRPI,ISPIN)+1),VRT(IRPA,ISPIN))
      ENDIF
  270 CONTINUE
  280 CONTINUE
  290 CONTINUE
  300 CONTINUE
C
C     IPASS = 1  :  T(Ae,Mn) * <Ie||Mn>
C     IPASS = 2  :  T(Ea,Nm) * <Ei||Nm>
C
      DO   400 IPASS=1,IUHF+1
C
      IF(IPASS.EQ.1)THEN
      ISPIN1 = 1
      ISPIN2 = 2
      ELSE
      ISPIN1 = 2
      ISPIN2 = 1
      ENDIF
C
      LISTV = 11 - ISPIN1
      LISTT = 46
C
      DO   390 IRPMN=1,NIRREP
      IF(IRPDPD(IRPMN,13).EQ.0.OR.IRPDPD(IRPMN,13-ISPIN1).EQ.0.OR.
     1   IRPDPD(IRPMN,14).EQ.0)                           GOTO 390
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPMN,14) * IRPDPD(IRPMN,13-ISPIN1)
      I020 = I010 +    MAX(
     1              IRPDPD(IRPMN,14) , IRPDPD(IRPMN,13))
C
      NEED = IINTFP * I020
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
      CALL ERREX
      ENDIF
C
C     Put V(Ie,Mn)/V(Ei,Nm) at CORE(I000).
C
      CALL GETTRN(CORE(I000),CORE(I010),IRPDPD(IRPMN,14),
     1                                  IRPDPD(IRPMN,13-ISPIN1),
     1            2,IRPMN,LISTV)
C
      DO   380 MN=1,IRPDPD(IRPMN,14)
C
      CALL GETLST(CORE(I010),MN,1,1,IRPMN,LISTT)
C
      DO   370 IRPE=1,NIRREP
      IRPI = DIRPRD(IRPE,IRPMN)
      IRPA = IRPI
C
      IF(VRT(IRPA,ISPIN1).EQ.0.OR.POP(IRPI,ISPIN1).EQ.0.OR.
     1                            VRT(IRPE,ISPIN2).EQ.0) GOTO 370
      IF(IPASS.EQ.1)THEN
      IOFFT = I010 + IOFFVV(IRPE,IRPMN,5)
      IOFFV = I000 + (MN-1)*IRPDPD(IRPMN,13-ISPIN1) +
     1        IOFFOV(IRPE,IRPMN,3)
      CALL XGEMM('N','T',
     1           VRT(IRPA,ISPIN1),POP(IRPI,ISPIN1),VRT(IRPE,ISPIN2),
     1           -1.0D+00,
     1           CORE(IOFFT),VRT(IRPA,ISPIN1),
     1           CORE(IOFFV),POP(IRPI,ISPIN1),1.0D+00,
     1            S1A(IOFFT1(IRPI,IPASS)+1),VRT(IRPA,ISPIN1))
      ELSE
      IOFFT = I010 + IOFFVV(IRPA,IRPMN,5)
      IOFFV = I000 + (MN-1) * IRPDPD(IRPMN,13-ISPIN1) + 
     1        IOFFVO(IRPI,IRPMN,4)
      CALL XGEMM('T','N',
     1           VRT(IRPA,ISPIN1),POP(IRPI,ISPIN1),VRT(IRPE,ISPIN2),
     1           -1.0D+00,
     1           CORE(IOFFT),VRT(IRPE,ISPIN2),
     1           CORE(IOFFV),VRT(IRPE,ISPIN2),1.0D+00,
     1            S1B(IOFFT1(IRPI,IPASS)+1),VRT(IRPA,ISPIN1))
      ENDIF
  370 CONTINUE
  380 CONTINUE
  390 CONTINUE
C
  400 CONTINUE
C
 1000 CONTINUE
C
c YAU : old
c     IF(IUHF.EQ.0) CALL ICOPY(IINTFP*NTAA,S1A,1,S1B,1)
c YAU : new
      IF(IUHF.EQ.0) CALL DCOPY(NTAA,S1A,1,S1B,1)
c YAU : end
C
c      DO 500 K=1,ntaa
c      write(6,*) s1a(k),s1b(k)
c  500 Continue
c      CALL SUMBLK(S1A,NTAA)
c      CALL SUMBLK(S1B,NTBB)
      RETURN
      END
