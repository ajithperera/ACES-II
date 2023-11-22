      SUBROUTINE S1S214N(S1A,S1B,CORE,W,IADW,ISPIN,I,J,K,
     1                  IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,IMODE)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION S1A(1),S1B(1),CORE(1),W(1)
      DIMENSION IADW(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /INFO/   NOCCO(2),NVRTO(2)
C
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
      INDEX(I) = I*(I-1)/2
C
C     Generic subroutine for forming, in combination with SYMCONTW,
C     disconnected triple amplitudes for AAA or BBB cases. This routine
C     computes
C
C      ABC       A  BC          A  BC           A  BC
C     S      =  S  S      -    S  S       +    S  S
C      IJK       I  JK          J  IK           K  IJ
C     
C     For example,
C
C                   A      A         BC
C                  S   =  T    ;    S    =  <BC//JK>
C                   I      I         JK
C
C     IMODE = 1 : S1 is T1       , S2 is integral.
C     IMODE = 2 : S1 is F        , S2 is T2
C     IMODE = 3 : S1 is passed in, S2 is T2
C     IMODE = 4 : S1 is T1       , S2 is T2
C
      IF(IMODE.EQ.1)THEN
       LISTS2 = 13 + ISPIN
       LISTS1 = 90
       PARTS1 =      ISPIN
      ENDIF
      IF(IMODE.EQ.2)THEN
       LISTS2 = 43 + ISPIN
       LISTS1 = 93
       PARTS1 =  2 + ISPIN
      ENDIF
      IF(IMODE.EQ.3)THEN
       LISTS2 = 43 + ISPIN
      ENDIF
      IF(IMODE.EQ.4)THEN
       LISTS2 = 43 + ISPIN
       LISTS1 = 90
       PARTS1 =      ISPIN
      ENDIF
      IF(IMODE.NE.1.AND.IMODE.NE.2.AND.IMODE.NE.3.AND.IMODE.NE.4)THEN
       WRITE(6,*) ' @S1S214-I, Invalid value of IMODE. '
       CALL ERREX
      ENDIF
C
      IF(IRPI.EQ.IRPJ)THEN
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN) + INDEX(J-1) + I
      ELSE
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN) + (J-1)*POP(IRPI,ISPIN) + I
      ENDIF
      IF(IRPI.EQ.IRPK)THEN
      IK = IOFFOO(IRPK,IRPIK,ISPIN) + INDEX(K-1) + I
      ELSE
      IK = IOFFOO(IRPK,IRPIK,ISPIN) + (K-1)*POP(IRPI,ISPIN) + I
      ENDIF
      IF(IRPJ.EQ.IRPK)THEN
      JK = IOFFOO(IRPK,IRPJK,ISPIN) + INDEX(K-1) + J
      ELSE
      JK = IOFFOO(IRPK,IRPJK,ISPIN) + (K-1)*POP(IRPJ,ISPIN) + J
      ENDIF
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPIJ,ISPIN)
      I020 = I010 + IRPDPD(IRPIK,ISPIN)
      I030 = I020 + IRPDPD(IRPJK,ISPIN)
      CALL GETLST(CORE(I000),IJ,1,1,IRPIJ,LISTS2)
      CALL GETLST(CORE(I010),IK,1,1,IRPIK,LISTS2)
      CALL GETLST(CORE(I020),JK,1,1,IRPJK,LISTS2)
C
      IF(ISPIN.EQ.1)THEN
      NS1 = NTAA
      ELSE
      NS1 = NTBB
      ENDIF
      I040 = I030 + NS1
      IF(IMODE.EQ.1.OR.IMODE.EQ.2.OR.IMODE.EQ.4)THEN
      CALL GETLST(CORE(I030),1,1,1,PARTS1,LISTS1)
      ENDIF
      IF(IMODE.EQ.3)THEN
c YAU : old
c      IF(ISPIN.EQ.1) CALL ICOPY(IINTFP*NS1,S1A,1,CORE(I030),1)
c      IF(ISPIN.EQ.2) CALL ICOPY(IINTFP*NS1,S1B,1,CORE(I030),1)
c YAU : new
       IF(ISPIN.EQ.1) CALL DCOPY(NS1,S1A,1,CORE(I030),1)
       IF(ISPIN.EQ.2) CALL DCOPY(NS1,S1B,1,CORE(I030),1)
c YAU : end
      ENDIF
C
      IRPABC = IRPIJK
C
      DO  100 IRPA=1,NIRREP
      IF(VRT(IRPA,ISPIN).EQ.0) GOTO 100
C
      IRPBC = DIRPRD(IRPA,IRPABC)
      IF(IRPDPD(IRPBC,ISPIN).EQ.0) GOTO 100
C
      IF(IRPA.EQ.IRPI)THEN
C
      DO   30  A=1,VRT(IRPA,ISPIN)
      DO   20 BC=1,IRPDPD(IRPBC,ISPIN)
C
       W(IADW(IRPA) + (A-1)*IRPDPD(IRPBC,ISPIN) + BC - 1) =
     1 W(IADW(IRPA) + (A-1)*IRPDPD(IRPBC,ISPIN) + BC - 1) +
     1 CORE(I020 + BC - 1) * 
     1 CORE(I030 - 1 + IOFFVO(IRPI,1,ISPIN)  + 
     1                 (I-1)*VRT(IRPA,ISPIN) + A)
   20 CONTINUE
   30 CONTINUE
C
      ENDIF
C
      IF(IRPA.EQ.IRPJ)THEN
C
      DO   50  A=1,VRT(IRPA,ISPIN)
      DO   40 BC=1,IRPDPD(IRPBC,ISPIN)
C
       W(IADW(IRPA) + (A-1)*IRPDPD(IRPBC,ISPIN) + BC - 1) =
     1 W(IADW(IRPA) + (A-1)*IRPDPD(IRPBC,ISPIN) + BC - 1) -
     1 CORE(I010 + BC - 1) *
     1 CORE(I030 - 1 + IOFFVO(IRPJ,1,ISPIN)  +
     1                 (J-1)*VRT(IRPA,ISPIN) + A)
   40 CONTINUE
   50 CONTINUE
C
      ENDIF
C
      IF(IRPA.EQ.IRPK)THEN
C
      DO   70  A=1,VRT(IRPA,ISPIN)
      DO   60 BC=1,IRPDPD(IRPBC,ISPIN)
C
       W(IADW(IRPA) + (A-1)*IRPDPD(IRPBC,ISPIN) + BC - 1) =
     1 W(IADW(IRPA) + (A-1)*IRPDPD(IRPBC,ISPIN) + BC - 1) +
     1 CORE(I000 + BC - 1) *
     1 CORE(I030 - 1 + IOFFVO(IRPK,1,ISPIN)  +
     1                 (K-1)*VRT(IRPA,ISPIN) + A)
   60 CONTINUE
   70 CONTINUE
C
      ENDIF
C
      IF(IRPA.NE.IRPI.AND.IRPA.NE.IRPJ.AND.IRPA.NE.IRPK)THEN
C
      DO   90  A=1,VRT(IRPA,ISPIN)
      DO   80 BC=1,IRPDPD(IRPBC,ISPIN)
C
       W(IADW(IRPA) + (A-1)*IRPDPD(IRPBC,ISPIN) + BC - 1) = 0.0D+00

   80 CONTINUE
   90 CONTINUE
C
      ENDIF
C
  100 CONTINUE
      RETURN
      END
