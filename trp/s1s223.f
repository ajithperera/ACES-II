      SUBROUTINE S1S223(CORE,T3D,W,IADT3,IADW,ISPIN1,ISPIN2,I,J,K,
     1                  IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,IUHF,
     1                  IMODE,LCCSD_T)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION CORE(1),T3D(1),W(1)
      DIMENSION IADT3(8),IADW(8)
      LOGICAL LCCSD_T
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
C     Generic subroutine for forming, in combination with EXPSC2/3,
C     disconnected triple amplitudes for AAB/BBA spin cases. This routine
C     computes
C
C      ABC       A  BC          B  AC   A  AB      B BA
C     S      =  S  S      -    S  S
C      IJK       I  JK          I  JK   A  AB      B BA
C
C                A  BC          B  AC   A  AB      B BA
C            -  S  S      +    S  S
C                J  IK          J  IK   A  AB      B BA
C
C                C  AB                  B  AA      A BB
C            +  S  S
C                K  IJ                  B  AA      A BB
C
C     
C     For example,
C
C                   A      A         BC
C                  S   =  T    ;    S    =  <BC//JK>
C                   I      I         JK
C
C     AAB triples :
C
C     ISPIN1 = 1
C     ISPIN2 = 2
C
C     BBA triples :
C
C     ISPIN1 = 2
C     ISPIN2 = 1
C
C     IMODE = 1 : S1 is T1, S2 is integral
C     IMODE = 2 : S1 is F , S2 is T2
C
       IF (LCCSD_T) THEN
         LISTOFF=100
       ELSE
         LISTOFF=0
       END IF

      IF(IMODE.EQ.1)THEN
       LISTS2A = 13 + ISPIN1
       LISTS2B = 16
       LISTS1=90+LISTOFF
      
       PARTS1A =      ISPIN2
       PARTS1B =      ISPIN1
      ELSE
      IF(IMODE.EQ.2)THEN
          LISTS1 = 93
       LISTS2A=43+ISPIN1+LISTOFF
       LISTS2B=46+LISTOFF
          
       PARTS1A =  2 + ISPIN2
       PARTS1B =  2 + ISPIN1
      ELSE
      WRITE(6,*) ' @S1S223-I, Invalid value of IMODE. '
      CALL ERREX
      ENDIF
      ENDIF
C
      IF(ISPIN1.EQ.1)THEN
      NS1A = NTBB
      NS1B = NTAA
      ELSE
      NS1A = NTAA
      NS1B = NTBB
      ENDIF
C
      IF(IRPI.EQ.IRPJ)THEN
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN1) + INDEX(J-1) + I
      ELSE
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN1) + (J-1)*POP(IRPI,ISPIN1) + I
      ENDIF
C
C     To make work for case 3 as well, we will probably have to eliminate
C     the number 5 here.
C
      IF(ISPIN1.EQ.1)THEN
      IK = IOFFOO(IRPK,IRPIK,5) + (K-1)*POP(IRPI,ISPIN1) + I
      JK = IOFFOO(IRPK,IRPJK,5) + (K-1)*POP(IRPJ,ISPIN1) + J
      ELSE
      IK = IOFFOO(IRPI,IRPIK,5) + (I-1)*POP(IRPK,ISPIN2) + K
      JK = IOFFOO(IRPJ,IRPJK,5) + (J-1)*POP(IRPK,ISPIN2) + K
      ENDIF
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPIJ,ISPIN1)
      I020 = I010 + IRPDPD(IRPIK,13)
      I030 = I020 + IRPDPD(IRPJK,13)
      CALL GETLST(CORE(I000),IJ,1,1,IRPIJ,LISTS2A)
      CALL GETLST(CORE(I010),IK,1,1,IRPIK,LISTS2B)
      CALL GETLST(CORE(I020),JK,1,1,IRPJK,LISTS2B)
C
C      IF(ISPIN1.EQ.2.AND.IMODE.EQ.2)THEN
C      CALL SYMTRA2(IRPIK,VRT(1,2),VRT(1,1),IRPDPD(IRPIK,13),1,
C     1             CORE(I010),CORE(I030))
C      CALL ICOPY(IRPDPD(IRPIK,13)*IINTFP,CORE(I030),1,CORE(I010),1)
C      CALL SYMTRA2(IRPJK,VRT(1,2),VRT(1,1),IRPDPD(IRPJK,13),1,
C     1             CORE(I020),CORE(I030))
C      CALL ICOPY(IRPDPD(IRPJK,13)*IINTFP,CORE(I030),1,CORE(I020),1)
C      ENDIF
C
      I040 = I030 + NS1B
      CALL GETLST(CORE(I030),1,1,1,PARTS1B,LISTS1)
C
      IRPABC = IRPIJK
C
      DO  100 IRPBC=1,NIRREP
      IF(IRPDPD(IRPBC,13).EQ.0) GOTO 100
C
      IRPA = DIRPRD(IRPBC,IRPABC)
      IF(VRT(IRPA,ISPIN1).EQ.0) GOTO 100
C
      IF(IRPA.EQ.IRPI)THEN
C
      DO   30 BC=1,IRPDPD(IRPBC,13)
      DO   20  A=1,VRT(IRPA,ISPIN1)
C
       W(IADW(IRPBC) + (BC-1)*VRT(IRPA,ISPIN1) + A - 1) =
     1 W(IADW(IRPBC) + (BC-1)*VRT(IRPA,ISPIN1) + A - 1) +
     1 CORE(I020 + BC - 1) * 
     1 CORE(I030 - 1 + IOFFVO(IRPI,1,ISPIN1)  +
     1                 (I-1)*VRT(IRPA,ISPIN1) + A)
   20 CONTINUE
   30 CONTINUE
C
      ENDIF
C
      IF(IRPA.EQ.IRPJ)THEN
C
      DO   50 BC=1,IRPDPD(IRPBC,13)
      DO   40  A=1,VRT(IRPA,ISPIN1)
C
       W(IADW(IRPBC) + (BC-1)*VRT(IRPA,ISPIN1) + A - 1) =
     1 W(IADW(IRPBC) + (BC-1)*VRT(IRPA,ISPIN1) + A - 1) -
     1 CORE(I010 + BC - 1) *
     1 CORE(I030 - 1 + IOFFVO(IRPJ,1,ISPIN1)  +
     1                 (J-1)*VRT(IRPA,ISPIN1) + A)
   40 CONTINUE
   50 CONTINUE
C
      ENDIF
C
      IF(IRPA.NE.IRPI.AND.IRPA.NE.IRPJ)THEN
C
      DO   90 BC=1,IRPDPD(IRPBC,13)
      DO   80  A=1,VRT(IRPA,ISPIN1)
C
       W(IADW(IRPBC) + (BC-1)*VRT(IRPA,ISPIN1) + A - 1) = 0.0D+00
C
   80 CONTINUE
   90 CONTINUE
C
      ENDIF
C
  100 CONTINUE
C
      I040 = I030 + NS1A
      IF(IUHF.EQ.0)THEN
      CALL GETLST(CORE(I030),1,1,1,PARTS1B,LISTS1)
      ELSE
      CALL GETLST(CORE(I030),1,1,1,PARTS1A,LISTS1)
      ENDIF
C
      DO  200 IRPC=1,NIRREP
      IF(VRT(IRPC,ISPIN2).EQ.0) GOTO 200
C
      IRPAB = DIRPRD(IRPC,IRPABC)
      IF(IRPDPD(IRPAB,ISPIN1).EQ.0) GOTO 200
C
      IF(IRPC.EQ.IRPK)THEN
C
      DO  130  C=1,VRT(IRPC,ISPIN2)
      DO  120 AB=1,IRPDPD(IRPAB,ISPIN1)
C
       T3D(IADT3(IRPC) + (C-1)*IRPDPD(IRPAB,ISPIN1) + AB - 1) =
     1 T3D(IADT3(IRPC) + (C-1)*IRPDPD(IRPAB,ISPIN1) + AB - 1) +
     1 CORE(I000 + AB - 1) *
     1 CORE(I030 - 1 + IOFFVO(IRPK,1,ISPIN2)  +
     1                 (K-1)*VRT(IRPC,ISPIN2) + C)
  120 CONTINUE
  130 CONTINUE
C
      ELSE
C
      CALL ZERO(T3D(IADT3(IRPC)),IRPDPD(IRPAB,ISPIN1)*VRT(IRPC,ISPIN2))
C
      ENDIF
C
  200 CONTINUE
C
      IF(ISPIN1.EQ.1)THEN
      CALL EXPSC2(T3D,W,IADT3,IADW,IRPIJK)
      ELSE
      CALL EXPSC3(T3D,W,IADT3,IADW,IRPIJK)
      ENDIF
      RETURN
      END
