
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
C     IMODE = 1 : S1 is T1, S2 is integral.
C     IMODE = 2 : S1 is F , S2 is T2

      SUBROUTINE S1S214(CORE,W,IADW,ISPIN,I,J,K,
     1                  IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,
     &                  IMODE, LMODE)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION CORE(1),W(1)
      LOGICAL LMODE
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

C ----------------------------------------------------------------------

      IF (IMODE.EQ.1) THEN
         LISTS2 = 13 + ISPIN
         IF (LMODE) THEN
            LISTS1 = 190
         ELSE
            LISTS1 = 90
         END IF
         PARTS1 =      ISPIN
      ELSE
         IF (IMODE.EQ.2) THEN
            IF (LMODE) THEN
               LISTS2 = 143 + ISPIN
            ELSE
               LISTS2 = 43 + ISPIN
               LISTS1 = 93
            END IF 
            LISTS1 = 93 
            PARTS1 =  2 + ISPIN
         ELSE
            WRITE(*,*) '@S1S214: Invalid value of IMODE.'
            CALL ERREX
         END IF
      END IF
C
      IF (IRPI.EQ.IRPJ) THEN
         IJ = IOFFOO(IRPJ,IRPIJ,ISPIN) + INDEX(J-1) + I
      ELSE
         IJ = IOFFOO(IRPJ,IRPIJ,ISPIN) + (J-1)*POP(IRPI,ISPIN) + I
      END IF
      IF (IRPI.EQ.IRPK) THEN
         IK = IOFFOO(IRPK,IRPIK,ISPIN) + INDEX(K-1) + I
      ELSE
         IK = IOFFOO(IRPK,IRPIK,ISPIN) + (K-1)*POP(IRPI,ISPIN) + I
      END IF
      IF (IRPJ.EQ.IRPK) THEN
         JK = IOFFOO(IRPK,IRPJK,ISPIN) + INDEX(K-1) + J
      ELSE
         JK = IOFFOO(IRPK,IRPJK,ISPIN) + (K-1)*POP(IRPJ,ISPIN) + J
      END IF
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPIJ,ISPIN)
      I020 = I010 + IRPDPD(IRPIK,ISPIN)
      I030 = I020 + IRPDPD(IRPJK,ISPIN)
      CALL GETLST(CORE(I000),IJ,1,1,IRPIJ,LISTS2)
      CALL GETLST(CORE(I010),IK,1,1,IRPIK,LISTS2)
      CALL GETLST(CORE(I020),JK,1,1,IRPJK,LISTS2)
C
      IF (ISPIN.EQ.1) THEN
         NS1 = NTAA
      ELSE
         NS1 = NTBB
      END IF
      I040 = I030 + NS1
      CALL GETLST(CORE(I030),1,1,1,PARTS1,LISTS1)
C
      IRPABC = IRPIJK
C
      DO IRPA  = 1, NIRREP
         IRPBC = DIRPRD(IRPA,IRPABC)
         MAX_A  = VRT(IRPA,ISPIN)
         MAX_BC = IRPDPD(IRPBC,ISPIN)
      IF ((MAX_A.NE.0).AND.(MAX_BC.NE.0)) THEN
C
      IF (IRPA.EQ.IRPI) THEN
         DO A = 1, MAX_A
            DO BC = 1, MAX_BC
         W(IADW(IRPA) + (A-1)*MAX_BC + BC - 1) =
     1   W(IADW(IRPA) + (A-1)*MAX_BC + BC - 1)
     1 + (   CORE(I020 + BC - 1)
     1     * CORE(I030 - 1 + IOFFVO(IRPI,1,ISPIN) + (I-1)*MAX_A + A) )
            END DO
         END DO
      END IF
C
      IF (IRPA.EQ.IRPJ) THEN
         DO A = 1, MAX_A
            DO BC = 1, MAX_BC
         W(IADW(IRPA) + (A-1)*MAX_BC + BC - 1) =
     1   W(IADW(IRPA) + (A-1)*MAX_BC + BC - 1)
     1 - (   CORE(I010 + BC - 1)
     1     * CORE(I030 - 1 + IOFFVO(IRPJ,1,ISPIN) + (J-1)*MAX_A + A) )
            END DO
         END DO
      END IF
C
      IF (IRPA.EQ.IRPK) THEN
         DO A = 1, MAX_A
            DO BC = 1, MAX_BC
         W(IADW(IRPA) + (A-1)*MAX_BC + BC - 1) =
     1   W(IADW(IRPA) + (A-1)*MAX_BC + BC - 1)
     1 + (   CORE(I000 + BC - 1)
     1     * CORE(I030 - 1 + IOFFVO(IRPK,1,ISPIN) + (K-1)*MAX_A + A) )
            END DO
         END DO
      END IF
C
      IF ((IRPA.NE.IRPI).AND.(IRPA.NE.IRPJ).AND.(IRPA.NE.IRPK)) THEN
         DO A = 1, MAX_A
            DO BC = 1, MAX_BC
               W(IADW(IRPA) + (A-1)*MAX_BC + BC - 1) = 0.0D0
            END DO
         END DO
      END IF
C
C     END IF ((MAX_A.NE.0).AND.(MAX_BC.NE.0))
      END IF
C     END DO IRPA = 1, NIRREP
      END DO

      RETURN
      END

