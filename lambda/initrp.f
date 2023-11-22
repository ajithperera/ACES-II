      SUBROUTINE INITRP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /FLAGS/  IFLAGS(100)
      EQUIVALENCE(ICLLVL,IFLAGS( 2))
      EQUIVALENCE(IDRLVL,IFLAGS( 3))
      EQUIVALENCE(IREFNC,IFLAGS(11))
      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
     1                LWIC15,LWIC16,LWIC17,LWIC18,
     1                LWIC21,LWIC22,LWIC23,LWIC24,
     1                LWIC25,LWIC26,LWIC27,LWIC28,
     1                LWIC31,LWIC32,LWIC33,
     1                LWIC34,LWIC35,LWIC36,
     1                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
C
      WRITE(6,1000)
 1000 FORMAT(' @INITRP-I, Initializing W intermediate lists. ')
C
C            --- basic D3T3 = WT2 contraction ---
C
C     --- MBPT(4),CCSD+T(CCSD),CCSD(T),QCISD(T),CCSDT-1A,CCSDT-1B,UCC(4) ---
C
      IF(ICLLVL.EQ. 4.OR.ICLLVL.EQ.11.OR.ICLLVL.EQ.13.OR.ICLLVL.EQ.9.OR.
     1   ICLLVL.EQ.14.OR.ICLLVL.EQ.21.OR.ICLLVL.EQ.22)THEN
C
      LWIC11 =   7 + 300 - 300
      LWIC12 =   8 + 300 - 300
      LWIC13 =   9 + 300 - 300
      LWIC14 =  10 + 300 - 300
C
      LWIC15 =  27 + 300 - 300
      LWIC16 =  28 + 300 - 300
      LWIC17 =  29 + 300 - 300
      LWIC18 =  30 + 300 - 300
C
      ENDIF
C
C     --- CCSDT-2, CCSDT-3, CCSDT-4, CCSDT ---
C
      IF(ICLLVL.EQ.15.OR.ICLLVL.EQ.16.OR.ICLLVL.EQ.17.OR.
     1   ICLLVL.EQ.18)THEN
C
      LWIC11 =   7 + 100
      LWIC12 =   8 + 100
      LWIC13 =   9 + 100
      LWIC14 =  10 + 100
C
      LWIC15 =  27 + 100
      LWIC16 =  28 + 100
      LWIC17 =  29 + 100
      LWIC18 =  30 + 100
C
      ENDIF
C
C              --- T3 inclusion in T2 equation ---
C
C     --- MBPT(4),CCSD+T(CCSD),CCSD(T),QCISD(T) gradients;CCSDT-1A,UCC(4) ---
C
      IF(ICLLVL.EQ. 4.OR.ICLLVL.EQ.11.OR.ICLLVL.EQ.13.OR.ICLLVL.EQ.9.OR.
     1   ICLLVL.EQ.21.OR.ICLLVL.EQ.22)THEN
C
      LWIC21 =   7 + 300 - 300
      LWIC22 =   8 + 300 - 300
      LWIC23 =   9 + 300 - 300
      LWIC24 =  10 + 300 - 300
C
      LWIC25 =  27 + 300 - 300
      LWIC26 =  28 + 300 - 300
      LWIC27 =  29 + 300 - 300
      LWIC28 =  30 + 300 - 300
C
      ENDIF
C
C     --- CCSDT-1B, CCSDT-2, CCSDT-3, CCSDT-4, CCSDT ---
C
      IF(ICLLVL.EQ.14.OR.ICLLVL.EQ.15.OR.ICLLVL.EQ.16.OR.
     1   ICLLVL.EQ.17.OR.ICLLVL.EQ.18)THEN
C
      LWIC21 =   7 + 300
      LWIC22 =   8 + 300
      LWIC23 =   9 + 300
      LWIC24 =  10 + 300
C
      LWIC25 =  27 + 300
      LWIC26 =  28 + 300
      LWIC27 =  29 + 300
      LWIC28 =  30 + 300
C
      ENDIF
C
C         --- D3T3 = WT3 contraction ---
C
C     --- CCSDT-4, CCSDT (wishful thinking) ---
C
      IF(ICLLVL.EQ.17.OR.ICLLVL.EQ.18)THEN
      ENDIF
      RETURN
      END
