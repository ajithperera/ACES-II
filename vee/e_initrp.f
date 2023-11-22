      SUBROUTINE E_INITRP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /FLAGS/  IFLAGS(100)
      EQUIVALENCE(ICLLVL,IFLAGS( 2))
      EQUIVALENCE(IDRLVL,IFLAGS( 3))
      EQUIVALENCE(IREFNC,IFLAGS(11))
      COMMON /FLAGS2/ IFLAGS2(500)
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
      IF(ICLLVL.NE.10.AND.ICLLVL.NE.13.AND.ICLLVL.NE.14.AND.
     &   ICLLVL.NE.16.AND.ICLLVL.NE.5.AND.ICLLVL.NE.6.AND.
     &   ICLLVL.NE.22.AND.ICLLVL.NE.33.AND.ICLLVL.NE.34)THEN
        WRITE(6,*) ' @INITRP-F, EOM not coded for method ',ICLLVL
        CALL ERREX
      ENDIF
C
C     Special CCSD EOM
C
      IF(ICLLVL.EQ.10.OR.ICLLVL.EQ.5.OR.ICLLVL.EQ.6)THEN
       LWIC11 = 7   +   300
       LWIC12 = 8   +   300 
       LWIC13 = 9   +   300 
       LWIC14 = 10   +   300 
       LWIC15 = 27   +   300 
       LWIC16 = 28   +   300 
       LWIC17 = 29   +   300 
       LWIC18 = 30   +   300 
       LWIC21 = 7   +   300
       LWIC22 = 8   +   300
       LWIC23 = 9   +   300
       LWIC24 = 10   +   300
       LWIC25 = 27   +   300
       LWIC26 = 28   +   300
       LWIC27 = 29   +   300
       LWIC28 = 30   +   300
      ENDIF
C            --- basic D3T3 = WT2 contraction ---
C
      IF(ICLLVL.EQ.13.OR.ICLLVL.EQ.14.OR.ICLLVL.EQ.16.OR.
     &   ICLLVL.EQ.22.OR.ICLLVL.EQ.33.OR.ICLLVL.EQ.34)THEN
C
      LWIC11 =   7 + 300
      LWIC12 =   8 + 300
      LWIC13 =   9 + 300
      LWIC14 =  10 + 300
C
      LWIC15 =  27 + 300
      LWIC16 =  28 + 300
      LWIC17 =  29 + 300
      LWIC18 =  30 + 300
C
      ENDIF
C
C              --- T3 inclusion in T2 equation ---
C
      IF( (ICLLVL.EQ.22.AND. (IFLAGS2(124).EQ.1 .OR. IFLAGS2(124).EQ.2))
     &    .OR.  ICLLVL.EQ.13 )THEN
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
      IF( (ICLLVL.EQ.22 .AND. IFLAGS2(124).GE.3) .OR.
     &     ICLLVL.EQ.14 .OR.  ICLLVL.EQ.16 .OR.  ICLLVL.EQ.18 .OR.
     &     ICLLVL.EQ.33 .OR.  ICLLVL.EQ.34)THEN
C
      LWIC21 =   7
      LWIC22 =   8
      LWIC23 =   9
      LWIC24 =  10
C
      LWIC25 =  27
      LWIC26 =  28
      LWIC27 =  29
      LWIC28 =  30
C
      ENDIF
C
      RETURN
      END
