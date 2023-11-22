      SUBROUTINE SETLWI(METHOD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
     1                LWIC15,LWIC16,LWIC17,LWIC18,
     1                LWIC21,LWIC22,LWIC23,LWIC24,
     1                LWIC25,LWIC26,LWIC27,LWIC28,
     1                LWIC31,LWIC32,LWIC33,
     1                LWIC34,LWIC35,LWIC36,
     1                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
      IF(METHOD.EQ.12.OR.(METHOD.GE.26.AND.METHOD.LE.31))THEN
      LWIC11 = 107
      LWIC12 = 108
      LWIC13 = 109
      LWIC14 = 110
      LWIC15 = 127
      LWIC16 = 128
      LWIC17 = 129
      LWIC18 = 130
C     This is a fudge so TSPABCI does not do anything.
      LWIC21 = 107
      LWIC22 = 108
      LWIC23 = 109
      LWIC24 = 110
      LWIC25 = 127
      LWIC26 = 128
      LWIC27 = 129
      LWIC28 = 130
      ENDIF
C
      IF(METHOD.EQ.22)THEN
      LWIC11 =   7
      LWIC12 =   8
      LWIC13 =   9
      LWIC14 =  10
      LWIC15 =  27
      LWIC16 =  28
      LWIC17 =  29
      LWIC18 =  30
C     This is a fudge so TSPABCI does not do anything.
      LWIC21 =   7
      LWIC22 =   8
      LWIC23 =   9
      LWIC24 =  10
      LWIC25 =  27
      LWIC26 =  28
      LWIC27 =  29
      LWIC28 =  30
      ENDIF
      RETURN
      END
