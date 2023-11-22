      SUBROUTINE MKOFVO
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION LENVO(8,8,4)
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
C     ROUTINE TO COMPUTE VARIOUS VO OFFSETS.
C
C     #     V,O
C
C     1     A,A
C     2     B,B
C     3     B,A
C     4     A,B
C
      DO   30 ISPIN=1,2
      DO   20 IRPVO=1,NIRREP
      DO   10 IRPO =1,NIRREP
      IRPV = DIRPRD(IRPVO,IRPO)
C
      LENVO(IRPO,IRPVO,ISPIN) =  POP(IRPO,ISPIN) * VRT(IRPV,ISPIN)
C
      IF(IRPO.EQ.1)THEN
      IOFFVO(IRPO,IRPVO,ISPIN) = 0
      ELSE
      IOFFVO(IRPO,IRPVO,ISPIN) = IOFFVO(IRPO-1,IRPVO,ISPIN) + 
     1                            LENVO(IRPO-1,IRPVO,ISPIN)
      ENDIF
   10 CONTINUE
   20 CONTINUE
   30 CONTINUE
C
      DO  120 IRPVO=1,NIRREP
      DO  110 IRPO =1,NIRREP
      IRPV = DIRPRD(IRPVO,IRPO)
C
      LENVO(IRPO,IRPVO,3) =  POP(IRPO,1) * VRT(IRPV,2)
C
      IF(IRPO.EQ.1)THEN
      IOFFVO(IRPO,IRPVO,3) = 0
      ELSE
      IOFFVO(IRPO,IRPVO,3) = IOFFVO(IRPO-1,IRPVO,3) + 
     1                        LENVO(IRPO-1,IRPVO,3)
      ENDIF
  110 CONTINUE
  120 CONTINUE
C
      DO  220 IRPVO=1,NIRREP
      DO  210 IRPO =1,NIRREP
      IRPV = DIRPRD(IRPVO,IRPO)
C
      LENVO(IRPO,IRPVO,4) =  POP(IRPO,2) * VRT(IRPV,1)
C
      IF(IRPO.EQ.1)THEN
      IOFFVO(IRPO,IRPVO,4) = 0
      ELSE
      IOFFVO(IRPO,IRPVO,4) = IOFFVO(IRPO-1,IRPVO,4) + 
     1                        LENVO(IRPO-1,IRPVO,4)
      ENDIF
  210 CONTINUE
  220 CONTINUE
      RETURN
      END
