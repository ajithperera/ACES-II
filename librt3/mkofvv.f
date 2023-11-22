      SUBROUTINE MKOFVV
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION LENVV(8,8)
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
C     ROUTINE TO COMPUTE VARIOUS VV OFFSETS.
C
C     #    SPIN   ORDER
C
C     1     AA     A<B
C     2     BB     A<B
C     3     AA     A,B
C     4     BB     A,B
C     5     AB     A,B
C     6     BA     A,B
C
      DO   40 ISPIN=1,2
      DO   10 IRPBC=1,NIRREP
      DO   10 IRPC =1,NIRREP
       LENVV(IRPC,IRPBC)       = 0
      IOFFVV(IRPC,IRPBC,ISPIN) = 0
   10 CONTINUE
C
      DO   30 IRPBC=1,NIRREP
      DO   20 IRPC =1,NIRREP
      IRPB = DIRPRD(IRPBC,IRPC)
C
      IF(IRPB.GT.IRPC)THEN
      LENVV(IRPC,IRPBC) = 0
      ENDIF
C
      IF(IRPB.EQ.IRPC)THEN
      LENVV(IRPC,IRPBC) = 
     1             (VRT(IRPC,ISPIN)-1)*VRT(IRPC,ISPIN)/2
      ENDIF
      IF(IRPB.LT.IRPC)THEN
      LENVV(IRPC,IRPBC) = 
     1              VRT(IRPC,ISPIN)   *VRT(IRPB,ISPIN)
      ENDIF
C
      IF(IRPC.EQ.1)THEN
      IOFFVV(IRPC,IRPBC,ISPIN) = 0
      ELSE
      IOFFVV(IRPC,IRPBC,ISPIN) = IOFFVV(IRPC-1,IRPBC,ISPIN) + 
     1                            LENVV(IRPC-1,IRPBC)
      ENDIF
   20 CONTINUE
   30 CONTINUE
   40 CONTINUE
C
      DO  140 ISPIN=3,4
      DO  110 IRPBC=1,NIRREP
      DO  110 IRPC =1,NIRREP
       LENVV(IRPC,IRPBC)       = 0
      IOFFVV(IRPC,IRPBC,ISPIN) = 0
  110 CONTINUE
C
      DO  130 IRPBC=1,NIRREP
      DO  120 IRPC =1,NIRREP
      IRPB = DIRPRD(IRPBC,IRPC)
C
      LENVV(IRPC,IRPBC) = VRT(IRPB,ISPIN-2) *
     1                    VRT(IRPC,ISPIN-2)
C
      IF(IRPC.EQ.1)THEN
      IOFFVV(IRPC,IRPBC,ISPIN) = 0
      ELSE
      IOFFVV(IRPC,IRPBC,ISPIN) = IOFFVV(IRPC-1,IRPBC,ISPIN) + 
     1                            LENVV(IRPC-1,IRPBC)
      ENDIF
  120 CONTINUE
  130 CONTINUE
  140 CONTINUE
C
      DO  210 IRPAB=1,NIRREP
      DO  210 IRPB =1,NIRREP
      IRPA = DIRPRD(IRPAB,IRPB)
C
      LENVV(IRPB,IRPAB) =  VRT(IRPA,1) * VRT(IRPB,2)
C
      IF(IRPB.EQ.1)THEN
      IOFFVV(IRPB,IRPAB,5) = 0
      ELSE
      IOFFVV(IRPB,IRPAB,5) = IOFFVV(IRPB-1,IRPAB,5) + 
     1                        LENVV(IRPB-1,IRPAB)
      ENDIF
  210 CONTINUE
C
      DO  310 IRPAB=1,NIRREP
      DO  310 IRPB =1,NIRREP
      IRPA = DIRPRD(IRPAB,IRPB)
C
      LENVV(IRPB,IRPAB) =  VRT(IRPA,2) * VRT(IRPB,1)
C
      IF(IRPB.EQ.1)THEN
      IOFFVV(IRPB,IRPAB,6) = 0
      ELSE
      IOFFVV(IRPB,IRPAB,6) = IOFFVV(IRPB-1,IRPAB,6) + 
     1                        LENVV(IRPB-1,IRPAB)
      ENDIF
  310 CONTINUE
C
      RETURN
      END
