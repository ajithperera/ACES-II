      SUBROUTINE W5T1ABCD(CORE,MAXCOR,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,G,AG,AM,DIRPRD,VRT,POP
      DIMENSION CORE(1)
      DIMENSION IOFFVV(8,8,3),IOFFVO(8,8,3),IOFFOV(8,8)
      DIMENSION  LENVV(8,8,3), LENVO(8,8,3), LENOV(8,8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYM   / POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,
     &                NF2BB
C
      INDEX(I) = (I*(I-1))/2
C
C     This is a straightforward, slow, clear (?!) routine to perform the 
C     contraction
C
C     W(ef,am) = t(g,m) * <ef||ag>
C
C     This term must be added to the W(ef,am) formed in lambda in order to
C     complete the Hbar element.
C
      CALL IZERO(IOFFVV,8*8*3)
      CALL IZERO(IOFFVO,8*8*3)
      CALL IZERO(IOFFOV,8*8  )
      CALL IZERO( LENVV,8*8*3)
      CALL IZERO( LENVO,8*8*3)
      CALL IZERO( LENOV,8*8  )
C
      DO    20 IRPPQ=1,NIRREP
      DO    10 IRPQ =1,NIRREP
C
      IRPP = DIRPRD(IRPQ,IRPPQ)
C
      IF(IRPP.LE.IRPQ)THEN
       IF(IRPP.EQ.IRPQ)THEN
        LENVV(IRPQ,IRPPQ,1) = (VRT(IRPP,1) * (VRT(IRPP,1) - 1))/2
        LENVV(IRPQ,IRPPQ,2) = (VRT(IRPP,2) * (VRT(IRPP,2) - 1))/2
       ELSE
        LENVV(IRPQ,IRPPQ,1) =  VRT(IRPP,1) *  VRT(IRPQ,1)
        LENVV(IRPQ,IRPPQ,2) =  VRT(IRPP,2) *  VRT(IRPQ,2)
       ENDIF
      ENDIF
C
      LENVV(IRPQ,IRPPQ,3)  =  VRT(IRPP,1) *  VRT(IRPQ,2)
C
      LENVO(IRPQ,IRPPQ,1)  =  VRT(IRPP,1) *  POP(IRPQ,1)
      LENVO(IRPQ,IRPPQ,2)  =  VRT(IRPP,2) *  POP(IRPQ,2)
      LENVO(IRPQ,IRPPQ,3)  =  VRT(IRPP,1) *  POP(IRPQ,2)
C
      LENOV(IRPQ,IRPPQ  )  =  POP(IRPP,1) *  VRT(IRPQ,2)
   10 CONTINUE
   20 CONTINUE
C
      DO    40 IRPPQ=1,NIRREP
      DO    30 IRPQ =1,NIRREP
C
      IRPP = DIRPRD(IRPQ,IRPPQ)
C
      IF(IRPQ.EQ.1)THEN
       IOFFVV(IRPQ,IRPPQ,3) = 0
       IOFFVO(IRPQ,IRPPQ,1) = 0
       IOFFVO(IRPQ,IRPPQ,2) = 0
       IOFFVO(IRPQ,IRPPQ,3) = 0
       IOFFOV(IRPQ,IRPPQ  ) = 0
      ELSE
       IOFFVV(IRPQ,IRPPQ,3) = 
     &                    IOFFVV(IRPQ-1,IRPPQ,3) + LENVV(IRPQ-1,IRPPQ,3)
       IOFFVO(IRPQ,IRPPQ,1) =
     &                    IOFFVO(IRPQ-1,IRPPQ,1) + LENVO(IRPQ-1,IRPPQ,1)
       IOFFVO(IRPQ,IRPPQ,2) =
     &                    IOFFVO(IRPQ-1,IRPPQ,2) + LENVO(IRPQ-1,IRPPQ,2)
       IOFFVO(IRPQ,IRPPQ,3) =
     &                    IOFFVO(IRPQ-1,IRPPQ,3) + LENVO(IRPQ-1,IRPPQ,3)
       IOFFOV(IRPQ,IRPPQ  ) =
     &                    IOFFOV(IRPQ-1,IRPPQ  ) + LENOV(IRPQ-1,IRPPQ  )
      ENDIF
C
      IF(IRPQ.EQ.1)THEN
       IOFFVV(IRPQ,IRPPQ,1) = 0
       IOFFVV(IRPQ,IRPPQ,2) = 0
      ELSE
       IOFFVV(IRPQ,IRPPQ,1) = 
     &                 IOFFVV(IRPQ-1,IRPPQ,1) + LENVV(IRPQ-1,IRPPQ,1)
       IOFFVV(IRPQ,IRPPQ,2) = 
     &                 IOFFVV(IRPQ-1,IRPPQ,2) + LENVV(IRPQ-1,IRPPQ,2)
      ENDIF
C
   30 CONTINUE
   40 CONTINUE
C
      I000    =   1
      LW5OFF  = 126
      LINTOFF = 230
C
C
C     ----- T(G,M) * <EF || AG> , T(g,m) * <ef || ag> -----
C
      IF(IUHF.NE.0)THEN
C
      DO  170 ISPIN=1,IUHF+1
C
      ISPIN1  = ISPIN
      ISPIN2  = ISPIN
C
      CALL GETLST(CORE(I000),1,1,1,ISPIN,90)
C
      DO  160 IRPAM=1,NIRREP
C
      IRPEF = IRPAM
      IRPAG = IRPAM
C
      IF(IRPDPD(IRPEF,ISPIN).EQ.0.OR.IRPDPD(IRPAM,8+ISPIN).EQ.0) 
     &                                                      GOTO 160
C
      I010 = I000 + IRPDPD(    1,8 + ISPIN)
      I020 = I010 + IRPDPD(IRPEF,    ISPIN)
      I030 = I020 + IRPDPD(IRPEF,    ISPIN)
      NEED = IINTFP * I030
C
      IF(NEED.GT.MAXCOR)THEN
       CALL INSMEM('W5T1ABCD',NEED,MAXCOR)
       CALL ERREX
      ENDIF
C
      DO  150 IRPM =1,NIRREP
C
      IRPG = IRPM
      IRPA = DIRPRD(IRPM,IRPAM)
C
      IF(POP(IRPM,ISPIN2).EQ.0.OR.VRT(IRPA,ISPIN1).EQ.0.OR.
     &                            VRT(IRPG,ISPIN2).EQ.0   ) GOTO 150
C
      IF(IRPAG.EQ.1.AND.VRT(IRPA,ISPIN1).LT.2)              GOTO 150
C
      DO  140    G=1,VRT(IRPG,ISPIN2)
      DO  130    A=1,VRT(IRPA,ISPIN1)
C
      IF(IRPAG.EQ.1.AND.A.EQ.G)                             GOTO 130
C
      IF(IRPA.EQ.IRPG)THEN
       AG = IOFFVV(IRPG,IRPAG,ISPIN) + INDEX(MAX(A,G)-1) + MIN(A,G)
      ENDIF
      IF(IRPA.LT.IRPG)THEN
       AG = IOFFVV(IRPG,IRPAG,ISPIN) + (G-1)*VRT(IRPA,ISPIN1) + A
      ENDIF
      IF(IRPA.GT.IRPG)THEN
       AG = IOFFVV(IRPA,IRPAG,ISPIN) + (A-1)*VRT(IRPG,ISPIN2) + G
      ENDIF
C
      CALL GETLST(CORE(I010),AG,1,2,IRPAG,LINTOFF + ISPIN)
C
      IF((IRPAG.EQ.1.AND.A.GT.G).OR.IRPA.GT.IRPG)THEN
       CALL VMINUS(CORE(I010),IRPDPD(IRPEF,ISPIN))
      ENDIF
C
      DO  120    M=1,POP(IRPM,ISPIN2)
C
      AM = IOFFVO(IRPM,IRPAM,ISPIN) + (M-1)*VRT(IRPA,ISPIN1) + A
      CALL GETLST(CORE(I020),AM,1,2,IRPAM,LW5OFF + ISPIN)
C
      FACT = CORE(I000 - 1 + IOFFVO(IRPM,1,ISPIN2) + 
     &                    (M-1)*VRT(IRPG,ISPIN2) + G)
C
      CALL VADD(CORE(I020),CORE(I020),CORE(I010),IRPDPD(IRPEF,ISPIN),
     &          FACT)
C
      CALL PUTLST(CORE(I020),AM,1,2,IRPAM,LW5OFF + ISPIN)
C
  120 CONTINUE
  130 CONTINUE
  140 CONTINUE
C
  150 CONTINUE
C
  160 CONTINUE
C
  170 CONTINUE
C
C     ----- T(G,M) * <Ef | Ga> -----
C
      ISPIN1  = 2
      ISPIN2  = 1
C
      CALL GETLST(CORE(I000),1,1,1,ISPIN2,90)
C
      DO  260 IRPAM=1,NIRREP
C
      IRPEF = IRPAM
      IRPAG = IRPAM
C
      IF(IRPDPD(IRPEF,13).EQ.0.OR.IRPDPD(IRPAM,12).EQ.0)    GOTO 260
C
      I010 = I000 + IRPDPD(    1,8 + ISPIN2)
      I020 = I010 + IRPDPD(IRPEF,        13)
      I030 = I020 + IRPDPD(IRPEF,        13)
      NEED = IINTFP * I030
C
      IF(NEED.GT.MAXCOR)THEN
       CALL INSMEM('W5T1ABCD',NEED,MAXCOR)
       CALL ERREX
      ENDIF
C
      DO  250 IRPM =1,NIRREP
C
      IRPG = IRPM
      IRPA = DIRPRD(IRPM,IRPAM)
C
      IF(POP(IRPM,ISPIN2).EQ.0.OR.VRT(IRPA,ISPIN1).EQ.0.OR.
     &                            VRT(IRPG,ISPIN2).EQ.0   ) GOTO 250
C
      DO  240    G=1,VRT(IRPG,ISPIN2)
      DO  230    A=1,VRT(IRPA,ISPIN1)
C
      AG = IOFFVV(IRPA,IRPAG,3) + (A-1)*VRT(IRPG,ISPIN2) + G
C
      CALL GETLST(CORE(I010),AG,1,2,IRPAG,LINTOFF + 3)
C
      DO  220    M=1,POP(IRPM,ISPIN2)
C
      AM = IOFFOV(IRPA,IRPAM) + (A-1)*POP(IRPM,ISPIN2) + M
      CALL GETLST(CORE(I020),AM,1,2,IRPAM,LW5OFF + 3)
C
      FACT = CORE(I000 - 1 + IOFFVO(IRPM,1,ISPIN2) + 
     &                    (M-1)*VRT(IRPG,ISPIN2) + G)
C
      CALL VADD(CORE(I020),CORE(I020),CORE(I010),IRPDPD(IRPEF,13),
     &          FACT)
C
      CALL PUTLST(CORE(I020),AM,1,2,IRPAM,LW5OFF + 3)
C
  220 CONTINUE
  230 CONTINUE
  240 CONTINUE
C
  250 CONTINUE
C
  260 CONTINUE
  270 continue
C
      ENDIF
C
C     ----- T(g,m) * <Ef | Ag> -----
C
      ISPIN1  = 1
      ISPIN2  = 2
C
      IF(IUHF.EQ.0)THEN
       CALL GETLST(CORE(I000),1,1,1,ISPIN1,90)
      ELSE
       CALL GETLST(CORE(I000),1,1,1,ISPIN2,90)
      ENDIF
C
      DO  360 IRPAM=1,NIRREP
C
      IRPEF = IRPAM
      IRPAG = IRPAM
C
      IF(IRPDPD(IRPEF,13).EQ.0.OR.IRPDPD(IRPAM,11).EQ.0)    GOTO 360
C
      I010 = I000 + IRPDPD(    1,8 + ISPIN2)
      I020 = I010 + IRPDPD(IRPEF,        13)
      I030 = I020 + IRPDPD(IRPEF,        13)
      NEED = IINTFP * I030
C
      IF(NEED.GT.MAXCOR)THEN
       CALL INSMEM('W5T1ABCD',NEED,MAXCOR)
       CALL ERREX
      ENDIF
C
      DO  350 IRPM =1,NIRREP
C
      IRPG = IRPM
      IRPA = DIRPRD(IRPM,IRPAM)
C
      IF(POP(IRPM,ISPIN2).EQ.0.OR.VRT(IRPA,ISPIN1).EQ.0.OR.
     &                            VRT(IRPG,ISPIN2).EQ.0   ) GOTO 350
C
      DO  340    G=1,VRT(IRPG,ISPIN2)
      DO  330    A=1,VRT(IRPA,ISPIN1)
C
      AG = IOFFVV(IRPG,IRPAG,3) + (G-1)*VRT(IRPA,ISPIN1) + A
C
      CALL GETLST(CORE(I010),AG,1,2,IRPAG,LINTOFF + 3)
C
      DO  320    M=1,POP(IRPM,ISPIN2)
C
      AM = IOFFVO(IRPM,IRPAM,3) + (M-1)*VRT(IRPA,ISPIN1) + A
      CALL GETLST(CORE(I020),AM,1,2,IRPAM,LW5OFF + 4)
C
      FACT = CORE(I000 - 1 + IOFFVO(IRPM,1,ISPIN2) + 
     &                    (M-1)*VRT(IRPG,ISPIN2) + G)
C
      CALL VADD(CORE(I020),CORE(I020),CORE(I010),IRPDPD(IRPEF,13),
     &          FACT)
C
      CALL PUTLST(CORE(I020),AM,1,2,IRPAM,LW5OFF + 4)
C
  320 CONTINUE
  330 CONTINUE
  340 CONTINUE
C
  350 CONTINUE
C
  360 CONTINUE
C
      RETURN
      END
