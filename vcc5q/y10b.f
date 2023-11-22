      SUBROUTINE Y10B(CORE,MAXCOR,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER B,D,E,F,BD,EF,BF,ED,DF,EB,BDEF,BFED,DFEB
      INTEGER DISSIZ,DISSIZAA,DISSIZBB,DISSIZAB,DISSIZT
      INTEGER POP,VRT,DIRPRD
      DIMENSION CORE(1)
      DIMENSION DISSIZ(8),DISSIZAA(8),DISSIZBB(8),DISSIZAB(8),
     1          DISSIZT(8)
      DIMENSION NDIST(8),NSUM(8)
      DIMENSION IOFFY5(8),IOFFY5AA(8),IOFFY5BB(8),IOFFY5AB(8),
     1          IOFFT(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/    POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /INFO/   NOCA,NOCB,NVRTA,NVRTB
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
      INDEX(I) = I * (I-1)/2
C
      WRITE(6,1000)
 1000 FORMAT(' @Y10B-I, Y5 * Y5 contributions. ')
c1000 FORMAT(' @Y10B-I, Welcome to one of the monsters ! ')
C
C     AAAA/BBBB
C
      IF(IUHF.GT.0)THEN
C
      DO  400 ISPIN=1,2
C
      LEN = 0
      DO   20 IRREP=1,NIRREP
      DISSIZ(IRREP)  =  IRPDPD(IRREP,18+ISPIN)
      LEN = LEN + DISSIZ(IRREP) * DISSIZ(IRREP)
      IF(IRREP.EQ.1)THEN
      IOFFY5(IRREP) = 0
      ELSE
      IOFFY5(IRREP) = IOFFY5(IRREP-1) + DISSIZ(IRREP-1)*DISSIZ(IRREP-1)
      ENDIF
   20 CONTINUE
C
      I000 = 1
      I010 = I000 + LEN
      I020 = I010 + LEN
      I030 = I020 + 2 * (NVRTA * NVRTB + NVRTA * NVRTB)
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
 1020 FORMAT(' @Y10B-I, Insufficient memory. Need ',I15,' . Got ',I15)
      STOP 'Y10B'
      ENDIF
C
C     Read abcf intermediates and expand each symmetry block in turn to
C     (a,b,c,f).
C
      DO   30 IRREP=1,NIRREP
      CALL GETLIST(CORE(I010+IOFFY5(IRREP)),
     1             1,IRPDPD(IRREP,ISPIN),2,IRREP,ISPIN)
C
C     (a<b,c<f) ---> (a,b,c<f)
      CALL SYMEXP2(IRREP,VRT(1,ISPIN),DISSIZ(IRREP),
     1             IRPDPD(IRREP,ISPIN),IRPDPD(IRREP,ISPIN),
     1             CORE(I010+IOFFY5(IRREP)),CORE(I010+IOFFY5(IRREP)))
C     (a,b,c<f) ---> (a,b,c,f)
      CALL  SYMEXP(IRREP,VRT(1,ISPIN),DISSIZ(IRREP),
     1             CORE(I010+IOFFY5(IRREP)))
C
   30 CONTINUE
C
C     Now reorder (a,b,c,f) to (a,c,b,f). This is also effectively
C     reordering (a,e,c,d) to (a,c,e,d).
C
      CALL SSTGEN(CORE(I010),CORE(I000),LEN,
     1            VRT(1,ISPIN),VRT(1,ISPIN),VRT(1,ISPIN),VRT(1,ISPIN),
     1            CORE(I020),1,'1324')
C
C     Now perform massive matrix multiply (ac,bf) * (ac,ed) = (bf,ed)
C
C     Y5 intermediates are now at I000. Target can be at I010. Allocation
C     is same.
C
      DO   40 IRPAC=1,NIRREP
C
      CALL XGEMM('T','N',DISSIZ(IRPAC),DISSIZ(IRPAC),DISSIZ(IRPAC),
     1            1.0D+00,
     1           CORE(I000+IOFFY5(IRPAC)),DISSIZ(IRPAC),
     1           CORE(I000+IOFFY5(IRPAC)),DISSIZ(IRPAC),
     1            0.0D+00,
     1           CORE(I010+IOFFY5(IRPAC)),DISSIZ(IRPAC))
C
   40 CONTINUE
C
C     Now loop over symmetry blocks of Y10 and insert appropriate pieces.
C
      DO  190 IRPBD=1,NIRREP
C
      I030 = I020 + IRPDPD(IRPBD,ISPIN) * IRPDPD(IRPBD,ISPIN)
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10B'
      ENDIF
C
      CALL GETLIST(CORE(I020),1,IRPDPD(IRPBD,ISPIN),2,IRPBD,3+ISPIN)
C
      IRPEF = IRPBD
C
      DO  180 IRPD =1,NIRREP
C
      IRPB = DIRPRD(IRPD,IRPBD)
      IF(IRPB.GT.IRPD) GOTO 180
C
      DO  170 IRPF =1,NIRREP
C
      IRPE = DIRPRD(IRPF,IRPEF)
      IF(IRPE.GT.IRPF) GOTO 170
C
      IRPED = DIRPRD(IRPE,IRPD)
      IRPBF = DIRPRD(IRPB,IRPF)
      IRPDF = DIRPRD(IRPF,IRPD)
      IRPEB = DIRPRD(IRPE,IRPB)
C
      IF(IRPBD.EQ.1)THEN
C
      IF(VRT(IRPD,ISPIN).GE.2.AND.VRT(IRPF,ISPIN).GE.2)THEN
C
      DO   80 F=2,VRT(IRPF,ISPIN)
      DO   70 E=1,F-1
      DO   60 D=2,VRT(IRPD,ISPIN)
      DO   50 B=1,D-1
C
      BD = IOFFVV(IRPD,IRPBD,ISPIN) + INDEX(D-1) + B
      EF = IOFFVV(IRPF,IRPEF,ISPIN) + INDEX(F-1) + E
C
      BDEF = (EF-1) * IRPDPD(IRPBD,ISPIN) + BD
C
      BF = IOFFVV(IRPF,IRPBF,2+ISPIN) + (F-1)*VRT(IRPB,ISPIN) + B
      ED = IOFFVV(IRPD,IRPED,2+ISPIN) + (D-1)*VRT(IRPE,ISPIN) + E
      DF = IOFFVV(IRPF,IRPDF,2+ISPIN) + (F-1)*VRT(IRPD,ISPIN) + D
      EB = IOFFVV(IRPB,IRPEB,2+ISPIN) + (B-1)*VRT(IRPE,ISPIN) + E
C
      BFED = I010 + IOFFY5(IRPBF) + (ED-1)*DISSIZ(IRPBF) + BF - 1
      DFEB = I010 + IOFFY5(IRPDF) + (EB-1)*DISSIZ(IRPDF) + DF - 1
C
      CORE(I020 + BDEF - 1) = CORE(I020 + BDEF - 1)
     1                      + CORE(BFED) - CORE(DFEB)
C
   50 CONTINUE
   60 CONTINUE
   70 CONTINUE
   80 CONTINUE
C
      ENDIF
C
      ELSE
C
      IF(VRT(IRPB,ISPIN).GE.1.AND.VRT(IRPD,ISPIN).GE.1.AND.
     1   VRT(IRPE,ISPIN).GE.1.AND.VRT(IRPF,ISPIN).GE.1)THEN
C
      DO  120 F=1,VRT(IRPF,ISPIN)
      DO  110 E=1,VRT(IRPE,ISPIN)
      DO  100 D=1,VRT(IRPD,ISPIN)
      DO   90 B=1,VRT(IRPB,ISPIN)
C
      BD = IOFFVV(IRPD,IRPBD,ISPIN) + (D-1)*VRT(IRPB,ISPIN) + B
      EF = IOFFVV(IRPF,IRPEF,ISPIN) + (F-1)*VRT(IRPE,ISPIN) + E
C
      BDEF = (EF-1) * IRPDPD(IRPBD,ISPIN) + BD
C
      BF = IOFFVV(IRPF,IRPBF,2+ISPIN) + (F-1)*VRT(IRPB,ISPIN) + B
      ED = IOFFVV(IRPD,IRPED,2+ISPIN) + (D-1)*VRT(IRPE,ISPIN) + E
      DF = IOFFVV(IRPF,IRPDF,2+ISPIN) + (F-1)*VRT(IRPD,ISPIN) + D
      EB = IOFFVV(IRPB,IRPEB,2+ISPIN) + (B-1)*VRT(IRPE,ISPIN) + E
C
      BFED = I010 + IOFFY5(IRPBF) + (ED-1)*DISSIZ(IRPBF) + BF - 1
      DFEB = I010 + IOFFY5(IRPDF) + (EB-1)*DISSIZ(IRPDF) + DF - 1
C
      CORE(I020 + BDEF - 1) = CORE(I020 + BDEF - 1)
     1                      + CORE(BFED) - CORE(DFEB)
C
   90 CONTINUE
  100 CONTINUE
  110 CONTINUE
  120 CONTINUE
C
      ENDIF
C
      ENDIF
C
  170 CONTINUE
  180 CONTINUE
C
      CALL PUTLIST(CORE(I020),1,IRPDPD(IRPBD,ISPIN),2,IRPBD,3+ISPIN)
C
  190 CONTINUE
C
C     alpha-beta contribution.
C
      IF(ISPIN.EQ.1)THEN
      ISPIN1 = 1
      ISPIN2 = 2
      ELSE
      ISPIN1 = 2
      ISPIN2 = 1
      ENDIF
C
      LEN  = 0
      LEN2 = 0
      LENT = 0
      DO  220 IRREP=1,NIRREP
      DISSIZAB(IRREP) = IRPDPD(IRREP,13)
       DISSIZT(IRREP) = IRPDPD(IRREP,18+ISPIN)
         NDIST(IRREP) = IRPDPD(IRREP,18+ISPIN)
          NSUM(IRREP) = IRPDPD(IRREP,21-ISPIN)
      LEN  = LEN  + DISSIZAB(IRREP) * DISSIZAB(IRREP)
      LEN2 = LEN2 +  DISSIZT(IRREP) *     NSUM(IRREP)
      LENT = LENT +  DISSIZT(IRREP) *    NDIST(IRREP)
      IF(IRREP.EQ.1)THEN
      IOFFY5AB(IRREP)  = 0
      IOFFY5(IRREP)    = 0
      IOFFT(IRREP)     = 0
      ELSE
      IOFFY5AB(IRREP)  = IOFFY5AB(IRREP-1) + 
     1                   DISSIZAB(IRREP-1)*DISSIZAB(IRREP-1)
      IOFFY5(IRREP)    = IOFFY5(IRREP-1) + 
     1                    DISSIZT(IRREP-1)*    NSUM(IRREP-1)
      IOFFT(IRREP)     = IOFFT(IRREP-1) +
     1                    DISSIZT(IRREP-1)*   NDIST(IRREP-1)
      ENDIF
  220 CONTINUE
C
      IF(LEN.NE.LEN2)THEN
      WRITE(6,1030)
 1030 FORMAT(' @Y10B-F, A screwup... No way anything can be right. ')
      CALL ERREX
      ENDIF
C
      I000 = 1
      I010 = I000 + LEN
      I020 = I010 + LEN
      I030 = I020 + 2 * (NVRTA * NVRTB + NVRTA * NVRTB)
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10B'
      ENDIF
C
C     Read abcf intermediates.
C
      DO  230 IRREP=1,NIRREP
      CALL GETLIST(CORE(I010+IOFFY5AB(IRREP)),
     1             1,IRPDPD(IRREP,13),2,IRREP,3)
  230 CONTINUE
C
C     ISPIN = 1 : (B,a,F,c), (E,a,D,c) ---> (B,F,a,c), (E,D,a,c)
C     ISPIN = 2 : (A,b,C,f), (A,e,C,d) ---> (A,C,b,f), (A,C,e,d)
C
      CALL SSTGEN(CORE(I010),CORE(I000),LEN,
     1            VRT(1,1),VRT(1,2),VRT(1,1),VRT(1,2),
     1            CORE(I020),1,'1324')
C
C     Now perform massive matrix multiply :
C
C     ISPIN = 1 : (BF,ac) * (ED,ac) = (BF,ED)
C     ISPIN = 2 : (AC,bf) * (AC,ed) = (bf,ed)
C
C     Y5 intermediates are now at I000. Target can be at I010.
C
      I020 = I010 + LENT
C
      NEED = IINTFP * I020
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10B'
      ENDIF
C
      DO  240 IRPAC=1,NIRREP
C
      IF(ISPIN.EQ.1)THEN
      CALL XGEMM('N','T',DISSIZT(IRPAC),NDIST(IRPAC),NSUM(IRPAC),
     1            1.0D+00,
     1           CORE(I000+IOFFY5(IRPAC)),DISSIZT(IRPAC),
     1           CORE(I000+IOFFY5(IRPAC)),  NDIST(IRPAC),
     1            0.0D+00,
     1           CORE(I010+IOFFT(IRPAC)),DISSIZT(IRPAC))
      ELSE
      CALL XGEMM('T','N',DISSIZT(IRPAC),NDIST(IRPAC),NSUM(IRPAC),
     1            1.0D+00,
     1           CORE(I000+IOFFY5(IRPAC)),NSUM(IRPAC),
     1           CORE(I000+IOFFY5(IRPAC)),NSUM(IRPAC),
     1            0.0D+00,
     1           CORE(I010+IOFFT(IRPAC)),DISSIZT(IRPAC))
      ENDIF
C
  240 CONTINUE
C
C     Now loop over symmetry blocks of Y10 and insert appropriate pieces.
C
      DO  390 IRPBD=1,NIRREP
C
      I030 = I020 + IRPDPD(IRPBD,ISPIN) * IRPDPD(IRPBD,ISPIN)
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10B'
      ENDIF
C
      CALL GETLIST(CORE(I020),1,IRPDPD(IRPBD,ISPIN),2,IRPBD,3+ISPIN)
C
      IRPEF = IRPBD
C
      DO  380 IRPD =1,NIRREP
C
      IRPB = DIRPRD(IRPD,IRPBD)
      IF(IRPB.GT.IRPD) GOTO 380
C
      DO  370 IRPF =1,NIRREP
C
      IRPE = DIRPRD(IRPF,IRPEF)
      IF(IRPE.GT.IRPF) GOTO 370
C
      IRPED = DIRPRD(IRPE,IRPD)
      IRPBF = DIRPRD(IRPB,IRPF)
      IRPDF = DIRPRD(IRPF,IRPD)
      IRPEB = DIRPRD(IRPE,IRPB)
C
      IF(IRPBD.EQ.1)THEN
C
      IF(VRT(IRPD,ISPIN).GE.2.AND.VRT(IRPF,ISPIN).GE.2)THEN
C
      DO  280 F=2,VRT(IRPF,ISPIN)
      DO  270 E=1,F-1
      DO  260 D=2,VRT(IRPD,ISPIN)
      DO  250 B=1,D-1
C
      BD = IOFFVV(IRPD,IRPBD,ISPIN) + INDEX(D-1) + B
      EF = IOFFVV(IRPF,IRPEF,ISPIN) + INDEX(F-1) + E
C
      BDEF = (EF-1) * IRPDPD(IRPBD,ISPIN) + BD
C
      BF = IOFFVV(IRPF,IRPBF,2+ISPIN) + (F-1)*VRT(IRPB,ISPIN) + B
      ED = IOFFVV(IRPD,IRPED,2+ISPIN) + (D-1)*VRT(IRPE,ISPIN) + E
      DF = IOFFVV(IRPF,IRPDF,2+ISPIN) + (F-1)*VRT(IRPD,ISPIN) + D
      EB = IOFFVV(IRPB,IRPEB,2+ISPIN) + (B-1)*VRT(IRPE,ISPIN) + E
C
c      BFED = I010 + IOFFY5(IRPBF) + (ED-1)*DISSIZ(IRPBF) + BF - 1
c      DFEB = I010 + IOFFY5(IRPDF) + (EB-1)*DISSIZ(IRPDF) + DF - 1
      BFED = I010 + IOFFT(IRPBF) + (ED-1)*DISSIZ(IRPBF) + BF - 1
      DFEB = I010 + IOFFT(IRPDF) + (EB-1)*DISSIZ(IRPDF) + DF - 1
C
      CORE(I020 + BDEF - 1) = CORE(I020 + BDEF - 1)
     1                      + CORE(BFED) - CORE(DFEB)
C
  250 CONTINUE
  260 CONTINUE
  270 CONTINUE
  280 CONTINUE
C
      ENDIF
C
      ELSE
C
      IF(VRT(IRPB,ISPIN).GE.1.AND.VRT(IRPD,ISPIN).GE.1.AND.
     1   VRT(IRPE,ISPIN).GE.1.AND.VRT(IRPF,ISPIN).GE.1)THEN
C
      DO  320 F=1,VRT(IRPF,ISPIN)
      DO  310 E=1,VRT(IRPE,ISPIN)
      DO  300 D=1,VRT(IRPD,ISPIN)
      DO  290 B=1,VRT(IRPB,ISPIN)
C
      BD = IOFFVV(IRPD,IRPBD,ISPIN) + (D-1)*VRT(IRPB,ISPIN) + B
      EF = IOFFVV(IRPF,IRPEF,ISPIN) + (F-1)*VRT(IRPE,ISPIN) + E
C
      BDEF = (EF-1) * IRPDPD(IRPBD,ISPIN) + BD
C
      BF = IOFFVV(IRPF,IRPBF,2+ISPIN) + (F-1)*VRT(IRPB,ISPIN) + B
      ED = IOFFVV(IRPD,IRPED,2+ISPIN) + (D-1)*VRT(IRPE,ISPIN) + E
      DF = IOFFVV(IRPF,IRPDF,2+ISPIN) + (F-1)*VRT(IRPD,ISPIN) + D
      EB = IOFFVV(IRPB,IRPEB,2+ISPIN) + (B-1)*VRT(IRPE,ISPIN) + E
C
c      BFED = I010 + IOFFY5(IRPBF) + (ED-1)*DISSIZ(IRPBF) + BF - 1
c      DFEB = I010 + IOFFY5(IRPDF) + (EB-1)*DISSIZ(IRPDF) + DF - 1
      BFED = I010 + IOFFT(IRPBF) + (ED-1)*DISSIZ(IRPBF) + BF - 1
      DFEB = I010 + IOFFT(IRPDF) + (EB-1)*DISSIZ(IRPDF) + DF - 1
C
      CORE(I020 + BDEF - 1) = CORE(I020 + BDEF - 1)
     1                      + CORE(BFED) - CORE(DFEB)
C
  290 CONTINUE
  300 CONTINUE
  310 CONTINUE
  320 CONTINUE
C
      ENDIF
C
      ENDIF
C
  370 CONTINUE
  380 CONTINUE
C
      CALL PUTLIST(CORE(I020),1,IRPDPD(IRPBD,ISPIN),2,IRPBD,3+ISPIN)
C
  390 CONTINUE
C
  400 CONTINUE
C
      ENDIF
C
C     ABAB
C
C     BdEf = aBCf * aECd
C     (BaCf * EaCd ---> BfCa * EdCa = BfEd ---> BdEf)
C
      LEN  = 0
      LENT = 0
      DO  420 IRREP=1,NIRREP
      DISSIZAB(IRREP) = IRPDPD(IRREP,13)
       DISSIZT(IRREP) = IRPDPD(IRREP,13)
         NDIST(IRREP) = IRPDPD(IRREP,13)
          NSUM(IRREP) = IRPDPD(IRREP,13)
      LEN  = LEN  + DISSIZAB(IRREP) * DISSIZAB(IRREP)
      LENT = LENT +  DISSIZT(IRREP) *    NDIST(IRREP)
      IF(IRREP.EQ.1)THEN
      IOFFY5AB(IRREP)  = 0
      IOFFT(IRREP)     = 0
      ELSE
      IOFFY5AB(IRREP)  = IOFFY5AB(IRREP-1) + 
     1                   DISSIZAB(IRREP-1)*DISSIZAB(IRREP-1)
      IOFFT(IRREP)     = IOFFT(IRREP-1)    +
     1                   DISSIZT(IRREP-1) *NDIST(IRREP-1)
      ENDIF
  420 CONTINUE
C
      I000 = 1
      I010 = I000 + LEN
      I020 = I010 + LEN
      I030 = I020 + 2 * (NVRTA * NVRTB + NVRTA * NVRTB)
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10B'
      ENDIF
C
C     Read BaCf/EaCd intermediates.
C
      DO  430 IRREP=1,NIRREP
      CALL GETLIST(CORE(I010+IOFFY5AB(IRREP)),
     1             1,IRPDPD(IRREP,13),2,IRREP,3)
  430 CONTINUE
C
C     (B,a,C,f), (E,a,C,d) ---> (B,f,C,a), (E,d,C,a)
C
      CALL SSTGEN(CORE(I010),CORE(I000),LEN,
     1            VRT(1,1),VRT(1,2),VRT(1,1),VRT(1,2),
     1            CORE(I020),1,'1432')
C
C     Now perform massive matrix multiply :
C
C     (Bf,Ca) * (Ed,Ca) = (Bf,Ed)
C
C     Y5 intermediates are now at I000. Target can be at I010.
C
      I020 = I010 + LENT
C
      NEED = IINTFP * I020
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10B'
      ENDIF
C
      DO  440 IRPAC=1,NIRREP
C
      CALL XGEMM('N','T',DISSIZT(IRPAC),NDIST(IRPAC),NSUM(IRPAC),
     1            1.0D+00,
     1           CORE(I000+IOFFY5AB(IRPAC)),DISSIZT(IRPAC),
     1           CORE(I000+IOFFY5AB(IRPAC)),NDIST(IRPAC),
     1            0.0D+00,
     1           CORE(I010+IOFFT(IRPAC)),DISSIZT(IRPAC))
C
  440 CONTINUE
C
C     We have matrix (Bf,Ed). We want (Bd,Ef). Apply SSTGEN and sum
C     directly to symmetry blocks of Y10.
C
C     (B,f,E,d) ---> (B,d,E,f)
C
      CALL SSTGEN(CORE(I010),CORE(I000),LEN,
     1            VRT(1,1),VRT(1,2),VRT(1,1),VRT(1,2),
     1            CORE(I020),1,'1432')
C
      DO  590 IRPBD=1,NIRREP
C
      CALL GETLIST(CORE(I010),1,IRPDPD(IRPBD,13),2,IRPBD, 6)
      CALL VADD(CORE(I010),CORE(I010),CORE(I000 + IOFFT(IRPBD)),
     1          IRPDPD(IRPBD,13)*IRPDPD(IRPBD,13), 1.0D+00)
      CALL PUTLIST(CORE(I010),1,IRPDPD(IRPBD,13),2,IRPBD, 6)
C
  590 CONTINUE
C
C     BdEf = -BcEa * adcf = + BcEa * dacf
C     (BcEa * dacf ---> BEca * dfca = BEdf ---> BdEf
C
      LENAA = 0
      LENBB = 0
      LENAB = 0
      LENT  = 0
      DO  620 IRREP=1,NIRREP
      DISSIZAA(IRREP) = IRPDPD(IRREP,19)
      DISSIZBB(IRREP) = IRPDPD(IRREP,20)
      DISSIZAB(IRREP) = IRPDPD(IRREP,13)
      DISSIZT(IRREP)  = IRPDPD(IRREP,19)
      NDIST(IRREP)    = IRPDPD(IRREP,20)
      NSUM(IRREP)     = IRPDPD(IRREP,20)
      LENAA = LENAA + DISSIZAA(IRREP)  * DISSIZAA(IRREP)
      LENBB = LENBB + DISSIZBB(IRREP)  * DISSIZBB(IRREP)
      LENAB = LENAB + DISSIZAB(IRREP)  * DISSIZAB(IRREP)
      LENT  = LENT  + DISSIZT(IRREP)   * NDIST(IRREP)
      IF(IRREP.EQ.1)THEN
      IOFFY5AA(IRREP)  = 0
      IOFFY5BB(IRREP)  = 0
      IOFFY5AB(IRREP)  = 0
      IOFFT(IRREP)     = 0
      ELSE
      IOFFY5AA(IRREP)  = IOFFY5AA(IRREP-1) + 
     1                   DISSIZAA(IRREP-1)*DISSIZAA(IRREP-1)
      IOFFY5BB(IRREP)  = IOFFY5BB(IRREP-1) + 
     1                   DISSIZBB(IRREP-1)*DISSIZBB(IRREP-1)
      IOFFY5AB(IRREP)  = IOFFY5AB(IRREP-1) + 
     1                   DISSIZAB(IRREP-1)*DISSIZAB(IRREP-1)
      IOFFT(IRREP)     = IOFFT(IRREP-1) +
     1                   DISSIZT(IRREP-1) *NDIST(IRREP-1)
      ENDIF
  620 CONTINUE
C
      I000 = 1
      I010 = I000 + LENAB
      I020 = I010 + LENAB
      I030 = I020 + 2 * (NVRTA * NVRTB + NVRTA * NVRTB)
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10B'
      ENDIF
C
C     Read BcEa intermediates.
C
      DO  630 IRREP=1,NIRREP
      CALL GETLIST(CORE(I010+IOFFY5AB(IRREP)),
     1             1,IRPDPD(IRREP,13),2,IRREP,3)
  630 CONTINUE
C
C     (B,c,E,a) ---> (B,E,c,a)
C
      CALL SSTGEN(CORE(I010),CORE(I000),LENAB,
     1            VRT(1,1),VRT(1,2),VRT(1,1),VRT(1,2),
     1            CORE(I020),1,'1324')
C
C     (B,E,c,a) at I000. Read dacf at I020 and form dfca at I010.
C
      I020 = I010 + LENBB
      I030 = I020 + LENBB
      I040 = I030 + 2 * (NVRTA*NVRTB + NVRTA*NVRTB)
C
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10B'
      ENDIF
C
      DO  640 IRREP=1,NIRREP
      CALL GETLIST(CORE(I020+IOFFY5BB(IRREP)),
     1             1,IRPDPD(IRREP,2),2,IRREP,2+IUHF-1)
C
C     (d<a,c<f) ---> (d,a,c<f)
      CALL SYMEXP2(IRREP,VRT(1,2),DISSIZBB(IRREP),
     1             IRPDPD(IRREP,2),IRPDPD(IRREP,2),
     1             CORE(I020+IOFFY5BB(IRREP)),
     1             CORE(I020+IOFFY5BB(IRREP)))
C     (d,a,c<f) ---> (d,a,c,f)
c      write(6,*) irrep,vrt(1,2),dissizbb(irrep),ioffy5bb(irrep)
c      write(6,*) i020
      CALL  SYMEXP(IRREP,VRT(1,2),DISSIZBB(IRREP),
     1             CORE(I020+IOFFY5BB(IRREP)))
C
  640 CONTINUE
C
C     (d,a,c,f) ---> (d,f,c,a)
C
      CALL SSTGEN(CORE(I020),CORE(I010),LENBB,
     1            VRT(1,2),VRT(1,2),VRT(1,2),VRT(1,2),
     1            CORE(I030),1,'1432')
C
C     Now perform massive matrix multiply :
C
C     (BE,ca) * (df,ca) = (BE,df)
C
C     Y5 intermediates are now at I000 and I010. Target can be at I020.

      I030 = I020 + LENT
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10B'
      ENDIF
C
      DO  650 IRPAC=1,NIRREP
C
      CALL XGEMM('N','T',DISSIZT(IRPAC),NDIST(IRPAC),NSUM(IRPAC),
     1            1.0D+00,
     1           CORE(I000+   IOFFT(IRPAC)),DISSIZT(IRPAC),
     1           CORE(I010+IOFFY5BB(IRPAC)),NDIST(IRPAC),
     1            0.0D+00,
     1           CORE(I020+IOFFT(IRPAC)),DISSIZT(IRPAC))
C
  650 CONTINUE
C
C     We have matrix (BE,df). Form (Bd,Ef).
C
      CALL SSTGEN(CORE(I020),CORE(I000),LENT,
     1            VRT(1,1),VRT(1,1),VRT(1,2),VRT(1,2),
     1            CORE(I030),1,'1324')
C
C     Now loop over symmetry blocks of Y10 and insert appropriate pieces.
C
      DO  790 IRPBD=1,NIRREP
C
      CALL GETLIST(CORE(I010),1,IRPDPD(IRPBD,13),2,IRPBD, 6)
      CALL VADD(CORE(I010),CORE(I010),CORE(I000 + IOFFY5AB(IRPBD)),
     1          IRPDPD(IRPBD,13)*IRPDPD(IRPBD,13), 1.0D+00)
      CALL PUTLIST(CORE(I010),1,IRPDPD(IRPBD,13),2,IRPBD, 6)
C
  790 CONTINUE
C
C     BdEf = -CBAE * AdCf = + BCAE * AdCf
C     (BCAE * AdCf ---> BEAC * ACdf = BEdf ---> BdEf
C
      LENAA = 0
      LENBB = 0
      LENAB = 0
      LENT  = 0
      DO  820 IRREP=1,NIRREP
      DISSIZAA(IRREP) = IRPDPD(IRREP,19)
      DISSIZBB(IRREP) = IRPDPD(IRREP,20)
      DISSIZAB(IRREP) = IRPDPD(IRREP,13)
      DISSIZT(IRREP)  = IRPDPD(IRREP,19)
      NDIST(IRREP)    = IRPDPD(IRREP,20)
      NSUM(IRREP)     = IRPDPD(IRREP,19)
      LENAA = LENAA + DISSIZAA(IRREP)  * DISSIZAA(IRREP)
      LENBB = LENBB + DISSIZBB(IRREP)  * DISSIZBB(IRREP)
      LENAB = LENAB + DISSIZAB(IRREP)  * DISSIZAB(IRREP)
      LENT  = LENT  + DISSIZT(IRREP)   * NDIST(IRREP)
      IF(IRREP.EQ.1)THEN
      IOFFY5AA(IRREP)  = 0
      IOFFY5BB(IRREP)  = 0
      IOFFY5AB(IRREP)  = 0
      IOFFT(IRREP)     = 0
      ELSE
      IOFFY5AA(IRREP)  = IOFFY5AA(IRREP-1) + 
     1                   DISSIZAA(IRREP-1)*DISSIZAA(IRREP-1)
      IOFFY5BB(IRREP)  = IOFFY5BB(IRREP-1) + 
     1                   DISSIZBB(IRREP-1)*DISSIZBB(IRREP-1)
      IOFFY5AB(IRREP)  = IOFFY5AB(IRREP-1) + 
     1                   DISSIZAB(IRREP-1)*DISSIZAB(IRREP-1)
      IOFFT(IRREP)     = IOFFT(IRREP-1) +
     1                   DISSIZT(IRREP-1) *NDIST(IRREP-1)
      ENDIF
  820 CONTINUE
C
      I000 = 1
      I010 = I000 + LENAB
      I020 = I010 + LENAB
      I030 = I020 + 2 * (NVRTA * NVRTB + NVRTA * NVRTB)
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10B'
      ENDIF
C
C     Read AdCf intermediates.
C
      DO  830 IRREP=1,NIRREP
      CALL GETLIST(CORE(I010+IOFFY5AB(IRREP)),
     1             1,IRPDPD(IRREP,13),2,IRREP,3)
  830 CONTINUE
C
C     (A,d,C,f) ---> (A,C,d,f)
C
      CALL SSTGEN(CORE(I010),CORE(I000),LENAB,
     1            VRT(1,1),VRT(1,2),VRT(1,1),VRT(1,2),
     1            CORE(I020),1,'1324')
C
C     (A,C,d,f) at I000. Read BCAE at I020 and form BEAC at I010.
C
      I020 = I010 + LENAA
      I030 = I020 + LENAA
      I040 = I030 + 2 * (NVRTA*NVRTB + NVRTA*NVRTB)
C
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10B'
      ENDIF
C
      DO  840 IRREP=1,NIRREP
      CALL GETLIST(CORE(I020+IOFFY5AA(IRREP)),
     1             1,IRPDPD(IRREP,1),2,IRREP,1)
C
C     (B<C,A<E) ---> (B,C,A<E)
      CALL SYMEXP2(IRREP,VRT(1,1),DISSIZAA(IRREP),
     1             IRPDPD(IRREP,1),IRPDPD(IRREP,1),
     1             CORE(I020+IOFFY5AA(IRREP)),
     1             CORE(I020+IOFFY5AA(IRREP)))
C     (B,C,A<E) ---> (B,C,A,E)
      CALL  SYMEXP(IRREP,VRT(1,1),DISSIZAA(IRREP),
     1             CORE(I020+IOFFY5AA(IRREP)))
C
  840 CONTINUE
C
C     (B,C,A,E) ---> (B,E,A,C)
C
      CALL SSTGEN(CORE(I020),CORE(I010),LENBB,
     1            VRT(1,1),VRT(1,1),VRT(1,1),VRT(1,1),
     1            CORE(I030),1,'1432')
C
C     Now perform massive matrix multiply :
C
C     (BE,AC) * (AC,df) = (BE,df)
C
C     Y5 intermediates are now at I000 and I010. Target can be at I020.

      I030 = I020 + LENT
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10B'
      ENDIF
C
      DO  850 IRPAC=1,NIRREP
C
      CALL XGEMM('N','N',DISSIZT(IRPAC),NDIST(IRPAC),NSUM(IRPAC),
     1            1.0D+00,
     1           CORE(I010+IOFFY5AA(IRPAC)),DISSIZT(IRPAC),
     1           CORE(I000+   IOFFT(IRPAC)),NSUM(IRPAC),
     1            0.0D+00,
     1           CORE(I020+IOFFT(IRPAC)),DISSIZT(IRPAC))
C
  850 CONTINUE
C
C     We have matrix (BE,df). Form (Bd,Ef).
C
      CALL SSTGEN(CORE(I020),CORE(I000),LENT,
     1            VRT(1,1),VRT(1,1),VRT(1,2),VRT(1,2),
     1            CORE(I030),1,'1324')
C
C     Now loop over symmetry blocks of Y10 and insert appropriate pieces.
C
      DO  990 IRPBD=1,NIRREP
C
      CALL GETLIST(CORE(I010),1,IRPDPD(IRPBD,13),2,IRPBD, 6)
      CALL VADD(CORE(I010),CORE(I010),CORE(I000 + IOFFY5AB(IRPBD)),
     1          IRPDPD(IRPBD,13)*IRPDPD(IRPBD,13), 1.0D+00)
      CALL PUTLIST(CORE(I010),1,IRPDPD(IRPBD,13),2,IRPBD, 6)
C
  990 CONTINUE
C
      RETURN
      END
