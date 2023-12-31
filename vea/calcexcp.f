      SUBROUTINE CALCEXCP(KSPIN, IUHF, SCR, MAXCOR, STSPLIT,
     $   DROPOPEN, H0COLUMN, LIST)
C
C AN EXCITATION PATTERN IS FORMED IN THE SAME FORMAT AS THE S-VECTORS ARE
C STORED IN ROUTINE EADAVID (SEE LOADS). THE EXCITATION PATTERN CONTAINS A
C ONE (1.0) FOR THOSE OPERATORS (CONFIGURATIONS) THAT ARE OF INTEREST,
C A ZERO (0.0)  OTHERWHISE. THE PATTERN IS PUT ON COLUMN H0COLUMN OF 
C LISTH0
C
C IF EXCICORE THEN
C
C   THE CONTENTS OF THE EXCITATION PATTERN DEPENDS ON DROPOPEN, STSPLIT, EXCICORE,
C   SINGONLY.
C   IF (EXCICORE = .TRUE.) WE HAVE THE MOST OPTIONS
C
C   DEFINE OCC = OCCO(A) + OCCO(B) + OCCO(I)
C      WHERE
C      OCCO(P) = 0.0 UNLESS
C      P = OCCUPIED CORE-ORBITAL : OCCO(I) = -3.0
C      P = VIRTUAL  CORE-ORBITAL : OCCO(i) =  1.0
C
C      OCC CAN TAKE THE VALUES -3, -2, -1, 0, 1
C         OCC = 0    : NO OPEN SHELL ORBITALS
C         OCC = 1    : VIRT. CORE-ORBITAL
C         OCC = -3   : OCC.  CORE-ORBITAL
C         OCC = -1,-2:  SPIN FLIP OPERATOR      
C
C
C      DROPOPEN  STSPLIT             COLUMNH0(i) = 1.0 if OCC equals
C
C       FALSE     FALSE                 ALL
C       FALSE     TRUE               OCC=0,1,-3 (NOT -1,-2)
C       TRUE      FALSE             OCC=0,-1,-2 (NOT 1,-3)
C       TRUE      TRUE                 OCC=0
C
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION SCR, OCCI, OCCBI, OCC
      DIMENSION SCR(MAXCOR), IOFFVRT(8,2), IOFFPOP(8,2), NDUMS(8)
      LOGICAL LEFTHAND, EXCICORE, SINGONLY, DROPOPEN, STSPLIT, print
C
      COMMON/SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/SINFO/NS(8), SIRREP
      COMMON/EACALC/LEFTHAND, EXCICORE, SINGONLY, DROPCORE
      COMMON/COREINFO/IREPCORE, SPINCORE, IORBCORE, IORBOCC
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
C
      print = .false.      
      NBASA = NOCCO(1) + NVRTO(1)
      IF (IUHF.NE.0) THEN
         NBASB = NOCCO(2)  + NVRTO(2)
      ELSE
         NBASB = 0
      ENDIF
      I000 = 1
      I010 = I000 + NBASA
      I020 = I010 + NBASB
      I030 = I020 + NSIZEC
      I040 = I030 + NSIZEC
C
      IF (.NOT. EXCICORE) THEN
         DO 1 I=I020,I020+NSIZEC -1
            SCR(I) = 1.0D0
 1       CONTINUE
         IF (SINGONLY) THEN
            DO 2 I=I020+VRT(SIRREP,KSPIN), I020+NSIZEC-1
               SCR(I) = 0.0D0
 2          CONTINUE
         ENDIF
         CALL PUTLST(SCR(I020), H0COLUMN,1,1,1,LIST)
      ENDIF
C
      IF (EXCICORE) THEN
C
C WE ASSIGN AN 'OCCO-NUMBER TO EACH ORBITAL'
C 
C DOUBLY OCCUPIED ORBITALS : OCCO = 0
C OCCUPIED CORE ORBITAL    : OCCO = - 3
C VIRTUAL  CORE ORBITAL    : OCCO = 1
C
C NEXT WE EVALUATE 
C    OCC(A) = OCCO(A)
C    OCCO(AB,I) = OCCO(A) + OCCO(B) + OCCO(I)
C    OCCO(aB,i) = OCCO(a) + OCCO(B) + OCCO(i)
C    Where B has kspin (input)
C   
      DO 5 I = 1, NBASA
         SCR(I000+I-1) = 0.0D0
 5    CONTINUE
      IF (IUHF. NE. 0) THEN
         DO 7 I = 1, NBASB
            SCR(I010+I-1) = 0.0D0
 7       CONTINUE
      ENDIF
C
C  CALCULATE OFFSETS IN INTEREST NUMBER ARRAYS
C
      IOFFPOP(1, 1) = I000 - 1
      IOFFPOP(1, 2) = I010 - 1
      DO 10 IRREP = 2, NIRREP
         DO 20 ISPIN = 1, 1 + IUHF
            IOFFPOP(IRREP, ISPIN) = IOFFPOP(IRREP-1, ISPIN) +
     $         POP(IRREP-1, ISPIN)
 20      CONTINUE
 10   CONTINUE
      IOFFVRT(1, 1) = IOFFPOP(NIRREP, 1) + POP(NIRREP,1)
      IF (IUHF .NE. 0) THEN
         IOFFVRT(1, 2) = IOFFPOP(NIRREP, 2) + POP(NIRREP,2)
      ENDIF
      DO 30 IRREP = 2, NIRREP
         DO 40 ISPIN = 1, 1 + IUHF
            IOFFVRT(IRREP, ISPIN) = IOFFVRT(IRREP-1, ISPIN) +
     $         VRT(IRREP-1, ISPIN)
 40      CONTINUE
 30   CONTINUE
C
         SCR(IOFFPOP(IREPCORE,3-SPINCORE) + IORBOCC) = -3.0D0
         SCR(IOFFVRT(IREPCORE,SPINCORE) + IORBCORE)  =  1.0D0
C
C  FILL SINGLES
C
      ICOUNT = I020
      DO 50 A = 1, VRT(SIRREP, KSPIN)
         SCR(ICOUNT) = SCR(IOFFVRT(SIRREP, KSPIN)+A)
         ICOUNT = ICOUNT + 1
 50   CONTINUE
      IF (IUHF .NE. 0) THEN
C
C   FILL AA-PART OF H0
C
         ICOUNT0 = ICOUNT
         MSPIN = KSPIN
         CALL IZERO(NDUMS, 8)
         NDUMS(SIRREP) = 1
         CALL GETLEN2(LENAA, IRPDPD(1,KSPIN),POP(1,KSPIN), NDUMS)
         ISCR = ICOUNT + LENAA
         DO 60 XIRREP = 1, NIRREP
            ICOUNT = ISCR
            MIRREP = DIRPRD(XIRREP, SIRREP)
            DO 70 I = 1, POP(MIRREP, MSPIN)
               OCCI = SCR(IOFFPOP(MIRREP,MSPIN)+I)
               DO 80 BIRREP = 1, NIRREP
                  AIRREP = DIRPRD(BIRREP, XIRREP)
                  DO 90 B=1, VRT(BIRREP, KSPIN)
                     OCCBI = SCR(IOFFVRT(BIRREP, KSPIN)+B) + OCCI
                     DO 100 A = 1, VRT(AIRREP, MSPIN)
                        OCC = SCR(IOFFVRT(AIRREP, MSPIN)+A) + OCCBI
                        SCR(ICOUNT) = OCC
                        ICOUNT = ICOUNT + 1
 100                 CONTINUE
 90               CONTINUE
 80            CONTINUE
 70         CONTINUE
         CALL SQSYM(XIRREP, VRT(1, MSPIN), IRPDPD(XIRREP,MSPIN),
     $         IRPDPD(XIRREP, 18+MSPIN), POP(MIRREP, MSPIN),
     $          SCR(ICOUNT0), SCR(ISCR))
         ICOUNT0 = ICOUNT0 + IRPDPD(XIRREP, KSPIN) * POP(MIRREP,MSPIN)
 60      CONTINUE
         ICOUNT = ICOUNT0
         IF (ICOUNT. NE. ISCR) THEN
            WRITE(6,*)' SOMETHING WRONG IN CALCEXCP', ICOUNT, ISCR
         ENDIF
      ENDIF
C
C  FILL AB PART OF H0
C
      MSPIN = 2+IUHF-KSPIN
      DO 160 XIRREP = 1, NIRREP
         MIRREP = DIRPRD(XIRREP, SIRREP)
         DO 170 I = 1, POP(MIRREP, MSPIN)
            OCCI = SCR(IOFFPOP(MIRREP,MSPIN)+I)
            DO 180 BIRREP = 1, NIRREP
               AIRREP = DIRPRD(BIRREP, XIRREP)
               DO 190 B=1, VRT(BIRREP, KSPIN)
                  OCCBI = SCR(IOFFVRT(BIRREP, KSPIN)+B) + OCCI
                  DO 200 A = 1, VRT(AIRREP, MSPIN)
                     OCC = SCR(IOFFVRT(AIRREP, MSPIN)+A) + OCCBI
                     SCR(ICOUNT) = OCC
                     ICOUNT = ICOUNT + 1
 200              CONTINUE
 190           CONTINUE
 180        CONTINUE
 170     CONTINUE
 160  CONTINUE
C
      IF (PRINT) THEN
         WRITE(6,*) ' OCCUPATION NUMBER ARRAY IN CALCEXCP'
         CALL OUTPUT(SCR(I020),1,1,1,NSIZEC,1,NSIZEC,1)
      ENDIF
C
      IF (DROPOPEN) THEN
         IF (STSPLIT) THEN
            WRITE(6,*) ' ALL OPEN SHELL OPERATORS ARE EXCLUDED'
         ELSE
            WRITE(6,*) ' ALL OPEN SHELL, EXCEPT PURE-SPIN-FLIP,',
     $         ' OPERATORS ARE EXCLUDED'
         ENDIF
      ELSE
         IF (STSPLIT) THEN
            WRITE(6,*) ' SPIN FLIP OPERATORS ARE EXCLUDED'
         ELSE
            WRITE(6,*) ' ALL OPEN SHELL OPERATORS ARE INCLUDED'
         ENDIF
      ENDIF
C
C      DROPOPEN  STSPLIT             EXCP(i) = 1.0 if OCC equals, 0 OTHERWHISE
C
C       FALSE     FALSE                 ALL
C       FALSE     TRUE               OCC=0,1,-3 (NOT -1,-2)
C       TRUE      FALSE             OCC=0,-1,-2 (NOT 1,-3)
C       TRUE      TRUE                 OCC=0
C
      IF (DROPOPEN) THEN
         IF (STSPLIT) THEN
C
C       TRUE      TRUE          OCC=0 
C
            DO 320 I = 1, NSIZEC
               IF (ABS(SCR(I020 + I -1)) .LT. 0.1) THEN
                  SCR(I030+I-1) = 1.0D0
               ELSE
                  SCR(I030+I-1) = 0.0D0
               ENDIF
 320        CONTINUE
         ELSE
C
C       TRUE      FALSE         OCC=0,-1,-2 (NOT 1,-3)
C
            DO 330 I = 1, NSIZEC
               IF ((SCR(I020 + I -1) .LT. -2.5) .OR.
     $            (SCR(I020+I-1) .GT. 0.5)) THEN
                  SCR(I030+I-1) = 0.0D0
               ELSE
                  SCR(I030+I-1) = 1.0D0
               ENDIF
 330        CONTINUE
         ENDIF
      ELSE
         IF (STSPLIT) THEN
C
C       FALSE     TRUE          OCC=0,1,-3 (NOT -1,-2)
C
            DO 340 I = 1, NSIZEC
               IF ((SCR(I020 + I -1) .GT. -2.5) .AND.
     $            (SCR(I020+I-1) .LT. -0.5)) THEN
                  SCR(I030+I-1) = 0.0D0
               ELSE
                  SCR(I030+I-1) = 1.0D0
               ENDIF
 340        CONTINUE
         ELSE
C
C       FALSE     FALSE             ALL
C
            DO 350 I = 1, NSIZEC
               SCR(I030+I-1) = 1.0D0
 350        CONTINUE
         ENDIF
      ENDIF
C
      IF (PRINT) THEN
         WRITE(6,*) ' CALCEXCP, COLUMN ', H0COLUMN
         CALL OUTPUT(SCR(I030),1,1,1,NSIZEC,1,NSIZEC,1)
      ENDIF
C
      CALL PUTLST(SCR(I030),H0COLUMN,1,1,1,LIST)
C
      ENDIF
C
      RETURN
      END
