      SUBROUTINE CALCH0(KSPIN, IUHF, SCR, MAXCOR)
C
C THE DENOMINATOR IS FORMED IN THE SAME FORMAT AS THE S-VECTORS ARE
C STORED IN ROUTINE EADAVID (SEE LOADS)
C
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION SCR, EI, EBI, EA, FACTOR
      DIMENSION SCR(MAXCOR), IOFFVRT(8,2), IOFFPOP(8,2), NDUMS(8)
      logical print, EXCICORE
C
      COMMON/SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/LISTDAV/LISTC, LISTHC, LISTH0
      COMMON/SINFO/NS(8), SIRREP
      COMMON/EACALC/LEFTHAND, EXCICORE, SINGONLY, DROPCORE
      COMMON/COREINFO/IREPCORE, SPINCORE, IORBCORE, IORBOCC
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
C
      print =.false.
      NBASA = NOCCO(1) + NVRTO(1)
      IF (IUHF.NE.0) THEN
         NBASB = NOCCO(2)  + NVRTO(2)
      ELSE
         NBASB = 0
      ENDIF
      I000 = 1
      I010 = I000 + NBASA
      I020 = I010 + NBASB
      CALL GETREC(20, 'JOBARC', 'SCFEVALA', IINTFP*NBASA, SCR(I000))
      if (print) then
         write(6,*) ' Hartree Fock orbital energies : alfa'
         call output(scr(i000), 1, 1, 1, nbasa, 1, nbasa, 1)
      endif
      IF (IUHF. NE. 0) THEN
         CALL GETREC(20, 'JOBARC', 'SCFEVALB', IINTFP*NBASB, SCR(I010))
      if (print) then
         write(6,*) ' Hartree Fock orbital energies : beta'
         call output(scr(i010), 1, 1, 1, nbasb, 1, nbasb, 1)
      endif
      ENDIF
C
C  CALCULATE OFFSETS IN ORBITAL ENERGY ARRAYS
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
      if (print) then
         write(6,*) 'occupied offsets in calch0 '
         do 12 irrep = 1, 1, nirrep
            write(6,*) (ioffpop(irrep, ispin), pop(irrep, ispin),
     $         ispin = 1, 1+iuhf)
 12      continue
         write(6,*) 'virtual offsets in calch0 '
         do 14 irrep = 1, 1, nirrep
            write(6,*) (ioffvrt(irrep, ispin), vrt(irrep, ispin),
     $         ispin = 1, 1+iuhf)
 14      continue
      endif
C
      ICOUNT = I020 
      DO 50 A = 1, VRT(SIRREP, KSPIN)
         SCR(ICOUNT) = SCR(IOFFVRT(SIRREP, KSPIN) + A)
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
               EI = SCR(IOFFPOP(MIRREP,MSPIN)+I)
               DO 80 BIRREP = 1, NIRREP
                  AIRREP = DIRPRD(BIRREP, XIRREP)
                  DO 90 B=1, VRT(BIRREP, KSPIN)
                     EBI = SCR(IOFFVRT(BIRREP, KSPIN)+B) - EI
                     DO 100 A = 1, VRT(AIRREP, MSPIN)
                        EA = SCR(IOFFVRT(AIRREP, MSPIN)+A)
                        SCR(ICOUNT) = EA + EBI
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
            WRITE(6,*)' SOMETHING WRONG IN CALCH0', ICOUNT, ISCR
         ENDIF
      ENDIF
C
C  FILL AB PART OF H0
C
      MSPIN = 2+IUHF-KSPIN
      DO 160 XIRREP = 1, NIRREP
         MIRREP = DIRPRD(XIRREP, SIRREP)
         DO 170 I = 1, POP(MIRREP, MSPIN)
            EI = SCR(IOFFPOP(MIRREP,MSPIN)+I)
            DO 180 BIRREP = 1, NIRREP
               AIRREP = DIRPRD(BIRREP, XIRREP)
               DO 190 B=1, VRT(BIRREP, KSPIN)
                  EBI = SCR(IOFFVRT(BIRREP, KSPIN)+B) - EI
                  DO 200 A = 1, VRT(AIRREP, MSPIN)
                     EA = SCR(IOFFVRT(AIRREP, MSPIN)+A)
                     SCR(ICOUNT) = EA + EBI
                     ICOUNT = ICOUNT + 1
 200              CONTINUE
 190           CONTINUE
 180        CONTINUE
 170     CONTINUE
 160  CONTINUE
      nea = icount - i020
      print = nea .le. 10
      if (print) then
         write(6,*)' zeroth order eigenvalues'
         call output(scr(i020), 1, 1, 1, nea, 1, nea, 1)
      endif
C
C  WRITE H0 TO LISTH0 = 472
C
      CALL PUTLST(SCR(I020), 1,1,1,1,LISTH0)
C
C MAKE GUESSLIST OF INTERESTING EXCITATIONS, BASED ON ZEROTH ORDER ENERGIES.
C (ON COLUMN 4 OF LIST H0)
C
      I000 = 1
      I010 = I000 + NSIZEC
      I020 = I010 + NSIZEC
      I030 = I020 + NSIZEC
C
      CALL GETLST(SCR(I000),1,1,1,1,LISTH0)
      CALL GETLST(SCR(I010),5,1,1,1,LISTH0)
      DO 210 I=I010, I010+NSIZEC-1
         SCR(I) = 1.0D0-SCR(I)
 210  CONTINUE
C
      FACTOR = 1.0D+00
C
      CALL VADD(SCR(I020),SCR(I000),SCR(I010),NSIZEC,FACTOR)
      CALL PUTLST(SCR(I020),4,1,1,1,LISTH0)
C
      RETURN
      END
