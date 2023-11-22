      SUBROUTINE CALCH0(KSPIN, IUHF, SCR, MAXCOR)
C
C THE DENOMINATOR IS FORMED IN THE SAME FORMAT AS THE S-VECTORS ARE
C STORED IN ROUTINE EADAVID (SEE LOADS)
C
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION SCR,EI,EBI,EA,FACTOR
      DIMENSION SCR(MAXCOR),IOFFVRT(8,2),IOFFPOP(8,2),NDUMS(8)
      LOGICAL PRINT
C
#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "info.com"
#
      COMMON/LISTDAV/LISTC, LISTHC, LISTH0
      COMMON/SINFO/NS(8), SIRREP
      COMMON/IPCALC/LEFTHAND,SINGONLY,DROPCORE
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
C
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
#ifdef _DEBUG_LVL0
         write(6,*) ' Hartree Fock orbital energies : Alfa'
         call output(scr(i000), 1, 1, 1, nbasa, 1, nbasa, 1)
#endif 
      IF (IUHF. NE. 0) THEN
         CALL GETREC(20, 'JOBARC', 'SCFEVALB', IINTFP*NBASB, SCR(I010))
#ifdef _DEBUG_LVL0
         write(6,*) ' Hartree Fock orbital energies : Beta'
         call output(scr(i010), 1, 1, 1, nbasb, 1, nbasb, 1)
#endif 
      ENDIF
C
C CALCULATE OFFSETS IN ORBITAL ENERGY ARRAYS
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

#ifdef _DEBUG_LVL0
         Write(6,*)
         write(6,"(a)") 'occupied offsets in calch0 '
         do 12 irrep = 1, nirrep
            write(6,*) (ioffpop(irrep, ispin), pop(irrep, ispin),
     $         ispin = 1, 1+iuhf)
 12      continue
         Write(6,*)
         write(6,"(a)") 'virtual offsets in calch0 '
         do 14 irrep = 1, nirrep
            write(6,*) (ioffvrt(irrep, ispin), vrt(irrep, ispin),
     $         ispin = 1, 1+iuhf)
 14      continue
#endif 
C
      ICOUNT = I020 
      DO 50 A = 1, POP(SIRREP, KSPIN)
         SCR(ICOUNT) = -SCR(IOFFPOP(SIRREP,KSPIN) + A)
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
         CALL GETLEN_HHP2(LENAA,IRPDPD(1,KSPIN+2),VRT(1,KSPIN),NDUMS)
         ISCR = ICOUNT + LENAA
         DO 60 XIRREP = 1, NIRREP
            ICOUNT = ISCR
            AIRREP = DIRPRD(XIRREP, SIRREP)
            DO 70 A = 1, VRT(AIRREP, MSPIN)
               EA = SCR(IOFFVRT(AIRREP,MSPIN)+A)
               DO 80 JIRREP = 1, NIRREP
                  IIRREP = DIRPRD(JIRREP, XIRREP)
                  DO 90 J=1, POP(JIRREP, MSPIN)
                     EBJ = -SCR(IOFFPOP(JIRREP, KSPIN)+J) + EA
                     DO 100 I = 1, POP(IIRREP, KSPIN)
                        EAI = -SCR(IOFFPOP(IIRREP, KSPIN)+I)
                        SCR(ICOUNT) = EAI + EBJ
                        ICOUNT = ICOUNT + 1
 100                 CONTINUE
 90               CONTINUE
 80            CONTINUE
 70         CONTINUE

          CALL SQSYM(XIRREP,POP(1,MSPIN),IRPDPD(XIRREP,MSPIN+2),
     $              IRPDPD(XIRREP,20+MSPIN),VRT(MIRREP,MSPIN),
     $              SCR(ICOUNT0),SCR(ISCR))
         ICOUNT0 = ICOUNT0 + IRPDPD(XIRREP,KSPIN+2) * VRT(MIRREP,MSPIN)
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
         DO 170 I = 1, VRT(MIRREP, MSPIN)
            EA = SCR(IOFFVRT(MIRREP,MSPIN)+I)
            DO 180 BIRREP = 1, NIRREP
               AIRREP = DIRPRD(BIRREP, XIRREP)
               DO 190 B=1, POP(BIRREP, KSPIN)
                  EBI = -SCR(IOFFPOP(BIRREP, KSPIN)+B) + EA
                  DO 200 A = 1, POP(AIRREP, MSPIN)
                     EAJ = -SCR(IOFFPOP(AIRREP, MSPIN)+A)
                     SCR(ICOUNT) = EAJ + EBI
                     ICOUNT = ICOUNT + 1
 200              CONTINUE
 190           CONTINUE
 180        CONTINUE
 170     CONTINUE
 160  CONTINUE

#ifdef _DEBUG_LVLM
      nip = icount - i020
      print = nip .le. 10
      if (print) then
         write(6,*)' zeroth order eigenvalues'
         call output(scr(i020), 1, 1, 1, nip, 1, nip, 1)
      endif
#endif 
C
C WRITE H0 TO LISTH0 = 472
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
      CALL VADD(SCR(I020),SCR(I000),SCR(I010),NSIZEC,FACTOR)
      CALL PUTLST(SCR(I000),4,1,1,1,LISTH0)
C
      RETURN
      END
