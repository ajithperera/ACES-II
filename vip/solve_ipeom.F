#include "flags.h"
C
C ***************************************************************
C FOLLOWING ARE A LARGE NUMBER OF SUBROUTINES THAT RELATE TO THE DIRECT 
C DIAGONALIZATION PROCEDURE, BASED ON DAVIDSON'S ALGORITHM
C ***************************************************************
C 
      SUBROUTINE SOLVE_IPEOM(ICORE, MAXCOR, IUHF)
C
C THIS ROUTINE DRIVES THE CALCULATION OF THE IOINIZATION POTENTIALS 
C BY DIAGONALIZATION OF THE IP-EOM MATRIX BY THE DAVIDSON PROCEDURE
C
      IMPLICIT INTEGER (A-Z)
      LOGICAL CONVRG, EXCICORE, LEFTHAND, SINGONLY, DROPCORE
      DOUBLE PRECISION EIGVAL, OSCSTR
      CHARACTER *5 SPIN(2)
C
      DIMENSION ICORE(MAXCOR)
C
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "sympop.com"
#include "flags.com"
#include "flags2.com"

      COMMON /EXTINF/NDIMR,IOLDEST
      COMMON /EXTINF2/ROOT
      COMMON /EXTINF3/IROOT,LOCROOT,ITROOT
      COMMON /EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON /IPINFO/NUMROOT(8,3)
      COMMON /IPCALC/LEFTHAND,SINGONLY,DROPCORE
      COMMON /SINFO/ NS(8), SIRREP
      COMMON /ROOTS/EIGVAL(100,8,3),OSCSTR(100,8,3)
      COMMON /STINFO/ITOTALS, STMODE, LISTST, STCALC, NSIZEST
C
      DATA SPIN /'ALPHA', 'BETA '/
      PARAMETER (MAXORD = 100)
C
      DMAXCOR = MAXCOR/IINTFP
C
      DO 100 IRREP = 1, NIRREP
         DO 50 ISPIN = 1, 1 + IUHF

         ITOTALS = 1
         IF(NUMROOT(IRREP, ISPIN) .GT. 0) THEN

C
C FIRST FIND THE DIMENSION OF THE VECTORS
C
         CALL IZERO(NS,8)
         NS(IRREP) = 1

         SIRREP = IRREP

         write(6,"(a)") '*********************************************'
         WRITE (6,"(a,i2)") '  SYMMETRY-BLOCK: ', SIRREP
         IF (IUHF .EQ. 0) THEN
            WRITE(6,"(2a)") '  SPIN STATE : DOUBLET ', SPIN(ISPIN)
         ELSE
            WRITE(6,"(2a)") '  SPIN STATE : UNKNOWN ', SPIN(ISPIN)
         ENDIF

         CALL GETNIP(NIP, SIRREP, ISPIN, IUHF)
         CALL GETNIP(NSIZEST, SIRREP, 3-ISPIN, IUHF)

         write(6,1001) NIP
 1001    format('  DIMENSION OF CURRENT IP VECTOR : ', i10)
         WRITE(6,1002) NUMROOT(IRREP,ISPIN),IRREP
 1002    FORMAT(I3, ' ROOTS WILL BE SEARCHED IN IRREP',1X,I2) 
         write(6,"(a)") '*********************************************'
         WRITE(6,*)
C
         IF (NIP.EQ.0) GO TO 50
C
C CONSTRUCT NEW LIST POSITIONS        
C
         ITYP = 2
         CALL IP_NEWLIST(SIRREP, 1, IUHF, ISPIN, ITYP)
         CALL INIT_IPEOMDAV(NIP,SIRREP,ISPIN,IUHF,ICORE,DMAXCOR)

         ISIDE = 1 
         IF (LEFTHAND) ISIDE   = 2 
         TOTITER = 0
         ITER    = 1
         NROOT   = MIN(NUMROOT(SIRREP,ISPIN), NIP)
         ITMAX   = IFLAGS2(h_IFLAGS2_eom_maxcyc)

         IF (ISIDE .EQ. 1) WRITE(6,"(2A)") " THE RIGHT EOM-IP-CC",
     &                                     " EQUATIONS ARE SOLVED"
         IF (ISIDE .EQ. 2) WRITE(6,"(2A)") " THE LEFT  EOM-IP-CC",
     &                                     " EQUATIONS ARE SOLVED"
         WRITE(6,*)
C
 1       IF ((IROOT.LT. NROOT) .AND. (NDIMR .LE. NIP)) THEN
C
C  FORM C_new using Hbar(IP_EOM)C_old = C_new 
C
            CALL HC_MULT(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)
C
C CALL DAVIDSON EXTRAPOLATOR
C            
            I000 = 1
            I010 = I000 + IINTFP*(MAXORD * (2*MAXORD + 3) + 5)
            I020 = I010 + IINTFP * MAXORD * MAXORD

            CALL IPDAVID(ICORE(I000),ICORE(I010),ICORE(I020),
     $                  (MAXCOR - I020 + 1)/IINTFP,CONVRG,NROOT,
     $                  ISIDE,ISPIN,IUHF,ITER,ITMAX)
C
            IF (CONVRG) THEN
               TOTITER = TOTITER + ITER
               ITER = 0
            ENDIF

            ITER = ITER + 1
            NDIMR = NDIMR + 1
            IF (ITER .NE. ITMAX) THEN
               GO TO 1
            ELSE
                WRITE(6,"(a,a)")'TOO MANY ITERATIONS:',
     +                          ' GIVING UP ON THIS SYMMETRY'
            ENDIF
C
         ENDIF
C
         NUMROOT(SIRREP, ISPIN) = IROOT
C
         ENDIF
C
 50   CONTINUE
 100  CONTINUE
C
      CALL TABLEIP(IUHF)
      RETURN
 1000 FORMAT(/)
      END
