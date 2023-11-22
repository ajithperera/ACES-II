      SUBROUTINE INITNUMR_IP(SCR, MAXCOR, IUHF)
C
C IF THE NUMBER OF ROOTS TO BE DETERMINED IS NOT GIVEN (I.E. ALL ARE ZERO)
C THEN THE NUMBER OF ROOTS TO BE DETERMINED IN THE IP EOM CALCULATION IS 
C DETERMINED FROM INSPECTION OF THE SCF ENERGY EIGENVALUES. POSSIBLE 
C DEGENERACIES BETWEEN BLOCKS IS ACCOUNTED FOR.
C
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION SCR, IP(8,2), IPMIN, DIFFE
      DIMENSION SCR(MAXCOR),IOFFVRT(8,2), IOFFPOP(8,2)
C
#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"
#include "info.com"

      COMMON/IPCALC/LEFTHAND,SINGONLY, DROPCORE
      COMMON/IpINFO/NUMROOT(8,3)
C
      IROOT = 0
      DO IRREP = 1, NIRREP
         DO ISPIN = 1, 1+IUHF
            IROOT = IROOT + NUMROOT(IRREP, ISPIN)
         ENDDO
      ENDDO
      IF (IROOT .GT. 0) THEN
         RETURN
      ELSE
         WRITE(6,*)
         WRITE(6,"(a,a)") ' NUMBER OF DESIRED ROOTS',
     $                    ' IS ESTIMATED FROM ORBITAL EIGENVALUES'
         WRITE(6,*)
      ENDIF
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
         write(6,*) ' Hartree Fock orbital energies : alfa'
         call output(scr(i000), 1, 1, 1, nbasa, 1, nbasa, 1)
#endif 
      IF (IUHF. NE. 0) THEN
         CALL GETREC(20, 'JOBARC', 'SCFEVALB', IINTFP*NBASB, SCR(I010))

#ifdef _DEBUG_LVL0
         write(6,*) ' Hartree Fock orbital energies : beta'
         call output(scr(i010), 1, 1, 1, nbasb, 1, nbasb, 1)
#endif 
      ENDIF
C
C  CALCULATE OFFSETS IN ORBITAL ENERGY ARRAYS
C
      IOFFPOP(1,1) = I000
      IOFFPOP(1,2) = I010

      DO IRREP = 2, NIRREP
         DO ISPIN = 1, 1 + IUHF
            IOFFPOP(IRREP,ISPIN) = IOFFPOP(IRREP-1,ISPIN) +
     $                             POP(IRREP-1, ISPIN)
         ENDDO 
      ENDDO 

      IOFFVRT(1,1) = IOFFPOP(NIRREP,1) + POP(NIRREP,1)

      IF (IUHF .NE. 0) THEN
         IOFFVRT(1,2) = IOFFPOP(NIRREP,2) + POP(NIRREP,2)
      ENDIF

      DO IRREP = 2, NIRREP
         DO ISPIN = 1, 1 + IUHF
            IOFFVRT(IRREP,ISPIN) = IOFFVRT(IRREP-1, ISPIN) +
     $                             VRT(IRREP-1, ISPIN)
         ENDDO 
      ENDDO 
C
C DETERMINE MINIMUM SCF IP PER IRREP
C
      DO IRREP = 1, NIRREP
         DO ISPIN = 1, 1 + IUHF
            IF (POP(IRREP,ISPIN) .NE. 0) THEN
            IP(IRREP,ISPIN) = DABS(SCR(IOFFPOP(IRREP,ISPIN)+
     &                                         POP(IRREP,ISPIN)-1))
            ENDIF 
         ENDDO 
      ENDDO 

#ifdef _DEBUG_LVL0
         write(6,*)
         WRITE(6,*) '  Minumm SCF Ionization potentials'
         CALL OUTPUT(IP,1,nirrep,1,2,8,2,1)
#endif 
C
C DETERMINE minimum Ionization Potential.
C
      IPMIN = IP(1,1)
      DO  ISPIN = 1, 1 + IUHF
         DO  IRREP = 1, NIRREP
            IF (POP(IRREP,ISPIN).GT.0) IPMIN=MIN(IPMIN,IP(IRREP,ISPIN))
         ENDDO 
      ENDDO
     
#ifdef _DEBUG_LVL0
      write(6,*)
      Write(6,"(a,F10.6)") " IPMIN =",IPMIN
#endif 
C      
      DO  ISPIN =1, 1+IUHF
         DO  IRREP = 1, NIRREP
             NUMROOT(IRREP,ISPIN) = 0
             IF (POP(IRREP,ISPIN).GT.0) THEN
                 IF (ABS(IPMIN-IP(IRREP,ISPIN)).LT.1.0D0) THEN
                     NUMROOT(IRREP,ISPIN) = 1
                 ENDIF
            ENDIF
         ENDDO
      ENDDO
C
C CHECK FOR DEGENERACIES
C
      DO ISPIN = 1, 1+IUHF
         DO IRREP = 1, NIRREP - 1
            IF (NUMROOT(IRREP,ISPIN) .NE. 0) THEN
               DO JRREP = IRREP + 1, NIRREP
                  IF ((NUMROOT(JRREP,ISPIN) .NE. 0) .AND.
     +               (POP(IRREP,ISPIN) .EQ. POP(JRREP,ISPIN)) .AND.
     +               (VRT(IRREP,ISPIN) .EQ. VRT(JRREP,ISPIN))) THEN
                     DIFFE = 0.0D0
                     DO I = 0, POP(IRREP,ISPIN) - 1
                     DIFFE = DIFFE + ABS(SCR(IOFFPOP(IRREP,ISPIN) + I) -
     +                                   SCR(IOFFPOP(JRREP,ISPIN) + I))
                     ENDDO
                     DO I = 0, VRT(IRREP,ISPIN) - 1
                     DIFFE = DIFFE + ABS(SCR(IOFFVRT(IRREP,ISPIN) + I) -
     $                               SCR(IOFFVRT(JRREP,ISPIN) + I))
                     ENDDO
                     IF (ABS(DIFFE) .LT. 0.0001)  THEN
                        NUMROOT(JRREP,ISPIN) = 0
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
C
      RETURN
      END
