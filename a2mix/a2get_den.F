
      SUBROUTINE A2GET_DEN(WORK, MAXCOR, NMBR_OF_PERTS, NAOBFNS,
     &                     IBEGIN_P_DENS, IBEGIN_P_OPRT, 
     &                     IBEGIN_AO_OVRLP, IBEGIN_MO_OVRLP,
     &                     IBEGIN_MO_VECTS, IBEGIN_MEM_LEFT)
C     
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL SPHERICAL, SPIN_D
C
      CHARACTER*7 DENSITY_TYPE
      CHARACTER*8 LABELPSCF, LABELDSCF, LABELDENS, LABELSDEN,
     &            RECNAMEA, RECNAMEB
      DIMENSION WORK(MAXCOR), PRDUTINT(NAOBFNS, NAOBFNS)
C
      DATA IZR0, IONE /0, 1/
C
#include "machsp.com"
#include "syminf.com"
#include "jodaflags.com"
C
      SPHERICAL = (IFLAGS(62) .EQ. 1)
      NAOBFNS2  = NAOBFNS*NAOBFNS
C
      IRWND         = IZR0
      IERR          = IZR0
      IOFF          = INOE
      IBEGIN_P_DENS = IONE
      ISCF_TDEN     = IONE
      ISCF_DDEN     = IONE
      ICOR_TDEN     = IONE
      ICOR_DDEN     = IONE
C
      MAX_MEM_NEED = 5*NAOBFNS2 + 2*NAOBFNS2*NMBR_OF_PERTS
      IF (MAXCOR .LE. MAX_MEM_NEED) CALL INSMEM("A2GET_DEN", 
     &                                                MAX_MEM_NEED, 
     &                                                MAXCOR)
C
      IKEEP_ONE_DENA = IONE 
      IKEEP_ONE_DENB = IKEEP_ONE_DENA + NAOBFNS2*NMBR_OF_PERTS
      IKEEP_MOVECS_A = IKEEP_ONE_DENB + NAOBFNS2
      IKEEP_MOVECS_B = IKEEP_MOVECS_A + NAOBFNS2
      IKEEP_MO_OVRLP = IKEEP_MOVECS_B + NAOBFNS2
      IKEEP_P_TDEN_N = IKEEP_MO_OVRLP + NAOBFNS2

      IBEGIN_P_DENS  = IKEEP_P_TDEN_N  

      DO IPERT = 1, NMBR_OF_PERTS
C
         IKEEP_P_DDEN_N =  IKEEP_P_TDEN_N + NAOBFNS2
         INEXT          =  IKEEP_P_DDEN_N
C
#ifdef _DEBUG_LVLM1
      Print*, "Offsets for 1-order density" 
      Write(*,'(a,1x,i5)')  "The number of pert.", NMBR_OF_PERTS
      Write(*,'(6(1x,i5))') IKEEP_ONE_DENA, IKEEP_ONE_DENB, 
     &                      IKEEP_MOVECS_A,IKEEP_MOVECS_B,
     &                      IKEEP_MO_OVRLP, IKEEP_P_TDEN_N
      Write(6,*) 
#endif
        WRITE(RECNAMEA, "(A,I2)") "PTDENA", IPERT
        CALL GETREC(-20, "JOBARC", RECNAMEA, NAOBFNS2*IINTFP, 
     &               WORK(IKEEP_ONE_DENA))
C
        WRITE(RECNAMEB, "(A,I2)") "PTDENB", IPERT
        CALL GETREC(-20, "JOBARC", RECNAMEB, NAOBFNS2*IINTFP,
     &               WORK(IKEEP_ONE_DENB))
C
#ifdef _DEBUG_LVLM1
      Write(6,'(1x,a,a,1x,a)') "The record labels ", RECNAMEA, 
     &      RECNAMEB
      Write(6,*)
      Write(6,*) "The  1-st order alpha density matrices in AO basis"
      CALL OUTPUT(WORK(IKEEP_ONE_DENA), 1, NAOBFNS, 1, NAOBFNS,
     &            NAOBFNS, NAOBFNS, 1)
C   
      Write(6,*)
      Write(6,*) "The 1-st order beta density matrices in AO basis"
      CALL OUTPUT(WORK(IKEEP_ONE_DENB), 1, NAOBFNS, 1, NAOBFNS,
     &            NAOBFNS, NAOBFNS, 1)
           
#endif   
        CALL ZERO(WORK(IKEEP_P_TDEN_N), NAOBFNS2)
C
C Build perturbed total density for n perturbation 
C
        CALL SAXPY(NAOBFNS2, 1.0D0, WORK(IKEEP_ONE_DENA), 1, 
     &             WORK(IKEEP_P_TDEN_N), 1) 
        CALL SAXPY(NAOBFNS2, -1.0D0, WORK(IKEEP_ONE_DENB), 1, 
     &             WORK(IKEEP_P_TDEN_N), 1) 
C
#ifdef _DEBUG_LVLM1
      Write(6,*) IKEEP_P_TDEN_N
      Write(6,*) "The  1-st order tot. density matrices in AO basis"
CSSS        CALL GETREC(20, 'JOBARC', "HFDENSTY", NAOBFNS2*IINTFP, 
CSSS     &                WORK(IKEEP_P_TDEN_N))
      CALL OUTPUT(WORK(IKEEP_P_TDEN_N), 1, NAOBFNS, 1, NAOBFNS, 
     &            NAOBFNS, NAOBFNS, 1)
C    
#endif
        IKEEP_P_TDEN_N =  INEXT
C
      ENDDO
C
      IOFF          = IKEEP_ONE_DENA
      IBEGIN_P_OPRT = IOFF
C
      OPEN (UNIT=30, FILE='VPOUT', FORM='UNFORMATTED',
     &      STATUS='OLD')
      
      DO JPERT = 1, NMBR_OF_PERTS

         CALL SEEKLB ('   DEN  ', IERR, IRWND, 30)
         IF (IERR .NE. 0) CALL ERREX
C
         CALL LOAD_DEFINED(WORK(IOFF), NAOBFNS*NAOBFNS, NAOBFNS)

#ifdef _DEBUG_LVLM1
         CALL HEADER ('Matrix elements of the defined operator', 1,
     &                6)
         CALL OUTPUT(WORK(IOFF), 1, NAOBFNS, 1, NAOBFNS, NAOBFNS,
     &               NAOBFNS, 1)
#endif
C
         IOFF  = IOFF + NAOBFNS*NAOBFNS
         IRWND = 1
      ENDDO
      CLOSE(30)
C
C Read the atomic overlap integral; Unfortuantely none of these 
C are really needed (these were done when it was not clear to
C what was going on)
C
      IBEGIN_AO_OVRLP = IKEEP_ONE_DENB
C
      CALL GETREC(20,'JOBARC','AOOVRLAP',NAOBFNS*NAOBFNS*IINTFP,
     &            WORK(IBEGIN_AO_OVRLP))

#ifdef _DEBUG_LVLM1
         CALL HEADER ('The Atomic Overlap Matrix', 1,
     &                6)
         CALL OUTPUT(WORK(IBEGIN_AO_OVRLP), 1, NAOBFNS, 1, NAOBFNS,
     &               NAOBFNS, NAOBFNS, 1)
#endif
      IBEGIN_MO_OVRLP = IKEEP_MO_OVRLP
      IBEGIN_MO_VECTS = IKEEP_MOVECS_A
      IBEGIN_MEM_LEFT = MAX_MEM_NEED
C                                                                               C Get the eigenvectors from the JOBARC file.
C 
         CALL GETREC(20, 'JOBARC', 'SCFEVCA0',  NAOBFNS2*IINTFP, 
     &               WORK(IKEEP_MOVECS_A))
         CALL GETREC(20, 'JOBARC', 'SCFEVCB0',  NAOBFNS2*IINTFP, 
     &               WORk(IKEEP_MOVECS_B))
#ifdef _DEBUG_LVLM1
         CALL HEADER ('The MO Overlap Matrix', 1, 6)
         CALL OUTPUT(WORK(IKEEP_MO_OVRLP), 1, NAOBFNS, 1, NAOBFNS,
     &               NAOBFNS, NAOBFNS, 1)
#endif
C
C At this point we need to transform the density matrices to Cartesian
C basis in order to directly contract with the Cartesian AO integrals.
C 
C$$$      IF (SPHERICAL) THEN
C$$$         CALL GETREC(20, 'JOBARC','CMP2CART', NBAS*NBASP*IINTFP, 
C$$$     &               SPH2CART)

C$$$      Write(6,*) "The spherical to Cartesian transformation matrix"
C$$$      CALL OUTPUT(SPH2CART, 1, NBAS, 1, NBASP, NBAS, NBASP, 1)
C
C$$$         CALL TRANS_SPH2CART(SCFDEN, SPH2CART, TMP, NBAS, NBASP)
C$$$         CALL TRANS_SPH2CART(RELDEN, SPH2CART, TMP, NBAS, NBASP)
C$$$      ENDIF
C
C
      RETURN
      END


