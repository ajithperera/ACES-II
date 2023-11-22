#include "flags.h"
      SUBROUTINE Template_orbs(ONEH,DENS,FOCK,EVAL,EVEC,LDIM1,
     &                         LDIM2,NBAS,IUHF,IOS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "maxbasfn.par"
#include "machsp.com"
#include "flags.com"
#include "flags2.com"
#include "symm2.com"
C
      DIMENSION ONEH(LDIM1),DENS((IUHF+1)*LDIM1),FOCK((IUHF+1)*LDIM1)
      DIMENSION EVEC((IUHF+1)*LDIM2),EVAL((IUHF+1)*NBAS)
C
      DIMENSION IDUMMY(MAXBASFN),IDUMMY2(MAXBASFN),ILOCATE(MAXBASFN)
      DIMENSION IDUMMY3(MAXBASFN), NSUM(8, 2)
      DIMENSION DOCC(MAXBASFN*2)
      LOGICAL REMOVE, ADD, OCCNUMS_FILE, DROP_ORB
      CHARACTER*80 Blank, FNAME
C
      COMMON /POPUL/ NOCC(8,2)

      DATA ONE /1.0/
      DATA TWO /2.0/
C
      INDX2(I,J,N)=I+(J-1)*N

      IUNIT=5
C
C Template orbital convention. The electrons are addeded to the 
C alpha virtuals and removed from beta occupied.
c
      REMOVE   = .FALSE. 
      ADD      = .FALSE.
      DROP_ORB = .FALSE.

      CALL GFNAME('OCCNUMS',FNAME,ILENGTH)
      INQUIRE(FILE=FNAME(1:7), EXIST=OCCNUMS_FILE)

      IF (OCCNUMS_FILE) THEN
         OPEN(UNIT=IUNIT, FILE="OCCNUMS", FORM="FORMATTED")
         REWIND(IUNIT)
         READ(IUNIT, "(80a)") Blank
      ELSE
C
C Template orbitals without OCCNUM file should not do anything here.
C
         Write(6,"(a,a)") " Warning: Template orbital requested but no",
     &                    " OCCNMUS file or SCF_MAXCYX=0,"
         Write(6,"(a)")   " Calculation will proceed!"
         Write(6,*) 
         RETURN
      ENDIF
C
      IF (BLANK(1:3) .EQ. "ADD")    ADD    = .TRUE.
      IF (BLANK(1:6) .EQ. "REMOVE") REMOVE = .TRUE.

      READ(IUNIT, "(80a)") Blank
      CLOSE(IUNIT)
    
      IF (Blank(1:4) .EQ. "DROP")  DROP_ORB = .TRUE.

      CALL GETREC(20,"JOABRC","ORB_OCCA",NBAS,DOCC)
      IF (IUHF .NE.0) CALL GETREC(20,"JOABRC","ORB_OCCB",NBAS,
     &                            DOCC(NBAS+1))
C
      IF (REMOVE) THEN

      IOFF = 0
      ILOC = 1
      DO IRREP = 1, NIRREP
         NOCCA_FOR_IRREP = NOCC(IRREP, 2)
         IPOS            = NOCCA_FOR_IRREP
         DO NOCCA = 1, NOCCA_FOR_IRREP
            IF (DOCC(NBAS+IOFF+NOCCA) .EQ. 0.0D0) THEN
                IDUMMY(ILOC)  = IRREP
                IDUMMY2(ILOC) = IPOS
                IDUMMY3(ILOC) = NOCCA 
                ILOC = ILOC + 1
                IPOS = IPOS - 1 
            ENDIF 
         ENDDO
         IOFF = IOFF + NBFIRR(IRREP)
      ENDDO

      NMODIFY = ILOC - 1

      ENDIF 

#ifdef _DEBUG_LVLM
      if (remove) then
      Write(6,"(a)") "Printing from the remove block"
      Write(6,"(a)") "Irrep,begning and ending orbital indices"
      Write(6,"(6i4)") (Idummy(i),i=1,iloc-1)
      Write(6,"(6i4)") (Idummy2(i),i=1,iloc-1)
      Write(6,"(6i4)") (Idummy3(i),i=1,iloc-1)
      Write(6,*)
      endif 
#endif 
      IF (ADD) THEN

C There are no-holes in beta spin (compared to the alpha spin). It
C may be possible that particles are added to the alpha space. Lets
C determine this.
     
      IOFF = 1
      ILOC = 1
      DO IRREP = 1, NIRREP
         NOCCA_FOR_IRREP = NOCC(IRREP, 1)
         NOCCB_FOR_IRREP = NOCC(IRREP, 2)
         NOCCS           = NOCCA_FOR_IRREP + NOCCB_FOR_IRREP
         IPOS1           = NOCCB_FOR_IRREP+1
         IPOS2           = NOCCS
         IOFF            = IRPOFF(IRREP) + NOCCB_FOR_IRREP + 1
         DO NOCCA =IOFF, IRPOFF(IRREP+1)
            IF (DOCC(NOCCA) .EQ. 1.0D0) THEN
                IDUMMY(ILOC)  = IRREP
                IDUMMY2(ILOC) = IPOS1
                IDUMMY3(ILOC) = IPOS2 
                ILOC = ILOC + 1
                IPOS1 = IPOS1 + 1
            ENDIF
            IPOS2 = IPOS2 + 1
         ENDDO
      ENDDO

      NMODIFY = ILOC - 1

      ENDIF 

#ifdef _DEBUG_LVLM
      If (.not. remove) then
      Write(6,"(a)") "Printing from the add block"
      Write(6,"(a)") "Irrep,begning and ending orbital indices"
      Write(6,"(6i4)") (Idummy(i),i=1,iloc-1)
      Write(6,"(6i4)") (Idummy2(i),i=1,iloc-1)
      Write(6,"(6i4)") (Idummy3(i),i=1,iloc-1)
      Write(6,*)
      endif 
#endif 
      IF (NMODIFY .EQ. 0) Then
         Write(6,*)
         Write(6,"(a,a)") " Warning: Entered template_orbs but there",
     &                    " are no orbitals to swicth." 
         Write(6,"(a)")  " Calculation will proceed as requested!" 
         Write(6,*)
         RETURN
      Endif 
C
      DO I = 1, NMODIFY
         IRRP  = IDUMMY(I)
         IBGN  = IDUMMY2(I)
         IEND  = IDUMMY3(I)
         IF (REMOVE) THEN
            ISPN  = 2
         ELSEIF (ADD) THEN
            ISPN  = 1
         ENDIF 
         CALL SWITCH_ORBS(IRRP,IBGN,IEND,NBAS,EVAL,EVEC,LDIM2,
     &                    ISPN,IOS)

      ENDDO
C
C This should eventually be automated, but lets say you want to drop
C the orbital that became a new virtual after the electron is removed.
C The easiest thing to do is to set the eigenvalue and vector to zero.
C 
      IF (DROP_ORB) THEN
         DO I = 1, NMODIFY
            IEND = IDUMMY2(I)
            DO ISPIN = 2, 2
               IOFF1 = (ISPIN-1)*LDIM2+ISQROF(IRRP)-1+
     &                  INDX2(1,IEND,NBFIRR(IRRP))
               CALL DZERO(EVEC(IOFF1),NBFIRR(IRRP))
   
               IOFF2 = IRPOFF(IRRP) + IEND
               EVAL((ISPIN-1)*NBAS+IOFF2) = 0.0D0
            ENDDO
         ENDDO 
      ENDIF
C
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The Template orbital Alpha and Beta Eigenvalues"
      Write(6, "(6F12.6)") (EVAL(i),i=1,2*Nbas)
      Write(6,"(a)") "The Template orbital Alpha and Beta Eigenvectors"
      Write(6,*)
      ioff = 1 
      do ispin = 1, 2
         Write(6, "(a,i2)") "The eigenvectors of spin : ", ISPIN
      do irrp =1, nirrep 
         nbf = nbfirr(irrp)
         Write(6, "(a,i2)") "The eigenvectors of irrep: ", Irrp
         call output(EVEC(ioff), 1, nbf, 1, nbf, nbf, nbf, 1)
         ioff = ioff + nbfirr(irrp) * nbfirr(irrp)
      enddo
      enddo 
#endif 
C
      RETURN
      END
