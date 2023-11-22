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
      DIMENSION IDUMMY3(MAXBASFN), NSUM(8, 2),NTEST(8,2)
      DIMENSION DOCC(MAXBASFN*2),NOC_ORB(8,2)
      LOGICAL REMOVE, ADD, OCCNUMS_FILE, DROP_ORB,CLOSED_SHELL_ORBS
      LOGICAL OLDMOS_PRESENT, EXCITE
      CHARACTER*80 Blank, FNAME
C
      COMMON /POPUL/ NOCC(8,2)

      DATA ONE /1.0/
      DATA TWO /2.0/
C
      INDX2(I,J,N)=I+(J-1)*N

      IUNIT=5
      LUGSS=71

      INQUIRE(FILE='OLDMOS',EXIST=OLDMOS_PRESENT)
C 
      IF (OLDMOS_PRESENT) THEN
         OPEN(LUGSS,FILE='OLDMOS',STATUS='OLD',ACCESS='SEQUENTIAL',
     &        FORM='FORMATTED')
         READ(LUGSS,"(8(1x,I5))") (NTEST(i,1), i=1, NIRREP)
         READ(LUGSS,"(8(1x,I5))") (NTEST(i,2), i=1, NIRREP)
  
         CLOSED_SHELL_ORBS= .TRUE.
         DO I = 1, 8
           IF (NTEST(I,1) .NE.NTEST(I,2)) CLOSED_SHELL_ORBS = .FALSE.
         ENDDO 
         CLOSE(LUGSS)

         IF (.NOT. CLOSED_SHELL_ORBS) Then
            Write(6,"(a,a,a)") " The template orbitals must be from",
     &                         " closed shell state; Try a closest",
     &                         " closed shell state to generate",
     &                         " template orbitals."
            Return 
         Endif 
      ELSE 
        RETURN
      ENDIF 
C
C Template orbital convention. The electrons are addeded to the 
C alpha virtuals and removed from beta occupied.
c
      REMOVE   = .FALSE. 
      ADD      = .FALSE.
      DROP_ORB = .FALSE.
      EXCITE   = .FALSE.

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
     &                    " OCCNMUS file. Calculation will proceed!"
         Write(6,*) 
         RETURN
      ENDIF
C
      IF (BLANK(1:3) .EQ. "ADD")    ADD    = .TRUE.
      IF (BLANK(1:6) .EQ. "REMOVE") REMOVE = .TRUE.
      IF (BLANK(1:6) .EQ. "EXCITE") EXCITE = .TRUE.

      WRITE(6,"(2a,a6,1x,a)") " The template orbital calculation",
     &                      " is performed with the ",BLANK, "option."
      WRITE(6,"(a,a)") " EMPTY request means no electrons were added", 
     &               " or removed"
      WRITE(6,*) 
      READ(IUNIT, "(80a)") Blank

      READ(IUNIT,10,END=19) (NOC_ORB(i, 1), i=1, NIRREP)
      READ(IUNIT,10,END=19) (NOC_ORB(i, 2), i=1, NIRREP)
 
      READ(IUNIT, "(80a)") Blank

 10   FORMAT(16I5) 
 19   CLOSE(IUNIT)

#ifdef _DEBUG_LVL0
         Write(6,"(a)") "The Alpha and Beta no. of occupied orbitals"
         Write(6,*) (NOC_ORB(i,1), i=1, NIRREP)
         Write(6,*)
         Write(6,*) (NOC_ORB(i,2), i=1, NIRREP)
#endif 
    
      IF (Blank(1:4) .EQ. "DROP")  DROP_ORB = .TRUE.

      CALL GETREC(20,"JOABRC","ORB_OCCA",NBAS,DOCC)
      IF (IUHF .NE.0) CALL GETREC(20,"JOABRC","ORB_OCCB",NBAS,
     &                            DOCC(NBAS+1))
C
      MAXDIFF = 1
      DO IRREP = 1, NIRREP
         DO ISPIN = 1, 2
            NOCC_DIFF = NOC_ORB(IRREP,ISPIN) - NOCC(IRREP,ISPIN)
            MAXDIFF  = MAX(NOCC_DIFF,MAXDIFF)
         ENDDO
      ENDDO 
      IF (MAXDIFF .NE. 1) THEN
         WRITE(6,"(a,a)") " More than one electron transfer per irrep",
     +                  " is not currently allowed." 
         CALL ERREX 
      ENDIF 

      IF (REMOVE) THEN

C This assume that the electrons are removed from the beta orbitals. 

      IOFF = 0
      ILOC = 1
      DO IRREP = 1, NIRREP
         NOCCA_FOR_IRREP = NOC_ORB(IRREP, 2)
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

#ifdef _DEBUG_LVL0
      if (remove) then
      Write(6,"(a)") "Printing from the remove block"
      Write(6,"(a)") "Irrep,begning and ending orbital indices"
      Write(6,"(6i4)") (Idummy(i),i=1,iloc-1)
      Write(6,"(6i4)") (Idummy3(i),i=1,iloc-1)
      Write(6,"(6i4)") (Idummy2(i),i=1,iloc-1)
      Write(6,*)
      endif
#endif
      ENDIF 

      IF (ADD) THEN

C This assume that we add electrons to the Alpha spin.

      IOFF = 1
      ILOC = 1
      DO IRREP = 1, NIRREP
         NOCCA_FOR_IRREP = NOC_ORB(IRREP, 1)
         NOCCB_FOR_IRREP = NOC_ORB(IRREP, 2)
         NOCCS           = NOCCA_FOR_IRREP + NOCCB_FOR_IRREP
         IPOS1           = NOCCB_FOR_IRREP 
         IPOS2           = NOCCB_FOR_IRREP + 1
         IOFF            = IRPOFF(IRREP) + NOCCB_FOR_IRREP + 1
         DO NOCCA =IOFF, IRPOFF(IRREP+1)
            IF (DOCC(NOCCA) .EQ. 1.0D0) THEN
                IDUMMY(ILOC)  = IRREP
                IDUMMY2(ILOC) = IPOS1
                IDUMMY3(ILOC) = IPOS2
                ILOC  = ILOC + 1
                IPOS1 = IPOS1 + 1
            ENDIF
            IPOS2 = IPOS2 + 1
         ENDDO
      ENDDO

      NMODIFY = ILOC - 1

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Printing from the add block"
      Write(6,"(a)") "Irrep,begning and ending orbital indices"
      Write(6,"(6i4)") (Idummy(i),i=1,iloc-1)
      Write(6,"(6i4)") (Idummy3(i),i=1,iloc-1)
      Write(6,"(6i4)") (Idummy2(i),i=1,iloc-1)
      Write(6,*)
#endif
     
      ENDIF 

      IF (ADD .OR. REMOVE .AND. (NMODIFY .EQ. 0)) Then
         Write(6,*)
         Write(6,"(a,a)") " Warning: Entered template_orbs with ADD",
     &                    "/REMOVE options but there are no"
         Write(6,"(a)")  " orbitals to switch."
         Write(6,"(a)")  " Calculation will proceed as requested!"
         Write(6,*)
         RETURN
      Endif

      IF (REMOVE) THEN
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
     &                       ISPN,IOS)

         ENDDO
      ENDIF 

      IF (ADD) THEN
         DO I = 1, NMODIFY
            IRRP  = IDUMMY(I)
            IEND  = IDUMMY2(I)
            IBGN  = IDUMMY3(I)
            IF (REMOVE) THEN
               ISPN  = 2
            ELSEIF (ADD) THEN
               ISPN  = 1
            ENDIF
            CALL SWITCH_ORBS(IRRP,IBGN,IEND,NBAS,EVAL,EVEC,LDIM2,
     &                       ISPN,IOS)

         ENDDO
      ENDIF 
C
      IF (EXCITE) THEN

C Alpha block

      ILOC1 = 1
      ILOC2 = 1
      IOFF  = 1
      DO IRREP = 1, NIRREP
         NOCCA_FOR_IRREP_N = NOC_ORB(IRREP,1)
         NOCCA_FOR_IRREP_O = NOCC(IRREP,1)
#ifdef _DEBUG_LVL0
      Write(6,"(a,a,1x,3I3)")" IRREP,NOCCA_FOR_IRREP_N",
     +                       "NOCCA_FOR_IRREP_O:",IRREP,
     +                        NOCCA_FOR_IRREP_N,NOCCA_FOR_IRREP_O
#endif
 
         IF (NOCCA_FOR_IRREP_N .GT. NOCCA_FOR_IRREP_O) THEN
            DO NOCCA =IOFF,  NOCCA_FOR_IRREP_O + (IOFF-1)

#ifdef _DEBUG_LVL0
      Write(6,"(a,1x,i3,1x,f3.1,1x,i2)") " Blk-1:",nocca,
     +         docc(nocca),NOCCA_FOR_IRREP_O +(Ioff-1)
#endif
               IF (DOCC(NOCCA) .EQ. 0.0D0 .AND. 
     +             NOCCA .LE. NOCCA_FOR_IRREP_O + (IOFF-1)) THEN
                   IDUMMY(ILOC1)  = IRREP
                   IDUMMY3(ILOC1) = NOCCA
                   ILOC1  = ILOC1 + 1
               ENDIF
            ENDDO

            IBGN  = IRPOFF(IRREP) + NOCCA_FOR_IRREP_N
            DO NOCCA =IBGN,  NOCCA_FOR_IRREP_N + (IBGN-1)

#ifdef _DEBUG_LVL0
      Write(6,"(a,1x,i3,1x,f3.1,1x,i2)") " Blk-2:",nocca,
     +         docc(nocca),NOCCA_FOR_IRREP_N
#endif
               IF (DOCC(NOCCA) .EQ. 1.0D0) THEN
                   IDUMMY(ILOC2)  = IRREP
                   IDUMMY2(ILOC2) = NOCCA
                   ILOC2  = ILOC2 + 1
               ENDIF
            ENDDO
         ENDIF 
         IOFF = IRPOFF(IRREP+1) + 1
      ENDDO

      NMODIFY_A = ILOC2 - 1

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Printing from the excite alpha block"
      Write(6,"(a)") "Irrep,begning and ending orbital indices"
      Write(6,"(6i4)") (Idummy(i),i=1,nmodify_a)
      Write(6,"(6i4)") (Idummy3(i),i=1,nmodify_a)
      Write(6,"(6i4)") (Idummy2(i),i=1,nmodify_a)
      Write(6,*)
#endif
      DO I = 1, NMODIFY_A
         IRRP  = IDUMMY(I)
         IBGN  = IDUMMY2(I)
         IEND  = IDUMMY3(I)
         ISPN  = 1
         CALL SWITCH_ORBS(IRRP,IBGN,IEND,NBAS,EVAL,EVEC,LDIM2,
     &                    ISPN,IOS)
      ENDDO

C Beta block 

      ILOC1 = 1
      ILOC2 = 1
      IOFF  = 1
      DO IRREP = 1, NIRREP
         NOCCB_FOR_IRREP_N = NOC_ORB(IRREP,2)
         NOCCB_FOR_IRREP_O = NOCC(IRREP,2)
#ifdef _DEBUG_LVL0
      Write(6,"(a,a,1x,3I3)") "IRREP,NOCCB_FOR_IRREP_N",
     +                        "NOCCB_FOR_IRREP_O:",IRREP,
     +                         NOCCB_FOR_IRREP_N,NOCCB_FOR_IRREP_O
#endif
         IF (NOCCB_FOR_IRREP_N .GT. NOCCB_FOR_IRREP_O) THEN
            DO NOCCB =IOFF,  NOCCB_FOR_IRREP_O + (IOFF-1)
#ifdef _DEBUG_LVL0
      Write(6,"(a,1x,i3,1x,f3.1,1x,i2)") " Blk-1:",noccb,
     +         docc(noccb),NOCCB_FOR_IRREP_O +(Ioff-1)
#endif
               IF (DOCC(NBAS+NOCCB) .EQ. 0.0D0 .AND.
     +             NOCCB .LE. NOCCB_FOR_IRREP_O + (IOFF-1)) THEN
                   IDUMMY(ILOC1)  = IRREP
                   IDUMMY3(ILOC1) = NOCCB
                   ILOC1  = ILOC1 + 1
               ENDIF
            ENDDO

            IBGN  = IRPOFF(IRREP) + NOCCB_FOR_IRREP_N
            DO NOCCB =IBGN, NOCCB_FOR_IRREP_N + (IBGN-1)

#ifdef _DEBUG_LVL0
      Write(6,"(a,1x,i3,1x,f3.1,1x,i2)") " Blk-2:",noccb,
     +         docc(nbas+noccb),NOCCA_FOR_IRREP_N
#endif
               IF (DOCC(NBAS+NOCCB) .EQ. 1.0D0) THEN
                   IDUMMY(ILOC2)  = IRREP
                   IDUMMY2(ILOC2) = NOCCB
                   ILOC2  = ILOC2 + 1
               ENDIF
            ENDDO
        ENDIF 
        IOFF = IRPOFF(IRREP+1) + 1

      ENDDO

      NMODIFY_B = ILOC2 - 1
      IF (EXCITE .AND. (NMODIFY_A .EQ.0 .AND. NMODIFY_B .EQ.0)) Then
         Write(6,*)
         Write(6,"(a,a)") " Warning: Entered template_orbs with",
     &                    " EXCITE options but there are no"
         Write(6,"(a)")  "  orbitals to switch."
         Write(6,"(a)")  " Calculation will proceed as requested!"
         Write(6,*)
         RETURN
      Endif

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Printing from the excite beta block"
      Write(6,"(a)") "Irrep,begning and ending orbital indices"
      Write(6,"(6i4)") (Idummy(i),i=1,nmodify_b)
      Write(6,"(6i4)") (Idummy3(i),i=1,nmodify_b)
      Write(6,"(6i4)") (Idummy2(i),i=1,nmodify_b)
      Write(6,*)
#endif
      DO I = 1, NMODIFY_B
         IRRP  = IDUMMY(I)
         IBGN  = IDUMMY2(I)
         IEND  = IDUMMY3(I)
         ISPN  = 2
         CALL SWITCH_ORBS(IRRP,IBGN,IEND,NBAS,EVAL,EVEC,LDIM2,
     &                    ISPN,IOS)
      ENDDO

      ENDIF
C
C This should eventually be automated, but lets say you want to drop
C the orbital that became a new virtual after the electron is removed.
C The easiest thing to do is to set the eigenvalue and vector to zero.
C 
      IF (DROP_ORB .AND. REMOVE) THEN
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
