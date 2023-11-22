      SUBROUTINE READ_ACT1AMPS(IUNIT,ILIST,ISPIN,MODE,ACTISFULL)

      IMPLICIT INTEGER (A-Z)
      LOGICAL ACTISFULL
      CHARACTER*4 MODE

#include <active_space.par>
#include <active_spacet12.com>
C
      COMMON /SYMINF/ NSTART, NIRREP, IRREPS(255,2), DIRPRD(8,8)
      COMMON /OFF4ASPACE/ IACT_POP_OFF(8), IACT_VRT_OFF(8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)

CSSS      READ(IUNIT, "(8I4)") (IACT_POP_OFF(IRREP), IRREP=1, NIRREP)
CSSS      READ(IUNIT, "(8I4)") (IACT_VRT_OFF(IRREP), IRREP=1, NIRREP)
C
      READ(IUNIT,*) NT1SIZE

      DO I = 1, NT1SIZE
         READ(IUNIT,9) I1(I),IR1(I),A1(I),AR1(I),T1AS(I)
      ENDDO
    
 9    FORMAT(3X,I4,3X,I1,3X,I4,3X, I1, 3X, F15.12)

#ifdef _DEBUG_LVLM
      Write(6,*) 
      Write(6,"(a)") "Active space amplitudes read in READ_ACT1AMPS"
      Write(6, "(a)") "Active space T1 amps."
      do i=1, nt1size
      write(6,9) I1(I),IR1(I),A1(I),AR1(I),T1AS(I)
      enddo
      Write(6,*)
#endif
C
      NT1_INT_SIZ=IRPDPD(1,8+ISPIN)
      ACTISFULL = .FALSE.
C
      IF (NT1_INT_SIZ .EQ. NT1SIZE) THEN
#ifdef _DEBUG_LVLM
          Write(6,*) "Active and full space are identical"
#endif
          ACTISFULL = .TRUE.
          IF(MODE .EQ. "AMPS") CALL PUTLST(T1AS,1,1,1,ISPIN,ILIST)
          IF(MODE .EQ. "RESD") CALL PUTLST(T1AS,1,1,1,ISPIN+2,ILIST)
      ENDIF 

C
      RETURN
      END
