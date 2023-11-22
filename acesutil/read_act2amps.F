      SUBROUTINE READ_ACT2AMPS(IUNIT, ILIST, ISPIN, ACT, ACTISFULL)

      IMPLICIT INTEGER (A-Z)
      LOGICAL ACTISFULL
      CHARACTER*4 ACT

#include <active_space.par>
#include <active_spacet12.com>

      READ(IUNIT,*) NT2SIZE
    
      IF (MAX_ACT_SPACE**4 .LT. NT2SIZE) THEN
        Write(6,"(1x,a,a)")"The current MAX_ACT_SPACE value of 100 is"
     &                    ," too small." 
        CALL ERREX
      ENDIF

      DO I = 1, NT2SIZE
          READ(IUNIT,99) I2(I),IR2(I),J2(I),JR2(I),A2(I),AR2(I),
     &                     B2(I),BR2(I),T2AS(I)
      ENDDO

 99   FORMAT(3X, I4, 3X, I1,3X, I4, 3X, I1, 3X, I4, 3X, I1, 3X, 
     &       I4, 3X, I1, 3X, F15.12)

#ifdef _DEBUG_LVLM
      Write(6,*) 
      Write(6,"(a)") "Active space amplitudes read in READ_ACT2AMPS"
      Write(6,"(a)") "Active space T2 amps."
      do i=1, nt2size
      write(6,99) I2(I),IR2(I),J2(I),JR2(I),A2(I),AR2(I),B2(I),
     &            BR2(I),T2AS(I)
      enddo
      Write(6,*)
#endif
      IF(ISPIN.EQ.1)THEN
        NT2_INT_SIZ=ISYMSZ(ISPIN,ISPIN+2)
      ELSEIF(ISPIN.EQ.2)THEN
        NT2_INT_SIZ=ISYMSZ(ISPIN,ISPIN+2)
      ELSEIF(ISPIN.EQ.3)THEN
        NT2_INT_SIZ=ISYMSZ(13,14)
      ENDIF
      ACTISFULL = .FALSE. 
C
      IF (NT2_INT_SIZ .EQ. NT2SIZE) THEN
#ifdef _DEBUG_LVLM
          Write(6,*) "Active and full space are identical"
#endif
         ACTISFULL = .TRUE.
         IF (ACT .EQ. "READ") CALL PUTALL(T2AS,NT2SIZ,1,ILIST)
      ENDIF

      RETURN
      END
