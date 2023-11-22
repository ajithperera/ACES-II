      SUBROUTINE READ_ACTAMPS(IUNIT)

      IMPLICIT INTEGER (A-Z)

#include <active_space.par>
#include <active_space.com>

      READ(IUNIT,*) NT1SIZE

      DO I = 1, NT1SIZE
         READ(IUNIT,9) I1(I),IR1(I),A1(I),AR1(I),T1AS(I)
      ENDDO
      
      READ(IUNIT,*) NT2SIZE
      Write(6,*) NT2SIZE
    
      DO I = 1, NT2SIZE
          READ(IUNIT,99) I2(I),IR2(I),J2(I),JR2(I),A2(I),AR2(I),
     &                     B2(I),BR2(I),T2AS(I)
      ENDDO

 9    FORMAT(3X,I4,3X,I1,3X,I4,3X, I1, 3X, F15.12)
 99   FORMAT(3X, I4, 3X, I1,3X, I4, 3X, I1, 3X, I4, 3X, I1, 3X, 
     &       I4, 3X, I1, 3X, F15.12)

#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,"(a)") "Active space amplitudes read in READ_ACTAMPS"
      Write(6, "(a)") "Active space T1 amps."
      do i=1, nt1size
      write(6,9) I1(I),IR1(I),A1(I),AR1(I),T1AS(I)
      enddo
      Write(6,"(a)") "Active space T2 amps."
      do i=1, nt2size
      write(6,99) I2(I),IR2(I),J2(I),JR2(I),A2(I),AR2(I),B2(I),
     &            BR2(I),T2AS(I)
      enddo
      Write(6,*)
#endif
      RETURN
      END
