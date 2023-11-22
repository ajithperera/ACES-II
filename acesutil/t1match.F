      SUBROUTINE T1MATCH(IF, IFR, AF, AFR, T)

      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION T

#include <active_space.par>
#include <active_spacet12.com>

      DO I = 1, NT1SIZE
            
         IA  = I1 (I)
         IAR = IR1(I)
         AA  = A1 (I)
         AAR = AR1(I)
#ifdef _DEBUG_LVLM1
         Write(6,*)
         Write(6,"(a,1x,i4, 1x, i1, 1x, i4, 1x, i1)") "Full space",
     &         IF, IFR, AF, AFR
         Write(6,"(a,1x,i4, 1x, i1, 1x, i4, 1x, i1)") "Act. space",
     &         IA, IAR, AA, AAR
#endif

         IF (IA .EQ. IF .AND.  AA .EQ. AF .AND. IAR .EQ. IFR 
     &       .AND. AAR .EQ. AFR) THEN

#ifdef _DEBUG_LVLM1
              Write(6,"(a,1x,i4,1x,i1,1x,i4,1x,i1)")"Matched indices",
     &         IF, IFR, AF, AFR
             Write(6,"(a,F15.12,a,F15.12)") "T1 full ", T, 
     &                " is replaced by tailored T1 ", T1AS(I)
#endif
             T = T1AS(I)

             RETURN

          ENDIF
         
      ENDDO
   
      RETURN
      END
      
