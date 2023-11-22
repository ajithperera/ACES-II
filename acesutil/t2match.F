      SUBROUTINE T2MATCH(IF, IFR, JF, JFR, AF, AFR, BF, BFR, T)

      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION T
C
#include <active_space.par>
#include <active_spacet12.com>

      DO I = 1, NT2SIZE
            
         IA  = I2 (I)
         IAR = IR2(I)
         JA  = J2 (I)
         JAR = JR2(I)
         AA  = A2 (I)
         AAR = AR2 (I)
         BA  = B2 (I)
         BAR = BR2 (I)
#ifdef _DEBUG_LVLM1
         Write(6,*)
         Write(6,"(a, 1x, i4, 1x, i1, 1x, i4, 1x, i1, 1x, i4, 1x, i1, 
     &             1x, i4, 1x, i1)") "Full space", IF, IFR, JF, JFR,
     &             AF, AFR, BF, BFR
         Write(6,"(a, 1x, i4, 1x, i1, 1x, i4, 1x, i1, 1x, i4, 1x, i1, 
     &             1x, i4, 1x, i1)") "Act. space", IA, IAR, JA, JAR,
     &             AA, AAR, BA, BAR

#endif
         IF (IA .EQ. IF .AND. JA .EQ. JF .AND. AA .EQ. AF .AND.
     &       BA .EQ. BF .AND. IAR .EQ. IFR .AND. JAR .EQ. JFR
     &      .AND. AAR .EQ. AFR .AND. BAR .EQ. BFR) THEN

#ifdef _DEBUG_LVLM1
        Write(6,"(a, 1x, i4, 1x, i1, 1x, i4, 1x, i1, 1x, i4, 1x, i1,
     &             1x, i4, 1x, i1)") "Matched indices", IF, IFR, JF,
     &             JFR, AF, AFR, BF, BFR
             Write(6,"(a,F13.7,a,F13.7)") "T2 full ", T,
     &                " is replaced by tailored T2 ", T2AS(I)
#endif

             T = T2AS(I)
             RETURN 
         ENDIF
         
      ENDDO
   
      RETURN
      END
      
