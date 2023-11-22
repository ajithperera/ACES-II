      SUBROUTINE Mod_WTWTW(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE FORMS W(MBEJ) "TWIDLE-TWIDLE", WHICH IS CONSTRUCTED
C  FROM W(MBEJ) "TWIDLE" AND AN INTERMEDIATE QUANTITY THAT WAS WRITTEN
C  OUT IN T1RING
C       =           _
C       W(MBEJ) = 2*W-WX
C
C
C  WHERE
C =
C W(mbej)=<mb||ej> + SUM t(fb,jn) <mn||ef>
C
C _
C W(mbej)=<mb||ej> + SUM t(f,j)<mb||ef> - SUM t(b,n)<mn||ej> 
C  
C         + [t(fb,jn) + t(f,j)*t(b,n)] <mn||ef>
C 
C WX(mbej)=<mb||ej> + 2 * SUM t(f,j)<mb||ef> - 2 * SUM t(b,n)<mn||ej> 
C  
C         + [t(fb,jn) + 2 * t(f,j)*t(b,n)] <mn||ef>
C
C 
C ON THE CALL, THE INTERMEDIATES ARE STORED ON THE FOLLOWING LISTS
C                                 _
C     TYPE            WX          W
C
C     MBEJ           123         254
C     mbej           124         255  (UHF only)
C     MbEj           118         256
C     mBeJ           117         257  (UHF only)
C     mBEj           125         258
C     MbeJ           126         259  (UHF only)
C
C   =
C   W IS WRITTEN TO THE LISTS WHICH HOLD WX
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,TWO
      DIMENSION ICORE(MAXCOR)
      LOGICAL CC2, HBAR_4LCCSD, ADC2
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS2/IFLAGS2(500)
      DATA ONE /1.0/
      DATA ONEM /-1.0/
      DATA TWO /2.0/

      CC2 = IFLAGS(2) .EQ. 47
      ADC2        = .FALSE.
      HBAR_4LCCSD = .FALSE.  
      HBAR_4LCCSD = (IFLAGS(2) .EQ. 6 .OR. IFLAGS2(117) .EQ. 7)
      ADC2 =  (IFLAGS(2) .EQ. 1 .AND. IFLAGS(87) .EQ. 3 .AND.
     &         IFLAGS2(117) .EQ. 10)
      IF (ADC2) HBAR_4LCCSD = .TRUE. 

      CALL SAXLST(ICORE,MAXCOR,254,123,123,TWO,ONEM)
      CALL SAXLST(ICORE,MAXCOR,258,125,125,TWO,ONEM)

      IF(IUHF.NE.0)THEN
         CALL SAXLST(ICORE,MAXCOR,256,118,118,TWO,ONEM)
         CALL SAXLST(ICORE,MAXCOR,255,124,124,TWO,ONEM)
         CALL SAXLST(ICORE,MAXCOR,257,117,117,TWO,ONEM)
         CALL SAXLST(ICORE,MAXCOR,259,126,126,TWO,ONEM)
      ELSE
         CALL MODF_QUIKAB2(ICORE,MAXCOR,IUHF)
      ENDIF
C
      RETURN
      END
