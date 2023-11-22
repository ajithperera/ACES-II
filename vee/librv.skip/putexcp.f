      SUBROUTINE PUTEXCP(IUHF, SCR, MAXCOR, IRREPX, 
     &   IPROJECT, LISTPTRN, COLPTRN, NSIZEC, IOPT)
C
C  THIS SUBROUTINE PUTS CALCULATED EXCITATION PATTERNS
C  ON THE APPROPRIATE COLUMNS OF LISTH0
C  
C  IPROJECT INDICATES THE USE OF THE PROJECTION:
C
C     IPROJECT = 0:  NO EXCITATION PATTERN USED
C     IPROJECT = 1:  USED TO DETERMINE EIGENVECTORS OF INTEREST
C     IPROJECT = 2:  EOMCC MATRIX IS REALLY PROJECTED ON SUBSPACE
C     IPROJECT = 3:  CONVERGENCE DISTURBING ELEMENTS ARE PROJECTED OUT
C
C  IOPT SPECIFIES IF COMPONENTS OF INTEREST ARE GENERATED (IOPT = 1)
C                 OR ACTUAL PROJECTION IS DETERMINED  (IOPT = 2)
c
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION SCR(MAXCOR), E, ZILCH, ONE, HUGE, TOL
      LOGICAL PRINT
C
      DATA ZILCH, ONE /0.0D0, 1.0D0/
C
      IF (IPROJECT .EQ. 0) RETURN
C
        CALL LOADVEC1(IRREPX,SCR,MAXCOR,IUHF,490,0,443,NSIZEC,
     &              .FALSE.)
C
C PUT PROJECTION VECTOR ON APPROPRIATE LIST
C
        IF (IOPT .EQ. 1) THEN
          CALL PUTLST(SCR,COLPTRN,1,1,1, LISTPTRN)
        ELSEIF (IOPT .EQ. 2) THEN
          IF (IPROJECT .EQ. 2) THEN
            CALL PUTLST(SCR,COLPTRN+1,1,1,1, LISTPTRN)
          ENDIF
          IF (IPROJECT .EQ. 3) THEN
            I000 = 1
            I010 = I000 + NSIZEC
            I020 = I010 + NSIZEC
C
C READ IN CURRENT EXCITATION PATTERN
C          
            CALL GETLST(SCR(I000), COLPTRN,1,1,1,LISTPTRN)
C
C  READ IN DENOMINATORS
C
            CALL GETLST(SCR(I010), 1,1,1,1,472)
C
C  PUT INTERESTING DENOMINATORS IN SCR(I020), SET OTHERS TO HUGE
C          
            HUGE = 1.0D8
            DO I = 0, NSIZEC - 1
              SCR(I020+I) = SCR(I010+I) + (ONE - SCR(I000+I))*HUGE
            ENDDO
C
C  FIND MINIMUM AMONG INTERSTING DENOMINATORS
C
            INDEX=ISAMIN(NSIZEC,SCR(I020),1)
            E = SCR(I010 + INDEX - 1)
            PRINT = .TRUE.
            IF (PRINT) THEN
              WRITE(6,*) ' REFERENCE ENERGY IN PUTEXCP ', E
            ENDIF
C
C  NOW DETERMINE TRUE PROJECTION MANIFOLD
C
            DO I = 0, NSIZEC - 1
              SCR(I020+I) = ONE
            ENDDO
            TOL = 1.0D0
C
C  PRESENT SCHEME IS TAILORED TO CORE-EXCITATIONS
C
            DO I = 0, NSIZEC - 1
              IF ( (SCR(I000+I) .LT. 0.75) .AND.
     &           (ABS(E-SCR(I010+I)) .LT. TOL) ) THEN
                SCR(I020+I) = ZILCH
              ENDIF
            ENDDO
C
            CALL PUTLST(SCR(I020),COLPTRN+1,1,1,1, LISTPTRN)
C
          ENDIF
        ELSE
          WRITE(6,*) ' UNKNOWN OPTION IN PUTEXCP ' , IOPT
          CALL ERREX
        ENDIF
C
      RETURN
      END