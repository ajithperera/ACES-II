      SUBROUTINE MODHBAR(ICORE, MAXCOR, IUHF, TRANABCD, TRANABCI)
C
C  1: IF (TRANABCD) THE HBARABCD INTEGRALS ARE CALCULATED.
C   IT IS ASSUMED THAT THE ABCD INTEGRALS ARE STORED IN FULL FORM. 
C
C  2: IF (TRANABCI) THE ABCI INTERMEDIATES ARE UPDATED: 
C
C      W(Ab,ci) =+ V(ab,cd) * t(d,i)

C  FOR THE VARIOUS SPINCASES. THIS REQUIRES bare abcd integrals
C
      INTEGER MAXCOR, ICORE, IUHF
      LOGICAL TRANABCD, TRANABCI
C
      DIMENSION ICORE(MAXCOR)
C
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /FLAGS2/IFLAGS2(500)
C
C  CALCULATE HBARABCI INTEGRALS
C      
      IF (TRANABCI) THEN
        WRITE(6,"(a)") ' ABCD contribution to HBARABCI is included.'
        CALL W5T1ABCD(ICORE,MAXCOR,IUHF)
        IFLAGS2(123) = 2
      ENDIF
C
C  CALCULATE HBARABCD INTEGRALS
C
      IF (TRANABCD) THEN
        IF (ISYTYP(1, 233) .EQ. 5) THEN
          WRITE(6,*)
          WRITE(6,*)' COMPRESSED ABCD INTEGRALS'
          WRITE(6,*)' NOT SUPPORTED IN CURRENT VERSION OF MODHBAR'
          CALL ERREX
        ELSE
          WRITE(6,"(a)") ' Full HBARABCD is explicitly constructed'
          CALL MODAIBC(ICORE, MAXCOR, IUHF, ONEM)
          CALL FORMW1(ICORE, MAXCOR, IUHF, .TRUE.)
          CALL MODAIBC(ICORE, MAXCOR, IUHF, ONE)
          CALL HBARABCD(ICORE, MAXCOR, IUHF, .TRUE.)
          IFLAGS2(122) = 2
        ENDIF
      ENDIF
C
      RETURN
      END
