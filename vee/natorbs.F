      SUBROUTINE NATORBS(EDEN,SCR1,SCR2,NMO,IUHF,ISPIN)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION EDEN(NMO,NMO),SCR1(NMO,NMO),SCR2(NMO*NMO)
      CHARACTER*8 LABSCF(2)

      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)

      DATA JUNK,ONEM,DONE,TWO /0,-1.0D0,1.0D0,2.0D0/
      DATA LABSCF/'SCFDENSA','SCFDENSB'/

      CALL DCOPY(NMO*NMO,EDEN,1,SCR1,1)

C Add in the reference density matrix contribution. Simply ones and
C zeros for the occupied orbital.

      SCFOCC = TWO
      IF (IUHF .EQ. 1) SCFOCC = DONE 
      DO IORB = 1, NOCCO(ISPIN)
         SCR1(IORB,IORB) = SCR1(IORB,IORB) + SCFOCC
      ENDDO
      CALL EIG(SCR1,SCR2,JUNK,NMO,ONEM)
      CALL DCOPY(NMO*NMO,SCR1,1,SCR2,1)

      WRITE(6,1000)
      WRITE(6,1001)
      WRITE(6,1000)

      IF (ISPIN.EQ.1) THEN
       WRITE(6,1002)
       WRITE(6,1000)
      ELSE
       WRITE(6,1003)
       WRITE(6,1000)
      ENDIF 
     
      TRACE = SSUM(NMO,SCR2,NMO+1)
      WRITE(6,2001) TRACE
      WRITE(6,2000) (SCR2(I),I=1,NMO*NMO,NMO+1)

1000  FORMAT(T3,70('-'))
1001  FORMAT(T18,'Natural orbital occupation numbers')
1002  FORMAT(T29,'Alpha spin')
1003  FORMAT(T29,'Beta spin')
2000  FORMAT((8(2X,F7.5)))
2001  FORMAT(T3,'Trace of density matrix : ',F14.10,'.')

      RETURN
      END

      


