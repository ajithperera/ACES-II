
      SUBROUTINE MKDDEN_4IRREP(COEF,DDENS,SCOEF,NOCC,DOCC,
     &                         NBF,ISPIN,IUHF)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION COEF(NBF,NOCC), DDENS(NBF,NBF), DOCC(NBF),
     &          SCOEF(NBF, NOCC)
      DATA AZERO, ONE, TWO/0.D0, 1.0D0, 2.0D0/
C
#ifdef _DEBUG_LVLM
      call output(COEF, 1, nbf, 1, nocc, nbf, nocc, 1)
      write(6,*)
      write(6,"(6(1x,f10.5))") (docc(i), i=1,nbf)
      Write(6,"(a,i5)") "The NOCC: ", NOCC
      call dzero(ddens,nbf*nbf)
#endif
      CALL DCOPY(NBF*NOCC, COEF, 1, SCOEF, 1)
      DO JBF = 1, NOCC
         DO IBF = 1, NBF
            SCOEF(IBF, JBF) = DOCC(JBF)*COEF(IBF, JBF)
         ENDDO
      ENDDO
 
      CALL XGEMM('N','T',nbf,nbf,NOCC,ONE,COEF,nbf,SCOEF,nbf,
     &           AZERO,DDENS,nbf)
C
      IF (IUHF.EQ.0) CALL SSCAL(NBF*NBF,TWO,DDENS,1)
#ifdef _DEBUG_LVLM
      Write(6,*)
      call output(ddens, 1, nbf, 1, nbf, nbf, nbf, 1)
#endif


      RETURN
      END

