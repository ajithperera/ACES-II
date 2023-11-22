      SUBROUTINE SCALEBLOCK(WMAT, NOCCA, NOCCB, NVRTA, NVRTB, NBAS,
     &                     OCC_OCC, VRT_VRT, VRT_OCC, OCC_VRT, IUHF,
     &                     FACTOR) 
C
C In order to zero-out a block of a Matrix
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL OCC_OCC, VRT_VRT, VRT_OCC, OCC_VRT
      DIMENSION WMAT(NBAS, NBAS)
C
      IF (IUHF .EQ. 1) THEN
          NOCC = NOCCA
          NVRT = NVRTA
      ELSE IF (IUHF .EQ. 2) THEN
          NOCC = NOCCB
          NVRT = NVRTB
      END IF

      IF (OCC_OCC .AND. VRT_VRT .AND. VRT_OCC) THEN
         DO I = 1, NOCC
           DO J = 1, NOCC
                WMAT(I, J) = FACTOR*WMAT(I, J)
           ENDDo
         ENDDO
         DO I = 1+NOCC, NBAS
             DO J = 1+NOCC, NBAS
                WMAT(I, J) = FACTOR*WMAT(I, J)
             ENDDo
         ENDDo

         DO I = 1+NOCC, NBAS
            DO J = 1, NOCC
                 WMAT(I, J) = FACTOR*WMAT(I, J)
            ENDDO
         ENDDO
 
      ELSE IF (OCC_OCC .AND. VRT_VRT) THEN
         DO I = 1, NOCC
             DO J = 1, NOCC
                WMAT(I, J) = FACTOR*WMAT(I, J)
             ENDDo
          ENDDO
          DO I = 1+NOCC, NBAS
              DO J = 1+NOCC, NBAS
                 WMAT(I, J) = FACTOR*WMAT(I, J)
              ENDDo
           ENDDo
      ELSE IF (VRT_OCC .AND. OCC_VRT) THEN
C
            DO I = 1+NOCC, NBAS
               DO J = 1, NOCC
                  WMAT(I, J) = FACTOR*WMAT(I, J)
               ENDDO
            ENDDO
C
            DO I = 1, NOCC
               DO J =1+NOCC, NBAS
                  WMAT(I, J) = FACTOR*WMAT(I, J)
                ENDDo
            ENDDO
   
      ELSE IF (OCC_OCC) THEN 
C
          DO I = 1, NOCC
             DO J = 1, NOCC
                WMAT(I, J) = FACTOR*WMAT(I, J)
             ENDDo
          ENDDO
C
      ELSE IF (VRT_VRT) THEN

           DO I = 1+NOCC, NBAS
              DO J = 1+NOCC, NBAS
                 WMAT(I, J) = FACTOR*WMAT(I, J)
              ENDDo
           ENDDo
C
      ElSE IF (VRT_OCC) THEN
   
            DO I = 1+NOCC, NBAS
               DO J = 1, NOCC 
                  WMAT(I, J) = FACTOR*WMAT(I, J)
               ENDDO
            ENDDO
C
      ELSE IF (OCC_VRT) THEN
C
             DO I = 1, NOCC
                DO J =1+NOCC, NBAS
                  WMAT(I, J) = FACTOR*WMAT(I, J)
                ENDDo
             ENDDO
C
      END IF
C 
#ifdef _DEBUG_LVL0
      Print*, "Matrix after zeroing choosen blocks"
      call output(wmat, 1, nbas, 1, nbas, nbas, nbas, 1)
#endif
          
      RETURN
      END

