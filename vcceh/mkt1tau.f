C
      SUBROUTINE MKT1TAU(T2, T1A, T1B, DISSYT, NUMSYT, POP1, POP2,
     &                   VRT1, VRT2, IRREP, IRREPXA, IRREPXB,
     &                   ISPIN, FACT)
C
C This subroutine forms the symmetry packed TAU(AB,IJ),
C
C  TAU(AB,IJ) =  T1(A,I)*T1(B,J)*FACT - T1(A,J)*T1(B,I)*FACT
C
C For ISPIN = 1 (AAAA case) and ISPIN = 2 (BBBB case)
C the equation given applies directly. For ISPIN = 3
C (ABAB case) it reduces to 
C
C  TAU(Ab,Ij) = T1(A,I)*T1(b,j)*FACT
C 
C Note that the symmetry packing is used and that the symmetry
C information is also used in order to decide if there are any 
C single contribution or not.
C For the ABAB spin case, there are only contributions when the 
C irrep of A is equal to the irrep of I (similarly, 
C IRREPB equal to IRREPJ), which is forced by the requirement 
C that the T2 amplitudes are total symmetric in the AAAA
C and BBBB spin cases. There are contributions either if
C IRREPA equals IRREPI (first term) or IRREPA equals IRREPJ
C (second term).
C
C Code by JULY/90 JG and modified for the quadratic term by Ajith 07/94.
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER POP1,POP2,VRT1,VRT2,DIRPRD,DISSYT
      DIMENSION T2(DISSYT,NUMSYT),T1A(1),T1B(1)
      DIMENSION POP1(8),POP2(8),VRT1(8),VRT2(8)
      DIMENSION IOFFT2O(8),IOFFT2V(8),IOFFT1A(8),
     &          IOFFT1B(8)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
C
      NNM1O2(I) = ((I-1)*I)/2
C     
C Lets zero-out the space for T2 at the begining since
C every call to MKT1TAU in the calculation of quadratic
C term assume to have zero T2 vector. This is just to
C be safer, if somebody forget to pass the zero T2 vector from
C the calling routine. Notice the differnce from the CC code.
C
      CALL ZERO(T2, DISSYT*NUMSYT)
C
      IF (ISPIN .LT. 3) THEN      
C
C AAAA or BBBB spin cases
C
C Determine offset for T2 and T1
C
         IOFFT2O(1) = 0
         IOFFT2V(1) = 0
         IOFFT1A(1) = 0
         IOFFT1B(1) = 0
C
         DO 1000 IRREPO = 1, NIRREP - 1
            IRREPV = DIRPRD(IRREPO, IRREPXA)
            IOFFT1A(IRREPO + 1) = IOFFT1A(IRREPO) + POP1(IRREPO)*
     &                            VRT1(IRREPV)
 1000    CONTINUE
C
         DO 2000 IRREPO = 1, NIRREP - 1
            IRREPV = DIRPRD(IRREPO, IRREPXB)
            IOFFT1B(IRREPO + 1) = IOFFT1B(IRREPO) + POP1(IRREPO)*
     &                            VRT1(IRREPV)
 2000    CONTINUE
C
         IF(IRREP .EQ. 1) THEN
C
            DO 2 IRREPJ = 1, NIRREP - 1
               IOFFT2O(IRREPJ + 1) = IOFFT2O(IRREPJ) +
     &                               NNM1O2(POP1(IRREPJ))
               IOFFT2V(IRREPJ + 1) = IOFFT2V(IRREPJ)+
     &                               NNM1O2(VRT1(IRREPJ))
 2          CONTINUE
         ELSE
            DO 3 IRREPJ = 1, NIRREP - 1
               IRREPI = DIRPRD(IRREP, IRREPJ)
               IF(IRREPI .LT. IRREPJ) THEN
                  IOFFT2O(IRREPJ + 1) = IOFFT2O(IRREPJ) + 
     &                                  POP1(IRREPI)*POP1(IRREPJ)
                  IOFFT2V(IRREPJ + 1) = IOFFT2V(IRREPJ) + 
     &                                  VRT1(IRREPI)*VRT1(IRREPJ)
                  ELSE
                  IOFFT2O(IRREPJ + 1) = IOFFT2O(IRREPJ)
                  IOFFT2V(IRREPJ + 1) = IOFFT2V(IRREPJ)
               ENDIF
 3          CONTINUE
         ENDIF
C     
         DO 100 IRREPJ = 1, NIRREP
C     
            IRREPI  = DIRPRD(IRREP, IRREPJ)
            IRREPAO = DIRPRD(IRREPI, IRREPXA)
            IRREPBO = DIRPRD(IRREPJ, IRREPXB)
C
            NUMJ   = POP1(IRREPJ)
            INDIJ  = IOFFT2O(IRREPJ)
            INDJ0  = IOFFT1B(IRREPJ)
            INDI0  = IOFFT1A(IRREPI)
C     
            IF (IRREPI .EQ. IRREPJ) THEN
C
               IF (IRREPAO .NE. IRREPBO) THEN
                  WRITE(6, *) ' @-MKT1TAU SYMMETRY VIOLATION '
                  CALL ERREX
               ENDIF
C     
               IF(NNM1O2(NUMJ) .EQ. 0) GO TO 100
C     
               NUMB   = VRT1(IRREPBO)
               NUMA   = VRT1(IRREPAO)
C
               INDAB0 = IOFFT2V(IRREPBO)
C     
               DO 10 J = 2, NUMJ
                  INDJ = INDJ0 + (J-1)*NUMB
                  DO 20 I = 1, J-1
                     INDI = INDI0 + (I-1)*NUMA
C     
                     INDIJ = INDIJ + 1
                     INDAB = INDAB0
C 
                     DO 30 IB = 2, NUMB
                        INDBI = INDI + IB
                        INDBJ = INDJ + IB
                        DO 40 IA = 1, IB - 1
                           INDAI = INDI + IA
                           INDAJ = INDJ + IA
                           INDAB = INDAB + 1
C     
                           T2(INDAB,INDIJ) = T2(INDAB,INDIJ) +  
     &                                       FACT*T1A(INDAI)*T1B(INDBJ)
     &                                      -FACT*T1B(INDAJ)*T1A(INDBI)
C                           
 40                     CONTINUE
 30                  CONTINUE
 20               CONTINUE
 10            CONTINUE
C     
            ELSE IF (IRREPI .LT. IRREPJ) THEN

               IRREPAO = DIRPRD(IRREPI, IRREPXA)
               IRREPBO = DIRPRD(IRREPJ, IRREPXB)

               IF (IRREPAO .LT. IRREPBO) THEN
C
C The normal case
C
                  NUMI   = POP1(IRREPI)
                  INDIJ  = IOFFT2O(IRREPJ)
                  INDAB0 = IOFFT2V(IRREPBO)
                  NUMB   = VRT1(IRREPBO)
                  NUMA   = VRT1(IRREPAO)
C
                  DO 50 J = 1, NUMJ
C
                     INDJ = IOFFT1B(IRREPJ) + (J-1)*NUMB
C     
                     DO 60 I = 1, NUMI
C     
                        INDI  = IOFFT1A(IRREPI) + (I-1)*NUMA
                        INDAB = INDAB0
                        INDIJ = INDIJ + 1         
C     
                        DO 70 IB = 1, NUMB
C
                           INDBJ = INDJ + IB
C
                           DO 80 IA = 1, NUMA
C
                              INDAI = INDI + IA
                              INDAB = INDAB + 1
C
                              T2(INDAB,INDIJ) = T2(INDAB,INDIJ) +
     &                                          FACT*T1A(INDAI)*
     &                                          T1B(INDBJ)
 80                        CONTINUE
 70                     CONTINUE
 60                  CONTINUE
 50               CONTINUE
C
               ENDIF
C
               IRREPAO = DIRPRD(IRREPJ, IRREPXA)
               IRREPBO = DIRPRD(IRREPI, IRREPXB)

               IF (IRREPAO .LT. IRREPBO) THEN
C
C Calculate T(AB,IJ) = - T(B,I)*T(A,J)
C
                  NUMI   = POP1(IRREPI)
                  INDIJ  = IOFFT2O(IRREPJ)
                  INDAB0 = IOFFT2V(IRREPBO)
                  NUMB   = VRT1(IRREPBO)
                  NUMA   = VRT1(IRREPAO)
C
                  DO 150 J = 1, NUMJ
C
                     INDJ = IOFFT1B(IRREPJ) + (J-1)*NUMA
C     
                     DO 160 I = 1, NUMI
C     
                        INDI  = IOFFT1A(IRREPI) + (I-1)*NUMB
                        INDAB = INDAB0
                        INDIJ = INDIJ + 1         
C
                        DO 170 IB = 1, NUMB
C     
                           INDBI = INDI + IB
C
                           DO 180 IA = 1, NUMA
C
                              INDAJ = INDJ + IA
C
                              INDAB = INDAB + 1
                              
                              T2(INDAB,INDIJ) = T2(INDAB,INDIJ) -
     &                                          FACT*T1A(INDBI)*
     &                                          T1B(INDAJ)
 180                       CONTINUE
 170                    CONTINUE
 160                 CONTINUE
 150              CONTINUE
C
               ENDIF
C
            ENDIF
C
 100     CONTINUE
C
      ELSE
C     
C ABAB spin case
C     
C First get the offsets of T1 and T2
C     
         IOFFT2O(1) = 0
         IOFFT2V(1) = 0
         IOFFT1A(1) = 0
         IOFFT1B(1) = 0
C     
         DO 101 IRREPJ = 1, NIRREP - 1
C
            IRREPI = DIRPRD(IRREPJ, IRREP)
C
            IOFFT2O(IRREPJ + 1) = IOFFT2O(IRREPJ) +
     &                            POP2(IRREPJ)*POP1(IRREPI)
            IOFFT2V(IRREPJ + 1) = IOFFT2V(IRREPJ) + 
     &                            VRT2(IRREPJ)*VRT1(IRREPI)
 101     CONTINUE
C
         DO 102 IRREPO = 1, NIRREP - 1
            IRREPV = DIRPRD(IRREPO, IRREPXA)
            IOFFT1A(IRREPO + 1) = IOFFT1A(IRREPO) + POP1(IRREPO)*
     &                            VRT1(IRREPV)
 102     CONTINUE
C
         DO 103 IRREPO = 1, NIRREP - 1
            IRREPV = DIRPRD(IRREPO, IRREPXB)
            IOFFT1B(IRREPO + 1) = IOFFT1B(IRREPO) + POP2(IRREPO)*
     &                            VRT2(IRREPV)
 103     CONTINUE
C     
          DO 200 IRREPJ = 1, NIRREP
C     
            IRREPI  = DIRPRD(IRREP, IRREPJ)
            IRREPAO = DIRPRD(IRREPI, IRREPXA)
            IRREPBO = DIRPRD(IRREPJ, IRREPXB)
C
            NUMJ   = POP2(IRREPJ)
            NUMI   = POP1(IRREPI)
            INDIJ  = IOFFT2O(IRREPJ)
            INDAB0 = IOFFT2V(IRREPBO)
            NUMA   = VRT1(IRREPAO)           
            NUMB   = VRT2(IRREPBO)
C     
            DO 210 J = 1, NUMJ
C
               INDJ = IOFFT1B(IRREPJ) + (J-1)*NUMB
C
               DO 220 I = 1, NUMI
C
                  INDI = IOFFT1A(IRREPI) + (I-1)*NUMA
                  INDIJ = INDIJ + 1
                  INDAB = INDAB0
C     
                  DO 230 IB = 1, NUMB
C
                     INDBJ = INDJ + IB
C
                     DO 240 IA = 1, NUMA
C
                        INDAB = INDAB + 1
                        INDAI = INDI + IA
C
                        T2(INDAB, INDIJ) = T2(INDAB,INDIJ) + 
     &                                     FACT*T1A(INDAI)*T1B(INDBJ)

 240                 CONTINUE
 230              CONTINUE
 220           CONTINUE
 210        CONTINUE
 200     CONTINUE
C
      ENDIF
C
      RETURN
      END
