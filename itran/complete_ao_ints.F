      SUBROUTINE COMPLETE_AO_INTS(AOINT,BUF1,KLPAIRS_4IRREP,NBAS,
     &                            NBASIS,IOFF_IRREP_4AO,
     &                            MAX_AOS_IRREP,NFIRST)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL RHF,DOALL

      INTEGER DIRPRD,DISZAO,DISZTRAO
C
      DIMENSION AOINT(1),BUF1(1)
      DIMENSION NBAS(8),KLPAIRS_4IRREP(8),IOFF_IRREP_4AO(8),
     &          NFIRST(8)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NDUMMY,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
      DATA AZERO,ONE,FOUR/0.D0, 1.D0, 4.0D0/
C
      Indx(i, j) = j+(i*(i-1))/2
#ifdef _DEBUG_LVL0
      Write(6,*)
      Print*, "The Nbas Nfirst and Isize Ioff_4IRREP_AO arrays"
      Print*, (Nbas(I), i=1, 8)
      Print*, (Ioff_Irrep_4AO(I), i=1, 8)
      Print*, "Number of basis functions", nbasis
#endif
C
       IOFF_MUNU_LAMSIG = 1
        IOFF_MUP_LAMSIG = 1
        IBGN_MUP_LAMSIG = 0
         IOFF_4AOINTS = 1
         IOFF_4DAOINT = MAX_AOS_IRREP*MAX_AOS_IRREP + 1
C
       DO IRREP=1,NIRREP
C
            NAOI = NBAS(IRREP)
          DISZAO = NAOI*(NAOI+1)/2
          NTRNI  = NFIRST(IRREP)
C          
C NAOI and is the # of AOs in a given irrep. The DISZAO is the lengths 
C of the two index AO or AO triangular arrays.
C
          KLPAIRS_READ = KLPAIRS_4IRREP(IRREP)
          IOFF_4IRREP  = IOFF_IRREP_4AO(IRREP)

          CALL DCOPY(KLPAIRS_READ*DISZAO, AOINT(IOFF_4AOINTS), 1, 
     &               BUF1(IOFF_4DAOINT), 1)
C
          DO KL = 1, KLPAIRS_READ
             KTmp_INDEX = INT(0.5d0 + 0.5d0*SQRT(DBLE(1 + 8*(KL-1))))
             LTmp_INDEX = KL - KTmp_INDEX*(KTmp_INDEX - 1)/2
c
#ifdef _DEBUG_LVLM1
        Write(*, '(a, 3(2X, I4))') "The KL, KTmp and LTmp",
     &  KL, KTmp_INDEX, LTmp_INDEX
#endif

                   DO MU_INDEX = 1, NAOI
                      DO NU_INDEX = 1, MU_INDEX
C 
                         I_INDEX = KTmp_INDEX
                         J_INDEX = LTmp_INDEX
                         K_INDEX = MU_INDEX
                         L_INDEX = NU_INDEX

                         IJ_INDEX = INDX(I_INDEX, J_INDEX)
                         KL_INDEX = INDX(K_INDEX, L_INDEX)
C
#ifdef _DEBUG_LVLM1
        If (KL_INDEX .NE. KL) 
     &Write(*, '(a, 3(2X, I4))') "The modf. IJ and KL indices 
     &with corrspn. KLPAIR and MUNU", KL,IJ_INDEX, KL_INDEX
#endif      
                         IF (KL_INDEX .NE. KL) THEN
                            IOFF_4KLINT = (KL_INDEX - 1)*DISZAO
     &                                    + IJ_INDEX + IOFF_4IRREP
                            IOFF_4IJINT = (IJ_INDEX - 1)*DISZAO
     &                                    + KL_INDEX + IOFF_4DAOINT
     &                                    - 1
C
#ifdef _DEBUG_LVLM1
      Write(*, '(a, 2(2X, I4))') "Offsets for AO update IOFF_4KL&IJ",
     & IOFF_4KLINT, IOFF_4IJINT 
#endif
                            AOINT(IOFF_4KLINT) = AOINT(IOFF_4KLINT) +
     &                                           BUF1(IOFF_4IJINT)
#ifdef _DEBUG_LVLM1
      Write(*, '(a, 2(2X, F12.5))') "Updated AOINTS",
     & AOINT(IOFF_4KLINT), AOINT(IOFF_4IJINT)
#endif
C
                         ENDIF
C
                     END DO
               END DO
C
          END DO

C end do for irrep
C
         IOFF_4AOINTS = IOFF_4AOINTS + DISZAO*KLPAIRS_READ
         IOFF_4DAOINT = IOFF_4DAOINT + NTRNI*NAOI*KLPAIRS_READ
C
       END DO
C
#ifdef _DEBUG_LVLM1
       Write(*, '(a)') "Building AO integrals is complete"
       Write(6,*)
#endif
       RETURN
       END

