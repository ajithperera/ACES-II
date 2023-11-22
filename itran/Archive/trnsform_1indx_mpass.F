      SUBROUTINE TRNSFORM_1INDX_MPASS(EVEC,AOINT,BUF1,MUP_TMP_KEEP,
     &                                KLPAIRS_READ,NFIRST,NSTART,
     &                                NEND,ISPIN,NBAS,NMO,NORB,NBASIS)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL RHF,DOALL

      INTEGER DIRPRD,DISZAO,DISZMO,DISZTRAO
      DOUBLE PRECISION MUP_TMP_KEEP(NBASIS*NBASIS)
C
      DIMENSION EVEC(1),AOINT(1),BUF1(1)
      DIMENSION NBAS(8),NMO(8),NFIRST(8),NSTART(8),NEND(8)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/AOOFST/INDOCC(8,2)
      COMMON/SYMINF/NDUMMY,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
      DATA AZERO,ONE,FOUR/0.D0, 1.D0, 4.0D0/
C
      Indx(i, j) = j+(i*(i-1))/2
#ifdef _DEBUG_LVLM1
      Print*, "The nstartx and nendx arrays"
      Print*, (nStart(I), i=1, 8)
      Print*, (nEnd(I), i=1, 8)
      Write(6,*)
      Print*, "The Nbas Nfirst and Isize arrays"
      Print*, (Nbas(I), i=1, 8)
      Print*, (NFirst(I), i=1, 8)
      Write(6,*)
#endif
C
       IOFF_MUNU_LAMSIG = 1
        IOFF_MUP_LAMSIG = 1
        IBGN_MUP_LAMSIG = 0
C
       DO IRREP=1,NIRREP
C
            NAOI = NBAS(IRREP)
            NMOI = NMO(IRREP)
          DISZAO = NAOI*(NAOI+1)/2
          DISZMO = NAOI*(NMOI+1)/2
C          
C NAOI and NMOI are the # of AOs and and the total coorelated MOs
C in a given irrep. The DISZAO and DISZMO are the lengths of the 
C two index AO or AO-MO triangular arrays.
C
            NTRNI  = NFIRST(IRREP)
          DISZTRAO = NAOI*(NRNI+1)/2
C
          NSTARTI  = NSTART(IRREP)
          NENDI    = NEND(IRREP)
C
          CALL DCOPY(KLPAIRS_READ*NAOI*NAOI, AOINT, 1, BUF1, 1)
C
          DO KL = 1, KLPAIRS_READ
             KTmp_INDEX = INT(0.5d0 + 0.5d0*SQRT(DBLE(1 + 8*(KL-1))))
             LTmp_INDEX = KL - KTmp_INDEX*(KTmp_INDEX - 1)/2
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
     &                                    + IJ_INDEX
                            IOFF_4IJINT = (IJ_INDEX - 1)*DISZAO
     &                                    + KL_INDEX

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

C
C NFIRSTI keep the range of MOS that need to be transformed. It is 
C identical to NMOI when no orbitals are dropped. Depending on 
C whether NOABCD set NFIRSTI has the full or occupide MOs set. 
C
C NSTARTI and NENDI are the starting and ending addresses of the range
C of MOs that is C being transformed. NSTARTI is the starting address 
C of the range if the entire set can not be handled in a one pass.
C 
C
          IF (NENDI .GE. NSTARTI) Then 
             NUMDIS=NMOI*(NENDI+1-NSTARTI)
C
C DISZAO and DISZTRAO are the lengths of the two index AO or AO-MO triangular 
C arrays.
C
C The first index can be over all the orbitals or occupied orbitals
C if ABCDTYPE is set to AOBASIS (the NMO is set accordingly).
C The AOINT(mu,nu,lam,sig) is ordered as for given lam > sig all 
C mu > nu. First expand the AOINT(mu>nu|lam>sig> to 
C <mu,nu|lam>sig> (assuming that we have enough memory to read at
C least one KL pair).

             DO KLPAIRS = 1, KLPAIRS_READ
                CALL EXPND2(AOINT(IOFF_MUNU_LAMSIG), BUF1, NAOI)
C
#ifdef _DEBUG_LVL0
         Print*, "The AO ints for a KL PAIRS:", "KL = ", KLPAIRS
         Call output(BUF1, 1, NAOI, 1, NAOI, NAOI, NAOI, 1)
#endif
C
C The transformed matrix is DISZMO <= DISZAO. After copying we
C can zero out the DISZMO from AOINT so that the transformed ones can
C be written to it as they get formed.
         
                CALL ZERO(MUP_TMP_KEEP, NAOI*NTRNI)
C
                IOFF_MO = INDOCC(IRREP,ISPIN)

#ifdef _DEBUG_LVL0
         Print*, "The MO vectors, occ-first, vrt-second"
         Call output(EVEC(IOFF_MO), 1, NAOI, 1, NTRNI, NAOI, NTRNI, 1)
#endif
C
                DO IFIRST=NSTARTI, NENDI
C
                   CALL XGEMM('N','N',NAOI,NTRNI,NAOI,ONE,
     &                        BUF1,NAOI,EVEC(IOFF_MO),NAOI,
     &                        AZERO,MUP_TMP_KEEP,NAOI)
                END DO
#ifdef _DEBUG_LVL0
         Print*, "The 1st index transformed AO ints"
         Call output(MUP_TMP_KEEP, 1, NAOI, 1, NTRNI, NAOI, NTRNI, 1)
#endif
C 
C Copy the DISZTRAO elements to BUF1. Eventually we need to copy this
C back to AOINT at the moment lets do it in two steps. First copy 
C it to BUF1 (starting address is nbasis*nbasis + 1 since the first
C nbasis*nbasis elements is left for temporary space for the 
C transformation. In principle this tmp space shoud only be 
C (Max(NBAS(IRREP;IRREP=1,8)))**2, but the difference is two small
C to worry too much). 
C
                IOFF_MUP_LAMSIG = NBASIS*NBASIS +IBGN_MUP_LAMSIG + 1

                DO INDEX = 1, NTRNI
                   CALL DCOPY(NAOI, MUP_TMP_KEEP(NAOI*(INDEX-1)+1), 1, 
     &                        BUF1(IOFF_MUP_LAMSIG), 1)

#ifdef _DEBUG_LVLM1
         Write(6,*)
         Print*, "The 1st index transformed (mu,p,lam,sig) AO ints"
         Write(*, '(a, 3(2X, I4))'), "The starting address for the 
     &   KL pair and MO index", KLPAIRS, INDEX, IOFF_MUP_LAMSIG
         Call output(BUF1(IOFF_MUP_LAMSIG), 1, NAOI, 1, NTRNI, 
     &               NAOI, NTRNI, 1)
#endif
                  IOFF_MUP_LAMSIG  = IOFF_MUP_LAMSIG  + 
     &                               (KLPAIRS_READ-1)*NAOI + NAOI
                END DO
C     
                IOFF_MUNU_LAMSIG = IOFF_MUNU_LAMSIG + DISZAO
                IBGN_MUP_LAMSIG  = IBGN_MUP_LAMSIG + NAOI
C
             END DO
C
C Endif for nendi .ge. nstarti
C
         END IF
C
C end do for irrep
C
       END DO
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Print*, "The transformed integrals in trnsform_.F"
      Ioff = NBASIS*NBASIS
      Do Irrep = 1, Nirrep
         Nao = NBAS(IrreP)
         NTRNI = NFirst(IrreP)
         KLPAIRS = NAO*(NAO+1)/2
         IJPAIRS = NAO*NTRNI
         print*, IJPAIRS
         Do KL = 1, KLPAIRS

            Ioff = NBASIS*NBASIS + (KL -1)*IJPAIRS 
            Print*, "The KL Pair = ", KL, Ioff
            Write(*,'(4(1X,F12.7)))'), (BUF1(Ioff + I), I=1, IJPAIRS)

         End Do
      End Do
#endif
C
       RETURN
       END

