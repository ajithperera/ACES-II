      SUBROUTINE TRNSFORM_1INDX_PONCE(EVEC,AOINT,BUF1,MUP_TMP_KEEP,
     &                                KLPAIRS_4IRREP,NFIRST,NSTART,
     &                                NEND,ISPIN,NBAS,NMO,NMOS_HDLE,
     &                                NPASS,NBASIS,IOFF_IRREP_4AO,
     &                                MAX_AOS_IRREP, IRREPX)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL RHF,DOALL

      INTEGER DIRPRD,DISZAO,DISZMO,DISZTRAO
      DOUBLE PRECISION MUP_TMP_KEEP(NBASIS*NBASIS)
C
      DIMENSION EVEC(1),AOINT(1),BUF1(1)
      DIMENSION NBAS(8),NMO(8),NFIRST(8),NSTART(8),NEND(8),
     &          KLPAIRS_4IRREP(8),IOFF_IRREP_4AO(8)
C
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/AOOFST/INDOCC(8,2)
      COMMON/SYMINF/NDUMMY,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
      DATA AZERO,ONE,FOUR/0.D0, 1.D0, 4.0D0/
C
      Indx(i, j) = j+(i*(i-1))/2
#ifdef _DEBUG_LVL0
      Print*, "The npass and irrepx", npass, irrepx
      Print*, "The nstartx and nendx arrays"
      Print*, (nStart(I), i=1, 8)
      Print*, (nEnd(I), i=1, 8)
      Print*, ((indocc(i,j), i=1,8), j=1,2)
      Write(6,*)
      Print*, "The Nbas Nfirst and Isize arrays"
      Print*, (Nbas(I), i=1, 8)
      Print*, (Ioff_Irrep_4AO(I), i=1, 8)
      Print*, (NFirst(I), i=1, 8)
      Write(6,*)
#endif
C
       IOFF_MUNU_LAMSIG = 1
        IOFF_MUP_LAMSIG = 1
           IOFF_4AOINTS = 1
           IOFF_4DAOINT = MAX_AOS_IRREP*MAX_AOS_IRREP + 1
C
C
       DO IRREP=1,NIRREP
C
          IF (IRREP .EQ. IRREPX) THEN
              NAOI = NBAS(IRREP)
              NMOI = NMO(IRREP)
            DISZAO = NAOI*(NAOI+1)/2
            DISZMO = NAOI*(NMOI+1)/2
C
C NAOI and NMOI are the # of AOs and and the total coorelated MOs
C in a given irrep. The DISZAO and DISZMO are the lengths of the 
C two index AO or AO-MO triangular arrays.
C
            NTRNI    = NMOS_HDLE
            DISZTRAO = NAOI*NMOS_HDLE
C
            NSTARTI  = NSTART(IRREP)
            NENDI    = NEND(IRREP)
C  
            KLPAIRS_READ      = KLPAIRS_4IRREP(IRREP)
            IOFF_MUNU_LAMSIG  = IOFF_IRREP_4AO(IRREP) + 1
C
C NFIRSTI keep the range of MOS that need to be transformed. It is 
C identical to NMOI when no orbitals are dropped. Depending on 
C whether NOABCD set NFIRSTI has the full or occupide MOs set. 
C
C NSTARTI and NENDI are the starting and ending addresses of the range
C of MOs that is C being transformed. NSTARTI is the starting address 
C of the range if the entire set can not be handled in a one pass.
C
           IBGN_MUP_LAMSIG = 0
           IF (NENDI .GE. NMOS_HDLE) Then 
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
C
CSSS             IBGN_MUP_LAMSIG = (NPASS - 1)*NAOI*NMOS_HDLE*
CSSS     &                          KLPAIRS_READ

             DO KLPAIRS = 1, KLPAIRS_READ
                CALL EXPND2(AOINT(IOFF_MUNU_LAMSIG), BUF1, NAOI)
C
#ifdef _DEBUG_LVL0
         Print*, "The AO ints for a KL PAIRS:", "KL = ", KLPAIRS
         If ((npass .eq. 1 .or. npass .eq.2) .and. klpairs .eq.1) 
     &       Call output(BUF1, 1, NAOI, 1, NAOI, NAOI, NAOI, 1)
#endif
C
C The transformed matrix is DISZMO <= DISZAO. After copying we
C can zero out the DISZMO from AOINT so that the transformed ones can
C be written to it as they get formed.
         
                CALL ZERO(MUP_TMP_KEEP, NAOI*NTRNI)
C
                IOFF_MO = INDOCC(IRREP,ISPIN) + (NSTART(IRREP) - 1)
     &                    *NAOI
C
#ifdef _DEBUG_LVL0
         Print*, "The MO vectors, occ-first, vrt-second"
         If (klpairs .eq. 1 .and. npass .eq. 1) 
     &   Call output(EVEC(IOFF_MO), 1, NAOI, 1, NTRNI, NAOI, NTRNI, 1)
#endif
C
CSSS                DO IFIRST=NSTARTI, NENDI
C
                   CALL XGEMM('N','N',NAOI,NTRNI,NAOI,ONE,
     &                        BUF1,NAOI,EVEC(IOFF_MO),NAOI,
     &                        AZERO,MUP_TMP_KEEP,NAOI)
CSSS                END DO
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
                IOFF_MUP_LAMSIG = IBGN_MUP_LAMSIG + IOFF_4DAOINT
C
                DO INDEX = 1, NTRNI
#ifdef _DEBUG_LVLM1
         Write(6,*)
         Print*, "The 1st index transformed (mu,p,lam,sig) AO ints"
         Write(*, '(a, 3(2X, I10))'), "The starting address for the
     &   KL pair and MO index", KLPAIRS, INDEX, IOFF_MUP_LAMSIG
         Call output(BUF1(IOFF_MUP_LAMSIG), 1, NAOI, 1, NTRNI,
     &               NAOI, NTRNI, 1)
#endif
                   CALL DCOPY(NAOI, MUP_TMP_KEEP(NAOI*(INDEX-1)+1), 1, 
     &                        BUF1(IOFF_MUP_LAMSIG), 1)

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

           IOFF_4DAOINT = IOFF_4DAOINT + NTRNI*NAOI*KLPAIRS_READ
C
C End if for irrep=irrepx
C
         ENDIF

C
C End if for irrep
C
       END DO
C
C Copy the entire first index transformed integral to the location of
C the origianal AOINTs. This will makesure that the INTRN1 will
C procced correctly. The BUF1 location is uesd in INTRN1 as a scratch
C array.
C
       IBGN_AO_MO   = 0
       IOFF_4AOINTS = 0
       DO IRREP = 1, NIRREP
C
          IF (IRREP .EQ. IRREPX) Then
             NAO     = NBAS(IRREP)
             KLPAIRS = NAO*(NAO+1)/2
             NTRNI   = NMOS_HDLE
             IJPAIRS = NAO*NTRNI
C
             IF (NAO .GT. 0) THEN
                 DO KL = 1, KLPAIRS 
                    IOFF_MO = IBGN_AO_MO + MAX_AOS_IRREP*
     &                        MAX_AOS_IRREP + (KL -1)*IJPAIRS
                    IOFF_AO = IBGN_AO_MO + (KL -1)*IJPAIRS + 1
                    CALL DCOPY(IJPAIRS, BUF1(IOFF_MO+1), 1,
     &                         AOINT(Ioff_AO), 1)
                 ENDDO
                 IBGN_AO_MO = KLPAIRS*IJPAIRS
             END IF
            IOFF_4AOINTS = IOFF_4AOINTS + NAO*NTRNI*
     &                                    KLPAIRS_4IRREP(IRREP)
            IBGN_AO_MO   = IOFF_4AOINTS 
          END IF 
C
       ENDDO
C    
#ifdef _DEBUG_LVL0
      Write(6,*)     
      Print*, nbasis
      Print*, "The transformed integrals in trnsform_.F"
CSSS      Call checksum2("AOLAD1", W, KLPAIRS*IJPAIRS)
      IBGN_AO_MO   = 0
      IOFF_4AOINTS = 0
      Do Irrep = 1, Nirrep
         If (irrep .eq. irrepx) Then
            Nao = NBAS(IrreP)
            NTRNI =  NMOS_HDLE
            KLPAIRS = NAO*(NAO+1)/2
            IJPAIRS = NAO*NTRNI
#ifdef _DEBUG_LVL0
            IF (NAO .GT. 0) THEN
               Do KL = 1, KLPAIRS
                  Ioff = IBGN_AO_MO + (KL -1)*IJPAIRS
                  Print*, "The KL Pair = ", KL, Ioff
                  Write(*,'(4(1X,F12.7)))'), (AOINT(Ioff + I), 
     &                                        I=1, IJPAIRS)
CSSS            Write(*,'(4(1X,F12.7)))'), (Buf1(Ioff + I), I=1, IJPAIRS)
CSSS            Call checksum2("AOLAD1", BUf1(Ioff), IJPAIRS)
              End Do
#endif
                IBGN_AO_MO = KLPAIRS*IJPAIRS
            ENDIF  
           IOFF_4AOINTS = IOFF_4AOINTS + NAO*NTRNI*
     &                                   KLPAIRS_4IRREP(IRREP)
           IBGN_AO_MO   = IOFF_4AOINTS
         ENd IF
C 
      End Do
#endif
C
       RETURN
       END

