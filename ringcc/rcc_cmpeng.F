      SUBROUTINE RCC_CMPENG(ICORE,MAXCOR,NLIST2,NLIST1,ECORR,
     &                      ETOT,ETOTT2,IUHF,IPRINT)
      IMPLICIT INTEGER (A-Z)
      LOGICAL RHF
      CHARACTER*2 SPCASE(3)
      DOUBLE PRECISION E,ETOT,ECORR(3),ESPIN,ET2,ETOTT2,
     &                 ESING,SDOT,E_flip,E_no_flip
      DOUBLE PRECISION FACTOR_AA,FACTOR_AB
      DOUBLE PRECISION E_S, E_T
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/ IFLAGS(100)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD

      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD

      COMMON /NHFREF/NONHF
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C
      DIMENSION I0T(2),I0F(2)
C
      EQUIVALENCE (IFLAGS(2),METHOD)
C
      DATA SPCASE /'AA','BB','AB'/

      MXCOR=MAXCOR
      MXCOR=MAXCOR - MOD(MAXCOR,2)
      IFREE=1

      ETOT=0.D0
      ETOTT2=0.D0
      RHF = .FALSE.
      IF (IUHF .EQ.0) RHF = .TRUE.

c Natural factor of CCD AAAA and BBBB diagrams (FACTOR_AA) are 1/4 but 
c since for  antisymmetric cases we only store A<B;I<J the overall factor
c that is used in the code is 1.0. For symmetric case we have AB;IJ and 
c the factor is 1/4 for rCCD  and 1/2 for drCCD (based on plasmon formula(?)). 
C For the ABAB block for CC based methods, the factor (FACTOR_AB) is one and
c that make sense. However, In order to match my RPA code (and also the 
c litrature) I need factor 1.0 for DRCCD and 1/2 for RCCD. This is bit 
c puzzling. It seems both RCC and DRCC needs an extra factor half for AB 
c block (instead of 1.0 as for CC based methods) and I have already 
c included extra factor 2.0 in front (I am not sure about this!!!).
c See below for the comments about the factor for RHF spin-adapted cases.

      IF (DRCCD) THEN
         FACTOR_AA=0.50D0
         FACTOR_AB=1.00D0
      ELSE IF (RCCD) THEN
         FACTOR_AA=0.25D0
         FACTOR_AB=0.50D0
      ENDIF 

      IF (.NOT. RHF) THEN

      DO 10 ISPIN=1,IUHF+1
       LISTT=43+ISPIN
       ESPIN=0.D0
C
C  THE T2(IJ,AB) <IJ|AB> CONTRIBUTION TO THE ENERGY (AA AND BB SPIN
C  CASES
C
       DO 100 IRREP=1,NIRREP
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
        IF(MIN(NUMSYT,DISSYT).NE.0) THEN
        I001=IFREE
        I002=I001+IINTFP*NUMSYT*DISSYT
        I003=I002+IINTFP*NUMSYT*DISSYT
        I004=I003+NUMSYT
        IF(I004.LT.MXCOR) THEN
        IF (DRCCD) THEN
         CALL DRCC_TENER(NLIST2,ET2,E,NUMSYT,DISSYT,ICORE(I001),
     &                  ICORE(I002),ICORE(I0T(ISPIN)),ICORE(I0T(ISPIN)),
     &                  ISPIN,TAU,IRREP,POP(1,ISPIN),POP(1,ISPIN),
     &                  VRT(1,ISPIN),VRT(1,ISPIN),ICORE(I003))

         ELSE IF (RCCD) THEN
         LISTW =13+ISPIN
         DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
         NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
         CALL GETLST(ICORE(I002),1,NUMSYW,1,IRREP,LISTW)
C W(A<B;I<J) -> W(A<B;IJ)
         CALL SYMEXP(IRREP,POP(1,ISPIN),DISSYW,ICORE(I002))
C W(A<B;IJ) -> W(AB,IJ)
         CALL SYMEXP2(IRREP,VRT(1,ISPIN),DISSYT,DISSYW,
     +                NUMSYT,ICORE(I002),ICORE(I002))
         
         CALL RCC_TENER(NLIST2,ET2,E,NUMSYT,DISSYT,ICORE(I001),
     &                  ICORE(I002),ICORE(I0T(ISPIN)),ICORE(I0T(ISPIN)),
     &                  ISPIN,TAU,IRREP,POP(1,ISPIN),POP(1,ISPIN),
     &                  VRT(1,ISPIN),VRT(1,ISPIN),ICORE(I003))
         ENDIF 
        ELSE
         CALL INSMEM('CMPENG',I004,MXCOR)
        ENDIF

        ETOTT2 = ETOTT2+ FACTOR_AA*ET2
        ETOT   = ETOT  + FACTOR_AA*E
        ESPIN  = ESPIN + FACTOR_AA*E
        ENDIF
100    CONTINUE

       ECORR(ISPIN)=ESPIN
       IF(IPRINT.NE.0)WRITE(*,80)SPCASE(ISPIN),ESPIN
80     FORMAT(T3,' The ',A2,' contribution to the correlation ',
     &        'energy is: ',F12.7,' a.u.')
10    CONTINUE
      ENDIF

      IF (RHF .AND. RCCD) THEN
          ESPIN = 0.0D0
          Write(6,*) "Entering RCC_RHF_CMPENEG"
          CALL RCC_RHF_CMPENEG(ICORE,MAXCOR,IUHF,NLIST2,E_S,E_T)

          IF(IPRINT.NE.0)WRITE(*,91) E_S
          IF(IPRINT.NE.0)WRITE(*,92) E_T

          E_flip     = 0.25D0*(E_t + E_s)
          E_no_flip  = 0.75D0*E_t + 0.25D0*E_s 
          ETOT  = E_flip
          ESPIN = E_flip

 91       FORMAT(T3,' The singlet spin adapted contribution to',
     &           ' the correlation energy is: ',F12.7,' a.u.')
 92       FORMAT(T3,' The triplet spin adapted contribution to',
     &           ' the correlation energy is: ',F12.7,' a.u.')
          ECORR(3)=ESPIN
          Write(6,"(a,a,F15.12)") " The rCCD correlation energy" ,
     &                            " with spin-flip excitations    is ", 
     &                              E_flip

          Write(6,"(a,a,F15.12)") " The rCCD correlation energy" ,
     &                            " without spin-flip excitations is ", 
     &                              E_no_flip
          Write(6,*)
          Write(6,"(a,a)") " The rCCD energy printed in the table",
     &           " correspond to rCCD with spin-flip"
          Write(6,"(a,a,a)") " excitations. This energy is identical",
     &            " to the energy computed with the"
          Write(6,"(a)") " spin-orbital formulation."
          Write(6,*)
          RETURN
      ENDIF 

C THE T2(Ij,Ab) <Ij|Ab> CONTRIBUTION TO THE ENERGY (SPIN CASE AB)
C The external multiplication factor is 1.0 here.
C For rCCD this is the singlet-spin adapated part.

      ESPIN=0.D0
      DO 200 IRREP=1,NIRREP
       DISSYT=IRPDPD(IRREP,ISYTYP(1,46))
       NUMSYT=IRPDPD(IRREP,ISYTYP(2,46))
       IF(MIN(NUMSYT,DISSYT).NE.0) THEN
        I001=IFREE
        I002=I001+IINTFP*NUMSYT*DISSYT
        I003=I002+IINTFP*NUMSYT*DISSYT
        I004=I003+NUMSYT
        IF(I004.LT.MXCOR) THEN
         CALL RCC_TENER(NLIST2,ET2,E,NUMSYT,DISSYT,ICORE(I001),
     &                  ICORE(I002),ICORE(I0T(1)),ICORE(I0T(2)),3,TAU,
     &                  IRREP,POP(1,1),POP(1,2),VRT(1,1),VRT(1,2),
     &                  ICORE(I003))
        ELSE
         CALL INSMEM('CMPENG',I004,MXCOR)
        ENDIF
        ETOTT2= ETOTT2 + FACTOR_AB*ET2
        ETOT  = ETOT   + FACTOR_AB*E
        ESPIN = ESPIN  + FACTOR_AB*E
       ENDIF
200   CONTINUE
      IF(IPRINT.NE.0)WRITE(*,80)SPCASE(3),ESPIN

c This multiplication is only for spin-adapted RHF. In otherwords,
c this accounts for the contributions from AAAA and BBBB. 

      IF (DRCCD .AND. RHF) THEN
         ESPIN = 2.0D0 * ESPIN
         ETOT  = 2.0D0 * ETOT
      ENDIF 

      ECORR(3)=ESPIN
      IF(IPRINT.NE.0)WRITE(*,81)ETOT
81    FORMAT(T3,' The total correlation energy is ',F15.12,' a.u.')
     
      RETURN
      END
