      SUBROUTINE ROHFPT(ICORE,MAXCOR,IUHF)
C
C DRIVER FOR ROHF-MBPT CALCULATIONS.  
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ENERGY,ET1,ET2,ETOT
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &        QCISD,ROHF4,ITRFLG,TRIPNI,TRIPNI1,TRIPIT,T3STOR,UCC,
     &        BRUECK
      DIMENSION ICORE(MAXCOR),ECORR(3)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC
      COMMON /TRIPLES/ TRIPNI,TRIPNI1,TRIPIT,T3STOR
      COMMON /ENERGY/ ENERGY(500,2),IXTRLE(500)
      COMMON /ROHF/ ROHF4,ITRFLG
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /BRUECKNER/BRUECK
#include "syminf.com" /* for nirrep */
C
C-----------------------------------------------------------------------
CJDW 3/98. Attempts have been made to convert this routine and its
C          dependents from being only an ROHF-MBPT open-shell code to
C          a non-HF MBPT code for REF=RHF,UHF.
C-----------------------------------------------------------------------
C
      WRITE(*,1000)
1000  FORMAT(T3,' Calculating MBPT energy for ROHF reference function.')
C
C  SINCE THIS IS NON-HF, WE NEED TO DEAL WITH T1, HENCE SING1 IS TRUE.
C  ITRFLG IS A FLAG WHICH TELLS CERTAIN ROUTINES WHETHER TO PICK UP
C  THEIR NORMAL LISTS (TRUE) OR TO USE DIFFERENT LISTS FOR ITERATIVE
C  SDQ- AND SDTQ-MBPT[4].
C
      SING1=.TRUE.
C
C INITIALIZATION ROUTINES
C
      CALL INITSN(ICORE,MAXCOR,IUHF)
C
C  IF NON-SEMICANONICAL ORBITALS ARE USED, THEN WE NEED TO ITERATE
C  IN ORDER TO GET T2[1] AND T1[1].
C
      IF(IFLAGS(2).LE.4.AND.IFLAGS(39).EQ.0)THEN
       CALL T1ITER2(ICORE,MAXCOR,IUHF,ET1)
       CALL T2ITER2(ICORE,MAXCOR,IUHF,ET2)
      ENDIF
      CALL ROHFENG(ICORE,MAXCOR,ENERGY(1,1),ENERGY(1,2),43,0,2,IUHF)
C
      IF(IUHF.NE.0) THEN
       WRITE(*,4000)
 4000  FORMAT(/,T3,'MBPT(2) spin multiplicity: ',/)
       CALL S2PRJIT(ICORE,MAXCOR,IUHF,SING1,2)
       WRITE(*,4001)
 4001  FORMAT(/)
      ENDIF
C
C THE FOLLOWING CODE IS EXPERIMENTAL AND FOR A MBPT BASED
C ON BRUECKNER ORBITALS
C
      IF(BRUECK) THEN
C
C RESORT T2 AMPLITUDES
C
       CALL RNABIJ(ICORE,MAXCOR,IUHF,'T')
C
C LOOP OVER SPIN CASES  
C
       DO 10 ISPIN=1,IUHF+1
C
C COPY f(m,e) TO F(m,e)
C
        CALL GETLST(ICORE,1,1,1,ISPIN+2,93)
        CALL PUTLST(ICORE,1,1,1,ISPIN,93)
C
10     CONTINUE
C
C LOOP OVER SPIN CASES  
C
       DO 15 ISPIN=1,IUHF+1
C
C CALCULATE F*T2 INTO T1
C
        CALL FMECONT(ICORE,MAXCOR,IUHF,ISPIN)
C
15     CONTINUE
C
C CALCULATE W*T2 INTO T1
C
       CALL E3S(ICORE,MAXCOR,IUHF,EDUMMY)
C
C DENOMINATOR WEIGHT T1 INCREMENTS
C
       CALL NEWT1(ICORE,MAXCOR,IUHF,1)
C
C COPY TOTAL T1 AMPLITUDES TO LISTS 1,90 AND 2,90 FOR BRUECKNER UPDATE
C
       DO 120 ISPIN=1,IUHF+1
        CALL GETLST(ICORE,1,1,1,ISPIN+2,90)
        CALL PUTLST(ICORE,1,1,1,ISPIN,90)
120    CONTINUE
C
      ENDIF

      IF(IFLAGS(2).EQ.1) RETURN
C
      CALL SETLST(ICORE,MAXCOR,IUHF)
      CALL RNABIJ(ICORE,MAXCOR,IUHF,'T')
C
C GENERATE W AND F INTERMEDIATES FOR MBPT(2), MBPT(3) AND LINEAR MBPT(4).
C THESE ARE OF COURSE ONLY THE LEADING TERMS (LINEAR CONTRIBUTIONS)
C
      CALL GENINT(ICORE,MAXCOR,IUHF,3)
C
C THIRD ORDER STUFF.
C
C
C 1. DO LINEAR INCLUSION OF F-T2[1] TERM INTO T1[2].
C
      CALL FMECONT(ICORE,MAXCOR,IUHF,1)
      IF(IUHF.NE.0) CALL FMECONT(ICORE,MAXCOR,IUHF,2)
C
C 2. DO LINEAR INCLUSION OF W-T1[1] INTO T2[2]
C
      CALL T1INT2A(ICORE,MAXCOR,IUHF)
      CALL T1INT2B(ICORE,MAXCOR,IUHF)
C
C 3. DO LINEAR INCLUSION OF W-T1[1] INTO T1[2]
C
      CALL T1INT1(ICORE,MAXCOR,IUHF,1)
      IF(IUHF.NE.0) CALL T1INT1(ICORE,MAXCOR,IUHF,2)
C
C 4. DO LINEAR INCLUSION OF W-T2[1] TO T2[2]
C
      CALL DRE3EN(ICORE,MAXCOR,IUHF,0)
C
C 5. DO LINEAR INCLUSION OF T2[1] TO T1[2].
C
      CALL E3S(ICORE,MAXCOR,IUHF,EDUMMY)
C
C 6. DENOMINATOR WEIGHT THESE INCREMENTS TO FORM THE REAL T1[2]
C    AND T2[2].  AFTER THIS:
C
C           QUANTITY                   LOCATION 
C            T1[1]                      1,90 2,90
C            T1[2]                      3,90 4,90
C            T2[1]                     LISTS 44-46
C            T2[2]                     LISTS 61-63
C
C
      IF(IFLAGS(39).EQ.1) THEN
        CALL NEWT1(ICORE,MAXCOR,IUHF,0)
        CALL NEWT2(ICORE,MAXCOR,IUHF)
      ELSE
        CALL T2ITER3(ICORE,MAXCOR,IUHF,3)
        CALL T1ITER3(ICORE,MAXCOR,IUHF)
      ENDIF
C
C COMPUTE THE THIRD ORDER CONTRIBUTION TO THE ENERGY FROM
C
C     E(3) =   SUM <ab||ij> TAU[2](ij,ab) + SUM f(ai) T1[2](a,i)
C            a<b,c<d                       a,i
C
C  (NOTE THAT THE SECOND ORDER TAU IS MADE UP OF SECOND ORDER T2 
C   AND THE FIRST ORDER T1)
C
      CALL ROHFENG(ICORE,MAXCOR,ENERGY(2,1),ENERGY(2,2),60,2,3,IUHF)
      IF(IUHF.NE.0)THEN
       WRITE(*,4010)
 4010  FORMAT(/,T3,'MBPT(3) spin multiplicity: ',/)
       CALL S2PRJIT(ICORE,MAXCOR,IUHF,SING1,3)
       WRITE(*,4001)
      ENDIF
C
      IF(IFLAGS(2).EQ.2) RETURN
C
C FOURTH ORDER !
C
      ROHF4=.TRUE.
C
C
C FIRST DO THE LINEAR FOURTH-ORDER TERM AND THE T1[1]*T1[2]*W
C  CONTRIBUTION IF THIS CALCULATION USES SEMICANONICAL ORBITALS.
C
      IF(IFLAGS(39).EQ.1) THEN
        CALL E4LNHF(ICORE,MAXCOR,IUHF,ENERGY(3,1),ENERGY(4,2))
      ELSE
C
C  Move the T1[2] and T2[2] pieces to new lists.
C  After this move we'll have:
C
C       Iterate			List
C	 T1[1]                  1,90;  2,90
C        T1[2]                  9,90; 10,90
C        T2[1]                  44-46
C        T2[2]                  94-96
C
C        T1[3]                  3,90; 4,90
C        T2[3]                  61-63
C
C
C  The 94-96 lists were created and zeroed in t2iter3.
C  We need to create and zero the 9,90 and 10,90 lists.
C
        NSIZE=ISYMSZ(ISYTYP(1,61),ISYTYP(2,61))
        CALL GETALL(ICORE,NSIZE,1,61)
        CALL PUTALL(ICORE,NSIZE,1,94)
        IF(IUHF.NE.0)THEN
         NSIZE=ISYMSZ(ISYTYP(1,62),ISYTYP(2,62))
         CALL GETALL(ICORE,NSIZE,1,62)
         CALL PUTALL(ICORE,NSIZE,1,95)
        ENDIF
        NSIZE=ISYMSZ(ISYTYP(1,63),ISYTYP(2,63))
        CALL GETALL(ICORE,NSIZE,1,63)
        CALL PUTALL(ICORE,NSIZE,1,96)
C
        CALL UPDMOI(1,NT(1),9,90,0,0)
        IF(IUHF.NE.0) CALL UPDMOI(1,NT(2),10,90,0,0)
        CALL GETLST(ICORE,1,1,1,3,90)
        CALL PUTLST(ICORE,1,1,1,9,90)
        IF(IUHF.NE.0)THEN
         CALL GETLST(ICORE,1,1,1,4,90)
         CALL PUTLST(ICORE,1,1,1,10,90)
        ENDIF
      ENDIF
C
      CALL ZEROLIST(ICORE,MAXCOR,61)
      IF(IUHF.NE.0) CALL ZEROLIST(ICORE,MAXCOR,62)
      CALL ZEROLIST(ICORE,MAXCOR,63)
      IF(IUHF.NE.0) CALL ZEROLIST(ICORE,MAXCOR,40)
      CALL ZEROLIST(ICORE,MAXCOR,41)
      CALL ZEROLIST(ICORE,MAXCOR,42)
      CALL ZEROLIST(ICORE,MAXCOR,43)
C
C  This call effectively zeros out lists 3,90 and 4,90.
C
      CALL INITSN(ICORE,MAXCOR,IUHF)
C
C FIRST CALL GENINT WITH A FLAG WHICH COMPUTES ONLY THE SECOND-ORDER
C  PARTS OF EACH INTERMEDIATE
C
      CALL GENINT(ICORE,MAXCOR,IUHF,1)
C
C CALL FMECONT TO GET THE F(ME) CONTRIBUTION TO T1(3)
C
C  This piece is:
C
C      SUM T(im,ae)[1]*F(me)[2]
C       me
C
      CALL FMECONT(ICORE,MAXCOR,IUHF,1)
      IF(IUHF.NE.0) CALL FMECONT(ICORE,MAXCOR,IUHF,2)
C
C NOW COPY THE f(AI) [FOCK MATRIX] ONTO THE F(ME) INTERMEDIATE LIST 
C  [FIRST-ORDER F(ME)]
C
       CALL GETLST(ICORE,1,1,1,3,93)
       CALL PUTLST(ICORE,1,1,1,1,93)
       IF(IUHF.NE.0)THEN
        CALL GETLST(ICORE,1,1,1,4,93)
        CALL PUTLST(ICORE,1,1,1,2,93)
       ENDIF
C
C NOW CALL FEACONT AND FMICONT TO GET THE REST OF T1(3) AND PARTS OF 
C  T2(3)
C
C  These pieces are (feacont):
C
C     P_(AB) SUM T(ij,ae)[1]*F(be)[2]
C             e
C
C
C    - 1/2 P_(AB) SUM T(ij,ae)[1] * SUM T(m,b)[1]*F(me)[1]
C             e                 m
C
C
C     SUM T(i,e)[1]*F(ae)[2]
C      e
C
C
C  These pieces are (fmicont):
C
C
C    -P_(IJ) SUM T(im,ef)[1]*F(mj)[2]
C             m
C
C    - 1/2 P_(ij) SUM T(im,ef)[1] * SUM T(j,e)[1]*F(me)[1]
C                  m                 e
C
C
C     SUM T(m,a)[1]*F(mi)[2]
C      m
C
      CALL FEACONT(ICORE,MAXCOR,IUHF)
      CALL FMICONT(ICORE,MAXCOR,IUHF)
C
C  If this calculation uses standard orbitals, we need to deal with
C  the terms outlined below.
C
      IF(IFLAGS(39).EQ.0) THEN
C
C  The DRLAD calls in DRE3EN handle:
C
C     1/2 SUM tau(mn,ab)[1]*W(mn,ij)[2]
C          mn
C
C  This term is handled computationally with a scale factor, since
C  it is equal to another quantity already calculated.  It is never
C  calculated explicitly.
C
C        1/2 SUM tau(ij,ef)[1]*W(ab,ef)[2]
C             ef
C
C  The T12INT2 call in DRE3EN handles:
C
C     -P_(ij)P_(ab) SUM T(i,e)[1]*T(m,a)[1]*<mb||ej>
C                    me
C
C
C  The DRRNG call in DRE3EN handles:
C
C     P_(ij)P_(ab) SUM T(im,ae)[1]*W(mb,ej)[2]
C                   me
C
C
C
        CALL DRE3EN(ICORE,MAXCOR,IUHF,1)
C
C
C  Now we handled all the terms which require second order intermediates.
C  We now generate the first order W intermediates ONLY.
C
        CALL GENINT(ICORE,MAXCOR,IUHF,0)
C
C  This call resorts the T2[2] elements from lists 94-96 to lists
C  34-39.
C
        CALL RNABIJ(ICORE,MAXCOR,IUHF,'C')
C
C  Clean up the stuff left on the increment lists from the DRRNG call
C  in DRE3EN above.
C
        IF(IUHF.NE.0) CALL ZEROLIST(ICORE,MAXCOR,40)
        CALL ZEROLIST(ICORE,MAXCOR,41)
        CALL ZEROLIST(ICORE,MAXCOR,42)
        CALL ZEROLIST(ICORE,MAXCOR,43)
C
C  The DRLAD calls in DRE3EN handle:
C
C     1/2 SUM tau(mn,ab)[2]*W(mn,ij)[1]
C          mn
C
C
C     1/2 SUM tau(ij,ef)[2]*W(ab,ef)[1]
C          ef
C
C
C  The DRRNG call in DRE3EN handles:
C
C     P_(ij)P_(ab) SUM T(im,ae)[2]*W(mb,ej)[1]
C                   me
C
C
C  ITRFLG=false tells DRE3EN not to include the t1**2 piece;
C  this was accounted for in the immediately previous call to DRE3EN.
C  It also tells DRLAD to use the correct lists for formation of tau
C  and to use the integrals directly (1st order W).
C
C
         ITRFLG=.FALSE.
         CALL DRE3EN(ICORE,MAXCOR,IUHF,0)
C
C  The following two calls calculate
C
C     -P_(ab) SUM T(m,a)[2]*<mb||ij>
C              m
C
C
C      P_(ij) SUM T(i,e)[2]*<ab||ej>
C              e
C
C
        CALL T1INT2A(ICORE,MAXCOR,IUHF)
        CALL T1INT2B(ICORE,MAXCOR,IUHF)

C
C   The routine FMECONT calculate
C
C       SUM T(im,ae)[2]*F(me)[1]
C        me
C
        CALL FMECONT(ICORE,MAXCOR,IUHF,1)
        IF(IUHF.NE.0) CALL FMECONT(ICORE,MAXCOR,IUHF,2)
C
C  The routine T1INT1 calculates
C
C    - SUM T(n,f)[2]*<na||if>
C       nf
C
        CALL T1INT1(ICORE,MAXCOR,IUHF,1)
        IF(IUHF.NE.0) CALL T1INT1(ICORE,MAXCOR,IUHF,2)
C
C  The routine E3S calculates
C
C     -1/2 SUM T(im,ef)[2]*<ma||ef>
C          mef
C
C
C     -1/2 SUM T(mn,ae)[2]*<nm||ei>
C          men
C
        IFLAGS(39)=99
CJDW 3/98  Fudge to force routines called by E3S to read right lists.
        IFL11OLD=IFLAGS(11)
        IFLAGS(11)=2
        CALL E3S(ICORE,MAXCOR,IUHF,EDUMMY)
        IFLAGS(39)=0
        IFLAGS(11)=IFL11OLD
C
C  We now have the full T1[3] amplitudes on 3,90 and 4,90
C  We now have the full T2[3] amplitudes on 61-63
C
      ENDIF
C
      IF(IFLAGS(39).EQ.1) THEN
C
C  For semicanonical orbitals:
C
C THE "THIRD-ORDER" T1 AMPLITUDES ARE NOW STORED ON LISTS 3,90 AND 4,90.
C
C COMPUTE MOST PARTS OF T2(3)
C
        CALL DRE3EN(ICORE,MAXCOR,IUHF,1)
        iwalt=999
        if(iwalt.eq.999) then
C
C  Generate first order W intermediates (just the integrals)
C
        CALL GENINT(ICORE,MAXCOR,IUHF,0)
C
C  T(m,a)[1]*T(n,b)[1]*W[1]  contraction, performed in two steps.
C
c        IMODE=1
        call aces_io_remove(51,'GAMLAM')
        IMODE=0
        IF(IUHF.NE.0)THEN
         CALL INIPCK(1,ISYTYP(1,7),ISYTYP(2,7),107,IMODE,0,1)
         CALL INIPCK(1,ISYTYP(1,8),ISYTYP(2,8),108,IMODE,0,1)
         do iFam = 107, 108
         do iGrp = 1, nirrep
            call aces_list_memset(iGrp,iFam,0)
         end do
         end do
        ENDIF
        CALL INIPCK(1,ISYTYP(1,10),ISYTYP(2,10),110,IMODE,0,1)
        do iGrp = 1, nirrep
           call aces_list_memset(iGrp,110,0)
        end do
        IF(IUHF.NE.0)THEN
         CALL W4AA1(ICORE,MAXCOR,IUHF,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
        ENDIF
        CALL W4AB2(ICORE,MAXCOR,IUHF,.FALSE.,.TRUE.,.FALSE.,.FALSE.)
        CALL YT1INT2A(ICORE,MAXCOR,IUHF,1)
C
        IF(IUHF.NE.0)THEN
         CALL INIPCK(1,ISYTYP(1,27),ISYTYP(2,27),127,IMODE,0,1)
         CALL INIPCK(1,ISYTYP(1,28),ISYTYP(2,28),128,IMODE,0,1)
         do iFam = 127, 128
         do iGrp = 1, nirrep
            call aces_list_memset(iGrp,iFam,0)
         end do
         end do
        ENDIF
        CALL INIPCK(1,ISYTYP(1,30),ISYTYP(2,30),130,IMODE,0,1)
        do iGrp = 1, nirrep
           call aces_list_memset(iGrp,130,0)
        end do
        IF(IUHF.NE.0)THEN
         CALL W5AA1(ICORE,MAXCOR,IUHF,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
        ENDIF
        CALL W5AB2(ICORE,MAXCOR,IUHF,.FALSE.,.TRUE.,.FALSE.,.FALSE.)
        CALL YT1INT2B(ICORE,MAXCOR,IUHF,1)
        endif
C
        CALL NEWT1(ICORE,MAXCOR,IUHF,0)
        CALL NEWT2(ICORE,MAXCOR,IUHF)
        CALL ROHFENG(ICORE,MAXCOR,ENERGY(4,1),ENERGY(3,2),60,2,4,IUHF)
        ENERGY(4,1)=ENERGY(4,1)+ENERGY(3,2)
        ENERGY(3,2)=0.D0
C
C  At this point, we're done with the lists on GAMLAM, so blow GAMLAM
C  away.  This will prevent any interference with gradient logic which
C  comes after ROHFPT.
C
C        CALL ACES_IO_REMOVE(51,'GAMLAM')
      ELSE
C
C  ROHFT1SQ does the W*T1[1]*T1[2] part.
C  Then we do the iterative parts of T1[3] and T2[3].
C
        CALL ROHFT1SQ(ICORE,MAXCOR,8,IUHF,ENERGY(3,2))
        CALL T1ITER3(ICORE,MAXCOR,IUHF)
        CALL T2ITER3(ICORE,MAXCOR,IUHF,4)
        CALL ROHFENG(ICORE,MAXCOR,ENERGY(5,1),ENERGY(5,2),60,2,4,IUHF)
        IF(IUHF.NE.0)THEN
         WRITE(*,4020)
 4020    FORMAT(/,T3,'SDQ-MBPT(4) spin multiplicity: ',/)
         CALL S2PRJIT(ICORE,MAXCOR,IUHF,SING1,4)
C        CALL S2PRJIT(ICORE,MAXCOR,IUHF,.false.,4)
         WRITE(*,4001)
        ENDIF
      ENDIF
C
C  If we are doing full 4th order, we need the triples.
C
      IF(TRIPNI) THEN
        CALL TRPS(ICORE,MAXCOR,IUHF,ENERGY(6,1))
        IF(IFLAGS(39).EQ.0) THEN
          WRITE(*,4021)
 4021     FORMAT(/,T3,'SDTQ-MBPT(4) spin multiplicity: ',/)
          CALL S2PRJIT(ICORE,MAXCOR,IUHF,SING1,5)
C          CALL S2PRJIT(ICORE,MAXCOR,IUHF,.false.,5)
          WRITE(*,4001)
        ELSE
          WRITE(*,4025)
 4025     FORMAT(/,T3,'@ROHFPT-I, Full third order wavefunction not ',
     &              'formed.',/,
     &           T14,'Multiplicity not computed.',/)
        ENDIF
      ENDIF
C
      RETURN
      END
