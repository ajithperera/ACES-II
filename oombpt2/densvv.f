       SUBROUTINE DENSVV(DVV,ICORE,MAXCOR,IUHF,SINGLES)
C
C
C THIS ROUTINE CALCULATES THE VIRTUAL-VIRTUAL BLOCK
C OF THE RELAXED DENSITY MATRIX IN CORRELATION METHODS
C
C THE FORMULAS ARE
C
C D(A,B) = 1/2 SUM M,N,E T[1](MN,AE) T[1](MN,BE)              (MBPT(2))
C
C D(A,B) = 1/2 SUM M,N,E T[1](MN,AE) T[1](MN,BE)
C
C          + SUM M T[1](M,A) T[1](M,B)                        (ROHF-MBPT(2))
C
C D(A,B) = 1/2 P(AB) SUM M,N,E T[1](MN,AE) (1/2 T[1](MN,BE)+T2(MN,BE))
C
C                                                             (MBPT(3))
C
C D(A,B) =  1/4 P(AB) SUM M,N,E T[1](IM,EF) (L[2](JM,EF)+T[2](JM,EF))
C
C           1/2 SUM M,N,E T[1](IM,EF) T[1](JM,EF))
C
C           1/2 P(AB) SUM M T[1](M,A) (T[1](I,E) + L[2](M,B)+ T[2](M,B))
C
C                                                         (ROHF-MBPT(3))
C
C D(A,B) = SUM M T[2](M,A) T[2](M,B)
C
C          + 1/2 SUM M,N,E {T[1](MN,AE)+T[2](MN,AE)}{T[2](MN,BE)+T[1](MN,BE)}
C
C          + 1/2 P(AB) SUM M,N,E T[1](MN,AE) ( T[3](MN,BE)
C
C                                            + 1/2 X(MN,BE))
C
C          + 1/12 SUM M,N,O,E,F T[2](MNO,AEF) T[2](MNO,BEF)     (MBPT(4))
C
C D(A,B) = 1/4 P(AB) SUM M,N,E T(MN,AE) LAMBDA(MN,BE)       (CCD)
C
C D(A,B) = 1/4 P(AB) SUM M,N,E T(MN,AE) LAMNDA(MN,BE)
C
C          + 1/2 P(AB) SUM N T(M,A) LAMBDA(M,B)             (CCSD, QCISD)
C
C D(A,B) = 1/4 P(AB) SUM M,N,E T(MN,AE) T(MN,BE)
C
C          + 1/2 P(AB) SUM N T(M,A) T(M,B)                    (UCC)
C
C
C THERE ARE THE FOLLOWING SPIN TYPES TO CONSIDER
C
C          D(AB)                T(MN,AE),....
C
C          AA                   AAAA, ABAB
C
C          BB                   BBBB, BABA
C
C THIS SUBROUTINE USES EXPLICITELY SYMMETRY
C
C THE TRIPLES CONTRIBUTION IS CALCULATED IN TRPS AND HERE SIMPLY
C ADDED TO THE REMAINING TERMS.
C
C IN THE RHF CASE EXPLICIT SPIN ADAPTED CODE IS USED
C
CEND
C
C CODED   JULY/90    JG
C
C EXTENDED FOR ROHF-MBPT(3), DECEMBER 92
C
      implicit none
C Common blocks
      integer pop(8,2),vrt(8,2),ntaa,ntbb,nf1(2),nf2(2)
      common/sym/pop,vrt,ntaa,ntbb,nf1,nf2
      integer nocco(2),nvrto(2)
      common/info/nocco,nvrto
      integer iintln,ifltln,iintfp,ialone,ibitwd
      common/machsp/iintln,ifltln,iintfp,ialone,ibitwd
      integer iflags(100)
      common/flags/iflags
      integer nstart,nirrep,irrepa(255),irrepb(255),dirprd(8,8)
      common/syminf/nstart,nirrep,irrepa,irrepb,dirprd
      integer irpdpd(8,22),isytyp(2,500),ntot(18)
      common/sympop/irpdpd,isytyp,ntot
C input variables
      integer maxcor,iuhf
      logical singles
C output variables
      double precision dvv(nf2(1)+iuhf*nf2(2))
C pre-allocated local variables
      integer icore(maxcor)
C Local variables
      integer mxcor,i0t1a,i0t1b,i001,i002,i003,i004,i005,dissyt,iflag,
     &        ioff,iofft,irrep,ioffd,irrepj,ispin,listt1,listt2,listt3,
     &        nocc,numsyt,nvrt2sq,nvrt
      double precision one,two,onem,half,fact
      DATA ONE,TWO,ONEM,HALF /1.0D0,2.0D0,-1.0D0,0.5D0/
C
      MXCOR=MAXCOR
C
C CONTRIBUTIONS DUE TO SINGLE EXCITATIONS
C
C
C ALLOCATE MEMORY FOR T1 AMPLITUDES
C
      if (singles) then
      I0T1A=MXCOR+1-NTAA*IINTFP
      MXCOR=MXCOR-NTAA*IINTFP
      CALL GETLST(ICORE(I0T1A),1,1,1,1,90)
C
      IF(IUHF.EQ.0) THEN
C
       I0T1B=I0T1A
C
      ELSE
C
       I0T1B=I0T1A-NTBB*IINTFP
       MXCOR=MXCOR-NTBB*IINTFP
       CALL GETLST(ICORE(I0T1B),1,1,1,2,90)
C
      ENDIF
C
C PERFORM MULTIPLICATION:    SUM M   T(M,A) L(M,B)
C
C   FACT IS HERE ALWAYS ONE
C
      FACT=ONE
C
      DO 300 ISPIN=1,IUHF+1
C
       IF(ISPIN.EQ.1) THEN
        IOFFT=I0T1A
        IOFFD=1
       ELSE
        IOFFT=I0T1B
        IOFFD=1+NF2(1)
       ENDIF
C
       DO 250 IRREP=1,NIRREP
C
        NOCC=POP(IRREP,ISPIN)
        NVRT=VRT(IRREP,ISPIN)
        IF(MIN(NVRT,NOCC).GT.0) THEN
         CALL XGEMM('N','T',NVRT,NVRT,NOCC,FACT,ICORE(IOFFT),NVRT,
     &              ICORE(IOFFT),NVRT,ONE,DVV(IOFFD),NVRT)
        ENDIF
        IOFFT=IOFFT+NOCC*NVRT*IINTFP
        IOFFD=IOFFD+NVRT*NVRT
250    CONTINUE
300   CONTINUE
      endif
C
C CONTRIBUTION DUE TO DOUBLES
C
C    RESET MXCOR
C
       MXCOR=MAXCOR
C
C LOOP OVER SPIN CASES (ISPIN=1: ALPHA DENSITY MATRIX, =2: BETA DENSITY MATRIX)
C
       DO 1000 ISPIN=1,IUHF+1
C
       IF(ISPIN.EQ.1) THEN
        IOFF=1
       ELSE
        IOFF=NF2(1)+1
       ENDIF
       IF(IUHF.EQ.1) THEN
C
C  CONTRIBUTION DUE TO AAAA OR BBBB AMPLITUDES (IS SKIPPED FOR RHF)
C
C  IN MBPT3 :
C            LISTT1 CONTAINS THE FIRST ORDER AMPLITUDES
C            LISTT2 CONTAINS THE SECOND ORDER AMPLITUDES
C  IN MBPT4 :
C            LISTT1 CONTAINS THE SECOND ORDER AMPLITUDES
C            LISTT2 CONTAINS THE SECOND ORDER AMPLITUDES
C            LISTT3 CONTAINS THE FIRST ORDER AMPLITUDES
C            LISTT4 CONTAINS THE THIRD ORDER DELTA AMPITUDES
C            LISTT5 CONTAINS THE X-CONTRIBUTION DUE TO QUADS
C  IN CC  :
C            LISTT1 CONTAINS THE CC AMPLITUDES
C            LISTT2 CONTAINS THE LAMBDA AMPLITUDES
C
C NOTE SOME FURTHER LOGICAL STUFF HERE
C
C MBPT2 LISTT1 IS EQUAL TO LISTT2 (DON'T READ IT TWICE)
C MBPT4 THERE IS AN ADDITIONAL STEP REQUIRED HERE (LISTT4+5 TIMES LISTT3)
C
C THE SYMMETRIZATION OF THE DENSITY MATRIX
C
C         D(A,B) = (D(A,B) + D(B,A))/2
C
C HAS TO BE CARRIED OUT FOR ALL METHODS EXCEPT MBPT2
C
         LISTT1=43+ISPIN
         LISTT2=43+ISPIN
         FACT=ONE
         IFLAG=1
C
C
C LOOP OVER IR REPS OF MN BLOCK (THE SAME IRREPS AS THE AF AND EF BLOCKS C HAVE
C
       DO 100 IRREP=1,NIRREP
C
C DETERMINE LENGTH OF EXPANDED VIRTUAL-VIRTUAL BLOCK
C
        NVRT2SQ=0
        DO 110 IRREPJ=1,NIRREP
         NVRT2SQ=NVRT2SQ+VRT(IRREPJ,ISPIN)*
     &                   VRT(DIRPRD(IRREPJ,IRREP),ISPIN)
110     CONTINUE
C
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT1))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT1))
        I001=1
        I002=I001+IINTFP*NUMSYT*NVRT2SQ
        I003=I002+IINTFP*NUMSYT*NVRT2SQ
        I004=I003+IINTFP*NUMSYT*DISSYT
        IF(MIN(NUMSYT,DISSYT).NE.0) THEN
         IF(I004.LT.MXCOR) THEN
C
C         IN CORE VERSION
C
          CALL DVVAA(ICORE(I001),ICORE(I002),ICORE(I003),DVV(IOFF),
     &               FACT,ISPIN,POP(1,ISPIN),VRT(1,ISPIN),
     &               DISSYT,NUMSYT,LISTT1,LISTT2,LISTT3,IRREP,IFLAG)
         ELSE
          CALL INSMEM('DVVAA',I004,MXCOR)
         ENDIF
        ELSE
        ENDIF
100    CONTINUE
       ENDIF
C
C       AB SPIN CASE
C
        LISTT1=46
        LISTT2=46
        FACT=ONE
        IFLAG=1


C      LOOP OVER IRREPS.
C
       DO 200 IRREP=1,NIRREP
C
        DISSYT=IRPDPD(IRREP,ISYTYP(1,46))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,46))
        I001=1
        I002=I001+IINTFP*NUMSYT*DISSYT
        I003=I002+IINTFP*NUMSYT*DISSYT
C
C WE NEED ONLY TWO TIMES THE LENGTH OF T2 FOR MBPT(2)
C
         I004=I003
        IF(MIN(NUMSYT,DISSYT).NE.0) THEN
         I005=I004+3*IINTFP*MAX(DISSYT,NUMSYT)
         IF(I005.LE.MXCOR) THEN
C
C         IN CORE VERSION
C
          CALL DVVAB(ICORE(I001),ICORE(I002),ICORE(I003),DVV(IOFF),
     &               FACT,ISPIN,POP(1,ISPIN),
     &               POP(1,3-ISPIN),VRT(1,ISPIN),VRT(1,3-ISPIN),
     &               DISSYT,NUMSYT,LISTT1,LISTT2,LISTT3,IRREP,
     &               ICORE(I004),IUHF,IFLAG)

         ELSE
          CALL INSMEM('DVVAB',I005,MXCOR)
         ENDIF
        ELSE
C
C
        ENDIF
200    CONTINUE
C
C  FOR PERTURBED CANONICAL ORBITALS, SET OFF-DIAGONAL
C  ELEMENTS TO ZERO
C
CSS        CALL DENCAN(DVV(IOFF),ICORE,NF2(ISPIN),NIRREP,VRT(1,ISPIN))

C
1000  CONTINUE
      RETURN
      END
