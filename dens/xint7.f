      SUBROUTINE XINT7(XIA,ICORE,MAXCOR,IUHF,bRedundant)
C
C THIS ROUTINE COMPUTES THE CONTRACTION
C
C    Z(a,i) =   SUM  <am||ne> Gamma(im,ne)
C               mne
C
C FOR THE FOLLOWING SPIN CASES:
C
C    Z(A,I) =   SUM  <AM||NE> Gamma(IM,NE) + SUM <Am|Ne> Gamma(Im,Ne)
C               MNE       (AAA)              mNe       (BAB)
C
C             + SUM  <Am|En>  Gamma(Im,En)   [RHF and UHF]
C               mEn   (AAB)
C    
C               
C    Z(a,i) =   SUM  <am||ne> Gamma(im,ne) + SUM <aM|nE> Gamma(iM,nE)
C               mne       (BBB)              MnE       (ABA)
C
C             + SUM  <aM|eN>  Gamma(iM,eN)   [UHF only]
C               MeN   (BBA)
C
C  SPIN ADAPTION HAS BEEN ADDED
C
C USING BOTH IN-CORE AND OLD-OF-CORE ALGORITHMS.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ALPHA,BETA,XIA,ONEM,TWO
      DIMENSION IOFFOV(8),IOFFSM(8,4),XIA(1)
      LOGICAL INCORE,RHF,bRedundant
      DIMENSION ICORE(MAXCOR),FULLMN(8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP2/ IRPDPD(8,22)
      COMMON /SYMPOP/ IRP_DM(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM2/ POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
      common /dropgeo/ ndrgeo
      COMMON /SHIFT/ ISHIFT 
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      DATA ONE   /1.0/
      DATA ONEM  /-1.0/
      DATA TWO/2.0D0/
C
C  FOR UHF RESORT SOME INTEGRALS
C
      IF(IUHF.EQ.1) CALL NWRNGAA(ICORE,MAXCOR,IUHF)
C
      RHF=.FALSE.
      IF(IUHF.EQ.0)RHF=.TRUE.
      XIAOFF=1
      DO 10 ISPIN=1,1+IUHF
       TARSIZ=NT(ISPIN)
       CALL IZERO(ICORE,TARSIZ*IINTFP)
       IOFFOV(1)=1
       DO 11 IRREP=1,NIRREP-1
        IOFFOV(IRREP+1)=IOFFOV(IRREP)+POP(IRREP,ISPIN)*
     &                    VRT(IRREP,ISPIN)*IINTFP
11     CONTINUE
C
C THIS PART CAN BE SKIPPED FOR RHF
C
       IF(.NOT.RHF) THEN
C
C SPIN CASES AAA AND BBB
C
C
C FILL OUT FULLMN VECTOR 
C
       CALL IZERO(FULLMN,NIRREP)
       DO 1000 IRREP=1,NIRREP
        DO 1001 IRREP1=1,NIRREP
         IRREP2=DIRPRD(IRREP1,IRREP)
         FULLMN(IRREP)=FULLMN(IRREP)+POP(IRREP1,ISPIN)*
     &                 POP(IRREP2,ISPIN) 
1001    CONTINUE
1000   CONTINUE
C
C SET THINGS UP
C
       LISTI=22+ISPIN + ISHIFT 
       LISTG=106+ISPIN
       DO 20 IRREPDO=1,NIRREP
        GAMDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTG))
        GAMDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTG)) 
        INTDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTI))
        INTDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTI))
        FULLSZ=FULLMN(IRREPDO)
C
C SEE IF THERE IS SUFFICIENT CORE TO DO THIS WITH AN INCORE ALGORITHM.
C  THE ALGORITHM USED HERE IS A GENERALIZATION OF A PUBLISHED METHOD.
C  [SEE J.F. STANTON, J. GAUSS, J.D. WATTS AND R.J. BARTLETT,
C                  JCP XX, XXXX (1991)].
C
        INEED=INTDIS*INTDSZ+GAMDIS*FULLSZ+TARSIZ
     &        +3*MAX(GAMDSZ,INTDSZ)
        IF(INEED.LE.MAXCOR/IINTFP)THEN
         INCORE=.TRUE.
        ELSE
         INCORE=.FALSE.
        ENDIF
C
C ALLOCATE CORE FOR TARGET AND INTEGRAL MATRICES.
C
        I000=1
        I010=I000+TARSIZ*IINTFP
        I020=I010+INTDSZ*INTDIS*IINTFP
C
C DO IN-CORE ALGORITHM IF POSSIBLE
C
        IF(INCORE)THEN 
         ALPHA=ONEM
         BETA=ONE
         I030=I020+FULLSZ*GAMDIS*IINTFP
         I040=I030+GAMDSZ
         if (ndrgeo.eq.0) then
           CALL GETTRN(ICORE(I020),ICORE(I030),GAMDSZ,GAMDIS,
     &                 2,IRREPDO,LISTG) 
         else  
           I040=I030+GAMDIS
           if (i040.le.MAXCOR) then
           CALL GETGX7TU(ICORE(I020),ICORE(I030),GAMDSZ,
     &                 GAMDIS,2,IRREPDO,LISTG,ispin) 
           else
            write(6,9013) i040,maxcor
 9013       format(4x,'---in xint7-1  i040=',i9,5x,'maxcor=',i9) 
            call errex 
           endif 
         endif 
C
C MATRIX RETURNED FROM THIS CALL IS Gamma(NE;I<M) [Gamma(ne;i<m)], 
C    WE NEED TO EXPAND THIS TO I,M [i,m].  
C     
         CALL SYMEXP(IRREPDO,POP(1,ISPIN),GAMDIS,ICORE(I020))
C
C NOW REDEFINE GAMDIS AND GAMDSZ TO REFER TO THE MATRIX WHICH IS HELD
C   (Gamma(NE;IM) [Gamma(ne;im)] 
C
         GAMDSZ=GAMDIS
         GAMDIS=FULLSZ
C
C NOW TRANSPOSE KET SIDE TO MAKE I THE SLOWEST INDEX
C
         I040=I030+GAMDSZ*IINTFP
         I050=I040+GAMDSZ*IINTFP
         I060=I050+GAMDSZ*IINTFP
         CALL SYMTR1(IRREPDO,POP(1,ISPIN),POP(1,ISPIN),GAMDSZ,
     &               ICORE(I020),ICORE(I030),ICORE(I040),
     &               ICORE(I050))
C
C NOW PLAY EXACTLY THE SAME GAME WITH THE INTEGRALS
C
         CALL GETLST(ICORE(I010),1,INTDIS,2,IRREPDO,LISTI)
         I040=I030+MAX(INTDIS,INTDSZ)*IINTFP
         I050=I040+MAX(INTDIS,INTDSZ)*IINTFP
         I060=I050+MAX(INTDIS,INTDSZ)*IINTFP
         CALL SYMTR1(IRREPDO,VRT(1,ISPIN),POP(1,ISPIN),INTDSZ,
     &               ICORE(I010),ICORE(I030),ICORE(I040),
     &               ICORE(I050))
         CALL SYMTR3(IRREPDO,VRT(1,ISPIN),POP(1,ISPIN),INTDSZ,
     &               INTDIS,ICORE(I010),ICORE(I030),
     &               ICORE(I040),ICORE(I050))
C
C NOW WE HAVE SOMETHING NICE.  THE INTEGRALS ARE ORDERED
C 
C     <NE||MA>   N,E ; M,A   [<en||ma> n,e ; m,a]
C
C  AND THE GAMMAS
C
C     G(NE,MI)   N,E ; M,I   [G(ne,mi) n,e ; m,i]
C
C AND WE CAN IMMEDIATELY GENERATE THE TARGET Z(A,I) [Z(a,i)] WITH NIRREP
C  MATRIX MULTIPLIES.  
C
        IOFFZ=I000
        IOFFI=I010
        IOFFG=I020
        DO 30 IRREPAI=1,NIRREP
C
C FIRST FIGURE OUT THE DIMENSIONS OF THE
C  MATRICES G(NEM;I) AND I(NEM;A) FOR EACH IRREP.
C
         NROWG=GAMDSZ*POP(DIRPRD(IRREPAI,IRREPDO),ISPIN)
         NROWI=NROWG
         NCOLG=POP(IRREPAI,ISPIN)
         NCOLI=VRT(IRREPAI,ISPIN)
C
C DO MATRIX MULTIPLY AND ACCUMULATE
C
         IF(MIN(NCOLI,NCOLG,NROWG).NE.0)THEN
          CALL XGEMM('T','N',NCOLI,NCOLG,NROWG,ALPHA,ICORE(IOFFI),
     &               NROWI,ICORE(IOFFG),NROWG,BETA,ICORE(IOFFZ),
     &               NCOLI)
         ENDIF
         IOFFI=IOFFI+NCOLI*NROWI*IINTFP
         IOFFG=IOFFG+NCOLG*NROWG*IINTFP
         IOFFZ=IOFFZ+NCOLG*NCOLI*IINTFP
30      CONTINUE
       ELSEIF(.NOT.INCORE)THEN
C
C OUT-OF-CORE ALGORITHM FOR THIS IRREP.  FIRST READ IN INTEGRALS
C  AND TRANSPOSE KET INDICES.
C
        CALL GETLST(ICORE(I010),1,INTDIS,2,IRREPDO,LISTI)
        I030=I020+INTDSZ*IINTFP
        I040=I030+INTDSZ*IINTFP
        I050=I040+INTDSZ*IINTFP
        CALL SYMTR1(IRREPDO,VRT(1,ISPIN),POP(1,ISPIN),INTDSZ,
     &              ICORE(I010),ICORE(I020),ICORE(I030),
     &              ICORE(I040))
C----------- for drop--mo ----------
C Make a space -- INDXG in GETINO4U -- for drop-mo business.
C Make another space -- NRDISG in GETINO4U -- for drop-mo business.
C
        IF (NDRGEO.NE.0) THEN 
          IINDXG = I020 
          ITMP = IINDXG + GAMDIS + 1 
          I020 = ITMP + IINTFP*FULLSZ + 1
        ENDIF 
C-----------------------------------
C
C DETERMINE HOW MANY LOGICAL RECORDS OF GAMMAS CAN BE HELD IN CORE
C  SIMULTANEOUSLY.
C
        ILEFT=MAXCOR-I020
        NUMIN=ILEFT/(FULLSZ*IINTFP)
C
C NOW DETERMINE HOW MANY PASSES MUST BE MADE TO PROCESS ALL INTEGRALS.
C
        NPASS=GAMDIS/NUMIN
        NLAST=GAMDIS
        IF(GAMDIS.NE.NPASS*NUMIN)THEN
         NLAST=GAMDIS-NPASS*NUMIN
         NPASS=NPASS+1
        ENDIF
c---------- for drop--mo -----------
c    pre-requisit for the out-of-core jobs in drop-mo case
c
      IF (NDRGEO.NE.0) THEN
        CALL GETINO6U(ICORE(IINDXG),GAMDIS,GAMDSZ,LISTG,IRREPDO,ISPIN)
        NIKI = 1
        IFIRDR = 1
      ENDIF 
c---------- for drop--mo -----------
C
C LOOP OVER PASSES
C
        DO 40 IPASS=1,NPASS
         IOFFZ=0  
         IOFFIL=0
         IOFFIR=0
         IF(IPASS.NE.NPASS)THEN
          NUMGET=NUMIN
         ELSE
          NUMGET=NLAST
         ENDIF
         IFIRST=1+(IPASS-1)*NUMIN
C
C PICK UP A LOAD OF GAMMAS AND EXPAND THEM TO
C
C         G(IM,NE) I,M ; N,E
C
         if (ndrgeo.eq.0) then 
           CALL GETLST(ICORE(I020),IFIRST,NUMGET,2,IRREPDO,LISTG)
         else  
           CALL GETGO6UO(ICORE(I020),ICORE(ITMP),IFIRDR,NUMGET,2,
     x                    IRREPDO,LISTG,
     x                    ispin,listg,gamdsz,ICORE(IINDXG),NIKI)
         endif 
         CALL SYMEXP2(IRREPDO,POP(1,ISPIN),FULLSZ,GAMDSZ,NUMGET,
     &                ICORE(I020),ICORE(I020))
C
C NOW PROCESS ALL N,E PAIRS WHICH ARE IN CORE
C
         ITHRU=0
         DO 50 INUMNE=IFIRST,IFIRST+NUMGET-1
          ITHRU=ITHRU+1
          IOFFIR=I010+IINTFP*INTDSZ*(INUMNE-1)
          IOFFGR=I020+IINTFP*FULLSZ*(ITHRU-1)
          IOFFIL=0
          IOFFGL=0
          DO 51 IRREPM=1,NIRREP
           IRREPI=DIRPRD(IRREPM,IRREPDO)
           IRREPA=IRREPI
           IOFFI =IOFFIR+IOFFIL
           IOFFG =IOFFGR+IOFFGL
           IOFFZ =IOFFOV(IRREPI)
           NROW  =VRT(IRREPA,ISPIN)
           NCOL  =POP(IRREPI,ISPIN)
           NSUM  =POP(IRREPM,ISPIN)
           ALPHA =ONE
           BETA  =ONE
           IF(MIN(NROW,NCOL,NSUM).GT.0)THEN
            CALL XGEMM('N','T',NROW,NCOL,NSUM,ALPHA,ICORE(IOFFI),NROW,
     &                 ICORE(IOFFG),NCOL,BETA,ICORE(IOFFZ),NROW)
           ENDIF
           IOFFIL=IOFFIL+NROW*NSUM*IINTFP
           IOFFGL=IOFFGL+NCOL*NSUM*IINTFP
51        CONTINUE
50       CONTINUE
40      CONTINUE
       ENDIF
20    CONTINUE
      ENDIF
C
C SPIN CASES BAB AND ABA
C
C    Z(A,I) =  Z(A,I) + SUM <Am|Ne> Gamma(Im,Ne)
C                       mNe       
C
      LISTI=23-ISPIN + ISHIFT 
      LISTI2=25 + ISHIFT 
      IF(IUHF.EQ.0)LISTI=21 + ISHIFT 
      LISTG=111-ISPIN

      DO 120 IRREPDO=1,NIRREP
       GAMDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTG))
       GAMDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTG))
       INTDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTI))
       INTDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTI))
C
C SEE IF THERE IS SUFFICIENT CORE FOR INCORE ALGORITHM.
C
       INEED=INTDIS*INTDSZ+GAMDIS*GAMDSZ+TARSIZ
       IF(INEED.LE.MAXCOR)THEN
        INCORE=.TRUE.
       ELSE
        INCORE=.FALSE.
       ENDIF
C
C ALLOCATE CORE FOR THE TARGET AND INTEGRAL MATRICES
C
       I000=1
       I010=I000+TARSIZ*IINTFP
       I020=I010+INTDSZ*INTDIS*IINTFP
C
C DO IN-CORE ALGORITHM
C
       IF(INCORE)THEN
C
C NOW READ THE INTEGRALS
C
       IF(bRedundant) THEN
          CALL GETLST(ICORE(I010),1,INTDIS,2,IRREPDO,LISTI)
       ELSE
          CALL GETLST_NR(ICORE(I010),ICORE(I020),MAXCOR-I020,
     &                  LISTI,IRREPDO)
       ENDIF

C
C SPIN ADAPT FOR RHF 
C
        IF(RHF) THEN 
         CALL SSCAL(INTDIS*INTDSZ,TWO,ICORE(I010),1)
         I030=I020+IINTFP*INTDIS*INTDSZ
         IF(I030.GE.MAXCOR) STOP 'XINT7'
         CALL GETLST(ICORE(I020),1,INTDIS,2,IRREPDO,LISTI2)
         CALL SAXPY(INTDIS*INTDSZ,ONEM,ICORE(I020),1,ICORE(I010),1)
        ENDIF
C
        I030=I020+GAMDSZ*GAMDIS*IINTFP
        I040=I030+GAMDSZ
         if (ndrgeo.eq.0) then
           CALL GETTRN (ICORE(I020),ICORE(I030),GAMDSZ,
     &                             GAMDIS,2,IRREPDO,LISTG)
         else
          I040=I030+GAMDIS
          if (i040.le.MAXCOR) then
          if (iuhf.eq.1 .and. ispin.eq.2) then 
           CALL GETGX7TR(ICORE(I020),ICORE(I030),GAMDSZ,
     &                             GAMDIS,2,IRREPDO,LISTG,ispin) 
          else
           CALL GETGX7T(ICORE(I020),ICORE(I030),GAMDSZ,
     &                             GAMDIS,2,IRREPDO,LISTG,ispin) 
          endif 
          else
            write(6,9014) i040,maxcor
 9014       format(4x,'---in xint7-2   i040=',i9,5x,'maxcor=',i9) 
            call errex 
          endif 
         endif 
C
C MATRIX RETURNED FROM THIS CALL IS :
C
C         Gamma(Ne,Im) N,e ; I,m (ISPIN=1)
C         Gamma(En,mI) E,n ; M,i (ISPIN=2)
C
C REDEFINE GAMDIS AND GAMDSZ TO REFER TO THIS STORAGE MODE.
C
        ITMP=GAMDIS
        GAMDIS=GAMDSZ
        GAMDSZ=ITMP
C
C TRANSPOSE KET SIDE TO MAKE I THE SLOWEST INDEX IF ISPIN = 1.
C  THIS GIVES
C
C         Gamma(Ne,mI) N,e ; m,I (ISPIN=1)
C 
C
        IF(ISPIN.EQ.1)THEN
         I040=I030+GAMDSZ*IINTFP
         I050=I040+GAMDSZ*IINTFP
         I060=I050+GAMDSZ*IINTFP
         CALL SYMTR1(IRREPDO,POP(1,1),POP(1,2),GAMDSZ,ICORE(I020),
     &               ICORE(I030),ICORE(I040),ICORE(I050))
        ENDIF
C
C INTEGRALS HAVE BEEN ALREADY READ IN.
C MATRIX RETURNED FROM THIS CALL IS :
C
C         I(eN;Am) (ISPIN=1)
C         I(En;aM) (ISPIN=2)
C
C TRANSPOSE KET SIDE TO MAKE A THE SLOWEST INDEX.
C
        I040=I030+INTDSZ*IINTFP
        I050=I040+INTDSZ*IINTFP
        I060=I050+INTDSZ*IINTFP
        CALL SYMTR1(IRREPDO,VRT(1,ISPIN),POP(1,3-ISPIN),INTDSZ,
     &              ICORE(I010),ICORE(I030),ICORE(I040),
     &              ICORE(I050))
C
C NOW WE HAVE:
C
C        I(eN;mA) AND Gamma(Ne;Im) (ISPIN=1)
C        I(En;Ma) AND Gamma(En;Mi) (ISPIN=2)
C
C FOR ISPIN=1 WE NEED TO TRANSPOSE THE BRA SIDE OF THE INTEGRALS.
C
        IF(ISPIN.EQ.1)THEN
         I040=I030+INTDIS*IINTFP
         I050=I040+INTDIS*IINTFP
         CALL SYMTR3(IRREPDO,VRT(1,2),POP(1,1),INTDSZ,INTDIS,
     &               ICORE(I010),ICORE(I030),ICORE(I040),
     &               ICORE(I050))
        ENDIF 
C
C NOW WE HAVE
C        
C        I(Ne;mA) AND Gamma(Ne;mI) (ISPIN=1)
C        I(En;Ma) AND Gamma(En;Mi) (ISPIN=2)
C       
C AND THE PRODUCT CAN BE FORMED.
C
        IOFFZ=I000
        IOFFI=I010
        IOFFG=I020
        ALPHA=ONE
        BETA =ONE
        DO 130 IRREPAI=1,NIRREP
C
C FIRST FIGURE OUT THE DIMENSIONS OF THE
C  MATRICES G(NEM;I) AND I(NEM;A) FOR EACH IRREP.
C
         NROWG=GAMDSZ*POP(DIRPRD(IRREPAI,IRREPDO),3-ISPIN)
         NROWI=NROWG
         NCOLG=POP(IRREPAI,ISPIN)
         NCOLI=VRT(IRREPAI,ISPIN)
C
C DO MATRIX MULTIPLY AND ACCUMULATE
C
         CALL XGEMM('T','N',NCOLI,NCOLG,NROWG,ALPHA,ICORE(IOFFI),
     &              NROWI,ICORE(IOFFG),NROWG,BETA,ICORE(IOFFZ),
     &              NCOLI)
         IOFFI=IOFFI+NCOLI*NROWI*IINTFP
         IOFFG=IOFFG+NCOLG*NROWG*IINTFP
         IOFFZ=IOFFZ+NCOLG*NCOLI*IINTFP
130     CONTINUE
       ELSEIF(.NOT.INCORE)THEN
C
C OUT-OF-CORE ALGORITHM FOR THIS IRREP.  
C
C
C GET OFFSET INFORMATION
C
        CALL SYMOFF(IOFFSM,IRREPDO)
C
C READ IN INTEGRALS.
C
        LISTI=20+ISPIN + ISHIFT 
        INTDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTI))
        INTDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTI))
C
C NOW READ THE INTEGRALS
C
        CALL GETLST(ICORE(I010),1,INTDIS,2,IRREPDO,LISTI)
C
C SPIN ADAPT FOR RHF 
C
        IF(RHF) THEN 
         CALL SSCAL(INTDIS*INTDSZ,TWO,ICORE(I010),1)
         I030=I020+IINTFP*INTDIS*INTDSZ
         IF(I030.GE.MAXCOR) STOP 'XINT7'
         CALL GETLST(ICORE(I020),1,INTDIS,2,IRREPDO,LISTI2)
         CALL SAXPY(INTDIS*INTDSZ,ONEM,ICORE(I020),1,ICORE(I010),1)
        ENDIF
C
C INTEGRALS ARE
C
C         I(Am;eN) (ISPIN=1)
C         I(aM;En) (ISPIN=2)
C
C AND KET PART MUST BE TRANSPOSED IF ISPIN=1.
C
        IF(ISPIN.EQ.1)THEN
         I030=I020+INTDSZ*IINTFP
         I040=I030+INTDSZ*IINTFP
         I050=I040+INTDSZ*IINTFP
         CALL SYMTR1(IRREPDO,VRT(1,2),POP(1,1),INTDSZ,
     &               ICORE(I010),ICORE(I030),ICORE(I040),
     &               ICORE(I050))
        ENDIF
C-------------- for drop--mo -------------
C Make a space -- INDXG in GETINO4(R) -- for drop-mo business.
C Make another space -- NRDISG in GETINO4(R) -- for drop-mo business.
C
        IF (NDRGEO.NE.0) THEN 
          IINDXG = I020 
          ITMP = IINDXG + GAMDIS + 1 
          I020 = ITMP + IINTFP*GAMDSZ + 1 
        ENDIF 
c----------------------------------------
C
C DETERMINE HOW MANY LOGICAL RECORDS OF GAMMAS CAN BE HELD IN CORE
C  SIMULTANEOUSLY.
C
        ILEFT=MAXCOR-I020
        NUMIN=ILEFT/(GAMDSZ*IINTFP)
C
C NOW DETERMINE HOW MANY PASSES MUST BE MADE TO PROCESS ALL INTEGRALS.
C
        NPASS=GAMDIS/NUMIN
        NLAST=GAMDIS
        IF(GAMDIS.NE.NPASS*NUMIN)THEN
         NLAST=GAMDIS-NPASS*NUMIN
         NPASS=NPASS+1
        ENDIF
C-------------- for drop--mo -------------
c    pre-requisit for the out-of-core jobs in drop-mo case
c
      IF (NDRGEO.NE.0) THEN
        if (iuhf.eq.1 .and. ispin.eq.2) then 
         CALL GETINO6R(ICORE(IINDXG),GAMDIS,GAMDSZ,LISTG,IRREPDO)
        else 
         CALL GETINO6 (ICORE(IINDXG),GAMDIS,GAMDSZ,LISTG,IRREPDO)
        endif 
        NIKI = 1
        IFIRDR = 1
      ENDIF 
c----------------------------------------
C
C LOOP OVER PASSES
C
        DO 140 IPASS=1,NPASS
         IOFFIL=0
         IOFFIR=0
         IF(IPASS.NE.NPASS)THEN
          NUMGET=NUMIN
         ELSE
          NUMGET=NLAST
         ENDIF
         IFIRST=1+(IPASS-1)*NUMIN
C
C PICK UP A LOAD OF GAMMAS.
C
C         Gamma(Im;Ne) (ISPIN=1)
C         Gamma(Mi;En) (ISPIN=2)
C
         if (ndrgeo.eq.0) then 
           CALL GETLST(ICORE(I020),IFIRST,NUMGET,2,IRREPDO,LISTG)
         else  
          CALL GETGO6O(ICORE(I020),ICORE(ITMP),IFIRDR,NUMGET,2,
     x                  IRREPDO,LISTG,
     x                  ispin,listg,gamdsz,ICORE(IINDXG),NIKI)
         endif 
C
C WE HAVE
C
C        I(Am;Ne) AND Gamma(Im;Ne) (ISPIN=1)
C        I(aM;En) AND Gamma(Mi,En) (ISPIN=2)
C
C NOW PROCESS ALL N,e [E,n] PAIRS WHICH ARE IN CORE
C
         ITHRU=0
         ALPHA=ONE
         BETA =ONE
         DO 150 INUMNE=IFIRST,IFIRST+NUMGET-1
          ITHRU=ITHRU+1
          IOFFIR=I010+IINTFP*INTDSZ*(INUMNE-1)
          IOFFGR=I020+IINTFP*GAMDSZ*(ITHRU-1)
          IOFFGL=0
          IF(ISPIN.EQ.1)THEN
           DO 151 IRREPM=1,NIRREP
            IRREPI=DIRPRD(IRREPM,IRREPDO)
            IRREPA=IRREPI
            IOFFI =IOFFIR+IOFFSM(IRREPM,4)*IINTFP
            IOFFG =IOFFGR+IOFFGL
            IOFFZ =IOFFOV(IRREPI)
            NROW  =VRT(IRREPA,ISPIN)
            NCOL  =POP(IRREPI,ISPIN)
            NSUM  =POP(IRREPM,3-ISPIN)
            IF(MIN(NROW,NCOL,NSUM).GT.0)THEN
             CALL XGEMM('N','T',NROW,NCOL,NSUM,ALPHA,ICORE(IOFFI),NROW,
     &                  ICORE(IOFFG),NCOL,BETA,ICORE(IOFFZ),NROW)
            ENDIF
            IOFFGL=IOFFGL+NCOL*NSUM*IINTFP
151        CONTINUE
          ELSE
           DO 152 IRREPI=1,NIRREP
            IRREPM=DIRPRD(IRREPI,IRREPDO)
            IRREPA=IRREPI
            IOFFI =IOFFIR+IOFFSM(IRREPM,3)*IINTFP 
            IOFFG =IOFFGR+IOFFGL
            IOFFZ =IOFFOV(IRREPI)
            NROW  =VRT(IRREPA,ISPIN)
            NCOL  =POP(IRREPI,ISPIN)
            NSUM  =POP(IRREPM,3-ISPIN)
            IF(MIN(NROW,NCOL,NSUM).GT.0)THEN
             CALL XGEMM('N','N',NROW,NCOL,NSUM,ALPHA,ICORE(IOFFI),NROW,
     &                  ICORE(IOFFG),NSUM,BETA,ICORE(IOFFZ),NROW)
            ENDIF
            IOFFGL=IOFFGL+NCOL*NSUM*IINTFP
152        CONTINUE
          ENDIF
150      CONTINUE
140     CONTINUE
       ENDIF
120   CONTINUE
C
C SPIN CASES AAB AND BBA
C
C             - SUM  <Am|En>  Gamma(Im,En)   (ISPIN=1)
C               mEn   
C
C             - SUM  <aM|eN>  Gamma(iM,eN)   (ISPIN=2 OR RHF)
C               MeN   
C
      LISTI=24+ISPIN + ISHIFT 
      LISTI2=21 + ISHIFT 
      LISTG=108+ISPIN
      IF(RHF)LISTG=110

      DO 220 IRREPDO=1,NIRREP
       GAMDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTG))
       GAMDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTG))
       INTDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTI))
       INTDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTI))
C
C SEE IF THERE IS SUFFICIENT CORE FOR INCORE ALGORITHM.
C
       INEED=INTDIS*INTDSZ+GAMDIS*GAMDSZ+TARSIZ
       IF(INEED.LE.MAXCOR)THEN
        INCORE=.TRUE.
       ELSE
        INCORE=.FALSE.
       ENDIF
C
C ALLOCATE CORE FOR THE TARGET AND INTEGRAL MATRICES
C
       I000=1
       I010=I000+TARSIZ*IINTFP
       I020=I010+INTDSZ*INTDIS*IINTFP
C
C NOW READ THE INTEGRALS
C
        CALL GETLST(ICORE(I010),1,INTDIS,2,IRREPDO,LISTI)

C
C SPIN ADAPT FOR RHF 
C
        IF(RHF) THEN 
         CALL SSCAL(INTDIS*INTDSZ,TWO,ICORE(I010),1)
         I030=I020+IINTFP*INTDIS*INTDSZ
         IF(I030.GE.MAXCOR) STOP 'XINT7'
         if(bRedundant) then
           CALL GETLST(ICORE(I020),1,INTDIS,2,IRREPDO,LISTI2)
         else
           CALL GETLST_NR(ICORE(I020),ICORE(I030),MAXCOR-I030,
     &                   LISTI2,IRREPDO)
         endif

         CALL SAXPY(INTDIS*INTDSZ,ONEM,ICORE(I020),1,ICORE(I010),1)
        ENDIF
C
C
C DO IN-CORE ALGORITHM
C
       IF(INCORE)THEN
        I030=I020+GAMDSZ*GAMDIS*IINTFP
        I040=I030+GAMDSZ
         if (ndrgeo.eq.0) then
           CALL GETTRN (ICORE(I020),ICORE(I030),GAMDSZ,
     &                             GAMDIS,2,IRREPDO,LISTG)
         else
          if (iuhf.eq.1 .and.ispin.eq.1) then 
           CALL GETGX7TR(ICORE(I020),ICORE(I030),GAMDSZ,
     &                             GAMDIS,2,IRREPDO,LISTG,ispin) 
          else
           CALL GETGX7T(ICORE(I020),ICORE(I030),GAMDSZ,
     &                             GAMDIS,2,IRREPDO,LISTG,ispin) 
          endif 
         endif 
C
C MATRIX RETURNED FROM THIS CALL IS :
C
C         Gamma(En;Im)  (ISPIN=1)
C         Gamma(Ne;Mi)  (ISPIN=2 OR RHF)
C
C REDEFINE GAMDIS AND GAMDSZ TO REFER TO THIS STORAGE MODE.
C
        ITMP=GAMDIS
        GAMDIS=GAMDSZ
        GAMDSZ=ITMP
C
C TRANSPOSE KET SIDE TO MAKE I THE SLOWEST INDEX IF ISPIN = 1.
C  THIS GIVES
C
C         Gamma(En;mI) 
C 
C
        IF(ISPIN.EQ.1.AND..NOT.RHF)THEN
         I040=I030+GAMDSZ*IINTFP
         I050=I040+GAMDSZ*IINTFP
         I060=I050+GAMDSZ*IINTFP
         CALL SYMTR1(IRREPDO,POP(1,1),POP(1,2),GAMDSZ,ICORE(I020),
     &               ICORE(I030),ICORE(I040),ICORE(I050))
        ENDIF
C
C TRANSPOSE BRA SIDE IF ISPIN=2 OR RHF.
C
        IF(ISPIN.EQ.2..OR.RHF)THEN
         I040=I030+GAMDIS*IINTFP
         I050=I040+GAMDIS*IINTFP
         CALL SYMTR3(IRREPDO,POP(1,1),VRT(1,2),GAMDSZ,GAMDIS,
     &               ICORE(I020),ICORE(I030),ICORE(I040),
     &               ICORE(I050))
        ENDIF
C
C INTEGRALS ARE ALREADY READ IN
C MATRIX RETURNED FROM THIS CALL IS :
C
C         I(En;Am) (ISPIN=1)
C         I(eN;aM) (ISPIN=2 OR RHF)
C
C TRANSPOSE KET SIDE TO MAKE A THE SLOWEST INDEX.
C
        I040=I030+INTDSZ*IINTFP
        I050=I040+INTDSZ*IINTFP
        I060=I050+INTDSZ*IINTFP
        CALL SYMTR1(IRREPDO,VRT(1,ISPIN),POP(1,3-ISPIN),INTDSZ,
     &              ICORE(I010),ICORE(I030),ICORE(I040),
     &              ICORE(I050))
C
C NOW WE HAVE:
C
C        I(En;mA) AND Gamma(En;mI) (ISPIN=1)
C        I(eN;Ma) AND Gamma(eN;Mi) (ISPIN=2 OR RHF)
C
C AND THE PRODUCT CAN BE FORMED.
C
        IOFFZ=I000
        IOFFI=I010
        IOFFG=I020
        ALPHA=ONE
        BETA =ONE
        DO 230 IRREPAI=1,NIRREP
C
C FIRST FIGURE OUT THE DIMENSIONS OF THE
C  MATRICES G(NEM;I) AND I(NEM;A) FOR EACH IRREP.
C
         NROWG=GAMDSZ*POP(DIRPRD(IRREPAI,IRREPDO),3-ISPIN)
         NROWI=NROWG
         NCOLG=POP(IRREPAI,ISPIN)
         NCOLI=VRT(IRREPAI,ISPIN)
C
C DO MATRIX MULTIPLY AND ACCUMULATE
C
         CALL XGEMM('T','N',NCOLI,NCOLG,NROWG,ALPHA,ICORE(IOFFI),
     &              NROWI,ICORE(IOFFG),NROWG,BETA,ICORE(IOFFZ),
     &              NCOLI)
         IOFFI=IOFFI+NCOLI*NROWI*IINTFP
         IOFFG=IOFFG+NCOLG*NROWG*IINTFP
         IOFFZ=IOFFZ+NCOLG*NCOLI*IINTFP
230     CONTINUE
       ELSEIF(.NOT.INCORE)THEN
C
C GET OFFSET INFORMATION
C
        CALL SYMOFF(IOFFSM,IRREPDO)
C
C OUT-OF-CORE ALGORITHM FOR THIS IRREP.  INTEGRALS ARE ALREADY READ IN.
C INTEGRALS ARE
C
C         I(Am;En) (ISPIN=1)
C         I(aM;eN) (ISPIN=2 OR RHF)
C
C AND KET PART MUST BE TRANSPOSED IF ISPIN=2.
C
        IF(ISPIN.EQ.2.OR.RHF)THEN
         I030=I020+INTDSZ*IINTFP
         I040=I030+INTDSZ*IINTFP
         I050=I040+INTDSZ*IINTFP
         CALL SYMTR1(IRREPDO,VRT(1,ISPIN),POP(1,3-ISPIN),INTDSZ,
     &               ICORE(I010),ICORE(I030),ICORE(I040),
     &               ICORE(I050))
        ENDIF
C----------  for drop--mo  -----------
C Make a space -- INDXG in GETINO4(R) -- for drop-mo business.
C Make another space -- NRDISG in GETINO4(R) -- for drop-mo business.
C
        IF (NDRGEO.NE.0) THEN 
          IINDXG = I020 
          ITMP = IINDXG + GAMDIS + 1 
          I020 = ITMP + IINTFP*GAMDSZ + 1 
        ENDIF 
c-------------------------------------
C
C DETERMINE HOW MANY LOGICAL RECORDS OF GAMMAS CAN BE HELD IN CORE
C  SIMULTANEOUSLY.
C
        ILEFT=MAXCOR-I020
        NUMIN=ILEFT/(GAMDSZ*IINTFP)
C
C NOW DETERMINE HOW MANY PASSES MUST BE MADE TO PROCESS ALL INTEGRALS.
C
        NPASS=GAMDIS/NUMIN
        NLAST=GAMDIS
        IF(GAMDIS.NE.NPASS*NUMIN)THEN
         NLAST=GAMDIS-NPASS*NUMIN
         NPASS=NPASS+1
        ENDIF
C----------  for drop--mo  -----------
c
c    pre-requisit for the out-of-core jobs in drop-mo case
c
      IF (NDRGEO.NE.0) THEN
        if (iuhf.eq.1 .and. ispin.eq.1) then 
         CALL GETINO6R(ICORE(IINDXG),GAMDIS,GAMDSZ,LISTG,IRREPDO)
        else 
         CALL GETINO6 (ICORE(IINDXG),GAMDIS,GAMDSZ,LISTG,IRREPDO)
        endif 
        NIKI = 1
        IFIRDR = 1
      ENDIF 
c-------------------------------------
C
C LOOP OVER PASSES
C
        DO 240 IPASS=1,NPASS
         IOFFZ=0  
         IOFFIL=0
         IOFFIR=0
         IF(IPASS.NE.NPASS)THEN
          NUMGET=NUMIN
         ELSE
          NUMGET=NLAST
         ENDIF
         IFIRST=1+(IPASS-1)*NUMIN
C
C PICK UP A LOAD OF GAMMAS.
C
C         Gamma(Im;En) (ISPIN=1)
C         Gamma(Mi;Ne) (ISPIN=2 OR RHF)
C
         if (ndrgeo.eq.0) then 
           CALL GETLST(ICORE(I020),IFIRST,NUMGET,2,IRREPDO,LISTG)
         else 
           CALL GETGO6O(ICORE(I020),ICORE(ITMP),IFIRDR,NUMGET,2,
     x                  IRREPDO,LISTG,
     x                  ispin,listg,gamdsz,ICORE(IINDXG),NIKI)
         endif 
C
C WE HAVE
C
C        I(Am;En) AND Gamma(Im;En) (ISPIN=1)
C        I(aM;Ne) AND Gamma(Mi,Ne) (ISPIN=2 OR RHF)
C
C NOW PROCESS ALL N,e [E,n] PAIRS WHICH ARE IN CORE
C
         ITHRU=0
         ALPHA =ONE
         BETA  =ONE
         DO 250 INUMNE=IFIRST,IFIRST+NUMGET-1
          ITHRU=ITHRU+1
          IOFFIR=I010+IINTFP*INTDSZ*(INUMNE-1)
          IOFFGR=I020+IINTFP*GAMDSZ*(ITHRU-1)
          IOFFGL=0
          IF(ISPIN.EQ.1.AND..NOT.RHF)THEN
           DO 251 IRREPM=1,NIRREP
            IRREPI=DIRPRD(IRREPM,IRREPDO)
            IRREPA=IRREPI
            IOFFI =IOFFIR+IOFFSM(IRREPM,4)*IINTFP
            IOFFG =IOFFGR+IOFFGL
            IOFFZ =IOFFOV(IRREPI)
            NROW  =VRT(IRREPA,ISPIN)
            NCOL  =POP(IRREPI,ISPIN)
            NSUM  =POP(IRREPM,3-ISPIN)
            IF(MIN(NROW,NCOL,NSUM).GT.0)THEN
             CALL XGEMM('N','T',NROW,NCOL,NSUM,ALPHA,ICORE(IOFFI),NROW,
     &                  ICORE(IOFFG),NCOL,BETA,ICORE(IOFFZ),NROW)
            ENDIF
            IOFFGL=IOFFGL+NCOL*NSUM*IINTFP
251        CONTINUE
          ELSE
           DO 252 IRREPI=1,NIRREP
            IRREPM=DIRPRD(IRREPI,IRREPDO)
            IRREPA=IRREPI
            IOFFI =IOFFIR+IOFFSM(IRREPM,3)*IINTFP
            IOFFG =IOFFGR+IOFFGL
            IOFFZ =IOFFOV(IRREPI)
            NROW  =VRT(IRREPA,ISPIN)
            NCOL  =POP(IRREPI,ISPIN)
            NSUM  =POP(IRREPM,3-ISPIN)
            IF(MIN(NROW,NCOL,NSUM).GT.0)THEN
             CALL XGEMM('N','N',NROW,NCOL,NSUM,ALPHA,ICORE(IOFFI),NROW,
     &                  ICORE(IOFFG),NSUM,BETA,ICORE(IOFFZ),NROW)
            ENDIF
            IOFFGL=IOFFGL+NCOL*NSUM*IINTFP
252        CONTINUE
          ENDIF
250      CONTINUE
240     CONTINUE
       ENDIF
220   CONTINUE
C
C NOW INCREMENT XIA.
C
      CALL SAXPY(TARSIZ,ONE,ICORE(I000),1,XIA(XIAOFF),1)
      XIAOFF=XIAOFF+TARSIZ
10    CONTINUE
C
C    FOR UHF RESORT SOME INTEGRALS 
C
      IF(IUHF.EQ.1) CALL NWRNGAA(ICORE,MAXCOR,IUHF)
C 
      RETURN
      END
