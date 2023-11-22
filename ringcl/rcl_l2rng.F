C
C DRIVER FOR RING CONTRIBUTIONS TO L2 AMPLITUDES.
C
C THIS CONTRIBUTION IS:
C
C        ab              ab
C       Z   =  P(ij|ab) Q
C        ij              ij
C
C WHERE P(ij|ab) IS THE USUAL PERMUTATION OPERATOR, AND THE QUANTITY
C  Q IS DEFINED BY 
C
C        ab           ae
C       Q    =  SUM  L   W(ejmb)
C        ij     m,e   im                                         
c
C 
C  WHERE W IS A TWO-PARTICLE INTERMEDIATE WHICH IS DIFFERENT FOR
C  DIFFERENT METHODS.  
C
C
C   METHOD                          DEFINITION OF W(mbej)
C   ------                          ---------------------
C                                              ef
C    RCCD                     <ej||mb> +  SUM  t  <jn||bf>
C                                               mn        
C
C                                        =  2 W(ejmb) - <ej||mb>
C
C  SINCE THE INTEGRALS REPRESENTING THE LEADING TERM OF THE EXPANSION ARE
C  USUALLY STORED <bm||ej>, THIS ROUTINE ACTUALLY COMPUTES THE *NEGATIVE*
C  OF Z.  THE INCREMENTS COMPUTED IN THE VARIOUS CALLS TO LRNGDRV ARE ADDED
C  TOGETHER AND THEN NEGATED IN ROUTINE SUMRNG.
C  IN THE CODES TO COMPUTE THE VARIOUS PIECES OF W(mbej) THE SIGNS ARE 
C  CONSISTENT WITH THE EQUATION ABOVE.
C
C FOR THE VARIOUS SPIN CASES, Z IS DEFINED BY:
C
C       AB    AB  BA  AB  BA
C      Z  =  Q  +Q  -Q  -Q     (AAAA)
C       IJ    IJ  JI  IJ  JI
C
C       Ab    Ab  bA  Ab  bA
C      Z  =  Q  +Q  +Q  +Q     (ABAB)
C       Ij    Ij  jI  Ij  jI
C
C       aB    aB  Ba  aB  Ba
C      Z  =  Q  +Q  +Q  +Q     (ABBA)
C       Ij    Ij  jI  Ij  jI
C    
C
C  THE BBBB, BABA AND BAAB SPIN CASES CAN BE OBTAINED FROM THESE EQUATIONS
C   BY CHANGING ALL ALPHA LABELS TO BETA AND ALL BETA LABELS TO ALPHA.
C
C
C  THE SPIN CASES FOR Q ARE GIVEN BY:
C
C  AB             AE           Ae                
C Q  = Q(AAAA) = t  W(EJMB) + t  W(eJmB)
C  IJ             IM           Im                   
C
C  ab             ae           aE                
C Q  = Q(BBBB) = t  W(ejmb) + t  W(EjMb)
C  ij             im           iM                   
C
C  Ab             AE           Ae                
C Q  = Q(ABAB) = t  W(EjMb) + t  W(ejmb)
C  Ij             IM           Im                   
C
C  aB             ae           aE
C Q  = Q(BABA) = t  W(eJmB) + t  W(EJMB)
C  iJ             im           iM
C
C  Ab               eA
C Q  = Q(BAAB) = - t  W(eJMb)
C  iJ               iM      
C
C  aB               Ae
C Q  = Q(ABBA) = - t  W(EjmB)
C  Ij               Im
C
C (ALL QUANTITIES ARE SUMMED OVER m AND e).
C
C
C HAVING THE Q QUANTITIES, THE RING CONTRIBUTION TO THE L2 AMPLITUDES
C  (THE Z(ij,ab)) CAN BE COMPUTED FROM:
C
C                 AB             AB
C SPIN CASE AA : Z   = P(IJ|AB) Q
C                 IJ             IJ
C
C                 ab             ab
C SPIN CASE BB : Z   = P(ij|ab) Q
C                 ij             ij
C
C                 Ab    Ab    bA    bA    Ab
C SPIN CASE AB : Z   = Q   + Q   - Q   - Q  
C                 Ij    Ij    jI    Ij    jI
C
C
C    OR Z(AB)=Q(ABAB)+Q(BABA)-Q(ABBA)-Q(BAAB)
C
C
C THE COMMENTS PRECEDING EACH CALL TO LRNGDRV INCLUDE ONLY THE
C  LEADING TERM OF THE W(ejmb) INTERMEDIATE WHICH IS USED.  IT SHOULD
C  BE UNDERSTOOD THAT THE ACTUAL W IS BEING USED.
C
      SUBROUTINE RCL_L2RNG(ICORE,MAXCOR,ISPIN,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD
      CHARACTER*4 SPCASE
      LOGICAL MBPT4
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      DIMENSION ICORE(MAXCOR)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /METH/MBPT2,MBPt3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /FLAGS/ IFLAGS(100)
      DATA ONE /1.0/
      NNP1O2(I)=I*(I+1)/2 
C
      INCREM=0
      NOCCA=NOCCO(1)
      NOCCB=NOCCO(2)
      NVRTA=NVRTO(1)
      NVRTB=NVRTO(2)
      IF(ISPIN.LT.3)THEN
       IF(ISPIN.EQ.1)SPCASE='AAAA'
       IF(ISPIN.EQ.2)SPCASE='BBBB'
C***********************************************************************
C
C   Q(AAAA) AND Q(BBBB).  COMMENTS BELOW REFER TO THE AAAA CASE ONLY.
C                            FLIPPING THE SENSE OF A AND B WILL GIVE
C                            THE EQUATIONS FOR THE BBBB CASE.
C                           
C
C***********************************************************************
C               AB         AE
C   SOLVE FOR  Q   = SUM  T    <JE||MB>  [FIRST PART OF Q(AAAA)]
C               IJ   M,E   IM
C***********************************************************************
       LISTT=133+ISPIN
       LISTW=53+ISPIN
       CALL RCL_LRNGDRV_U(ICORE,MAXCOR,LISTW,LISTT,'TxW',MAXSIZ,
     &                     INCREM,39+ISPIN)
C***********************************************************************
C               AB         Ae
C   SOLVE FOR  Q   = SUM  T    <Jm|Be>  [SECOND PART OF Q(AAAA)].
C               IJ   m,e   Im
C***********************************************************************
       LISTW=55+ISPIN
       LISTT=138-ISPIN
       CALL RCL_LRNGDRV_U(ICORE,MAXCOR,LISTW,LISTT,'TxW',MAXSIZ,
     &                    1,39+ISPIN)

       NOCC=NOCCO(ISPIN)
       NVRT=NVRTO(ISPIN)
       NVSIZ=(NVRT*(NVRT-1))/2
       NOSIZ=(NOCC*(NOCC-1))/2
       NTOTSZ=ISYMSZ(ISYTYP(1,39+ISPIN),ISYTYP(2,39+ISPIN))
       NTOTSZ2=ISYMSZ(ISYTYP(1,60+ISPIN),ISYTYP(2,60+ISPIN))
       ISCRSZ=NVSIZ+NOSIZ+NOCC*NVRT
       I000=1
       I010=I000+MAX(NTOTSZ,NTOTSZ2)*IINTFP
       I020=I010+MAX(NTOTSZ,NTOTSZ2)*IINTFP
       I030=I020+ISCRSZ
       IEND=I030
       IF (IEND .GE. MAXCOR) CALL INSMEM("RCC_RNG",IEND,MAXCOR)
       MEMLEFT = MAXCOR - IEND
       CALL GETALL(ICORE(I000),NTOTSZ,1,39+ISPIN)
       CALL RCC_SST03I(ICORE(I000),ICORE(I010),ICORE(IEND),MEMLEFT,
     &                 NTOTSZ,NTOTSZ2,ICORE(I020),SPCASE,1)
C
       CALL GETALL(ICORE(I000),NTOTSZ2,1,60+ISPIN)
       CALL SAXPY(NTOTSZ2,ONE,ICORE(I000),1,ICORE(I010),1)
#ifdef _DEBUG_LVL0
       If (ispin .eq.1) call checksum("rcc_drng,61 :",ICORE(I010),
     &                                 NTOTSZ2)
       If (ispin .eq.2) call checksum("rcc_drng,62 :",ICORE(I010),
     &                                 NTOTSZ2)
#endif
       CALL PUTALL(ICORE(I010),NTOTSZ2,1,60+ISPIN)
C 
      ELSEIF(ISPIN.EQ.3)THEN
C 
C UHF AND ROHF CODE
C
       IF(IUHF.NE.0) THEN
C
C***********************************************************************
C
C  SPIN CASES  ABAB, BABA
C
C***********************************************************************
       MAXSIZ=0
C***********************************************************************
C               Ab           AE
C   SOLVE FOR  Q   =   SUM  T    <Mj|Eb>  [PART 1 OF Q(ABAB)]
C               Ij     M,E   IM
C***********************************************************************
        LISTT=134
        LISTW=57
        CALL RCL_LRNGDRV_U(ICORE,MAXCOR,LISTW,LISTT,'TxW',MAXSIZ,
     &                     INCREM,42)
C***********************************************************************
C               Ab         Ae
C   SOLVE FOR  Q   = SUM  T    <mb||je> [PART 2 OF Q(ABAB)]
C               Ij   m,e   Im
C***********************************************************************
C
C USE AA LIST FOR RHF REFERENCE SINCE IT IS IDENTICAL WITH THE (NOT STOR
C  BB LIST.
C
       LISTT=137
       LISTW=55
       CALL RCL_LRNGDRV_U(ICORE,MAXCOR,LISTW,LISTT,'TxW',MAXSIZ,
     &                    1,42)
C***********************************************************************
C               aB         aE
C   SOLVE FOR  Q   = SUM  T    <MB||JE> [PART 1 OF Q(BABA)].
C               iJ   M,E   iM
C***********************************************************************
C
C NEED TO CHANGE LISTW FOR W(mbej) INTERMEDIATE PICKUP!
C
       LISTT=136
       LISTW=54
       CALL RCL_LRNGDRV_U(ICORE,MAXCOR,LISTW,LISTT,'WxT',MAXSIZ,
     &                    1,42)
C***********************************************************************
C               aB           ae
C   SOLVE FOR  Q   =   SUM  T    <Bm|Je> [PART 2 OF Q(BABA)].
C               iJ     m,e   im
C***********************************************************************
C
C NEED TO CHANGE LISTW FOR W(mbej) INTERMEDIATE PICKUP!
C
       LISTT=135
       LISTW=56
       CALL RCL_LRNGDRV_U(ICORE,MAXCOR,LISTW,LISTT,'WxT',MAXSIZ,
     &                    1,42)
C***********************************************************************
C SPIN-ADAPTED RHF CODE
C
       ELSE IF(IUHF.EQ.0) THEN
C
C FIRST TERM  0.5 (2 T(Im,Ae) - T(Im,Ea))* (2 W(MbEj)+ W(MbeJ))
C
        LISTT=137
        LISTW=56

        CALL RCL_LRNGDRV_R(ICORE,MAXCOR,LISTW,LISTT,'WxT',MAXSIZ,0,
     &                     42)
       ENDIF
C
C SUM TOGETHER THE (AI,bj) AND (Aj,bI) RING INCREMENTS AND FORM THE
C  FULL Z(AB) PIECE.  
C
C  ***THE INCREMENTS COMPUTED THUS FAR ARE NEGATED IN THE ROUTINE SUMRNG***
C
       NSZAB=ISYMSZ(ISYTYP(1,42),ISYTYP(2,42))
       ISCRSZ=(NOCCA+NOCCB)*(NVRTA+NVRTB)
       I000=1
       I010=I000+NSZAB*IINTFP
       I020=I010+NSZAB*IINTFP
       I030=I020+ISCRSZ
       IF(I030.GT.MAXCOR)CALL INSMEM('RCL_RNG',I030,MAXCOR)
       CALL RCC_SUMRNG(ICORE(I000),ICORE(I010),ICORE(I020),NSZAB,NSZAB,
     &                 ISCRSZ)
C
C NOW CONVERT THESE TO Ab-Ij INCREMENTS AND AUGMENT THE T2 INCREMENT LIST.
C
       I030=I020+ISCRSZ
       IF(I030.GT.MAXCOR)CALL INSMEM('L2RNG',I030,MAXCOR)
       CALL SST02I(ICORE(I000),ICORE(I010),NSZAB,NSZAB,ICORE(I020),
     &             'AABB')
C
C NOW SUM THIS WITH THE EXISTING AB INCREMENT AND OVERWRITE IT.
C
#ifdef _DEBUG_LVL0
       call checksum("rcc_drng,63 :",ICORE(I010),nszab)
#endif
       CALL SUMSYM(ICORE(I010),ICORE(I000),NSZAB,63)
      ENDIF
      RETURN
      END
