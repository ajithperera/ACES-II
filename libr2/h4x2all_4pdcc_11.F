      SUBROUTINE H4X2ALL_4PDCC_11(ICORE,MAXCOR,ISPIN,IUHF,TAU,
     &                           T3L)
C
C DRIVER FOR RING CONTRIBUTIONS TO X(IJ,AB) AMPLITUDES.
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
C       Q    =  SUM  t   H4(mejb)
C        ij     m,e   im                                         
c
C 
C  WHERE H4 IS A TWO-PARTICLE INTERMEDIATE WHICH IS DIFFERENT FOR
C   DIFFERENT METHODS.  
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
C Q  = Q(AAAA) = t  H4(EMJB) + t  H4(meJB)
C  IJ             IM           Im                   
C
C  ab             ae           aE                
C Q  = Q(BBBB) = t  H4(mejb) + t  H4(MEjb)
C  ij             im           iM                   
C
C  Ab             AE           Ae                
C Q  = Q(ABAB) = t  H4(MEjb) + t  H4(mejb)
C  Ij             IM           Im                   
C
C  aB             ae           aE
C Q  = Q(BABA) = t  H4(meJB) + t  H4(MEJB)
C  iJ             im           iM
C
C  Ab               eA
C Q  = Q(BAAB) = - t  H4(MeJb)
C  iJ               iM      
C
C  aB               Ae
C Q  = Q(ABBA) = - t  H4(mEjB)
C  Ij               Im
C
C (ALL QUANTITIES ARE SUMMED OVER m AND e).
C
C
C HAVING THE Q QUANTITIES, THE RING CONTRIBUTION TO THE X(IJ,AB) AMPLITUDES
C  (THE Z(ij,ab)) CAN BE COMPUTED FROM:
C
C                 AB             AB
C SPIN CASE AA : X   = P(IJ|AB) Q
C                 IJ             IJ
C
C                 ab             ab
C SPIN CASE BB : X   = P(ij|ab) Q
C                 ij             ij
C
C                 Ab    Ab    bA    bA    Ab
C SPIN CASE AB : X   = Q   + Q   - Q   - Q  
C                 Ij    Ij    jI    Ij    jI
C
C
C    OR X(AB)=Q(ABAB)+Q(BABA)-Q(ABBA)-Q(BAAB)
C
C NOTE THAT WE CALCULATE FOR THE SPIN CASES AAAA AND BBBB THE POSITIVE
C Q(IA,JB) WHILE FOR ALL OTHER SPIN CASES THE NEGATIVE Q(IA,JB) IS     
C EVALUATED.
C
C FOR CCSD IN ADDITION A L1*T1**3 TERM IS ADDED WHICH IS CALCULATED IN
C ROUTINE L1T13. ALSO, WE NEED THIS FUNNY TAU : T(IJ,AB)-2*T(I,B)*T(J,A)
C FOR THE FISRT TWO SPIN CASES. THE FACTOR IS TWO AND THIS CONTRIBUTION
C IS FORMED IN F2TAU.
C 
C THE FACTORS ARE CHOSEN SO THAT HERE THE OVERALL FACTOR IS ONE OR ONEM,
C THE FACTOR OF HALF IS INCLUDED ALREADY IN H4(MEJB). 
C
CEND 
C
C  CODED AUGUST/90  JG
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TAU,T3L
      INTEGER DIRPRD,POP,VRT,DISSYQ
      CHARACTER*4 SPCASE
      DIMENSION ICORE(MAXCOR),IOFFT1(2),IOFFL1(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON/INFO/NOCCO(2),NVRTO(2) 
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)

#include "parcc.com"
C
      DATA ONE /1.0D0/,TWO/2.0D0/,THREE/3.D0/,THREEM /-3.D0/
C
      NOCCA=NOCCO(1)
      NOCCB=NOCCO(2)
      NVRTA=NVRTO(1)
      NVRTB=NVRTO(2)
C 
C   FOR CCSD (TAU = TRUE) GET T1 AND L1 AMPLITUDES
C
      MXCOR=MAXCOR 
      IF(TAU) THEN
       IOFFT1(1)=MXCOR+1-IINTFP*NT(1)
       IOFFL1(1)=IOFFT1(1)-IINTFP*NT(1)
       MXCOR=MXCOR-2*IINTFP*NT(1)
       CALL GETLST(ICORE(IOFFT1(1)),1,1,1,1,90)
       CALL GETLST(ICORE(IOFFL1(1)),1,1,2,1,190)
       IF(IUHF.EQ.0) THEN
        IOFFT1(2)=IOFFT1(1)
        IOFFL1(2)=IOFFL1(1)
       ELSE
        IOFFT1(2)=IOFFL1(1)-IINTFP*NT(2)
        IOFFL1(2)=IOFFT1(2)-IINTFP*NT(2)
        MXCOR=MXCOR-2*IINTFP*NT(2)
        CALL GETLST(ICORE(IOFFT1(2)),1,1,1,2,90)
        CALL GETLST(ICORE(IOFFL1(2)),1,1,2,2,190)
       ENDIF
      ENDIF
      
      INCREM=0
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
C               AB         EA
C   SOLVE FOR  Q   = SUM  T   *[- H4(MEJB)]  [FIRST PART OF Q(AAAA)]
C               IJ   M,E   IM
C***********************************************************************
       LISTT=33+ISPIN
       LISTW=253+ISPIN
       CALL HRNGDRV_1(ICORE,MXCOR,LISTW,LISTT,'TxW',MAXSIZ,
     &             INCREM,39+ISPIN,TAU,ICORE(IOFFT1(ISPIN)),
     &             ICORE(IOFFT1(ISPIN)),POP(1,ISPIN),POP(1,ISPIN),
     &             VRT(1,ISPIN),VRT(1,ISPIN),TWO,SPCASE,IUHF)
C***********************************************************************
C               AB         Ae
C   SOLVE FOR  Q   = SUM  T    H(meJB)  [SECOND PART OF Q(AAAA)].
C               IJ   m,e   Im
C**********************************************************************
       LISTW=258-ISPIN
       LISTT=38-ISPIN
       CALL HRNGDRV_1(ICORE,MXCOR,LISTW,LISTT,'TxW',MAXSIZ,
     &             1,39+ISPIN,.FALSE.,ICORE(IOFFT1(ISPIN)),
     &             ICORE(IOFFT1(ISPIN)),POP(1,ISPIN),POP(1,ISPIN),
     &             VRT(1,ISPIN),VRT(1,ISPIN),TWO,SPCASE,IUHF) 
C
C  FOR CCSD, THE L1*T1*T1*T1 CONTRIBUTION MUST BE ADDED HERE
C
       IF(T3L) THEN
        NUMSYQ=IRPDPD(1,ISYTYP(2,39+ISPIN))
        DISSYQ=IRPDPD(1,ISYTYP(1,39+ISPIN))
        I000=1
        I001=I000+IINTFP*NUMSYQ*DISSYQ
        CALL L1T13(ICORE(I000),DISSYQ,NUMSYQ,ICORE(I001),
     &             ICORE(IOFFT1(ISPIN)),ICORE(IOFFT1(ISPIN)),
     &             ICORE(IOFFL1(ISPIN)),THREE,ISPIN,ISPIN,39+ISPIN,'N')
       ENDIF
C
C RESORT AND ANTISYMMETRIZE X2 INCREMENT, AND THEN AUGMENT THE 
C  T2 INCREMENT LIST. 
C
       NOCC=NOCCO(ISPIN)
       NVRT=NVRTO(ISPIN)
       NVSIZ=(NVRT*(NVRT-1))/2
       NOSIZ=(NOCC*(NOCC-1))/2
       NTOTSZ=ISYMSZ(ISYTYP(1,39+ISPIN),ISYTYP(2,39+ISPIN))
       NTOTSZ2=ISYMSZ(ISYTYP(1,13+ISPIN),ISYTYP(2,13+ISPIN))
       ISCRSZ=NVSIZ+NOSIZ+NOCC*NVRT
       I000=1
       I010=I000+MAX(NTOTSZ,NTOTSZ2)*IINTFP
       I020=I010+MAX(NTOTSZ,NTOTSZ2)*IINTFP
       I030=I020+ISCRSZ
       CALL GETALL(ICORE(I000),NTOTSZ,1,39+ISPIN)
       CALL SST03I(ICORE(I000),ICORE(I010),NTOTSZ,NTOTSZ2,ICORE(I020),
     &             SPCASE,1)
       CALL GETALL(ICORE(I000),NTOTSZ2,1,113+ISPIN)
       CALL SAXPY(NTOTSZ2,ONE,ICORE(I000),1,ICORE(I010),1)

#ifdef _DEBUG_LVL0
       Call checksum("@-H4X2ALL_4PDCC_1_AA:",ICORE(I010),NTOTSZ2)
#endif 
CSSS       IF (ISPAR .AND. COULOMB) THEN
           CALL PUTALL(ICORE(I010),NTOTSZ2,1,213+ISPIN)
           CALL PUTALL(ICORE(I010),NTOTSZ2,1,113+ISPIN)
CSSS       ELSE 
CSSS           CALL PUTALL(ICORE(I010),NTOTSZ2,1,113+ISPIN)
CSSS       ENDIF 
      ELSEIF(ISPIN.EQ.3)THEN
C
C***********************************************************************
C
C  SPIN CASES  ABAB, BABA, ABBA AND BAAB.
C
C***********************************************************************
       MAXSIZ=0
C***********************************************************************
C               Ab           EA
C   SOLVE FOR  Q   =  SUM  T    H(MEjb)  [PART 1 OF Q(ABAB)]
C               Ij     M,E   IM
C***********************************************************************
        LISTT=34
        LISTW=256
        CALL HRNGDRV_1(ICORE,MXCOR,LISTW,LISTT,'TxW',MAXSIZ,INCREM,42,
     &               TAU,ICORE(IOFFT1(1)),ICORE(IOFFT1(1)),      
     &               POP(1,1),POP(1,1),VRT(1,1),VRT(1,1),TWO,
     &               'AAAA',IUHF)
C***********************************************************************
C               Ab         Ae
C   SOLVE FOR  Q   = SUM  T   *[-H(mejb)] [PART 2 OF Q(ABAB)]
C               Ij   m,e   Im
C***********************************************************************
C
C USE AA LIST FOR RHF REFERENCE SINCE IT IS IDENTICAL WITH THE (NOT STOR
C  BB LIST.
C
       LISTW=254+IUHF
       LISTT=37
       CALL HRNGDRV_1(ICORE,MXCOR,LISTW,LISTT,'TxW',MAXSIZ,1,42,
     &              .FALSE.,ICORE(IOFFT1(1)),ICORE(IOFFT1(1)),
     &               POP(1,1),POP(1,1),VRT(1,1),VRT(1,1),TWO,
     &              'AAAA',IUHF)
C
C  FOR CCSD, THE L1*T1*T1*T1 CONTRIBUTION MUST BE ADDED HERE
C
       IF(T3L) THEN
        NUMSYQ=IRPDPD(1,ISYTYP(2,42))
        DISSYQ=IRPDPD(1,ISYTYP(1,42))
        I000=1
        I001=I000+NUMSYQ*DISSYQ*IINTFP
        CALL L1T13(ICORE(I000),DISSYQ,NUMSYQ,ICORE(I001),
     &             ICORE(IOFFT1(1)),ICORE(IOFFT1(2)),
     &             ICORE(IOFFL1(2)),THREEM,1,2,42,'N')
       ENDIF
C
C NEXT TWO TERMS DO NOT HAVE TO BE COMPUTED FOR RHF REFERENCE SINCE
C  THEY ARE SIMPLY THE TRANSPOSE OF THE Q(ABAB) EVALUATED ABOVE.
C
      IF(IUHF.EQ.0)GOTO 3001
C***********************************************************************
C               aB         aE
C   SOLVE FOR  Q   = SUM  T   *[-H(MEJB)] [PART 1 OF Q(BABA)].
C               iJ   M,E   iM
C***********************************************************************
C
C NEED TO CHANGE LISTW FOR W(mbej) INTERMEDIATE PICKUP!
C
       LISTT=36
       LISTW=254
       CALL HRNGDRV_1(ICORE,MXCOR,LISTW,LISTT,'WxT',MAXSIZ,1,42,
     &              .FALSE.,ICORE(IOFFT1(2)),ICORE(IOFFT1(2)),
     &               POP(1,2),POP(1,2),VRT(1,2),VRT(1,2),TWO,
     &              'BBBB',IUHF)
C***********************************************************************
C               aB           ea
C   SOLVE FOR  Q   =   SUM  T    H(meJB) [PART 2 OF Q(BABA)].
C               iJ     m,e   im
C***********************************************************************
C
C NEED TO CHANGE LISTW FOR W(mbej) INTERMEDIATE PICKUP!
C
       LISTT=35
       LISTW=257
       CALL HRNGDRV_1(ICORE,MXCOR,LISTW,LISTT,'WxT',MAXSIZ,1,42,
     &               TAU,ICORE(IOFFT1(2)),ICORE(IOFFT1(2)),      
     &               POP(1,2),POP(1,2),VRT(1,2),VRT(1,2),TWO,
     &               'BBBB',IUHF)
C
C  FOR CCSD, THE L1*T1*T1*T1 CONTRIBUTION MUST BE ADDED HERE
C
       IF(T3L) THEN
        NUMSYQ=IRPDPD(1,ISYTYP(2,42))
        DISSYQ=IRPDPD(1,ISYTYP(1,42))
        I000=1
        I001=I000+NUMSYQ*DISSYQ*IINTFP
        CALL L1T13(ICORE,DISSYQ,NUMSYQ,ICORE(I001),
     &             ICORE(IOFFT1(2)),ICORE(IOFFT1(1)),
     &             ICORE(IOFFL1(1)),THREEM,2,1,42,'T')
       ENDIF
C***********************************************************************
C
C RHF REENTRY POINT - FORM Q(ABAB) + Q(ABAB) (t) IF RHF.
C
C***********************************************************************
3001   IF(IUHF.EQ.0)THEN
        I000=1
        I010=I000+MAXSIZ*IINTFP
        CALL RHFQ(ICORE,MAXSIZ,42)
       ELSE
        CONTINUE
       ENDIF

#ifdef _KEEP258
C***********************************************************************
C               Ab         eA
C   SOLVE FOR  Q   = SUM  T   *[- H(MeJb)]  [Q(BAAB)]
C               iJ   M,e   iM
C***********************************************************************
       LISTW=258+IUHF
       LISTT=39
       CALL HRNGDRV_1(ICORE,MXCOR,LISTW,LISTT,'TxW',MAXSIZ,INCREM,43,
     &              TAU,ICORE(IOFFT1(1)),ICORE(IOFFT1(2)),      
     &              POP(1,1),POP(1,2),VRT(1,1),VRT(1,2),TWO,
     &              'BABA',IUHF)
C
C SKIP NEXT TERM FOR RHF
C
      IF(IUHF.EQ.0)GOTO 3002
C***********************************************************************
C               aB         Ea
C   SOLVE FOR  Q   = SUM  T    *[-H(mEjB]  [Q(ABBA)]
C               Ij   m,E   Im
C***********************************************************************
C
C NEED TO CHANGE LISTW FOR W INTERMEDIATE PICKUP!
C
       LISTT=38
       LISTW=258
       CALL HRNGDRV_1(ICORE,MXCOR,LISTW,LISTT,'WxT',MAXSIZ,1,43,
     &              TAU,ICORE(IOFFT1(2)),ICORE(IOFFT1(1)),      
     &              POP(1,2),POP(1,1),VRT(1,2),VRT(1,1),TWO,
     &              'ABAB',IUHF)
C
C RHF REENTRY POINT #2. FORM Q(BAAB) + Q(BAAB)(t).
C
3002   IF(IUHF.EQ.0)THEN
        I000=1
        I010=I000+MAXSIZ*IINTFP
        CALL RHFQ(ICORE,MAXSIZ,43)
       ELSE
        CONTINUE
       ENDIF
#endif 
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
       IF(I030.GT.MAXCOR)CALL INSMEM('H4X2ALL',I030,MAXCOR)
       CALL SUMRNG_11(ICORE(I000),ICORE(I010),ICORE(I020),NSZAB,NSZAB,
     &                ISCRSZ)
C
C NOW CONVERT THESE TO Ab-Ij INCREMENTS AND AUGMENT THE T2 INCREMENT LIST.
C
       I030=I020+ISCRSZ
       IF(I030.GT.MAXCOR)CALL INSMEM('H4X2ALL',I030,MAXCOR)
       CALL SST02I(ICORE(I000),ICORE(I010),NSZAB,NSZAB,ICORE(I020),
     &             'AABB')
C
C NOW SUM THIS WITH THE EXISTING AB INCREMENT AND OVERWRITE IT.
C
#ifdef _DEBUG_LVL0
CSSS       call dscal(NSZAB,0.0D0,ICORE(I010),1)
       Call checksum("@-H4X2ALL_4PDCC_1:",ICORE(I010),NSZAB)
#endif 
CSSS       IF (ISPAR .AND. COULOMB) THEN
           CALL PUTALL(ICORE(I010),NSZAB,1,216)
           CALL SUMSYM(ICORE(I010),ICORE(I000),NSZAB,116)
CSSS       ELSE
CSSS          CALL SUMSYM(ICORE(I010),ICORE(I000),NSZAB,116)
CSSS       ENDIF 

      ENDIF
      RETURN
      END
