      SUBROUTINE RCC_SST03I(WIN,WOUT,TMP,MEMLEFT,NSIZIN,NSIZOT,
     +                      ISCR,SPTYPE)
C
C THE INCREMENT LIST (A,J;B,I) IS PRESUMED TO BE PACKED AJ-BI.  THIS 
C   ROUTINE RETURNS THE SYMMETRIZED INCREMENT LIST AS (A<=B;I<=J) STORED
C   AS AB-IJ.  
C
C *** IMPORTANT *** : THE SPACE ALLOCATED FOR BOTH WIN AND WOUT IN THE
C                     CALLING ROUTINE MUST BE EQUAL TO THE SIZE OF
C                     *WIN*
C
C INPUT: 
C           WIN  - THE SYMMETRY PACKED AB-IJ LIST.
C         NSIZIN - THE TOTAL SIZE OF THE SYM. PACKED INPUT VECTOR.
C         NSIZOT - THE TOTAL SIZE OF THE SYM. PACKED OUTPUT VECTOR.
C         SPTYPE - THE SPIN TYPE FOR THE INPUT LIST
C
C                        'AAAA' FOR (AJ-BI)
C                        'BBBB' FOR (aj-bi)
C
C OUTPUT: 
C          WOUT  - THE SYMMETRY PACKED A<=B-I<=J OR a<=b-i<=j LIST.
C       
C SCRATCH:
C         ISCR   - SCRATCH AREA TO HOLD THE SYMMETRY VECTORS AND INVERSE
C                   SYMMETRY VECTORS WHICH ARE NEEDED. 
C                   (SIZE: NVRT*(NVRT-1)/2+NOCC*(NOCC-1)/2+NOCC*NVRT)
C         
CEND
      IMPLICIT INTEGER (A-Z)
      CHARACTER*4 REORTP,SPTYPE,RERDTP
      DOUBLE PRECISION WIN(NSIZIN),WOUT(NSIZOT),ISCR(*)
      DOUBLE PRECISION TMP(MEMLEFT)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
C
      IRREPX=1
C
      IF(SPTYPE.EQ.'AAAA')THEN
       ISPIN=1
      ELSEIF(SPTYPE.EQ.'BBBB')THEN
       ISPIN=2
      ELSE
       WRITE(6,1000)SPTYPE
1000   FORMAT(T3,'@SST03I-F, Unknown spin type ',A4,' passed in.')
       CALL ERREX
      ENDIF

C Reorder (AI,BJ) --> (AB,IJ)

      CALL SSTGEN(WIN,WOUT,NSIZIN,VRT(1,ISPIN),POP(1,ISPIN),
     &            VRT(1,ISPIN),POP(1,ISPIN),ISCR,IRREPX,'1324')
C
C#ifdef _NOSKIP
      IOFFT=1
      IOFFF=1
      DO 10 IRREPR=1,NIRREP
       IRREPL=DIRPRD(IRREPR,IRREPX)
       DISSZT=IRPDPD(IRREPL,ISPIN)
       NUMDST=IRPDPD(IRREPR,2+ISPIN)
       DISSZF=IRPDPD(IRREPL,18+ISPIN)
       NUMDSF=IRPDPD(IRREPR,20+ISPIN)
       I000 = 1
       I010 = I000 + DISSZF * NUMDSF
       I020 = I010 + DISSZF * NUMDSF
       I020 = IEND 
       IF (IEND .GE. MEMLEFT) CALL INSMEM("RCC_SST03I",IEND,MEMLEFT)

C T2(AB,IJ) ---> T2(AB,JI)
       CALL SYMTRA(IRREPR,POP(1,ISPIN),POP(1,ISPIN),DISSZF,
     &             WOUT(IOFFF),TMP(I000))

C T2(AB,JI) ---> T2(BA,JI)
       CALL SYMTRA2(IRREPL,VRT(1,ISPIN),VRT(1,ISPIN),DISSZF,
     &              NUMDSF,TMP(I000),TMP(I010))

C T2(AB,IJ) + T2(BA,IJ)
       CALL DAXPY(DISSZF*NUMDSF,1.0D0,TMP(I010),1,WOUT(IOFFF),1)
#ifdef _DEBUG_LVLM
       Write(6,*) "T2(AAAA/BBBB) increment in sst03i"
       call output(wout(iofff),1,DISSZF,1,NUMDSF,DISSZF,NUMDSF,1)
#endif 

C Convert T2(AB,IJ) --> T2(<AB;I<=J)
CSSS
CSSS       CALL SQ3SYM(IRREPR,POP(1,ISPIN),DISSZF,TMP(I000),WIN(IOFFF))
CSSS
C Convert <AB,I<=J> --> <A<=B;I<=J>
C
CSSS       CALL SQ4SYM(IRREPL,VRT(1,ISPIN),DISSZT,DISSZF,NUMDST,
CSSS     +             WOUT(IOFFF),TMP(I000))
CSSS
       IOFFF=IOFFF+NUMDSF*DISSZF
       IOFFT=IOFFT+NUMDST*DISSZT
10    CONTINUE
C#endif 
C
      RETURN
      END
