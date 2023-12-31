      SUBROUTINE RCC_SST003(WIN,WOUT,TMP,MAXCOR,NSIZIN,NSIZOT,ISCR,
     +                      SPTYPE,RERDTP,LISTIN)
C
C THIS ROUTINE ACCEPTS A SYMMETRY PACKED FOUR-INDEX LIST AND RETURNS
C   THE SAME LIST BUT WITH AN ALTERNATIVE SCHEME FOR SYMMETRY PACKING.
C
C THE LIST (A<=B;I<=J) IS PRESUMED TO BE PACKED AB-IJ.  THIS ROUTINE RETURNS
C   THE LIST AS (A,I;B,J) PACKED AS AI-BJ OR AJ-BI (THESE POSSIBILITIES
C   DIFFER BY A SIGN).
C
C *** IMPORTANT *** : THE SPACE ALLOCATED FOR BOTH WIN AND WOUT IN THE
C                     CALLING ROUTINE MUST BE EQUAL TO THE SIZE OF
C                     *WOUT*
C
C INPUT: 
C           WIN  - THE SYMMETRY PACKED AB-IJ LIST.
C         NSIZIN - THE TOTAL SIZE OF THE SYM. PACKED INPUT VECTOR.
C         NSIZOT - THE TOTAL SIZE OF THE SYM. PACKED OUTPUT VECTOR.
C         SPTYPE - THE SPIN TYPE FOR THE INPUT LIST
C
C                        'AAAA' FOR (AB-IJ)
C                        'BBBB' FOR (ab-ij)
C         RERDTP - 'AIBJ' FOR AI-BJ POSITIVE, 'AJBI' FOR AJ-BI POSITIVE.
C OUTPUT: 
C          WOUT  - THE SYMMETRY PACKED AI-bj OR bj-AI LIST.
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
      DOUBLE PRECISION TMP(MAXCOR)
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
1000   FORMAT(T3,'@SST003-F, Unknown spin type ',A4,' passed in.')
       CALL ERREX
      ENDIF
C
      IF(RERDTP.EQ.'AIBJ')THEN
       REORTP='1324'
      ELSEIF(RERDTP.EQ.'AJBI')THEN
       REORTP='1324'
      ELSE
       WRITE(6,1001)RERDTP
1001   FORMAT(T3,'@SST003-F, Unknown reorder type ',A4,' passed in.')
       CALL ERREX
      ENDIF
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "------Printing from RCC_SST003-------"
      Write(6,*) 
#endif
      CALL SSTGEN(WIN,WOUT,NSIZOT,VRT(1,ISPIN),VRT(1,ISPIN),
     &            POP(1,ISPIN),POP(1,ISPIN),ISCR,IRREPX,REORTP)

      RETURN
      END

#ifdef _NOSKIP
C
C EXPAND LIST FROM A<=B,I<=J TO AB,IJ
C
      IOFFT=1
      IOFFF=1
      DO 10 IRREPR=1,NIRREP
       IRREPL=DIRPRD(IRREPX,IRREPR)
       DISSZT=IRPDPD(IRREPL,4+ISPIN)
       NUMDST=IRPDPD(IRREPR,6+ISPIN)
       DISSZF=IRPDPD(IRREPL,18+ISPIN)
       NUMDSF=IRPDPD(IRREPR,20+ISPIN)
       I000 =1 
       I010 = I000 + DISSZF * NUMDSF
       I020 = I010 + DISSZF * NUMDSF
       I030 = I020 + DISSZF * NUMDSF
       IEND = I030
       IF (IEND .GT. MAXCOR) CALL INSMEM("RCC_SST003",IEND,MAXCOR)

#ifdef _DEBUG_LVLM
       Write(6,"(a,1x,2i2)") "T2<A<=B;I<=J> for irrep and and spin:",
     +                        irrepr,ispin
       call output(win(iofft),1,disszt,1,numdst,disszt,numdst,1)
#endif
C <A<=B;I<=J> --> <AB;I<=J>  
       CALL SYMEXP6(IRREPL,VRT(1,ISPIN),POP(1,ISPIN),DISSZF,DISSZT,
     &              NUMDST,TMP(I000),WIN(IOFFT),TMP(I020))
#ifdef _DEBUG_LVLM
       Write(6,*) "T2<AB;I<=J>"
       call output(TMP(I000),1,disszf,1,numdst,disszf,numdst,1)
#endif
C  <AB;I<=J>  --> <IJ,A<=B>

       CALL TRANSP(TMP(I000),TMP(I010),NUMDST,DISSZF)

C  <IJ;A<=B>  --> <IJ,AB>

       CALL SYMEXP4(IRREPR,POP(1,ISPIN),VRT(1,ISPIN),NUMDSF,NUMDST,
     &              DISSZF,TMP(I000),TMP(I010),TMP(I020))

C  <I<=J,A<=B>--> <A<=B;I<=J>
       CALL TRANSP(TMP(I000),WOUT(IOFFF),DISSZF,NUMDSF)

#ifdef _DEBUG_LVLM
       Write(6,*) "T2<AB|IJ>"
       call output(WOUT(IOFFF),1,disszf,1,numdsf,disszf,numdsf,1)
#endif 
       IOFFF=IOFFF+NUMDSF*DISSZF
       IOFFT=IOFFT+NUMDST*DISSZT
10    CONTINUE
      ITOTSZF=IOFFF-1
      CALL SCOPY(ITOTSZF,WOUT,1,WIN,1)

#ifdef _DEBUG_LVLM
      call rcc_check_sym(win,itotszf,1,44,"pphh","noread")
#endif
#endif
