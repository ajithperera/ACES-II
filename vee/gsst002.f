       SUBROUTINE GSST002(WIN,WOUT,NSIZIN,NSIZOT,ISCR,SPTYPE,IRREPX)
C
C THIS ROUTINE ACCEPTS A SYMMETRY PACKED FOUR-INDEX LIST AND RETURNS
C   THE SAME LIST BUT WITH AN ALTERNATIVE SCHEME FOR SYMMETRY PACKING.
C
C THE LIST (A,b;I,j) IS PRESUMED TO BE PACKED Ab-Ij.  THIS ROUTINE RETURNS
C   THE LIST PACKED IN TWO POSSIBLE WAYS (SEE BELOW).
C
C INPUT: 
C           WIN  - THE SYMMETRY PACKED AB-IJ LIST.
C         NSIZIN - THE TOTAL SIZE OF THE SYM. PACKED INPUT VECTOR.
C         NSIZOT - THE TOTAL SIZE OF THE SYM. PACKED OUTPUT VECTOR.
C         SPTYPE - THE SPIN TYPE FOR THE INPUT LIST
C
C                        'AABB' FOR (AI-bj) RETURNED
C                        'BBAA' FOR (bj-AI) RETURNED
C
C OUTPUT: 
C          WOUT  - THE SYMMETRY PACKED AI-bj OR bj-AI LIST.
C       
C SCRATCH:
C         ISCR   - SCRATCH AREA TO HOLD THE SYMMETRY VECTORS AND INVERSE
C                   SYMMETRY VECTORS WHICH ARE NEEDED. 
C                   (SIZE: NVRTA*NVRTB+NOCCA*NOCCB+NVRTA*NOCCB+
C                          NVRTA*NOCCA+NVRTB*NOCCB)
C         
CEND
      IMPLICIT INTEGER (A-Z)
      CHARACTER*4 REORTP,SPTYPE
      DOUBLE PRECISION WIN(NSIZIN),WOUT(NSIZOT),ISCR(*)
      DIMENSION IOFFTAR(8)
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
C
      REORTP='1324'
      CALL SSTGEN(WIN,WOUT,NSIZIN,VRT(1,1),VRT(1,2),POP(1,1),POP(1,2),
     &            ISCR,IRREPX,REORTP)
      IF(SPTYPE.EQ.'BBAA')THEN
C
C CALCULATE SOME OFFSETS
C
       IOFFTAR(1)=1
       DO 5 IRREPR=1,NIRREP-1
        IRREPL=DIRPRD(IRREPX,IRREPR)
        DISSIZ=IRPDPD(IRREPL,10)
        NUMDIS=IRPDPD(IRREPR,9)
        IOFFTAR(IRREPR+1)=IOFFTAR(IRREPR)+DISSIZ*NUMDIS
5      CONTINUE
       IOFF1=1
       DO 10 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPX,IRREPR)
        NUMDIS=IRPDPD(IRREPR,10)
        DISSIZ=IRPDPD(IRREPL,9)
        CALL TRANSP(WOUT(IOFF1),WIN(IOFFTAR(IRREPL)),NUMDIS,DISSIZ)
        IOFF1=IOFF1+NUMDIS*DISSIZ
10     CONTINUE
       TOTSIZ=IOFF1-1
       CALL SCOPY(TOTSIZ,WIN,1,WOUT,1)
      ENDIF
C
      RETURN
      END
