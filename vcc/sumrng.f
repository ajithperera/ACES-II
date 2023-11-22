
      SUBROUTINE SUMRNG(T2AIBJ,T2AJBI,ISCR,NSZ1,NSZ2,ISCRSZ)
C
C FORMS THE AB T2 RING INCREMENT IN THE ORDER AI-bj.
C
C DURING THE CALCULATION, LIST 42 HOLDS ONE CONTRIBUTION TO 
C  THE INCREMENT WHICH IS ORDERED AI-bj, WHILE LIST 43 HOLDS
C  ANOTHER CONTRIBUTION WHICH IS ORDERED Aj-bI.  THIS ROUTINE
C  RESORTS THE LATTER AND OVERWRITES LIST 42 WITH THE SUM OF 
C  THE TWO LISTS.
C
C    INPUT:  
C
C         NSZ1   - THE TOTAL SIZE OF AN AI-bj SYMMETRY PACKED LIST.
C         NSZ2   - THE TOTAL SIZE OF AN Aj-bI SYMMETRY PACKED LIST.
C         T2AIBJ - THIS ARRAY HOLDS THE LIST 43 INCREMENT AFTER
C                  IT HAS BEEN REORDERED TO AI-bj.  
C                  (SIZE: NSZ1).
C         T2AJBI - THIS ARRAY IS USED TWICE, FOR TWO DIFFERENT PURPOSES.
C                  FIRST, IT HOLDS THE Aj-bI INCREMENTS FROM LIST 43 
C                  AND IS LATER USED TO HOLD THE AI-bj INCREMENTS FROM
C                  LIST 42.  CONSEQUENTLY, THE SPACE ALLOCATED TO THIS
C                  ARRAY IN THE CALLING ROUTINE MUST BE MAX(NSZ1,NSZ2).
C         ISCR   - SCRATCH ROUTINE USED TO HOLD SYMMETRY VECTORS IN
C                  SSTRNG.  THE SIZE ALLOCATED MUST BE 
C                  (NVRTA+NVRTB)*(NOCCA+NOCCB).
C         ISCRSZ - (NVRTA+NVRTB)*(NOCCA+NOCCB).
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION T2AIBJ(NSZ1),T2AJBI(NSZ2),ONE
      DIMENSION ISCR(ISCRSZ)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      DATA ONE /1.0/
      CALL GETALL(T2AJBI,NSZ2,1,43)
      CALL SSTRNG(T2AJBI,T2AIBJ,NSZ2,NSZ1,ISCR,'ABBA')
      CALL GETALL(T2AJBI,NSZ1,1,42)
      CALL SAXPY(NSZ1,ONE,T2AJBI,1,T2AIBJ,1)
      CALL VMINUS(T2AIBJ,NSZ1)
      CALL PUTALL(T2AIBJ,NSZ1,1,42)
      RETURN
      END