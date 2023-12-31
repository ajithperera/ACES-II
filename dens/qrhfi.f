




      SUBROUTINE QRHFI(AIOOA,AIVVA,AIOVA,DROOA,DRVVA,DROVA,
     &                 AIOOB,AIVVB,AIOVB,DROOB,DRVVB,DROVB,
     &                 DSOOA,DSVVA,DSOVA,DSOOB,DSVVB,DSOVB,
     &                 EVAL,SCR,MAXCOR,AIVOA,AIVOB)
C
C THIS ROUTINE FIRST CONSTRUCTS I" ACCORDING TO EQNS (17-18) OF
C    THE QRHF GRADIENT PAPER. IT SIMPLY COPIES THE UPPER TRIANGLE
C    INTO LOWER ONE
C THEN IT FORMS THE FOLLOWING CONTRIBUTION TO THE I
C    INTERMEDIATES FOR QRHF-CCSD GRADIENT CALCULATIONS.
C
C PART I.
C
C      I(P,Q) = I(P,Q) - D(P,Q) * EPSILON(Q,RHF)
C
C WHERE P,Q RUN OVER ALL INDICES WHICH ARE CONSIDERED IN THE
C Z-VECTOR EQUATIONS (THAT IS, ALL AI WHERE THE DESIGNATIONS
C REFER TO THE RHF OCCUPATIONS.
C
C THIS IS MOST CONVENIENTLY HANDLED WITH THE FULL MATRIX EXPANSION
C TRICK.
C
C PART II.
C
C WE NEED THE CONTRIBUTION :
C
C    2 * D(pq) * <pm||qn>  
C
C WHERE m,n ARE OCCUPIED IN THE RHF REFERENCE AND p,q ARE THE
C ORBITALS INVOLVED IN THE Z-VECTOR EQUATION.  THE TRICK HERE 
C IS TO DEAL WITH THE SPIN-SYMMETRIZED DENSITY MATRIX AND TO 
C ABANDON PROPER SPIN PARTITIONING OF THE I INTERMEDIATE.
C
C
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT,POPRHF,VRTRHF,POPDOC,VRTDOC
c&line delted
      DIMENSION AIOOA(*),AIVVA(*),AIOVA(*),DROOA(*),DRVVA(*),DROVA(*)
      DIMENSION AIOOB(*),AIVVB(*),AIOVB(*),DROOB(*),DRVVB(*),DROVB(*)
      DIMENSION DSOOA(*),DSVVA(*),DSOVA(*) 
      DIMENSION DSOOB(*),DSVVB(*),DSOVB(*)
      DIMENSION AIVOA(*),AIVOB(*)
      DIMENSION EVAL(*),SCR(*)
C
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/FLAGS/IFLAGS(100)
      COMMON/QRHFINF/POPRHF(8),VRTRHF(8),NOSH1(8),NOSH2(8),
     &               POPDOC(8),VRTDOC(8),NAI,N1I,NA2,
c&line mod
     &               NUMISCF,NUMASCF,ISPINP,ISPINM,IQRHF
C
      DATA ONE,ONEM /1.0D0,-1.D0/
      DATA HALF,TWO /0.5D0,2.D0/
c
      NBAS=NOCCO(1)+NVRTO(1)
c&lines delted
C
C PART I.
C
C      I(P,Q) = I(P,Q) - D(P,Q) * EPSILON(Q,RHF)
C
C WHERE P,Q RUN OVER ALL INDICES WHICH ARE CONSIDERED IN THE
C Z-VECTOR EQUATIONS (THAT IS, ALL AI WHERE THE DESIGNATIONS
C REFER TO THE RHF OCCUPATIONS.
C
C THIS IS MOST CONVENIENTLY HANDLED WITH THE FULL MATRIX EXPANSION
C TRICK.
C
C
C LOOP OVER SYMMETRY BLOCKS
C
      IOFFAI=1
      IOFF1I=1
      IOFFA2=1
      IOFFOOA=1
      IOFFOOB=1
      IOFFVOA=1
      IOFFVOB=1
      IOFFVVA=1
      IOFFVVB=1
      I000=1
      CALL GETREC(20,'JOBARC','RHFEVAL ',NBAS*IINTFP,SCR(I000))
      IOFFOCC=I000
      IOFFVRT=I000+NUMISCF
C
      DO 10 IRREP=1,NIRREP
C
       NOCCA  =POP(IRREP,1)
       NOCCB  =POP(IRREP,2)
       NVRTA  =VRT(IRREP,1)
       NVRTB  =VRT(IRREP,2)
       NOCCRHF=POPRHF(IRREP)
       NVRTRHF=VRTRHF(IRREP)
       NDOCC  =POPDOC(IRREP)
       NDVRT  =VRTDOC(IRREP)
c&lines deleted
       NSIZE  =NOCCA+NVRTA
       I010=I000+NBAS
       I020=I010+NSIZE*NSIZE
       I030=I020+NSIZE*NSIZE
       I040=I030+NSIZE*NSIZE
       I050=I040+NSIZE*NSIZE
       I060=I050+NSIZE*NSIZE
       I070=I060+NSIZE*NSIZE
C
C FORM SQUARED-OUT RELAXATION DENSITY AND I MATRIX FOR BOTH SPINS.
C ALSO FORM SPIN-SYMMETRIZED DENSITIES.
C
C    DR(alpha) HELD IN SCR(I010)
C    DR(beta)  HELD IN SCR(I020)
C    I (alpha) HELD IN SCR(I030)
C    I (beta)  HELD IN SCR(I040)
C    HALF*(DAA+DBB) HELD IN SCR(I060)
C
       CALL TRANSP(DROVA(IOFFVOA),SCR(I050),NOCCA,NVRTA)
       CALL SQROUT(DROOA(IOFFOOA),DROVA(IOFFVOA),SCR(I050),
     &             DRVVA(IOFFVVA),NOCCA,NVRTA,SCR(I010),1)
       CALL TRANSP(DROVB(IOFFVOB),SCR(I050),NOCCB,NVRTB)
       CALL SQROUT(DROOB(IOFFOOB),DROVB(IOFFVOB),SCR(I050),
     &             DRVVB(IOFFVVB),NOCCB,NVRTB,SCR(I020),1)
       CALL VADD  (SCR(I060),SCR(I010),SCR(I020),NSIZE*NSIZE,ONE)
       CALL SSCAL (NSIZE*NSIZE,HALF,SCR(I060),1)
C
       CALL VMINUS(AIVOA(IOFFVOA),NOCCA*NVRTA)
       CALL TRANSP(AIOVA(IOFFVOA),SCR(I050),NOCCA,NVRTA)
       CALL SQROUT(AIOOA(IOFFOOA),AIVOA(IOFFVOA),SCR(I050),
     &             AIVVA(IOFFVVA),NOCCA,NVRTA,SCR(I030),1)
       CALL VMINUS(AIVOA(IOFFVOA),NOCCA*NVRTA)
C
       CALL VMINUS(AIVOB(IOFFVOB),NOCCB*NVRTB)
       CALL TRANSP(AIOVB(IOFFVOB),SCR(I050),NOCCB,NVRTB)
       CALL SQROUT(AIOOB(IOFFOOB),AIVOB(IOFFVOB),SCR(I050),
     &             AIVVB(IOFFVVB),NOCCB,NVRTB,SCR(I040),1)
       CALL VMINUS(AIVOB(IOFFVOB),NOCCB*NVRTB)
C
C   NOW COPY UPPER TRIANGLE INTO LOWER
C   THIS IS NECESSARY FOR ORBITALS FOR WHICH ORBITAL RELAXATION
C   HAS BEEN CALCULATED. FOR OTHERS IT DOES NOT DO ANYTHING SINCE
C   THOSE PARTS OF THE MATRIX IS ZERO
C
       DO 100 ISIZE1=1,NSIZE
          DO 101 ISIZE2=ISIZE1+1,NSIZE
             IOFF2=(ISIZE1-1)*NSIZE+ISIZE2-1
             IOFF1=(ISIZE2-1)*NSIZE+ISIZE1-1
             SCR(I030+IOFF2)=SCR(I030+IOFF1)
             SCR(I040+IOFF2)=SCR(I040+IOFF1)
 101      CONTINUE
100    CONTINUE
C
C NOW FILL "FULL" EIGENVALUE VECTOR FOR THIS IRREP
C
       CALL SCOPY(NOCCRHF,SCR(IOFFOCC),1,SCR(I050),1)
       CALL SCOPY(NVRTRHF,SCR(IOFFVRT),1,SCR(I050+NOCCRHF),1)
C
C NOW FORM I(P,Q) = I(P,Q) - D(P,Q)*E(Q) FOR BOTH SPIN CASES
C
       IOFF=0
       DO 200 I=1,NSIZE
        CALL SSCAL(NSIZE,SCR(I050-1+I),SCR(I010+IOFF),1)
        CALL SSCAL(NSIZE,SCR(I050-1+I),SCR(I020+IOFF),1)
        IOFF=IOFF+NSIZE
200    CONTINUE
       CALL SAXPY(NSIZE*NSIZE,ONEM,SCR(I010),1,SCR(I030),1)
       CALL SAXPY(NSIZE*NSIZE,ONEM,SCR(I020),1,SCR(I040),1)
C
C SYMMETRIZE
C
       CALL SQUEZ2(SCR(I030),SCR(I050),NSIZE)
       CALL EXPND2(SCR(I050),SCR(I030),NSIZE)
       CALL SQUEZ2(SCR(I040),SCR(I050),NSIZE)
       CALL EXPND2(SCR(I050),SCR(I040),NSIZE)
C
C NOW COPY I MATRIX PIECES BACK OUT OF FULL MATRIX
C
       CALL SQROUT(AIOOA(IOFFOOA),AIOVA(IOFFVOA),SCR(I050),
     &             AIVVA(IOFFVVA),NOCCA,NVRTA,SCR(I030),0)
       CALL SQROUT(AIOOB(IOFFOOB),AIOVB(IOFFVOB),SCR(I050),
     &             AIVVB(IOFFVVB),NOCCB,NVRTB,SCR(I040),0)
C
C NOW COPY SYMMETRIZED DENSITY TO TARGET POSITIONS
C
      CALL SQROUT(DSOOA(IOFFOOA),DSOVA(IOFFVOA),SCR(I050),
     &            DSVVA(IOFFVVA),NOCCA,NVRTA,SCR(I060),0) 
      CALL SQROUT(DSOOB(IOFFOOB),DSOVB(IOFFVOB),SCR(I050),
     &            DSVVB(IOFFVVB),NOCCB,NVRTB,SCR(I060),0) 
C
      IOFFOCC=IOFFOCC+POPRHF(IRREP)
      IOFFVRT=IOFFVRT+VRTRHF(IRREP)
      IOFFOOA=IOFFOOA+NOCCA*NOCCA
      IOFFOOB=IOFFOOB+NOCCB*NOCCB
      IOFFVOA=IOFFVOA+NVRTA*NOCCA
      IOFFVOB=IOFFVOB+NVRTB*NOCCB
      IOFFVVA=IOFFVVA+NVRTA*NVRTA
      IOFFVVB=IOFFVVB+NVRTB*NVRTB
10    CONTINUE
C
C PART II.  A FANCY TRICK.
C
C WE NEED THE CONTRIBUTION :
C
C    2 * D(pq) * <pm||qn>  
C
C WHERE m,n ARE OCCUPIED IN THE RHF REFERENCE AND p,q ARE THE
C ORBITALS INVOLVED IN THE Z-VECTOR EQUATION.  THE TRICK HERE 
C IS TO DEAL WITH THE SPIN-SYMMETRIZED DENSITY MATRIX AND TO 
C ABANDON PROPER SPIN PARTITIONING OF THE I INTERMEDIATE.
C
C THE CONTRIBUTIONS ARE :
C
C  SUM D(EM) <EI||MJ> + SUM D(MN) <MI||NJ>   [QRHF_M]
C  E,M                  M,N
C
C  SUM D(EM) <EI||MJ> + SUM D(MN) <MI||NJ>   [QRHF_P]
C  E,M                  M,N
C
C
C  THE SPIN ORBITAL EQUATIONS FOR THE TWO TERMS ARE
C
C I(IJ) = I(IJ) + D(EM) <EI||MJ> + D(em) <eI|mJ>
C I(ij) = I(ij) + D(EM) <Ei|Mj> + D(em) <EI||MJ>
C
C I(IJ) = I(IJ) + D(MN) <MI||NJ> + D(mn) <MI|nJ>
C I(ij) = I(ij) + D(MN) <Mi|Nj>  + D(mn) <Mi||nj>
C
C HERE, IJ ARE OCCUPIED IN THE RHF REFERENCE, AND ALL OTHER LABELS ARE
C EXPRESSED IN THE SPIN-ORBITAL BASIS.  SINCE THE DENSITY IS SYMMETRIZED,
C BOTH ALPHA AND BETA CONTRIBUTIONS TO THE I INTERMEDIATES ARE THE SAME, AND
C WE WILL DEAL ONLY WITH ONE SPIN CASE, SINCE THE OCCUPATION SCHEME FOR
C BETA OR ALPHA HAS TO BE THE SAME AS THE RHF CASE.
C
C CALCULATE SOME OFFSETS. ADDRESSING IN SAXPY CALLS ASSUMES 
C THAT ALPHA AND BETA SPIN CASES OF I INTERMEDIATES ARE CONTIGUOUS
C
c&lines modif
       IOFF=(IQRHF-1)*NF1(1)
       NLEN=NF1(IQRHF)
c&end modif
C
C DO FIRST TERM
C
       I000=1
c&line corrected
       I010=I000+NF1(1)+NF1(2)
       I020=I010+NT(1)
       I030=I020+NT(2)
       CALL ZERO  (SCR(I000),NF1(1)+NF1(2))
       CALL SYMTRA(1,VRT(1,1),POP(1,1),1,DSOVA,SCR(I010))
       CALL SYMTRA(1,VRT(1,2),POP(1,2),1,DSOVB,SCR(I020))
       CALL VOINOO(SCR(I000),SCR(I010),SCR(I030),MAXCOR,1)
       CALL SAXPY (NLEN,TWO,SCR(I000+IOFF),1,AIOOA(IOFF+1),1)
C
C DO SECOND TERM
C
       CALL ZERO  (SCR(I000),NF1(1)+NF1(2))
       CALL OOINOO(SCR(I000),DSOOA,SCR(I010),MAXCOR,1)
       CALL SAXPY (NLEN,TWO,SCR(I000+IOFF),1,AIOOA(IOFF+1),1)
C
C DO THIRD TERM
C
       CALL ZERO  (SCR(I000),NF1(1)+NF1(2))
       CALL VVINOO(SCR(I000),DSVVA,SCR(I010),MAXCOR,1)
       CALL SAXPY (NLEN,TWO,SCR(I000+IOFF),1,AIOOA(IOFF+1),1)
C
C ALL DONE. WE HAVE CORRECT TOTAL I, WRONG SPIN PARTITIONING, BUT WHO CARES
C
C
c&lines deleted
      RETURN
      END
