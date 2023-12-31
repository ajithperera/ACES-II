       SUBROUTINE DDENSVV(DDVV,ICORE,MAXCOR,IUHF,ANTI)
C
C THIS ROUTINE CALCULATES THE VIRTUAL-VIRTUAL BLOCK
C OF THE DERIVATIVE OF THE RELAXED DENSITY MATRIX IN
C CORRELATION METHODS.
C
C THE FORMULAS ARE
C
C MBPT(2):
C
C DD(A,B) = 1/2 P(A,B) SUM M,N,E (DELTA T1(MN,AE)/DELTA X) T1(MN,BE)          
C
C FOR IMAGINARY PERTURBATIONS, THE FORMULA IS
C
C DD(A,B) = 1/2 SUM M,N,E (DELTA T1*(MN,AE)/DELTA X) T1(MN,BE)
C
C           + 1/2 SUM M,N,E T1*(MN,AE) (DELTA T1(MN,BE)/DELTA X)
C
C         = - 1/2 P_(A,B) SUM M,N,E (DELTA T1(MN,AE)/DELTA X) T1(MN,BE)
C
C ROHF-MBPT(2):
C
C DD(A,B) = 1/2 P(A,B) SUM M,N,E (DELTA T1(MN,AE)/DELTA X) T1(MN,BE)            
C
C          + P(A,B) SUM M (DELTA T(M,A)/DELTA X) T(M,B)                     
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
C IN THE RHF CASE EXPLICIT SPIN ADAPTED CODE IS USED
C
CEND
C
C CODED   APRIL/91 JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYL,DISSYR,POP,VRT
      LOGICAL MBPT4,CC
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      LOGICAL DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,ANTI
C
      DIMENSION ICORE(MAXCOR),DDVV(1)
C
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1(2),NF2(2)
      COMMON/DSYM/IRREPX,IPERT,NDTAA,NDTBB,NDF1(2),NDF2(2),
     &            IOFFIJ(8,2),IOFFAB(8,2),IOFFAI(8,2) 
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /METH/ MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      COMMON/DERIV/DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C
      DATA ONE,TWO,ONEM,TWOM /1.0D0,2.0D0,-1.0D0,-2.0D0/
C
      MBPT4=M4DQ.OR.M4SDQ.OR.M4SDTQ
      CC=CCD.OR.CCSD.OR.QCISD
      MXCOR=MAXCOR
C
      IF(CCSD.OR.QCISD.OR.M4SDQ.OR.M4SDTQ.OR.ROHF) THEN
C
C    ALLOCATE MEMORY FOR T1 AMPLITUDES
C
      I0T1A=MXCOR+1-NDTAA*IINTFP
      MXCOR=MXCOR-NDTAA*IINTFP
      CALL GETLST(ICORE(I0T1A),1,1,1,1,490)
      IF(CC) THEN
       LISTL1=190
      ELSE
       LISTL1=90
      ENDIF
      I0T2A=I0T1A-NTAA*IINTFP
      MXCOR=MXCOR-NTAA*IINTFP
      CALL GETLST(ICORE(I0T2A),1,1,2,1,LISTL1)
      IF(IUHF.EQ.0) THEN
       I0T1B=I0T1A
       I0T2B=I0T2A
      ELSE
       I0T1B=I0T2A-NDTBB*IINTFP
       MXCOR=MXCOR-NDTBB*IINTFP
       CALL GETLST(ICORE(I0T1B),1,1,1,2,490)
       I0T2B=I0T1B-NTBB*IINTFP
       MXCOR=MXCOR-NTBB*IINTFP
       CALL GETLST(ICORE(I0T2B),1,1,2,2,LISTL1)
      ENDIF
C
C   NOW PERFORM MULTIPLICATION    SUM M   T(M,A) L(M,B)
C
C   FACT IS HERE ALWAYS ONE
C
      FACT=TWO
      IF(ANTI) FACT=TWOM
C
      DO 300 ISPIN=1,IUHF+1
C
       IF(ISPIN.EQ.1) THEN
        IOFFT1=I0T1A
        IOFFT2=I0T2A
        IOFFD=1
       ELSE
        IOFFT1=I0T1B
        IOFFT2=I0T2B
        IOFFD=1+NDF2(1)
       ENDIF
C 
       DO 250 IRREPR=1,NIRREP
C
        IRREPL=DIRPRD(IRREPR,IRREPX)
        NOCCR=POP(IRREPR,ISPIN)
        NVRTR=VRT(IRREPR,ISPIN)
        NVRTL=VRT(IRREPL,ISPIN)
        IF(MIN(NVRTR,NVRTL,NOCCR).GT.0) THEN
         CALL XGEMM('N','T',NVRTL,NVRTR,NOCCR,FACT,ICORE(IOFFT1),NVRTL, 
     &              ICORE(IOFFT2),NVRTR,ONE,DDVV(IOFFD),NVRTL)
        ENDIF
        IOFFT1=IOFFT1+NOCCR*NVRTL*IINTFP
        IOFFT2=IOFFT2+NOCCR*NVRTR*IINTFP
        IOFFD=IOFFD+NVRTR*NVRTL
250    CONTINUE
300   CONTINUE
C
      ENDIF
C
       MXCOR=MAXCOR 

       DO 1000 ISPIN=1,IUHF+1
C
C      AA AND BB SPIN CASES
C
       IF(ISPIN.EQ.1) THEN
        IOFF=1
       ELSE
        IOFF=NDF2(1)+1
       ENDIF
       IF(IUHF.EQ.1) THEN
C
C       AA AND BB SPIN CASES
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
C         DD(A,B) = (DD(A,B) + DD(B,A))/2
C
C HAS TO BE CARRIED OUT FOR ALL METHODS.
C
        IF(MBPT2) THEN
         LISTT1=443+ISPIN
         LISTT2=43+ISPIN
C
C THE FACTOR IS ALWAYS TWO
C
         FACT=TWO
         IF(ANTI) FACT=TWOM
         IFLAG=2
        ELSE IF(MBPT3) THEN
c         LISTT1=43+ISPIN
c         LISTT2=60+ISPIN
         FACT=ONE
         IF(ANTI) FACT=TWOM
c         IFLAG=3
        ELSE IF(MBPT4) THEN
c         LISTT1=143+ISPIN
c         LISTT2=143+ISPIN
c         LISTT3=43+ISPIN
c         LISTT4=60+ISPIN
c         LISTT5=113+ISPIN
         FACT=ONE
         IF(ANTI) FACT=TWOM
c         IFLAG=1
        ELSE IF(CC) THEN
         LISTT1=443+ISPIN
         LISTT2=143+ISPIN
         FACT=ONE
         IF(ANTI) FACT=TWOM
         IFLAG=2
        ENDIF
C
C
C LOOP OVER IR REPS OF MN BLOCK (THE SAME IRREPS AS THE AF AND EF BLOCKS C HAVE
C
       DO 100 IRREPR=1,NIRREP
C
        IRREPL=DIRPRD(IRREPR,IRREPX)
C
C DETERMINE LENGTH OF EXPANDED VIRTUAL-VIRTUAL BLOCK
C
        NV2SQR=0
        NV2SQL=0
        DO 110 IRREPJ=1,NIRREP
         NV2SQR=NV2SQR+VRT(IRREPJ,ISPIN)*
     &                   VRT(DIRPRD(IRREPJ,IRREPR),ISPIN)
         NV2SQL=NV2SQL+VRT(IRREPJ,ISPIN)*
     &                   VRT(DIRPRD(IRREPJ,IRREPL),ISPIN)
110     CONTINUE
C
        DISSYR=IRPDPD(IRREPR,ISYTYP(1,43+ISPIN))
        NUMSYR=IRPDPD(IRREPR,ISYTYP(2,43+ISPIN)) 
        DISSYL=IRPDPD(IRREPL,ISYTYP(1,43+ISPIN))
        NUMSYL=NUMSYR
        I001=1
        I002=I001+IINTFP*NUMSYL*NV2SQL
        I003=I002+IINTFP*NUMSYR*NV2SQR
        I004=I003+IINTFP*MAX(NUMSYR,DISSYR,NUMSYL,DISSYL)
        IF(MIN(NUMSYR,DISSYL,DISSYR).NE.0) THEN
         IF(I004.LT.MXCOR) THEN
C  
C         IN CORE VERSION
C
          CALL DDVVAA(ICORE(I001),ICORE(I002),DDVV(IOFF),
     &               FACT,ISPIN,POP(1,ISPIN),VRT(1,ISPIN),DISSYL,
     &               NUMSYL,DISSYR,NUMSYR,LISTT1,LISTT2,LISTT3,
     &               IRREPL,IRREPR,IFLAG,ICORE(I003))
          IF(MBPT4) THEN
c          CALL DVVAA(ICORE(I001),ICORE(I002),ICORE(I003),DDVV(IOFF),
c     &               FACT,ISPIN,POP(1,ISPIN),VRT(1,ISPIN),
c     &               DISSYT,NUMSYT,LISTT4,LISTT5,LISTT3,IRREP,4)
          ENDIF
         ELSE
          CALL INSMEM('DDVVAA',I004,MXCOR)
         ENDIF
        ELSE
        ENDIF      
100    CONTINUE
       call checksum('ddvvaa  ',ddvv(ioff),ndf2(ispin))
       ENDIF
C
C       AB SPIN CASE
C
       IF(MBPT2) THEN
        LISTT1=446
        LISTT2=46
C
C FACTOR IS ALWAYS TWO      
C
        FACT=TWO
        IF(ANTI) FACT=TWOM
        IFLAG=2
       ELSE IF(MBPT3) THEN
c        LISTT1=46
c        LISTT2=63
        FACT=ONE
        IF(ANTI) FACT=TWOM
c        IFLAG=3
c       ELSE IF(MBPT4) THEN
c        LISTT1=146
c        LISTT2=146
c        LISTT3=46
c        LISTT4=63
c        LISTT5=116
        FACT=ONE
        IF(ANTI) FACT=TWOM
c        IFLAG=1
       ELSE IF(CC) THEN
        LISTT1=446
        LISTT2=146
        FACT=ONE
        IF(ANTI) FACT=TWOM
        IFLAG=2
       ENDIF
C
C      LOOP OVER IRREPS.
C
       DO 200 IRREPR=1,NIRREP
C
        IRREPL=DIRPRD(IRREPX,IRREPR)
C
        DISSYR=IRPDPD(IRREPR,ISYTYP(1,46))
        NUMSYR=IRPDPD(IRREPR,ISYTYP(2,46))
        DISSYL=IRPDPD(IRREPL,ISYTYP(1,46))
        NUMSYL=NUMSYR
        I001=1
        I002=I001+IINTFP*NUMSYL*DISSYL
        I003=I002+IINTFP*NUMSYR*DISSYR
        IF(MIN(NUMSYR,DISSYR,DISSYL).NE.0) THEN
         I004=I003+3*IINTFP*MAX(DISSYL,DISSYR,NUMSYR)
         IF(I004.LE.MXCOR) THEN
C
C         IN CORE VERSION
C
          CALL DDVVAB(ICORE(I001),ICORE(I002),DDVV(IOFF),
     &                FACT,ISPIN,POP(1,ISPIN),
     &                POP(1,3-ISPIN),VRT(1,ISPIN),VRT(1,3-ISPIN), 
     &                DISSYL,NUMSYL,DISSYR,NUMSYR,LISTT1,LISTT2,
     &                LISTT3,IRREPL,IRREPR,
     &                ICORE(I004),IUHF,IFLAG)  
          IF(MBPT4) THEN
c           CALL DVVAB(ICORE(I001),ICORE(I002),ICORE(I003),DDVV(IOFF),
c     &                FACT,ISPIN,POP(1,ISPIN),
c     &                POP(1,3-ISPIN),VRT(1,ISPIN),VRT(1,3-ISPIN), 
c     &                DISSYT,NUMSYT,LISTT4,LISTT5,LISTT3,IRREP,
c     &                ICORE(I004),IUHF,4)  
          ENDIF
         ELSE
          CALL INSMEM('DDVVAB',I004,MXCOR)
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
c       IF(CANON) THEN
c        CALL DENCAN(DVV(IOFF),ICORE,NF2(ISPIN),NIRREP,VRT(1,ISPIN))
c       ENDIF
C
C  ADD HERE TRIPLE CONTRIBUTIONS IN MBPT(4)
C
       IF(M4SDTQ) THEN
        CALL ERREX
       ENDIF
C
C   SYMMETRIZE DERIVATIVE OF THE DENSITY MATRIX
C    
c      call checksum('bdsymmet',ddvv(ioff),ndf2(ispin))
      CALL DSYMMET(IRREPX,DDVV(IOFF),ICORE,VRT(1,ISPIN),
     &             IOFFAB(1,ISPIN),ANTI)
c      call checksum('adsymmet',ddvv(ioff),ndf2(ispin))
C
      call checksum('ddvvab  ',ddvv(ioff),ndf2(ispin))
C
1000  CONTINUE
      RETURN
      END
