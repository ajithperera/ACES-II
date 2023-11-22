      PROGRAM VTRAN
C
C INTEGRAL TRANSFORMATION PROGRAM
C
C  THIS PROGRAM TRANSFORMS SYMMETRY-ADPATED AO INTEGRALS
C  TO THE MO BASIS EXPLOITING ABELIAN SYMMETRY. BESIDE
C  A PARTIAL TRANSFORMATION FOR MBPT(2) METHODS IT
C  INCLUDES A FULL INTEGRAL TRANSFORMATION AS IT IS
C  CURRENTLY REQUIRED FOR CC CALCULATIONS.
C
CEND
C
C  WRITTEN BY J. GAUSS,  UNIVERITY KARLSRUHE
C
C  EXPERIMENTAL VERSION MARCH/92   
C  DROPPED CORE IMPLEMENTED, SEPTEMBER/92
C  OPEN-SHELL IMPLEMENTED, NOVEMBER/92
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      LOGICAL DOALL,ABIJ,NOABCD,IGNORE_DROP_MO
C
#include "flags.h"
C
      COMMON//ICORE(1)
      COMMON/ISTART/I0,ICRSIZ
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS2/IFLAGS2(500)
      COMMON/SYMINF/NDUMMY,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/INFO2/NBASIS,NBAS(8),NMO(8),POP(8,2),VRT(8,2)
      COMMON/AOOFST/INDOCC(8,2)
      COMMON/HF2FIL/LUHF2,LUHFA,LUHFB
      COMMON/INTTOL/TOL
      COMMON/DOINTS/DOALL,ABIJ,NOABCD
      COMMON/VTINFO/NPASS1,NPASS2,NPASS3,NPASS4,
     &              NLOAD1,NLOAD2,NLOAD3,NLOAD4,
     &              NWRIT1,NWRIT2,NWRIT3,NWRIT4,
     &              NWRIT1A,NWRIT2A,NWRIT3A,NWRIT4A,
     &              NWRIT1B,NWRIT2B,NWRIT3B,NWRIT4B
      COMMON/DANGRD/ndrop0 
C
      DATA TEN/10.D0/
C
C SET UNIT NUMBER FOR HF2 FILES
C
      LUHF2=15
      LUHFA=16
      LUHFB=17
C
      NPASS1=0
      NPASS2=0
      NPASS3=0
      NPASS4=0
      NLOAD1=0
      NLOAD2=0
      NLOAD3=0
      NLOAD4=0
      NWRIT1=0
      NWRIT2=0
      NWRIT3=0
      NWRIT4=0
      NWRIT1A=0
      NWRIT2A=0
      NWRIT3A=0
      NWRIT4A=0
      NWRIT1B=0
      NWRIT2B=0
      NWRIT3B=0
      NWRIT4B=0
C
C INITIALIZE THE PROGRAM BY CALLING CRAPSI, GET IUHF FLAG
C
      CALL CRAPSI(ICORE,IUHF,0)
C
C MAXIMUM CORE AVAILABLE FOR THIS CALCULATION
C
      MAXCOR=ICRSIZ
C
      CALL SETMET(IUHF)
      IGNORE_DROP_MO = (Iflags(h_IOPPAR_calclevel) .EQ. 0 .OR. 
     &                  Iflags(h_IOPPAR_vibration) .EQ. 1 .OR.
     &                  Iflags(h_IOPPAR_eval_hess) .GT. 0)
C
C CUTOFF FOR TRANSFORMED INTEGRALS WRITTEN TO HF2 (LATER TO SORT FILES)
C
      IEXP=-IFLAGS(6)
      TOL=TEN**IEXP
C
C NUMBER OF BASIS FUNCTIONS
C
      NBASIS=NOCCO(1)+NVRTO(1)
c----------------------------------------------------------------------
c---- Get the iflag for the drop-mo analytical gradient  ---  KB ------
c---- ndrgeo = 0 :  not a drop-mo analytic gradient  ------------------
c----        = 1 :  first  run for the drop-mo analytic gradient ------
c----        = 2 :  second run for the drop-mo analytic gradient ------
c----------------------------------------------------------------------
      CALL GETREC(20,'JOBARC','NDROPGEO',1,NDRGEO)  
c
C
C GET POPULATION VECTORS FOR ALL IRREPS
C
      CALL GETREC(20,'JOBARC','SYMPOPOA',NIRREP,POP(1,1))
      CALL GETREC(20,'JOBARC','SYMPOPVA',NIRREP,VRT(1,1))
C
      IF(IUHF.EQ.0) THEN
       DO 2 IRREP=1,NIRREP
        POP(IRREP,2)=POP(IRREP,1)
        VRT(IRREP,2)=VRT(IRREP,1)
2      CONTINUE
      ELSE
       CALL GETREC(20,'JOBARC','SYMPOPOB',NIRREP,POP(1,2))
       CALL GETREC(20,'JOBARC','SYMPOPVB',NIRREP,VRT(1,2))
      ENDIF
C
C GET NUMBER OF BASIS FUNCTIONS PER IRREP
C
      DO 1 IRREP=1,NIRREP
       NBAS(IRREP)=POP(IRREP,1)+VRT(IRREP,1)
1     CONTINUE
C
C REORDER EIGENVECTORS, ...
C
C   THE SCF CODE PRODUCES THE ORDERING AS FOLLOWS :
C
C     'SCF' ORDERING        :   ALL OCCUPIED FIRST, THEN ALL VIRTUALS
C     'CORRELATED' ORDERING :   WITHIN ONE IRREP ALL OCCUPIED FIRST,
C                               THEN ALL VIRTUALS
C
      CALL GTNRB2(ICORE(I0),MAXCOR,IUHF,NBASIS,NBASIS,1)
C 
C DROP MOS AND FORM CMO AND EVAL LISTS FOR CORRELATED 
C CALCULATION. The following block of code was added
C to fix a long existed problem. For cases where the 
C nirrep is more than the twice the number basis functions,
C the code was erroneous. To alleviate this, we have now
C defined a new variable MAXSIZE and that is used for 
C memory allocation and to define the size of arrays in    
C dropvc. 02/2001 Ajith Perera   
C
      IF (2*NBASIS .LE. NIRREP) THEN
         MAXSIZE = 2*NIRREP
      ELSE
         MAXSIZE = 2*NBASIS
      ENDIF 
C  
      I000=I0
      I010=I000+NBASIS*NBASIS*IINTFP
C 
C Note the correction 
cSSS  I020=I010+MAX(2*NBASIS,NIRREP) 
C 
      I020=I010+MAXSIZE 
      I030=I020+NBASIS
      I040=I030+NBASIS
      I050=I040+NBASIS
      I060=I050+NBASIS
C
      IF(I060-I0.GE.MAXCOR) CALL INSMEM('VTRAN1',I060,MAXCOR)
C
c-----------------------------------------------------------------------
c--  For the Drop-MO analytical gradient  ----------------- KB ---------
c-----------------------------------------------------------------------
      if (ndrgeo.ne.0) then 
        CALL GETREC(20,'JOBARC','NUMDROPA',1,ndrop0)  
        if (ndrop0.ne.0) then
         ijzero=0
         CALL PUTREC(20,'JOBARC','NUMDROPA',1,ijzero)   
         if (iuhf.ne.0) CALL PUTREC(20,'JOBARC','NUMDROPB',1,ijzero)   
        endif 
      endif 
c-----------------------------------------------------------------------
c Ignore all drop core processing for scf only calculations. The drop
c core (correctly should be called as post-SCF frozen core approximation)
c should only apply for post scf electronic structure methods (scf 
c analytic second derivatives s not a post scf calculations, but needs
c transformed integrals). Ajith Perera 06/2006. 
c
      If (IGNORE_DROP_MO) Then
         CALL DROPVC(NBASIS,MAXSIZE,NDROP,IUHF,ICORE(I000),
     &               ICORE(I010), ICORE(I020),ICORE(I030),
     &               ICORE(I040),ICORE(I050),ICORE(I060))
          NDROP = 0
         NDRGEO = 0  
         CALL PUTREC(20, 'JOBARC', 'NDROPGEO', 1, 0)
         CALL PUTREC(20, 'JOBARC', 'NCOMPA  ', 1, NBASIS)
         CALL PUTREC(20, 'JOBARC', 'NBASCOMP', 1, NBASIS)
         CALL PUTREC(20, 'JOBARC', 'NDROTPOP', 1, 0)
         CALL PUTREC(20, 'JOBARC', 'NDROTVRT', 1, 0)
         CALL IZERO(ICORE(I000), NBASIS)
         CALL PUTREC(20, 'JOBARC', 'IDROPA  ', NBASIS, ICORE(I000))
         CALL GETREC(20, 'JOBARC', 'NDROPGEO', 1, NDRGEO)
c
         IF (IUHF .NE. 0) THEN
            CALL PUTREC(20, 'JOBARC', 'NCOMPB  ', 1, NBASIS)
            CALL PUTREC(20, 'JOBARC', 'IDROPB  ', NBASIS, 
     &                  ICORE(I000))
         ENDIF 
      Else
         CALL DROPVC(NBASIS,MAXSIZE,NDROP,IUHF,ICORE(I000),
     &               ICORE(I010), ICORE(I020),ICORE(I030),
     &               ICORE(I040),ICORE(I050),ICORE(I060))
      Endif
C
C NUMBER OF CORRELATED MOS
C
      NCOMP=NBASIS-NDROP
C
C FORM SYMMETRY VECTORS
C
      CALL SYMDRV(ICORE(I0),MAXCOR,NCOMP,IUHF)
c-----------------------------------------------------------------------
c--  For the second run in Drop-MO analytical gradient  ------ KB ------
c-----------------------------------------------------------------------
      if (ndrgeo.eq.2) then
        call aces_fin
        stop 
      endif
c-----------------------------------------------------------------------
C
C IF NDROP NE 0, ADJUST POP AND VRT
C
      IF(NDROP.NE.0) THEN
C
       CALL GETREC(20,'JOBARC','SYMPOPOA',NIRREP,POP(1,1))
       CALL GETREC(20,'JOBARC','SYMPOPVA',NIRREP,VRT(1,1))
       IF(IUHF.EQ.0) THEN
        DO 3 IRREP=1,NIRREP
         POP(IRREP,2)=POP(IRREP,1)
         VRT(IRREP,2)=VRT(IRREP,1)
3       CONTINUE
       ELSE
        CALL GETREC(20,'JOBARC','SYMPOPOB',NIRREP,POP(1,2))
        CALL GETREC(20,'JOBARC','SYMPOPVB',NIRREP,VRT(1,2))
       ENDIF
C
      ENDIF
C
C  FILL NMO VECTOR WHICH CONTAINS TOTAL NUMBER OF CORRELATED ORBITALS
C
      DO 4 IRREP=1,NIRREP
       NMO(IRREP)=POP(IRREP,1)+VRT(IRREP,1)
4     CONTINUE
C
C CURRENT VERSION ONLY, CREATE HF2 FILE AND WRITE PREAMBLE
C
      CALL INIHF2(ICORE(I0),NBASIS,IUHF)
C
C ALLOCATE MEMORY FOR EIGEN VECTORS AND GET THEM FROM DISK
C
      ICMO=I0
      JSTART=ICMO+NBASIS*NCOMP*IINTFP*(IUHF+1)
      MXCOR=MAXCOR-IINTFP*NBASIS*NCOMP*(IUHF+1)
C
#ifdef _DEBUG
      Write(6,*)
      Write(6,*) "@-VTRAN The value of JSTART at start: ", JSTART
#endif

      IF(JSTART-I0 .GE.MAXCOR) CALL INSMEM('VTRAN2',JSTART,MAXCOR)
C
      CALL PROCEVEC(ICORE(ICMO),MAXCOR,IUHF,NBASIS,NCOMP)

      CALL TRIIII(ICORE(ICMO),ICORE(JSTART),MXCOR,IUHF)
      IF (NIRREP.GT.1) THEN
         IF (NOABCD) THEN
            CALL TRIIJJ(ICORE(ICMO),ICORE(JSTART),MXCOR,IUHF)
            CALL TRIJIJ(ICORE(ICMO),ICORE(JSTART),MXCOR,IUHF)
            IF (NIRREP.GT.2) THEN
               CALL TRIJKL(ICORE(ICMO),ICORE(JSTART),MXCOR,IUHF)
            END IF
         ELSE
            CALL TFIIJJ(ICORE(ICMO),ICORE(JSTART),MXCOR,IUHF)
            CALL TFIJIJ(ICORE(ICMO),ICORE(JSTART),MXCOR,IUHF)
            IF (NIRREP.GT.2) THEN
               CALL TFIJKL(ICORE(ICMO),ICORE(JSTART),MXCOR,IUHF)
            END IF
         END IF
      END IF
C
C WRITE FINAL BUFFER TO HF2 FILE
C
#ifdef _DEBUG
      Write(6,*) 
      Write(6,*) "@-VTRAN The value of JSTART at end: ", JSTART
#endif
      ILNBUF=600
      IBUF1=JSTART
      IBUF2=JSTART+ILNBUF*IINTFP 
      INUT=-1
C
      CALL FINAL(LUHF2,ICORE(IBUF1),ICORE(IBUF2),ILNBUF,INUT)
      IF(IUHF.EQ.1) THEN
       CALL FINAL(LUHFA,ICORE(IBUF1),ICORE(IBUF2),ILNBUF,INUT)
       CALL FINAL(LUHFB,ICORE(IBUF1),ICORE(IBUF2),ILNBUF,INUT)
      ENDIF
C
C CLOSE HF2 FILE
C
      CLOSE(UNIT=LUHF2,STATUS='KEEP')
      IF(IUHF.NE.0) THEN
       CLOSE(UNIT=LUHFA,STATUS='KEEP')
       CLOSE(UNIT=LUHFB,STATUS='KEEP')
      ENDIF
C
C WRITE TRANSFORMATION INFO
C
      IF(IFLAGS(1).GE.10)THEN
C
       write(*,9100)
9100   format('  Transformation of IIII integrals :')
       if(npass1.eq.1) then
        write(*,9101) npass1
       else if(npass1.le.9) then
        write(*,9111) npass1
       else 
        write(*,9121) npass1
       endif
9101   format (i3,' pass through the ',
     &         'AO integral file was required.')
9111   format (i3,' passes through the ',
     &         'AO integral file were required.')
9121   format (i4,' passes through the ',
     &         'AO integral file were required.')
       write(*,9102) nload1
9102   format ('  ',I10,' AO integrals were read.')
       if(iuhf.eq.0) then
        write(*,9103) nwrit1
       else
        write(*,9104) nwrit1a
        write(*,9105) nwrit1b
        write(*,9106) nwrit1
       endif
9103   format ('  ',I10,' MO integrals were written to HF2.')
9104   format ('  ',I10,' MO integrals (Spin case AAAA) were',
     &         ' written to HF2AA.')
9105   format ('  ',I10,' MO integrals (Spin case BBBB) were',
     &         ' written to HF2BB.')
9106   format ('  ',I10,' MO integrals (Spin case AABB) were',
     &         ' written to HF2AB.')
       if(nirrep.gt.1) then
        write(*,9200)
9200    format('  Transformation of IIJJ integrals :')
        if(npass2.eq.1) then
         write(*,9101) npass2
        else if(npass2.le.9) then
         write(*,9111) npass2
        else 
         write(*,9121) npass2
        endif
        write(*,9202) nload2
9202    format ('  ',I10,' AO integrals were read.')
        if(iuhf.eq.0) then
         write(*,9103) nwrit2
        else
         write(*,9104) nwrit2a
         write(*,9105) nwrit2b
         write(*,9106) nwrit2
        endif
        write(*,9300)
9300    format('  Transformation of IJIJ integrals :')
        if(npass3.eq.1) then
         write(*,9101) npass3
        else if(npass3.le.9) then
         write(*,9111) npass3
        else 
         write(*,9121) npass3
        endif
        write(*,9302) nload3
9302    format ('  ',I10,' AO integrals were read.')
        if(iuhf.eq.0) then
         write(*,9103) nwrit3
        else
         write(*,9104) nwrit3a
         write(*,9105) nwrit3b
         write(*,9106) nwrit3
        endif
        if(nirrep.gt.2) then
         write(*,9400)
9400     format('  Transformation of IJKL integrals :')
         if(npass4.eq.1) then
          write(*,9101) npass4
         else if(npass4.le.9) then
          write(*,9111) npass4
         else 
          write(*,9121) npass4
         endif
         write(*,9402) nload4
9402     format ('  ',I10,' AO integrals were read.')
         if(iuhf.eq.0) then
          write(*,9103) nwrit4
         else
          write(*,9104) nwrit4a
          write(*,9105) nwrit4b
          write(*,9106) nwrit4
         endif
        endif
       endif
      ENDIF
C
      I000=I0
      I010=I000+IINTFP*NBASIS
      IF(IFLAGS(1).GE.10)THEN
       CALL ORBSUM(ICORE(I000),ICORE(I010),NCOMP,IUHF)
      ENDIF
C
C CALL ACES_FIN
C
c-----------------------------------------------------------------------
c---  For Drop-MO analytical gradient  ------  Mar. 1995  ----   KB ----
c-----------------------------------------------------------------------
c---  In the first run of vtran for drop-mo analytic gradient ----------
c---  HF2 file is made for all-mo index,  ------------------------------
c---  But, others are made for drop-mo space, here  --------------------
c-----------------------------------------------------------------------
      if (ndrgeo.eq.1) then
       CALL PUTREC(20,'JOBARC','NUMDROPA',1,ndrop0)  
       if (iuhf.ne.0) CALL PUTREC(20,'JOBARC','NUMDROPB',1,ndrop0)  
       CALL DROPVC(NBASIS,MAXSIZE,NDROP,IUHF,ICORE(I000),ICORE(I010),
     &            ICORE(I020),ICORE(I030),ICORE(I040),ICORE(I050),
     &            ICORE(I060))
       NCOMP=NBASIS-NDROP
       CALL SYMDRV(ICORE(I0),MAXCOR,NCOMP,IUHF)
      endif 
c-----------------------------------------------------------------------
      call aces_fin   
C
      STOP
C        
      END
