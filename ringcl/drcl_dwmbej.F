      SUBROUTINE DRCL_DWMBEJ(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE FORMS THE RING INTERMEDIATES [W(mbej)] WHICH
C  ARE APPROPRIATE FOR THE CCD LAMBDA  EQUATIONS.
C
C          W(mbej) = 2 W'(mbej) - <mb||ej>
C
C WHERE W'(mbej) IS THE INTERMEDIATE USED IN THE ENERGY CALCULATION.
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL bRedundant
      LOGICAL HBAR_4LCCSD,ADC2
      DOUBLE PRECISION ONE,TWO,ONEM
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD,RLE,EOM_SDXRCCD,EOM_SFDXRCCD
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      COMMON/EXCITE/EOM,CIS,RPA,EOM_SRCCD,EOM_SDRCCD,EOM_SFRCCD,
     &              EOM_SFDRCCD,EOM_SDXRCCD,EOM_SFDXRCCD
      COMMON /FLAGS2/IFLAGS2(500)
      COMMON /FLAGS/IFLAGS(100)
      DATA ONE /1.0/
      DATA TWO /2.0/
      DATA ONEM /-1.0/
C
C DO SPIN CASE AAAA
C
      IF (IUHF.EQ.1) THEN
      ISIZE=ISYMSZ(ISYTYP(1,54),ISYTYP(2,54))
      I000=1
      I010=I000+IINTFP*ISIZE
      I020=I010+IINTFP*ISIZE
      IF (I020.GT.MAXCOR+1) CALL INSMEM("FORMWL",(I020-1),MAXCOR)
      CALL GETALL(ICORE,ISIZE,1,54)

      CALL SSCAL(ISIZE,TWO,ICORE,1)
      IF (DRCCD) THEN
         CALL GETALL(ICORE(I010),ISIZE,1,119)
      ELSE 
         CALL GETALL(ICORE(I010),ISIZE,1,23)
      ENDIF 
      CALL SAXPY(ISIZE,ONEM,ICORE(I010),1,ICORE(I000),1)
#ifdef _DEBUG_LVL0
      call checksum("List-54: ",icore(i000),isize)
C      nrow=irpdpd(1,isytyp(1,54))
C      ncol=irpdpd(1,isytyp(2,54))
C      call output(icore(i000),1,nrow,1,ncol,nrow,ncol,1)
#endif 
      CALL PUTALL(ICORE,ISIZE,1,54)

      ENDIF 
C
C DO SPIN CASE ABAB
C
      ISIZE=ISYMSZ(ISYTYP(1,56),ISYTYP(2,56))
      I000=1
      I010=I000+IINTFP*ISIZE
      I020=I010+IINTFP*ISIZE
      IF (I020.GT.MAXCOR+1) CALL INSMEM("FORMWL",(I020-1),MAXCOR)
      CALL GETALL(ICORE,ISIZE,1,56)
      CALL SSCAL(ISIZE,TWO,ICORE,1)
      CALL GETALL(ICORE(I010),ISIZE,1,18)

C Spin-adapted Coulomb integrals for DRCCD and RHF. I am not
C sure why I needed to multiply the entire term (2W(mb,ej)-2<mb|ej>)
C by a factor of 2, but that is the only way to maintain the identity of
C closed shell RHF and UHF calcs. Well, I this multiplication no longer
C necessary. This factor must be accounted when constructing 
C L2*W(mb,ej)-> L2. This is now done in rcl_lrngdrv_r.F

      IF (DRCCD .AND. IUHF .EQ.0) THEN
          CALL SSCAL(ISIZE,2.0D0,ICORE(I010),1)
          CALL SAXPY(ISIZE,ONEM,ICORE(I010),1,ICORE(I000),1)
CSSS         CALL SSCAL(ISIZE,2.0D0,ICORE(I000),1)
      ELSE
         CALL SAXPY(ISIZE,ONEM,ICORE(I010),1,ICORE(I000),1)
      ENDIF 

#ifdef _DEBUG_LVL0
      call checksum("List-56: ",icore(i000),isize)
C      nrow=irpdpd(1,isytyp(1,56))
C      ncol=irpdpd(1,isytyp(2,56))
C      call output(icore(i000),1,nrow,1,ncol,nrow,ncol,1)
#endif 
      CALL PUTALL(ICORE,ISIZE,1,56)
C
C REMAINDER OF ROUTINE FOR UHF ONLY
C
      IF(IUHF.EQ.1)THEN
C
C DO SPIN CASE BBBB
C
      ISIZE=ISYMSZ(ISYTYP(1,55),ISYTYP(2,55))
      I000=1
      I010=I000+IINTFP*ISIZE
      I020=I010+IINTFP*ISIZE
      IF (I020.GT.MAXCOR+1) CALL INSMEM("FORMWL",(I020-1),MAXCOR)
      CALL GETALL(ICORE,ISIZE,1,55)
      CALL SSCAL(ISIZE,TWO,ICORE,1)
      IF (DRCCD) THEN
         CALL GETALL(ICORE(I010),ISIZE,1,120)
      ELSE
         CALL GETALL(ICORE(I010),ISIZE,1,24)
      ENDIF 
      CALL SAXPY(ISIZE,ONEM,ICORE(I010),1,ICORE(I000),1)
#ifdef _DEBUG_LVL0
      call checksum("List-55: ",icore(i000),isize)
#endif 
      CALL PUTALL(ICORE,ISIZE,1,55)
C
C DO SPIN CASE BABA
C
      ISIZE=ISYMSZ(ISYTYP(1,57),ISYTYP(2,57))
      I000=1
      I010=I000+IINTFP*ISIZE
      I020=I010+IINTFP*ISIZE
      IF (I020.GT.MAXCOR+1) CALL INSMEM("FORMWL",(I020-1),MAXCOR)
      CALL GETALL(ICORE,ISIZE,1,57)
      CALL SSCAL(ISIZE,TWO,ICORE,1)
      CALL GETALL(ICORE(I010),ISIZE,1,17)
      CALL SAXPY(ISIZE,ONEM,ICORE(I010),1,ICORE(I000),1)
#ifdef _DEBUG_LVL0
      call checksum("List-57: ",icore(i000),isize)
C      nrow=irpdpd(1,isytyp(1,57))
C      ncol=irpdpd(1,isytyp(2,57))
C      call output(icore(i000),1,nrow,1,ncol,nrow,ncol,1)
#endif 
      CALL PUTALL(ICORE,ISIZE,1,57)
C
      ENDIF

      RETURN
      END
