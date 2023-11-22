C
C THIS ROUTINE CALCULATES THE SECOND TERM OF THE VIRTUAL-VIRTUAL BLOCK
C OF THE INTERMEDIATE I
C
C THE TERM IS GIVEN BY
C
C   - 2  SUM E,F,G G(EF,BG) <EF//AG>
C
C  THE SPIN CASES ARE
C
C   AA       AAAA AAAA
C            ABAB ABAB
C   BB       BBBB BBBB
C            BABA BABA
C
C NOTE THAT 4 G(AB,CD) HAS BEEN STORED, SO THAT
C THE ACTUAL FACTOR IN THE CALCULATION IS ONEM
C 
C FOR RHF A SPIN ADAPTED VERSION IS USED SINCE
C THE AAAA INTEGRALS AND GAMMA INTERMEDIATES ARE
C NOT AVAILABLE
C
C THE FORMULA IS THERE
C
C  -2 SUM E,F,G (2 G(Ef,Bg) <Ef//Ag>
C
C               + 2 G(Ef,Bg) <Ef//Ag> - 2 G(Ef,Bg) <Fe//Ag>)
C
C = -4 SUM E,F,G G(Ef,Bg) ( 2 <Ef//Ag> - <Fe//Ag>)
C
C THIS SUBROUTINE USES EXPLICITELY SYMMETRY
C
C CODED   JULY/90    JG
C
      SUBROUTINE INTV2(AIVV,ICORE,MAXCOR,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,DISSYW,POP,VRT
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      DIMENSION ICORE(MAXCOR),AIVV(1)
      COMMON /SYM2/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON/SYMPOP2/IRPDPD(8,22)
      COMMON/SYMPOP/IRP_DM(8,22),ISYTYP(2,500),NTOT(18) 
      COMMON /CONTROL/ IPRNT,IXXX,IXXX2
      COMMON /SHIFT/ ISHIFT 
      COMMON /dropgeo/ ndrgeo 
C
      DATA ONE,ONEM /1.0D+0,-1.0D+0/
C
      MXCOR=MAXCOR
C
       DO 1000 ISPIN=1,IUHF+1
C
C      AA AND BB SPIN CASES
C
       IF(ISPIN.EQ.1) THEN
        IOFF=1
       ELSE
        IOFF=NF2AA+1
       ENDIF
       IF(IUHF.EQ.1) THEN
C
C       AA AND BB SPIN CASES
C
C
C  LISTG:   GAMMA AMPLITUDES
C  LISTW:   INTEGRALS
C  FACT:    PREFACTOR
C
        LISTW=230+ISPIN
        LISTG=130+ISPIN 
        FACT=ONEM
C
C
C LOOP OVER IR REPS OF MN BLOCK (THE SAME IRREPS AS THE AF AND EF BLOCKS C HAVE
C
       DO 100 IRREP=1,NIRREP
C
C DETERMINE LENGTH OF EXPANDED VIRTUAL-VIRTUAL BLOCK
C
        NVRT2SQ=0
        DO 110 IRREPJ=1,NIRREP
         NVRT2SQ=NVRT2SQ+VRT(IRREPJ,ISPIN)*
     &                   VRT(DIRPRD(IRREPJ,IRREP),ISPIN)
110     CONTINUE
C
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTW))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTW)) 
        IF(MIN(NUMSYT,DISSYT,NVRT2SQ,DISSYT).NE.0) THEN
         if (ndrgeo.eq.0) then 
          MAXDIS=MAXCOR/(IINTFP*(2*NVRT2SQ+DISSYT))
         else
          MAXDIS=(MAXCOR-MAX(NUMSYT,DISSYT)-1)/
     x                           (IINTFP*(2*NVRT2SQ+DISSYT))
         endif 
         I001=1
C
CJDW 10/2/95. I think the next three lines are not quite right,
C             particularly for small cases. Thus, in those cases when
C             we can fit a lot more than we need in core, MAXDIS
C             is greater than NUMSYT. As a result we will use more memory
C             than we need and we cannot add arrays.
C
c        I002=I001+IINTFP*MAXDIS*NVRT2SQ
c        I003=I002+IINTFP*MAXDIS*NVRT2SQ
c        I004=I003+IINTFP*MAXDIS*DISSYT
C
         I002=I001+IINTFP*MIN(MAXDIS,NUMSYT)*NVRT2SQ
         I003=I002+IINTFP*MIN(MAXDIS,NUMSYT)*NVRT2SQ
         I004=I003+IINTFP*MIN(MAXDIS,NUMSYT)*DISSYT
         I005=I004
C
CJDW KKB stuff
         if (ndrgeo.ne.0) I005=I004+MAX(NUMSYT,DISSYT)+1 
     &     +mod((MAX(NUMSYT,DISSYT)+1),2)
CJDW END
C
         IF(MAXDIS.NE.0 .AND. I005.LE.MAXCOR) THEN
C  
C         IN CORE VERSION
C
C  THE FULL IN CORE VERSION IS IDENTICAL TO IV1AA, SHOULD BE
C  REPLACED VERY SOON BY A CLEVER OUT CORE VERSION
C
          CALL IV2AA(ICORE(I001),ICORE(I002),ICORE(I003),MAXDIS,
     &               AIVV(IOFF),FACT,ISPIN,POP(1,ISPIN),VRT(1,ISPIN),
     &               DISSYT,NUMSYT,LISTG,LISTW,IRREP,ICORE(I004))
         ELSE
          CALL INSMEM('IV2AA',MAXDIS,0)
         ENDIF
        ELSE
        ENDIF      
100    CONTINUE
       ENDIF
C
C       AB SPIN CASE
C
       LISTW=233
       LISTG=133
       FACT=ONEM
C
C      LOOP OVER IRREPS.
C
       DO 200 IRREP=1,NIRREP
C
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTW))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTW))
        NUMSYW=NUMSYT
        DISSYW=DISSYT
        IF(MIN(NUMSYT,DISSYT).NE.0) THEN
C
         if (ndrgeo.eq.0) then
          MAXDIS=(MXCOR-3*IINTFP*MAX(NUMSYT,DISSYT)-1)/(3*IINTFP*NUMSYT)
         else 
          MAXDIS=(MXCOR 
     x     -3*IINTFP*MAX(NUMSYT,DISSYT)-MAX(NUMSYT,DISSYT)-2) 
     x             /(3*IINTFP*NUMSYT)
         endif
C
         I001=1
CJDW 10/2/95. See above comments.
c        I002=I001+IINTFP*NUMSYT*MAXDIS
c        I003=I002+IINTFP*NUMSYT*MAXDIS
c        I004=I003+IINTFP*NUMSYT*MAXDIS
         I002=I001+IINTFP*NUMSYT*MIN(MAXDIS,NUMSYT)
         I003=I002+IINTFP*NUMSYT*MIN(MAXDIS,NUMSYT)
         I004=I003+IINTFP*NUMSYT*MIN(MAXDIS,NUMSYT)
         I005=I004+3*IINTFP*MAX(DISSYT,NUMSYT)
         if (ndrgeo.ne.0) I005=I005+MAX(DISSYT,NUMSYT)+1 
     &     +mod((MAX(DISSYT,NUMSYT)+1),2)
         IF(I005.LE.MXCOR) THEN
C
C         IN CORE VERSION
C
          CALL IV2AB(ICORE(I001),ICORE(I002),ICORE(I003),MAXDIS,
     &               AIVV(IOFF),FACT,ISPIN,POP(1,ISPIN),
     &               POP(1,3-ISPIN),VRT(1,ISPIN),VRT(1,3-ISPIN),
     &               DISSYT,NUMSYT,DISSYW,NUMSYW,LISTG,LISTW,
     &               IRREP,ICORE(I004),IUHF)  
         ELSE
          CALL INSMEM('IV2AB',I005,MXCOR)
         ENDIF
        ELSE
C
C
        ENDIF
200   CONTINUE
1000  CONTINUE
      RETURN
      END