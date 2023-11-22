      SUBROUTINE TPDIJKA4(ICORE,MAXCOR,IUHF,LISTL2,LISTR1,LISTR1OFF,
     &   LISTR2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,DISSYG,DISSYR,DISSYL,POP,VRT
      DIMENSION ICORE(MAXCOR)
      DIMENSION I0T(2),I0R(2),I0G(2)
C
C CALCULATION OF THE FOURTH IJKA CONTRIBUTION TO
C THE EOM-CCSD TWO-PARTICLE DENSITY MATRIX
C
C  G(IJ,KA) = - 1/4 P(IJ) { R(IM,EF) L(EF,MK) 
C
C ALSO, ONLY THE SUM OF G(IJ,KA) AND G(IA,JK) IS STORED ON DISK
C
CEND 
C
C CODED JG SEPTEMBER/93
C
      COMMON/STATSYM/IRREPX
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA ONE /1.0D0/

C READ IN T1 AND R1 AMPLITUDES
C
      I0T(1)=1
      I0T(2)=I0T(1)+IRPDPD(1,9)*IINTFP*IUHF
      I0R(1)=I0T(2)+IRPDPD(1,10)*IINTFP
      I0R(2)=I0R(1)+IRPDPD(IRREPX,9)*IINTFP*IUHF
      I0G(1)=I0R(2)+IRPDPD(IRREPX,10)*IINTFP
      I0G(2)=I0G(1)+IRPDPD(1,21)*IINTFP*IUHF
      ISTART=I0G(2)+IRPDPD(1,22)*IINTFP
C
c      CALL GETLST(ICORE(I0R(1)),1,1,1,3,490)
      CALL GETLST(ICORE(I0R(1)),1,1,1,LISTR1OFF+1,LISTR1)
      CALL GETLST(ICORE(I0T(1)),1,1,1,1,90)
      CALL ZERO(ICORE(I0G(1)),IRPDPD(1,21))
      IF(IUHF.EQ.1) THEN
c       CALL GETLST(ICORE(I0R(2)),1,1,1,4,490)
       CALL GETLST(ICORE(I0R(2)),1,1,1,LISTR1OFF+2,LISTR1)
       CALL GETLST(ICORE(I0T(2)),1,1,1,2,90)
       CALL ZERO(ICORE(I0G(2)),IRPDPD(1,22))
      ENDIF
C
C LOOP OVER SPIN CASES
C
      DO 1000 ISPIN=1,1+IUHF
C
C AA AND BB SPIN CASES
C
       I000=I0G(ISPIN)
C
       IF(IUHF.EQ.1) THEN
C
c        LISTL=443+ISPIN
        LISTL=LISTL2-1+ISPIN
c        LISTR=460+ISPIN
        LISTR=LISTR2-1+ISPIN
C
C LOOP OVER IRREPS
C
        DO 100 IRREPR=1,NIRREP
C
         IRREPL=DIRPRD(IRREPR,IRREPX) 
C
         NUMSYR=IRPDPD(IRREPR,ISYTYP(2,LISTR))
         DISSYR=IRPDPD(IRREPL,ISYTYP(1,LISTR))
         NUMSYL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
         DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
         NOCCSQ=IRPDPD(IRREPR,20+ISPIN)
C
C  ALLOCATE MEMORY
C
         I010=ISTART
         I020=I010+IINTFP*NOCCSQ*DISSYR
         IEND=I020+IINTFP*NOCCSQ*DISSYL
         IF(IEND.GE.MAXCOR) CALL INSMEM('TPDIJKA4',IEND,MAXCOR)
C
C GET R AMPLITUDE
C
         CALL GETLST(ICORE(I010),1,NUMSYR,1,IRREPR,LISTR)
C
C GET L AMPLITUDES
C
         CALL GETLST(ICORE(I020),1,NUMSYL,1,IRREPR,LISTL)
C
C EXPAND RIGHT-HANDSIDE OF R AND L
C
         CALL SYMEXP(IRREPR,POP(1,ISPIN),DISSYR,ICORE(I010))
         CALL SYMEXP(IRREPR,POP(1,ISPIN),DISSYL,ICORE(I020))
C
         IOFFG=0
         IOFFR=0
         IOFFL=0
C
         DO 70 IRREPRR=1,NIRREP
C
          NOCCI=POP(IRREPRR,ISPIN)
          NOCCK=POP(IRREPRR,ISPIN)
          IRREPRL=DIRPRD(IRREPR,IRREPRR)
          NOCCM=POP(IRREPRL,ISPIN)
C
C PERFORM MULTIPLICATION
C
          CALL XGEMM('T','N',NOCCI,NOCCK,NOCCM*DISSYR,ONE,
     &               ICORE(I010+IOFFR),DISSYR*NOCCM,
     &               ICORE(I020+IOFFL),DISSYL*NOCCM,
     &               ONE,ICORE(I000+IOFFG),NOCCI)
C
C UPDATE OFFSETS
C
          IOFFG=IOFFG+IINTFP*NOCCI*NOCCK
          IOFFR=IOFFR+IINTFP*DISSYR*NOCCI*NOCCM
          IOFFL=IOFFL+IINTFP*DISSYL*NOCCM*NOCCK
C
70       CONTINUE
C
100     CONTINUE
C
       ENDIF
C
c       LISTL=446
       LISTL=LISTL2+2
c       LISTR=463
       LISTR=LISTR2+2
C
C LOOP OVER IRREPS
C
       DO 1100 IRREPR=1,NIRREP
C
        IRREPL=DIRPRD(IRREPR,IRREPX)
C
        NUMSYR=IRPDPD(IRREPR,ISYTYP(2,LISTR))
        DISSYR=IRPDPD(IRREPL,ISYTYP(1,LISTR))
        NUMSYL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
        DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
C
C  ALLOCATE MEMORY
C
        I010=ISTART
        I020=I010+IINTFP*DISSYR*NUMSYR
        ITMP=I020+IINTFP*DISSYL*NUMSYL
        IEND=ITMP+3*IINTFP*MAX(NUMSYR,DISSYR)
        IF(IEND.GE.MAXCOR) CALL INSMEM('TPDIJKA4',IEND,MAXCOR)
C
C GET R AMPLITUDE
C
        CALL GETLST(ICORE(I010),1,NUMSYR,1,IRREPR,LISTR)
C
C FOR RHF, SPIN ADAPT R AMPLITUDES
C
        IF(IUHF.EQ.0) CALL SPINAD1(IRREPR,POP(1,1),DISSYR,
     &                             ICORE(I010),ICORE(ITMP),
     &                             ICORE(ITMP+IINTFP*DISSYR))
C
C GET L AMPLITUDES
C
         CALL GETLST(ICORE(I020),1,NUMSYL,1,IRREPR,LISTL)
C
C FOR ISPIN = 1 AND UHF, INTERCHANGE THE TWO SLOWEST INDICES
C 
         IF(IUHF.EQ.1.AND.ISPIN.EQ.1) THEN
C
          CALL SYMTR1(IRREPR,POP(1,1),POP(1,2),DISSYR,ICORE(I010),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYR),
     &                ICORE(ITMP+2*IINTFP*DISSYR))
C
          CALL SYMTR1(IRREPR,POP(1,1),POP(1,2),DISSYL,ICORE(I020),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYL),
     &                ICORE(ITMP+2*IINTFP*DISSYL))
C
         ENDIF

         IOFFG=0
         IOFFR=0
         IOFFL=0
C
         DO 1070 IRREPRR=1,NIRREP
C
          NOCCI=POP(IRREPRR,ISPIN)
          NOCCK=POP(IRREPRR,ISPIN)
          IRREPRL=DIRPRD(IRREPR,IRREPRR)
          NOCCM=POP(IRREPRL,3-ISPIN)
C
C PERFORM MULTIPLICATION
C
          CALL XGEMM('T','N',NOCCI,NOCCK,NOCCM*DISSYR,ONE,
     &               ICORE(I010+IOFFR),DISSYR*NOCCM,
     &               ICORE(I020+IOFFL),DISSYL*NOCCM,
     &               ONE,ICORE(I000+IOFFG),NOCCI)
C
C UPDATE OFFSETS
C
          IOFFG=IOFFG+IINTFP*NOCCI*NOCCK
          IOFFR=IOFFR+IINTFP*DISSYR*NOCCI*NOCCM
          IOFFL=IOFFL+IINTFP*DISSYL*NOCCM*NOCCK
C
1070     CONTINUE
C
1100    CONTINUE
C
       call checksum('tpdijka4',icore(i0g(ispin)),nfmi(ispin))
1000   CONTINUE
C
C FORM G*T TERMS
C
       IF(IUHF.EQ.1) THEN
C
C AAAA AND BBBB SPIN CASES
C
        DO 2000 ISPIN=1,2
C
         LISTG=106+ISPIN
C
         DO 2100 IRREP=1,NIRREP
C
          NUMSYG=IRPDPD(IRREP,ISYTYP(2,LISTG))
          DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
C
          I010=ISTART
          IEND=I010+IINTFP*NUMSYG*DISSYG
          IF(IEND.GE.MAXCOR) CALL INSMEM('TPDIJKA4',IEND,MAXCOR)
          CALL GETLST(ICORE(I010),1,NUMSYG,1,IRREP,LISTG)
          CALL GTTAU(ICORE(I010),ICORE(I0G(ISPIN)),ICORE(I0T(ISPIN)),
     &              DISSYG,NUMSYG,POP(1,ISPIN),POP(1,ISPIN),
     &              VRT(1,ISPIN),1,IRREP,ISPIN)
C
          CALL PUTLST(ICORE(I010),1,NUMSYG,1,IRREP,LISTG)
C
2100     CONTINUE
C
2000    CONTINUE
C
       ENDIF
C
       DO 3000 ISPIN=1,1+IUHF
C
        LISTG=111-ISPIN
C
        DO 3100 IRREP=1,NIRREP
C
         NUMSYG=IRPDPD(IRREP,ISYTYP(2,LISTG))
         DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
C
         I010=ISTART
         IEND=I010+IINTFP*NUMSYG*DISSYG
         IF(IEND.GE.MAXCOR) CALL INSMEM('TPDIJKA4',IEND,MAXCOR)
         CALL GETLST(ICORE(I010),1,NUMSYG,1,IRREP,LISTG)
         CALL GTTAU(ICORE(I010),ICORE(I0G(ISPIN)),
     &             ICORE(I0T(3-ISPIN)),DISSYG,NUMSYG,
     &             POP(1,ISPIN),POP(1,3-ISPIN),VRT(1,3-ISPIN),
     &             1,IRREP,ISPIN+2)
C
         call checksum('tpdijka4',icore(i010),numsyg*dissyg)
         CALL PUTLST(ICORE(I010),1,NUMSYG,1,IRREP,LISTG)
C
3100    CONTINUE
C
3000   CONTINUE
C
C ALL DONE, RETURN
C
      TWO=0.D0
      if(iuhf.eq.0) then
      call checkgam1(icore,10,110,two,iuhf,2,pop)
      endif
      IF(IUHF.EQ.1) THEN
      CALL CHECKGAM(ICORE,10,110,TWO)
       CALL CHECKGAM(ICORE,7,107,TWO)
       CALL CHECKGAM(ICORE,8,108,TWO)
       CALL CHECKGAM(ICORE,9,109,TWO)
      ENDIF
C
      RETURN
      END
