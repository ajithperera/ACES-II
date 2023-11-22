
      SUBROUTINE TPDIJKA6(ICORE,MAXCOR,IUHF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,DISSYG,DISSYH,DISSYH1,DISSYH2,POP,VRT
      DIMENSION ICORE(MAXCOR)
      DIMENSION I0T(2)
C
C CALCULATION OF THE SIXTH IJKA CONTRIBUTION TO
C THE EOM-CCSD TWO-PARTICLE DENSITY MATRIX
C
C  G(IJ,KA) = - 1/2 P(IJ) H(KF,IA) T(J,F)
C
CEND 
C
C CODED JG SEPTEMBER/93
C
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA AZERO,HALFM,ONE /0.D0,-0.5D0,1.0D0/

C READ IN T1 AND R1 AMPLITUDES
C
      I0T(1)=1
      I0T(2)=I0T(1)+IRPDPD(1,9)*IINTFP*IUHF
      ISTART=I0T(2)+IRPDPD(1,10)*IINTFP
C
      CALL GETLST(ICORE(I0T(1)),1,1,1,1,90)
      IF(IUHF.EQ.1) THEN
       CALL GETLST(ICORE(I0T(2)),1,1,1,2,90)
      ENDIF
C
C LOOP OVER SPIN CASES
C
      IF(IUHF.NE.0) THEN
C
       DO 1000 ISPIN=1,1+IUHF
C
C AAAA AND BBBB SPIN CASES
C
        LISTH=53+ISPIN
        LISTG=106+ISPIN
C
        DO 100 IRREP=1,NIRREP
C
         NUMSYH=IRPDPD(IRREP,ISYTYP(2,LISTH))
         DISSYH=IRPDPD(IRREP,ISYTYP(1,LISTH))
         NUMSYG=IRPDPD(IRREP,ISYTYP(2,LISTG))
         DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
         NOCCSQ=IRPDPD(IRREP,20+ISPIN)
C
         I000=ISTART
         I010=I000+IINTFP*MAX(NUMSYG*NOCCSQ,NUMSYG*DISSYG)
         ITMP=I010+IINTFP*MAX(NUMSYH*DISSYH,NUMSYG*DISSYG)
         IEND=ITMP+3*IINTFP*MAX(NUMSYH,DISSYH,NUMSYG,DISSYG)
         IF(IEND.GE.MAXCOR) CALL INSMEM('TPDIJKA6',IEND,MAXCOR)
C
C GET TRANSPOSED H INTERMEDIATES

         CALL GETTRN(ICORE(I010),ICORE(ITMP),DISSYH,
     &                NUMSYH,1,IRREP,LISTH)
C
C ZERO G ARRAY
C
         CALL ZERO(ICORE(I000),NUMSYG*NOCCSQ)
C
C TRANSPOSE THE LAST TWO INDICES
C
         CALL SYMTR1(IRREP,VRT(1,ISPIN),POP(1,ISPIN),
     &               NUMSYH,ICORE(I010),ICORE(ITMP),
     &               ICORE(ITMP+IINTFP*NUMSYH),
     &               ICORE(ITMP+2*IINTFP*NUMSYH))
C
C  LOOP OVER IRREPS OF F AND MULTIPLY
C
         IOFFH=0
         IOFFG=0
         IOFFT=0
C
         DO 101 IRREPJ=1,NIRREP
C
C  GET POPULATIONS
C
          NVRTJ=VRT(IRREPJ,ISPIN)
          NOCCJ=POP(IRREPJ,ISPIN)
          IRREPI=DIRPRD(IRREPJ,IRREP)
          NOCCI=POP(IRREPI,ISPIN)
C
C IF ONE OF THE POPULATIONS IS ZERO, SKIP MULTIPLICATION
C
          IF(MIN(NVRTJ,NOCCI,NOCCJ).NE.0) THEN
C
           CALL XGEMM('N','N',NUMSYH*NOCCI,NOCCJ,NVRTJ,HALFM,
     &                ICORE(I010+IOFFH),NUMSYH*NOCCI,
     &                ICORE(I0T(ISPIN)+IOFFT),NVRTJ,AZERO,
     &                ICORE(I000+IOFFG),NUMSYG*NOCCI)
C
          ENDIF
C
C UPDATE OFFSETS
C
          IOFFT=IOFFT+IINTFP*NOCCJ*NVRTJ
          IOFFG=IOFFG+IINTFP*NUMSYG*NOCCI*NOCCJ
          IOFFH=IOFFH+IINTFP*NUMSYH*NOCCI*NVRTJ
C
101      CONTINUE
C
C ANTISYMMETRIZE G IN THE LAST TWO INDICES    : A,K ; I<J
C
         CALL ASSYM(IRREP,POP(1,ISPIN),NUMSYG,NUMSYG,ICORE(I010),
     &              ICORE(I000))
C
C TRANSPOSE THE WHOLE G MATRIX : --> I<J ; A,K
C
         CALL TRANSP(ICORE(I010),ICORE(I000),DISSYG,NUMSYG)
C
C  TRANPOSE LAST TWO INDICES IN G --> I<J ; K,A
C
         CALL SYMTR1(IRREP,VRT(1,ISPIN),POP(1,ISPIN),DISSYG,
     &               ICORE(I000),ICORE(ITMP),
     &               ICORE(ITMP+IINTFP*DISSYG),
     &               ICORE(ITMP+2*IINTFP*DISSYG))
C
C  SAVE G ON DISK
C
          call checksum('tpdijka6',icore(i010),numsyg*dissyg)
         CALL GETLST(ICORE(I010),1,NUMSYG,1,IRREP,LISTG)
         CALL SAXPY (NUMSYG*DISSYG,ONE,ICORE(I010),1,ICORE(I000),1)
         CALL PUTLST(ICORE(I000),1,NUMSYG,1,IRREP,LISTG)
C
100     CONTINUE
C
1000   CONTINUE
C
      ENDIF
C
C ABAB AND BABA SPIN CASES
C
      DO 2000 ISPIN=1,1+IUHF
C
       LISTG=111-ISPIN
       LISTH1=55+ISPIN
       LISTH2=59-ISPIN+IUHF
C
       DO 1100 IRREP=1,NIRREP
C
        NUMSYH1=IRPDPD(IRREP,ISYTYP(2,LISTH1))
        DISSYH1=IRPDPD(IRREP,ISYTYP(1,LISTH1))
        NUMSYH2=IRPDPD(IRREP,ISYTYP(2,LISTH2))
        DISSYH2=IRPDPD(IRREP,ISYTYP(1,LISTH2))
        NUMSYG=IRPDPD(IRREP,ISYTYP(2,LISTG))
        DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
C
        I000=ISTART
        I010=I000+IINTFP*NUMSYG*DISSYG
        ITMP=I010+IINTFP*MAX(DISSYG*NUMSYG,DISSYH1*NUMSYH1,
     &                       DISSYH2*NUMSYH2)
        IEND=ITMP+3*IINTFP*MAX(DISSYG,NUMSYG,DISSYH1,DISSYH2,
     &                         NUMSYH1,NUMSYH2)
C
        IF(IEND.GE.MAXCOR) CALL INSMEM('TPDIJKA6',IEND,MAXCOR)
C
        CALL GETTRN(ICORE(I010),ICORE(ITMP),DISSYH1,NUMSYH1,
     &              1,IRREP,LISTH1)

        CALL ZERO(ICORE(I000),NUMSYG*DISSYG)

        CALL SYMTR1(IRREP,VRT(1,ISPIN),POP(1,3-ISPIN),
     &              NUMSYH1,ICORE(I010),ICORE(ITMP),
     &              ICORE(ITMP+IINTFP*NUMSYH1),
     &              ICORE(ITMP+2*IINTFP*NUMSYH1))
C
C  LOOP OVER IRREPS OF F AND MULTIPLY WITH T1A
C
        IOFFH=0
        IOFFG=0
        IOFFT=0
C
        DO 1101 IRREPJ=1,NIRREP
C
C  GET POPULATIONS
C
         NVRTJ=VRT(IRREPJ,ISPIN)
         NOCCJ=POP(IRREPJ,ISPIN)
         IRREPI=DIRPRD(IRREPJ,IRREP)
         NOCCI=POP(IRREPI,3-ISPIN)
C
C IF ONE OF THE POPULATIONS IS ZERO, SKIP MULTIPLICATION
C
         IF(MIN(NVRTJ,NOCCI,NOCCJ).NE.0) THEN
C
          CALL XGEMM('N','N',NUMSYH1*NOCCI,NOCCJ,NVRTJ,HALFM,
     &               ICORE(I010+IOFFH),NUMSYH1*NOCCI,
     &               ICORE(I0T(ISPIN)+IOFFT),NVRTJ,AZERO,
     &               ICORE(I000+IOFFG),NUMSYG*NOCCI)
         ENDIF
C
C UPDATE OFFSETS
C
         IOFFT=IOFFT+IINTFP*NOCCJ*NVRTJ
         IOFFG=IOFFG+IINTFP*NUMSYG*NOCCI*NOCCJ
         IOFFH=IOFFH+IINTFP*NUMSYH1*NOCCI*NVRTJ
C
1101    CONTINUE
        call checksum('tp6 in g',icore(i000),numsyg*dissyg)
        call checksum('tp6 inh1',icore(i010),numsyh1*dissyh1)
C
C TRANSPOSE THE G LIST FOR SECOND CONTRACTION :    a,K ; j,I --> a,K ; I,j
C
        CALL SYMTR1(IRREP,POP(1,3-ISPIN),POP(1,ISPIN),NUMSYG,
     &              ICORE(I000),ICORE(ITMP),
     &              ICORE(ITMP+IINTFP*NUMSYG),
     &              ICORE(ITMP+2*IINTFP*NUMSYG))
C
        CALL GETTRN(ICORE(I010),ICORE(ITMP),DISSYH2,NUMSYH2,
     &              1,IRREP,LISTH2)
C
C TRANPOSE THE LAST TWO INDICES OF H4 : a,K ; f,I  ---> a,K ; I,f
C
        CALL SYMTR1(IRREP,VRT(1,3-ISPIN),POP(1,ISPIN),NUMSYH2,
     &              ICORE(I010),ICORE(ITMP),
     &              ICORE(ITMP+IINTFP*NUMSYH2),
     &              ICORE(ITMP+2*IINTFP*NUMSYH2))
C
C PERFORM MULTIPLICATION WITH T1B
C
C
        IOFFG=0
        IOFFH=0
        IOFFT=0
C
        DO 1200 IRREPJ=1,NIRREP
C
C  GET POPULATIONS FOR MULTIPLICATION
C
         NOCCJ=POP(IRREPJ,3-ISPIN)
         NVRTJ=VRT(IRREPJ,3-ISPIN)
         IRREPI=DIRPRD(IRREP,IRREPJ)
         NOCCI=POP(IRREPI,ISPIN)
C
C  IF ONE OF THE POPULATIONS IS ZERO, SKIP THE MULTIPLICATION
C
         IF(MIN(NOCCI,NOCCJ,NVRTJ).NE.0) THEN
C
          CALL XGEMM('N','N',NUMSYH2*NOCCI,NOCCJ,NVRTJ,HALFM,
     &               ICORE(I010+IOFFH),NUMSYH2*NOCCI,
     &               ICORE(I0T(3-ISPIN)+IOFFT),NVRTJ,
     &               ONE,ICORE(I000+IOFFG),NUMSYG*NOCCI)
C
         ENDIF
C
C   UPDATE OFFSETS
C
         IOFFT=IOFFT+IINTFP*NOCCJ*NVRTJ
         IOFFG=IOFFG+IINTFP*NUMSYG*NOCCI*NOCCJ
         IOFFH=IOFFH+IINTFP*NUMSYH2*NOCCI*NVRTJ
C
1200    CONTINUE
        call checksum('tp6 in g',icore(i000),numsyg*dissyg)
        call checksum('tp6 inh2',icore(i010),numsyh2*dissyh2)
C
C TRANSPOSE THE WHOLE G MATRIX : a,K ; I,j  -- --> I,j ; a,K
C
        IF(ISPIN.EQ.2) THEN
         CALL SYMTR1(IRREP,POP(1,ISPIN),POP(1,3-ISPIN),
     &               NUMSYG,ICORE(I000),ICORE(ITMP),
     &               ICORE(ITMP+IINTFP*NUMSYG),
     &               ICORE(ITMP+2*IINTFP*NUMSYG))
        ENDIF
C  
        CALL TRANSP(ICORE(I000),ICORE(I010),DISSYG,NUMSYG)
C
C TRANPOSE LAST TWO INDICES IN G --> I,j ; K,a, BUT ONLY FOR ISPIN=1
C
        IF(ISPIN.EQ.1) THEN
         CALL SYMTR1(IRREP,VRT(1,3-ISPIN),POP(1,ISPIN),DISSYG,
     &               ICORE(I010),ICORE(ITMP),
     &               ICORE(ITMP+IINTFP*DISSYG),
     &               ICORE(ITMP+2*IINTFP*DISSYG))
        ENDIF
C
C SAVE G ON DISK
C
        call checksum('tpdijka6',icore(i010),numsyg*dissyg)
        CALL GETLST(ICORE(I000),1,NUMSYG,1,IRREP,LISTG)
        CALL SAXPY (NUMSYG*DISSYG,ONE,ICORE(I000),1,ICORE(I010),1)
        CALL PUTLST(ICORE(I010),1,NUMSYG,1,IRREP,LISTG)
C
1100   CONTINUE
C
2000  CONTINUE
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
