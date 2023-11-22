      SUBROUTINE MODF(ICORE,MAXCOR,IUHF,IRREPT,IRREPF,
     &                LISTT1,IOFFT1,LISTF3,IOFFF3,LISTF1, 
     &                IOFFF1,LISTF2,IOFFF2)
C
C ADDS TO THE GENERAL F(MI) AND F(AE) INTERMEDIATES THE
C FOLLOWING CONTRIBUTIONS
C
C F(a,e) = F(a,e) - 1/2 T(a,m) F(m,e)
C
C F(M,i) = F(m,i) + 1/2 T(e,i) F(m,e)
C
C THE CODE IS GENERAL IN THE SENSE THAT T AND F ARE NOT 
C RESTRICTED TO BE TOTALLY SYMMETRIC, BUT RATHER TRANSFORM
C AS IRREPT (T(a,i)), IRREPF (F(m,e)), AND IRREPX=DIRPRD
C (IRREPT,IRREPF) (F(a,e) and F(m,i)).
C
C NO CHANGE FOR IMAGINARY PERTURBATIONS REQUIRED.
C
CEND
C
C CODED JULY/94 JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
C
      DIMENSION ICORE(MAXCOR)
      DIMENSION I0T(2),I0F1(2),I0F2(2),I0F3(2)
C
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
C
      DATA ZILCH,HALF,HALFM,ONE/0.D0,0.5D0,-0.5D0,1.D0/
C
      write(*,*)
      write(*,*) ' Full Fbar is constructed for CCSD second ',
     &           'derivatives'
      write(*,*)
      IRREPX=DIRPRD(IRREPF,IRREPT)
C
C ALLOCATE MEMORY AND GET T1 AMPLITUDES AND ALL F INTERMEDIATES
C
      I0T(1)=1
      I0T(2)=I0T(1)+IINTFP*IUHF*IRPDPD(IRREPT,9)
      ISTART=I0T(2)+IINTFP*IRPDPD(IRREPT,10)
      CALL GETLST(ICORE(I0T(1)),1,1,1,1+IOFFT1,LISTT1)
      IF(IUHF.EQ.1) CALL GETLST(ICORE(I0T(2)),1,1,1,2+IOFFT1,LISTT1)
C
      I0F1(1)=ISTART
      I0F1(2)=I0F1(1)+IINTFP*IUHF*IRPDPD(IRREPX,21)
      ISTART=I0F1(2)+IINTFP*IRPDPD(IRREPX,22)
      CALL GETLST(ICORE(I0F1(1)),1,1,1,1+IOFFF1,LISTF1)
      IF(IUHF.EQ.1) CALL GETLST(ICORE(I0F1(2)),1,1,1,2+IOFFF1,LISTF1)
C
      I0F2(1)=ISTART
      I0F2(2)=I0F2(1)+IINTFP*IUHF*IRPDPD(IRREPX,19)
      ISTART=I0F2(2)+IINTFP*IRPDPD(IRREPX,20)
      CALL GETLST(ICORE(I0F2(1)),1,1,1,1+IOFFF2,LISTF2)
      IF(IUHF.EQ.1) CALL GETLST(ICORE(I0F2(2)),1,1,1,2+IOFFF2,LISTF2)
C
      I0F3(1)=ISTART
      I0F3(2)=I0F3(1)+IINTFP*IUHF*IRPDPD(IRREPF,9)
      ISTART=I0F3(2)+IINTFP*IRPDPD(IRREPF,10)
      CALL GETLST(ICORE(I0F3(1)),1,1,1,1+IOFFF3,LISTF3)
      IF(IUHF.EQ.1) CALL GETLST(ICORE(I0F3(2)),1,1,1,2+IOFFF3,LISTF3)
C
      DO 1000 ISPIN=1,1+IUHF
C
       DO 100 IRREPA=1,NIRREP
C
        IRREPE=DIRPRD(IRREPX,IRREPA)
        IRREPM=DIRPRD(IRREPT,IRREPA)
C
        NVRTA=VRT(IRREPA,ISPIN)
        NVRTE=VRT(IRREPE,ISPIN)
        NOCCM=POP(IRREPM,ISPIN)
C
        IOFFT=I0T(ISPIN)+IINTFP*(ISYMOFF(IRREPM,IRREPT,8+ISPIN)-1)
        IOFFF=I0F3(ISPIN)+IINTFP*(ISYMOFF(IRREPM,IRREPF,8+ISPIN)-1)
        IOFFZ=I0F2(ISPIN)+IINTFP*(ISYMOFF(IRREPA,IRREPX,18+ISPIN)-1)
C
        CALL XGEMM('N','T',NVRTE,NVRTA,NOCCM,HALFM,ICORE(IOFFF),NVRTE,
     &             ICORE(IOFFT),NVRTA,ONE,ICORE(IOFFZ),NVRTE)
C
100    CONTINUE
C
       DO 200 IRREPI=1,NIRREP
C
        IRREPM=DIRPRD(IRREPX,IRREPI)
        IRREPE=DIRPRD(IRREPT,IRREPI)
C
        NOCCI=POP(IRREPI,ISPIN)
        NOCCM=POP(IRREPM,ISPIN)
        NVRTE=VRT(IRREPE,ISPIN)
C
        IOFFT=I0T(ISPIN)+IINTFP*(ISYMOFF(IRREPI,IRREPT,8+ISPIN)-1)
        IOFFF=I0F3(ISPIN)+IINTFP*(ISYMOFF(IRREPM,IRREPF,8+ISPIN)-1)
        IOFFZ=I0F1(ISPIN)+IINTFP*(ISYMOFF(IRREPI,IRREPX,20+ISPIN)-1)
C
        CALL XGEMM('T','N',NOCCM,NOCCI,NVRTE,HALF,ICORE(IOFFF),NVRTE,
     &             ICORE(IOFFT),NVRTE,ONE,ICORE(IOFFZ),NOCCM)
C
200    CONTINUE
C
1000  CONTINUE
C
C SAVE MODIFIED F INTERMEDIATES ON DISK
C
      CALL PUTLST(ICORE(I0F1(1)),1,1,1,1+IOFFF1,LISTF1)
      IF(IUHF.EQ.1) CALL PUTLST(ICORE(I0F1(2)),1,1,1,2+IOFFF1,LISTF1)
      CALL PUTLST(ICORE(I0F2(1)),1,1,1,1+IOFFF2,LISTF2)
      IF(IUHF.EQ.1) CALL PUTLST(ICORE(I0F2(2)),1,1,1,2+IOFFF2,LISTF2)
C
C ALL DONE
C
      RETURN
      END