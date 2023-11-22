      SUBROUTINE DT1INT2A(ICORE,MAXCOR,IUHF,IRROMEGA,LISTW0,LISTZ0,
     &                    LISTT1)
C
C THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF T1 TO
C   T2 (T2<-T1).  THIS CODE IS GENERAL IN THE SENSE THAT
C   IT DOES NOT ASSUME THAT T (AND Z) ARE TOTALLY SYMMETRIC
C   QUANTITIES, BUT RATHER TRANSFORM AS IRROMEGA.
C
C CONTRACTION 1:
C
C     Z(ab,ij) = SUM  T(a,m) * <ij||mb>   [AAAA and BBBB]
C                 e
C
C     Z(Ab,Ij) = SUM  T(A,M) * <Ij||Mb>   [ABAB]
C                 e
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      LOGICAL RHF
C
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2),IOFFW(8),IOFFZ(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      DATA ONE  /1.0D0/
      DATA ONEM /-1.0D0/
      DATA ZILCH/0.0D0/
c
C
      CALL GETDT1(ICORE,MAXCOR,MXCOR,IUHF,IOFFT1,IRROMEGA,LISTT1)
      RHF=.FALSE.
      IF(IUHF.EQ.0)RHF=.TRUE.
C
C
C DO ALPHA-BETA SPIN CASE.  CONTAINS BOTH GENERAL AB AND 
C SPIN-ADAPTED RHF CODE.  THIS BLOCK OF CODE ALWAYS RUNS.
C
      DO 100 IRREPZR=1,NIRREP
C
C LOOP OVER KET IRREPS OF *TARGET*.  THIS IS NOT THE SAME AS THE
C  IRREPS OF THE INTEGRALS AND AMPLITUDES UNLESS IRROMEGA=1.
C
       IRREPZL=DIRPRD(IRREPZR,IRROMEGA)
       IRREPW=IRREPZR
       LISTW=LISTW0+3
       LISTZ=LISTZ0+2
       DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW))
       NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW)) 
       DISSYZ=IRPDPD(IRREPZL,ISYTYP(1,LISTZ))
       NUMDSZ=IRPDPD(IRREPZR,ISYTYP(2,LISTZ))
       MAXW=MAX(NUMDSW,DISSYW)
       MAXZ=MAX(DISSYZ,NUMDSZ)
       MAXALL=MAX(MAXW,MAXZ)
       I000=1
       I010=I000+IINTFP*MAX(DISSYZ*NUMDSZ,3*MAXALL)
       I020=I010+IINTFP*MAX(DISSYW*NUMDSW,NUMDSZ*DISSYZ,3*MAXALL)
       IF(I020.LT.MXCOR)THEN
C
C DO IN-CORE ALGORITHM
C
C
C READ W INTO W(Ij,Mb) AND REORDER TO W(Ij,bM)
C
        CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPW,LISTW)
        ITMP1=I000
        ITMP2=ITMP1+IINTFP*MAXW
        ITMP3=ITMP2+IINTFP*MAXW
        ITMP4=ITMP3+IINTFP*MAXW
        CALL SYMTR1 (IRREPW,POP(1,1),VRT(1,2),DISSYW,ICORE(I010),
     &               ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
C
C COMPUTE OFFSETS FOR W AND Z VECTORS
C
        IW=0
        IZ=0
        DO 110 IRR2=1,NIRREP
         IOFFW(IRR2)=I010+IW
         IOFFZ(IRR2)=I000+IZ
         IRR1W=DIRPRD(IRR2,IRREPW)
         IRR1Z=DIRPRD(IRR2,IRREPZL)
         IW=IW+DISSYW*VRT(IRR1W,2)*POP(IRR2,1)*IINTFP
         IZ=IZ+NUMDSZ*VRT(IRR1Z,2)*VRT(IRR2,1)*IINTFP
110     CONTINUE
C
C PERFORM MATRIX MULTIPLICATION
C
C                                         +
C         Z(Ij,bA) = SUM W(Ijb,M) * T(A,M)
C                     M
C  
        DO 120 IRREPM=1,NIRREP
         IRREPA=DIRPRD(IRREPM,IRROMEGA)
         IRREPB=DIRPRD(IRREPM,IRREPW)
         NROW=DISSYW*VRT(IRREPB,2)
         NCOL=VRT(IRREPA,1)
         NSUM=POP(IRREPM,1)
         IZ=IOFFZ(IRREPA)
         IW=IOFFW(IRREPM)
         IT=IOFFT1(IRREPM,1) 
         CALL XGEMM('N','T',NROW,NCOL,NSUM,ONEM,ICORE(IW),NROW,
     &              ICORE(IT),NCOL,ZILCH,ICORE(IZ),NROW)
120     CONTINUE
C
       ELSE
C
C OUT-OF-CORE ALGORITHM
C
        WRITE(6,*)' out-of-core AB not coded '
        call errex
C
       ENDIF
C
C Z(Ij,bA) --> Z(Ij,Ab)
C
       ITMP1=I010
       ITMP2=ITMP1+IINTFP*MAXZ
       ITMP3=ITMP2+IINTFP*MAXZ
       ITMP4=ITMP3+IINTFP*MAXZ
       CALL SYMTR1 (IRREPZL,VRT(1,2),VRT(1,1),NUMDSZ,ICORE(I000),
     &              ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
C
       IF(RHF)THEN
C
C REORDER TARGET BACK TO Z(Ab,Ij)
C
        CALL TRANSP(ICORE(I000),ICORE(I010),DISSYZ,NUMDSZ)
c YAU : old
c       CALL ICOPY(IINTFP*NUMDSZ*DISSYZ,ICORE(I010),1,ICORE(I000),1)
c YAU : new
        CALL DCOPY(NUMDSZ*DISSYZ,ICORE(I010),1,ICORE(I000),1)
c YAU : end
C
C IF RHF, SYMMETRIZE TARGET, DUMP TO DISK AND RETURN
C
C           Z(ab,ij) = Z(ab,ij) + Z(ba,ji)
C
        CALL SYMRHF3(IRREPZL,IRREPZR,VRT(1,1),POP(1,1),DISSYZ,
     &               ICORE(I000),ICORE(ITMP1),ICORE(ITMP2),
     &               ICORE(ITMP3))
        CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPZR,LISTZ)
        CALL SAXPY (NUMDSZ*DISSYZ,ONE,ICORE(I010),1,ICORE(I000),1)
        CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPZR,LISTZ)
C
       ELSE
C
C IF UHF, HAVE TO DO ANOTHER CONTRACTION FOR AB CASE
C
C        Z(Ij,Ab) = W(Ij,Am) * T(bm)
C
C
C READ W INTO W(Ab,Je)
C
        LISTW=LISTW0+2
        DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW)) 
        NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW)) 
        I020=I010+IINTFP*MAX(DISSYW*NUMDSW,3*MAXZ,NUMDSZ*DISSYZ)
        IF(I020.LT.MXCOR)THEN
         CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPW,LISTW)
C
C COMPUTE OFFSETS FOR W AND Z VECTORS
C
         IW=0
         IZ=0
         DO 140 IRR2=1,NIRREP
          IOFFW(IRR2)=I010+IW
          IOFFZ(IRR2)=I000+IZ
          IRR1W=DIRPRD(IRR2,IRREPW)
          IRR1Z=DIRPRD(IRR2,IRREPZL)
          IW=IW+DISSYW*VRT(IRR1W,1)*POP(IRR2,2)*IINTFP
          IZ=IZ+NUMDSZ*VRT(IRR1Z,1)*VRT(IRR2,2)*IINTFP
140      CONTINUE
C
C PERFORM MATRIX MULTIPLICATION
C                                         +
C         Z(Ij,Ab) = SUM W(Ij,Am) * T(b,m)
C                     e
C  
         DO 150 IRREPM=1,NIRREP
          IRREPB=DIRPRD(IRREPM,IRROMEGA)
          IRREPA=DIRPRD(IRREPM,IRREPW)
          NROW=DISSYW*VRT(IRREPA,1)
          NCOL=VRT(IRREPB,2)
          NSUM=POP(IRREPM,2)
          IZ=IOFFZ(IRREPB)
          IW=IOFFW(IRREPM)
          IT=IOFFT1(IRREPM,2) 
          CALL XGEMM('N','T',NROW,NCOL,NSUM,ONEM,ICORE(IW),NROW,
     &               ICORE(IT),NCOL,ONE,ICORE(IZ),NROW)
150      CONTINUE
C
C REORDER TARGET BACK TO Z(Ab,Ij)
C
         CALL TRANSP(ICORE(I000),ICORE(I010),DISSYZ,NUMDSZ)
c YAU : old
c        CALL ICOPY(IINTFP*NUMDSZ*DISSYZ,ICORE(I010),1,ICORE(I000),1)
c YAU : new
         CALL DCOPY(NUMDSZ*DISSYZ,ICORE(I010),1,ICORE(I000),1)
c YAU : end
C
         CALL SUMSYM2(ICORE(I000),ICORE(I010),NUMDSZ*DISSYZ,1,
     &                IRREPZR,LISTZ)
C
        ELSE
C
         WRITE(6,*)' out-of-core not coded yet '
         call errex
C
        ENDIF
C
       ENDIF
C
100   CONTINUE
C
      IF(RHF)RETURN
C
C NOW DO AAAA AND BBBB SPIN CASES IF WE HAVEN'T RETURNED YET.
C
      DO 200 ISPIN=1,1+IUHF
       DO 210 IRREPZR=1,NIRREP
C
C LOOP OVER KET IRREPS OF *TARGET*.  THIS IS NOT THE SAME AS THE
C  IRREPS OF THE INTEGRALS AND AMPLITUDES UNLESS IRROMEGA=1.
C
        IRREPZL=DIRPRD(IRREPZR,IRROMEGA)
        IRREPW=IRREPZR
        LISTW=LISTW0-1+ISPIN
        LISTZ=LISTZ0-1+ISPIN
        DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW)) 
        NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW)) 
        DISSYZ=IRPDPD(IRREPZL,ISYTYP(1,LISTZ))
        NUMDSZ=IRPDPD(IRREPZR,ISYTYP(2,LISTZ))
        DISSYZ0=IRPDPD(IRREPZL,18+ISPIN)
        MAXW=MAX(DISSYW,NUMDSW)
        MAXZ=MAX(DISSYZ,NUMDSZ,DISSYZ0)
        I000=1
        I010=I000+IINTFP*MAX(DISSYZ0*NUMDSZ,3*MAXW)
        I020=I010+IINTFP*MAX(DISSYW*NUMDSW,DISSYZ*NUMDSZ)
        IF(I020.LT.MXCOR)THEN
C
C DO IN-CORE ALGORITHM
C
C
C READ W INTO W(i<j,mb) AND TRANSPOSE KET INDICES TO W(i<j,bm)
C
         CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPW,LISTW)
         ITMP1=I000
         ITMP2=ITMP1+IINTFP*MAXW
         ITMP3=ITMP2+IINTFP*MAXW
         ITMP4=ITMP3+IINTFP*MAXW
         CALL SYMTR1 (IRREPW,POP(1,ISPIN),VRT(1,ISPIN),DISSYW,
     &                ICORE(I010),ICORE(ITMP1),ICORE(ITMP2),
     &                ICORE(ITMP3))
C
C COMPUTE OFFSETS FOR W AND Z VECTORS
C
         IW=0
         IZ=0
         DO 220 IRR2=1,NIRREP
          IOFFW(IRR2)=I010+IW
          IOFFZ(IRR2)=I000+IZ
          IRR1W=DIRPRD(IRR2,IRREPW)
          IRR1Z=DIRPRD(IRR2,IRREPZL)
          IW=IW+DISSYW*VRT(IRR1W,ISPIN)*POP(IRR2,ISPIN)*IINTFP
          IZ=IZ+NUMDSZ*VRT(IRR1Z,ISPIN)*VRT(IRR2,ISPIN)*IINTFP
220      CONTINUE
C
C PERFORM MATRIX MULTIPLICATION
C
C         Z(i<j,ba) = SUM W(i<jb,m) * T(ma)
C                      m
C  
         DO 230 IRREPM=1,NIRREP
          IRREPA=DIRPRD(IRREPM,IRROMEGA)
          IRREPB=DIRPRD(IRREPM,IRREPW)
          NROW=DISSYW*VRT(IRREPB,ISPIN)
          NCOL=VRT(IRREPA,ISPIN)
          NSUM=POP(IRREPM,ISPIN)
          IZ=IOFFZ(IRREPA)
          IW=IOFFW(IRREPM)
          IT=IOFFT1(IRREPM,ISPIN) 
          CALL XGEMM('N','T',NROW,NCOL,NSUM,ONEM,ICORE(IW),NROW,
     &               ICORE(IT),NCOL,ZILCH,ICORE(IZ),NROW)
          call vminus(icore(iz),nrow*ncol)
230      CONTINUE
C
C NOW ANTISYMMETRIZE TARGET TO Z(i<j,a<b) AND DUMP IT TO INCREMENTS
C
         CALL ASSYM2(IRREPZL,VRT(1,ISPIN),NUMDSZ,ICORE(I000))
         CALL TRANSP(ICORE(I000),ICORE(I010),DISSYZ,NUMDSZ)
c YAU : old
c        CALL ICOPY(IINTFP*NUMDSZ*DISSYZ,ICORE(I010),1,ICORE(I000),1)
c YAU : new
         CALL DCOPY(NUMDSZ*DISSYZ,ICORE(I010),1,ICORE(I000),1)
c YAU : end
         CALL SUMSYM2(ICORE(I000),ICORE(I010),NUMDSZ*DISSYZ,1,
     &                IRREPZR,LISTZ)
C
        ELSE
C
C OUT-OF-CORE ALGORITHM
C
         WRITE(6,*)' out-of-core not coded yet '
         call errex
C
        ENDIF
C
210    CONTINUE
200   CONTINUE
C
      RETURN
      END
