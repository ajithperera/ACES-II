      SUBROUTINE DT1INT2B(ICORE,MAXCOR,IUHF,IRREPX,LISTW0,LISTZ0,
     &                     LISTT1)
C
C THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF T1 TO
C   T2 (T2<-T1).  THIS CODE IS GENERAL IN THE SENSE THAT
C   IT DOES NOT ASSUME THAT T (AND Z) ARE TOTALLY SYMMETRIC
C   QUANTITIES, BUT RATHER TRANSFORM AS IRREPX.  THIS ALSO
C   INCLUDES OUT-OF-CORE SOLUTIONS.
C
C CONTRACTION 1:
C
C     Z(ab,ij) = SUM  T(e,i) * <ab||ej>   [AAAA and BBBB]
C                 e
C
C     Z(Ab,Ij) = SUM  T(E,I) * <Ab||Ej> + T(e,j) * <Ab||Ie> [ABAB]
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
      COMMON /SYMLOC/ ISYMOFF(8,8,25)
C
      DATA ONE  /1.0D0/
      DATA ONEM /-1.0D0/
      DATA ZILCH/0.0D0/
C
      CALL GETDT1(ICORE,MAXCOR,MXCOR,IUHF,IOFFT1,IRREPX,LISTT1)
      RHF=IUHF.EQ.0
C
C
C DO ALPHA-BETA SPIN CASE.  CONTAINS BOTH GENERAL AB AND 
C SPIN-ADAPTED RHF CODE.  THIS BLOCK OF CODE ALWAYS RUNS.
C
      DO 100 IRREPZR=1,NIRREP
C
C LOOP OVER KET IRREPS OF *TARGET*.  THIS IS NOT THE SAME AS THE
C  IRREPS OF THE INTEGRALS AND AMPLITUDES UNLESS IRREPX=1.
C
       IRREPZL=DIRPRD(IRREPZR,IRREPX)
       IRREPW=IRREPZL
       LISTW=LISTW0+3
       LISTZ=LISTZ0+2
       DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW))
       NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW)) 
       DISSYZ=IRPDPD(IRREPZL,ISYTYP(1,LISTZ))
       NUMDSZ=IRPDPD(IRREPZR,ISYTYP(2,LISTZ))
       MAXBUF=MAX(NUMDSW,DISSYW,DISSYZ,NUMDSZ)
       I000=1
       I010=I000+IINTFP*MAX(DISSYZ*NUMDSZ,3*MAXBUF)
       I020=I010+IINTFP*MAX(DISSYW*NUMDSW,3*MAXBUF,NUMDSZ*DISSYZ)
       IF(I020.LT.MXCOR)THEN
C
C DO IN-CORE ALGORITHM
C
C
C READ W INTO W(Ab,Ej) AND TRANSPOSE KET INDICES TO W(Ab,jE).  
C
        ITMP1=I000
        ITMP2=ITMP1+IINTFP*MAXBUF
        ITMP3=ITMP2+IINTFP*MAXBUF
        ITMP4=ITMP3+IINTFP*MAXBUF
        CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPW,LISTW)
        CALL SYMTR1 (IRREPW,VRT(1,1),POP(1,2),DISSYW,ICORE(I010),
     &               ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
C
C PERFORM MATRIX MULTIPLICATION
C
C         Z(Ab,jI) = SUM W(Abj,E) * T(E,I)
C                     e
C  
        DO 120 IRREPE=1,NIRREP
         IRREPI=DIRPRD(IRREPE,IRREPX)
         IRREPJ=DIRPRD(IRREPE,IRREPW)
         NROW=DISSYW*POP(IRREPJ,2)
         NCOL=POP(IRREPI,1)
         NSUM=VRT(IRREPE,1)
         IZ=I000+IINTFP*DISSYZ*(ISYMOFF(IRREPI,IRREPZR,24)-1)
         IW=I010+IINTFP*DISSYW*(ISYMOFF(IRREPE,IRREPW,25) -1)
         IT=IOFFT1(IRREPI,1) 
         IF (NSUM .NE.0) THEN 
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IW),NROW,
     &              ICORE(IT),NSUM,ZILCH,ICORE(IZ),NROW)
         ELSE
         CALL ZERO(ICORE(IZ),NROW*NCOL)
         ENDIF 
120     CONTINUE
C
C TRANSPOSE KET INDICES BACK TO Z(Ab,Ij)
C
        ITMP1=I010
        ITMP2=ITMP1+IINTFP*MAXBUF
        ITMP3=ITMP2+IINTFP*MAXBUF
        ITMP4=ITMP3+IINTFP*MAXBUF
        CALL SYMTR1 (IRREPZR,POP(1,2),POP(1,1),DISSYZ,ICORE(I000),
     &               ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
       ELSE
C
C OUT-OF-CORE ALGORITHM
C
        IZ=I000
        IDISW1=1
        DO 130 IRREPJ=1,NIRREP
         IRREPE=DIRPRD(IRREPJ,IRREPW)
         IRREPI=DIRPRD(IRREPJ,IRREPZR)
         NUMJ=POP(IRREPJ,2)
         NUMI=POP(IRREPI,1)
         NUME=VRT(IRREPE,1)
         ITOP=I010+IINTFP*MAX(DISSYZ*NUMDSZ,DISSYW*NUME)
         IF(ITOP.GT.MXCOR)THEN
          CALL INSMEM('DT1INT2B2',ITOP,MXCOR)
         ENDIF
C
C LOOP OVER J, READ IN ALL W(Ab,E) AND FORM
C
C   Z(Ab,I) = W(Ab,E) * T(EI) FOR THIS j
C
         DO 131 J=1,NUMJ
          CALL GETLST(ICORE(I010),IDISW1,NUME,1,IRREPW,LISTW)
          NROW=DISSYW
          NCOL=NUMI
          NSUM=NUME
          IT=IOFFT1(IRREPI,1) 
          IF (NSUM .NE. 0) THEN
          CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(I010),NROW,
     &               ICORE(IT),NSUM,ZILCH,ICORE(IZ),NROW)
          ELSE
          CALL ZERO(ICORE(IZ),NROW*NCOL)
          ENDIF
          IZ=IZ+IINTFP*DISSYZ*NUMI
          IDISW1=IDISW1+NUME
131      CONTINUE
130     CONTINUE
C
       ENDIF
C
       IF(RHF)THEN
C
C IF RHF, SYMMETRIZE TARGET, DUMP TO DISK AND RETURN
C
C           Z(ab,ij) = Z(ab,ij) + Z(ba,ji)
C     
        ITMP1=I010
        ITMP2=ITMP1+IINTFP*MAXBUF
        ITMP3=ITMP2+IINTFP*MAXBUF
        ITMP4=ITMP3+IINTFP*MAXBUF
        CALL SYMRHF3(IRREPZL,IRREPZR,VRT(1,1),POP(1,1),DISSYZ,
     &               ICORE(I000),ICORE(ITMP1),ICORE(ITMP2),
     &               ICORE(ITMP3))
        CALL GETLST (ICORE(I010),1,NUMDSZ,1,IRREPZR,LISTZ)
        CALL SAXPY  (NUMDSZ*DISSYZ,ONE,ICORE(I010),1,ICORE(I000),1)
        CALL PUTLST (ICORE(I000),1,NUMDSZ,1,IRREPZR,LISTZ)
C
       ELSE
C
C IF UHF, HAVE TO DO ANOTHER CONTRACTION FOR AB CASE
C
C        Z(Ab,Ij) = W(Ab,Je) * T(e,i)
C
C
C READ W INTO W(Ab,Je)
C
        LISTW=LISTW0+2
        DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW)) 
        NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW)) 
        I020=I010+IINTFP*MAX(DISSYW*NUMDSW,3*MAXBUF,NUMDSZ*DISSYZ)
        IF(I020.LT.MXCOR)THEN
C
C IN-CORE ALGORITHM
C
         CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPW,LISTW)
C
C PERFORM MATRIX MULTIPLICATION
C
C         Z(Ab,Ji) = SUM W(Ab,Je) * T(e,i)
C                     e
C  
         DO 150 IRREPE=1,NIRREP
          IRREPI=DIRPRD(IRREPE,IRREPX)
          IRREPJ=DIRPRD(IRREPE,IRREPW)
          NROW=DISSYW*POP(IRREPJ,1)
          NCOL=POP(IRREPI,2)
          NSUM=VRT(IRREPE,2)
          IZ=I000+IINTFP*DISSYZ*(ISYMOFF(IRREPI,IRREPZR,14)-1)
          IW=I010+IINTFP*DISSYW*(ISYMOFF(IRREPE,IRREPW ,18)-1)
          IT=IOFFT1(IRREPI,2) 
          IF (NSUM .NE. 0) THEN
          CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IW),NROW,
     &               ICORE(IT),NSUM,ONE,ICORE(IZ),NROW)
          ENDIF 
150      CONTINUE
C
        ELSE
C
C OUT-OF-CORE ALGORITHM
C
         IDISW1=1
         DO 151 IRREPE=1,NIRREP
          IRREPI=DIRPRD(IRREPE,IRREPX)
          IRREPJ=DIRPRD(IRREPI,IRREPZR)
          IZ=I000+IINTFP*DISSYZ*(ISYMOFF(IRREPI,IRREPZR,14)-1)
          NUMJ=POP(IRREPJ,1)
          NUMI=POP(IRREPI,2)
          NUME=VRT(IRREPE,2)
          ITOP=I010+IINTFP*MAX(DISSYZ*NUMDSZ,DISSYW*NUMJ)
          IF(ITOP.GT.MXCOR)THEN
           CALL INSMEM('DT1INT2B2',ITOP,MXCOR)
          ENDIF
C
C LOOP OVER e, READ IN ALL W(Ab,J) AND FORM
C
C   Z(Ab,Ji) = W(Ab,J) * T(e,i) FOR THIS e
C
          DO 152 E=1,NUME
           CALL GETLST(ICORE(I010),IDISW1,NUMJ,1,IRREPW,LISTW)
           NROW=DISSYZ*NUMJ
           NCOL=NUMI
           NSUM=1
           IT=IOFFT1(IRREPI,2)+IINTFP*(E-1) 
           CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(I010),NROW,
     &                ICORE(IT),NUME,ONE,ICORE(IZ),NROW)
           IDISW1=IDISW1+NUMJ
152       CONTINUE
151      CONTINUE
C
        ENDIF
C
        CALL SUMSYM2(ICORE(I000),ICORE(I010),NUMDSZ*DISSYZ,1,
     &                IRREPZR,LISTZ)
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
C  IRREPS OF THE INTEGRALS AND AMPLITUDES UNLESS IRREPX=1.
C
        IRREPZL=DIRPRD(IRREPZR,IRREPX)
        IRREPW=IRREPZL
        LISTW=LISTW0-1+ISPIN
        LISTZ=LISTZ0-1+ISPIN
        DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW)) 
        NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW)) 
        DISSYZ=IRPDPD(IRREPZL,ISYTYP(1,LISTZ))
        NUMDSZ=IRPDPD(IRREPZR,ISYTYP(2,LISTZ))
        NUMDSZ0=IRPDPD(IRREPZR,20+ISPIN)
        MAXBUF=MAX(DISSYZ,NUMDSZ,NUMDSZ0,DISSYW,NUMDSW)
        I000=1
        I010=I000+IINTFP*MAX(DISSYZ*NUMDSZ0,3*MAXBUF)
        I020=I010+IINTFP*MAX(DISSYW*NUMDSW,DISSYZ*NUMDSZ)
        IF(I020.LT.MXCOR)THEN
C
C DO IN-CORE ALGORITHM
C
C
C READ W INTO W(a<b,ej) AND TRANSPOSE KET INDICES TO W(a<b,je)
C
         CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPW,LISTW)
         ITMP1=I000
         ITMP2=ITMP1+IINTFP*MAXBUF
         ITMP3=ITMP2+IINTFP*MAXBUF
         ITMP4=ITMP3+IINTFP*MAXBUF
         CALL SYMTR1 (IRREPW,VRT(1,ISPIN),POP(1,ISPIN),DISSYW,
     &                ICORE(I010),ICORE(ITMP1),ICORE(ITMP2),
     &                ICORE(ITMP3))
C
C PERFORM MATRIX MULTIPLICATION
C
C         Z(a<b,ji) = SUM W(a<bj,e) * T(e,i)
C                     e
C  
         DO 230 IRREPE=1,NIRREP
          IRREPI=DIRPRD(IRREPE,IRREPX)
          IRREPJ=DIRPRD(IRREPE,IRREPW)
          NROW=DISSYW*POP(IRREPJ,ISPIN)
          NCOL=POP(IRREPI,ISPIN)
          NSUM=VRT(IRREPE,ISPIN)
          IZ=I000+IINTFP*DISSYZ*(ISYMOFF(IRREPI,IRREPZR,20+ISPIN)-1)
          IW=I010+IINTFP*DISSYW*(ISYMOFF(IRREPE,IRREPW,15+ISPIN)-1)
          IT=IOFFT1(IRREPI,ISPIN) 
          IF (NSUM .NE. 0) THEN
          CALL XGEMM('N','N',NROW,NCOL,NSUM,ONEM,ICORE(IW),NROW,
     &               ICORE(IT),NSUM,ZILCH,ICORE(IZ),NROW)
          ELSE
          CALL ZERO(ICORE(IZ),NROW*NCOL)
          ENDIF 
230      CONTINUE
C
        ELSE
C
C OUT-OF-CORE ALGORITHM
C
         IZ=I000
         IDISW1=1
         DO 231 IRREPJ=1,NIRREP
          IRREPE=DIRPRD(IRREPJ,IRREPW)
          IRREPI=DIRPRD(IRREPJ,IRREPZR)
          NUMJ=POP(IRREPJ,ISPIN)
          NUMI=POP(IRREPI,ISPIN)
          NUME=VRT(IRREPE,ISPIN)
          ITOP=I010+IINTFP*MAX(DISSYZ*NUMDSZ,DISSYW*NUME)
          IF(ITOP.GT.MXCOR)THEN
           CALL INSMEM('DT1INT2B2',ITOP,MXCOR)
          ENDIF
C
C LOOP OVER j, READ IN ALL W(a<b,e) AND FORM
C
C   Z(Ab,I) = W(a<be) * T(ei) FOR THIS j
C
          DO 232 J=1,NUMJ
           CALL GETLST(ICORE(I010),IDISW1,NUME,1,IRREPW,LISTW)
           NROW=DISSYW
           NCOL=NUMI
           NSUM=NUME
           IT=IOFFT1(IRREPI,ISPIN) 
           IF (NSUM .NE. 0) THEN
           CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(I010),NROW,
     &                ICORE(IT),NSUM,ZILCH,ICORE(IZ),NROW)
           ELSE
           CALL ZERO(ICORE(IZ),NROW*NCOL)
           ENDIF 
           IZ=IZ+IINTFP*DISSYZ*NUMI
           IDISW1=IDISW1+NUME
232       CONTINUE
231      CONTINUE
C
        ENDIF
C
C
C NOW ANTISYMMETRIZE TARGET TO Z(a<b,i<j) AND DUMP IT TO INCREMENTS
C
        CALL ASSYM2 (IRREPZR,POP(1,ISPIN),DISSYZ,ICORE(I000))
        CALL SUMSYM2(ICORE(I000),ICORE(I010),NUMDSZ*DISSYZ,1,
     &               IRREPZR,LISTZ)
210    CONTINUE
200   CONTINUE
C
      RETURN
      END
