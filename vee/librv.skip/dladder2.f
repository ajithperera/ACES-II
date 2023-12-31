      SUBROUTINE DLADDER2(ICORE,MAXCOR,IUHF,IRREPX,ITYPE,
     &                   LISTT0,LISTW0,LISTZ0,ISIDE)
C
C THIS SUBROUTINE CALCULATES THE PART-PART LADDER CONTRIBUTION
C USING COMPRESSED A<=B,CD INTEGRAL LISTS.  IT ALSO EVALUATES 
C AN ADDITIONAL CONTRIBUTION ARISING FORMALLY FROM THE ABCI 
C EFFECTIVE HAMILTONIAN MATRIX ELEMENTS IN A RATHER INDIRECT
C AND SNEAKY WAY (LEFT-HAND EIGENVECTOR ONLY).
C
C
C     Z(ab,ij) = SUM W(ab,ef) * T(ef,ij)  [ITYPE=6]
C                e,f
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ZILCH,FACT
C
      DIMENSION ICORE(MAXCOR)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      DATA ONE  /1.0D0/
      DATA ZILCH/0.0D0/
C
      ISPIN=3
C
C ALPHA-BETA SPIN CASE ONLY (THIS RUNS ONLY FOR RHF!)
C
      DO 100 IRREPZR=1,NIRREP
C
C LOOP OVER KET IRREPS OF *TARGET*.  THIS IS NOT THE SAME AS THE
C  IRREPS OF THE INTEGRALS AND AMPLITUDES UNLESS IRREPX=1.
C
       IRREPZL=DIRPRD(IRREPZR,IRREPX)
       LISTW=LISTW0+ISPIN-1
       LISTT=LISTT0+ISPIN-1
       LISTZ=LISTZ0+ISPIN-1
       IRREPW=IRREPZL
       DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW))
       NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW)) 
       DISSYZ=IRPDPD(IRREPZL,ISYTYP(1,LISTZ))
       NUMDSZ=IRPDPD(IRREPZR,ISYTYP(2,LISTZ))
       DISSYT=DISSYZ
       NUMDST=NUMDSZ
       I000=1
       I010=I000+IINTFP*DISSYZ*NUMDSZ
       I020=I010+IINTFP*DISSYT*NUMDST
C
C USE GENERAL ALGORITHM ALLOWING BOTH IN-CORE AND OUT-OF-CORE
C SOLUTIONS
C
       FACT=ZILCH
       CORLFT=MAXCOR-I020+1
       IF(DISSYW.NE.0)THEN
        NINCOR=CORLFT/(DISSYW*IINTFP)
        IF(NINCOR.EQ.0)THEN
         WRITE(6,1000)
1000     FORMAT(T3,'@DLADDER-F, Not enough memory for ladders.')
         CALL INSMEM('DLADDER',DISSYW*IINTFP,CORLFT)
        ENDIF
       ELSE
        NINCOR=DISSYW
       ENDIF
       IOFFT=I010
       IOFFZ=I000
       NLEFT=NUMDSW
       NFIRST=1
       CALL GETLST(ICORE(I010),1,NUMDST,1,IRREPZR,LISTT)
1      NREAD =MIN(NLEFT,NINCOR)
       CALL GETLST(ICORE(I020),NFIRST,NREAD,1,IRREPW,LISTW)
C
       CALL XGEMM ('N','N',DISSYW,NUMDSZ,NREAD,ONE,ICORE(I020),
     &             DISSYW,ICORE(IOFFT),DISSYT,FACT,ICORE(IOFFZ),
     &             DISSYW)
       IOFFT=IOFFT+IINTFP*NREAD
       FACT =ONE
       NFIRST=NFIRST+NREAD
       NLEFT =NLEFT -NREAD
       IF(NLEFT.NE.0)GOTO 1
C
       CALL SYMEXP5(IRREPZL,IRREPZR,VRT(1,1),POP(1,1),DISSYT,DISSYW,
     &              NUMDST,ICORE(I000),ICORE(I000),ICORE(I010))
C
       IF(ISIDE.EQ.2)THEN
        CALL SCOPY(DISSYZ*NUMDSZ,ICORE(I000),1,ICORE(I010),1)
        ITMP1=I020
        ITMP2=ITMP1+IINTFP*MAX(DISSYZ,NUMDSZ,IRPDPD(IRREPX,9),NT(1))
        ITMP3=ITMP2+IINTFP*MAX(DISSYZ,NUMDSZ,IRPDPD(IRREPX,9),NT(1))
        CALL SPINAD1(IRREPZR,POP(1,1),DISSYZ,ICORE(I010),ICORE(ITMP1),
     &               ICORE(ITMP2))
        CALL GETLST(ICORE(ITMP1),1,1,1,3,490)
        CALL GETLST(ICORE(ITMP2),1,1,1,1,90)
        CALL DDDOT24(IRREPX,1,IRREPZL,IRREPZR,ICORE(ITMP1),ICORE(ITMP2),
     &               ICORE(I010),ICORE(ITMP3),DISSYZ,VRT(1,1),POP(1,1),
     &               VRT(1,1),VRT(1,1),POP(1,1),POP(1,1),'STST')
        CALL PUTLST(ICORE(ITMP1),1,1,1,3,490)
       ENDIF
C
       CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPZR,LISTZ)
C         
       CALL SAXPY (NUMDSZ*DISSYZ,ONE,ICORE(I010),1,ICORE(I000),1)
       CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPZR,LISTZ)
100    CONTINUE
50    CONTINUE
C
      RETURN
      END
