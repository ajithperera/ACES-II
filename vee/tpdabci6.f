      SUBROUTINE TPDABCI6(ICORE,MAXCOR,IUHF)
C
C CALCULATION OF THE SIXTH ABCI CONTRIBUTION TO
C THE EOM-CCSD TWO-PARTICLE DENSITY MATRIX
C
C  G(AB,CI) = 1/2 P(AB) H(KB,CI) T(KA)
C
CEND 
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ZILCH,HALF,ONE,TWO
      DIMENSION ICORE(MAXCOR)
      DIMENSION I0T(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMLOC/ISYMOFF(8,8,25)
C
      DATA ZILCH,HALF,ONE /0.D0,0.5D0,1.0D0/

C READ IN T1 AMPLITUDES
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
       DO 10 ISPIN=1,1+IUHF
        LISTH=53+ISPIN
        LISTG=126+ISPIN
        DO 100 IRREP=1,NIRREP
         DISSYH=IRPDPD(IRREP,ISYTYP(1,LISTH))
         NUMDSH=IRPDPD(IRREP,ISYTYP(2,LISTH))
         DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
         NUMDSG=IRPDPD(IRREP,ISYTYP(2,LISTG))
         NVRTSQ=IRPDPD(IRREP,18+ISPIN)
         I000=ISTART
         I010=I000+IINTFP*NUMDSG*NVRTSQ
         IEND=I010+IINTFP*MAX(NUMDSH*DISSYH,NUMDSG*DISSYG)
         IF(IEND.LE.MAXCOR)THEN
C
C IN-CORE ALGORITHM
C
          CALL GETLST(ICORE(I010),1,NUMDSH,1,IRREP,LISTH)
C
C READ H(CI,BK) FROM DISK
C
C
C FORM PRODUCT  G(CI,BA) = H(CI,BK)*T(AK)
C
          DO 110 IRREPK=1,NIRREP
           IRREPA=IRREPK
           IRREPB=DIRPRD(IRREPK,IRREP)
           NUMK=POP(IRREPK,ISPIN)
           NUMB=VRT(IRREPB,ISPIN)
           NUMA=VRT(IRREPA,ISPIN)
           NROW=NUMDSG*NUMB
           NCOL=NUMA
           NSUM=NUMK
           IOFFH=I010+(ISYMOFF(IRREPK,IRREP,8+ISPIN)-1)*DISSYH*IINTFP
           IOFFG=I000+(ISYMOFF(IRREPA,IRREP,18+ISPIN)-1)*NUMDSG*IINTFP
           IOFFT=I0T(ISPIN)+(ISYMOFF(IRREPK,1,8+ISPIN)-1)*IINTFP
           CALL XGEMM('N','T',NROW,NCOL,NSUM,HALF,ICORE(IOFFH),
     &                NROW,ICORE(IOFFT),NCOL,ZILCH,ICORE(IOFFG),
     &                NROW)
110       CONTINUE
          CALL ASSYM2(IRREP,VRT(1,ISPIN),NUMDSG,ICORE(I000))
          CALL TRANSP(ICORE(I000),ICORE(I010),DISSYG,NUMDSG)
          CALL GETLST(ICORE(I000),1,NUMDSG,1,IRREP,LISTG)
          CALL SAXPY (NUMDSG*DISSYG,ONE,ICORE(I010),1,ICORE(I000),1)
          CALL PUTLST(ICORE(I000),1,NUMDSG,1,IRREP,LISTG)
C
         ELSE
C
C OUT-OF-CORE ALGORITHM
C
C   G(CI,BA) = H(CI,BK)*R(AK) , ONE CI AT A TIME
C
          I010=I000+IINTFP*DISSYH*NUMDSH
          I020=I010+IINTFP*MAX(NUMDSH,DISSYH)
          I030=I020+IINTFP*NVRTSQ
          CALL GETLST(ICORE(I000),1,NUMDSH,1,IRREP,LISTH)
          I0H=I000
          DO 111 IPASS=1,NUMDSG
           CALL SCOPY(NUMDSH,ICORE(I0H),DISSYH,ICORE(I010),1)
           I0H=I0H+IINTFP
           DO 112 IRREPK=1,NIRREP
            IRREPA=IRREPK
            IRREPB=DIRPRD(IRREPK,IRREP)
            NUMK=POP(IRREPK,ISPIN)
            NUMB=VRT(IRREPB,ISPIN)
            NUMA=VRT(IRREPA,ISPIN)
            NROW=NUMB
            NCOL=NUMA
            NSUM=NUMK
            IOFFH=I010+(ISYMOFF(IRREPK,IRREP,8+ISPIN)-1)*IINTFP
            IOFFT=I0T(ISPIN)+(ISYMOFF(IRREPK,1,8+ISPIN)-1)*IINTFP
            IOFFG=I020+(ISYMOFF(IRREPA,IRREP,18+ISPIN)-1)*IINTFP
            CALL XGEMM('N','T',NROW,NCOL,NSUM,HALF,ICORE(IOFFH),
     &                 NROW,ICORE(IOFFT),NCOL,ZILCH,ICORE(IOFFG),
     &                 NROW)
112        CONTINUE
           CALL ASSYM2(IRREP,VRT(1,ISPIN),1,ICORE(I020))
           CALL GETLST(ICORE(I030),IPASS,1,1,IRREP,LISTG)
           CALL SAXPY (DISSYG,ONE,ICORE(I030),1,ICORE(I020),1)
           CALL PUTLST(ICORE(I020),IPASS,1,1,IRREP,LISTG)
111       CONTINUE
         ENDIF
100     CONTINUE
10     CONTINUE
      ENDIF
C
C ABAB AND BABA SPIN CASES
C
      DO 20 ISPIN=1,1+IUHF
       LISTG=131-ISPIN
       LISTH1=55+ISPIN
       LISTH2=57+ISPIN
       DO 200 IRREP=1,NIRREP
        DISSYH1=IRPDPD(IRREP,ISYTYP(1,LISTH1))
        NUMDSH1=IRPDPD(IRREP,ISYTYP(2,LISTH1))
        DISSYH2=IRPDPD(IRREP,ISYTYP(1,LISTH2))
        NUMDSH2=IRPDPD(IRREP,ISYTYP(2,LISTH2))
        DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
        NUMDSG=IRPDPD(IRREP,ISYTYP(2,LISTG))
        MAXT=MAX(NUMDSG,DISSYG,NUMDSH2,DISSYH2,NUMDSH1,DISSYH1)
        I000=ISTART
        I010=I000+IINTFP*NUMDSG*DISSYG
        ITMP1=I010+IINTFP*MAX(DISSYG*NUMDSG,DISSYH1*NUMDSH1,
     &                       DISSYH2*NUMDSH2)
        ITMP2=ITMP1+IINTFP*MAXT
        ITMP3=ITMP2+IINTFP*MAXT
        IEND=ITMP3 +IINTFP*MAXT
C
        IF(IEND.LE.MAXCOR)THEN
C
C IN-CORE ALGORITHM
C
         CALL ZERO(ICORE(I000),DISSYG*NUMDSG)
         CALL GETLST(ICORE(I010),1,NUMDSH1,1,IRREP,LISTH1)
C
C FORM PRODUCT:
C                 G(Ci,bA) <= H(Ci,bK) * T(KA) [ISPIN=1]
C                 G(cI,Ab) <= H(cI,Ak) * T(kb) [ISPIN=2]
C
         DO 210 IRREPK=1,NIRREP
          IRREPA=IRREPK
          IRREPB=DIRPRD(IRREPK,IRREP)
          NUMK=POP(IRREPK,ISPIN)
          NUMA=VRT(IRREPA,ISPIN)
          NUMB=VRT(IRREPB,3-ISPIN)
          NROW=NUMDSG*NUMB
          NCOL=NUMA
          NSUM=NUMK
          ITH=13-ISPIN
          ITG=23-10*(ISPIN-1)
          IOFFH=I010+(ISYMOFF(IRREPK,IRREP,ITH)-1)*DISSYH1*IINTFP
          IOFFG=I000+(ISYMOFF(IRREPA,IRREP,ITG)-1)*NUMDSG*IINTFP
          IOFFT=I0T(ISPIN)+(ISYMOFF(IRREPK,1,8+ISPIN)-1)*IINTFP
          CALL XGEMM('N','T',NROW,NCOL,NSUM,HALF,ICORE(IOFFH),
     &               NROW,ICORE(IOFFT),NCOL,ZILCH,ICORE(IOFFG),
     &               NROW)
210      CONTINUE
C
CSSS         call checksum('tp6 in g',icore(i000),numdsg*dissyg,s)
CSSS         call checksum('tp6 inh1',icore(i010),numdsh1*dissyh1,s)
C
C TRANSPOSE LAST TWO INDICES OF G FOR SECOND CONTRACTION
C
         CALL SYMTR1(IRREP,VRT(1,3-ISPIN),VRT(1,ISPIN),NUMDSG,
     &               ICORE(I000),ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
C
         CALL GETLST(ICORE(I010),1,NUMDSH2,1,IRREP,LISTH2)
C
C FORM SECOND PRODUCT:
C
C                 G(Ci,Ab) <= H(Ci,Ak)*T(kb)  [ISPIN=1]
C                 G(cI,bA) <= H(cI,bK)*T(KA)  [ISPIN=2]
C
         DO 220 IRREPK=1,NIRREP
          IRREPB=IRREPK
          IRREPA=DIRPRD(IRREPK,IRREP)
          NUMK=POP(IRREPK,3-ISPIN)
          NUMA=VRT(IRREPA,ISPIN)
          NUMB=VRT(IRREPB,3-ISPIN)
          NROW=NUMDSG*NUMA
          NCOL=NUMB
          NSUM=NUMK
          ITH=10+ISPIN
          ITG=13+(ISPIN-1)*10
          IOFFH=I010+(ISYMOFF(IRREPK,IRREP,ITH)-1)*DISSYH2*IINTFP
          IOFFG=I000+(ISYMOFF(IRREPB,IRREP,ITG)-1)*NUMDSG*IINTFP
          IOFFT=I0T(3-ISPIN)+(ISYMOFF(IRREPK,1,11-ISPIN)-1)*IINTFP
          CALL XGEMM('N','T',NROW,NCOL,NSUM,HALF,ICORE(IOFFH),
     &               NROW,ICORE(IOFFT),NCOL,ONE,ICORE(IOFFG),
     &               NROW)
220      CONTINUE
         IF(ISPIN.EQ.2) THEN
          CALL SYMTR1(IRREP,VRT(1,2),VRT(1,1),NUMDSG,ICORE(I000),
     &                ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
          CALL SYMTR3(IRREP,VRT(1,2),POP(1,1),NUMDSG,DISSYG,ICORE(I000),
     &                ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
         ENDIF
         CALL TRANSP(ICORE(I000),ICORE(I010),DISSYG,NUMDSG)
C
C SAVE G ON DISK
C
CSSS         call checksum('tpdabci6',icore(i010),numdsg*dissyg,s)
         CALL GETLST(ICORE(I000),1,NUMDSG,1,IRREP,LISTG)
         CALL SAXPY (NUMDSG*DISSYG,ONE,ICORE(I010),1,ICORE(I000),1)
         CALL PUTLST(ICORE(I000),1,NUMDSG,1,IRREP,LISTG)
        ELSE
C
C OUT-OF-CORE ALGORITHM
C
         I010=I000+IINTFP*DISSYH1*NUMDSH1
         I020=I010+IINTFP*DISSYH2*NUMDSH2
         I030=I020+IINTFP*MAXT
         I040=I030+IINTFP*MAXT
         I050=I040+IINTFP*MAXT
         I060=I050+IINTFP*MAXT
         CALL GETLST(ICORE(I000),1,NUMDSH1,1,IRREP,LISTH1)
         CALL GETLST(ICORE(I010),1,NUMDSH2,1,IRREP,LISTH2)
         IF(ISPIN.EQ.2)THEN
          CALL SYMTR3(IRREP,VRT(1,2),POP(1,1),DISSYH1,NUMDSH1,
     &                ICORE(I000),ICORE(I020),ICORE(I030),ICORE(I040))
          CALL SYMTR3(IRREP,VRT(1,2),POP(1,1),DISSYH2,NUMDSH2,
     &                ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040))
         ENDIF
C
         I0H1=I000
         I0H2=I010
         DO 230 IPASS=1,NUMDSG
          CALL SCOPY(NUMDSH1,ICORE(I0H1),DISSYH1,ICORE(I020),1)
          CALL SCOPY(NUMDSH2,ICORE(I0H2),DISSYH2,ICORE(I030),1)
          I0H1=I0H1+IINTFP
          I0H2=I0H2+IINTFP
C
C FORM PRODUCT:
C                 G(Ci,bA) <= H(Ci,bK) * R(KA) [ISPIN=1]
C                 G(Ic,Ab) <= H(Ic,Ak) * R(kb) [ISPIN=2]
C
          DO 240 IRREPK=1,NIRREP
           IRREPA=IRREPK
           IRREPB=DIRPRD(IRREPK,IRREP)
           NUMK=POP(IRREPK,ISPIN)
           NUMA=VRT(IRREPA,ISPIN)
           NUMB=VRT(IRREPB,3-ISPIN)
           NROW=NUMB
           NCOL=NUMA
           NSUM=NUMK
           ITH=13-ISPIN
           ITG=23-10*(ISPIN-1)
           IOFFH=I020+(ISYMOFF(IRREPK,IRREP,ITH)-1)*IINTFP
           IOFFG=I040+(ISYMOFF(IRREPA,IRREP,ITG)-1)*IINTFP
           IOFFT=I0T(ISPIN)+(ISYMOFF(IRREPK,1,8+ISPIN)-1)*IINTFP
           CALL XGEMM('N','T',NROW,NCOL,NSUM,HALF,ICORE(IOFFH),
     &                NROW,ICORE(IOFFT),NCOL,ZILCH,ICORE(IOFFG),
     &                NROW)
240       CONTINUE
          CALL SYMTRA(IRREP,VRT(1,3-ISPIN),VRT(1,ISPIN),1,
     &                ICORE(I040),ICORE(I050))
C
C FORM SECOND PRODUCT:
C
C                 G(Ci,Ab) <= H(Ci,Ak)*R(kb)  [ISPIN=1]
C                 G(cI,bA) <= H(cI,bK)*R(KA)  [ISPIN=2]
C
          DO 250 IRREPK=1,NIRREP
           IRREPB=DIRPRD(IRREPK,1)
           IRREPA=DIRPRD(IRREPK,IRREP)
           NUMK=POP(IRREPK,3-ISPIN)
           NUMA=VRT(IRREPA,ISPIN)
           NUMB=VRT(IRREPB,3-ISPIN)
           NROW=NUMA
           NCOL=NUMB
           NSUM=NUMK
           ITH=10+ISPIN
           ITG=13+(ISPIN-1)*10
           IOFFH=I030+(ISYMOFF(IRREPK,IRREP,ITH)-1)*IINTFP
           IOFFG=I050+(ISYMOFF(IRREPB,IRREP,ITG)-1)*IINTFP
           IOFFT=I0T(3-ISPIN)+(ISYMOFF(IRREPK,1,11-ISPIN)-1)*IINTFP
           CALL XGEMM('N','T',NROW,NCOL,NSUM,HALF,ICORE(IOFFH),
     &                NROW,ICORE(IOFFT),NCOL,ONE,ICORE(IOFFG),
     &                NROW)
250       CONTINUE
          IF(ISPIN.EQ.2) THEN
           CALL SYMTRA(IRREP,VRT(1,2),VRT(1,1),1,ICORE(I050),
     &                 ICORE(I040))
           CALL SCOPY (DISSYG,ICORE(I040),1,ICORE(I050),1)
          ENDIF
C
C SAVE G ON DISK
C
          CALL GETLST(ICORE(I040),IPASS,1,1,IRREP,LISTG)
          CALL SAXPY (DISSYG,ONE,ICORE(I050),1,ICORE(I040),1)
          CALL PUTLST(ICORE(I040),IPASS,1,1,IRREP,LISTG)
230      CONTINUE
C
        ENDIF
C
200    CONTINUE
20    CONTINUE
C
C ALL DONE, RETURN
C
      TWO=0.D0
      if(iuhf.eq.0) then
       call checkgam1(icore,30,130,two,iuhf,2,vrt,s)
      endif
      IF(IUHF.EQ.1) THEN
       CALL CHECKGAM(ICORE,30,130,TWO)
       CALL CHECKGAM(ICORE,27,127,TWO)
       CALL CHECKGAM(ICORE,28,128,TWO)
       CALL CHECKGAM(ICORE,29,129,TWO)
      ENDIF
C
      RETURN
      END
