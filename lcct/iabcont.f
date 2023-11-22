      SUBROUTINE IABCONT(ICORE,MAXCOR,IUHF,LISTL2,LISTR1,LISTR1OFF,
     &   LISTR2)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,HALF,HALFM,FACT
      LOGICAL TAU2
C
      DIMENSION ICORE(MAXCOR)
      DIMENSION I0T(2),I0R(2),I0I(2)
C
C CALCULATES CONTRIBUTION OF G(AB,CD) TO THE VIRTUAL-OCCUPIED
C INTERMEDIATE WITHOUT EXPLICIT CONSTRUCTION OR STORAGE OF
C G(AB,CD).
C
C  I(ai) <= G(af,eg)*<bf||eg>
C
C WHERE
C
C            = [RMOD(af,mn)*L(eg,mn)+RMOD(eg,mn)*L(af,mn)]*<bf||eg>
C
C
C (RMOD IS THE USUAL R VECTOR PLUS A DTAU CONTRIBUTION
C
CEND 
C
C CODED JFS DECEMBER/93 DURING JG's VISIT!
C
      COMMON/STATSYM/IRREPX
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/FLAGS/IFLAGS(100)
      COMMON/EIGPROB/ISIDE
C
      DATA ONE,HALF,HALFM/1.D0,0.5D0,-0.5D0/
C
C READ IN T1 AND R1 AMPLITUDES AND MAKE SPACE FOR THE I(AB)
C INTERMEDIATES
C
      I0T(1)=1
      I0T(2)=I0T(1)+IRPDPD(1,9)*IINTFP*IUHF
      I0R(1)=I0T(2)+IRPDPD(1,10)*IINTFP
      I0R(2)=I0R(1)+IRPDPD(IRREPX,9)*IINTFP*IUHF
      I0I(1)=I0R(2)+IRPDPD(IRREPX,10)*IINTFP
      I0I(2)=I0I(1)+IRPDPD(1,19)*IINTFP*IUHF
      ISTART=I0I(2)+IRPDPD(1,20)*IINTFP
C
      MXCOR=MAXCOR-ISTART
C
      CALL GETLST(ICORE(I0T(1)),1,1,1,1,90)
c      CALL GETLST(ICORE(I0R(1)),1,1,1,3,490)
      CALL GETLST(ICORE(I0R(1)),1,1,1,LISTR1OFF+1,LISTR1)
      CALL ZERO  (ICORE(I0I(1)),NFEA(1))
      IF(IUHF.EQ.1) THEN
       CALL GETLST(ICORE(I0T(2)),1,1,1,2,90)
c       CALL GETLST(ICORE(I0R(2)),1,1,1,4,490)
       CALL GETLST(ICORE(I0R(2)),1,1,1,LISTR1OFF+2,LISTR1)
       CALL ZERO  (ICORE(I0I(2)),NFEA(2))
      ENDIF
C
C AO-BASED ALGORITHM
C
      IF(IFLAGS(93).EQ.2) THEN
C
       ISIDE=1
       NPASS=2
       FACT=HALFM
       TAU2=.TRUE.
C
       CALL AOLADLST(IUHF,0,IRREPX)
C
       DO 10000 IPASS=1,NPASS
C
        IF(IPASS.EQ.1) THEN
C
c           CALL T2TOAO(ICORE(ISTART),MXCOR,IUHF,.FALSE.,TAU2,
c     &        460,213,IRREPX)
           CALL T2TOAO(ICORE(ISTART),MXCOR,IUHF,.FALSE.,TAU2,
     &        LISTR2-1,213,IRREPX)
C
        ELSE
C
c         CALL T2TOAO(ICORE(ISTART),MXCOR,IUHF,.FALSE.,.FALSE.,
c     &               443,213,IRREPX)
         CALL T2TOAO(ICORE(ISTART),MXCOR,IUHF,.FALSE.,.FALSE.,
     &               LISTL2-1,213,IRREPX)
C
        ENDIF
C
        IF(IFLAGS(95).EQ.1) THEN
C
c         CALL AOLAD2(ICORE(ISTART),MXCOR,IUHF,.FALSE.,IRREPX,213,463)
           CALL AOLAD2(ICORE(ISTART),MXCOR,IUHF,.FALSE.,IRREPX,213,
     &        LISTR2+2)
C
        ELSE
C
c         CALL AOLAD3(ICORE(ISTART),MXCOR,IUHF,.FALSE.,IRREPX,213,463)
           CALL AOLAD3(ICORE(ISTART),MXCOR,IUHF,.FALSE.,IRREPX,213,
     &        LISTR2+2)
C
        ENDIF
C
c        CALL Z2TOMO(ICORE(ISTART),MXCOR,IUHF,.FALSE.,IRREPX,280,280,
c     &              463,.TRUE.)
        CALL Z2TOMO(ICORE(ISTART),MXCOR,IUHF,.FALSE.,IRREPX,280,280,
     &              LISTR2+2,.TRUE.)
C
        DO 10010 ISPIN=3,3-2*IUHF,-1
C
         IF(IPASS.EQ.1) THEN
c          LIST=443+ISPIN
          LIST=LISTL2-1+ISPIN
          TAU2=.FALSE. 
         ELSE
c          LIST=460+ISPIN
          LIST=LISTR2-1+ISPIN
          TAU2=.TRUE.
         ENDIF
C
         LISTQ=280+ISPIN
C
         DO 10020 IRREP1=1,NIRREP
C
          IRREP2=DIRPRD(IRREP1,IRREPX)
          DISDSZ=IRPDPD(IRREP2,ISYTYP(1,LIST))
          NUMDSZ=IRPDPD(IRREP1,ISYTYP(2,LIST))
C
          IF(ISPIN.LE.2) THEN
C
           NVRTFUL=IRPDPD(IRREP2,18+ISPIN)
C
          ELSE
C
           NVRTFUL=DISDSZ
C
          ENDIF
C
          I010=ISTART
          I020=I010+IINTFP*NVRTFUL*NUMDSZ
          I030=I020+IINTFP*NVRTFUL*NUMDSZ
          IEND1=I030+IINTFP*DISDSZ*NUMDSZ
          ITMP1=I030
          ITMP2=ITMP1+IINTFP*MAX(NUMDSZ,DISDSZ)
          ITMP3=ITMP2+IINTFP*MAX(NUMDSZ,DISDSZ)
          IEND2=ITMP3+IINTFP*MAX(NUMDSZ,DISDSZ)
          IEND=MAX(IEND1,IEND2)
          IF(IEND.GE.MAXCOR) CALL INSMEM('IABCONT',IEND,MAXCOR)
C
          CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREP1,LIST)
C
          IF(TAU2) THEN
C
           IF(ISPIN.LE.2) THEN
C
            CALL DTAU(IRREP2,IRREP1,1,IRREPX,ICORE(I010),
     &                ICORE(I0T(ISPIN)),ICORE(I0T(ISPIN)),
     &                ICORE(I0R(ISPIN)),ICORE(I0R(ISPIN)),
     &                ISPIN,ONE)
C
           ELSE
C
            CALL DTAU(IRREP2,IRREP1,1,IRREPX,ICORE(I010),
     &                ICORE(I0T(1)),ICORE(I0T(2)),
     &                ICORE(I0R(1)),ICORE(I0R(2)),3,ONE)
C
           ENDIF
          ENDIF
C
          CALL GETTRN(ICORE(I020),ICORE(ITMP1),DISDSZ,NUMDSZ,1,
     &                IRREP1,LISTQ)
C
          IF(ISPIN.LE.2) THEN
C
           CALL SYMEXP(IRREP2,VRT(1,ISPIN),NUMDSZ,ICORE(I020))
C
          ENDIF
C
C SPIN ADAPT AMPLITUDES
C
          IF(IUHF.EQ.0) THEN
C
           CALL SPINAD1(IRREP1,POP(1,1),DISDSZ,ICORE(I010),
     &                  ICORE(ITMP1),ICORE(ITMP2))
C
          ENDIF
C
          CALL TRANSP(ICORE(I010),ICORE(I030),NUMDSZ,DISDSZ)
          CALL SCOPY(NUMDSZ*DISDSZ,ICORE(I030),1,ICORE(I010),1)
C
          IF(ISPIN.LT.3) THEN
C
           CALL SYMEXP(IRREP2,VRT(1,ISPIN),NUMDSZ,ICORE(I010))
C
           ISPIN1=ISPIN
           ISPIN2=ISPIN
           IOFFI=I0I(ISPIN)
C
          ELSE
C
           ISPIN1=1
           ISPIN2=2
c           IOFFI=I0I(1)+IINTFP*IUHF*NFEA(1)
           IOFFI=I0I(1)+IINTFP*IUHF*NFEA(1)
C
          ENDIF
C
          IOFF=0
C
          DO 10110 IRREPA=1,NIRREP
C
           IRREPE=DIRPRD(IRREPA,IRREP2)
           IRREPB=IRREPA
C
           NUME=VRT(IRREPE,ISPIN1)
           NUMA=VRT(IRREPA,ISPIN2)
           NUMB=VRT(IRREPB,ISPIN2)
C
           CALL XGEMM('T','N',NUMA,NUMB,NUMDSZ*NUME,
     &                FACT,ICORE(I020+IOFF),NUMDSZ*NUME,
     &                ICORE(I010+IOFF),NUMDSZ*NUME,ONE,
     &                ICORE(IOFFI),NUMA) 
C
           IOFF=IOFF+IINTFP*NUMDSZ*NUME*NUMA
           IOFFI=IOFFI+IINTFP*NUMA*NUMB
C
10110     CONTINUE
C 
          IF(IUHF.EQ.1.AND.ISPIN.EQ.3) THEN
C
           CALL SYMTR1(IRREP2,VRT(1,1),VRT(1,2),NUMDSZ,
     &                 ICORE(I010),ICORE(ITMP1),ICORE(ITMP2),
     &                 ICORE(ITMP3))
           CALL SYMTR1(IRREP2,VRT(1,1),VRT(1,2),NUMDSZ,
     &                 ICORE(I020),ICORE(ITMP1),ICORE(ITMP2),
     &                 ICORE(ITMP3))
C
           ISPIN1=2
           ISPIN2=1
           IOFFI=I0I(1)
C
           IOFF=0
C
           DO 10120 IRREPA=1,NIRREP
C
            IRREPE=DIRPRD(IRREPA,IRREP2)
            IRREPB=IRREPA
C
            NUME=VRT(IRREPE,ISPIN1)
            NUMA=VRT(IRREPA,ISPIN2)
            NUMB=VRT(IRREPB,ISPIN2)
C
            CALL XGEMM('T','N',NUMA,NUMB,NUMDSZ*NUME,
     &                 FACT,ICORE(I020+IOFF),NUMDSZ*NUME,
     &                 ICORE(I010+IOFF),NUMDSZ*NUME,ONE,
     &                 ICORE(IOFFI),NUMA) 
C
            IOFF=IOFF+IINTFP*NUMDSZ*NUME*NUMA
            IOFFI=IOFFI+IINTFP*NUMA*NUMB
C
10120      CONTINUE
C 
          ENDIF
C
10020    CONTINUE
C
10010   CONTINUE
C
10000  CONTINUE
C
       CALL PUTLST(ICORE(I0I(1)),1,1,1,1,92)
       IF(IUHF.EQ.1) CALL PUTLST(ICORE(I0I(2)),1,1,1,2,92)
C
      ELSE
C
C LOOP OVER SPIN CASES
C
      DO 10 ISPIN=3,3-2*IUHF,-1
c       LISTL=443+ISPIN 
       LISTL=LISTL2-1+ISPIN 
c       LISTR=460+ISPIN
       LISTR=LISTR2-1+ISPIN
       LISTW=230+ISPIN
       DO 100 IRREP1=1,NIRREP
        IRREP2=DIRPRD(IRREP1,IRREPX)
        IRREPW=IRREP2
        DISSYL=IRPDPD(IRREP2,ISYTYP(1,LISTL))
        NUMDSL=IRPDPD(IRREP1,ISYTYP(2,LISTL))
        DISSYR=IRPDPD(IRREP2,ISYTYP(1,LISTR))
        NUMDSR=IRPDPD(IRREP1,ISYTYP(2,LISTR))
        DISSYW=IRPDPD(IRREP2,ISYTYP(1,LISTW))
        NUMDSW=IRPDPD(IRREP2,ISYTYP(2,LISTW))
        DISSYQ=IRPDPD(IRREP1,ISYTYP(2,43+ISPIN))
        NUMDSQ=IRPDPD(IRREP2,ISYTYP(1,43+ISPIN))
        IF(ISPIN.LE.2)THEN
         NVRTFUL=IRPDPD(IRREP2,18+ISPIN)
        ELSE
         NVRTFUL=DISSYL
        ENDIF
        MAXT=MAX(DISSYL,DISSYR,DISSYQ,DISSYW,
     &           NUMDSL,NUMDSR,NUMDSQ,NUMDSW,
     &           NVRTFUL)
        I0Q1=ISTART
        I0Q2=I0Q1+IINTFP*DISSYQ*MAX(NVRTFUL,NUMDSQ)
        I000=I0Q2+IINTFP*DISSYQ*MAX(NVRTFUL,NUMDSQ)
        I010=I000+IINTFP*MAX(NVRTFUL,DISSYL)*NUMDSL
        I020=I010+IINTFP*MAX(NVRTFUL,DISSYR)*NUMDSR
        ITMP1=I020
        ITMP2=ITMP1+IINTFP*MAXT
        ITMP3=ITMP2+IINTFP*MAXT
        CALL GETLST(ICORE(I000),1,NUMDSL,1,IRREP1,LISTL)
        IF(IUHF.EQ.0)THEN
         CALL SPINAD1(IRREP1,POP(1,1),DISSYL,ICORE(I000),
     &                ICORE(ITMP1),ICORE(ITMP2))
        ENDIF 
        CALL GETLST(ICORE(I010),1,NUMDSR,1,IRREP1,LISTR)
        IF(ISPIN.LE.2)THEN
         CALL DTAU(IRREP2,IRREP1,1,IRREPX,ICORE(I010),
     &             ICORE(I0T(ISPIN)),ICORE(I0T(ISPIN)),
     &             ICORE(I0R(ISPIN)),ICORE(I0R(ISPIN)),ISPIN,ONE)
        ELSE
         CALL DTAU(IRREP2,IRREP1,1,IRREPX,ICORE(I010),
     &             ICORE(I0T(1)),ICORE(I0T(2)),
     &             ICORE(I0R(1)),ICORE(I0R(2)),3,ONE)
        ENDIF
        IFIRST =I000
        ISECOND=I010 
C
C FIRST CONTRACT X(E<G,M<N) WITH W(E<G,B<F)
C
C
C         Q1(mn,bf) = L+(eg,mn) *W+(bf,eg)
C         Q2(mn,bf) = R+(eg,mn) *W+(bf,eg) 
C
        CORLFT=MAXCOR-I020+1
        IF(DISSYW.NE.0)THEN
         NINCOR=CORLFT/(DISSYW*IINTFP)
        ELSE
         NINCOR=NUMDSW
        ENDIF
        NLEFT =NUMDSW
        NFIRST=1
        IOFFQ1=I0Q1
        IOFFQ2=I0Q2
        IOFF1 =IFIRST
        IOFF2 =ISECOND
        CALL ZERO(ICORE(I0Q1),DISSYQ*NUMDSQ)
        CALL ZERO(ICORE(I0Q2),DISSYQ*NUMDSQ)
1       NBATCH=MIN(NLEFT,NINCOR)
        CALL GETLST(ICORE(I020),NFIRST,NBATCH,1,IRREPW,LISTW)
        CALL XGEMM('T','T',NUMDSL,DISSYW,NBATCH,ONE,ICORE(IOFF1),
     &             DISSYL,ICORE(I020),DISSYW,ONE,ICORE(I0Q1),
     &             NUMDSL)
        CALL XGEMM('T','T',NUMDSL,DISSYW,NBATCH,ONE,ICORE(IOFF2),
     &             DISSYL,ICORE(I020),DISSYW,ONE,ICORE(I0Q2),
     &             NUMDSL)
        IOFF1=IOFF1+IINTFP*NBATCH
        IOFF2=IOFF2+IINTFP*NBATCH
        NLEFT=NLEFT-NBATCH
        NFIRST=NFIRST+NBATCH
        IF(NLEFT.NE.0)GOTO 1
C
C        I(ab) <= R(af,mn)*Q1(mn,bf) + L(af,mn)*Q2(mn,bf) [ISPIN=1,2]
C 
        CALL TRANSP(ICORE(I000),ICORE(I020),NUMDSL,DISSYL)
c YAU : old
c       CALL ICOPY(NUMDSL*DISSYL*IINTFP,ICORE(I020),1,ICORE(I000),1)
c YAU : new
        CALL DCOPY(NUMDSL*DISSYL,ICORE(I020),1,ICORE(I000),1)
c YAU : end
        CALL TRANSP(ICORE(I010),ICORE(I020),NUMDSL,DISSYL)
c YAU : old
c       CALL ICOPY(NUMDSL*DISSYL*IINTFP,ICORE(I020),1,ICORE(I010),1)
c YAU : new
        CALL DCOPY(NUMDSL*DISSYL,ICORE(I020),1,ICORE(I010),1)
c YAU : end
        IFIRST =I010
        ISECOND=I000
C
        IF(ISPIN.LE.2)THEN
C
C EXPAND RHS OF EVERYTHING 
C
         CALL SYMEXP(IRREP2,VRT(1,ISPIN),DISSYQ,ICORE(I0Q1))
         CALL SYMEXP(IRREP2,VRT(1,ISPIN),DISSYQ,ICORE(I0Q2))
         CALL SYMEXP(IRREP2,VRT(1,ISPIN),DISSYQ,ICORE(I000))
         CALL SYMEXP(IRREP2,VRT(1,ISPIN),DISSYQ,ICORE(I010))
C
C DO CONTRACTIONS 
C                           
C      I(AB) <= X1(mn,fa)*Q1(mn,fb) + X2(mn,fa)*Q2(mn,fb)
C
         IOFFQ1=I0Q1
         IOFFQ2=I0Q2
         IOFF1 =IFIRST
         IOFF2 =ISECOND
         IOFFI=I0I(ISPIN) 
         DO 110 IRREPB=1,NIRREP
          IRREPA=IRREPB
          IRREPF=DIRPRD(IRREPA,IRREP2)
          NUMA=VRT(IRREPA,ISPIN)
          NUMB=VRT(IRREPB,ISPIN)
          NUMF=VRT(IRREPF,ISPIN)
          NROW=NUMA
          NCOL=NUMB
          NSUM=NUMDSR*NUMF
          CALL XGEMM('T','N',NROW,NCOL,NSUM,-HALF,ICORE(IOFF1),NSUM,
     &               ICORE(IOFFQ1),NSUM,ONE,ICORE(IOFFI),NROW)
          CALL XGEMM('T','N',NROW,NCOL,NSUM,-HALF,ICORE(IOFF2),NSUM,
     &               ICORE(IOFFQ2),NSUM,ONE,ICORE(IOFFI),NROW)
          IOFF1=IOFF1+NROW*NSUM*IINTFP
          IOFF2=IOFF2+NROW*NSUM*IINTFP
          IOFFQ1=IOFFQ1+NCOL*NSUM*IINTFP
          IOFFQ2=IOFFQ2+NCOL*NSUM*IINTFP
          IOFFI=IOFFI+NROW*NCOL*IINTFP
110      CONTINUE
C
        ELSE
C
C DO CONTRACTION 
C                           
C      I(ab) <= X1(Mn,Fa)*Q1(Mn,Fb) + X2(Mn,Fa)*Q2(Mn,Fb) [ICASE=1]
C      I(AB) <= X1(Nm,fA)*Q1(Nm,fB) + X2(Nm,fA)*Q2(Nm,fB) [ICASE=2]
C
         DO 120 ICASE=1,1+IUHF
          IF(ICASE.EQ.2)THEN
           CALL SYMTR1(IRREP2,VRT(1,1),VRT(1,2),DISSYQ,ICORE(I0Q1),
     &                 ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
           CALL SYMTR1(IRREP2,VRT(1,1),VRT(1,2),DISSYQ,ICORE(I0Q2),
     &                 ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
           CALL SYMTR1(IRREP2,VRT(1,1),VRT(1,2),DISSYQ,ICORE(I000),
     &                 ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
           CALL SYMTR1(IRREP2,VRT(1,1),VRT(1,2),DISSYQ,ICORE(I010),
     &                 ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
          ENDIF
          IOFFQ1=I0Q1
          IOFFQ2=I0Q2
          IOFF1 =IFIRST
          IOFF2 =ISECOND
          IOFFI =I0I(3-ICASE)
          DO 210 IRREPB=1,NIRREP
           IRREPA=IRREPB
           IRREPF=DIRPRD(IRREPA,IRREP2)
           NUMA=VRT(IRREPA,3-ICASE)
           NUMB=VRT(IRREPB,3-ICASE)
           NUMF=VRT(IRREPF,ICASE)
           NROW=NUMA
           NCOL=NUMB
           NSUM=DISSYQ*NUMF
           CALL XGEMM('T','N',NROW,NCOL,NSUM,-HALF,ICORE(IOFF1),NSUM,
     &                ICORE(IOFFQ1),NSUM,ONE,ICORE(IOFFI),NROW)
           CALL XGEMM('T','N',NROW,NCOL,NSUM,-HALF,ICORE(IOFF2),NSUM,
     &                ICORE(IOFFQ2),NSUM,ONE,ICORE(IOFFI),NROW)
           IOFF1=IOFF1+NROW*NSUM*IINTFP
           IOFF2=IOFF2+NROW*NSUM*IINTFP
           IOFFQ1=IOFFQ1+NCOL*NSUM*IINTFP
           IOFFQ2=IOFFQ2+NCOL*NSUM*IINTFP
           IOFFI=IOFFI+NROW*NCOL*IINTFP
210       CONTINUE
C
120      CONTINUE
C
        ENDIF
C
100    CONTINUE
C
10    CONTINUE
C
C ALL DONE, DO A TRANSPOSITION OF I(A,B) AND SAVE IT ON DISK
C
      
       DO 300 ISPIN=1,1+IUHF
        IOFF=I0I(ISPIN)
        DO 301 IRREP=1,NIRREP
         NVRT=VRT(IRREP,ISPIN)
         CALL MTRAN2(ICORE(IOFF),NVRT)
         IOFF=IOFF+NVRT*NVRT*IINTFP
301     CONTINUE
        CALL PUTLST(ICORE(I0I(ISPIN)),1,1,1,ISPIN,92)
300    CONTINUE
C
      ENDIF
C
      RETURN
      END