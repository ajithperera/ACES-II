
      SUBROUTINE ORBANAL(SCR,MAXCOR,IUHF)
C
C DRIVER FOR FORMATION OF EXCITED STATE DENSITY MATRICES
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP,VRT,DIRPRD
      CHARACTER*2 CART(6)
      CHARACTER*5 LABSPN(2)
      CHARACTER*8 LABEL2(6),LABEL(2)
      DIMENSION SCR(MAXCOR),ILOC(6),SECMOM(6)
      COMMON/FLAGS/IFLAGS(100)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON/EXTINF3/IROOT,LOCROOT,ITROOT
C
      DATA ONE,ONEM,ZILCH,HALF /1.0D0,-1.0D0,0.0D0,0.5D0/
      DATA LABEL/'SCFEVALA','SCFEVALB'/
      DATA LABSPN/'alpha','beta '/
      DATA LABEL2/'2NDMO_XX','2NDMO_YY','2NDMO_ZZ',
     &            '2NDMO_XY','2NDMO_XZ','2NDMO_YZ'/
      DATA CART /'XX','YY','ZZ','XY','XZ','YZ'/
C
      NNP1O2(I)=(I*(I+1))/2
C
      IONE=1
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NAO)
      NMO=NOCCO(1)+NVRTO(1)
      I000=1
      I010=I000+NAO*NAO
      I020=I010+NAO*NAO
      I030=I020+NAO*NAO
      I040=I030+NAO*NAO
      ILOC(1)=I040+NAO*NAO
      ILOC(2)=ILOC(1)+NAO*NAO
      ILOC(3)=ILOC(2)+NAO*NAO
      ILOC(4)=ILOC(3)+NAO*NAO
      ILOC(5)=ILOC(4)+NAO*NAO
      ILOC(6)=ILOC(5)+NAO*NAO
      ITOP=ILOC(6)+NAO*NAO
      LENGTH=NNP1O2(NAO)
      DO 20 IXYZ=1,6
       CALL GETREC(20,'JOBARC',LABEL2(IXYZ),LENGTH*IINTFP,
     &             SCR(ILOC(IXYZ)))
       CALL EXPND2(SCR(ILOC(IXYZ)),SCR(I030),NAO)
       CALL SCOPY (NAO*NAO,SCR(I030),1,SCR(ILOC(IXYZ)),1)
20    CONTINUE
C 
      DO 9 ISPIN=1,1+IUHF
       WRITE(6,1005)LABSPN(ISPIN)
       WRITE(6,*)
1005   FORMAT(T3,'Summary of active ',A5,' molecular orbitals:')
       IORB=0
       CALL GETREC(20,'JOBARC',LABEL(ISPIN),NMO*IINTFP,SCR(ITOP))
       WRITE(6,1001)
       WRITE(6,1002)
       WRITE(6,1001)
       DO 11 isym=1,nirrep
        DO 12 j=1,pop(isym,ispin) 
         IORB=IORB+1
         CALL ZERO(SCR(I000),NMO*NMO)
         SCR(IORB+(IORB-1)*NMO)=ONE
         CALL MO2AO3(SCR(I000),SCR(I010),SCR(I020),SCR(I030),
     &               NAO,NMO,ISPIN)
         DO 13 IXYZ=1,6
          SECMOM(IXYZ)=
     &            SDOT(NAO*NAO,SCR(ILOC(IXYZ)),1,SCR(I010),1)
13       CONTINUE
         WRITE(6,1000)IORB,ISYM,SCR(ITOP-1+IORB),(SECMOM(K),K=1,6)
12      CONTINUE
11     CONTINUE
       DO 111 ISYM=1,NIRREP
        DO 112 J=1,VRT(ISYM,ISPIN)
         IORB=IORB+1
         CALL ZERO(SCR(I000),NMO*NMO)
         SCR(IORB+(IORB-1)*NMO)=ONE
         CALL MO2AO3(SCR(I000),SCR(I010),SCR(I020),SCR(I030),
     &               NAO,NMO,ISPIN)
         DO 113 IXYZ=1,6
          SECMOM(IXYZ)=
     &            SDOT(NAO*NAO,SCR(ILOC(IXYZ)),1,SCR(I010),1)
113      CONTINUE
         WRITE(6,1000)IORB,ISYM,SCR(ITOP-1+IORB),(SECMOM(K),K=1,6)
112     CONTINUE
111    CONTINUE
C
9     CONTINUE
      WRITE(6,1001)
C
      RETURN
1002  FORMAT('Orbital',T9,'Symm.',T15,'Eigenvalue',T27,'<xx>',
     &       T36,'<yy>',T45,'<zz>',T54,'<xy>',T63,'<xz>',T72,'<yz>')
1001  FORMAT(78('-'))
1000  FORMAT(T3,I3,T11,I1,T14,F11.4,T25,6(F8.3,1X))
      END 
