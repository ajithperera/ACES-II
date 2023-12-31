      SUBROUTINE GT1INL1(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE COMPUTES THE CONTRIBUTION OF LINEAR L1 IN THE L1
C  EQUATION.
C
C L1INC(A,I) = [ SUM F G(FE) T(M,F) - SUM N G(MN) T(N,E) ] <MI||EA> 
C
C              [ SUM f G(fe) T(m,f) - SUM n G(mn) T(n,e) ] <mI||eA>  (AA)
C
C L1INC(a,i) = [ SUM f G(f,e) T(m,f) - SUM n G(mn) T(n,e) ] <mi||ea> 
C
C            + [ SUM F G(FE) T(M,F) - SUM N G(MN) T(N,E) ] <Mi|Ea> (BB)
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP,VRT,DIRPRD
      LOGICAL BREDUNDANT
      DIMENSION ICORE(MAXCOR),I0T(2),I0L(2),I0GT(2),I0GMN(2),
     &          I0GFE(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON /FLAGS2/ IFLAGS2(500)
C
      DATA AZERO,ONE,ONEM /0.D0,1.D0,-1.D0/
C
      MXCOR=MAXCOR

      BREDUNDANT=iflags2(155).eq.0
C
C  ALLOCATE FIRST MEMORY FOR T1, G, L1 AND THE PRODUCT G*T
C
      I0T(1)=MXCOR+1-IINTFP*NT(1)
      I0L(1)=I0T(1)-IINTFP*NT(1)
      I0GT(1)=I0L(1)-IINTFP*NT(1)
      I0GMN(1)=I0GT(1)-IINTFP*NF1(1)
      I0GFE(1)=I0GMN(1)-IINTFP*NF2(1)
      MXCOR=MXCOR-IINTFP*(3*NT(1)+NF1(1)+NF2(1))
      CALL GETLST(ICORE(I0T(1)),1,1,1,1,90)
      CALL GETLST(ICORE(I0L(1)),1,1,1,3,90)

      CALL GETLST(ICORE(I0GMN(1)),1,1,2,1,191)
      CALL GETLST(ICORE(I0GFE(1)),1,1,2,1,192)

      IF(IUHF.EQ.0) THEN
       I0T(2)=I0T(1)
       I0L(2)=I0L(1)
       I0GT(2)=I0GT(1)
       I0GMN(2)=I0GMN(1)
       I0GFE(2)=I0GFE(1)
      ELSE
       I0T(2)=I0GFE(1)-IINTFP*NT(2)
       I0L(2)=I0T(2)-IINTFP*NT(2)
       I0GT(2)=I0L(2)-IINTFP*NT(2)
       I0GMN(2)=I0GT(2)-IINTFP*NF1(2)
       I0GFE(2)=I0GMN(2)-IINTFP*NF2(2)
       MXCOR=MXCOR-IINTFP*(3*NT(2)+NF1(2)+NF2(2))
       CALL GETLST(ICORE(I0T(2)),1,1,1,2,90)
       CALL GETLST(ICORE(I0L(2)),1,1,1,4,90)
       CALL GETLST(ICORE(I0GMN(2)),1,1,2,2,191)
       CALL GETLST(ICORE(I0GFE(2)),1,1,2,2,192)
      ENDIF
C
C   FORM FIRST THE PRODUCTS G*T
C
       DO 1000 ISPIN=1,IUHF+1
C
        IOFFGMN=I0GMN(ISPIN)
        IOFFGFE=I0GFE(ISPIN)
        IOFFT=I0T(ISPIN)
        IOFFGT=I0GT(ISPIN)
        CALL IZERO(ICORE(I0GT(ISPIN)),NT(ISPIN)*IINTFP)
C
        DO 100 IRREP=1,NIRREP
         NOCC=POP(IRREP,ISPIN)
         NVRT=VRT(IRREP,ISPIN)
         IF(MIN(NVRT,NOCC).NE.0) THEN
          CALL XGEMM('N','N',NVRT,NOCC,NVRT,ONE,ICORE(IOFFGFE),
     &               NVRT,ICORE(IOFFT),NVRT,AZERO,ICORE(IOFFGT),
     &               NVRT)
          CALL XGEMM('N','T',NVRT,NOCC,NOCC,ONEM,ICORE(IOFFT),
     &               NVRT,ICORE(IOFFGMN),NOCC,ONE,ICORE(IOFFGT),
     &               NVRT)
         ENDIF
         IOFFGMN=IOFFGMN+IINTFP*NOCC*NOCC
         IOFFGFE=IOFFGFE+IINTFP*NVRT*NVRT
         IOFFT=IOFFT+IINTFP*NOCC*NVRT
         IOFFGT=IOFFGT+IINTFP*NOCC*NVRT
100     CONTINUE
1000   CONTINUE

C
C   NOW THE SECOND PART MULTIPLICATION WITH THE INTEGRALS
C
      DO 2000 ISPIN=1,IUHF+1
C
C  SET INTEGRAL LIST NUMBERS
C
       LISTW1=18+ISPIN
       LISTW2=16+ISPIN
       IF(IUHF.EQ.0) LISTW2=18
C
C  WE ONLY NEED TO CONSIDER ONE DPD IRREP -- THE TOTALLY SYMMETRIC ONE
C
       NTOTAR=NT(ISPIN)
       NSIZ1=NT(ISPIN)
       NSIZ2=NT(3-ISPIN)
C
C   I000 HOLDS THE INTEGRALS
C
       I000=1 
       I001=I000+IINTFP*NSIZ1*MAX(NSIZ1,NSIZ2)
       IF(MXCOR.LE.I001) STOP 'GT1INL1'
C
C  GET INTEGRALS FOR LISTW1 FROM LIST
C

       IF(BREDUNDANT) THEN
          CALL GETLST(ICORE(I000),1,NTOTAR,2,1,LISTW1)
       ELSE
          CALL GETLST_NR(ICORE(I000),ICORE(I001),MXCOR-I001,
     &                  LISTW1,1)
       ENDIF
c       call chksums("gt1inl1 listw1 1: ",icore(i000),ntotar*nsiz1)

       CALL XGEMM('N','N',1,NTOTAR,NSIZ1,ONE,ICORE(I0GT(ISPIN)),
     &            1,ICORE(I000),NSIZ1,ONE,ICORE(I0L(ISPIN)),1)

       IF(BREDUNDANT) THEN
          CALL GETLST(ICORE(I000),1,NTOTAR,2,1,LISTW2)
       ELSE
          CALL GETLST_NR(ICORE(I000),ICORE(I001),MXCOR-I001,
     &                  LISTW2,1)
       ENDIF
c       call chksums("gt1inl1 listw2 1: ",icore(i000),ntotar*nsiz1)

       CALL XGEMM('N','N',1,NTOTAR,NSIZ2,ONE,ICORE(I0GT(3-ISPIN)),
     &            1,ICORE(I000),NSIZ2,ONE,ICORE(I0L(ISPIN)),1)

       CALL PUTLST(ICORE(I0L(ISPIN)),1,1,1,2+ISPIN,90)

2000  CONTINUE
C
      RETURN
      END
