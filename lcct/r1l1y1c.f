      SUBROUTINE R1L1Y1C(L1,Y1,ICORE,MAXCOR,IUHF)
C
C Y1(a,i) = R(me)*L(e,n)*W(mi,na) w/ the VV component too!
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,TWO,HALF
      CHARACTER*4 SPCASE(2)
      DIMENSION ICORE(MAXCOR),L1(*),Y1(*)
      COMMON/STATSYM/IRREPX
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA ONE,ONEM,TWO,ZILCH/1.0D0,-1.0D0,2.0D0,0.0D0/
      DATA HALF/0.5D0/
      DATA SPCASE/'AAAA','BBBB'/
C
      IF(IUHF.EQ.0)THEN
C
C SPIN ADAPTED RHF CODE.
C
C  Y(b,j) = L(EM)*[2*X(Eb,Mj)-X(Eb,Jm)]
C
C FIRST FORM 2*W(Eb,Mj) - W(Eb,Jm) AND STORE IN A W(EM,bj) MATRIX
C
       LISTW1=428
       LISTW2=426
       DISSYW=IRPDPD(IRREPX,ISYTYP(1,LISTW1))
       NUMDSW=IRPDPD(IRREPX,ISYTYP(2,LISTW1))
       I010=I000+IINTFP*DISSYW*NUMDSW
       I020=I010+IINTFP*DISSYW*NUMDSW
       IF(I020.GE.MAXCOR) CALL INSMEM('R1L0Y1',I020,MAXCOR)
       CALL GETLST(ICORE(I000),1,NUMDSW,1,IRREPX,LISTW1)
       CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPX,LISTW2)
       CALL SSCAL (NUMDSW*DISSYW,TWO,ICORE(I000),1)
       CALL SAXPY (NUMDSW*DISSYW,ONEM,ICORE(I010),1,ICORE(I000),1)
C
C CALCULATE L(EM)*W(EM,bj)=Y(bj)
C
       CALL XGEMM('N','N',1,NUMDSW,DISSYW,ONE,L1,
     &            1,ICORE(I000),DISSYW,ONE,Y1,1)
C
      ELSE
C
C UHF CODE, FIRST THE AAAA AND BBBB SPIN CASES
C
C  Y(A,I) = R(EM) * W(EA,MI) + R(em) * W(eA,mI)
C
       DO 100 ISPIN=1,2
C
        IOFFY = 1 + (ISPIN-1)*NT(1)*IINTFP
        LISTW1=423+ISPIN
        LISTW2=425+ISPIN
        DISSYW1=IRPDPD(IRREPX,ISYTYP(1,LISTW1))
        NUMDSW1=IRPDPD(IRREPX,ISYTYP(2,LISTW1))
        DISSYW2=IRPDPD(IRREPX,ISYTYP(1,LISTW2))
        NUMDSW2=IRPDPD(IRREPX,ISYTYP(2,LISTW2))

        IEND=I000+IINTFP*MAX(NUMDSW1*DISSYW1,NUMDSW2*DISSYW2)
        IF(IEND.GE.MAXCOR) CALL INSMEM('R1L0Y1',IEND,MAXCOR)
C
C CALCULATE L(EM)*W(EM,BJ) = Y(BJ)
C
        CALL GETLST(ICORE(I000),1,NUMDSW1,1,IRREPX,LISTW1)
C
        IOFFL = 1 + (ISPIN-1)*NT(1)*IINTFP
        CALL XGEMM('N','N',1,NUMDSW1,DISSYW1,ONE,L1(IOFFL),1,
     &             ICORE(I000),DISSYW1,ONE,Y1(IOFFY),1)
C
C CALCULATE L(em)*W(em,BJ) => Y(BJ)
C
        CALL GETLST(ICORE(I000),1,NUMDSW2,1,IRREPX,LISTW2)
C
        IOFFL = 1 + (2-ISPIN)*NT(1)*IINTFP
        CALL XGEMM('N','N',1,NUMDSW2,DISSYW2,ONE,L1(IOFFL),1,
     &             ICORE(I000),DISSYW2,ONE,Y1(IOFFY),1)

100    CONTINUE
C
      ENDIF
      RETURN
      END