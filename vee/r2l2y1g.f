      SUBROUTINE R2L2Y1G(Y1,ICORE,MAXCOR,IUHF)
C
C Y1(ai) = - 1/2 R(mn,ef)*W(mn,of)*L(io,ae)
C
C        = L(io,ae) * Q(eo)
C
C WHERE Q(eo) = R(mn,ef)*W(mn,of)
C  
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,TWO
      DIMENSION ICORE(MAXCOR),Y1(*),IAI(2)
      COMMON/STATSYM/IRREPX
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA ONE,ONEM,TWO/1.0D0,-1.0D0,2.0D0/
C
      I0G=1
      I000=I0G+IINTFP*(IRPDPD(IRREPX,9)+IUHF*IRPDPD(IRREPX,10))
C
      IOFFG=I0G
      DO 10 ISPIN=1,1+IUHF
       CALL ZERO(ICORE(IOFFG),IRPDPD(IRREPX,8+ISPIN))
       CALL PUTLST(ICORE(IOFFG),1,1,1,2+ISPIN,490)
       IOFFG=IOFFG+IINTFP*IRPDPD(IRREPX,8+ISPIN)
10    CONTINUE
      CALL DT2INT1B(ICORE(I0G),MAXCOR,IUHF,IRREPX,7,461,490)
      IOFFG=I0G
      DO 11 ISPIN=1,1+IUHF
       CALL GETLST(ICORE(IOFFG),1,1,1,2+ISPIN,490)
       IOFFG=IOFFG+IINTFP*IRPDPD(IRREPX,8+ISPIN)
11    CONTINUE
C
      IF(IUHF.EQ.0)THEN
C
C SPIN ADAPTED RHF CODE.
C
C   Y(ai) = [2 L(Io,Ae) - L(Oi,Ae)]*G(eo)
C
C
C FORM SPIN ADAPTED 2*L(Io,Ae)-L(Oi,Ae) AND PUT IN L(AI,eo) MATRIX
C
       LISTL1=437
       LISTL2=439
       NUMDSL=IRPDPD(IRREPX,ISYTYP(2,LISTL1))
       DISSYL=IRPDPD(1,ISYTYP(1,LISTL1))
       I010=I000+IINTFP*NUMDSL*DISSYL
       I020=I010+IINTFP*NUMDSL*DISSYL
       CALL GETLST(ICORE(I000),1,NUMDSL,1,IRREPX,LISTL1)
       CALL GETLST(ICORE(I010),1,NUMDSL,1,IRREPX,LISTL2)
       CALL SSCAL (NUMDSL*DISSYL,TWO,ICORE(I000),1)
       CALL SAXPY (NUMDSL*DISSYL,ONEM,ICORE(I010),1,ICORE(I000),1)
C
C EVALUTE PRODUCT Y(AI) = L(AI,eo)*Q(eo)
C
       CALL XGEMM('N','N',DISSYL,1,NUMDSL,ONE,ICORE(I000),DISSYL,
     &            ICORE(I0G),NUMDSL,ONE,Y1,DISSYL)
C
      ELSE
C
C  Y(AI) <= L(IO,AE) * Q(EO)
C
       IAI(1)=I0G
       IAI(2)=I0G+IINTFP*IRPDPD(IRREPX,9)
       IOFFY=1
       DO 12 ISPIN=1,2
        LISTL1=433+ISPIN
        NUMDSL=IRPDPD(IRREPX,ISYTYP(2,LISTL1))
        DISSYL=IRPDPD(1,ISYTYP(1,LISTL1))
        I010=I000+IINTFP*NUMDSL*DISSYL
        CALL GETLST(ICORE(I000),1,NUMDSL,1,IRREPX,LISTL1)
        CALL XGEMM('N','N',DISSYL,1,NUMDSL,ONEM,ICORE(I000),DISSYL,
     &             ICORE(IAI(ISPIN)),NUMDSL,ONE,Y1(IOFFY),DISSYL)
C
C  Y(AI) <= L(AI,eo) * Q(eo)
C
        LISTL1=438-ISPIN
        NUMDSL=IRPDPD(IRREPX,ISYTYP(2,LISTL1))
        DISSYL=IRPDPD(1,ISYTYP(1,LISTL1))
        I010=I000+IINTFP*NUMDSL*DISSYL
        CALL GETLST(ICORE(I000),1,NUMDSL,1,IRREPX,LISTL1)
        CALL XGEMM('N','N',DISSYL,1,NUMDSL,ONE,ICORE(I000),DISSYL,
     &             ICORE(IAI(3-ISPIN)),NUMDSL,ONE,Y1(IOFFY),DISSYL)
C
        IOFFY=IOFFY+IINTFP*NT(ISPIN)
C
12     CONTINUE
C
      ENDIF
C
      RETURN
      END
