      SUBROUTINE R1L2Y2E(R1,ICORE,MAXCOR,IUHF,LISTL2)
C
C Y1(ab,ij) = P(ab) L(af,ij)*R(me)*W(be,fm)
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ZILCH
      DIMENSION ICORE(MAXCOR),R1(*),I0X(2)
      COMMON/STATSYM/IRREPX
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMLOC/ISYMOFF(8,8,25)
C
      DATA ONE,ZILCH/1.0D0,0.0D0/
C
C FORM INTERMEDIATE X(bf) = R(em)*W(be,fm)
C
      I0X(1)=1
      IF(IUHF.NE.0)THEN
       I0X(2)=I0X(1)+IINTFP*IRPDPD(IRREPX,19)
      ELSE
       I0X(2)=I0X(1)
      ENDIF
      I000=I0X(2)+IINTFP*IRPDPD(IRREPX,20)
      MXCOR = MAXCOR - I000 + 1
      LENAB=IRPDPD(IRREPX,19)+IUHF*IRPDPD(IRREPX,20)
      CALL ZERO(ICORE(I0X(1)),LENAB) 
      CALL DVOINVV(IRREPX,ICORE(I0X(1)),R1,ICORE(I000),MXCOR,IUHF,27)
C
      IF(IUHF.EQ.0)THEN
C
C SPIN ADAPTED RHF CODE.
C
C   Y(Ab,Ij) = Q(Ab,Ij) + Q(Ba,Ji)
C
C WHERE
C
C   Q(ab,ij) = L(Af,Ij) * R(me)* [2W(Be,Fm) - W(Eb,Fm)]
C
C            = L(Af,Ij) * X(bf)
C
C FORM PRODUCT X(bf) = R(em) * W(Be,Fm)
C
C
       LISTL=LISTL2+2
       LISTZ=63
       DO 10 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPX)
        DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
        NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
        DISSYZ=IRPDPD(IRREPR,ISYTYP(1,LISTZ))
        NUMDSZ=IRPDPD(IRREPR,ISYTYP(2,LISTZ))
        MAXSIZ=MAX(DISSYL,NUMDSL,DISSYZ,NUMDSZ)
        I010=I000+IINTFP*MAX(DISSYL*NUMDSL,DISSYZ*NUMDSZ)
        I020=I010+IINTFP*MAX(DISSYL*NUMDSL,DISSYZ*NUMDSZ)
        ITMP1=I010
        ITMP2=ITMP1+IINTFP*MAXSIZ
        ITMP3=ITMP2+IINTFP*MAXSIZ
        CALL GETLST(ICORE(I010),1,NUMDSL,1,IRREPR,LISTL)
        CALL TRANSP(ICORE(I010),ICORE(I000),NUMDSL,DISSYL)
C
C                            +
C Y(Ij,Ab) = L(Ij,Af) * X(bf)
C
        DO 20 IRREPF=1,NIRREP
         IRREPA=DIRPRD(IRREPF,IRREPL)
         IRREPB=DIRPRD(IRREPF,IRREPX)
         NUMF=VRT(IRREPF,1)
         NUMA=VRT(IRREPA,1)
         NUMB=VRT(IRREPB,1)
         NROW=NUMDSL*NUMA
         NCOL=NUMB
         NSUM=NUMF
         IOFFL=I000+IINTFP*(ISYMOFF(IRREPF,IRREPL,19)-1)*NUMDSL
         IOFFZ=I010+IINTFP*(ISYMOFF(IRREPB,IRREPR,19)-1)*NUMDSL
         IOFFX=I0X(1)+IINTFP*(ISYMOFF(IRREPF,IRREPX,19)-1)
         CALL XGEMM('N','T',NROW,NCOL,NSUM,ONE,ICORE(IOFFL),NROW,
     &              ICORE(IOFFX),NCOL,ZILCH,ICORE(IOFFZ),NROW)
20      CONTINUE
C
        CALL TRANSP(ICORE(I010),ICORE(I000),DISSYZ,NUMDSZ)
        CALL SYMRHF(IRREPR,VRT(1,1),POP(1,1),DISSYZ,ICORE(I000),
     &              ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
        CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPR,LISTZ)
        CALL SAXPY (NUMDSZ*DISSYZ,ONE,ICORE(I010),1,ICORE(I000),1)
        CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPR,LISTZ)
C
10     CONTINUE
C
      ELSE
C
C UHF:
C 
C Y1(AB,IJ) = L(AF,IJ)*X(BF) + ANTISYMMETRIZATION
C Y1(Ab,Ij) = L(Af,Ij)*X(bf) + L(Fb,Ij)*X(AF)
C
C
C AAAA AND BBBB:
C
       DO 100 ISPIN=1,2
        DO 110 IRREPR=1,NIRREP
         LISTL=LISTL2-1+ISPIN
         LISTZ=60+ISPIN
         IRREPL=DIRPRD(IRREPR,IRREPX)
         DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
         DISSYLX=IRPDPD(IRREPL,18+ISPIN)
         NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
         DISSYZ=IRPDPD(IRREPR,ISYTYP(1,LISTZ))
         DISSYZX=IRPDPD(IRREPR,18+ISPIN)
         NUMDSZ=IRPDPD(IRREPR,ISYTYP(2,LISTZ))
         I010=I000+IINTFP*MAX(DISSYLX*NUMDSL,DISSYZX*NUMDSZ)
         I020=I010+IINTFP*MAX(DISSYLX*NUMDSL,DISSYZX*NUMDSZ)
         CALL GETLST(ICORE(I010),1,NUMDSL,1,IRREPR,LISTL)
         CALL TRANSP(ICORE(I010),ICORE(I000),NUMDSL,DISSYL)
         CALL SYMEXP(IRREPL,VRT(1,ISPIN),NUMDSL,ICORE(I000))
C                           
C Y(IJ,AB) = L(IJ,AF) * X(BF)
C
         DO 120 IRREPF=1,NIRREP
          IRREPA=DIRPRD(IRREPF,IRREPL)
          IRREPB=DIRPRD(IRREPF,IRREPX)
          NUMF=VRT(IRREPF,ISPIN)
          NUMA=VRT(IRREPA,ISPIN)
          NUMB=VRT(IRREPB,ISPIN)
          NROW=NUMDSL*NUMA
          NCOL=NUMB
          NSUM=NUMF
          IT=18+ISPIN
          IOFFL=I000+IINTFP*(ISYMOFF(IRREPF,IRREPL,IT)-1)*NUMDSL
          IOFFZ=I010+IINTFP*(ISYMOFF(IRREPB,IRREPR,IT)-1)*NUMDSL
          IOFFX=I0X(ISPIN)+IINTFP*(ISYMOFF(IRREPF,IRREPX,IT)-1)
          CALL XGEMM('N','T',NROW,NCOL,NSUM,ONE,ICORE(IOFFL),NROW,
     &               ICORE(IOFFX),NCOL,ZILCH,ICORE(IOFFZ),NROW)
120      CONTINUE
C
         CALL ASSYM2(IRREPR,VRT(1,ISPIN),NUMDSZ,ICORE(I010))
         CALL TRANSP(ICORE(I010),ICORE(I000),DISSYZ,NUMDSZ)
         CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPR,LISTZ)
         CALL SAXPY (NUMDSZ*DISSYZ,ONE,ICORE(I010),1,ICORE(I000),1)
         CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPR,LISTZ)
C
110     CONTINUE
C
100    CONTINUE
C
C ABAB
C
C Y1(Ab,Ij) = L(Af,Ij)*X(bf) + L(Fb,Ij)*X(AF)
C
       LISTL=LISTL2+2
       LISTZ=63
       DO 210 IRREPR=1,NIRREP
C
C Y1(Ab,Ij) <= L(Af,Ij)*X(bf)
C
        IRREPL=DIRPRD(IRREPR,IRREPX)
        DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
        NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
        DISSYZ=IRPDPD(IRREPR,ISYTYP(1,LISTZ))
        NUMDSZ=IRPDPD(IRREPR,ISYTYP(2,LISTZ))
        MAXSIZ=MAX(DISSYL,NUMDSL,DISSYZ,NUMDSZ)
        I010=I000+IINTFP*MAX(DISSYL*NUMDSL,DISSYZ*NUMDSZ)
        I020=I010+IINTFP*MAX(DISSYL*NUMDSL,DISSYZ*NUMDSZ)
        ITMP1=I020
        ITMP2=ITMP1+IINTFP*MAXSIZ
        ITMP3=ITMP2+IINTFP*MAXSIZ
        CALL GETLST(ICORE(I010),1,NUMDSL,1,IRREPR,LISTL)
        CALL TRANSP(ICORE(I010),ICORE(I000),NUMDSL,DISSYL)
C
C Y1(Ij,Ab) <= L(Ij,Af)*X(bf)
C
        DO 220 IRREPF=1,NIRREP
         IRREPA=DIRPRD(IRREPF,IRREPL)
         IRREPB=DIRPRD(IRREPF,IRREPX)
         NUMF=VRT(IRREPF,2)
         NUMA=VRT(IRREPA,1)
         NUMB=VRT(IRREPB,2)
         NROW=NUMDSL*NUMA
         NCOL=NUMB
         NSUM=NUMF
         IOFFL=I000+IINTFP*(ISYMOFF(IRREPF,IRREPL,13)-1)*NUMDSL
         IOFFZ=I010+IINTFP*(ISYMOFF(IRREPB,IRREPR,13)-1)*NUMDSL
         IOFFX=I0X(2)+IINTFP*(ISYMOFF(IRREPF,IRREPX,20)-1)
         CALL XGEMM('N','T',NROW,NCOL,NSUM,ONE,ICORE(IOFFL),NROW,
     &              ICORE(IOFFX),NCOL,ZILCH,ICORE(IOFFZ),NROW)
220     CONTINUE
C
C Y1(Ij,Ab) -> Y1(Ij,bA)
C L(Ij,Fb) -> L(Ij,bF)
C
        CALL SYMTR1(IRREPR,VRT(1,1),VRT(1,2),NUMDSZ,ICORE(I010),
     &              ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
        CALL SYMTR1(IRREPL,VRT(1,1),VRT(1,2),NUMDSL,ICORE(I000),
     &              ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
C
C Y1(Ij,bA) <= L(Ij,bF)*X(AF)
C
        DO 230 IRREPF=1,NIRREP
         IRREPB=DIRPRD(IRREPF,IRREPL)
         IRREPA=DIRPRD(IRREPF,IRREPX)
         NUMF=VRT(IRREPF,1)
         NUMA=VRT(IRREPA,1)
         NUMB=VRT(IRREPB,2)
         NROW=NUMDSL*NUMB
         NCOL=NUMA
         NSUM=NUMF
         IOFFL=I000+IINTFP*(ISYMOFF(IRREPF,IRREPL,23)-1)*NUMDSL
         IOFFZ=I010+IINTFP*(ISYMOFF(IRREPA,IRREPR,23)-1)*NUMDSL
         IOFFX=I0X(1)+IINTFP*(ISYMOFF(IRREPF,IRREPX,19)-1)
         CALL XGEMM('N','T',NROW,NCOL,NSUM,ONE,ICORE(IOFFL),NROW,
     &              ICORE(IOFFX),NCOL,ONE,ICORE(IOFFZ),NROW)
230     CONTINUE
C
C Y1(Ij,bA) -> Y1(Ij,Ab) -> Y1(Ab,Ij)
C
        CALL SYMTR1(IRREPR,VRT(1,2),VRT(1,1),NUMDSZ,ICORE(I010),
     &              ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
        CALL TRANSP(ICORE(I010),ICORE(I000),DISSYZ,NUMDSZ)
        CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPR,LISTZ)
        CALL SAXPY (NUMDSZ*DISSYZ,ONE,ICORE(I010),1,ICORE(I000),1)
        CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPR,LISTZ)
C
210    CONTINUE
C
      ENDIF
C
      RETURN
      END    