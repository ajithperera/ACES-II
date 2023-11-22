      SUBROUTINE R1L1Y1F(F1,L1,R1,Y1,ICORE,MAXCOR,IUHF)
C
C Y1(ab,ij) = R(me)*F(ma)*L(e,i)
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      DIMENSION ICORE(MAXCOR),R1(*),F1(*),I0X(2),I0F(2),I0R(2),L1(*),
     &          Y1(*)
      COMMON/STATSYM/IRREPX
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMLOC/ISYMOFF(8,8,25)
C
      DATA ONE,ONEM,ZILCH/1.0D0,-1.0D0,0.0D0/
C
C CALCULATE X(ea) = R(me)*F(ma) INTERMEDIATES
C
      I0X(1)=1
      I0F(1)=1
      I0R(1)=1
      IF(IUHF.NE.0)THEN
       I0X(2)=I0X(1)+IINTFP*IRPDPD(IRREPX,19)
       I0F(2)=I0F(1)+IINTFP*NT(1)
       I0R(2)=I0R(1)+IINTFP*IRPDPD(IRREPX,9)
      ELSE
       I0X(2)=I0X(1)
       I0F(2)=I0F(1)
       I0R(2)=I0R(1)
      ENDIF
      I000=I0X(2)+IINTFP*IRPDPD(IRREPX,20)
      DO 5 ISPIN=1,1+IUHF
       DO 6 IRREPA=1,NIRREP
        IRREPE=DIRPRD(IRREPA,IRREPX)
        IRREPM=IRREPA
        NUMA=VRT(IRREPA,ISPIN)
        NUME=VRT(IRREPE,ISPIN)
        NUMM=POP(IRREPM,ISPIN)
        NROW=NUME
        NCOL=NUMA
        NSUM=NUMM
        IOFFR=I0R(ISPIN)+IINTFP*(ISYMOFF(IRREPM,IRREPX,8+ISPIN)-1)
        IOFFF=I0F(ISPIN)+IINTFP*(ISYMOFF(IRREPM,1,8+ISPIN)-1)
        IOFFX=I0X(ISPIN)+IINTFP*(ISYMOFF(IRREPA,IRREPX,18+ISPIN)-1)
        CALL XGEMM('N','T',NROW,NCOL,NSUM,ONE,R1(IOFFR),NROW,
     &             F1(IOFFF),NCOL,ZILCH,ICORE(IOFFX),NROW)
6      CONTINUE
5     CONTINUE

      IOFFL=1
      IOFFY=1
      DO 10 ISPIN=1,1+IUHF
       IOFFX = I0X(ISPIN)
C
C FORM PRODUCT Y1(ai) = -X(ea)*L(ei)
C
       DO 11 IRREP=1,NIRREP
        NUMA=VRT(IRREP,ISPIN)
        NUMI=POP(IRREP,ISPIN)
        CALL XGEMM('T','N',NUMA,NUMI,NUMA,ONEM,ICORE(IOFFX),NUMA,
     &             L1(IOFFL),NUMA,ONE,Y1(IOFFY),NUMA)
        IOFFL=IOFFL+IINTFP*NUMA*NUMI
        IOFFX=IOFFX+IINTFP*NUMA*NUMA
        IOFFY=IOFFY+IINTFP*NUMA*NUMI
11     CONTINUE
C
10    CONTINUE
C
      RETURN
      END
