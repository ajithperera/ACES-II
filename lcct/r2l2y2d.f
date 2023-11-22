      SUBROUTINE R2L2Y2D(ICORE,MAXCOR,IUHF,LISTL2)
C
C Y1(ab,ij) = - 1/2 P(a,b) R(mn,ef)*W(mn,bf)*L(ae,ij)
C
C  RHF : SPIN ADAPTED CODE
C  UHF : ALL SPIN CASES ARE EXPLICITELY CONSIDERED
C
CEND
C
C CODED SEPTEMBER/93 JFS AND JG
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ZILCH
C
      DIMENSION ICORE(MAXCOR)
      DIMENSION I00G(2)
C
      COMMON/STATSYM/IRREPX
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMLOC/ISYMOFF(8,8,25)
C
      DATA ONE,ZILCH/1.0D0,0.0D0/
C
      IF(IUHF.EQ.0)THEN
       LISTL=LISTL2+2
       LISTZ=63
       I0G=1
       I000=I0G+IINTFP*IRPDPD(IRREPX,19)
C
C SPIN ADAPTED RHF CODE.
C
C  Y1(Ab,Ij) = Z(Ab,Ij) + Z(Ba,Ji)
C
C WHERE
C
C   Z(Ab,Ij) = [2 R(Mn,Ef)-R(Nm,Ef)]*W(Mn,Af)*L(Eb,Ij)
C
C
C FIRST PICK UP G(EA) = -W(Mn,Af)*[2*R(Mn,Ef)-R(Nm,Ef)] FROM LIST 492
C
       CALL GETLST(ICORE(I0G),1,1,1,1,492)
C
C FORM PRODUCT Z(Ab,Ij)=L(Eb,Ij)*G(EA)
C
       DO 10 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPX)
        DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
        NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
        DISSYZ=IRPDPD(IRREPR,ISYTYP(1,LISTZ))
        NUMDSZ=IRPDPD(IRREPR,ISYTYP(2,LISTZ))
        MAXSIZ=MAX(DISSYZ*NUMDSZ,DISSYL*NUMDSL)
        MAXT  =MAX(DISSYZ,NUMDSZ,DISSYL,NUMDSL)
        I010=I000+IINTFP*MAXSIZ
        I020=I010+IINTFP*MAXSIZ
        ITMP1=I020
        ITMP2=ITMP1+IINTFP*MAXT
        ITMP3=ITMP2+IINTFP*MAXT
        ITMP4=ITMP3+IINTFP*MAXT
        CALL GETLST(ICORE(I010),1,NUMDSL,1,IRREPR,LISTL)
        CALL TRANSP(ICORE(I010),ICORE(I000),NUMDSL,DISSYL)
        CALL SYMTR1(IRREPL,VRT(1,1),VRT(1,1),NUMDSL,ICORE(I000),
     &              ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
C
C NOW FORM Z(Ij,bA) = L(Ij,bE)*G(EA)
C
        DO 20 IRREPE=1,NIRREP
         IRREPB=DIRPRD(IRREPE,IRREPL)
         IRREPA=DIRPRD(IRREPE,IRREPX)
         NUME=VRT(IRREPE,1)
         NUMB=VRT(IRREPB,1)
         NUMA=VRT(IRREPA,1)
         NROW=NUMDSL*NUMB
         NCOL=NUMA
         NSUM=NUME
         IOFFL=I000+(ISYMOFF(IRREPE,IRREPL,19)-1)*NUMDSL*IINTFP
         IOFFZ=I010+(ISYMOFF(IRREPA,IRREPR,19)-1)*NUMDSZ*IINTFP
         IOFFG=I0G+(ISYMOFF(IRREPA,IRREPX,19)-1)*IINTFP
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFL),NROW,
     &              ICORE(IOFFG),NSUM,ZILCH,ICORE(IOFFZ),NROW)
20      CONTINUE
        CALL SYMTR1(IRREPR,VRT(1,1),VRT(1,1),NUMDSZ,ICORE(I010),
     &              ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
        CALL TRANSP(ICORE(I010),ICORE(I000),DISSYZ,NUMDSZ)
C
        CALL SYMRHF(IRREPR,VRT(1,1),POP(1,1),DISSYZ,ICORE(I000),
     &              ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
C
        CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPR,LISTZ)
        CALL SAXPY (NUMDSZ*DISSYZ,ONE,ICORE(I000),1,ICORE(I010),1) 
        CALL PUTLST(ICORE(I010),1,NUMDSZ,1,IRREPR,LISTZ)
C
10     CONTINUE
C
      ELSE
C
       I00G(1)=1
       I00G(2)=I00G(1)+IINTFP*IRPDPD(IRREPX,19)
       ISTART=I00G(2)+IINTFP*IRPDPD(IRREPX,20)
       CALL GETLST(ICORE(I00G(1)),1,1,1,1,492)
       CALL GETLST(ICORE(I00G(2)),1,1,1,2,492)
C
C DO FIRST SPIN CASES AAAA AND BBBB 
C
       DO 1000 ISPIN=1,2
C
        LISTL=LISTL2-1+ISPIN
        LISTZ=60+ISPIN
C
C LOOP OVER IRREPS
C 
        DO 100 IRREPR=1,NIRREP
C
         IRREPL=DIRPRD(IRREPR,IRREPX)
         NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
         DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
         NUMDSZ=IRPDPD(IRREPR,ISYTYP(2,LISTZ))
         DISSYZ=IRPDPD(IRREPR,ISYTYP(1,LISTZ))
         NVRTSQL=IRPDPD(IRREPL,18+ISPIN)
         NVRTSQR=IRPDPD(IRREPR,18+ISPIN)
C
         I000=ISTART
         I010=I000+IINTFP*MAX(NUMDSL*NVRTSQL,NUMDSZ*DISSYZ)
         ITMP=I010+IINTFP*NUMDSZ*NVRTSQR
         IEND=ITMP+IINTFP*3*MAX(NUMDSZ,NUMDSL,NVRTSQL,NVRTSQR)
         IF(IEND.GE.MAXCOR) CALL INSMEM('R2L2Y2D',IEND,MAXCOR)
C
         CALL GETTRN(ICORE(I000),ICORE(ITMP),DISSYL,NUMDSL,1,
     &               IRREPR,LISTL)
C
         CALL SYMEXP(IRREPL,VRT(1,ISPIN),NUMDSL,ICORE(I000))
C
         DO 120 IRREPE=1,NIRREP
C
          IRREPA=DIRPRD(IRREPE,IRREPL)
          IRREPB=DIRPRD(IRREPX,IRREPE) 
          NUME=VRT(IRREPE,ISPIN)
          NUMA=VRT(IRREPA,ISPIN)
          NUMB=VRT(IRREPB,ISPIN)
C
          NROW=NUMDSL*NUMA
          NCOL=NUMB
          NSUM=NUME
C
          IOFFL=I000+IINTFP*NUMDSL*
     &               (ISYMOFF(IRREPE,IRREPL,18+ISPIN)-1)
          IOFFZ=I010+IINTFP*NUMDSZ*
     &               (ISYMOFF(IRREPB,IRREPR,18+ISPIN)-1)
          IOFFG=I00G(ISPIN)+IINTFP*
     &               (ISYMOFF(IRREPB,IRREPX,18+ISPIN)-1)
C
          CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFL),
     &               NROW,ICORE(IOFFG),NSUM,ZILCH,ICORE(IOFFZ),
     &               NROW)
C
120      CONTINUE
C
         CALL ASSYM(IRREPR,VRT(1,ISPIN),NUMDSZ,NUMDSZ,
     &              ICORE(I000),ICORE(I010))
C    
         CALL TRANSP(ICORE(I000),ICORE(I010),DISSYZ,NUMDSZ)
C
         CALL GETLST(ICORE(I000),1,NUMDSZ,1,IRREPR,LISTZ)
         CALL SAXPY(NUMDSZ*DISSYZ,ONE,ICORE(I010),1,ICORE(I000),1)
         CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPR,LISTZ)
C
100     CONTINUE
C
1000   CONTINUE
C
C DO ABAB SPIN CASE
C
       LISTL=LISTL2+2
       LISTZ=63
C
C LOOP OVER IRREPS
C 
       DO 1100 IRREPR=1,NIRREP
C
        IRREPL=DIRPRD(IRREPR,IRREPX)
        NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
        DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
        NUMDSZ=IRPDPD(IRREPR,ISYTYP(2,LISTZ))
        DISSYZ=IRPDPD(IRREPR,ISYTYP(1,LISTZ))
C
        I000=ISTART
        I010=I000+IINTFP*MAX(NUMDSL*DISSYL,NUMDSZ*DISSYZ)
        ITMP=I010+IINTFP*NUMDSZ*DISSYZ
        IEND=ITMP+IINTFP*3*MAX(NUMDSZ,NUMDSL,DISSYL,DISSYZ)
        IF(IEND.GE.MAXCOR) CALL INSMEM('R2L2Y2D',IEND,MAXCOR)
C
        CALL GETTRN(ICORE(I000),ICORE(ITMP),DISSYL,NUMDSL,1,
     &              IRREPR,LISTL)
C
        DO 1120 IRREPE=1,NIRREP
C
         IRREPA=DIRPRD(IRREPE,IRREPL)
         IRREPB=DIRPRD(IRREPX,IRREPE) 
         NUME=VRT(IRREPE,2)
         NUMA=VRT(IRREPA,1)
         NUMB=VRT(IRREPB,2)
C
         NROW=NUMDSL*NUMA
         NCOL=NUMB
         NSUM=NUME
C
         IOFFL=I000+IINTFP*NUMDSL*
     &              (ISYMOFF(IRREPE,IRREPL,13)-1)
         IOFFZ=I010+IINTFP*NUMDSZ*
     &              (ISYMOFF(IRREPB,IRREPR,13)-1)
         IOFFG=I00G(2)+IINTFP*
     &              (ISYMOFF(IRREPB,IRREPX,20)-1)
C
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFL),
     &              NROW,ICORE(IOFFG),NSUM,ZILCH,ICORE(IOFFZ),
     &              NROW)
C
1120    CONTINUE
C
        CALL SYMTR1(IRREPL,VRT(1,1),VRT(1,2),NUMDSL,ICORE(I000),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*NUMDSL),
     &              ICORE(ITMP+2*IINTFP*NUMDSL))
        CALL SYMTR1(IRREPR,VRT(1,1),VRT(1,2),NUMDSZ,ICORE(I010),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*NUMDSZ),
     &              ICORE(ITMP+2*IINTFP*NUMDSZ))

        DO 1130 IRREPE=1,NIRREP
C
         IRREPB=DIRPRD(IRREPE,IRREPL)
         IRREPA=DIRPRD(IRREPX,IRREPE) 
         NUME=VRT(IRREPE,1)
         NUMA=VRT(IRREPA,1)
         NUMB=VRT(IRREPB,2)
C
         NROW=NUMDSL*NUMB
         NCOL=NUMA
         NSUM=NUME
C
         IOFFL=I000+IINTFP*NUMDSL*
     &              (ISYMOFF(IRREPE,IRREPL,23)-1)
         IOFFZ=I010+IINTFP*NUMDSZ*
     &              (ISYMOFF(IRREPA,IRREPR,23)-1)
         IOFFG=I00G(1)+IINTFP*
     &              (ISYMOFF(IRREPA,IRREPX,19)-1)
C
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFL),
     &              NROW,ICORE(IOFFG),NSUM,ONE,ICORE(IOFFZ),
     &              NROW)
C
1130    CONTINUE
C
        CALL SYMTR1(IRREPR,VRT(1,2),VRT(1,1),NUMDSZ,ICORE(I010),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*NUMDSZ),
     &              ICORE(ITMP+2*IINTFP*NUMDSZ))
        CALL TRANSP(ICORE(I010),ICORE(I000),DISSYZ,NUMDSZ)
C
        CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPR,LISTZ)
        CALL SAXPY(NUMDSZ*DISSYZ,ONE,ICORE(I000),1,ICORE(I010),1)
        CALL PUTLST(ICORE(I010),1,NUMDSZ,1,IRREPR,LISTZ)
C
1100   CONTINUE
C
      ENDIF
C
C ALL DONE, RETURN
C
      RETURN
      END  