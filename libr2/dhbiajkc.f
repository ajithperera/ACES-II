      SUBROUTINE DHBIAJKC(ICORE,MAXCOR,IUHF,IRREPT,IRREPW,
     &                    LISTT0,LISTW0,LISTZ0)
C
C CALCULATES ABAB SPIN CASE OF "THE CONTRACTION FROM HELL" FOR
C UHF CASES:
C
C  Z(iA,jK) = T(Ae,Km)*W(im,je) + T(AE,KM)*W(iM,jE) 
C            -T(eA,jM)*W(iM,eK)
C
C FOR EITHER DIFFERENTIATED OR UNDIFFERENTIATED INTEGRALS
C
C CONTRACTIONS:
C
CEND
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR),I0T(2)
      DIMENSION IW(8)
      logical bRedundant
      CHARACTER*4 SPCASE
      DOUBLE PRECISION ONE,ONEM,ZILCH,HALF,TWO,THIRD
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /FLAGS2/IFLAGS2(500)
C
      DATA ONE,ONEM,ZILCH /1.0D0,-1.0D0,0.0D0/
C
      bRedundant = iflags2(155).eq.0
      IRREPZ=DIRPRD(IRREPW,IRREPT)
C
C CALCULATE MAXIMUM SIZE FOR W ARRAY
C
      LISTBBBB=LISTW0+2
      LISTABAB=LISTW0+4
      LISTBABA=LISTW0+3
      LISTZ   =LISTZ0+3
      BBBBSIZW=IDSYMSZ(IRREPW,22,ISYTYP(2,LISTBBBB))
      ABABSIZW=IDSYMSZ(IRREPW,ISYTYP(1,LISTABAB),ISYTYP(2,LISTABAB))
      BABASIZW=IDSYMSZ(IRREPW,ISYTYP(1,LISTBABA),ISYTYP(2,LISTBABA))
      ISIZEZ  =IDSYMSZ(IRREPZ,ISYTYP(1,LISTZ),ISYTYP(2,LISTZ))
      MAXSIZE =MAX(BBBBSIZW,ABABSIZW,BABASIZW,ISIZEZ)
      IFUL1=1
      IFUL2=IFUL1+IINTFP*MAXSIZE
      IFUL3=IFUL2+IINTFP*MAXSIZE
      IFUL4=IFUL3+IINTFP*MAXSIZE
      I000=IFUL3
C
C DO BBBB CONTRACTION FIRST.  READ IN ENTIRE LIST OF INTEGRALS,
C SYMMETRY EXPAND AND RESORT THEM
C
      IOFF=IFUL2
      DO 100 IRREPR=1,NIRREP
       IRREPL=DIRPRD(IRREPW,IRREPR)
       NUMDIS=IRPDPD(IRREPR,ISYTYP(2,LISTBBBB))
       DISSIZ=IRPDPD(IRREPL,22)
       DISPCK=IRPDPD(IRREPL,ISYTYP(1,LISTBBBB))
       CALL GETLST(ICORE(IOFF),1,NUMDIS,1,IRREPR,LISTBBBB)
       CALL SYMEXP2(IRREPL,POP(1,2),DISSIZ,DISPCK,NUMDIS,
     &              ICORE(IOFF),ICORE(IOFF))
       IW(IRREPR)=IOFF-IINTFP*MAXSIZE
       IOFF=IOFF+IINTFP*NUMDIS*DISSIZ 
100   CONTINUE
C
C W(im,je) => W(ij,em)
C
      CALL SSTGEN(ICORE(IFUL2),ICORE(IFUL1),BBBBSIZW,POP(1,2),
     &            POP(1,2),POP(1,2),VRT(1,2),
     &            ICORE(IFUL3),IRREPW,'1342')
C
C NOW LOOP OVER IRREPS AND FORM PRODUCT Z(ij,AK) = W(ij,em)*T(em,AK)
C
      IOFFZ=IFUL2
      LISTT=LISTT0+3
      DO 200 IRREPZR=1,NIRREP
       IRREPZL=DIRPRD(IRREPZR,IRREPZ)
       IRREPTR=IRREPZR
       IRREPTL=DIRPRD(IRREPTR,IRREPT)
       IRREPWL=IRREPZL
       IRREPWR=IRREPTL
       IOFFW=IW(IRREPWR)
       DISSYT=IRPDPD(IRREPTL,ISYTYP(1,LISTT))
       NUMDST=IRPDPD(IRREPTR,ISYTYP(2,LISTT))
       DISSYW=IRPDPD(IRREPWL,22)
       NUMDSW=DISSYT
       DISSYZ=DISSYW
       NUMDSZ=NUMDST
       IF(bRedundant) THEN
          CALL GETLST(ICORE(I000),1,NUMDST,1,IRREPTR,LISTT)
       ELSE
           I010=I000+NUMDST*DISSYT*IINTFP
           CALL GETLST_NR(ICORE(I000),ICORE(I010),MAXCOR-I010,
     &                    LISTT,IRREPTR)
       ENDIF
       CALL XGEMM ('N','N',DISSYZ,NUMDSZ,NUMDSW,ONEM,ICORE(IOFFW),
     &             DISSYW,ICORE(I000),DISSYT,ZILCH,ICORE(IOFFZ),
     &             DISSYZ)
       IOFFZ=IOFFZ+IINTFP*DISSYZ*NUMDSZ
200   CONTINUE
C
C NOW READ   W(Mi,Ej) => W(ij,EM)
C
      CALL GETALL(ICORE(IFUL3),BABASIZW,IRREPW,LISTBABA)
      CALL SSTGEN(ICORE(IFUL3),ICORE(IFUL1),BABASIZW,POP(1,1),
     &            POP(1,2),VRT(1,1),POP(1,2),
     &            ICORE(IFUL4),IRREPW,'2431')
C
C CALCULATE W OFFSETS
C
      IOFF=IFUL1
      DO 201 IRREPR=1,NIRREP
       IRREPL=DIRPRD(IRREPR,IRREPW)
       DIS=IRPDPD(IRREPL,22)
       NUM=IRPDPD(IRREPR,9)
       IW(IRREPR)=IOFF
       IOFF=IOFF+IINTFP*DIS*NUM
201   CONTINUE
C
C NOW LOOP OVER IRREPS AND FORM PRODUCT
C
C    Z(ij,AK) <= W(ij,EM)*T(EM,AK)
C
      IOFFZ=IFUL2
      LISTT=LISTT0+1
      DO 300 IRREPZR=1,NIRREP
       IRREPZL=DIRPRD(IRREPZR,IRREPZ)
       IRREPTR=IRREPZR
       IRREPTL=DIRPRD(IRREPTR,IRREPT)
       IRREPWL=IRREPZL
       IRREPWR=IRREPTL
       IOFFW=IW(IRREPWR)
       DISSYT=IRPDPD(IRREPTL,ISYTYP(1,LISTT))
       NUMDST=IRPDPD(IRREPTR,ISYTYP(2,LISTT))
       DISSYW=IRPDPD(IRREPWL,22)
       NUMDSW=DISSYT
       DISSYZ=DISSYW
       NUMDSZ=NUMDST
       IF(bRedundant) THEN
          CALL GETLST(ICORE(I000),1,NUMDST,1,IRREPTR,LISTT)
       ELSE
           I010=I000+NUMDST*DISSYT*IINTFP
           CALL GETLST_NR(ICORE(I000),ICORE(I010),MAXCOR-I010,
     &                    LISTT,IRREPTR)
       ENDIF
       CALL XGEMM ('N','N',DISSYZ,NUMDSZ,NUMDSW,ONE,ICORE(IOFFW),
     &             DISSYW,ICORE(I000),DISSYT,ONE,ICORE(IOFFZ),
     &             DISSYZ)
       IOFFZ=IOFFZ+IINTFP*DISSYZ*NUMDSZ
300   CONTINUE
C
C NOW WE NEED TO REORDER THE TARGET
C
C   Z(ij,AK) => Z(iK,Aj)
C
      CALL SSTGEN(ICORE(IFUL2),ICORE(IFUL1),ISIZEZ,POP(1,2),
     &            POP(1,2),VRT(1,1),POP(1,1),
     &            ICORE(IFUL4),IRREPZ,'1432')
      CALL SCOPY (ISIZEZ,ICORE(IFUL1),1,ICORE(IFUL2),1)
C
C NOW THE EXCHANGE CONTRIBUTION
C
C NOW READ   W(Mi,Ke) => W(ie,KM) => W(iK,eM)
C
      CALL GETALL(ICORE(IFUL1),ABABSIZW,IRREPW,LISTABAB)
      CALL SSTGEN(ICORE(IFUL1),ICORE(IFUL3),ABABSIZW,POP(1,1),
     &            POP(1,2),POP(1,1),VRT(1,2),
     &            ICORE(IFUL4),IRREPW,'2431')
      CALL SSTGEN(ICORE(IFUL3),ICORE(IFUL1),ABABSIZW,POP(1,2),
     &            VRT(1,2),POP(1,1),POP(1,1),
     &            ICORE(IFUL4),IRREPW,'1324')
C
C CALCULATE W OFFSETS
C
      IOFF=IFUL1
      DO 301 IRREPR=1,NIRREP
       IRREPL=DIRPRD(IRREPR,IRREPW)
       DIS=IRPDPD(IRREPL,14)
       NUM=IRPDPD(IRREPR,12)
       IW(IRREPR)=IOFF
       IOFF=IOFF+IINTFP*DIS*NUM
301   CONTINUE
C
C NOW LOOP OVER IRREPS AND FORM PRODUCT
C
C    Z(iK,Aj) <= W(iK,eM)*T(eM,Aj)
C
      IOFFZ=IFUL2
      LISTT=LISTT0+5
      DO 400 IRREPZR=1,NIRREP
       IRREPZL=DIRPRD(IRREPZR,IRREPZ)
       IRREPTR=IRREPZR
       IRREPTL=DIRPRD(IRREPTR,IRREPT)
       IRREPWL=IRREPZL
       IRREPWR=IRREPTL
       IOFFW=IW(IRREPWR)
       DISSYT=IRPDPD(IRREPTL,ISYTYP(1,LISTT))
       NUMDST=IRPDPD(IRREPTR,ISYTYP(2,LISTT))
       DISSYW=IRPDPD(IRREPWL,14)
       NUMDSW=DISSYT
       DISSYZ=DISSYW
       NUMDSZ=NUMDST
       IF(bRedundant) THEN
          CALL GETLST(ICORE(I000),1,NUMDST,1,IRREPTR,LISTT)
       ELSE
           I010=I000+NUMDST*DISSYT*IINTFP
           CALL GETLST_NR(ICORE(I000),ICORE(I010),MAXCOR-I010,
     &                    LISTT,IRREPTR)
       ENDIF
       CALL XGEMM ('N','N',DISSYZ,NUMDSZ,NUMDSW,ONE,ICORE(IOFFW),
     &             DISSYW,ICORE(I000),DISSYT,ONE,ICORE(IOFFZ),
     &             DISSYZ)
       IOFFZ=IOFFZ+IINTFP*DISSYZ*NUMDSZ
400   CONTINUE
C
C NOW WE NEED TO REORDER THE TARGET
C
C   Z(iK,Aj) => Z(Kj,Ai)
C
      CALL SSTGEN(ICORE(IFUL2),ICORE(IFUL1),ISIZEZ,POP(1,2),
     &            POP(1,1),VRT(1,1),POP(1,2),
     &            ICORE(IFUL4),IRREPZ,'2431')
C
C ALL DONE
C
      CALL GETALL(ICORE(IFUL2),ISIZEZ,IRREPZ,LISTZ)
      CALL SAXPY (ISIZEZ,ONEM,ICORE(IFUL1),1,ICORE(IFUL2),1)
      CALL PUTALL(ICORE(IFUL2),ISIZEZ,IRREPZ,LISTZ)
      RETURN
      END
