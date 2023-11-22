      SUBROUTINE DHBIAJKA(ICORE,MAXCOR,IUHF,IRREPT,IRREPW,
     &                    LISTT0,LISTW0,LISTZ0)
C
C CALCULATES AAAA SPIN CASE OF "THE CONTRACTION FROM HELL" 
C
C    Z(IA,JK) =  T(AE,KM)*W(IM,JE) + T(Ae,Km)*W(Im,Je)
C
C FOR EITHER DIFFERENTIATED OR UNDIFFERENTIATED INTEGRALS
C
C CONTRACTIONS:
C
CEND
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR),I0T(2)
      DIMENSION IW(8)
      LOGICAL bRedundant
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
      IRREPZ=DIRPRD(IRREPW,IRREPT)

      bRedundant = IFLAGS2(155).EQ.0
C
      DO 10 ISPIN=1,1+IUHF
C
C CALCULATE MAXIMUM SIZE FOR W ARRAY
C
       LISTAAAA=LISTW0+ISPIN
       LISTABAB=LISTW0+5-ISPIN
       LISTZ   =LISTZ0+ISPIN
       AAAASIZW=IDSYMSZ(IRREPW,20+ISPIN,ISYTYP(2,LISTAAAA))
       ABABSIZW=IDSYMSZ(IRREPW,ISYTYP(1,LISTABAB),ISYTYP(2,LISTABAB))
       AAAASIZZ=IDSYMSZ(IRREPZ,20+ISPIN,ISYTYP(2,LISTZ))
       ISIZEZ  =IDSYMSZ(IRREPZ,ISYTYP(1,LISTZ),ISYTYP(2,LISTZ))
       MAXSIZE =MAX(AAAASIZW,AAAASIZZ,ABABSIZW,ISIZEZ)
       IFUL1=1
       IFUL2=IFUL1+IINTFP*MAXSIZE
       IFUL3=IFUL2+IINTFP*MAXSIZE
       IFUL4=IFUL3+IINTFP*MAXSIZE
       I000=IFUL3
C
C DO AAAA CONTRACTION FIRST.  READ IN ENTIRE LIST OF INTEGRALS,
C SYMMETRY EXPAND AND RESORT THEM
C
       IOFF=IFUL2
       DO 100 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPW,IRREPR)
        NUMDIS=IRPDPD(IRREPR,ISYTYP(2,LISTAAAA))
        DISSIZ=IRPDPD(IRREPL,20+ISPIN)
        DISPCK=IRPDPD(IRREPL,ISYTYP(1,LISTAAAA))
        CALL GETLST(ICORE(IOFF),1,NUMDIS,1,IRREPR,LISTAAAA)
        CALL SYMEXP2(IRREPL,POP(1,ISPIN),DISSIZ,DISPCK,NUMDIS,
     &               ICORE(IOFF),ICORE(IOFF))
        IW(IRREPR)=IOFF-IINTFP*MAXSIZE
        IOFF=IOFF+IINTFP*NUMDIS*DISSIZ 
100    CONTINUE
C
C W(IM,JE) => W(IJ,EM)
C
       CALL SSTGEN(ICORE(IFUL2),ICORE(IFUL1),AAAASIZW,POP(1,ISPIN),
     &             POP(1,ISPIN),POP(1,ISPIN),VRT(1,ISPIN),
     &             ICORE(IFUL3),IRREPW,'1342')
C
C NOW LOOP OVER IRREPS AND FORM PRODUCT Z(IJ,AK) = W(IJ,EM)*T(EM,AK)
C
       IOFFZ=IFUL2
       LISTT=LISTT0+ISPIN

       DO 200 IRREPZR=1,NIRREP
        IRREPZL=DIRPRD(IRREPZR,IRREPZ)
        IRREPTR=IRREPZR
        IRREPTL=DIRPRD(IRREPTR,IRREPT)
        IRREPWL=IRREPZL
        IRREPWR=IRREPTL
        IOFFW=IW(IRREPWR)
        DISSYT=IRPDPD(IRREPTL,ISYTYP(1,LISTT))
        NUMDST=IRPDPD(IRREPTR,ISYTYP(2,LISTT))
        DISSYW=IRPDPD(IRREPWL,20+ISPIN)
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
c        call checksum("DHBIAJKA 1 ",icore(i000),numdst*dissyt)

        CALL XGEMM ('N','N',DISSYZ,NUMDSZ,NUMDSW,ONE,ICORE(IOFFW),
     &              DISSYW,ICORE(I000),DISSYT,ZILCH,ICORE(IOFFZ),
     &              DISSYZ)
        IOFFZ=IOFFZ+IINTFP*DISSYZ*NUMDSZ
200    CONTINUE
C
C NOW READ   W(Im,Je) => W(IJ,em) [ISPIN=1]
C            W(Mi,Ej) => W(im,EM) [ISPIN=2]
C
       CALL GETALL(ICORE(IFUL3),ABABSIZW,IRREPW,LISTABAB)
       IF(ISPIN.EQ.1)THEN
        CALL SSTGEN(ICORE(IFUL3),ICORE(IFUL1),ABABSIZW,POP(1,1),
     &              POP(1,2),POP(1,1),VRT(1,2),
     &              ICORE(IFUL4),IRREPW,'1342')
       ELSE
        CALL SSTGEN(ICORE(IFUL3),ICORE(IFUL1),ABABSIZW,POP(1,1),
     &              POP(1,2),VRT(1,1),POP(1,2),
     &              ICORE(IFUL4),IRREPW,'2431')
       ENDIF
C
C CALCULATE W OFFSETS
C
       IOFF=IFUL1
       DO 201 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPW)
        DIS=IRPDPD(IRREPL,20+ISPIN)
        NUM=IRPDPD(IRREPR,11-ISPIN)
        IW(IRREPR)=IOFF
        IOFF=IOFF+IINTFP*DIS*NUM
201    CONTINUE
C
C NOW LOOP OVER IRREPS AND FORM PRODUCT
C
C    Z(IJ,AK) <= W(IJ,em)*T(em,AK)
C    Z(ij,ak) <= W(ij,EM)*T(EM,ak)
C
       IOFFZ=IFUL2
       LISTT=LISTT0+2+ISPIN
       DO 300 IRREPZR=1,NIRREP
        IRREPZL=DIRPRD(IRREPZR,IRREPZ)
        IRREPTR=IRREPZR
        IRREPTL=DIRPRD(IRREPTR,IRREPT)
        IRREPWL=IRREPZL
        IRREPWR=IRREPTL
        IOFFW=IW(IRREPWR)
        DISSYT=IRPDPD(IRREPTL,ISYTYP(1,LISTT))
        NUMDST=IRPDPD(IRREPTR,ISYTYP(2,LISTT))
        DISSYW=IRPDPD(IRREPWL,20+ISPIN)
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
c       call checksum("DHBIAJKA T ",icore(i000),numdst*dissyt)
c       call checksum("DHBIAJKA W ",icore(ioffw),DISSYW*NUMDSW)

        CALL XGEMM ('N','N',DISSYZ,NUMDSZ,NUMDSW,ONEM,ICORE(IOFFW),
     &              DISSYW,ICORE(I000),DISSYT,ONE,ICORE(IOFFZ),
     &              DISSYZ)

        IOFFZ=IOFFZ+IINTFP*DISSYZ*NUMDSZ
300    CONTINUE
C
C NOW PUT TARGET BACK TO Z(JK,IA) ORDERING AND WRITE TO DISK
C
       CALL SSTGEN(ICORE(IFUL2),ICORE(IFUL1),AAAASIZZ,POP(1,ISPIN),
     &             POP(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),
     &             ICORE(IFUL3),IRREPZ,'2413')
C
C ANTISYMMETRIZE ON LHS INDICES
C
       IOFF=IFUL1
       DO 400 IRREPZR=1,NIRREP
        IRREPZL=DIRPRD(IRREPZ,IRREPZR)
        DISSYZX=IRPDPD(IRREPZL,20+ISPIN)
        DISSYZ =IRPDPD(IRREPZL,ISYTYP(1,LISTZ))
        NUMDSZ =IRPDPD(IRREPZR,ISYTYP(2,LISTZ))
        CALL TRANSP(ICORE(IOFF),ICORE(IFUL2),NUMDSZ,DISSYZX)
        CALL ASSYM2(IRREPZL,POP(1,ISPIN),NUMDSZ,ICORE(IFUL2))
        CALL TRANSP(ICORE(IFUL2),ICORE(IOFF),DISSYZ,NUMDSZ)
        CALL GETLST(ICORE(IFUL2),1,NUMDSZ,1,IRREPZR,LISTZ)
        CALL SAXPY (NUMDSZ*DISSYZ,ONEM,ICORE(IOFF),1,ICORE(IFUL2),1)
        CALL PUTLST(ICORE(IFUL2),1,NUMDSZ,1,IRREPZR,LISTZ)

        IOFF=IOFF+IINTFP*NUMDSZ*DISSYZX
400    CONTINUE
C
10    CONTINUE
C
      RETURN
      END
