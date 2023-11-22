      SUBROUTINE S2HHPPT2S2(ICORE,MAXCOR,IUHF, ISIDE, ISPIN)
C
C   THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF S2 TO
C   S2 INVOLVING THE THREE-PARTICLE EFFECTIVE INTEGRALS W.
C   THESE CONTRIBUTIONS ARE EVALUATED IN TWO STEPS.
C
C CONTRACTIONS :
C
C Q(F,P) =  SUM <NM||EF> *S2(NM,PE) =  <NM|EF> * S2(NM,EP)
C           NME
C Z(IJ,AP) = -SUM [Q(E,P) T2(EA,IJ) ]
C                  E
C
C   IT IS SEEN THAT IRREPNMF=IRREPE=SIRREP; NF has ISPIN and M has MIXSPIN.
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
C
      DIMENSION ICORE(MAXCOR)
      DIMENSION LISTW0(2,2),NI(8),NJ(8),NN(8),NM(8)
      LOGICAL PRINT
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SINFO/ NS(8), SIRREP
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
C
      DATA ONE  /1.0/
      DATA ZILCH/0.0/
      DATA ONEM/-1.0/

      Write(6,"(a)") "-----------------S2HHPPT2S2-------------------"
C
C 16 : V(Ab,Ij) =  W(Ij, Ab) ordered as Ab;Ij
C 46 : T(Ij; Ab) ordered as Ab;Ij
C
      LISTW0(1,1) = 14
      LISTW0(1,2) = 16
      LISTW0(2,1) = 16
      LISTW0(2,2) = 15
C
C  THE SECOND SPINLABEL OF LISTW0 EQUALS ISPIN, THE FIRST LABEL DETERMINES
C  THE SPIN-TYPE OF THE S-VECTOR (IMIXSPIN)
C  LISTT0 = 30 + LISTW0
C  LIST 16 (46) CONTAINS INTEGRALS V(Ab, Ij) and T2(Ab,Ij), WHICH MEANS THAT THE ORDERING IS
C  CONVENIENT FOR ISPIN = 2, BUT FOR ISPIN = 1 BOTH THE ROW AND COLUMNS HAVE TO
C  BE SUBJECTED TO SYMTR.
C
       LENTAR= NS(SIRREP) * VRT(SIRREP,ISPIN)
       I000=1
       I010=I000+LENTAR*IINTFP
       CALL ZERO(ICORE(I000),LENTAR)

       DO 50 IMIXSPIN = 1, 1+IUHF
          LISTW = LISTW0(IMIXSPIN, ISPIN+1-IUHF) + (ISIDE - 1) * 30
          LISTS2IN = LS2IN(ISPIN,IMIXSPIN+1-IUHF)
          MSPIN = IMIXSPIN
          IF (((ISPIN .EQ. 2) .AND.
     $       (MSPIN .EQ.1)) .OR. (IUHF.EQ.0)) THEN
C
C RHF: S2(AB,A) W(ABAB,16) (IUHF=0)
C UHF: S2(AA,A) W(AAAA,14) (ISPIN=2)
C
             CASE = 1

          ELSEIF (ISPIN.EQ.MSPIN) THEN
C
C UHF: S2(AAA) W(AAAA,14), S2(BBB) W(BBBB)
C
             CASE = 2
          ELSE
C
C UHF: S2(AB,A) W(ABAB,16)
C
             CASE  = 3
          ENDIF
C
C LOOP OVER IRREPS OF RHS OF S VECTOR S(NM,EP) and W(NM,EF)
C
       DO 100 XIRREP=1,NIRREP
          IRREPE=DIRPRD(XIRREP,SIRREP)
          IRREPW=XIRREP 
          IRREPN=IRREPE
          IRREPM=SIRREP
          IRREPF=IRREPM
          DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW)) 
          NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW))
          IF (CASE .EQ. 2) THEN
              DISSYW = IRPDPD(XIRREP,18+ISPIN)
          ENDIF
          DISSYS=IRPDPD(XIRREP,ISYTYP(1,LISTS2IN))
          NUMDSS=VRT(IRREPE,MSPIN)*NS(SIRREP)
          MAXW=MAX(NUMDSW,DISSYW)
          I020=I010+IINTFP*DISSYW*NUMDSW
          I030=I020+IINTFP*MAX(DISSYS*NUMDSS,3*MAXW)

          IF(I030.GT.MAXCOR)THEN
C
C OUT-OF-CORE ALGORITHM
C   
             CALL INSMEM("S2IHHPPT2S2",I030,MAXCOR)
          ENDIF 

C DO IN-CORE ALGORITHM
C
C
C TREAT IMIXSPIN = 2, ISPIN =1 AS SPECIAL CASE (CASE=3)
C IMIXSPIN = ISPIN, UHF  AS SPECIAL CASE (CASE=2)
C
       IF (CASE.EQ.1) THEN
C
C READ W and Spin-adapt for RHF.
C     
C W(NM,EF) = W(EF,NM) =  2 W(Ef,NM) - W(Fe,Nm)   [ RHF ] W is 
C ORDERED AS FE,NM. 

           CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPW,LISTW)
           IF (IUHF.EQ.0) THEN
              ITMP1=I020
              ITMP2=ITMP1+IINTFP*MAXW
              IEND=ITMP2+IINTFP*MAXW
              IF (IEND .GT. MAXCOR) CALL INSMEM("@-S2IHHPPT2S2",
     &                                        IEND,MAXCOR)
              CALL SPINAD3(IRREPW,VRT(1,1),DISSYW,NUMDSW,
     &                     ICORE(I010),ICORE(ITMP1),ICORE(ITMP2))
           ENDIF
#ifdef _DEBUG_LVL0
c      write(6,"(a,4i4)") "case,listw,numdsw,dissyw:",case,listw,
c     &                    numdsw,dissyw
      call checksum("->w22   :",icore(i010),numdsw*dissyw)
#endif 
C
       ELSE IF (CASE.EQ.2) THEN

            DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW)) 
            CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPW,LISTW)
            CALL DCOPY(NUMDSW*DISSYW,ICORE(I010),1,ICORE(I020),1)
            CALL TRANSP(ICORE(I020),ICORE(I010),NUMDSW,DISSYW)
            CALL SYMEXP(IRREPW,VRT(1,ISPIN),NUMDSW,ICORE(I010))
            DISSYW = IRPDPD(IRREPW,18+ISPIN)

#ifdef _DEBUG_LVL0
      write(6,"(a,4i4)") "case,listw,numdsw,dissyw:",case,listw,
     &                    numdsw,dissyw
      call checksum("->w22   :",icore(i010),numdsw*dissyw)
#endif 
        ELSE         
C
C TREAT SPECIAL CASE ISPIN =2, IMIXSPIN = 1, UHF CASE
C W(Nm|Ef) is needed. Read in the W(Ef,Nm) and permute both
C bra and ket indices to obtain W(fE,mN)=W(mN,fE)
C UHF: S2(AA,B) W(ABAB,21)

           CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPW,LISTW)
           IF (ISPIN .EQ. 1 .AND. MSPIN .EQ. 1) THEN
C
C NOW PERMUTE W(Ef|Mn) INTO W(fE|nM), ROW-INDICES FIRST
C
           NUMW = IRPDPD(IRREPW,14)
           NUMW = POP(IRREPM,ISPIN)*POP(IRREPN,MSPIN)
           MAXW = MAX(NUMW, DISSYW)
           I030 = I020 + MAXW * IINTFP
           I040 = I030 + MAXW * IINTFP
           IEND = I040 + MAXW * IINTFP
           IF (IEND .GT. MAXCOR) CALL INSMEM("S2IHHPPT2S2",I030,MAXCOR)
           CALL SYMTR3(IRREPW,VRT(1,MSPIN),VRT(1,ISPIN),DISSYW,NUMW,
     $                 ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040))
C
C TO PERMUTE COLUMN-INDICES WE NEED POPULATION VECTORS (fE,Nm)->(fE,mN)
C
         CALL IZERO(NN,8)
         CALL IZERO(NM,8)
         NN(IRREPN) = POP(IRREPN,MSPIN)
         NM(IRREPM) = POP(IRREPM,ISPIN)
           CALL SYMTR1(IRREPW,NN,NM,DISSYW,ICORE(I010),ICORE(I020),
     $                 ICORE(I030),ICORE(I040))
         ENDIF 
#ifdef _DEBUG_LVL0
      write(6,"(a,4i4)") "case,listw,numdsw,dissyw:",case,listw,
     &                    numdsw,dissyw
      call checksum("->w22   :",icore(i010),numdsw*dissyw)
#endif 
        ENDIF

        I030 = I020 + NUMDSW*DISSYW*IINTFP
        IEND = I030 + NUMDSW*DISSYW*IINTFP
        IF (IEND .GT. MAXCOR) CALL INSMEM("@-S2IHHPPT2S2",
     &                                        IEND,MAXCOR)
        IF (CASE .NE. 2) THEN
           CALL DCOPY(NUMDSW*DISSYW,ICORE(I010),1,ICORE(I020),1)
           CALL TRANSP(ICORE(I020),ICORE(I010),NUMDSW,DISSYW)
        ENDIF 
C
C  S2in(NM,EP)
C
        CALL GETLST(ICORE(I020),1,NUMDSS,1,XIRREP,LISTS2IN)

#ifdef _DEBUG_LVL0
      write(6,"(a,5(1x,i4))") "@-W:nfirst,nblock,numdsw,dissyw,numdss:",
     +                    nfirst,nblock,numdsw,dissyw,numdss
      write(6,"(a,2(1x,i4))") "@-W:dissys,numdss:",dissys,numdss
c      call output(icore(i010),1,numdsw,1,dissyw,numdsw,dissyw,1)
      call checksum("->s2    :",icore(i020),numdss*dissys)
c      call output(icore(i020),1,dissys,1,numdss,dissys,numdss,1)
#endif
        CALL CONTRACT_w2S2(ICORE(I010),ICORE(I020),ICORE(I000),NUMDSW,
     &                     DISSYW,DISSYS,NUMDSS,LENTAR,XIRREP,ISPIN,
     &                     MSPIN)
#ifdef _DEBUG_LVL0
      call checksum("q(f,p)I :",icore(i000),lentar)
#endif

 100  CONTINUE
 50   CONTINUE

#ifdef _DEBUG_LVL0
      call checksum("q(f,p)F :",icore(i000),lentar)
C      WRITE(6,*)' Q-INTERMEDIATES IN S2CON3S2'
C      CALL OUTPUT(ICORE(I000), 1, NROW, 1, NCOL, NROW, NCOL, 1)
#endif
C
C NOW CONTRACT Q WITH T
C
      DO 75 IMIXSPIN = 1, 1+IUHF
         LISTT    = LISTW0(IMIXSPIN,ISPIN+1-IUHF) + (2-ISIDE) * 30
         LISTS2EX = LS2OUT(ISPIN,IMIXSPIN+1-IUHF)
         MSPIN    = IMIXSPIN 

         IF (((ISPIN .EQ. 2) .AND.
     $        (MSPIN .EQ.1)) .OR. (IUHF.EQ.0)) THEN
C
C RHF: S2(AA,A) W(BABA,12) (IUHF=0)
C UHF: S2(BB,A) W(BABA,12) (ISPIN=2)
C
              CASE = 1

         ELSEIF (ISPIN.EQ.MSPIN) THEN
C
C UHF: S2(AAA) W(AAAA), S2(BBB) W(BBBB)
C
              CASE = 2
         ELSE
C
C UHF: S2(AA,B) W(ABAB,21)
C
              CASE  = 3
         ENDIF
C
C Z(IJ,AP) = - SUM [T2(AF,IJ)*Q(F,P) ] Note that IRREPF=SIRREP (MSPIN)
C
         DO 200 XIRREP = 1, NIRREP
            IRREPA=DIRPRD(XIRREP,SIRREP)
            IRREPT=XIRREP
            IRREPI=IRREPA
            IRREPJ=SIRREP
            IRREPF=IRREPJ
            DISSYT=IRPDPD(XIRREP,ISYTYP(1,LISTT)) 
            NUMDST=IRPDPD(XIRREP,ISYTYP(2,LISTT)) 
            IF (CASE.EQ.2) THEN
               DISSYT = IRPDPD(IRREPT,18+ISPIN)
            ENDIF
            MAXT = MAX(DISSYT,NUMDST)
            DISSYS=IRPDPD(XIRREP,ISYTYP(1,LISTS2EX))
            NUMDSS=VRT(IRREPA,MSPIN)*NS(SIRREP)

            IF (NUMDST * NUMDSS .GT. 0) THEN
               I020=I010+IINTFP*DISSYT*NUMDST          
               IEND=I020+MAX(DISSYS*NUMDSS, 3*MAXT)*IINTFP
               IF (IEND .GT. MAXCOR) CALL INSMEM("@-S2IHHPPT2S2",
     &                                            IEND,MAXCOR)
C 
            IF (CASE.EQ.1) THEN

            CALL GETLST(ICORE(I010),1,NUMDST,1,IRREPT,LISTT)

            ELSE IF (CASE.EQ.2) THEN

            DISSYT = IRPDPD(IRREPT,ISYTYP(1,LISTT))
            CALL GETLST(ICORE(I010),1,NUMDST,1,IRREPT,LISTT)
            CALL DCOPY(NUMDST*DISSYT,ICORE(I010),1,ICORE(I020),1)
            CALL TRANSP(ICORE(I020),ICORE(I010),NUMDST,DISSYT)
            CALL SYMEXP(IRREPT,VRT(1,ISPIN),NUMDST,ICORE(I010))
            DISSYT = IRPDPD(IRREPT,18+ISPIN)

            ELSE         

C TREAT SPECIAL CASE ISPIN =2, MSPIN = 1, UHF CASE
C Read T2 as T2(Af,Ji) and transpose to both bra and Ket 
C to iobtain T2(fA|iJ).
C
C ONLY THE KIRREP=SIRREP, AND MIRREP INTEGRALS ARE NEEDED.
C THEREFORE CALULATE OFFSET
C
         CALL GETLST(ICORE(I010),1,NUMDST,1,IRREPT,LISTT)
C
C NOW PERMUTE T(Af|Ji) INTO T(fA|Ji), ROW-INDICES FIRST
C
         IF (ISPIN .EQ. 1 .AND. MSPIN .EQ. 1) THEN
         NUMT = IRPDPD(XIRREP,ISYTYP(2,LISTT))
         NUMT = POP(IRREPI,MSPIN)*POP(IRREPJ,ISPIN)
         MAXW = MAX(NUMT,DISSYT)
         I030 = I020 + MAXW * IINTFP
         I040 = I030 + MAXW * IINTFP
         IEND = I040 + MAXW * IINTFP
         IF (IEND .GT. MAXCOR) CALL INSMEM("@-S2IHHPPT2S2",
     &                                      IEND,MAXCOR)
         CALL SYMTR3(IRREPT,VRT(1,ISPIN),VRT(1,MIXPIN),
     $               DISSYT,NUMT,ICORE(I010),ICORE(I020),
     $               ICORE(I030),ICORE(I040))
C
C TO PERMUTE COLUMN-INDICES WE NEED POPULATION VECTORS 
C T2(fA|Ji) -> T2(fA,iJ)
C
         CALL IZERO(NI,8)
         CALL IZERO(NJ,8)
         NI(IRREPI) = POP(IRREPI,MSPIN)
         NJ(IRREPJ) = POP(IRREPJ,ISPIN)
         CALL SYMTR1(IRREPT,NI,NJ,DISSYT,
     $               ICORE(I010),ICORE(I020),ICORE(I030),
     $               ICORE(I040))
         ENDIF 
         ENDIF 
C            
         I030 = I020 + NUMDST*DISSYT*IINTFP
         IEND = I030 + NUMDST*DISSYT*IINTFP
         IF (IEND .GT. MAXCOR) CALL INSMEM("@-S2IHHPPT2S2",
     &                                     IEND,MAXCOR)
         IF (CASE .NE. 2) THEN
            CALL DCOPY(NUMDST*DISSYT,ICORE(I010),1,ICORE(I020),1)
            CALL TRANSP(ICORE(I020),ICORE(I010),NUMDST,DISSYT)
         ENDIF 

         CALL GETLST(ICORE(I020),1,NUMDSS,1,XIRREP,LISTS2EX)

#ifdef _DEBUG_LVL0
      call checksum("->hs2   :",icore(i020),numdss*dissys)
#endif 
C
C Z(IJ,AP) = - T2(IJA,F) Q(F,P) 
C
        CALL CONTRACT_T2Q(ICORE(I010),ICORE(I000),ICORE(I020),NUMDST,
     &                    DISSYT,DISSYS,NUMDSS,LENTAR,XIRREP,ISPIN,
     &                    MSPIN)
         
         CALL PUTLST(ICORE(I020),1,NUMDSS,1,XIRREP,LISTS2EX)

#ifdef _DEBUG_LVL0
      call checksum("s2-fs2  :",icore(i020),numdss*dissys)
#endif 
       ENDIF
 200   CONTINUE
 75   CONTINUE
C
      RETURN
      END
