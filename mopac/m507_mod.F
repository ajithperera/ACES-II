      SUBROUTINE EXTPAR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
************************************************************************
      CHARACTER NUMBRS(0:9)*1, PARTYP(28)*5, FILES*64, DUMMY(2)*50,     IR0494 
     1          KEYWRD*80, TEXT*50, TXTNEW*50, ELEMNT(107)*2 
      COMMON /ATHEAT/ ATHEAT
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /ATOMIC/ EISOL(107),EHEAT(107)
      COMMON /KEYWRD/ KEYWRD
C
C     COMMON BLOCKS FOR SRP (Ivan Rossi - April 94)                    
C
      COMMON /BETAS / BETAS(107),BETAP(107),BETAD(107)
      COMMON /SRPI/ IBTPTR(107), NATPTR(MXATSP), NATSP
     *       /SRPL/ ISSRP
     *       /SRPR/ BETSS(MAXBET), BETSP(MXATSP,MXATSP), BETPP(MAXBET)
C
      DIMENSION  IJPARS(5,500), PARSIJ(500), IBLNK(2)
         SAVE                                                           GL0892
      LOGICAL ISSRP                                                     IR0494

      DATA NUMBRS/' ','1','2','3','4','5','6','7','8','9'/
      DATA PARTYP/'USS  ','UPP  ','UDD  ','ZS   ','ZP   ','ZD   ',
     1    'BETAS','BETAP','BETAD','GSS  ','GSP  ','GPP  ','GP2  ',
     2    'HSP  ','AM1  ','EXPC ','GAUSS','ALP  ','GSD  ','GPD  ',
     3    'GDD  ','FN1  ','FN2  ','FN3  ','ORB  ','BETPP','BETSP',
     4    'BETSS'/                                                      IR0494
      DATA (ELEMNT(I),I=1,107)/'H ','HE',
     1 'LI','BE','B ','C ','N ','O ','F ','NE',
     2 'NA','MG','AL','SI','P ','S ','CL','AR',
     3 'K ','CA','SC','TI','V ','CR','MN','FE','CO','NI','CU',
     4 'ZN','GA','GE','AS','SE','BR','KR',
     5 'RB','SR','Y ','ZR','NB','MO','TC','RU','RH','PD','AG',
     6 'CD','IN','SN','SB','TE','I ','XE',
     7 'CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB','DY',
     8 'HO','ER','TM','YB','LU','HF','TA','W ','RE','OS','IR','PT',
     9 'AU','HG','TL','PB','BI','PO','AT','RN',
     1 'FR','RA','AC','TH','PA','U ','NP','PU','AM','CM','BK','CF','XX',
     2 'FM','MD','CB','++','+','--','-','TV'/
      DATA (IBLNK(i),i=1,2) / 32, 160/                                  IR0494 
      ISSRP=.FALSE.                                                     IR0494
      I=INDEX(KEYWRD,'EXTERNAL=')+9
      J=INDEX(KEYWRD(I:),' ')+I-1
      FILES=KEYWRD(I:J)
      WRITE(6,'(//5X,'' PARAMETER TYPE      ELEMENT    PARAMETER'')')
      OPEN(14,STATUS='OLD',FILE=FILES)
      I=0
   10 READ(14,'(A40)',ERR=90,IOSTAT=IOS)TEXT
      IF(TEXT.EQ.' ')GOTO 90
      IF(INDEX(TEXT,'END').NE.0)GOTO 90
      ILOWA = ICHAR('a')
      ILOWZ = ICHAR('z')
      ICAPA = ICHAR('A')
************************************************************************
      DO 20 I=1,50
         ILINE=ICHAR(TEXT(I:I))
         IF(ILINE.GE.ILOWA.AND.ILINE.LE.ILOWZ) THEN
            TEXT(I:I)=CHAR(ILINE+ICAPA-ILOWA)
         ENDIF
   20 CONTINUE
************************************************************************
      IF(INDEX(TEXT,'END') .NE. 0) GOTO 90
      DO 30 J=1,25
         IF(J.GT.21) THEN
            IT=INDEX(TEXT,'FN')
            TXTNEW = TEXT(1:IT+2)
            IF(INDEX(TXTNEW,PARTYP(J)) .NE. 0) GOTO 40
         ENDIF
         IF(INDEX(TEXT,PARTYP(J)) .NE. 0) GOTO 40
   30 CONTINUE
C
C SPECIAL SRP BETA initialization (Ivan Rossi - April '94)
C
      DO 31 J=26,28
         IF(INDEX(TEXT,PARTYP(J)) .NE. 0) THEN 
            IF(.NOT. ISSRP) THEN 
               ISSRP=.TRUE.                          
C Check if one of the parameter specified before is a BETA and update
               DO 32 IPAR=1,LPARS
                 IF (IJPARS(3,IPAR).eq. 7 .or.IJPARS(3,IPAR).eq.8) THEN 
                    CALL UPDATE( IJPARS(3,IPAR), IJPARS(2,IPAR),
     1                           PARSIJ(IPAR), 1, IJPARS(1,IPAR) )  
                 ENDIF
  32           CONTINUE
               CALL INIBET                          
            ENDIF
            GOTO 40
         ENDIF                                  
  31  CONTINUE
C
C End SRP BETA init
C
      WRITE(6,'(''  FAULTY LINE:'',A)')TXTNEW
      WRITE(6,'(''  FAULTY LINE:'',A)')TEXT
      WRITE(6,'(''   NAME NOT FOUND'')')
      STOP
   40 IPARAM=J
      IF(IPARAM.GT.21 .AND. J.LT.26) THEN                               IR0494
         I=INDEX(TEXT,'FN')
         KFN=READA(TEXT,I+3)
      ELSE
         KFN=0
         I=INDEX(TEXT,PARTYP(J))
      ENDIF
C     K=INDEX(TEXT(I:),' ')+1
C     DUMMY=TEXT(K:)
C     TEXT=DUMMY
C Quick and dirty parsing for the first element in the string           IR0494
      i=i+3
   42 NASC=ICHAR(text(i:i))
      IF(NASC .eq. IBLNK(1) .or. NASC .eq. IBLNK(2)) goto 43
      i=i+1
      goto 42
   43 NASC=ICHAR(TEXT(I:I))
      IF(NASC .ne. IBLNK(1) .and. NASC .ne. IBLNK(2)) goto 44
      i=i+1
      goto 43
   44 DUMMY(1)=TEXT(i-1 : i+1)
      DUMMY(2)=TEXT(i:)
      TEXT=DUMMY(2)
      DO 50 J=1,107
   50 IF(INDEX(DUMMY(1), ' '//ELEMNT(J)) .NE. 0) GOTO 60
   51 WRITE(6,'('' ELEMENT NOT FOUND '')')
      WRITE(6,*)' FAULTY LINE: "'//TEXT//'"'
      STOP
   60 IELMNT=J
C
C READ second element for SRP BETA and CHECK (Ivan Rossi - April '94)
C
      IF (ISSRP) THEN
         DUMMY(2)=TEXT(3:) 
         TEXT=DUMMY(2)
         DO 61 J=1,107
  61        IF(INDEX(TEXT,' '//ELEMNT(J)) .NE. 0) goto 62
C ...if second element is wrong then error
  62     IF(J .gt. 107 ) GOTO 51
         IELMN2=j
      ENDIF
C
C END READ second element
C
      PARAM=READA(TEXT,INDEX(TEXT,ELEMNT(J)))
C
C SET special SRP BETA (Ivan Rossi - April '94)
C
      IF(ISSRP) THEN
         CALL SETBET( IELMNT, IELMN2, PARTYP(IPARAM), PARAM)
         WRITE(6,'(10x,"SRP",A6,8X,A2,1X,A2,F14.6)') 
     *            PARTYP(IPARAM),ELEMNT(IELMNT),ELEMNT(IELMN2),PARAM
      ENDIF
C
C end SET 
C
      DO 70 I=1,LPARS
         IF(IJPARS(1,I).EQ.KFN.AND.IJPARS(2,I).EQ.IELMNT.AND.
     1IJPARS(3,I).EQ.IPARAM) GOTO 80
   70 CONTINUE
      LPARS=LPARS+1
      I=LPARS
   80 IJPARS(1,I)=KFN
      IJPARS(2,I)=IELMNT
      IJPARS(3,I)=IPARAM
      PARSIJ(I)=PARAM
      GOTO 10
   90 CONTINUE
      CLOSE(14)
      DO 120 J=1,107
         DO 110 K=1,25
            DO 100 I=1,LPARS
               IPARAM=IJPARS(3,I)
               KFN=IJPARS(1,I)
               IELMNT=IJPARS(2,I)
               IF(IPARAM.NE.K) GOTO 100
               IF(IELMNT.NE.J) GOTO 100
               PARAM=PARSIJ(I)
               IF(KFN.NE.0)THEN
                  WRITE(6,'(10X,A6,11X,A2,F17.6)')
     1PARTYP(IPARAM)(:3)//NUMBRS(KFN)//'  ',
     2ELEMNT(IELMNT),PARAM
               ELSE
                  WRITE(6,'(10X,A6,11X,A2,F17.6)')
     1PARTYP(IPARAM)//NUMBRS(KFN),
     2ELEMNT(IELMNT),PARAM
               ENDIF
               CALL UPDATE(IPARAM,IELMNT,PARAM,1,KFN)
  100       CONTINUE
  110    CONTINUE
  120 CONTINUE
      CALL MOLDAT(1)                                                    0304WH93
      CALL CALPAR
      ATHEAT=0.D0
      ETH=0.D0
      DO 130 I=1,NUMAT
         NI=NAT(I)
         ATHEAT=ATHEAT+EHEAT(NI)
  130 ETH=ETH+EISOL(NI)
      ATHEAT=ATHEAT-ETH*23.061D0
      RETURN
      END
C
      SUBROUTINE ANALYT(PSUM,PALPHA,PBETA,COORD,NAT,JJA,JJD,
     1IIA,IID,NORB,ENG,KREP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION COORD(3,*),ENG(3), PSUM(*), PALPHA(*), PBETA(*),NAT(*)
************************************************************************
*                                                                      *
*         CALCULATION OF ANALYTICAL DERIVATIVES                        *
*                                                                      *
************************************************************************
C
C COMMON BLOCKS 'OWNED' BY REST OF PROGRAM.
C
      COMMON /CORE  / CORE(107)
     1       /BETAS / BETAS(107),BETAP(107),BETAD(107)
     2       /EXPONT/ ZS(107),ZP(107),ZD(107)
     3       /ALPHA / ALPA(107)
      COMMON /TWOEL3/ F03(107)
      COMMON /NATORB/ NATORB(107)
      COMMON /ALPHA3/ ALP3(153)
      COMMON /IDEAS / FN1(107,10),FN2(107,10),FN3(107,10)
      COMMON /WMATRX/ W(N2ELEC), DUMMY(N2ELEC)
      COMMON /NATYPE/ NZTYPE(107),MTYPE(30),LTYPE
      COMMON /BETA3 / BETA3(153)
      COMMON /VSIPS / VS(107),VP(107),VD(107)
      COMMON /KEYWRD/ KEYWRD
C
C COMMON BLOCKS 'OWNED' BY ANT
C
      COMMON /DERIVS/ DS(16),DG(22),DR(100),TDX(3),TDY(3),TDZ(3)
      COMMON /EXTRA/  G(22), TXYZ(9)
      COMMON /NUMCAL/ NUMCAL
C
C ON RETURN, ENG HOLDS ANALYTICAL DERIVATIVES
C
      COMMON /FORCE3/ IDMY(5),I3N,IX
      COMMON /SRPL/ ISSRP                                               IR0494
      DIMENSION EAA(3),EAB(3),ENUC(3), BI(4), BJ(4),TRMBET(4,4)         IR0494
      CHARACTER*80 KEYWRD
      LOGICAL AM1, MINDO3, ISSRP                                        IR0494
         SAVE                                                           GL0892
      DATA ICALCN /0/
      IF (ICALCN .NE. NUMCAL) THEN
         AM1=(INDEX(KEYWRD,'AM1')+INDEX(KEYWRD,'PM3').NE.0)             IR0492
         MINDO3=(INDEX(KEYWRD,'MINDO').NE.0)
         ICALCN = NUMCAL
      ENDIF
      A0=0.529167D0
      JD=JJD-JJA+1
      JA=1
      ID=IID-IIA+1+JD
      IA=JD+1
      DO 10 J=1,3
         EAA(J)=0.0D0
         EAB(J)=0.0D0
         ENUC(J)=0.0D0
         ENG(J)=0.0D0
   10 CONTINUE
      I=2
      NI=NAT(I)
      ISTART=NZTYPE(NI)*4-3
      J=1
      NJ=NAT(J)
      JSTART=NZTYPE(NJ)*4-3
      R2=(COORD(1,I)-COORD(1,J))**2+(COORD(2,I)-COORD(2,J))**2
     1   + (COORD(3,I)-COORD(3,J))**2
      RIJ=SQRT(R2)
      R0=RIJ/A0
      RR=R2/(A0*A0)
      DO 150 IX=1,3
         DEL1=COORD(IX,I)-COORD(IX,J)
         TERMAA=0.0D0
         TERMAB=0.0D0
         ISP=0
         IOL=0
C   THE FIRST DERIVATIVES OF OVERLAP INTEGRALS
         DO 30 K=IA,ID
            KA=K-IA
            KG=ISTART+KA
            DO 30 L=JA,JD
               LA=L-JA
               LG=JSTART+LA
               IOL=IOL+1
               DS(IOL)=0.0D0
               IF(KA.EQ.0.AND.LA.EQ.0) THEN
C   (S/S) TERM
                  IF(ABS(DEL1).LE.1.0D-6) GO TO 30
                  IS=1
               ELSEIF(KA.EQ.0.AND.LA.GT.0) THEN
C   (S/P) TERM
                  IS=3
                  IF(IX.EQ.LA) GO TO 20
                  IF(ABS(DEL1).LE.1.0D-6) GO TO 30
                  IS=2
                  DEL2=COORD(LA,I)-COORD(LA,J)
               ELSEIF(KA.GT.0.AND.LA.EQ.0) THEN
C   (P/S) TERM
                  IS=5
                  IF(IX.EQ.KA) GO TO 20
                  IF(ABS(DEL1).LE.1.0D-6) GO TO 30
                  IS=4
                  DEL2=COORD(KA,I)-COORD(KA,J)
               ELSE
C   (P/P) TERM
                  IF(KA.EQ.LA) THEN
C    P/P
                     IS=9
                     IF(IX.EQ.KA) GO TO 20
                     IF(ABS(DEL1).LE.1.0D-6) GO TO 30
C    P'/P'
                     IS=8
                     DEL2=COORD(KA,I)-COORD(KA,J)
                  ELSEIF(IX.NE.KA.AND.IX.NE.LA) THEN
C    P'/P"
                     IF(ABS(DEL1).LE.1.0D-6) GO TO 30
                     IS=7
                     DEL2=COORD(KA,I)-COORD(KA,J)
                     DEL3=COORD(LA,I)-COORD(LA,J)
                  ELSE
C    P/P' OR P'/P
                     DEL2=COORD(KA+LA-IX,I)-COORD(KA+LA-IX,J)
                     IS=6
                  ENDIF
               ENDIF
C
C        CALCULATE OVERLAP DERIVATIVES, STORE RESULTS IN DS
C
   20          CALL DERS(KG,LG,RR,DEL1,DEL2,DEL3,IS,IOL)
   30    CONTINUE
         IF((IX.EQ.1).and.(.not.mindo3)) READ (2) (G(I22),I22=1,22)
         IF(.NOT.MINDO3) CALL DELRI(DG,NI,NJ,R0,DEL1)
         CALL DELMOL(COORD,I,J,NI,NJ,IA,ID,JA,JD,IX,RIJ,DEL1,ISP)
C
C   THE FIRST DERIVATIVE OF NUCLEAR REPULSION TERM
         IF(MINDO3)THEN
            II=MAX(NI,NJ)
            NBOND=(II*(II-1))/2+NI+NJ-II
            ALPHA=0
            IF(NBOND.LT.154)THEN
               ALPHA=ALP3(NBOND)
            ELSE
               IF(NATORB(NI).EQ.0)ALPHA=ALPA(NI)
               IF(NATORB(NJ).EQ.0)ALPHA=ALPHA+ALPA(NJ)
            ENDIF
            C2=(7.1995D0/F03(NI)+7.1995D0/F03(NJ))**2
            C1=DEL1/RIJ*CORE(NI)*CORE(NJ)*14.399D0
            IF(NBOND.EQ.22.OR.NBOND.EQ.29)THEN
               TERMNC=-C1*ALPHA*(1.D0/RIJ**2 - RIJ*(RIJ**2+C2)**(-1.5D0)
     1 +  1.D0/RIJ - 1.D0/SQRT(RIJ**2+C2)) * EXP(-RIJ) -
     2C1*RIJ*(RIJ**2+C2)**(-1.5D0)
            ELSEIF(RIJ.LT.1.D0.AND.ALPHA.NE.0.D0)THEN
               TERMNC=0.D0
            ELSE
               TERMNC=-C1*(1.D0/RIJ**2 - RIJ*(RIJ**2+C2)**(-1.5D0) +
     1ALPHA/RIJ - ALPHA/SQRT(RIJ**2+C2)) * EXP(-ALPHA*RIJ) -
     2C1*RIJ*(RIJ**2+C2)**(-1.5D0)
            ENDIF
            DR1=DEL1/RIJ*14.399D0*RIJ*(RIJ**2+C2)**(-1.5D0)
         ELSE
C
C      CORE-CORE TERMS, MNDO AND AM1
C
C
C  SPECIAL TREATMENT FOR N-H AND O-H TERMS
C
            IF(RIJ.LT.1.D0.AND.NATORB(NI)*NATORB(NJ).EQ.0)THEN
               TERMNC=0.D0
               GOTO 50
            ENDIF
            IF(NI.EQ.1.AND.(NJ.EQ.7.OR.NJ.EQ.8)) THEN
               F3=1.0D0+EXP(-ALPA(1)*RIJ)+RIJ*EXP(-ALPA(NJ)*RIJ)
               DD=DG(1)*F3-G(1)*(DEL1/RIJ)*(ALPA(1)*EXP(-ALPA(1)*RIJ)
     1 +(ALPA(NJ)*RIJ-1.0D0)*EXP(-ALPA(NJ)*RIJ))
            ELSEIF((NI.EQ.7.OR.NI.EQ.8).AND.NJ.EQ.1) THEN
               F3=1.0D0+EXP(-ALPA(1)*RIJ)+RIJ*EXP(-ALPA(NI)*RIJ)
               DD=DG(1)*F3-G(1)*(DEL1/RIJ)*(ALPA(1)*EXP(-ALPA(1)*RIJ)
     1 +(ALPA(NI)*RIJ-1.0D0)*EXP(-ALPA(NI)*RIJ))
            ELSE
               F3=1.0D0+EXP(-ALPA(NI)*RIJ)+EXP(-ALPA(NJ)*RIJ)
               DD=DG(1)*F3-G(1)*(DEL1/RIJ)*(ALPA(NI)*EXP(-ALPA(NI)*RI
     1J) +ALPA(NJ)*EXP(-ALPA(NJ)*RIJ))
            ENDIF
            TERMNC=CORE(NI)*CORE(NJ)*DD
         ENDIF
C
C   ****   START OF THE AM1 SPECIFIC DERIVATIVE CODE   ***
C
C      ANALYT=-A*(1/(R*R)+2.D0*B*(R-C)/R)*EXP(-B*(R-C)**2)
C
C    ANALYTICAL DERIVATIVES
C
         IF( AM1 )THEN
            ANAM1=0.D0
            DO 40 IG=1,10
               IF(ABS(FN1(NI,IG)).GT.0.D0)
     1ANAM1=ANAM1+FN1(NI,IG)*
     2(1.D0/(RIJ*RIJ)+2.D0*FN2(NI,IG)*(RIJ-FN3(NI,IG))/RIJ)*
     3EXP(MAX(-30.D0,-FN2(NI,IG)*(RIJ-FN3(NI,IG))**2))
               IF(ABS(FN1(NJ,IG)).GT.0.D0)
     1ANAM1=ANAM1+FN1(NJ,IG)*
     2(1.D0/(RIJ*RIJ)+2.D0*FN2(NJ,IG)*(RIJ-FN3(NJ,IG))/RIJ)*
     3EXP(MAX(-30.D0,-FN2(NJ,IG)*(RIJ-FN3(NJ,IG))**2))
   40       CONTINUE
            ANAM1=ANAM1*CORE(NI)*CORE(NJ)
            TERMNC=TERMNC-ANAM1*DEL1/RIJ
         ENDIF
C
C   ****   END OF THE AM1 SPECIFIC DERIVATIVE CODE   ***
C
   50    CONTINUE
C
C   COMBINE TOGETHER THE OVERLAP DERIVATIVE PARTS
C
         IF(MINDO3)THEN
            II=MAX(NI,NJ)
            NBOND=(II*(II-1))/2+NI+NJ-II
            IF(NBOND.GT.153)GOTO 60
            BI(1)=BETA3(NBOND)*VS(NI)*2.D0
            BI(2)=BETA3(NBOND)*VP(NI)*2.D0
            BI(3)=BI(2)
            BI(4)=BI(2)
            BJ(1)=BETA3(NBOND)*VS(NJ)*2.D0
            BJ(2)=BETA3(NBOND)*VP(NJ)*2.D0
            BJ(3)=BJ(2)
            BJ(4)=BJ(2)
   60       CONTINUE
         ELSE
            BI(1)=BETAS(NI)
            BI(2)=BETAP(NI)
            BI(3)=BI(2)
            BI(4)=BI(2)
            BJ(1)=BETAS(NJ)
            BJ(2)=BETAP(NJ)
            BJ(3)=BJ(2)
            BJ(4)=BJ(2)
         ENDIF
C
C       CODE COMMON TO MINDO/3, MNDO, AND AM1
C
* 
* ....first generate pairwise BETA                                      IR0494
         IF(ISSRP) THEN                                                 IR0494
            TRMBET(1,1)=2.0d0 *GETBET( NI, NJ, 'BETSS')                 IR0494
            BSP=2.0d0 *GETBET(NJ, NI, 'BETSP')                          IR0494
            DO 62 IBET=2,4                                              IR0494
  62           TRMBET(IBET,1)=BSP                                       IR0494
            BSP=2.0d0 *GETBET(NI, NJ, 'BETSP')                          IR0494
            DO 64 JBET=2,4                                              IR0494
  64           TRMBET(1,JBET)=BSP                                       IR0494
            BPP=2.0d0 *GETBET(NI, NJ, 'BETPP')                          IR0494
            DO 66 JBET=2,4                                              IR0494
               DO 66 IBET=2,4                                           IR0494
  66           TRMBET(IBET,JBET)=BPP                                    IR0494
         ELSE                                                           IR0494
            DO 68 JBET=1,4                                              IR0494
               DO 68 IBET=1,4                                           IR0494
  68              TRMBET(IBET,JBET)=BI(IBET)+BJ(JBET)                   IR0494
         ENDIF                                                          IR0494
         IF (index(keywrd,' DEBUG') .ne. 0) THEN                        IR0494
            Write(6,'("BETA", I4, I4)') NI, NJ                          IR0494
            do 69 ibet=1,4                                              IR0494
  69           Write(6,'(4F10.6)') (TRMBET(IBET,JBET),JBET=1,4)         IR0494
         ENDIF                                                          IR0494
         IOL=0
         DO 70 K=IA,ID
            DO 70 L=JA,JD
               LK=L+K*(K-1)/2
               IOL=IOL+1
               TERMAB=TERMAB+TRMBET((K-IA+1),(L-JA+1))*PSUM(LK)*DS(IOL) IR0494
   70    CONTINUE
         IF(MINDO3)THEN
C
C        FIRST, CORE-ELECTRON ATTRACTION DERIVATIVES (MINDO/3)
C
C          ATOM CORE I AFFECTING A.O.S ON J
            DO 80 M=JA,JD
               MN=(M*(M+1))/2
   80       TERMAB=TERMAB+CORE(NI)*PSUM(MN)*DR1
C          ATOM CORE J AFFECTING A.O.S ON I
            DO 90 M=IA,ID
               MN=(M*(M+1))/2
   90       TERMAB=TERMAB+CORE(NJ)*PSUM(MN)*DR1
C
C   NOW FOR COULOMB AND EXCHANGE TERMS (MINDO/3)
C
            DO 100 I1=IA,ID
               II=(I1*(I1+1))/2
               DO 100 J1=JA,JD
                  JJ=(J1*(J1+1))/2
                  IJ=J1+(I1*(I1-1))/2
C
C           COULOMB TERM
C
                  TERMAA=TERMAA-PSUM(II)*DR1*PSUM(JJ)
C
C           EXCHANGE TERM
C
                  TERMAA=TERMAA+(PALPHA(IJ)*PALPHA(IJ)+PBETA(IJ)*PBETA(I
     1J))*DR1
  100       CONTINUE
         ELSE
C
C        FIRST, CORE-ELECTRON ATTRACTION DERIVATIVES (MNDO AND AM1)
C
C          ATOM CORE I AFFECTING A.O.S ON J
            ISP=0
            DO 110 M=JA,JD
               BB=1.D0
               DO 110 N=M,JD
                  MN=M+N*(N-1)/2
                  ISP=ISP+1
                  TERMAB=TERMAB-BB*CORE(NI)*PSUM(MN)*DR(ISP)
  110       BB=2.D0
C          ATOM CORE J AFFECTING A.O.S ON I
            K=MAX(JD-JA+1,1)
            K=(K*(K+1))/2
            ISP=-K+1
            DO 120 M=IA,ID
               BB=1.D0
               DO 120 N=M,ID
                  MN=M+N*(N-1)/2
                  ISP=ISP+K
                  TERMAB=TERMAB-BB*CORE(NJ)*PSUM(MN)*DR(ISP)
  120       BB=2.D0
            ISP=0
C
C   NOW FOR COULOMB AND EXCHANGE TERMS (MNDO AND AM1)
C
            DO 140 K=IA,ID
               AA=1.D0
               KK=(K*(K-1))/2
               DO 140 L=K,ID
                  LL=(L*(L-1))/2
                  DO 130 M=JA,JD
                     BB=1.D0
                     DO 130 N=M,JD
                        ISP=ISP+1
                        KL=K+LL
                        MN=M+N*(N-1)/2
C
C    COULOMB TERM
C
                        TERMAA=TERMAA+AA*BB*PSUM(KL)*PSUM(MN)*DR(ISP)
                        MK=M+KK
                        NK=N+KK
                        ML=M+LL
                        NL=N+LL
C
C    EXCHANGE TERM
C
                        TERMAA= TERMAA-0.5D0*AA*BB*(PALPHA(M
     1K)*PALPHA(NL)+PALPHA(NK)*PALPHA(ML)+PBETA(MK)*PBETA(NL)+PBETA(NK)*
     2PBETA(ML))*DR(ISP)
  130             BB=2.D0
  140       AA=2.D0
C           END OF MNDO AND AM1 SPECIFIC CODE
         ENDIF
         EAA(IX)=EAA(IX)+TERMAA
         EAB(IX)=EAB(IX)+TERMAB
         ENUC(IX)=ENUC(IX)+TERMNC
  150 CONTINUE
C#            WRITE(6,*)EAA,EAB,ENUC,NAT(1),NAT(2)
  160 CONTINUE
  170 CONTINUE
      DO 180 J=1,3
         ENG(J)=EAA(J)+EAB(J)+ENUC(J)
         ENG(J) = -ENG(J)*2.D0*23.061D0
  180 CONTINUE
      RETURN
      END
      BLOCK DATA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C   This blockdata file has been updated so that the Aluminum, Phosphorus,
C   and Sulfur AM1 parameters are from MOPAC-version 6.0, as well as the
C   following PM3 parameters: Beryllium, Magnesium, Zinc, Gallium, Germanium,
C   Arsenic, Selenium, Cadmium, Indium, Tin, Antimony, Tellurium, Mercury, 
C   Thallium, Lead, and Bismuth.  
C
      COMMON /NATORB/ NATORB(107)
***********************************************************************
*
*     COMMON BLOCKS FOR AM1
*
***********************************************************************
      COMMON /ELEMTS/ ELEMNT(107)
     1       /ALPHA / ALP(107)
     2       /CORE  / CORE(107)
     3       /MULTIP/ DD(107),QQ(107),AM(107),AD(107),AQ(107)
     4       /EXPONT/ ZS(107),ZP(107),ZD(107)
     5       /ONELEC/ USS(107),UPP(107),UDD(107)
     6       /BETAS / BETAS(107),BETAP(107),BETAD(107)
     7       /TWOELE/ GSS(107),GSP(107),GPP(107),GP2(107),HSP(107),
     8                GSD(107),GPD(107),GDD(107)
     9       /ATOMIC/ EISOL(107),EHEAT(107)
     1       /VSIPS / VS(107),VP(107),VD(107)
     2       /ISTOPE/ AMS(107)
     3       /IDEAS / GUESS1(107,10),GUESS2(107,10),GUESS3(107,10)
     4       /IDEAP / GUESP1(107,10),GUESP2(107,10),GUESP3(107,10)
     5       /GAUSS / FN1(107),FN2(107)
***********************************************************************
*
*     COMMON BLOCKS FOR MNDO
*
***********************************************************************
      COMMON /MNDO/  USSM(107), UPPM(107), UDDM(107), ZSM(107),
     1ZPM(107), ZDM(107), BETASM(107), BETAPM(107), BETADM(107),
     2ALPM(107), EISOLM(107), DDM(107), QQM(107), AMM(107), ADM(107),
     3AQM(107) ,GSSM(107), GSPM(107), GPPM(107), GP2M(107), HSPM(107),
     4POLVOM(107)
***********************************************************************
*
*     COMMON BLOCKS FOR PM3
*
***********************************************************************
      COMMON /PM3 /  USSPM3(107), UPPPM3(107), UDDPM3(107), ZSPM3(107),
     1ZPPM3(107), ZDPM3(107), BETASP(107), BETAPP(107), BETADP(107),
     2ALPPM3(107), EISOLP(107), DDPM3(107), QQPM3(107), AMPM3(107),
     3ADPM3(107), AQPM3(107) ,GSSPM3(107), GSPPM3(107), GPPPM3(107),
     4GP2PM3(107), HSPPM3(107),POLVOP(107)
     5       /REFS/ REFMN(107), REFM3(107), REFAM(107), REFPM3(107)
*
*  COMMON BLOCKS FOR MINDO/3
*
      COMMON /ONELE3 /  USS3(18),UPP3(18)
     1       /TWOEL3 /  F03(107)
     2       /ATOMI3 /  EISOL3(18),EHEAT3(18)
     3       /BETA3  /  BETA3(153)
     4       /ALPHA3 /  ALP3(153)
     5       /EXPON3 /  ZS3(18),ZP3(18)
*
*  END OF MINDO/3 COMMON BLOCKS
*
C
C ELECTRIC FIELD OPTIONS FOR POLARIZABILITY
C
      COMMON /FIELD/ EFIELD(3)
      CHARACTER ELEMNT*2, REFMN*80, REFM3*80, REFAM*80, REFPM3*80
      DATA EFIELD/0.0D00,0.0D00,0.0D00/
      DATA ELEMNT/' H','He',
     1 'Li','Be',' B',' C',' N',' O',' F','Ne',
     2 'Na','Mg','Al','Si',' P',' S','Cl','Ar',
     3 ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu',
     4 'Zn','Ga','Ge','As','Se','Br','Kr',
     5 'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag',
     6 'Cd','In','Sn','Sb','Te',' I','Xe',
     7 'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',
     8 'Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re','Os','Ir','Pt',
     9 'Au','Hg','Tl','Pb','Bi','Po','At','Rn',
     1 'Fr','Ra','Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','XX',
     2 'Fm','Md','Cb','++',' +','--',' -','Tv'/
C
C   NATORB IS THE NUMBER OF ATOMIC ORBITALS PER ATOM.
C
      DATA NATORB/2*1, 4, 7*4, 0, 7*4, 0, 4, 9*9, 7*4,
     12*4, 9*9, 7*4, 2*2, 14*8, 9*9, 7*4, 15*0,1,5*0/
***********************************************************************
*                      VALENCE SHELLS ARE DEFINED AS                  *
*  PQN   VALENCE SHELLS                                               *
*                 P-GROUP              F-GROUP    TRANSITION METALS   *
*   1       1S                                                        *
*   2       2S 2P                                                     *
*   3       3S 3P  OR  3S 3P 3D                                       *
*   4       4S 4P                                    4S 4P 3D         *
*   5       5S 5P                                    5S 5P 4D         *
*   6       6S 6P                       6S 4F        6S 6P 5D         *
*   7  NOT ASSIGNED YET  ****DO  NOT  USE****                         *
***********************************************************************
      DATA      POLVOM(1) /0.2287D0/
      DATA      POLVOM(6) /0.2647D0/
      DATA      POLVOM(7) /0.3584D0/
      DATA      POLVOM(8) /0.2324D0/
      DATA      POLVOM(9) /0.1982D0/
      DATA      POLVOM(17)/1.3236D0/
      DATA      POLVOM(35)/2.2583D0/
      DATA      POLVOM(53)/4.0930D0/
C
C                STANDARD ATOMIC MASSES
C
      DATA  AMS /  1.00790D0,  4.00260D0,  6.94000D0,  9.01218D0,
     110.81000D0, 12.01100D0, 14.00670D0, 15.99940D0, 18.99840D0,
     220.17900D0, 22.98977D0, 24.30500D0, 26.98154D0, 28.08550D0,
     330.97376D0, 32.06000D0, 35.45300D0, 39.94800D0, 39.09830D0,
     440.08000D0, 44.95590D0, 47.90000D0, 50.94150D0, 51.99600D0,
     554.93800D0, 55.84700D0, 58.93320D0, 58.71000D0, 63.54600D0,
     665.38000D0, 69.73500D0, 72.59000D0, 74.92160D0, 78.96000D0,
     779.90400D0, 83.80000D0, 85.46780D0, 87.62000D0, 88.90590D0,
     891.22000D0, 92.90640D0, 95.94000D0, 98.90620D0, 101.0700D0,
     9102.9055D0, 106.4000D0, 107.8680D0, 112.4100D0, 114.8200D0,
     1118.6900D0, 121.7500D0, 127.6000D0, 126.9045D0, 131.3000D0,
     2132.9054D0, 137.3300D0, 15*0.000D0, 178.4900D0, 180.9479D0,
     3183.8500D0, 186.2070D0, 190.2000D0, 192.2200D0, 195.0900D0,
     4196.9665D0, 200.5900D0, 204.3700D0, 207.2000D0, 208.9804D0,
     518*0.000D0,   1.0079D0,  5*0.000D0/
C
C   CORE IS THE CHARGE ON THE ATOM AS SEEN BY THE ELECTRONS
C
      DATA CORE/1.D0,0.D0,
     1 1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,0.D0,
     2 1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,0.D0,
     3 1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,9.D0,10.D0,11.D0,2.D0,
     4 3.D0,4.D0,5.D0,6.D0,7.D0,0.D0,
     5 1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,9.D0,10.D0,11.D0,2.D0,
     6 3.D0,4.D0,5.D0,6.D0,7.D0,0.D0,
     7 1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,9.D0,10.D0,
     8 11.D0,12.D0,13.D0,14.D0,15.D0,16.D0,
     9 3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,9.D0,10.D0,11.D0,2.D0,
     1 3.D0,4.D0,5.D0,6.D0,7.D0,0.D0,
     2  15*0.D0,1.D0,2.D0,1.D0,-2.D0,-1.D0,0.D0/
C
C     ENTHALPIES OF FORMATION OF GASEOUS ATOMS ARE TAKEN FROM \ANNUAL
C     REPORTS,1974,71B,P 117\  THERE ARE SOME SIGNIFICANT DIFFERENCES
C     BETWEEN THE VALUES REPORTED THERE AND THE VALUES PREVIOUSLY IN
C     THE BLOCK DATA OF THIS PROGRAM.  ONLY THE THIRD  ROW ELEMENTS
C     HAVE BEEN UPDATED.
C
* ALL THE OTHER ELEMENTS ARE TAKEN FROM CRC HANDBOOK 1981-1982
      DATA EHEAT(1)  / 52.102D0/
      DATA EHEAT(2)  /  0.000D0/
C
      DATA EHEAT(3)  / 38.410D0/
      DATA EHEAT(4)  / 76.960D0/
      DATA EHEAT(5)  /135.700D0/
      DATA EHEAT(6)  /170.890D0/
      DATA EHEAT(7)  /113.000D0/
      DATA EHEAT(8)  / 59.559D0/
      DATA EHEAT(9)  / 18.890D0/
      DATA EHEAT(10) /  0.000D0/
C
C#      DATA EHEAT(11) / 25.850D0/
      DATA EHEAT(12) / 35.000D0/
      DATA EHEAT(13) / 79.490D0/
      DATA EHEAT(14) /108.390D0/
      DATA EHEAT(15) / 75.570D0/
      DATA EHEAT(16) / 66.400D0/
      DATA EHEAT(17) / 28.990D0/
      DATA EHEAT(18) /  0.000D0/
C
C#      DATA EHEAT(19) / 21.420D0/
      DATA EHEAT(20) / 42.600D0/
      DATA EHEAT(21) / 90.300D0/
      DATA EHEAT(22) /112.300D0/
      DATA EHEAT(23) /122.900D0/
      DATA EHEAT(24) / 95.000D0/
      DATA EHEAT(25) / 67.700D0/
      DATA EHEAT(26) / 99.300D0/
      DATA EHEAT(27) /102.400D0/
      DATA EHEAT(28) /102.800D0/
      DATA EHEAT(29) / 80.700D0/
      DATA EHEAT(30) / 31.170D0/
      DATA EHEAT(31) / 65.400D0/
      DATA EHEAT(32) / 89.500D0/
      DATA EHEAT(33) / 72.300D0/
      DATA EHEAT(34) / 54.300D0/
      DATA EHEAT(35) / 26.740D0/
      DATA EHEAT(36) /  0.000D0/
C
      DATA EHEAT(37) / 19.600D0/
      DATA EHEAT(38) / 39.100D0/
      DATA EHEAT(39) /101.500D0/
      DATA EHEAT(40) /145.500D0/
      DATA EHEAT(41) /172.400D0/
      DATA EHEAT(42) /157.300D0/
      DATA EHEAT(44) /155.500D0/
      DATA EHEAT(45) /133.000D0/
      DATA EHEAT(46) / 90.000D0/
      DATA EHEAT(47) / 68.100D0/
      DATA EHEAT(48) / 26.720D0/
      DATA EHEAT(49) / 58.000D0/
      DATA EHEAT(50) / 72.200D0/
      DATA EHEAT(51) / 63.200D0/
      DATA EHEAT(52) / 47.000D0/
      DATA EHEAT(53) / 25.517D0/
      DATA EHEAT(54) /  0.000D0/
C
      DATA EHEAT(55) / 18.700D0/
      DATA EHEAT(56) / 42.500D0/
      DATA EHEAT(58) /101.300D0/
      DATA EHEAT(62) / 49.400D0/
      DATA EHEAT(68) / 75.800D0/
      DATA EHEAT(70) / 36.350D0/
      DATA EHEAT(72) /148.000D0/
      DATA EHEAT(73) /186.900D0/
      DATA EHEAT(74) /203.100D0/
      DATA EHEAT(75) /185.000D0/
      DATA EHEAT(76) /188.000D0/
      DATA EHEAT(77) /160.000D0/
      DATA EHEAT(78) /135.200D0/
      DATA EHEAT(79) / 88.000D0/
      DATA EHEAT(80) / 14.690D0/
      DATA EHEAT(81) / 43.550D0/
      DATA EHEAT(82) / 46.620D0/
      DATA EHEAT(83) / 50.100D0/
      DATA EHEAT(86) /  0.000D0/
      DATA EHEAT(102)  / 207.0D0/
C
      DATA VS(1) /  -13.605  /
      DATA VS(5)/-15.16D00/
      DATA VS(6)/-21.34D00/
      DATA VS(7)/-27.51D00/
      DATA VS(8)/-35.30D00/
      DATA VS(9)/-43.70D00/
      DATA VS(14)/-17.82D00/
      DATA VS(15)/-21.10D00/
      DATA VS(16)/-23.84D00/
      DATA VS(17)/-25.26D00/
      DATA VP(1)  /  0.0D00  /
      DATA VP(5)/-8.52D00/
      DATA VP(6)/-11.54D00/
      DATA VP(7)/-14.34D00/
      DATA VP(8)/-17.91D00/
      DATA VP(9)/-20.89D00/
      DATA VP(14)/-8.51D00/
      DATA VP(15)/-10.29D00/
      DATA VP(16)/-12.41D00/
      DATA VP(17)/-15.09D00/
C      DATA NPQ/1,1, 2,2,2,2,2,2,2,2, 3,3,3,3,3,3,3,3, 4,4,4,4,4,4,4,4,
C     +4,4,4,4,4,4,4,4,4,4, 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5/
C
C *** ONE CENTER REPULSION INTEGRALS
C     GSS ::= (SS,SS)
C     GPP ::= (PP,PP)
C     GSP ::= (SS,PP)
C     GP2 ::= (PP,P*P*)
C     HSP ::= (SP,SP)
************************************************************************
      DATA GSSM(1) / 12.848D00 /
      DATA GSSM(4)/9.00D00/
      DATA GSSM(5)/10.59D00/
      DATA GSSM(6) / 12.23D00 /
      DATA GSSM(7)/13.59D00/
      DATA GSSM(8)/15.42D00/
      DATA GSSM(9)/16.92D00/
      DATA GSSM(13)/8.09D00/
      DATA GSSM(14)/9.82D00/
      DATA GSSM(15)/11.56D00/
      DATA GSSM(16)/12.88D00/
      DATA GSSM(17)/15.03D00/
      DATA GSSM(35)/15.03643948D0/
      DATA GSSM(53)/15.04044855D0/
      DATA GPPM(4)/6.97D00/
      DATA GPPM(5)/8.86D00/
      DATA GPPM(6) / 11.08D00 /
      DATA GPPM(7)/12.98D00/
      DATA GPPM(8)/14.52D00/
      DATA GPPM(9)/16.71D00/
      DATA GPPM(13)/5.98D00/
      DATA GPPM(14)/7.31D00/
      DATA GPPM(15)/8.64D00/
      DATA GPPM(16)/9.90D00/
      DATA GPPM(17)/11.30D00/
      DATA GPPM(35)/11.27632539D0/
      DATA GPPM(53)/11.14778369D0/
      DATA GSPM(4)/7.43D00/
      DATA GSPM(5)/9.56D00/
      DATA GSPM(6) / 11.47D00 /
      DATA GSPM(7)/12.66D00/
      DATA GSPM(8)/14.48D00/
      DATA GSPM(9)/17.25D00/
      DATA GSPM(13)/6.63D00/
      DATA GSPM(14)/8.36D00/
      DATA GSPM(15)/10.08D00/
      DATA GSPM(16)/11.26D00/
      DATA GSPM(17)/13.16D00/
      DATA GSPM(35)/13.03468242D0/
      DATA GSPM(53)/13.05655798D0/
      DATA GP2M(4)/6.22D00/
      DATA GP2M(5)/7.86D00/
      DATA GP2M(6) / 9.84D00 /
      DATA GP2M(7)/11.59D00/
      DATA GP2M(8)/12.98D00/
      DATA GP2M(9)/14.91D00/
      DATA GP2M(13)/5.40D00/
      DATA GP2M(14)/6.54D00/
      DATA GP2M(15)/7.68D00/
      DATA GP2M(16)/8.83D00/
      DATA GP2M(17)/9.97D00/
      DATA GP2M(35)/9.85442552D0/
      DATA GP2M(53)/9.91409071D0/
      DATA HSPM(4)/1.28D00/
      DATA HSPM(5)/1.81D00/
      DATA HSPM(6) / 2.43D00 /
      DATA HSPM(7)/3.14D00/
      DATA HSPM(8)/3.94D00/
      DATA HSPM(9)/4.83D00/
      DATA HSPM(13)/0.70D00/
      DATA HSPM(14)/1.32D00/
      DATA HSPM(15)/1.92D00/
      DATA HSPM(16)/2.26D00/
      DATA HSPM(17)/2.42D00/
      DATA HSPM(35)/2.45586832D0/
      DATA HSPM(53)/2.45638202D0/
C
C     THE MONOCENTRIC INTEGRALS HSP AND GSP FOR ALUMINIUM ARE ONLY
C     ESTIMATES. A VALUE OF G1 FOR AL IS NEEDED TO RESOLVE OLEARIS
C     INTEGRALS.
C
C     OPTIMIZED MNDO PARAMETERS FOR H, BE, B, C, N, O, F
C                                                     CL
C     ESTIMATED MNDO PARAMETERS FOR       AL,SI, P, S
C
C     ELEMENTS H, C, N, O WERE PARAMETERIZED BY WALTER THIEL
C     ELEMENTS B,SI,P,S   WERE      ..          MICHAEL MCKEE
C     ELEMENTS BE,F,AL,CL WERE      ..          HENRY RZEPA
C
***********************************************************************
*
*    START OF MINDO/3 PARAMETERS
*
***********************************************************************
C *** F03 IS THE ONE CENTER AVERAGED REPULSION INTEGRAL FOR USE IN THE
C        TWO CENTER ELECTRONIC REPULSION INTEGRAL EVALUATION.
      DATA REFM3  ( 1)/'  H: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.CHEM.SO
     1C. 97,1285,1294,1302,1307 (1975)'/
      DATA REFM3  ( 5)/'  B: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.CHEM.SO
     1C. 97,1285,1294,1302,1307 (1975)'/
      DATA REFM3  ( 6)/'  C: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.CHEM.SO
     1C. 97,1285,1294,1302,1307 (1975)'/
      DATA REFM3  ( 7)/'  N: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.CHEM.SO
     1C. 97,1285,1294,1302,1307 (1975)'/
      DATA REFM3  ( 8)/'  O: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.CHEM.SO
     1C. 97,1285,1294,1302,1307 (1975)'/
      DATA REFM3  ( 9)/'  F: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.CHEM.SO
     1C. 97,1285,1294,1302,1307 (1975)'/
      DATA REFM3  (14)/' Si: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.CHEM.SO
     1C. 97,1285,1294,1302,1307 (1975)'/
      DATA REFM3  (15)/'  P: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.CHEM.SO
     1C. 97,1285,1294,1302,1307 (1975)'/
      DATA REFM3  (16)/'  S: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.CHEM.SO
     1C. 97,1285,1294,1302,1307 (1975)'/
      DATA REFM3  (17)/' Cl: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.CHEM.SO
     1C. 97,1285,1294,1302,1307 (1975)'/
      DATA F03              /  12.848D0, 10.0D0, 10.0D0, 0.0D0,
     1  8.958D0, 10.833D0, 12.377D0, 13.985D0, 16.250D0,
     2         10.000D0, 10.000D0, 0.000D0, 0.000D0,7.57D0 ,  9.00D0 ,
     3         10.20D0 , 11.73,10.0D0,35*0.D0,10.D0,53*10.D0/
C *** USS AND UPP ARE THE ONE-CENTER CORE ELECTRON ATTRACTION AND KINETI
C     ENERGY INTEGRALS FOR S AND P ELECTRONS RESPECTIVELY IN E.V.
      DATA USS3             / -12.505D0, 0.000D0, 0.000D0, 0.000D0,
     1                       -33.61D0, -51.79D0, -66.06D0, -91.73D0 ,
     2                       -129.86D0,
     3                        0.0000D0 , 0.000 D0 ,0.000D0 , 0.000D0 ,
     4          -39.82D0 , -56.23D0 , -73.39D0 , -98.99D0 ,.0D0/
      DATA UPP3             /   0.0D0, 0.0D0, 0.0D0, 0.0D0,
     1     -25.11D0 , -39.18D0 , -56.40D0 , -78.80D0 , -105.93D0 ,
     2                        0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 ,
     3         -29.15D0 , -42.31D0 , -57.25D0 , -76.43D0 ,.0D0/
C *** EISOL3 AND EHEAT3 ARE THE GS ELECTRONIC ENERGY OF THE NEUTRAL ATOM
C     (IN E.V.) AND THE HEAT OF FORMATION IF THE FREE ATOM (IN KCAL/MOL)
      DATA EISOL3             /-12.505D0 , 0.0D0 , 0.0D0 ,0.0D0 ,
     1        -61.70D0 ,-119.47D0 , -187.51D0 , -307.07D0 , -475.00D0 ,
     2                         0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 ,
     3          -90.98D0 , -150.81D0 , -229.15D0 , -345.93D0 , 0.0D0/
      DATA EHEAT3             / 52.102D0 , 0.0D0 , 0.0D0 , 0.0D0 ,
     1     135.7 D0 , 170.89D0 ,  113.0 D0 ,  59.559D0 ,  18.86D0 ,
     2                         0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 ,
     3     106.0D0 ,   79.8D0 ,  65.65D0 ,  28.95D0 , 0.0D0 /
C *** BETA3 AND ALP3 ARE THE BOND PARAMETERS USED IN THE
C     RESONANCE INTEGRAL AND THE CORE CORE REPULSION INTEGRAL RESPECTIVE
C     THAT IS ACCORDING TO THE FOLLOWING CONVENTION
C
C     HERE IS THE
C     BOND TYPE DESIGNATION
C
C
C         H   B   C   N   O   F  SI   P   S  CL
C       -----------------------------------------
C      H  1  11  16  22  29  37  92 106 121 137
C      B     15  20  26  33  41
C      C         21  27  34  42  97 111 126 142
C      N             28  35  43         127 143
C      O                 36  44         128 144
C      F                     45         129
C     SI                        105
C      P                            120     151
C      S                                136 152
C     CL                                    153
      DATA BETA3(1),ALP3(1)   /  0.244770D0 ,  1.489450D0 /
      DATA BETA3(11),ALP3(11)   /  0.185347D0 ,  2.090352D0 /
      DATA BETA3(15),ALP3(15)   /  0.151324D0 ,  2.280544D0 /
      DATA BETA3(16),ALP3(16)   /  0.315011D0 ,  1.475836D0 /
      DATA BETA3(20),ALP3(20)   /  0.250031D0 ,  2.138291D0 /
      DATA BETA3(21),ALP3(21)   /  0.419907D0 ,  1.371208D0 /
      DATA BETA3(22),ALP3(22)   /  0.360776D0 ,  0.589380D0 /
      DATA BETA3(26),ALP3(26)   /  0.310959D0 ,  1.909763D0 /
      DATA BETA3(27),ALP3(27)   /  0.410886D0 ,  1.635259D0 /
      DATA BETA3(28),ALP3(28) /  0.377342D0 ,  2.029618D0 /
      DATA BETA3(29),ALP3(29) /  0.417759D0 ,  0.478901D0 /
      DATA BETA3(33),ALP3(33) /  0.349745D0 ,  2.484827D0 /
      DATA BETA3(34),ALP3(34) /  0.464514D0 ,  1.820975D0 /
      DATA BETA3(35),ALP3(35) /  0.458110D0 ,  1.873859D0 /
      DATA BETA3(36),ALP3(36) /  0.659407D0 ,  1.537190D0 /
      DATA BETA3(37),ALP3(37) /  0.195242D0 ,  3.771362D0 /
      DATA BETA3(41),ALP3(41) /  0.219591D0 ,  2.862183D0 /
      DATA BETA3(42),ALP3(42) /  0.247494D0 ,  2.725913D0 /
      DATA BETA3(43),ALP3(43) /  0.205347D0 ,  2.861667D0 /
      DATA BETA3(44),ALP3(44) /  0.334044D0 ,  2.266949D0 /
      DATA BETA3(45),ALP3(45) /  0.197464D0 ,  3.864997D0 /
      DATA BETA3(92),ALP3(92) /  0.289647D0 ,  0.940789D0 /
      DATA BETA3(97),ALP3(97) /  0.411377D0 ,  1.101382D0 /
      DATA BETA3(105),ALP3(105) /  0.291703D0 ,  0.918432D0 /
      DATA BETA3(106),ALP3(106) /  0.320118D0 ,  0.923170D0 /
      DATA BETA3(111),ALP3(111) /  0.457816D0 ,  1.029693D0 /
      DATA BETA3(120),ALP3(120) /  0.311790D0 ,  1.186652D0 /
      DATA BETA3(121),ALP3(121) /  0.220654D0 ,  1.700698D0 /
      DATA BETA3(126),ALP3(126) /  0.284620D0 ,  1.761370D0 /
      DATA BETA3(127),ALP3(127) /  0.313170D0 ,  1.878176D0/
      DATA BETA3(128),ALP3(128) /  0.422890D0 ,  2.077240D0 /
      DATA BETA3(129),ALP3(129)  /  0.000000D0 ,  0.000000D0 /
      DATA BETA3(136),ALP3(136) /  0.202489D0 ,  1.751617D0 /
      DATA BETA3(137),ALP3(137) /  0.231653D0 ,  2.089404D0 /
      DATA BETA3(142),ALP3(142) /  0.315480D0 ,  1.676222D0 /
      DATA BETA3(143),ALP3(143) /  0.302298D0 ,  1.817064D0 /
      DATA BETA3(144),ALP3(144) /  0.000000D0 ,  0.000000D0 /
      DATA BETA3(151),ALP3(151) /  0.277322D0 ,  1.543720D0 /
      DATA BETA3(152),ALP3(152) /  0.221764D0 ,  1.950318D0 /
      DATA BETA3(153),ALP3(153) /  0.258969D0 ,  1.792125D0 /
C *** HERE COMES THE OPTIMIZED SLATER_S EXPONENTS FOR THE EVALUATION
C     OF THE OVERLAP INTEGRALS AND MOLECULAR DIPOLE MOMENTS.
      DATA ZS3(1),ZP3(1)      /  1.3D0       ,  0.0D0      /
      DATA ZS3(5),ZP3(5)      /  1.211156D0 ,  0.972826D0 /
      DATA ZS3(6),ZP3(6)      /  1.739391D0 ,  1.709645D0 /
      DATA ZS3(7),ZP3(7)      /  2.704546D0 ,  1.870839D0 /
      DATA ZS3(8),ZP3(8)      /  3.640575D0 ,  2.168448D0 /
      DATA ZS3(9),ZP3(9)      /  3.111270D0 ,  1.41986D0 /
      DATA ZS3(14),ZP3(14)    /  1.629173D0 ,  1.381721D0 /
      DATA ZS3(15),ZP3(15)    /  1.926108D0 ,  1.590665D0 /
      DATA ZS3(16),ZP3(16)    /  1.719480D0 ,  1.403205D0 /
      DATA ZS3(17),ZP3(17)    /  3.430887D0 ,  1.627017D0 /
*************************************************************
*                                                           *
*               DATA FOR THE SPARKLES                       *
*                                                           *
*************************************************************
*                               DATA FOR THE " ++ " SPARKLE
      DATA EHEAT(103)    / 0.0D0/
      DATA VS(103)       /10.0D0/
      DATA ALP(103)      / 1.5D0/
      DATA EISOL(103)    / 0.0D0/
      DATA AM(103)       / 0.5D0/
      DATA ALPM(103)      / 1.5D0/
      DATA EISOLM(103)    / 0.0D0/
      DATA AMM(103)       / 0.5D0/
      DATA ALPPM3(103)      / 1.5D0/
      DATA EISOLP(103)    / 0.0D0/
      DATA AMPM3(103)       / 0.5D0/
*                               DATA FOR THE " + " SPARKLE
      DATA EHEAT(104)    / 0.0D0/
      DATA VS(104)       /10.0D0/
      DATA ALP(104)      / 1.5D0/
      DATA EISOL(104)    / 0.0D0/
      DATA AM(104)       / 0.5D0/
      DATA ALPM(104)      / 1.5D0/
      DATA EISOLM(104)    / 0.0D0/
      DATA AMM(104)       / 0.5D0/
      DATA ALPPM3(104)      / 1.5D0/
      DATA EISOLP(104)    / 0.0D0/
      DATA AMPM3(104)       / 0.5D0/
*                               DATA FOR THE " -- " SPARKLE
      DATA EHEAT(105)    / 0.0D0/
      DATA VS(105)       /10.0D0/
      DATA ALP(105)      / 1.5D0/
      DATA EISOL(105)    / 0.0D0/
      DATA AM(105)       / 0.5D0/
      DATA ALPM(105)      / 1.5D0/
      DATA EISOLM(105)    / 0.0D0/
      DATA AMM(105)       / 0.5D0/
      DATA ALPPM3(105)      / 1.5D0/
      DATA EISOLP(105)    / 0.0D0/
      DATA AMPM3(105)       / 0.5D0/
*                               DATA FOR THE " - " SPARKLE
      DATA EHEAT(106)    / 0.0D0/
      DATA VS(106)       /10.0D0/
      DATA ALP(106)      / 1.5D0/
      DATA EISOL(106)    / 0.0D0/
      DATA AM(106)       / 0.5D0/
      DATA ALPM(106)      / 1.5D0/
      DATA EISOLM(106)    / 0.0D0/
      DATA AMM(106)       / 0.5D0/
      DATA ALPPM3(106)      / 1.5D0/
      DATA EISOLP(106)    / 0.0D0/
      DATA AMPM3(106)       / 0.5D0/
***********************************************************************
*
*    START OF MNDO PARAMETERS
*
***********************************************************************
C                    DATA FOR ELEMENT  1        HYDROGEN
      DATA REFMN  ( 1)/'  H: (MNDO):  M.J.S. DEWAR, W. THIEL, J. AM. CHE
     1M. SOC., 99, 4899, (1977)       '/
      DATA USSM   ( 1)/     -11.9062760D0/
      DATA BETASM ( 1)/      -6.9890640D0/
      DATA ZSM    ( 1)/       1.3319670D0/
      DATA ALPM   ( 1)/       2.5441341D0/
      DATA EISOLM ( 1)/     -11.9062760D0/
      DATA AMM    ( 1)/       0.4721793D0/
      DATA ADM    ( 1)/       0.4721793D0/
      DATA AQM    ( 1)/       0.4721793D0/
C                    DATA FOR ELEMENT  3        LITHIUM
      DATA REFMN  ( 3)/' Li: (MNDO):  TAKEN FROM MNDOC BY W.THIEL,
     1QCPE NO.438, V. 2, P.63, (1982).'/
      DATA USSM   (  3)/      -5.1280000D0/
      DATA UPPM   (  3)/      -2.7212000D0/
      DATA BETASM (  3)/      -1.3500400D0/
      DATA BETAPM (  3)/      -1.3500400D0/
      DATA ZSM    (  3)/       0.7023800D0/
      DATA ZPM    (  3)/       0.7023800D0/
      DATA ALPM   (  3)/       1.2501400D0/
      DATA EISOLM (  3)/      -5.1280000D0/
      DATA GSSM   (  3)/       7.3000000D0/
      DATA GSPM   (  3)/       5.4200000D0/
      DATA GPPM   (  3)/       5.0000000D0/
      DATA GP2M   (  3)/       4.5200000D0/
      DATA HSPM   (  3)/       0.8300000D0/
      DATA DDM    (  3)/       2.0549783D0/
      DATA QQM    (  3)/       1.7437069D0/
      DATA AMM    (  3)/       0.2682837D0/
      DATA ADM    (  3)/       0.2269793D0/
      DATA AQM    (  3)/       0.2614581D0/
C                    DATA FOR ELEMENT  4        BERYLLIUM
      DATA REFMN  ( 4)/' Be: (MNDO):  M.J.S. DEWAR, H.S. RZEPA, J. AM. C
     1HEM. SOC., 100, 777, (1978)     '/
      DATA USSM   ( 4)/     -16.6023780D0/
      DATA UPPM   ( 4)/     -10.7037710D0/
      DATA BETASM ( 4)/      -4.0170960D0/
      DATA BETAPM ( 4)/      -4.0170960D0/
      DATA ZSM    ( 4)/       1.0042100D0/
      DATA ZPM    ( 4)/       1.0042100D0/
      DATA ALPM   ( 4)/       1.6694340D0/
      DATA EISOLM ( 4)/     -24.2047560D0/
      DATA DDM    ( 4)/       1.4373245D0/
      DATA QQM    ( 4)/       1.2196103D0/
      DATA AMM    ( 4)/       0.3307607D0/
      DATA ADM    ( 4)/       0.3356142D0/
      DATA AQM    ( 4)/       0.3846373D0/
C                    DATA FOR ELEMENT  5        BORON
      DATA REFMN  ( 5)/'  B: (MNDO):  M.J.S. DEWAR, M.L. MCKEE, J. AM. C
     1HEM. SOC., 99, 5231, (1977)     '/
      DATA USSM   ( 5)/     -34.5471300D0/
      DATA UPPM   ( 5)/     -23.1216900D0/
      DATA BETASM ( 5)/      -8.2520540D0/
      DATA BETAPM ( 5)/      -8.2520540D0/
      DATA ZSM    ( 5)/       1.5068010D0/
      DATA ZPM    ( 5)/       1.5068010D0/
      DATA ALPM   ( 5)/       2.1349930D0/
      DATA EISOLM ( 5)/     -64.3159500D0/
      DATA DDM    ( 5)/       0.9579073D0/
      DATA QQM    ( 5)/       0.8128113D0/
      DATA AMM    ( 5)/       0.3891951D0/
      DATA ADM    ( 5)/       0.4904730D0/
      DATA AQM    ( 5)/       0.5556979D0/
C                    DATA FOR ELEMENT  6        CARBON
      DATA REFMN  ( 6)/'  C: (MNDO):  M.J.S. DEWAR, W. THIEL, J. AM. CHE
     1M. SOC., 99, 4899, (1977)       '/
      DATA USSM   ( 6)/     -52.2797450D0/
      DATA UPPM   ( 6)/     -39.2055580D0/
      DATA BETASM ( 6)/     -18.9850440D0/
      DATA BETAPM ( 6)/      -7.9341220D0/
      DATA ZSM    ( 6)/       1.7875370D0/
      DATA ZPM    ( 6)/       1.7875370D0/
      DATA ALPM   ( 6)/       2.5463800D0/
      DATA EISOLM ( 6)/    -120.5006060D0/
      DATA DDM    ( 6)/       0.8074662D0/
      DATA QQM    ( 6)/       0.6851578D0/
      DATA AMM    ( 6)/       0.4494671D0/
      DATA ADM    ( 6)/       0.6149474D0/
      DATA AQM    ( 6)/       0.6685897D0/
C                    DATA FOR ELEMENT  7        NITROGEN
      DATA REFMN  ( 7)/'  N: (MNDO):  M.J.S. DEWAR, W. THIEL, J. AM. CHE
     1M. SOC., 99, 4899, (1977)       '/
      DATA USSM   ( 7)/     -71.9321220D0/
      DATA UPPM   ( 7)/     -57.1723190D0/
      DATA BETASM ( 7)/     -20.4957580D0/
      DATA BETAPM ( 7)/     -20.4957580D0/
      DATA ZSM    ( 7)/       2.2556140D0/
      DATA ZPM    ( 7)/       2.2556140D0/
      DATA ALPM   ( 7)/       2.8613420D0/
      DATA EISOLM ( 7)/    -202.5662010D0/
      DATA DDM    ( 7)/       0.6399037D0/
      DATA QQM    ( 7)/       0.5429763D0/
      DATA AMM    ( 7)/       0.4994487D0/
      DATA ADM    ( 7)/       0.7843643D0/
      DATA AQM    ( 7)/       0.8126445D0/
C                    DATA FOR ELEMENT  8        OXYGEN
      DATA REFMN  ( 8)/'  O: (MNDO):  M.J.S. DEWAR, W. THIEL, J. AM. CHE
     1M. SOC., 99, 4899, (1977)       '/
      DATA USSM   ( 8)/     -99.6443090D0/
      DATA UPPM   ( 8)/     -77.7974720D0/
      DATA BETASM ( 8)/     -32.6880820D0/
      DATA BETAPM ( 8)/     -32.6880820D0/
      DATA ZSM    ( 8)/       2.6999050D0/
      DATA ZPM    ( 8)/       2.6999050D0/
      DATA ALPM   ( 8)/       3.1606040D0/
      DATA EISOLM ( 8)/    -317.8685060D0/
      DATA DDM    ( 8)/       0.5346024D0/
      DATA QQM    ( 8)/       0.4536252D0/
      DATA AMM    ( 8)/       0.5667034D0/
      DATA ADM    ( 8)/       0.9592562D0/
      DATA AQM    ( 8)/       0.9495934D0/
C                    DATA FOR ELEMENT  9        FLUORINE
      DATA REFMN  ( 9)/'  F: (MNDO):  M.J.S. DEWAR, H.S. RZEPA, J. AM. C
     1HEM. SOC., 100, 777, (1978)     '/
      DATA USSM   ( 9)/    -131.0715480D0/
      DATA UPPM   ( 9)/    -105.7821370D0/
      DATA BETASM ( 9)/     -48.2904660D0/
      DATA BETAPM ( 9)/     -36.5085400D0/
      DATA ZSM    ( 9)/       2.8484870D0/
      DATA ZPM    ( 9)/       2.8484870D0/
      DATA ALPM   ( 9)/       3.4196606D0/
      DATA EISOLM ( 9)/    -476.6837810D0/
      DATA DDM    ( 9)/       0.5067166D0/
      DATA QQM    ( 9)/       0.4299633D0/
      DATA AMM    ( 9)/       0.6218302D0/
      DATA ADM    ( 9)/       1.0850301D0/
      DATA AQM    ( 9)/       1.0343643D0/
*                               DATA FOR THE SODIUM-LIKE SPARKLE
      DATA REFMN  (11)/' Na: (MNDO):  SODIUM-LIKE SPARKLE.   USE WITH CA
     1RE.                             '/
      DATA EHEAT(11)    / 0.0D0/
      DATA VS(11)       /10.0D0/
      DATA ALP(11)      / 1.32D0/
      DATA EISOL(11)    / 0.0D0/
      DATA AM(11)       / 0.5D0/
      DATA ALPM(11)      / 1.32D0/
      DATA EISOLM(11)    / 0.0D0/
      DATA AMM(11)       / 0.5D0/
C                    DATA FOR ELEMENT 13        ALUMINUM
      DATA REFMN  (13)/' Al: (MNDO):  L.P. DAVIS, ET.AL.  J. COMP. CHEM.
     1, 2, 433, (1981) SEE MANUAL.    '/
      DATA USSM   (13)/     -23.8070970D0/
      DATA UPPM   (13)/     -17.5198780D0/
      DATA BETASM (13)/      -2.6702840D0/
      DATA BETAPM (13)/      -2.6702840D0/
      DATA ZSM    (13)/       1.4441610D0/
      DATA ZPM    (13)/       1.4441610D0/
      DATA ZDM    (13)/       1.0000000D0/
      DATA ALPM   (13)/       1.8688394D0/
      DATA EISOLM (13)/     -44.4840720D0/
      DATA DDM    (13)/       1.3992387D0/
      DATA QQM    (13)/       1.1586797D0/
      DATA AMM    (13)/       0.2973172D0/
      DATA ADM    (13)/       0.2635574D0/
      DATA AQM    (13)/       0.3673560D0/
C                    DATA FOR ELEMENT 14          SILICON
      DATA REFMN  (14)/' Si: (MNDO): M.J.S.DEWAR, ET. AL. ORGANOMETALLIC
     1S  5, 375 (1986)                '/
      DATA USS   M(14)/     -37.0375330D0/
      DATA UPP   M(14)/     -27.7696780D0/
      DATA BETAS M(14)/      -9.0868040D0/
      DATA BETAP M(14)/      -1.0758270D0/
      DATA ZS    M(14)/       1.3159860D0/
      DATA ZP    M(14)/       1.7099430D0/
      DATA ZD    M(14)/       1.0000000D0/
      DATA ALP   M(14)/       2.2053160D0/
      DATA EISOL M(14)/     -82.8394220D0/
      DATA DD    M(14)/       1.2580349D0/
      DATA QQ    M(14)/       0.9785824D0/
      DATA AM    M(14)/       0.3608967D0/
      DATA AD    M(14)/       0.3664244D0/
      DATA AQ    M(14)/       0.4506740D0/
C                    DATA FOR ELEMENT 15        PHOSPHORUS
      DATA REFMN  (15)/'  P: (MNDO): M.J.S.DEWAR, M.L.MCKEE, H.S.RZEPA,
     1J. AM. CHEM. SOC., 100 3607 1978'/
      DATA USSM   (15)/     -56.1433600D0/
      DATA UPPM   (15)/     -42.8510800D0/
      DATA BETASM (15)/      -6.7916000D0/
      DATA BETAPM (15)/      -6.7916000D0/
      DATA ZSM    (15)/       2.1087200D0/
      DATA ZPM    (15)/       1.7858100D0/
      DATA ZDM    (15)/       1.0000000D0/
      DATA ALPM   (15)/       2.4152800D0/
      DATA EISOLM (15)/    -152.9599600D0/
      DATA DDM    (15)/       1.0129699D0/
      DATA QQM    (15)/       0.9370090D0/
      DATA AMM    (15)/       0.4248438D0/
      DATA ADM    (15)/       0.4882420D0/
      DATA AQM    (15)/       0.4979406D0/
C                    DATA FOR ELEMENT 16        SULFUR
C
      DATA REFMN  (16)/'  S: (MNDO): M.J.S.DEWAR, C.H. REYNOLDS, J. COM
     1P. CHEM. 7, 140-143 (1986)      '/
      DATA USSM   (16)/     -72.2422810D0/
      DATA UPPM   (16)/     -56.9732070D0/
      DATA BETASM (16)/     -10.7616700D0/
      DATA BETAPM (16)/     -10.1084330D0/
      DATA ZSM    (16)/       2.3129620D0/
      DATA ZPM    (16)/       2.0091460D0/
      DATA ZDM    (16)/       1.0000000D0/
      DATA ALPM   (16)/       2.4780260D0/
      DATA EISOLM (16)/    -226.0123900D0/
      DATA DDM    (16)/       0.9189935D0/
      DATA QQM    (16)/       0.8328514D0/
      DATA AMM    (16)/       0.4733554D0/
      DATA ADM    (16)/       0.5544502D0/
      DATA AQM    (16)/       0.5585244D0/
C                    DATA FOR ELEMENT 17        CHLORINE
      DATA REFMN  (17)/' Cl: (MNDO): M.J.S.DEWAR, H.S.RZEPA, J. COMP. CH
     1EM., 4, 158, (1983)             '/
      DATA USSM   (17)/    -100.2271660D0/
      DATA UPPM   (17)/     -77.3786670D0/
      DATA BETASM (17)/     -14.2623200D0/
      DATA BETAPM (17)/     -14.2623200D0/
      DATA ZSM    (17)/       3.7846450D0/
      DATA ZPM    (17)/       2.0362630D0/
      DATA ZDM    (17)/       1.0000000D0/
      DATA ALPM   (17)/       2.5422010D0/
      DATA EISOLM (17)/    -353.1176670D0/
      DATA DDM    (17)/       0.4986870D0/
      DATA QQM    (17)/       0.8217603D0/
      DATA AMM    (17)/       0.5523705D0/
      DATA ADM    (17)/       0.8061220D0/
      DATA AQM    (17)/       0.6053435D0/
*                               DATA FOR THE POTASSIUM-LIKE SPARKLE
      DATA REFMN  (19)/'  K: (MNDO):  POTASSIUM-LIKE SPARKLE.  USE WITH
     1CARE.                           '/
*                               DATA FOR THE POTASSIUM-LIKE SPARKLE
      DATA REFAM  (19)/' K:  (AM1):  POTASSIUM-LIKE SPARKLE.   USE WITH
     1 CARE.                          '/
      DATA EHEAT(19)    / 0.0D0/
      DATA VS(19)       /10.0D0/
      DATA ALP(19)      / 1.16D0/
      DATA EISOL(19)    / 0.0D0/
      DATA AM(19)       / 0.5D0/
      DATA ALPM(19)      / 1.16D0/
      DATA EISOLM(19)    / 0.0D0/
      DATA AMM(19)       / 0.5D0/
C                    DATA FOR ELEMENT 24  CHROMIUM
      DATA REFMN  (24)/' Cr: (MNDO):  M.J.S. DEWAR, E.F. HEALY, J.J.P.
     1STEWART (IN PREPARATION)        '/
      DATA USSM   (24)/     -17.5170270D0/
      DATA UPPM  (24)/     -12.5337290D0/
      DATA UDDM  (24)/     -44.1249280D0/
      DATA BETASM (24)/      -0.1000000D0/
      DATA BETAPM (24)/      -0.1000000D0/
      DATA BETADM (24)/      -8.7766360D0/
      DATA ZSM    (24)/       1.5000000D0/
      DATA ZPM    (24)/       1.5000000D0/
      DATA ZDM    (24)/       2.8845490D0/
      DATA ALPM   (24)/       3.0683070D0/
      DATA EISOLM (24)/    -134.8187920D0/
      DATA GSSM   (24)/       6.0000000D0/
      DATA GSPM   (24)/       4.1500000D0/
      DATA GPPM   (24)/       5.0000000D0/
      DATA GP2M   (24)/       3.5000000D0/
      DATA HSPM   (24)/       1.0000000D0/
      DATA GSD    (24)/       2.8746410D0/
      DATA GPD    (24)/       3.0000000D0/
      DATA GDD    (24)/       8.8949670D0/
      DATA DDM    (24)/       1.7320508D0/
      DATA QQM    (24)/       1.4142136D0/
      DATA AMM    (24)/       0.2205072D0/
      DATA ADM    (24)/       0.2711332D0/
      DATA AQM    (24)/       0.4464656D0/
C                    DATA FOR ELEMENT 30        ZINC
      DATA REFMN  (30)/' Zn: (MNDO):  M.J.S. DEWAR, K.M. MERZ, ORGANOMET
     1ALLICS, 5, 1494-1496 (1986)     '/
C                    DATA FOR ELEMENT 30
      DATA USSM  ( 30)/     -20.8397160D0/
      DATA UPPM  ( 30)/     -19.6252240D0/
      DATA BETASM( 30)/      -1.0000000D0/
      DATA BETAPM( 30)/      -2.0000000D0/
      DATA ZSM   ( 30)/       2.0473590D0/
      DATA ZPM   ( 30)/       1.4609460D0/
      DATA ZDM   ( 30)/       1.0000000D0/
      DATA ALPM  ( 30)/       1.5064570D0/
      DATA EISOLM( 30)/     -29.8794320D0/
      DATA GSSM  ( 30)/      11.8000000D0/
      DATA GSPM  ( 30)/      11.1820180D0/
      DATA GPPM  ( 30)/      13.3000000D0/
      DATA GP2M  ( 30)/      12.9305200D0/
      DATA HSPM  ( 30)/       0.4846060D0/
      DATA DDM   ( 30)/       1.3037826D0/
      DATA QQM   ( 30)/       1.4520183D0/
      DATA AMM   ( 30)/       0.4336641D0/
      DATA ADM   ( 30)/       0.2375912D0/
      DATA AQM   ( 30)/       0.2738858D0/
C                    DATA FOR ELEMENT 32        GERMANIUM
      DATA REFMN  (32)/' Ge: (MNDO): M.J.S.DEWAR, G.L.GRADY, E.F.HEALY,O
     1RGANOMETALLICS 6 186-189, (1987)'/
      DATA USSM  ( 32)/     -33.9493670D0/
      DATA UPPM  ( 32)/     -27.4251050D0/
      DATA BETASM( 32)/      -4.5164790D0/
      DATA BETAPM( 32)/      -1.7555170D0/
      DATA ZSM   ( 32)/       1.2931800D0/
      DATA ZPM   ( 32)/       2.0205640D0/
      DATA ALPM  ( 32)/       1.9784980D0/
      DATA EISOLM( 32)/     -76.2489440D0/
      DATA GSSM  ( 32)/       9.8000000D0/
      DATA GSPM  ( 32)/       8.3000000D0/
      DATA GPPM  ( 32)/       7.3000000D0/
      DATA GP2M  ( 32)/       6.5000000D0/
      DATA HSPM  ( 32)/       1.3000000D0/
      DATA DDM   ( 32)/       1.2556091D0/
      DATA QQM   ( 32)/       1.0498655D0/
      DATA AMM   ( 32)/       0.3601617D0/
      DATA ADM   ( 32)/       0.3643722D0/
      DATA AQM   ( 32)/       0.4347337D0/
C                    DATA FOR ELEMENT 35        BROMINE
      DATA REFMN  (35)/' Br: (MNDO): M.J.S.DEWAR, E.F. HEALY, J. COMP. C
     1HEM., 4, 542, (1983)            '/
      DATA USSM   (35)/     -99.9864405D0/
      DATA UPPM   (35)/     -75.6713075D0/
      DATA BETASM (35)/      -8.9171070D0/
      DATA BETAPM (35)/      -9.9437400D0/
      DATA ZSM    (35)/       3.8543019D0/
      DATA ZPM    (35)/       2.1992091D0/
      DATA ZDM    (35)/       1.0000000D0/
      DATA ALPM   (35)/       2.4457051D0/
      DATA EISOLM (35)/    -346.6812500D0/
      DATA DDM    (35)/       0.6051074D0/
      DATA QQM    (35)/       0.9645873D0/
      DATA AMM    (35)/       0.5526068D0/
      DATA ADM    (35)/       0.7258330D0/
      DATA AQM    (35)/       0.5574589D0/
C                    DATA FOR ELEMENT 50        TIN
      DATA REFMN  (50)/' Sn: (MNDO): M.J.S.DEWAR,G.L.GRADY,J.J.P.STEWART
     1, J.AM.CHEM.SOC.,106 6771 (1984)'/
      DATA USSM  (50)/     -40.8518020D0/
      DATA UPPM   (50)/     -28.5602490D0/
      DATA BETASM (50)/      -3.2351470D0/
      DATA BETAPM (50)/      -4.2904160D0/
      DATA ZSM    (50)/       2.0803800D0/
      DATA ZPM   (50)/       1.9371060D0/
      DATA ALPM   (50)/       1.8008140D0/
      DATA EISOLM (50)/     -92.3241020D0/
      DATA GSSM   (50)/       9.8000000D0/
      DATA GSPM   (50)/       8.3000000D0/
      DATA GPPM   (50)/       7.3000000D0/
      DATA GP2M   (50)/       6.5000000D0/
      DATA HSPM   (50)/       1.3000000D0/
      DATA DDM    (50)/       1.5697766D0/
      DATA QQM    (50)/       1.3262292D0/
      DATA AMM    (50)/       0.3601617D0/
      DATA ADM    (50)/       0.3219998D0/
      DATA AQM    (50)/       0.3713827D0/
C                    DATA FOR ELEMENT 53        IODINE
      DATA REFMN  (53)/'  I: (MNDO): M.J.S.DEWAR, E.F. HEALY, J.J.P. STE
     1WART, J.COMP.CHEM., 5,358,(1984)'/
      DATA USSM   (53)/    -100.0030538D0/
      DATA UPPM   (53)/     -74.6114692D0/
      DATA BETASM (53)/      -7.4144510D0/
      DATA BETAPM (53)/      -6.1967810D0/
      DATA ZSM    (53)/       2.2729610D0/
      DATA ZPM    (53)/       2.1694980D0/
      DATA ZDM    (53)/       1.0000000D0/
      DATA ALPM   (53)/       2.2073200D0/
      DATA EISOLM (53)/    -340.5983600D0/
      DATA DDM    (53)/       1.4253233D0/
      DATA QQM    (53)/       1.1841707D0/
      DATA AMM    (53)/       0.5527541D0/
      DATA ADM    (53)/       0.4593451D0/
      DATA AQM    (53)/       0.4585376D0/
C                    DATA FOR ELEMENT 80        MERCURY
      DATA REFMN  (80)/' Hg: (MNDO): M.J.S.DEWAR,  ET. AL. ORGANOMETALLI
     1CS 4, 1964, (1985) SEE MANUAL   '/
      DATA USSM   ( 80)/     -19.8095740D0/
      DATA UPPM   ( 80)/     -13.1025300D0/
      DATA BETASM ( 80)/      -0.4045250D0/
      DATA BETAPM ( 80)/      -6.2066830D0/
      DATA ZSM    ( 80)/       2.2181840D0/
      DATA ZPM    ( 80)/       2.0650380D0/
      DATA ALPM   ( 80)/       1.3356410D0/
      DATA EISOLM ( 80)/     -28.8191480D0/
      DATA GSSM   ( 80)/      10.8000000D0/
      DATA GSPM   ( 80)/       9.3000000D0/
      DATA GPPM   ( 80)/      14.3000000D0/
      DATA GP2M   ( 80)/      13.5000000D0/
      DATA HSPM   ( 80)/       1.3000000D0/
      DATA DDM    ( 80)/       1.7378048D0/
      DATA QQM    ( 80)/       1.4608064D0/
      DATA AMM    ( 80)/       0.3969129D0/
      DATA ADM    ( 80)/       0.3047694D0/
      DATA AQM    ( 80)/       0.3483102D0/
C                    DATA FOR ELEMENT 82        LEAD
      DATA REFMN  (82)/' Pb: (MNDO): M.J.S.DEWAR, ET.AL ORGANOMETALLICS
     14 1973-1980 (1985)              '/
      DATA USSM   ( 82)/     -47.3196920D0/
      DATA UPPM   ( 82)/     -28.8475600D0/
      DATA BETASM ( 82)/      -8.0423870D0/
      DATA BETAPM ( 82)/      -3.0000000D0/
      DATA ZSM    ( 82)/       2.4982860D0/
      DATA ZPM    ( 82)/       2.0820710D0/
      DATA ALPM   ( 82)/       1.7283330D0/
      DATA EISOLM ( 82)/    -105.8345040D0/
      DATA GSSM   ( 82)/       9.8000000D0/
      DATA GSPM   ( 82)/       8.3000000D0/
      DATA GPPM   ( 82)/       7.3000000D0/
      DATA GP2M   ( 82)/       6.5000000D0/
      DATA HSPM   ( 82)/       1.3000000D0/
      DATA DDM    ( 82)/       1.5526624D0/
      DATA QQM    ( 82)/       1.4488558D0/
      DATA AMM    ( 82)/       0.3601617D0/
      DATA ADM    ( 82)/       0.3239309D0/
      DATA AQM    ( 82)/       0.3502057D0/
C
C     START OF "OLD" ELEMENTS: THESE ARE OLD PARAMETERS WHICH
C     CAN BE USED, IF DESIRED, BY SPECIFYING "<CHEMICAL SYMBOL>YEAR"
C     AS IN SI1978 OR  S1983.
C
C                    DATA FOR ELEMENT 90        SILICON
      DATA REFMN  (90)/' Si: (MNDO): M.J.S.DEWAR, M.L.MCKEE, H.S.RZEPA,
     1J. AM. CHEM. SOC., 100 3607 1978'/
      DATA USSM   (90)/     -40.5682920D0/
      DATA UPPM   (90)/     -28.0891870D0/
      DATA BETASM (90)/      -4.2562180D0/
      DATA BETAPM (90)/      -4.2562180D0/
      DATA ZSM    (90)/       1.4353060D0/
      DATA ZPM    (90)/       1.4353060D0/
      DATA ZDM    (90)/       1.0000000D0/
      DATA ALPM   (90)/       2.1961078D0/
      DATA EISOLM (90)/     -90.5399580D0/
      DATA DDM    (90)/       1.4078712D0/
      DATA QQM    (90)/       1.1658281D0/
      DATA AMM    (90)/       0.3608967D0/
      DATA ADM    (90)/       0.3441817D0/
      DATA AQM    (90)/       0.3999442D0/
      DATA HSPM(90)/1.32D00/
      DATA GP2M(90)/6.54D00/
      DATA GPPM(90)/7.31D00/
      DATA GSPM(90)/8.36D00/
      DATA GSSM(90)/9.82D00/
      DATA REFMN  (91)/'  S: (MNDO): M.J.S.DEWAR, H.S. RZEPA, M.L.MCKEE,
     1 J.AM.CHEM.SOC.100, 3607 (1978).'/
      DATA USSM   (91)/     -75.2391520D0/
      DATA UPPM   (91)/     -57.8320130D0/
      DATA BETASM (91)/     -11.1422310D0/
      DATA BETAPM (91)/     -11.1422310D0/
      DATA ZSM    (91)/       2.6135910D0/
      DATA ZPM    (91)/       2.0343930D0/
      DATA ZDM    (91)/       1.0000000D0/
      DATA ALPM   (91)/       2.4916445D0/
      DATA EISOLM (91)/    -235.4413560D0/
      DATA DDM    (91)/       0.8231596D0/
      DATA QQM    (91)/       0.8225156D0/
      DATA AMM    (91)/       0.4733554D0/
      DATA ADM    (91)/       0.5889395D0/
      DATA AQM    (91)/       0.5632724D0/
      DATA REFMN (102)/' Cb: (MNDO):  Capped Bond  (Hydrogen-like, takes
     1 on a  zero charge.)            '/
      DATA USS  M(102)/     -11.9062760D0/
      DATA BETASM(102)/-9999999.0000000D0/
      DATA ZS   M(102)/       4.0000000D0/
      DATA ZP   M(102)/       0.3000000D0/
      DATA ZD   M(102)/       0.3000000D0/
      DATA ALP  M(102)/       2.5441341D0/
      DATA EISOLM(102)/       4.0000000D0/
      DATA GSS  M(102)/      12.8480000D0/
      DATA HSP  M(102)/       0.1000000D0/
      DATA DD   M(102)/       0.0684105D0/
      DATA QQ   M(102)/       1.0540926D0/
      DATA AM   M(102)/       0.4721793D0/
      DATA AD   M(102)/       0.9262742D0/
      DATA AQ   M(102)/       0.2909059D0/
***********************************************************************
*
*    START OF AM1 PARAMETERS
*
***********************************************************************
C                    DATA FOR ELEMENT  1       AM1:   HYDROGEN
      DATA REFAM  ( 1)/'  H: (AM1): M.J.S. DEWAR ET AL, J. AM. CHEM. SOC
     1. 107 3902-3909 (1985)          '/
      DATA USS   ( 1)/     -11.3964270D0/
      DATA BETAS ( 1)/      -6.1737870D0/
      DATA ZS    ( 1)/       1.1880780D0/
      DATA ALP   ( 1)/       2.8823240D0/
      DATA EISOL ( 1)/     -11.3964270D0/
      DATA GSS   ( 1)/      12.8480000D0/
      DATA AM    ( 1)/       0.4721793D0/
      DATA AD    ( 1)/       0.4721793D0/
      DATA AQ    ( 1)/       0.4721793D0/
      DATA GUESS1( 1,1)/       0.1227960D0/
      DATA GUESS2( 1,1)/       5.0000000D0/
      DATA GUESS3( 1,1)/       1.2000000D0/
      DATA GUESS1( 1,2)/       0.0050900D0/
      DATA GUESS2( 1,2)/       5.0000000D0/
      DATA GUESS3( 1,2)/       1.8000000D0/
      DATA GUESS1( 1,3)/      -0.0183360D0/
      DATA GUESS2( 1,3)/       2.0000000D0/
      DATA GUESS3( 1,3)/       2.1000000D0/
C                    DATA FOR ELEMENT  3       AM1:   LITHIUM    *
      DATA REFAM  ( 3)/' Li: (MNDO):  TAKEN FROM MNDOC BY W.THIEL,
     1QCPE NO.438, V. 2, P.63, (1982).'/
      DATA USS   (  3)/      -5.1280000D0/
      DATA UPP   (  3)/      -2.7212000D0/
      DATA BETAS (  3)/      -1.3500400D0/
      DATA BETAP (  3)/      -1.3500400D0/
      DATA ZS    (  3)/       0.7023800D0/
      DATA ZP    (  3)/       0.7023800D0/
      DATA ALP   (  3)/       1.2501400D0/
      DATA EISOL (  3)/      -5.1280000D0/
      DATA GSS   (  3)/       7.3000000D0/
      DATA GSP   (  3)/       5.4200000D0/
      DATA GPP   (  3)/       5.0000000D0/
      DATA GP2   (  3)/       4.5200000D0/
      DATA HSP   (  3)/       0.8300000D0/
      DATA DD    (  3)/       2.0549783D0/
      DATA QQ    (  3)/       1.7437069D0/
      DATA AM    (  3)/       0.2682837D0/
      DATA AD    (  3)/       0.2269793D0/
      DATA AQ    (  3)/       0.2614581D0/
C                    DATA FOR ELEMENT  4       AM1:   BERYLLIUM  *
      DATA REFAM  ( 4)/' Be: (MNDO):  M.J.S. DEWAR, H.S. RZEPA, J. AM. C
     1HEM. SOC., 100, 777, (1978)     '/
      DATA USS   ( 4)/     -16.6023780D0/
      DATA UPP   ( 4)/     -10.7037710D0/
      DATA BETAS ( 4)/      -4.0170960D0/
      DATA BETAP ( 4)/      -4.0170960D0/
      DATA ZS    ( 4)/       1.0042100D0/
      DATA ZP    ( 4)/       1.0042100D0/
      DATA ALP   ( 4)/       1.6694340D0/
      DATA EISOL ( 4)/     -24.2047560D0/
      DATA GSS   ( 4)/       9.0000000D0/
      DATA GSP   ( 4)/       7.4300000D0/
      DATA GPP   ( 4)/       6.9700000D0/
      DATA GP2   ( 4)/       6.2200000D0/
      DATA HSP   ( 4)/       1.2800000D0/
      DATA DD    ( 4)/       1.4373245D0/
      DATA QQ    ( 4)/       1.2196103D0/
      DATA AM    ( 4)/       0.3307607D0/
      DATA AD    ( 4)/       0.3356142D0/
      DATA AQ    ( 4)/       0.3846373D0/
C                    DATA FOR ELEMENT  5       AM1:   BORON  *
      DATA REFAM  ( 5)/'  B: (AM1):  M.J.S. DEWAR, C. JIE, E. G. ZOEBISC
     1H ORGANOMETALLICS 7, 513 (1988) '/
C                    DATA FOR ELEMENT  5
      DATA USS   (  5)/     -34.4928700D0/
      DATA UPP   (  5)/     -22.6315250D0/
      DATA BETAS (  5)/      -9.5991140D0/
      DATA BETAP (  5)/      -6.2737570D0/
      DATA ZS    (  5)/       1.6117090D0/
      DATA ZP    (  5)/       1.5553850D0/
      DATA ALP   (  5)/       2.4469090D0/
      DATA EISOL (  5)/     -63.7172650D0/
      DATA GSS   (  5)/      10.5900000D0/
      DATA GSP   (  5)/       9.5600000D0/
      DATA GPP   (  5)/       8.8600000D0/
      DATA GP2   (  5)/       7.8600000D0/
      DATA HSP   (  5)/       1.8100000D0/
      DATA DD    (  5)/       0.9107622D0/
      DATA QQ    (  5)/       0.7874223D0/
      DATA AM    (  5)/       0.3891951D0/
      DATA AD    (  5)/       0.5045152D0/
      DATA AQ    (  5)/       0.5678856D0/
C                    DATA FOR ELEMENT  6       AM1:   CARBON
      DATA REFAM  ( 6)/'  C: (AM1): M.J.S. DEWAR ET AL, J. AM. CHEM. SOC
     1. 107 3902-3909 (1985)          '/
      DATA USS   ( 6)/     -52.0286580D0/
      DATA UPP   ( 6)/     -39.6142390D0/
      DATA BETAS ( 6)/     -15.7157830D0/
      DATA BETAP ( 6)/      -7.7192830D0/
      DATA ZS    ( 6)/       1.8086650D0/
      DATA ZP    ( 6)/       1.6851160D0/
      DATA ALP   ( 6)/       2.6482740D0/
      DATA EISOL ( 6)/    -120.8157940D0/
      DATA GSS   ( 6)/      12.2300000D0/
      DATA GSP   ( 6)/      11.4700000D0/
      DATA GPP   ( 6)/      11.0800000D0/
      DATA GP2   ( 6)/       9.8400000D0/
      DATA HSP   ( 6)/       2.4300000D0/
      DATA DD    ( 6)/       0.8236736D0/
      DATA QQ    ( 6)/       0.7268015D0/
      DATA AM    ( 6)/       0.4494671D0/
      DATA AD    ( 6)/       0.6082946D0/
      DATA AQ    ( 6)/       0.6423492D0/
      DATA GUESS1( 6,1)/       0.0113550D0/
      DATA GUESS2( 6,1)/       5.0000000D0/
      DATA GUESS3( 6,1)/       1.6000000D0/
      DATA GUESS1( 6,2)/       0.0459240D0/
      DATA GUESS2( 6,2)/       5.0000000D0/
      DATA GUESS3( 6,2)/       1.8500000D0/
      DATA GUESS1( 6,3)/      -0.0200610D0/
      DATA GUESS2( 6,3)/       5.0000000D0/
      DATA GUESS3( 6,3)/       2.0500000D0/
      DATA GUESS1( 6,4)/      -0.0012600D0/
      DATA GUESS2( 6,4)/       5.0000000D0/
      DATA GUESS3( 6,4)/       2.6500000D0/
C                    DATA FOR ELEMENT  7       AM1:   NITROGEN
      DATA REFAM  ( 7)/'  N: (AM1): M.J.S. DEWAR ET AL, J. AM. CHEM. SOC
     1. 107 3902-3909 (1985)          '/
      DATA USS   ( 7)/     -71.8600000D0/
      DATA UPP   ( 7)/     -57.1675810D0/
      DATA BETAS ( 7)/     -20.2991100D0/
      DATA BETAP ( 7)/     -18.2386660D0/
      DATA ZS    ( 7)/       2.3154100D0/
      DATA ZP    ( 7)/       2.1579400D0/
      DATA ALP   ( 7)/       2.9472860D0/
      DATA EISOL ( 7)/    -202.4077430D0/
      DATA GSS   ( 7)/      13.5900000D0/
      DATA GSP   ( 7)/      12.6600000D0/
      DATA GPP   ( 7)/      12.9800000D0/
      DATA GP2   ( 7)/      11.5900000D0/
      DATA HSP   ( 7)/       3.1400000D0/
      DATA DD    ( 7)/       0.6433247D0/
      DATA QQ    ( 7)/       0.5675528D0/
      DATA AM    ( 7)/       0.4994487D0/
      DATA AD    ( 7)/       0.7820840D0/
      DATA AQ    ( 7)/       0.7883498D0/
      DATA GUESS1( 7,1)/       0.0252510D0/
      DATA GUESS2( 7,1)/       5.0000000D0/
      DATA GUESS3( 7,1)/       1.5000000D0/
      DATA GUESS1( 7,2)/       0.0289530D0/
      DATA GUESS2( 7,2)/       5.0000000D0/
      DATA GUESS3( 7,2)/       2.1000000D0/
      DATA GUESS1( 7,3)/      -0.0058060D0/
      DATA GUESS2( 7,3)/       2.0000000D0/
      DATA GUESS3( 7,3)/       2.4000000D0/
C                    DATA FOR ELEMENT  8       AM1:   OXYGEN
      DATA REFAM  ( 8)/'  O: (AM1): M.J.S. DEWAR ET AL, J. AM. CHEM. SOC
     1. 107 3902-3909 (1985)          '/
      DATA USS   ( 8)/     -97.8300000D0/
      DATA UPP   ( 8)/     -78.2623800D0/
      DATA BETAS ( 8)/     -29.2727730D0/
      DATA BETAP ( 8)/     -29.2727730D0/
      DATA ZS    ( 8)/       3.1080320D0/
      DATA ZP    ( 8)/       2.5240390D0/
      DATA ALP   ( 8)/       4.4553710D0/
      DATA EISOL ( 8)/    -316.0995200D0/
      DATA GSS   ( 8)/      15.4200000D0/
      DATA GSP   ( 8)/      14.4800000D0/
      DATA GPP   ( 8)/      14.5200000D0/
      DATA GP2   ( 8)/      12.9800000D0/
      DATA HSP   ( 8)/       3.9400000D0/
      DATA DD    ( 8)/       0.4988896D0/
      DATA QQ    ( 8)/       0.4852322D0/
      DATA AM    ( 8)/       0.5667034D0/
      DATA AD    ( 8)/       0.9961066D0/
      DATA AQ    ( 8)/       0.9065223D0/
      DATA GUESS1( 8,1)/       0.2809620D0/
      DATA GUESS2( 8,1)/       5.0000000D0/
      DATA GUESS3( 8,1)/       0.8479180D0/
      DATA GUESS1( 8,2)/       0.0814300D0/
      DATA GUESS2( 8,2)/       7.0000000D0/
      DATA GUESS3( 8,2)/       1.4450710D0/
C                    DATA FOR ELEMENT  9       AM1:   FLUORINE  *
C
C THE FOLLOWING NUMBERS ARE PARAMETRIZED AT THE AM1 LEVEL, BUT
C ARE NOT YET PUBLISHED.  WHEN PUBLISHED UNCOMMENT THE FOLLOWING
C BLOCK, AND REMOVE THE PRECEDING SET OF PARAMETERS.
      DATA REFAM  ( 9)/'  F: (AM1): DEWAR AND ZOEBISCH, UNPUBLISHED WORK
     1                                '/
      DATA USS   ( 9)/    -136.1055790D0/
      DATA UPP   ( 9)/    -104.8898850D0/
      DATA BETAS ( 9)/     -69.5902770D0/
      DATA BETAP ( 9)/     -27.9223600D0/
      DATA ZS    ( 9)/       3.7700820D0/
      DATA ZP    ( 9)/       2.4946700D0/
      DATA ALP   ( 9)/       5.5178000D0/
      DATA EISOL ( 9)/    -482.2905830D0/
      DATA GSS   ( 9)/      16.9200000D0/
      DATA GSP   ( 9)/      17.2500000D0/
      DATA GPP   ( 9)/      16.7100000D0/
      DATA GP2   ( 9)/      14.9100000D0/
      DATA HSP   ( 9)/       4.8300000D0/
      DATA DD    ( 9)/       0.4145203D0/
      DATA QQ    ( 9)/       0.4909446D0/
      DATA AM    ( 9)/       0.6218302D0/
      DATA AD    ( 9)/       1.2088792D0/
      DATA AQ    ( 9)/       0.9449355D0/
      DATA GUESS1( 9,1)/       0.2420790D0/
      DATA GUESS2( 9,1)/       4.8000000D0/
      DATA GUESS3( 9,1)/       0.9300000D0/
      DATA GUESS1( 9,2)/       0.0036070D0/
      DATA GUESS2( 9,2)/       4.6000000D0/
      DATA GUESS3( 9,2)/       1.6600000D0/
*                               DATA FOR THE SODIUM-LIKE SPARKLE
      DATA REFAM  (11)/' Na: (MNDO):  SODIUM-LIKE SPARKLE.   USE WITH CA
     1RE.                             '/
C                    DATA FOR ELEMENT 13       AM1:   ALUMINUM  *
      DATA REFAM  (13)/' Al: (AM1):  M. J. S. Dewar, A. J. Holder, Organ
     1ometallics, 9, 508-511 (1990).  '/
      DATA USS   ( 13)/     -24.3535850D0/
      DATA UPP   ( 13)/     -18.3636450D0/
      DATA BETAS ( 13)/      -3.8668220D0/
      DATA BETAP ( 13)/      -2.3171460D0/
      DATA ZS    ( 13)/       1.5165930D0/
      DATA ZP    ( 13)/       1.3063470D0/
      DATA ZD    ( 13)/       1.0000000D0/
      DATA ALP   ( 13)/       1.9765860D0/
      DATA EISOL ( 13)/     -46.4208150D0/
      DATA GSS   ( 13)/       8.0900000D0/
      DATA GSP   ( 13)/       6.6300000D0/
      DATA GPP   ( 13)/       5.9800000D0/
      DATA GP2   ( 13)/       5.4000000D0/
      DATA HSP   ( 13)/       0.7000000D0/
      DATA DD    ( 13)/       1.4040443D0/
      DATA QQ    ( 13)/       1.2809154D0/
      DATA AM    ( 13)/       0.2973172D0/
      DATA AD    ( 13)/       0.2630229D0/
      DATA AQ    ( 13)/       0.3427832D0/
      DATA GUESS1( 13,1)/       0.0900000D0/
      DATA GUESS2( 13,1)/      12.3924430D0/
      DATA GUESS3( 13,1)/       2.0503940D0/
C                    DATA FOR ELEMENT 14       AM1:   SILICON  *
      DATA REFAM  (14)/' Si: (AM1): M.J.S.DEWAR, C. JIE, ORGANOMETALLICS
     1, 6, 1486-1490 (1987).          '/
      DATA USS   (14)/     -33.9536220D0/
      DATA UPP   (14)/     -28.9347490D0/
      DATA BETAS (14)/      -3.784852D0/
      DATA BETAP (14)/      -1.968123D0/
      DATA ZS    (14)/       1.830697D0/
      DATA ZP    (14)/       1.2849530D0/
      DATA ZD    (14)/       1.0000000D0/
      DATA ALP   (14)/       2.257816D0/
      DATA EISOL (14)/     -79.0017420D0/
      DATA GSS   (14)/       9.8200000D0/
      DATA GSP   (14)/       8.3600000D0/
      DATA GPP   (14)/       7.3100000D0/
      DATA GP2   (14)/       6.5400000D0/
      DATA HSP   (14)/       1.3200000D0/
      DATA DD    (14)/       1.1631107D0/
      DATA QQ    (14)/       1.3022422D0/
      DATA AM    (14)/       0.3608967D0/
      DATA AD    (14)/       0.3829813D0/
      DATA AQ    (14)/       0.3712106D0/
      DATA GUESS1(14,1)/       0.25D0/
      DATA GUESS2(14,1)/       9.000D0/
      DATA GUESS3(14,1)/       0.911453D0/
      DATA GUESS1(14,2)/       0.061513D0/
      DATA GUESS2(14,2)/       5.00D0/
      DATA GUESS3(14,2)/       1.995569D0/
      DATA GUESS1(14,3)/       0.0207890D0/
      DATA GUESS2(14,3)/       5.00D0/
      DATA GUESS3(14,3)/       2.990610D0/
C                    DATA FOR ELEMENT 15        PHOSPHORUS
      DATA REFAM  (15)/'  P: (AM1): M.J.S.DEWAR, JIE, C, THEOCHEM, 187,
     11 (1989)                        '/
      DATA USS   ( 15)/     -42.0298630D0/
      DATA UPP   ( 15)/     -34.0307090D0/
      DATA BETAS ( 15)/      -6.3537640D0/
      DATA BETAP ( 15)/      -6.5907090D0/
      DATA ZS    ( 15)/       1.9812800D0/
      DATA ZP    ( 15)/       1.8751500D0/
      DATA ZD    ( 15)/       1.0000000D0/
      DATA ALP   ( 15)/       2.4553220D0/
      DATA EISOL ( 15)/    -124.4368355D0/
      DATA GSS   ( 15)/      11.5600050D0/
      DATA GSP   ( 15)/       5.2374490D0/
      DATA GPP   ( 15)/       7.8775890D0/
      DATA GP2   ( 15)/       7.3076480D0/
      DATA HSP   ( 15)/       0.7792380D0/
      DATA DD    ( 15)/       1.0452022D0/
      DATA QQ    ( 15)/       0.8923660D0/
      DATA AM    ( 15)/       0.4248440D0/
      DATA AD    ( 15)/       0.3275319D0/
      DATA AQ    ( 15)/       0.4386854D0/
      DATA GUESS1( 15,1)/      -0.0318270D0/
      DATA GUESS2( 15,1)/       6.0000000D0/
      DATA GUESS3( 15,1)/       1.4743230D0/
      DATA GUESS1( 15,2)/       0.0184700D0/
      DATA GUESS2( 15,2)/       7.0000000D0/
      DATA GUESS3( 15,2)/       1.7793540D0/
      DATA GUESS1( 15,3)/       0.0332900D0/
      DATA GUESS2( 15,3)/       9.0000000D0/
      DATA GUESS3( 15,3)/       3.0065760D0/
C                    DATA FOR ELEMENT 16       AM1:   SULFUR  *
C
      DATA REFAM  (16)/'  S: (AM1): M.J.S.DEWAR, Y-C YUAN, THEOCHEM, IN
     1 PRESS                          '/
      DATA USS   (16)/     -56.6940560D0/
      DATA UPP   (16)/     -48.7170490D0/
      DATA BETAS (16)/      -3.9205660D0/
      DATA BETAP (16)/      -7.9052780D0/
      DATA ZS    (16)/       2.3665150D0/
      DATA ZP    (16)/       1.6672630D0/
      DATA ZD    (16)/       1.0000000D0/
      DATA ALP   (16)/       2.4616480D0/
      DATA EISOL (16)/    -191.7321930D0/
      DATA GSS   (16)/      11.7863290D0/
      DATA GSP   (16)/       8.6631270D0/
      DATA GPP   (16)/      10.0393080D0/
      DATA GP2   (16)/       7.7816880D0/
      DATA HSP   (16)/       2.5321370D0/
      DATA DD    (16)/       0.9004265D0/
      DATA QQ    (16)/       1.0036329D0/
      DATA AM    (16)/       0.4331617D0/
      DATA AD    (16)/       0.5907115D0/
      DATA AQ    (16)/       0.6454943D0/
      DATA GUESS1(16,1)/      -0.5091950D0/
      DATA GUESS2(16,1)/       4.5936910D0/
      DATA GUESS3(16,1)/       0.7706650D0/
      DATA GUESS1(16,2)/      -0.0118630D0/
      DATA GUESS2(16,2)/       5.8657310D0/
      DATA GUESS3(16,2)/       1.5033130D0/
      DATA GUESS1(16,3)/       0.0123340D0/
      DATA GUESS2(16,3)/      13.5573360D0/
      DATA GUESS3(16,3)/       2.0091730D0/
C                    DATA FOR ELEMENT 17       AM1:   CHLORINE  *
C
C THE FOLLOWING NUMBERS ARE PARAMETRIZED AT THE AM1 LEVEL, BUT
C ARE NOT YET PUBLISHED.  WHEN PUBLISHED UNCOMMENT THE FOLLOWING
C BLOCK, AND REMOVE THE PRECEDING SET OF PARAMETERS.
      DATA REFAM  (17)/' Cl: (AM1): DEWAR AND ZOEBISCH, UNPUBLISHED WORK
     1                                '/
      DATA USS   (17)/    -111.6139480D0/
      DATA UPP   (17)/     -76.6401070D0/
      DATA BETAS (17)/     -24.5946700D0/
      DATA BETAP (17)/     -14.6372160D0/
      DATA ZS    (17)/       3.6313760D0/
      DATA ZP    (17)/       2.0767990D0/
      DATA ZD    (17)/       1.0000000D0/
      DATA ALP   (17)/       2.9193680D0/
      DATA EISOL (17)/    -372.1984310D0/
      DATA GSS   (17)/      15.0300000D0/
      DATA GSP   (17)/      13.1600000D0/
      DATA GPP   (17)/      11.3000000D0/
      DATA GP2   (17)/       9.9700000D0/
      DATA HSP   (17)/       2.4200000D0/
      DATA DD    (17)/       0.5406286D0/
      DATA QQ    (17)/       0.8057208D0/
      DATA AM    (17)/       0.5523705D0/
      DATA AD    (17)/       0.7693200D0/
      DATA AQ    (17)/       0.6133369D0/
      DATA GUESS1(17,1)/       0.0942430D0/
      DATA GUESS2(17,1)/       4.0000000D0/
      DATA GUESS3(17,1)/       1.3000000D0/
      DATA GUESS1(17,2)/       0.0271680D0/
      DATA GUESS2(17,2)/       4.0000000D0/
      DATA GUESS3(17,2)/       2.1000000D0/
C                    DATA FOR ELEMENT 30        ZINC
      DATA REFAM  (30)/' Zn: (MNDO):  M.J.S. DEWAR, K.M. MERZ, ORGANOMET
     1ALLICS, 5, 1494-1496 (1986)     '/
      DATA USS   ( 30)/     -20.8397160D0/
      DATA UPP   ( 30)/     -19.6252240D0/
      DATA BETAS ( 30)/      -1.0000000D0/
      DATA BETAP ( 30)/      -2.0000000D0/
      DATA ZS    ( 30)/       2.0473590D0/
      DATA ZP    ( 30)/       1.4609460D0/
      DATA ZD    ( 30)/       1.0000000D0/
      DATA ALP   ( 30)/       1.5064570D0/
      DATA EISOL ( 30)/     -29.8794320D0/
      DATA GSS   ( 30)/      11.8000000D0/
      DATA GSP   ( 30)/      11.1820180D0/
      DATA GPP   ( 30)/      13.3000000D0/
      DATA GP2   ( 30)/      12.9305200D0/
      DATA HSP   ( 30)/       0.4846060D0/
      DATA DD    ( 30)/       1.3037826D0/
      DATA QQ    ( 30)/       1.4520183D0/
      DATA AM    ( 30)/       0.4336641D0/
      DATA AD    ( 30)/       0.2375912D0/
      DATA AQ    ( 30)/       0.2738858D0/
C                    DATA FOR ELEMENT 35       AM1:   BROMINE  *
C
C THE FOLLOWING NUMBERS ARE PARAMETRIZED AT THE AM1 LEVEL, BUT
C ARE NOT YET PUBLISHED.  WHEN PUBLISHED UNCOMMENT THE FOLLOWING
C BLOCK, AND REMOVE THE PRECEDING SET OF PARAMETERS.
      DATA REFAM  (35)/' Br: (AM1): DEWAR AND ZOEBISCH, UNPUBLISHED WORK
     1                                '/
      DATA USS   (35)/    -104.6560630D0/
      DATA UPP   (35)/     -74.9300520D0/
      DATA BETAS (35)/     -19.3998800D0/
      DATA BETAP (35)/      -8.9571950D0/
      DATA ZS    (35)/       3.0641330D0/
      DATA ZP    (35)/       2.0383330D0/
      DATA ZD    (35)/       1.0000000D0/
      DATA ALP   (35)/       2.5765460D0/
      DATA EISOL (35)/    -352.3142087D0/
      DATA GSS   (35)/      15.0364395D0/
      DATA GSP   (35)/      13.0346824D0/
      DATA GPP   (35)/      11.2763254D0/
      DATA GP2   (35)/       9.8544255D0/
      DATA HSP   (35)/       2.4558683D0/
      DATA DD    (35)/       0.8458104D0/
      DATA QQ    (35)/       1.0407133D0/
      DATA AM    (35)/       0.5526071D0/
      DATA AD    (35)/       0.6024598D0/
      DATA AQ    (35)/       0.5307555D0/
      DATA GUESS1(35,1)/       0.0666850D0/
      DATA GUESS2(35,1)/       4.0000000D0/
      DATA GUESS3(35,1)/       1.5000000D0/
      DATA GUESS1(35,2)/       0.0255680D0/
      DATA GUESS2(35,2)/       4.0000000D0/
      DATA GUESS3(35,2)/       2.3000000D0/
C                    DATA FOR ELEMENT 53       AM1:   IODINE  *
C
C THE FOLLOWING NUMBERS ARE PARAMETRIZED AT THE AM1 LEVEL, BUT
C ARE NOT YET PUBLISHED.  WHEN PUBLISHED UNCOMMENT THE FOLLOWING
C BLOCK, AND REMOVE THE PRECEDING SET OF PARAMETERS.
      DATA REFAM  (53)/'  I: (AM1): DEWAR AND ZOEBISCH, UNPUBLISHED WORK
     1                                '/
      DATA USS   (53)/    -103.5896630D0/
      DATA UPP   (53)/     -74.4299970D0/
      DATA BETAS (53)/      -8.4433270D0/
      DATA BETAP (53)/      -6.3234050D0/
      DATA ZS    (53)/       2.1028580D0/
      DATA ZP    (53)/       2.1611530D0/
      DATA ZD    (53)/       1.0000000D0/
      DATA ALP   (53)/       2.2994240D0/
      DATA EISOL (53)/    -346.8642857D0/
      DATA GSS   (53)/      15.0404486D0/
      DATA GSP   (53)/      13.0565580D0/
      DATA GPP   (53)/      11.1477837D0/
      DATA GP2   (53)/       9.9140907D0/
      DATA HSP   (53)/       2.4563820D0/
      DATA DD    (53)/       1.4878778D0/
      DATA QQ    (53)/       1.1887388D0/
      DATA AM    (53)/       0.5527544D0/
      DATA AD    (53)/       0.4497523D0/
      DATA AQ    (53)/       0.4631775D0/
      DATA GUESS1(53,1)/       0.0043610D0/
      DATA GUESS2(53,1)/       2.3000000D0/
      DATA GUESS3(53,1)/       1.8000000D0/
      DATA GUESS1(53,2)/       0.0157060D0/
      DATA GUESS2(53,2)/       3.0000000D0/
      DATA GUESS3(53,2)/       2.2400000D0/
C
C     START OF "OLD" ELEMENTS: THESE ARE OLD PARAMETERS WHICH
C     CAN BE USED, IF DESIRED, BY SPECIFYING "<CHEMICAL SYMBOL>YEAR"
C     AS IN SI1978 OR  S1983.
C
C                    DATA FOR ELEMENT 90        SILICON
      DATA REFAM  (90)/' Si: (MNDO): M.J.S.DEWAR, M.L.MCKEE, H.S.RZEPA,
     1J. AM. CHEM. SOC., 100 3607 1978'/
      DATA USS   (90)/     -40.5682920D0/
      DATA UPP   (90)/     -28.0891870D0/
      DATA BETAS (90)/      -4.2562180D0/
      DATA BETAP (90)/      -4.2562180D0/
      DATA ZS    (90)/       1.4353060D0/
      DATA ZP    (90)/       1.4353060D0/
      DATA ZD    (90)/       1.0000000D0/
      DATA ALP   (90)/       2.1961078D0/
      DATA EISOL (90)/     -90.5399580D0/
      DATA DD    (90)/       1.4078712D0/
      DATA QQ    (90)/       1.1658281D0/
      DATA AM    (90)/       0.3608967D0/
      DATA AD    (90)/       0.3441817D0/
      DATA AQ    (90)/       0.3999442D0/
      DATA HSP(90)/1.32D00/
      DATA GP2(90)/6.54D00/
      DATA GPP(90)/7.31D00/
      DATA GSP(90)/8.36D00/
      DATA GSS(90)/9.82D00/
      DATA REFAM  (91)/'  S: (MNDO): M.J.S.DEWAR, H.S. RZEPA, M.L.MCKEE,
     1 J.AM.CHEM.SOC.100, 3607 (1978).'/
      DATA USS   (91)/     -75.2391520D0/
      DATA UPP   (91)/     -57.8320130D0/
      DATA BETAS (91)/     -11.1422310D0/
      DATA BETAP (91)/     -11.1422310D0/
      DATA ZS    (91)/       2.6135910D0/
      DATA ZP    (91)/       2.0343930D0/
      DATA ZD    (91)/       1.0000000D0/
      DATA ALP   (91)/       2.4916445D0/
      DATA EISOL (91)/    -235.4413560D0/
      DATA GSS   (91)/      12.8800000D0/
      DATA GSP   (91)/      11.2600000D0/
      DATA GPP   (91)/       9.9000000D0/
      DATA GP2   (91)/       8.8300000D0/
      DATA HSP   (91)/       2.2600000D0/
      DATA DD    (91)/       0.8231596D0/
      DATA QQ    (91)/       0.8225156D0/
      DATA AM    (91)/       0.4733554D0/
      DATA AD    (91)/       0.5889395D0/
      DATA AQ    (91)/       0.5632724D0/
      DATA REFAM (102)/' Cb: (AM1):  Capped Bond  (Hydrogen-like, takes
     1 on zero charge.)               '/
C                    DATA FOR ELEMENT102
      DATA USS   (102)/     -11.9062760D0/
      DATA BETAS (102)/-9999999.0000000D0/
      DATA ZS    (102)/       4.0000000D0/
      DATA ZP    (102)/       0.3000000D0/
      DATA ZD    (102)/       0.3000000D0/
      DATA ALP   (102)/       2.5441341D0/
      DATA EISOL (102)/       4.0000000D0/
      DATA GSS   (102)/      12.8480000D0/
      DATA HSP   (102)/       0.1000000D0/
      DATA DD    (102)/       0.0684105D0/
      DATA QQ    (102)/       1.0540926D0/
      DATA AM    (102)/       0.4721793D0/
      DATA AD    (102)/       0.9262742D0/
      DATA AQ    (102)/       0.2909059D0/
C
C     START OF MNDO-PM3 PARAMETER SET
C
C                    DATA FOR ELEMENT  1        HYDROGEN
      DATA REFPM3 ( 1)/'  H: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1 (IN PRESS) (1988)              '/
      DATA USSPM3(  1)/     -13.0733210D0/
      DATA BETASP(  1)/      -5.6265120D0/
      DATA ZSPM3 (  1)/       0.9678070D0/
      DATA ALPPM3(  1)/       3.3563860D0/
      DATA EISOLP(  1)/     -13.0733210D0/
      DATA GSSPM3(  1)/      14.7942080D0/
      DATA AMPM3 (  1)/       0.5437048D0/
      DATA ADPM3 (  1)/       0.5437048D0/
      DATA AQPM3 (  1)/       0.5437048D0/
      DATA GUESP1(  1,1)/       1.1287500D0/
      DATA GUESP2(  1,1)/       5.0962820D0/
      DATA GUESP3(  1,1)/       1.5374650D0/
      DATA GUESP1(  1,2)/      -1.0603290D0/
      DATA GUESP2(  1,2)/       6.0037880D0/
      DATA GUESP3(  1,2)/       1.5701890D0/
C                    DATA FOR ELEMENT  4        BERYLLIUM
      DATA REFPM3( 4)/ ' Be: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1(ACCEPTED)                      '/
      DATA USSPM3(  4)/     -17.2647520D0/
      DATA UPPPM3(  4)/     -11.3042430D0/
      DATA BETASP(  4)/      -3.9620530D0/
      DATA BETAPP(  4)/      -2.7806840D0/
      DATA ZSPM3 (  4)/       0.8774390D0/
      DATA ZPPM3 (  4)/       1.5087550D0/
      DATA ALPPM3(  4)/       1.5935360D0/
      DATA EISOLP(  4)/     -25.5166530D0/
      DATA GSSPM3(  4)/       9.0128510D0/
      DATA GSPPM3(  4)/       6.5761990D0/
      DATA GPPPM3(  4)/       6.0571820D0/
      DATA GP2PM3(  4)/       9.0052190D0/
      DATA HSPPM3(  4)/       0.5446790D0/
      DATA DDPM3 (  4)/       1.0090531D0/
      DATA QQPM3 (  4)/       0.8117586D0/
      DATA AMPM3 (  4)/       0.3312330D0/
      DATA ADPM3 (  4)/       0.2908996D0/
      DATA AQPM3 (  4)/       0.3530008D0/
      DATA GUESP1(  4,1)/       1.6315720D0/
      DATA GUESP2(  4,1)/       2.6729620D0/
      DATA GUESP3(  4,1)/       1.7916860D0/
      DATA GUESP1(  4,2)/      -2.1109590D0/
      DATA GUESP2(  4,2)/       1.9685940D0/
      DATA GUESP3(  4,2)/       1.7558710D0/
C                    DATA FOR ELEMENT  6      CARBON
      DATA REFPM3 ( 6)/'  C: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1 (IN PRESS) (1988)              '/
      DATA USSPM3(  6)/     -47.2703200D0/
      DATA UPPPM3(  6)/     -36.2669180D0/
      DATA BETASP(  6)/     -11.9100150D0/
      DATA BETAPP(  6)/      -9.8027550D0/
      DATA ZSPM3 (  6)/       1.5650850D0/
      DATA ZPPM3 (  6)/       1.8423450D0/
      DATA ALPPM3(  6)/       2.7078070D0/
      DATA EISOLP(  6)/    -111.2299170D0/
      DATA GSSPM3(  6)/      11.2007080D0/
      DATA GSPPM3(  6)/      10.2650270D0/
      DATA GPPPM3(  6)/      10.7962920D0/
      DATA GP2PM3(  6)/       9.0425660D0/
      DATA HSPPM3(  6)/       2.2909800D0/
      DATA DDPM3 (  6)/       0.8332396D0/
      DATA QQPM3 (  6)/       0.6647750D0/
      DATA AMPM3 (  6)/       0.4116394D0/
      DATA ADPM3 (  6)/       0.5885862D0/
      DATA AQPM3 (  6)/       0.7647667D0/
      DATA GUESP1(  6,1)/       0.0501070D0/
      DATA GUESP2(  6,1)/       6.0031650D0/
      DATA GUESP3(  6,1)/       1.6422140D0/
      DATA GUESP1(  6,2)/       0.0507330D0/
      DATA GUESP2(  6,2)/       6.0029790D0/
      DATA GUESP3(  6,2)/       0.8924880D0/
C                    DATA FOR ELEMENT  7      NITROGEN
      DATA REFPM3 ( 7)/'  N: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1 (IN PRESS) (1988)              '/
      DATA USSPM3(  7)/     -49.3356720D0/
      DATA UPPPM3(  7)/     -47.5097360D0/
      DATA BETASP(  7)/     -14.0625210D0/
      DATA BETAPP(  7)/     -20.0438480D0/
      DATA ZSPM3 (  7)/       2.0280940D0/
      DATA ZPPM3 (  7)/       2.3137280D0/
      DATA ALPPM3(  7)/       2.8305450D0/
      DATA EISOLP(  7)/    -157.6137755D0/
      DATA GSSPM3(  7)/      11.9047870D0/
      DATA GSPPM3(  7)/       7.3485650D0/
      DATA GPPPM3(  7)/      11.7546720D0/
      DATA GP2PM3(  7)/      10.8072770D0/
      DATA HSPPM3(  7)/       1.1367130D0/
      DATA DDPM3 (  7)/       0.6577006D0/
      DATA QQPM3 (  7)/       0.5293383D0/
      DATA AMPM3 (  7)/       0.4375151D0/
      DATA ADPM3 (  7)/       0.5030995D0/
      DATA AQPM3 (  7)/       0.7364933D0/
      DATA GUESP1(  7,1)/       1.5016740D0/
      DATA GUESP2(  7,1)/       5.9011480D0/
      DATA GUESP3(  7,1)/       1.7107400D0/
      DATA GUESP1(  7,2)/      -1.5057720D0/
      DATA GUESP2(  7,2)/       6.0046580D0/
      DATA GUESP3(  7,2)/       1.7161490D0/
C                    DATA FOR ELEMENT  8      OXYGEN
      DATA REFPM3 ( 8)/'  O: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1 (IN PRESS) (1988)              '/
      DATA USSPM3(  8)/     -86.9930020D0/
      DATA UPPPM3(  8)/     -71.8795800D0/
      DATA BETASP(  8)/     -45.2026510D0/
      DATA BETAPP(  8)/     -24.7525150D0/
      DATA ZSPM3 (  8)/       3.7965440D0/
      DATA ZPPM3 (  8)/       2.3894020D0/
      DATA ALPPM3(  8)/       3.2171020D0/
      DATA EISOLP(  8)/    -289.3422065D0/
      DATA GSSPM3(  8)/      15.7557600D0/
      DATA GSPPM3(  8)/      10.6211600D0/
      DATA GPPPM3(  8)/      13.6540160D0/
      DATA GP2PM3(  8)/      12.4060950D0/
      DATA HSPPM3(  8)/       0.5938830D0/
      DATA DDPM3 (  8)/       0.4086173D0/
      DATA QQPM3 (  8)/       0.5125738D0/
      DATA AMPM3 (  8)/       0.5790430D0/
      DATA ADPM3 (  8)/       0.5299630D0/
      DATA AQPM3 (  8)/       0.8179630D0/
      DATA GUESP1(  8,1)/      -1.1311280D0/
      DATA GUESP2(  8,1)/       6.0024770D0/
      DATA GUESP3(  8,1)/       1.6073110D0/
      DATA GUESP1(  8,2)/       1.1378910D0/
      DATA GUESP2(  8,2)/       5.9505120D0/
      DATA GUESP3(  8,2)/       1.5983950D0/
C                    DATA FOR ELEMENT  9      FLUORINE
      DATA REFPM3 ( 9)/'  F: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1 (IN PRESS) (1988)              '/
      DATA USSPM3(  9)/    -110.4353030D0/
      DATA UPPPM3(  9)/    -105.6850470D0/
      DATA BETASP(  9)/     -48.4059390D0/
      DATA BETAPP(  9)/     -27.7446600D0/
      DATA ZSPM3 (  9)/       4.7085550D0/
      DATA ZPPM3 (  9)/       2.4911780D0/
      DATA ALPPM3(  9)/       3.3589210D0/
      DATA EISOLP(  9)/    -437.5171690D0/
      DATA GSSPM3(  9)/      10.4966670D0/
      DATA GSPPM3(  9)/      16.0736890D0/
      DATA GPPPM3(  9)/      14.8172560D0/
      DATA GP2PM3(  9)/      14.4183930D0/
      DATA HSPPM3(  9)/       0.7277630D0/
      DATA DDPM3 (  9)/       0.3125302D0/
      DATA QQPM3 (  9)/       0.4916328D0/
      DATA AMPM3 (  9)/       0.3857650D0/
      DATA ADPM3 (  9)/       0.6768503D0/
      DATA AQPM3 (  9)/       0.6120047D0/
      DATA GUESP1(  9,1)/      -0.0121660D0/
      DATA GUESP2(  9,1)/       6.0235740D0/
      DATA GUESP3(  9,1)/       1.8568590D0/
      DATA GUESP1(  9,2)/      -0.0028520D0/
      DATA GUESP2(  9,2)/       6.0037170D0/
      DATA GUESP3(  9,2)/       2.6361580D0/
C                    DATA FOR ELEMENT 12        MAGNESIUM
      DATA REFPM3(12)/ ' Mg: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1(ACCEPTED)                      '/
      DATA USSPM3( 12)/     -14.6236880D0/
      DATA UPPPM3( 12)/     -14.1734600D0/
      DATA BETASP( 12)/      -2.0716910D0/
      DATA BETAPP( 12)/      -0.5695810D0/
      DATA ZSPM3 ( 12)/       0.6985520D0/
      DATA ZPPM3 ( 12)/       1.4834530D0/
      DATA ALPPM3( 12)/       1.3291470D0/
      DATA EISOLP( 12)/     -22.5530760D0/
      DATA GSSPM3( 12)/       6.6943000D0/
      DATA GSPPM3( 12)/       6.7939950D0/
      DATA GPPPM3( 12)/       6.9104460D0/
      DATA GP2PM3( 12)/       7.0908230D0/
      DATA HSPPM3( 12)/       0.5433000D0/
      DATA DDPM3 ( 12)/       1.1403950D0/
      DATA QQPM3 ( 12)/       1.1279899D0/
      DATA AMPM3 ( 12)/       0.2460235D0/
      DATA ADPM3 ( 12)/       0.2695751D0/
      DATA AQPM3 ( 12)/       0.2767522D0/
      DATA GUESP1( 12,1)/       2.1170500D0/
      DATA GUESP2( 12,1)/       6.0094770D0/
      DATA GUESP3( 12,1)/       2.0844060D0/
      DATA GUESP1( 12,2)/      -2.5477670D0/
      DATA GUESP2( 12,2)/       4.3953700D0/
      DATA GUESP3( 12,2)/       2.0636740D0/
C                    DATA FOR ELEMENT 13      ALUMINUM
      DATA REFPM3 (13)/' Al: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1 (IN PRESS) (1988)              '/
      DATA USSPM3( 13)/     -24.8454040D0/
      DATA UPPPM3( 13)/     -22.2641590D0/
      DATA BETASP( 13)/      -0.5943010D0/
      DATA BETAPP( 13)/      -0.9565500D0/
      DATA ZSPM3 ( 13)/       1.7028880D0/
      DATA ZPPM3 ( 13)/       1.0736290D0/
      DATA ZDPM3 ( 13)/       1.0000000D0/
      DATA ALPPM3( 13)/       1.5217030D0/
      DATA EISOLP( 13)/     -46.8647630D0/
      DATA GSSPM3( 13)/       5.7767370D0/
      DATA GSPPM3( 13)/      11.6598560D0/
      DATA GPPPM3( 13)/       6.3477900D0/
      DATA GP2PM3( 13)/       6.1210770D0/
      DATA HSPPM3( 13)/       4.0062450D0/
      DATA DDPM3 ( 13)/       1.2102799D0/
      DATA QQPM3 ( 13)/       1.5585645D0/
      DATA AMPM3 ( 13)/       0.2123020D0/
      DATA ADPM3 ( 13)/       0.6418584D0/
      DATA AQPM3 ( 13)/       0.2262838D0/
      DATA GUESP1( 13,1)/      -0.4730900D0/
      DATA GUESP2( 13,1)/       1.9158250D0/
      DATA GUESP3( 13,1)/       1.4517280D0/
      DATA GUESP1( 13,2)/      -0.1540510D0/
      DATA GUESP2( 13,2)/       6.0050860D0/
      DATA GUESP3( 13,2)/       2.5199970D0/
C                    DATA FOR ELEMENT 14      SILICON
      DATA REFPM3 (14)/' Si: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1 (IN PRESS) (1988)              '/
      DATA USSPM3( 14)/     -26.7634830D0/
      DATA UPPPM3( 14)/     -22.8136350D0/
      DATA BETASP( 14)/      -2.8621450D0/
      DATA BETAPP( 14)/      -3.9331480D0/
      DATA ZSPM3 ( 14)/       1.6350750D0/
      DATA ZPPM3 ( 14)/       1.3130880D0/
      DATA ZDPM3 ( 14)/       1.0000000D0/
      DATA ALPPM3( 14)/       2.1358090D0/
      DATA EISOLP( 14)/     -67.7882140D0/
      DATA GSSPM3( 14)/       5.0471960D0/
      DATA GSPPM3( 14)/       5.9490570D0/
      DATA GPPPM3( 14)/       6.7593670D0/
      DATA GP2PM3( 14)/       5.1612970D0/
      DATA HSPPM3( 14)/       0.9198320D0/
      DATA DDPM3 ( 14)/       1.3144550D0/
      DATA QQPM3 ( 14)/       1.2743396D0/
      DATA AMPM3 ( 14)/       0.1854905D0/
      DATA ADPM3 ( 14)/       0.3060715D0/
      DATA AQPM3 ( 14)/       0.4877432D0/
      DATA GUESP1( 14,1)/      -0.3906000D0/
      DATA GUESP2( 14,1)/       6.0000540D0/
      DATA GUESP3( 14,1)/       0.6322620D0/
      DATA GUESP1( 14,2)/       0.0572590D0/
      DATA GUESP2( 14,2)/       6.0071830D0/
      DATA GUESP3( 14,2)/       2.0199870D0/
      DATA GUESP1( 14,3)/       0.0207890D0/
      DATA GUESP2( 14,3)/       5.0000000D0/
      DATA GUESP3( 14,3)/       2.9906100D0/
C                    DATA FOR ELEMENT 15      PHOSPHORUS
      DATA REFPM3 (15)/'  P: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1 (IN PRESS) (1988)              '/
      DATA USSPM3( 15)/     -40.4130960D0/
      DATA UPPPM3( 15)/     -29.5930520D0/
      DATA BETASP( 15)/     -12.6158790D0/
      DATA BETAPP( 15)/      -4.1600400D0/
      DATA ZSPM3 ( 15)/       2.0175630D0/
      DATA ZPPM3 ( 15)/       1.5047320D0/
      DATA ZDPM3 ( 15)/       1.0000000D0/
      DATA ALPPM3( 15)/       1.9405340D0/
      DATA EISOLP( 15)/    -117.9591740D0/
      DATA GSSPM3( 15)/       7.8016150D0/
      DATA GSPPM3( 15)/       5.1869490D0/
      DATA GPPPM3( 15)/       6.6184780D0/
      DATA GP2PM3( 15)/       6.0620020D0/
      DATA HSPPM3( 15)/       1.5428090D0/
      DATA DDPM3 ( 15)/       1.0644947D0/
      DATA QQPM3 ( 15)/       1.1120386D0/
      DATA AMPM3 ( 15)/       0.2867187D0/
      DATA ADPM3 ( 15)/       0.4309446D0/
      DATA AQPM3 ( 15)/       0.3732517D0/
      DATA GUESP1( 15,1)/      -0.6114210D0/
      DATA GUESP2( 15,1)/       1.9972720D0/
      DATA GUESP3( 15,1)/       0.7946240D0/
      DATA GUESP1( 15,2)/      -0.0939350D0/
      DATA GUESP2( 15,2)/       1.9983600D0/
      DATA GUESP3( 15,2)/       1.9106770D0/
C                    DATA FOR ELEMENT 16      SULFUR
      DATA REFPM3 (16)/'  S: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1 (IN PRESS) (1988)              '/
      DATA USSPM3( 16)/     -49.8953710D0/
      DATA UPPPM3( 16)/     -44.3925830D0/
      DATA BETASP( 16)/      -8.8274650D0/
      DATA BETAPP( 16)/      -8.0914150D0/
      DATA ZSPM3 ( 16)/       1.8911850D0/
      DATA ZPPM3 ( 16)/       1.6589720D0/
      DATA ZDPM3 ( 16)/       1.0000000D0/
      DATA ALPPM3( 16)/       2.2697060D0/
      DATA EISOLP( 16)/    -183.4537395D0/
      DATA GSSPM3( 16)/       8.9646670D0/
      DATA GSPPM3( 16)/       6.7859360D0/
      DATA GPPPM3( 16)/       9.9681640D0/
      DATA GP2PM3( 16)/       7.9702470D0/
      DATA HSPPM3( 16)/       4.0418360D0/
      DATA DDPM3 ( 16)/       1.1214313D0/
      DATA QQPM3 ( 16)/       1.0086488D0/
      DATA AMPM3 ( 16)/       0.3294622D0/
      DATA ADPM3 ( 16)/       0.6679118D0/
      DATA AQPM3 ( 16)/       0.6137472D0/
      DATA GUESP1( 16,1)/      -0.3991910D0/
      DATA GUESP2( 16,1)/       6.0006690D0/
      DATA GUESP3( 16,1)/       0.9621230D0/
      DATA GUESP1( 16,2)/      -0.0548990D0/
      DATA GUESP2( 16,2)/       6.0018450D0/
      DATA GUESP3( 16,2)/       1.5799440D0/
C                    DATA FOR ELEMENT 17      CHLORINE
      DATA REFPM3 (17)/' Cl: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1 (IN PRESS) (1988)              '/
      DATA USSPM3( 17)/    -100.6267470D0/
      DATA UPPPM3( 17)/     -53.6143960D0/
      DATA BETASP( 17)/     -27.5285600D0/
      DATA BETAPP( 17)/     -11.5939220D0/
      DATA ZSPM3 ( 17)/       2.2462100D0/
      DATA ZPPM3 ( 17)/       2.1510100D0/
      DATA ZDPM3 ( 17)/       1.0000000D0/
      DATA ALPPM3( 17)/       2.5172960D0/
      DATA EISOLP( 17)/    -315.1949480D0/
      DATA GSSPM3( 17)/      16.0136010D0/
      DATA GSPPM3( 17)/       8.0481150D0/
      DATA GPPPM3( 17)/       7.5222150D0/
      DATA GP2PM3( 17)/       7.5041540D0/
      DATA HSPPM3( 17)/       3.4811530D0/
      DATA DDPM3 ( 17)/       0.9175856D0/
      DATA QQPM3 ( 17)/       0.7779230D0/
      DATA AMPM3 ( 17)/       0.5885190D0/
      DATA ADPM3 ( 17)/       0.6814522D0/
      DATA AQPM3 ( 17)/       0.3643694D0/
      DATA GUESP1( 17,1)/      -0.1715910D0/
      DATA GUESP2( 17,1)/       6.0008020D0/
      DATA GUESP3( 17,1)/       1.0875020D0/
      DATA GUESP1( 17,2)/      -0.0134580D0/
      DATA GUESP2( 17,2)/       1.9666180D0/
      DATA GUESP3( 17,2)/       2.2928910D0/
C                    DATA FOR ELEMENT 30        ZINC
      DATA REFPM3(30)/ ' Zn: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1(ACCEPTED)                      '/
      DATA USSPM3( 30)/     -18.5321980D0/
      DATA UPPPM3( 30)/     -11.0474090D0/
      DATA BETASP( 30)/      -0.7155780D0/
      DATA BETAPP( 30)/      -6.3518640D0/
      DATA ZSPM3 ( 30)/       1.8199890D0/
      DATA ZPPM3 ( 30)/       1.5069220D0/
      DATA ZDPM3 ( 30)/       1.0000000D0/
      DATA ALPPM3( 30)/       1.3501260D0/
      DATA EISOLP( 30)/     -27.3872000D0/
      DATA GSSPM3( 30)/       9.6771960D0/
      DATA GSPPM3( 30)/       7.7362040D0/
      DATA GPPPM3( 30)/       4.9801740D0/
      DATA GP2PM3( 30)/       4.6696560D0/
      DATA HSPPM3( 30)/       0.6004130D0/
      DATA DDPM3 ( 30)/       1.5005758D0/
      DATA QQPM3 ( 30)/       1.4077174D0/
      DATA AMPM3 ( 30)/       0.3556485D0/
      DATA ADPM3 ( 30)/       0.2375689D0/
      DATA AQPM3 ( 30)/       0.2661069D0/
      DATA GUESP1( 30,1)/      -0.1112340D0/
      DATA GUESP2( 30,1)/       6.0014780D0/
      DATA GUESP3( 30,1)/       1.5160320D0/
      DATA GUESP1( 30,2)/      -0.1323700D0/
      DATA GUESP2( 30,2)/       1.9958390D0/
      DATA GUESP3( 30,2)/       2.5196420D0/
C                    DATA FOR ELEMENT 31        GALLIUM
      DATA REFPM3(31)/ ' Ga: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1(ACCEPTED)                      '/
      DATA USSPM3( 31)/     -29.8555930D0/
      DATA UPPPM3( 31)/     -21.8753710D0/
      DATA BETASP( 31)/      -4.9456180D0/
      DATA BETAPP( 31)/      -0.4070530D0/
      DATA ZSPM3 ( 31)/       1.8470400D0/
      DATA ZPPM3 ( 31)/       0.8394110D0/
      DATA ALPPM3( 31)/       1.6051150D0/
      DATA EISOLP( 31)/     -57.3280250D0/
      DATA GSSPM3( 31)/       8.4585540D0/
      DATA GSPPM3( 31)/       8.9256190D0/
      DATA GPPPM3( 31)/       5.0868550D0/
      DATA GP2PM3( 31)/       4.9830450D0/
      DATA HSPPM3( 31)/       2.0512600D0/
      DATA DDPM3 ( 31)/       0.9776692D0/
      DATA QQPM3 ( 31)/       2.5271534D0/
      DATA AMPM3 ( 31)/       0.3108620D0/
      DATA ADPM3 ( 31)/       0.5129360D0/
      DATA AQPM3 ( 31)/       0.1546208D0/
      DATA GUESP1( 31,1)/      -0.5601790D0/
      DATA GUESP2( 31,1)/       5.6232730D0/
      DATA GUESP3( 31,1)/       1.5317800D0/
      DATA GUESP1( 31,2)/      -0.2727310D0/
      DATA GUESP2( 31,2)/       1.9918430D0/
      DATA GUESP3( 31,2)/       2.1838640D0/
C                    DATA FOR ELEMENT 32        GERMANIUM
      DATA REFPM3(32)/ ' Ge: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1(ACCEPTED)                      '/
      DATA USSPM3( 32)/     -35.4671955D0/
      DATA UPPPM3( 32)/     -31.5863583D0/
      DATA BETASP( 32)/      -5.3250024D0/
      DATA BETAPP( 32)/      -2.2501567D0/
      DATA ZSPM3 ( 32)/       2.2373526D0/
      DATA ZPPM3 ( 32)/       1.5924319D0/
      DATA ALPPM3( 32)/       1.9723370D0/
      DATA EISOLP( 32)/     -84.0156006D0/
      DATA GSSPM3( 32)/       5.3769635D0/
      DATA GSPPM3( 32)/      10.2095293D0/
      DATA GPPPM3( 32)/       7.6718647D0/
      DATA GP2PM3( 32)/       6.9242663D0/
      DATA HSPPM3( 32)/       1.3370204D0/
      DATA DDPM3 ( 32)/       1.1920304D0/
      DATA QQPM3 ( 32)/       1.3321263D0/
      DATA AMPM3 ( 32)/       0.1976098D0/
      DATA ADPM3 ( 32)/       0.3798182D0/
      DATA AQPM3 ( 32)/       0.3620669D0/
      DATA GUESP1( 32,1)/       0.9631726D0/
      DATA GUESP2( 32,1)/       6.0120134D0/
      DATA GUESP3( 32,1)/       2.1633655D0/
      DATA GUESP1( 32,2)/      -0.9593891D0/
      DATA GUESP2( 32,2)/       5.7491802D0/
      DATA GUESP3( 32,2)/       2.1693724D0/
C                    DATA FOR ELEMENT 33        ARSENIC
      DATA REFPM3(33)/ ' As: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1(ACCEPTED)                      '/
      DATA USSPM3( 33)/     -38.5074240D0/
      DATA UPPPM3( 33)/     -35.1524150D0/
      DATA BETASP( 33)/      -8.2321650D0/
      DATA BETAPP( 33)/      -5.0173860D0/
      DATA ZSPM3 ( 33)/       2.6361770D0/
      DATA ZPPM3 ( 33)/       1.7038890D0/
      DATA ALPPM3( 33)/       1.7944770D0/
      DATA EISOLP( 33)/    -122.6326140D0/
      DATA GSSPM3( 33)/       8.7890010D0/
      DATA GSPPM3( 33)/       5.3979830D0/
      DATA GPPPM3( 33)/       8.2872500D0/
      DATA GP2PM3( 33)/       8.2103460D0/
      DATA HSPPM3( 33)/       1.9510340D0/
      DATA DDPM3 ( 33)/       0.9679655D0/
      DATA QQPM3 ( 33)/       1.2449874D0/
      DATA AMPM3 ( 33)/       0.3230063D0/
      DATA ADPM3 ( 33)/       0.5042239D0/
      DATA AQPM3 ( 33)/       0.2574219D0/
      DATA GUESP1( 33,1)/      -0.4600950D0/
      DATA GUESP2( 33,1)/       1.9831150D0/
      DATA GUESP3( 33,1)/       1.0867930D0/
      DATA GUESP1( 33,2)/      -0.0889960D0/
      DATA GUESP2( 33,2)/       1.9929440D0/
      DATA GUESP3( 33,2)/       2.1400580D0/
C                    DATA FOR ELEMENT 34        SELENIUM
      DATA REFPM3(34)/ ' Se: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1(ACCEPTED)                      '/
      DATA USSPM3( 34)/     -55.3781350D0/
      DATA UPPPM3( 34)/     -49.8230760D0/
      DATA BETASP( 34)/      -6.1578220D0/
      DATA BETAPP( 34)/      -5.4930390D0/
      DATA ZSPM3 ( 34)/       2.8280510D0/
      DATA ZPPM3 ( 34)/       1.7325360D0/
      DATA ALPPM3( 34)/       3.0439570D0/
      DATA EISOLP( 34)/    -192.7748115D0/
      DATA GSSPM3( 34)/       7.4325910D0/
      DATA GSPPM3( 34)/      10.0604610D0/
      DATA GPPPM3( 34)/       9.5683260D0/
      DATA GP2PM3( 34)/       7.7242890D0/
      DATA HSPPM3( 34)/       4.0165580D0/
      DATA DDPM3 ( 34)/       0.8719813D0/
      DATA QQPM3 ( 34)/       1.2244019D0/
      DATA AMPM3 ( 34)/       0.2731566D0/
      DATA ADPM3 ( 34)/       0.7509697D0/
      DATA AQPM3 ( 34)/       0.5283737D0/
      DATA GUESP1( 34,1)/       0.0478730D0/
      DATA GUESP2( 34,1)/       6.0074000D0/
      DATA GUESP3( 34,1)/       2.0817170D0/
      DATA GUESP1( 34,2)/       0.1147200D0/
      DATA GUESP2( 34,2)/       6.0086720D0/
      DATA GUESP3( 34,2)/       1.5164230D0/
C                    DATA FOR ELEMENT 35      BROMINE
      DATA REFPM3 (35)/' Br: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1 (IN PRESS) (1988)              '/
      DATA USSPM3( 35)/    -116.6193110D0/
      DATA UPPPM3( 35)/     -74.2271290D0/
      DATA BETASP( 35)/     -31.1713420D0/
      DATA BETAPP( 35)/      -6.8140130D0/
      DATA ZSPM3 ( 35)/       5.3484570D0/
      DATA ZPPM3 ( 35)/       2.1275900D0/
      DATA ZDPM3 ( 35)/       1.0000000D0/
      DATA ALPPM3( 35)/       2.5118420D0/
      DATA EISOLP( 35)/    -352.5398970D0/
      DATA GSSPM3( 35)/      15.9434250D0/
      DATA GSPPM3( 35)/      16.0616800D0/
      DATA GPPPM3( 35)/       8.2827630D0/
      DATA GP2PM3( 35)/       7.8168490D0/
      DATA HSPPM3( 35)/       0.5788690D0/
      DATA DDPM3 ( 35)/       0.2759025D0/
      DATA QQPM3 ( 35)/       0.9970532D0/
      DATA AMPM3 ( 35)/       0.5859399D0/
      DATA ADPM3 ( 35)/       0.6755383D0/
      DATA AQPM3 ( 35)/       0.3823719D0/
      DATA GUESP1( 35,1)/       0.9604580D0/
      DATA GUESP2( 35,1)/       5.9765080D0/
      DATA GUESP3( 35,1)/       2.3216540D0/
      DATA GUESP1( 35,2)/      -0.9549160D0/
      DATA GUESP2( 35,2)/       5.9447030D0/
      DATA GUESP3( 35,2)/       2.3281420D0/
C                    DATA FOR ELEMENT 48        CADMIUM
      DATA REFPM3(48)/ ' Cd: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1(ACCEPTED)                      '/
      DATA USSPM3( 48)/     -15.8285840D0/
      DATA UPPPM3( 48)/       8.7497950D0/
      DATA BETASP( 48)/      -8.5819440D0/
      DATA BETAPP( 48)/      -0.6010340D0/
      DATA ZSPM3 ( 48)/       1.6793510D0/
      DATA ZPPM3 ( 48)/       2.0664120D0/
      DATA ALPPM3( 48)/       1.5253820D0/
      DATA EISOLP( 48)/     -22.4502080D0/
      DATA GSSPM3( 48)/       9.2069600D0/
      DATA GSPPM3( 48)/       8.2315390D0/
      DATA GPPPM3( 48)/       4.9481040D0/
      DATA GP2PM3( 48)/       4.6696560D0/
      DATA HSPPM3( 48)/       1.6562340D0/
      DATA DDPM3 ( 48)/       1.5982681D0/
      DATA QQPM3 ( 48)/       1.2432402D0/
      DATA AMPM3 ( 48)/       0.3383668D0/
      DATA ADPM3 ( 48)/       0.3570290D0/
      DATA AQPM3 ( 48)/       0.2820582D0/
C                    DATA FOR ELEMENT 49        INDIUM
      DATA REFPM3(49)/ ' In: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1(ACCEPTED)                      '/
      DATA USSPM3( 49)/     -26.1762050D0/
      DATA UPPPM3( 49)/     -20.0058220D0/
      DATA BETASP( 49)/      -2.9933190D0/
      DATA BETAPP( 49)/      -1.8289080D0/
      DATA ZSPM3 ( 49)/       2.0161160D0/
      DATA ZPPM3 ( 49)/       1.4453500D0/
      DATA ALPPM3( 49)/       1.4183850D0/
      DATA EISOLP( 49)/     -51.9750470D0/
      DATA GSSPM3( 49)/       6.5549000D0/
      DATA GSPPM3( 49)/       8.2298730D0/
      DATA GPPPM3( 49)/       6.2992690D0/
      DATA GP2PM3( 49)/       4.9842110D0/
      DATA HSPPM3( 49)/       2.6314610D0/
      DATA DDPM3 ( 49)/       1.5766241D0/
      DATA QQPM3 ( 49)/       1.7774563D0/
      DATA AMPM3 ( 49)/       0.2409004D0/
      DATA ADPM3 ( 49)/       0.4532655D0/
      DATA AQPM3 ( 49)/       0.3689812D0/
      DATA GUESP1( 49,1)/      -0.3431380D0/
      DATA GUESP2( 49,1)/       1.9940340D0/
      DATA GUESP3( 49,1)/       1.6255160D0/
      DATA GUESP1( 49,2)/      -0.1095320D0/
      DATA GUESP2( 49,2)/       5.6832170D0/
      DATA GUESP3( 49,2)/       2.8670090D0/
C                    DATA FOR ELEMENT 50        TIN
      DATA REFPM3(50)/ ' Sn: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1(ACCEPTED)                      '/
      DATA USSPM3( 50)/     -34.5501920D0/
      DATA UPPPM3( 50)/     -25.8944190D0/
      DATA BETASP( 50)/      -2.7858020D0/
      DATA BETAPP( 50)/      -2.0059990D0/
      DATA ZSPM3 ( 50)/       2.3733280D0/
      DATA ZPPM3 ( 50)/       1.6382330D0/
      DATA ALPPM3( 50)/       1.6996500D0/
      DATA EISOLP( 50)/     -78.8877790D0/
      DATA GSSPM3( 50)/      10.1900330D0/
      DATA GSPPM3( 50)/       7.2353270D0/
      DATA GPPPM3( 50)/       5.6738100D0/
      DATA GP2PM3( 50)/       5.1822140D0/
      DATA HSPPM3( 50)/       1.0331570D0/
      DATA DDPM3 ( 50)/       1.3120038D0/
      DATA QQPM3 ( 50)/       1.5681814D0/
      DATA AMPM3 ( 50)/       0.3744959D0/
      DATA ADPM3 ( 50)/       0.3218163D0/
      DATA AQPM3 ( 50)/       0.2832529D0/
      DATA GUESP1( 50,1)/      -0.1503530D0/
      DATA GUESP2( 50,1)/       6.0056940D0/
      DATA GUESP3( 50,1)/       1.7046420D0/
      DATA GUESP1( 50,2)/      -0.0444170D0/
      DATA GUESP2( 50,2)/       2.2573810D0/
      DATA GUESP3( 50,2)/       2.4698690D0/
C                    DATA FOR ELEMENT 51        ANTIMONY
      DATA REFPM3(51)/ ' Sb: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1(ACCEPTED)                      '/
      DATA USSPM3( 51)/     -56.4321960D0/
      DATA UPPPM3( 51)/     -29.4349540D0/
      DATA BETASP( 51)/     -14.7942170D0/
      DATA BETAPP( 51)/      -2.8179480D0/
      DATA ZSPM3 ( 51)/       2.3430390D0/
      DATA ZPPM3 ( 51)/       1.8999920D0/
      DATA ALPPM3( 51)/       2.0343010D0/
      DATA EISOLP( 51)/    -148.9382890D0/
      DATA GSSPM3( 51)/       9.2382770D0/
      DATA GSPPM3( 51)/       5.2776800D0/
      DATA GPPPM3( 51)/       6.3500000D0/
      DATA GP2PM3( 51)/       6.2500000D0/
      DATA HSPPM3( 51)/       2.4244640D0/
      DATA DDPM3 ( 51)/       1.4091903D0/
      DATA QQPM3 ( 51)/       1.3521354D0/
      DATA AMPM3 ( 51)/       0.3395177D0/
      DATA ADPM3 ( 51)/       0.4589010D0/
      DATA AQPM3 ( 51)/       0.2423472D0/
      DATA GUESP1( 51,1)/       3.0020280D0/
      DATA GUESP2( 51,1)/       6.0053420D0/
      DATA GUESP3( 51,1)/       0.8530600D0/
      DATA GUESP1( 51,2)/      -0.0188920D0/
      DATA GUESP2( 51,2)/       6.0114780D0/
      DATA GUESP3( 51,2)/       2.7933110D0/
C                    DATA FOR ELEMENT 52        TELLURIUM
      DATA REFPM3(52)/ ' Te: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1(ACCEPTED)                      '/
      DATA USSPM3( 52)/     -44.9380360D0/
      DATA UPPPM3( 52)/     -46.3140990D0/
      DATA BETASP( 52)/      -2.6651460D0/
      DATA BETAPP( 52)/      -3.8954300D0/
      DATA ZSPM3 ( 52)/       4.1654920D0/
      DATA ZPPM3 ( 52)/       1.6475550D0/
      DATA ALPPM3( 52)/       2.4850190D0/
      DATA EISOLP( 52)/    -168.0945925D0/
      DATA GSSPM3( 52)/      10.2550730D0/
      DATA GSPPM3( 52)/       8.1691450D0/
      DATA GPPPM3( 52)/       7.7775920D0/
      DATA GP2PM3( 52)/       7.7551210D0/
      DATA HSPPM3( 52)/       3.7724620D0/
      DATA DDPM3 ( 52)/       0.3484177D0/
      DATA QQPM3 ( 52)/       1.5593085D0/
      DATA AMPM3 ( 52)/       0.3768862D0/
      DATA ADPM3 ( 52)/       1.1960743D0/
      DATA AQPM3 ( 52)/       0.2184786D0/
      DATA GUESP1( 52,1)/       0.0333910D0/
      DATA GUESP2( 52,1)/       5.9563790D0/
      DATA GUESP3( 52,1)/       2.2775750D0/
      DATA GUESP1( 52,2)/      -1.9218670D0/
      DATA GUESP2( 52,2)/       4.9732190D0/
      DATA GUESP3( 52,2)/       0.5242430D0/
C                    DATA FOR ELEMENT 53      IODINE
      DATA REFPM3 (53)/'  I: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1 (IN PRESS) (1988)              '/
      DATA USSPM3( 53)/     -96.4540370D0/
      DATA UPPPM3( 53)/     -61.0915820D0/
      DATA BETASP( 53)/     -14.4942340D0/
      DATA BETAPP( 53)/      -5.8947030D0/
      DATA ZSPM3 ( 53)/       7.0010130D0/
      DATA ZPPM3 ( 53)/       2.4543540D0/
      DATA ZDPM3 ( 53)/       1.0000000D0/
      DATA ALPPM3( 53)/       1.9901850D0/
      DATA EISOLP( 53)/    -288.3160860D0/
      DATA GSSPM3( 53)/      13.6319430D0/
      DATA GSPPM3( 53)/      14.9904060D0/
      DATA GPPPM3( 53)/       7.2883300D0/
      DATA GP2PM3( 53)/       5.9664070D0/
      DATA HSPPM3( 53)/       2.6300350D0/
      DATA DDPM3 ( 53)/       0.1581469D0/
      DATA QQPM3 ( 53)/       1.0467302D0/
      DATA AMPM3 ( 53)/       0.5009902D0/
      DATA ADPM3 ( 53)/       1.6699104D0/
      DATA AQPM3 ( 53)/       0.5153082D0/
      DATA GUESP1( 53,1)/      -0.1314810D0/
      DATA GUESP2( 53,1)/       5.2064170D0/
      DATA GUESP3( 53,1)/       1.7488240D0/
      DATA GUESP1( 53,2)/      -0.0368970D0/
      DATA GUESP2( 53,2)/       6.0101170D0/
      DATA GUESP3( 53,2)/       2.7103730D0/
C                    DATA FOR ELEMENT 80        MERCURY
      DATA REFPM3(80)/ ' Hg: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1(ACCEPTED)                      '/
      DATA USSPM3( 80)/     -17.7622290D0/
      DATA UPPPM3( 80)/     -18.3307510D0/
      DATA BETASP( 80)/      -3.1013650D0/
      DATA BETAPP( 80)/      -3.4640310D0/
      DATA ZSPM3 ( 80)/       1.4768850D0/
      DATA ZPPM3 ( 80)/       2.4799510D0/
      DATA ALPPM3( 80)/       1.5293770D0/
      DATA EISOLP( 80)/     -28.8997380D0/
      DATA GSSPM3( 80)/       6.6247200D0/
      DATA GSPPM3( 80)/      10.6392970D0/
      DATA GPPPM3( 80)/      14.7092830D0/
      DATA GP2PM3( 80)/      16.0007400D0/
      DATA HSPPM3( 80)/       2.0363110D0/
      DATA DDPM3 ( 80)/       1.2317811D0/
      DATA QQPM3 ( 80)/       1.2164033D0/
      DATA AMPM3 ( 80)/       0.2434664D0/
      DATA ADPM3 ( 80)/       0.4515472D0/
      DATA AQPM3 ( 80)/       0.2618394D0/
      DATA GUESP1( 80,1)/       1.0827200D0/
      DATA GUESP2( 80,1)/       6.4965980D0/
      DATA GUESP3( 80,1)/       1.1951460D0/
      DATA GUESP1( 80,2)/      -0.0965530D0/
      DATA GUESP2( 80,2)/       3.9262810D0/
      DATA GUESP3( 80,2)/       2.6271600D0/
C                    DATA FOR ELEMENT 81        THALLIUM
      DATA REFPM3(81)/ ' Tl: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1(ACCEPTED)                      '/
      DATA USSPM3( 81)/     -30.0531700D0/
      DATA UPPPM3( 81)/     -26.9206370D0/
      DATA BETASP( 81)/      -1.0844950D0/
      DATA BETAPP( 81)/      -7.9467990D0/
      DATA ZSPM3 ( 81)/       6.8679210D0/
      DATA ZPPM3 ( 81)/       1.9694450D0/
      DATA ALPPM3( 81)/       1.3409510D0/
      DATA EISOLP( 81)/     -56.6492050D0/
      DATA GSSPM3( 81)/      10.4604120D0/
      DATA GSPPM3( 81)/      11.2238830D0/
      DATA GPPPM3( 81)/       4.9927850D0/
      DATA GP2PM3( 81)/       8.9627270D0/
      DATA HSPPM3( 81)/       2.5304060D0/
      DATA DDPM3 ( 81)/       0.0781362D0/
      DATA QQPM3 ( 81)/       1.5317110D0/
      DATA AMPM3 ( 81)/       0.3844326D0/
      DATA ADPM3 ( 81)/       2.5741815D0/
      DATA AQPM3 ( 81)/       0.2213264D0/
      DATA GUESP1( 81,1)/      -1.3613990D0/
      DATA GUESP2( 81,1)/       3.5572260D0/
      DATA GUESP3( 81,1)/       1.0928020D0/
      DATA GUESP1( 81,2)/      -0.0454010D0/
      DATA GUESP2( 81,2)/       2.3069950D0/
      DATA GUESP3( 81,2)/       2.9650290D0/
C                    DATA FOR ELEMENT 82        LEAD
      DATA REFPM3(82)/ ' Pb: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1(ACCEPTED)                      '/
      DATA USSPM3( 82)/     -30.3227560D0/
      DATA UPPPM3( 82)/     -24.4258340D0/
      DATA BETASP( 82)/      -6.1260240D0/
      DATA BETAPP( 82)/      -1.3954300D0/
      DATA ZSPM3 ( 82)/       3.1412890D0/
      DATA ZPPM3 ( 82)/       1.8924180D0/
      DATA ALPPM3( 82)/       1.6200450D0/
      DATA EISOLP( 82)/     -73.4660775D0/
      DATA GSSPM3( 82)/       7.0119920D0/
      DATA GSPPM3( 82)/       6.7937820D0/
      DATA GPPPM3( 82)/       5.1837800D0/
      DATA GP2PM3( 82)/       5.0456510D0/
      DATA HSPPM3( 82)/       1.5663020D0/
      DATA DDPM3 ( 82)/       0.9866290D0/
      DATA QQPM3 ( 82)/       1.5940562D0/
      DATA AMPM3 ( 82)/       0.2576991D0/
      DATA ADPM3 ( 82)/       0.4527678D0/
      DATA AQPM3 ( 82)/       0.2150175D0/
      DATA GUESP1( 82,1)/      -0.1225760D0/
      DATA GUESP2( 82,1)/       6.0030620D0/
      DATA GUESP3( 82,1)/       1.9015970D0/
      DATA GUESP1( 82,2)/      -0.0566480D0/
      DATA GUESP2( 82,2)/       4.7437050D0/
      DATA GUESP3( 82,2)/       2.8618790D0/
C                    DATA FOR ELEMENT 83        BISMUTH
      DATA REFPM3(83)/ ' Bi: (PM3): J. J. P. STEWART, J. COMP. CHEM.
     1(ACCEPTED)                      '/
      DATA USSPM3( 83)/     -33.4959380D0/
      DATA UPPPM3( 83)/     -35.5210260D0/
      DATA BETASP( 83)/      -5.6072830D0/
      DATA BETAPP( 83)/      -5.8001520D0/
      DATA ZSPM3 ( 83)/       4.9164510D0/
      DATA ZPPM3 ( 83)/       1.9349350D0/
      DATA ALPPM3( 83)/       1.8574310D0/
      DATA EISOLP( 83)/    -109.2774910D0/
      DATA GSSPM3( 83)/       4.9894800D0/
      DATA GSPPM3( 83)/       6.1033080D0/
      DATA GPPPM3( 83)/       8.6960070D0/
      DATA GP2PM3( 83)/       8.3354470D0/
      DATA HSPPM3( 83)/       0.5991220D0/
      DATA DDPM3 ( 83)/       0.2798609D0/
      DATA QQPM3 ( 83)/       1.5590294D0/
      DATA AMPM3 ( 83)/       0.1833693D0/
      DATA ADPM3 ( 83)/       0.6776013D0/
      DATA AQPM3 ( 83)/       0.2586520D0/
      DATA GUESP1( 83,1)/       2.5816930D0/
      DATA GUESP2( 83,1)/       5.0940220D0/
      DATA GUESP3( 83,1)/       0.4997870D0/
      DATA GUESP1( 83,2)/       0.0603200D0/
      DATA GUESP2( 83,2)/       6.0015380D0/
      DATA GUESP3( 83,2)/       2.4279700D0/
      DATA REFPM3 (102)/' Cb: (PM3):  Capped Bond  (Hydrogen-like, takes
     1 on a  zero charge.)             '/
      DATA USSPM3(102)/     -11.9062760D0/
      DATA BETASP(102)/-9999999.0000000D0/
      DATA ZSPM3 (102)/       4.0000000D0/
      DATA ZPPM3 (102)/       0.3000000D0/
      DATA ZDPM3 (102)/       0.3000000D0/
      DATA ALPPM3(102)/       2.5441341D0/
      DATA EISOLP(102)/       4.0000000D0/
      DATA GSSPM3(102)/      12.8480000D0/
      DATA HSPPM3(102)/       0.1000000D0/
      DATA DDPM3 (102)/       0.0684105D0/
      DATA QQPM3 (102)/       1.0540926D0/
      DATA AMPM3 (102)/       0.4721793D0/
      DATA ADPM3 (102)/       0.9262742D0/
      DATA AQPM3 (102)/       0.2909059D0/
      END
      SUBROUTINE BONDS(P)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      PARAMETER (NATMS2=MAXPAR*MAXPAR-MAXORB*MAXORB)
      DIMENSION P(*),B(MAXORB,MAXORB)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      LOGICAL PRINT
      COMMON /BONDCM/ BONDAB(NATMS2),PRINT
C*TNT      COMMON /SCRACH/ B(MAXORB,MAXORB), BONDAB(NATMS2)
C***********************************************************************
C
C   CALCULATES, AND PRINTS, THE BOND INDICES AND VALENCIES OF ATOMS
C
C  FOR REFERENCE, SEE "BOND INDICES AND VALENCY", J. C. S. DALTON,
C  ARMSTRONG, D.R., PERKINS, P.G., AND STEWART, J.J.P., 838 (1973)
C
C   ON INPUT
C            P = DENSITY MATRIX, LOWER HALF TRIANGLE, PACKED.
C            P   IS NOT ALTERED BY BONDS.
C
C***********************************************************************
         SAVE                                                           GL0892
      IF(PRINT) THEN
        WRITE(6,10)
   10   FORMAT(//20X,'BOND ORDERS AND VALENCIES',//)
      END IF
      K=0
      DO 20 I=1,NORBS
         DO 20 J=1,I
            K=K+1
            B(I,J)=P(K)
   20 B(J,I)=P(K)
      IJ = 0
      DO 60 I=1,NUMAT
         L=NFIRST(I)
         LL=NLAST(I)
         DO 40 J=1,I
            IJ = IJ + 1
            K=NFIRST(J)
            KK=NLAST(J)
            X=0.0
            DO 30 IL=L,LL
               DO 30 IH=K,KK
   30       X=X+B(IL,IH)*B(IL,IH)
   40    BONDAB(IJ)=X
         X=-BONDAB(IJ)
         DO 50 J=L,LL
   50    X=X+2.D0*B(J,J)
         BONDAB(IJ)=X
   60 CONTINUE
      IF (PRINT) THEN
      CALL VECPRT( BONDAB, NUMAT)
      END IF
      RETURN
      END
      SUBROUTINE CNVG(PNEW, P, P1,NORBS, NITER, PL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION P1(*), P(*), PNEW(*)
      LOGICAL EXTRAP
C***********************************************************************
C
C  CNVG IS A TWO-POINT INTERPOLATION ROUTINE FOR SPEEDING CONVERGENCE
C       OF THE DENSITY MATRIX.
C
C ON OUTPUT P      = NEW DENSITY MATRIX
C           P1     = DIAGONAL OF OLD DENSITY MATRIX
C           PL     = LARGEST DIFFERENCE BETWEEN OLD AND NEW DENSITY
C                    MATRICES
C***********************************************************************
      COMMON/KEYWRD/ KEYWRD
      COMMON /NUMCAL/ NUMCAL
      CHARACTER*80 KEYWRD
         SAVE                                                           GL0892
      DATA ICALCN /0/
      IF (ICALCN .NE. NUMCAL) THEN
         ICALCN = NUMCAL
         IF(INDEX(KEYWRD,'UHF').NE.0)THEN
            RHFUHF=1.D0
         ELSE
            RHFUHF=2.D0
         ENDIF
      ENDIF
      PL=0.0D00
      FACA=0.0D00
      DAMP=1.D10
      IF(NITER.GT.3)DAMP=0.05D0
      FACB=0.0D00
      FAC=0.0D00
      II=MOD(NITER,3)
      EXTRAP=II.NE.0
      SUM1=0.D0
      K=0
      DO 20 I=1,NORBS
         K=K+I
         A=PNEW(K)
         SUM1=SUM1+A
         SA=ABS(A-P(K))
         IF (SA.GT.PL) PL=SA
         IF (EXTRAP) GO TO 10
         FACA=FACA+SA**2
         FACB=FACB+(A-2.D00*P(K)+P1(I))**2
   10    P1(I)=P(K)
   20 P(K)=A
      IF (FACB.LE.0.0D00) GO TO 30
      IF (FACA.LT.(100.D00*FACB)) FAC=SQRT(FACA/FACB)
   30 IE=0
      SUM2=0.D0
      DO 50 I=1,NORBS
         II=I-1
         DO 40 J=1,II
            IE=IE+1
            A=PNEW(IE)
            P(IE)=A+FAC*(A-P(IE))
            PNEW(IE)=P(IE)
   40    CONTINUE
         IE=IE+1
         IF(ABS(P(IE)-P1(I)) .GT. DAMP) THEN
            P(IE)=P1(I)+SIGN(DAMP,P(IE)-P1(I))
         ELSE
            P(IE)=P(IE)+FAC*(P(IE)-P1(I))
         ENDIF
         P(IE)=MIN(RHFUHF,MAX(P(IE),0.D0))
         SUM2=SUM2+P(IE)
   50 PNEW(IE)=P(IE)
C
C   RE-NORMALIZE IF ANY DENSITY MATRIX ELEMENTS HAVE BEEN TRUNCATED
C
      SUM0=SUM1
   60 SUM=SUM1/SUM2
      SUM1=SUM0
      IF(ABS(SUM-1.D0).GT.1.D-5)THEN
C#      WRITE(6,'(6F12.6)')(P((I*(I+1))/2),I=1,NORBS)
         SUM2=0.D0
         DO 70 I=1,NORBS
            J=(I*(I+1))/2
C
C   ADD ON A SMALL NUMBER IN CASE AN OCCUPANCY IS EXACTLY ZERO
C
            P(J)=P(J)*SUM+1.D-20
            P(J)=MAX(P(J),0.D0)
C
C  SET UP RENORMALIZATION OVER PARTLY OCCUPIED M.O.'S ONLY.  FULL M.O.'S
C  CAN'T BE FILLED ANY MORE
C
            IF(P(J).GT.RHFUHF)THEN
               P(J)=RHFUHF
               SUM1=SUM1-RHFUHF
            ELSE
               SUM2=SUM2+P(J)
            ENDIF
   70    PNEW(J)=P(J)
         GOTO 60
      ENDIF
      RETURN
      END
      SUBROUTINE COMPFG(XPARAM,INT,ESCF,FULSCF,GRAD,LGRAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION XPARAM(MAXPAR),GRAD(MAXPAR)
      LOGICAL LGRAD, FULSCF
      COMMON /GEOVAR/ XDUMMY(MAXPAR),NVAR,LOC(2,MAXPAR)                 IR0394
      COMMON /GEOSYM/ NDEP,LOCPAR(MAXPAR),IDEPFN(MAXPAR),LOCDEP(MAXPAR)
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /ATHEAT/ ATHEAT
      COMMON /WMATRX/ WJ(N2ELEC), WK(N2ELEC)
      COMMON /ENUCLR/ ENUCLR
      COMMON /NATYPE/ NZTYPE(107),MTYPE(30),LTYPE
      COMMON /ELECT / ELECT
      PARAMETER (MDUMY=MAXPAR**2-MPACK)
      COMMON /SCRACH/ RXYZ(MPACK), XDUMY(MDUMY)
      COMMON /HMATRX/ H(MPACK)
      COMMON /MOLMEC/ HTYPE(4),NHCO(4,20),NNHCO,ITYPE,USEMM
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
C***********************************************************************
C
C   COMPFG CALCULATES (A) THE HEAT OF FORMATION OF THE SYSTEM, AND
C                     (B) THE GRADIENTS, IF LGRAD IS .TRUE.
C
C   ON INPUT  XPARAM = ARRAY OF PARAMETERS TO BE USED IN INTERNAL COORDS
C             LGRAD  = .TRUE. IF GRADIENTS ARE NEEDED, .FALSE. OTHERWISE
C             INT    = NEVER USED, RESERVED FOR FUTURE USE
C             FULSCF = .TRUE. IF FULL SCF TO BE DONE, .FALSE. OTHERWISE.
C
C   ON OUTPUT ESCF  = HEAT OF FORMATION.
C             GRAD   = ARRAY OF GRADIENTS, IF LGRAD = .TRUE.
C
C***********************************************************************
      COMMON /NUMCAL/ NUMCAL
      COMMON /MESAGE/ IFLEPO, IERR                                      IR0694
      COMMON /KEYWRD/KEYWRD
      CHARACTER*80 KEYWRD
      LOGICAL DEBUG, INT, PRINT, ANALYT, LARGE, USEMM
      DIMENSION COORD(3,NUMATM), W(N2ELEC), DEGREE(3)
      EQUIVALENCE (W,WJ)
         SAVE                                                           GL0892
      DATA ICALCN /0/
C                 MNDO     AM1      PM3      MINDO/3
C
C On IBM you can initialize common only inside BLOCK DATA               IR0394
C     DATA HTYPE/ 6.1737D0,  3.3191D0,  7.1853D0,  1.7712D0/
      htype(1)=6.1737D0
      htype(2)=3.3191D0
      htype(3)=7.1853D0
      htype(4)=1.7712D0
C
      IF (ICALCN .NE. NUMCAL) THEN
         ICALCN = NUMCAL
         LTYPE=0
         DO 30 I=1,NUMAT
            IF(NAT(I).LT.99)THEN
               DO 10 J=1,LTYPE
   10          IF(NAT(I).EQ.MTYPE(J)) GOTO 20
               LTYPE=LTYPE+1
               MTYPE(LTYPE)=NAT(I)
               NZTYPE(NAT(I))=LTYPE
C
C       LTYPE = NUMBER OF TYPES OF REAL ATOM PRESENT
C       MTYPE = TYPES OF REAL ATOMS PRESENT
               J=LTYPE
   20          CONTINUE
            ENDIF
   30    CONTINUE
         ANALYT=(INDEX(KEYWRD,'ANALYT').NE.0)
         IF(ANALYT)CALL SETUPG
         DEGREE(1)=1.D0
         IF(INDEX(KEYWRD,' XYZ').NE.0)THEN
            DEGREE(2)=1.D0
         ELSE
            DEGREE(2)=180.D0/3.141592652589D0
         ENDIF
         DEGREE(3)=DEGREE(2)
         LARGE=(INDEX(KEYWRD,'LARGE') .NE. 0)
         PRINT=(INDEX(KEYWRD,'COMPFG') .NE. 0)
         DEBUG=(INDEX(KEYWRD,'DEBUG') .NE. 0 .AND. PRINT)
      ENDIF
C
C SET UP COORDINATES FOR CURRENT CALCULATION
C
C       PLACE THE NEW VALUES OF THE VARIABLES IN THE ARRAY GEO.
C       MAKE CHANGES IN THE GEOMETRY.
      DO 40 I=1,NVAR
         K=LOC(1,I)
         L=LOC(2,I)
   40 GEO(L,K)=XPARAM(I)
C#      WRITE(6,'(3F18.11)')(XPARAM(I),I=1,3)
C      IMPOSE THE SYMMETRY CONDITIONS + COMPUTE THE DEPENDENT-PARAMETERS
      IF(NDEP.NE.0) CALL SYMTRY
C      NOW COMPUTE THE ATOMIC COORDINATES.
      IF( DEBUG ) THEN
         IF( LARGE ) THEN
            K=NUMAT
         ELSE
            K=MIN(5,NUMAT)
         ENDIF
         WRITE(6,FMT='('' INTERNAL COORDS'',/100(/,3F12.6))')
     1            ((GEO(J,I)*DEGREE(J),J=1,3),I=1,K)
      END IF
      CALL GMETRY(GEO,COORD)
      IF( DEBUG ) THEN
         IF( LARGE ) THEN
            K=NUMAT
         ELSE
            K=MIN(5,NUMAT)
         ENDIF
         WRITE(6,FMT='('' CARTESIAN COORDS'',/100(/,3F12.6))')
     1            ((COORD(J,I),J=1,3),I=1,K)
      END IF
      IF(ANALYT)REWIND 2
      CALL HCORE(COORD, H, W, WJ, WK, ENUCLR)
C
C COMPUTE THE HEAT OF FORMATION.
C
      IF(NORBS.GT.0.AND.NELECS.GT.0) THEN
         CALL ITER(H, W, WJ, WK, ELECT, FULSCF,.TRUE.)
C GA-FRIENDLY : uncomment the next line if linked with GA code
C        IF(IERR .eq. 2) RETURN                                           GA-ON
      ELSE
         ELECT=0.D0
      ENDIF
      ESCF=(ELECT+ENUCLR)*23.061D0+ATHEAT
      DO 50 I=1,NNHCO
         CALL DIHED(COORD,NHCO(1,I),NHCO(2,I),NHCO(3,I),NHCO(4,I),ANGLE)
         ESCF=ESCF+HTYPE(ITYPE)*SIN(ANGLE)**2
   50 CONTINUE
      IF(PRINT) WRITE(6,'(/10X,'' HEAT OF FORMATION'',G30.17)')ESCF
C
C FIND DERIVATIVES IF DESIRED
C
      IF(LGRAD) CALL DERIV(GEO,GRAD)
      RETURN
      END
      SUBROUTINE DCART (COORD,DXYZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION COORD(3,*), DXYZ(3,*)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /DENSTY/ P(MPACK), PA(MPACK), PB(MPACK)
C***********************************************************************
C
C    DCART CALCULATES THE DERIVATIVES OF THE ENERGY WITH RESPECT TO THE
C          CARTESIAN COORDINATES. THIS IS DONE BY FINITE DIFFERENCES.
C
C    THE MAIN ARRAYS IN DCART ARE:
C        DXYZ   ON EXIT CONTAINS THE CARTESIAN DERIVATIVES.
C
C***********************************************************************
      COMMON /KEYWRD/ KEYWRD
      COMMON /EULER / TVEC(3,3), ID
      COMMON /MOLMEC/ HTYPE(4),NHCO(4,20),NNHCO,ITYPE,USEMM
      COMMON /UCELL / L1L,L2L,L3L,L1U,L2U,L3U
      COMMON /DCARTC/ K1L,K2L,K3L,K1U,K2U,K3U
      COMMON /NUMCAL/ NUMCAL
      CHARACTER*80 KEYWRD
      DIMENSION PDI(171),PADI(171),PBDI(171),
     1CDI(3,2),NDI(2),LSTOR1(6), LSTOR2(6), ENG(3)
      LOGICAL DEBUG, FORCE, MAKEP, ANADER, USEMM
      EQUIVALENCE (LSTOR1(1),L1L), (LSTOR2(1), K1L)
         SAVE                                                           GL0892
      DATA CHNGE,CHNGE2 /1.D-4,5.D-5/
*
* CHNGE IS A MACHINE-PRECISION DEPENDENT CONSTANT
* CHNGE2=CHNGE/2
*
      DATA ICALCN /0/
      IF (ICALCN .NE. NUMCAL) THEN
         ANADER= (INDEX(KEYWRD,'ANALYT') .NE. 0)
         DEBUG = (INDEX(KEYWRD,'DCART') .NE. 0)
         FORCE = (INDEX(KEYWRD,'PRECISE')+INDEX(KEYWRD,'FORCE') .NE. 0)
         ICALCN = NUMCAL
      ENDIF
      NCELLS=(L1U-L1L+1)*(L2U-L2L+1)*(L3U-L3L+1)
      DO 10 I=1,6
         LSTOR2(I)=LSTOR1(I)
   10 LSTOR1(I)=0
      IOFSET=(NCELLS+1)/2
      NUMTOT=NUMAT*NCELLS
      DO 20 I=1,NUMTOT
         DO 20 J=1,3
   20 DXYZ(J,I)=0.D0
      IF(ANADER) REWIND 2
      KREP=0
      DO 130 II=1,NUMAT
         III=NCELLS*(II-1)+IOFSET
         IM1=II
         IF=NFIRST(II)
         IM=NMIDLE(II)
         IL=NLAST(II)
         NDI(2)=NAT(II)
         DO 30 I=1,3
   30    CDI(I,2)=COORD(I,II)
         DO 130 JJ=1,IM1
            JJJ=NCELLS*(JJ-1)
C  FORM DIATOMIC MATRICES
            JF=NFIRST(JJ)
            JM=NMIDLE(JJ)
            JL=NLAST(JJ)
C   GET FIRST ATOM
            NDI(1)=NAT(JJ)
            MAKEP=.TRUE.
            DO 120 IK=K1L,K1U
               DO 120 JK=K2L,K2U
                  DO 120 KL=K3L,K3U
                     JJJ=JJJ+1
                     DO 40 L=1,3
   40                CDI(L,1)=COORD(L,JJ)+TVEC(L,1)*IK+TVEC(L,2)*JK+TVEC
     1(L,3)*KL
                     IF(.NOT. MAKEP) GOTO 90
                     MAKEP=.FALSE.
                     IJ=0
                     DO 50 I=JF,JL
                        K=I*(I-1)/2+JF-1
                        DO 50 J=JF,I
                           IJ=IJ+1
                           K=K+1
                           PADI(IJ)=PA(K)
                           PBDI(IJ)=PB(K)
   50                PDI(IJ)=P(K)
C GET SECOND ATOM FIRST ATOM INTERSECTION
                     DO 80 I=IF,IL
                        L=I*(I-1)/2
                        K=L+JF-1
                        DO 60 J=JF,JL
                           IJ=IJ+1
                           K=K+1
                           PADI(IJ)=PA(K)
                           PBDI(IJ)=PB(K)
   60                   PDI(IJ)=P(K)
                        K=L+IF-1
                        DO 70 L=IF,I
                           K=K+1
                           IJ=IJ+1
                           PADI(IJ)=PA(K)
                           PBDI(IJ)=PB(K)
   70                   PDI(IJ)=P(K)
   80                CONTINUE
   90                CONTINUE
                     IF(II.EQ.JJ) GOTO  120
                     IF(ANADER)THEN
                        CALL ANALYT(PDI,PADI,PBDI,CDI,NDI,JF,JL,IF,IL
     1,                 NORBS,ENG,KREP)
                        DO 100 K=1,3
                           DXYZ(K,III)=DXYZ(K,III)+ENG(K)
  100                   DXYZ(K,JJJ)=DXYZ(K,JJJ)-ENG(K)
                     ELSE
                        IF( .NOT. FORCE) THEN
                           CDI(1,1)=CDI(1,1)+CHNGE2
                           CDI(2,1)=CDI(2,1)+CHNGE2
                           CDI(3,1)=CDI(3,1)+CHNGE2
                           CALL DHC(PDI,PADI,PBDI,CDI,NDI,JF,JM,JL,IF,IM
     1,IL,                 NORBS,AA)
                        ENDIF
                        DO 110 K=1,3
                           IF( FORCE )THEN
                              CDI(K,2)=CDI(K,2)-CHNGE2
                              CALL DHC(PDI,PADI,PBDI,CDI,NDI,JF,JM,JL,IF
     1,IM,IL,                 NORBS,AA)
                           ENDIF
                           CDI(K,2)=CDI(K,2)+CHNGE
                           CALL DHC(PDI,PADI,PBDI,CDI,NDI,JF,JM,JL,IF,IM
     1,IL,                 NORBS,EE)
                           CDI(K,2)=CDI(K,2)-CHNGE2
                           IF( .NOT. FORCE) CDI(K,2)=CDI(K,2)-CHNGE2
                           DERIV=(AA-EE)*46.122D0/CHNGE
                           DXYZ(K,III)=DXYZ(K,III)+DERIV
                           DXYZ(K,JJJ)=DXYZ(K,JJJ)-DERIV
  110                   CONTINUE
                     ENDIF
  120       CONTINUE
  130 CONTINUE
      IF(USEMM)THEN
C
C   NOW ADD IN MOLECULAR-MECHANICS CORRECTION TO THE H-N-C=O TORSION
C
         DEL=1.D-5
         DO 160 I=1,NNHCO
            DO 150 J=1,4
               DO 140 K=1,3
                  COORD(K,NHCO(J,I))=COORD(K,NHCO(J,I))-DEL
                  CALL DIHED(COORD,NHCO(1,I),NHCO(2,I),NHCO(3,I),NHCO(4,
     1I),ANGLE)
                  REFH=HTYPE(ITYPE)*SIN(ANGLE)**2
                  COORD(K,NHCO(J,I))=COORD(K,NHCO(J,I))+DEL*2.D0
                  CALL DIHED(COORD,NHCO(1,I),NHCO(2,I),NHCO(3,I),NHCO(4,
     1I),ANGLE)
                  COORD(K,NHCO(J,I))=COORD(K,NHCO(J,I))-DEL
                  HEAT=HTYPE(ITYPE)*SIN(ANGLE)**2
                  SUM=(REFH-HEAT)/DEL
                  DXYZ(K,NHCO(J,I))=DXYZ(K,NHCO(J,I))+SUM
  140          CONTINUE
  150       CONTINUE
  160    CONTINUE
      ENDIF
      DO 170 I=1,6
  170 LSTOR1(I)=LSTOR2(I)
      IF (  .NOT. DEBUG) RETURN
      WRITE(6,'(//10X,''CARTESIAN COORDINATE DERIVATIVES'',//3X,
     1''ATOM  AT. NO.'',5X,''X'',12X,''Y'',12X,''Z'',/)')
      WRITE(6,'(2I6,F13.6,2F13.6)')
     1 (I,NAT((I-1)/NCELLS+1),(DXYZ(J,I),J=1,3),I=1,NUMTOT)
      IF (ANADER) REWIND 2
      RETURN
      END
      SUBROUTINE DELRI(DG,NI,NJ,RR,DEL1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DG(22)
************************************************************************
*                                                                      *
*    ON INPUT NI = ATOMIC NUMBER OF FIRST ATOM                         *
*             NJ = ATOMIC NUMBER OF SECOND ATOM                        *
*             RR = INTERATOMIC DISTANCE IN BOHRS                       *
*                                                                      *
************************************************************************
      COMMON  /MULTIP/  DD(107),QQ(107),BDD(107,3)
      COMMON /TWOEL3/ F03(107)
      COMMON /NATORB/ NATORB(107)
      COMMON /ALPHA3/ ALP3(153)
      COMMON /KEYWRD/ KEYWRD
      COMMON /NUMCAL/ NUMCAL
      CHARACTER*80 KEYWRD
      LOGICAL MINDO3
         SAVE                                                           GL0892
      DATA ICALCN /0/
      IF (ICALCN .NE. NUMCAL) THEN
         ICALCN = NUMCAL
         MINDO3=(INDEX(KEYWRD,'MINDO').NE.0)
      ENDIF
      A0=0.529167D0
      TERM=(27.21D0*DEL1)/(RR*A0*A0)
      DA=DD(NI)
      DB=DD(NJ)
      QA=QQ(NI)
      QB=QQ(NJ)
C   HYDROGEN-HYDROGEN
      AEE=0.25D0*(1.0D0/BDD(NI,1)+1.0D0/BDD(NJ,1))**2
      EE    =-RR/(SQRT(RR**2+AEE))**3
      DG(1)=TERM*EE
      IF(NATORB(NI).LE.2.AND.NATORB(NJ).LE.2) RETURN
      IF(NATORB(NI).LE.2) GO TO 10
C   HEAVY ATOM-HYDROGEN
      ADE=0.25D0*(1.0D0/BDD(NI,2)+1.0D0/BDD(NJ,1))**2
      AQE=0.25D0*(1.0D0/BDD(NI,3)+1.0D0/BDD(NJ,1))**2
      DZE   = (RR+DA)/(SQRT((RR+DA)**2+ADE))**3
     1       -(RR-DA)/(SQRT((RR-DA)**2+ADE))**3
      QZZE  =-(RR+2.0D0*QA)/(SQRT((RR+2.0D0*QA)**2+AQE))**3
     1       -(RR-2.0D0*QA)/(SQRT((RR-2.0D0*QA)**2+AQE))**3
     2       +(2.0D0*RR)/(SQRT(RR**2+AQE))**3
      QXXE  =-(2.0D0*RR)/(SQRT(RR**2+4.0D0*QA**2+AQE))**3
     1       +(2.0D0*RR)/(SQRT(RR**2+AQE))**3
      DG(2)=-(TERM*DZE)/2.0D0
      DG(3)=TERM*(EE+QZZE/4.0D0)
      DG(4)=TERM*(EE+QXXE/4.0D0)
      IF(NATORB(NJ).LE.2) RETURN
C   HYDROGEN-HEAVY ATOM
   10 AED=0.25D0*(1.0D0/BDD(NI,1)+1.0D0/BDD(NJ,2))**2
      AEQ=0.25D0*(1.0D0/BDD(NI,1)+1.0D0/BDD(NJ,3))**2
      EDZ   = (RR-DB)/(SQRT((RR-DB)**2+AED))**3
     1       -(RR+DB)/(SQRT((RR+DB)**2+AED))**3
      EQZZ  =-(RR-2.0D0*QB)/(SQRT((RR-2.0D0*QB)**2+AEQ))**3
     1       -(RR+2.0D0*QB)/(SQRT((RR+2.0D0*QB)**2+AEQ))**3
     2       +(2.0D0*RR)/(SQRT(RR**2+AEQ))**3
      EQXX  =-(2.0D0*RR)/(SQRT(RR**2+4.0D0*QB**2+AEQ))**3
     1       +(2.0D0*RR)/(SQRT(RR**2+AEQ))**3
      DG(5)=-(TERM*EDZ)/2.0D0
      DG(11)=TERM*(EE+EQZZ/4.0D0)
      DG(12)=TERM*(EE+EQXX/4.0D0)
      IF(NATORB(NI).LE.2) RETURN
C   HEAVY ATOM-HEAVY ATOM
      ADD=0.25D0*(1.D0/BDD(NI,2)+1.D0/BDD(NJ,2))**2
      ADQ=0.25D0*(1.D0/BDD(NI,2)+1.D0/BDD(NJ,3))**2
      AQD=0.25D0*(1.D0/BDD(NI,3)+1.D0/BDD(NJ,2))**2
      AQQ=0.25D0*(1.D0/BDD(NI,3)+1.D0/BDD(NJ,3))**2
      DXDX  =-(2.D0*RR)/(SQRT(RR**2+(DA-DB)**2+ADD))**3
     1       +(2.D0*RR)/(SQRT(RR**2+(DA+DB)**2+ADD))**3
      DZDZ  =-(RR+DA-DB)/(SQRT((RR+DA-DB)**2+ADD))**3
     1       -(RR-DA+DB)/(SQRT((RR-DA+DB)**2+ADD))**3
     2       +(RR-DA-DB)/(SQRT((RR-DA-DB)**2+ADD))**3
     3       +(RR+DA+DB)/(SQRT((RR+DA+DB)**2+ADD))**3
      DZQXX = 2.D0*(RR+DA)/(SQRT((RR+DA)**2+4.D0*QB**2+ADQ))**3
     1       -2.D0*(RR-DA)/(SQRT((RR-DA)**2+4.D0*QB**2+ADQ))**3
     2       -2.D0*(RR+DA)/(SQRT((RR+DA)**2+ADQ))**3
     3       +2.D0*(RR-DA)/(SQRT((RR-DA)**2+ADQ))**3
      QXXDZ = 2.D0*(RR-DB)/(SQRT((RR-DB)**2+4.D0*QA**2+AQD))**3
     1       -2.D0*(RR+DB)/(SQRT((RR+DB)**2+4.D0*QA**2+AQD))**3
     2       -2.D0*(RR-DB)/(SQRT((RR-DB)**2+AQD))**3
     3       +2.D0*(RR+DB)/(SQRT((RR+DB)**2+AQD))**3
      DZQZZ = (RR+DA-2.D0*QB)/(SQRT((RR+DA-2.D0*QB)**2+ADQ))**3
     1       -(RR-DA-2.D0*QB)/(SQRT((RR-DA-2.D0*QB)**2+ADQ))**3
     2       +(RR+DA+2.D0*QB)/(SQRT((RR+DA+2.D0*QB)**2+ADQ))**3
     3       -(RR-DA+2.D0*QB)/(SQRT((RR-DA+2.D0*QB)**2+ADQ))**3
     4       +2.D0*(RR-DA)/(SQRT((RR-DA)**2+ADQ))**3
     5       -2.D0*(RR+DA)/(SQRT((RR+DA)**2+ADQ))**3
      QZZDZ = (RR+2.D0*QA-DB)/(SQRT((RR+2.D0*QA-DB)**2+AQD))**3
     1       -(RR+2.D0*QA+DB)/(SQRT((RR+2.D0*QA+DB)**2+AQD))**3
     2       +(RR-2.D0*QA-DB)/(SQRT((RR-2.D0*QA-DB)**2+AQD))**3
     3       -(RR-2.D0*QA+DB)/(SQRT((RR-2.D0*QA+DB)**2+AQD))**3
     4       -2.D0*(RR-DB)/(SQRT((RR-DB)**2+AQD))**3
     5       +2.D0*(RR+DB)/(SQRT((RR+DB)**2+AQD))**3
      QXXQXX=-(2.D0*RR)/(SQRT(RR**2+4.D0*(QA-QB)**2+AQQ))**3
     1       -(2.D0*RR)/(SQRT(RR**2+4.D0*(QA+QB)**2+AQQ))**3
     2       +(4.D0*RR)/(SQRT(RR**2+4.D0*QA**2+AQQ))**3
     3       +(4.D0*RR)/(SQRT(RR**2+4.D0*QB**2+AQQ))**3
     4       -(4.D0*RR)/(SQRT(RR**2+AQQ))**3
      QXXQYY=-(4.D0*RR)/(SQRT(RR**2+4.D0*QA**2+4.D0*QB**2+AQQ))**3
     1       +(4.D0*RR)/(SQRT(RR**2+4.D0*QA**2+AQQ))**3
     2       +(4.D0*RR)/(SQRT(RR**2+4.D0*QB**2+AQQ))**3
     3       -(4.D0*RR)/(SQRT(RR**2+AQQ))**3
      QXXQZZ=
     1     -2.D0*(RR-2.D0*QB)/(SQRT((RR-2.D0*QB)**2+4.D0*QA**2+AQQ))**3
     2     -2.D0*(RR+2.D0*QB)/(SQRT((RR+2.D0*QB)**2+4.D0*QA**2+AQQ))**3
     3       +2.D0*(RR-2.D0*QB)/(SQRT((RR-2.D0*QB)**2+AQQ))**3
     4       +2.D0*(RR+2.D0*QB)/(SQRT((RR+2.D0*QB)**2+AQQ))**3
     5       +(4.D0*RR)/(SQRT(RR**2+4.D0*QA**2+AQQ))**3
     6       -(4.D0*RR)/(SQRT(RR**2+AQQ))**3
      QZZQXX=
     1     -2.D0*(RR+2.D0*QA)/(SQRT((RR+2.D0*QA)**2+4.D0*QB**2+AQQ))**3
     2     -2.D0*(RR-2.D0*QA)/(SQRT((RR-2.D0*QA)**2+4.D0*QB**2+AQQ))**3
     3       +2.D0*(RR+2.D0*QA)/(SQRT((RR+2.D0*QA)**2+AQQ))**3
     4       +2.D0*(RR-2.D0*QA)/(SQRT((RR-2.D0*QA)**2+AQQ))**3
     5       +(4.D0*RR)/(SQRT(RR**2+4.D0*QB**2+AQQ))**3
     6       -(4.D0*RR)/(SQRT(RR**2+AQQ))**3
      QZZQZZ=
     1     -(RR+2.D0*QA-2.D0*QB)/(SQRT((RR+2.D0*QA-2.D0*QB)**2+AQQ))**3
     2     -(RR+2.D0*QA+2.D0*QB)/(SQRT((RR+2.D0*QA+2.D0*QB)**2+AQQ))**3
     3     -(RR-2.D0*QA-2.D0*QB)/(SQRT((RR-2.D0*QA-2.D0*QB)**2+AQQ))**3
     4     -(RR-2.D0*QA+2.D0*QB)/(SQRT((RR-2.D0*QA+2.D0*QB)**2+AQQ))**3
     5       +2.D0*(RR-2.D0*QA)/(SQRT((RR-2.D0*QA)**2+AQQ))**3
     6       +2.D0*(RR+2.D0*QA)/(SQRT((RR+2.D0*QA)**2+AQQ))**3
     7       +2.D0*(RR-2.D0*QB)/(SQRT((RR-2.D0*QB)**2+AQQ))**3
     8       +2.D0*(RR+2.D0*QB)/(SQRT((RR+2.D0*QB)**2+AQQ))**3
     9       -(4.D0*RR)/(SQRT(RR**2+AQQ))**3
      DXQXZ = 2.D0*(RR-QB)/(SQRT((RR-QB)**2+(DA-QB)**2+ADQ))**3
     1       -2.D0*(RR+QB)/(SQRT((RR+QB)**2+(DA-QB)**2+ADQ))**3
     2       -2.D0*(RR-QB)/(SQRT((RR-QB)**2+(DA+QB)**2+ADQ))**3
     3       +2.D0*(RR+QB)/(SQRT((RR+QB)**2+(DA+QB)**2+ADQ))**3
      QXZDX = 2.D0*(RR+QA)/(SQRT((RR+QA)**2+(QA-DB)**2+AQD))**3
     1       -2.D0*(RR-QA)/(SQRT((RR-QA)**2+(QA-DB)**2+AQD))**3
     2       -2.D0*(RR+QA)/(SQRT((RR+QA)**2+(QA+DB)**2+AQD))**3
     3       +2.D0*(RR-QA)/(SQRT((RR-QA)**2+(QA+DB)**2+AQD))**3
      QXYQXY=-(4.D0*RR)/(SQRT(RR**2+2.D0*(QA-QB)**2+AQQ))**3
     1       -(4.D0*RR)/(SQRT(RR**2+2.D0*(QA+QB)**2+AQQ))**3
     2       +(8.D0*RR)/(SQRT(RR**2+2.D0*(QA**2+QB**2)+AQQ))**3
      QXZQXZ=-2.D0*(RR+QA-QB)/(SQRT((RR+QA-QB)**2+(QA-QB)**2+AQQ))**3
     1       +2.D0*(RR+QA+QB)/(SQRT((RR+QA+QB)**2+(QA-QB)**2+AQQ))**3
     2       +2.D0*(RR-QA-QB)/(SQRT((RR-QA-QB)**2+(QA-QB)**2+AQQ))**3
     3       -2.D0*(RR-QA+QB)/(SQRT((RR-QA+QB)**2+(QA-QB)**2+AQQ))**3
     4       +2.D0*(RR+QA-QB)/(SQRT((RR+QA-QB)**2+(QA+QB)**2+AQQ))**3
     5       -2.D0*(RR+QA+QB)/(SQRT((RR+QA+QB)**2+(QA+QB)**2+AQQ))**3
     6       -2.D0*(RR-QA-QB)/(SQRT((RR-QA-QB)**2+(QA+QB)**2+AQQ))**3
     7       +2.D0*(RR-QA+QB)/(SQRT((RR-QA+QB)**2+(QA+QB)**2+AQQ))**3
      DG(6)=(TERM*DZDZ)/4.0D0
      DG(7)=(TERM*DXDX)/4.0D0
      DG(8)=-TERM*(EDZ/2.0D0+QZZDZ/8.0D0)
      DG(9)=-TERM*(EDZ/2.0D0+QXXDZ/8.0D0)
      DG(10)=-(TERM*QXZDX)/8.0D0
      DG(13)=-TERM*(DZE/2.0D0+DZQZZ/8.0D0)
      DG(14)=-TERM*(DZE/2.0D0+DZQXX/8.0D0)
      DG(15)=-(TERM*DXQXZ)/8.0D0
      DG(16)=TERM*(EE+EQZZ/4.0D0+QZZE/4.0D0+QZZQZZ/16.0D0)
      DG(17)=TERM*(EE+EQZZ/4.0D0+QXXE/4.0D0+QXXQZZ/16.0D0)
      DG(18)=TERM*(EE+EQXX/4.0D0+QZZE/4.0D0+QZZQXX/16.0D0)
      DG(19)=TERM*(EE+EQXX/4.0D0+QXXE/4.0D0+QXXQXX/16.0D0)
      DG(20)=(TERM*QXZQXZ)/16.0D0
      DG(21)=TERM*(EE+EQXX/4.0D0+QXXE/4.0D0+QXXQYY/16.0D0)
      DG(22)=TERM*(QXXQXX-QXXQYY)/32.0D0
      RETURN
      END
      SUBROUTINE DEPVAR (A,I,W,L)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(3,*)
C***********************************************************************
C
C  IN SUBROUTINE HADDON WHEN M, THE SYMMETRY OPERATION, IS 18 DEPVAR IS
C  CALLED. DEPVAR SHOULD THEN CONTAIN A USER-WRITTEN SYMMETRY OPERATION.
C  SEE HADDON TO GET THE IDEA ON HOW TO WRITE DEPVAR.
C
C ON INPUT:
C           A = ARRAY OF INTERNAL COORDINATES
C           I = ADDRESS OF REFERENCE ATOM
C ON OUTPUT:
C           L = 1 (IF A BOND-LENGTH IS THE DEPENDENT FUNCTION)
C             = 2 (IF AN ANGLE IS THE DEPENDENT FUNCTION)
C             = 3 (IF A DIHEDRAL ANGLE IS THE DEPENDENT FUNCTION)
C           W = VALUE OF THE FUNCTION
C
C  NOTE:  IT IS THE WRITER'S RESPONSIBILITY TO MAKE CERTAIN THAT THE
C         SUBROUTINE DOES NOT CONTAIN ANY ERRORS!
C***********************************************************************
      COMMON /KEYWRD/ KEYWRD
      COMMON /NUMCAL/ NUMCAL
      CHARACTER*80 KEYWRD
         SAVE                                                           GL0892
      DATA ICALCN /0/
      IF (ICALCN .NE. NUMCAL) THEN
         ICALCN = NUMCAL
         FACT=READA(KEYWRD,INDEX(KEYWRD,'DEPVAR'))
         WRITE(6,'(''  UNIT CELL LENGTH ='',F14.7,
     1'' TIMES BOND LENGTH'')')FACT
      ENDIF
      W=A(1,I)*FACT
      L=1
      RETURN
      END
      SUBROUTINE DFPSAV(TOTIME,XPARAM,GD,XLAST,FUNCT1,MDFP,XDFP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
************************************************************************
      DIMENSION XPARAM(*), GD(*), XLAST(*), MDFP(9),XDFP(9)
**********************************************************************
*
* DFPSAV STORES AND RESTORES DATA USED IN THE D-F-P GEOMETRY
*        OPTIMISATION.
*
*  ON INPUT TOTIME = TOTAL CPU TIME ELAPSED DURING THE CALCULATION.
*           XPARAM = CURRENT VALUE OF PARAMETERS.
*           GD     = OLD GRADIENT.
*           XLAST  = OLD VALUE OF PARAMETERS.
*           FUNCT1 = CURRENT VALUE OF HEAT OF FORMATION.
*           MDFP   = INTEGER CONSTANTS USED IN D-F-P.
*           XDFP   = REAL CONSTANTS USED IN D-F-P.
*           MDFP(9)= 1 FOR DUMP, 0 FOR RESTORE.
**********************************************************************
      COMMON /KEYWRD/ KEYWRD
      COMMON /TITLES/ KOMENT,TITLE
      COMMON /GRADNT/ GRAD(MAXPAR),GNORM
      COMMON /GEOVAR/ XDUMMY(MAXPAR),NVAR,LOC(2,MAXPAR)                 IR0394
      COMMON /DENSTY/ P(MPACK), PA(MPACK), PB(MPACK)
      COMMON /ALPARM/ ALPARM(3,MAXPAR),X0, X1, X2, ILOOP
      COMMON /REACTN/ STEP, GEOA(3,NUMATM), GEOVEC(3,NUMATM),CALCST
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1                NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /ELEMTS/ ELEMNT(107)
C    changed common path for portability  (IR)
      COMMON /PATHI / LATOM,LPARAM
      COMMON /PATHR / REACT(200)
C
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /FMATRX/ HESINV(MAXHES)
      COMMON /GEOSYM/ NDEP,LOCPAR(MAXPAR),IDEPFN(MAXPAR),
     1                     LOCDEP(MAXPAR)
      COMMON /ERRFN / ERRFN(MAXPAR)
      COMMON /NUMCAL/ NUMCAL
      DIMENSION IEL1(3),Q(3), COORD(3,NUMATM)
      CHARACTER ELEMNT*2, KEYWRD*80,KOMENT*80, TITLE*80
      LOGICAL INTXYZ
         SAVE                                                           GL0892
      DATA ICALCN /0/
      OPEN(UNIT=9,FILE='FOR009',STATUS='UNKNOWN',FORM='UNFORMATTED')
      REWIND 9
      OPEN(UNIT=10,FILE='FOR010',STATUS='UNKNOWN',FORM='UNFORMATTED')
      REWIND 10
      DEGREE=57.29577951D0
      IR=9
      IF(MDFP(9) .NE. 0) THEN
         IF(MDFP(9) .EQ. 1) THEN
            WRITE(6,'(//10X,''- - - - - - - TIME UP - - - - - - -'',//)'
     1)
            IF(INDEX(KEYWRD,'SADDLE') .NE. 0) THEN
               WRITE(6,'(//10X,'' NO RESTART EXISTS FOR SADDLE'',//
     1  10X,'' HERE IS A DATA-FILE FILES THAT MIGHT BE SUITABLE'',/
     2  10X,'' FOR RESTARTING THE CALCULATION'',///)')
               WRITE(6,'(A)')KEYWRD,KOMENT,TITLE
               INTXYZ=(NA(1).EQ.0)
               DO 60 ILOOP=1,2
                  IF(INTXYZ)THEN
                     GEO(2,1)=0.D0
                     GEO(3,1)=0.D0
                     GEO(1,1)=0.D0
                     GEO(2,2)=0.D0
                     GEO(3,2)=0.D0
                     GEO(3,3)=0.D0
                     DO 10 I=1,NATOMS
                        DO 10 J=1,3
   10                COORD(J,I)=GEO(J,I)
                  ELSE
                     CALL XYZINT(GEO,NUMAT,NA,NB,NC,1.D0,COORD)
                  ENDIF
                  IVAR=1
                  NA(1)=0
                  DO 40 I=1,NATOMS
                     DO 20 J=1,3
   20                IEL1(J)=0
   30                CONTINUE
                     IF(LOC(1,IVAR).EQ.I) THEN
                        IEL1(LOC(2,IVAR))=1
                        IVAR=IVAR+1
                        GOTO 30
                     ENDIF
                     IF(I.LT.4) THEN
                        IEL1(3)=0
                        IF(I.LT.3) THEN
                           IEL1(2)=0
                           IF(I.LT.2) THEN
                              IEL1(1)=0
                           ENDIF
                        ENDIF
                     ENDIF
                     IF(I.EQ.LATOM)IEL1(LPARAM)=-1
                     Q(1)=COORD(1,I)
                     Q(2)=COORD(2,I)*DEGREE
                     Q(3)=COORD(3,I)*DEGREE
   40             WRITE(6,'(2X,A2,3(F12.6,I3),I4,2I3)')
     1    ELEMNT(LABELS(I)),(Q(K),IEL1(K),K=1,3),NA(I),NB(I),NC(I)
                  I=0
                  X=0.D0
                  WRITE(6,'(I4,3(F12.6,I3),I4,2I3)')
     1    I,X,I,X,I,X,I,I,I,I
                  DO 50 I=1,NATOMS
                     DO 50 J=1,3
   50             GEO(J,I)=GEOA(J,I)
                  NA(1)=99
   60          CONTINUE
               WRITE(6,'(///10X,''CALCULATION TERMINATED HERE'')')
               STOP
            ENDIF
            WRITE(6,'(//10X,'' - THE CALCULATION IS BEING DUMPED TO DISK
     1'',  /10X,''   RESTART IT USING THE MAGIC WORD "RESTART"'')')
            WRITE(6,'(//10X,''CURRENT VALUE OF HEAT OF FORMATION =''
     1  ,F12.6)')FUNCT1
         ENDIF
         IF(MDFP(9) .EQ. 1)THEN
            IF(NA(1) .EQ. 99) THEN
C
C  CONVERT FROM CARTESIAN COORDINATES TO INTERNAL
C
               DO 70 I=1,NATOMS
                  DO 70 J=1,3
   70          COORD(J,I)=GEO(J,I)
               CALL XYZINT(COORD,NUMAT,NA,NB,NC,1.D0,GEO)
            ENDIF
            GEO(2,1)=0.D0
            GEO(3,1)=0.D0
            GEO(1,1)=0.D0
            GEO(2,2)=0.D0
            GEO(3,2)=0.D0
            GEO(3,3)=0.D0
            IVAR=1
            NA(1)=0
            WRITE(6,'(A)')KEYWRD,KOMENT,TITLE
            DO 100 I=1,NATOMS
               DO 80 J=1,3
   80          IEL1(J)=0
   90          CONTINUE
               IF(LOC(1,IVAR).EQ.I) THEN
                  IEL1(LOC(2,IVAR))=1
                  IVAR=IVAR+1
                  GOTO 90
               ENDIF
               IF(I.LT.4) THEN
                  IEL1(3)=0
                  IF(I.LT.3) THEN
                     IEL1(2)=0
                     IF(I.LT.2) THEN
                        IEL1(1)=0
                     ENDIF
                  ENDIF
               ENDIF
               IF(I.EQ.LATOM)IEL1(LPARAM)=-1
               Q(1)=GEO(1,I)
               Q(2)=GEO(2,I)*DEGREE
               Q(3)=GEO(3,I)*DEGREE
  100       WRITE(6,'(2X,A2,3(F12.6,I3),I4,2I3)')
     1ELEMNT(LABELS(I)),(Q(K),IEL1(K),K=1,3),NA(I),NB(I),NC(I)
            I=0
            X=0.D0
            WRITE(6,'(I4,3(F12.6,I3),I4,2I3)')
     1I,X,I,X,I,X,I,I,I,I
            IF(NDEP.NE.0)THEN
               DO 110 I=1,NDEP
  110          WRITE(6,'(3(I4,'',''))')LOCPAR(I),IDEPFN(I),LOCDEP(I)
               WRITE(6,*)
            ENDIF
         ENDIF
         WRITE(IR)MDFP,XDFP,TOTIME,FUNCT1
         WRITE(IR)(XPARAM(I),I=1,NVAR),(GD(I),I=1,NVAR)
         WRITE(IR)(XLAST(I),I=1,NVAR),(GRAD(I),I=1,NVAR)
         LINEAR=(NVAR*(NVAR+1))/2
         WRITE(IR)(HESINV(I),I=1,LINEAR)
         LINEAR=(NORBS*(NORBS+1))/2
         WRITE(10)(PA(I),I=1,LINEAR)
         IF(NALPHA.NE.0)WRITE(10)(PB(I),I=1,LINEAR)
         IF(LATOM .NE. 0) THEN
            WRITE(IR)((ALPARM(J,I),J=1,3),I=1,NVAR)
            WRITE(IR)ILOOP,X0, X1, X2
         ENDIF
         WRITE(IR)(ERRFN(I),I=1,NVAR)
         CLOSE (9)
         CLOSE (10)
      ELSE
         IF (ICALCN .NE. NUMCAL) THEN
           WRITE(6,'(//10X,'' RESTORING DATA FROM DISK''/)')
         END IF
         READ(IR,END=130,ERR=130)MDFP,XDFP,TOTIME,FUNCT1
         IF (ICALCN .NE. NUMCAL) THEN
           WRITE(6,'(10X,''FUNCTION ='',F13.6//)')FUNCT1
         END IF
         READ(IR)(XPARAM(I),I=1,NVAR),(GD(I),I=1,NVAR)
         READ(IR)(XLAST(I),I=1,NVAR),(GRAD(I),I=1,NVAR)
         LINEAR=(NVAR*(NVAR+1))/2
         READ(IR)(HESINV(I),I=1,LINEAR)
         LINEAR=(NORBS*(NORBS+1))/2
         READ(10)(PA(I),I=1,LINEAR)
         IF(NALPHA.NE.0)READ(10)(PB(I),I=1,LINEAR)
         IF(LATOM.NE.0) THEN
            READ(IR)((ALPARM(J,I),J=1,3),I=1,NVAR)
            READ(IR)ILOOP,X0, X1, X2
         ENDIF
         READ(IR)(ERRFN(I),I=1,NVAR)
  120    ICALCN = NUMCAL
         RETURN
  130    WRITE(6,'(//10X,''NO RESTART FILE EXISTS!'')')
         STOP
      ENDIF
      END
      FUNCTION DIPOLE (P,Q,COORD,DIPVEC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /KEYWRD/ KEYWRD
      COMMON /ISTOPE/ AMS(107)
      COMMON /MULTIP/ DD(107), QQ(107), AM(107), AD(107), AQ(107)
      COMMON /NUMCAL/ NUMCAL
      DIMENSION P(*),Q(*),COORD(3,*),DIPVEC(3), CENTER(3)
      CHARACTER*80 KEYWRD
C***********************************************************************
C     DIPOLE CALCULATES DIPOLE MOMENTS
C
C  ON INPUT P     = DENSITY MATRIX
C           Q     = TOTAL ATOMIC CHARGES, (NUCLEAR + ELECTRONIC)
C           NUMAT = NUMBER OF ATOMS IN MOLECULE
C           NAT   = ATOMIC NUMBERS OF ATOMS
C           NFIRST= START OF ATOM ORBITAL COUNTERS
C           COORD = COORDINATES OF ATOMS
C
C  OUTPUT  DIPOLE = DIPOLE MOMENT
C***********************************************************************
C
C     IN THE ZDO APPROXIMATION, ONLY TWO TERMS ARE RETAINED IN THE
C     CALCULATION OF DIPOLE MOMENTS.
C     1. THE POINT CHARGE TERM (INDEPENDENT OF PARAMETERIZATION).
C     2. THE ONE-CENTER HYBRIDIZATION TERM, WHICH ARISES FROM MATRIX
C     ELEMENTS OF THE FORM <NS/R/NP>. THIS TERM IS A FUNCTION OF
C     THE SLATER EXPONENTS (ZS,ZP) AND IS THUS DEPENDENT ON PARAMETER-
C     IZATION. THE HYBRIDIZATION FACTORS (HYF(I)) USED IN THIS SUB-
C     ROUTINE ARE CALCULATED FROM THE FOLLOWING FORMULAE.
C     FOR SECOND ROW ELEMENTS <2S/R/2P>
C     HYF(I)= 469.56193322*(SQRT(((ZS(I)**5)*(ZP(I)**5)))/
C           ((ZS(I) + ZP(I))**6))
C     FOR THIRD ROW ELEMENTS <3S/R/3P>
C     HYF(I)=2629.107682607*(SQRT(((ZS(I)**7)*(ZP(I)**7)))/
C           ((ZS(I) + ZP(I))**8))
C     FOR FOURTH ROW ELEMENTS AND UP :
C     HYF(I)=2*(2.10716)*DD(I)
C     WHERE DD(I) IS THE CHARGE SEPARATION IN ATOMIC UNITS
C
C
C     REFERENCES:
C     J.A.POPLE & D.L.BEVERIDGE: APPROXIMATE M.O. THEORY
C     S.P.MCGLYNN, ET AL: APPLIED QUANTUM CHEMISTRY
C
      DIMENSION DIP(4,3)
      DIMENSION HYF(107,2)
      LOGICAL FORCE, CHARGD
         SAVE                                                           GL0892
      DATA HYF(1,1)     / 0.0D00           /
      DATA   HYF(1,2) /0.0D0     /
      DATA   HYF(5,2) /6.520587D0/
      DATA   HYF(6,2) /4.253676D0/
      DATA   HYF(7,2) /2.947501D0/
      DATA   HYF(8,2) /2.139793D0/
      DATA   HYF(9,2) /2.2210719D0/
      DATA   HYF(14,2)/6.663059D0/
      DATA   HYF(15,2)/5.657623D0/
      DATA   HYF(16,2)/6.345552D0/
      DATA   HYF(17,2)/2.522964D0/
      DATA ICALCN /0/
      IF (ICALCN .NE. NUMCAL) THEN
         DO 10 I=2,107
   10    HYF(I,1)= 5.0832*DD(I)
         WTMOL=0.D0
         SUM=0.D0
         DO 20 I=1,NUMAT
            WTMOL=WTMOL+AMS(NAT(I))
   20    SUM=SUM+Q(I)
         CHARGD=(ABS(SUM).GT.0.5D0)
         ICALCN =NUMCAL
         FORCE=(INDEX(KEYWRD,'FORCE') +INDEX(KEYWRD,'IRC').NE. 0)
         ITYPE=1
         IF(INDEX(KEYWRD,'MINDO') .NE. 0)ITYPE=2
      ENDIF
      IF(CHARGD)THEN
C
C   NEED TO RESET ION'S POSITION SO THAT THE CENTER OF MASS IS AT THE
C   ORIGIN.
C
         DO 30 I=1,3
   30    CENTER(I)=0.D0
         DO 40 I=1,3
            DO 40 J=1,NUMAT
   40    CENTER(I)=CENTER(I)+AMS(NAT(J))*COORD(I,J)
         DO 50 I=1,3
   50    CENTER(I)=CENTER(I)/WTMOL
         DO 60 I=1,3
            DO 60 J=1,NUMAT
   60    COORD(I,J)=COORD(I,J)-CENTER(I)
      ENDIF
      DO 70 I=1,4
         DO 70 J=1,3
   70 DIP(I,J)=0.0D00
      DO 90 I=1,NUMAT
         NI=NAT(I)
         IA=NFIRST(I)
         L=NLAST(I)-IA
         DO 80 J=1,L
            K=((IA+J)*(IA+J-1))/2+IA
   80    DIP(J,2)=DIP(J,2)-HYF(NI,ITYPE)*P(K)
         DO 90 J=1,3
   90 DIP(J,1)=DIP(J,1)+4.803D00*Q(I)*COORD(J,I)
      DO 100 J=1,3
  100 DIP(J,3)=DIP(J,2)+DIP(J,1)
      DO 110 J=1,3
  110 DIP(4,J)=SQRT(DIP(1,J)**2+DIP(2,J)**2+DIP(3,J)**2)
      IF( FORCE) THEN
         DIPVEC(1)=DIP(1,3)
         DIPVEC(2)=DIP(2,3)
         DIPVEC(3)=DIP(3,3)
      ELSE
         WRITE (6,120) ((DIP(I,J),I=1,4),J=1,3)
      ENDIF
      DIPOLE = DIP(4,3)
      RETURN
C
  120 FORMAT (' DIPOLE',11X,2HX ,8X,2HY ,8X,2HZ ,6X,'TOTAL',/,
     1' POINT-CHG.',4F10.3/,' HYBRID',4X,4F10.3/,' SUM',7X,4F10.3)
C
      END
*************************************************************
*                                                           *
*    DRCOUT PRINTS THE GEOMETRY, ETC. FOR A DRC AT A        *
*    POSITION DETERMINED BY FRACT.                          *
*    ON INPUT XYZ3  = QUADRATIC EXPRESSION FOR THE GEOMETRY *
*             VEL3  = QUADRATIC EXPRESSION FOR THE VELOCITY *
*             ESCF3 = QUADRATIC EXPRESSION FOR THE P.E.     *
*             EKIN3 = QUADRATIC EXPRESSION FOR THE K.E.     *
*                                                           *
*************************************************************
      SUBROUTINE DRCOUT(XYZ3,GEO3,VEL3,NVAR,TIME,ESCF3,EKIN3,
     1GTOT3,ETOT3,XTOT3,ILOOP,CHARGE,FRACT,TEXT1,TEXT2,II,JLOOP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      CHARACTER TEXT1*3,TEXT2*2
      DIMENSION XYZ3(3,NVAR), VEL3(3,NVAR), ESCF3(3), EKIN3(3)
      DIMENSION XYZ(3,NUMATM), VEL(3,NUMATM), CHARGE(NUMATM),
     1 GEO3(3,NUMATM), ETOT3(3), GTOT3(3), XTOT3(3)
      COMMON /KEYWRD/ KEYWRD
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1                NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /ELEMTS/ ELEMNT(107)
      COMMON /GEOVAR/ XPARAM(MAXPAR),IDUM, LOC(2,MAXPAR)                IR0394
      COMMON /TITLES/ KOMENT,TITLE
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /NUMCAL/ NUMCAL
      DIMENSION  IEL1(3), GG(3)
      CHARACTER*80 KEYWRD,KOMENT,TITLE, ALPHA*2, ELEMNT*2
      INTEGER PRTKOM, PRTITL, PRTKEY
      LOGICAL LARGE, DRC
         SAVE                                                           GL0892
      DATA ICALCN /0/
      IF (ICALCN  .NE. NUMCAL) THEN
         ICALCN = NUMCAL
         IF( INDEX(KEYWRD,'REST').EQ.0.OR.INDEX(KEYWRD,'IRC=').NE.0)THEN
            JLOOP=0
         ENDIF
         DO 10 I=80,1,-1
   10    IF(KEYWRD(I:I).NE.' ')GOTO 20
         I=1
   20    PRTKEY=I
         DO 30 I=80,1,-1
   30    IF(KOMENT(I:I).NE.' ')GOTO 40
         I=1
   40    PRTKOM=I
         DO 50 I=80,1,-1
   50    IF(TITLE(I:I).NE.' ')GOTO 60
         I=1
   60    PRTITL=I
         DRC=(INDEX(KEYWRD,'DRC').NE.0)
         LARGE=(INDEX(KEYWRD,'LARGE').NE.0)
         IF(DRC) THEN
            WRITE(6,'(//,'' FEMTOSECONDS  POINT  POTENTIAL +''
     1,'' KINETIC  =  TOTAL     ERROR    REF%'')')
         ELSE
            WRITE(6,'(//,''     POINT   POTENTIAL  +  ''
     1,''ENERGY LOST   =   TOTAL      ERROR    REF%'')')
         ENDIF
      ELSE
         IF(DRC) THEN
            WRITE(6,'(//,'' FEMTOSECONDS  POINT  POTENTIAL +''
     1,'' KINETIC  =  TOTAL     ERROR    REF'')')
         ELSE
            WRITE(6,'(//,''     POINT   POTENTIAL  +  ''
     1,''ENERGY LOST   =   TOTAL      ERROR    REF'')')
         ENDIF
      ENDIF
      JLOOP=JLOOP+1
C#      FRACT=0.D0
      ESCF=ESCF3(1)+ESCF3(2)*FRACT+ESCF3(3)*FRACT**2
      EKIN=EKIN3(1)+EKIN3(2)*FRACT+EKIN3(3)*FRACT**2
      ETOT=ETOT3(1)+ETOT3(2)*FRACT+ETOT3(3)*FRACT**2
      GTOT=GTOT3(1)+GTOT3(2)*FRACT+GTOT3(3)*FRACT**2
      XTOT=XTOT3(1)+XTOT3(2)*FRACT+XTOT3(3)*FRACT**2
      IF(II.NE.0)THEN
         IF(DRC) THEN
            WRITE(6,'(F10.3,I8,F12.5,F11.5,F11.5,
     1F10.5,'' '',I5,3X,''%'',A,A,I3)')TIME, ILOOP-2, ESCF, EKIN,
     2 ESCF+EKIN,ESCF+EKIN-ETOT,JLOOP,TEXT1,TEXT2,II
            WRITE(6,'(9X,A,F9.4,A)')' MOVEMENT FROM START =',XTOT,' ANGS
     1TROMS'
         ELSE
            WRITE(6,'(I8,F14.5,F13.5,F17.5,
     1F10.5,'' '',I5,3X,''%'',A,A,I3)') ILOOP-2, ESCF, EKIN,
     2 ESCF+EKIN,ESCF+EKIN-ETOT,JLOOP,TEXT1,TEXT2,II
C#      WRITE(6,'(24X,'' INTEGRATED GRADIENT ERROR ='',G10.3,
C#     1'' KCALS/ANGSTROM'')')GTOT
            WRITE(6,'(9X,A,F9.4,A)')' MOVEMENT FROM START =',XTOT,' ANGS
     1TROMS'
         ENDIF
      ELSE
         IF(DRC) THEN
            WRITE(6,'(F10.3,I8,F12.5,F11.5,F11.5,
     1F10.5,'' '',I5,3X,''%'',A,A,I3)')TIME, ILOOP-2, ESCF, EKIN,
     2 ESCF+EKIN,ESCF+EKIN-ETOT,JLOOP,TEXT1,TEXT2
            WRITE(6,'(9X,A,F9.4,A)')' MOVEMENT FROM START =',XTOT,' ANGS
     1TROMS'
         ELSE
            WRITE(6,'(I8,F14.5,F13.5,F17.5,
     1F10.5,'' '',I5,3X,''%'',A,A,I3)') ILOOP-2, ESCF, EKIN,
     2 ESCF+EKIN,ESCF+EKIN-ETOT,JLOOP,TEXT1,TEXT2
C#      WRITE(6,'(24X,'' INTEGRATED GRADIENT ERROR ='',F10.5,
C#     1'' KCALS/ANGSTROM'')')GTOT
            WRITE(6,'(9X,A,F9.4,A)')' MOVEMENT FROM START =',XTOT,' ANGS
     1TROMS'
         ENDIF
      ENDIF
      NATOMS=NVAR/3
      L=0
      DO 80 I=1,NATOMS
         DO 70 J=1,3
            L=L+1
            VEL(J,I)=VEL3(1,L)+VEL3(2,L)*FRACT+VEL3(3,L)*FRACT**2
   70    XYZ(J,I)=XYZ3(1,L)+XYZ3(2,L)*FRACT+XYZ3(3,L)*FRACT**2
   80 CONTINUE
      IF(LARGE)THEN
         WRITE(6,*)'                CARTESIAN GEOMETRY           '
     1//'VELOCITY (IN CM/SEC)'
         WRITE(6,*)'  ATOM        X          Y          Z'
     1      //'                X          Y          Z'
         DO 90 I=1,NUMAT
            LL=(I-1)*3+1
            LU=LL+2
            WRITE(6,'(I4,3X,A2,3F11.5,2X,3F11.1)')
     1I, ELEMNT(NAT(I)),(XYZ(J,I),J=1,3),(VEL(J,I),J=1,3)
   90    CONTINUE
      ENDIF
      IVAR=1
      NA(1)=0
      L=0
      WRITE(6,'(//10X,''FINAL GEOMETRY OBTAINED'',33X,''CHARGE'')'
     1)
      WRITE(6,'(A)')KEYWRD(:PRTKEY),KOMENT(:PRTKOM),TITLE(:PRTITL)
      LMAX=L
      L=0
      DO 120 I=1,NUMAT
         J=I/26
         ALPHA(1:1)=CHAR(ICHAR('A')+J)
         J=I-J*26
         ALPHA(2:2)=CHAR(ICHAR('A')+J-1)
         DO 100 J=1,3
  100    IEL1(J)=0
  110    CONTINUE
         IF(LOC(1,IVAR).EQ.I) THEN
            IEL1(LOC(2,IVAR))=1
            IVAR=IVAR+1
            GOTO 110
         ENDIF
         IF(I.LT.4) THEN
            IEL1(3)=0
            IF(I.LT.3) THEN
               IEL1(2)=0
               IF(I.LT.2) THEN
                  IEL1(1)=0
               ENDIF
            ENDIF
         ENDIF
         IF(LABELS(I).LT.99)THEN
            L=L+1
            GG(1)=GEO3(1,I*3-2)+GEO3(2,I*3-2)*FRACT+GEO3(3,I*3-2)*FRACT*
     1*2
            GG(2)=GEO3(1,I*3-1)+GEO3(2,I*3-1)*FRACT+GEO3(3,I*3-1)*FRACT*
     1*2
            GG(3)=GEO3(1,I*3  )+GEO3(2,I*3  )*FRACT+GEO3(3,I*3  )*FRACT*
     1*2
            WRITE(6,'(2X,A2,3(F12.6,I3),I4,2I3,F13.4,I5,A)')
     1    ELEMNT(LABELS(I)),(GG(K),IEL1(K),K=1,3),
     2    NA(I),NB(I),NC(I),CHARGE(L),JLOOP,ALPHA//'*'
         ELSE
            WRITE(6,'(2X,A2,3(F12.6,I3),I4,2I3,13X,I5,A)')
     1    ELEMNT(LABELS(I)),(GG(K),IEL1(K),K=1,3),
     2NA(I),NB(I),NC(I),JLOOP,ALPHA//'%'
         ENDIF
  120 CONTINUE
      NA(1)=99
      RETURN
      END
      SUBROUTINE FOCK1(F, PTOT, PA, PB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      COMMON /NUMCAL/ NUMCAL
      DIMENSION F(*), PTOT(*), PA(*), PB(*)
C *********************************************************************
C
C *** COMPUTE THE REMAINING CONTRIBUTIONS TO THE ONE-CENTRE ELEMENTS.
C
C *********************************************************************
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /GAUSS / FN1(107),FN2(107)
     1       /MOLORB/ USPD(MAXORB),DUMY(MAXORB)
      COMMON /TWOELE/ GSS(107),GSP(107),GPP(107),GP2(107),HSP(107)
     1                ,GSD(107),GPD(107),GDD(107)
      DIMENSION QTOT(NUMATM), QA(NUMATM), QB(NUMATM)
      COMMON /KEYWRD/ KEYWRD
      CHARACTER*80 KEYWRD
         SAVE                                                           GL0892
      DATA ICALCN /0/
      IF (ICALCN .NE. NUMCAL) THEN
         ICALCN = NUMCAL
      ENDIF
      CALL CHRGE(PTOT,QTOT)
      CALL CHRGE(PA,QA)
      DO 10 I=1,NUMAT
   10 QB(I)=QTOT(I)-QA(I)
      DO 100 II=1,NUMAT
         IA=NFIRST(II)
         IB=NMIDLE(II)
         IC=NLAST(II)
         NI=NAT(II)
         DTPOP=0.D0
         DAPOP=0.D0
         PTPOP=0.D0
         PAPOP=0.D0
         GOTO (100,40,30,30,30,20,20,20,20,20)IC-IA+2
   20    DTPOP=PTOT((IC*(IC+1))/2)+PTOT(((IC-1)*(IC))/2)
     1        +PTOT(((IC-2)*(IC-1))/2)+PTOT(((IC-3)*(IC-2))/2)
     2        +PTOT(((IC-4)*(IC-3))/2)
         DAPOP=PA((IC*(IC+1))/2)+PA(((IC-1)*(IC))/2)
     1        +PA(((IC-2)*(IC-1))/2)+PA(((IC-3)*(IC-2))/2)
     2        +PA(((IC-4)*(IC-3))/2)
   30    PTPOP=PTOT((IB*(IB+1))/2)+PTOT(((IB-1)*(IB))/2)
     1        +PTOT(((IB-2)*(IB-1))/2)
         PAPOP=PA((IB*(IB+1))/2)+PA(((IB-1)*(IB))/2)
     1        +PA(((IB-2)*(IB-1))/2)
   40    DBPOP=DTPOP-DAPOP
         PBPOP=PTPOP-PAPOP
         IF(NI.EQ.1)THEN
            SUM=0.D0
         ELSE
            SUM2=0.D0
            SUM1=0.D0
            DO 60 I=IA,IB
               IM1=I-1
               DO 50 J=IA,IM1
   50          SUM1=SUM1+PTOT(J+(I*(I-1))/2)**2
   60       SUM2=SUM2+PTOT((I*(I+1))/2)**2
            SUM=SUM1*2.D0+SUM2
            SUM=SQRT(SUM)-QTOT(II)*0.5D0
         ENDIF
         SUM=SUM*FN1(NI)
C
C     F(S,S)
C
         KA=(IA*(IA+1))/2
         F(KA)=F(KA)+PB(KA)*GSS(NI)+PTPOP*GSP(NI)
     1         -PAPOP*HSP(NI) + DTPOP*GSD(NI)
         IF (NI.LT.3) GO TO 100
         IPLUS=IA+1
         L=KA
         DO 70 J=IPLUS,IB
            M=L+IA
            L=L+J
C
C     F(P,P)
C
            F(L)=F(L)+PTOT(KA)*GSP(NI)-PA(KA)*HSP(NI)+
     1      PB(L)*GPP(NI)+(PTPOP-PTOT(L))*GP2(NI)
     2      -0.5D0*(PAPOP-PA(L))*(GPP(NI)-GP2(NI))
     3      +DTPOP*GPD(NI)
C
C     F(S,P)
C
   70    F(M)=F(M)+2.D0*PTOT(M)*HSP(NI)-PA(M)*(HSP(NI)+GSP(NI))
C
C     F(P,P*)
C
         IMINUS=IB-1
         DO 80 J=IPLUS,IMINUS
            ICC=J+1
            DO 80 L=ICC,IB
               M=(L*(L-1))/2+J
   80    F(M)=F(M)+PTOT(M)*(GPP(NI)-GP2(NI))
     1      -0.5D0*PA  (M)*(GPP(NI)+GP2(NI))
         DO 90 J=IB+1,IC
            M=(J*(J+1))/2
   90    F(M)=F(M)+PTOT(KA)*GSD(NI)
     1         +PTPOP*GPD(NI)
     2         +(DTPOP-PA(M))*GDD(NI)
  100 CONTINUE
      RETURN
      END
      SUBROUTINE GETGEO(IREAD,LABELS,GEO,LOPT,NA,NB,NC,AMS,NATOMS,INT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION GEO(3,*),NA(*),NB(*),NC(*),AMS(*), LOPT(3,*)
     1,LABELS(*)
      LOGICAL INT
************************************************************************
*
*   GETGEO READS IN THE GEOMETRY. THE ELEMENT IS SPECIFIED BY IT'S
*          CHEMICAL SYMBOL, OR, OPTIONALLY, BY IT'S ATOMIC NUMBER.
*
*  ON INPUT   IREAD  = CHANNEL NUMBER FOR READ, NORMALLY 5
*             AMS    = DEFAULT ATOMIC MASSES.
*
* ON OUTPUT LABELS = ATOMIC NUMBERS OF ALL ATOMS, INCLUDING DUMMIES.
*           GEO    = INTERNAL COORDINATES, IN ANGSTROMS, AND DEGREES.
*           LOPT   = INTEGER ARRAY, A '1' MEANS OPTIMIZE THIS PARAMETER,
*                    '0' MEANS DO NOT OPTIMIZE, AND A '-1' LABELS THE
*                    REACTION COORDINATE.
*           NA     = INTEGER ARRAY OF ATOMS (SEE DATA INPUT)
*           NB     = INTEGER ARRAY OF ATOMS (SEE DATA INPUT)
*           NC     = INTEGER ARRAY OF ATOMS (SEE DATA INPUT)
*           ATMASS = ATOMIC MASSES OF ATOMS.
************************************************************************
      COMMON /PATHR / REACT(3,66), DUMM1,DUMM2                          IR0394
      COMMON /ATMASS/ ATMASS(NUMATM)
      COMMON /KEYWRD/ KEYWRD
      DIMENSION ISTART(40), XYZ(3,NUMATM), VALUE(4)
      LOGICAL LEADSP, IRCDRC
      CHARACTER*80 KEYWRD
      CHARACTER ELEMNT(107)*2, LINE*80, SPACE*1, NINE*1,ZERO*1,
     1TAB*1, COMMA*1, STRING*80, ELE*2, TURN*1
         SAVE                                                           GL0892
      DATA (ELEMNT(I),I=1,107)/'H','HE',
     1 'LI','BE','B','C','N','O','F','NE',
     2 'NA','MG','AL','SI','P','S','CL','AR',
     3 'K','CA','SC','TI','V','CR','MN','FE','CO','NI','CU',
     4 'ZN','GA','GE','AS','SE','BR','KR',
     5 'RB','SR','Y','ZR','NB','MO','TC','RU','RH','PD','AG',
     6 'CD','IN','SN','SB','TE','I','XE',
     7 'CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB','DY',
     8 'HO','ER','TM','YB','LU','HF','TA','W','RE','OS','IR','PT',
     9 'AU','HG','TL','PB','BI','PO','AT','RN',
     1 'FR','RA','AC','TH','PA','U','NP','PU','AM','CM','BK','CF','XX',
     2 'FM','MD','CB','++','+','--','-','TV'/
      DATA COMMA,SPACE,NINE,ZERO/',',' ','9','0'/
      TAB=CHAR(9)
      IRCDRC=(INDEX(KEYWRD,'IRC')+INDEX(KEYWRD,'DRC') .NE.0)
      ILOWA = ICHAR('a')
      ILOWZ = ICHAR('z')
      ICAPA = ICHAR('A')
      NATOMS=0
      NUMAT=0
   10 READ(IREAD,'(A)',END=90,ERR=180)LINE
      IF(LINE.EQ.' ') GO TO 90
      NATOMS=NATOMS+1
      IF(NATOMS.GT.NUMATM)THEN
         WRITE(6,'(//10X,''****  MAX. NUMBER OF ATOMS ALLOWED:'',I4)')
     1NUMATM
         STOP
      ENDIF
*   CLEAN THE INPUT DATA
************************************************************************
      DO 20 I=1,80
         ILINE=ICHAR(LINE(I:I))
         IF(ILINE.GE.ILOWA.AND.ILINE.LE.ILOWZ) THEN
            LINE(I:I)=CHAR(ILINE+ICAPA-ILOWA)
         ENDIF
   20 CONTINUE
************************************************************************
      ICOMMA=ICHAR(COMMA)
      ITAB=ICHAR(TAB)
      DO 30 I=1,80
      KHAR=ICHAR(LINE(I:I))
      IF(KHAR.EQ.ICOMMA.OR.KHAR.EQ.ITAB)LINE(I:I)=SPACE
  30  CONTINUE
*
*   INITIALIZE ISTART TO INTERPRET BLANKS AS ZERO'S
      DO 40 I=1,10
   40 ISTART(I)=80
* FIND INITIAL DIGIT OF ALL NUMBERS, CHECK FOR LEADING SPACES FOLLOWED
*     BY A CHARACTER AND STORE IN ISTART
      LEADSP=.TRUE.
      NVALUE=0
      DO 50 I=1,80
         IF (LEADSP.AND.LINE(I:I).NE.SPACE) THEN
            NVALUE=NVALUE+1
            ISTART(NVALUE)=I
         END IF
         LEADSP=(LINE(I:I).EQ.SPACE)
   50 CONTINUE
*
* ESTABLISH THE ELEMENT'S NAME AND ISOTOPE, CHECK FOR ERRORS OR E.O.DATA
*
      WEIGHT=0.D0
      STRING=LINE(ISTART(1):ISTART(2)-1)
      IF( STRING(1:1) .GE. ZERO .AND. STRING(1:1) .LE. NINE) THEN
*  ATOMIC NUMBER USED: NO ISOTOPE ALLOWED
         LABEL=READA(STRING,1)
         IF (LABEL.EQ.0) GO TO 80
         IF (LABEL.LT.0.OR.LABEL.GT.107) THEN
            WRITE(6,'(''  ILLEGAL ATOMIC NUMBER'')')
            GO TO 190
         END IF
         GO TO 70
      END IF
*  ATOMIC SYMBOL USED
      REAL=ABS(READA(STRING,1))
      IF (REAL.LT.1.D-15) THEN
*   NO ISOTOPE
         ELE=STRING(1:2)
      ELSE
         WEIGHT=REAL
         IF( STRING(2:2) .GE. ZERO .AND. STRING(2:2) .LE. NINE) THEN
            ELE=STRING(1:1)
         ELSE
            ELE=STRING(1:2)
         END IF
      END IF
*   CHECK FOR ERROR IN ATOMIC SYMBOL
      IF(ELE(1:1).EQ.'-'.AND.ELE(2:2).NE.'-')ELE(2:2)=' '
      DO 60 I=1,107
         IF(ELE.EQ.ELEMNT(I)) THEN
            LABEL=I
            GO TO 70
         END IF
   60 CONTINUE
      WRITE(6,'(''  UNRECOGNIZED ELEMENT NAME: ('',A,'')'')')ELE
      GOTO 190
*
* ALL O.K.
*
      WRITE(6,*) "AJith get HERE"
   70 IF (LABEL.NE.99) NUMAT=NUMAT+1
      IF(WEIGHT.NE.0.D0)THEN
         WRITE(6,'('' FOR ATOM'',I4,''  ISOTOPIC MASS:''
     1    ,F15.5)')NATOMS, WEIGHT
         ATMASS(NUMAT)=WEIGHT
      ELSE
         IF(LABEL .NE. 99)  ATMASS(NUMAT)=AMS(LABEL)
      ENDIF
      LABELS(NATOMS)   =LABEL
      GEO(1,NATOMS)    =READA(LINE,ISTART(2))
      GEO(2,NATOMS)    =READA(LINE,ISTART(4))
      GEO(3,NATOMS)    =READA(LINE,ISTART(6))
      IF(IRCDRC)THEN
         TURN=LINE(ISTART(3):ISTART(3))
         IF(TURN.EQ.'T')THEN
            LOPT(1,NATOMS)=1
            IF(NATOMS.EQ.1)WRITE(6,'(A)')' IN DRC MONITOR POTENTIAL ENER
     1GY'//' TURNING POINTS'
         ELSE
            LOPT(1,NATOMS)=0
         ENDIF
         TURN=LINE(ISTART(5):ISTART(5))
         IF(TURN.EQ.'T')THEN
            LOPT(2,NATOMS)=1
         ELSE
            LOPT(2,NATOMS)=0
         ENDIF
         TURN=LINE(ISTART(7):ISTART(7))
         IF(TURN.EQ.'T')THEN
            LOPT(3,NATOMS)=1
         ELSE
            LOPT(3,NATOMS)=0
         ENDIF
      ELSE
         LOPT(1,NATOMS)   =READA(LINE,ISTART(3))
         LOPT(2,NATOMS)   =READA(LINE,ISTART(5))
         LOPT(3,NATOMS)   =READA(LINE,ISTART(7))
      ENDIF
      NA(NATOMS)       =READA(LINE,ISTART(8))
      NB(NATOMS)       =READA(LINE,ISTART(9))
      NC(NATOMS)       =READA(LINE,ISTART(10))
      GOTO 10
*
* ALL DATA READ IN, CLEAN UP AND RETURN
*
   80 NATOMS=NATOMS-1
   90 NA(2)=1
      IF(NATOMS.GT.3)THEN
         INT=(NA(4).NE.0)
      ELSE
         IF(GEO(2,3).LT.10.AND.NATOMS.EQ.3)
     1WRITE(6,'(//10X,'' WARNING: INTERNAL COORDINATES ARE ASSUMED -'',/
     210X,'' FOR THREE-ATOM SYSTEMS '',//)')
         INT=.TRUE.
      ENDIF
C
C     READ IN VELOCITY VECTOR, IF PRESENT
C
      IF(INDEX(KEYWRD,'VELO').GT.0)THEN
         IF(INT)THEN
            WRITE(6,'(A)')' COORDINATES MUST BE CARTESIAN WHEN VELOCITY
     1VECTOR IS USED.'
            STOP
         ENDIF
C#      WRITE(6,'(/10X,A)')'INITIAL VELOCITY VECTOR FOR DRC'
         DO 110 I=1,NATOMS
            READ(5,'(A)') LINE
            CALL NUCHAR(LINE,VALUE,NDMY)
            IF(NDMY.NE.3)THEN
               WRITE(6,'(/10X,A)')
     1'  THERE MUST BE EXACTLY THREE VELOCITY DATA PER LINE'
               STOP
            ENDIF
            DO 100 J=1,3
  100       REACT(J,I+2)=VALUE(J)
C#      WRITE(6,'(2X,A2,2X,3F13.5)')ELEMNT(LABELS(I)),(VALUE(J),J=1,3)
  110    CONTINUE
         DO 120 I=1,3
            DO 120 J=1,2
  120    REACT(I,J)=GEO(I,J+1)-GEO(I,1)
C
C  NOW TO ROTATE VELOCITY VECTOR TO SUIT INTERNAL COORDINATE DEFINITION
C
C
C   ROTATE AROUND THE 1-2 X-AXIS TO AS TO ELIMINATE REACT(3,2)
C   (PUT ATOM 2 IN X-Y PLANE)
         SA=REACT(3,1)/SQRT(REACT(2,1)**2+REACT(3,1)**2+1.D-20)
         CA=SIGN(SQRT(1.D0-SA**2),REACT(2,1))
C#      LABELS(NATOMS+1)=1
C#      LABELS(NATOMS+2)=1
C#      WRITE(6,*)' FIRST ROTATION, ABOUT 1-2 X-AXIS'
         DO 130 I=1,NATOMS+2
            TEMP1= REACT(2,I)*CA+REACT(3,I)*SA
            TEMP2=-REACT(2,I)*SA+REACT(3,I)*CA
            REACT(2,I)=TEMP1
            REACT(3,I)=TEMP2
C#      WRITE(6,'(2X,A2,2X,3F13.5)')ELEMNT(LABELS(I)),(REACT(J,I),J=1,3)
  130    CONTINUE
C   ROTATE AROUND THE 1-2 Z-AXIS TO AS TO ELIMINATE REACT(2,2)
C   (PUT ATOM 2 ON X AXIS)
         CA=REACT(1,1)/SQRT(REACT(2,1)**2+REACT(1,1)**2+1.D-20)
         SA=SIGN(SQRT(1.D0-CA**2),REACT(2,1))
C#      WRITE(6,*)' SECOND ROTATION, ABOUT 1-2 Z-AXIS'
         DO 140 I=1,NATOMS+2
            TEMP1= REACT(1,I)*CA+REACT(2,I)*SA
            TEMP2=-REACT(1,I)*SA+REACT(2,I)*CA
            REACT(1,I)=TEMP1
            REACT(2,I)=TEMP2
C#      WRITE(6,'(2X,A2,2X,3F13.5)')ELEMNT(LABELS(I)),(REACT(J,I),J=1,3)
  140    CONTINUE
C   ROTATE AROUND THE 2-3 X-AXIS TO AS TO ELIMINATE REACT(3,3)
C   (PUT ATOM 3 ON X-Y PLANE)
         SA=REACT(3,2)/SQRT(REACT(2,2)**2+REACT(3,2)**2+1.D-20)
         CA=SIGN(SQRT(1.D0-SA**2),REACT(2,2))
C#      WRITE(6,*)' THIRD ROTATION, ABOUT 2-3 X-AXIS'
         DO 150 I=1,NATOMS+2
            TEMP1= REACT(2,I)*CA+REACT(3,I)*SA
            TEMP2=-REACT(2,I)*SA+REACT(3,I)*CA
            REACT(2,I)=TEMP1
            REACT(3,I)=TEMP2
C#      WRITE(6,'(2X,A2,2X,3F13.5)')ELEMNT(LABELS(I)),(REACT(J,I),J=1,3)
  150    CONTINUE
C
C  STRIP OFF FIRST TWO COORDINATES; THESE WERE THE COORDINATE AXIS
C  DEFINITIONS
C
         DO 160 I=1,NATOMS
            DO 160 J=1,3
  160    REACT(J,I)=REACT(J,I+2)
      ENDIF
      IF(  .NOT. INT ) THEN
         DO 170 I=1,NATOMS
            DO 170 J=1,3
  170    XYZ(J,I)=GEO(J,I)
         DEGREE=90.D0/ASIN(1.D0)
         CALL XYZINT(XYZ,NATOMS,NA,NB,NC,DEGREE,GEO)
         IF(ABS(GEO(2,3)-180.D0).LT.1.D-4.OR.ABS(GEO(2,3)).LT.1.D-4)
     1THEN
            WRITE(6,'(A)')' DUE TO PROGRAM BUG, THE FIRST THREE ATOMS MU
     1ST NOT LIE IN A STRAIGHT LINE.'
            STOP
         ENDIF
      ELSEIF (.NOT.IRCDRC) THEN
         IF(LOPT(1,1)+LOPT(2,1)+LOPT(3,1)+LOPT(2,2)+LOPT(3,2)+
     1        LOPT(3,3) .GT. 0)THEN
            LOPT(1,1)=0
            LOPT(2,1)=0
            LOPT(3,1)=0
            LOPT(2,2)=0
            LOPT(3,2)=0
            LOPT(3,3)=0
            WRITE(6,'(//10X,'' AN UNOPTIMIZABLE GEOMETRIC PARAMETER HAS'
     1',/10X,'' BEEN MARKED FOR OPTIMIZATION. THIS IS A NON-FATAL ''
     2,''ERROR'')')
         ENDIF
      ENDIF
      IF(NA(3).EQ.0) THEN
         NB(3)=1
         NA(3)=2
      ENDIF
      RETURN
* ERROR CONDITIONS
  180 IF(IREAD.EQ.5) THEN
         WRITE(6,'( '' ERROR DURING READ AT ATOM NUMBER '', I3 )')NATOMS
      ELSE
         NATOMS=0
         RETURN
      ENDIF
  190 J=NATOMS-1
      WRITE(6,'('' DATA CURRENTLY READ IN ARE '')')
      DO 200 K=1,J
  200 WRITE(6,210)LABELS(K),(GEO(JJ,K),LOPT(JJ,K),JJ=1,3),
     1NA(K),NB(K),NC(K)
  210 FORMAT(I4,2X,3(F10.5,2X,I2,2X),3(I2,1X))
      STOP
      END
      SUBROUTINE GETSYM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      COMMON /KEYWRD/ KEYWRD
      COMMON /GEOSYM/ NDEP, LOCPAR(MAXPAR), IDEPFN(MAXPAR),
     1                LOCDEP(MAXPAR)
      COMMON /IOCM/ IREAD
C***********************************************************************
C
C   GETSYM READS IN THE SYMMETRY DEPENDENCE RELATIONSHIPS.
C
C   ON EXIT     NDEP    = NUMBER OF SYMMETRY RELATIONS.
C               LOCPAR  = ARRAY OF REFERENCE FUNCTION INDICES.
C               IDEPFN  = ARRAY OF REFERENCE ATOM LOCATIONS.
C               LOCDEP  = ARRAY OF DEPENDENT ATOM LOCATIONS.
C
C***********************************************************************
C
C     LOCDEP IS THE ATOM WHOSE COORDINATES DEPEND ON THE COORDINATES OF
C     LOCPAR.
C     LOCPAR IS THE ATOM WHOSE COORDINATES ARE USED TO CALCULATE THOSE
C     OF LOCDEP
C     IDEPFN POINTS TO THE PARTICULAR FUNCTION TO BE USED (SEE NDDO)
C
C***********************************************************************
      DIMENSION IVALUE(40),VALUE(40)
      CHARACTER  TEXT(18)*60, KEYWRD*80, LINE*80
         SAVE                                                           GL0892
      DATA TEXT/
     1' BOND LENGTH    IS SET EQUAL TO THE REFERENCE BOND LENGTH   ',
     2' BOND ANGLE     IS SET EQUAL TO THE REFERENCE BOND ANGLE    ',
     3' DIHEDRAL ANGLE IS SET EQUAL TO THE REFERENCE DIHEDRAL ANGLE',
     4' DIHEDRAL ANGLE VARIES AS  90 DEGREES - REFERENCE DIHEDRAL  ',
     5' DIHEDRAL ANGLE VARIES AS  90 DEGREES + REFERENCE DIHEDRAL  ',
     6' DIHEDRAL ANGLE VARIES AS 120 DEGREES - REFERENCE DIHEDRAL  ',
     7' DIHEDRAL ANGLE VARIES AS 120 DEGREES + REFERENCE DIHEDRAL  ',
     8' DIHEDRAL ANGLE VARIES AS 180 DEGREES - REFERENCE DIHEDRAL  ',
     9' DIHEDRAL ANGLE VARIES AS 180 DEGREES + REFERENCE DIHEDRAL  ',
     1' DIHEDRAL ANGLE VARIES AS 240 DEGREES - REFERENCE DIHEDRAL  ',
     2' DIHEDRAL ANGLE VARIES AS 240 DEGREES + REFERENCE DIHEDRAL  ',
     3' DIHEDRAL ANGLE VARIES AS 270 DEGREES - REFERENCE DIHEDRAL  ',
     4' DIHEDRAL ANGLE VARIES AS 270 DEGREES - REFERENCE DIHEDRAL  ',
     5' DIHEDRAL ANGLE VARIES AS - REFERENCE DIHEDRAL              ',
     6' BOND LENGTH VARIES AS HALF THE REFERENCE BOND LENGTH       ',
     7' BOND ANGLE VARIES AS HALF THE REFERENCE BOND ANGLE         ',
     8' BOND ANGLE VARIES AS 180 DEGREES - REFERENCE BOND ANGLE    ',
     9' THE USER HAS TO SUPPLY THIS FUNCTION IN DEPVAR             '/
C
C TITLE OUTPUT
      WRITE (6,10)
   10 FORMAT (///5X,25HPARAMETER DEPENDENCE DATA//
     1'        REFERENCE ATOM      FUNCTION NO.    DEPENDENT ATOM(S)')
C
C INPUT SYMMETRY : FUNCTION, REFERANCE PARAMETER, AND DEPENDENT ATOMS
C
      NDEP=0
   20 READ(IREAD,'(A)',END=70) LINE
      CALL NUCHAR(LINE,VALUE,NVALUE)
C   INTEGER VALUES
      DO 30 I=1,NVALUE
   30 IVALUE(I)=VALUE(I)
C   FILL THE LOCDEP ARRAY
      IF(NVALUE.EQ.0.OR.IVALUE(3).EQ.0) GO TO 70
      DO 40 I=3,NVALUE
         IF(IVALUE(I).EQ.0) GOTO 50
         NDEP=NDEP+1
         LOCDEP(NDEP)=IVALUE(I)
         LOCPAR(NDEP)=IVALUE(1)
         IDEPFN(NDEP)=IVALUE(2)
   40 CONTINUE
   50 LL=I-1
      WRITE(6,60)IVALUE(1),IVALUE(2),(IVALUE(J),J=3,LL)
   60 FORMAT(I13,I19,I14,20I3)
      GO TO 20
C
C CLEAN UP
   70 CONTINUE
      WRITE(6,80)
   80 FORMAT(/10X,'   DESCRIPTIONS OF THE FUNCTIONS USED',/)
      DO 120 J=1,18
         DO 90 I=1,NDEP
            IF(IDEPFN(I).EQ.J) GOTO 100
   90    CONTINUE
         GOTO 120
  100    WRITE(6,110)J,TEXT(J)
  110    FORMAT(I4,5X,A)
  120 CONTINUE
      RETURN
      END
      SUBROUTINE GMETRY(GEO,COORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1NA(NUMATM),NB(NUMATM),NC(NUMATM)
     2       /EULER / TVEC(3,3), ID
      COMMON /REACTN/ STEP, GEOA(3,NUMATM), GEOVEC(3,NUMATM),COLCST
      COMMON /GEOOK/ IGEOOK
      COMMON /NUMCAL/ NUMCAL
      DIMENSION GEO(3,*),COORD(3,*)
      CHARACTER *15 NDIMEN(4)
      LOGICAL GEOOK
         SAVE                                                           GL0892
      DATA NDIMEN/' MOLECULE     ',' POLYMER       ',
     1'LAYER STRUCTURE',' SOLID         '/
      DATA ICALCN /0/
C***********************************************************************
C
C    GMETRY  COMPUTES COORDINATES FROM BOND-ANGLES AND LENGTHS.
C *** IT IS ADAPTED FROM THE PROGRAM WRITTEN BY M.J.S. DEWAR.
C
C     THREE SEPARATE OPTIONS EXIST WITHIN GMETRY. THESE ARE:
C    (A) IF NA(1) IS EQUAL TO 99 (IMPOSSIBLE UNDER NORMAL CIRCUMSTANCES)
C        THEN GEO IS ASSUMED TO BE IN CARTESIAN RATHER THAN INTERNAL
C        COORDINATES, AND COORD IS THEN SET EQUAL TO GEO.
C    (B) IF STEP IS NON-ZERO (THIS IS THE CASE WHEN "SADDLE" IS USED)
C        THEN GEO IS FIRST MODIFIED BY SHIFTING THE INTERNAL COORDINATES
C        ALONG A RADIUS FROM GEOA TO PLACE GEO AT ADISTANCESTEPFROMGEOA.
C    (C) NORMAL CONVERSION FROM INTERNAL TO CARTESIAN COORDINATESISDONE.
C
C  ON INPUT:
C         GEO    = ARRAY OF INTERNAL COORDINATES.
C         NATOMS = NUMBER OF ATOMS, INCLUDING DUMMIES.
C         NA     = ARRAY OF ATOM LABELS FOR BOND LENGTHS.
C
C  ON OUTPUT:
C         COORD  = ARRAY OF CARTESIAN COORDINATES
C
C***********************************************************************
C                                     OPTION (B)
      GEOOK=(IGEOOK.EQ.99)
      IF(ABS(STEP) .GT. 1.D-4) THEN
         SUM=0.D0
         DO 10 I=1,NATOMS
            DO 10 J=1,3
               GEOVEC(J,I)=GEO(J,I)-GEOA(J,I)
   10    SUM=SUM+GEOVEC(J,I)**2
         SUM=SQRT(SUM)
         ERROR=(SUM-STEP)/SUM
      ELSE
         ERROR=0.D0
      ENDIF
      DO 20 I=1,NATOMS
         DO 20 J=1,3
   20 GEO(J,I)=GEO(J,I)-ERROR*GEOVEC(J,I)
C                                     OPTION (A)
      IF(NA(1).EQ.99) THEN
         DO 30 I=1,3
            DO 30 J=1,NATOMS
   30    COORD(I,J)=GEO(I,J)
         GOTO 110
      ENDIF
C                                     OPTION (C)
      COORD(1,1)=0.0D00
      COORD(2,1)=0.0D00
      COORD(3,1)=0.0D00
      COORD(1,2)=GEO(1,2)
      COORD(2,2)=0.0D00
      COORD(3,2)=0.0D00
      IF(NATOMS.EQ.2) GOTO 110
      CCOS=COS(GEO(2,3))
      IF(NA(3).EQ.1)THEN
         COORD(1,3)=COORD(1,1)+GEO(1,3)*CCOS
      ELSE
         COORD(1,3)=COORD(1,2)-GEO(1,3)*CCOS
      ENDIF
      COORD(2,3)=GEO(1,3)*SIN(GEO(2,3))
      COORD(3,3)=0.0D00
      DO 100 I=4,NATOMS
         COSA=COS(GEO(2,I))
         MB=NB(I)
         MC=NA(I)
         XB=COORD(1,MB)-COORD(1,MC)
         YB=COORD(2,MB)-COORD(2,MC)
         ZB=COORD(3,MB)-COORD(3,MC)
         RBC=1.0D00/SQRT(XB*XB+YB*YB+ZB*ZB)
         IF (ABS(COSA).LT.1.D0-1.D-12) GO TO 40
C
C     ATOMS MC, MB, AND (I) ARE COLLINEAR
C
         RBC=GEO(1,I)*RBC*COSA
         COORD(1,I)=COORD(1,MC)+XB*RBC
         COORD(2,I)=COORD(2,MC)+YB*RBC
         COORD(3,I)=COORD(3,MC)+ZB*RBC
         GO TO 100
C
C     THE ATOMS ARE NOT COLLINEAR
C
   40    MA=NC(I)
         XA=COORD(1,MA)-COORD(1,MC)
         YA=COORD(2,MA)-COORD(2,MC)
         ZA=COORD(3,MA)-COORD(3,MC)
C
C     ROTATE ABOUT THE Z-AXIS TO MAKE YB=0, AND XB POSITIVE.  IF XYB IS
C     TOO SMALL, FIRST ROTATE THE Y-AXIS BY 90 DEGREES.
C
         XYB=SQRT(XB*XB+YB*YB)
         K=-1
         IF (XYB.GT.0.1D00) GO TO 50
         XPA=ZA
         ZA=-XA
         XA=XPA
         XPB=ZB
         ZB=-XB
         XB=XPB
         XYB=SQRT(XB*XB+YB*YB)
         K=+1
C
C     ROTATE ABOUT THE Y-AXIS TO MAKE ZB VANISH
C
   50    COSTH=XB/XYB
         SINTH=YB/XYB
         XPA=XA*COSTH+YA*SINTH
         YPA=YA*COSTH-XA*SINTH
         SINPH=ZB*RBC
         COSPH=SQRT(ABS(1.D00-SINPH*SINPH))
         XQA=XPA*COSPH+ZA*SINPH
         ZQA=ZA*COSPH-XPA*SINPH
C
C     ROTATE ABOUT THE X-AXIS TO MAKE ZA=0, AND YA POSITIVE.
C
         YZA=SQRT(YPA**2+ZQA**2)
         IF(YZA.LT.1.D-4)GOTO 70
         IF(YZA.LT.2.D-2 .AND. .NOT.GEOOK)THEN
            WRITE(6,'(//20X,'' CALCULATION ABANDONED AT THIS POINT'')')
            WRITE(6,'(//10X,'' THREE ATOMS BEING USED TO DEFINE THE'',/
     110X,'' COORDINATES OF A FOURTH ATOM, WHOSE BOND-ANGLE IS'')')
            WRITE(6,'(10X,'' NOT ZERO OR 180 DEGREEES, ARE '',
     1''IN AN ALMOST STRAIGHT'')')
            WRITE(6,'(10X,'' LINE.  THERE IS A HIGH PROBABILITY THAT THE
     1'',/10X,'' COORDINATES OF THE ATOM WILL BE INCORRECT.'')')
            WRITE(6,'(//20X,''THE FAULTY ATOM IS ATOM NUMBER'',I4)')I
            CALL GEOUT
            WRITE(6,'(//20X,''CARTESIAN COORDINATES UP TO FAULTY ATOM'')
     1')
            WRITE(6,'(//5X,''I'',12X,''X'',12X,''Y'',12X,''Z'')')
            DO 60 J=1,I
   60       WRITE(6,'(I6,F16.5,2F13.5)')J,(COORD(K,J),K=1,3)
            WRITE(6,'(//6X,'' ATOMS'',I3,'','',I3,'', AND'',I3,
     1'' ARE WITHIN'',F7.4,'' ANGSTROMS OF A STRAIGHT LINE'')')
     2MC,MB,MA,YZA
            STOP
         ENDIF
         COSKH=YPA/YZA
         SINKH=ZQA/YZA
         GOTO 80
   70    CONTINUE
C
C   ANGLE TOO SMALL TO BE IMPORTANT
C
         COSKH=1.D0
         SINKH=0.D0
   80    CONTINUE
C
C     COORDINATES :-   A=(XQA,YZA,0),   B=(RBC,0,0),  C=(0,0,0)
C     NONE ARE NEGATIVE.
C     THE COORDINATES OF I ARE EVALUATED IN THE NEW FRAME.
C
         SINA=SIN(GEO(2,I))
         SIND=-SIN(GEO(3,I))
         COSD=COS(GEO(3,I))
         XD=GEO(1,I)*COSA
         YD=GEO(1,I)*SINA*COSD
         ZD=GEO(1,I)*SINA*SIND
C
C     TRANSFORM THE COORDINATES BACK TO THE ORIGINAL SYSTEM.
C
         YPD=YD*COSKH-ZD*SINKH
         ZPD=ZD*COSKH+YD*SINKH
         XPD=XD*COSPH-ZPD*SINPH
         ZQD=ZPD*COSPH+XD*SINPH
         XQD=XPD*COSTH-YPD*SINTH
         YQD=YPD*COSTH+XPD*SINTH
         IF (K.LT.1) GO TO 90
         XRD=-ZQD
         ZQD=XQD
         XQD=XRD
   90    COORD(1,I)=XQD+COORD(1,MC)
         COORD(2,I)=YQD+COORD(2,MC)
         COORD(3,I)=ZQD+COORD(3,MC)
  100 CONTINUE
C
C *** NOW REMOVE THE TRANSLATION VECTORS, IF ANY, FROM THE ARRAY COOR
C
  110 CONTINUE
      K=NATOMS
  120 IF(LABELS(K).NE.107) GOTO 130
      K=K-1
      GOTO 120
  130 K=K+1
      IF(K.GT.NATOMS) GOTO 180
C
C   SYSTEM IS A SOLID, OF DIMENSION NATOMS+1-K
C
      L=0
      DO 140 I=K,NATOMS
         L=L+1
         MC=NA(I)
         DO 140 J=1,3
            TVEC(J,L)=COORD(J,I)-COORD(J,MC)
  140 CONTINUE
      ID=L
      IF (ICALCN .NE. NUMCAL) THEN
         ICALCN = NUMCAL
         WRITE(6,150)NDIMEN(ID+1)
  150    FORMAT(/10X,'    THE SYSTEM IS A ',A15,/)
         IF(ID.EQ.0) GOTO 180
         WRITE(6,160)
         WRITE(6,170)(I,(TVEC(J,I),J=1,3),I=1,ID)
  160    FORMAT(/,'                UNIT CELL TRANSLATION VECTORS',/
     1/,'              X              Y              Z')
  170    FORMAT('    T',I1,' = ',F11.7,'    ',F11.7,'    ',F11.7)
      ENDIF
  180 CONTINUE
      J=0
      DO 200 I=1,NATOMS
         IF (LABELS(I).EQ.99.OR.LABELS(I).EQ.107) GO TO 200
         J=J+1
         DO 190 K=1,3
  190    COORD(K,J)=COORD(K,I)
  200 CONTINUE
      CUTOFF=200.D0
      RETURN
      END
      SUBROUTINE H1ELEC(NI,NJ,XI,XJ,SMAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XI(3),XJ(3),SMAT(9,9), BI(9), BJ(9)
C***********************************************************************
C
C  H1ELEC FORMS THE ONE-ELECTRON MATRIX BETWEEN TWO ATOMS.
C
C   ON INPUT    NI   = ATOMIC NO. OF FIRST ATOM.
C               NJ   = ATOMIC NO. OF SECOND ATOM.
C               XI   = COORDINATES OF FIRST ATOM.
C               XJ   = COORDINATES OF SECOND ATOM.
C
C   ON OUTPUT   SMAT = MATRIX OF ONE-ELECTRON INTERACTIONS.
C
C***********************************************************************
      COMMON /BETAS / BETAS(107),BETAP(107),BETAD(107)
      COMMON /BETA3 / BETA3(153)
      COMMON /KEYWRD/ KEYWRD
      COMMON /EULER / TVEC(3,3), ID
      COMMON /VSIPS / VS(107),VP(107),VD(107)
      COMMON /NATORB/ NATORB(107)
      COMMON /UCELL / L1L,L2L,L3L,L1U,L2U,L3U
      COMMON /NUMCAL/ NUMCAL
      COMMON /SRPL/ ISSRP                                               IR0494
      DIMENSION SBITS(9,9), LIMS(3,2), XJUC(3)
      CHARACTER*80 KEYWRD
      LOGICAL ISSRP                                                     IR0494
      EQUIVALENCE (L1L,LIMS(1,1))
         SAVE                                                           GL0892
      DATA ITYPE /1/
      DATA ICALCN /0/
      IF(NI.EQ.102.OR.NJ.EQ.102)THEN
         IF(SQRT((XI(1)-XJ(1))**2+
     1        (XI(2)-XJ(2))**2+
     2        (XI(3)-XJ(3))**2) .GT.1.8)THEN
            DO 10 I=1,9
               DO 10 J=1,9
   10       SMAT(I,J)=0.D0
            RETURN
         ENDIF
      ENDIF
      IF(ID.EQ.0) THEN
         CALL DIAT(NI,NJ,XI,XJ,SMAT)
      ELSE
      IF (ICALCN .NE. NUMCAL) THEN
            ICALCN = NUMCAL
            DO 20 I=1,ID
               LIMS(I,1)=-1
   20       LIMS(I,2)= 1
            DO 30 I=ID+1,3
               LIMS(I,1)=0
   30       LIMS(I,2)=0
         ENDIF
         DO 40 I=1,9
            DO 40 J=1,9
   40    SMAT(I,J)=0
         DO 70 I=L1L,L1U
            DO 70 J=L2L,L2U
               DO 70 K=L3L,L3U
                  DO 50 L=1,3
   50             XJUC(L)=XJ(L)+TVEC(L,1)*I+TVEC(L,2)*J+TVEC(L,3)*K
                  CALL DIAT(NI,NJ,XI,XJUC,SBITS)
                  DO 60 L=1,9
                     DO 60 M=1,9
   60             SMAT(L,M)=SMAT(L,M)+SBITS(L,M)
   70    CONTINUE
      ENDIF
   80 GOTO (90,100,110) ITYPE
   90 IF(INDEX(KEYWRD,'MINDO') .NE. 0) THEN
         ITYPE=2
      ELSE
         ITYPE=3
      ENDIF
      GOTO 80
  100 CONTINUE
      II=MAX(NI,NJ)
      NBOND=(II*(II-1))/2+NI+NJ-II
      IF(NBOND.GT.153)GOTO 120
      BI(1)=BETA3(NBOND)*VS(NI)
      BI(2)=BETA3(NBOND)*VP(NI)
      BI(3)=BI(2)
      BI(4)=BI(2)
      BJ(1)=BETA3(NBOND)*VS(NJ)
      BJ(2)=BETA3(NBOND)*VP(NJ)
      BJ(3)=BJ(2)
      BJ(4)=BJ(2)
      GOTO 120
  110 CONTINUE
      BI(1)=BETAS(NI)*0.5D0
      BI(2)=BETAP(NI)*0.5D0
      BI(3)=BI(2)
      BI(4)=BI(2)
      BI(5)=BETAD(NI)*0.5D0
      BI(6)=BI(5)
      BI(7)=BI(5)
      BI(8)=BI(5)
      BI(9)=BI(5)
      BJ(1)=BETAS(NJ)*0.5D0
      BJ(2)=BETAP(NJ)*0.5D0
      BJ(3)=BJ(2)
      BJ(4)=BJ(2)
      BJ(5)=BETAD(NJ)*0.5D0
      BJ(6)=BJ(5)
      BJ(7)=BJ(5)
      BJ(8)=BJ(5)
      BJ(9)=BJ(5)
  120 CONTINUE
      NORBI=NATORB(NI)
      NORBJ=NATORB(NJ)
      IF(NORBI.EQ.9.OR.NORBJ.EQ.9) THEN
C
C    IN THE CALCULATION OF THE ONE-ELECTRON TERMS THE GEOMETRIC MEAN
C    OF THE TWO BETA VALUES IS BEING USED IF ONE OF THE ATOMS
C    CONTAINS D-ORBITALS.
         DO 130 J=1,NORBJ
            DO 130 I=1,NORBI
  130    SMAT(I,J)=-2.0D0*SMAT(I,J)*SQRT(BI(I)*BJ(J))
      ELSE IF (ISSRP) THEN                                              IR0494
C
C ...USE SPECIAL SRP BETAs   (Ivan Rossi - Aprile '94)
C
         IF(NORBI .gt. 1 .and. norbj .gt. 1) THEN
           BPP=GETBET( NI, NJ, 'BETPP')
           DO 132 JBET=2,4
             DO 132 IBET=2,4
  132          SMAT(IBET,JBET)=SMAT(IBET,JBET)*BPP
         ENDIF
         IF (NORBI .gt. 1 ) THEN
            BSP=GETBET( NJ, NI, 'BETSP')
            DO 134 IBET=2,4
  134          SMAT(IBET,1)=SMAT(IBET,1)*BSP
         ENDIF
         IF (NORBJ .gt. 1 ) THEN
            BSP=GETBET( NI, NJ, 'BETSP')
            DO 136 JBET=2,4
  136          SMAT(1,JBET)=SMAT(1,JBET)*BSP
         ENDIF
         SMAT(1,1)=SMAT(1,1)*GETBET( NJ, NI, 'BETSS')
C end of SRP beta
      ELSE
         DO 140 J=1,NORBJ
            DO 140 I=1,NORBI
  140    SMAT(I,J)=SMAT(I,J)*(BI(I)+BJ(J))
      ENDIF
      RETURN
      END
      SUBROUTINE HCORE (COORD,H,W, WJ,WK,ENUCLR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      LOGICAL FLDON
      DIMENSION COORD(3,*),H(*), WJ(*), WK(*), W(*)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /MOLORB/ USPD(MAXORB),DUMY(MAXORB)
      COMMON /KEYWRD/ KEYWRD
      COMMON /EULER / TVEC(3,3), ID
      COMMON /MULTIP/ DD(107),QQ(107),AM(107),AD(107),AQ(107)
      COMMON /CORE  / CORE(107)
      COMMON /FIELD / EFIELD(3)
      COMMON /NUMCAL/ NUMCAL
************************************************************************
C
C   HCORE GENERATES THE ONE-ELECTRON MATRIX AND TWO ELECTRON INTEGRALS
C         FOR A GIVEN MOLECULE WHOSE GEOMETRY IS DEFINED IN CARTESIAN
C         COORDINATES.
C
C  ON INPUT  COORD   = COORDINATES OF THE MOLECULE.
C
C  ON OUTPUT  H      = ONE-ELECTRON MATRIX.
C             W      = TWO-ELECTRON INTEGRALS.
C             ENUCLR = NUCLEAR ENERGY
************************************************************************
      CHARACTER*80 KEYWRD
      LOGICAL DEBUG
      DIMENSION E1B(10),E2A(10),DI(9,9), WJD(100), WKD(100)
         SAVE                                                           GL0892
      DATA ICALCN /0/
      IF (ICALCN .NE. NUMCAL) THEN
         IONE=1
         CUTOFF=1.D10
         IF(ID.NE.0)CUTOFF=60.D0
         IF(ID.NE.0)IONE=0
         ICALCN = NUMCAL
         DEBUG=(INDEX(KEYWRD,'HCORE') .NE. 0)
      ENDIF
      FLDON = .FALSE.
      IF ((EFIELD(1).NE.0.0D00).OR.(EFIELD(2).NE.0.0D00).OR.
     1    (EFIELD(3).NE.0.0D00)) THEN
         FLDCON = 51.4257D00
         FLDON = .TRUE.
      ENDIF
      DO 10 I=1,(NORBS*(NORBS+1))/2
   10 H(I)=0
      ENUCLR=0.D0
      KR=1
      DO 110 I=1,NUMAT
         IA=NFIRST(I)
         IB=NLAST(I)
         IC=NMIDLE(I)
         NI=NAT(I)
C
C FIRST WE FILL THE DIAGONALS, AND OFF-DIAGONALS ON THE SAME ATOM
C
         DO 30 I1=IA,IB
            I2=I1*(I1-1)/2+IA-1
            DO 20 J1=IA,I1
               I2=I2+1
               H(I2)=0.D0
               IF (FLDON) THEN
                  IO1 = I1 - IA
                  JO1 = J1 - IA
                  IF ((JO1.EQ.0).AND.(IO1.EQ.1)) THEN
                     HTERME = -0.529177D00*DD(NI)*EFIELD(1)*FLDCON
                     H(I2) = HTERME
                  ENDIF
                  IF ((JO1.EQ.0).AND.(IO1.EQ.2)) THEN
                     HTERME = -0.529177D00*DD(NI)*EFIELD(2)*FLDCON
                     H(I2) = HTERME
                  ENDIF
                  IF ((JO1.EQ.0).AND.(IO1.EQ.3)) THEN
                     HTERME = -0.529177D00*DD(NI)*EFIELD(3)*FLDCON
                     H(I2) = HTERME
                  ENDIF
               ENDIF
   20       CONTINUE
            H(I2) = USPD(I1)
            IF (FLDON) THEN
               FNUC = -(EFIELD(1)*COORD(1,I) + EFIELD(2)*COORD(2,I) +
     1              EFIELD(3)*COORD(3,I))*FLDCON
               H(I2) = H(I2) + FNUC
            ENDIF
   30    CONTINUE
C
C   FILL THE ATOM-OTHER ATOM ONE-ELECTRON MATRIX<PSI(LAMBDA)|PSI(SIGMA)>
C
         IM1=I-IONE
         DO 100 J=1,IM1
            HALF=1.D0
            IF(I.EQ.J)HALF=0.5D0
            JA=NFIRST(J)
            JB=NLAST(J)
            JC=NMIDLE(J)
            NJ=NAT(J)
            CALL H1ELEC(NI,NJ,COORD(1,I),COORD(1,J),DI)
            I2=0
            DO 40 I1=IA,IB
               II=I1*(I1-1)/2+JA-1
               I2=I2+1
               J2=0
               JJ=MIN(I1,JB)
               DO 40 J1=JA,JJ
                  II=II+1
                  J2=J2+1
   40       H(II)=H(II)+DI(I2,J2)
C
C   CALCULATE THE TWO-ELECTRON INTEGRALS, W; THE ELECTRON NUCLEAR TERMS
C   E1B AND E2A; AND THE NUCLEAR-NUCLEAR TERM ENUC.
C
            IF(ID.EQ.0) THEN
               CALL ROTATE(NI,NJ,COORD(1,I),COORD(1,J),
     1 W(KR), KR,E1B,E2A,ENUC,CUTOFF)
            ELSE
               KRO=KR
               CALL SOLROT(NI,NJ,COORD(1,I),COORD(1,J),
     1                WJD, WKD,KR,E1B,E2A,ENUC,CUTOFF)
               JJ=0
               DO 50 II=KRO,KR-1
                  JJ=JJ+1
                  WJ(II)=WJD(JJ)
   50          WK(II)=WKD(JJ)
            ENDIF
            ENUCLR = ENUCLR + ENUC
C
C   ADD ON THE ELECTRON-NUCLEAR ATTRACTION TERM FOR ATOM I.
C
            I2=0
            DO 60 I1=IA,IC
               II=I1*(I1-1)/2+IA-1
               DO 60 J1=IA,I1
                  II=II+1
                  I2=I2+1
   60       H(II)=H(II)+E1B(I2)*HALF
            DO  70 I1=IC+1,IB
               II=(I1*(I1+1))/2
   70       H(II)=H(II)+E1B(1)*HALF
C
C   ADD ON THE ELECTRON-NUCLEAR ATTRACTION TERM FOR ATOM J.
C
            I2=0
            DO 80 I1=JA,JC
               II=I1*(I1-1)/2+JA-1
               DO 80 J1=JA,I1
                  II=II+1
                  I2=I2+1
   80       H(II)=H(II)+E2A(I2)*HALF
            DO 90 I1=JC+1,JB
               II=(I1*(I1+1))/2
   90       H(II)=H(II)+E2A(1)*HALF
  100    CONTINUE
  110 CONTINUE
      IF( .NOT. DEBUG) RETURN
      WRITE(6,'(//10X,''ONE-ELECTRON MATRIX FROM HCORE'')')
      CALL VECPRT(H,NORBS)
      J=MIN(400,KR)
      IF(ID.EQ.0) THEN
         WRITE(6,'(//10X,''TWO-ELECTRON MATRIX IN HCORE''/)')
         WRITE(6,120)(W(I),I=1,J)
      ELSE
         WRITE(6,'(//10X,''TWO-ELECTRON J MATRIX IN HCORE''/)')
         WRITE(6,120)(WJ(I),I=1,J)
         WRITE(6,'(//10X,''TWO-ELECTRON K MATRIX IN HCORE''/)')
         WRITE(6,120)(WK(I),I=1,J)
      ENDIF
  120 FORMAT(10F8.4)
      RETURN
      END
      SUBROUTINE INTERP(N,NP,NQ,MODE,E,FP,CP,VEC,FOCK,P,H,VECL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION FP(MPACK), CP(N,N)
      DIMENSION VEC(N,N), FOCK(N,N),
     1          P(N,N), H(N*N), VECL(N*N)
**********************************************************************
*
* INTERP: AN INTERPOLATION PROCEDURE FOR FORCING SCF CONVERGANCE
*         ORIGINAL THEORY AND FORTRAN WRITTEN BY R.N. CAMP AND
*         H.F. KING, J. CHEM. PHYS. 75, 268 (1981)
**********************************************************************
*
* ON INPUT N     = NUMBER OF ORBITALS
*          NP    = NUMBER OF FILLED LEVELS
*          NQ    = NUMBER OF EMPTY LEVELS
*          MODE  = 1, DO NOT RESET.
*          E     = ENERGY
*          FP    = FOCK MATRIX, AS LOWER HALF TRIANGLE, PACKED
*          CP    = EIGENVECTORS OF FOCK MATRIX OF ITERATION -1
*                  AS PACKED ARRAY OF N*N COEFFICIENTS
*
* ON OUTPUT CP   = BEST GUESSED SET OF EIGENVECTORS
*           MODE = 2 OR 3 - USED BY CALLING PROGRAM
**********************************************************************
      DIMENSION THETA(MAXORB), IA(MAXORB)
      COMMON /KEYWRD/ KEYWRD
C common FIT splitted for portability (Ivan Rossi 0394 &8) )
      COMMON /FITPTS/ NPNTS
     1       /FITDAT/ XLOW,XHIGH,XMIN,EMIN,DEMIN,X(12),F(12),DF(12)
      COMMON /NUMCAL/ NUMCAL
      LOGICAL DEBUG, DEBUG1
      CHARACTER*80 KEYWRD
         SAVE                                                           GL0892
      DATA ICALCN /0/
      DATA ZERO,FF,RADMAX/0.0,0.9,1.5708/
      IF (ICALCN .NE. NUMCAL) THEN
         ICALCN = NUMCAL
         DEBUG=(INDEX(KEYWRD,'INTERP').NE.0)
         DEBUG1=(INDEX(KEYWRD,'DEBUG').NE.0.AND.DEBUG)
         DO 10 I=1,MAXORB
   10    IA(I)=(I*I-I)/2
      ENDIF
C
C     RADMAX=MAXIMUM ROTATION ANGLE (RADIANS).  1.5708 = 90 DEGREES.
C         FF=FACTOR FOR CONVERGENCE TEST FOR 1D SEARCH.
C
      MINPQ=MIN0(NP,NQ)
      MAXPQ=MAX0(NP,NQ)
      NP1=NP+1
      NP2=MAX0(1,NP/2)
      IF(MODE.EQ.2) GO TO 110
C
C     (MODE=1 OR 3 ENTRY)
C     TRANSFORM FOCK MATRIX TO CURRENT MO BASIS.
C     ONLY THE OFF DIAGONAL OCC-VIRT BLOCK IS COMPUTED.
C     STORE IN FOCK ARRAY
C
      II=0
      DO 50 I=1,N
         I1=I+1
         DO 40 J=1,NQ
            DUM=ZERO
            DO 20 K=1,I
   20       DUM=DUM+FP(II+K)*CP(K,J+NP)
            IF(I.EQ.N) GO TO 40
            IK=II+I+I
            DO 30 K=I1,N
               DUM=DUM+FP(IK)*CP(K,J+NP)
   30       IK=IK+K
   40    P(I,J)=DUM
   50 II=II+I
      DO 80 I=1,NP
         DO 70 J=1,NQ
            DUM=ZERO
            DO 60 K=1,N
   60       DUM=DUM+CP(K,I)*P(K,J)
            FOCK(I,J)=DUM
   70    CONTINUE
   80 CONTINUE
      IF(MODE.EQ.3) GO TO 100
C
C     CURRENT POINT BECOMES OLD POINT (MODE=1 ENTRY)
C
      DO 90 I=1,N
         DO 90 J=1,N
   90 VEC(I,J)=CP(I,J)
      EOLD=E
      XOLD=1.0
      MODE=2
      RETURN
C
C     (MODE=3 ENTRY)
C     FOCK CORRESPONDS TO CURRENT POINT IN CORRESPONDING REPRESENTATION.
C     VEC DOES NOT HOLD CURRENT VECTORS. VEC SET IN LAST MODE=2 ENTRY.
C
  100 NPNTS=NPNTS+1
      IF(DEBUG)WRITE(6,'(''   INTERPOLATED ENERGY:'',F13.6)')E*23.061D0
      IPOINT=NPNTS
      GO TO 500
C
C    (MODE=2 ENTRY) CALCULATE THETA, AND U, V, W MATRICES.
C                   U ROTATES CURRENT INTO OLD MO.
C                   V ROTATES CURRENT INTO CORRESPONDING CURRENT MO.
C                   W ROTATES OLD INTO CORRESPONDING OLD MO.
C
  110 J1=1
      DO 140 I=1,N
         IF(I.EQ.NP1) J1=NP1
         DO 130 J=J1,N
            P(I,J)=ZERO
            DO 120 K=1,N
  120       P(I,J)=P(I,J)+CP(K,I)*VEC(K,J)
  130    CONTINUE
  140 CONTINUE
C
C     U = CP(DAGGER)*VEC IS NOW IN P ARRAY.
C     VEC IS NOW AVAILABLE FOR TEMPORARY STORAGE.
C
      IJ=0
      DO 170 I=1,NP
         DO 160 J=1,I
            IJ=IJ+1
            H(IJ)=0.D0
            DO 150 K=NP1,N
  150       H(IJ)=H(IJ)+P(I,K)*P(J,K)
  160    CONTINUE
  170 CONTINUE
      CALL HQRII(H,NP,NP,THETA,VECL)
      DO 180 I=NP,1,-1
         IL=I*NP-NP
         DO 180 J=NP,1,-1
  180 VEC(J,I)=VECL(J+IL)
      DO 200 I=1,NP2
         DUM=THETA(NP1-I)
         THETA(NP1-I)=THETA(I)
         THETA(I)=DUM
         DO 190 J=1,NP
            DUM=VEC(J,NP1-I)
            VEC(J,NP1-I)=VEC(J,I)
  190    VEC(J,I)=DUM
  200 CONTINUE
      DO 210 I=1,MINPQ
         THETA(I)=MAX(THETA(I),ZERO)
         THETA(I)=MIN(THETA(I),1.D0)
  210 THETA(I)=ASIN(SQRT(THETA(I)))
C
C     THETA MATRIX HAS NOW BEEN CALCULATED, ALSO UNITARY VP MATRIX
C     HAS BEEN CALCULATED AND STORED IN FIRST NP COLUMNS OF VEC MATRIX.
C     NOW COMPUTE WQ
C
      DO 240 I=1,NQ
         DO 230 J=1,MINPQ
            VEC(I,NP+J)=ZERO
            DO 220 K=1,NP
  220       VEC(I,NP+J)=VEC(I,NP+J)+P(K,NP+I)*VEC(K,J)
  230    CONTINUE
  240 CONTINUE
      CALL SCHMIT(VEC(1,NP1),NQ,N)
C
C     UNITARY WQ MATRIX NOW IN LAST NQ COLUMNS OF VEC MATRIX.
C     TRANSPOSE NP BY NP BLOCK OF U STORED IN P
C
      DO 260 I=1,NP
         DO 250 J=1,I
            DUM=P(I,J)
            P(I,J)=P(J,I)
  250    P(J,I)=DUM
  260 CONTINUE
C
C     CALCULATE WP MATRIX AND HOLD IN FIRST NP COLUMNS OF P
C
      DO 300 I=1,NP
         DO 270 K=1,NP
  270    H(K)=P(I,K)
         DO 290 J=1,NP
            P(I,J)=ZERO
            DO 280 K=1,NP
  280       P(I,J)=P(I,J)+H(K)*VEC(K,J)
  290    CONTINUE
  300 CONTINUE
      CALL SCHMIB(P,NP,N)
C
C     CALCULATE VQ MATRIX AND HOLD IN LAST NQ COLUMNS OF P MATRIX.
C
      DO 340 I=1,NQ
         DO 310 K=1,NQ
  310    H(K)=P(NP+I,NP+K)
         DO 330 J=NP1,N
            P(I,J)=ZERO
            DO 320 K=1,NQ
  320       P(I,J)=P(I,J)+H(K)*VEC(K,J)
  330    CONTINUE
  340 CONTINUE
      CALL SCHMIB(P(1,NP1),NQ,N)
C
C     CALCULATE (DE/DX) AT OLD POINT
C
      DEDX=ZERO
      DO 370 I=1,NP
         DO 360 J=1,NQ
            DUM=ZERO
            DO 350 K=1,MINPQ
  350       DUM=DUM+THETA(K)*P(I,K)*VEC(J,NP+K)
  360    DEDX=DEDX+DUM*FOCK(I,J)
  370 CONTINUE
C
C     STORE OLD POINT INFORMATION FOR SPLINE FIT
C
      DEOLD=-4.0*DEDX
      X(2)=XOLD
      F(2)=EOLD
      DF(2)=DEOLD
C
C     MOVE VP OUT OF VEC ARRAY INTO FIRST NP COLUMNS OF P MATRIX.
C
      DO 380 I=1,NP
         DO 380 J=1,NP
  380 P(I,J)=VEC(I,J)
      K1=0
      K2=NP
      DO 410 J=1,N
         IF(J.EQ.NP1) K1=NP
         IF(J.EQ.NP1) K2=NQ
         DO 400 I=1,N
            DUM=ZERO
            DO 390 K=1,K2
  390       DUM=DUM+CP(I,K1+K)*P(K,J)
  400    VEC(I,J)=DUM
  410 CONTINUE
C
C     CORRESPONDING CURRENT MO VECTORS NOW HELD IN VEC.
C     COMPUTE VEC(DAGGER)*FP*VEC
C     STORE OFF-DIAGONAL BLOCK IN FOCK ARRAY.
C
  420 II=0
      DO 460 I=1,N
         I1=I+1
         DO 450 J=1,NQ
            DUM=ZERO
            DO 430 K=1,I
  430       DUM=DUM+FP(II+K)*VEC(K,J+NP)
            IF(I.EQ.N) GO TO 450
            IK=II+I+I
            DO 440 K=I1,N
               DUM=DUM+FP(IK)*VEC(K,J+NP)
  440       IK=IK+K
  450    P(I,J)=DUM
  460 II=II+I
      DO 490 I=1,NP
         DO 480 J=1,NQ
            DUM=ZERO
            DO 470 K=1,N
  470       DUM=DUM+VEC(K,I)*P(K,J)
            FOCK(I,J)=DUM
  480    CONTINUE
  490 CONTINUE
C
C     SET LIMITS ON RANGE OF 1-D SEARCH
C
      NPNTS=2
      IPOINT=1
      XNOW=ZERO
      XHIGH=RADMAX/THETA(1)
      XLOW=-0.5*XHIGH
C
C     CALCULATE (DE/DX) AT CURRENT POINT AND
C     STORE INFORMATION FOR SPLINE FIT
C     ***** JUMP POINT FOR MODE=3 ENTRY *****
C
  500 DEDX=ZERO
      DO 510 K=1,MINPQ
  510 DEDX=DEDX+THETA(K)*FOCK(K,K)
      DENOW=-4.0*DEDX
      ENOW=E
      IF(IPOINT.LE.12) GO TO 530
  520 FORMAT(//34H EXCESSIVE DATA PNTS FOR SPLINE./
     1,9H IPOINT =,I3,15H MAXIMUM IS 12.)
C
C     PERFORM 1-D SEARCH AND DETERMINE EXIT MODE.
C
  530 X(IPOINT)=XNOW
      F(IPOINT)=ENOW
      DF(IPOINT)=DENOW
      CALL SPLINE
      IF((EOLD-ENOW).GT.FF*(EOLD-EMIN).OR.IPOINT.GT.10) GO TO 560
C
C     (MODE=3 EXIT) RECOMPUTE CP VECTORS AT PREDICTED MINIMUM.
C
      XNOW=XMIN
      DO 550 K=1,MINPQ
         CK=COS(XNOW*THETA(K))
         SK=SIN(XNOW*THETA(K))
         IF(DEBUG)WRITE(6,'('' ROTATION ANGLE:'',F12.4)')SK*57.29578
         DO 540 I=1,N
            CP(I,K)   =CK*VEC(I,K)-SK*VEC(I,NP+K)
  540    CP(I,NP+K)=SK*VEC(I,K)+CK*VEC(I,NP+K)
  550 CONTINUE
      MODE=3
      RETURN
C
C     (MODE=2 EXIT) CURRENT VECTORS GIVE SATISFACTORY ENERGY IMPROVEMENT
C     CURRENT POINT BECOMES OLD POINT FOR THE NEXT 1-D SEARCH.
C
  560 IF(MODE.EQ.2) GO TO 580
      DO 570 I=1,N
         DO 570 J=1,N
  570 VEC(I,J)=CP(I,J)
      MODE=2
  580 ROLD=XOLD*THETA(1)*57.29578
      RNOW=XNOW*THETA(1)*57.29578
      RMIN=XMIN*THETA(1)*57.29578
      IF(DEBUG)WRITE(6,600) XOLD,EOLD*23.061,DEOLD,ROLD
     1,             XNOW,ENOW*23.061,DENOW,RNOW
     2,             XMIN,EMIN*23.061,DEMIN,RMIN
      EOLD=ENOW
      IF(NPNTS.LE.200) RETURN
      WRITE(6,610)
      DO 590 K=1,NPNTS
  590 WRITE(6,620) K,X(K),F(K),DF(K)
      WRITE(6,630)
      RETURN
  600 FORMAT(
     1/14X,3H X ,10X,6H F(X) ,9X,7H DF/DX ,21H   ROTATION (DEGREES),
     2/10H      OLD ,F10.5,3F15.10,
     3/10H  CURRENT ,F10.5,3F15.10,
     4/10H PREDICTED,F10.5,3F15.10/)
  610 FORMAT(3H  K,10H     X(K) ,15H       F(K)    ,10H     DF(K))
  620 FORMAT(I3,F10.5,2F15.10)
  630 FORMAT(10X)
      END
      FUNCTION MECI(EIGS,COEFF,COEFFS,EIGA,N,NMOS,IRESET,FINISH)
***********************************************************************
*
*                 PROGRAM MECI
*
*   A MULTI-ELECTRON CONFIGURATION INTERACTION CALCULATION
*
*   WRITTEN BY JAMES J. P. STEWART, AT THE
*              FRANK J. SEILER RESEARCH LABORATORY
*              USAFA, COLORADO SPRINGS, CO 80840
*
*              1985
*
***********************************************************************
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DOUBLE PRECISION MECI
C
C   MATRICES FOR SEC.DET., VECTORS, AND EIGENVALUES.
C
      DIMENSION CIMAT(NMECI**4), CONF(NMECI**4),
     1 EIG(NMECI**2), DIAG(2*NMECI**3)
C
C   MATRICES TO HOLD ELECTRON CONFIGURATIONS
C
      DIMENSION MICROA(NMECI,2*NMECI**3), MICROB(NMECI,2*NMECI**3),
     1IOCCA1(NMECI), IOCCA2(NMECI), IOCCB1(NMECI),
     2IOCCB2(NMECI), NALPHA(NMECI**2)
C
C   MATRICES FOR PERMUTATION WORK
C
      DIMENSION NFA(2*NMECI), NPERMA(NMECI,6*NMECI),
     1NPERMB(NMECI,6*NMECI)
C
C   MATRICES FOR ONE AND TWO ELECTRON INTEGRALS, AND M.O.S
C
      DIMENSION EIGA(MAXORB), EIGS(MAXORB), RJKAA(NMECI,NMECI),
     1          RJKAB(NMECI,NMECI), COEFF(N,N),  COEFFS(N,N)
C
C   SPIN MATRICES
C
      DIMENSION SPIN(NMECI**2)
      LOGICAL DEBUG,  LARGE, PRNT, LSPIN, LSPIN1, FINISH,
     1 BIGPRT, SING, DOUB, TRIP, QUAR, QUIN, SEXT, LAST1,
     2 PRNT2, FORCE, GEOOK
      CHARACTER KEYWRD*80, TSPIN(7)*8, LINE*80
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NDUMMY(2),NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /SPQR/ ISPQR(NMECI**2,NMECI),IS,I,K
      COMMON /DENSTY/ P(MPACK), PA(MPACK), PB(MPACK)
      COMMON /KEYWRD/ KEYWRD
      COMMON /BASEOC/ OCCA(NMECI)
      COMMON /XYIJKL/ XY(NMECI,NMECI,NMECI,NMECI)
      COMMON /NUMCAL/ NUMCAL
      COMMON /IOCM/ IREAD
         SAVE                                                           GL0892
      DATA TSPIN/'SINGLET ','DOUBLET ','TRIPLET ','QUARTET ','QUINTET ',
     1'SEXTET  ','SEPTET  '/
      DATA ICALCN /0/, ICALC1 /0/
         IF (ICALCN .NE. NUMCAL) THEN
            ICALCN = NUMCAL
         GEOOK=(INDEX(KEYWRD,'GEO-OK').NE.0)
         FORCE=(INDEX(KEYWRD,'FORCE').NE.0)
         LAST1=(INDEX(KEYWRD,'1SCF').NE.0)
         LSPIN1=(INDEX(KEYWRD,'ESR').NE.0)
         MDIM=NMECI**2
         DEBUG=(INDEX(KEYWRD,'DEBUG').NE.0)
         PRNT2=(INDEX(KEYWRD,'MECI').NE.0)
         DEBUG=(DEBUG.AND.PRNT2)
         LARGE=(INDEX(KEYWRD,'LARGE').NE.0)
C C.I.=(m,n) keyword                                                    IR0494
         NDOUBL=99
         IF(INDEX(KEYWRD,'C.I.=(').NE.0)THEN
            NDOUBL=READA(KEYWRD,INDEX(KEYWRD,'C.I.=(')+7)
            NMOS=READA(KEYWRD,INDEX(KEYWRD,'C.I.=(')+5)
         ELSEIF (INDEX(KEYWRD,'C.I.=').NE.0)THEN
            NMOS=READA(KEYWRD,INDEX(KEYWRD,'C.I.=')+5)
         ELSE
            NMOS=NOPEN-NCLOSE
         ENDIF
C END C.I.=(m,n) keyword
         LROOT=1
         IF(INDEX(KEYWRD,'EXCI').NE.0)LROOT=2
         I=INDEX(KEYWRD,'ROOT')
         IF(I.NE.0)LROOT=READA(KEYWRD,I)
C C.I.=(m,n) keyword                                                    IR0494
         IF(NDOUBL.EQ.99)THEN
            J=MAX(MIN((NCLOSE+NOPEN+1)/2-(NMOS-1)/2,NORBS-NMOS+1),1)
         ELSE
            J=NCLOSE-NDOUBL+1
            IF(FRACT.GT.1.99D0)J=J+1
         ENDIF
C END C.I.=(m,n) keyword
         L=0
         DO 10 I=J,NCLOSE
            L=L+1
   10    OCCA(L)=1
         DO 20 I=NCLOSE+1,NOPEN
            L=L+1
   20    OCCA(L)=FRACT*0.5D0
         DO 30 I=NOPEN+1,J+NMOS-1
            L=L+1
   30    OCCA(L)=0.D0
C#         WRITE(6,'('' INITIAL ORBITAL OCCUPANCIES'')')
C#         WRITE(6,'(6F12.6)')(OCCA(L),L=1,NMOS)
         SING=(INDEX(KEYWRD,'SINGL')+
     1         INDEX(KEYWRD,'EXCI')+
     2         INDEX(KEYWRD,'BIRAD').NE.0)
         DOUB=(INDEX(KEYWRD,'DOUBL').NE.0)
         TRIP=(INDEX(KEYWRD,'TRIPL').NE.0)
         QUAR=(INDEX(KEYWRD,'QUART').NE.0)
         QUIN=(INDEX(KEYWRD,'QUINT').NE.0)
         SEXT=(INDEX(KEYWRD,'SEXTE').NE.0)
         SMULT=-.5D0
         IF(SING) SMULT=0.00D0
         IF(DOUB) SMULT=0.75D0
         IF(TRIP) SMULT=2.00D0
         IF(QUAR) SMULT=3.75D0
         IF(QUIN) SMULT=6.00D0
         IF(SEXT) SMULT=8.75D0
C#      WRITE(6,'('' ORBITAL COUNTERS, PER ATOM, FIRST LAST'')')
C#      WRITE(6,'(I40,I6)')(NFIRST(I),NLAST(I),I=1,NUMAT)
         X=0.D0
         DO 40 J=1,NMOS
   40    X=X+OCCA(J)
         XX=X+X
         NE=XX+0.5
         NELEC=(NELECS-NE+1)/2
         NLEFT=NORBS-NMOS-NELEC
      ENDIF
      IF(IRESET.EQ.1)THEN
         ONE=-1.D0
         GOTO 520
      ELSE
         ONE=1.D0
      ENDIF
      PRNT=(DEBUG.OR.FINISH.AND.PRNT2)
      BIGPRT=(PRNT.AND.LARGE)
      LAST1=(LAST1.OR.FINISH)
C
C    TEST TO SEE IF THE SET OF ENERGY LEVELS USED IN MECI IS COMPLETE,
C    I.E., ALL COMPONENTS OF DEGENERATE IRREDUCIBLE REPRESENTATIONS
C    ARE USED.  IF NOT, THEN RESULTS WILL BE NONSENSE.  GIVE USERS A
C    CHANCE TO REALLY FOUL THINGS UP BY ALLOWING JOB TO CONTINUE IF
C    'GEO-OK' IS SPECIFIED.
C
      DO 60 I=1,NMOS
         IN=I+NELEC
         DO 50 J=1,NORBS
   50    COEFFS(J,I)=COEFF(J,IN)
   60 EIGA(I)=EIGS(IN)
      LSPIN=(LSPIN1.AND. FINISH)
      IF(BIGPRT)THEN
         WRITE(6,'(''  INITIAL EIGENVALUES'')')
         WRITE(6,'(5F12.6)')(EIGA(I),I=1,NMOS)
         WRITE(6,'(//10X,''NUMBER OF ELECTRONS IN C.I. ='',F5.1)')XX
      ENDIF
      IF(.NOT.GEOOK.AND.NELEC.GT.0)THEN
         IF(ABS(EIGS(NELEC+1)-EIGS(NELEC)).LT.1.D-1.OR.
     1ABS(EIGS(NELEC+1+NMOS)-EIGS(NELEC+NMOS)).LT.1.D-1)THEN
            WRITE(6,'(///10X,A)')'DEGENERATE ENERGY LEVELS DETECTED IN M
     1ECI'
            WRITE(6,'(10X,A)')'SOME OF THESE LEVELS WOULD BE TREATED BY'
     1//'MECI,'
            WRITE(6,'(10X,A)')'WHILE OTHERS WOULD NOT.  THIS WOULD RESUL
     1T IN'
            WRITE(6,'(10X,A)')'NON-REPRODUCIBLE ELECTRONIC ENERGIES.'
            WRITE(6,'(10X,A)')'  JOB STOPPED.  TO CONTINUE, SPECIFY "GEO
     1-OK"'
            STOP
         ENDIF
      ENDIF
      KALPHA=(NE+1)/2
      KBETA=NE-KALPHA
C
C   FOR DEBUGGING PURPOSES, UNCOMMENT THE FOLLOWING LINES
C   THEY LOCK THE PHASE-FACTOR
C
C#      DO 55 I=1,NMOS
C#      SUM=0.D0
C#      DO 56 J=1,NORBS
C#  56  IF(ABS(SUM).LT.ABS(COEFFS(J,I)))SUM=COEFFS(J,I)
C#      IF(SUM.LT.0.D0)THEN
C#      DO 57 J=1,NORBS
C#  57  COEFFS(J,I)=-COEFFS(J,I)
C#      ENDIF
C#  55  CONTINUE
C
C   FOR DEBUGGING PURPOSES, UNCOMMENT THE PRECEEDING LINES
C
      IF( BIGPRT ) THEN
         WRITE(6,'(//10X,''EIGENVECTORS'',/)')
         DO 70 I=1,NORBS
   70    WRITE(6,'(6F12.6)')(COEFFS(I,J),J=1,NMOS)
      ENDIF
      DO 80 I=1,NMOS
         DO 80 J=1,NMOS
            DO 80 K=1,NMOS
               DO 80 L=1,NMOS
   80 XY(I,J,K,L)=100.0D0
      NFA(2)=1
      NFA(1)=1
      DO 90 I=3,NMECI+1
   90 NFA(I)=NFA(I-1)*(I-1)
      DO 100 I=1,NMOS
         DO 100 J=1,NMOS
            DO 100 K=1,NMOS
               DO 100 L=1,NMOS
  100 XY(I,J,K,L)=100.0D0
      DO 110 I=1,NMOS
         I1=I
         DO 110 J=1,I
            J1=J
            DO 110 K=1,NMOS
               K1=K
               DO 110 L=1,K
                  L1=L
  110 CALL IJKL(I1,K1,J1,L1,X,COEFFS,NORBS)
      DO 120 I=1,NMOS
         DO 120 J=1,NMOS
            RJKAA(I,J)=XY(I,I,J,J)-XY(I,J,I,J)
  120 RJKAB(I,J)=XY(I,I,J,J)
      DO 140 I=1,NMOS
         X=0.0
         DO 130 J=1,NMOS
            X=X+(RJKAA(I,J)+RJKAB(I,J))*OCCA(J)
  130    CONTINUE
         EIGA(I)=EIGA(I)-X
C#      IF(ABS(OCCA(I)-0.5).LT.1.D-4)EIGA(I)=EIGA(I)+XY(I,I,I,I)*0.25D0
  140 CONTINUE
      IF(BIGPRT) THEN
         WRITE(6,150)
  150    FORMAT(/,5X,'EIGENVALUES AFTER REMOVAL OF INTER-ELECTRONIC INTE
     1RACTIONS',/)
         WRITE(6,'(6F12.6)')(EIGA(I),I=1,NMOS)
         WRITE(6,'(///10X,''TWO-ELECTRON J-INTEGRALS'',/)')
         DO 160 I1=1,NMOS
  160    WRITE(6,'(10F10.4)')(RJKAB(I1,J1),J1=1,NMOS)
         WRITE(6,'(///10X,''TWO-ELECTRON K-INTEGRALS'',/)')
         DO 170 I1=1,NMOS
  170    WRITE(6,'(10F10.4)')(RJKAB(I1,J1)-RJKAA(I1,J1),J1=1,NMOS)
      ENDIF
      NATOMS=NUMAT
      DO 180 I=1,NMOS
         DO 180 J=1,NMOS
            RJKAA(I,J)=RJKAA(I,J)*0.5D0
  180 CONTINUE
         IF (ICALC1 .NE. NUMCAL) THEN
         ICALC1 = NUMCAL
         I=INDEX(KEYWRD,'MICROS')
         IF(I.NE.0)THEN
            K=READA(KEYWRD,I)
            LAB=K
            IF(PRNT)WRITE(6,'(''    MICROSTATES READ IN'')')
            NTOT=XX+0.5
            OPEN(UNIT=IREAD,FILE='FOR004',STATUS='OLD',BLANK='ZERO')
            REWIND IREAD
            DO 190 I=1,3
  190       READ(IREAD,'(A)')LINE
            DO 200 I=1,1000
               READ(IREAD,'(A)')LINE
  200       IF(INDEX(LINE,'MICRO').NE.0)GOTO 210
  210       DO 240 I=1,LAB
               READ(IREAD,'(A)')LINE
               IZERO=MAX(0,MIN(INDEX(LINE,'0'),INDEX(LINE,'1'))-1)
               DO 220 J=1,NMOS
                  IF(LINE(J+IZERO:J+IZERO).NE.'1')
     1            LINE(J+IZERO:J+IZERO)='0'
                  IF(LINE(J+NMOS+IZERO:J+NMOS+IZERO).NE.'1')
     1            LINE(J+NMOS+IZERO:J+NMOS+IZERO)='0'
                  MICROA(J,I)=ICHAR(LINE(J+IZERO:J+IZERO))-
     1          ICHAR('0')
                  MICROB(J,I)=ICHAR(LINE(J+NMOS+IZERO:J+NMOS+IZERO))-
     1          ICHAR('0')
  220          CONTINUE
               IF(PRNT)WRITE(6,'(20I6)')(MICROA(J,I),J=1,NMOS),
     1        (MICROB(J,I),J=1,NMOS)
               K=0
               DO 230 J=1,NMOS
  230          K=K+MICROA(J,I)+MICROB(J,I)
               IF(K.NE.NTOT)THEN
                  NTOT=K
                  XX=K
                  WRITE(6,'(/,''NUMBER OF ELECTRONS IN C.I. REDEFINED TO
     1:'',I4,/)')K
               ENDIF
  240       CONTINUE
            GOTO 310
         ENDIF
         NUPP=KALPHA
         NDOWN=KBETA
         AMS=(NUPP-NDOWN)*0.5D0
         IF(PRNT)WRITE(6,250) AMS
  250    FORMAT(10X,'COMPONENT OF SPIN  = ',F4.1)
         IF(NUPP*NDOWN.GE.0) GOTO 270
         WRITE(6,260)
  260    FORMAT(/10X,28H IMPOSSIBLE VALUE OF DELTA S/)
         STOP
  270    LIMA=NFA(NMOS+1)/(NFA(NUPP+1)*NFA(NMOS-NUPP+1))
         LIMB=NFA(NMOS+1)/(NFA(NDOWN+1)*NFA(NMOS-NDOWN+1))
         LAB=LIMA*LIMB
         IF(PRNT)WRITE(6,280) LAB
  280    FORMAT(//10X,35H NO OF CONFIGURATIONS CONSIDERED = ,I4)
C#      IF(LAB.LT.101) GOTO 240
C#      WRITE(6,230)
C#  230 FORMAT(10X,24H TOO MANY CONFIGURATIONS/)
C#      GOTO 160
C#  240 CONTINUE
         CALL PERM(NPERMA, NUPP, NMOS, NMECI, LIMA)
         CALL PERM(NPERMB, NDOWN, NMOS, NMECI, LIMB)
         K=0
         DO 290 I=1,LIMA
            DO 290 J=1,LIMB
               K=K+1
               DO 290 L=1,NMOS
                  MICROA(L,K)=NPERMA(L,I)
  290    MICROB(L,K)=NPERMB(L,J)
  300    FORMAT(10I1)
  310    CONTINUE
         LIMA=LAB
         LIMB=LAB
      ENDIF
      GSE=0.0D0
      DO 320 I=1,NMOS
C#         IF(ABS(OCCA(I)-0.5).LT.0.01)GSE=GSE-0.25D0*XY(I,I,I,I)
         GSE=GSE+EIGA(I)*OCCA(I)*2.D0
         GSE=GSE+XY(I,I,I,I)*OCCA(I)*OCCA(I)
         DO 320 J=I+1,NMOS
  320 GSE=GSE+2.D0*(2.D0*XY(I,I,J,J) - XY(I,J,I,J))*OCCA(I)*OCCA(J)
C#    IF(PRNT)WRITE(6,'('' GROUND STATE ENERGY:'',F13.6,'' E.V.'')')GSE
C     ..........
      IF( PRNT )
     1WRITE(6,'(//10X,''CONFIGURATIONS CONSIDERED IN C.I.'',//)')
      J=0
      DO 330 I=1,LAB
         DIAG(I)=DIAGI(MICROA(1,I),MICROB(1,I),EIGA,XY,NMOS)-GSE
  330 CONTINUE
  340 CONTINUE
      IF(LAB.LE.MDIM) GOTO 380
      X=-100.D0
      DO 350 I=1,LAB
         IF(DIAG(I).GT.X)THEN
            X=DIAG(I)
            J=I
         ENDIF
  350 CONTINUE
      IF(J.NE.LAB) THEN
         DO 370 I=J,LAB
            I1=I+1
            DO 360 K=1,NMOS
               MICROA(K,I)=MICROA(K,I1)
  360       MICROB(K,I)=MICROB(K,I1)
  370    DIAG(I)=DIAG(I1)
      ENDIF
      LAB=LAB-1
      GOTO 340
  380 CONTINUE
C
C   BEFORE STARTING, CHECK THAT THE ROOT WANTED CAN EXIST
C
      IF(LAB.LT.LROOT)THEN
         WRITE(6,'(//10X,''C.I. IS OF SIZE LESS THAN ROOT SPECIFIED'')')
         WRITE(6,'(10X,''MODIFY SIZE OF C.I. OR ROOT NUMBER'')')
         STOP
      ENDIF
C
C  MAIN LOOP TO FILL SECULAR DETERMINANT
C
      IK=0
      DO 430 I=1,LAB
         K=0
         DO 390 J=1,NMOS
            IOCCB1(J)=MICROB(J,I)
            IOCCA1(J)=MICROA(J,I)
  390    K=K+IOCCA1(J)
         NALPHA(I)=K
         IF(PRNT)  THEN
            WRITE(6,'(/10X,I4,6X,10I4)') I,(IOCCA1(K),K=1,NMOS)
            WRITE(6,'(20X,10I4)')(IOCCB1(K),K=1,NMOS)
         ENDIF
         IS=2
C
C   INNER LOOP TO FILL SECULAR DETERMINANT
C
         DO 420 K=1,I
            IK=IK+1
            CIMAT(IK)=I*1.D-15+J*.1D-17
            IX=0
            IY=0
            DO 400 J=1,NMOS
               IOCCB2(J)=MICROB(J,K)
               IOCCA2(J)=MICROA(J,K)
               IX=IX+ABS(IOCCA1(J)-IOCCA2(J))
  400       IY=IY+ABS(IOCCB1(J)-IOCCB2(J))
C
C                              CHECK IF MATRIX ELEMENT HAS TO BE ZERO
C
            IF(IX+IY.GT.4 .OR. NALPHA(I).NE.NALPHA(K))  GOTO 420
            IF(IX+IY.EQ.4) THEN
               IF(IX.EQ.0)THEN
                  CIMAT(IK)=BABBCD(IOCCA1, IOCCB1, IOCCA2, IOCCB2, NMOS)
               ELSEIF(IX.EQ.2)THEN
                  CIMAT(IK)=AABBCD(IOCCA1, IOCCB1, IOCCA2, IOCCB2, NMOS)
               ELSE
                  CIMAT(IK)=AABACD(IOCCA1, IOCCB1, IOCCA2, IOCCB2, NMOS)
               ENDIF
            ELSEIF(IX.EQ.2)THEN
               CIMAT(IK)=AABABC(IOCCA1, IOCCB1, IOCCA2, IOCCB2, NMOS)
            ELSEIF(IY.EQ.2)THEN
               CIMAT(IK)=BABBBC(IOCCA1, IOCCB1, IOCCA2, IOCCB2, NMOS)
            ELSE
               CIMAT(IK)=DIAG(I)
               X=0.0D0
               DO 410 J=1,NMOS
  410          X=X+IOCCA1(J)*IOCCB1(J)
               SPIN(I)=(-(XX-2*NALPHA(I))*(XX-2*NALPHA(I))+X*4.D0)
            ENDIF
  420    CONTINUE
         ISPQR(I,1)=IS-1
  430 CONTINUE
      IF(BIGPRT)THEN
         WRITE(6,'(//,'' C.I. MATRIX'')')
         CALL VECPRT(CIMAT,-LAB)
      ELSE
         IF(PRNT)WRITE(6,'(//,'' DIAGONAL OF C.I. MATRIX'')')
         IF(PRNT)WRITE(6,'(5F13.6)')(CIMAT((I*(I+1))/2),I=1,LAB)
      ENDIF
      CALL RSP(CIMAT,LAB,LAB,EIG,CONF)
C
C   DECIDE WHICH ROOT TO EXTRACT
C
      KROOT=0
      IF(SMULT.LT.-0.1D0)THEN
         MECI=EIG(LROOT)
         KROOT=LROOT
      ENDIF
      IF(PRNT)  THEN                                                    IR0594
         WRITE(6,'(//20X,''STATE VECTORS'',//)')
         CALL MATOUT(CONF,EIG,LAB,LAB,LAB)
      ENDIF
      IF(PRNT)THEN
         WRITE(6,440)
  440    FORMAT(///,' STATE ENERGIES '
     1,' EXPECTATION VALUE OF S**2  S FROM S**2=S(S+1)',//)
      ENDIF
      IROOT=0
      DO 450 I=1,9
  450 CIMAT(I)=0.1D0
      DO 490 I=1,LAB
         X=0.5D0*XX
         II=(I-1)*LAB
         DO 480 J=1,LAB
            JI=J+II
            X=X-CONF(JI)*CONF(JI)*SPIN(J)*0.25D0
            K=ISPQR(J,1)
            IF(K.EQ.1)  GOTO  470
            DO 460 K=2,K
               LI=ISPQR(J,K)+II
  460       X=X+CONF(JI)*CONF(LI)*2.D0
  470       CONTINUE
  480    CONTINUE
         Y=(-1.D0+SQRT(1.D0+4.D0*X))*0.5D0
         IF(ABS(SMULT-X).LT.0.01)THEN
            IROOT=IROOT+1
            IF(IROOT.EQ.LROOT) THEN
               KROOT=I
               MECI=EIG(I)
            ENDIF
         ENDIF
         J=Y*2.D0+1.5D0
         CIMAT(J)=CIMAT(J)+1
  490 IF(PRNT)WRITE(6,510) I,EIG(I),TSPIN(J),X,Y
      IF(KROOT.EQ.0)THEN
         WRITE(6,'(//10X,''THE STATE REQUIRED IS NOT PRESENT IN THE'')')
         WRITE(6,'(10X,  ''    SET OF CONFIGURATIONS AVAILABLE'')')
         WRITE(6,'(/ 4X,''NUMBER OF STATES ACCESSIBLE USING CURRENT KEY-
     1WORDS'',/)')
         DO 500 I=1,7
  500    IF(CIMAT(I).GT.0.5D0)
     1WRITE(6,'((24X,A8,I4))')TSPIN(I),NINT(CIMAT(I))
         STOP
      ENDIF
  510 FORMAT(I5,F12.6,3X,A8,F15.5,F10.5)
C#      M=0
C#      DO 440 I=1,NMOS
C#         WRITE(6,*)
C#         DO 440 J=1,NMOS
C#            WRITE(6,*)
C#            DO 440 K=1,NMOS
C#  440 WRITE(6,'(4I2,8F12.6)')I,J,K,M,(XY(I,J,K,L),L=1,NMOS)
  520 CONTINUE
      IF(FORCE.OR.LAST1)THEN
C
C   REFORM DENSITY MATRIX
C
         K=(KROOT-1)*LAB
         DO 540 I=1,NMOS
            SUM=0.D0
            DO 530 J=1,LAB
  530       SUM=SUM+(MICROA(I,J)+MICROB(I,J))*CONF(J+K)**2
  540    EIGA(I)=SUM-OCCA(I)*2.D0
         L=0
         DO 560 I=1,NORBS
            DO 560 J=1,I
               SUM=0.D0
               DO 550 K=1,NMOS
  550          SUM=SUM+EIGA(K)*COEFFS(I,K)*COEFFS(J,K)
               L=L+1
  560    P(L)=P(L)+SUM*ONE
      ENDIF
      MAXVEC=0
      IF(LSPIN)MAXVEC=MIN(4,LAB)
      IF(LSPIN.AND.(NE/2)*2.EQ.NE) THEN
         WRITE(6,'(''   ESR SPECIFIED FOR AN EVEN-ELECTRON SYSTEM'')')
      ENDIF
      DO 570 I=1,NMOS
         DO 570 J=1,NORBS
  570 COEFFS(J,I)=COEFFS(J,I)**2
      DO 640 IUJ=1,MAXVEC
         IOFSET=(IUJ-1)*LAB
         WRITE(6,'(//,''      MICROSTATE CONTRIBUTIONS TO '',
     1''STATE EIGENFUNCTION'',I3)')IUJ
         WRITE(6,'(5F13.6)')(CONF(I+IOFSET),I=1,LAB)
         DO 580 I=1,LAB
  580    CONF(I)=CONF(I+IOFSET)**2
C                                             SECOND VECTOR!
         DO 600 I=1,NMOS
            SUM=0.D0
            DO 590 J=1,LAB
  590       SUM=SUM+(MICROA(I,J)-MICROB(I,J))*CONF(J)
  600    EIGA(I)=SUM
         WRITE(6,'(/,''    SPIN DENSITIES FROM EACH M.O., ENERGY:''
     1,F7.3)')EIG(IUJ)
         WRITE(6,'(5F12.6)') (EIGA(I),I=1,NMOS)
         WRITE(6,*)
         WRITE(6,*)'     SPIN DENSITIES FROM EACH ATOMIC ORBITAL'
         WRITE(6,*)'                              S        PX        '//
     1'PY        PZ        TOTAL'
         DO 630 I=1,NATOMS
            IL=NFIRST(I)
            IU=NLAST(I)
            L=0
            SUMM=0.D0
            DO 620 K=IL,IU
               L=L+1
               SUM=0.D0
               DO 610 J=1,NMOS
  610          SUM=SUM+COEFFS(K,J)*EIGA(J)
               SUMM=SUMM+SUM
  620       EIGS(L)=SUM
            IF(L.EQ.4)THEN
               WRITE(6,'(''  ATOM'',I4,''    SPIN DENSITY  '',5F10.7)')
     1I,(EIGS(K),K=1,L),SUMM
            ELSE
               WRITE(6,'(''  ATOM'',I4,''    SPIN DENSITY  '',F10.7,30X,
     1F10.7)')I,EIGS(1),SUMM
            ENDIF
  630    CONTINUE
  640 CONTINUE
      RETURN
      END
C
      SUBROUTINE MOLDAT(IP)                                             0304WH93
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1                NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /MOLMEC/ HTYPE(4),NHCO(4,20),NNHCO,ITYPE,USEMM
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
     3       /MOLKSR/ FRACT
     4       /KEYWRD/ KEYWRD
     5       /NATORB/ NATORB(107)
      COMMON /CORE  / CORE(107)
     1       /BETAS / BETAS(107),BETAP(107),BETAD(107)
     2       /MOLORB/ USPD(MAXORB),PSPD(MAXORB)
     3       /VSIPS / VS(107),VP(107),VD(107)
     4       /ONELEC/ USS(107),UPP(107),UDD(107)
      COMMON /ATHEAT/ ATHEAT
     1       /POLVOL/ POLVOL(107)
     2       /MULTIP/ DD(107),QQ(107),AM(107),AD(107),AQ(107)
     3       /TWOELE/ GSS(107),GSP(107),GPP(107),GP2(107),HSP(107)
     4                ,GSD(107),GPD(107),GDD(107)
     5       /IDEAS / GUESS1(107,10),GUESS2(107,10),GUESS3(107,10)
     6       /IDEAP / GUESP1(107,10),GUESP2(107,10),GUESP3(107,10)
      COMMON /ALPHA / ALP(107)
     1       /REFS/ ALLREF(107,4)
      COMMON /MNDO/  USSM(107), UPPM(107), UDDM(107), ZSM(107),
     1ZPM(107), ZDM(107), BETASM(107), BETAPM(107), BETADM(107),
     2ALPM(107), EISOLM(107), DDM(107), QQM(107), AMM(107),
     3ADM(107), AQM(107), GSSM(107), GSPM(107), GPPM(107),
     4GP2M(107), HSPM(107), POLVOM(107)
      COMMON /PM3 /  USSPM3(107), UPPPM3(107), UDDPM3(107), ZSPM3(107),
     1ZPPM3(107), ZDPM3(107), BETASP(107), BETAPP(107), BETADP(107),
     2ALPPM3(107), EISOLP(107), DDPM3(107), QQPM3(107), AMPM3(107),
     3ADPM3(107), AQPM3(107) ,GSSPM3(107), GSPPM3(107), GPPPM3(107),
     4GP2PM3(107), HSPPM3(107),POLVOP(107)
      COMMON /GEOM  / GEO(3,NUMATM)
      PARAMETER (MDUMY=MAXPAR**2-MPACK)
      COMMON /SCRACH/ RXYZ(MPACK), XDUMY(MDUMY)
*
*  COMMON BLOCKS FOR MINDO/3
*
      COMMON /ONELE3 /  USS3(18),UPP3(18)
     1       /ATOMI3 /  EISOL3(18),EHEAT3(18)
     2       /EXPON3 /  ZS3(18),ZP3(18)
*
*  END OF MINDO/3 COMMON BLOCKS
*
      COMMON /EXPONT/ ZS(107),ZP(107),ZD(107)
      COMMON /ATOMIC/ EISOL(107),EHEAT(107)
      DIMENSION COORD(3,NUMATM), ISWAP(2,20)
      character allrf2
      CHARACTER*80 KEYWRD, OLDE(20)*6, ALLREF
      LOGICAL DEBUG, UHF,EXCI, TRIP, MINDO3, BIRAD, AM1, OPEN, LPM3,
     1        USEMM,SLOW,HALFE                                          IR0494
         SAVE                                                           GL0892
C
      IF (IP .EQ. 1) GOTO 35                                            0304WH93
C     
      DEBUG = (INDEX(KEYWRD,'MOLDAT').NE.0)
      LPM3  = (INDEX(KEYWRD,'PM3').NE.0)
      MINDO3= (INDEX(KEYWRD,'MINDO').NE.0)
      UHF=(INDEX(KEYWRD,'UHF') .NE. 0)
      AM1= (INDEX(KEYWRD,'AM1').NE.0)
      KHARGE=0
      I=INDEX(KEYWRD,'CHARGE')
      IF(I.NE.0) KHARGE=READA(KEYWRD,I)
      ELECS=-KHARGE
      NDORBS=0
      ATHEAT=0.D0
      EAT=0.D0
      NUMAT=0
      IF (   .NOT.   AM1  .AND. .NOT. LPM3 ) THEN
*
*    SWITCH IN MNDO PARAMETERS
*
C
C       ZERO OUT GAUSSIAN 1 FOR CARBON.  THIS WILL BE USED IN
C       ROTATE TO DECIDE WHETHER OR NOT TO USE AM1-TYPE GAUSSIANS
C
         GUESS1(6,1)=0.D0
         DO 10 I=1,107
            IF(.NOT.MINDO3) POLVOL(I)=POLVOM(I)
            ZS(I)=ZSM(I)
            ZP(I)=ZPM(I)
            ZD(I)=ZDM(I)
            USS(I)=USSM(I)
            UPP(I)=UPPM(I)
            UDD(I)=UDDM(I)
            BETAS(I)=BETASM(I)
            BETAP(I)=BETAPM(I)
            BETAD(I)=BETADM(I)
            ALP(I)=ALPM(I)
            EISOL(I)=EISOLM(I)
            DD(I)=DDM(I)
            QQ(I)=QQM(I)
            AM(I)=AMM(I)
            AD(I)=ADM(I)
            AQ(I)=AQM(I)
            GSS(I)=GSSM(I)
            GPP(I)=GPPM(I)
            GSP(I)=GSPM(I)
            GP2(I)=GP2M(I)
            HSP(I)=HSPM(I)
   10    CONTINUE
      ELSEIF( .NOT. AM1 .AND. LPM3) THEN
*
*    SWITCH IN MNDO-PM3 PARAMETERS
*
         DO 30 I=1,107
            DO 20 J=1,10
               GUESS1(I,J)=GUESP1(I,J)
               GUESS2(I,J)=GUESP2(I,J)
   20       GUESS3(I,J)=GUESP3(I,J)
            POLVOL(I)=POLVOP(I)
            ZS(I)=ZSPM3(I)
            ZP(I)=ZPPM3(I)
            ZD(I)=ZDPM3(I)
            USS(I)=USSPM3(I)
            UPP(I)=UPPPM3(I)
            UDD(I)=UDDPM3(I)
            BETAS(I)=BETASP(I)
            BETAP(I)=BETAPP(I)
            BETAD(I)=BETADP(I)
            ALP(I)=ALPPM3(I)
            EISOL(I)=EISOLP(I)
            DD(I)=DDPM3(I)
            QQ(I)=QQPM3(I)
            AM(I)=AMPM3(I)
            AD(I)=ADPM3(I)
            AQ(I)=AQPM3(I)
            GSS(I)=GSSPM3(I)
            GPP(I)=GPPPM3(I)
            GSP(I)=GSPPM3(I)
            GP2(I)=GP2PM3(I)
            HSP(I)=HSPPM3(I)
   30    CONTINUE
      ENDIF
C
      IF (IP .EQ. 0) RETURN                                             0304WH93
C
C        SWAP IN OLD PARAMETERS FOR ELEMENTS.  OLDE CONTAINS THE
C        CHARACTER NAME OF THE ELEMENT, AND ISWAP(1,1:NEWELE) CONTAINS
C        THE ATOMIC NUMBER OF THE ELEMENT. ISWAP(2,1:NEWELE) CONTAINS
C        THE STORAGE ADDRESS OF THE OLD SET OF PARAMETERS.
C
35    NEWELE=2                                                          0304WH93
      OLDE(1)=' S1978'
      ISWAP(1,1)=16
      ISWAP(2,1)=91
      OLDE(2)='SI1978'
      ISWAP(1,2)=14
      ISWAP(2,2)=90
      DO 40 K=1,NEWELE
         IF(INDEX(KEYWRD,OLDE(K)).NE.0)THEN
            I=ISWAP(1,K)
            J=ISWAP(2,K)
      allrf2 = allref(j,1)
      allref(i,3) = allrf2
      allref(i,1) = allrf2
            ZS(I)=ZS(J)
            ZP(I)=ZP(J)
            ZD(I)=ZD(J)
            USS(I)=USS(J)
            UPP(I)=UPP(J)
            UDD(I)=UDD(J)
            BETAS(I)=BETAS(J)
            BETAP(I)=BETAP(J)
            BETAD(I)=BETAD(J)
            ALP(I)=ALP(J)
            EISOL(I)=EISOL(J)
            DD(I)=DD(J)
            QQ(I)=QQ(J)
            AM(I)=AM(J)
            AD(I)=AD(J)
            AQ(I)=AQ(J)
            IF(GSS(J).NE.0)GSS(I)=GSS(J)
            IF(GPP(J).NE.0)GPP(I)=GPP(J)
            IF(GSP(J).NE.0)GSP(I)=GSP(J)
            IF(GP2(J).NE.0)GP2(I)=GP2(J)
            IF(HSP(J).NE.0)HSP(I)=HSP(J)
         ENDIF
   40 CONTINUE
      IF( MINDO3 ) THEN
         DO 50 I=1,17
            IF(I.EQ.2.OR.I.EQ.10)GOTO 50
            USS(I)=USS3(I)
            UPP(I)=UPP3(I)
            EISOL(I)=EISOL3(I)
            EHEAT(I)=EHEAT3(I)
            ZS(I)=ZS3(I)
            ZP(I)=ZP3(I)
   50    CONTINUE
      ENDIF
      IF(USS(1) .GT. -1.D0) THEN
         WRITE(6,'(''  THE HAMILTONIAN REQUESTED IS NOT AVAILABLE IN''
     1,'' THIS PROGRAM'')')
         STOP
      ENDIF
      IA=1
      IB=0
      NHEAVY=0
      DO 100 II=1,NATOMS
         IF(LABELS(II).EQ.99.OR.LABELS(II).EQ.107) GOTO 100
         NUMAT=NUMAT+1
         NAT(NUMAT)=LABELS(II)
         NFIRST(NUMAT)=IA
         NI=NAT(NUMAT)
         ATHEAT=ATHEAT+EHEAT(NI)
         EAT   =EAT   +EISOL(NI)
         ELECS=ELECS+CORE(NI)
         IB=IA+NATORB(NI)-1
         NMIDLE(NUMAT)=IB
         IF(NATORB(NI).EQ.9)NDORBS=NDORBS+5
         IF(NATORB(NI).EQ.9)NMIDLE(NUMAT)=IA+3
         NLAST(NUMAT)=IB
         IF(IA.GT.MAXORB) GOTO 240
         USPD(IA)=USS(NI)
         IF(IA.EQ.IB) GOTO 90
         K=IA+1
         K1=IA+3
         DO 60 J=K,K1
            IF(J.GT.MAXORB) GOTO 240
            USPD(J)=UPP(NI)
   60    CONTINUE
         NHEAVY=NHEAVY+1
   70    IF(K1.EQ.IB)GOTO 90
         K=K1+1
         DO 80 J=K,IB
   80    USPD(J)=UDD(NI)
   90    CONTINUE
  100 IA=IB+1
      CALL REFER
      ATHEAT=ATHEAT-EAT*23.061D0
      NORBS=NLAST(NUMAT)
      IF(NORBS.GT.MAXORB)THEN
         WRITE(6,'(//10X,''**** MAX. NUMBER OF ORBITALS:'',I4,/
     1            10X,''NUMBER OF ORBITALS IN SYSTEM:'',I4)')
     2MAXORB,NORBS
         STOP
      ENDIF
      NLIGHT=NUMAT-NHEAVY
      N2EL=50*NHEAVY*(NHEAVY-1)+10*NHEAVY*NLIGHT+(NLIGHT*(NLIGHT-1))/2
      IF(N2EL.GT.N2ELEC)THEN
         WRITE(6,'(//10X,''**** MAX. NUMBER OF TWO-ELECTRON INTEGRALS:''
     1,I8,/
     2            10X,''NUMBER OF TWO ELECTRON INTEGRALS IN SYSTEM:'',
     3I8)')
     4N2ELEC,N2EL
         STOP
      ENDIF
C
C   NOW TO CALCULATE THE NUMBER OF LEVELS OCCUPIED
      TRIP=(INDEX(KEYWRD,'TRIPLET').NE.0)
      EXCI=(INDEX(KEYWRD,'EXCITED').NE.0)
      BIRAD=(EXCI.OR.INDEX(KEYWRD,'BIRAD').NE.0)
      IF(INDEX(KEYWRD,'C.I.') .NE. 0 .AND. UHF ) THEN
         WRITE(6,'(//10X,''C.I. NOT ALLOWED WITH UHF '')')
         STOP
      ENDIF
C
C NOW TO WORK OUT HOW MANY ELECTRONS ARE IN EACH TYPE OF SHELL
C
      NALPHA=0
      NBETA=0
C
C      PROTECT DUMB USERS FROM DUMB ERRORS!
C
      NELECS=MAX(ELECS+0.5D0,0.D0)
      NELECS=MIN(2*NORBS,NELECS)
      NCLOSE=0
      NOPEN=0
      IF( UHF ) THEN
         FRACT=1.D0
         NBETA=NELECS/2
         IF( TRIP ) THEN
            IF(NBETA*2 .NE. NELECS) THEN
               WRITE(6,'(//10X,''TRIPLET SPECIFIED WITH ODD NUMBER'',
     1            '' OF ELECTRONS, CORRECT FAULT '')')
               STOP
            ELSE
               WRITE(6,'(//'' TRIPLET STATE CALCULATION'')')
               NBETA=NBETA-1
            ENDIF
         ENDIF
         IF(INDEX(KEYWRD,'QUART').NE.0) THEN
            IF(NBETA*2 .EQ. NELECS) THEN
               WRITE(6,'(//10X,''QUARTET SPECIFIED WITH EVEN NUMBER'',
     1            '' OF ELECTRONS, CORRECT FAULT '')')
               STOP
            ELSE
               WRITE(6,'(//'' QUARTET STATE CALCULATION'')')
               NBETA=NBETA-1
            ENDIF
         ENDIF
         IF(INDEX(KEYWRD,'QUINT').NE.0) THEN
            IF(NBETA*2 .NE. NELECS) THEN
               WRITE(6,'(//10X,''QUINTET SPECIFIED WITH ODD NUMBER'',
     1            '' OF ELECTRONS, CORRECT FAULT '')')
               STOP
            ELSE
               WRITE(6,'(//'' QUINTET STATE CALCULATION'')')
               NBETA=NBETA-2
            ENDIF
         ENDIF
         IF(INDEX(KEYWRD,'SEXT').NE.0) THEN
            IF(NBETA*2 .EQ. NELECS) THEN
               WRITE(6,'(//10X,''SEXTET SPECIFIED WITH EVEN NUMBER'',
     1            '' OF ELECTRONS, CORRECT FAULT '')')
               STOP
            ELSE
               WRITE(6,'(//'' SEXTET STATE CALCULATION'')')
               NBETA=NBETA-1
            ENDIF
         ENDIF
         NALPHA=NELECS-NBETA
         WRITE(6,'(//10X,''UHF CALCULATION, NO. OF ALPHA ELECTRONS ='',I
     13,/27X,''NO. OF BETA  ELECTRONS ='',I3)')NALPHA,NBETA
      ELSE
C
C   NOW TO DETERMINE OPEN AND CLOSED SHELLS
C
         OPEN=.FALSE.
         IELEC=0
         ILEVEL=0
         IF( TRIP .OR. EXCI .OR. BIRAD ) THEN
            IF( (NELECS/2)*2 .NE. NELECS) THEN
               WRITE(6,'(//10X,''SYSTEM SPECIFIED WITH ODD NUMBER'',
     1            '' OF ELECTRONS, CORRECT FAULT '')')
               STOP
            ENDIF
            IF(BIRAD)WRITE(6,'(//'' SYSTEM IS A BIRADICAL'')')
            IF(TRIP )WRITE(6,'(//'' TRIPLET STATE CALCULATION'')')
            IF(EXCI )WRITE(6,'(//'' EXCITED STATE CALCULATION'')')
            IELEC=2
            ILEVEL=2
         ELSEIF((NELECS/2)*2.NE.NELECS) THEN
            IELEC=1
            ILEVEL=1
         ENDIF
         IF(INDEX(KEYWRD,'QUART').NE.0) THEN
            WRITE(6,'(//'' QUARTET STATE CALCULATION'')')
            IELEC=3
            ILEVEL=3
         ENDIF
         IF(INDEX(KEYWRD,'QUINT').NE.0) THEN
            WRITE(6,'(//'' QUINTET STATE CALCULATION'')')
            IELEC=4
            ILEVEL=4
         ENDIF
         IF(INDEX(KEYWRD,'SEXT').NE.0) THEN
            WRITE(6,'(//'' SEXTET STATE CALCULATION'')')
            IELEC=5
            ILEVEL=5
         ENDIF
         I=INDEX(KEYWRD,'OPEN(')
         IF(I.NE.0)THEN
            IELEC=READA(KEYWRD,I)
            ILEVEL=READA(KEYWRD,I+7)
         ENDIF
         NCLOSE=NELECS/2
         NOPEN = NELECS-NCLOSE*2
         IF( IELEC.NE.0 )THEN
            IF((NELECS/2)*2.EQ.NELECS .NEQV.
     1                  (IELEC/2)*2.EQ.IELEC) THEN
               WRITE(6,'('' IMPOSSIBLE NUMBER OF OPEN SHELL ELECTR
     1ONS'')')
               STOP
            ENDIF
            NCLOSE=NCLOSE-IELEC/2
            NOPEN=ILEVEL
            FRACT=IELEC*1.D0/ILEVEL
            WRITE(6,'('' THERE ARE'',I3,'' DOUBLY FILLED LEVELS'')
     1')NCLOSE
         ENDIF
         WRITE(6,'(//10X,''RHF CALCULATION, NO. OF '',
     1''DOUBLY OCCUPIED LEVELS ='',I3)')NCLOSE
         IF(NOPEN.NE.0.AND.ABS(FRACT-1.D0).LT.1.D-4)
     1WRITE(6,'(/27X,''NO. OF SINGLY OCCUPIED LEVELS ='',I3)')NOPEN
         IF(NOPEN.NE.0.AND.ABS(FRACT-1.D0).GT.1.D-4)
     1WRITE(6,'(/27X,''NO. OF LEVELS WITH OCCUPANCY'',F6.3,''  ='',I3)')
     2FRACT,NOPEN
C
C C.I=(n,m) modifications                                               IR0494
C
         IF(INDEX(KEYWRD,'C.I.=(').NE.0) THEN
            I=READA(KEYWRD,INDEX(KEYWRD,'C.I.=(')+5)-
     1      READA(KEYWRD,INDEX(KEYWRD,'C.I.=(')+7)
            IF(NOPEN.GT.I)THEN
               WRITE(6,'(//,'' NUMBER OF OPEN-SHELLS ALLOWED IN C.I. IS
     1LESS ''            /''    THAN THAT SPECIFIED BY OTHER KEYWORDS'')
     2')
               STOP
            ENDIF
         ENDIF
         IF(INDEX(KEYWRD,'C.I.').NE.0.AND.NOPEN.EQ.0)THEN
            NOPEN=1
            NCLOSE=NCLOSE-1
            FRACT=2.D0
         ENDIF
C
C END of  C.I=(n,m) modifications
C
         NOPEN=NOPEN+NCLOSE
      ENDIF
C other C.I.=(n,m) checks                                               IR0494
C
C  WORK OUT IF DEFINED SPIN-STATE ALLOWED
C
      MSDEL=INDEX(KEYWRD,' MS')
      IF(MSDEL.NE.0)THEN
         MSDEL=1.0001D0*READA(KEYWRD,INDEX(KEYWRD,' MS'))
      ELSE
         IF(INDEX(KEYWRD,'TRIP')+INDEX(KEYWRD,'QUAR').GT.0)MSDEL=1
         IF(INDEX(KEYWRD,'QUIN')+INDEX(KEYWRD,'SEXT').GT.0)MSDEL=2
      ENDIF
      IF(MSDEL.NE.0.AND..NOT.UHF)THEN
C
C   MSDEL = NUMBER OF ALPHA ELECTRONS - NUMBER OF BETA ELECTRONS
C
         NDOUBL=99
         IF(INDEX(KEYWRD,'C.I.=(').NE.0)THEN
            NDOUBL=READA(KEYWRD,INDEX(KEYWRD,'C.I.=(')+7)
            NMOS=READA(KEYWRD,INDEX(KEYWRD,'C.I.=(')+5)
         ELSEIF (INDEX(KEYWRD,'C.I.=').NE.0)THEN
            NMOS=READA(KEYWRD,INDEX(KEYWRD,'C.I.=')+5)
         ELSE
            NMOS=NOPEN-NCLOSE
         ENDIF
         IF(NDOUBL.EQ.99)THEN
            J=MAX(MIN((NCLOSE+NOPEN+1)/2-(NMOS-1)/2,NORBS-NMOS+1),1)
         ELSE
            J=NCLOSE-NDOUBL+1
         ENDIF
         NE=MAX(0.D0,(NCLOSE-J+1.D0))*2.D0+
     1     MAX(0.D0,(NOPEN-NCLOSE)*FRACT) + 0.5D0
         NUPP=(NE+1)/2+MSDEL
         NDOWN=NE-NUPP
C
C  NUPP  = NUMBER OF ALPHA ELECTRONS IN ACTIVE SPACE
C  NDOWN = NUMBER OF BETA  ELECTRONS IN ACTIVE SPACE
C
         IF(NUPP*NDOWN.LT.0.OR.NUPP.GT.NMOS.OR.NDOWN.GT.NMOS)THEN
            WRITE(6,'(A)')
     1' SPECIFIED SPIN COMPONENT NOT SPANNED BY ACTIVE SPACE'
            STOP
         ENDIF
      ENDIF
C#      WRITE(6,'(''  NOPEN,NCLOSE,NALPHA,NBETA,FRACT'',4I4,F12.5)')
C#     1 NOPEN, NCLOSE, NLAPHA, NBETA, FRACT
C
C   MAKE SURE ANALYT IS NOT USED WITH ANALYTICAL C.I. DERIVATIVES
C
      HALFE = (NOPEN.GT.NCLOSE.AND.FRACT.NE.2.D0.AND.FRACT.NE.0.D0
     1         .OR. (INDEX(KEYWRD,'C.I.').NE.0))                        IR0494
      SLOW=(INDEX(KEYWRD,'EXCI').NE.0.OR.
     1INDEX(KEYWRD,'ROOT').NE.0.AND.INDEX(KEYWRD,'ROOT=1').EQ.0)
      IF(HALFE)HALFE=(.NOT.SLOW)
      IF(INDEX(KEYWRD,'ANALYT').NE.0.AND.HALFE)THEN
         WRITE(6,*)
         WRITE(6,'(A)')' KEYWORD ''ANALYT'' CANNOT BE USED HERE: ',
     1' ANALYICAL C.I. DERIVATIVES MUST USE FINITE DIFFERENCES',
     2' TO CORRECT, REMOVE KEYWORD ''ANALYT'''
         STOP
      ENDIF
C#      WRITE(6,'(''  NOPEN,NCLOSE,NALPHA,NBETA,FRACT'',4I4,F12.5)')
C#     1 NOPEN, NCLOSE, NLAPHA, NBETA, FRACT
C  END of other C.I.=(n,m) checks
      YY=FLOAT(KHARGE)/(NORBS+1.D-10)
      L=0
      DO 130 I=1,NUMAT
         NI=NAT(I)
         XX=1.D0/(NLAST(I)-NFIRST(I)+1+1.D-10)
         W=CORE(NI)*XX-YY
         IA=NFIRST(I)
         IC=NMIDLE(I)
         IB=NLAST(I)
         DO 110 J=IA,IC
            L=L+1
  110    PSPD(L)=W
         DO 120 J=IC+1,IB
            L=L+1
  120    PSPD(L)=0.D0
  130 CONTINUE
C
C   WRITE OUT THE INTERATOMIC DISTANCES
C
      CALL GMETRY(GEO,COORD)
      RMIN=100.D0
      L=0
      DO 140 I=1,NUMAT
         DO 140 J=1,I
            L=L+1
            RXYZ(L)=SQRT((COORD(1,I)-COORD(1,J))**2+
     1                     (COORD(2,I)-COORD(2,J))**2+
     2                     (COORD(3,I)-COORD(3,J))**2)
            IF(RMIN.GT.RXYZ(L) .AND. I .NE. J .AND.
     1 (NAT(I).LT.103 .OR. NAT(J).LT.103)) THEN
               IMINR=I
               JMINR=J
               RMIN=RXYZ(L)
            ENDIF
  140 CONTINUE
      NNHCO=0
C
C   SET UP MOLECULAR-MECHANICS CORRECTION TO -(C=O)-(NH)- LINKAGE
C   THIS WILL BE USED IF MMOK HAS BEEN SPECIFIED.
C
         ITYPE=1
         IF(INDEX(KEYWRD,'AM1').NE.0)ITYPE=2
         IF(INDEX(KEYWRD,'PM3').NE.0)ITYPE=3
         IF(INDEX(KEYWRD,'MINDO').NE.0)ITYPE=4
C
C   IDENTIFY O=C-N-H SYSTEMS VIA THE INTERATOMIC DISTANCES MATRIX
         DO 190 I=1,NUMAT
            IF(NAT(I).NE.8) GOTO 190
            DO 180 J=1,NUMAT
               IF(NAT(J).NE.6) GOTO 180
               IJ=MAX(I,J)
               JI=I+J-IJ
               IF(RXYZ((IJ*(IJ-1))/2+JI).GT.1.3)GOTO 180
               DO 170 K=1,NUMAT
                  IF(NAT(K).NE.7) GOTO 170
                  JK=MAX(J,K)
                  KJ=J+K-JK
                  IF(RXYZ((JK*(JK-1))/2+KJ).GT.1.6)GOTO 180
                  DO 160 L=1,NUMAT
                     IF(NAT(L).NE.1) GOTO 160
                     KL=MAX(K,L)
                     LK=K+L-KL
                     IF(RXYZ((KL*(KL-1))/2+LK).GT.1.3)GOTO 180
C
C   WE HAVE A H-N-C=O SYSTEM.  THE ATOM NUMBERS ARE L-K-J-I
C   NOW SEARCH OUT ATOM ATTACHED TO NITROGEN, THIS SPECIFIES
C   THE SYSTEM X-N-C=O
C
                     DO 150 M=1,NUMAT
                        IF(M.EQ.K.OR.M.EQ.L.OR.M.EQ.J) GOTO 150
                        MK=MAX(M,K)
                        KM=M+K-MK
                        IF(RXYZ((MK*(MK-1))/2+KM).GT.1.7)GOTO 150
                        NNHCO=NNHCO+1
                        NHCO(1,NNHCO)=I
                        NHCO(2,NNHCO)=J
                        NHCO(3,NNHCO)=K
                        NHCO(4,NNHCO)=M
                        NNHCO=NNHCO+1
                        NHCO(1,NNHCO)=I
                        NHCO(2,NNHCO)=J
                        NHCO(3,NNHCO)=K
                        NHCO(4,NNHCO)=L
                        GOTO 160
  150                CONTINUE
  160             CONTINUE
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
      IF(NNHCO.NE.0)THEN
      IF(INDEX(KEYWRD,'MMOK').NE.0) THEN
      WRITE(6,'(A)')' MOLECULAR MECHANICS CORRECTION APPLIED TO PEPTIDE
     +LINKAGE'
      ELSEIF(INDEX(KEYWRD,'NOMM').NE.0)THEN
      WRITE(6,'(A,I2,A)')' THERE ARE ',NNHCO/2,' PEPTIDE LINKAGES
     1IDENTIFIED IN THIS SYSTEM'
      WRITE(6,'(A)')' IF YOU WANT MM CORRECTION TO THE CONH BARRIER, ADD
     + THE KEY-WORD "MMOK"'
      ELSE
      WRITE(6,'(A)')' THIS SYSTEM CONTAINS -HNCO- GROUPS.'
      WRITE(6,'(A)')' YOU MUST SPECIFY "NOMM" OR "MMOK" REGARDING MOLECU
     +LAR MECHANICS CORRECTION'
      STOP
      ENDIF
      ENDIF
      IF (INDEX(KEYWRD,'NOINTER') .EQ. 0) THEN
         WRITE(6,'(//10X,''  INTERATOMIC DISTANCES'')')
         CALL VECPRT(RXYZ,NUMAT)
      ENDIF
      IF(RMIN.LT.0.8D0.AND.INDEX(KEYWRD,'GEO-OK') .EQ.0) THEN
         WRITE(6,200)IMINR,JMINR,RMIN
  200    FORMAT(//,'   ATOMS',I3,' AND',I3,' ARE SEPARATED BY',F8.4,
     1' ANGSTROMS.',/'   TO CONTINUE CALCULATION SPECIFY "GEO-OK"')
         STOP
      ENDIF
      IF(.NOT. DEBUG) RETURN
      WRITE(6,210)NUMAT,NORBS,NDORBS,NATOMS
  210 FORMAT('   NUMBER OF REAL ATOMS:',I4,/
     1      ,'   NUMBER OF ORBITALS:  ',I4,/
     2      ,'   NUMBER OF D ORBITALS:',I4,/
     3      ,'   TOTAL NO. OF ATOMS:  ',I4)
      WRITE(6,220)(USPD(I),I=1,NORBS)
  220 FORMAT('   ONE-ELECTRON DIAGONAL TERMS',/,10(/,10F8.3))
      WRITE(6,230)(PSPD(I),I=1,NORBS)
  230 FORMAT('   INITIAL P FOR ALL ATOMIC ORBITALS',/,10(/,10F8.3))
      RETURN
  240 WRITE(6,'(//10X,'' MAXIMUM NUMBER OF ATOMIC ORBITALS EXCEEDED'')')
      WRITE(6,'(  10X,'' MAXIMUM ALLOWED ='',I4)')MAXORB
      STOP
      END
      SUBROUTINE READMO
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "mopcvar.par"
C
C
C
C MODULE TO READ IN GEOMETRY FILE, OUTPUT IT TO THE USER,
C AND CHECK THE DATA TO SEE IF IT IS REASONABLE.
C EXIT IF NECESSARY.
C
C
C
C  ON EXIT NATOMS    = NUMBER OF ATOMS PLUS DUMMY ATOMS (IF ANY).
C          KEYWRD    = KEYWORDS TO CONTROL CALCULATION
C          KOMENT    = COMMENT CARD
C          TITLE     = TITLE CARD
C          LABELS    = ARRAY OF ATOMIC LABELS INCLUDING DUMMY ATOMS.
C          GEO       = ARRAY OF INTERNAL COORDINATES.
C          LOPT      = FLAGS FOR OPTIMIZATION OF MOLECULE
C          NA        = ARRAY OF LABELS OF ATOMS, BOND LENGTHS.
C          NB        = ARRAY OF LABELS OF ATOMS, BOND ANGLES.
C          NC        = ARRAY OF LABELS OF ATOMS, DIHEDRAL ANGLES.
C          LATOM     = LABEL OF ATOM OF REACTION COORDINATE.
C          LPARAM    = RC: 1 FOR LENGTH, 2 FOR ANGLE, AND 3 FOR DIHEDRAL
C          REACT(200)= REACTION COORDINATE PARAMETERS
C          LOC(1,I)  = LABEL OF ATOM TO BE OPTIMIZED.
C          LOC(2,I)  = 1 FOR LENGTH, 2 FOR ANGLE, AND 3 FOR DIHEDRAL.
C          NVAR      = NUMBER OF PARAMETERS TO BE OPTIMIZED.
C          XPARAM    = STARTING VALUE OF PARAMETERS TO BE OPTIMIZED.
C
************************************************************************
C *** INPUT THE TRIAL GEOMETRY  \IE.  KGEOM=0\
C   LABEL(I) = THE ATOMIC NUMBER OF ATOM\I\.
C            = 99, THEN THE I-TH ATOM IS A DUMMY ATOM USED ONLY TO
C              SIMPLIFY THE DEFINITION OF THE MOLECULAR GEOMETRY.
C   GEO(1,I) = THE INTERNUCLEAR SEPARATION \IN ANGSTROMS\ BETWEEN ATOMS
C              NA(I) AND (I).
C   GEO(2,I) = THE ANGLE NB(I):NA(I):(I) INPUT IN DEGREES; STORED IN
C              RADIANS.
C   GEO(3,I) = THE ANGLE BETWEEN THE VECTORS NC(I):NB(I) AND NA(I):(I)
C              INPUT IN DEGREES - STORED IN RADIANS.
C  LOPT(J,I) = -1 IF GEO(J,I) IS THE REACTION COORDINATE.
C            = +1 IF GEO(J,I) IS A PARAMETER TO BE OPTIMIZED
C            =  0 OTHERWISE.
C *** NOTE:    MUCH OF THIS DATA IS NOT INCLUDED FOR THE FIRST 3 ATOMS.
C     ATOM1  INPUT LABELS(1) ONLY.
C     ATOM2  INPUT LABELS(2) AND GEO(1,2) SEPARATION BETWEEN ATOMS 1+2
C     ATOM3  INPUT LABELS(3), GEO(1,3)    SEPARATION BETWEEN ATOMS 2+3
C              AND GEO(2,3)              ANGLE ATOM1 : ATOM2 : ATOM3
C
************************************************************************
C
      DIMENSION LOPT(3,NUMATM)
      CHARACTER*80 KEYWRD,KOMENT,TITLE,LINE, BANNER
      CHARACTER KEYS(80)*1, SPACE*1, SPACE2*2, CH*1, CH2*2
      CHARACTER ELEMNT*2, IDATE*24
      COMMON /KEYWRD/ KEYWRD
      COMMON /TITLES/ KOMENT,TITLE
      COMMON /GEOVAR/ XPARAM(MAXPAR), NVAR, LOC(2,MAXPAR)               IR0394
C    changed common path for portability  (IR)
      COMMON /PATHI / LATOM,LPARAM
      COMMON /PATHR / REACT(200)
      COMMON /MESH  / LATOM1, LPARA1, LATOM2, LPARA2
      COMMON /ELEMTS/ ELEMNT(107)
      COMMON /ISTOPE/ AMS(107)
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /GEOSYM/ NDEP, LOCPAR(MAXPAR), IDEPFN(MAXPAR),
     1                      LOCDEP(MAXPAR)
      COMMON /IOCM/ IREAD
      LOGICAL INT
      DIMENSION COORD(3,NUMATM),VALUE(40)
      EQUIVALENCE (KEYS(1),KEYWRD)
         SAVE                                                           GL0892
      DATA SPACE, SPACE2/' ','  '/
C
c      WRITE(6,*) "I GET TO READMO"
      READ(IREAD,'(A)')KEYWRD,KOMENT,TITLE
      ILOWA = ICHAR('a')
      ILOWZ = ICHAR('z')
      ICAPA = ICHAR('A')
************************************************************************
      DO 10 I=1,80
         ILINE=ICHAR(KEYWRD(I:I))
         IF(ILINE.GE.ILOWA.AND.ILINE.LE.ILOWZ) THEN
            KEYWRD(I:I)=CHAR(ILINE+ICAPA-ILOWA)
         ENDIF
   10 CONTINUE
************************************************************************
      IF(INDEX(KEYWRD,'ECHO').NE.0)THEN
         REWIND IREAD
         DO 30 I=1,1000
            READ(IREAD,'(A)',END=40)KEYWRD
            DO 20 J=80,2,-1
   20       IF(KEYWRD(J:J).NE.' ')GOTO 30
            J=1
   30    WRITE(6,'(1X,A)')KEYWRD(1:J)
      ENDIF
   40 REWIND IREAD
      IF(INDEX(KEYWRD,'ECHO').NE.0)WRITE(6,'(''1'')')
      READ(IREAD,'(A)')KEYWRD,KOMENT,TITLE
************************************************************************
      DO 50 I=1,80
         ILINE=ICHAR(KEYWRD(I:I))
         IF(ILINE.GE.ILOWA.AND.ILINE.LE.ILOWZ) THEN
            KEYWRD(I:I)=CHAR(ILINE+ICAPA-ILOWA)
         ENDIF
   50 CONTINUE
************************************************************************
      IF(KEYWRD(1:1) .NE. SPACE) THEN
         CH=KEYWRD(1:1)
         KEYWRD(1:1)=SPACE
         DO 60 I=2,78
            CH2=KEYWRD(I:I)
            KEYWRD(I:I)=CH
            CH=CH2
            IF(KEYWRD(I+1:I+2) .EQ. SPACE2) THEN
               KEYWRD(I+1:I+1)=CH
               GOTO 70
            ENDIF
   60    CONTINUE
         IF(KEYWRD(80:80) .NE. SPACE)
     1WRITE(6,'(//10X,''WARNING-THE LAST KEYWORD MAY BE CORRUPTED'')')
         KEYWRD(80:80)=CH
   70    CONTINUE
      ENDIF
      IF(KOMENT(1:1) .NE. SPACE) THEN
         CH=KOMENT(1:1)
         KOMENT(1:1)=SPACE
         DO 80 I=2,78
            CH2=KOMENT(I:I)
            KOMENT(I:I)=CH
            CH=CH2
            IF(KOMENT(I+1:I+2) .EQ. SPACE2) THEN
               KOMENT(I+1:I+1)=CH
               GOTO 90
            ENDIF
   80    CONTINUE
         IF(KOMENT(80:80) .NE. SPACE)
     1WRITE(6,'(//1X,A)')'WARNING-THE LAST WORD IN FIRST COMMENT LINE'//
     2' MAY BE CORRUPTED'
         KOMENT(80:80)=CH
   90    CONTINUE
      ENDIF
      IF(TITLE(1:1) .NE. SPACE) THEN
         CH=TITLE(1:1)
         TITLE(1:1)=SPACE
         DO 100 I=2,78
            CH2=TITLE(I:I)
            TITLE(I:I)=CH
            CH=CH2
            IF(TITLE(I+1:I+2) .EQ. SPACE2) THEN
               TITLE(I+1:I+1)=CH
               GOTO 110
            ENDIF
  100    CONTINUE
         IF(TITLE(80:80) .NE. SPACE)
     1WRITE(6,'(//1X,A)')'WARNING-THE LAST WORD IN SECOND COMMENT'//
     2' LINE MAY BE CORRUPTED'
         TITLE(80:80)=CH
  110    CONTINUE
      ENDIF
      CALL GETGEO(IREAD,LABELS,GEO,LOPT,NA,NB,NC,AMS,NATOMS,INT)
C
C
C OUTPUT FILE TO UNIT 6
C
C    WRITE HEADER
      IDATE=' '
      CALL DATE(IDATE)
C          WRITE(6,'(//)')
      WRITE(6,'(1X,15(''*****''),''****'')')
C
C     CHANGE THE FOLLOWING LINE TO SUIT LOCAL ENVIRONMENT, IF DESIRED
C
      BANNER=' **  Welcome to MOPAC 505 - A version of MOPAC '//
     1'for Direct Molecular Dynamics  **'
           WRITE(6,'(A)')BANNER
C
C    THE BANNER DOES NOT APPEAR ANYWHERE ELSE.
C
      WRITE(6,'(1X,79(''*''))')
      LINE='   MNDO'
      IF(INDEX(KEYWRD,'MINDO') .NE. 0) LINE='MINDO/3'
      IF(INDEX(KEYWRD,'AM1')   .NE. 0) LINE='    AM1'
      IF(INDEX(KEYWRD,'PM3')   .NE. 0) LINE='    PM3'
      WRITE(6,'(/29X,A,'' CALCULATION RESULTS'',28X,///1X,15(''*****'')
     1,''****'' )')LINE(:7)
      WRITE(6,'('' *'',10X,''MOPAC:  VERSION '',F5.2,
     115X,''CALC''''D. '',A24)') VERSON, IDATE
C
C CHECK DATA
C
      DO 120 I=1,NATOMS
         IF (LABELS(I) .LE. 0 ) THEN
            WRITE(6,'('' ATOMIC NUMBER OF '',I3,'' ?'')') LABELS(I)
            STOP
         ENDIF
         IF (  NA(I).GE.I.OR. NB(I).GE.I.OR. NC(I).GE.I
     1  .OR. (NA(I).EQ.NB(I))   .AND. I.GT.1
     2  .OR. (NA(I).EQ.NC(I).OR.NB(I).EQ.NC(I))  .AND. I.GT.2
     3  .OR.  NA(I)*NB(I)*NC(I).EQ.0  .AND. I.GT.3) THEN
            WRITE(6,'('' ATOM NUMBER '',I3,'' IS ILLDEFINED'')') I
            STOP
         ENDIF
  120 CONTINUE
C
C WRITE KEYWORDS BACK TO USER AS FEEDBACK
      CALL WRTKEY(KEYWRD)
      WRITE(6,'(1X,14(''*****''),''*'',I3.3,''BY'',I3.3)')MAXHEV,MAXLIT
C
C CONVERT ANGLES TO RADIANS
      DO 130 I=1,NATOMS
         DO 130 J=2,3
            GEO(J,I) = GEO(J,I) * 2.D0*ASIN(1.D0)/180.D0
  130 CONTINUE
C
C FILL IN GEO MATRIX IF NEEDED
      NDEP=0
      IF( INDEX(KEYWRD,'SYM') .NE. 0) CALL GETSYM
      IF(NDEP.NE.0) CALL SYMTRY
C
C INITIALIZE FLAGS FOR OPTIMIZE AND PATH
      IFLAG = 0
      NVAR  = 0
      LATOM = 0
      NUMAT=0
      DO 160 I=1,NATOMS
         IF(LABELS(I).NE.99.AND.LABELS(I).NE.107)NUMAT=NUMAT+1
         DO 160 J=1,3
            IF (LOPT(J,I) ) 140, 160, 150
C    FLAG FOR PATH
  140       CONVRT=1.D0
            IF ( IFLAG .NE. 0 ) THEN
               IF(INDEX(KEYWRD,'STEP1').NE.0)THEN
                  LPARA1=LPARAM
                  LATOM1=LATOM
                  LPARA2=J
                  LATOM2=I
                  LATOM=0
                  IFLAG=0
                  GOTO 160
               ELSE
                  WRITE(6,'('' ONLY ONE REACTION COORDINATE PERMITTED'')
     1')
                  STOP
               ENDIF
            ENDIF
            LATOM  = I
            LPARAM = J
            IF(J.GT.1) CONVRT=0.01745329252D00
            REACT(1)  = GEO(J,I)
            IREACT=1
            IFLAG = 1
            GO TO 160
C    FLAG FOR OPTIMIZE
  150       NVAR = NVAR + 1
            LOC(1,NVAR) = I
            LOC(2,NVAR) = J
            XPARAM(NVAR)   = GEO(J,I)
  160 CONTINUE
C READ IN PATH VALUES
      IF(IFLAG.EQ.0) GO TO 200
      IF(INDEX(KEYWRD,'NLLSQ').NE.0)THEN
         WRITE(6,'(A)')' NLLSQ USED WITH REACTION PATH; '//
     1'THIS OPTION IS NOT ALLOWED'
         STOP
      ENDIF
      IF(INDEX(KEYWRD,'SIGMA').NE.0)THEN
         WRITE(6,'(A)')' SIGMA USED WITH REACTION PATH; '//
     1'THIS OPTION IS NOT ALLOWED'
         STOP
      ENDIF
  170 READ(IREAD,'(A)',END=190) LINE
      CALL NUCHAR(LINE,VALUE,NREACT)
      DO 180 I=1,NREACT
         IJ=IREACT+I
         IF(IJ.GT.200)THEN
            WRITE(6,'(///,''    ONLY TWO HUNDRED POINTS ALLOWED IN REACT
     1ION'','' COORDINATE'')')
            STOP
         ENDIF
         REACT(IJ)=VALUE(I)*CONVRT
         IF(ABS(REACT(IJ)-REACT(IJ-1)).LT.1.D-5)THEN
            DUM1 = REACT(IJ)/CONVRT
            DUM2 = REACT(IJ-1)/CONVRT
            WRITE(6,'(///,'' TWO ADJACENT POINTS ARE IDENTICAL:  '',
     1 F7.3,2X,F7.3,/,'' THIS IS NOT ALLOWED IN A PATH CALCULATION'')')
     2 DUM1,DUM2
            STOP
         ENDIF
  180 CONTINUE
      IREACT=IREACT+NREACT
      GO TO 170
  190 CONTINUE
      DEGREE=1.D0
      IF(LPARAM.GT.1)DEGREE=90.D0/ASIN(1.D0)
      IF(IREACT.LE.1) THEN
         WRITE(6,'(//10X,'' NO POINTS SUPPLIED FOR REACTION PATH'')')
         WRITE(6,'(//10X,'' GEOMETRY AS READ IN IS AS FOLLOWS'')')
         CALL GEOUT
         STOP
      ELSE
         WRITE(6,'(//10X,'' POINTS ON REACTION COORDINATE'')')
         WRITE(6,'(10X,8F8.2)')(REACT(I)*DEGREE,I=1,IREACT)
      ENDIF
      IEND=IREACT+1
      REACT(IEND)=-1.D12
C
C OUTPUT GEOMETRY AS FEEDBACK
C
  200 WRITE(6,'(A)')KEYWRD,KOMENT,TITLE
      CALL GEOUT
      CALL GMETRY(GEO,COORD)
      IF (INDEX(KEYWRD,'NOXYZ') .EQ. 0) THEN
         WRITE(6,'(//10X,''CARTESIAN COORDINATES '',/)')
         WRITE(6,'(4X,''NO.'',7X,''ATOM'',9X,''X'',
     1  9X,''Y'',9X,''Z'',/)')
         L=0
         DO 210 I=1,NATOMS
            IF(LABELS(I) .EQ. 99.OR.LABELS(I).EQ.107) GOTO 210
            L=L+1
            WRITE(6,'(I6,8X,A2,4X,3F10.4)')
     1  L,ELEMNT(LABELS(I)),(COORD(J,L),J=1,3)
  210    CONTINUE
      ENDIF
      IF(   INDEX(KEYWRD,' XYZ') .NE. 0 )THEN
         IF( NVAR .NE. 0 .AND.
     1 INT.AND.(NDEP .NE. 0 .OR.  NVAR.LT.3*NUMAT-6)) THEN
            IF(NDEP.NE.0)
     1WRITE(6,'(//10X,'' INTERNAL COORDINATES READ IN, AND SYMMETRY''
     2,/10X,'' SPECIFIED, BUT CALCULATION TO BE RUN IN CARTESIAN ''
     3,''COORDINATES'')')
            IF(NVAR.LT.3*NUMAT-6)
     1WRITE(6,'(//10X,'' INTERNAL COORDINATES READ IN, AND'',
     2'' CALCULATION '',/10X,''TO BE RUN IN CARTESIAN COORDINATES, '',
     3/10X,''BUT NOT ALL COORDINATES MARKED FOR OPTIMISATION'')')
            WRITE(6,'(//10X,'' THIS INVOLVES A LOGICALLLY ABSURD CHOICE'
     1',/10X,'' SO THE CALCULATION IS TERMINATED AT THIS POINT'')')
            STOP
         ENDIF
         SUMX=0.D0
         SUMY=0.D0
         SUMZ=0.D0
         DO 220 J=1,NUMAT
            SUMX=SUMX+COORD(1,J)
            SUMY=SUMY+COORD(2,J)
  220    SUMZ=SUMZ+COORD(3,J)
         SUMX=SUMX/NUMAT
         SUMY=SUMY/NUMAT
         SUMZ=SUMZ/NUMAT
         DO 230 J=1,NUMAT
            GEO(1,J)=COORD(1,J)-SUMX
            GEO(2,J)=COORD(2,J)-SUMY
  230    GEO(3,J)=COORD(3,J)-SUMZ
         NA(1)=99
         J=0
         MVAR=1
         NVAR=1
         DO 260 I=1,NATOMS
            IF(LABELS(I).NE.99)THEN
               J=J+1
  240          IF(LOC(1,MVAR) .EQ. I) THEN
                  XPARAM(NVAR)=GEO(LOC(2,NVAR),J)
                  LOC(1,NVAR)=J
                  MVAR=MVAR+1
                  NVAR=NVAR+1
                  GOTO 240
               ENDIF
               LABELS(J)=LABELS(I)
            ELSE
  250          IF(LOC(1,MVAR) .EQ. I)THEN
                  MVAR=MVAR+1
                  GOTO 250
               ENDIF
            ENDIF
  260    CONTINUE
         NVAR=NVAR-1
         NATOMS=NUMAT
      ELSE
         IF(NVAR.EQ.0) RETURN
         IF( .NOT. INT.AND.(NDEP .NE. 0 .OR.  NVAR.LT.3*NUMAT-6)) THEN
            IF(NDEP.NE.0)
     1WRITE(6,'(//10X,'' CARTESIAN COORDINATES READ IN, AND SYMMETRY''
     2,/10X,'' SPECIFIED, BUT CALCULATION TO BE RUN IN INTERNAL ''
     3,''COORDINATES'')')
            IF(NVAR.LT.3*NUMAT-6)
     1WRITE(6,'(//10X,'' CARTESIAN COORDINATES READ IN, AND'',
     2'' CALCULATION '',/10X,''TO BE RUN IN INTERNAL COORDINATES, '',
     3/10X,''BUT NOT ALL COORDINATES MARKED FOR OPTIMISATION'')')
            WRITE(6,'(//10X,''MOPAC, BY DEFAULT, USES INTERNAL COORDINAT
     1ES'',/10X,''TO SPECIFY CARTESIAN COORDINATES USE KEY-WORD :XYZ:'')
     2')
            WRITE(6,'(10X,''YOUR CURRENT CHOICE OF KEY-WORDS INVOLVES A
     1'',''LOGICALLLY '',
     2/10X,''ABSURD CHOICE SO THE CALCULATION IS TERMINATED AT THIS ''
     3,''POINT'')')
            STOP
         ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE SOLROT (NI,NJ,XI,XJ,WJ,WK,KR,E1B,E2A,ENUC,CUTOFF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XI(3), XJ(3), WJ(*), WK(*), E1B(10), E2A(10)
************************************************************************
*
*   SOLROT FORMS THE TWO-ELECTRON TWO-ATOM J AND K INTEGRAL STRINGS.
*          ON EXIT WJ = "J"-TYPE INTEGRALS
*                  WK = "K"-TYPE INTEGRALS
*
*      FOR MOLECULES, WJ = WK.
************************************************************************
      COMMON /EULER / TVEC(3,3), ID
      COMMON /UCELL / L1L,L2L,L3L,L1U,L2U,L3U
      COMMON /NUMCAL/ NUMCAL
      DIMENSION WSUM(100), WBITS(100), LIMS(3,2), XJUC(3), E1BITS(10),
     1E2BITS(10), WMAX(100)
      EQUIVALENCE (L1L,LIMS(1,1))
         SAVE                                                           GL0892
      DATA ICALCN /0/
      IF (ICALCN .NE. NUMCAL) THEN
         ICALCN = NUMCAL
         DO 10 I=1,ID
            LIMS(I,1)=-1
   10    LIMS(I,2)= 1
         DO 20 I=ID+1,3
            LIMS(I,1)=0
   20    LIMS(I,2)=0
      ENDIF
      ONE=1.D0
      IF(XI(1).EQ.XJ(1) .AND. XI(2).EQ.XJ(2) .AND. XI(3).EQ. XJ(3))
     1ONE=0.5D0
      DO 30 I=1,100
         WMAX(I)=0.D0
         WSUM(I)=0.D0
   30 WBITS(I)=0.D0
      NEQUAL=1
      DO 40 I=1,10
         E1B(I)=0.D0
   40 E2A(I)=0.D0
      ENUC=0.D0
      DO 90 I=L1L,L1U
         DO 90 J=L2L,L2U
            DO 90 K=L3L,L3U
               DO 50 L=1,3
   50          XJUC(L)=XJ(L)+TVEC(L,1)*I+TVEC(L,2)*J+TVEC(L,3)*K
               KB=1
               CALL ROTATE (NI,NJ,XI,XJUC,WBITS,KB,E1BITS,E2BITS,
     1ENUBIT,CUTOFF)
               KB=KB-1
               DO 60 II=1,KB
   60          WSUM(II)=WSUM(II)+WBITS(II)
               IF(WMAX(1).LT.WBITS(1))THEN
                  DO 70 II=1,KB
   70             WMAX(II)=WBITS(II)
               ENDIF
               DO 80 II=1,10
                  E1B(II)=E1B(II)+E1BITS(II)
   80          E2A(II)=E2A(II)+E2BITS(II)
               ENUC=ENUC+ENUBIT*ONE
   90 CONTINUE
      IF(ONE.LT.0.9D0) THEN
         DO 100 I=1,KB
  100    WMAX(I)=0.D0
      ENDIF
      DO 110 I=1,KB
         WK(I)=WMAX(I)
  110 WJ(I)=WSUM(I)
      KR=KB+KR
      RETURN
      END
      DOUBLE PRECISION FUNCTION SS(NA,NB,LA1,LB1,M1,UA,UB,R1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /NUMCAL/ NUMCAL
      DIMENSION FA(0:13),AFF(0:2,0:2,0:2),AF(0:19),BF(0:19),
     1BI(0:12,0:12)
         SAVE                                                           GL0892
      DATA AFF/27*0. D0/
      DATA FA/1.D0,1.D0,2.D0,6.D0,24.D0,120.D0,720.D0,5040.D0,40320.D0,
     1362880.D0,3628800.D0,39916800.D0,479001600.D0,6227020800.D0/
      DATA ICALCN /0/
      M=M1-1
      LB=LB1-1
      LA=LA1-1
      R=R1/0.529167D0
      IF (ICALCN .NE. NUMCAL) THEN
        ICALCN = NUMCAL
C
C           INITIALISE SOME CONSTANTS
C
C                  BINOMIALS
C
         DO 10 I=0,12
            BI(I,0)=1.D0
            BI(I,I)=1.D0
   10    CONTINUE
         DO 20 I=0,11
            I1=I-1
            DO 20 J=0,I1
               BI(I+1,J+1)=BI(I,J+1)+BI(I,J)
   20    CONTINUE
         AFF(0,0,0)=1.D0
         AFF(1,0,0)=1.D0
         AFF(1,1,0)=1.D0
         AFF(2,0,0)=1.5D0
         AFF(2,1,0)=1.73205D0
         AFF(2,2,0)=1.224745D0
         AFF(2,0,2)=-0.5D0
      ENDIF
      P=(UA+UB)*R*0.5D0
      B=(UA-UB)*R*0.5D0
      EX=EXP(B)
      QUO=1/P
      AF(0)=QUO*EXP(-P)
      DO 30 N=1,19
         AF(N)=N*QUO*AF(N-1)+AF(0)
   30 CONTINUE
      CALL BFN(B,BF)
      SUM=0.D0
      LAM1=LA-M
      LBM1=LB-M
C
C          START OF OVERLAP CALCULATION PROPER
C
      DO 50 I=0,LAM1,2
         IA=NA+I-LA
         IC=LA-I-M
         DO 50 J=0,LBM1,2
            IB=NB+J-LB
            ID=LB-J-M
            SUM1=0.D0
            IAB=IA+IB
            DO 40 K1=0,IA
               DO 40 K2=0,IB
                  DO 40 K3=0,IC
                     DO 40 K4=0,ID
                        DO 40 K5=0,M
                           IAF=IAB-K1-K2+K3+K4+2*K5
                           DO 40 K6=0,M
                              IBF=K1+K2+K3+K4+2*K6
                              JX=(-1)**(M+K2+K4+K5+K6)
                              SUM1=SUM1+BI(ID,K4)*
     1BI(M,K5)*BI(IC,K3)*BI(IB,K2)*BI(IA,K1)*
     2BI(M,K6)*JX*AF(IAF)*BF(IBF)
   40       CONTINUE
            SUM=SUM+SUM1*AFF(LA,M,I)*AFF(LB,M,J)
   50 CONTINUE
      SS=SUM*R**(NA+NB+1)*UA**NA*UB**NB/(2.D0**(M+1))*
     1SQRT(UA*UB/(FA(NA+NA)*FA(NB+NB))*((LA+LA+1)*(LB+LB+1)))
      RETURN
      END
      SUBROUTINE WRTKEY(KEYWRD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      CHARACTER*80 KEYWRD, ALLKEY
***********************************************************************
*
*  WRTKEY CHECKS ALL KEY-WORDS AND PRINTS THOSE IT RECOGNIZES.  IF IT
*  FINDS A WORD IT DOES NOT RECOGNIZE THE PROGRAM WILL BE STOPPED.
*
***********************************************************************
      LOGICAL UHF, TRIP, BIRAD, EXCI, CI, WORD
      LOGICAL AM1, MNDO, MINDO3, PM3
      CHARACTER SPACE*1, DOT*1, ZERO*1, NINE*1, CH*1, CHRONO*7
         SAVE                                                           GL0892
      DATA SPACE,DOT,ZERO,NINE /' ','.','0','9'/
      DATA AM1, MNDO, MINDO3, PM3 /4*.FALSE./
      ALLKEY=KEYWRD
      IF(WORD(ALLKEY,'AUTHOR')) THEN
         WRITE(6,'(10X,'' MOPAC - A GENERAL MOLECULAR ORBITAL PACKAGE'',
     1/         ,10X,''   ORIGINAL VERSION WRITTEN IN 1983'')')
         WRITE(6,'(10X,''     BY JAMES J. P. STEWART AT THE'',/
     1         ,10X,''     UNIVERSITY OF TEXAS AT AUSTIN'',/
     2         ,10X,''          AUSTIN, TEXAS, 78712'')')
      ENDIF
      IF (WORD(ALLKEY,'VECT') ) WRITE(6,70)
      IF (WORD(ALLKEY,' EXTE') ) THEN
         I=INDEX(KEYWRD,' EXTE')
         J=INDEX(KEYWRD(I:),'=')+I
         I=INDEX(KEYWRD(J:),' ')+J-1
         WRITE(6,80)KEYWRD(J:I)
      ENDIF
      MAXGEO=0
      IF (WORD(ALLKEY,' DENS') ) WRITE(6,90)
      IF (WORD(ALLKEY,'SPIN') ) WRITE(6,100)
      IF (WORD(ALLKEY,' DEPVAR') )
     1WRITE(6,110)READA(KEYWRD,INDEX(KEYWRD,'DEPVAR'))
      IF (WORD(ALLKEY,' DEP ') )WRITE(6,120)
      IF (WORD(ALLKEY,'VELO') )WRITE(6,130)
      IF (WORD(ALLKEY,'TIMES') )WRITE(6,140)
      IF (WORD(ALLKEY,'PARASOK') ) WRITE(6,150)
      IF (WORD(ALLKEY,'FLEPO') ) WRITE(6,160)
      IF (WORD(ALLKEY,'BONDS') ) WRITE(6,170)
      IF (WORD(ALLKEY,'GEO-OK') ) WRITE(6,180)
      IF (WORD(ALLKEY,'FOCK') ) WRITE(6,190)
      IF (WORD(ALLKEY,'LARGE') ) WRITE(6,200)
      IF (WORD(ALLKEY,' S1978') ) WRITE(6,210)
      IF (WORD(ALLKEY,' SI1978') ) WRITE(6,220)
      IF (WORD(ALLKEY,' GRAP') ) WRITE(6,230)
      IF (WORD(ALLKEY,'1ELEC') ) WRITE(6,240)
      IF (WORD(ALLKEY,' NOMM') ) WRITE(6,260)
      IF (WORD(ALLKEY,' MMOK') ) WRITE(6,265)
      IF (WORD(ALLKEY,'INTERP') ) WRITE(6,270)
      IF (WORD(ALLKEY,' ESR') ) WRITE(6,250)
      IF (WORD(ALLKEY,'DFP') ) WRITE(6,280)
      IF (WORD(ALLKEY,'ANALYT') ) WRITE(6,290)
      IF (WORD(ALLKEY,' MECI') ) WRITE(6,300)
      IF (WORD(ALLKEY,'LOCAL') ) WRITE(6,340)
      IF (WORD(ALLKEY,'MULLIK') ) WRITE(6,350)
      IF (WORD(ALLKEY,' XYZ') ) WRITE(6,360)
      IF (WORD(ALLKEY,' PI') ) WRITE(6,370)
      IF (WORD(ALLKEY,'ECHO') ) WRITE(6,380)
      IF (WORD(ALLKEY, 'SING') ) WRITE(6,540)
      IF (WORD(ALLKEY, 'DOUB') ) WRITE(6,550)
      IF (WORD(ALLKEY, 'QUAR') ) WRITE(6,570)
      IF (WORD(ALLKEY, 'QUIN') ) WRITE(6,580)
      IF (WORD(ALLKEY, 'SEXT') ) WRITE(6,590)
      IF (WORD(ALLKEY,'H-PRIO') ) WRITE(6,390)
      IF (WORD(ALLKEY,'X-PRIO') ) WRITE(6,400)
      IF (WORD(ALLKEY,'T-PRIO') ) WRITE(6,410)
      IF (WORD(ALLKEY,'POWSQ') ) WRITE(6,420)
      IF (WORD(ALLKEY,'POLAR') ) WRITE(6,430)
      IF (WORD(ALLKEY,'DEBUG ') ) WRITE(6,440)
      IF (WORD(ALLKEY,'RESTART') ) WRITE(6,450)
      IF (WORD(ALLKEY,'IRC=') ) THEN
         MAXGEO=1
         WRITE(6,460)NINT(READA(KEYWRD,INDEX(KEYWRD,'IRC=')))
      ELSEIF (WORD(ALLKEY,'IRC') ) THEN
         MAXGEO=1
         WRITE(6,470)
      ENDIF
      IF (WORD(ALLKEY,'CHARGE') )
     1 WRITE(6,480)NINT(READA(KEYWRD,INDEX(KEYWRD,'CHARGE')))
      IF (WORD(ALLKEY,'GRAD') ) WRITE(6,490)
      UHF=(WORD(ALLKEY,'UHF') )
      IF(UHF)WRITE(6,500)
      BIRAD=(WORD(ALLKEY,'BIRAD') )
      IF(BIRAD)WRITE(6,520)
      EXCI=(WORD(ALLKEY,'EXCITED') )
      IF(EXCI) WRITE(6,530)
      TRIP=(WORD(ALLKEY,'TRIPLET') )
      IF(TRIP)WRITE(6,560)
      IF (WORD(ALLKEY,'SYM') ) WRITE(6,600)
      IF(WORD(ALLKEY,'OPEN('))THEN
         I=INDEX(KEYWRD,'OPEN(')
         IELEC=READA(KEYWRD,I)
         ILEVEL=READA(KEYWRD,I+7)
         WRITE(6,620)IELEC,ILEVEL
      ENDIF
      IF(WORD(ALLKEY,'MICROS'))
     1WRITE(6,610)INT(READA(KEYWRD,INDEX(KEYWRD,'MICROS')))
      IF(WORD(ALLKEY,'DRC='))THEN
         MAXGEO=1
         WRITE(6,320)READA(KEYWRD,INDEX(KEYWRD,'DRC='))
      ELSEIF (WORD(ALLKEY,' DRC') ) THEN
         MAXGEO=1
         WRITE(6,310)
      ENDIF
      IF(WORD(ALLKEY,'KINE'))
     1WRITE(6,330)READA(KEYWRD,INDEX(KEYWRD,'KINE'))
      CHRONO='SECONDS'
      IF(WORD(ALLKEY,' T=')) THEN
         I=INDEX(KEYWRD,' T=')
         TIME=READA(KEYWRD,I)
         DO 10 J=I+3,80
            IF( KEYWRD(J+1:J+1).EQ.' ') THEN
               CH=KEYWRD(J:J)
               IF( CH .EQ. 'M') CHRONO='MINUTES'
               IF( CH .EQ. 'H') CHRONO='HOURS'
               IF( CH .EQ. 'D') CHRONO='DAYS'
               GOTO 20
            ENDIF
   10    CONTINUE
   20    CONTINUE
         IF(TIME.LT.99999.9D0)THEN
            WRITE(6,630)TIME,CHRONO
         ELSE
            WRITE(6,640)TIME,CHRONO
         ENDIF
      ELSE
         WRITE(6,630)MAXTIM*1.D0,CHRONO
      ENDIF
               CHRONO='SECONDS'
      IF(WORD(ALLKEY,' DUMP')) THEN
         TDUMP=READA(KEYWRD,INDEX(KEYWRD,' DUMP'))
         DO 30 J=I+6,80
            IF( KEYWRD(J+1:J+1).EQ.' ') THEN
               CH=KEYWRD(J:J)
               IF( CH .EQ. 'M') CHRONO='MINUTES'
               IF( CH .EQ. 'H') CHRONO='HOURS'
               IF( CH .EQ. 'D') CHRONO='DAYS'
               GOTO 40
            ENDIF
   30    CONTINUE
   40    CONTINUE
         IF(TDUMP.LT.99999.9D0)THEN
            WRITE(6,650)TDUMP,CHRONO
         ELSE
            WRITE(6,660)TDUMP,CHRONO
         ENDIF
      ELSE
         WRITE(6,650)MAXDMP*1.D0,CHRONO
      ENDIF
      IF (WORD(ALLKEY,'1SCF') ) THEN
         WRITE(6,670)
         IF(INDEX(KEYWRD,'RESTART').EQ.0)MAXGEO=MAXGEO+1
      ENDIF
      CI=WORD(ALLKEY,'C.I.')
C C.I.=(n,m) keyword                                                    IR0494
      IF (CI) THEN
         J=INDEX(KEYWRD,' C.I.=(')
         IF(J.NE.0)THEN
            WRITE(6,681)INT(READA(KEYWRD,INDEX(KEYWRD,'C.I.=(')+7)),
     1      INT(READA(KEYWRD,INDEX(KEYWRD,'C.I.=(')+5))
         ELSE
            WRITE(6,680)INT(READA(KEYWRD,INDEX(KEYWRD,'C.I.')+5))
         ENDIF
      ENDIF
C End of C.I.=(n,m)
      IF (WORD(ALLKEY,'DFORCE') ) THEN
         WRITE(6,690)
c        MAXGEO=MAXGEO+1
      ELSEIF (WORD(ALLKEY,'FORCE') ) THEN
         WRITE(6,700)
c        MAXGEO=MAXGEO+1
      ENDIF
      METHOD=0
      IF (WORD(ALLKEY,'MINDO') ) THEN
         WRITE(6,710)
         MINDO3=.TRUE.
         METHOD=1
      ENDIF
      IF (WORD(ALLKEY,'AM1') ) THEN
         WRITE(6,720)
         AM1=.TRUE.
         METHOD=METHOD+1
      ENDIF
      IF (WORD(ALLKEY,'PM3') ) THEN
         WRITE(6,730)
         PM3=.TRUE.
         METHOD=METHOD+1
      ENDIF
      IF (WORD(ALLKEY,'MNDO') ) THEN
         MNDO=.TRUE.
         METHOD=METHOD+1
      ENDIF
      IF (WORD(ALLKEY,'PREC') ) WRITE(6,740)
      IF (WORD(ALLKEY,'NOINTER') ) WRITE(6,750)
      IF (WORD(ALLKEY,'ISOTOPE') ) WRITE(6,760)
      IF (WORD(ALLKEY,'DENOUT') ) WRITE(6,770)
      IF (WORD(ALLKEY,'SHIFT') ) WRITE(6,780)
     1 READA(KEYWRD,INDEX(KEYWRD,'SHIFT'))
      IF (WORD(ALLKEY,'OLDENS') ) WRITE(6,790)
      IF (WORD(ALLKEY,'SCFCRT') ) WRITE(6,800)
     1 READA(KEYWRD,INDEX(KEYWRD,'SCFCRT'))
      IF (WORD(ALLKEY,'ENPART') ) WRITE(6,810)
      IF (WORD(ALLKEY,'NOXYZ') ) WRITE(6,820)
      IF (WORD(ALLKEY,'SIGMA') ) THEN
         WRITE(6,830)
         MAXGEO=MAXGEO+1
      ENDIF
      IF (WORD(ALLKEY,'NLLSQ') ) THEN
         WRITE(6,840)
         MAXGEO=MAXGEO+1
      ENDIF
      IF (WORD(ALLKEY,'ROOT') ) WRITE(6,850)
     1 NINT(READA(KEYWRD,INDEX(KEYWRD,'ROOT')))
      IF (WORD(ALLKEY,'TRANS=') ) THEN
         WRITE(6,870)NINT(READA(KEYWRD,INDEX(KEYWRD,'TRANS=')))
      ELSEIF (WORD(ALLKEY,'TRANS') ) THEN
         WRITE(6,860)
      ENDIF
      IF (WORD(ALLKEY,'SADDLE') ) THEN
         WRITE(6,880)
         MAXGEO=MAXGEO+1
      ENDIF
C
C	Start of EF section						IR0494
C
      IF (WORD(ALLKEY,'EF') ) THEN
         WRITE(6,1190)
         MAXGEO=MAXGEO+1
      ENDIF
      IF (WORD(ALLKEY,'TS') ) THEN
         WRITE(6,1200)
         MAXGEO=MAXGEO+1
      ENDIF
      IF (WORD(ALLKEY,' IUPD')) THEN
         I=NINT(READA(KEYWRD,INDEX(KEYWRD,' IUPD=')))
         IF (I.EQ.0) WRITE(6,1210)
         IF (I.EQ.1) WRITE(6,1220)
         IF (I.EQ.2) WRITE(6,1230)
      ENDIF
      IF (WORD(ALLKEY,' NONR') )WRITE(6,1240)
      IF (WORD(ALLKEY,' HESS')) THEN
         I=NINT(READA(KEYWRD,INDEX(KEYWRD,' HESS=')))
         IF (I.EQ.0) WRITE(6,1250)
         IF (I.EQ.1) WRITE(6,1260)
         IF (I.EQ.2) WRITE(6,1270)
         IF (I.EQ.3) WRITE(6,1280)
      ENDIF
      IF (WORD(ALLKEY,' MODE')) WRITE(6,1290)
     1 NINT(READA(KEYWRD,INDEX(KEYWRD,'MODE=')))
      IF (WORD(ALLKEY,' RSCAL') )WRITE(6,1300)
      IF (WORD(ALLKEY,' RECALC')) WRITE(6,1310)
     1 NINT(READA(KEYWRD,INDEX(KEYWRD,'RECALC')))
      IF (WORD(ALLKEY,' DMAX')) WRITE(6,1320)
     1 READA(KEYWRD,INDEX(KEYWRD,'DMAX='))
      IF (WORD(ALLKEY,' OMIN')) WRITE(6,1330)
     1 READA(KEYWRD,INDEX(KEYWRD,'OMIN='))
      IF (WORD(ALLKEY,' RMIN')) WRITE(6,1340)
     1 READA(KEYWRD,INDEX(KEYWRD,' RMIN='))
      IF (WORD(ALLKEY,' RMAX')) WRITE(6,1350)
     1 READA(KEYWRD,INDEX(KEYWRD,' RMAX='))
      IF (WORD(ALLKEY,' DDMIN')) WRITE(6,1360)
     1 READA(KEYWRD,INDEX(KEYWRD,' DDMIN='))
      IF (WORD(ALLKEY,' DDMAX')) WRITE(6,1370)
     1 READA(KEYWRD,INDEX(KEYWRD,' DDMAX='))
      IF (WORD(ALLKEY,' PRNT')) WRITE(6,1380)
C
C	END of EF section						IR0494
C
      IF (WORD(ALLKEY,' LET') ) WRITE(6,890)
      IF (WORD(ALLKEY,'COMPFG') ) WRITE(6,900)
      IF (WORD(ALLKEY,'DERIV') ) WRITE(6,910)
      IF (WORD(ALLKEY,'FULSCF') ) WRITE(6,920)
      IF (WORD(ALLKEY,'DCART') ) WRITE(6,930)
      IF (WORD(ALLKEY,'GNORM') ) WRITE(6,940)
     1 READA(KEYWRD,INDEX(KEYWRD,'GNORM'))
      IF (WORD(ALLKEY,'FMAT') ) WRITE(6,950)
      IF (WORD(ALLKEY,'HCORE') ) WRITE(6,960)
      IF (WORD(ALLKEY,'ITER') ) WRITE(6,970)
      IF (WORD(ALLKEY,'PULAY') ) WRITE(6,980)
      IF (WORD(ALLKEY,'LINMIN') ) WRITE(6,990)
      IF (WORD(ALLKEY,'LOCMIN') ) THEN
         WRITE(6,1000)
         I=INDEX(KEYWRD,'LOCMIN')
         KEYWRD(I:I+5)='LINMIN'
      ENDIF
      IF (WORD(ALLKEY,' STEP1')  )WRITE(6,1010)
     1 READA(KEYWRD,INDEX(KEYWRD,'STEP1')+6)
      IF (WORD(ALLKEY,' STEP2')  )WRITE(6,1020)
     1 READA(KEYWRD,INDEX(KEYWRD,'STEP2')+6)
      IF (WORD(ALLKEY,'BAR') ) WRITE(6,1030)
     1 READA(KEYWRD,INDEX(KEYWRD,'BAR'))
      IF (WORD(ALLKEY,'DEBUGPULAY') ) WRITE(6,1040)
      IF (WORD(ALLKEY,'CAMP') ) WRITE(6,1050)
      IF (WORD(ALLKEY,'KING') ) WRITE(6,1050)
      IF (WORD(ALLKEY,'EIGS') ) WRITE(6,1060)
      IF (WORD(ALLKEY,'MOLDAT') ) WRITE(6,1070)
      IF (WORD(ALLKEY,'HYPERF') ) WRITE(6,1080)
      IF (WORD(ALLKEY,'OPCI') ) WRITE(6,1090)
      IF (WORD(ALLKEY,' PL') ) WRITE(6,1100)
      IF (WORD(ALLKEY,'SEARCH') ) THEN
         WRITE(6,1110)
         I=INDEX(KEYWRD,'SEARCH')
         KEYWRD(I:I+5)='LINMIN'
      ENDIF
      IF (WORD(ALLKEY,'CYCLES') ) WRITE(6,1130)
     1 NINT(READA(KEYWRD,INDEX(KEYWRD,'CYCLES')))
      IF (WORD(ALLKEY,'FILL') ) WRITE(6,1120)
     1 NINT(READA(KEYWRD,INDEX(KEYWRD,'FILL')))
      IF (WORD(ALLKEY,'ITRY') ) WRITE(6,1160)
     1 NINT(READA(KEYWRD,INDEX(KEYWRD,'ITRY')))
      IF (WORD(ALLKEY,'0SCF') ) WRITE(6,1180)
      IF(UHF)THEN
         IF(BIRAD.OR.EXCI.OR.CI)THEN
            WRITE(6,'(//10X,
     1'' UHF USED WITH EITHER BIRAD, EXCITED OR C.I. '')')
            WRITE(6,1170)
            GOTO 60
         ENDIF
      ELSE
         IF(EXCI.AND. TRIP) THEN
            WRITE(6,'(//10X,'' EXCITED USED WITH TRIPLET'')')
            WRITE(6,1170)
            GOTO 60
         ENDIF
      ENDIF
      IF (INDEX(KEYWRD,'T-PRIO').NE.0.AND.
     1INDEX(KEYWRD,'DRC').EQ.0) THEN
         WRITE (6,'(//10X,''T-PRIO AND NO DRC'')')
         WRITE (6,1170)
         GOTO 60
      ENDIF
      IF ( METHOD .GT. 1) THEN
         WRITE(6,'(//10X,
     1'' ONLY ONE OF MINDO, MNDO, AM1 AND PM3 ALLOWED'')')
         WRITE (6,1170)
         GOTO 60
      ENDIF
      IF (WORD(ALLKEY,'THERMO') )THEN
         WRITE(6,1140)
         IF(WORD(ALLKEY,' ROT')) THEN
            WRITE(6,1150)NINT(READA(KEYWRD,INDEX(KEYWRD,' ROT')))
         ELSE
            WRITE(6,'
     1    (//10X,'' YOU MUST SUPPLY THE SYMMETRY NUMBER "ROT"'')')
            STOP
         ENDIF
      ENDIF
      IF(MAXGEO.GT.1)THEN
         WRITE(6,'(//10X,''MORE THAN ONE GEOMETRY OPTION HAS BEEN '',
     1''SPECIFIED'',/10X,
     2''CONFLICT MUST BE RESOLVED BEFORE JOB WILL RUN'')')
         STOP
      ENDIF
      IF(ALLKEY.NE.' ')THEN
C
C  AN UNRECOGNIZED KEY-WORD HAS BEEN USED. FIND IT AND THEN ABANDON THE
C  RUN
         J=0
         DO 50 I=1,79
            IF(ALLKEY(I:I).NE.' '.OR.ALLKEY(I:I+1).NE.'  ')THEN
               J=J+1
               KEYWRD(J:J)=ALLKEY(I:I)
            ENDIF
   50    CONTINUE
         IF(ALLKEY(80:80).NE.' ')THEN
            J=J+1
            KEYWRD(J:J)=ALLKEY(80:80)
         ENDIF
         J=MAX(1,J)
         WRITE(6,'(///10X,''UNRECOGNIZED KEY-WORDS: ('',A,'')'')')
     1KEYWRD(:J)
         WRITE(6,'(///10X,''CALCULATION STOPPED TO AVOID WASTING TIME.''
     1)')
         STOP
      ENDIF
      RETURN
   60 WRITE(6,'(//10X,'' CALCULATION ABANDONED, SORRY!'')')
      STOP
   70 FORMAT(' *  VECTORS  - FINAL EIGENVECTORS TO BE PRINTED')
   80 FORMAT(' *  EXTERNAL - USE ATOMIC PARAMETERS FROM THE FOLLOWING '
     1,'FILE',/15X,A)
   90 FORMAT(' *  DENSITY  - FINAL DENSITY MATRIX TO BE PRINTED')
  100 FORMAT(' *  SPIN     - FINAL UHF SPIN MATRIX TO BE PRINTED')
  110 FORMAT(' *  DEPVAR=N - SPECIFIED DISTANCE IS',F7.4,
     1' TIMES BOND LENGTH')
  120 FORMAT(' *  DEP      - OUTPUT FORTRAN CODE FOR BLOCK-DATA')
  130 FORMAT(' *  VELOCITY - INPUT STARTING VELOCITIES FOR DRC')
  140 FORMAT(' *  TIMES    - TIMES OF VARIOUS STAGES TO BE PRINTED')
  150 FORMAT(' *  PARASOK  - USE SOME MNDO PARAMETERS IN AN AM1 CALCULA'
     1,'TION')
  160 FORMAT(' *  FLEPO    - PRINT DETAILS OF GEOMETRY OPTIMISATION')
  170 FORMAT(' *  BONDS    - FINAL BOND-ORDER MATRIX TO BE PRINTED')
  180 FORMAT(' *  GEO-OK   - OVERRIDE INTERATOMIC DISTANCE CHECK')
  190 FORMAT(' *  FOCK     - LAST FOCK MATRIX TO BE PRINTED')
  200 FORMAT(' *  LARGE    - EXPANDED OUTPUT TO BE PRINTED')
  210 FORMAT(' *  S1978    - 1978 SULFUR PARAMETERS TO BE USED')
  220 FORMAT(' *  SI1978   - 1978 SILICON PARAMETERS TO BE USED')
  230 FORMAT(' *  GRAPH    - GENERATE FILE FOR GRAPHICS')
  240 FORMAT(' *  1ELECTRON- FINAL ONE-ELECTRON MATRIX TO BE PRINTED')
  250 FORMAT(' *  ESR      - RHF SPIN DENSITY CALCULATION REQUESTED')
  260 FORMAT(' *  NOMM     - DO NOT MAKE MM CORRECTION TO CONH BARRIER')
  265 FORMAT(' *  MMOK     - APPLY MM CORRECTION TO CONH BARRIER')
  270 FORMAT(' *  INTERP   - PRINT DETAILS OF CAMP-KING CONVERGER')
  280 FORMAT(' *  DFP      - USE DAVIDON FLETCHER POWELL OPTIMIZER')
  290 FORMAT(' *  ANALYT   - USE ANALYTIC DERIVATIVES ')
  300 FORMAT(' *  MECI     - M.E.C.I. WORKING TO BE PRINTED')
  310 FORMAT(' *  DRC      - DYNAMIC REACTION COORDINATE CALCULATION')
  320 FORMAT(' *  DRC=     - HALF-LIFE FOR KINETIC ENERGY LOSS =',F9.2,
     1' * 10**(-14) SECONDS')
  330 FORMAT(' *  KINETIC= - ',F7.3,' KCAL KINETIC ENERGY ADDED TO DRC')
  340 FORMAT(' *  LOCALIZE - LOCALIZED ORBITALS TO BE PRINTED')
  350 FORMAT(' *  MULLIK   - THE MULLIKEN ANALYSIS TO BE PERFORMED')
  360 FORMAT(' *   XYZ     - CARTESIAN COORDINATE SYSTEM TO BE USED')
  370 FORMAT(' *   PI      - BONDS MATRIX, SPLIT INTO SIGMA-PI-DELL',
     1' COMPONENTS, TO BE PRINTED')
  380 FORMAT(' *  ECHO     - ALL INPUT DATA TO BE ECHOED BEFORE RUN')
  390 FORMAT(' *  H-PRIOR  - HEAT OF FORMATION TAKES PRIORITY IN DRC')
  400 FORMAT(' *  X-PRIOR  - GEOMETRY CHANGES TAKE PRIORITY IN DRC')
  410 FORMAT(' *  T-PRIOR  - TIME TAKES PRIORITY IN DRC')
  420 FORMAT(' *  POWSQ    - PRINT DETAILS OF WORKING IN POWSQ')
  430 FORMAT(' *  POLAR    - CALCULATE FIRST, SECOND AND THIRD-ORDER'
     1,' POLARIZABILITIES')
  440 FORMAT(' *  DEBUG    - DEBUG OPTION TURNED ON')
  450 FORMAT(' *  RESTART  - CALCULATION RESTARTED')
  460 FORMAT(' *  IRC=N    - INTRINSIC REACTION COORDINATE',I3,
     1' DEFINED')
  470 FORMAT(' *  IRC      - INTRINSIC REACTION COORDINATE CALCULATION')
  480 FORMAT(3(' *',/),' *',15X,'  CHARGE ON SYSTEM =',I3,3(/,' *'))
  490 FORMAT(' *  GRADIENTS- ALL GRADIENTS TO BE PRINTED')
  500 FORMAT(' *  UHF      - UNRESTRICTED HARTREE-FOCK CALCULATION')
  510 FORMAT(' *  SINGLET  - STATE REQUIRED MUST BE A SINGLET')
  520 FORMAT(' *  BIRADICAL- SYSTEM HAS TWO UNPAIRED ELECTRONS')
  530 FORMAT(' *  EXCITED  - FIRST EXCITED STATE IS TO BE OPTIMIZED')
  540 FORMAT(' *  SINGLET  - SPIN STATE DEFINED AS A SINGLET')
  550 FORMAT(' *  DOUBLET  - SPIN STATE DEFINED AS A DOUBLET')
  560 FORMAT(' *  TRIPLET  - SPIN STATE DEFINED AS A TRIPLET')
  570 FORMAT(' *  QUARTET  - SPIN STATE DEFINED AS A QUARTET')
  580 FORMAT(' *  QUINTET  - SPIN STATE DEFINED AS A QUINTET')
  590 FORMAT(' *  SEXTET   - SPIN STATE DEFINED AS A SEXTET')
  600 FORMAT(' *  SYMMETRY - SYMMETRY CONDITIONS TO BE IMPOSED')
  610 FORMAT(' *  MICROS=N -',I4,' MICROSTATES TO BE SUPPLIED FOR C.I.')
  620 FORMAT(' *  OPEN(N,N)- THERE ARE',I2,' ELECTRONS IN',I2,' LEVELS')
  630 FORMAT(' *   T=      - A TIME OF',F8.1,' ',A7,' REQUESTED')
  640 FORMAT(' *   T=      - A TIME OF',G11.3,' ',A7,' REQUESTED')
  650 FORMAT(' *  DUMP=N   - RESTART FILE WRITTEN EVERY',F8.1,
     1' ',A7)
  660 FORMAT(' *  DUMP=N   - RESTART FILE WRITTEN EVERY',G11.3,
     1' ',A7)
  670 FORMAT(' *  1SCF     - DO 1 SCF AND THEN STOP ')
  680 FORMAT(' *  C.I.=N   -',I2,' M.O.S TO BE USED IN C.I.')
  681 FORMAT(' * C.I.=(N,M)-',I2,' DOUBLY FILLED LEVELS USED IN A ',/
     1,      ' *             C.I. INVOLVING ',I2,' M.O.''S')
  690 FORMAT(' *  DFORCE   - PRINT HESSIAN MATRIX IN FORCE')
  700 FORMAT(' *  FORCE    - FORCE CALCULATION SPECIFIED')
  710 FORMAT(' *  MINDO/3  - THE MINDO/3 HAMILTONIAN TO BE USED')
  720 FORMAT(' *  AM1      - THE AM1 HAMILTONIAN TO BE USED')
  730 FORMAT(' *  PM3      - THE PM3 HAMILTONIAN TO BE USED')
  740 FORMAT(' *  PRECISE  - CRITERIA TO BE INCREASED BY 100 TIMES')
  750 FORMAT(' *  NOINTER  - INTERATOMIC DISTANCES NOT TO BE PRINTED')
  760 FORMAT(' *  ISOTOPE  - FORCE MATRIX WRITTEN TO DISK (CHAN. 9 )')
  770 FORMAT(' *  DENOUT   - DENSITY MATRIX OUTPUT ON CHANNEL 10')
  780 FORMAT(' *  SHIFT    - A DAMPING FACTOR OF',F8.2,' DEFINED')
  790 FORMAT(' *  OLDENS   - INITIAL DENSITY MATRIX READ OF DISK')
  800 FORMAT(' *  SCFCRT   - DEFAULT SCF CRITERION REPLACED BY',G12.3)
  810 FORMAT(' *  ENPART   - ENERGY TO BE PARTITIONED INTO COMPONENTS')
  820 FORMAT(' *  NOXYZ    - CARTESIAN COORDINATES NOT TO BE PRINTED')
  830 FORMAT(' *  SIGMA    - GEOMETRY TO BE OPTIMIZED USING SIGMA.')
  840 FORMAT(' *  NLLSQ    - GRADIENTS TO BE MINIMIZED USING NLLSQ.')
  850 FORMAT(' *  ROOT     - IN A C.I. CALCULATION, ROOT',I2,
     1                       ' TO BE OPTIMIZED.')
  860 FORMAT(' *  TRANS    - THE REACTION VIBRATION TO BE DELETED FROM',
     1' THE THERMO CALCULATION')
  870 FORMAT(' *  TRANS=   - ',I4,' VIBRATIONS ARE TO BE DELETED FROM',
     1' THE THERMO CALCULATION')
  880 FORMAT(' *  SADDLE   - TRANSITION STATE TO BE OPTIMIZED')
  890 FORMAT(' *   LET     - OVERRIDE SOME SAFETY CHECKS')
  900 FORMAT(' *  COMPFG   - PRINT HEAT OF FORMATION CALC''D IN COMPFG')
  910 FORMAT(' *  DERIV    - PRINT PART OF WORKING IN SUB. DERIV')
  920 FORMAT(' *  FULSCF   - IN SEARCHES, FULL SCF CALCN''S TO BE DONE')
  930 FORMAT(' *  DCART    - PRINT DETAILS OF WORKING IN DCART')
  940 FORMAT(' *  GNORM=   - EXIT WHEN GRADIENT NORM DROPS BELOW ',G8.3)
  950 FORMAT(' *  FMAT     - PRINT DETAILS OF WORKING IN FMAT')
  960 FORMAT(' *  HCORE    - PRINT DETAILS OF WORKING IN HCORE')
  970 FORMAT(' *  ITER     - PRINT DETAILS OF WORKING IN ITER')
  980 FORMAT(' *  PULAY    - PULAY''S METHOD TO BE USED IN SCF')
  990 FORMAT(' *  LINMIN   - PRINT DETAILS OF LINE MINIMIZATION')
 1000 FORMAT(' *  LOCMIN   - USE LINMIN INSTEAD OF THIS KEYWORD')
 1010 FORMAT(' *  STEP1    - FIRST  STEP-SIZE IN GRID =',F7.2)
 1020 FORMAT(' *  STEP2    - SECOND STEP-SIZE IN GRID =',F7.2)
 1030 FORMAT(' *  BAR=     - REDUCE BAR LENGTH BY A MAX. OF',F7.2)
 1040 FORMAT(' *  DEBUGPULAY-PRINT DETAILS OF WORKING IN PULAY')
 1050 FORMAT(' *  CAMP,KING- THE CAMP-KING CONVERGER TO BE USED')
 1060 FORMAT(' *  EIGS     - PRINT ALL EIGENVALUES IN ITER')
 1070 FORMAT(' *  MOLDAT   - PRINT DETAILS OF WORKING IN MOLDAT')
 1080 FORMAT(' *  HYPERFINE- HYPERFINE COUPLING CONSTANTS TO BE'
     1,' PRINTED')
 1090 FORMAT(' *  OPCI     - PRINT DETAILS OF WORKING IN OPCI')
 1100 FORMAT(' *   PL      - MONITOR CONVERGANCE IN DENSITY MATRIX')
 1110 FORMAT(' *  SEARCH   - USE LINMIN INSTEAD OF THIS KEYWORD')
 1120 FORMAT(' *  FILL=    - IN RHF CLOSED SHELL, FORCE M.O.',I3,' TO BE
     1 FILLED')
 1130 FORMAT(' *  CYCLES=  - DO A MAXIMUM OF ',I4,' STEPS')
 1140 FORMAT(' *  THERMO   - THERMODYNAMIC QUANTITIES TO BE CALCULATED')
 1150 FORMAT(' *  ROT      - SYMMETRY NUMBER OF',I3,' SPECIFIED')
 1160 FORMAT(' *  ITRY=    - DO A MAXIMUM OF',I6,' ITERATIONS FOR SCF')
 1170 FORMAT( //10X,' IMPOSSIBLE OPTION REQUESTED,')
 1180 FORMAT(' *  0SCF     - AFTER READING AND PRINTING DATA, STOP')
C
C	EF-specific FORMATs                                             IR0494
C
 1190 FORMAT(' *  EF       - USE EF ROUTINE FOR MINIMUM SEARCH')
 1200 FORMAT(' *  TS       - USE EF ROUTINE FOR TS SEARCH')
 1210 FORMAT(' *  IUPD=0   - EF: HESSIAN WILL NOT BE UPDATED')
 1220 FORMAT(' *  IUPD=1   - EF: HESSIAN UPDATED USING POWELL')
 1230 FORMAT(' *  IUPD=2   - EF: HESSIAN UPDATED USING BFGS')
 1240 FORMAT(' *  NONR     - EF: DO NOT USE NEWTON-RAPHSON STEP')
 1250 FORMAT(' *  HESS=0   - EF: DIAGONAL HESSIAN AS INITIAL GUESS')
 1260 FORMAT(' *  HESS=1   - EF: HESSIAN CALCULATED NUMERICALLY')
 1270 FORMAT(' *  HESS=2   - EF: INITIAL HESSIAN READ FROM DISK')
 1280 FORMAT(' *  HESS=3   - EF: HESSIAN CALCULATED DOUBLE NUMERICALLY')    
 1290 FORMAT(' *  MODE=    - EF: FOLLOW HESSIAN MODE',I3,' TOWARD TS')
 1300 FORMAT(' *  RSCAL    - EF: SCALE P-RFO STEP TO TRUST RADIUS')
 1310 FORMAT(' *  RECALC=  - EF: HESSIAN RECALC EVERY',I4,' CYCLES ')
 1320 FORMAT(' *  DMAX=    - EF: STARTING TRUST RADIUS',F7.3,' ANG/RAD')   
 1330 FORMAT(' *  OMIN=    - EF: MINIMUM EIGENVECTOR OVERLAP ',F7.3)
 1340 FORMAT(' *  RMIN=    - EF: MIN. CALC./PRED. ENERGY STEP ',F7.3)
 1350 FORMAT(' *  RMAX=    - EF: MAX. CALC./PRED. ENERGY STEP ',F7.3)
 1360 FORMAT(' *  DDMIN=   - EF: MINIMUM TRUST RADIUS',F7.3,' ANG/RAD')
 1370 FORMAT(' *  DDMAX=   - EF: MAXIMUM TRUST RADIUS',F7.3,' ANG/RAD')
 1380 FORMAT(' *  PRNT     - EF: EXTRA PRINTING')
C
C	End of  EF-specific FORMATs
C
      END
