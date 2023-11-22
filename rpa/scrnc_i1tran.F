      SUBROUTINE SCRNC_I1TRAN(ICASE,I2MO,I2MIX,EVEC1,EVEC2,
     &                        SCR,NAO,NMO,ISPIN,IUHF,
     &                        IRREPR,IRREPL,MOTOAO,SIDE)
C
C THIS ROUTINE TRANSFORMS INTEGRAL VECTORS BETWEEN THE MO AND MIXED
C AO/MO REPRESENTATIONS:

C        I2(AI,BJ) ==> I2(XX,XX) [ICASE=1,2]
C        I2(XX,XX) ==> I2(IJ,XX) [ICASE=3]
C        I2(XX,XX) ==> I2(AB,XX) [ICASE=4]
C        I2(XX,XX) ==> I2(IA,XX) [ICASE=5]
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      LOGICAL AOTOMO,MOTOAO
      CHARACTER*4 SIDE 
      CHARACTER*8 LABELX(2),LABELY(2)
      DOUBLE PRECISION I2MO(*),I2MIX(*)
      DIMENSION EVEC1(NAO,NMO),EVEC2(NAO,NMO),SCR(*)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)

      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      DATA LABELX /'SCFEVECA','SCFEVECB'/
      DATA LABELY /'SCFEVECA','SCFEVECB'/
      DATA ZILCH,ONE /0.0D0,1.0D0/
C
      NNM1O2(I)=(I*(I-1))/2
C
      IF (ICASE .EQ. 1) THEN
         NUMBJDIS=IRPDPD(IRREPR,ISYTYP(2,4+ISPIN))
      ELSE IF (ICASE .EQ. 2 .OR. ICASE .EQ. 3 .OR.
     +         ICASE .EQ. 4 .OR. ICASE .EQ. 5) THEN
         NUMBJDIS=IRPDPDAO(IRREPR)
      ENDIF 
#ifdef _DEBUG_LVLM
      Write(6,"(a,I5)") " The number of distributions:", numbjdis 
#endif 
      IOFFMO =1
      IOFFMIX=1
      IOFFXX =1
      IOFFOO =1
      IOFFVV =1
      IOFFVO =1
C
      IF(IUHF .EQ. 0)THEN
       CALL GETREC(20,'JOBARC',LABELX(ISPIN),IINTFP*NAO*NMO,EVEC1)
       ISPIN1=ISPIN
       ISPIN2=ISPIN
      ELSE
       CALL GETREC(20,'JOBARC',LABELX(1),IINTFP*NAO*NMO,EVEC1)
       IF (ISPIN .EQ. 2) CALL GETREC(20,'JOBARC',LABELX(1+IUHF),
     &                               IINTFP*NAO*NMO,EVEC1)
       ISPIN1=1
       ISPIN2=2
      ENDIF 

#ifdef _DEBUG_LVLM
       write(6,*) "The (C) matrix"
       CALL output(EVEC1,1,nao,1,nao,nao,nao,1)
#endif 
       IF (MOTOAO .AND. SIDE .EQ. "LEFT ") THEN
       IOVL1 = 1
       IOVL2 = IOVL1 + NAO * NAO
       IOVL3 = IOVL2 + NAO * NAO
       IF (IEND .GT. MAXCOR) CALL INSMEM("@-SCRNC_I1TRAN",IEND,MAXCOR)

       CALL GETREC(20,"JOBARC","AOOVRLAP",NAO*NAO*IINTFP,SCR(IOVL1))

#ifdef _DEBUG_LVLM
       write(6,*) "The S matrix"
       CALL output(scr(iovl1),1,nao,1,nao,nao,nao,1)
#endif 
       CALL XGEMM("N","N",NMO,NAO,NAO,ONE,SCR(IOVL1),NAO,EVEC1,NAO,
     &             ZILCH,SCR(IOVL2),NMO)
       CALL DCOPY(NAO*NMO,SCR(IOVL2),1,EVEC1,1)

#ifdef _DEBUG_LVLM
       write(6,*) "The SC matrix"
       CALL output(evec1,1,nao,1,nao,nao,nao,1)
       CALL GETREC(20,"JOBARC",LABELX(1),NAO*NAO*IINTFP,SCR(IOVL1))
       CALL XGEMM("T","N",NMO,NAO,NAO,ONE,SCR(IOVL1),NAO,EVEC1,NAO,
     &             ZILCH,SCR(IOVL2),NMO)
       CALL output(scr(iovl2),1,nao,1,nao,nao,nao,1)
#endif 
       ENDIF 
C
      IF(ICASE.EQ.1 .OR. ICASE .EQ. 2)THEN
C
C LOOP OVER DISTRIBUTIONS
C
       DO 10 IBJ=1,NUMBJDIS
C
C LOOP OVER SYMMETRIES OF I(AI,BJ)->(XX,BJ)
C
        IRREPBJ = IRREPR
        IRREPAI = IRREPL
       
        DO 20 IRREPI=1,NIRREP
         IRREPA=DIRPRD(IRREPI,IRREPAI)
         NUMA=VRT(IRREPA,ISPIN)
         NUMI=POP(IRREPI,ISPIN)
         NAO1=IAOPOP(IRREPA)
         NAO2=IAOPOP(IRREPI)
         IOFFR1=IOFFV(IRREPA,ISPIN)
         IOFFL1=IOFFAO(IRREPA)
         IOFFR2=IOFFO(IRREPI,ISPIN)
         IOFFL2=IOFFAO(IRREPI)

         IF (MIN(NAO1,NUMI,NUMA).NE.0) THEN
            CALL XGEMM('N','N',NAO1,NUMI,NUMA,ONE,EVEC1(IOFFL1,IOFFR1),
     &                 NAO,I2MO(IOFFMO),NUMA,ZILCH,SCR,NAO1)
         ELSE
            CALL DZERO(SCR,NAO1*NUMI)
         END IF

#ifdef _DEBUG_LVLM
      if (ispin.eq.1) then
      Write(6,"(a,a,3i2)") " I(AI,BJ)->(XI,BJ);"," irrepa,irrepi,Nbj:",
     +irrepa,irrepi,ibj
      write(6,*) "CMO1"
      call output(EVEC1(IOFFL1,IOFFR1),1,NAO1,1,NUMA,NAO1,NUMA,1)
      write(6,*) "I(AI)"
      call output(I2MO(IOFFMO),1,NUMA,1,NUMI,NUMA,NUMI,1)
      write(6,*) "I(XI)"
      call output(SCR,1,NAO1,1,NUMI,NAO1,NUMI,1) 
      endif 
#endif 
        IF (MIN(NAO1,NAO2,NUMI).NE.0) THEN
            CALL XGEMM('N','T',NAO1,NAO2,NUMI,ONE,SCR,NAO1,
     &                 EVEC1(IOFFL2,IOFFR2),NAO,
     &                 ZILCH,I2MIX(IOFFMIX),NAO1)
         ELSE
            CALL DZERO(I2MIX(IOFFMIX),NAO1*NAO2)
         END IF

#ifdef _DEBUG_LVLM
      if (ispin .eq. 1) then
      write(6,*) "nao1 and nao2", nao1,nao2
      write(6,*) "CMO2"
      call output(EVEC1(IOFFL2,IOFFR2),1,NAO2,1,NUMI,NAO2,NUMI,1)
      write(6,*) "I(XX)"
      call output(I2MIX(IOFFMIX),1,NAO1,1,NAO2,NAO1,NAO2,1) 
      endif 
#endif 
         IOFFMO =IOFFMO +NUMA*NUMI
         IOFFMIX=IOFFMIX+NAO1*NAO2
20      CONTINUE
10     CONTINUE 

      ELSE IF (ICASE .EQ. 3) THEN 

C LOOP OVER SYMMETRIES OF I(XX,XX)->(IJ,XX)

       IRREPXX = IRREPR
       IRREPIJ = IRREPL

       DO 30 IBJ=1,NUMBJDIS
#ifdef _DEBUG_LVLM
      If (ispin.eq.1) then
      write(6,*)
      Write(6,"(a,i4)") "The XX distribution number:", ibj
      endif 
#endif
        DO 40 IRREPJ=1,NIRREP 
         IRREPI=DIRPRD(IRREPJ,IRREPIJ)
         NUMI=POP(IRREPI,ISPIN)
         NUMJ=POP(IRREPJ,ISPIN)
         NAO1=IAOPOP(IRREPI)
         NAO2=IAOPOP(IRREPJ)

         IOFFR1=IOFFO(IRREPI,ISPIN)
         IOFFL1=IOFFAO(IRREPI)

         IOFFR2=IOFFO(IRREPJ,ISPIN)
         IOFFL2=IOFFAO(IRREPJ)

#ifdef _DEBUG_LVLM
      if (ispin .eq.1) then
      Write(6,"(a,a,2I2)") "I(XX,XX)->(IX,XX)","Irrepl,irrepr:",
     +Irrepj,irrepi
      write(6,*) "CMO"
      call output(EVEC1(IOFFL1,IOFFR1),1,NAO1,1,NUMI,NAO1,NUMI,1)
      write(6,*) "I2(XX)"
      call output(I2MIX(IOFFXX),1,NAO1,1,NAO2,NAO1,NAO2,1)
      endif 
#endif 
         IF (MIN(NUMI,NAO2,NAO1).NE.0) THEN
            CALL XGEMM('T','N',NUMI,NAO2,NAO1,ONE,EVEC1(IOFFL1,IOFFR1),
     &                 NAO,I2MIX(IOFFXX),NAO1,ZILCH,SCR,NUMI)
         ELSE
            CALL DZERO(SCR,NAO2*NUMI)
         END IF

#ifdef _DEBUG_LVLM
      if (ispin .eq.1) then 
      Write(6,"(a)") "I(IX,XX)->(IJ,XX)"
      write(6,*) "I2MO(IX)"
      call output(SCR,1,NUMI,1,NAO2,NUMI,NAO2,1)
      write(6,*) "CMO"
      call output(EVEC1(IOFFL2,IOFFR2),1,NAO2,1,NUMJ,NAO2,NUMJ,1)
      endif
#endif
         IF (MIN(NUMI,NUMJ,NAO2).NE.0) THEN
            CALL XGEMM('N','N',NUMI,NUMJ,NAO2,ONE,SCR,NUMI,
     &                 EVEC1(IOFFL2,IOFFR2),NAO,ZILCH,I2MO(IOFFOO),
     &                 NUMI)
         ELSE
            CALL DZERO(I2MO(IOFFOO),NUMI*NUMJ)
         END IF

#ifdef _DEBUG_LVLM
      write(6,"(a,i2)") "I2MO(IJ),Ispin:",ispin
      call output(I2MO(IOFFOO),1,NUMI,1,NUMJ,NUMI,NUMJ,1)
#endif 
         IOFFXX=IOFFXX+NAO1*NAO2
         IOFFOO=IOFFOO+NUMI*NUMJ
40      CONTINUE
30     CONTINUE
C
      ELSE IF (ICASE .EQ. 4) THEN

C LOOP OVER SYMMETRIES OF I(XX,XX)->(AB,XX)

       IRREPXX = IRREPR
       IRREPAB = IRREPL

       DO 50 IBJ=1,NUMBJDIS
#ifdef _DEBUG_LVLM
      If (ispin.eq.1) then
      write(6,*)
      Write(6,"(a,i4,i2)") "The XX distribution number:", ibj,ispin 
      endif
#endif
        DO 60 IRREPB=1,NIRREP
         IRREPA=DIRPRD(IRREPB,IRREPAB)
         NUMA=VRT(IRREPA,ISPIN)
         NUMB=VRT(IRREPB,ISPIN)
         NAO1=IAOPOP(IRREPA)
         NAO2=IAOPOP(IRREPB)

         IOFFR1=IOFFV(IRREPA,ISPIN)
         IOFFL1=IOFFAO(IRREPA)

         IOFFR2=IOFFV(IRREPB,ISPIN)
         IOFFL2=IOFFAO(IRREPB)

#ifdef _DEBUG_LVLM
c      if (ispin .eq.1) then
      Write(6,"(a,i2)") "I(XX,XX)->(AX,XX)",irrepl,irrepr:",
     +irrepa,irrepb
      write(6,*) "CMO"
      call output(EVEC1(IOFFL1,IOFFR1),1,NAO1,1,NUMA,NAO1,NUMA,1)
      write(6,*) "I2(XX)"
      call output(I2MIX(IOFFXX),1,NAO1,1,NAO2,NAO1,NAO2,1)
c      endif 
#endif 
         IF (MIN(NUMA,NAO2,NAO1).NE.0) THEN
            CALL XGEMM('T','N',NUMA,NAO2,NAO1,ONE,EVEC1(IOFFL1,IOFFR1),
     &                 NAO,I2MIX(IOFFXX),NAO1,ZILCH,SCR,NUMA)
         ELSE
            CALL DZERO(SCR,NAO2*NUMA)
         END IF

#ifdef _DEBUG_LVLM
c      if (ispin .eq. 1) then 
      Write(6,"(a)") "I(AX,XX)->(AB,XX)"
      write(6,*) "CMO"
      call output(EVEC1(IOFFL2,IOFFR2),1,NAO2,1,NUMB,NAO1,NUMB,1)
      write(6,*) "I2MO(AX)"
      call output(SCR,1,NUMA,1,NAO2,NUMA,NAO2,1)
c      endif
#endif
         IF (MIN(NUMA,NUMB,NAO2).NE.0) THEN
            CALL XGEMM('N','N',NUMA,NUMB,NAO2,ONE,SCR,NUMA,
     &                 EVEC1(IOFFL2,IOFFR2),NAO,ZILCH,I2MO(IOFFVV),
     &                 NUMA)
         ELSE
            CALL DZERO(I2MO(IOFFVV),NUMA*NUMB)
         END IF

#ifdef _DEBUG_LVLM
      write(6,"(a,i2)") "I2MO(AB),Ispin:",ispin 
      call output(I2MO(IOFFVV),1,NUMA,1,NUMB,NUMA,NUMB,1)
#endif
         IOFFXX=IOFFXX+NAO1*NAO2
         IOFFVV=IOFFVV+NUMA*NUMB
60      CONTINUE
50     CONTINUE

C LOOP OVER SYMMETRIES OF I(XX,XX)->(AI,XX)

      ELSE IF (ICASE .EQ. 5) THEN

       IRREPXX = IRREPR
       IRREPAI = IRREPL

       DO 70 IBJ=1,NUMBJDIS
#ifdef _DEBUG_LVLM
      If (ispin.eq.1) then
      write(6,*)
      Write(6,"(a,i4)") "The XX distribution number:", ibj
      endif
#endif
        DO 80 IRREPI=1,NIRREP
         IRREPA=DIRPRD(IRREPI,IRREPAI)
         NUMA=VRT(IRREPA,ISPIN)
         NUMI=POP(IRREPI,ISPIN)
         NAO1=IAOPOP(IRREPA)
         NAO2=IAOPOP(IRREPI)

         IOFFR1=IOFFV(IRREPA,ISPIN)
         IOFFL1=IOFFAO(IRREPA)

         IOFFR2=IOFFO(IRREPI,ISPIN)
         IOFFL2=IOFFAO(IRREPI)

#ifdef _DEBUG_LVLM
      if (ispin .eq.1) then
      Write(6,"(a,a,2I2)") "I(XX,XX)->(AX,XX)","Irrepl,irrepr:",
     +Irrepa,irrepi
      write(6,*) "CMO"
      call output(EVEC1(IOFFL1,IOFFR1),1,NAO1,1,NUMA,NAO1,NUMA,1)
      write(6,*) "I2(XX)"
      call output(I2MIX(IOFFXX),1,NAO1,1,NAO2,NAO1,NAO2,1)
      endif 
#endif
         IF (MIN(NUMA,NAO2,NAO1).NE.0) THEN
            CALL XGEMM('T','N',NUMA,NAO2,NAO1,ONE,EVEC1(IOFFL1,IOFFR1),
     &                 NAO,I2MIX(IOFFXX),NAO1,ZILCH,SCR,NUMA)
         ELSE
            CALL DZERO(SCR,NAO2*NUMA)
         END IF
#ifdef _DEBUG_LVLM
      if (ispin .eq.1) then
      Write(6,"(a)") "I(IX,XX)->(AI,XX)"
      write(6,*) "I2MO(AX)"
      call output(SCR,1,NUMA,1,NAO2,NUMA,NAO2,1)
      write(6,*) "CMO"
      call output(EVEC1(IOFFL2,IOFFR2),1,NAO2,1,NUMI,NAO2,NUMI,1)
      endif
#endif
         IF (MIN(NUMA,NUMI,NAO2).NE.0) THEN
            CALL XGEMM('N','N',NUMA,NUMI,NAO2,ONE,SCR,NUMA,
     &                  EVEC1(IOFFL2,IOFFR2),NAO,ZILCH,I2MO(IOFFVO),
     &                  NUMA)
         ELSE
            CALL DZERO(I2MO(IOFFVO),NUMI*NUMA)
         END IF

#ifdef _DEBUG_LVLM
      write(6,"(a,i2)") "I2MO(AI),Ispin:"ispin
      call output(I2MO(IOFFVO),1,NUMA,1,NUMI,NUMA,NUMI,1)
#endif
         IOFFXX=IOFFXX+NAO1*NAO2
         IOFFVO=IOFFVO+NUMI*NUMA
80      CONTINUE
70     CONTINUE

      ENDIF 

      RETURN
      END      

