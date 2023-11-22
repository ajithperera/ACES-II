      SUBROUTINE DG6TRAN(G6MO,G6MIX,EVEC1,EVEC2,SCR1,SCR2,
     &                   NAO,NMO,NUMABDIS,NUMAIDIS,
     &                   NUMXXDIS,NUMXIDIS,ISPIN,IRREPAI,
     &                   IRREPAB,IUHF)
C
C THIS ROUTINE TRANSFORMS G6 VECTORS BETWEEN THE MO AND MIXED
C AO/MO REPRESENTATIONS:
C
C G6(XX,XI) <== G6(AB,CI) 
C
C INPUT:
C
C   G6MO...........G(AB,CI) IN FULL MO REPRESENTATION
C   EVEC1..........SCRATCH ARRAY FOR SCF EIGENVECTORS (SPIN CASE 1)
C   EVEC2..........SCRACTH ARRAY FOR SCF EIGENVECTORS (SPIN CASE 2)
C   SCR1, SCR2.....SCRATCH ARRAYS
C   NAO............NUMBER OF ATOMIC ORBITALS
C   NMO............NUMBER OF MOLECULAR ORBITALS
C   NUMABDIS.......DISTRIBUTION LENGTH FOR G(AB,CI)
C   NUMAIDIS.......NUMBER OF DISTRIBUTIONS FOR G(AB,CI)
C   NUMXXDIS.......DISTRIBUTION LENGTH FOR G(XX,XI)
C   NUMXIDIS.......NUMBER OF DISTRIBUTIONS FOR G(XX,XI)
C   ISPIN..........SPIN CASE
C   IUHF...........UHF FLAG
C
C OUTPUT:
C
C   G6MIX..........G(XX,XI) IN MIXED AO/MO REPRESENTATION
C
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      CHARACTER*8 LABELX(2)
      LOGICAL DONE,STERM,ANTI,MOEQAO
C
      DIMENSION G6MO(*),G6MIX(*)
      DIMENSION EVEC1(NAO,NMO),EVEC2(NAO,NMO),SCR1(*),SCR2(*)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &             IRPDPDAO(8),IRPDPDAOS(8),ISTART(8,8),ISTARTMO(8,3)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/G6OOC/IPOPS(8,8),IPOPE(8,8),NDONE(8,8),NSTARTM,
     &             ISTARTI,NSZIG,IOFFSET,DONE
      DATA LABELX /'SCFEVECA','SCFEVECB'/
      DATA ZILCH,ONE,ONEM /0.0D0,1.0D0,-1.D0/
      DATA IONE/1/
C
      IOFF    =1
      IOFFMIX0=1
C     
      IF (ISPIN.LE.2) THEN
         CALL GETREC(20,'JOBARC',LABELX(ISPIN),IINTFP*NAO*NMO,EVEC1) 
         CALL SCOPY (NAO*NMO,EVEC1,1,EVEC2,1)
         ISPIN1=ISPIN
         ISPIN2=ISPIN
         ISPIN3=ISPIN
         ISPIN4=ISPIN
      ELSE 
         CALL GETREC(20,'JOBARC',LABELX(1),IINTFP*NAO*NMO,EVEC1)
         IF(IUHF.EQ.0) THEN
            CALL SCOPY (NAO*NMO,EVEC1,1,EVEC2,1)
         ELSE
            CALL GETREC(20,'JOBARC',LABELX(2),IINTFP*NAO*NMO,EVEC2)
         ENDIF
         IF(ISPIN.EQ.3) THEN
            ISPIN1=2
            ISPIN2=1
            ISPIN3=2
            ISPIN4=1
         ELSE
            ISPIN1=1
            ISPIN2=2
            ISPIN3=1
            ISPIN4=2
         ENDIF
      ENDIF
C     
      DO 9 IRREPI=1,NIRREP
C     
         IRREPC=DIRPRD(IRREPI,IRREPAI)
         NUMIE=IPOPE(IRREPI,IRREPAI)
         NUMIS=IPOPS(IRREPI,IRREPAI)
         NUMC=VRT(IRREPC,ISPIN3)
C     
         DO 10 II=NUMIS,NUMIE
C     
            NAO3=IAOPOP(IRREPC)
            IOFFL3=IOFFAO(IRREPC)
            IOFFR3=IOFFV(IRREPC,ISPIN3)
C
C Transform G(AB,CI) to G(AB,lambdaI)
C
            IF(ISPIN3.EQ.1) THEN
               CALL XGEMM('N','T',NUMABDIS,NAO3,NUMC,ONE,G6MO(IOFF),
     &                     NUMABDIS,EVEC1(IOFFL3,IOFFR3),NAO,ZILCH,
     &                     SCR1,NUMABDIS)
            ELSE
               CALL XGEMM('N','T',NUMABDIS,NAO3,NUMC,ONE,G6MO(IOFF),
     &                     NUMABDIS,EVEC2(IOFFL3,IOFFR3),NAO,ZILCH,
     &                     SCR1,NUMABDIS)
            ENDIF
C     
            IOFF=IOFF+NUMC*NUMABDIS
C     
            IOFFMO=1
C     
            DO 11 IC=1,NAO3
C     
               IOFFMIX=IOFFMIX0
C     
               DO 20 IRREPB=1,NIRREP
C     
                  IRREPA=DIRPRD(IRREPB,IRREPAB)
                  NUMA=VRT(IRREPA,ISPIN1)
                  NUMB=VRT(IRREPB,ISPIN2)
                  NAO1=IAOPOP(IRREPA)
                  NAO2=IAOPOP(IRREPB)
                  IOFFR1=IOFFV(IRREPA,ISPIN1)
                  IOFFR2=IOFFV(IRREPB,ISPIN2)
                  IOFFL1=IOFFAO(IRREPA)
                  IOFFL2=IOFFAO(IRREPB)
C
C Transform G(AB|lambdaI) to G(mu,nu,lambda,I)
C
                  IF(ISPIN1.EQ.1) THEN
                     CALL XGEMM('N','N',NAO1,NUMB,NUMA,ONE,EVEC1
     &                          (IOFFL1,IOFFR1),NAO,SCR1(IOFFMO),
     &                           NUMA,ZILCH,SCR2,NAO1)
                  ELSE
                     CALL XGEMM('N','N',NAO1,NUMB,NUMA,ONE,EVEC2
     &                          (IOFFL1,IOFFR1),NAO,SCR1(IOFFMO),
     &                           NUMA,ZILCH,SCR2,NAO1)
                  ENDIF
                  IF(ISPIN2.EQ.1) THEN
                     CALL XGEMM('N','T',NAO1,NAO2,NUMB,ONE,SCR2,NAO1,
     &                           EVEC1(IOFFL2,IOFFR2),NAO,ZILCH,
     &                           G6MIX(IOFFMIX),NAO1)
                  ELSE
                     CALL XGEMM('N','T',NAO1,NAO2,NUMB,ONE,SCR2,NAO1,
     &                           EVEC2(IOFFL2,IOFFR2),NAO,ZILCH,
     &                           G6MIX(IOFFMIX),NAO1)
                  ENDIF
C     
                  IOFFMIX=IOFFMIX+NAO1*NAO2
                  IOFFMO=IOFFMO+NUMA*NUMB
C     
 20            CONTINUE
C     
               IOFFMIX0=IOFFMIX0+NUMXXDIS
C     
 11         CONTINUE
C     
 10      CONTINUE
C     
 9    CONTINUE
C     
      RETURN
      END
