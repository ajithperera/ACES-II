      SUBROUTINE TDENS(IRREPX,SCR,MAXCOR,IUHF,ISIDE,FACT,R0,ROOT,
     &                 CORE_SEARCH,DOUBLE_EES)
C
C DRIVER FOR FORMATION OF TRANSITION DENSITY
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MAXCENT = 30)
C
C Hard coded limit on no. of centers. Only apply for NMR interpretative 
C tools and in no way effected to other aspects of the code. This hard
C coded limit (very reasonable) is imposed to avoid changes in several 
C different places in the code. 04/96 S. Ajith perera
C
      INTEGER POP,VRT
      LOGICAL PRINT,PRINT2,VPROP,CORE_SEARCH,DOUBLE_EES
      CHARACTER*8 LABEL(3),STRNG(2), STRNGNMR(2)
      DIMENSION SCR(MAXCOR),LENVV(2),LENOO(2),LENVO(2),TM(3),
     &          TMNMR (6*MAXCENT)
      COMMON/FLAGS/IFLAGS(100)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /VDINTLEN/ LENS(8),LENT(8)
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &             IRPDPDAO(8),IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
      COMMON /VDINTPRT/NTPERT,NPERT(8),KPERT(8),IDIPSM(3),
     &                 IYZPERT,IXZPERT,IXYPERT,ITRANSX,
     &                 ITRANSY,ITRANSZ,NUCIND
      COMMON /INTPROG/ VPROP
      COMMON /FILES/ LUOUT, MOINTS
C
      DATA ONE,ZILCH /1.0D0,0.0D0/
      DATA LABEL /'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/
      DATA STRNG /'TMRIGHT ','TMLEFT  '/
      DATA STRNGNMR /'NMRTMRHT','NMRTMLHT'/
C
      NNP1O2(I)=(I*(I+1))/2
C
      PRINT=IFLAGS(1).GE.2
      PRINT2=IFLAGS(1).GE.30
C
      CALL ZERO(TM,3)
      CALL ZERO (TMNMR, MAXCENT*6)
C
      CALL IZERO(LENVV,2)
      CALL IZERO(LENOO,2)
      CALL IZERO(LENVO,2)
      DO 10 ISPIN=1,1+IUHF
       LENVV(ISPIN)=IRPDPD(IRREPX,18+ISPIN)
       LENOO(ISPIN)=IRPDPD(IRREPX,20+ISPIN)
       LENVO(ISPIN)=IRPDPD(IRREPX, 8+ISPIN)
10    CONTINUE
C
      IDOO=1
      IDVV=IDOO+LENOO(1)+IUHF*LENOO(2)
      IDVO=IDVV+LENVV(1)+IUHF*LENVV(2)
      IDOV=IDVO+LENVO(1)+IUHF*LENVO(2)
      IDOOG=IDOV+LENVO(1)+IUHF*LENVO(2)
      IDVVG=IDOOG+NFMI(1)+IUHF*NFMI(2)
      IDVOG=IDVVG+NFEA(1)+IUHF*NFEA(2)
      IDOVG=IDVOG+NT(1)+IUHF*NT(2)
      ITOP =IDOVG+NT(1)+IUHF*NT(2)
      MXCOR=MAXCOR-ITOP+1
C
      IONE=1
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NAO)
      NMO=NOCCO(1)+NVRTO(1)
      IDFLMO=ITOP
      IDPROP=IDFLMO+(IUHF+1)*NAO*NAO
      ITMP1 =IDPROP+NAO*NAO
      ITMP2 =ITMP1 +NAO*NAO
      ITMP3 =ITMP2 +NAO*NAO
C
      MXCOR=MAXCOR-ITOP+1
C
      IF(ISIDE.EQ.2)THEN
C
C LEFT HAND DENSITY
C
       CALL RESORT  (SCR,MAXCOR,IUHF,IRREPX,444,434)
       lstr1=-1
       lstl1=490
       lstr1off=0
       lstl1off=0
       listl2=444
       listr2=-1
       lstl2rs=434
       lstr2rs=-1
       lstgrl =-1
       lstgtl =400
       lstgrlof=0
       lstgtlof=0
       lsttmp  =490
       lsttmpof=2
       lstt1 = 90
       lstt1off=0
       call gdens(irrepx,1,irrepx,scr(idoo),scr(idvv),scr(idvo),
     &            scr(idov),scr(itop),mxcor,iuhf,1.0d0,one,zilch,zilch,
     &            lstgrl,lstgtl,lstgrlof,lstgtlof,lsttmp,lsttmpof,
     &            lstr1,lstl1,lstr1off,lstl1off,listr2,listl2,
     &            lstr2rs,lstl2rs,lstt1, lstt1off)
c       CALL GFORMG  (IRREPX,1,444,44,400,SCR,MAXCOR,0,ONE,ONE,IUHF)
c       CALL ZERO(SCR,ITOP-1)
c       CALL TDENSLOO(IRREPX,SCR(IDOO),SCR(ITOP),MXCOR,IUHF)
c       CALL TDENSLVV(IRREPX,SCR(IDVV),SCR(ITOP),MXCOR,IUHF)
c       CALL TDENSLVO(IRREPX,SCR(IDVO),SCR(ITOP),MXCOR,IUHF)
c       CALL GETLST  (SCR(IDOV),1,1,1,1,490)
c       IF(IUHF.NE.0)THEN
c        CALL GETLST(SCR(IDOV+LENVO(1)),1,1,1,2,490)
c       ENDIF
      ELSE
C
C RIGHT HAND DENSITY
C
       CALL RESORT  (SCR,MAXCOR,IUHF,IRREPX,444,434)
       lstr1=490
       lstl1=190
       lstr1off=0
       lstl1off=0
       listl2=144
       listr2=444
       lstl2rs=134
       lstr2rs=434
       lstgrl =400
       lstgtl =0 
       lstgrlof=0
       lstgtlof=2
       lsttmp  =490
       lsttmpof=2
       lstt1 = 90
       lstt1off=0
       call gdens(1,irrepx,irrepx,scr(idoo),scr(idvv),scr(idvo),
     &            scr(idov),scr(itop),mxcor,iuhf,1.0d0,r0,one,zilch,
     &            lstgrl,lstgtl,lstgrlof,lstgtlof,lsttmp,lsttmpof,
     &            lstr1,lstl1,lstr1off,lstl1off,listr2,listl2,
     &            lstr2rs,lstl2rs,lstt1,lstt1off)
 
c       CALL GFORMG  (1,IRREPX,144,444,400,SCR,MAXCOR,0,ONE,ONE,IUHF)
c       CALL ZERO    (SCR,ITOP-1)
c       CALL TDENSROO(IRREPX,SCR(IDOO),SCR(IDOV),SCR(ITOP),MXCOR,IUHF)
c       CALL TDENSRVV(IRREPX,SCR(IDVV),SCR(IDOV),SCR(ITOP),MXCOR,IUHF)
c        CALL GSDEN(SCR(IDOOG),SCR(IDVVG),SCR(IDVOG),SCR(IDOVG),
c     &            0,190,0,90,144,44,ZILCH,SCR(ITOP),MXCOR,IUHF)
c       CALL TDENSRVO(IRREPX,SCR(IDVO),SCR(IDOV),SCR(IDOOG),SCR(IDVVG),
c     &               SCR(ITOP),MXCOR,IUHF)
C
C NOW ADD IN S0 TIMES THE GROUND STATE DENSITY
C
c       IF(IRREPX.EQ.1)THEN
c        CALL SAXPY (LENOO(1)+IUHF*LENOO(2),R0,SCR(IDOOG),1,SCR(IDOO),1)
c        CALL SAXPY (LENVV(1)+IUHF*LENVV(2),R0,SCR(IDVVG),1,SCR(IDVV),1)
c        CALL SAXPY (LENVO(1)+IUHF*LENVO(2),R0,SCR(IDVOG),1,SCR(IDVO),1)
c        CALL SAXPY (LENVO(1)+IUHF*LENVO(2),R0,SCR(IDOVG),1,SCR(IDOV),1)
c       ENDIF
C
      ENDIF
C      
      DO 100 ISPIN=1,1+IUHF
       IOFFOO=IDOO+(ISPIN-1)*LENOO(1)
       IOFFVV=IDVV+(ISPIN-1)*LENVV(1)
       IOFFVO=IDVO+(ISPIN-1)*LENVO(1)
       IOFFOV=IDOV+(ISPIN-1)*LENVO(1)
C
C CALCULATE FULL UNPACKED DENSITY MATRIX 
C
       CALL EXPDEN(SCR(IOFFOO),SCR(IOFFVV),SCR(IOFFVO),SCR(IOFFOV),
     &             SCR(IDFLMO),NMO,IRREPX,ISPIN,.FALSE.)
C
       IF (PRINT) THEN
         WRITE(6,5100)
 5100    FORMAT(T3,'@TDENS-I, Largest elements of transition density ')
         N=MIN(NMO*NMO,IFLAGS(14))
         CALL LARGE(N,NMO,SCR(IDFLMO),SCR(ITMP1),SCR(ITMP2)) 
       ENDIF
       IF(PRINT2)THEN
        WRITE(6,5101)
5101    FORMAT(T3,'@TDENS-I, Transition density matrix : ')
        CALL OUTPUT(SCR(IDFLMO),1,NMO,1,NMO,NMO,NMO,1)
       ENDIF
       CALL MO2AO3(SCR(IDFLMO),SCR(ITMP1),SCR(ITMP2),SCR(ITMP3),
     &             NAO,NMO,ISPIN)
C
C READ IN THE DIPOLE INTEGRALS AND TRANSFORM THEM TO THE MO BASIS
C          
       DO 105 IXYZ=1,3
          IF(VPROP)THEN
             LENGTH=NNP1O2(NAO)
         CALL GETREC(20,'JOBARC',LABEL(IXYZ),LENGTH*IINTFP,
     &               SCR(ITMP2))
         CALL EXPND2(SCR(ITMP2),SCR(IDPROP),NAO) 
        ELSE
         IRRPRT=IDIPSM(IXYZ)
         LENGTH=LENT(IRRPRT)
         CALL GETREC(20,'JOBARC',LABEL(IXYZ),LENGTH*IINTFP,
     &               SCR(ITMP2))
         CALL VMINUS(SCR(ITMP2),LENGTH)
         CALL ZERO(SCR(IDPROP),NAO*NAO)
         CALL MATEXP(IRRPRT,IAOPOP,SCR(ITMP2),SCR(IDPROP))
         CALL MATEXP2(IRRPRT,SCR(IDPROP),SCR(ITMP2),NAO)
         CALL SCOPY  (NAO*NAO,SCR(ITMP2),1,SCR(IDPROP),1)
        ENDIF
        TM(IXYZ)=TM(IXYZ)+SDOT(NAO*NAO,SCR(IDPROP),1,SCR(ITMP1),1)*FACT
105    CONTINUE
C Extend the transition moments beyound the dipole approximation. 
C Ajith Perera, 02/2020.

        IF (CORE_SEARCH .OR. DOUBLE_EES) THEN
           ITOP   = ITMP3 + NAO*NAO
           MAXCOR = MAXCOR-ITOP+1
           CALL BEYOND_DIPOLE(SCR(ITMP1),SCR(ITMP2),SCR(ITMP3),
     &                        TM,SCR,MAXCOR,NAO,ISIDE)
        ENDIF 
C
C Compute the effectve one-particle spin-orbit coupling contributions
C
       MXCOR = MAXCOR - ITMP2 + IONE
CSSS       Call Spin_orbit(Scr(Itmp1), Scr(Itmp3), Mxcor, Nao)

       IF (IFLAGS(1) .GE. 40) THEN
          CALL HEADER ('@-TDENS Transition density matrix', 0, 6) 
          CALL TAB (6, SCR(ITMP1), NAO, NAO, NAO, NAO)
       ENDIF
       ITHR=3*IINTFP
       CALL PUTREC(20,'JOBARC',STRNG(ISIDE),ITHR,TM)
C
C    &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C    &                                                            &
C    & In a wild chase to develop interpretative tools for NMR    &
C    & Marcel and I decided to take a look at excitation energies &
C    & (denominators in the sum-over-state expression) and        &
C    & the transition moments over FC, SD and PSO operators.      &
C    & (numerators in the sum-over-sate expression)               &
C    & Implemented by S. Ajith Perera 04/01/96                    &
C    &                                                            &
C    &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C Get some important information from the JOBARC file.
C      
       IF (IFLAGS(18) .EQ. 10 .OR. IFLAGS(18) .EQ. 9 .OR.
     &    IFLAGS(18) .EQ. 8) THEN
C
         CALL GETREC (20, 'JOBARC', 'NREALATM', IONE, NATOM)
C
         IF (IFLAGS(18) .EQ. 8) THEN
           NTPERT = 3*NATOM
         ELSE IF (IFLAGS(18) .EQ. 9) THEN
           NTPERT = NATOM
         ELSE IF (IFLAGS(18) .EQ. 10) THEN
           NTPERT = 6*NATOM
         ENDIF
C
         MXCOR = MAXCOR - ITMP2 + IONE
         CALL INTPRTNMR (SCR(ITMP1), SCR(ITMP2), TMNMR, MXCOR, NAO,
     &      NMO, ISIDE, ISPIN, NTPERT, MAXCENT, NATOM, 
     &      FACT)
C
       ENDIF

 100  CONTINUE
C
      IF (IFLAGS(18) .EQ. 10 .OR. IFLAGS(18) .EQ. 9 .OR.
     &   IFLAGS(18) .EQ. 8) THEN
C
        IF (IFLAGS(1) .GT. 40) THEN
C
          IF (ISIDE .EQ. 1) THEN
            WRITE (LUOUT, 150)
            WRITE (LUOUT, *)
            DO 60 IPERT = 1, NTPERT
              WRITE (LUOUT, 250)  TMNMR (IPERT)
 60         CONTINUE
C
          ELSE IF (ISIDE .EQ. 2) THEN
            WRITE (LUOUT, 200)
            WRITE (LUOUT, *)
            DO 70 IPERT = 1, NTPERT
              WRITE (LUOUT, 250)  TMNMR (IPERT)
 70         CONTINUE
C
          ENDIF
        ENDIF
C
 150    FORMAT(T9, 'RIGHT transition moment over NMR operators')
 200    FORMAT(T9, 'LEFT transition moment over NMR operators')
 250    FORMAT(T15, F12.8)
C
C Write left and right transition moments over NMR operators
C to the JOBARC file
C
        LENGTH = NTPERT*IINTFP
        CALL PUTREC (20, 'JOBARC', STRNGNMR(ISIDE), LENGTH, TMNMR)
C
        IF (ISIDE .EQ. 2) THEN
          ITOP = ITMP2
          ITOP1 = ITOP +  NTPERT*NTPERT
          ITOP2 = ITOP1 + NTPERT
          ITOP3 = ITOP2 + NTPERT
C
          CALL PRTDPSNMR (SCR(ITOP), SCR(ITOP1), SCR(ITOP2),
     &       ROOT, NTPERT, NATOM)
        ENDIF
      ENDIF
C
      RETURN
      END

