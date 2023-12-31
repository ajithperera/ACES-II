      SUBROUTINE EDENS_4GTDENS(IRREPX,SCR,MAXCOR,IUHF,ISIDE,FACT2,
     +                         R0,IRPUT,NSIZEC)
C
C DRIVER FOR FORMATION OF EXCITED STATE DENSITY MATRICES
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD
      LOGICAL PRINT,PRINT2,VPROP,LAMTHERE,MBPT2,CC,CCD,RCCD,DRCCD
      LOGICAL LCCD,LCCSD,CC2
      CHARACTER*8 LABEL(3),LABEL2(6),LABSCF(2)
      DIMENSION SCR(MAXCOR),LENVV(2),LENOO(2),LENVO(2),MAP(6)
      DIMENSION SECMOM(6),GSECMOM(6),RSECMOM(6),TM(3),GM(3),RM(3)
      COMMON/FLAGS/IFLAGS(100)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/REFTYPE/MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/LAMSTATE/LAMTHERE
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &             IRPDPDAO(8),IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
      COMMON /VDINTLEN/ LENS(8),LENT(8)
      COMMON /VDINTPRT/NTPERT,NPERT(8),KPERT(8),IDIPSM(3),
     &                 IYZPERT,IXZPERT,IXYPERT,ITRANSX,
     &                 ITRANSY,ITRANSZ,NUCIND
      COMMON /INTPROG/ VPROP
C
      DATA MAP /1,4,2,5,6,3/
      DATA ONE,ONEM,ZILCH,HALF /1.0D0,-1.0D0,0.0D0,0.5D0/
      DATA LABEL /'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/
      DATA LABSCF/'SCFDENSA','SCFDENSB'/
      DATA LABEL2/'2NDMO_XX','2NDMO_YY','2NDMO_ZZ',
     &            '2NDMO_XZ','2NDMO_YZ','2NDMO_XZ'/
C
      NNP1O2(I)=(I*(I+1))/2
C
      FACT = FACT2
      PRINT=IFLAGS(1).GE.2
      PRINT2=IFLAGS(1).GE.20
C
      CALL GETREC(20,'JOBARC','NUCDIPOL',3*IINTFP,TM)
      CALL GETREC(20,'JOBARC','NUCDIPOL',3*IINTFP,RM)
      CALL GETREC(20,'JOBARC','NUCDIPOL',3*IINTFP,GM)
      CALL ZERO(SECMOM,6)
      CALL ZERO(GSECMOM,6)
      CALL ZERO(RSECMOM,6)
C
C DO SOME ACCOUNTING AND SET UP ADDRESSES
C
      CALL IZERO(LENVV,2)
      CALL IZERO(LENOO,2)
      CALL IZERO(LENVO,2)
      DO 10 ISPIN=1,1+IUHF
       LENVV(ISPIN)=IRPDPD(1,18+ISPIN)
       LENOO(ISPIN)=IRPDPD(1,20+ISPIN)
       LENVO(ISPIN)=IRPDPD(1, 8+ISPIN)
10    CONTINUE
C
      IDOO=1
      IDVV=IDOO+LENOO(1)+IUHF*LENOO(2)
      IDVO=IDVV+LENVV(1)+IUHF*LENVV(2)
      IDOV=IDVO+LENVO(1)+IUHF*LENVO(2)
      IDOOG=IDOV+LENVO(1)+IUHF*LENVO(2)
      IDVVG=IDOOG+LENOO(1)+IUHF*LENOO(2)
      IDVOG=IDVVG+LENVV(1)+IUHF*LENVV(2)
      IDOVG=IDVOG+LENVO(1)+IUHF*LENVO(2)
      IDOOR=IDOVG+LENVO(1)+IUHF*LENVO(2)
      IDVVR=IDOOR+LENOO(1)+IUHF*LENOO(2)
      IDVOR=IDVVR+LENVV(1)+IUHF*LENVV(2)
      IDOVR=IDVOR+LENVO(1)+IUHF*LENVO(2)
      ITOP=IDOVR+LENVO(1)+IUHF*LENVO(2)
      MXCOR=MAXCOR-ITOP+1
      IONE=1
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NAO)
      NMO=NOCCO(1)+NVRTO(1)
      IDFLMO=ITOP
      IDFLMOG=ITOP+NAO*NAO
      IDFLMOR=IDFLMOG+NAO*NAO
      IDPROP=IDFLMOR+(IUHF+1)*NAO*NAO
      ITMP1 =IDPROP+NAO*NAO
      ITMP2 =ITMP1 +NAO*NAO
      ITMP3 =ITMP2 +NAO*MAX(NAO,6)
      ITMP4 =ITMP3 +NAO*MAX(NAO,6)
      ITMP5 =ITMP4 +NAO*MAX(NAO,6)
C
C PUT RIGHT-HAND VECTOR ON LISTS 461-463 AND 3,490 - 4,490
C PUT LEFT-HAND VECTOR ON LISTS 444-446 AND 1,490 - 2,490
C
C     PREPARE R AND L VECTORS HERE    !!!!!!!!!
C
      CALL GETLST(SCR(ITOP),2,1,1,1,472)
#ifdef _DEBUG_LVL0
         Write(6,"(a,1x,i2)")"Writing the right vector of state :",irput
         Call checksum("list-472      :",scr(itop),nsizec,s)
         Write(6,*)
#endif
      CALL PUTLST(SCR(ITOP),IRPUT,1,1,3,472)
      CALL UPDATES(IRREPX,SCR(ITOP),461,2,490,IUHF)
      CALL GETLST(SCR(ITOP),3,1,1,1,472)
CSS      Write(6,"(a,1x,i2)")"Writing the left vector of state :",irput
CSS      Call checksum("list-472      :",scr(itop),nsizec,s)
      CALL UPDATES(IRREPX,SCR(ITOP),444,0,490,IUHF)
C
C PUT RESORTED RIGHT AND LEFT VECTORS ON DISK
C
       CALL RESORT(SCR(ITOP),IINTFP*MAXCOR,IUHF,IRREPX,444,434)
       CALL RESORT(SCR(ITOP),IINTFP*MAXCOR,IUHF,IRREPX,461,454)
C
      CALL WFANAL(SCR(ITOP),MAXCOR,IUHF,44,461,444,34,454,434,
     &            90,490,490,0,2,0,IRREPX,R0,ZILCH)
C
C
C CALCULATE THE EXCITED STATE DENSITY MATRIX
C
      LSTR1=490
      LSTL1=490
      LSTR1OFF=2
      LSTL1OFF=0
      LISTL2=444
      LISTR2=461
      LSTL2RS=434
      LSTR2RS=454
      LSTGRL =0
      LSTGTL =400
      LSTGRLOF=2
      LSTGTLOF=0
      LSTTMP  =90
      LSTTMPOF=2
      LSTT1   =90
      LSTT1OFF=0
      IF(MBPT2 .OR. LCCD .OR. LCCSD .OR. RCCD .OR. DRCCD .OR.CC2)
     &   R0=0.0d0
      CALL GDENS(IRREPX,IRREPX,1,SCR(IDOO),SCR(IDVV),SCR(IDVO),
     &           SCR(IDOV),SCR(ITOP),MXCOR,IUHF,FACT,R0,ZILCH,ONE,
     &           LSTGRL,LSTGTL,LSTGRLOF,LSTGTLOF,LSTTMP,LSTTMPOF,
     &           LSTR1,LSTL1,LSTR1OFF,LSTL1OFF,LISTR2,LISTL2,
     &           LSTR2RS,LSTL2RS,LSTT1,LSTT1OFF)
C
C EVALUATE GROUND STATE REDUCED DENSITY MATRIX (IF POSSIBLE)
C
      IF(LAMTHERE)THEN
       CALL GSDEN(SCR(IDOOG),SCR(IDVVG),SCR(IDVOG),SCR(IDOVG),
     &            0,190,0,90,144,44,ONE,SCR(ITOP),MXCOR,IUHF)
      ENDIF
C
C CALCULATE ONE-ELECTRON ENERGY
C
      IF (PRINT) THEN
        IF(CC)THEN
          CALL CALCONEL(SCR(IDOO),SCR(IDVV),SCR(IDVO),SCR(IDOV),
     &       SCR(IDFLMO),IUHF,ENONE)
        ENDIF
      ENDIF
C
C IF WE ARE CALCULATING THE T-AMPLITUDE RESPONSE, WE MUST EVALUATE
C THE ZETA AMPLITUDES.  DO THAT NOW. 
C
      IF(IFLAGS(91).GE.2)THEN
C
C PUT RESORTED RIGHT VECTOR ON LISTS 454-459
C
       CALL CALCXI(SCR(ITOP),IINTFP*MXCOR,IUHF)
C
C                         + 
C NOW CALCULATE <0|ZETA {p q exp(T)}|0> CONTRIBUTION TO THE DENSITY.
C
       IF(CC)THEN
        CALL GSDEN(SCR(IDOOR),SCR(IDVVR),SCR(IDVOR),SCR(IDOVR),
     &             0,190,0,90,144,44,ZILCH,SCR(ITOP),MXCOR,IUHF)
       ELSEIF(MBPT2 .OR. LCCD .OR. LCCSD .OR. RCCD .OR. DRCCD)THEN
        CALL GSDEN(SCR(IDOOR),SCR(IDVVR),SCR(IDVOR),SCR(IDOVR),
     &             0,190,0,90,61,44,ZILCH,SCR(ITOP),MXCOR,IUHF)
       ENDIF
       CALL SAXPY(LENOO(1)+IUHF*LENOO(2),ONE,SCR(IDOO),1,
     &            SCR(IDOOR),1)
       CALL SAXPY(LENVV(1)+IUHF*LENVV(2),ONE,SCR(IDVV),1,
     &            SCR(IDVVR),1)
       CALL SAXPY(LENVO(1)+IUHF*LENVO(2),ONE,SCR(IDVO),1,
     &            SCR(IDVOR),1)
       CALL SAXPY(LENVO(1)+IUHF*LENVO(2),ONE,SCR(IDOV),1,
     &            SCR(IDOVR),1)
      ENDIF
      IF (PRINT) THEN
      CALL CALCONEL(SCR(IDOOr),SCR(IDVVr),SCR(IDVOr),SCR(IDOVr),
     &     SCR(IDFLMO),IUHF,ENONE)
      ENDIF
C 
      FACT=2.0D0-DFLOAT(IUHF)
      DO 100 ISPIN=1,1+IUHF
       IOFFOO=IDOO+(ISPIN-1)*LENOO(1)
       IOFFVV=IDVV+(ISPIN-1)*LENVV(1)
       IOFFVO=IDVO+(ISPIN-1)*LENVO(1)
       IOFFOV=IDOV+(ISPIN-1)*LENVO(1)
       IOFFOOG=IDOOG+(ISPIN-1)*LENOO(1)
       IOFFVVG=IDVVG+(ISPIN-1)*LENVV(1)
       IOFFVOG=IDVOG+(ISPIN-1)*LENVO(1)
       IOFFOVG=IDOVG+(ISPIN-1)*LENVO(1)
       IOFFOOR=IDOOR+(ISPIN-1)*LENOO(1)
       IOFFVVR=IDVVR+(ISPIN-1)*LENVV(1)
       IOFFVOR=IDVOR+(ISPIN-1)*LENVO(1)
       IOFFOVR=IDOVR+(ISPIN-1)*LENVO(1)
C
C CALCULATE FULL UNPACKED EXCITED STATE DENSITY IN AO BASIS
C
       CALL EXPDEN(SCR(IOFFOO),SCR(IOFFVV),SCR(IOFFVO),SCR(IOFFOV),
     &             SCR(IDFLMO),NMO,1,ISPIN,.FALSE.)
C
C PUT REDUCED DENSITY ON DISK
C
       CALL WRITDG(SCR(IOFFOO),SCR(IOFFVV),SCR(IOFFVO),SCR(IOFFOV),
     &             LENOO(ISPIN),LENVV(ISPIN),LENVO(ISPIN),ISPIN,2)
       IF(PRINT2)THEN
        WRITE(6,5100)
5100    FORMAT(T3,'@EDENS-I, Excited state density matrix in ',
     &            'MO basis: ')
        CALL SSCAL(NMO*NMO,FACT,SCR(IDFLMO),1)
        CALL OUTPUT(SCR(IDFLMO),1,NMO,1,NMO,NMO,NMO,1)
        CALL SSCAL(NMO*NMO,ONE/FACT,SCR(IDFLMO),1)
       ENDIF
       CALL SYMMET2(SCR(IDFLMO),NMO)
       CALL NATORBS(SCR(IDFLMO),SCR(ITMP1),SCR(ITMP2),NMO,IUHF,ISPIN)
       CALL MO2AO3(SCR(IDFLMO),SCR(ITMP1),SCR(ITMP2),SCR(ITMP3),
     &             NAO,NMO,ISPIN)
       CALL GETREC(20,'JOBARC',LABSCF(ISPIN),NAO*NAO*IINTFP,SCR(ITMP3))
       CALL SAXPY (NAO*NAO,ONE/FACT,SCR(ITMP3),1,SCR(ITMP1),1) 
C
CSSS       call output(SCR(ITMP1), 1,nao,1,nao,nao,nao,1)
C CALCULATE FULL UNPACKED GROUND STATE DENSITY IN AO BASIS
C
       IF(LAMTHERE)THEN
        CALL EXPDEN(SCR(IOFFOOG),SCR(IOFFVVG),SCR(IOFFVOG),SCR(IOFFOVG),
     &              SCR(IDFLMOG),NMO,1,ISPIN,.FALSE.)
        IF(PRINT2)THEN
         WRITE(6,5101)
5101     FORMAT(T3,'@EDENS-I, Ground state density matrix in ',
     &             'MO basis: ')
         CALL SSCAL(NMO*NMO,FACT,SCR(IDFLMOG),1)
         CALL OUTPUT(SCR(IDFLMOG),1,NMO,1,NMO,NMO,NMO,1)
         CALL SSCAL(NMO*NMO,ONE/FACT,SCR(IDFLMOG),1)
        ENDIF
        CALL SYMMET2(SCR(IDFLMOG),NMO)
        CALL MO2AO3(SCR(IDFLMOG),SCR(ITMP2),SCR(ITMP3),SCR(ITMP4),
     &              NAO,NMO,ISPIN)
        CALL GETREC(20,'JOBARC',LABSCF(ISPIN),NAO*NAO*IINTFP,SCR(ITMP3))
        CALL SAXPY (NAO*NAO,ONE/FACT,SCR(ITMP3),1,SCR(ITMP2),1) 
       ENDIF
C
C CALCULATE FULL UNPACKED RESPONSE EXCITED STATE DENSITY IN AO BASIS
C
       IF(IFLAGS(91).GE.2)THEN
        CALL EXPDEN(SCR(IOFFOOR),SCR(IOFFVVR),SCR(IOFFVOR),SCR(IOFFOVR),
     &             SCR(IDFLMOR),NMO,1,ISPIN,.FALSE.)
C
C PUT RESPONSE DENSITY ON DISK
C
        CALL WRITDG(SCR(IOFFOOR),SCR(IOFFVVR),SCR(IOFFVOR),SCR(IOFFOVR),
     &              LENOO(ISPIN),LENVV(ISPIN),LENVO(ISPIN),ISPIN,1)
        IF(PRINT2)THEN
         WRITE(6,5102)
5102     FORMAT(T3,'@EDENS-I, Response excited state density matrix ',
     &             'in MO basis: ')
         CALL SSCAL(NMO*NMO,FACT,SCR(IDFLMOR),1)
         CALL OUTPUT(SCR(IDFLMOR),1,NMO,1,NMO,NMO,NMO,1)
         CALL SSCAL(NMO*NMO,ONE/FACT,SCR(IDFLMOR),1)
        ENDIF
        CALL SYMMET2(SCR(IDFLMOR),NMO)
        CALL MO2AO3(SCR(IDFLMOR),SCR(ITMP3),SCR(ITMP4),SCR(ITMP5),
     &              NAO,NMO,ISPIN)
        CALL GETREC(20,'JOBARC',LABSCF(ISPIN),NAO*NAO*IINTFP,SCR(ITMP4))
        CALL SAXPY (NAO*NAO,ONE/FACT,SCR(ITMP4),1,SCR(ITMP3),1) 
CSS        call output(SCR(ITMP3), 1,nao,1,nao,nao,nao,1)
       ENDIF
C
C READ IN THE DIPOLE INTEGRALS AND CALCULATE DIPOLE MOMENTS:
C
       DO 105 IXYZ=1,3
        IF(VPROP)THEN
         LENGTH=NNP1O2(NAO)
         CALL GETREC(20,'JOBARC',LABEL(IXYZ),LENGTH*IINTFP,
     &               SCR(ITMP4))
         CALL EXPND2(SCR(ITMP4),SCR(IDPROP),NAO)
        ELSE
         IRRPRT=IDIPSM(IXYZ)
         LENGTH=LENT(IRRPRT)
         CALL GETREC(20,'JOBARC',LABEL(IXYZ),LENGTH*IINTFP,
     &               SCR(ITMP4))
         CALL ZERO(SCR(IDPROP),NAO*NAO)
         CALL MATEXP(IRRPRT,IAOPOP,SCR(ITMP4),SCR(IDPROP))
         CALL MATEXP2(IRRPRT,SCR(IDPROP),SCR(ITMP4),NAO)
         CALL SCOPY  (NAO*NAO,SCR(ITMP4),1,SCR(IDPROP),1)
        ENDIF
        TM(IXYZ)=TM(IXYZ)+SDOT(NAO*NAO,SCR(IDPROP),1,SCR(ITMP1),1)*FACT
        GM(IXYZ)=GM(IXYZ)+SDOT(NAO*NAO,SCR(IDPROP),1,SCR(ITMP2),1)*FACT
        IF(IFLAGS(91).GE.2)THEN
         RM(IXYZ)=RM(IXYZ)+SDOT(NAO*NAO,SCR(IDPROP),1,SCR(ITMP3),1)*FACT
        ENDIF
105    CONTINUE
C
C READ IN THE SECOND MOMENT INTEGRALS AND CALCULATE DIFFERENCE IN
C EXTENT OF ELECTRONIC DISTRIBUTION 
C
       IXYZ=0
       DO 106 IXYZ1=1,3
        DO 107 IXYZ2=1,IXYZ1
         IXYZ=IXYZ+1
         IF(VPROP)THEN
          LENGTH=NNP1O2(NAO) 
          CALL GETREC(20,'JOBARC',LABEL2(MAP(IXYZ)),LENGTH*IINTFP,
     &                SCR(ITMP4))
          CALL EXPND2(SCR(ITMP4),SCR(IDPROP),NAO)
         ELSE
          IRRPRT=DIRPRD(IDIPSM(IXYZ1),IDIPSM(IXYZ2))
          LENGTH=LENT(IRRPRT)
          CALL GETREC(20,'JOBARC',LABEL2(MAP(IXYZ)),LENGTH*IINTFP,
     &                SCR(ITMP4))
          CALL VMINUS(SCR(ITMP4),LENGTH)
          CALL ZERO(SCR(IDPROP),NAO*NAO)
          CALL MATEXP(IRRPRT,IAOPOP,SCR(ITMP4),SCR(IDPROP))
          CALL MATEXP2(IRRPRT,SCR(IDPROP),SCR(ITMP4),NAO)
          CALL SCOPY  (NAO*NAO,SCR(ITMP4),1,SCR(IDPROP),1)
         ENDIF
         SECMOM(MAP(IXYZ))=SECMOM(MAP(IXYZ))+
     &           SDOT(NAO*NAO,SCR(IDPROP),1,SCR(ITMP1),1)*FACT
         GSECMOM(MAP(IXYZ))=GSECMOM(MAP(IXYZ))+
     &           SDOT(NAO*NAO,SCR(IDPROP),1,SCR(ITMP2),1)*FACT
         IF(IFLAGS(91).GE.2)THEN
          RSECMOM(MAP(IXYZ))=RSECMOM(MAP(IXYZ))+
     &            SDOT(NAO*NAO,SCR(IDPROP),1,SCR(ITMP3),1)*FACT
         ENDIF
107     CONTINUE
106    CONTINUE
C
C PROCESS EXCITED STATE DENSITY DIFFERENCE IF POSSIBLE
C
       IF(LAMTHERE)THEN
         IF (PRINT) THEN
           CALL SCOPY  (NMO*NMO,SCR(IDFLMO),1,SCR(ITMP1),1)
           CALL EDENPRC(SCR(ITMP1),SCR(ITMP2),SCR(ITMP4),NMO,NAO,
     &        IUHF,ISPIN,2) 
         ENDIF
        CALL SAXPY  (NMO*NMO,ONEM,SCR(IDFLMOG),1,SCR(IDFLMO),1)
        CALL EDENPRC(SCR(IDFLMO),SCR(ITMP2),SCR(ITMP4),NMO,NAO,
     &               IUHF,ISPIN,1) 
C
        CALL EXPDEN(SCR(IOFFOO),SCR(IOFFVV),SCR(IOFFVO),SCR(IOFFOV),
     &              SCR(IDFLMO),NMO,1,ISPIN,.FALSE.)
        QEL=FACT*HALF*SASUM(NMO,SCR(IDFLMO),NMO+1)
C
        CALL EXPDEN(SCR(IOFFOO),SCR(IOFFVV),SCR(IOFFVO),SCR(IOFFOV),
     &              SCR(IDFLMO),NMO,1,ISPIN,.TRUE.)
        CALL EXPDEN(SCR(IOFFOOG),SCR(IOFFVVG),SCR(IOFFVOG),SCR(IOFFOVG),
     &             SCR(IDFLMOG),NMO,1,ISPIN,.TRUE.)
        CALL SAXPY  (NMO*NMO,ONEM,SCR(IDFLMOG),1,SCR(IDFLMO),1)
        WRITE(6,5000)QEL
        IF(PRINT) CALL POPANAL(NAO,NMO,SCR(IDFLMO),SCR(ITMP1),
     &     SCR(ITMP2),MAXCOR-ITMP3,ISPIN)
       ENDIF
C
100   CONTINUE
      CALL PUTREC(20,'JOBARC','DENSDIFF',IINTFP*NMO*NMO,SCR(IDFLMO))
C
      IF(LAMTHERE)THEN
       WRITE(6,3000)
       WRITE(6,4000)
       WRITE(6,3000)
       WRITE(6,5001)GM
       WRITE(6,5002)GSECMOM,(GSECMOM(1)+GSECMOM(2)+GSECMOM(3))/3.0D0
      ENDIF
      WRITE(6,3000)
      WRITE(6,4001)
      WRITE(6,3000)
      WRITE(6,5001)TM
      WRITE(6,5002)SECMOM,(SECMOM(1)+SECMOM(2)+SECMOM(3))/3.0D0
      IF(IFLAGS(91).GE.2)THEN
       WRITE(6,3000)
       WRITE(6,4002)
       WRITE(6,3000)
       WRITE(6,5001)RM
       WRITE(6,5002)RSECMOM,(RSECMOM(1)+RSECMOM(2)+RSECMOM(3))/3.0D0
      ENDIF
      CALL GETLST(SCR(ITOP),2,1,1,1,472)
      CALL UPDATES(IRREPX,SCR(ITOP),461,2,490,IUHF)
      CALL GETLST(SCR(ITOP),3,1,1,1,472)
      CALL UPDATES(IRREPX,SCR(ITOP),444,0,490,IUHF)
      IF (PRINT) WRITE(6,6000)ENONE
C
C PUT RESORTED RIGHT AND LEFT VECTORS ON DISK
C
       CALL RESORT(SCR(ITOP),IINTFP*MAXCOR,IUHF,IRREPX,444,434)
       CALL RESORT(SCR(ITOP),IINTFP*MAXCOR,IUHF,IRREPX,461,454)
      RETURN
3000  FORMAT(72('-'))
4000  FORMAT(T25,' Ground state properties ')
4001  FORMAT(T25,' Excited state properties ')
4002  FORMAT(T25,' Response Excited state properties ')
5000  FORMAT(T20,' Approximate excitation level (AEL) : ',F10.5,'.')
5003  FORMAT(T3,'Quasiparticle excitation level (2)',F10.5,'.')
5001  FORMAT(T3,'Dipole moment ',3(F16.12,1X))
5002  FORMAT(T3,'Second moment of charge distribution :',
     &       /,T3,'XX : ',F12.6,' YY : ',F12.6,' ZZ : ',F12.6,
     &       /,T3,'XY : ',F12.6,' YZ : ',F12.6,' YZ : ',F12.6,
     &       /,T3,'Spherical average : ',F12.6)
6000  FORMAT(T3,'Total one-electron energy ',F15.10,' a.u.')
      END
