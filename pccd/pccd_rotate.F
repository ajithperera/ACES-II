      SUBROUTINE PCCD_ROTATE(U,ICORE,MAXCOR,ISPIN,NBAS,ICYCLE)
      IMPLICIT INTEGER (A-Z)
      PARAMETER(MINISIZ=50)
      LOGICAL PRINT,PRINT2
      DOUBLE PRECISION ONE,ONEM,ZILCH,SNRM2,X,ERROR,SDOT,DIISCOEF
      CHARACTER*1 ISP(2)
      DIMENSION ICORE(MAXCOR),error(minisiz,minisiz),DIISCOEF(MINISIZ)
      DIMENSION U(NBAS,NBAS)
      DIMENSION X(2)
      COMMON /SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF(4)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON/EXTRAPO/DIIS

      DATA ONE  /1.0D0/
      DATA ZILCH/0.0D0/
      DATA ONEM /-1.0D0/
      DATA ISP /'A','B'/
      PRINT = IFLAGS(1) .GE. 10
      PRINT2= IFLAGS(1) .GE. 20
      IONE=1
      NBAS2=NBAS*NBAS
CSSSS      DO 10 ISPIN=1,1+IUHF
       NOCC=NOCCO(ISPIN)
       NVRT=NVRTO(ISPIN)
       IMOOAO=1
       IT1PACKED=IMOOAO+NBAS2*IINTFP
       IT1FULL=IT1PACKED+NBAS2*IINTFP
       IR1MO=IT1FULL+NBAS2*IINTFP
       IR1OAO=IR1MO+NBAS2*IINTFP
       IJUNK1=IR1OAO+NBAS2*IINTFP
       IJUNK2=IJUNK1+NBAS2*IINTFP
       IJUNK3=IJUNK2+NBAS2*IINTFP
       IMOOAONEW=IJUNK3+NBAS2*IINTFP
       IAOOVRLAP=IMOOAONEW+NBAS2*IINTFP
       IMOAO=IAOOVRLAP+NBAS2*IINTFP
       ISHALF=IMOAO+NBAS2*IINTFP
       ISMHALF=ISHALF+NBAS2*IINTFP
       IMOAONEW=ISMHALF+NBAS2*IINTFP
       IERROR=IMOAONEW+NBAS2*IINTFP
       IDIISMAT=IERROR+NBAS2*MINISIZ*IINTFP
       iMASTERR1=IDIISMAT+2*MINISIZ*MINISIZ*IINTFP
       CALL GETREC(20,'JOBARC','SMHALF  ',NBAS2*IINTFP,
     &             ICORE(ISMHALF))
       CALL PUTREC(20,'JOBARC','SMHALFX ',NBAS2*IINTFP,
     &             ICORE(ISMHALF))
       CALL MINV(ICORE(ISMHALF),NBAS,NBAS,ICORE(IJUNK1),X,ZILCH,0,0)
       CALL PUTREC(20,'JOBARC','SHALF   ',NBAS2*IINTFP,
     &             ICORE(ISMHALF))
       CALL GETREC(20,'JOBARC','SCFEVEC'//ISP(ISPIN),NBAS2*IINTFP,
     &             ICORE(IJUNK1))
       CALL GETREC(-1,'JOBARC','BMBPT2IT',IONE,ICYCLE)
       if(icycle.ne.0)then
       CALL GETREC(20,'JOBARC','SAVEVEC'//ISP(ISPIN),NBAS2*IINTFP,
     &             ICORE(IJUNK1))
       endif
       CALL XGEMM('T','N',NBAS,NBAS,NBAS,ONE,ICORE(ISMHALF),NBAS,
     &            ICORE(IJUNK1),NBAS,ZILCH,ICORE(IMOOAO),NBAS)
       CALL PUTREC(20,'JOBARC','EVECOAO'//ISP(ISPIN),NBAS2*IINTFP,
     &             ICORE(IMOOAO))
       IF(ICYCLE.EQ.0)THEN
        CALL PUTREC(20,'JOBARC','EMATRIX'//ISP(ISPIN),
     &              NBAS2*MINISIZ*IINTFP,ICORE)
        CALL PUTREC(20,'JOBARC','MASTER1'//ISP(ISPIN),
     &              NBAS2*MINISIZ*IINTFP,
     &              ICORE(ierror))
        CALL PUTREC(20,'JOBARC','EVCOAOX'//ISP(ISPIN),NBAS2*IINTFP,
     &              ICORE(IMOOAO))
        CALL SCOPY(NBAS2,ICORE(IMOOAO),1,ICORE(IJUNK1),1)
        CALL MINV(ICORE(IJUNK1),NBAS,NBAS,ICORE(IJUNK2),X,ZILCH,0,0)
        CALL PUTREC(20,'JOBARC','EVCOAOI'//ISP(ISPIN),NBAS2*IINTFP,
     &              ICORE(IJUNK1))
       ENDIF
       ICYCLE=ICYCLE+1
       CALL MO2OAO(U,ICORE(IR1OAO),ICORE(IJUNK1),ICORE(IJUNK2),
     &             NBAS,ISPIN)

#ifdef _DEBUG_LVL0
       IF(PRINT)THEN
        write(6,*)' OAO rotation matrix determined by current gradient'
        CALL OUTPUT(ICORE(IR1OAO),1,NBAS,1,NBAS,NBAS,NBAS,1)
       ENDIF
#endif 
       CALL XGEMM('N','N',NBAS,NBAS,NBAS,ONE,ICORE(IR1OAO),
     &            NBAS,ICORE(IMOOAO),NBAS,ZILCH,ICORE(IMOOAONEW),NBAS)

#ifdef _DEBUG_LVL0
       IF(PRINT)THEN
        WRITE(6,*) 'New OAO orbitals:'
        CALL OUTPUT(ICORE(IMOOAONEW),1,NBAS,1,NBAS,NBAS,NBAS,1)
       ENDIF
#endif 

       CALL GETREC(20,'JOBARC','SMHALF  ',NBAS2*IINTFP,ICORE(ISMHALF))
       CALL XGEMM ('N','N',NBAS,NBAS,NBAS,ONE,ICORE(ISMHALF),NBAS,
     &             ICORE(IMOOAONEW),NBAS,ZILCH,ICORE(IMOAONEW),NBAS)
       CALL ORTHOGMO(ICORE(IMOAONEW),ICORE(IJUNK1),ICORE(IAOOVRLAP),
     &               ICORE(IJUNK2),ICORE(IJUNK3),ICORE(ISMHALF),NBAS,
     &               ISPIN)
       CALL GETREC(20,'JOBARC','SHALF  ',IINTFP*NBAS2,ICORE(ISMHALF))
       CALL XGEMM('T','N',NBAS,NBAS,NBAS,ONE,ICORE(ISMHALF),NBAS,
     &            ICORE(IJUNK1),NBAS,ZILCH,ICORE(IJUNK2),NBAS)

#ifdef _DEBUG_LVL0
       IF(PRINT)THEN
        Write(6,*)
        WRITE(6,*)' Orthogonalized new vectors in OAO basis '
        CALL OUTPUT(ICORE(IJUNK2),1,NBAS,1,NBAS,NBAS,NBAS,1)
       ENDIF
#endif 
       IF (DIIS) THEN
       CALL GETREC(20,'JOBARC','EVCOAOI'//ISP(ISPIN),NBAS2*IINTFP,
     &              ICORE(IMOOAO))
       CALL XGEMM ('N','N',NBAS,NBAS,NBAS,ONE,ICORE(IJUNK2),NBAS,
     &             ICORE(IMOOAO),NBAS,ZILCH,ICORE(IJUNK3),NBAS)
      CALL GETREC(20,'JOBARC','MASTER1'//ISP(ISPIN),
     &             IINTFP*NBAS2*(ICYCLE-1),ICORE(IERROR))
       IOFF=IERROR+IINTFP*NBAS2*(ICYCLE-1)
       CALL SCOPY(NBAS2,ICORE(IJUNK3),1,ICORE(IOFF),1) 
       CALL PUTREC(20,'JOBARC','MASTER1'//ISP(ISPIN), 
     &             IINTFP*NBAS2*ICYCLE,ICORE(IERROR))
       CALL GETREC(20,'JOBARC','EMATRIX'//ISP(ISPIN),
     &             IINTFP*NBAS2*(ICYCLE-1),ICORE(IERROR))
       IOFF=IERROR+IINTFP*NBAS2*(ICYCLE-1)
       CALL SCOPY(NBAS2,ICORE(IR1OAO),1,ICORE(IOFF),1)
       CALL PUTREC(20,'JOBARC','EMATRIX'//ISP(ISPIN),
     &             IINTFP*NBAS2*ICYCLE,ICORE(IERROR))
      ELSE

       CALL PUTREC(20,'JOBARC','SCFEVC'//ISP(ISPIN)//'0',
     &             NBAS2*IINTFP,ICORE(IJUNK2))
       IF(IUHF.EQ.0)THEN
        CALL PUTREC(20,'JOBARC','SCFEVCB0',
     &             NBAS2*IINTFP,ICORE(IJUNK2))
       ENDIF 

      ENDIF

CSSS10    CONTINUE

#ifdef _NOSKIP

       call zero(error,icycle*icycle)
       do 11 ispin=1,1+iuhf
       call getrec(20,'JOBARC','EMATRIX'//ISP(ISPIN),
     &             iintfp*nbas2*icycle,
     &             icore(ierror))
       do 30 i=1,icycle
        do 31 j=1,i
         ioffi=ierror+iintfp*nbas2*(i-1)
         ioffj=ierror+iintfp*nbas2*(j-1)
         error(i,j)=sdot(nbas2,icore(ioffi),1,icore(ioffj),1)-
     &              dfloat(nbas)+error(i,j)
         error(j,i)=error(i,j)
31      continue
30     continue
11     continue

#ifdef _DEBUG_LVL0
       IF(PRINT2)THEN
        write(6,*)' current error matrix '
        CALL OUTPUT(error,1,minisiz,1,minisiz,minisiz,minisiz,1)
       ENDIF
#endif
       call pccd_dodiis(error,icore(idiismat),icore(ijunk3),minisiz,
     &                  icycle)

#ifdef _DEBUG_LVL0
       IF(PRINT)THEN
        write(6,*)
        write(6,*)' current DIIS coefficients '
        call prvecr(icore(ijunk3),icycle)
       ENDIF
#endif 
       call scopy(icycle,icore(ijunk3),1,diiscoef,1)
       do 12 ispin=1,1+iuhf
       call getrec(20,'JOBARC','MASTER1'//ISP(ISPIN),
     &              iintfp*nbas2*icycle,
     &             icore(imasterr1))
       CALL ZERO(ICORE(IJUNK2),NBAS2)
       DO 100 I=1,ICYCLE
        IOFFR1=IMASTERR1+IINTFP*NBAS2*(I-1)
        CALL SAXPY(NBAS2,DIISCOEF(I),ICORE(IOFFR1),1,ICORE(IJUNK2),1)
100    CONTINUE

#ifdef _DEBUG_LVL0
       IF(PRINT)THEN
        write(6,*)
        write(6,*)' extrapolated R1 matrix '
        CALL OUTPUT(ICORE(Ijunk2),1,NBAS,1,NBAS,NBAS,NBAS,1)
       ENDIF
#endif 

       call getrec(20,'JOBARC','EVCOAOX'//ISP(ispin),nbas2*iintfp,
     &             icore(ismhalf))       
       call xgemm('n','n',nbas,nbas,nbas,one,icore(ijunk2),nbas,
     &            icore(ismhalf),nbas,zilch,icore(ijunk3),nbas)

#ifdef _DEBUG_LVL0
       IF(PRINT)THEN
        write(6,*)
        write(6,*)' extrapolated vectors in original OAO basis '
        CALL OUTPUT(ICORE(IJUNK3),1,NBAS,1,NBAS,NBAS,NBAS,1)
       ENDIF
#endif 

       call scopy(nbas2,icore(ijunk3),1,icore(imoaonew),1)
       CALL ORTHOGMO(ICORE(IMOAONEW),ICORE(IJUNK1),ICORE(IAOOVRLAP),
     &               ICORE(IJUNK2),ICORE(IJUNK3),ICORE(ISMHALF),NBAS,
     &               ISPIN)
       CALL SCOPY(NBAS2,ICORE(IJUNK1),1,ICORE(IMOAONEW),1)

#ifdef _DEBUG_LVL0
       IF(PRINT)THEN
        write(6,*)
        write(6,*)' extrapolated vectors in original AO basis '
        CALL OUTPUT(ICORE(Imoaonew),1,NBAS,1,NBAS,NBAS,NBAS,1)
       ENDIF
#endif 

       call scopy(nbas2,icore(imoaonew),1,icore(ijunk1)
     &    ,1)
9999   continue
       CALL PUTREC(20,'JOBARC','SCFEVC'//ISP(ISPIN)//'0',
     &             NBAS2*IINTFP,ICORE(IJUNK1))
       CALL PUTREC(20,'JOBARC','SAVEVEC'//ISP(ISPIN),
     &             NBAS2*IINTFP,ICORE(IJUNK1))
       IF(IUHF.EQ.0)THEN
        CALL PUTREC(20,'JOBARC','SCFEVCB0',
     &             NBAS2*IINTFP,ICORE(IJUNK1))
       CALL PUTREC(20,'JOBARC','SAVEVECB',
     &             NBAS2*IINTFP,ICORE(IJUNK1))
       ENDIF
12     continue     
      CALL GETREC(20,'JOBARC','OCCUPYA ',NIRREP,ICORE)
      CALL PUTREC(20,'JOBARC','OCCUPYA0',NIRREP,ICORE)
      IF(IUHF.NE.0)THEN
       CALL GETREC(20,'JOBARC','OCCUPYB ',NIRREP,ICORE)
       CALL PUTREC(20,'JOBARC','OCCUPYB0',NIRREP,ICORE)
      ENDIF
      CALL PUTREC(20,'JOBARC','BMBPT2IT',IONE,ICYCLE)
#endif 

      RETURN
      END
