      Subroutine Pccd_rotate_aux(Icore,Maxcor,Nbas,Iuhf,icycle)

      Implicit Integer (A-Z)
      Double Precision One,Onem,Zilch,Snrm2,X,Error,Sdot,Diiscoef
      Logical Oldmos_exist 
      Character*1 Isp(2)

      Parameter(Minisiz=50)
      Dimension Error(Minisiz,Minisiz),Diiscoef(Minisiz)
      Dimension Icore(Maxcor)

#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "flags2.com"
 
      Data One,Zilch,Onem,Ione/1.0D0,0.0D0,-1.0D0,1/
      Data Isp /"A","B"/

      Nbas2 = Nbas*Nbas 

      Ijunk1    = Ione
      Ijunk2    = Ijunk1 + Nbas2*Iintfp
      Ijunk3    = Ijunk2 + Nbas2*Iintfp
      Ierror    = Ijunk3 + Nbas2 *Iintfp
      Idissmat  = Ierror + Nbas2*Iintfp
      Imasterr1 = Idissmat + Minsiz*Minsiz*Iintfp
      Ismhalf   = Imasterr1 + Nbas2*Iintfp
      Imoaonew  = Ismhalf  + Nbas2*Iintfp
      Iend      = Imoaonew + Nbas2*Iintfp
      If (Iend .Gt. Maxcor) Call Insmem("pccd_rotate_aux",Iend,Maxcor)

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
31    continue
30    continue
11    continue

#ifdef _DEBUG_LVL0
      write(6,*)' current error matrix '
      call output(error,1,Icycle,1,Icycle,minisiz,minisiz,1)
#endif
      call dodiis(error,icore(idiismat),icore(ijunk3),minisiz,
     &            icycle)

#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,*)' current DIIS coefficients '
      call prvecr(icore(ijunk3),icycle)
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
 100  CONTINUE

#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,*)' extrapolated R1 matrix '
      CALL OUTPUT(ICORE(Ijunk2),1,NBAS,1,NBAS,NBAS,NBAS,1)
#endif
      call getrec(20,'JOBARC','EVCOAOX'//ISP(ispin),nbas2*iintfp,
     &             icore(ismhalf))
      call xgemm('n','n',nbas,nbas,nbas,one,icore(ijunk2),nbas,
     &            icore(ismhalf),nbas,zilch,icore(ijunk3),nbas)
#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,*)' extrapolated vectors in original OAO basis '
      CALL OUTPUT(ICORE(IJUNK3),1,NBAS,1,NBAS,NBAS,NBAS,1)
#endif
      call scopy(nbas2,icore(ijunk3),1,icore(imoaonew),1)
      CALL ORTHOGMO(ICORE(IMOAONEW),ICORE(IJUNK1),ICORE(IAOOVRLAP),
     &               ICORE(IJUNK2),ICORE(IJUNK3),ICORE(ISMHALF),NBAS,
     &               ISPIN)
      CALL SCOPY(NBAS2,ICORE(IJUNK1),1,ICORE(IMOAONEW),1)

      Inquire(File="OLDMOS",Exist=Oldmos_exist)
      If (Oldmos_exist) Then
          Call runit("rm OLDMOS")
      Endif

       If (Icycle .Eq. 0) Then
           Direction = Onem
           Write(15,*) E_new, Direction
       Endif

#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,*)' extrapolated vectors in original AO basis '
      CALL OUTPUT(ICORE(Imoaonew),1,NBAS,1,NBAS,NBAS,NBAS,1)
#endif
      call scopy(nbas2,icore(imoaonew),1,icore(ijunk1),1)

9999   continue
      CALL PUTREC(20,'JOBARC','SCFEVC'//ISP(ISPIN)//'0',
     &            NBAS2*IINTFP,ICORE(IJUNK1))
      CALL PUTREC(20,'JOBARC','SAVEVEC'//ISP(ISPIN),
     &            NBAS2*IINTFP,ICORE(IJUNK1))
      IF(IUHF.EQ.0)THEN
      CALL PUTREC(20,'JOBARC','SCFEVCB0',
     &            NBAS2*IINTFP,ICORE(IJUNK1))
      CALL PUTREC(20,'JOBARC','SAVEVECB',
     &            NBAS2*IINTFP,ICORE(IJUNK1))
      ENDIF
12    continue

      CALL GETREC(20,'JOBARC','OCCUPYA ',NIRREP,ICORE)
      CALL PUTREC(20,'JOBARC','OCCUPYA0',NIRREP,ICORE)
      IF(IUHF.NE.0)THEN
       CALL GETREC(20,'JOBARC','OCCUPYB ',NIRREP,ICORE)
       CALL PUTREC(20,'JOBARC','OCCUPYB0',NIRREP,ICORE)
      ENDIF
      CALL PUTREC(20,'JOBARC','BMBPT2IT',IONE,ICYCLE)

      RETURN
      END



