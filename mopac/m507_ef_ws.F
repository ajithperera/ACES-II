      SUBROUTINE EF(XPARAM, NVAR, FUNCT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION LAMDA,LAMDA0
C
#include "mopcvar.par"
C
      DIMENSION XPARAM(MAXPAR)
**********************************************************************
*
*   EF IS A QUASI NEWTON RAPHSON OPTIMIZATION ROUTINE BASED ON
*      Jacs Simons P-RFO algorithm as implemented by Jon Baker
*      (J.COMP.CHEM. 7, 385). Step scaling to keep length within
*      trust radius is taken from Culot et al. (Theo. Chim. Acta 82, 189)
*      The trust radius can be updated dynamically according to Fletcher.
*      Safeguards on valid step for TS searches based on actual/predicted
*      function change and change in TS mode are own modifications
*
*  ON ENTRY XPARAM = VALUES OF PARAMETERS TO BE OPTIMISED.
*           NVAR   = NUMBER OF PARAMETERS TO BE OPTIMISED.
*
*  ON EXIT  XPARAM = OPTIMISED PARAMETERS.
*           FUNCT  = HEAT OF FORMATION IN KCAL/MOL.
*
*  Current version implementing combined NR, P-RFO and QA algorithm
*      together with thrust radius update and step rejection was
*      made october 1992 by F.Jensen, Odense, DK
*
**********************************************************************
C
      COMMON /MESAGE/ IFLEPO,ISCF
      COMMON /GEOVAR/ XDUMMY(MAXPAR), NDUM, LOC(2,MAXPAR)               IR0494
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /GEOSYM/ NDEP,LOCPAR(MAXPAR),IDEPFN(MAXPAR),LOCDEP(MAXPAR)
      COMMON /ISTOPE/ AMS(107)
      COMMON /LAST  / LAST
      COMMON /KEYWRD/ KEYWRD
      COMMON /TIMER / TIME0
      COMMON /GRADNT/ GRAD(MAXPAR),GNFINA
C Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
     3       /MOLKSR/ FRACT
      COMMON /NUMCAL/ NUMCAL
      COMMON /TIMDMP/ TLEFT, TDUMP
      COMMON /SIGMA2/ GNEXT1(MAXPAR), GMIN1(MAXPAR)
      COMMON /NLLCOM/ HESS(MAXPAR,MAXPAR),BMAT(MAXPAR,MAXPAR),
     1                PMAT(MAXPAR,MAXPAR)                                IR0494
      COMMON/OPTEF/OLDF(MAXPAR),D(MAXPAR),VMODE(MAXPAR),
     $U(MAXPAR,MAXPAR),DD,rmin,rmax,omin,xlamd,xlamd0,skal,
     $MODE,NSTEP,NEGREQ,IPRNT

      DIMENSION IPOW(9), EIGVAL(MAXPAR),TVEC(MAXPAR),SVEC(MAXPAR),
     1FX(MAXPAR),HESSC(MAXHES),UC(MAXPAR**2),oldfx(maxpar),
     1oldeig(maxpar),
     $oldhss(maxpar,maxpar),oldu(maxpar,maxpar),ooldf(maxpar)
      DIMENSION BB(MAXPAR,MAXPAR) 

      LOGICAL RESTRT,SCF1,LIMSCF
      LOGICAL LUPD,lts,lrjk,lorjk,rrscal,donr,gnmin                      IR0794
      CHARACTER*80 KEYWRD
      EQUIVALENCE(IPOW(1),IHESS)
      DATA  ICALCN,ZERO,ONE,TWO    /0,0.D0,1.D0,2.D0/
      DATA tmone /1.0d-1/, TMTWO/1.0D-2/, TMSIX/1.0D-06/
      data three/3.0d0/, four/4.0d0/, 
     1pt25/0.25d0/, pt5/0.50d0/, pt75/0.75d0/
      data demin/2.0d-2/, gmin/5.0d0/
C
C       Variables check                                                 IR0494
c
      if(nvar.gt.(3*NUMAT-6)) then
        igeok=index(keywrd,' GEO-OK')
        if ((igeok.eq.0).and.(nvar.ne.(3*NUMAT))) then
             write(6,'(/,5x,"EF: Redundant variables. ",
     1       "Eigenvector Following NOT recommended",/,5x,
     2       "Specify GEO-OK if you want to continue",/)')
             stop
        endif
      elseif(nvar.eq.(3*NUMAT)) then
          isXYZ=index(keywrd,' XYZ')
          if((isXYZ.eq.0).and.(igeok.eq.0)) then
           write(6,'(/,5x,"EF: Exactly 3*N variables.",/,5x,
     1      "Specif XYZ if you are using CARTESIAN coordinates",/,5x,
     2      "If you use INTERNAL coordinate,  you can specify GEO-OK,",/,
     3      5x,"but Eigenvector Following is NOT recommended",/)')
           stop
          endif
      endif
C
C      Zero out some variables (for MD-type runs)                       IR0894
C
      do 5 i=1,MAXPAR
  5      EIGVAL(i)=0.0d0
      do 6 i=1,MAXPAR
  6      FX(i)=0.0d0
      odd=0.0d0
      odmax=0.0d0
      oolde=0.0d0
C
C     GET ALL INITIALIZATION DATA
C
      IF(ICALCN.NE.NUMCAL) 
     1CALL EFSTR(XPARAM,FUNCT,IHESS,NTIME,ILOOP,IGTHES,MXSTEP,IRECLC,    IR0494
     *           IUPD,DMAX,DDMAX,dmin,TOL2,TOTIME,nvar,SCF1,LUPD,ldump,
     *           rrscal,donr,gnmin)
      lts=.false.
      if (negreq.eq.1) lts=.true.
      lorjk=.false.
c     osmin is smallest step for which a ts-mode overlap less than omin
c     will be rejected. for updated hessians there is little hope of
c     better overlap by reducing the step below 0.005. for exact hessian
c     the overlap should go toward one as the step become smaller, but
c     don't allow very small steps 
      osmin=0.005d0
      if(ireclc.eq.1)osmin=0.001d0
      IF (SCF1) THEN
         GNFINA=SQRT(DDOT(NVAR,GRAD,1,GRAD,1))                          IR0494
         IFLEPO=1
         RETURN
      ENDIF
C     CHECK THAT GEOMETRY IS NOT ALREADY OPTIMIZED
         RMX=SQRT(DDOT(NVAR,GRAD,1,GRAD,1))                             IR0494
         IF (RMX.LT.TOL2) THEN
            IFLEPO=2
            LAST=1
            RETURN
         ENDIF
C     GET INITIAL HESSIAN. IF ILOOP IS .LE.0 THIS IS AN OPTIMIZATION RESTART 
C     AND HESSIAN SHOULD ALREADY BE AVAILABLE                                   
      IF (ILOOP .GE. 0) CALL GETHES(XPARAM,IGTHES,NVAR,iloop,TOTIME)     IR0594

C     START OF MAIN LOOP
C     WE NOW HAVE GRADIENTS AND A HESSIAN. IF THIS IS THE FIRST
C     TIME THROUGH DON'T UPDATE THE HESSIAN. FOR LATER LOOPS ALSO
C     CHECK IF WE NEED TO RECALCULATE THE HESSIAN
      IFLEPO=0
      itime=0
   10 CONTINUE
c     store various things for possibly omin rejection
      do 30 i=1,nvar
         oldfx(i)=fx(i)
         ooldf(i)=oldf(i)
         oldeig(i)=eigval(i)
         do 20 j=1,nvar
         oldhss(i,j)=hess(i,j)
           oldu(i,j)=u(i,j)
20       continue
30    continue
      IF (IHESS.GE.IRECLC.AND.IFLEPO.NE.15) THEN
         ILOOP=1
         IHESS=0
         if (igthes.ne.3)IGTHES=1
         CALL GETHES(XPARAM,IGTHES,NVAR,iloop,TOTIME)
      ENDIF
         IF (IHESS.GT.0) CALL UPDHES(SVEC,TVEC,NVAR,IUPD)
         IF(IPRNT.GE.2) CALL GEOUT                                      IR0494
         IF(IPRNT.GE.2) THEN
            WRITE(6,'('' XPARAM '')')
            WRITE(6,'(5(2I3,F10.4))')(LOC(1,I),LOC(2,I),
     1                                XPARAM(I),I=1,NVAR)
            WRITE(6,'('' GRADIENTS'')')
            WRITE(6,'(3X,8F9.3)')(GRAD(I),I=1,NVAR)
         ENDIF
C
C        PRINT RESULTS IN CYCLE
         GNFINA=SQRT(DDOT(NVAR,GRAD,1,GRAD,1))                          IR0494
cmgc      TIME2=SECOND()
      if (itime.eq.0) time1=time0
      TSTEP=TIME2-TIME1
      IF (TSTEP.LT.ZERO)TSTEP=ZERO
      TLEFT=TOTIME-TIME2                                                IR0494
      TIME1=TIME2
      itime=itime+1
      WRITE(6,40) NSTEP+1, TSTEP, TLEFT, GNFINA, FUNCT
   40 FORMAT(' CYCLE:',I4,' TIME:',F7.2,' TIME LEFT:',F9.1,
     1       ' GRAD.:',F10.3,' HEAT:',G13.7)
      IF (TLEFT .LT. TSTEP*TWO) GOTO 280
      IHESS=IHESS+1
      NSTEP=NSTEP+1
C
C        TEST FOR CONVERGENCE
C
      RMX=SQRT(DDOT(NVAR,GRAD,1,GRAD,1))                                IR0494
      IF (RMX.LT.TOL2)GOTO 250
      OLDE  = FUNCT
      oldgn = rmx
      DO 60 I=1,NVAR
         OLDF(I)=GRAD(I)
60    CONTINUE
C
C     if the optimization is in cartesian coordinates, we should remove
C     translation and rotation modes. Possible problem if run is in
C     internal but with exactly 3*natoms variable (i.e. dummy atoms
C     are also optimized).
      if ( isXYZ.ne.0 ) then                                            IR0494
        write(6,'(10x,"Optiminzing using CARTESIAN coordinates..")')    IR0494
        call prjfc(hess,xparam,nvar)
      endif
      IJ=0
      DO 80 I=1,NVAR
         DO 80 J=1,I
            IJ=IJ+1
            HESSC(IJ)=HESS(J,I)
   80 CONTINUE
      CALL RSP(HESSC,NVAR,NVAR,EIGVAL,UC)
      IJ=0
      DO 90 I=1,NVAR
         IF (ABS(EIGVAL(I)).LT.TMSIX) EIGVAL(I)=ZERO
         DO 90 J=1,NVAR
            IJ=IJ+1
            U(J,I)=UC(IJ)
   90 CONTINUE
      IF (IPRNT.GE.3) CALL PRTHES(EIGVAL,NVAR)
      IF (MXSTEP.EQ.0) nstep=0
      IF (MXSTEP.EQ.0) GOTO 280                                                 

      NEG=0                                                                     
      DO 100 I=1,NVAR                                                           
         IF (EIGVAL(I) .LT. ZERO)NEG=NEG+1                                     
  100 CONTINUE                                                                  
      IF (IPRNT.GE.1)WRITE(6,110)NEG,(eigval(i),i=1,neg)
  110 FORMAT(/,10X,'HESSIAN HAS',I3,' NEGATIVE EIGENVALUE(S)',6f7.1,/)
c     if an eigenvalue has been zero out it is probably one of the T,R modes
c     in a cartesian optimization. zero corresponding fx to allow formation
c     of step without these contributions. a more safe criteria for deciding
c     whether this actually is a cartesian optimization should be put in 
c     some day...
      DO 120 I=1,NVAR                                                           
         FX(I)=DDOT(NVAR,U(1,I),1,GRAD,1)                               IR0494
         if (abs(eigval(i)).eq.zero) fx(i)=zero
  120 CONTINUE                                                                  
c     form geometry step d
130   CALL FORMD(EIGVAL,FX,NVAR,DMAX,OSMIN,lts,LRJK,lorjk,rrscal,donr)  
c     if lorjk is true, then ts mode overlap is less than omin, reject prev step
      if (lorjk) then
         if (iprnt.ge.1)write(6,*)'      Now undoing previous step'
         dmax=odmax
         dd=odd
         olde=oolde
         do 131 i=1,nvar
            fx(i)=oldfx(i)
                   oldf(i)=ooldf(i)
            eigval(i)=oldeig(i)
            do 132 j=1,nvar
               hess(i,j)=oldhss(i,j)
               u(i,j)=oldu(i,j)
132         CONTINUE
131      CONTINUE
         DO 140 I=1,NVAR
            XPARAM(I)=XPARAM(I)-D(I)
            K=LOC(1,I)
            L=LOC(2,I)
            GEO(L,K)=XPARAM(I)
140      CONTINUE
         IF(NDEP.NE.0) CALL SYMTRY
         dmax=min(dmax,dd)/two
         odmax=dmax
         odd=dd
         nstep=nstep-1
         if (dmax.lt.dmin) goto 230
      if (iprnt.ge.1)write(6,*)
     1'      Finish undoing, now going for new step'
         goto 130
      endif
C
C  FORM NEW TRIAL XPARAM AND GEO
C
      DO 150 I=1,NVAR
         XPARAM(I)=XPARAM(I)+D(I)
         K=LOC(1,I)
         L=LOC(2,I)
         GEO(L,K)=XPARAM(I)
  150 CONTINUE
      IF(NDEP.NE.0) CALL SYMTRY
C
C     COMPARE PREDICTED E-CHANGE WITH ACTUAL 
C
      depre=zero
      imode=1
      if (mode.ne.0)imode=mode
      do 160 i=1,nvar
         xtmp=xlamd
         if (lts .and. i.eq.imode) xtmp=xlamd0
         if (abs(xtmp-eigval(i)).lt.tmtwo) then
         ss=zero
         else
         ss=skal*fx(i)/(xtmp-eigval(i))
         endif
         frodo=ss*fx(i) + pt5*ss*ss*eigval(i)
c        write(6,88)i,fx(i),ss,xtmp,eigval(i),frodo
         depre=depre+frodo
160   continue
c88   format(i3,f10.3,f10.6,f10.3,4f10.6)
C
C     GET GRADIENT FOR NEW GEOMETRY 
C
      CALL COMPFG(XPARAM, .TRUE., FUNCT, .TRUE., GRAD, .TRUE.)
      IF(GNMIN) GNTEST=SQRT(DDOT(NVAR,GRAD,1,GRAD,1))                   IR0494
      DEACT = FUNCT-OLDE
      RATIO = DEACT/DEPRE
      if(iprnt.ge.1)WRITE(6,170)DEACT,DEPRE,RATIO       
  170 FORMAT(5X,'ACTUAL, PREDICTED ENERGY CHANGE, RATIO',2F10.3,F10.5)
      lrjk=.false.
C     if this is a minimum search, don't allow the energy to raise
      if (.not.lts .and. funct.gt.olde) then
         if (iprnt.ge.1)write(6,180)funct,min(dmax,dd)/two
180      format(1x,'energy raises ',f10.4,' rejecting step, ',
     $             'reducing dmax to',f7.4)
         lrjk=.true.
      endif
      if (gnmin .and. gntest.gt.oldgn) then
         if (iprnt.ge.1)write(6,181)gntest,min(dmax,dd)/two
181      format(1x,'gradient norm raises ',f10.4,' rejecting step, ',
     $             'reducing dmax to',f7.4)
         lrjk=.true.
      endif
      if (lts .and. (ratio.lt.rmin .or. ratio.gt.rmax) .and.
     $(abs(depre).gt.demin .or. abs(deact).gt.demin)) then
         if (iprnt.ge.1)write(6,190)min(dmax,dd)/two
190   format(1x,'unacceptable ratio,',
     $          ' rejecting step, reducing dmax to',f7.4)
         lrjk=.true.
      endif
      if (lrjk) then
         DO 200 I=1,NVAR
            XPARAM(I)=XPARAM(I)-D(I)
            K=LOC(1,I)
            L=LOC(2,I)
            GEO(L,K)=XPARAM(I)
200      CONTINUE
         IF(NDEP.NE.0) CALL SYMTRY
         dmax=min(dmax,dd)/two
         if (dmax.lt.dmin) goto 230
         goto 130
      endif
      IF(IPRNT.GE.1)WRITE(6,210)DD
  210 FORMAT(5X,'STEPSIZE USED IS',F9.5)
      IF(IPRNT.GE.2) THEN
         WRITE(6,'('' CALCULATED STEP'')')
         WRITE(6,'(3X,8F9.5)')(D(I),I=1,NVAR)
      ENDIF
C
C     POSSIBLE USE DYNAMICAL TRUST RADIUS
      odmax=dmax
      odd=dd
      oolde=olde
      IF (LUPD .and. ( (RMX.gt.gmin) .or.
     $    (abs(depre).gt.demin .or. abs(deact).gt.demin) ) ) THEN
c     Fletcher recommend dmax=dmax/4 and dmax=dmax*2
c     these are are a little more conservative since hessian is being updated
c     don't reduce trust radius due to ratio for min searches
      if (lts .and. ratio.le.tmone .or. ratio.ge.three)
     $    dmax=min(dmax,dd)/two
      if (lts .and. ratio.ge.pt75 .and. ratio.le.(four/three) 
     $                  .and. dd.gt.(dmax-tmsix)) 
     $   dmax=dmax*sqrt(two)
c     allow wider limits for increasing trust radius for min searches
      if (.not.lts .and. ratio.ge.pt5 
     $                  .and. dd.gt.(dmax-tmsix)) 
     $   dmax=dmax*sqrt(two)
c     be brave if  0.90 < ratio < 1.10 ...
      if (abs(ratio-one).lt.tmone) dmax=dmax*sqrt(two)
      dmax=max(dmax,dmin-tmsix)
      dmax=min(dmax,ddmax)
      ENDIF
c     allow stepsize up to 0.1 in the end-game where changes are less 
c     than demin and gradient is less than gmin
      IF (LUPD .and. RMX.lt.gmin .and.
     $   (abs(depre).lt.demin .and. abs(deact).lt.demin) )
     $    dmax=max(dmax,tmone)
      if(iprnt.ge.1)WRITE(6,220)DMAX
 220  FORMAT(5X,'CURRENT TRUST RADIUS = ',F7.5)                  
230   if (dmax.lt.dmin) then
         write(6,240)dmin
240      format(/,5x,'TRUST RADIUS NOW LESS THAN ',F7.5,' OPTIMIZATION',
     $   ' TERMINATING',/,5X,
     1' GEOMETRY MAY NOT BE COMPLETELY OPTIMIZED')
         goto 270
      endif

C     CHECK STEPS AND ENOUGH TIME FOR ANOTHER PASS
      if (nstep.ge.mxstep) goto 280
C     IN USER UNFRIENDLY ENVIROMENT, SAVE RESULTS EVERY 1 CPU HRS
      ITTEST=AINT((TIME2-TIME0)/TDUMP)
      IF (ITTEST.GT.NTIME) THEN
         LDUMP=1
         NTIME=MAX(ITTEST,(NTIME+1))
         IPOW(9)=2
cmgc         TT0=SECOND()-TIME0
         CALL EFSAV(TT0,HESS,FUNCT,GRAD,XPARAM,PMAT,-NSTEP,NSTEP,BMAT,
     1              IPOW)
      ELSE
         LDUMP=0
      ENDIF
C     RETURN FOR ANOTHER CYCLE
      GOTO 10                                                                  
C
C     ****** OPTIMIZATION TERMINATION ******
C
  250 CONTINUE
      WRITE(6,260)RMX,TOL2
  260 FORMAT(/,5X,'RMS GRADIENT =',F9.5,'  IS LESS THAN CUTOFF =',
     1F9.5,//)
  270 IFLEPO=15
      LAST=1
C     SAVE HESSIAN ON FILE 9
      IPOW(9)=2
cmgc      TT0=SECOND()-TIME0
      CALL EFSAV(TT0,HESS,FUNCT,GRAD,XPARAM,PMAT,-NSTEP,NSTEP,BMAT,
     1           IPOW)
C     CALL COMPFG TO CALCULATE ENERGY FOR FIXING MO-VECTOR BUG
      CALL COMPFG(XPARAM, .TRUE., FUNCT, .TRUE., GRAD, .FALSE.)
      RETURN
  280 CONTINUE
C     WE RAN OUT OF TIME or too many iterations. DUMP RESULTS
      IF (TLEFT .LT. TSTEP*TWO) THEN
         WRITE(6,290)
  290    FORMAT(/,5X,'NOT ENOUGH TIME FOR ANOTHER CYCLE')
      ENDIF
      IF (nstep.ge.mxstep) THEN
         WRITE(6,300)
  300    FORMAT(/,5X,'EXCESS NUMBER OF OPTIMIZATION CYCLES')
      ENDIF
      IPOW(9)=1
cmgc      TT0=SECOND()-TIME0
      CALL EFSAV(TT0,HESS,FUNCT,GRAD,XPARAM,PMAT,-NSTEP,NSTEP,BMAT,
     1           IPOW)
      RETURN
      END
      SUBROUTINE EFSAV(TT0,HESS,FUNCT,GRAD,XPARAM,PMAT,IL,JL,BMAT,IPOW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
      CHARACTER ELEMNT*2, KEYWRD*80
      DIMENSION HESS(MAXPAR,*),GRAD(*),BMAT(MAXPAR,*),IPOW(9),
     1 XPARAM(*), PMAT(*)
**********************************************************************
*
* EFSAV STORES AND RETRIEVE DATA USED IN THE EF GEOMETRY
*        OPTIMISATION. VERY SIMILAR TO POWSAV.
*
*  ON INPUT HESS   = HESSIAN MATRIX, PARTIAL OR WHOLE.
*           GRAD   = GRADIENTS.
*           XPARAM = CURRENT STATE OF PARAMETERS.
*           IL     = INDEX OF HESSIAN,
*           JL     = CYCLE NUMBER REACHED SO-FAR.
*           BMAT   = "B" MATRIX!
*           IPOW   = INDICES AND FLAGS.
*           IPOW(9)= 0 FOR RESTORE, 1 FOR DUMP, 2 FOR SILENT DUMP
*
**********************************************************************
      COMMON /GEOVAR/ XDUMMY(MAXPAR), NVAR, LOC(2,MAXPAR)               IR0494
      COMMON /ELEMTS/ ELEMNT(107)
      COMMON /GEOSYM/ NDEP,LOCPAR(MAXPAR),IDEPFN(MAXPAR),
     1                     LOCDEP(MAXPAR)
      COMMON/OPTEF/OLDF(MAXPAR),D(MAXPAR),VMODE(MAXPAR),
     $U(MAXPAR,MAXPAR),DD,rmin,rmax,omin,xlamd,xlamd0,skal,
     $MODE,NSTEP,NEGREQ,IPRNT
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1                NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /LOCVAR/ LOCVAR(2,MAXPAR)
      COMMON /NUMSCF/ NSCF
      COMMON /KEYWRD/ KEYWRD
      COMMON /VALVAR/ VALVAR(MAXPAR),NUMVAR
      COMMON /DENSTY/ P(MPACK), PA(MPACK), PB(MPACK)
      COMMON /ALPARM/ ALPARM(3,MAXPAR),X0, X1, X2, JLOOP
C Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
     3       /MOLKSR/ FRACT
      COMMON /PATHI / LATOM,LPARAM                                      IR0494
      OPEN(UNIT=9,FILE='FOR009',STATUS='UNKNOWN',FORM='UNFORMATTED')
      REWIND 9
      OPEN(UNIT=10,FILE='FOR010',STATUS='UNKNOWN',FORM='UNFORMATTED')
      REWIND 10
      IR=9
      IF(IPOW(9) .EQ. 1 .OR. IPOW(9) .EQ. 2) THEN
         FUNCT1=SQRT(DDOT(NVAR,GRAD,1,GRAD,1))                          IR0494
         IF(IPOW(9).EQ.1)THEN
            WRITE(6,'(//10X,''CURRENT VALUE OF GRADIENT NORM =''
     1  ,F12.6)')FUNCT1
            WRITE(6,'(/10X,''CURRENT VALUE OF GEOMETRY'',/)')
            CALL GEOUT                                                  IR0494
         ENDIF
C
C  IPOW(1) AND IPOW(9) ARE USED ALREADY, THE REST ARE FREE FOR USE
C  
         IPOW(8)=NSCF
         WRITE(IR)IPOW,IL,JL,FUNCT,TT0
         WRITE(IR)(XPARAM(I),I=1,NVAR)
         WRITE(IR)(  GRAD(I),I=1,NVAR)
         WRITE(IR)((HESS(J,I),J=1,NVAR),I=1,NVAR)
         WRITE(IR)((BMAT(J,I),J=1,NVAR),I=1,NVAR)
         WRITE(IR)(OLDF(I),I=1,NVAR),(D(I),I=1,NVAR),(VMODE(I),I=1,NVAR)
         WRITE(IR)DD,MODE,NSTEP,NEGREQ
         LINEAR=(NVAR*(NVAR+1))/2
         WRITE(IR)(PMAT(I),I=1,LINEAR)
         LINEAR=(NORBS*(NORBS+1))/2
         WRITE(10)(PA(I),I=1,LINEAR)
         IF(NALPHA.NE.0)WRITE(10)(PB(I),I=1,LINEAR)
         IF(LATOM .NE. 0) THEN
            WRITE(IR)((ALPARM(J,I),J=1,3),I=1,NVAR)
            WRITE(IR)JLOOP,X0, X1, X2
         ENDIF
         CLOSE(9)
         CLOSE(10)
         RETURN
      ELSE
C#         WRITE(6,'(//10X,'' READING DATA FROM DISK''/)')
         READ(IR,END=10,ERR=10)IPOW,IL,JL,FUNCT,TT0
         NSCF=IPOW(8)
         I=TT0/1000000
         TT0=TT0-I*1000000
         WRITE(6,'(//10X,''TOTAL TIME USED SO FAR:'',
     1    F13.2,'' SECONDS'')')TT0
         WRITE(6,'(  10X,''              FUNCTION:'',F17.6)')FUNCT
         READ(IR)(XPARAM(I),I=1,NVAR)
         READ(IR)(  GRAD(I),I=1,NVAR)
         READ(IR)((HESS(J,I),J=1,NVAR),I=1,NVAR)
         READ(IR)((BMAT(J,I),J=1,NVAR),I=1,NVAR)
         READ(IR)(OLDF(I),I=1,NVAR),(D(I),I=1,NVAR),(VMODE(I),I=1,NVAR)
         READ(IR)DD,MODE,NSTEP,NEGREQ
         LINEAR=(NVAR*(NVAR+1))/2
         READ(IR)(PMAT(I),I=1,LINEAR)
         LINEAR=(NORBS*(NORBS+1))/2
C        READ DENSITY MATRIX
         READ(10)(PA(I),I=1,LINEAR)
         IF(NALPHA.NE.0)READ(10)(PB(I),I=1,LINEAR)
         IF(LATOM.NE.0) THEN
            READ(IR)((ALPARM(J,I),J=1,3),I=1,NVAR)
            READ(IR)JLOOP,X0, X1, X2
            IL=IL+1
         ENDIF
         CLOSE(9)
         CLOSE(10)
         RETURN
   10    WRITE(6,'(//10X,''NO RESTART FILE EXISTS!'')')
         STOP
      ENDIF
      END
      SUBROUTINE EFSTR(XPARAM,FUNCT,IHESS,NTIME,ILOOP,IGTHES,MXSTEP,
     $                 IRECLC,IUPD,DMAX,DDMAX,dmin,TOL2,TOTIME,nvar,
     $                 SCF1,LUPD,ldump,rrscal,donr,gnmin)                IR0494
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C
#include "mopcvar.par"
C
      DIMENSION XPARAM(*)                                                       
C                                                                               
      COMMON /ISTOPE/ AMS(107)                                                  
      COMMON /LAST  / LAST                                                      
      COMMON /KEYWRD/ KEYWRD                                                    
      COMMON /TIMER/ TIME0                                                     
      COMMON /GRADNT/ GRAD(MAXPAR),GNFINA                                       
C Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
     3       /MOLKSR/ FRACT
      COMMON /NUMCAL/ NUMCAL                                                    
      COMMON /NLLCOM/ HESS(MAXPAR,MAXPAR),BMAT(MAXPAR,MAXPAR),                  
     1                PMAT(MAXPAR,MAXPAR)                                IR0494
      COMMON/OPTEF/OLDF(MAXPAR),D(MAXPAR),VMODE(MAXPAR),
     $U(MAXPAR,MAXPAR),DD,rmin,rmax,omin,xlamd,xlamd0,skal,
     $MODE,NSTEP,NEGREQ,IPRNT
      DIMENSION IPOW(9)               
      LOGICAL RESTRT,SCF1,LDUM,LUPD,rrscal,donr,gnmin 
      CHARACTER*80 KEYWRD                                               IR0494
      CHARACTER CHDOT*1,ZERO*1,NINE*1,CH*1
      DATA CHDOT,ZERO,NINE   /'.','0','9'/                     
      DATA  ICALCN,ZZERO  /0,0.D0/                          
C     GET ALL INITIALIZATION DATA                                               
         NVAR=ABS(NVAR)
         LDUMP=0
         ICALCN=NUMCAL
         LUPD=(INDEX(KEYWRD,' IUPD=0') .EQ. 0)                          IR0494
         RESTRT=(INDEX(KEYWRD,'RESTART') .NE. 0)
         SCF1=(INDEX(KEYWRD,'1SCF') .NE. 0)
         NSTEP=0
         IHESS=0
         LAST=0
         NTIME=0
         ILOOP=1
         IMIN=INDEX(KEYWRD,' EF')
         IF(IMIN.NE.0) THEN
            MODE=0
            IGTHES=0
            IUPD  =2
            NEGREQ=0
            ddmax=0.5d0
         ENDIF
         ITS=INDEX(KEYWRD,' TS')
         IF(ITS.NE.0) THEN
            MODE=1
            IGTHES=1
            IUPD  =1
            NEGREQ=1
            rmin=0.0d0
            rmax=4.0d0
            omin=0.8d0
            ddmax=0.3d0
         ENDIF
         rrscal=.false.
         I=INDEX(KEYWRD,' RSCAL') 
         IF(I.NE.0) rrscal=.true.
         donr=.true.
         I=INDEX(KEYWRD,' NONR') 
         IF(I.NE.0) donr=.false.
         gnmin=.false.
         I=INDEX(KEYWRD,' GNMIN') 
         IF(I.NE.0) gnmin=.true.
         IPRNT=0
         IP=INDEX(KEYWRD,' PRNT=') 
         IF(IP.NE.0) IPRNT=READA(KEYWRD,IP) 
         IF(IPRNT.GT.5)IPRNT=5                                                  
         IF(IPRNT.LT.0)IPRNT=0                                                  
         MXSTEP=100                                                             
         I=INDEX(KEYWRD,' CYCLES=')                                             
         IF(I.NE.0) MXSTEP=READA(KEYWRD,I)                                      
         IF (I.NE.0 .AND. MXSTEP.EQ.0 .AND. IP.EQ.0) IPRNT=3
         IRECLC=999999
         I=INDEX(KEYWRD,' RECALC=')
         IF(I.NE.0) IRECLC=READA(KEYWRD,I)
         I=INDEX(KEYWRD,' IUPD=')
         IF(I.NE.0) IUPD=READA(KEYWRD,I)
         I=INDEX(KEYWRD,' MODE=')
         IF(I.NE.0) MODE=READA(KEYWRD,I)
         DMIN=1.0D-3
         I=INDEX(KEYWRD,' DDMIN=')
         IF(I.NE.0) DMIN=READA(KEYWRD,I)
         DMAX=0.2D0
         I=INDEX(KEYWRD,' DMAX=')
         IF(I.NE.0) DMAX=READA(KEYWRD,I)
         I=INDEX(KEYWRD,' DDMAX=')
         IF(I.NE.0) DDMAX=READA(KEYWRD,I)
         TOL2=1.D+0
         IF(INDEX(KEYWRD,' PREC') .NE. 0) TOL2=5.D-2
         I=INDEX(KEYWRD,' GNORM=')
         IF(I.NE.0) TOL2=READA(KEYWRD,I)
         IF(INDEX(KEYWRD,' LET').EQ.0.AND.TOL2.LT.0.01D0)THEN
            WRITE(6,'(/,A)')'  GNORM HAS BEEN SET TOO LOW, RESET TO 0
     1.01. SPECIFY LET AS KEYWORD TO ALLOW GNORM LESS THAN 0.01'
            TOL2=0.01D0
         ENDIF
         I=INDEX(KEYWRD,' HESS=')
         IF(I.NE.0) IGTHES=READA(KEYWRD,I)
         I=INDEX(KEYWRD,' RMIN=')
         IF(I.NE.0) RMIN=READA(KEYWRD,I)
         I=INDEX(KEYWRD,' RMAX=')
         IF(I.NE.0) RMAX=READA(KEYWRD,I)
         I=INDEX(KEYWRD,' OMIN=')
         IF(I.NE.0) OMIN=READA(KEYWRD,I)
C fix timing problem                                                     IR0494
         I=INDEX(KEYWRD,' T=')
         IF(I.NE.0) THEN
            TOTIME=READA(KEYWRD,I)
         ELSE
            TOTIME=MAXTIM
         ENDIF
C   DONE WITH ALL INITIALIZING STUFF.
C   CHECK THAT OPTIONS REQUESTED ARE RESONABLE
         IF((ITS.NE.0).AND.(IUPD.EQ.2))THEN
            WRITE(6,*)' TS SEARCH AND BFGS UPDATE WILL NOT WORK'
            STOP
         ENDIF
         IF((ITS.NE.0).AND.(IGTHES.EQ.0))THEN
            WRITE(6,*)' TS SEARCH REQUIRE BETTER THAN DIAGONAL HESSIAN'
            STOP
         ENDIF
         IF((IGTHES.LT.0).OR.(IGTHES.GT.3))THEN
            WRITE(6,*)' UNRECOGNIZED HESS OPTION',IGTHES
            STOP
         ENDIF
         IF((OMIN.LT.0.d0).OR.(OMIN.GT.1.d0))THEN
            WRITE(6,*)' OMIN MUST BE BETWEEN 0 AND 1',OMIN
            STOP
         ENDIF
         IF (RESTRT) THEN
C
C   RESTORE DATA. I INDICATES (HESSIAN RESTART OR OPTIMIZATION
C   RESTART). IF I .GT. 0 THEN HESSIAN RESTART AND I IS LAST
C   STEP CALCULATED IN THE HESSIAN. IF I .LE. 0 THEN J (NSTEP)
C   IN AN OPTIMIZATION HAS BEEN DONE.
C
            IPOW(9)=0
            mtmp=mode
            CALL EFSAV(TT0,HESS,FUNCT,GRAD,XPARAM,PMAT,I,J,BMAT,IPOW)
            mode=mtmp
            ILOOP=I
            IF (I .GT. 0) THEN
               IGTHES=4
               NSTEP=J
               WRITE(6,'(10X,''RESTARTING HESSIAN AT POINT'',I4)')ILOOP
               IF(NSTEP.NE.0)WRITE(6,'(10X,''IN OPTIMIZATION STEP'',I4)'
     1)NSTEP
            ELSE
               NSTEP=J
               WRITE(6,'(//10X,''RESTARTING OPTIMIZATION AT STEP'',I4)')
     1NSTEP
               DO 26 I=1,NVAR
   26          GRAD(I)=ZZERO
               CALL COMPFG(XPARAM, .TRUE., FUNCT, .TRUE., GRAD, .TRUE.)
            ENDIF
         ELSE
C   NOT A RESTART, WE NEED TO GET THE GRADIENTS
            DO 30 I=1,NVAR
   30       GRAD(I)=ZZERO
            CALL COMPFG(XPARAM, .TRUE., FUNCT, .TRUE., GRAD, .TRUE.)
         ENDIF
      return
      end
      SUBROUTINE FORMD(EIGVAL,FX,NVAR,DMAX,
     1osmin,ts,lrjk,lorjk,rrscal,donr)
C     This version forms geometry step by either pure NR, P-RFO or QA
C     algorithm, under the condition that the steplength is less than dmax
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION LAMDA,lamda0
C
#include "mopcvar.par"
C
      logical ts,rscal,frodo1,frodo2,lrjk,lorjk,rrscal,donr              IR0794
      DIMENSION EIGVAL(MAXPAR),FX(MAXPAR)                                       
      COMMON/OPTEF/OLDF(MAXPAR),D(MAXPAR),VMODE(MAXPAR),
     $U(MAXPAR,MAXPAR),DD,rmin,rmax,omin,xlamd,xlamd0,skal,
     $MODE,NSTEP,NEGREQ,IPRNT
      DATA ZERO/0.0D0/, HALF/0.5D0/, TWO/2.0D+00/, TOLL/1.0D-8/         
      DATA STEP/5.0D-02/, TEN/1.0D+1/, ONE/1.0D+0/, BIG/1.0D+3/         
      DATA FOUR/4.0D+00/
      DATA TMTWO/1.0D-2/, TMSIX/1.0D-06/, SFIX/1.0D+01/, EPS/1.0D-12/ 
C                                                                       
      MAXIT=999                                                         
      NUMIT=0                                                           
      SKAL=ONE
      rscal=rrscal
      it=0
      jt=1
      if (ts) then
      IF(MODE.NE.0) THEN 
      CALL EFOVLP(dmax,osmin,NEWMOD,NVAR,lorjk)
      if (lorjk) return
C                                                                               
C  ON RETURN FROM EFOVLP, NEWMOD IS THE TS MODE
C                                                                               
      IF(NEWMOD.NE.MODE .and. iprnt.ge.1) WRITE(6,1000) MODE,NEWMOD
1000  FORMAT(5X,'WARNING! MODE SWITCHING. WAS FOLLOWING MODE ',I3,             
     $       ' NOW FOLLOWING MODE ',I3)                                         
      MODE=NEWMOD                                                               
      IT=MODE                                                                   
      ELSE
      IT=1
      ENDIF
      eigit=eigval(it)
      IF (IPRNT.GE.1) THEN                                                      
         WRITE(6,900)IT,EIGIT
         WRITE(6,910)(U(I,IT),I=1,NVAR)                                         
900      FORMAT(/,5X,'TS MODE IS NUMBER',I3,' WITH EIGENVALUE',F9.1,/,          
     *5X,'AND COMPONENTS',/)                                                    
910      FORMAT(5X,8F9.4)                                                       
      ENDIF                                                                     
      endif
      if (it.eq.1) jt=2
      eone=eigval(jt)                                                    
      ssmin=max(abs(eone)*eps,(ten*eps))
      ssmax=max(big,abs(eone))
      ssmax=ssmax*big
      sstoll=toll
      d2max=dmax*dmax                                                   
c     write(6,*)'from formd, eone, ssmin, ssmax, sstoll',
c    $eone,ssmin,ssmax,sstoll
      
C  SOLVE ITERATIVELY FOR LAMDA                                          
C  INITIAL GUESS FOR LAMDA IS ZERO EXCEPT NOTE THAT                     
C  LAMDA SHOULD BE LESS THAN EIGVAL(1)                                  
C  START BY BRACKETING ROOT, THEN HUNT IT DOWN WITH BRUTE FORCE BISECT. 
C                                                                       
         frodo1=.false.
         frodo2=.false.
         LAMDA=ZERO                                                     
         lamda0=zero
      if (ts .and. eigit.lt.zero .and. eone.ge.zero .and. donr) then
         if (iprnt.ge.1) then
         write(6,*)' ts search, correct hessian, trying pure NR step'
         endif
         goto 776
      endif
      if (.not.ts .and. eone.ge.zero .and. donr) then
         if (iprnt.ge.1) then
         write(6,*)' min search, correct hessian, trying pure NR step'
         endif
         goto 776
      endif
5     if (ts) then
         lamda0=eigval(it)+sqrt(eigval(it)**2+four*fx(it)**2)
         lamda0=lamda0*half
         if (iprnt.ge.1)WRITE(6,1030) LAMDA0 
      endif
         SSTEP = STEP                                                          
         IF(EONE.LE.ZERO) LAMDA=EONE-SSTEP                              
         IF(EONE.GT.ZERO) SSTEP=EONE                                           
         BL = LAMDA - SSTEP                                             
         BU = LAMDA + SSTEP*HALF                                        
20       FL = ZERO                                                      
         FU = ZERO                                                      
         DO 30 I = 1,NVAR                                               
            if (i.eq.it) goto 30
            FL   = FL + (FX(I)*FX(I))/(BL-EIGVAL(I))                    
            FU   = FU + (FX(I)*FX(I))/(BU-EIGVAL(I))                    
30       CONTINUE                                                       
         FL = FL - BL                                                   
         FU = FU - BU                                                   
c        write(6,*)'bl,bu,fl,fu from brack'                             
c        write(6,668)bl,bu,fl,fu                                        
c668     format(6f20.15)
         IF (FL*FU .LT. ZERO) GOTO 40                                   
         BL = BL - (EONE-BL)                                            
         BU = BU + HALF*(EONE-BU)                                       
         IF (BL.LE.-SSMAX) then
            BL = -SSMAX
            frodo1=.true.
         endif
         IF (abs(eone-bu).le.ssmin) then
            BU = EONE-SSMIN           
            frodo2=.true.
         endif
         IF (frodo1.and.frodo2) THEN              
            WRITE(6,*)'NUMERICAL PROBLEMS IN BRACKETING LAMDA',
     $                    EONE,BL,BU,FL,FU
            write(6,*)' going for fixed step size....'                       
            goto 450                                                           
         ENDIF                                                          
         GOTO 20                                                        
                                                                        
40       CONTINUE                                                       
         NCNT = 0                                                       
         XLAMDA = ZERO                                                  
50       CONTINUE                                                       
         FL = ZERO                                                      
         FU = ZERO                                                      
         FM = ZERO                                                      
         LAMDA = HALF*(BL+BU)                                                  
         DO 60 I = 1,NVAR                                               
            if (i.eq.it) goto 60
            FL   = FL + (FX(I)*FX(I))/(BL-EIGVAL(I))                    
            FU   = FU + (FX(I)*FX(I))/(BU-EIGVAL(I))                    
            FM   = FM + (FX(I)*FX(I))/(LAMDA-EIGVAL(I))                 
60       CONTINUE                                                       
         FL = FL - BL                                                   
         FU = FU - BU                                                   
         FM = FM - LAMDA                                                
c        write(6,*)'bl,bu,lamda,fl,fu,fm from search'                   
c        write(6,668)bl,bu,lamda,fl,fu,fm                               
         IF (ABS(XLAMDA-LAMDA).LT.sstoll) GOTO 776
         NCNT = NCNT + 1                                                
         IF (NCNT.GT.1000) THEN                                         
            WRITE(6,*)'TOO MANY ITERATIONS IN LAMDA BISECT',
     $                    BL,BU,LAMDA,FL,FU
            STOP                                                        
         ENDIF                                                          
         XLAMDA = LAMDA                                                 
         IF (FM*FU.LT.ZERO) BL = LAMDA                                  
         IF (FM*FL.LT.ZERO) BU = LAMDA                                  
         GOTO 50                                                        
C                                                                       
776   if (iprnt.ge.1) WRITE(6,1031) LAMDA 
C                                                                       
C  CALCULATE THE STEP                                                   
C                                                                       
      DO 310 I=1,NVAR                                                   
      D(I)=ZERO                                                         
310   CONTINUE                                                          
      DO 330 I=1,NVAR                                                   
      if (lamda.eq.zero .and. abs(eigval(i)).lt.tmtwo) then
      temp=zero
      else
      TEMP=FX(I)/(LAMDA-EIGVAL(I))                                      
      endif
      if (i.eq.it) then
      TEMP=FX(IT)/(LAMDA0-EIGVAL(IT)) 
      endif
      if (iprnt.ge.5) write(6,*)'formd, delta step',i,temp
      DO 320 J=1,NVAR                                                   
      D(J)=D(J)+TEMP*U(J,I)                                             
320   CONTINUE                                                          
330   CONTINUE                                                          
      DD=SQRT(DDOT(NVAR,D,1,D,1))                                       IR0494
      if(lamda.eq.zero .and. lamda0.eq.zero .and.iprnt.ge.1)
     1 write(6,777)dd
777   format(1x,'pure NR-step has length',f10.5)
      if(lamda.ne.zero .and. lamda0.ne.-lamda .and.iprnt.ge.1) 
     1write(6,778)dd
778   format(1x,'P-RFO-step   has length',f10.5)
      if (dd.lt.(dmax+tmsix)) then
         xlamd=lamda
         xlamd0=lamda0
         return
      endif
      if (lamda.eq.zero .and. lamda0.eq.zero) goto 5
      if (rscal) then
         SKAL=DMAX/DD
         DO 160 I=1,NVAR
            D(I)=D(I)*SKAL
160      CONTINUE
         DD=SQRT(DDOT(NVAR,D,1,D,1))                                     IR0494
         IF(IPRNT.GE.1)WRITE(6,170)SKAL
170      FORMAT(5X,'CALCULATED STEP SIZE TOO LARGE, SCALED WITH',F9.5)
         xlamd=lamda
         xlamd0=lamda0
         return
      endif

450      LAMDA=ZERO                                                     
         frodo1=.false.
         frodo2=.false.
         SSTEP = STEP                                                          
         IF(EONE.LE.ZERO) LAMDA=EONE-SSTEP                              
         if (ts .and. -eigit.lt.eone) lamda=-eigit-sstep
         IF(EONE.GT.ZERO) SSTEP=EONE                                           
         BL = LAMDA - SSTEP                                             
         BU = LAMDA + SSTEP*HALF                                        
520      FL = ZERO                                                      
         FU = ZERO                                                      
         DO 530 I = 1,NVAR                                              
            if (i.eq.it) goto 530
            FL   = FL + (FX(I)/(BL-EIGVAL(I)))**2                       
            FU   = FU + (FX(I)/(BU-EIGVAL(I)))**2                       
530      CONTINUE                                                       
         if (ts) then
            FL   = FL + (FX(IT)/(BL+EIGVAL(IT)))**2                       
            FU   = FU + (FX(IT)/(BU+EIGVAL(IT)))**2                       
         endif
         FL = FL - d2max                                                
         FU = FU - d2max                                                
c        write(6,*)'bl,bu,fl,fu from brack2'                            
c        write(6,668)bl,bu,fl,fu                                        
         IF (FL*FU .LT. ZERO) GOTO 540                                  
         BL = BL - (EONE-BL)                                            
         BU = BU + HALF*(EONE-BU)                                       
         IF (BL.LE.-SSMAX) then
            BL = -SSMAX
            frodo1=.true.
         endif
         IF (abs(eone-bu).le.ssmin) then
            BU = EONE-SSMIN           
            frodo2=.true.
         endif
         IF (frodo1.and.frodo2) THEN              
            WRITE(6,*)'NUMERICAL PROBLEMS IN BRACKETING LAMDA',
     $                    EONE,BL,BU,FL,FU
            write(6,*)' going for fixed level shifted NR step...'
c           both lamda searches failed, go for fixed level shifted nr    
c           this is unlikely to produce anything useful, but maybe we're lucky
            lamda=eone-sfix                                                    
            lamda0=eigit+sfix
            rscal=.true.                                                
            goto 776                                                    
         ENDIF                                                          
         GOTO 520                                                       
                                                                        
540      CONTINUE                                                       
         NCNT = 0                                                       
         XLAMDA = ZERO                                                  
550      CONTINUE                                                       
         FL = ZERO                                                      
         FU = ZERO                                                      
         FM = ZERO                                                      
         LAMDA = HALF*(BL+BU)                                                  
         DO 560 I = 1,NVAR                                              
            if (i.eq.it) goto 560
            FL   = FL + (FX(I)/(BL-EIGVAL(I)))**2                       
            FU   = FU + (FX(I)/(BU-EIGVAL(I)))**2                       
            FM   = FM + (FX(I)/(LAMDA-EIGVAL(I)))**2                    
560      CONTINUE                                                       
         if (ts) then
            FL   = FL + (FX(IT)/(BL+EIGVAL(IT)))**2                       
            FU   = FU + (FX(IT)/(BU+EIGVAL(IT)))**2                       
            FM   = FM + (FX(IT)/(LAMDA+EIGVAL(IT)))**2                    
         endif
         FL = FL - d2max                                                
         FU = FU - d2max                                                
         FM = FM - d2max                                                
c        write(6,*)'bl,bu,lamda,fl,fu,fm from search2'                  
c        write(6,668)bl,bu,lamda,fl,fu,fm                               
         IF (ABS(XLAMDA-LAMDA).LT.sstoll) GOTO 570                        
         NCNT = NCNT + 1                                                
         IF (NCNT.GT.1000) THEN                                         
            WRITE(6,*)'TOO MANY ITERATIONS IN LAMDA BISECT',
     $                    BL,BU,LAMDA,FL,FU
            STOP                                                        
         ENDIF                                                          
         XLAMDA = LAMDA                                                 
         IF (FM*FU.LT.ZERO) BL = LAMDA                                  
         IF (FM*FL.LT.ZERO) BU = LAMDA                                  
         GOTO 550                                                       
C                                                                       
570      CONTINUE                                                       
         lamda0=-lamda
         rscal=.true.                                                   
         goto 776                                                       
C                                                                       
1030  FORMAT(1X,'lamda that maximizes along ts modes =   ',F15.5)       
1031  FORMAT(1X,'lamda that minimizes along all modes =  ',F15.5)       
      END                                                               
      SUBROUTINE GETHES(XPARAM,IGTHES,NVAR,iloop,TOTIME) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C
#include "mopcvar.par"
C
C     GET THE HESSIAN. DEPENDING ON IGTHES WE GET IT FROM :                     
C                                                                               
C      0 : DIAGONAL MATRIX, DGHSX*I (DEFAULT FOR MIN-SEARCH)                    
C      1 : CALCULATE IT NUMERICALLY (DEFAULT FOR TS-SEARCH)                     
C      2 : READ IN FROM FTN009                                                  
C      3 : CALCULATE IT BY DOUBLE NUMERICAL DIFFERENTIATION
C      4 : READ IN FROM FTN009 (DURING RESTART, PARTLY OR WHOLE,                
C          ALREADY DONE AT THIS POINT)                                          
      COMMON /GEOVAR/ XDUMMY(MAXPAR), NDUM, LOC(2,MAXPAR)               IR0494
      COMMON /GEOM  / GEO(3,NUMATM)                                             
      COMMON /GEOSYM/ NDEP,LOCPAR(MAXPAR),IDEPFN(MAXPAR),LOCDEP(MAXPAR)         
      COMMON /LAST  / LAST                                                      
      COMMON /KEYWRD/ KEYWRD                                                    
      COMMON /TIMER/ TIME0                                                     
      COMMON /GRADNT/ GRAD(MAXPAR),GNFINA                                       
C Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
     3       /MOLKSR/ FRACT
      COMMON /NUMCAL/ NUMCAL                                                    
      COMMON /SIGMA2/ GNEXT1(MAXPAR), GMIN1(MAXPAR)                             
      COMMON /NLLCOM/ HESS(MAXPAR,MAXPAR),BMAT(MAXPAR,MAXPAR),                  
     1                PMAT(MAXPAR,MAXPAR)                                IR0494
      COMMON /TIMDMP/ TLEFT, TDUMP
      COMMON/OPTEF/OLDF(MAXPAR),D(MAXPAR),VMODE(MAXPAR),
     $U(MAXPAR,MAXPAR),DD,rmin,rmax,omin,xlamd,xlamd0,skal,
     $MODE,NSTEP,NEGREQ,IPRNT
      DIMENSION IPOW(9), EIGVAL(MAXPAR),TVEC(MAXPAR),SVEC(MAXPAR),              
     *FX(MAXPAR),HESSC(MAXHES),UC(MAXPAR**2)                                    
      DIMENSION XPARAM(*),tmp(150,150)
      LOGICAL RESTRT,SCF1,LDUM 
      CHARACTER*80 KEYWRD                                               IR0494
      CHARACTER CHDOT*1,ZERO*1,NINE*1,CH*1
      DATA CHDOT,ZERO,NINE  /'.','0','9'/                     
      DATA  ICALCN,ZZERO,ONE,TWO    /0,0.D0,1.D0,2.D0/                          
C
      DATA DGHSS,DGHSA,DGHSD /1000.d0,500.d0,200.d0/          
      DATA XINC /1.d-3/          
C     DGHSX IS HESSIAN DIAGONAL FOR IGTHES=0 (STRETCHING, ANGLE,
C     DIHEDRAL).  THE VALUES SHOULD BE 'OPTIMUM' FOR CYCLOHEXANONE
C     XINC IS STEPSIZE FOR HESSIAN CALCULATION. TESTS SHOWS THAT IT SHOULD
C     BE IN THE RANGE 10(-2) TO 10(-4). 10(-3) APPEARS TO BE 
C     A REASONABLE COMPROMISE BETWEEN ACCURACY AND NUMERICAL PROBLEMS
      IF (IGTHES.EQ.0) THEN
         WRITE(6,60)
   60    FORMAT(/,10X,'DIAGONAL MATRIX USED AS START HESSIAN',/)
         DO 70 I=1,NVAR
            DO 70 J=1,NVAR
               HESS(I,J)=ZZERO
   70    CONTINUE
         IJ=1
         DO 80 J=1,NUMATM
            DO 80 I=1,3
               IF (LOC(2,IJ).EQ.I.AND.LOC(1,IJ).EQ.J)THEN
                  IF (I.EQ.1)HESS(IJ,IJ)=DGHSS
                  IF (I.EQ.2)HESS(IJ,IJ)=DGHSA
                  IF (I.EQ.3)HESS(IJ,IJ)=DGHSD
                  IJ=IJ+1
               ENDIF
   80    CONTINUE
         IJ=IJ-1
         IF(IJ.NE.NVAR)WRITE(*,*)'ERROR IN IGTHES=0,IJ,NVAR',IJ,NVAR
      ENDIF
C
      IF (IGTHES.EQ.2) THEN
         WRITE(6,100)
  100    FORMAT(/,10X,'HESSIAN READ FROM DISK',/)
         IPOW(9)=0
C        USE DUMMY ARRAY FOR CALL EXCEPT FOR HESSIAN
C        TEMPORARY SET NALPHA = 0, THEN WE CAN READ HESSIAN FROM RHF
C        RUN FOR USE IN SAY UHF RUNS
C        ALSO SAVE MODE, TO ALLOW FOLLOWING A DIFFERENT MODE THAN THE ONE
C        CURRENTLY ON RESTART FILE
         nxxx=nalpha
         nalpha=0
         mtmp=mode
         CALL EFSAV(TDM,HESS,FDMY,GNEXT1,GMIN1,PMAT,IIDUM,J,BMAT,IPOW)
         nalpha=nxxx
         mode=mtmp
         nstep=0
      ENDIF
      IF((IGTHES.EQ.1).OR.(IGTHES.EQ.3).OR.(IGTHES.EQ.4))THEN
C       IF IGTHES IS .EQ. 4, THEN THIS IS A HESSIAN RESTART.
C       USE GNEXT1 AND DUMMY FOR CALLS TO COMPFG DURING HESSIAN
C       CALCULATION
         IF (IGTHES.EQ.1)WRITE(6,190)
  190    FORMAT(/,10X,'HESSIAN CALCULATED NUMERICALLY',/)
         IF (IGTHES.EQ.3)WRITE(6,191)
  191    FORMAT(/,10X,'HESSIAN CALCULATED DOUBLE NUMERICALLY',/)
            IF(IPRNT.GE.5)WRITE(6,'(I3,12(8F9.4,/3X))')
     1    0,(Grad(IF),IF=1,NVAR)
cmgc         TIME1=SECOND()
         TSTORE=TIME1
         DO 210 I=ILOOP,NVAR
            XPARAM(I)=XPARAM(I) + XINC
            CALL COMPFG(XPARAM, .TRUE., DUMMY, .TRUE., GNEXT1, .TRUE.)
            IF(IPRNT.GE.5)WRITE(6,'(I3,12(8F9.4,/3X))')
     1    I,(GNEXT1(IF),IF=1,NVAR)
            XPARAM(I)=XPARAM(I) - XINC
            if (igthes.eq.3) then
            XPARAM(I)=XPARAM(I) - XINC
            CALL COMPFG(XPARAM, .TRUE., DUMMY, .TRUE., GMIN1, .TRUE.)
            IF(IPRNT.GE.5)WRITE(6,'(I3,12(8F9.4,/3X))')
     1    -I,(GMIN1(IF),IF=1,NVAR)
            XPARAM(I)=XPARAM(I) + XINC
            DO 199 J=1,NVAR
  199       HESS(I,J)= (GNEXT1(J)-GMIN1(J))/(XINC+XINC)
            else
            DO 200 J=1,NVAR
  200       HESS(I,J)= (GNEXT1(J)-GRAD(J))/XINC
            endif
cmgc            TIME2=SECOND()
            TSTEP=TIME2-TIME1
            TLEFT=TOTIME-TIME2                                           IR0494
            TIME1=TIME2
            IF( TLEFT .LT. TSTEP*TWO) THEN
C
C  STORE PARTIAL HESSIAN PATRIX
C  STORE GRADIENTS FOR GEOMETRY AND ILOOP AS POSITIVE
               WRITE(6,'(A)')' NOT ENOUGH TIME TO COMPLETE HESSIAN'
               WRITE(6,'(A,I4)')' STOPPING IN HESSIAN AT COORDINATE:',I
               IPOW(9)=1
cmgc               TT0=SECOND()-TIME0
               CALL EFSAV(TT0,HESS,FUNCT,GRAD,XPARAM,PMAT,I,NSTEP,BMAT,
     1IPOW)
               STOP
            ENDIF
  210    CONTINUE
C     fix last entry in geo array, this is currently at value-xinc
         K=LOC(1,nvar)
         L=LOC(2,nvar)
         GEO(L,K)=XPARAM(nvar)
         IF(NDEP.NE.0) CALL SYMTRY
c        add all time used back to tleft, this will then be subtracted
c        again in main ef routine
cmgc         TIME2=SECOND()
         TSTEP=TIME2-TSTORE
         TLEFT=TOTIME-TIME2                                              IR0494
      ENDIF
C
C     SYMMETRIZE HESSIAN
      DO 220 I=1,NVAR
         DO 220 J=1,I-1
            HESS(I,J)=(HESS(I,J)+HESS(J,I))/TWO
            HESS(J,I)=HESS(I,J)
  220 CONTINUE
      RETURN
      END
      SUBROUTINE EFOVLP(dmax,osmin,NEWMOD,NVAR,lorjk)             
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mopcvar.par"
      COMMON /NUMCAL/ NUMCAL                                            IR1294
      COMMON/OPTEF/OLDF(MAXPAR),D(MAXPAR),VMODE(MAXPAR),
     $U(MAXPAR,MAXPAR),DD,rmin,rmax,omin,xlamd,xlamd0,skal,
     $MODE,NSTEP,NEGREQ,IPRNT
      dimension xo(maxpar)
      logical lorjk
      DATA ICALCN /0/                                                   IR1294
C
C  ON THE FIRST STEP SIMPLY DETERMINE WHICH MODE TO FOLLOW
C
c     IF(NSTEP.EQ.1) THEN
      IF(ICALCN.NE.NUMCAL) THEN                                         IR1294
         ICALCN=NUMCAL                                                  IR1294
         IF(MODE.GT.NVAR)THEN
            WRITE(6,*)'ERROR!! MODE IS LARGER THAN NVAR',MODE
            STOP
         ENDIF
         IT=MODE
         if (iprnt.ge.1) WRITE(6,40) MODE
   40 FORMAT(5X,'HESSIAN MODE FOLLOWING SWITCHED ON'/
     1     '     FOLLOWING MODE ',I3)
C
      ELSE
C
C  ON SUBSEQUENT STEPS DETERMINE WHICH HESSIAN EIGENVECTOR HAS
C  THE GREATEST OVERLAP WITH THE MODE WE ARE FOLLOWING
C
         IT=1
         lorjk=.false.
         TOVLP=DDOT(NVAR,U(1,1),1,VMODE,1)                              IR0494
         TOVLP=ABS(TOVLP)
c        xo(1)=tovlp
         DO 10 I=2,NVAR
            OVLP=DDOT(NVAR,U(1,I),1,VMODE,1)                            IR0494
            OVLP=ABS(OVLP)
c           xo(i)=ovlp
            IF(OVLP.GT.TOVLP) THEN
               TOVLP=OVLP
               IT=I
            ENDIF
   10    CONTINUE
C
         if (iprnt.ge.5) then
         do 11 j=1,5
         xxx=0.d0
         do 12 i=1,nvar
         if (xo(i).gt.xxx)ix=i
         if (xo(i).gt.xxx)xxx=xo(i)
12       continue
         xo(ix)=0.d0
         write(6,*)'overlaps',ix,xxx
11       continue
         endif

         if(iprnt.ge.1)WRITE(6,30) IT,TOVLP
         if (tovlp.lt.omin) then
            if (dmax.gt.osmin) then
            lorjk=.true.
            if (iprnt.ge.1)write(6,31)omin
            return
            else
            if (iprnt.ge.1)write(6,32)omin,dmax,osmin
            endif
         endif
      ENDIF
   30 FORMAT(5X,'OVERLAP OF CURRENT MODE',I3,' WITH PREVIOUS MODE IS ',
     $       F6.3)
   31 FORMAT(5X,'OVERLAP LESS THAN OMIN',
     1F6.3,' REJECTING PREVIOUS STEP')
   32 FORMAT(5X,'OVERLAP LESS THAN OMIN',F6.3,' BUT TRUST RADIUS',F6.3,
     $          ' IS LESS THAN',F6.3,/,5X,' ACCEPTING STEP')
C
C  SAVE THE EIGENVECTOR IN VMODE
C
      DO 20 I=1,NVAR
         VMODE(I)=U(I,IT)
   20 CONTINUE
C
      NEWMOD=IT
      RETURN
C
      END
      SUBROUTINE PRJFC(F,xparam,nvar) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#include "mopcvar.par"
C                                                                       
C  CALCULATES PROJECTED FORCE CONSTANT MATRIX (F).                      
C  THIS ROUTINE CAME ORIGINALLY FROM POLYRATE. IT IS USED BY PERMISSION
C  OF D. TRUHLAR. THE CURRENT VERSION IS LIFTED FROM GAMESS AND 
C  ADAPTED BY F.JENSEN, ODENSE, DK             
C  IF WE ARE AT A STATIONARY POINT (STPT=.T.), I.E. GNORM .LT. 10,      
C  THEN THE ROTATIONAL AND TRANSLATIONAL MODES ARE PROJECTED OUT        
C  AND THEIR FREQUENCIES BECOME IDENTICAL ZERO. IF NOT AT A STATIONARY  
C  POINT THEN THE MASS-WEIGHTED GRADIENT IS ALSO PROJECTED OUT AND      
C  THE CORRESPONDING FREQUENCY BECOME ZERO.                             
C ************************************************                      
C   X : MASS-WEIGHTED COORDINATE                                        
C   DX: NORMALIZED MASS-WEIGHTED GRADIENT VECTOR                        
C   F : MASS-WEIGHTED FORCE CONSTANT MATRIX                             
C   RM: INVERSION OF SQUARE ROOT OF MASS                                
C   P, COF: BUFFER                                                      
C                                                                       
      COMMON /ATMASS/ ATMASS(NUMATM)
      DIMENSION X(MAXPAR),RM(MAXPAR),F(MAXPAR,MAXPAR),                  
     *          P(MAXPAR,MAXPAR),COF(MAXPAR,MAXPAR)                     
      DIMENSION TENS(3,3,3),ROT(3,3),SCR(3,3),ISCR(6),CMASS(3)          
      dimension coord(3,numatm),dx(maxpar),xparam(maxpar)
      PARAMETER (ZERO=0.0d+00, ONE=1.0d+00, EPS=1.0d-14,                
     *           CUT5=1.0d-05, CUT8=1.0d-08)                            
C                                                                       
C TOTALLY ASYMMETRIC CARTESIAN TENSOR.                                  
      DATA TENS/ 0.0d+00,  0.0d+00,  0.0d+00,                           
     X           0.0d+00,  0.0d+00, -1.0d+00,                           
     X           0.0d+00,  1.0d+00,  0.0d+00,                           
     Y           0.0d+00,  0.0d+00,  1.0d+00,                           
     Y           0.0d+00,  0.0d+00,  0.0d+00,                           
     Y          -1.0d+00,  0.0d+00,  0.0d+00,                           
     Z           0.0d+00, -1.0d+00,  0.0d+00,                           
     Z           1.0d+00,  0.0d+00,  0.0d+00,                           
     Z           0.0d+00,  0.0d+00,  0.0d+00  /                         
C                                                                       
      natm=nvar/3                                                        
      nc1=nvar
      ij=1
      do 2 i=1,natm                                                    
         coord(1,i)=xparam(ij)
         coord(2,i)=xparam(ij+1)
         coord(3,i)=xparam(ij+2)
         ij=ij+3
2     continue                                                          
C     CALCULATE 1/SQRT(MASS)                                            
      L=0                                                               
      DO 3 I=1,NATM                                                     
         TMP=ONE/SQRT(ATMASS(I))                                        
         DO 3 J=1,3                                                     
            L=L+1                                                       
3           RM(L)=TMP                                                   
C     PREPARE GRADIENT                                                  
         DO 4 I=1,NC1                                                   
4        DX(I)=ZERO                                                     
C     FIND CMS AND CALCULATED MASS WEIGHTED COORDINATES                 
      totm=zero                                                         
      cmass(1)=zero                                                     
      cmass(2)=zero                                                     
      cmass(3)=zero                                                     
      DO 6 I=1,NATM                                                     
         TOTM=TOTM+ATMASS(I)                                            
         DO 6 J=1,3                                                     
            CMASS(J)=CMASS(J)+ATMASS(I)*COORD(J,I)                      
6     CONTINUE                                                          
      DO 7 J=1,3                                                        
7     CMASS(J)=CMASS(J)/TOTM                                            
      L=0                                                               
      DO 8 I=1,NATM                                                     
         DO 8 J=1,3                                                     
         TMP=SQRT(ATMASS(I))                                            
             L=L+1                                                      
             X(L)=TMP*(COORD(J,I)-CMASS(J))                             
8     CONTINUE                                                          
c     WRITE(6,9020)                                                     
c     CALL prsq(f,nc1,nc1,maxpar,1)                                          
c9020 FORMAT(/1X,'ENTER THE SUBROUTINE <PRJFC>'//                       
c    *        1X,'UNPROJECTED FORCE CONSTANT MATRIX (HARTREE/BOHR**2)') 
c     WRITE(6,*)' MASS-WEIGHTED COORDINATES AND CORRESPONDING GRADIENT' 
c     DO 9 I=1,NC1                                                      
c9       WRITE(6,*)X(I),DX(I)                                           
C                                                                       
C 2. COMPUTE INERTIA TENSOR.                                            
      DO 10 I=1,3                                                       
       DO 10 J=1,3                                                      
   10   ROT(I,J)=ZERO                                                   
      DO 20 I=1,NATM                                                    
       L=3*(I-1)+1                                                      
       ROT(1,1)=ROT(1,1)+X(L+1)**2+X(L+2)**2                            
       ROT(1,2)=ROT(1,2)-X(L)*X(L+1)                                    
       ROT(1,3)=ROT(1,3)-X(L)*X(L+2)                                    
       ROT(2,2)=ROT(2,2)+X(L)**2+X(L+2)**2                              
       ROT(2,3)=ROT(2,3)-X(L+1)*X(L+2)                                  
   20  ROT(3,3)=ROT(3,3)+X(L)**2+X(L+1)**2                              
      ROT(2,1)=ROT(1,2)                                                 
      ROT(3,1)=ROT(1,3)                                                 
      ROT(3,2)=ROT(2,3)                                                 
C                                                                       
CHECK THE INERTIA TENSOR.                                               
      CHK=ROT(1,1)*ROT(2,2)*ROT(3,3)                                    
      IF(ABS(CHK).GT.CUT8) GO TO 21                                     
c     WRITE(6,23)                                                       
c  23 FORMAT(/1X,'MATRIX OF INERTIA MOMENT')                            
c     CALL PRSQ(ROT,3,3,3,3)                                              
      IF(ABS(ROT(1,1)).GT.CUT8) GO TO 11                                
C X=0                                                                   
      IF(ABS(ROT(2,2)).GT.CUT8) GO TO 12                                
C X,Y=0                                                                 
      IF(ABS(ROT(3,3)).GT.CUT8) GO TO 13                                
      WRITE(6,14) ROT(1,1),ROT(2,2),ROT(3,3)                            
   14 FORMAT(1X,'EVERY DIAGONAL ELEMENTS ARE ZERO ?',3F20.10)           
      RETURN                                                            
C                                                                       
C* 1. X,Y=0 BUT Z.NE.0                                                  
   13 ROT(3,3)=ONE/ROT(3,3)                                             
      GO TO 22                                                          
C Y.NE.0                                                                
   12 IF(ABS(ROT(3,3)).GT.CUT8) GO TO 15                                
C* 2. X,Z=0 BUT Y.NE.0                                                  
      ROT(2,2)=ONE/ROT(2,2)                                             
      GO TO 22                                                          
C X.NE.0                                                                
   11 IF(ABS(ROT(2,2)).GT.CUT8) GO TO 16                                
      IF(ABS(ROT(3,3)).GT.CUT8) GO TO 17                                
C* 3. Y,Z=0 BUT X.NE.0                                                  
      ROT(1,1)=ONE/ROT(1,1)                                             
      GO TO 22                                                          
C* 4. X,Y.NE.0 BUT Z=0                                                  
   16 DET=ROT(1,1)*ROT(2,2)-ROT(1,2)*ROT(2,1)                           
      TRP=ROT(1,1)                                                      
      ROT(1,1)=ROT(2,2)/DET                                             
      ROT(2,2)=TRP/DET                                                  
      ROT(1,2)=-ROT(1,2)/DET                                            
      ROT(2,1)=-ROT(2,1)/DET                                            
      GO TO 22                                                          
C* 5. X,Z.NE.0 BUT Y=0                                                  
   17 DET=ROT(1,1)*ROT(3,3)-ROT(1,3)*ROT(3,1)                           
      TRP=ROT(1,1)                                                      
      ROT(1,1)=ROT(3,3)/DET                                             
      ROT(3,3)=TRP/DET                                                  
      ROT(1,3)=-ROT(1,3)/DET                                            
      ROT(3,1)=-ROT(3,1)/DET                                            
      GO TO 22                                                          
C* 6. Y,Z.NE.0 BUT X=0                                                  
   15 DET=ROT(3,3)*ROT(2,2)-ROT(3,2)*ROT(2,3)                           
      TRP=ROT(3,3)                                                      
      ROT(3,3)=ROT(2,2)/DET                                             
      ROT(2,2)=TRP/DET                                                  
      ROT(3,2)=-ROT(3,2)/DET                                            
      ROT(2,3)=-ROT(2,3)/DET                                            
      GO TO 22                                                          
   21 CONTINUE                                                          
C                                                                       
C.DEBUG.                                                                
c      CALL PRSQ(TENS(1,1,1),3,3,3,3)                                     
c      CALL PRSQ(TENS(1,1,2),3,3,3,3)                                     
c      CALL PRSQ(TENS(1,1,3),3,3,3,3)                                     
c      CALL PRSQ(ROT,3,3,3,3)                                             
C                                                                       
C 4. COMPUTE INVERSION MATRIX OF ROT.                                   
C     CALL MXLNEQ(ROT,3,3,DET,JRNK,EPS,SCR,+0)                          
C     IF(JRNK.LT.3) STOP 1                                              
      INFO=0                                                            
      CALL DGEFA(ROT,3,3,ISCR,INFO)                                     
      IF(INFO.NE.0) STOP                                                
      DET=ZERO                                                          
      CALL DGEDI(ROT,3,3,ISCR,DET,SCR,1)                                
C                                                                       
   22 CONTINUE                                                          
c     WRITE (6,702)                                                     
c 702 FORMAT(/1X,'INVERSE MATRIX OF MOMENT OF INERTIA.')                
c     CALL PRSQ(ROT,3,3,3,3)                                              
C                                                                       
C 5. TOTAL MASS ---> TOTM.                                              
C                                                                       
C 6. COMPUTE P MATRIX                                                   
C    ----------------                                                   
      DO 100 IP=1,NATM                                                  
       INDX=3*(IP-1)                                                    
       DO 100 JP=1,IP                                                   
        JNDX=3*(JP-1)                                                   
        DO 70 IC=1,3                                                    
         JEND=3                                                         
         IF(JP.EQ.IP) JEND=IC                                           
         DO 70 JC=1,JEND                                                
          SUM=ZERO                                                      
          DO 50 IA=1,3                                                  
           DO 50 IB=1,3                                                 
            IF(TENS(IA,IB,IC).EQ.0) GO TO 50                            
            DO 30 JA=1,3                                                
             DO 30 JB=1,3                                               
              IF(TENS(JA,JB,JC).EQ.0) GO TO 30                          
              SUM=SUM+TENS(IA,IB,IC)*TENS(JA,JB,JC)*ROT(IA,JA)*         
     &                X(INDX+IB)*X(JNDX+JB)                             
   30         CONTINUE                                                  
   50       CONTINUE                                                    
          II=INDX+IC                                                    
          JJ=JNDX+JC                                                    
          P(II,JJ)=SUM+DX(II)*DX(JJ)                                    
          IF(IC.EQ.JC) P(II,JJ)=P(II,JJ)+ONE/(RM(II)*RM(JJ)*TOTM)       
   70     CONTINUE                                                      
  100   CONTINUE                                                        
C                                                                       
C 7. COMPUTE DELTA(I,J)-P(I,J)                                          
      DO 110 I=1,NC1                                                    
       DO 110 J=1,I                                                     
        P(I,J)=-P(I,J)                                                  
        IF(I.EQ.J) P(I,J) = ONE +P(I,J)                                 
  110   CONTINUE                                                        
C                                                                       
C 8. NEGLECT SMALLER VALUES THAN 10**-8.                                
      DO 120 I=1,NC1                                                    
       DO 120 J=1,I                                                     
        IF(ABS(P(I,J)).LT.CUT8) P(I,J)=ZERO                             
        P(J,I)=P(I,J)                                                   
  120   CONTINUE                                                        
C                                                                       
C.DEBUG.                                                                
c     WRITE(6,703)                                                      
c 703 FORMAT(/1X,'PROJECTION MATRIX')                                   
c     CALL PRSQ(P,NC1,NC1,NC1)                                          
c     CALL PRSQ(P,NC1,NC1,maxpar,3)                                       
C                                                                       
C 10. POST AND PREMULTIPLY F BY P.                                      
C     USE COF FOR SCRATCH.                                              
      DO 150 I=1,NC1                                                    
       DO 150 J=1,NC1                                                   
        SUM=ZERO                                                        
        DO 140 K=1,NC1                                                  
  140    SUM=SUM+F(I,K)*P(K,J)                                          
  150   COF(I,J)=SUM                                                    
C                                                                       
C 11. COMPUTE P*F*P.                                                    
      DO 200 I=1,NC1                                                    
       DO 200 J=1,NC1                                                   
        SUM=ZERO                                                        
        DO 190 K=1,NC1                                                  
  190    SUM=SUM+P(I,K)*COF(K,J)                                        
  200   F(I,J)=SUM                                                      
C                                                                       
c     WRITE(6,9030)                                                     
c     CALL prsq(f,nc1,nc1,maxpar,1)                                          
c9030 FORMAT(/1X,'LEAVE THE SUBROUTINE <PRJFC>'//                       
c    *        1X,'PROJECTED FORCE CONSTANT MATRIX (HARTREE/BOHR**2)')   
      RETURN                                                            
      END                                                               
      SUBROUTINE PRTHES(EIGVAL,NVAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
#include "mopcvar.par"
      COMMON /NLLCOM/ HESS(MAXPAR,MAXPAR),BMAT(MAXPAR,MAXPAR),                  
     1                PMAT(MAXPAR,MAXPAR)                                IR0494
      COMMON/OPTEF/OLDF(MAXPAR),D(MAXPAR),VMODE(MAXPAR),                        
     $U(MAXPAR,MAXPAR),DD,rmin,rmax,omin,xlamd,xlamd0,skal,
     $MODE,NSTEP,NEGREQ,IPRNT
      DIMENSION EIGVAL(MAXPAR)
         IF (IPRNT.GE.4) THEN                                                   
         WRITE(6,*)' '                                                          
         WRITE(6,*)'              HESSIAN MATRIX'                               
         LOW=1                                                                  
         NUP=8                                                                  
540      NUP=MIN(NUP,NVAR)                                                      
         WRITE(6,1000) (I,I=LOW,NUP)                                            
         DO 550 I=1,NVAR                                                        
         WRITE(6,1010) I,(HESS(I,J),J=LOW,NUP)                                  
550      CONTINUE                                                               
         NUP=NUP+8                                                              
         LOW=LOW+8                                                              
         IF(LOW.LE.NVAR) GOTO 540                                               
         ENDIF                                                                  
         WRITE(6,*)' '                                                          
         WRITE(6,*)'              HESSIAN EIGENVALUES AND -VECTORS'             
         LOW=1                                                                  
         NUP=8                                                                  
560      NUP=MIN(NUP,NVAR)                                                      
         WRITE(6,1000) (I,I=LOW,NUP)                                            
         WRITE(6,1020) (EIGVAL(I),I=LOW,NUP)                                    
         DO 570 I=1,NVAR                                                        
         WRITE(6,1030) I,(U(I,J),J=LOW,NUP)                                     
570      CONTINUE                                                               
         NUP=NUP+8                                                              
         LOW=LOW+8                                                              
         IF(LOW.LE.NVAR) GOTO 560                                               
1000     FORMAT(/,3X,8I9)                                                       
1010     FORMAT(1X,I3,8F9.1)                                                    
1020     FORMAT(/,4X,8F9.1,/)                                                   
1030     FORMAT(1X,I3,8F9.4)                                                    
      RETURN
      END
      SUBROUTINE UPDHES(SVEC,TVEC,NVAR,IUPD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mopcvar.par"
      DIMENSION TVEC(*),SVEC(*)
      COMMON /NUMCAL/NUMCAL                                              IR1294
      COMMON/OPTEF/OLDF(MAXPAR),D(MAXPAR),VMODE(MAXPAR),
     $U(MAXPAR,MAXPAR),DD,rmin,rmax,omin,xlamd,xlamd0,skal,
     $MODE,NSTEP,NEGREQ,IPRNT
      COMMON /NLLCOM/ HESS(MAXPAR,MAXPAR),BMAT(MAXPAR,MAXPAR),                  
     1                PMAT(MAXPAR,MAXPAR)                                IR0494
      COMMON /GRADNT/ GRAD(MAXPAR),GNFINA
C
      DATA ZERO/0.0D0/
C
C  UPDATING OF THE HESSIAN
C  DEPENDS ON CURRENT GRADIENTS, OLD GRADIENTS AND THE
C  CORRECTION VECTOR USED ON THE LAST CYCLE
C  SVEC & TVEC ARE FOR TEMPORARY STORAGE
C
C  2 UPDATING PROCEDURES ARE POSSIBLE
C  (I)   THE POWELL UPDATE
C        THIS PRESERVES THE SYMMETRIC CHARACTER OF THE HESSIAN
C        WHILST ALLOWING ITS EIGENVALUE STRUCTURE TO CHANGE.
C        IT IS THE DEFAULT UPDATE FOR A TRANSITION STATE SEARCH
C  (II)  THE BFGS UPDATE
C        THIS UPDATE HAS THE IMPORTANT CHARACTERISTIC OF RETAINING
C        POSITIVE DEFINITENESS (NOTE: THIS IS NOT RIGOROUSLY
C        GUARANTEED, BUT CAN BE CHECKED FOR BY THE PROGRAM).
C        IT IS THE DEFAULT UPDATE FOR A MINIMUM SEARCH
C
C     SWITCH : IUPD
C       IUPD = 0  :  SKIP UPDATE
C       IUPD = 1  :  POWELL
C       IUPD = 2  :  BFGS
C
      IF(ICALCN.NE.NUMCAL) THEN                                         IR1294
         ICALCN=NUMCAL                                                  IR1294
         IF(IPRNT.GE.2) THEN
            IF (IUPD.EQ.0)WRITE(6,90)
            IF (IUPD.EQ.1)WRITE(6,80)
            IF (IUPD.EQ.2)WRITE(6,120)
         ENDIF
      ENDIF
      IF(IUPD.EQ.0) RETURN
      DO 5 I=1,NVAR
         TVEC(I)=ZERO
 5    CONTINUE
      DO 10 J=1,NVAR
         DO 10 I=1,NVAR
            TVEC(I)=TVEC(I) + HESS(I,J)*D(J)
   10 CONTINUE
C
      IF(IUPD.EQ.1) THEN
C
C   (I) POWELL UPDATE
C
         DO 20 I=1,NVAR
            TVEC(I)=GRAD(I)-OLDF(I)-TVEC(I)
            sVEC(I)=GRAD(I)-OLDF(I)
   20    CONTINUE
         DDS=DD*DD
         DDTD=DDOT(NVAR,TVEC,1,D,1)                                     IR0494
         DDTD=DDTD/DDS
C
         DO 40 I=2,NVAR
            DO 30 J=1,I-1
               TEMP=TVEC(I)*D(J) + D(I)*TVEC(J) - D(I)*DDTD*D(J)
               HESS(I,J)=HESS(I,J)+TEMP/DDS
               HESS(J,I)=HESS(I,J)
   30       CONTINUE
   40    CONTINUE
         DO 45 I=1,NVAR
            TEMP=D(I)*(2.0D0*TVEC(I) - D(I)*DDTD)
            HESS(I,I)=HESS(I,I)+TEMP/DDS
 45      CONTINUE
C
      ENDIF
      IF (IUPD.EQ.2) THEN
C
C  (II) BFGS UPDATE
C
         DO 50 I=1,NVAR
            SVEC(I)=GRAD(I)-OLDF(I)
   50    CONTINUE
         DDS=DDOT(NVAR,SVEC,1,D,1)                                      IR0494
C
C  IF DDS IS NEGATIVE, RETENTION OF POSITIVE DEFINITENESS IS NOT
C  GUARANTEED. PRINT A WARNING AND SKIP UPDATE THIS CYCLE.
C
cfrj With the current level shift technique I think the Hessian should
cfrj be allowed to aquire negative eigenvalues. Without updating the
cfrj optimization has the potential of stalling
cfrj     IF(DDS.LT.ZERO) THEN
cfrj        WRITE(6,100)
cfrj        WRITE(6,110)
cfrj        RETURN
cfrj     ENDIF
C
         DDTD=DDOT(NVAR,D,1,TVEC,1)                                     IR0494
C
         DO 70 I=2,NVAR
            DO 60 J=1,I-1
               TEMP= (SVEC(I)*SVEC(J))/DDS - (TVEC(I)*TVEC(J))/DDTD
               HESS(I,J)=HESS(I,J)+TEMP
               HESS(J,I)=HESS(I,J)
   60       CONTINUE
   70    CONTINUE
         DO 75 I=1,NVAR
            TEMP= (SVEC(I)*SVEC(I))/DDS - (TVEC(I)*TVEC(I))/DDTD
            HESS(I,I)=HESS(I,I)+TEMP
 75      CONTINUE
      ENDIF
C
      RETURN
C
   80 FORMAT(/,5X,'HESSIAN IS BEING UPDATED USING THE POWELL UPDATE',/)
   90 FORMAT(/,5X,'HESSIAN IS NOT BEING UPDATED',/)
c 100 FORMAT(5X,'WARNING! HEREDITARY POSITIVE DEFINITENESS ENDANGERED')
c 110 FORMAT(5X,'UPDATE SKIPPED THIS CYCLE')
  120 FORMAT(/,5X,'HESSIAN IS BEING UPDATED USING THE BFGS UPDATE',/)
      END
