      SUBROUTINE CCT4(NO,NU,O1,OEH,OEP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER A,B,E
      LOGICAL SDT1,SDT,Q1,Q2,Q2NHF,Q3,Q4,Q5,FQ,IREST,T32,Q1IT,print
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/NEWOPT/NOPT(6)
      COMMON/NT3T3/NTO3,NTT3,LT3
      COMMON/NVTAP/NTITER,NQ2,NVT
      COMMON/ITRAT/icycle,MAXIT,ICCNV
      COMMON/DMP/IDAMP,DAMP
      common/finerg/erg0,erg1
      common/flags/iflags(100)
      common/timeinfo/timein,timenow,timetot,timenew
      common/rstrt/irest,ccopt(10),t32
      common/activ/noa,nua
      COMMON/IQ1F/iq1f
      common/dift2/difmax
      COMMON /ENERGY/ ECORR(500,2),IXTRLE(500)
      common /rmat/ rmatrix(400)
      common/restart/nresf,nresl
      COMMON/NEWCCSD/NTT2
      COMMON/ENERGIES/ENRES(150)
      COMMON/ENRSDT/ECCSDT
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON/NEWTAPE/NT2T4,NT4INT
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      dimension o1(*),oeh(no),oep(nu)
      COMMON // ICORE(1)
      COMMON/SCFEN/ESCF
      COMMON /IOPOS/ ICRSIZ,ICHSIZ,IOFF(2),LENREC
      COMMON /ISTART/ I0
      DATA ZERO/0.0D+0/,TRESH/1.0D-12/,TWO/2.0D+00/,FOUR/4.0D+00/,
     *     HALF/0.5D+00/,ONE/1.0D+0/
      call ienter(4)
      print=iflags(1).ge.10
      call zclock('cct4    ',0)
      iocc=no-noa
      CALL getrec(20,'JOBARC','SCFENEG ',iintfp,ESCF)
      CLOSE(UNIT=NTITER,STATUS='KEEP')
      IT4=NOPT(5)
      IREST=NOPT(4).GT.0
      SDT1=NOPT(1).GT.0
      SDT=NOPT(1).GT.1
      Q1IT=.TRUE.
      T32=.FALSE.
      ESDTQTOT=ZERO
      CPUSTRT=CPTIME
      NO4U2=NO4*NU2
      NB1=NO+NU
      IF(.NOT.IREST)icycle=1
      call drgetfock(no,nu,o1,oeh,oep)
      call initt1(no,nu,o1)
      irlecyc=0
      istat=ishell('test -f inger')
      open(unit=49,file='inger',form='formatted',status='unknown')
      rewind 49
      if(istat.eq.0)then
         read(49,9001)maxit,idamp,iccnv,damp
	else
	damp=zero
	idamp=100
      write(49,9001)maxit,idamp,iccnv,damp
	call flush(6)
      CLOSE(UNIT=49,STATUS='KEEP')
      endif
      if(irest)call drrest(0,no,nu,o1)
 777  continue
      open(unit=49,file='inger',form='formatted',status='old')
      call timer(1)
      xxx0=timenow
      rewind 49
      read(49,9001)maxit,idamp,iccnv,damp
      CLOSE(UNIT=49,STATUS='KEEP')
 9001 format(3i3,f5.2)
      IF(icycle.LT.IT4.OR.NOPT(6).EQ.9) THEN
      Q1=.FALSE.
      Q2=.FALSE.
      Q2NHF=.FALSE.
      Q3=.FALSE.
      Q4=.FALSE.
      Q5=.FALSE.
      FQ=.FALSE.
      ELSE
      IF (icycle.EQ.IT4) THEN
      Q1=NOPT(1).GE.3
      Q2=NOPT(1).GE.4
      ELSE
      Q2NHF=NOPT(1).eq.5
      Q3=NOPT(1).GE.6
      Q4=NOPT(1).GE.7
      FQ=NOPT(1).EQ.8
      ENDIF
      ENDIF
      IF((IREST.OR.IT4.EQ.0).and.nopt(6).ne.9) THEN
      Q1=NOPT(1).GE.3
      Q2=NOPT(1).GE.4
      Q2NHF=NOPT(1).EQ.5
      Q3=NOPT(1).GE.6
      Q4=NOPT(1).GE.7
      FQ=NOPT(1).EQ.8
      ENDIF
      irlecyc=irlecyc+1
      call drtestd2(no,nu,o1,oeh,oep)
      CALL DRF1INT(NO,NU,O1)
      call zclock('f1int   ',1)
      CALL DRWDEXC(NO,NU,O1)
      call zclock('wdexc   ',1)
      if(Q1)CALL DRT2SEC(NO,NU,O1,OEH,OEP)
      call zclock('t2sec   ',1)
      call drintquat2(no,nu,o1)
      call zclock('intquat2',1)
      call drintrit2(no,nu,o1)
      call zclock('intrit2 ',1)
      call drcct1(no,nu,o1,oeh,oep)
      call zclock('cct1    ',1)
      call drcct2(no,nu,o1,oeh,oep)
      call zclock('cct2    ',1)
      call flush(6)
      IF(.NOT.SDT1)GOTO 111
      IF (Q2NHF)CALL drintqua(NO,NU,o1)
      IF(icycle.GT.1.AND.FQ.and.t32)then
      CALL drvt3int(NO,NU,O1)
      call zclock('vt3int  ',1)
      endif
      CPUTIM=CPTIME
      if(sdt)then
      call drintri(0,no,nu,o1,oeh,oep)
      call zclock('intri   ',1)
      endif
      IF (icycle.GT.1.AND.SDT)then
      CALL DRVEMVT3(0,no,nu,o1)
      call zclock('vemvt3  ',1)
      endif
      IF(SDT1)then
         CALL drt3(NO,NU,O1,OEH,OEP)
      call zclock('t3scr   ',1)
         IF (icycle.GT.1) THEN
            IF(SDT)CALL DRT3WT3(NO,NU,O1,OEH,OEP)
      call zclock('t3wt3   ',1)
            if(q4)call drt3wt4(no,nu,o1)
      call zclock('t3wt4   ',1)
         ENDIF
         IF(SDT)CALL DRT3DEN(NO,NU,o1,OEH,OEP)
      endif
      call flush(6)
      IF (icycle.GT.1.AND.(Q3.AND..NOT.FQ.OR.FQ.AND..NOT.T32))then
      call drvemvt3(1,no,nu,o1)
      call zclock('vemvt3  ',1)
      endif
      if(icycle.gt.1)then
      IF(Q2)CALL DRTRINT3(NO,NU,O1)
      call zclock('trint3  ',1)
      else
      call wrva0(no,nu,o1)
      endif
      call flush(6)
      IF (icycle.GT.1.AND.(Q4)) THEN
      call drvt4kl(no,nu,o1)
      call zclock('vt4kl   ',1)
      call drvt4cd(no,nu,o1)
      call zclock('vt4cd   ',1)
      ENDIF
      IF (Q1.AND.Q1IT) THEN
      call drinth(no,nu,o1)
      call zclock('inth    ',1)
      ENDIF
      if(q2nhf) then
      call drintri(1,no,nu,o1,oeh,oep)
      call zclock('intri2  ',1)
      endif
      call zeroma(o1,1,no2u2)
      IF (.NOT.Q1)GOTO 110
      if(q3)then
      call drvmadd(no,nu,o1)
      call zclock('vmadd   ',1)
      endif
      EQ=ZERO
      EQS=ZERO
      EQD=ZERO
      EQ2=ZERO
      EQ3=ZERO
      EQQ=ZERO
      KKIAS=0
 198  format('ijkl:',4i4,'  kk= ',i6)
      call zeroma(o1,1,no2u2)
      DO 100 I=1,NO
      EI=OEH(I)
      DO 100 J=1,I
      DIJ=EI+OEH(J)
      DO 100 K=1,J
      DIJK=DIJ+OEH(K)
      IF(I.EQ.K)GOTO 100
      DO 101 L=1,K
      DIJKL=DIJK+OEH(L)
      IF(J.EQ.L)GOTO 101
      call timer(1)
      xxx=timenow
      iit=no2u2+1
      call zeroma(o1(iit),1,nu4)
      KKIAS=KKIAS+1
CSSS      write(6,198)i,j,k,l,kkias
      if(nopt(6).lt.0.and.j.le.iocc)goto 95
      IF (icycle.GT.1.AND.Q4) THEN
      IF (FQ.AND.T32) THEN
      call drt4wt32(i,j,k,l,no,nu,o1)
      call zclock('vmadd   ',1)
      ENDIF
      call drt4fht4(i,j,k,l,no,nu,o1)
      call zclock('t4fht4  ',1)
      call drt4fpt4(i,j,k,l,no,nu,o1)
      call zclock('t4fpt4  ',1)
      call drt4ppt4(i,j,k,l,no,nu,o1)
      call zclock('t4ppt4  ',1)
      call drt4hht4(i,j,k,l,no,nu,o1)
      call zclock('t4hht4  ',1)
      call drt4hpt4(i,j,k,l,no,nu,o1)
      call zclock('t4hpt4  ',1)
      ENDIF
      IF(Q1IT)THEN
      call drt4ppt2(i,j,k,l,no,nu,o1)
      call zclock('t4ppt2  ',1)
      call drt4vot2(i,j,k,l,no,nu,o1)
      call zclock('t4vot2  ',1)
      if(icycle.gt.1)then
      call drt4ppt3(i,j,k,l,no,nu,o1)
      call zclock('t4ppt3  ',1)
      endif
      if(iq1f.ne.-3)CALL T4DEN(I,J,K,L,NU,DIJKL,o1(iit),OEP)
      endif
 95   continue
      IF(Q1IT)THEN
      IAS=KKIAS
      IF(icycle.EQ.1.AND.Q2.OR.Q4)
     *CALL WRT4(IAS,NU,o1(iit))
      call drdequa(i,j,k,l,no,nu,o1,oeh,oep)
      call zclock('dequa   ',1)
      ENDIF
      call timer(1)
      if(print)WRITE(6,1000) I,J,K,L,KKIAS,(timenow-xxx)
      OPEN(UNIT=NTITER,FILE='ITER',STATUS='OLD')
      WRITE(NTITER,1009)icycle,I,J,K,L,KKIAS,(CPTIM3-CPTIM2),ESDTQTOT
      CLOSE(UNIT=NTITER,STATUS='KEEP')
 101  CONTINUE
 100  CONTINUE
 110  CONTINUE
 111  CONTINUE
      call drenrcon(no,nu,o1,oeh,oep,esdtqtot,xxx0)
      call zclock('enrcon  ',1)
      ecorr(icycle+1,1)=esdtqtot
      zincrem=ecorr(icycle+1,1)-ecorr(icycle,1)
      IF (ABS(ZINCREM).LT.0.0001D+0)T32=.TRUE.
      IF (ABS(ZINCREM).LT.0.005D+0)Q1IT=.TRUE.
       call drnewrlet3(no,nu,o1)
      call zclock('newrlet3',1)
      if(sdt)then
       call movo3(no,nu,o1)
      endif
      icycle=icycle+1
      call movt1(no,nu,o1)
      call drmovo2(no,nu,o1)
      iiii=nt4
      nt4=not4
      not4=iiii
      CCCON=10.0D0**(-ICCNV)
      call drrest(1,no,nu,o1)
      IF(icycle.GT.MAXIT.OR.(icycle.GT.2.and.dabs(difmax).lt.CCCON).OR.
     *NOPT(6).EQ.8)then
      IF (NOPT(6).EQ.9) THEN
      NOPT(6)=8
      GOTO 777
      ENDIF
      call putrec(20,'JOBARC','RLESHIT ',1500,ecorr)
      call putrec(20,'JOBARC','RMATRIX ',400,rmatrix)
 1052 format('fin.corr:',f17.12,'    fin.tot:',f17.12)
      goto 9999
      endif
      OPEN(UNIT=NTITER,FILE='ITER',STATUS='OLD')
      WRITE(NTITER,1090)icycle,CPUIT,ESDTQTOT
      CLOSE(UNIT=NTITER,STATUS='KEEP')
      goto 777
 1200 FORMAT(1X,' Total cpu :',F15.3)
 1201 FORMAT(1X,' Cpu for one iteration :',F15.3)
 1000 FORMAT('i j k l kk:',5i3,' CPU TIME:',F14.2)
 1090 FORMAT('IT:',i3,' CPU:',F8.2,'  CCSDTQ =',F15.10)
 1009 FORMAT('IT,i j k l kk:',6i3,' CPU:',F8.2,'  CCSDTQ =',F15.10)
9999  continue
      call iexit(4)
      RETURN
      END
      subroutine mixt4(no,nu,o4,t4)
      implicit double precision (a-h,o-z)
      integer a,b,c,d
      dimension o4(nu,nu,nu,nu),t4(nu,nu,nu,nu)
      data damp/0.8d+0/,one/1.0d+0/
      do 10 i=1,no
      do 10 j=1,i
      do 10 k=1,j
      do 10 l=1,k
      if(i.eq.k.or.j.eq.l)goto 10
      kkk=it4(i,j,k,l)
      call rdt4(kkk,nu,o4)
      call rdt4n(kkk,nu,t4)
      do 5 a=1,nu
      do 5 b=1,nu
      do 5 c=1,nu
      do 5 d=1,nu
      t4(a,b,c,d)=damp*o4(a,b,c,d)+(one-damp)*t4(a,b,c,d)
 5    continue
      call wrt4(kkk,nu,t4)
 10   continue
      return
      end
