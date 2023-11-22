      PROGRAM CCSDTQ 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL ILBL,ILAB,DAMP
      LOGICAL OP,IREST
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      common/activ/noa,nua
      COMMON/ORBITALS/NO,NU
      COMMON/ITRAT/IT,MAXIT,ICCNV
      COMMON/DMP/IDAMP,DAMP
      COMMON // ICORE(1)
      COMMON/NN/NH2,NH3,NH4,NP2,NP3,NP4,NHP,NH2P,NH3P,NHP2,NHP3,NH2P2
      COMMON/DWORK/WORK(200)
      common /flags/ iflags(100)
      common/rstrt/irest,ccopt(10),t32
      COMMON/INFO/nocco(2),nvrto(2)
      COMMON /ISTART/ I0,icrsiz
      COMMON/NEWOPT/NOPT(6)
      COMMON/IQ1F/IQ1F
      COMMON/ITQ4/itq4
      COMMON/OPTSD/NOPTSD(6)
      COMMON/NT3T3/NO3,NT3,LT3
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      COMMON/NVTAP/NTITER,NQ2,NVT
      COMMON/NEWLEN/LT2T4,LT4INT,LOU2
      COMMON/NEWT3/NT3INT,LT3INT
      COMMON/LNT/lrecb,lrecw
      common/totmem/lw
      COMMON/PARAM/TRES,TREMO,N1,N2,N3,N4
      common/treshold/tsh
      DATA ZERO/0.0D+0/
      call ienter(1)
      call stamp
      call crapsi(icore,iuhf,0)
      lenwor=(icrsiz-iand(icrsiz,1))/2
      DAMP=ZERO
      IDAMP=0
      itq4=0
      no=nocco(1)
      nu=nvrto(1)
      call setsiz(no,nu)
      LOU2=NO*NU*NU*8
      LT3   =NU*NU*NU*8
      LT3INT=LT3
      LT4INT=LT3
      LL4   =LT3
      LT2T4 =NU*NU*NU*8
      LT4   =NU*NU*NU*NU*8
      do 1 i=1,6
         noptsd(i)=0
 1    continue
      lrecb=32768
       if(iintfp.eq.2)then
      lrecw=2728
      else
      lrecw=2048
      endif
 1189 FORMAT('LT3:',I7,'  LT2T4:',I7,'  LT4:',I8)
 1188 FORMAT('Restart amplitudes below ',f5.4,' will be dropped')
 1006 FORMAT('nopt:',6I3)

CSSS      nopt(1)=iflags(2)-10
      nopt(1)=iflags(2)-31
      nopt(2)=iflags(7)
      maxit=iflags(7)
      iccnv=iflags(4)
      nopt(3)=iflags(4)
      nopt(4)=0
      nopt(5)=0
      nopt(6)=0
      IQ1F=NOPT(1)
      IF (IQ1F.LT.0)NOPT(1)=-IQ1F
      if(nopt(6).eq.9.or.nopt(6).lt.0)then
         read(5,1006)noptsd
      if(nopt(6).lt.0)then
      noa=noptsd(3)
      nua=noptsd(4)
      endif
      endif
      if(nopt(4).eq.2)then
      read(5,*)tsh
      write(6,1188)tsh
      endif
      IREST=NOPT(4).GT.0
      IT=0
      IF(IREST) THEN
         IT=NOPT(4)
         KK=2*(IT/2)
         IF (IT.NE.KK) THEN
            II =NOT4
            NOT4=NT4
            NT4=II
         ENDIF
      ENDIF
      TRES =10.0D0**(-9)
      TREMO=10.0D0**(-9)
      MAXW=LW
      NFMX=NOMX
      call prtinf(no,nu,nopt)
      CALL OPENFI(1)
      LW=LENWOR
      LW=(LW/8)*8
      CALL DRIVER(icore(i0),LW)
      call iexit(1)
      call rcpu
      stop 'all done!'
      END
      BLOCK DATA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ACC/NACC
      COMMON/NEWTAPE/NT2T4,NT4INT
      COMMON/NT3T3/NO3,NT3,LT3
      COMMON/NEWCCSD/NTT2
      COMMON/NEWT3/NT3INT,LT3INT
      COMMON/NEWLEN/LT2T4,LT4INT,LOU2
      COMMON/NEWT4/NT4,NO4,LT4,NALL4,LL4
      COMMON/NVTAP/NTITER,NQ2,NVT
      COMMON/ARPO/NARP,nt1
      COMMON/PAK/INTG,INTR
      common/naniby/nv4
      common/unpak/ntt3,lnt3
c      COMMON/LNT/lrecb,lrecw
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/ETIM/TIME(150),RNAME(150),IENTER(150),tx(200)
      CHARACTER*8 RNAME
      DATA INTG,INTR,ntt3/1,2,3/
c      DATA lrecb,lrecw/32768,2728/
      DATA NT2T4, NT4INT,NT4,   NO4,   NALL4, NT3INT,NO3,   NT3,
     $     NTT2,  NTITER,NQ2,   NVT,   NARP,  NACC,  NV4
     $/    35,    36,    37,    38,    39,    40,    41,    42,
     $     43,    44,    45,    46,    47,    48,    31/
      DATA RNAME /
     *'CCSDTQ  ','CCT1    ','CCT2    ','CCT4    ','DEQUA   ','ENRCON  ',
     *'F1INT   ','INSITU  ','INTH2P4 ','INTH4P2 ','INTQUA  ','INTRI   ',
     *'INTRIH  ','INTRIP  ','INTRIT2 ','MATMUL  ','MOVT1   ','MOVO2   ',
     *'MTRANS  ','MTRSM   ','MTRSMN  ','RDIJKL  ','RDIJKM  ','RDIJMK  ',
     *'RDIJMN  ','RDILKM  ','RDIMJK  ','RDIMKN  ','RDIMNL  ','RDKJML  ',
     *'RDLJKM  ','RDMIJK  ','RDMJKN  ','RDMJNL  ','RDMNKL  ','RDVA    ',
     *'RDVT3   ','RDVT3I  ','RDWT3N  ','RDWT3O  ','RDWT3ONW','SYMT3   ',
     *'TRANMD  ','TRANSQ  ','TRANT3  ','TRINT3  ',
     *'TRINT3ME','TRINT3MN','TRINVTEF','TRINVTME','TRT3ALL ','T1WT3   ',
     *'T2FT3   ','T2SEC   ','T2SET   ','T2WT3   ','T3DEN   ','T3SCR   ',
     *'T3WT3   ','T3WT3HH ','T3WT3HP ','T3WT3PP ','T3WT4   ','T4DEN   ',
     *'T4FHT4  ','T4FPT4  ','T4HHT4  ','T4HPT4VE','T4PPT2IJ','T4PPT3IJ',
     *'T4PPT4  ','T4SQUA  ','T4VOT2IJ','T4WT32ME','T4WT32MN','VECCAD  ',
     *'VECCOP  ','VECMUL  ','VECSUB  ','VEMVT3  ','VMADD   ','VMINUS  ',
     *'VT3EF   ','VT3INT  ','VT3ME   ','VT4CD   ','VT4KL   ','WDEXC   ',
     *'WRGEN   ','WRT4    ','WRVA    ','WRVEM   ','WRVO    ','WRVOE   ',
     *'WRVOEP  ','WRVT3I  ','WRVT3IMN','WRVT3N  ','WRVT3O  ','WRWT3   ',
     *'IPERM   ','IPERM24 ','IPERM4  ','RANK1   ','T2PPT2SQ','T2HPT2SQ',
     *'T2HHT2SQ','T2WWT3  ','        ','110     ','111     ','112     ',
     *'113     ','114     ','115     ','116     ','117     ','118     ',
     *'119     ','120     ','ppppint ','hhhhint ','voeint  ','voint   ',
     *'store   ','scattr  ','mmf     ','zeroma  ','t1wt3ijk','t2wt3ijk',
     *'rdt4    ','        ','        ','        ','        ','        ',
     *'        ','        ','        ','        ','        ','        ',
     *'        ','        ','        ','        ','        ','        ',
     *'        ','        '/
      END
      SUBROUTINE DRIVER(WORK,LW)
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      DIMENSION WORK(1)
      COMMON/ORBITALS/NO,NU
      COMMON/PAK/INTG,INTR
      COMMON/LNT/lrecb,lrecw
      COMMON/ORBINF/NBASIS,NBCOMP
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/newio/nvv,nh4,npp,nvo,nt3,no1,nu1
      common/wpak/nfr(30),nsz(30)
      COMMON/DWORK/ORBEN(200)
      call getrec(20,'JOBARC','NBASTOT ',1,NBASIS)
      call getrec(20,'JOBARC','NBASCOMP',1,NBCOMP)
c      lrecb=32768
c       if(iintfp.eq.2)then
c      lrecw=2728
c      else
c      lrecw=2048
c      endif
      lto=nou+nu2+nu3+no2u2*2
      length=lw-lto
      nvv=33
      lr=8*nu2
      call getrec(20,'JOBARC','SCFEVALA',(NO+NU)*iintfp,ORBEN)
      OPEN(nvv,recl=lr,ACCESS='DIRECT',FORM='UNFORMATTED')
      l34=2*4096/iintfp
c      l34=2*8192/iintfp
c      l34=2*16384/iintfp
      open(unit=34,recl=l34,access='DIRECT',FORM='UNFORMATTED')
      io0=1
      io1=io0+nu2
      io2=io1+nu2
      io3=io2+nu2
      io4=io3+nu2
      io5=io4+nu2
      ito=io5+nu2
      call rdsymv(nu,work,work(io1),work(io2),work(io3),work(io4),
     *work(io5),233)
      inu2=nu2
      if(lrecw.gt.nu2)inu2=lrecw
      i0=1
      i1=i0+inu2
      i2=i1+inu2
      i3=i2+inu2
      i4=i3+inu2
      i5=i4+inu2
      write(6,*)'entering getv4'
      call flush(6)
      call getv4(nu,lrecw,work,work(i1),work(i2),work(i3),work(i4),
     *work(i5))
      write(6,*)'past getv4'
      call flush(6)
      i0=1
      i1=i0+nou3
      call  rhhhh(no,nu,work,work(i1))
      call rrghpp(no,nu,work,work(i1))
      call rldhpp(no,nu,work,work(i1))
      call  rhhhp(no,nu,work,work(i1))
      call  rhppp(no,nu,work,work(i1))
      i2=i1+nu2
      i3=i2+nu2
      call  rpppp(no,nu,work,work(i1),work(i2),work(i3))
      call ro2hppnew(no,nu,work,work(i1),orben,orben(NO+1))
      CALL CC(NO,NU,LW,ORBEN,ORBEN(NO+1),WORK)
      RETURN
      END
      SUBROUTINE OPENFI(ID)
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      COMMON/NN/NO2,Nno3,NnO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/NEWLEN/ LT2T4,LT4INT,LOU2
      COMMON/NEWTAPE/NT2T4,NT4INT
      COMMON/NEWT4/NT4,NO4,LT4,NALL4,LL4
      COMMON/NEWT3/NT3INT,LT3INT
      COMMON/NT3T3/NO3,NT3,LT3
      COMMON/NVTAP/NTITER,NQ2,NVT
      COMMON/NEWCCSD/NTT2
      COMMON/ARPO/NARP,nt1
      COMMON/PAK/INTG,INTR
      COMMON/UNPAK/ntt3,lnt3
      COMMON/LNT/lrecb,lrecw
      common/restart/nresf,nresl
      COMMON/NEWOPT/NOPT(6)
      COMMON/ACC/NACC
      write(6,*)'nt2t4,nt4int,nt4,no4,nall4,nt3int,no3,nt3,ntiter,nq2,
     *nvt,ntt2,narp,nt1,intg,intr,ntt3,nresf,nresl',
     *           nt2t4,nt4int,nt4,no4,nall4,nt3int,no3,nt3,ntiter,nq2,
     *nvt,ntt2,narp,nt1,intg,intr,ntt3,nresf,nresl
      call flush(6)
c      lrecb=8*32768/iintfp
c      lrecw=2*2728/iintfp
      write(6,*)'lrecb=',lrecb
      call flush(6)
      nresf=10
      nresl=16*16192/iintfp
      lnt3=8*nou2
      nnew3=3
      lnrle=8*(nou+no2u2)
      OPEN(unit=nresf,file='restart',RECL=nresl,ACCESS='DIRECT',
     *FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(ntt3,RECL=lnt3,ACCESS='DIRECT',FORM='UNFORMATTED',
     *STATUS='UNKNOWN')
      OPEN(NALL4,RECL=LL4,ACCESS='DIRECT',FORM='UNFORMATTED',
     *STATUS='UNKNOWN')
      OPEN(NTT2,RECL=LT3,ACCESS='DIRECT',FORM='UNFORMATTED')
      OPEN(NT2T4,RECL=LT2T4,ACCESS='DIRECT',FORM='UNFORMATTED')
      OPEN(NARP,RECL=LT3,ACCESS='DIRECT',FORM='UNFORMATTED')
      OPEN(UNIT=NTITER,FILE='ITER',ACCESS='SEQUENTIAL',FORM='FORMATTED')
      OPEN(NACC,RECL=LT3,ACCESS='DIRECT',FORM='UNFORMATTED',
     *     STATUS='UNKNOWN')
      OPEN(INTG,RECL=lrecb,ACCESS='DIRECT',FORM='UNFORMATTED',
     *     STATUS='UNKNOWN')
      OPEN(INTR,RECL=lrecb,ACCESS='DIRECT',FORM='UNFORMATTED',
     *     STATUS='UNKNOWN')
      OPEN(unit=11,FILE='RLE1',RECL=lnrle,ACCESS='DIRECT',
     *     FORM='UNFORMATTED',STATUS='UNKNOWN')
      IF (NOPT(1).EQ.0)RETURN
      OPEN(NT3INT,RECL=LT3INT,ACCESS='DIRECT',FORM='UNFORMATTED')
      OPEN(NT3,RECL=LT3,ACCESS='DIRECT',FORM='UNFORMATTED')
      OPEN(NO3,RECL=LT3,ACCESS='DIRECT',FORM='UNFORMATTED',
     *     STATUS='UNKNOWN')
      OPEN(unit=12,FILE='RLE2',RECL=LT3,ACCESS='DIRECT',
     *     FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(unit=13,FILE='RLE3',RECL=LT3,ACCESS='DIRECT',
     *     FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(unit=14,FILE='RLE4A',RECL=LT4,ACCESS='DIRECT',
     *     FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(unit=15,FILE='RLE4B',RECL=LT4,ACCESS='DIRECT',
     *     FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(NQ2,RECL=LT3,ACCESS='DIRECT',FORM='UNFORMATTED')
      OPEN(NT4INT,RECL=LT4INT,ACCESS='DIRECT',FORM='UNFORMATTED')
      OPEN(NT4,RECL=LT4,ACCESS='DIRECT',FORM='UNFORMATTED',
     *STATUS='UNKNOWN')
      OPEN(NO4,RECL=LT4,ACCESS='DIRECT',FORM='UNFORMATTED',
     *     STATUS='UNKNOWN')
      OPEN(NVT,RECL=LOU2,ACCESS='DIRECT',FORM='UNFORMATTED')
      nvtq4=4
      OPEN(nvtq4,RECL=LOU2,ACCESS='DIRECT',FORM='UNFORMATTED')
      RETURN
      END
      SUBROUTINE CC(NO,NU,NLT4,OEH,OEP,T4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER   A,B,C,D
      COMMON/ITRAT/IT,MAXIT,ICCNV
      common/wpak/nfr(30),nsz(30)
      COMMON/NEWOPT/NOPT(6)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION OEH(NO),OEP(NU),T4(NLT4)
      CPUTIM=CPTIME
      CALL ZEROMA(T4,1,NLT4)
      CALL CCT4(NO,NU,T4,oeh,oep)
      RETURN
      END
      SUBROUTINE PRTINF(NH,NP,NOPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*8 CCOPT
      COMMON/NN/NH2,NH3,NH4,NP2,NP3,NP4,NHP,NH2P,NH3P,NHP2,NHP3,NH2P2
      common/rstrt/irest,ccopt(10),t32
      DIMENSION NOPT(6)
      DATA CCOPT/'     CCD','  CCSD  ','  SDT-1 ',' FCCSDT ',' SDTQ-1 ',
     *' SDTQ-2 ','SDTQ2NHF',' SDTQ-3 ',' SDTQ-4 ','  FSDTQ '/
      NH2=NH*NH
      NH3=NH*NH2
      NH4=NH*NH3
      NP2=NP*NP
      NP3=NP*NP2
      NP4=NP*NP3
      NHP=NH*NP
      NH2P=NH2*NP
      NH3P=NH3*NP
      NHP2=NH*NP2
      NHP3=NH*NP3
      NH2P2=NH2*NP2
      WRITE(6,106)NH
      WRITE(6,107)NP
      K=NOPT(1)+2
      WRITE(6,108)CCOPT(K)
 101  FORMAT('NUMBER OF FROZEN CORE ORBITALS   ',I4) 
 102  FORMAT('FIRST CORRELATED OCCUPIED ORBITAL',I4)
 103  FORMAT('LAST  CORRELATED OCCUPIED ORBITAL',I4)
 104  FORMAT('FIRST CORRELATED VIRTUAL ORBITAL ',I4)
 105  FORMAT('LAST  CORRELATED VIRTUAL ORBITAL ',I4)
 106  FORMAT('NUMBER OF CORRELATED OCCUPIEDS   ',I4)
 107  FORMAT('NUMBER OF CORRELATED VIRTUALS    ',I4)
 108  FORMAT('CC OPTIONS:  ',A8)
      RETURN
      END
      SUBROUTINE STAMP 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      WRITE(6,*)'    **************************************************'
      WRITE(6,1000)
 1000 FORMAT(/T10,'COUPLED CLUSTER PROGRAM WITH FULL INCLUSION OF'/,
     *       T10,'SINGLES, DOUBLES, TRIPLES AND QUADRUPLES'//,
     *       T10,'  ******* CLOSED SHELL VERSION *******'////)
      write(6,1100)
 1100 FORMAT(T15,'Authors: Stanislaw A. Kucharski* and'/,
     *       T24,'Rodney J. Bartlett'//,
     *       T20,'Quantum Theory Project'/,
     *       T20,'University of Florida, Gainesville,USA'//)
     *
      write(6,1200)
      write(6,1300)
 1200 FORMAT(T5,'*Permanent address: Institute of Chemistry',/
     *           T25,'Silesian University,Katowice,Poland'//)
 1300 FORMAT(T5,'Relevant papers: S.A.Kucharski,R.J.Bartlett,'/,
     *       T22,'J.Chem.Phys.,97,4282(1992)',//T22,'S.A.Kucharski,R.J.Ba
     *rtlett,',/T22,'Theoret.Chim.Acta,80,387(1991)'//,T22,
     *'S.A.Kucharski, R.J.Bartlett,'/,T22,'Chem.Phys.Lett.,158,550(1989)
     *',//T22,'S.A.Kucharski,R.J.Bartlett,',/T22,'Chem.Phys.Lett.,206,574
     *(1993)')
      RETURN
      END
      SUBROUTINE rcpu
      IMPLICIT REAL*8 (a-h,o-z)
      COMMON /etim/tm(150),rutine(150),LUSE(150),ttim(200)
c      COMMON /timex/ttim(200)
      WRITE (6,100)
      WRITE (6,110)
      DO 1 K = 1,150
      IF (LUSE(K).NE.0) WRITE (6,120)k,rutine(K),LUSE(K),ttim(K)
   1  CONTINUE
      RETURN

 100  FORMAT (//1X,'CPU USED BY ROUTINES:')
 120  FORMAT (1X,I3,4X,A8,8X,I8,8X,F10.3)
 110  FORMAT (1X,'SUBROUTINE         No. of calls      ',' CPU',
     *       '(SECONDS)'/)
      END
      subroutine setsiz(nh,np)
      implicit double precision(a-h,o-z)
      COMMON/NN/NH2,NH3,NH4,NP2,NP3,NP4,NHP,NH2P,NH3P,NHP2,NHP3,NH2P2
      NH2=NH*NH
      NH3=NH*NH2
      NH4=NH*NH3
      NP2=NP*NP
      NP3=NP*NP2
      NP4=NP*NP3
      NHP=NH*NP
      NH2P=NH2*NP
      NH3P=NH3*NP
      NHP2=NH*NP2
      NHP3=NH*NP3
      NH2P2=NH2*NP2
      return
      end
