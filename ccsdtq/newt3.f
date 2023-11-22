      subroutine drintri(ic,no,nu,t,eh,ep)
      implicit double precision (a-h,o-z)
      logical print
      common/flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/totmem/nmem
      common/itrat/icycle,mx,icn
      dimension t(*),eh(no),ep(nu)
      print=iflags(1).gt.10
      i1=1         !ti
      i2=i1+nu3    !o2
      i3=i2+no2u2  !o2v
      i4=i3+no2u2  !vm
      i5=i4+no3u   !vm1
      i6=i5+no3u   !vm2
      i7=i6+no3u   !ve
      i8=i7+nu4   !ve1
      itot=i8+nou3
      if(icycle.eq.1.and.print)write(6,99)nmem,itot
 99   format('Space usage in intri: available - ',i8,'   used - ',i8)
      call intri(ic,no,nu,t,t(i2),t(i3),t(i4),t(i5),t(i6),t(i7),t(i8),
     &eh,ep)
      return
      end
      SUBROUTINE INTRI(IC,NO,NU,TI,O2,O2V,VM,VM1,VM2,VE,VE1,OEH,OEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TI(1),O2(1),O2V(1),VM(1),VM1(1),VM2(1),VE(1),VE1(1),
     *OEH(1),OEP(1)
      call ienter(12)
      CALL INTRIP(IC,NO,NU,TI,O2,O2V,VM,VE,VE1)
      CALL INTRIH(IC,NO,NU,TI,O2,O2V,VM,VM1,VE,VM2)
      call iexit(12)
      RETURN
      END
      SUBROUTINE INTRIH(IC,NO,NU,TI,FHP,T2HP,VM,VM1,VE,VM2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SDTQ2,SDTQ3
      COMMON/NEWOPT/NOPT(6)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION TI(NU,NU,NO),T2HP(NO,NU,NO,NU),FHP(NO,NU),
     *VE(NU,NU,NO,NU),VM(NO,NO,NO,NU),VM1(NO,NO,NO,NU),
     *VM2(NO,NU,NO,NO)
      DATA ZERO/0.0D+0/,TWO/2.0D+00/,HALF/0.5D+00/
      call ienter(13)
      SDTQ2=NOPT(1).EQ.4
      SDTQ3=NOPT(1).EQ.5
      I10=10
 10   CONTINUE
      IF (SDTQ2)THEN
      J10=I10
      I10=0
      IF(J10.EQ.0)I10=10
      ENDIF
      IF (SDTQ3.AND.IC.EQ.1)I10=0
      CALL RO2HPH(0,NO,NU,TI,T2HP)
      CALL RDVEM3(I10,NO,NU,TI,VE)
      CALL TRANMD(VE,NU,NU,NO,NU,14)
      I11=I10+1
      CALL RDVEM4(I11,NU,NO,TI,VM)
      CALL RDVEM4(1,NU,NO,TI,VM1)
      CALL VECMUL(T2HP,NO2U2,TWO)
      CALL MATMUL(VM,T2HP,VM1,NO2,NOU,NOU,0,0)
      CALL VECMUL(T2HP,NO2U2,HALF)
      CALL TRANMD(T2HP,NO,NU,NO,NU,24)
      CALL MATMUL(VM,T2HP,VM1,NO2,NOU,NOU,0,1)
      CALL TRANMD(T2HP,NO,NU,NO,NU,24)
      CALL TRANMD(VM,NO,NO,NO,NU,23)
      CALL MATMUL(VM,T2HP,VM1,NO2,NOU,NOU,0,1)
      CALL TRANMD(T2HP,NO,NU,NO,NU,24)
      CALL TRANMD(VM1,NO,NO,NO,NU,13)
      CALL MATMUL(VM,T2HP,VM1,NO2,NOU,NOU,0,1)
      CALL TRANMD(T2HP,NO,NU,NO,NU,24)
      CALL TRANMD(VM,NO,NO,NO,NU,23)
      CALL TRANMD(VM1,NO,NO,NO,NU,231)
      CALL RO2HHP(0,NO,NU,TI,T2HP)
      CALL MATMUL(T2HP,VE,VM1,NO2,NOU,NU2,0,0)
      CALL TRANMD(VM1,NO,NO,NO,NU,231)
      IF (SDTQ2.AND.(I10.EQ.0)) THEN
      call wrvemq2(1,nu,no,ti,vm1)
      GOTO 10
      ENDIF
      IF(SDTQ3.AND.IC.EQ.1) GOTO 11
      CALL RO2HPP(0,NO,NU,TI,T2HP)
      CALL TRANSQ(T2HP,NOU)
      CALL RDFPH(NO,NU,FHP,TI)
      CALL MATMUL(FHP,T2HP,VM1,NO,NO2U,NU,0,0)
 11   CONTINUE
      CALL RDT1PH(NO,NU,FHP,TI)
      IF(SDTQ3.AND.IC.EQ.1)THEN
      CALL RO2HPP(1,NO,NU,TI,T2HP)
      ELSE
      CALL ROVOEP(1,NO,NU,TI,T2HP)
      ENDIF
      CALL TRANSQ(T2HP,NOU)
      CALL TRANMD(VM1,NO,NO,NO,NU,12)
      CALL MATMUL(FHP,T2HP,VM1,NO,NO2U,NU,0,0)
      CALL TRANMD(VM1,NO,NO,NO,NU,231)
      IF(SDTQ3.AND.IC.EQ.1)THEN
      CALL RO2HPP(2,NO,NU,TI,T2HP)
      ELSE
      CALL ROVOEP(2,NO,NU,TI,T2HP)
      ENDIF
      CALL TRANSQ(T2HP,NOU)
      CALL MATMUL(FHP,T2HP,VM1,NO,NO2U,NU,0,0)
       CALL TRANMD(VM1,NO,NO,NO,NU,312)
      CALL RDOV4(1,NU,NO,TI,VM)
      CALL TRANMD(VM1,NO,NO,NO,NU,12)
      CALL MTRANS(VM,NO,7)
      CALL MATMUL(VM,FHP,VM1,NO3,NU,NO,0,1)
 200  CONTINUE
      IF(SDTQ3.AND.IC.EQ.1) GO TO 12
      CALL RO2HPP(0,NO,NU,TI,T2HP)
      CALL TRANSQ(T2HP,NOU)
      CALL RDFPH(NO,NU,FHP,TI)
      CALL MATMUL(FHP,T2HP,VM1,NO,NO2U,NU,0,1)
 12    CONTINUE
      IW11=11
      IF(SDTQ3.AND.IC.EQ.1)IW11=1
      call wrvemq2(iw11,nu,no,ti,vm1)
      call iexit(13)
      RETURN
      END
      SUBROUTINE INTRIP(IC,NO,NU,TI,FPH,T2PH,VM,VE,VE1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SDTQ2,SDTQ3
      COMMON/NEWOPT/NOPT(6)
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION TI(NU,NU,NO),FPH(NU,NO),T2PH(NU,NO,NU,NO),
     *VE(NU,NU,NU,NO),VM(NO,NO,NU,NO),VE1(NU,NU,NU,NO)
      DATA ZERO/0.0D+0/,TWO/2.0D+00/,HALF/0.5D+00/
      call ienter(14)
      SDTQ2=NOPT(1).EQ.4
      SDTQ3=NOPT(1).EQ.5
      I10=10
 10   CONTINUE
      IF (SDTQ2)THEN
      J10=I10
      I10=0
      IF(J10.EQ.0)I10=10
      ENDIF
      IF(IC.EQ.1.AND.SDTQ3)I10=0
      CALL RO2PHP(0,NO,NU,TI,T2PH)
      CALL RDVEM4(I10,NO,NU,TI,VE)
      CALL RDVEM4(0,NO,NU,TI,VE1)
      CALL VECMUL(T2PH,NO2U2,TWO)
      CALL MATMUL(VE,T2PH,VE1,NU2,NOU,NOU,0,0)
      CALL VECMUL(T2PH,NO2U2,HALF)
      CALL TRANMD(T2PH,NU,NO,NU,NO,24)
      CALL MATMUL(VE,T2PH,VE1,NU2,NOU,NOU,0,1)
      CALL TRANMD(T2PH,NU,NO,NU,NO,24)
      CALL TRANMD(VE,NU,NU,NU,NO,23)
      CALL MATMUL(VE,T2PH,VE1,NU2,NOU,NOU,0,1)
      CALL TRANMD(T2PH,NU,NO,NU,NO,24)
      CALL TRANMD(VE1,NU,NU,NU,NO,13)
      CALL MATMUL(VE,T2PH,VE1,NU2,NOU,NOU,0,1)
      CALL TRANMD(VE,NU,NU,NU,NO,23)
      CALL TRANMD(T2PH,NU,NO,NU,NO,24)
      CALL TRANMD(VE1,NU,NU,NU,NO,231)
      CALL RO2PPH(0,NO,NU,TI,T2PH)
      I11=I10+1
      CALL RDVEM3(I11,NU,NO,TI,VE)
      CALL TRANMD(VE,NO,NO,NU,NO,14)
      CALL MATMUL(T2PH,VE,VE1,NU2,NOU,NO2,0,0)
      CALL TRANMD(VE1,NU,NU,NU,NO,231)
      IF (SDTQ2.AND.(I10.EQ.0)) THEN
      call wrvemq2(0,no,nu,ti,ve1)
      GOTO 10
      ENDIF
      IF (SDTQ3.AND.IC.EQ.1) GOTO 11
      CALL RO2HPP(0,NO,NU,TI,T2PH)
      CALL RDFPH(NO,NU,TI,FPH)
      CALL MATMUL(FPH,T2PH,VE1,NU,NOU2,NO,0,1)
 11   CONTINUE
      CALL RDT1PH(NO,NU,TI,FPH)
      IF(SDTQ3.AND.IC.EQ.1) THEN
      CALL RO2HPP(1,NO,NU,TI,T2PH)
      ELSE
      CALL ROVOEP(1,NO,NU,TI,T2PH)
      ENDIF
      CALL TRANMD(T2PH,NO,NU,NU,NO,14)
      CALL TRANMD(T2PH,NO,NU,NU,NO,23)
      CALL TRANMD(VE1,NU,NU,NU,NO,12)
      CALL MATMUL(FPH,T2PH,VE1,NU,NOU2,NO,0,1)
      CALL TRANMD(VE1,NU,NU,NU,NO,231)
      IF(SDTQ3.AND.IC.EQ.1) THEN
      CALL RO2HPP(2,NO,NU,TI,T2PH)
      ELSE
      CALL ROVOEP(3,NO,NU,TI,T2PH)
      ENDIF
      CALL TRANMD(T2PH,NO,NU,NU,NO,23)
      CALL TRANMD(T2PH,NO,NU,NU,NO,14)
      CALL MATMUL(FPH,T2PH,VE1,NU,NOU2,NO,0,1)
      CALL TRANMD(VE1,NU,NU,NU,NO,312)
      CALL RDOV4(0,NO,NU,TI,VE)
      CALL TRANMD(VE1,NU,NU,NU,NO,12)
      CALL MATMUL(VE,FPH,VE1,NU3,NO,NU,0,0)
      IW10=10
      IF(SDTQ3.AND.IC.EQ.1)IW10=0
      call wrvemq2(iw10,no,nu,ti,ve1)
      call iexit(14)
      RETURN
      END
      subroutine drt3(nh,np,t,eh,ep)
      implicit double precision(a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/wpak/nfr(30),nsz(30)
      common/totmem/nlmem
      dimension eh(nh),ep(np),t(1)
      i1=1         
      i2=i1+nou    
      i3=i2+no2u2  
      i4=i3+no3u   
      i5=i4+nou3   
      i6=i5+nu3    
      i7=i6+nu3    
      i8=i7+no3u   
      i9=i8+no2u2  
      i10=i9+no2u2 
      i11=i10+nsz(6)/2+1  
      itot=i11+np
      if(itot.gt.nlmem)then
         write(6,*)'insuficient memory for t3 amplitudes'
         write(6,*)'required:',itot,'   available:',nlmem
         endif
      call t3scr(nh,np,t,t(i2),t(i3),t(i4),t(i5),t(i6),t(i7),t(i8),
     &t(i9),t(i10),t(i11),eh,ep)
      return
      end
      SUBROUTINE T3SCR(NO,NU,T1,T2PP,VM,VEpak,VT3M,t3,vm1,t2,voe,
     &ive,di,OEH,OEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,E,F
      LOGICAL SDT
      DIMENSION t1(nu,no),T2PP(NU,NU,NO,NO),VM(NO,NU,NO,NO),voe(1),
     &VT3M(nu,nu,nu),t3(1),t2(no,nu,nu,no),vepak(nu,nu,nu,no),ive(1),
     &di(1),OEH(NO),OEP(NU)
      COMMON/ITRAT/IT,MAXIT,ICCNV
      COMMON/NEWTAPE/NT2T4,NT4INT
      COMMON/RESLTS/CMP(30)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON/UNLIN/UN(15)
      common/wpak/nfr(30),nsz(30)
      common/isizes/isi(50)
      common/pak/intg,intr
      COMMON/unpak/nt3,lnt3
      COMMON/NEWCCSD/NTT2
      COMMON/NEWOPT/NOPT(6)
      DATA ZERO/0.0D+0/,HALF/0.5D+0/,TWO/2.0D+0/,FOUR/4.0D+0/,
     &EIGHT/8.0D+0/
      call ienter(58)
      SDT=NOPT(1).GT.1
      NO3U3=NO3*NU3
      call zeroma(t2,1,no2u2)
      call zeroma(t1,1,nou)
      call ro2pph(1,no,nu,vt3m,voe)
      IF ((.NOT.SDT).OR.(IT.EQ.1))THEN
      CALL RDVEM4(0,NO,NU,VT3M,VEpak)
      call tranmd(vepak,nu,nu,nu,no,23)
      CALL RDVEM2(1,NU,NO,VT3M,VM)
      call veccop(no3u,vm1,vm)
      call insitu(no,nu,no,no,vt3m,vm1,12)
      ELSE
      CALL RDVEM4A(10,NO,NU,VT3M,VEpak)
      call tranmd(vepak,nu,nu,nu,no,23)
      CALL RDVEM2A(11,NU,NO,VT3M,VM)
      ENDIF
      CALL TRANMD(VM,NO,NU,NO,NO,13)
      CALL TRANMD(VM1,Nu,No,No,No,23)
      CALL RO2PPH(0,NO,NU,VT3M,T2PP)
      CALL TRANMD(T2PP,NU,NU,NO,NO,34)
      E4T  =ZERO
      E5TDS=ZERO
      E6TDS=ZERO
      O4T2 =ZERO
      ETEMP=ZERO
      KK=0
         ioff=0
         nnve=nsz(6)
         if(sdt.and.it.gt.1)then
         nnve=nsz(16)
         ioff=no+1
         endif
      DO 351 I=1,NO
         i1=i-1
      DO 351 J=1,I1
         j1=j-1
      DO 351 K=1,J1
         kk=it3(i,j,k)
         call ienter(131)
      CALL MATMUL(T2PP(1,1,1,I),VM(1,1,K,J),VT3M,NU2,NU,NO,1,1)
      call trant3(vt3m,nu,1)
      CALL MATMUL(T2PP(1,1,1,I),VM(1,1,J,K),VT3M,NU2,NU,NO,0,1)
      call trant3(vt3m,nu,4)
      CALL MATMUL(T2PP(1,1,1,J),VM(1,1,K,I),VT3M,NU2,NU,NO,0,1)
      call trant3(vt3m,nu,1)
      CALL MATMUL(T2PP(1,1,1,J),VM(1,1,I,K),VT3M,NU2,NU,NO,0,1)
      call trant3(vt3m,nu,5)
      CALL MATMUL(T2PP(1,1,1,K),VM(1,1,J,I),VT3M,NU2,NU,NO,0,1)
      call trant3(vt3m,nu,1)
      CALL MATMUL(T2PP(1,1,1,K),VM(1,1,I,J),VT3M,NU2,NU,NO,0,1)
         call trant3(vt3m,nu,4)
         call iexit(131)
         call ienter(132)
      CALL MATMUL(T2PP(1,1,J,I),VEpak(1,1,1,K),VT3M,NU,NU2,NU,0,0)
      call trant3(vt3m,nu,1)
      CALL MATMUL(T2PP(1,1,K,I),VEpak(1,1,1,J),VT3M,NU,NU2,NU,0,0)
      call trant3(vt3m,nu,5)
      CALL MATMUL(T2PP(1,1,I,J),VEpak(1,1,1,K),VT3M,NU,NU2,NU,0,0)
      call trant3(vt3m,nu,1)
      CALL MATMUL(T2PP(1,1,K,J),VEpak(1,1,1,I),VT3M,NU,NU2,NU,0,0)
      call trant3(vt3m,nu,4)
      CALL MATMUL(T2PP(1,1,I,K),VEpak(1,1,1,J),VT3M,NU,NU2,NU,0,0)
      call trant3(vt3m,nu,1)
      CALL MATMUL(T2PP(1,1,J,K),VEpak(1,1,1,I),VT3M,NU,NU2,NU,0,0)
      call trant3(vt3m,nu,5)
      do 250 ll=1,nu
         vt3m(ll,ll,ll)=zero
 250  continue
      if(sdt)then
      CALL WRVT3N(KK,NU,VT3M)
      call iexit(132)
      goto 351
      endif
      deh=oeh(i)+oeh(j)+oeh(k)
      call adt3den(nu,deh,vt3m,oep)
      call tranmd(vepak,nu,nu,nu,no,23)
      call drt2wt3ijk(j,k,i,no,nu,t2,vt3m,vm1,t3,vepak,ive)
      call tranmd(vepak,nu,nu,nu,no,23)
      call drt1wt3ijk(i,j,k,no,nu,t1,voe,vt3m,t3)
 351  CONTINUE
      DO 352 I=1,NO
         i1=i-1
      DO 352 J=1,I1
      KK=it3(i,j,j)
      call ienter(131)
      CALL MATMUL(T2PP(1,1,1,I),VM(1,1,j,J),VT3M,NU2,NU,NO,1,1)
      call trant3(vt3m,nu,2)
      CALL MATMUL(T2PP(1,1,1,j),VM(1,1,j,i),vT3m,NU2,NU,NO,0,1)
      call trant3(vt3m,nu,1)
      CALL MATMUL(T2PP(1,1,1,J),VM(1,1,I,j),vT3m,NU2,NU,NO,0,1)
      call trant3(vt3m,nu,4)
      call iexit(131)
      call ienter(132)
      CALL MATMUL(T2PP(1,1,J,I),VEpak(1,1,1,j),VT3M,NU,NU2,NU,0,0)
      call trant3(vt3m,nu,3)
      CALL MATMUL(T2PP(1,1,I,J),VEpak(1,1,1,j),VT3M,NU,NU2,NU,0,0)
      call trant3(vt3m,nu,1)
      CALL MATMUL(T2PP(1,1,j,J),VEpak(1,1,1,I),VT3M,NU,NU2,NU,0,0)
         call iexit(132)
      call trant3(vt3m,nu,2)
      call symt311(vt3m,nu,23)
      if(sdt)then
      do 251 ll=1,nu
         vt3m(ll,ll,ll)=zero
 251  continue
      CALL WRVT3N(KK,NU,VT3M)
      goto 352
      endif
      deh=oeh(i)+oeh(j)+oeh(j)
      call adt3den(nu,deh,vt3m,oep)
      call tranmd(vepak,nu,nu,nu,no,23)
      call drt2wt3ij(i,j,no,nu,t2,vt3m,vm1,t3,vepak,ive)
      call tranmd(vepak,nu,nu,nu,no,23)
      call drt1wt3ij(i,j,no,nu,t1,voe,vt3m,t3)
 352  CONTINUE
      DO 353 J=1,no
         j1=j-1
      DO 353 K=1,J1
      KK=it3(j,j,k)
         call ienter(131)
      CALL MATMUL(T2PP(1,1,1,j),VM(1,1,K,J),VT3M,NU2,NU,NO,1,1)
      call trant3(vt3m,nu,1)
      CALL MATMUL(T2PP(1,1,1,j),VM(1,1,J,K),VT3M,NU2,NU,NO,0,1)
      call trant3(vt3m,nu,2)
      CALL MATMUL(T2PP(1,1,1,K),VM(1,1,J,j),VT3M,NU2,NU,NO,0,1)
      call trant3(vt3m,nu,2)
         call iexit(131)
         call ienter(132)
      CALL MATMUL(T2PP(1,1,J,j),VEpak(1,1,1,K),VT3M,NU,NU2,NU,0,0)
      call trant3(vt3m,nu,1)
      CALL MATMUL(T2PP(1,1,K,j),VEpak(1,1,1,J),VT3M,NU,NU2,NU,0,0)
      call trant3(vt3m,nu,3)
      CALL MATMUL(T2PP(1,1,j,K),VEpak(1,1,1,J),VT3M,NU,NU2,NU,0,0)
      call trant3(vt3m,nu,3)
      call iexit(132)
      call symt311(vt3m,nu,12)
      if(sdt)then
         do 252 ll=1,nu
            vt3m(ll,ll,ll)=zero
 252     continue
         CALL WRVT3N(KK,NU,VT3M)
         goto 353
      endif
      deh=oeh(j)+oeh(j)+oeh(k)
      call adt3den(nu,deh,vt3m,oep)
      call tranmd(vepak,nu,nu,nu,no,23)
      call drt2wt3jm(j,k,no,nu,t2,vt3m,vm1,t3,vepak,ive)
      call tranmd(vepak,nu,nu,nu,no,23)
      call drt1wt3jk(j,k,no,nu,t1,voe,vt3m,t3)
 353  CONTINUE
      if (sdt)goto 999
      call insitu(nu,nu,no,no,t3,t2,13)
      do 109 i=1,no
         call wr(nt3,i+1,nou2,t2(1,1,1,i))
 109  continue
      call trt1(nu,no,t3,t1)
      call wr(nt3,1,nou,t1)
 999  continue
      call iexit(58)
      RETURN
      END
      subroutine adt3den(nu,deh,t3,ep)
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      integer a,b,c
      dimension t3(nu,nu,nu),ep(nu)
      do 10 a=1,nu
         do 10 b=1,nu
            do 10 c=1,nu
               den=deh-ep(a)-ep(b)-ep(c)
               t3(a,b,c)=t3(a,b,c)/den
 10         continue
            return
            end
      subroutine drt3wt4(no,nu,t)
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      logical print
      common/flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension t(*)
      common/totmem/mem
      common/itrat/icycle,mx,icn
      print=iflags(1).gt.10
      nou4=no*nu4
      i1=1        !fph
      i2=i1+no3u  !vm
      i3=i2+nou3  !ve
      i4=i3+nou3  !t3
      i5=i4+nu3   !o4
      i6=i5+nou4  !t4n
      it=i6+nu4
      if(icycle.eq.2.and.print)write(6,99)mem,it
 99   format('Space usage in t3wt4: available - ',i8,'   used - ',i8)
      call t3wt4(no,nu,t(i1),t(i2),t(i3),t(i4),t(i5),t(i6))
      return
      end
      SUBROUTINE T3WT4(NO,NU,FPH,VM,VE,T3,O4,T4N)
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      LOGICAL IEJ,JEK,NE1,NE2
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION O4(NU,NU,NU,NU,NO),FPH(NU,NO),VE(NU,NU,NO,NU),
     *T3(NU,NU,NU), T4N(NU,NU,NU,NU),VM(NU,NO,NO,NO)
      data half/0.5d+0/
      call ienter(63)
      CALL RDFPH(NO,NU,T4N,FPH)
      CALL RDVEM3(10,NO,NU,T4N,VE)
      CALL TRANMD(VE,NU,NU,NO,NU,14)
      CALL TRANMD(VE,NU,NU,NO,NU,12)
      CALL RDVEM1(11,NU,NO,T4N,VM)
      CALL TRANMD(VM,NU,NO,NO,NO,34)
      CALL TRANMD(VM,NU,NO,NO,NO,23)
      KKK=0
      DO 1 I=1,NO
      DO 1 J=1,I
      IEJ=I.EQ.J
      DO 1 K=1,J
      JEK=J.EQ.K
      IF (I.EQ.K) GOTO 1
      KKK=KKK+1
      CALL RDVT3N(KKK,NU,T3)
      DO 5 M=1,NO
      IF (IEJ.AND.(J.EQ.K.OR.J.EQ.M).OR.JEK.AND.K.EQ.M)THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GO TO 5
      ENDIF
      CALL RDIJMNT3(I,J,K,M,NO,NU,T4N,O4(1,1,1,1,M))
 5    CONTINUE
      CALL MATMUL(O4,FPH,T3,NU3,1, NOU, 0,0)
      if(j.eq.k)call vecmul(t3,nu3,half)
      CALL MATMUL(O4,VE, T3,NU2,NU,NOU2,0,0)
      CALL MTRSM(O4,NO,NU,7)
      if(j.eq.k)then
         call symt311(t3,nu,23)
         goto 111
      endif
      CALL TRANT3(T3,NU,1)
      if(i.eq.j)call vecmul(t3,nu3,half)
      CALL MATMUL(O4,VE, T3,NU2,NU,NOU2,0,0)
 111  continue
      if(i.eq.j)then
         call symt311(t3,nu,13)
         goto 112
      endif
      CALL TRANT3(T3,NU,3)
      CALL MTRSM(O4,NO,NU,8)
      CALL MATMUL(O4,VE, T3,NU2,NU,NOU2,0,0)
 112  continue
      CALL TRANT3(T3,NU,4)
      if(j.eq.k)call vecmul(t3,nu3,half)
      DO 101 N=1,NO
      NE1=I.EQ.N
      NE2=J.EQ.N
      DO 10 M=1,NO
      IF (IEJ.AND.(J.EQ.M.OR.NE2).OR.M.EQ.N.AND.(NE1.OR.NE2))THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GO TO 10
      ENDIF
      CALL RDIJMNt3(I,J,M,N,NO,NU,T4N,O4(1,1,1,1,M))
 10   CONTINUE
      CALL MATMUL(O4,VM(1,1,N,K),T3,NU3,1,NOU,0,1)
 101  CONTINUE
      if(j.eq.k)then
         call symt311(t3,nu,23)
         goto 1111
      endif
      if(i.eq.j)call vecmul(t3,nu3,half)
      CALL TRANT3(T3,NU,1)
      DO 121 N=1,NO
      NE1=I.EQ.N
      NE2=K.EQ.N
      DO 12 M=1,NO
      IF(M.EQ.N.AND.(NE1.OR.NE2))THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GO TO 12
      ENDIF
      CALL RDIJMNT3(I,K,M,N,NO,NU,T4N,O4(1,1,1,1,M))
 12   CONTINUE
      CALL MATMUL(O4,VM(1,1,N,J),T3,NU3,1,NOU,0,1)
 121  CONTINUE
 1111 continue
      if(i.eq.j)then
         call symt311(t3,nu,13)
         goto 1112
      endif
      CALL TRANT3(T3,NU,3)
      DO 141 N=1,NO
      NE1=N.EQ.J
      NE2=N.EQ.K
      DO 14 M=1,NO
      IF (JEK.AND.(J.EQ.M.OR.NE1).OR.M.EQ.N.AND.(NE1.OR.NE2))THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GO TO 14
      ENDIF
      CALL RDIJMNT3(J,K,M,N,NO,NU,T4N,O4(1,1,1,1,M))
 14   CONTINUE
      CALL MATMUL(O4,VM(1,1,N,I),T3,NU3,1,NOU,0,1)
 141  CONTINUE
 1112 continue
      CALL TRANT3(T3,NU,4)
      CALL WRVT3N(KKK,NU,T3)
 1    CONTINUE
      call iexit(63)
      RETURN
      END
      subroutine drvmadd(no,nu,t)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension t(*)
      i1=1          !ti
      i2=i1+nou2    !t1
      i3=i2+nou     !t2
      i4=i3+no2u2   !vm
      it=i4+no3u
      call vmadd(no,nu,t(i1),t(i2),t(i3),t(i4))
      return
      end
      SUBROUTINE VMADD(NO,NU,TI,T1,T2,VM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T1(NO,NU),TI(NU,NU,NO),T2(NO,NU,NU,NO),VM(NO,NO,NO,NU)
      call ienter(81)
      NOU=NO*NU
      NO2U=NOU*NO
      CALL RDVEM4A(11,NU,NO,TI,VM)
      call rdfph(no,nu,T1,TI)
      call ro2hpp(0,no,nu,TI,T2)
      call transq(T2,NOU)
      CALL MATMUL(T1,T2,VM,NO,NO2U,NU,0,0)
      CALL WRVEMQ2(11,NU,NO,TI,VM)
      call iexit(81)
      RETURN
      END
      subroutine drvemvt3(isel,no,nu,t)
      implicit double precision (a-h,o-z)
      logical print
      common/flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/itrat/icycle,mx,icn
      common/totmem/nmem
      dimension t(*)
      print=iflags(1).gt.10
      i1=1         !ti
      i2=i1+nu3    !voe
      i3=i2+no2u2  !t3
      i4=i3+nou3    !vm
      i5=i4+no3u   !ve
      itot=i5+nou3
      if(icycle.eq.1.and.print)write(6,99)nmem,itot
 99   format('Space usage in vemvt3: available - ',i8,'   used - ',i8)
      call vemvt3(isel,no,nu,t,t(i2),t(i3),t(i4),t(i5))
      return
      end
      SUBROUTINE VEMVT3(IS,NO,NU,TI,VOE,T3,VM,VE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T3(1),VOE(1),TI(1),VM(1),VE(1)
      call ienter(80)
      call VEMVT3a(IS,NO,NU,TI,VOE,T3,VE)
      call VEMVT3b(IS,NO,NU,TI,VOE,T3,VM)
      call iexit(80)
      return
      end
      SUBROUTINE VEMVT3A(IS,NO,NU,TI,VOE,T3,VE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION T3(NU,NU,NU,NO),VOE(NU,NO,Nu,NO),TI(NU,NU,NU),
     *VE(NU,NU,NU,NO)
      IS1=1-IS
      CALL RO2PPH(1,NO,NU,TI,VOE)
      call insitu(nu,nu,no,no,ti,voe,23)
      CALL TRANMD(VOE,NU,NO,Nu,No,13)
      call rdvem4a(10,no,nu,ti,ve)
      call tranmd(ve,nu,nu,nu,no,312)
      do 21 j=1,no
      do 21 k=1,no
c
      do 501 m=1,no
      if(j.eq.k.and.j.eq.m)then
      call zeroma(t3(1,1,1,m),1,nu3)
      else
      call rdvt3onw(m,k,j,nu,t3(1,1,1,m))
      endif
501   continue
c
      call tranmd(t3,nu,nu,nu,no,23)
      do 22 i=1,no
      call symt3 (t3(1,1,1,i),nu,6)
 22   continue
      CALL MATMUL(T3,VOE(1,1,1,k),VE(1,1,1,j),NU2,NU,NOU,0,IS1)
21    continue
      CALL TRANMD(VE,NU,NU,NU,NO,231)
      call wrvemq2(10,no,nu,ti,ve)
      RETURN
      END
      SUBROUTINE VEMVT3B(IS,NO,NU,TI,VOE,T3,VM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,E,F
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION T3(NU,NU,NU,NO),VOE(1),TI(NU,NU,NU),VM(NU,NO,NO,NO)
      IS1=1-IS
      CALL RO2PPH(1,NO,NU,TI,VOE)
      call rdvem1a(11,nu,no,ti,vm)
      CALL TRANMD(vm,NU,NO,NO,NO,34)
      CALL TRANMD(vm,NU,NO,NO,NO,24)
      do 31 j=1,no
      do 31 k=1,no
c
      do 501 m=1,no
      if(j.eq.k.and.j.eq.m)then
      call zeroma(t3(1,1,1,m),1,nu3)
      else
      call rdvt3onw(m,j,k,nu,t3(1,1,1,m))
      endif
501   continue
c 
      call tranmd(t3,nu,nu,nu,no,231)
      do 32 i=1,no
      call symt3 (t3(1,1,1,i),nu,1)
32    continue
      CALL MATMUL(T3,VOE,VM(1,1,j,k),NU,No,NOU2,0,IS)
 31   continue
      CALL TRANMD(vm,NU,NO,NO,NO,24)
      CALL TRANMD(vm,NU,NO,NO,NO,34)
      call wrvem1q2(11,nu,no,ti,vm)
      RETURN
      END
      subroutine movo3(no,nu,t3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension t3(nu,nu,nu)
      do 10 i=1,no
      do 10 j=1,i
      j1=j
      if(i.eq.j)j1=j-1
      do 10 k=1,j1
      kkk=it3(i,j,k)
      call rdvt3n(kkk,nu,t3)
      call wrvt3o(kkk,nu,t3)
 10   continue
      return
      end
      SUBROUTINE SYMT311(V,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU)
      DATA TWO/2.0D+0/
      if (id.eq.23)goto 23
      if (id.eq.12)goto 12
      if (id.eq.13)goto 13
 23   CONTINUE
      DO 100 A=1,NU
      DO 100 B=1,nu
      DO 100 C=1,B
      x=v(a,b,c)+v(a,c,b)
      V(A,B,C)=x
      V(A,C,B)=x
  100 CONTINUE
      GO TO 1000
 12   CONTINUE
      DO 101 A=1,NU
      DO 101 B=1,a
      DO 101 C=1,nu
      x=v(a,b,c)+v(b,a,c)
      V(A,B,C)=x
      V(b,a,c)=x
 101  CONTINUE
      GO TO 1000
 13   CONTINUE
      DO 102 A=1,NU
      DO 102 B=1,nu
      DO 102 C=1,a
      x=v(a,b,c)+v(c,b,a)
      V(A,B,C)=x
      V(c,b,a)=x
 102  CONTINUE
      GO TO 1000
 1000 CONTINUE
      END
      SUBROUTINE GETL2(NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER A,B
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /NEWT4/NT4,NTO4,LT4,NALL4,LL4
      COMMON /NEWCCSD/NTT2
      DIMENSION TI(NU,NU,NO),T2(NO,NU,NU,NO)
      IOFF=1
      DO 1 I=1,NO
         CALL GETLST(TI,IOFF,NO,1,1,146)
         DO 10 J=1,NO
            DO 10 A=1,NU
               DO 10 B=1,NU
                  T2(I,A,B,J)=TI(B,A,J)
 10            continue
               IOFF=IOFF+NO
 1          CONTINUE
            RETURN
      END
