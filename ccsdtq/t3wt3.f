      SUBROUTINE DRT3WT3(NO,NU,T,EH,EP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/itrat/icycle,mx,icn
      DIMENSION T(1),EH(NO),EP(NU)
      call ienter(59)
      i1=1             !fhh,
      i2=i1+no2        !t3
      i3=i2+nou3       !ve
      i4=i3+nou3       !voe
      i5=i4+no2u2      !vo
      i6=i5+no2u2      !vtm
      itot=i6+nu3
      if(icycle.eq.1)write(6,*)'t3wt3hp-space:',itot
      call flush(6)
      call t3wt3hp(no,nu,t,t(i2),t(i3),t(i4),t(i5),t(i6))
      i1=1             !t3
      i2=1+nou3        !pz
      i3=i2+no4        !vtm
      itot=i3+nu3
      if(icycle.eq.1)write(6,*)'t3wt3hh-space:',itot
      call flush(6)
      call t3wt3hh(no,nu,t,t(i2),t(i3))
      i1=1       !fpp
      i2=1+nu2   !t3
      i3=i2+nu3  !vtm
      i4=i3+nu3  !v
      itot=i4+nu4
      if(icycle.eq.1)write(6,*)'t3wt3pp-space:',itot
      call flush(6)
      call t3wt3pp(no,nu,t,t(i2),t(i3),t(i4))
      call iexit(59)
      RETURN
      END
      SUBROUTINE T3WT3HH(NO,NU,T3,PZ,VTM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T3(NU,NU,NU,NO),PZ(NO,NO,NO,NO),vtm(*)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DATA ZERO/0.0D+0/,HALF/0.5D+00/,TWO/2.0D+00/,THREE/3.0D+00/
      call ienter(60)
      CALL RDOV4(11,NU,NO,VTM,PZ)
      DO 400 i=1,no
         im1=i-1
         ip1=i+1
      DO 300 M=1,NO
      do 200 n=1,no
         call rdvt3onw(n,m,i,nu,T3(1,1,1,n))
 200  continue
      DO 401 j=i+1,NO
      DO 401 k=1,j
         if(k.eq.i)goto 401
         if(k.gt.i)then
            kkk=it3(j,k,i)
            call rdvt3n(KKK,NU,VTM)
         else
            kkk=it3(j,i,k)
            CALL RDVT3N(KKK,NU,VTM)
            call trant3(vtm,nu,1)
         endif
         CALL  MATMUL(T3,PZ(1,M,J,k),VTM,NU3,1,NO,0,0)
         if(k.lt.i)call trant3(vtm,nu,1)
         CALL WRVT3N(KKK,NU,VTM)
 401  CONTINUE
      DO 1403 k=ip1,NO
         kkk=it3(k,i,i)
         CALL RDVT3N(KKK,NU,VTM)
         call vecmul(vtm,nu3,half)
         call trant3(vtm,nu,1)
         CALL  MATMUL(T3,PZ(1,M,k,i),VTM,NU3,1,NO,0,0)
         call symt311(vtm,nu,23)
         call trant3(vtm,nu,1)
         CALL WRVT3N(KKK,NU,VTM)
 1403 continue
      DO 1402 K=1,im1
         kkk=it3(i,i,k)
         CALL RDVT3N(KKK,NU,VTM)
         call vecmul(vtm,nu3,half)
         call trant3(vtm,nu,1)
         CALL  MATMUL(T3,PZ(1,M,i,K),VTM,NU3,1,NO,0,0)
         call trant3(vtm,nu,1)
         call symt311(vtm,nu,12)
         CALL WRVT3N(KKK,NU,VTM)
 1402 continue
      DO 2401 j=1,im1
         DO 2402 k=j,im1
            kkk=it3(i,k,j)
            CALL RDVT3N(KKK,NU,VTM)
            call trant3(vtm,nu,3)
            CALL  MATMUL(T3,PZ(1,M,j,k),VTM,NU3,1,NO,0,0)
            call trant3(vtm,nu,3)
            CALL WRVT3N(KKK,NU,VTM)
 2402    continue
 2401 continue
 300  continue
 400  CONTINUE
      call iexit(60)
      RETURN
      END
      SUBROUTINE T3WT3HP(NO,NU,FHH,T3,VE,VOE,VO,VTM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION T3(NU,NU,NU,NO),VOE(NU,NO,NU,NO),VO(NU,NO,NU,NO),
     *VTM(*),FHH(NO,NO),ve(nu,nu,nu,no)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DATA ZERO/0.0D+0/,HALF/0.5D+00/,TWO/2.0D+00/,THREE/3.0D+00/
      call ienter(61)
      READ(NALL4,REC=5*no+2*nu+5)FHH
      CALL RO2PHP(11,NO,NU,VTM,VOE)
      CALL RO2PHP(12,NO,NU,VTM,VO)
      CALL TRANSQ(VOE,NOU)
      CALL TRANSQ(VO,NOU)
      DO 402 I=1,NO
      DO 402 J=1,i
      do 401 m=1,no
      call rdvt3onw(m,i,j,nu,vtm)
      CALL VECCOP(NU3,T3(1,1,1,M),VTM)
401   continue
      call veccop(nou3,ve,t3)
      call tranmd(ve,nu,nu,nu,no,312)
      call symt3i(ve,no,nu)
      DO 403 K=1,no
         if(i.eq.j.and.i.eq.k)goto 403
      call rdvt3nnw(i,j,k,nu,vtm)
      if(k.eq.j.or.k.eq.i)call vecmul(vtm,nu3,half)
      call trant3(vtm,nu,4)
      CALL MATMUL(T3,FHH(1,K),VTM,NU3,1,NO,0,1)
      if(i.eq.j)call vecmul(vtm,nu3,half)
      CALL  MATMUL(T3,VO(1,1,1,K),VTM,NU2,NU,NOU,0,1)
      if(i.eq.j)call symt311(vtm,nu,23)
      CALL TRANMD(T3,NU,NU,NU,NO,23)      
      call trant3(vtm,nu,1)
      if(i.ne.j)CALL  MATMUL(T3,VO(1,1,1,K),VTM,NU2,NU,NOU,0,1)
      call trant3(vtm,nu,5)
      CALL TRANMD(T3,NU,NU,NU,NO,312)      
      CALL MATMUL(T3,VO(1,1,1,K),VTM,NU2,NU,NOU,0,1)
      CALL TRANMD(T3,NU,NU,NU,NO,13)      
      call trant3(vtm,nu,2)
      CALL  MATMUL(ve,VOE(1,1,1,K),VTM,NU2,NU,NOU,0,0)
      if(k.gt.j)then
         if(k.lt.i)then
            kkk=it3(i,k,j)
            call trant3(vtm,nu,1)
         else
            kkk=it3(k,i,j)
            call trant3(vtm,nu,4)
         endif
      else
         kkk=it3(i,j,k)
      endif
      if(k.eq.j)call symt311(vtm,nu,23)
      if(k.eq.i)call symt311(vtm,nu,12)
      CALL WRVT3N(KKK,NU,VTM)
 403  continue
 402  CONTINUE
      call iexit(61)
      RETURN
      END
      SUBROUTINE T3WT3PP(NO,NU,FPP,T3,VTM,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,E,F
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      COMMON/NEWOPT/NOPT(6)
      DIMENSION fpp(nu,nu),T3(*),vtm(*),V(*)
      COMMON/RESLTS/CMP(30)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DATA ZERO/0.0D+0/,HALF/0.5D+00/,TWO/2.0D+00/,THREE/3.0D+00/,
     *FOUR/4.0D+00/,SIX/6.0D+0/,EIGHT/8.0D+0/
      call ienter(62)
      NLAST=5*NO+2*NU
      IASV=NLAST+6
      READ(NALL4,REC=IASV)FPP
      call transq(fpp,nu)
      i10=10
      CALL RDOV4(I10,NO,NU,VTM,V)
      IF(I10.EQ.0)CALL MTRANS(V,NU,7)
      DO 402 I=1,NO
      DO 402 J=1,i
      DO 402 K=1,j
      IF (I.EQ.K) GO TO 402
      kkk=it3(i,j,k)
      CALL RDVT3N(KKK,NU,VTM)
      CALL RDVT3O(KKK,NU,T3)
      if(j.eq.k)call vecmul(vtm,nu3,half)
      call matmul(v,t3,vtm,nu2,nu,nu2,0,0)
      call trant3(t3,nu,4)
      call trant3(vtm,nu,4)
      call matmul(fpp,t3,vtm,nu,nu2,nu,0,0)
      if(i.eq.j)call vecmul(vtm,nu3,half)
      call trant3(t3,nu,3)
      call trant3(vtm,nu,3)
      if(j.eq.k)call symt311(vtm,nu,13)
      if(j.ne.k)call matmul(fpp,t3,vtm,nu,nu2,nu,0,0)
      call trant3(t3,nu,5)
      call trant3(vtm,nu,5)

      if(j.ne.k)call matmul(v,t3,vtm,nu2,nu,nu2,0,0)
      if(i.eq.j)call symt311(vtm,nu,13)
      if(i.ne.j)call matmul(fpp,t3,vtm,nu,nu2,nu,0,0)
      call trant3(t3,nu,5)
      call trant3(vtm,nu,5)
      if(i.ne.j)call matmul(v,t3,vtm,nu2,nu,nu2,0,0)
      call trant3(vtm,nu,3)
      CALL WRVT3N(KKK,NU,VTM)
 402  CONTINUE
      call iexit(62)
      RETURN
      END
