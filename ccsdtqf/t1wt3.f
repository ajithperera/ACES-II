      SUBROUTINE T1WT3(NO,NU,TI,T1,VOE,t3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D,E,F
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION TI(NU,NU,NU),VOE(Nu,NU,No,NO),T1(Nu,No),t3(nu,nu,nu)
      DATA ZERO/0.0D+0/,TWO/2.0D+00/,HALF/0.5D+00/,FOUR/4.0D+0/
      call ienter(52)
      call zeroma(t1,1,nou)
      call insitu(no,nu,nu,no,ti,voe,13)
      call tranmd(voe,nu,nu,no,no,12)
      DO 10 I=1,NO
         i1=i-1
      DO 10 j=1,i1
         j1=j-1
      DO 10 k=1,j1
      call drt1wt3ijk(i,j,k,NO,NU,T1,VOE,ti,t3)
 10   CONTINUE
      DO 11 I=1,NO
         i1=i-1
      DO 11 j=1,i1
      call drt1wt3ij(i,j,NO,NU,T1,VOE,ti,t3)
 11   continue
      DO 12 j=1,NO
         j1=j-1
         DO 12 k=1,j1
            call drt1wt3jk(j,k,NO,NU,T1,VOE,ti,t3)
 12      continue
      call tranmd(voe,nu,nu,no,no,12)
      call insitu(nu,nu,no,no,ti,voe,13)
      call trt1(nu,no,ti,t1)
      call iexit(52)
      RETURN
      END
      SUBROUTINE drt1wt3ijk(i,j,k,NO,NU,T1,VOE,ti,t3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION t1(1),VOE(1),t3(1),ti(1)
      CALL RDVT3ONW(i,j,k,NU,TI)
      call T1WT3ijk(i,j,k,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,1)
      call T1WT3ijk(i,k,j,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,4)
      call T1WT3ijk(j,i,k,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,1)
      call T1WT3ijk(j,k,i,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,5)
      call T1WT3ijk(k,i,j,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,1)
      call T1WT3ijk(k,j,i,NO,NU,T1,VOE,ti,t3)
      return
      end
      SUBROUTINE drt1wt3ij(i,j,NO,NU,T1,VOE,ti,t3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION t1(1),VOE(1),t3(1),ti(1)
      CALL RDVT3ONW(i,j,j,NU,TI)
      call T1WT3ijk(i,j,j,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,2)
      call T1WT3ijk(j,i,j,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,1)
      call T1WT3ijk(j,j,i,NO,NU,T1,VOE,ti,t3)
      return
      end
      SUBROUTINE drt1wt3jk(j,k,NO,NU,T1,VOE,ti,t3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION t1(1),VOE(1),t3(1),ti(1)
      CALL RDVT3ONW(j,j,k,NU,TI)
      call T1WT3ijk(j,j,k,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,1)
      call T1WT3ijk(j,k,j,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,2)
      call T1WT3ijk(k,j,j,NO,NU,T1,VOE,ti,t3)
      return
      end
      SUBROUTINE T1WT3ijk(i,j,k,NO,NU,T1,VOE,ti,t3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION t1(nu,no),VOE(Nu,NU,No,NO),t3(1),ti(1)
      call ienter(129)
      call smt3four(NU,T3,ti)
      call matmulsk(t3,voe(1,1,j,k),t1(1,i),nu,1,nu2,0,0)
      call iexit(129)
      return
      end
      SUBROUTINE smt3four(NU,T3,V3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D,E,F
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension t3(nu,nu,nu),v3(nu,nu,nu)
      data two/2.0d+0/
      do 1 a=1,nu
      do 1 b=1,nu
      do 1 c=1,nu
         t3(a,b,c)=(v3(a,b,c)-v3(b,a,c))*TWO-v3(a,c,b)+v3(b,c,a)
 1    continue
      return
      end
