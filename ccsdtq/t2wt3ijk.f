      subroutine drt2wt3ijk(i,j,m,no,nu,t2,t3,vm,tmp,vepak,ive)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION T2(1),T3(1),VM(1),tmp(1),vepak(1),ive(1)
      call veccop(nu3,tmp,t3)
      call t2wt3ijknew(i,j,m,no,nu,t2,tmp,vm,vepak,ive)
      call veccop(nu3,tmp,t3)
      call trant3(tmp,nu,3)
      call t2wt3ijknew(i,m,j,no,nu,t2,tmp,vm,vepak,ive)
      call veccop(nu3,tmp,t3)
      call trant3(tmp,nu,1)
      call t2wt3ijknew(j,i,m,no,nu,t2,tmp,vm,vepak,ive)
      call veccop(nu3,tmp,t3)
      call trant3(tmp,nu,5)
      call t2wt3ijknew(j,m,i,no,nu,t2,tmp,vm,vepak,ive)
      call veccop(nu3,tmp,t3)
      call trant3(tmp,nu,4)
      call t2wt3ijknew(m,i,j,no,nu,t2,tmp,vm,vepak,ive)
      call veccop(nu3,tmp,t3)
      call trant3(tmp,nu,2)
      call t2wt3ijknew(m,j,i,no,nu,t2,tmp,vm,vepak,ive)
      return
      end
      subroutine drt2wt3ij(i,j,no,nu,t2,t3,vm,tmp,vepak,ive)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION T2(1),T3(1),VE(1),VM(1),tmp(1),vepak(1),ive(1)
      call veccop(nu3,tmp,t3)
      call trant3(tmp,nu,2)
      call t2wt3ijknew(i,j,j,no,nu,t2,tmp,vm,vepak,ive)
      call veccop(nu3,tmp,t3)
      call trant3(tmp,nu,5)
      call t2wt3ijknew(j,i,j,no,nu,t2,tmp,vm,vepak,ive)
      call veccop(nu3,tmp,t3)
      call t2wt3ijknew(j,j,i,no,nu,t2,tmp,vm,vepak,ive)
      return
      end
      subroutine drt2wt3jm(j,m,no,nu,t2,t3,vm,tmp,vepak,ive)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION T2(1),T3(1),VE(1),VM(1),tmp(1),vepak(1),ive(1)
      call veccop(nu3,tmp,t3)
      call trant3(tmp,nu,3)
      call t2wt3ijknew(j,j,m,no,nu,t2,tmp,vm,vepak,ive)
      call veccop(nu3,tmp,t3)
      call t2wt3ijknew(j,m,j,no,nu,t2,tmp,vm,vepak,ive)
      call veccop(nu3,tmp,t3)
      call trant3(tmp,nu,5)
      call t2wt3ijknew(m,j,j,no,nu,t2,tmp,vm,vepak,ive)
      return
      end
      SUBROUTINE T2WT3ijknew(i,j,m,NO,NU,T2,T3,VM,vepak,ive)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical sdt
      common/itrat/iter,ittt,itt1
      common /newopt/nopt(6)
      common/wpak/nfr(30),nsz(30)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION T2(NU,NU,NO,NO),T3(1),VM(Nu,No,NO,NO),ive(1),
     *vepak(nu,nu,nu,no)
      call ienter(130)
      sdt=nopt(1).gt.1
      nnve=nsz(6)
      ioff=0
      if(sdt)then
         nnve=nsz(8)
         ioff=2*no+2
      endif
      CALL SYMT3(T3,NU,3)
      call trant3(t3,nu,4)
      call matmul(vepak(1,1,1,m),t3,t2(1,1,i,j),nu,nu,nu2,0,0)
      call trant3(t3,nu,1)
      CALL MATMUL(t3,VM(1,1,M,i),T2(1,1,1,J),NU2,no,NU,0,1)
      call iexit(130)
      RETURN
      END

