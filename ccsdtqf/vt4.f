      subroutine vt4(no,nu,ti,fpp,v,t4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension fpp(nu,nu),v(nu,nu,nu),t4(nu,nu),ti(1)
      do 301 i=1,nu
         call rpakv(35,nu,i,ti,ti(nu2+1),v)
         call trant3(v,nu,1)
         call symt21(v,nu,nu,nu,1,13)
         call matmul(fpp,v,t4(1,i),1,nu,nu2,1,1)
 301  continue
      return
      end
