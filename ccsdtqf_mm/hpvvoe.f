      subroutine hpvvoe(ipt,no,nu,v,voe,t2)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/vrinpi/ifirst(512)
      dimension v(1),voe(nu,no,no,nu),t2(nu,nu,no,no)
      do 302 i=1,nu
         irc=ifirst(i)
         call rdvv(46,nu3,irc,v)
         if (ipt.eq.1)goto 10
         call trant3(v,nu,3)
         if (ipt.eq.2)goto 10
         call symt21(v,nu,nu,nu,1,13)
 10      continue
         call matmul(v,voe(1,1,1,i),t2,nu2,no2,nu,0,1)
 302  continue
      return
      end
