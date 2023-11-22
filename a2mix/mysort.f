C This subroutine sorts the four indices of the two electron 
C integrals and calculates the symmetry index, because the two
C electron integrals are stored in a list not a 4 index array.
      Subroutine mysort(mu,nu,lam,sig,indtot)
      integer mu,lam,nu,sig,tmp1,tmp2,temp1,temp2,temp3,temp4
      indx(i,j)=j+(i*(i-1))/2
      temp1=mu
      temp2=nu
      temp3=lam
      temp4=sig
      If (nu.gt.mu) Then
         tmp1=nu
         nu=mu
         mu=tmp1
       End If
       If (sig.gt.lam) Then
           tmp1=sig
           sig=lam
           lam=tmp1
       End If
       If (lam.gt.mu) Then
           tmp1=lam
           tmp2=sig
           lam=mu
           sig=nu
           mu=tmp1
           nu=tmp2
        End If
        If ((lam.eq.mu).and.(sig.gt.nu)) Then
           tmp1=lam
           tmp2=sig
           lam=mu
           sig=nu
           mu=tmp1
           nu=tmp2
         End If
            indij=indx(mu,nu)
            indkl=indx(lam,sig)
            indtot=indx(indij,indkl)
      mu=temp1
      nu=temp2
      lam=temp3
      sig=temp4
      Return
      End
