      Subroutine Form_wpqrs_abab(AO_2ints,W_pxxx,W_pqxx,W_pqrx,
     +                          W_pqrs_abab,Cmo_a,Cmo_b,Nbasis,
     +                          Naobfns)

      implicit double precision (a-h,o-z)
      integer m, i, j, k,l,p,q,r,s

      Dimension AO_2ints(Naobfns,Naobfns,Naobfns,Naobfns)
      Dimension W_pqrs_abab(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pxxx(Nbasis,Naobfns,Naobfns,Naobfns)
      Dimension W_pqxx(Nbasis,Nbasis,Naobfns,Naobfns)
      Dimension W_pqrx(Nbasis,Nbasis,Nbasis,Naobfns)

      Dimension Cmo_a(Naobfns,Nbasis),Cmo_b(Naobfns,Nbasis)
    
      Do p=1,Nbasis
      Do j=1,Naobfns
      Do k=1,Naobfns
      Do l=1,Naobfns
      Sum =0.0d0
      
      do i=1,Naobfns

      Sum = Sum + (AO_2ints(i,j,k,l)*Cmo_a(i,p))
      
      End do
      W_pxxx(p,j,k,l)=Sum

      End do                           
      End do   
      End do
      
      End do

      Do p=1,Nbasis
      Do q=1,Nbasis
      Do k=1,Naobfns
      Do l=1,Naobfns
      Sum =0.0d0
      
      Do j=1,Naobfns
      Sum =sum + (W_pxxx(p,j,k,l)*Cmo_b(j,q))

      End do
      W_pqxx(p,q,k,l)=Sum  
      End do
      End do
      End do 
      End do
      
      Do p=1,Nbasis 
      Do q=1,Nbasis
      Do r=1,Nbasis
      Do l=1,Naobfns
      Sum=0.0d0
      do k=1,Naobfns

      Sum = sum +(W_pqxx(p,q,k,l)*Cmo_a(k,r))

      End do
      W_pqrx(p,q,r,l) = Sum
      End do
      End do
      End do
      End do
      
      Do p=1,Nbasis 
      Do q=1,Nbasis 
      Do r=1,Nbasis 
      Do s=1,Nbasis 
      Sum=0.0d0
      Do l=1,Naobfns

      Sum = sum +(W_pqrx(p,q,r,l)*Cmo_b(l,s))

      End do
      W_pqrs_abab(p,q,r,s) = Sum
      End do
      End do
      End do
      End do

      Return
      End 
