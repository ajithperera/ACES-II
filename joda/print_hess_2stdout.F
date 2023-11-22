      Subroutine Print_hess_2stdout(Hess,Ndim,Nreals)

      Implicit Double Precision(A-H,O-Z)
       
      Dimension Hess(Ndim,Ndim)
      Character*6 Names(Ndim)
      Character*6 String
      Character*1 V1
      Character*2 V2
      Character*2 V3

#include "mxatms.par"
#include "cbchar.com"
      Parameter (Kcol=6)

      Data Dnull /0.0D0/

C Here Ndim = 3*Nreals
C Construct the atom names (to appear like printing from vdint)

      Do I = 1, Nreals
         Ioff = (I-1)*3 + 1
         If (I .Lt. 10) Then
             Write(V1,"(I1)") I 
             V2 = "0"//V1
             Do J = 1, 3
                String = Zsym(I)
                If (J .Eq. 1) Names(Ioff)   = String(1:2)//V2//" x"
                If (J .Eq. 2) Names(Ioff+1) = String(1:2)//V2//" y"
                If (J .Eq. 3) Names(Ioff+2) = String(1:2)//V2//" z"
             Enddo
         Else if (I .Lt. 99) Then
             Write(V3,"(I2)") I 
             Do J = 1, 3
                String = Zsym(I)
                If (J .Eq. 1) Names(Ioff)   = String(1:2)//V3//" x"
                If (J .Eq. 2) Names(Ioff+1) = String(1:2)//V3//" y"
                If (J .Eq. 3) Names(Ioff+2) = String(1:2)//V3//" z"
             Enddo
         Else
             Write(6,"(a)") "@-print_hess-2stdout, More than 100",
     +                      " atoms is not allowed"
             Call Errex
         Endif 
      Enddo
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6, "(3(1x,A6))") (Names(J),J=1,Ndim)
#endif 
      Write(6,"(29x,a)") "Molecular hessian"
      Write(6,"(29x,a)") "-----------------"
      Write(6,"()")

      Ibegin = 1 
      Nrow   = 3*Nreals 
      Ilast  = Min(Nrow,Kcol)
 51   Continue 
      Write(6,1000) (Names(I),I=Ibegin,Ilast)
      Write(6,"()")
      Ncol = 1
      Do K = Ibegin, Nrow 
         Do I = 1, Ncol
C For vibrational averaging it is preferable to have zero rows printed.
C Otherwise processing the Hessian for displaced geos. become cumbersome.
C 
CSSS            If (Hess(K,(Ibegin-1)+I) .Ne. Dnull) Go to 401
            Go to 401
         Enddo
         Go to 301
  401    Continue 
         Write(6,2000) " ", Names(K), (Hess(K,(Ibegin-1)+J),J=1,Ncol)
         If (Mod(K,3) .Eq. 0) Write(6,"()")
  301    If (K .Lt. (Ibegin+Kcol-1)) Ncol = Ncol + 1
      Enddo

      Write(6,"()")
      Ilast  =  Min(Ilast+Kcol,Nrow)
      Ibegin = Ibegin + Ncol
      If (Ibegin .Le. Nrow) Goto 51

 1000 Format (8X,6(3X,A6,3X),(3X,A6,3X))
 2000 format (A1,A6,6F14.10)

      Return
      End

 
   
