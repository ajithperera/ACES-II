      Subroutine Make_Res_0(Lam, Mu, xph, Yph, Zph, Res)
      
      Implicit Double Precision (A-H, O-Z)
C
#include "ecp.par"
#include "ecp.com"
     
      JD = Lam*(Lam+1) - Mu + 1
     
      Imn = Lmf(Jd)
      Imx = Lmf(Jd+1) - 1
C
      Res = 0.0D0

#ifdef _DEBUG_LVLM1
      Write(6,"(a,10(1x,I3))") "At entry to make_res, Imn and Jmn ",
     &                          imn, Imx
#endif
      Do J = Jmn, Jmx
     
         Jx = Lmx(J)
         Jy = Lmy(J)
         Jz = Lmz(J)  

         If (Jx .EQ. 0 .AND. Jy .EQ. 0 .AND. Jz .EQ. 0) Then
             Res = Res + Zlm(J)
         Else if (Jx .EQ. 0) Then
             Res = Res + Zlm(J)*(Yph**Jy)*(Zph**Jz)
         Else if (Jy .EQ. 0) Then
             Res = Res + Zlm(J)*(Xph**Jx)*(Zph**Jz)
         Else if (Jz .EQ. 0) Then
             Res = Res + Zlm(J)*(Xph**Jx)*(Yph**Jy)
         Else
             Res = Res + Zlm(J)*(Xph**Jx)*(Yph**Jy)*(Zph**jz)
         Endif
C           
#ifdef _DEBUG_LVLM1
      Write(6,"(a,3(1x,i2),2(1x,F15.8))"), "The real Sph ",jx,jy,jz,
     &          Zlm(J),Res
#endif
      Enddo
     
      Return
      End
            
           
