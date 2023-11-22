      Subroutine Read_Lanczos_vars()

      Implicit Double Precision (A-H,O-Z)

      Logical Zmat_Present
      Character*80 Wrk

#include "lanczos_vars.com"

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " Entered --- Read_Lanczos_vars ---"
      Write(6,*)
#endif 
      Inquire(File="ZMAT",Exist=Zmat_present)
      
      If (Zmat_present) Then
         Open(Unit=4,File="ZMAT",Form="FORMATTED",Status="OLD")
         Rewind(4)
      Else 
        Write(6,"(1x,a,a)") "@-read_Lanczos_vars; ZMAT file not",
     +                     " present"
        Call Errex
      Endif 

 300  Read(4,"(a)",End=900) Wrk
      If (Wrk(1:13) .NE. "*LANCZOS_VARS") Goto 300 

C The maximum allowed number of look-ahead Lanczos steps

      Read(4,"(i4)",End=900) Nlim 

C Lead dimension of the tridiagonal matrix

      Read(4,"(i4)",End=900) Maxn

C Maximum length of the VW step

      Read(4,"(i4)",End=900) Maxvw

C Lead dimension of arrays iwk and dwk

      M = 2*Maxvw + 2

      Read(4,"(F4.1)",End=900) Norm

C The convergence monitoring (the two options are strict and light).

      Read(4,"(a)",End=900) Wrk
      Convrg = Wrk(1:6)

      Go to 999

 900  Write(6,901) 
 901  Format(T3,"@-read_lanczos_vars,LANCZOS_VARS namelist not found',
     + or incomplete.")
      Write(6,"(1x,a)") " The default values are being used"
      
C Set the default values.
#ifdef _DULAL
       Nlim  = 10
       Maxn  = 1000
       MaxvW = 10 
       M     = 22
       Norm  = 1.0D0
#else
       Maxn   = 10000
       Maxvw  = 1
       Convrg = "Light"
#endif 

 999  Continue 

#ifdef _DULAL 
      If (Maxn .Gt. 1000)  Then
         Write(6,"(a,a)") "The dimension of the tridiagonal-matrix",
     +                    " exceded the current limt of 1002." 
         Call Errex 
      Endif 

      If (Nlim .Gt. Maxn_int)  Then
         Write(6,"(a,a)") "The maximum allowed number of look-ahead",
     +                    " Lanczos steps exceed the current limit"
         Write(6,"(a)")   " of 10." 
         Call Errex 
      Endif 

      If (Maxvw .Gt. 10)  Then
         Write(6,"(a,a)") "The maximum allowed block size of the", 
     +                    " VW step exceed the current limit"
         Write(6,"(a)")   " of 10."
         Call Errex 
      Endif 

      If (M .Gt. 22)  Then
         Write(6,"(a,a)") "The lead dimension of arrays iwk and dwk", 
     +                    " exceed the current limit of 10." 
         Call Errex 
      Endif 
#else 
      If (Maxn .Gt. 10000)  Then
         Write(6,"(a,a)") "The dimension of the tridiagonal-matrix",
     +                    " exceded the current limt of 1002." 
         Call Errex 
      Endif 

      If (Maxn .Le. 0)  Then
         Write(6,"(a,a)") "The dimension of the tridiagonal-matrix",
     +                    " must be greater than zero." 
         Call Errex 
      Endif 
    
      Itest = 100
      If (Convrg(1:6) .EQ. "Light " .OR. Convrg(1:6) .EQ. "Strict") 
     +   Itest=0

      If (Itest .NE. 0) Then
         Write(6,"(a,a)") "Invalid option for the convergence",
     +                    " criteria. Only Strict or Light is allowed." 
         Call Errex 
         
      Endif 

#endif 
C Initialize the the full tridigonal matirx to zero.

      Call Dzero(Trd,Maxn*Maxn)

      Close(Unit=4,Status="KEEP")
#ifdef _DEBUG_LVL0
      Write(6,*)
#ifdef _DULAL
      Write(6,"(a)") " The settings for the Lanczos solve dulal"
      Write(6,"(a,a)") " Nlim,Maxn,Maxvw,M,Norm"
      Write(6,"(4(1x,i4),1x,F5.3)") Nlim,Maxn,Maxvw,M,norm 
#else
      Write(6,"(a)") " The settings for the Lanczos_light_v1"
      Write(6,"(a,a)") " Maxn and Maxvw"
      Write(6,"(2(1x,i4),1x,F5.3)") Maxn,Maxvw
      Write(6,"(1x,a)") Convrg 
#endif 
      Write(6,*)
      Write(6,*)
      Write(6,"(a)") "---Exit read_lanczos_vars ---"
      Write(6,*)
#endif
      Return
      End
