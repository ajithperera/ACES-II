      Subroutine form_modf_fbar(W,Maxcor,Iuhf,Imult)

      Implicit Double Precision (A-H,O-Z)

      Dimension W(Maxcor)

      Data Onem /-1.0D0/

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

C Modify F(a,b) and F(i,j) by removing the f(a,b) and f(i,j).
C This is relevant only for the non-HF methods and required 
C by doubles correction. 

      Do Ispin = 1, 1+Iuhf 
         Len_pq = Irpdpd(1,18+Ispin)

         I000 = 1 
         I010 = I000 + Len_pq
         Iend = I010
         If (Iend .Gt. Maxcor) Call Insmem("form_modf_fbar",Iend,Maxcor)

         If (Imult .EQ. 1) Then
           Call Getlst(W(I000),1,1,1,Ispin,92) 
         Else
           Call Getlst(W(I000),1,1,1,10,92) 
         Endif 

         Call Getlst(W(I010),1,1,1,Ispin+2,92) 
         Call Daxpy(Len_aa,ONEM,W(I010),1,W(I000),1) 

         If (Imult .EQ. 1) Then
            Call putlst(W(I000),1,1,1,10+Ispin,92)
         Else
            Call putlst(W(I000),1,1,1,13,92)
         Endif

         Len_pq = Irpdpd(1,20+Ispin) 

         I000 = 1
         I010 = I000 + Len_pq
         Iend = I010
         If (Iend .Gt. Maxcor) Call Insmem("form_modf_fbar",Iend,Maxcor)

         If (Imult .EQ. 1) Then
           Call Getlst(W(I000),1,1,1,Ispin,91) 
         Else
           Call Getlst(W(I000),1,1,1,10,91)   
         Endif

         Call Getlst(W(I010),1,1,1,Ispin+2,91)
         Call Daxpy(Len_aa,ONEM,W(I010),1,W(I000),1) 

         If (Imult .EQ. 1) Then
            Call putlst(W(I000),1,1,1,10+Ispin,91)
         Else
            Call putlst(W(I000),1,1,1,13,91)
         Endif

      Enddo 

      Return
      End 
