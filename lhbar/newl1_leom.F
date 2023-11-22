#include "flags.h"
      SUBROUTINE NEWL1_LEOM(Work,Length,Iuhf)

      Implicit Double Precision (A-H, O-Z)
      Dimension Work(Length)
      Logical Sing
      
#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"

       Sing   = .False.
       Sing   = (Iflags(h_IFLAGS_calc) .gt. 9)

       I000 = 1
       I010 = I000 + NT(1)
       I020 = I010 + NT(1)
       I030 = I020
       I040 = I020
       If (Iuhf .Gt. 0) Then
           I030 = I020 + NT(2)
           I040 = I030 + NT(2)
       Endif 

       If (I040 .ge. Length) Call Insmem("@-NEWL1_LEOM",I040,Length)

       If (Sing) Then
      
          Call Getlst(Work(I000), 1, 1, 1, 1, 93)
          Call Getlst(Work(I010), 1, 1, 1, 3, 90)
          Call Daxpy(NT(1), 1.0D0, Work(I000), 1, Work(I010), 1)
          Call Newl1(Work(I010), Work, Length, 1)
          Call Putlst(Work(I010), 1, 1, 1, 3, 90)

          If (Iuhf .ne. 0) Then

             Call Getlst(Work(I020), 1, 1, 1, 2, 93)
             Call Getlst(Work(I030), 1, 1, 1, 4, 90)
             Call Daxpy(NT(1), 1.0D0, Work(I020), 1, Work(I030), 1)
             Call Newl1(Work(I030), Work, Length, 1)
             Call Putlst(Work(I030), 1, 1, 1, 4, 90)

          Endif 
      Endif 
C
      Return
      End
