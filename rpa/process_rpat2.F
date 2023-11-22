      Subroutine Process_rpat2(Work,Length,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Length)
      Character*4 Spin(2)
      Dimension Ioff_row(9),Ioff_col(9)

#include "machsp.com"
#include "jodaflags.com"
#include "syminf.com"
#include "sympop.com"
#include "sym.com"
#include "info.com"

      Data Ione,One /1,1.0D0/
      Data Spin /"AAAA","BBBB"/

C The T2 vectors generated in the RPA code are ordered in ai,bj form 
C regardless of the spin. We need to reorder them into ab,ij form to
C interface with the rest of the code. The RPA vs rCC identity holds
C only when both singlet and triplets are included. Therefore RHF
C runs are flagged.

      If (Iuhf .EQ. 0) Then
         Write(6,"(a)") " Both singlet and triplet roots are needed.",
     +                  " Therfore switch to UHF reference to proceed."
         Call Errex 
      Else

C First AAAA and BBBB spin cases.

         Do Ispin = 1, IUHF+1  

            Nsize_s = 0
            Nsize_t = 0
            Do Irrepx = 1, Nirrep
               Nrow_s = Irpdpd(Irrepx,8+Ispin)
               Ncol_s = Irpdpd(Irrepx,8+Ispin)
               Nrow_t = Irpdpd(Irrepx,18+Ispin)
               Ncol_t = Irpdpd(Irrepx,20+Ispin)
               Nsize_s = Nsize_s + Nrow_s*Ncol_s
               Nsize_t = Nsize_t + Nrow_t*Ncol_t
            Enddo 

            If (Nsize_s .ne. Nsize_t) Then
             Write(6,"(a,a,a)") " The dimensions of the same spin",
     +                          " incoming and outgoing arrays",
     +                          " do not match."
               Call Errex
            Endif 
            It2_s = Ione
            It2_t = It2_s + Nsize_s
            Iend  = It2_t
            If (Iend .Ge. Length) Call Insmem("Process_rpat2",
     +                                         Iend,Maxcor)

C Read in T2(ai,bj) order and obtain T2(ab,ij). Note that this T2
C is not antisymmetric. 

            Ioff = It2_s
            Do Irrepx = 1, Nirrep 
               Nrow = Irpdpd(Irrepx,8+Ispin)
               Ncol = Irpdpd(Irrepx,8+Ispin)
               Call Getlst(Work(Ioff),1,Ncol,1,Irrepx,199+Ispin)
               Ioff = Ioff + Nrow*Ncol
            Enddo 

            Iscr   = It2_t + Nsize_t
            Iend   = Iscr  + Nsize_t
            If (Iend .Ge. Length) Call Insmem("Process_rpat2",
     +                                         Iend,Maxcor)
#ifdef _DEBUG_LVLM
       Write(6,*)
       Write(6,"(a,4a,a)") Spin(Ispin)," RPA T2 vectors"
       Write(6, "(5(F15.10,1x))") (Work(It2_s+i-1), i=1,Nsize_s)
#endif 
            Call Sstgen(Work(It2_s),Work(It2_t),Nsize_t,Vrt(1,Ispin),
     +                  Pop(1,Ispin),Vrt(1,Ispin),Pop(1,Ispin),
     +                  Work(Iscr),1,"1324")
#ifdef _DEBUG_LVLM
       Write(6,*)
       Write(6,"(a,4a,a)") " Reordered ",Spin(Ispin)," RPA T2 vectors"
       Write(6, "(5(F15.10,1x))") (Work(It2_t+i-1), i=1,Nsize_t)
       If (Ispin .Eq. 1) call checksum("T2-AAAA",Work(It2_t),Nsize_t)
       If (Ispin .Eq. 2) call checksum("T2-BBBB",Work(It2_t),Nsize_t)
#endif 
            Call Putall(Work(It2_t),Nsize_t,1,43+Ispin)
        Enddo 

C ABAB block

        Nsize_s = 0
        Nsize_t = 0
        Do Irrepx = 1, Nirrep
           Nrow_s = Irpdpd(Irrepx,9)
           Ncol_s = Irpdpd(Irrepx,10)
           Nrow_t = Irpdpd(Irrepx,15)
           Ncol_t = Irpdpd(Irrepx,14)
           Nsize_s = Nsize_s + Nrow_s*Ncol_s
           Nsize_t = Nsize_t + Nrow_t*Ncol_t
        Enddo 
        If (Nsize_s .ne. Nsize_t) Then
           Write(6,"(a,a,a)") " The dimensions of the opposite",
     +                        " spin incoming and outgoing arrays",
     +                        " do not match."
          Call Errex
        Endif 

        It2_s  = Ione
        It2_t  = It2_s + Nsize_s
        Iend   = It2_t
        If (Iend .Ge. Length) Call Insmem("Process_rpat2",
     +                                     Iend,Maxcor)

C Read in T2(AI,bj) order and obtain T2(Ab,Ij)

        Ioff = It2_s 
        Do Irrepx = 1, Nirrep
           Nrow = Irpdpd(Irrepx,9)
           Ncol = Irpdpd(Irrepx,10)
           Call Getlst(Work(Ioff),1,Ncol,1,Irrepx,199+Ispin)
           Ioff = Ioff + Nrow*Ncol
        Enddo 

#ifdef _DEBUG_LVLM
       Write(6,*)
       Write(6,"(a)") " ABAB RPA T2 vectors"
       Write(6, "(5(F15.10,1x))") (Work(It2_s+i-1), i=1,Nsize_s)
#endif
        Iscr   = It2_t + Nsize_t
        Iend   = Iscr  + Nsize_t
        If (Iend .Ge. Length) Call Insmem("Process_rpat2",
     +                                     Iend,Maxcor)

        Call Sstgen(Work(It2_s),Work(It2_t),Nsize_t,Vrt(1,1),Pop(1,1),
     +              Vrt(1,2),Pop(1,2),Work(Iscr),1,"1324")
#ifdef _DEBUG_LVLM
       Write(6,*)
       Write(6,"(a)") " Reordered ABAB RPA T2 vectors"
       Write(6, "(5(F15.10,1x))") (Work(It2_t+i-1), i=1,Nsize_t)
       call checksum("T2-ABAB",Work(It2_t),Nsize_t)
#endif 
        Call Putall(Work(It2_t),Nsize_t,1,46)

      Endif 
C Add a record to JOBARC to flag that the RPA T2 vectors are available
C on the list 44,45 and 46.

      Call Putrec(20,"JOBARC","RPAT2VEC",Ione,One)

      Return
      End 
