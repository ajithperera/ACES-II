      Subroutine form_new_c2(Work,Length,Iuhf,Irrepx,Listz0,Listd0)
      
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Length)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "sym.com"

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,"(a,1x,i2)") " IRREPX:",Irrepx
#endif 
C ABAB Block

#ifdef _DEBUG_LVL0
      e_ab = 0.0D0
#endif 
      Ioff = 1 
      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)

         Listz = Listz0 + 2
         Listd = Listd0 + 3
         Nrow  = Irpdpd(Irrepl,Isytyp(1,Listz))
         Ncol  = Irpdpd(Irrepr,Isytyp(2,Listz))
         Nsize = Nrow * Ncol

         I000  = 1
         I010  = I000  + Nsize
         I020  = I010  + Nsize
         IScr1 = I020  + Nsize
         Iscr2 = Iscr1 + Nsize
         Iend  = Iscr2 + Nsize 
         If (Iend .Ge. Length) Call Insmem("form_new_c2",
     +                                      Iend,Length)
         Call Getlst(Work(I000),1,Ncol,1,Irrepr,Listz)
         Call Getlst(Work(I010),1,Ncol,1,Irrepr,Listd)

#ifdef _DEBUG_LVL0
      write(6,"(a,3i4)") "List-z,Ncol,Nrow        :",Listz,Ncol,Nrow
      call checksum("DC2_to_C2(AB):",Work(I000),nrow*ncol,s)
      call checksum("DC2_to_D2(AB):",Work(I010),nrow*ncol,s)
      Call Vecprd(Work(I000),Work(I010),Work(I020),Nsize)
      e_ab = e_ab + Ddot(Nsize,Work(I000),1,Work(I020),1)
#endif
         Call Vecprd(Work(I000),Work(I010),Work(I000),Nsize)
C         If (Iuhf .Eq. 0) Then
C            Call Spinad1(Irrepr,Pop,Nrow,Work(Ioff),Work(Iscr1),
C     +                   Work(Iscr2))
C            Ioff = Ioff + Nrow*Ncol
C         Endif
         Call Putlst(Work(I000),1,Ncol,1,Irrepr,Listz)

      Enddo
#ifdef _DEBUG_LVL0
      e =  e_ab
      write(6,*)
      Write(6, "(a,1x,F15.9)") "(C1*W)^2/D     = ", e
      write(6,*)
#endif
      If (Iuhf .eq. 0) Return 

C AAAA and BBBB blocks

#ifdef _DEBUG_LVL0
         e_aa = 0.0D0
         e_bb = 0.0D0
#endif 
      Do Ispin = 1, 1+Iuhf
         Do Irrepr = 1, Nirrep
            Irrepl = Dirprd(Irrepr,Irrepx)

            Listz = Listz0 + Ispin - 1
            Listd = Listd0 + Ispin 
            Nrow  = Irpdpd(Irrepl,Isytyp(1,Listz))
            Ncol  = Irpdpd(Irrepr,Isytyp(2,Listz))
            Nsize = Nrow * Ncol

            I000 = 1
            I010 = I000 + Nsize 
            I020 = I010 + Nsize
            Iend = I020
            If (Iend .Ge. Length) Call Insmem("form_new_c2",
     +                                         Iend,Length)

            Call Getlst(Work(I000),1,Ncol,1,Irrepr,Listz)
            Call Getlst(Work(I010),1,Ncol,1,Irrepr,Listd)

#ifdef _DEBUG_LVL0
      write(6,"(a,3i4)") "List-z,Ncol,Nrow        :", Listz,Ncol,Nrow
      if (ispin.eq.1) call checksum("DC2_to_C2(AA):",Work(I000),
     +                               nrow*ncol,s)
      if (ispin.eq.1) call checksum("DC2_to_D2(AA):",Work(I010),
     +                               nrow*ncol,s)
      if (ispin.eq.2) call checksum("DC2_to_C2(BB):",Work(I000),
     +                               nrow*ncol,s)
      if (ispin.eq.2) call checksum("DC2_to_D2(BB):",Work(I010),
     +                               nrow*ncol,s)
      if (ispin .eq. 1) then
         Call Vecprd(Work(I000),Work(I010),Work(I020),Nsize)
         e_aa = e_aa + Ddot(Nsize,Work(I000),1,Work(I020),1)
      else 
         Call Vecprd(Work(I000),Work(I010),Work(I020),Nsize)
         e_bb = e_bb + Ddot(Nsize,Work(I000),1,Work(I020),1)
      endif 
#endif 
            Call Vecprd(Work(I000),Work(I010),Work(I000),Nsize)
            Call Putlst(Work(I000),1,Ncol,1,Irrepr,Listz)
         Enddo 
      Enddo 

#ifdef _DEBUG_LVL0
      e = e_aa + e_bb + e_ab
      write(6,*) 
      Write(6, "(a,3(1x,F15.9))") "e_ab,e_aa,e_bb = ", 
     +                              e_ab,e_aa,e_bb
      Write(6, "(a,1x,F15.9)") "(C1*W)^2/D     = ", e
      write(6,*)
#endif 

      Return
      End
