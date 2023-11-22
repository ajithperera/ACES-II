      Subroutine Add_t3in2t2_abij(Work,Length,Iuhf)

      Implicit Double Precision (A-H,O-Z)
      Logical Uhf

      Dimension Work(Length)

      Data Ione,onem,One /1,-1.0D0,1.0D0/

#include "sympop.com"
#include "sym.com"
#include "syminf.com"
#include "parcc.com"

#ifdef _DEBUG_LVL0
       Write(6,"(a)") "-----Entered add_t3int2_abij-----"
#endif 

      Irrepx = 1
      Uhf    = (Iuhf .NE. 0)
C Tai contributions are temporarily stored in 10,12 coulmns of list 
C 90. Add them to the current increment. 

      Lenhp_aa = Irpdpd(Irrepx,9)
      Lenhp_bb = Irpdpd(Irrepx,10)
      Maxln = Max(Lenhp_aa,Lenhp_bb)
      I000 = Ione
      I010 = I000 + Maxln
      Iend = I010 + Maxln
      If (Iend .Gt. Length)Call Insmem ("add_external_t3_to_h12",Maxln,
     +                                  Length)
      Call Getlst(Work(I000),1,1,1,5,93)
#ifdef _DEBUG_LVL0
      call checksum("T1-AA   :",Work(I000),lenhp_aa)
#endif 
      Call Getlst(Work(I010),1,1,1,3,90)
      Call Daxpy(Lenhp_aa,One,Work(I000),1,WOrk(I010),1)
      Call Putlst(Work(I010),1,1,1,3,90)

      If (Iuhf .Ne. 0) Then
         Call Getlst(Work(I000),1,1,1,6,93)
#ifdef _DEBUG_LVL0
      call checksum("T1-BB   :",Work(I000),lenhp_bb)
#endif 
         Call Getlst(Work(I010),1,1,1,4,90)
         Call Daxpy(Lenhp_bb,One,Work(I000),1,WOrk(I010),1)
         Call Putlst(Work(I010),1,1,1,4,90)
      Endif

C The <ab|ij> contributions from T3 are temporarliy stored in
C list 61-63 (T2 residual list). Add them to <ab||ij> list

      Length_61=Idsymsz(IRREPX,ISYTYP(1,61),ISYTYP(2,61))
      If (Iuhf .NE. 0) Length_62=Idsymsz(IRREPX,ISYTYP(1,62),
     +                                   ISYTYP(2,62))
      Length_63=Idsymsz(IRREPX,ISYTYP(1,63),ISYTYP(2,63))

      Maxln = Max(Length_61,Length_62,Length_63)
      I000 = Ione 
      I010 = I000 + Maxln
      Iend = I010 + Maxln
      If (Iend .Gt. Length) Call Insmem("add_external_t3_to_h12",
     +                                   Iend,Length)

      If (Uhf) Then
         Call Getall(Work(I000),Length_61,Irrepx,147)
#ifdef _DEBUG_LVL0
      call checksum("T2-AAAA :",Work(I000),Length_61)
#endif 
         Call Getall(Work(I010),Length_61,Irrepx,61)
         Call Daxpy(Length_61,One,Work(I000),1,Work(I010),1)
         Call Putall(Work(I010),Length_61,Irrepx,61)

         Call Getall(Work(I000),Length_62,Irrepx,148)
#ifdef _DEBUG_LVL0
      call checksum("T2-BBBB :",Work(I000),Length_62)
#endif 
         Call Getall(Work(I010),Length_62,Irrepx,62)
         Call Daxpy(Length_62,One,Work(I000),1,Work(I010),1)
         Call Putall(Work(I010),Length_62,Irrepx,62)
      Endif

      Call Getall(Work(I000),Length_63,Irrepx,149)
#ifdef _DEBUG_LVL0
      call checksum("T2-ABAB :",Work(I000),Length_63)
#endif 
      Call Getall(Work(I010),Length_63,Irrepx,63)
      Call Daxpy(Length_63,One,Work(I000),1,Work(I010),1)
      Call Putall(Work(I010),Length_63,Irrepx,63)

C This block is not needed since all the redundent lists are 
C generated from ABAB (list 63).

#ifdef _NOSKIP
      If (Iuhf .eq. 0) Then
         Memleft = Length - Iend 
         Call Form_t2aaaa(Work(I000),Work(Iend),Memleft,Length_63,
     &                    Length_61)
         Call Putall(Work(I000),Length_61,Irrepx,61)
      Endif 
#endif 

      Return
      End
