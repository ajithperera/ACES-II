      Subroutine Form_t2aaaa(T2ab,Work,Maxcor,Lenab,Lenaa)

      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Maxcor),T2ab(Lenab)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "sym.com"

      Data Ione,Inul /1,0/
      
      Ioff   = Ione 
      Joff   = Inul
      Irrepx = Ione

      call getall(T2ab,Lenab,Irrepx,63)

#ifdef _DEBUG_LVL0
      Call checksum("T2-ABAB :",T2ab,Lenab)
#endif 
      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepx,Irrepr)
 
         Ncol_isqj = Irpdpd(Irrepr,14)
         Nrow_asqb = Irpdpd(Irrepl,15)
         Ncol_itrJ = Irpdpd(Irrepr,3)
         Nrow_atrb = Irpdpd(Irrepl,1)

         I000 = Ione
         Iend = I000 + Lenaa
         If (Iend .Gt. Maxcor) Call Insmem("form_t2aaaa",Iend,
     +                                      Maxcor)
         Call Assym2(Irrepl,Pop(1,1),Nrow_asqb,T2ab(Ioff))

         Call Sqsym(Irrepl,Vrt(1,1),Nrow_atrb,Nrow_asqb,
     +              Ncol_itrJ,Work(I000+joff),T2ab(Ioff))

         Ioff = Ncol_isqj*Nrow_asqb + Ioff
         Joff = Ncol_itrJ*Nrow_atrb + Joff
      Enddo 

      Call Dcopy(Lenaa,Work(I000),1,T2ab,1)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " Forming T2AAAA from T2ABAB"
      Call checksum("T2-AAAA :",T2ab,Lenaa)
#endif 

      Return
      End

