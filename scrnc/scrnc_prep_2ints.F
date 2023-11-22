      Subroutine scrnc_prep_2ints(Work,Maxcor,Irrepx,List_xxxx)

      Implicit Double Precision (A-H,O-Z)

#include "maxbasfn.par"

      Dimension Work(Maxcor)
      Dimension Iord(Maxbasfn*Maxbasfn) 
      Double Precision Buf(600) 
      Integer Ibuf(600)

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
      Common/Sym_pqxx_ints/Irrep_aos(Maxbasfn),
     &                     Irrep_ao_pairs(Maxbasfn*Maxbasfn)

#include "syminf.com"
#include "sympop.com"
#include "sym.com"
    
      Data Ione /1/

      Call Getrec(20,"JOBARC","NBASTOT",Ione,Nbfns)
      Ndim4 = Nbfns**4
      Ndim2 = Nbfns**2 
   
      I000 = 1
      I010 = I000 + Ndim4 
      Iend = I010

      If (Iend .Gt. Maxcor) Call Insmem("scrnc_prep_2ints",Iend,
     +                                    Maxcor)
      Call Scrnc_get2ints(Work(I000),Buf,Ibuf,Nbfns)

#ifdef _DEBUG_LVL0
      Call output(Work(I000),1,Ndim2,1,Ndim2,Ndim2,Ndim2,1)
      Call checksum("I(XX,XX),Nosym:",work(i000),Ndim4)
#endif 
      
C Integrals are arranged in I(IJ,AB) format.First symmetry pack the IJ 
C block per each AB block.

      I000 = 1
      I010 = I000 + Ndim4
      I020 = I010 + Ndim4

      Ioff1  = I000
      Ioff2  = I010
      Maxcor = Maxcor - Ioff2 
      Irrepx = 1

      Ithru = 0
      Do Irrep = 1, Nirrep
         Do Iaos = 1, Iaopop(Irrep)
            Ithru = Ithru + 1
            Irrep_aos(Ithru) = Irrep 
         Enddo
      Enddo

#ifdef _DEBUG_LVLM
      write(6,"(a,10(1x,I4))")"Irreps of AOs      :",(Irrep_aos(i),i=1,
     +                       Nbfns)
#endif 
      Ithru = 0
      Do k = 1, Nbfns
         Irrep_k = Irrep_aos(k)
         Do l= 1, Nbfns 
            Irrep_l = Irrep_aos(l)
            Irrep = Dirprd(Irrep_l,Irrep_k)
            Ithru = Ithru + 1
            Irrep_ao_pairs(Ithru) = Irrep
         Enddo
      Enddo 
#ifdef _DEBUG_LVL0
      write(6,"(a,10(1x,I4))")"Irreps of AOs_pairs:",(Irrep_ao_pairs(i),
     +                       i=1,Nbfns*Nbfns)
      Write(6,*)
#endif 

      Call Scrnc_reord_cols(Work(Ioff1),Work(Ioff2),Nbfns,
     +                      Irrep_ao_pairs,Iord)
#ifdef _DEBUG_LVLM
      call scrnc_ao2mo(Work(Ioff1),Work(ioff2),Nbfns,Iuhf)
      Stop
#endif 

      Do Ipairs = 1, Ndim2

         Irrepr = Irrep_ao_pairs(Ipairs) 
         Call Scrnc_sympack_aoints(Work(Ioff1),Work(Ioff2),Nbfns,
     +                          Irrepr)
         Ioff1 = Ioff1 + Ndim2 
         Ioff2 = Ioff2 + Ndim2
      Enddo 
#ifdef _DEBUG_LVL0
      Write(6,"(a,a)") "@-scrnc_xxxx_pqxx,The <bra|ket> sym. packed",
     +           " integrals."
      Call output(Work(I010),1,Ndim2,1,Ndim2,Ndim2,Ndim2,1)
#endif 
      
      Ioff1 = I000
      Ioff2 = I010
      Do Ipairs = 1, Ndim2 
         Irrepr = Irrep_ao_pairs(Ipairs)
         Irrepl = Dirprd(Irrepr,Irrepx)
         Ndim = Irpdpdao(Irrepl) 
         Call Dcopy(Ndim,Work(Ioff2),1,Work(Ioff1),1)
         Ioff1 = Ioff1 + Ndim 
         Ioff2 = Ioff2 + Ndim2
      Enddo

      Iput = I000
      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)
         Ncol = Irpdpdao(Irrepr)
         Nrow = Irpdpdao(Irrepl)
         Call Putlst(Work(IPut),1,Ncol,1,Irrepr,List_xxxx)
#ifdef _DEBUG_LVL0
C         Write(6,"(a)") "@-scrnc_xxxx_pqxx,I(XX,XX)"
C         call output(work(iput),1,nrow,1,ncol,nrow,ncol,1)
         call checksum("I(XX,XX)Sym   :",work(iput),Nrow*Ncol)
#endif 
         Iput = Iput + Nrow * Ncol 
      Enddo 

      Return
      End
