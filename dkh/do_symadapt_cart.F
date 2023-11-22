      Subroutine Do_symadapt_cart(Nbfns, Ovlp, Ao2so, Tmp1, Tmp2,
     +                            Itriang_length)
  
      Implicit None

#include "machsp.com"  

      Integer Nbfns, Nbfirr(8), Irrep, Nirrep, Idex, I, J
      Integer Itriang_off(8)
      Integer Isquare_off(8)
      Integer Isquar2_off(8)
      Integer Ioff, Itriang_length

      Double Precision Ovlp(Nbfns,Nbfns), Tmp1(Nbfns,Nbfns),
     &                 Tmp2(Nbfns*Nbfns), Ao2So(Nbfns*Nbfns)

      Call Getrec (20, "JOBARC", "CMP2ZMAT", Nbfns*Nbfns*Iintfp,
     &             Ao2So)
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The AO to SO transformation matrix"
      Call output(Ao2so, 1, Nbfns, 1, Nbfns, Nbfns,
     &            Nbfns, 1)
C      Write(6,"(a)") "The Incomming PVP integrals (scalled)"
C      Call output(Ovlp, 1, Nbfns, 1, Nbfns, Nbfns, Nbfns,1)
#endif 
      Call Xgemm("N", "N", Nbfns, Nbfns, Nbfns, 1.0D0, Ovlp,
     &            Nbfns, Ao2so, Nbfns, 0.0D0, Tmp2, Nbfns)
      Call Xgemm("T", "N", Nbfns, Nbfns, Nbfns, 1.0D0, Ao2so,
     &            Nbfns, Tmp2, Nbfns, 0.0D0, Tmp1, Nbfns)

#ifdef _DEBUG_LVL0
      Write(6,"(a,a)") "Sym. adapted  PVP integrals",
     +                 " (scalled)"
      Call output(Tmp1, 1, Nbfns, 1, Nbfns, Nbfns, Nbfns,1)
#endif

C The following block of code has been obsoleted. 

#ifdef _NOSKIP
      Call Getrec (20, "JOBARC", "AO2SO   ", Nbfns*Nbfns*Iintfp,
     &             Ao2So)
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The AO to SO transformation matrix"
      Call output(Ao2so, 1, Nbfns, 1, Nbfns, Nbfns,
     &            Nbfns, 1)
      Write(6,"(a)") "The Incomming PVP integrals (scalled)"
      Call output(Ovlp, 1, Nbfns, 1, Nbfns, Nbfns, Nbfns,1)
#endif
      Call Xgemm("N", "N", Nbfns, Nbfns, Nbfns, 1.0D0, Ovlp,
     &            Nbfns, Ao2so, Nbfns, 0.0D0, Tmp2, Nbfns)

      Call Getrec (20, "JOBARC", "AO2SOINV", Nbfns*Nbfns*Iintfp,
     &             Ao2So)

      Call Xgemm("N", "N", Nbfns, Nbfns, Nbfns, 1.0D0, Ao2so,
     &            Nbfns, Tmp2, Nbfns, 0.0D0, Tmp1, Nbfns)
   
#ifdef _DEBUG_LVL0
      Write(6,"(a)")"Sym. adapted PVP integrals with contaminants"
      Call output(Tmp1, 1, Nbfns, 1, Nbfns, Nbfns,Nbfns, 1)
#endif
#endif 
      Call Getrec(20, "JOBARC", "NIRREP  ", 1, Nirrep)
      Call Getrec(20, "JOBARC", "NUMBASIR", Nirrep, Nbfirr)
C
      Idex = 0
      Do Irrep = 1, Nirrep
         Do I = 1 + Idex, Nbfirr(Irrep) + idex
              Do J = 1+ Idex, Nbfirr(Irrep) + idex
              Ovlp(j, i) = Tmp1(j,i)
           Enddo
        Enddo
           Idex = Idex + Nbfirr(Irrep)
      Enddo

      Itriang_off(1) = 1
      Isquare_off(1) = 1
      Isquar2_off(1) = 1
      Do Irrep = 1, Nirrep-1
         Itriang_off(Irrep+1) = Itriang_off(Irrep) + (Nbfirr(Irrep)+1)*
     &                          Nbfirr(Irrep)/2
         Isquare_off(Irrep+1) = Isquare_off(Irrep) + Nbfirr(Irrep) 
         Isquar2_off(Irrep+1) = Isquar2_off(Irrep) + Nbfirr(Irrep) *
     &                          Nbfirr(Irrep)
      Enddo 

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)")"Sym. Adapted PVP integrals sym. block basis"
      Do Irrep = 1, Nirrep
      If (Nbfirr(Irrep) .ne. 0) then
      Call output(Ovlp(Isquare_off(Irrep),Isquare_off(Irrep)), 
     &            1, Nbfirr(Irrep), 1, Nbfirr(Irrep), 
     &            Nbfns, Nbfns, 1)
      endif 
      Enddo 
#endif 
      Ioff = 1 
      Idex = 0
      Do Irrep = 1, Nirrep 
         Ioff = Isquar2_off(Irrep)
         Do I =1 +Idex, Nbfirr(Irrep) + Idex
            Call Dcopy(Nbfirr(Irrep),Ovlp(1+idex,i),1,Ao2so(Ioff),1) 
            Ioff = Ioff + Nbfirr(Irrep)
         Enddo
         Idex = Idex + Nbfirr(Irrep)
      Enddo 

      Do Irrep = 1, Nirrep
         Call Squez2(AO2so(Isquar2_off(Irrep)),
     +               Tmp2(Itriang_off(Irrep)),Nbfirr(Irrep))
      Enddo 

      Itriang_length = Itriang_off(Nirrep) + (Nbfirr(Nirrep)+1)*
     +                                        Nbfirr(Nirrep)/2
      Itriang_length = Itriang_length      -  1

      Return
      End
