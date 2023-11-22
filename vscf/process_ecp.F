      Subroutine process_ecp(Ecpint, Onehamil, Ldim, Nbas, Nirrep, 
     &                       Nbfirr)
C
      Implicit Double Precision (A-H, O-Z)

#include "machsp.com"
#include "baslims.par"
C
      Dimension Ecpint(Ldim), Onehamil(Ldim), Ao2so(Mxcbf*Mxcbf),
     &          Tmpmat1(Mxcbf*Mxcbf), Tmpmat2(Mxcbf*Mxcbf),
     &          Nbfirr(Nirrep)
 
      Call Getrec(20, 'JOBARC','NAOBASFN', 1, Naobfns)
      Length = Naobfns*(Naobfns + 1)/2
      
      Call Getrec(20, "JOBARC", "ECP1INTS", Length*Iintfp, Tmpmat2)
      Call Expnd2(Tmpmat2, Tmpmat1, Naobfns)

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "@-VSCF The ECP integrals"
      Call output(Tmpmat2, 1, Naobfns, 1, Naobfns, Naobfns, 
     &            Naobfns, 1)
#endif

      Call Getrec (20, "JOBARC", "CMP2ZMAT", Nbas*Naobfns*Iintfp,
     &             Ao2So)

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "@-VSCF The CMP2ZMAT"
      Call output(Ao2So, 1, Nbas, 1, Naobfns, Nbas, 
     &            Naobfns, 1)
#endif
      
      Call Xgemm("N", "N", Naobfns, Nbas, Naobfns, 1.0D0, Tmpmat1,
     &            Naobfns, Ao2So, Naobfns, 0.0D0, Tmpmat2, Naobfns)

      Call Xgemm("T", "N", Nbas, Nbas, Naobfns, 1.0D0, Ao2so,
     &            Naobfns, Tmpmat2, Naobfns, 0.0D0, Tmpmat1, Naobfns)

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "@-VSCF Sym. adapted  ECP integrals"
      Call output(Tmpmat1, 1, Naobfns, 1, Naobfns, Naobfns, 
     &            Naobfns, 1)
#endif
      Call sympack(Tmpmat1, Ecpint, Naobfns, Nirrep, Nbfirr)

#ifdef _DEBUG_LVM
      Write(6,"(a)") "The symmetry packed ECP integrals"
      Write(6, "(6(1x,F10.7))") (Ecpint(I), I=1, Ldim)
      Write(6,*) 
      Write(6,"(a)") "The Onehamil before ECP add"
      Write(6, "(6(1x,F10.7))") (Onehamil(I), I=1, Ldim)
#endif

      Call Daxpy(Ldim, 1.0D0, Ecpint, 1, Onehamil, 1)

#ifdef _DEBUG_LVLM
      Write(6,*) 
      Write(6,"(a)") "The Onehamil after ECP add"
      Write(6, "(6(1x,F10.7))") (Onehamil(I), I=1, Ldim)
#endif
    
      Return
      End
     
             
    
