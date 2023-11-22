#include "flags.h"
      Subroutine Ecp_dgrad(Dens_tao, Naofns, Ncoord, Nbfirr, ECP_grd)
C
      Implicit Double Precision(A-h, O-Z)
C
#include "machsp.com"
#include "mxatms.par"
#include "baslims.par"
#include "flags.com"
c
      Logical Spherical 

      Dimension Dens_tao(Naofns*(Naofns+1)/2), 
     &          Dens_fao(Mxcbf*Mxcbf), Grd_xyz(3,Mxatms),
     &          Mxpair(1), Map2zmat(Mxatms), Ecp_grd(3,Mxatms),
     &          Idexatom(Mxatms), Accord(3, Mxatms),
     &          Reorder_Mat(Mxcbf*Mxcbf), Tmpmat(Mxcbf*Mxcbf),
     &          Nbfirr(8)
    
      Common /Symminfo/ Ntotatoms, Idexatom, Acoord

      Call GetRec (20, "JOBARC", 'NIRREP  ', 1, Nirrep)
      Call Getrec(20, 'JOBARC','NATOMS  ', 1, Natoms) 
      call Getrec(20, 'JOBARC','MAP2ZMAT', Natoms, Map2zmat) 
      Call Getrec(20, 'JOBARC','NAOBASFN', 1, Naobfns)
      Call Getrec(20, "JOBARC", "NBASTOT ", 1, Nbfns)

#ifdef _ASSERT
      if (Naofns .Ne. Naobfns) then
         print *, '@REDUCE: Assertion failed.'
         print *, '   Nbas = ', Naofns
         print *, ' Naobf  = ', Naobfns
         call aces_exit(1)
      end if
#endif
#ifdef _DEBUG_LVLM1
      Write(6,"(a,2(1x,i4))") "Nbfns and NAOBASFN: ", Nbfns, Naobfns
      Write(6,*) "@-entry AO basis density matrix"
      Write(6, "(4(2x,F12.10))") (dens_tao(i), i=1, nbfns*(nbfns+1)/2)
#endif
      Call Dzero(Dens_fao, Naobfns*Naobfns)
      Call Sym_unpack(Dens_fao,Dens_tao,Naobfns,Nirrep,Nbfirr)

#ifdef _DEBUG_LVLM1
      Write(6,"(a)")"The Sym. Packed expanded AO basis density matrix"
      call debugout(Dens_fao, naobfns)
#endif

      Spherical = (iFlags(h_IFLAGS_spherical).Eq.1)
      If (Spherical) Then
         Call Getrec(20, "JOBARC", "CART2CMP", Nbfns*Naobfns*Iintfp,
     &                Reorder_Mat)
         Call Xgemm("N", "N", Nbfns, Naobfns, Naobfns, 1.0D0, 
     &               Reorder_Mat, Nbfns, Dens_fao, Naobfns, 
     &               0.0D0, Tmpmat, Nbfns)
         Call Xgemm("N", "T", Nbfns, Nbfns, Naobfns, 1.0D0, Tmpmat, 
     &               Nbfns, Reorder_Mat, Nbfns, 0.0D0, Dens_fao, 
     &               Naobfns)
      Endif 

#ifdef _DEBUG_LVL0
      Write(6,"(a)")"After transformation to spherical"
      call debugout(Dens_fao, nbfns)
#endif
  
      Call Getrec(20, "JOBARC", "CMP2ZMAT", Nbfns*Naobfns*Iintfp,
     &            Reorder_Mat)

      Call Xgemm("N", "N", Naobfns, Nbfns, Nbfns, 1.0D0, Reorder_Mat, 
     &            Naobfns, Dens_fao, Naobfns, 0.0D0, Tmpmat, Naobfns)
      Call Xgemm("N", "T", Naobfns, Naobfns, Nbfns, 1.0D0, Tmpmat, 
     &            Naobfns, Reorder_Mat, Naobfns, 0.0D0, Dens_fao, 
     &            Naobfns)

#ifdef _DEBUG_LVL0
      Write(6,"(a)")"The reordered AO basis density matrix"
      call debugout(Dens_fao, nbfns)
#endif
      Ntotatoms = Natoms
#ifdef _OLDECP
      call Getrec(20,'JOBARC','ATOMINDX', Natoms, Idexatom)
#endif
      call Getrec(20, 'JOBARC','COORD   ', 3*Ntotatoms*Iintfp, 
     &            Acoord) 
      Call Dzero(Grd_xyz, 3*Mxatms)
C
#ifdef _OLDECP    
      Call ecp_grdint(Mxpair, 1, Dens_fao, Naobfns, Grd_xyz)
#else
      Call ecp_grdint_main(Ntotatoms, Dens_fao, Naobfns, Grd_xyz)
#endif
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "Gradients returned from ecp_grdint_main"
      Write(6,"(3(1x,F10.7))") ((Grd_xyz(i,j),i=1,3),j=1,natoms)
#endif

      Do Iatms = 1, Natoms

         If (Map2zmat(Iatms) .NE. 0) Then

            Do Ixyz = 1, 3
           
               Ecp_grd(Ixyz, Iatms) = Grd_xyz(Ixyz, Map2zmat(Iatms))*
     &                                2.0D0
            Enddo

         Endif
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6, "(a)") "The gradient due scalar-ECP part"
      Do iatms = 1, Natoms
         Write(6, "(3(1x, F15.12))") (Ecp_grd(Ixyz,iatms), ixyz=1,3)
      Enddo 
#endif

      Call Tran_c2sgrad(Ecp_grd, Grd_xyz, Ncoord)
      Call Dcopy(Ncoord, Grd_xyz, 1, Ecp_grd, 1)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6, "(a)") "The sym. ada. gradient due scalar-ECP part"
      Do iatms = 1, Natoms
         Write(6, "(3(1x, F15.12))") (Ecp_grd(Ixyz,iatms), ixyz=1,3)
      Enddo 
#endif
      Return
      End
