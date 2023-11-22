#include "flags.h"
      Subroutine libint_driver(Dwork,Iwork,Dmaxcor,Imaxcor,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Integer Dmaxcor
 
#include "flags.com"
#include "libints.par"

      Dimension Dwork(Dmaxcor),Iwork(Imaxcor)
      Character*4 Atmnam(Max_centers)
      Character*4 Inp_name 
      Logical Spherical,Cartesian 
      Integer Dleft 
      Character*1 Kka(3,3)

      Data Ione /1/

C This call simply assign the number of atoms, total number of  shells (the 
C sum of maximum angular momentum per center) and the maximum number of 
C primitives of shells.
C In ACES II/III these data are in the MOL file (written from joda).

      Inp_name = "MOL"
      Call Libint_get_basics(Inp_name,Max_centers,Max_shells,
     &                       Max_Prims_ps,Ncenters,Nshells,
     &                       Nprims_ps,Nunique,Nsymgen,Kka,Iunits)

      If (Max_centers .Lt. Ncenters) Then 
         Write(6,"(a,I4,a,I4)") " The number of centers ",Ncenters, 
     &                        " is greater than the maximum allowed ", 
     &                          Max_centers
         Call Errex 
      Endif 

      If (Max_shells .Lt. Nshells) Then 
         Write(6,"(a,I4,a,I4)") " The number of shells ",Nshells, 
     &                        " is greater than the maximum allowed ", 
     &                          Max_shells 
         Call Errex 
      Endif 

      If (Max_prims .Lt. Nprims) Then 
         Write(6,"(2a,I4,a)") " The largest number of primitives",
     &                     " per shell ", Nprims_Ps, " is greater than"
         Write(6,"(a,I4)") " the maximum allowed ", Max_prims_ps
         Call Errex 
      Endif 

C Construct the transformation to form symmetry adapted AOs.

      Call Built_symtran(Iwork,Dwork,Dmaxcor,Imaxcor,Nunique,
     &                   Nsymgen,Kka)

      J000 = Ione
      I000 = Ione

      J010 = J000 + Ncenters
      J020 = J010 + Ncenters 
      J030 = J020 + Ncenters 
      J040 = J030 + Nshells 
      J050 = J040 + Nshells 
      J060 = J050 + Nshells 
      J070 = J060 + Nshells 
      J080 = J070 + Nshells 
      J090 = J080 + Ncenters*Nshells*Nshells 
      J100 = J090 + Ncenters*Nshells*Nshells
      J110 = J100 + Ncenters
      J120 = J110 + Ncenters*Nshells 
      J130 = J120 + Nshells 
      Jend = J130 + Nshells 

      If (Jend .Gt. Imaxcor) Call Insmem("libint_drive",Jend,Imaxcor)
 
      I010 = I000 + Ncenters 
      I020 = I010 + 3*Ncenters 
 
      Iecp       = 1
      Ispherical = 0
      Spherical  = .False.
      cartesian  = .False.
 
      If (IFLAGS(62) .EQ. 1) Then
          Ispherical = 1
          Spherical  = .True.
      Else
          Cartesian  = .True.
      Endif

      Call Libint_getbasis_data(Inp_name,Ncenters,Nshells,
     &                          Ncenters,Nshells,IWork(J000),
     &                          Cartesian,Max_prims_ps,Nprims,Lnp1,
     &                          Lnpo,Iwork(J010),Iwork(J020),Nbasis,
     &                          Naobasis,Iwork(J030),Iwork(J040),
     &                          Iwork(J050),Iwork(J060),
     &                          Iwork(J070),Vnn,Nshells,Nshells,Iecp,
     &                          Iwork(J080),Iwork(J090),Iwork(J100),
     &                          Iwork(J110),Iwork(J120),Iwork(J130),
     &                          Atmnam,Dwork(I000),Dwork(I010),
     &                          IWork(Jend))

      Nalpha   = 0
      Npcoef   = 0
      Max_cfps = 0
      Max_pfps = 0
      Do I = 1, Nshells
         Nalpha = Nalpha + Iwork(J040+I-1)
         Npcoef = Npcoef + Iwork(J040+I-1)*Iwork(J030+I-1)
         Max_cfps = Max(Max_cfps,Iwork(J030+I-1)) 
         Max_pfps = Max(Max_pfps,Iwork(J040+I-1)) 
      Enddo 

      J140 = Jend
      J150 = J140 + Ncenters
      J160 = J150 + Ncenters 
      J170 = J160 + Nshells
      J180 = J170 + Ncenters 
      J190 = J180 + Nshells 
      J200 = J190 + Nshells 
      J210 = J200 + Nshells 
      J220 = J210 + Nshells 
      Jend = J220 + Nshells 

      Do I = 1, Ncenters
         IWork(J140+I-1) = 1
         IWork(J150+I-1) = i
      Enddo 

      If (Jend .Gt. Imaxcor) Call Insmem("libint_drive",Jend,Imaxcor)

      I030 = I020 + Nalpha 
      I040 = I030 + NPcoef 
      I050 = I040 + 3*Nshells 
      I060 = I050 + 3*Ncenters 
      I070 = I060 + Nbasis 
      I080 = I070 + 4*Max_pfps
      IEnd = I080 + 4*Max_pfps*Max_cfps

      If (Iend .Gt. Imaxcor) Call Insmem("libint_driver",Iend,Dmaxcor)

      Call Set_atm_4Shell(Iwork(J160),Ncenters,Nshells,Iwork(J100),
     &                    Iwork(J110))

      Call Libint_readbasis_data(Inp_name,Ncenters,Ncenters,
     &                           Iwork(J140),Iwork(J150),Cartesian,
     &                           Nprims,Lnp1,Lnpo,Iwork(J010),Nbasis,
     &                           Dwork(I020),IWork(J180),DWork(I030),
     &                           IWork(J190),Max_centers,Atmnam,
     &                           Dwork(I040),Dwork(I050),Iwork(J050))

      Call Print_basis_info(Nshells,Ncenters,Nalpha,Npcoef,Naobasis,
     &                      Nbasis,Lnp1,Lnpo,Nprims,Iwork(J030),
     &                      Iwork(J040),Iwork(J100),Iwork(J080),
     &                      Iwork(J090),Iwork(J110),
     &                      Iwork(J120),Iwork(J130),IWork(J060),
     &                      Iwork(J180),Iwork(J190),Dwork(I020),
     &                      Dwork(I030),Dwork(I040),Dwork(I050),
     &                      Nshells,IecP)
      J230 = Jend 
      J240 = J230 + Nprims
      J250 = J240 + Nprims
      J260 = J250 + Nbasis 
      J270 = J260 + 4*Max_cfps
      Jend = J270 + 4*Max_cfps
     
      If (Jend .Gt. Imaxcor) Call Insmem("libint_drive",Jend,Imaxcor)

      Dleft = Dmaxcor - Iend 
      Ileft = Imaxcor - Jend 

      Call Libint_prep_4ints(Iwork(J060),Iwork(J030),Iwork(J040),
     &                       Dwork(I030),Iwork(J190),IWork(J200),
     &                       Iwork(J210),Iwork(J220),Iwork(J230),
     &                       Iwork(J240),Nshells,Npcoef,Nprims,
     &                       Cartesian,SPherical)

      Call Aces_to_erd(Iwork(J200),Iwork(J060),Nshells,Ispherical,
     &                 Iwork(J250),Dwork(I060))
      

      Dleft = Dmaxcor - Iend 
      Ileft = Imaxcor - Jend 

      Call Libint_ints_main(Nshells,Ncenters,Nalpha,Npcoef,Nbasis,
     &                      Nprims,Spherical,Cartesian,Iwork(J060),
     &                      Iwork(J030),Iwork(J040),Iwork(J180),
     &                      Iwork(J190),IWork(J200),Iwork(J210),
     &                      Iwork(J220),Iwork(J230),Iwork(J240),
     &                      Iwork(J250),Iwork(J160),Dwork(I020),
     &                      Dwork(I030),Dwork(I040),Dwork(I050),
     &                      Dwork(I060),Dwork(I070),Dwork(I080),
     &                      Iwork(J260),IWork(J270),Dwork(Iend),
     &                      Iwork(Jend),Dleft,Ileft,Max_cfps,
     &                      Max_pfps)

      Return
      End  
