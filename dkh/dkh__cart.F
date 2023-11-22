#include "flags.h"
      Subroutine Dkh__cart(Work,Mxdcor,Iuhf,Spherical,Cartesian,
     +                     Symmetry,Contract)
      
      Implicit Double Precision (A-H, O-Z)

#include "icore.com"
#include "istart.com"
#include "ecplims.par"
#include "cvwecp.h"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

      Parameter (Maxbfns = 1000,Mxicor= 1000000)
      Character*10 INP_Fname
      Character*4 ATMNAM(Max_centers), Namat(Max_centers)
      Logical Spherical, Cartesian, Sym_adapt 
      Logical Contract,Symmetry
      Integer Nelement,Dkh_order
    
      Integer CCbeg(Maxbfns), CCend(Maxbfns), Zmax, Zhrr, 
     &        End_nfps, Val1(4), Val2(4), Angmom, Atm_4shell,
     &        Erd_index(Maxbfns)

      DIMENSION WORK(Mxdcor), Iwork(Mxicor), Erd_scale(Maxbfns)

      DIMENSION NFCT(max_centers), NAOATM(max_centers),
     &          NUFCT(max_centers), NAOUATM(max_centers),
     &          iqmstr(max_centers),jcostr(max_centers,Maxjco),
     &          nucstr(max_centers,Maxang,Maxjco),
     &          nrcstr(max_centers,Maxang,Maxjco),
     &          Jstrt(max_shells),Jrs(max_shells),
     &          NCFpS(max_shells),Angmom(max_shells),
     &          NPFpS(max_shells), Atomic_label(max_shells),
     &          Charge(max_centers),Nspc(max_centers), 
     &          NPOP(max_centers),
     &          IREORDER(max_centers), 
     &          COORD(3,max_shells),Cord(3,max_centers),
     &          Ivshloff(max_shells),Nfps(max_shells),
     &          End_nfps(max_shells),Npfs(max_shells),
     &          End_nfps_uc(max_shells),
     &          ivAngMom(max_shells),Ixalpha(max_shells),
     &          Ixpcoef(max_shells),Indx_cc(max_shells),
     &          Ixshells(max_shells),Atm_4shell(max_shells),
     &          Temp(max_shells)
C
      Iecp = 1
      Ispherical = 0
      If (Spherical) Ispherical = 1
      INP_Fname = "MOL_ERDOED"
      Dkh_order =   IFLAGS2(h_IFLAGS2_dkh_order)
      
      Call SIMPLE_INSPECT_MOL(INP_Fname, max_centers, max_shells,
     &                        ncenters, nshells, nspc, cartesian, 
     &                        ITFCT, LNP1, lnpo, nfct, nufct, 
     &                        nbasis, NAOBASIS, nCFpS, 
     &                        nPFpS, NAOATM, angmom, atomic_label, 
     &                        vnn, Maxang, Maxjco, Iecp,
     &                        NUcstr, Nrcstr, Iqmstr, Jcostr, 
     &                        Jstrt, Jrs, Atmnam, Charge)

#ifdef _DEBUG_LVLM
      Write(6,"(a,a,7(1x,I4))") "natoms,nshells,ITFCT,LNP1,lnpo,",
     &                          "nbasis,NAOBASIS:",ncenters,nshells,
     &                                         ITFCT,LNP1,lnpo,
     &                                        nbasis,NAOBASIS
      Write(6,*) 
      Write(6,"(a)") "The number of con. function per shell"
      Write(6,"(6(1x,I4))") (nCFpS(i), i=1, nshells)
      Write(6,"(a)") "The number of prim. function per shell"
      Write(6,"(6(1x,I4))") (nPFpS(i), i=1, nshells)
#endif 

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The Iqmstr"
      Write(6,"(6(1x,I4))") (Iqmstr(i), i=1, ncenters)
      Write(6,"(a)") "The jcostr"
      Do i=1, Ncenters
         Write(6,"(4(1x,I4))") (jcostr(i,j),j=1,4)
      Enddo
      Write(6,"(a)") "The Nucstr"
      Do i=1, Ncenters
         Write(6,"(10(1x,I4))") ((Nucstr(i,j,k),j=1,4),k=1,ITFCT)
      Enddo
      Write(6,"(a)") "The Nrcstr"
      Do i=1, Ncenters
         Write(6,"(10(1x,I4))") ((Nrcstr(i,j,k),j=1,4),k=1,ITFCT)
      Enddo
      Write(6,"(a)") "The jstrt"
      Write(6,"(6(1x,I4))") (Jstrt(i), i=1, Nshells)
      Write(6,"(a)") "The jrs"
      Write(6,"(6(1x,I4))") (Jrs(i), i=1, Nshells)
#endif 

      nalpha    = 0
      npcoef    = 0
      npcoef_uc = 0 
      do i = 1, nshells
         nalpha    = nalpha + npfps(i)
         npcoef    = npcoef + npfps(i) * ncfps(i)
         npcoef_uc = npcoef_uc + npfps(i) * npfps(i) 
         Ivangmom(i) = Angmom(i)
      enddo

      IALPHA    = 1
      IPCOEF    = IALPHA + nalpha
      IPCOEF_UC = IPCOEF + npcoef
      INEXT     = IPCOEF_UC + npcoef_uc

      do i = 1, ncenters
         npop(i)     = 1
         ireorder(i) = i
         Namat(i) = ATMNAM(i)
      enddo
      Nshll = 0
      do iat=1,Ncenters
        do ilq=1,iqmstr(iat)
          do ijco=1,jcostr(iat,ilq)
             nshll=nshll+1
             Atm_4shell(nshll)=iat
C      write(*,*)'Atm(nshll)',Atm_4shell(nshll)
           Enddo
        Enddo
      Enddo
C
      Call READ_BASIS_INFO(INP_FNAME, Ncenters, Ncenters, NPOP,
     &                     IREORDER,
     &                     CARTESIAN, ITFCT, LNP1, LNPO,
     &                     NFCT, NBASIS, WORK(IALPHA), IXALPHA,
     &                     Work(IPCOEF), IXPCOEF, MAX_CENTERS,
     &                     ATMNAM, COORD, CORD, NAOATM)

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,"(a)") "The exponents"
      Write(6,"(6(1x,F12.8))") (Work(Ialpha+i), i=0, Nalpha-1)
      Write(6,*)
      Write(6,"(a)") "The contraction coefs"
      Write(6,"(6(1x,F12.8))") (Work(Ipcoef+i), i=0, Npcoef-1)
      Write(6,"(a)") "The jcostr"
      Do i=1, Ncenters
         Write(6,"(4(1x,I4))") (jcostr(i,j),j=1,4)
      Enddo

      Write(6,*)
      Write(6,*) "Atom Coordiantes of each shell"
      Do Iatm =1, Ncenters
         Write(6, "(6(1x,F12.8))") (Cord(i,Iatm), i=1,3)
      Enddo
#endif 

      Ixpcoef(1) = 1
      do i = 1, Nshells 
          Ixpcoef(i+1) = Npfps(i) * Npfps(i) + Ixpcoef(i)
      enddo 
C
C Here the number of contracted functions are the same as primitive 
C functions. 

      Call int_gen_init(Nshells, Ivshloff, Ivangmom, Nfps, Npfs, 
     &                  Npfps,
     &                  Ixshells, Ixalpha, Ixpcoef, End_nfps, Temp, 
     &                  Ispherical)

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,*) "The Ivangmom array" 
      Write(6,"(8(1x,I4))") (Ivangmom(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The Ixshells array" 
      Write(6,"(8(1x,I4))") (Ixshells(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The Nfps array" 
      Write(6,"(8(1x,I4))") (Nfps(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The Npfs array" 
      Write(6,"(8(1x,I4))") (Npfs(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The Npfps array" 
      Write(6,"(8(1x,I4))") (Npfps(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The IVshloff array" 
      Write(6,"(8(1x,I4))") (Ivshloff(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The Ixpcoef array" 
      Write(6,"(8(1x,I4))") (Ixpcoef(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The Ixalpha array" 
      Write(6,"(8(1x,I4))") (Ixalpha(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The End_nfps array" 
      Write(6,"(8(1x,I4))") (End_nfps(i), i=1, Nshells)
#endif 

      call dzero(work(ipcoef_uc), npcoef_uc)

      Index = 0
      Do ishell=1, Nshells 
         ishel_off = Ixpcoef(ishell) 
         iprim_off = npfps(ishell)
         index = 0
         Do i=1, iprim_off
            work(ipcoef_uc+ishel_off+index-1) = 1.0D0
            Index = Index + iprim_off + 1
         Enddo 
      Enddo 
#ifdef _DEBUG_LVLM
      write(6,*)
      Write(6,"(a)") "The primitive contraction coefs"
      Write(6,"(4(1x,F10.5))") (work(ipcoef_uc +i),i=0,npcoef_uc-1)
      Write(6,*)
#endif 
C
C Here the number of contracted functions are the same as primitive 
C functions. 
C
      Call Gen_oed_to_vmol(Nfps, Ivangmom, Nshells, Spherical, 
     &                 Erd_index,ERD_Scale)

      Call setup_ccbeg(WORK(Ialpha), IXalpha, Work(Ipcoef_uc),Ixpcoef,
     &                 Npfps, Npfps, Nshells, CCbeg, CCend, Indx_cc)
C
#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,*)"The CCbeg and CCend arrays" 
      Write(6,"(10(1x,I4))") (CCbeg(i), i=1, Itfct)
      Write(6,"(10(1x,I4))") (CCend(i), i=1, Itfct)
#endif 
 
      Call Erd_scratch_mem_calc(nshells, ivangmom, npfps, npfps,
     &                          atom, Coord, Work(ialpha),
     &                          Work(ipcoef_uc), 
     &                          ixalpha, ixpcoef,
     &                          Ccbeg, Ccend, indx_cc,
     &                          spherical, Ncenters, .true., 
     &                          intmax, zmax)

      Call Erd__memory_hrr_correction(ivAngMom, nshells,
     &                                spherical, ihrr, zhrr)
     
      Imem   = iflags(36)
      Intmax = Intmax + Ihrr
      Zmax   = Zmax   + Zhrr

#ifdef _DEBUG_LVLM
      Write(6,"(a,1x,I6,1x,I10)") "ERD memmory requirments Ibuf & Dbuf", 
     &                            Intmax, Zmax
#endif 
       
      NBFNS   = NAOBASIS
      NAOBFNS = NAOBASIS

      N_Cart_fns = Naobasis 
      N_Sphe_fns = Nbasis 

      Nprims = 0
      ncnfns = 0
      Do Ishell = 1, Nshells 
         If (Spherical) Then
            lvalue = Ivangmom(Ishell) 
            Nprims = Nprims + (2*lvalue+1) *
     &               npfps(ishell) 
            Ncnfns = ncnfns  + (2*lvalue+1) * ncfps(ishell)
         Else 
            lvalue = Ivangmom(Ishell) 
            Nprims = Nprims + (lvalue+1)*(lvalue+2)/2 * 
     &               npfps(ishell) 
            Ncnfns = Ncnfns + (lvalue+1)*(lvalue+2)/2 * 
     &               ncfps(ishell) 
         Endif 
     
      Enddo 
#ifdef _DEBUG_LVL0
      Write(6,"(a,1x,I4,1x,I4)") "Nprims,Ncnfns: ", Nprims,Ncnfns
#endif 
      NBFNS = Nprims 

      Ioed_dbuf     = INext
      Ioed_ovl      = Ioed_dbuf   + Zmax 
      Ioed_kin      = Ioed_ovl    + Nprims**2
      Ioed_nai      = Ioed_kin    + Nprims**2
      Ioed_pvp      = Ioed_nai    + Nprims**2
      Icoreham      = Ioed_pvp    + Nprims**2
      Idkh          = Icoreham    + Nprims**2 
      Ireord        = Idkh        + Nprims**2
      IEnd          = Ireord      + 40*Nprims**2 
      Ileft         = Mxdcor - Iend 
#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,*)"The Erd_index and Erd_scalae arrays"
      Write(6,"(10(1x,I4))") (Erd_index(i), i=1,Nbfns)
      Write(6,"(10(1x,F5.2))") (Erd_Scale(i), i=1,Nbfns)
      Write(6,*)
      Write(6,"(a,1x,i5,1x,i5)") "Iend and Mxdcor:", Iend, Mxdcor
#endif

      If (Iend .Ge. Mxdcor) Call insmem("Dkh:", Iend, Mxdcor)
      If (Intmax .Ge. Mxicor) Call insmem("Dkh:",Intmax,Mxicor)
C
C  There is one segment in serial work (not necessary, but lets work
C  with that in mind).
C
      Val1(1) = 1
      Val1(2) = 1
      Val1(3) = 1
      Val1(4) = 1
      Val2(1) = Nprims
      Val2(2) = Nprims 
      Val2(3) = Nprims 
      Val2(4) = Nprims 

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,*) "Atom Coordinates of each shell"
      Do Iatm =1, Ncenters
         Write(6, "(6(1x,F12.8))") (Cord(i,Iatm), i=1,3)
      Enddo

      Write(6,*)
      Write(6,*) "The atom label for each shell in main",Nshells
      Write(6,"(4(1x,I4))") (Atm_4shell(i), i=1, Nshells)
      Write(6,*)
#endif 
   
      Call compute_oed_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ioed_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ipcoef_uc),
     &                       Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, 
     &                       Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,Contract,dkh_order,
     &                       Work(Ioed_kin),
     &                       Work(Ioed_ovl),
     &                       Work(Ioed_nai),
     &                       Work(Ioed_pvp), 
     &                       Work(Icoreham),
     &                       Work(Ireord),
     &                       Erd_index,Erd_scale,
     &                       Nsend,
     &                       Nalpha,Nprims,Ncnfns,Npcoef,Npcoef_uc,
     &                       Ncenters,
     &                       charge,
     &                       max_centers,Work(Idkh),Ileft)
C
C Contracted DKH integrals are returned in Work(Ioed_pvp). First transform
C them to ACES II ordering (and normalizations). Reuse the memory locations.
C Also reset the Nbfns to the number of basis functions (it was set to
C number of primitives.

      If (Spherical) then
         Nbfns = N_Sphe_fns
      Else 
         Nbfns = N_Cart_fns
      Endif 

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Unscaled contracted DKH integrals"
      Call output(Work(Ioed_pvp), 1, Nbfns, 1, Nbfns, Nbfns,
     &            Nbfns, 1)
#endif 
       Call rescale_2ints_op1(Work(Ioed_pvp),Work(Ioed_kin),
     +                        N_Cart_fns,Erd_index,Erd_scale,
     +                        Symmetry)
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Scaled/reordered contracted DKH integrals"
      Call output(Work(Ioed_kin), 1, Nbfns, 1, Nbfns, Nbfns,
     &            Nbfns, 1)
#endif
C
C Do the symmetry adaptation. The locations Ioed_ovl,Ioed_nai,Icoreham
C is used to store the sym. transformation matrix and scratch arrays.

      Sym_adapt = .False. 
      If (Symmetry) Sym_adapt = .True.
C 
      If (Sym_adapt) Then

         Ioed_ovl      = Ioed_nai    + N_Cart_fns * N_Cart_fns
         Icoreham      = Ioed_ovl    + N_Cart_fns * N_Cart_fns
         Iend          = Icoreham    + N_Cart_fns * N_Cart_fns
         If (Iend .Ge. Mxdcor) Call insmem("Dkh:", Iend, Mxdcor)

         Call do_symadapt_cart(Nbfns,Work(Ioed_kin),
     &                         Work(Ioed_ovl),Work(Ioed_nai),
     &                         Work(Icoreham),Itriang_length)

        Call Putrec(20,"JOBARC","DKH_INTS",Itriang_length*IINTFP,
     &              Work(Icoreham))

      Else 

      Ndim = (Nbfns * (Nbfns+1))/2
#ifdef _DEBUG
       Ndim = Nbfns * Nbfns 
       Call Dzero(Work(Icoreham), Nbfns*Nbfns)
#endif 
        Call Squez2(Work(Ioed_kin), Work(Icoreham), Nbfns)
        Call Putrec(20,"JOBARC","DKH_INTS",Ndim*IINTFP,Work(Icoreham))

      Endif 
      
      Return
      End
