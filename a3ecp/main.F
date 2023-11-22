#include "flags.h"
      Program A3ecp
      
      Implicit Double Precision (A-H, O-Z)

#include "icore.com"
#include "istart.com"
C#include "ecp.par"
C#include "ecp.com"
#include "ecplims.par"
#include "cvwecp.h"
#include "flags.com"

      Parameter (Maxbfns = 1000, Mxdcor=80000000, Mxicor= 1000000)
      Character*4 INP_Fname
      Character*4 ATMNAM(Max_centers), Namat(Max_centers)
      Logical Spherical, Cartesian 
      Logical Pvp
    
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
     &          ivAngMom(max_shells),Ixalpha(max_shells),
     &          Ixpcoef(max_shells),Indx_cc(max_shells),
     &          Ixshells(max_shells),Atm_4shell(max_shells),
     &          Temp(max_shells),CordM(3,max_centers),
     &          CordP(3,max_centers) 

cPV 
       logical delta_int, twoCent_int,FourCent_soi, angmom_int,
     &         efield_int, Dshield_int, Nattrac, fourcent_dshield,
     &         fourcent_repl,Fourcent_spnspn_cpling,pvp_ints,
     &         secnd_der_nai,Two_cent_fld_grad
cendPV
C
      Call aces_init(icore, i0,icrsiz, iuhf, .true.) 
C
      If (IFLAGS(h_IFLAGS_ecp) .EQ. 1) Iecp = 1
      Ispherical = 0
      Spherical  = .False.
      cartesian  = .False.
      Iecp       = 1
      Eps        = 0.000001D0
      If (IFLAGS(h_IFLAGS_spherical) .EQ. 1) Then
          Ispherical = 1
          Spherical  = .True.
      Else
          Cartesian  = .True.
      Endif 
C
      INP_Fname = "MOL"
      Call SIMPLE_INSPECT_MOL(INP_Fname, max_centers, max_shells,
     &                        ncenters, nshells, nspc, cartesian, 
     &                        ITFCT, LNP1, lnpo, nfct, nufct, 
     &                        nbasis, NAOBASIS, nCFpS, 
     &                        nPFpS, NAOATM, angmom, atomic_label, 
     &                        vnn, Maxang, Maxjco, Iecp,
     &                        NUcstr, Nrcstr, Iqmstr, Jcostr, 
     &                        Jstrt, Jrs, Atmnam, Charge)

C---
      Write(6,"(a,a,7(1x,I4))") "natoms,nshells,ITFCT,LNP1,lnpo,",
     &                          "nbasis,NAOBASIS:",ncenters,nshells,
     &                                         ITFCT,LNP1,lnpo,
     &                                        nbasis,NAOBASIS
      Write(6,*) 
      Write(6,"(a)") "The number of con. function per shell"
      Write(6,"(6(1x,I4))") (nCFpS(i), i=1, nshells)
      Write(6,"(a)") "The number of prim. function per shell"
      Write(6,"(6(1x,I4))") (nPFpS(i), i=1, nshells)
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

      nalpha = 0
      npcoef = 0
      do i = 1, nshells
         nalpha = nalpha + npfps(i)
         npcoef = npcoef + npfps(i) * ncfps(i)
         Ivangmom(i) = Angmom(i)
      enddo

      IALPHA = 1
      IPCOEF = IALPHA + nalpha
      INEXT  = IPCOEF + npcoef

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

      Call int_gen_init(Nshells, Ivshloff, Ivangmom, Nfps, Npfs, Ncfps,
     &                  Ixshells, Ixalpha, Ixpcoef, End_nfps, Temp, 
     &                  Ispherical)
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
      Write(6,*) "The End_Nfps array" 
      Write(6,"(8(1x,I4))") (End_Nfps(i), i=1, Nshells)
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
      
      Call aces_to_erd(Nfps, Ivangmom, Nshells, ispherical, Erd_index,
     &                 ERD_Scale)

      Call setup_ccbeg(WORK(Ialpha), IXalpha, Work(Ipcoef), Ixpcoef,
     &                 Ncfps, Npfps, Nshells, CCbeg, CCend, Indx_cc)
C
      Write(6,*)
      Write(6,*)"The CCbeg and CCend arrays" 
      Write(6,"(10(1x,I4))") (CCbeg(i), i=1, Itfct)
      Write(6,"(10(1x,I4))") (CCend(i), i=1, Itfct)
      Write(6,*)"The Erd_index and Erd_scalae arrays" 
      Write(6,"(10(1x,I4))") (Erd_index(i), i=1,Nbasis)
      Write(6,"(10(1x,F5.2))") (Erd_Scale(i), i=1,Nbasis)
      Write(6,*)
 
      Call Erd_scratch_mem_calc(nshells, ivangmom, ncfps, npfps,
     &                          atom, Coord, Work(ialpha),
     &                          Work(ipcoef), ixalpha, ixpcoef,
     &                          Ccbeg, Ccend, indx_cc,
     &                          spherical, Ncenters, .true., 
     &                          intmax, zmax)

      Call Erd__memory_hrr_correction(ivAngMom, nshells,
     &                                spherical, ihrr, zhrr)
     
      Imem   = iflags(36)
      Intmax = Intmax + Ihrr
      Zmax   = Zmax   + Zhrr

      Write(6,"(a,1x,I6,1x,I10)") "ERD memmory requirments Ibuf & Dbuf", 
     &                            Intmax, Zmax
      NBFNS = NAOBASIS
      If (Spherical) NBFNS = Nbasis 

      Ierd_dbuf          = INext
      Ierd_d4cent_Jz     = Ierd_dbuf + Zmax
      Ierd_d4cent_Jx     = Ierd_d4cent_Jz + Nbfns**4
      Ierd_d4cent_Jy     = Ierd_d4cent_Jx + Nbfns**4
      Ioed_delta_int     = Ierd_d4cent_Jy + Nbfns**4
      Ioed_d2cent_intefx = Ioed_delta_int + (Nbfns**2)*Ncenters
      Ioed_d2cent_intefy = Ioed_d2cent_intefx + Nbfns**2
      Ioed_d2cent_intefz = Ioed_d2cent_intefy + Nbfns**2
      Ioed_d2cent_intlx  = Ioed_d2cent_intefz + Nbfns**2
      Ioed_d2cent_intly  = Ioed_d2cent_intlx  + Nbfns**2
      Ioed_d2cent_intlz  = Ioed_d2cent_intly  + Nbfns**2
      Iattrac            = Ioed_d2cent_intlz  + Nbfns**2
      Ioed_d2cent_intdxx = Iattrac            + Nbfns**2
      Ioed_d2cent_intdxy = Ioed_d2cent_intdxx + Nbfns**2
      Ioed_d2cent_intdxz = Ioed_d2cent_intdxy + Nbfns**2
      Ioed_d2cent_intdyx = Ioed_d2cent_intdxz + Nbfns**2
      Ioed_d2cent_intdyy = Ioed_d2cent_intdyx + Nbfns**2
      Ioed_d2cent_intdyz = Ioed_d2cent_intdyy + Nbfns**2
      Ioed_d2cent_intdzx = Ioed_d2cent_intdyz + Nbfns**2
      Ioed_d2cent_intdzy = Ioed_d2cent_intdzx + Nbfns**2
      Ioed_d2cent_intdzz = Ioed_d2cent_intdzy + Nbfns**2
      Ioed_d4cent_intdxx = Ioed_d2cent_intdzz + Nbfns**2 
      Ioed_d4cent_intdxy = Ioed_d4cent_intdxx + Nbfns**4
      Ioed_d4cent_intdxz = Ioed_d4cent_intdxy + Nbfns**4
      Ioed_d4cent_intdyx = Ioed_d4cent_intdxz + Nbfns**4
      Ioed_d4cent_intdyy = Ioed_d4cent_intdyx + Nbfns**4
      Ioed_d4cent_intdyz = Ioed_d4cent_intdyy + Nbfns**4
      Ioed_d4cent_intdzx = Ioed_d4cent_intdyz + Nbfns**4
      Ioed_d4cent_intdzy = Ioed_d4cent_intdzx + Nbfns**4
      Ioed_d4cent_intdzz = Ioed_d4cent_intdzy + Nbfns**4
      Ioed_4eri_int      = Ioed_d4cent_intdzz + Nbfns**4
      Ioed_4repls_int    = Ioed_4eri_int      + Nbfns**4
      Iaces_ord          = Ioed_4repls_int    + Nbfns**4
      Ioed_pvp_int       = Iaces_ord          + Nbfns**4
      Ioed_hess_int      = Ioed_pvp_int       + Nbfns**2
      Ioed_aces_ord      = Ioed_hess_int      + Nbfns**2 
      Ioed_nai_x1        = Ioed_aces_ord      + Nbfns**2
      Ioed_nai_y1        = Ioed_nai_x1        + Nbfns**2
      Ioed_nai_z1        = Ioed_nai_y1        + Nbfns**2
      Ioed_nai_x2        = Ioed_nai_z1        + Nbfns**2
      Ioed_nai_y2        = Ioed_nai_x2        + Nbfns**2
      Ioed_nai_z2        = Ioed_nai_y2        + Nbfns**2
      Ioed_nai           = Ioed_nai_z2        + Nbfns**2 
      Iend               = Ioed_nai           + Nbfns**2 
#ifdef _DEBUG_LVL0     
      Write(6,"(a,1x,i5,1x,i5)") "Iend and Mxdcor:", Iend, Mxdcor
#endif

      If (Iend .Ge. Mxdcor) Call insmem("Main dcore:", Iend, Mxdcor)
      If (Intmax .Ge. Mxicor) Call insmem("Main icore:",Intmax,Mxicor)
C
C  There is one segment in serial work (not necessary, but lets work
C  with that in mind).
C
      Val1(1) = 1
      Val1(2) = 1
      Val1(3) = 1
      Val1(4) = 1
      Val2(1) = Nbfns
      Val2(2) = Nbfns
      Val2(3) = Nbfns
      Val2(4) = Nbfns
      Write(6,*)
      Write(6,*) "Atom Coordinates of each shell"
      Do Iatm =1, Ncenters
         Write(6, "(6(1x,F12.8))") (Cord(i,Iatm), i=1,3)
      Enddo

      Write(6,*)
      Write(6,*) "The atom label for each shell in main",Nshells
      Write(6,"(4(1x,I4))") (Atm_4shell(i), i=1, Nshells)
      Write(6,*)
        
      Nattrac=.False. 
      delta_int=.false. 
      twoCent_int=.false.
      FourCent_soi=.false.
      FourCent_Repl=.true.
      angmom_int=.false.
      efield_int=.false.
      dshield_int=.false.
      Fourcent_dshield=.false.
      Fourcent_spnspn_cpling=.false. 
      pvp_ints=.false. 
      secnd_der_nai=.false.
      two_cent_fld_grad=.false.

      if(Nattrac) then

      Call compute_kin_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Iattrac),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,charge,
     &                       max_centers)

      Call compute_nattrac_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Iattrac),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,charge,
     &                       max_centers)

       end if

      if(FourCent_soi) then
      Call compute_4cent_integrals_SOI(val1(1),val2(1), val1(2),val2(2),
     &                       val1(3),val2(3),val1(4),val2(4),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef, 
     &                       Ispherical,
     &                       Work(Ierd_d4cent_Jz), 
     *                       work(Ierd_d4cent_Jx),
     *                       work(Ierd_d4cent_Jy),Nsend, 
     &                       Nalpha, Npcoef,ERD_index,ERD_scale,
     &                       Work(Iaces_ord))
      end if

c      Write(*,*) 'Time to calculate 1electron SOI'
c      write(*,*) 'delta_int=',delta_int
c      write(*,*) 'twoCent_int=',twoCent_int

      if(delta_int) then
      Call compute_delta_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &   Work(Ioed_delta_int),
     *   Nsend, Nalpha, Npcoef, Ncenters,  max_centers)
      end if

       if(angmom_int) then
      Call compute_angmom_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_d2cent_intlx), 
     &                       Work(Ioed_d2cent_intly), 
     &                       Work(Ioed_d2cent_intlz),
     &                       Nsend, 
     &                       Nalpha, Npcoef, Ncenters,
     &                       max_centers)
       end if

      if(pvp_ints) then
      Call compute_pvp_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_pvp_int),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,Charge, 
     &                       max_centers,ERD_index, 
     &                       ERD_scale,Work(Iaces_ord))
       end if

      if(secnd_der_nai) then

#ifdef _NOSKIP
      Do Iatm = 1, 2
      Call Dcopy(3*Ncenters,Cord,1,CordM,1)
      Do Ixyz  = 3, 3
         CordM(Ixyz,Iatm) = Cord(Ixyz,Iatm) - EPs
      Do Jatm = 1, Iatm 
      Call Dcopy(3*Ncenters,Cord,1,CordP,1)
      Do Jxyz  = 3, 3
         CordP(Ixyz,Jatm) = Cord(Ixyz,Jatm) + EPs

         Write(6,"(a,4(1x,i3))") "Iatm,jatm,Ixyz,Jxyz: ",Iatm,Jatm,
     &                            Ixyz,Jxyz
      Write(6,"(a)") "-------------------------------------"        

      Call compute_1dernai_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       CordM,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_nai_x1),Work(Ioed_nai_y1),
     &                       Work(Ioed_nai_z1),Work(Ioed_nai),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,Charge,
     &                       max_centers,ERD_index,
     &                       ERD_scale)

      Call compute_1dernai_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       CordP,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_nai_x2),Work(Ioed_nai_y2),
     &                       Work(Ioed_nai_z2),Work(Ioed_nai),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,Charge,
     &                       max_centers,ERD_index,
     &                       ERD_scale)

      Call Perform_numder(Work(Ioed_nai_x1),Work(Ioed_nai_y1),
     &                    Work(Ioed_nai_z1),Work(Ioed_nai_x2),
     &                    Work(Ioed_nai_y2),Work(Ioed_nai_z2),
     &                    Work(Ioed_nai),Nbfns,Eps,Ixyz,Jxyz)
      Enddo
      Enddo
      Enddo
      Enddo
#endif 
C#ifdef _NOSKIP
      Call compute_2dernai_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_pvp_int),Work(Ioed_hess_int),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,Charge,
     &                       max_centers,ERD_index,
     &                       ERD_scale)
C#endif 
      End if
C#ifdef _4CENTDSH
      if(efield_int) then
      Call compute_efield_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_d2cent_intefx), 
     &                       Work(Ioed_d2cent_intefy), 
     &                       Work(Ioed_d2cent_intefz),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,charge,
     &                       max_centers)

      Call compute_efield_integrals_test(val1(1),val2(1),
     &                       val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_d2cent_intefx),
     &                       Work(Ioed_d2cent_intefy),
     &                       Work(Ioed_d2cent_intefz),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,charge,
     &                       max_centers)


       end if
C#endif 

      if(dshield_int) then
      Call compute_dshield_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_d2cent_intdxx), 
     &                       Work(Ioed_d2cent_intdxy), 
     &                       Work(Ioed_d2cent_intdxz), 
     &                       Work(Ioed_d2cent_intdyx), 
     &                       Work(Ioed_d2cent_intdyy), 
     &                       Work(Ioed_d2cent_intdyz), 
     &                       Work(Ioed_d2cent_intdzx), 
     &                       Work(Ioed_d2cent_intdzy), 
     &                       Work(Ioed_d2cent_intdzz), 
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,charge,
     &                       max_centers)
       end if

      if(Fourcent_dshield) then
      Call compute_4cent_dshield_integrals(val1(1),val2(1),
     &                       val1(2),val2(2),
     &                       val1(3),val2(3),
     &                       val1(4),val2(4),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_d4cent_intdxx),
     &                       Work(Ioed_d4cent_intdxy),
     &                       Work(Ioed_d4cent_intdxz),
     &                       Work(Ioed_d4cent_intdyx),
     &                       Work(Ioed_d4cent_intdyy),
     &                       Work(Ioed_d4cent_intdyz),
     &                       Work(Ioed_d4cent_intdzx),
     &                       Work(Ioed_d4cent_intdzy),
     &                       Work(Ioed_d4cent_intdzz),
     &                       Work(Ioed_4eri_int), Nsend,
     &                       Nalpha, Npcoef, Ncenters,charge,
     &                       max_centers,ERD_index,ERD_scale,
     &                       Work(Iaces_ord))
       end if

      if(Fourcent_spnspn_cpling) then
      Call compute_4cent_spnspn_cpling_integrals(val1(1),val2(1),
     &                       val1(2),val2(2),
     &                       val1(3),val2(3),
     &                       val1(4),val2(4),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_d4cent_intdxx),
     &                       Work(Ioed_d4cent_intdxy),
     &                       Work(Ioed_d4cent_intdxz),
     &                       Work(Ioed_d4cent_intdyx),
     &                       Work(Ioed_d4cent_intdyy),
     &                       Work(Ioed_d4cent_intdyz),
     &                       Work(Ioed_d4cent_intdzx),
     &                       Work(Ioed_d4cent_intdzy),
     &                       Work(Ioed_d4cent_intdzz),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,charge,
     &                       max_centers,ERD_index,ERD_scale,
     &                       Work(Iaces_ord))
       end if

       if(fourCent_Repl) then
      Call compute_4cent_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       val1(3),val2(3),
     &                       val1(4),val2(4),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &   Work(Ioed_4repls_int),
     &   Nalpha, Npcoef)
       end if

       if(twoCent_int) then
      Call compute_2cent_integrals_SOI(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &   Work(Ioed_d2cent_intefx), Work(Ioed_d2cent_intefy),
     &   Work(Ioed_d2cent_intefz),Nsend,
     &   Nalpha, Npcoef, Ncenters,charge,max_centers)
       end if

      if(two_cent_fld_grad) then
      Call compute_2cent_fld_grad(val1(1),val2(1),
     &                       val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_d2cent_intdxx),
     &                       Work(Ioed_d2cent_intdxy),
     &                       Work(Ioed_d2cent_intdxz),
     &                       Work(Ioed_d2cent_intdyx),
     &                       Work(Ioed_d2cent_intdyy),
     &                       Work(Ioed_d2cent_intdyz),
     &                       Work(Ioed_d2cent_intdzx),
     &                       Work(Ioed_d2cent_intdzy),
     &                       Work(Ioed_d2cent_intdzz),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,charge,
     &                       max_centers,ERD_index,ERD_scale,
     &                       Work(Ioed_aces_ord))
       end if
C
      Call aces_fin
      
      Stop
      End
