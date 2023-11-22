#include "flags.h"
      Subroutine Tdee_form_nondiag_state(Work,Maxmem,Nondiag_component,
     +                                   Iuhf)

      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Maxmem)
      Double Precision Mubar_00, M_expect

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "flags.com"
#include "flags2.com"
#include "info.com"
#include "tdee_intgrt_vars.com"
#include "lanczos_vars.com"

      Integer Isympert(10)
      Character*8 Label(6)
      Character*8 Work_label 
      Character*1 DXYZ(3)
      Character*2 QXYZ(6)
      Logical Mu_bar   
      Logical Mubar_t
      Logical Source
      Logical Target
      Logical Sing
      Common /METH/MBPT2_DUMMY,MBPT3,M4DQ,M4SDQ,M4SDTQ,
     +              CCD_DUMMY,QCISD,CCSD,UCC
      Common /REFTYPE/MBPT2,CC,CCD,RPA,DRPA,LCCD,LCCSD,CC2
      Common /DRVHBAR/ SS, SD, DS, DD
      COMMON/LISTDENS/LDENS

      Data Ione, One /1, 1.0D0/
      Data DXYZ /"x", "y", "z"/
      Data QXYZ /"xx", "yy", "zz","xy","xz","yz"/

      Sing   = .False.
      Icontl = Iflags(h_IFLAGS_cc_conv)
      Iconvg = 1
      Ncycle = 0
      Sing   = (Iflags(h_IFLAGS_calc) .gt. 9)

      Mu_bar  = .True.
      Mubar_t = .True.
      Source  = .True.
      Target  = .True.

      Nbfns   = Nocco(1) + Nvrto(1)
      Call Getrec(20, "JOBARC","NBASTOT ", Ione, Naobfns)

      Call Tdee_get_pert_type(Work,Maxmem,Nbfns,Label,Isympert,
     +                        Iuhf)

      Do Ixyz = Nondiag_component,Nondiag_component

         Irrepx     = Isympert(Ixyz)
         Ipert      = Ixyz
         Work_label = Label(Ixyz)

#ifdef _DEBUG_LVL0
      Write(6,"(a,i2)") " The non-diagonal component: ",Ixyz
#endif
         Call Tdee_reset_lists(Work,Maxmem,Iuhf,Irrepx)

         Call Tdee_init_lists(Work,Maxmem,Iuhf,Irrepx,Source,
     +                        Target,Side)

         Length   = Naobfns * (Naobfns+1)/2
         Length21 = Nbfns * Nbfns
         Length22 = Nbfns * Naobfns
         Lenoo    = Irpdpd(Irrepx,21) + Iuhf * Irpdpd(Irrepx,22)
         Lenvv    = Irpdpd(Irrepx,19) + Iuhf * Irpdpd(Irrepx,20)
         Lenvo    = Irpdpd(Irrepx,9)  + Iuhf * Irpdpd(Irrepx,10)
         Lenvo2_aa= Irpdpd(Irrepx,9)  * Irpdpd(Irrepx,9)
         Lenvo2_bb= Irpdpd(Irrepx,10) * Irpdpd(Irrepx,10)
         Lenvo2_ab= Irpdpd(Irrepx,11) * Irpdpd(Irrepx,11)

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a,a)") "Length,Length21,Length22,Lenoo,Lenvv,Lenv0",
     +                   ",Lenvo2_aa,Lenvo2_bb,Lenvo2_ab: "
        Write(6,"(9(1x,I5))")Length,Length21,Length22,Lenoo,Lenvv,
     +                       Lenvo,Lenvo2_aa, Lenvo2_bb, Lenvo2_ab
        Write(6,*)
#endif

         Nsize = Irpdpd(Irrepx,9)
         If (Iuhf .EQ. 0) Then
            Nsize = Nsize  + Idsymsz(Irrepx,13,14)
         Else
           Nsize = Nsize + Irpdpd(Irrepx,10)
           Nsize = Nsize + Idsymsz(Irrepx,1,3)
           Nsize = Nsize + Idsymsz(Irrepx,2,4)
           Nsize = Nsize + Idsymsz(Irrepx,13,14)
         Endif

         Ibgn = Ione
         I000 = Ibgn
         I010 = I000 + Lenoo
         I020 = I010 + Lenvv
         I030 = I020 + Lenvo
         I040 = I030 + Length
         I050 = I040 + Length21
         I060 = I050 + Length22
         I070 = I060 + Length22
         Iend = I070
         Memleft = Maxmem - Iend

         If (Iend .Gt. Maxmem) Call Insmem("@-Tdee_driver",Iend,
     +                                      Memleft)

         If (D_Pole) Then
         Call Tdee_prep_dipole_ints(Work(I030),Work(I040),Work(I050),
     +                              Work(I060),Work(I000),Work(I010),
     +                              Work(I020),Work(Iend),
     +                              Memleft,Lenoo,Lenvv,Lenvo,Nbfns,
     +                              Naobfns,Iuhf,Ipert,Work_label,
     +                              Irrepx)
         Elseif (Q_pole) Then
         Call Tdee_prep_qdpole_ints(Work(I030),Work(I040),Work(I050),
     +                              Work(I060),Work(I000),Work(I010),
     +                              Work(I020),Work(Iend),
     +                              Memleft,Lenoo,Lenvv,Lenvo,Nbfns,
     +                              Naobfns,Iuhf,Ipert,Work_label,
     +                              Irrepx)
         Endif

C Built the  Mu_bar = exp(-T)Mexp(T)

         If (Mu_bar) Then

            Ibgn   = I030
            I040   = I030 + Lenoo
            I050   = I040 + Lenvv
            I060   = I050 + Lenvo
            I070   = I060 + Lenvo
            Iend   = I070
            Iside  = 1

            Memleft = Maxmem - Iend
            Call Dcopy(Lenoo,Work(I000),1,Work(I030),1)
            Call Dcopy(Lenvv,Work(I010),1,Work(I040),1)
            Call Dcopy(Lenvo,Work(I020),1,Work(I050),1)
            Call Dcopy(Lenvo,Work(I020),1,Work(I060),1)

            Call Tdee_form_mubar_00(Work(I000),Work(I010),Work(I020),
     +                             Work(Iend),Memleft,Iuhf,Irrepx,
     +                             M_expect,Mubar_00)

            Call Tdee_form_mubar_s(Work(I000),Work(I010),Work(I020),
     +                             Work(I030),Work(I040),Work(I050),
     +                             Work(I060),Work(Iend),Memleft,Iuhf,
     +                             Irrepx,Iside,Mubar_00)

            If (Iuhf .Eq. 0) Then

               Call Tdee_form_mubar_d_rhf(Work(I030),Work(I040),
     +                                    Work(I050),Work(Iend),
     +                                    Memleft,Iuhf,Irrepx,Iside,
     +                                    Mubar_00)
            Else


               Call Tdee_form_mubar_d_uhf(Work(I030),Work(I040),
     +                                    Work(I050),Work(Iend),
     +                                    Memleft,Iuhf,Irrepx,Iside,
     +                                    Mubar_00)
            Endif

         Endif

         If (Mubar_t) Then

            Ibgn   = I070
            I080   = I070 + Lenoo
            I090   = I080 + Lenvv
            I100   = I090 + Lenvo
            I110   = I100 + Lenvo
            Iend   = I110
            Iside  = 2

            Memleft = Maxmem - Iend
            Call Dcopy(Lenoo,Work(I000),1,Work(I070),1)
            Call Dcopy(Lenvv,Work(I010),1,Work(I080),1)
            Call Dcopy(Lenvo,Work(I020),1,Work(I090),1)
            Call Dcopy(Lenvo,Work(I020),1,Work(I100),1)

            Call Tdee_form_mubar_s(Work(I030),Work(I040),Work(I050),
     +                             Work(I070),Work(I080),Work(I090),
     +                             Work(I100),Work(Iend),Memleft,Iuhf,
     +                             Irrepx,Iside,Mubar_00-M_expect)

            If (Iuhf .Eq. 0) Then

               Call Tdee_form_mubar_d_rhf(Work(I030),Work(I040),
     +                                    Work(I060),Work(Iend),
     +                                    Memleft,Iuhf,Irrepx,Iside,
     +                                    Mubar_00-M_expect)
            Else

               Call Tdee_form_mubar_d_uhf(Work(I030),Work(I040),
     +                                    Work(I060),Work(Iend),
     +                                    Memleft,Iuhf,Irrepx,Iside,
     +                                    Mubar_00-M_expect)
            Endif

         Endif

         Memleft = Maxmem

C Save the M_bar_tilde(0) and M_bar(0) to built the autocorrelation
C function.
C Note that this is a safe measure on since any of the subsequent
C steps do not touch the list containing M_bar_tilde(0) and M_bar(0).

         Do Iside = 1, 2
            Call Tdee_save_initial_state(Work(I000),Memleft,Irrepx,
     +                                   Iuhf,Iside,Nsize)
         Enddo

      Enddo
      
      Return
      End

