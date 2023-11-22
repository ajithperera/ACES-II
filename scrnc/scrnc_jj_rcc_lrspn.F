      Subroutine Scrnc_jj_rcc_lrspn(Work,Maxcor,Iuhf,JJ_lrespn,Isympert,
     +                             Ioffpert,Nao_pair)

      Implicit Double Precision (A-H,O-Z)
      
      Double Precision M_expect,Mbar_00
      Double Precision JJ_lrespn(Nao_pair,Nao_pair)
      Logical Incore
      Logical ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,Drccd
      Logical Realfreq
      Integer Occ_pair,Occ_pair_4irrep
      Dimension Isympert(Nao_pair),Ioffpert(Nao_pair)
      Dimension Work(Maxcor)
      Character*8 Label(3),Work_label

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)

      Common /Extrap/Maxexp,Nreduce,Ntol,Nsizec
      Common /Eominfo/Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,
     +       Drccd
      Common /Freqinfo/Freq,Realfreq,Xfreq(200),No_freqs 
      Common /eomprop/Ieomprop

      Data Label /'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags2.com"
#include "flags.com"
#include "parcc.com"
#include "info.com"

      Irreprl     = 1
      Npert_count = 0

#ifdef _DEBUG_LVL0
      Write(6,"(1x,a)") "-----Entering scrnc_jj_rcc_lrspn------"
#endif
      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irreprl)
         Irrepx = Irrepl

         Nsizec = Irpdpd(Irrepx,9)
         If (Iuhf .NE. 0) Then
            Nsizec = Nsizec + Irpdpd(Irrepx,10)
         Endif

         Call Scrnc_reset_lrspn_rcc_lists(Work,Maxcor,Iuhf,Irrepx)

C List_ljki, List_abki and List_bjki are the OO, VV and OV integral
C list
         Lenoo = Irpdpd(Irrepx,21) + Iuhf*Irpdpd(Irrepx,22)
         Lenvv = Irpdpd(Irrepx,19) + Iuhf*Irpdpd(Irrepx,20)
         Lenvo = Irpdpd(Irrepx,9)  + Iuhf*Irpdpd(Irrepx,10)

         Nperts_4irrepr = Irpdpdao(Irrepr)
 
         Iao_pair_4irrep = 0
         Do Iao_pair= 1, Nperts_4irrepr

         I000 = 1
         I010 = I000 + Lenoo
         I020 = I010 + Lenvv
         I030 = I020 + Lenvo
         I040 = I030 + Lenoo
         I050 = I040 + Lenvv
         I060 = I050 + Lenvo
         I070 = I060 + Lenvo
         Iend = I070
         Mxcor = Maxcor - Iend 

         If (Iend .Gt. Maxcor) Call Insmem("Scrnc_jj_lrspn",Iend,
     +                                      Maxcor)
         Ioff_oo = I000
         Ioff_vv = I010 
         Ioff_vo = I020 

         Iao_pair_4irrep = Iao_pair_4irrep + 1 
         Npert_count     = Npert_count     + 1

         Call Scrnc_prep_culomb_ints(Work(Ioff_oo),Work(Ioff_vv),
     +                               Work(Ioff_vo),Iuhf,Lenoo,
     +                               Lenvv,Lenvo,Iao_pair_4irrep,
     +                               Irrepx,Irrepr)
         Call Dcopy(Lenoo,Work(I000),1,Work(I030),1)
         Call Dcopy(Lenvv,Work(I010),1,Work(I040),1)
         Call Dcopy(Lenvo,Work(I020),1,Work(I050),1)
         Call Dcopy(Lenvo,Work(I020),1,Work(I060),1)

C We construct e^-TMue(pq)^T, similarity transformed Mu(p,q).
            
         Iside = 1
         Call Scrnc_form_mubar_00(Work(Ioff_oo),Work(Ioff_vv),
     +                            Work(Ioff_vo),Work(Iend),
     +                            Mxcor,Iuhf,Irrepx,M_expect,
     +                            Mbar_00)

         Call Scrnc_form_mubar_s(Work(Ioff_oo),Work(Ioff_vv),
     +                           Work(Ioff_vo),Work(I030),
     +                           Work(I040),Work(I050),
     +                           Work(I060),Work(Iend),Mxcor, 
     +                           Iuhf,Irrepx,Iside,Mbar_00,
     +                           Ieomprop)

         Iside = 2
         Ibgn  = I070 
         I080  = I070 + Lenoo
         I090  = I080 + Lenvv
         I100  = I090 + Lenvo
         I110  = I100 + Lenvo
         Iend  = I110
         Mxcor = Maxcor - Iend 
         If (Iend .Gt. Maxcor) Call Insmem("Scrnc_lrspn_main",Iend,
     +                                     Maxcor)

         Call Dcopy(Lenoo,Work(I000),1,Work(I070),1)
         Call Dcopy(Lenvv,Work(I010),1,Work(I080),1)
         Call Dcopy(Lenvo,Work(I020),1,Work(I090),1)
         Call Dcopy(Lenvo,Work(I020),1,Work(I100),1)
 
         Call Scrnc_form_mubar_s(Work(I030),Work(I040),
     +                           Work(I050),Work(I070),
     +                           Work(I080),Work(I090),
     +                           Work(I100),Work(Iend),Mxcor,   
     +                           Iuhf,Irrepx,Iside,Mbar_00-
     +                           M_expect,Ieomprop)

         If (Iside .Eq. 2) Then 
             Mu_s    = 482
             Ioffset = 1
             Lenr1a = Irpdpd(Irrepx,9)
             Lenr1b = Irpdpd(Irrepx,10)

             I000 = 1
             I020 = I010 + Lenr1a
             Call Getlst(Work(I000),1,1,1,1,Mu_s)
             If (Iuhf .NE. 0) Then
                 I020 = I010 + Lenr1b
                 Call Getlst(Work(I010),1,1,1,2,Mu_s)
             Endif

             Call Putlst(Work(I010),1,1,0,Irrepx,374)

#ifdef _DEBUG_LVL0
      write(6,*)
      call checksum("Jbar_L<(S,D)|pq>:",Work(I010),Nsizec)
#endif
          Endif 

C Save the left amplitudeds.

          Call Getlst(Work(I000),Ioffset,1,1,Irrepx,372+Iside)
          Call Putlst(Work(I000),Npert_count,1,1,Irrepx,372+Iside)

         Enddo ! Enddo for Iao_pair 
      Enddo ! Enddo of irrepr

      Call Dzero(JJ_lrespn,Nao_pair*Nao_pair)

      I_pairs = Nao_pair
      J_pairs = Nao_pair

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " Symmetry of perturbations"
      Write(6,"(5(1x,I5))") (Isympert(i),i=1,Nao_pair)
      Write(6,"(a)") " Offset for perturbations"
      Write(6,"(5(1x,I5))") (Ioffpert(i),i=1,Nao_pair)
      Write(6,*)
#endif 
      Do Iao_pair = 1, I_pairs 

         Irrepx = Isympert(Iao_pair)

         Nsizec = Irpdpd(Irrepx,9)
         If (Iuhf .NE. 0) Then
            Nsizec = Nsizec + Irpdpd(Irrepx,10)
         Endif

         I000 = 1
         I010 = I000 + Nsizec 
         I020 = I010 + Nsizec 
         Iend = I020

         Memleft = Maxcor - Iend 
         If (Iend .Gt. Maxcor)  Call Insmem("Scrnc_jj_lrspn",Iend,
     +                                      Memleft)

         Call Getlst(Work(I000),Iao_pair,1,1,Irrepx,373)
#ifdef _DEBUG_LVL0
         Write(6,"(a)") " <0|(1+L)M_bar and first-order e^T|0>"
         call checksum("373",work(i000),nsizec)
#endif 
         If (Iuhf .eq. 0) Then
            Call Spntsing(Isympert(Iao_pair),Work(I000),Work(Iend),
     +                    Memleft) 
         Endif 

         Do Jao_pair = 1, J_pairs 

            If (Isympert(Iao_pair) .eq. Isympert(Jao_pair)) Then
 
               Call Getlst(Work(I010),Jao_pair,1,1,Irrepx,
     +                     374)
#ifdef _DEBUG_LVL0
               call checksum("374",work(i010),nsizec)
#endif 
               X = Ddot(Nsizec,Work(I010),1,Work(I000),1)

               JJ_lrespn(Iao_pair,Jao_pair) = 
     +                  JJ_lrespn(Iao_pair,Jao_pair) + X
            Endif 
              
         Enddo
      Enddo 
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " RPA JJ matrix"
      call output(JJ_lrespn,1,Nao_pair,1,Nao_pair,Nao_pair,
     +            Nao_pair,1)
      Write(6,*)
#endif 
      
      Return
      End 
