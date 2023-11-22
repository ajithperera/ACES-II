      Subroutine Scrnc_lrspn_cc_main(Work,Maxcor,Iuhf,Isympert,
     +                               Ioffpert,Nperts)

      Implicit Double Precision (A-H,O-Z)

#include "maxbasfn.par"
      
      Double Precision M_expect,Mbar_00
      Logical Incore
      Logical ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,Drccd
      Logical Realfreq
      Integer Occ_pair,Occ_pair_4irrep
      Dimension Isympert(Nperts),Ioffpert(Nperts) 
      Dimension Work(Maxcor)
      Character*8 Label(3),Work_label 

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
      Common/Sym_pqxx_ints/Irrep_aos(Maxbasfn),
     &                     Irrep_ao_pairs(Maxbasfn*Maxbasfn)
      Common /Extrap/Maxexp,Nreduce,Ntol,Nsizec
      Common /Eominfo/Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,Drccd
      Common /Freqinfo/Freq,Realfreq,Xfreq(200),No_freqs 
      Common /eomprop/Ieomprop

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags2.com"
#include "flags.com"
#include "parcc.com"
#include "info.com"

      Data Label /'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/
      Data Zilch /0.0D0/
     
C First step is to prepare the integrals. We have the repulsion
C integrals prepared in the form X(p(1)q(1),i(2)j(2)) where p and q can 
C be occ-occ,vrt-vrt or occ-vrt and i(2) and j(2) always kept 
C in AO basis. We process the p(1)q(1) integrals for a given pair 
C of X(2)X(2) C and lets designate them as Mu(pq).

      Irreprl     = 1
      Npert_count = 0

      Call Scrnc_form_lrspn_cc_lists(Work,Maxcor,Iuhf,Nperts)

      Do Irrep = 1, Nirrep
         Nperts_4irrep = Irpdpdao(Irrep)

         Iao_pair_4irrep = 0
         Do Iao_pair= 1, Nperts_4irrep
        
         Npert_count           = Npert_count     + 1
         Iao_pair_4irrep       = Iao_pair_4irrep + 1
         Irrep_ao_pair = Irrep_ao_pairs(Npert_count)
         Irrepr = Irrep_ao_pair
         Irrepl = Dirprd(Irrepr,Irreprl)
         Irrepx = Irrepl
         Isympert(Npert_count) = Irrepx
         Ioffpert(Npert_count) = Npert_count 

#ifdef _DEBUG
      Nbfns   = Nocco(1) + Nvrto(1)
      Call scrnc_get_pert_type(Work,Maxcor,Nbfns,Label,Isympert,
     +                        Iuhf)
      Do Occ_pair= 1, 3
      Irrepx     = Isympert(occ_pair)
      Ipert      = occ_pair
      Work_label = Label(occ_pair)
#endif 
         Nsizec = Irpdpd(Irrepx,9)
         If (Iuhf .EQ. 0) Then
            Nsizec = Nsizec  + Idsymsz(Irrepx,13,14)
         Else
            Nsizec = Nsizec + Irpdpd(Irrepx,10)
            Nsizec = Nsizec + Idsymsz(Irrepx,1,3)
            Nsizec = Nsizec + Idsymsz(Irrepx,2,4)
            Nsizec = Nsizec + Idsymsz(Irrepx,13,14)
         Endif

         Call Scrnc_reset_lrspn_cc_lists(Work,Maxcor,Iuhf,Irrepx)

C List_ijxx, List_abxx, List_aixx and List_iaxx are the OO, VV,
C VO and OV integral lists

         Lenoo = Irpdpd(Irrepx,21) + Iuhf*Irpdpd(Irrepx,22)
         Lenvv = Irpdpd(Irrepx,19) + Iuhf*Irpdpd(Irrepx,20)
         Lenvo = Irpdpd(Irrepx,9)  + Iuhf*Irpdpd(Irrepx,10)

         I000 = 1
         I010 = I000 + Lenoo
         I020 = I010 + Lenvv
         I030 = I020 + Lenvo
         I040 = I030 + Lenoo
         I050 = I040 + Lenvv
         I060 = I050 + Lenvo
         I070 = I060 + Lenvo
         Iend = I070
         Maxcor = Maxcor - Iend 

         If (Iend .Gt. Maxcor) Call Insmem("Scrnc_lrspn_main",Iend,
     +                                      Maxcor)
         Ioff_oo = I000
         Ioff_vv = I010
         Ioff_vo = I020

         Call Scrnc_prep_culomb_ints(Work(Ioff_oo),Work(Ioff_vv),
     +                               Work(Ioff_vo),Iuhf,Lenoo,
     +                                Lenvv,Lenvo,Iao_pair_4irrep,
     +                                Irrepx,Irrepr)
#ifdef _DEBUG
C      Nbfns   = Nocco(1) + Nvrto(1)
C      Call Getrec(20, "JOBARC","NBASTOT ", 1, Naobfns)
C      Length   = Naobfns * (Naobfns+1)/2
C      Length21 = Nbfns * Nbfns
C      Length22 = Nbfns * Naobfns
C      I080 = I070 + Length
C      I090 = I080 + Length21
C      I100 = I090 + Length22
C      I110 = I100 + Length22
C      Iend = I110
C      Memleft = Maxcor - Iend
C      If (Iend .Gt. Maxcor) Call Insmem("Scrnc_lrspn_main",Iend,
C     +                                      Memleft)
C      Call Scrnc_prep_dipole_ints(Work(I070),Work(I080),Work(I090),
C     +                           Work(I100),Work(I000),Work(I010),
C     +                           Work(I020),Work(Iend),
C     +                           Memleft,Lenoo,Lenvv,Lenvo,Nbfns,
C     +                           Naobfns,Iuhf,Ipert,Work_label,
C     +                           Irrepx)
C      Call Dcopy(Lenoo,Work(I000),1,Work(I030),1)
C      Call Dcopy(Lenvv,Work(I010),1,Work(I040),1)
C      Call Dcopy(Lenvo,Work(I020),1,Work(I050),1)
C      Call Dcopy(Lenvo,Work(I020),1,Work(I060),1)
#endif 
          Call Dcopy(Lenoo,Work(I000),1,Work(I030),1)
          Call Dcopy(Lenvv,Work(I010),1,Work(I040),1)
          Call Dcopy(Lenvo,Work(I020),1,Work(I050),1)
          Call Dcopy(Lenvo,Work(I020),1,Work(I060),1)

C We construct e^-TMue(pq)^T, similarity transformed Mu(p,q).

          Do ISide = 1, 1
               
             Call Scrnc_form_mubar_00(Work(Ioff_oo),Work(Ioff_vv),
     +                                Work(Ioff_vo),Work(Iend),
     +                                Maxcor,Iuhf,Irrepx,M_expect,
     +                                Mbar_00)

             Call Scrnc_form_mubar_s(Work(Ioff_oo),Work(Ioff_vv),
     +                               Work(Ioff_vo),Work(I030),
     +                               Work(I040),Work(I050),
     +                               Work(I060),Work(Iend),Maxcor, 
     +                               Iuhf,Irrepx,Iside,Mbar_00,
     +                                 Ieomprop)
             If (Iuhf .EQ. 0) Then 
                Call Scrnc_form_mubar_d_rhf(Work(I030),Work(I040),
     +                                      Work(I050),Work(Iend),
     +                                      Maxcor,Iuhf,Irrepx,
     +                                      Iside,Mbar_00,Ieomprop)
             Else
                Call Scrnc_form_mubar_d_uhf(Work(I030),Work(I040),
     +                                      Work(I050),Work(Iend),
     +                                      Maxcor,Iuhf,Irrepx,
     +                                      Iside,Mbar_00,Ieomprop)
             Endif 

             If (Iside .EQ. 1) Then
                 Mu_s   = 480
                 Mu_d   = 443
                 Ioff_s = 0
                 Ioff_d = 0
                 Call Scrnc_load_vec(Irrepx,Work(I000),Maxcor,Mu_s,
     +                               Ioff_s,Mu_d,Ioff_d,Iuhf,.False.)
                 Call Putlst(Work(I000),1,1,0,Irrepx,373)
                 Dmod_r=Ddot(Nsizec,Work(I000),1,Work(I000),1)
#ifdef _DEBUG_LVL0
      write(6,*)
      call checksum("Jbar_R<(S,D)|pq>:",Work(I000),Nsizec)
#endif
             Else
                 Mu_s   = 482
                 Mu_d   = 453
                 Ioff_s = 0
                 Ioff_d = 0
                 Call Scrnc_load_vec(Irrepx,Work(I000),Maxcor,Mu_s,
     +                               Ioff_s,Mu_d,Ioff_d,Iuhf,.False.)
                 Call Putlst(Work(I000),1,1,0,Irrepx,374)
                 Dmod_l=Ddot(Nsizec,Work(I000),1,Work(I000),1)
#ifdef _DEBUG_LVL0
      write(6,*)
      call checksum("Jbar_R<(S,D)|pq>:",Work(I000),Nsizec)
#endif
             Endif 

C If the entire perturbation is zero, then do not try to solve
C linear equations. Simply set the solution vector to zero. 

             If ((Iside .EQ. 1 .AND. Dmod_r .EQ. Zilch) .OR.
     +           (Iside .EQ. 2 .AND. Dmod_l .EQ. Zilch)) Then

                  Call Zero(Work(I000),Nsizec)
                  Call Putlst(Work(I000),Npert_count,1,1,Irrepx,
     +                   372+Iside)
             Else

             Incore = .True.
             If  (Nodavid) Then
                  Ndim    = Lenvo
                  Ndim_aa = Ndim 
                  Ndim    = Ndim + Irpdpd(Irrepx,10)
                  Kmax    = Min(Ndim,Maxexp-1)
            
                  I000 = 1
                  I010 = I000 + Kmax * Kmax
                  I020 = I010 + Kmax 
                  I030 = I020 + Kmax * Ndim
                  I040 = I030 + Kmax * Ndim
                  I050 = I040 + Ndim
                  I060 = I050 + Ndim
                  I070 = I060 + Ndim
                  I080 = I070 + 2 * Nsizec
                  Iend = I080

                  Mxcor = Maxcor - I070 
             
                  If (Iend .gt. Maxcor) Then
                      Call Insmem("Scrnc_lrspn_main",Iend,Maxcor)
                      Incore = .False.
                  Endif 
             Else 
                  Incore = .False.
             Endif 
                
             If (Incore) Then
                Ioffset  = 1
                Ipert    = Iao_pair
#ifdef _DEBUG_LVL0
      Write(6,"(1x,a)") "-----Entering Plineqy------"
#endif
                Call Plineqy(Work(I070),Mxcor,Iuhf,Work(I000),
     +                       Work(I010),Kmax,Irrepx,Ioffset,
     +                       Ipert,Freq,Iside,Realfreq,Work(I020), 
     +                       Work(I030),Work(I040),Work(I050),
     +                       Work(I060),Ndim_aa,Ndim)
             Else 

                Kmax = Maxexp-1 
                I000 = 1
                I010 = I000 + Kmax * Kmax
                I020 = I010 + Kmax
       
                Ioffset = 1
                Ipert   = Iao_pair
                Mxcor   = Maxcor - I020 
                Kmax    = Maxexp - 1
#ifdef _DEBUG_LVL0
      Write(6,"(1x,a)") "-----Entering Lineqy------"
#endif 
                Write(6,"(3a,1x,I5,I2)")" Linear equations are",
     +                           " solved from right side for" ,
     +                           " perturbation and symmetry:", 
     +                             Iao_pair,Irrepx 

                Call Lineqy(Work(I020),Mxcor,Iuhf,Work(I000),
     +                      Work(I010),Kmax,Irrepx,Ioffset,Ipert, 
     +                      Freq,Iside,Realfreq)

             Endif 

C Save the 1^st order perturbed amplitude for all perturbations

             Call Getlst(Work(I000),Ioffset,1,1,Irrepx,372+Iside)
  
             Call Putlst(Work(I000),Npert_count,1,1,Irrepx,
     +                   372+Iside)
#ifdef _DEBUG_LVL0
      write(6,*)
      call checksum("1^st order amps:",Work(I000),Nsizec)
#endif 
              Endif ! Endif for Dmod_r and Dmod_l

          Enddo ! Enddo for Iside
        Enddo ! Enddo for Iao_pair

      Enddo ! Enddo of irrepr

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,"(a)") " Symmetry of perturbations"
      Write(6,"(5(1x,I5))") (Isympert(i),i=1,Nperts)
      Write(6,"(a)") " Offset for perturbations"
      Write(6,"(5(1x,I5))") (Ioffpert(i),i=1,Nperts)
      stop
#endif
    
      Return
      End 
