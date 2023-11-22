      Subroutine Form_notrans_pt2(Work,Maxcor,Iwork,Imaxcor,Iuhf,
     +                            Irrepx,Nsizec,Ioffset,Npert,Iside)

      Implicit Double Precision(A-H,O-Z)
      Integer T2ln_aa,T2ln_bb,T2ln_ab,T2ln
      Integer T1ln,T1ln_aa,T1ln_bb
      Integer T2ln_aabb_ab
      Integer T2ln_aabb
      Logical Cct2,Cct1
      Logical Build_dens
      Logical Save_vecs
      Character*4 Type

      Dimension Work(Maxcor)
      Dimension Iwork(Imaxcor)
      Dimension Lenvv(2)
      Dimension Lenoo(2)
      Dimension Npert(8)

#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "flags.com"
#include "flags2.com"
#include "sympop.com"
#include "symoff.com"
#include "info.com"
#include "active_space.com"   

      Data Ione /1/

      T2ln_aa    = 0
      T2ln_bb    = 0
      T2ln_ab    = 0
      T1ln_aa    = 0
      T1ln_bb    = 0
      List2_feom = 372 + Iside
      List2_peom = 372 + Iside

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " No. Alpha OO and VV perturbed active:"
      Write(6,"(8(1x,i3))") (Pactive_oo(Ioffset,i,1),i=1,Nirrep)
      Write(6,"(8(1x,i3))") (Pactive_vv(Ioffset,i,1),i=1,Nirrep)
      Write(6,*)
      Write(6,"(a)") " No. Beta OO and VV perturbed active:"
      Write(6,"(8(1x,i3))") (Pactive_oo(Ioffset,i,2),i=1,Nirrep)
      Write(6,"(8(1x,i3))") (Pactive_vv(Ioffset,i,2),i=1,Nirrep)
#endif 

C Precompute the symmetry packed lengths of <AB|IJ> like distributions

      Do Ispin = 1, (Iuhf+1)
         If (Ispin .Eq. 1) Then
            Do Irrepr = 1, Nirrep
               Irrepl = Dirprd(Irrepr,Irrepx)
               Ncol = Irpdpd(Irrepr,2+Ispin)
               Nrow = Irpdpd(Irrepl,Ispin)
               T2ln_aa  = T2ln_aa + Nrow*Ncol
            Enddo
         Elseif (Ispin .Eq. 2) Then
            Do Irrepr = 1, Nirrep
               Irrepl = Dirprd(Irrepr,Irrepx)
               Ncol = Irpdpd(Irrepr,2+Ispin)
               Nrow = Irpdpd(Irrepl,Ispin)
               T2ln_bb  = T2ln_bb + Nrow*Ncol
            Enddo
         Endif
      Enddo

      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)
         Ncol = Irpdpd(Irrepr,14)
         Nrow = Irpdpd(Irrepl,13)
         T2ln_ab = T2ln_ab + Nrow*Ncol
      Enddo

      T1ln_aa = Irpdpd(Irrepx,9)
      If (Iuhf .Ne. 0) T1ln_bb = Irpdpd(Irrepx,10)

C Precompute the irrep offsets occ-occ and vrt-vrt blocks

      Ioff_oo(1,1) = Ione
      Ioff_oo(1,2) = Nfmi(1) + Ione
      Ioff_vv(1,1) = Ione
      Ioff_vv(1,2) = Nfea(1) + Ione
      Do Ispin = 1, Iuhf+1
         Do irrep = 1, Nirrep-1
            Ioff_oo(Irrep+1,Ispin) = Ioff_oo(Irrep,Ispin) +
     +                               Pop(Irrep,Ispin)*Pop(Irrep,Ispin)
            Ioff_vv(Irrep+1,Ispin) = Ioff_vv(Irrep,Ispin) +
     +                               Vrt(Irrep,Ispin)*Vrt(Irrep,Ispin)
        Enddo
      Enddo

      T2ln = T2ln_aa + Iuhf*T2ln_bb + T2ln_ab
      T1ln = T1ln_aa + Iuhf*T1ln_bb

C These quantities are totally symmetric.

      Do Ispin = 1, Iuhf+1
         Lenvv(Ispin) = Irpdpd(1,18+Ispin)
         Lenoo(Ispin) = Irpdpd(1,20+Ispin)
      Enddo

      If (T1ln + T2ln .Ne. Nsizec) Then 
         Write(6,"(2a)") " There is an inconsistency in the length of",
     +                   " perturbed T vector in"
         Write(6,"(a)")  " form_notrans_pt2." 
         Call Errex
      Endif 

      I000 = Ione
      I020 = I000 + Nsizec
      I030 = I020 + Lenoo(1)+Iuhf*Lenoo(2)
      Iend = I030 + Lenvv(1)+Iuhf*Lenvv(2)

      If (Iend .Gt. Maxcor) Call Insmem("form_notrans_pt2",Iend,
     +                                   Maxcor)
      Maxcor = Maxcor - Iend

      Call Getlst(Work(I020),Ioffset,1,1,Irrepx,380)
      Call Getlst(Work(I030),Ioffset,1,1,Irrepx,381)

      If (Iuhf .Ne. 0) Then
         Call Getlst(Work(I020+Lenoo(1)),Ioffset,1,1,Irrepx,382)
         Call Getlst(Work(I030+Lenvv(1)),Ioffset,1,1,Irrepx,383)
      Endif

#ifdef _DEBUG_LVLM
      Write(6,"(a)") " Check for unitary: O(t)O and V(t)V"
      I040 = Iend
      I050 = I040 + (Nfmi(1)+Iuhf*Nfmi(2))
      Iend = I050 + (Nfea(1)+Iuhf*Nfea(2))
      Ioff = Ione
      Joff = Ione
      do irrep =1,Nirrep
      Noo  = Pop(irrep,1)
      Nvv  = vrt(irrep,1)
      call Xgemm("T","N",Noo,Noo,Noo,1.0D0,Work(I020+ioff-1),Noo,
     +            Work(I020+ioff-1),Noo,0.0D0,Work(i040+ioff-1),Noo)
      call output(Work(i040+Ioff-1),1,noo,1,noo,noo,noo,1)
      Ioff = Ioff + noo*noo
      call Xgemm("T","N",Nvv,Nvv,Nvv,1.0D0,Work(I030+joff-1),Nvv,
     +            Work(I030+joff-1),Nvv,0.0D0,Work(i050+joff-1),Nvv)
      call output(Work(i050+joff-1),1,nvv,1,nvv,nvv,nvv,1)
      Joff = Joff + nvv*nvv
      enddo
#endif

C Transform the EOM-CCSD perturbed T amplitudes to PNO basis.
C These are now stored in 490 (T1x),444-446(T2,x)

      Call Getlst(Work(I000),Ioffset,1,1,Irrepx,List2_feom)

      List1 = 490
      List2 = 444

      Call Putlst(Work(I000),1,1,1,1,List1)
#ifdef _DEBUG_LVL0
      call checksum("T1aa :",Work(I000),T1ln_aa)
#endif 
      Ioff = T1ln_aa
      If (Iuhf .Ne. 0) Call Putlst(Work(I000+Ioff),1,1,1,2,List1)
#ifdef _DEBUG_LVL0
      call checksum("T1bb :",Work(I000+Ioff),T1ln_bb)
#endif 
      Ioff = T1ln_aa + Iuhf*T1ln_bb
      Call Putall(Work(I000+Ioff),T2ln_ab,Irrepx,List2+2)
#ifdef _DEBUG_LVL0
      call checksum("T2ab :",Work(I000+Ioff),T2ln_ab)
#endif 
      Ioff = T1ln_aa + Iuhf*T1ln_bb + T2ln_ab
      Call Putall(Work(I000+Ioff),T2ln_aa,Irrepx,List2)
      call checksum("T2aa :",Work(I000+Ioff),T2ln_aa)

      Ioff =  T1ln_aa + Iuhf*T1ln_bb + T2ln_ab + T2ln_aa
      If (Iuhf .Ne. 0) Call Putall(Work(I000+Ioff),T2ln_bb,
     +                             Irrepx,List2+1)
#ifdef _DEBUG_LVL0
      call checksum("T2bb :",Work(I000+Ioff),T2ln_bb)
      call checksum("Tx   :",Work(I000),Nsizec)
      Write(6,"(a)") "C2N"
#endif 
      Type  = "C2N"
      CCt2  = .True.
      CCt1  = .True.

      List1_h = 490
      List1_t = 493
      Call Transform_2index_pt(Work(I000),Work(I020),Work(I030),
     +                         Work(Iend),Maxcor,List1_h,List1_t,
     +                         T1ln,Iuhf,Irrepx,CCt1,Type)

      List2_h = 443 
      List2_t = 453
      Call Transform_4index_aa_pt(Work(I000),Work(I020),Work(I030),
     +                            Work(Iend),Maxcor,List2_h,List2_t,
     +                            T2ln,T2ln_aa,Iuhf,Irrepx,Cct2,Type)

      List2_h = 446 
      List2_t = 456
      Call Transform_4index_ab_pt(Work(I000),Work(I020),Work(I030),
     +                            Work(Iend),Maxcor,List2_h,List2_t,
     +                            T2ln,T2ln_aa,T2ln_bb,Iuhf,Irrepx,
     +                            Cct2,Type)

C Solve the PEOM-CCSD equations for the current perturbation. The NO 
C transformed EOM-CCSD vector for the current perturbation is at 
C Work(I000).

      Numpert    = 1
      Build_dens = .False.
      Save_vecs  = .False.

      Call Get_peom_nos(Work(Iend),Maxcor,Iwork,Imaxcor,Iuhf,Npert,
     +                  Ioffset,Numpert,Irrepx,Iside,Build_dens,
     +                  Save_vecs)

      Call Reorder_tx(Work(I000),Work(Iend),Maxcor,T1ln,T2ln,
     +                T1ln_aa,T1ln_bb,T2ln_aa,T2ln_bb,T2ln_ab,
     +                Nsizec,Iuhf)

      CCt2  = .False.
      CCt1  = .False.
      List1 = 0
      List2 = 0

      Call Transform_2index_pt(Work(I000),Work(I020),Work(I030),
     +                         Work(Iend),Maxcor,List1,List1,T1ln,
     +                         Iuhf,Irrepx,CCt1,Type)

      Ioff = T1ln_aa + Iuhf*T1ln_bb

      Call TranSform_4index_aa_pt(Work(I000+Ioff),Work(I020),
     +                            Work(I030),Work(Iend),Maxcor,
     +                            List2,List2,T2ln,T2ln_aa,Iuhf,
     +                            Irrepx,Cct2,Type)

      Call Transform_4index_ab_pt(Work(I000+Ioff),Work(I020),
     +                            Work(I030),Work(Iend),Maxcor,
     +                            List2,List2,T2ln,T2ln_aa,T2ln_bb,
     +                            Iuhf,Irrepx,Cct2,Type)

C Now as required by the next routine the PEOM-CCSD Tx are at Work(I000)
C and EOM-CCSD vectors are read from 190 and 244-246 lists.
    
      Call Form_new_pamps(Work(I000),Work(Iend),Maxcor,Iwork,Imaxcor,
     +                    T1ln,T1ln_aa,T1ln_bb,T2ln,T2ln_aa,
     +                    T2ln_bb,T2ln_ab,Nsizec,Ioffset,Irrepx,
     +                    Iuhf)

      List1 = 493
      List2 = 453

      Call Getlst(Work(I000),1,1,1,1,List1)
      Ioff = T1ln_aa
      If (Iuhf .Ne. 0) Call Getlst(Work(I000+Ioff),1,1,1,2,List1)

      Ioff = T1ln_aa + Iuhf*T1ln_bb
      Call Getall(Work(I000+Ioff),T2ln_aa,Irrepx,List2+1)
      Ioff = T1ln_aa + Iuhf*T1ln_bb + T2ln_aa

      If (Iuhf .Ne. 0) Call Getall(Work(I000+Ioff),T2ln_bb,
     +                             Irrepx,List2+2)

      Ioff = T1ln_aa + Iuhf*T1ln_bb + T2ln_aa + Iuhf*T2ln_bb
      Call Getall(Work(I000+Ioff),T2ln_ab,Irrepx,List2+3)

#ifdef _DEBUG_LVL0
      call checksum("T1aa :",Work(I000),T1ln_aa)
      Ioff = T1ln_aa
      If (Iuhf.eq.0) call checksum("T1aa :",Work(I000+ioff),T1ln_bb)
      Ioff = T1ln_aa + Iuhf*T1ln_bb
      call checksum("T2aa :",Work(I000+ioff),T2ln_aa)
      Ioff = T1ln_aa + Iuhf*T1ln_bb + T2ln_aa
      If (Iuhf.eq.0) call checksum("T1bb :",Work(I000+ioff),T2ln_bb)
      Ioff = T1ln_aa + Iuhf*T1ln_bb + T2ln_aa + Iuhf*T2ln_bb
      call checksum("T1ab :",Work(I000+ioff),T2ln_ab)
      Write(6,"(a)") "N2C"
#endif 

      Type  = "N2C"
      CCt2  = .False.
      CCt1  = .False.

      List1_h = 493
      List1_t = 493
      Call Transform_2index_pt(Work(I000),Work(I020),Work(I030),
     +                         Work(Iend),Maxcor,List1_h,List1_t,
     +                         T1ln,Iuhf,Irrepx,CCt1,Type)

      List2_h = 453
      List2_t = 453
      Ioff  = T1ln_aa + Iuhf*T1ln_bb

      Call TranSform_4index_aa_pt(Work(I000+Ioff),Work(I020),Work(I030),
     +                            Work(Iend),Maxcor,List2_h,List2_t,
     +                            T2ln,T2ln_aa,Iuhf,Irrepx,Cct2,Type)

      List2_h = 456
      List2_t = 456
      Call Transform_4index_ab_pt(Work(I000+Ioff),Work(I020),Work(I030),
     +                            Work(Iend),Maxcor,List2_h,List2_t,
     +                            T2ln,T2ln_aa,T2ln_bb,Iuhf,Irrepx,CCt2,
     +                            Type)

#ifdef _DEBUG_LVL0
      call checksum("T1aa :",Work(I000),T1ln_aa)
      Ioff = T1ln_aa
      If (Iuhf .Ne. 0) call checksum("T1bb :",Work(I000+Ioff),
     +                                T1ln_bb)
      Ioff = T1ln_aa + Iuhf*T1ln_bb
      call checksum("T2aa :",Work(I000+ioff),T2ln_aa)
      Ioff = T1ln_aa + Iuhf*T1ln_bb + T2ln_aa
      If (Iuhf .Ne. 0) call checksum("T2bb :",Work(I000+ioff),
     +                                T2ln_bb)
      Ioff =  T1ln_aa + Iuhf*T1ln_bb + T2ln_aa + Iuhf*T2ln_bb
      call checksum("T2ab :",Work(I000+ioff),T2ln_ab)
#endif 

C Reorder the vector from T1AA-T1BB-T2AAAA-T2BBBB-T2ABAB to 
C T1AA-T1BB-T2ABAB-T2AAAA-T2BBBB

      T2ln_aabb_ab = T2ln_aa + Iuhf*T2ln_bb + T2ln_ab 
      T2ln_aabb    = T2ln_aa + Iuhf*T2ln_bb 

      I000 = Ione 
      I010 = I000 + Nsizec 
      Iend = I010 + T2len_aabb_ab
      If (Iend .Gt. Maxcor) Call Insmem("form_notrans_pt2",Iend,
     +                                   Maxcor)

      Ioff = T1ln_aa + Iuhf*T1ln_bb 
      Call Dcopy(T2ln_aabb_ab,Work(I000+Ioff),1,Work(I010),1) 

      Ioff = T2ln_aa + Iuhf*T2ln_bb 
      Joff = T1ln_aa + Iuhf*T1ln_bb  
      Call Dcopy(T2ln_ab,Work(I010+Ioff),1,Work(I000+Joff),1)

      Joff = T1ln_aa + Iuhf*T1ln_bb + T2ln_ab
      If (Iuhf .Ne. 0) Call Dcopy(T2ln_aabb,Work(I010),1,
     +                            Work(I000+Joff),1)

      Call Putlst(Work(I000),Ioffset,1,1,Irrepx,List2_feom)
      call checksum("Tx   :",Work(I000),Nsizec)

      Return
      End
