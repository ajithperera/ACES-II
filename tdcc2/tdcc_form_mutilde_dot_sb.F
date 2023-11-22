      Subroutine Tdcc_form_mutilde_dot_sb(Work,Memleft,Irrepx,Iuhf)

      Implicit Double Precision (A-H, O-Z)
      Logical Rhf
C
      Dimension Work(Memleft)

      Dimension Ioff_hb_irrep(8),Ioff_mut_irrep(8)
      Dimension Ioff_mut_d_irrep(8)
      Integer Hbar_hphh_list
C
      Data One,Onem,Zero/1.0D0,-1.0D0,0.0D0/

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "symloc.com"
#include "sym.com"

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Entered  tdcc_form_mutilde_dot_sb---"
        Write(6,*)
#endif
      Rhf = .False.
      If (Iuhf .Eq. 0) Rhf = .True.
      Mut_d_list     = 396

C Z(A,I) = Hbar(IE,NM) * Mut(NM,AE) +  Hbar(Ie,Nm) * Mut(Nm,Ae)
C Z(a,i) = Hbar(ie,nm) * Mut(nm,ae) +  Hbar(iE,nM) * Mut(nM,aE)
  
      Do Ispin = 1, 1+Iuhf

         Ndim_ai = Irpdpd(Irrepx,8+Ispin)

C Hbar Read as Hbar(Nm,Ie)  Ispin=1
C              Hbar(Nm,Ei)  Ispin=2
C Mu_tilde read as Mut(Ae,Nm) 


         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Ndim_ai
         Call Dzero(Work(I000),Ndim_ai)

         Do Irrep_nm = 1, Nirrep 
            Hbar_hphh_list = 111 - Ispin
            Mut_list       = 326

            Irrep_ae = Dirprd(Irrep_nm,Irrepx)
            Irrep_ie = Irrep_nm

            Ndim_ae_z = Irpdpd(Irrep_ae,Isytyp(1,Mut_list))
            Ndim_nm_z = Irpdpd(Irrep_nm,Isytyp(2,Mut_list))
            Ndim_nm_h = Irpdpd(Irrep_nm,Isytyp(1,Hbar_hphh_list))
            Ndim_ie_h = Irpdpd(Irrep_ie,Isytyp(2,Hbar_hphh_list))
            Max_z     = Max(Ndim_ae_z,Ndim_Nm_z)
            Max_h     = Max(Ndim_nm_h,Ndim_ie_h)

            Ibgn = I010
            I020 = I010 + Ndim_ie_h * Ndim_nm_h
            I030 = I020 + Max(Ndim_ae_z*Ndim_nm_z,3*Max_h)
            I040 = I030 + Max(Ndim_ae_z*Ndim_nm_z,3*Max_z)
            Iend = I040 

           If (Iend .Lt. Memleft) Then 
               Call Getlst(Work(I010),1,Ndim_ie_h,1,Irrep_ie,
     +                     Hbar_hphh_list)
C Hbar(Nm,Ei) Ispin=2
            
               If (Ispin .EQ. 1) Then
       
C Hbar(Nm,Ie) -> Hbar(Nm,eI)

                  Ibgn  = I020 
                  Iptr1 = Ibgn
                  Iptr2 = Iptr1 + Max_h
                  Iptr3 = Iptr2 + Max_h
                  Iptr4 = Iptr3 + Max_h
                  Iend  = Iptr4
                  If (Iend .Gt. Memleft) Call Insmem(
     +                         "@-Tdcc_form_mutilde_dot_sb",
     +                          Memleft,Iend)
  
                  Call Symtr1(Irrep_ie,Pop(1,1),Vrt(1,2),Ndim_nm_h,
     +                        Work(I010),Work(Iptr1),Work(Iptr2),
     +                        Work(Iptr3))
               Endif 

C Mu(Ae,Nm) -> Mu(Nm,Ae) -> Mu(Nm,eA) Ispin=1
C Mu(aE,nM) -> Mu(nM,aE) -> Mu(Nm,Ae) Ispin=2 

               Call Getlst(Work(I030),1,Ndim_nm_z,1,Irrep_nm,
     +                     Mut_list)

               Call Transp(Work(I030),Work(I020),Ndim_nm_z,Ndim_ae_z)

               If (Ispin .EQ. 1) Then
                  Ibgn  = I030 
                  Iptr1 = Ibgn
                  Iptr2 = Iptr1 + Max_z
                  Iptr3 = Iptr2 + Max_z
                  Iptr4 = Iptr3 + Max_z
                  Iend  = Iptr4

                  If (Iend .Gt. Memleft) Call Insmem(
     +                         "@-Tdcc_form_mutilde_dot_sb",
     +                          Memleft,Iend)

                  Call Symtr1(Irrep_ae,Vrt(1,1),Vrt(1,2),Ndim_nm_z,
     +                         Work(I020),Work(IPtr1),Work(Iptr2),
     +                         Work(Iptr3))
                  If (RHF) Then
                     Call Spinad3(Irrep_nm,Pop(1,1),Ndim_nm_z,
     +                            Ndim_ae_z,Work(I020),Work(Iptr1),
     +                            Work(Iptr2))
                  Endif 
               Endif

C Z(A,I) =  Mu(Nm,eA)^t * Hbar(Nm,eI)  Ispin=1
C Z(a,i) =  Mu(Nm,aE)^t * Hbar(Nm,Ei)  Ispin=2

              Ioff_hb    = 0
              Ioff_mut   = 0
              Ioff_mut_d = 0

              Do Irrep_p = 1, Nirrep
                 Irrep_pz = Dirprd(Irrep_p,Irrepx)
                 Irrep_ph = Dirprd(Irrep_p,Irrep_ie)
                 Irrep_pm = Dirprd(Irrep_p,Irrep_ae)

                 Ioff_hb_irrep(Irrep_p)    = I010 + Ioff_hb
                 Ioff_mut_irrep(Irrep_p)   = I020 + Ioff_mut
                 Ioff_mut_d_irrep(Irrep_p) = I000 + Ioff_mut_d

                 Ioff_hb  = Ioff_hb +  Ndim_nm_h*
     +                                 Vrt(Irrep_ph,3-Ispin)
     +                                *Pop(Irrep_p,Ispin)
                 Ioff_mut = Ioff_mut + Ndim_nm_z*
     +                                 Vrt(Irrep_pm,3-Ispin)
     +                                *Vrt(Irrep_p,Ispin)
                 Ioff_mut_d = Ioff_mut_d + Vrt(Irrep_pz,Ispin)
     +                                    *Pop(Irrep_p,Ispin)  
              Enddo 
 
              Do Irrep_i = 1, Nirrep 
                 Irrep_e = Dirprd(Irrep_i,Irrep_ie)
                 Irrep_a = Dirprd(Irrep_e,Irrep_ae)

                 Nrow = Vrt(Irrep_a,Ispin)
                 Ncol = Pop(Irrep_i,Ispin)
                 Nsum = Ndim_nm_h * Vrt(Irrep_e,3-Ispin)
                 
                 Ioff_hb    = Ioff_hb_irrep(Irrep_i)
                 Ioff_mut   = Ioff_mut_irrep(Irrep_a)
                 Ioff_mut_d = Ioff_mut_d_irrep(Irrep_i)

                 Call Xgemm("T","N",Nrow,Ncol,Nsum,Onem,
     +                       Work(Ioff_mut),Nsum,
     +                       Work(Ioff_hb),Nsum,One,
     +                       Work(Ioff_mut_d),Nrow)
              Enddo 
#ifdef _DEBUG_LVLM
      call Checksum("Tdcc_form_mu_dot_sb,Mutilde_dot(1):",
     +               Work(I000),Ndim_ai,S)
#endif
          Else 

              Write(6,"(3a)") "Insufficient amount of memory", 
     +                        " to perform the contraction in",
     +                        " memory."
              Call Errex
          Endif 

C AAAA and BBBB contributions 
C Z(A,I) = Hbar(IE,NM) * Mut(NM,AE) 
C Z(a,i) = Hbar(ie,nm) * Mut(nm,ae) 

C Read Hbar(N<M,IE) and Hbar(n<m,ie) and Mut(A<E,N<M) and Mut(a<e,n<m)

           If (.NOT. RHF) Then

              Hbar_hphh_list = 106 + Ispin
              Mut_list       = 323 + Ispin

              Ndim_nm_c_h = Irpdpd(Irrep_nm,
     +                           Isytyp(1,Hbar_hphh_list))
              Ndim_ie_f_h = Irpdpd(Irrep_ie,
     +                           Isytyp(2,Hbar_hphh_list))

              Ndim_ae_c_z = Irpdpd(Irrep_ae,Isytyp(1,Mut_list))
              Ndim_nm_c_z = Irpdpd(Irrep_nm,Isytyp(2,Mut_list))

              Ndim_nm_f = Irpdpd(Irrep_nm,20+Ispin)
              Ndim_ae_f = Irpdpd(Irrep_ae,18+Ispin)

              Max_h = Max(Ndim_ie_f_h,Ndim_nm_c_h)
              Max_z = Max(Ndim_nm_c_z,Ndim_ae_c_z)

              Ibgn = I010
              I020 = I010 + Ndim_nm_c_h*Ndim_ie_f_h
              I030 = I020 + Max(Ndim_nm_c_h*Ndim_ae_f,
     +                          Ndim_nm_c_z*Ndim_ae_c_z,3*Max_h)
              I040 = I030 + Max(Ndim_nm_c_z*Ndim_ae_f,
     +                          3*Max_z)
              Iend = I040 

              If (Iend .Lt. Memleft) Then 
             
C Read Hbar(N<M,IE) (and Hbar(n<m,ie))
              
                 Call Getlst(Work(I010),1,Ndim_ie_f_h,1,Irrep_ie,
     +                       Hbar_hphh_list)

C Hbar(N<M,IE) -> Hbar(N<M,EI)  Ispin=1
C Hbar(n<m,ie) -> Hbar(n<m,ei)  Ispin=2

                  Ibgn  = I020
                  Iptr1 = Ibgn
                  Iptr2 = Iptr1 + Max_h
                  Iptr3 = Iptr2 + Max_h
                  Iptr4 = Iptr3 + Max_h
                  Iend  = Iptr4
                  If (Iend .Gt. Memleft) Call Insmem(
     +                         "@-Tdcc_form_mutilde_dot_sb",
     +                          Memleft,Iend)

                 Call Symtr1(Irrep_ie,Pop(1,Ispin),Vrt(1,Ispin),
     +                      Ndim_nm_c_h,Work(I010),Work(Iptr1),
     +                      Work(Iptr2),Work(Iptr3))

C Read Mut(A<E,N<M) Ispin=1
C Read Mut(a<e,n<m) Ispin=2

                 Call Getlst(Work(I030),1,Ndim_nm_c_z,1,Irrep_nm,
     +                       Mut_list)

C Mut(A<E,N<M) -> Mut(N<M,A<E); Mut(a<e,n<m) -> Mut(n<m,a<e) Ispin=1,2

                 Call Transp(Work(I030),Work(I020),Ndim_nm_c_z,
     +                       Ndim_ae_c_z)

C Mut(N<M,A<E) -> Mut(N<M,AE);  Mut(n<m,a<e) -> Mut(n<m,ae) Ispin=1,2

                 Call Symexp(Irrep_ae,Vrt(1,Ispin),Ndim_nm_c_z,
     +                       Work(I020))

C Z(A,I) = Mut(N<M,AE)^t * Hbar(N<M,EI)  Ispin=1
C Z(a,i) = Mut(n<m,ae)^t * Hbar(n<m,ei)  Ispin=2 
   
                 Ioff_hb    = 0
                 Ioff_mut   = 0
                 Ioff_mut_d = 0

                 Do Irrep_p = 1, Nirrep
                     Irrep_pz = Dirprd(Irrep_p,Irrepx)
                     Irrep_ph = Dirprd(Irrep_p,Irrep_ie)
                     Irrep_pm = Dirprd(Irrep_p,Irrep_ae)

                     Ioff_hb_irrep(Irrep_p)    = I010 + Ioff_hb
                     Ioff_mut_irrep(Irrep_p)   = I020 + Ioff_mut
                     Ioff_mut_d_irrep(Irrep_p) = I000 + Ioff_mut_d

                     Ioff_hb = Ioff_hb   + Ndim_nm_c_h*
     +                                     Vrt(Irrep_ph,Ispin)
     +                                    *Pop(Irrep_p,Ispin)
                     Ioff_mut = Ioff_mut + Ndim_nm_c_z*
     +                                     Vrt(Irrep_pm,Ispin)
     +                                    *Vrt(Irrep_p,Ispin)
                     Ioff_mut_d = Ioff_mut_d + 
     +                                    Pop(Irrep_p,Ispin)*
     +                                    vrt(Irrep_pz,Ispin)
                  Enddo
               
                 Do Irrep_i = 1, Nirrep
                    Irrep_e = Dirprd(Irrep_i,Irrep_ie)
                    Irrep_a = Dirprd(Irrep_e,Irrep_ae)

                    Nrow = Vrt(Irrep_a,Ispin)
                    Ncol = Pop(Irrep_i,Ispin)
                    Nsum = Ndim_nm_c_h*Vrt(Irrep_e,Ispin)
          
                    Ioff_hb    = Ioff_hb_irrep(Irrep_i)
                    Ioff_mut   = Ioff_mut_irrep(Irrep_a)
                    Ioff_mut_d = Ioff_mut_d_irrep(Irrep_i)

                    Call Xgemm("T","N",Nrow,Ncol,Nsum,One,
     +                         Work(Ioff_mut),
     +                         Nsum,Work(Ioff_hb),Nsum,One,
     +                         Work(Ioff_mut_d),Nrow)
                 Enddo
#ifdef _DEBUG_LVLM
      call Checksum("Tdcc_form_mu_dot_sb,Mutilde_dot(2):",
     +               Work(I000),Ndim_ai,S)
#endif
             Else

                 Write(6,"(3a)") "Insufficient amount of memory",
     +                           " to perform the contraction in",
     +                           " memory."
                 Call Errex
             Endif 

           Endif

         Enddo
#ifdef _DEBUG_LVL0
      call Checksum("Tdcc_form_mu_dot_sb,Mutilde_dot(f):",
     +               Work(I000),Ndim_ai,S)
#endif
         Call Getlst(Work(I010),1,1,1,Ispin,Mut_d_list)
         Call Daxpy(Ndim_ai,One,Work(I010),1,Work(I000),1)
         Call Putlst(Work(I000),1,1,1,Ispin,Mut_d_list)

      Enddo 

      Return
      End 

