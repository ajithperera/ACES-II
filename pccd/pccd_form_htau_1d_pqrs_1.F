      Subroutine Pccd_form_htau_1d_pqrs_1(Htau_pq,Htau_oo,Htau_vv,
     +                                    Htau_ov,Htau_vo,Ioo,Ivv,
     +                                    Work,Ioffo,Ioffv,Maxcor,
     +                                    Nocc,Nvrt,Nbas)

      Implicit Double Precision(A-H,O-Z)
      Double Precision Ioo,Ivv
      Integer a,b,c,d,an,bn,cn,dn
      Logical ONEP_ONLY
     
      Dimension Htau_pq(Nbas,Nbas)
      Dimension Htau_oo(Nocc,Nocc)
      Dimension Htau_vv(Nocc,Nocc)
      Dimension Htau_ov(Nocc,Nocc)
      Dimension Htau_vo(Nocc,Nocc)
      Dimension Ioo(Nocc,Nocc)
      Dimension Ivv(Nvrt,Nvrt)
      Dimension Work(Maxcor)
      Dimension Ioffo(8)
      Dimension Ioffv(8)
 
      COMMON/ORBR_HESS/ONEP_ONLY

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"

      Data One,Ione,Inull /1.0D0,1,0/

      Call Dzero(Ioo,Nocc*Nocc)
      Call Dzero(Ivv,Nvrt*Nvrt)

      Ioffo(1) = 0
      Ioffv(1) = Nocco(1)
       
      Do Irrep =2, Nirrep
         Ioffo(Irrep)=Ioffo(Irrep-1)+Pop(irrep-1,1)
         Ioffv(Irrep)=Ioffv(Irrep-1)+Vrt(irrep-1,1)
      Enddo 

      Do I = 1, Nocc
         Ioo(I,I) = One
      Enddo 
      Do I = 1, Nvrt
         Ivv(I,I) = One
      Enddo 
     
      Irrepx = Ione
      Ioff   = Ione
      Ispin  = Ione 

C OOOO block (H(ij,kl)=-delta(ik)(F(jl)+D(l,j)) and all permutations)

      List_h = 213 
      Do Irrep_kl = 1, Nirrep
         Irrep_ij = Dirprd(Irrep_kl,Irrepx)

         Nij = Irpdpd(Irrep_ij,14)
         Nkl = Irpdpd(Irrep_kl,14)
         I000 = Ione
         I010 = I000 + Nij*Nkl
         Iend = I010 + Nij*Nkl
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_1",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nij*Nkl)
         Icount = Inull
         Do Irrep_l = 1, Nirrep
            Irrep_k = Dirprd(Irrep_l,Irrep_kl)
            Do l = 1, Pop(Irrep_l,Ispin)
               Do k = 1, Pop(Irrep_k,Ispin)
                  Do Irrep_j = 1, Nirrep
                     Irrep_i = Dirprd(Irrep_j,Irrep_ij)
                     Do j = 1, Pop(Irrep_j,Ispin)
                        Do i = 1, Pop(Irrep_i,Ispin)
                           In = I + Ioffo(Irrep_i)
                           Jn = J + Ioffo(Irrep_j)
                           Kn = K + Ioffo(Irrep_k)
                           Ln = L + Ioffo(Irrep_l)
                           Work(Ioff+Icount) =  
     +                              -Ioo(i,j)*
     +                               (Htau_pq(kn,ln)+Htau_pq(ln,kn)) 
     +                              -Ioo(i,l)*
     +                               (Htau_pq(kn,jn)+Htau_pq(kn,jn))
     +                              +Ioo(j,k)*
     +                               (Htau_pq(in,ln)+Htau_pq(ln,in))
     +                              +Ioo(i,l)*
     +                               (Htau_pq(jn,kn)+Htau_pq(kn,jn))
                          Icount = Icount + 1
                        Enddo 
                      Enddo
                  Enddo                   
               Enddo
            Enddo 
          Enddo
          If (.Not. ONEP_ONLY) Then 
             Call Getlst(Work(I010),1,Nkl,1,Irrep_kl,List_h)
             Call Daxpy(Nij*Nkl,One,Work(I010),1,Work(Ioff),1)
          Endif
          Call Putlst(Work(Ioff),1,Nkl,1,Irrep_kl,List_h)
#ifdef _DEBUG_LVL0
          Call checksum("OOOO    :",Work(Ioff),Nij*Nkl)
#endif 
      Enddo 

C VVVV block (H(ab,cd)=-delta(ac)(F(bd)+F(db)) and all permutations)

      Ioff   = Ione
      List_h = 234
      Do Irrep_cd = 1, Nirrep
         Irrep_ab = Dirprd(Irrep_cd,Irrepx)

         Nab = Irpdpd(Irrep_ab,15)
         Ncd = Irpdpd(Irrep_cd,15)
         I000 = Ione
         I010 = I000 + Nab*Ncd
         Iend = I010 + Nab*Ncd
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_1",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nab*Ncd)
         Icount = Inull 
         Do Irrep_d = 1, Nirrep
            Irrep_c = Dirprd(Irrep_d,Irrep_cd)
            Do d = 1, Vrt(Irrep_d,Ispin)
               Do c = 1, Vrt(Irrep_c,Ispin)
                  Do Irrep_b = 1, Nirrep
                     Irrep_a = Dirprd(Irrep_b,Irrep_ab)
                     Do b = 1, Vrt(Irrep_b,Ispin)
                        Do a = 1, Vrt(Irrep_a,Ispin)
                           An = A + Ioffv(Irrep_a)
                           Bn = B + Ioffv(Irrep_b)
                           Cn = C + Ioffv(Irrep_c)
                           Dn = D + Ioffv(Irrep_d)
                           Work(Ioff+Icount) =
     +                              -Ivv(a,b)*
     +                               (Htau_pq(cn,dn)+Htau_pq(dn,cn))
     +                              -Ivv(c,d)*
     +                               (Htau_pq(an,bn)+Htau_pq(bn,an))
     +                              +Ivv(b,c)*
     +                               (Htau_pq(an,dn)+Htau_pq(dn,an))
     +                              +Ivv(a,d)*
     +                               (Htau_pq(bn,cn)+Htau_pq(cn,bn))
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
          If (.Not. ONEP_ONLY) Then
             Call Getlst(Work(I010),1,Ncd,1,Irrep_cd,List_h)
             Call Daxpy(Nab*Ncd,One,Work(I010),1,Work(Ioff),1)
          Endif 
          Call Putlst(Work(Ioff),1,Ncd,1,Irrep_cd,List_h)
#ifdef _DEBUG_LVL0
          Call checksum("VVVV    :",Work(Ioff),Nab*Ncd)
#endif 
      Enddo

C OOVV block (H(ij,ab)=-delta(ij)(F(ab)+F(ba)))

      Ioff   = Ione
      List_h = 217
      Do Irrep_ab = 1, Nirrep
         Irrep_ij = Dirprd(Irrep_ab,Irrepx)

         Nij = Irpdpd(Irrep_ij,14)
         Nab = Irpdpd(Irrep_ab,15)
         I000 = Ione
         I010 = I000 + Nab*Nij
         Iend = I010 + Nab*Nij
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_1",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nab*Nij)
         Icount = Inull 
         Do Irrep_b = 1, Nirrep
            Irrep_a = Dirprd(Irrep_b,Irrep_ab)
            Do b = 1, Vrt(Irrep_b,Ispin)
               Do a = 1, Vrt(Irrep_a,Ispin)
                  Do Irrep_j = 1, Nirrep
                     Irrep_i = Dirprd(Irrep_j,Irrep_ij)
                     Do j = 1, Pop(Irrep_j,Ispin)
                        Do i = 1, Pop(Irrep_i,Ispin)
                           In = I + Ioffo(Irrep_i)
                           Jn = J + Ioffo(Irrep_j)
                           An = A + Ioffv(Irrep_a)
                           Bn = B + Ioffv(Irrep_b)
                           Work(Ioff+Icount) =  
     +                              -Ioo(i,j)*
     +                               (Htau_pq(an,bn)+Htau_pq(bn,an))
     +                              -Ivv(a,b)*
     +                               (Htau_pq(in,jn)+Htau_pq(jn,in))
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
          If (.Not. ONEP_ONLY) Then
             Call Getlst(Work(I010),1,Nab,1,Irrep_ab,List_h)
             Call Daxpy(Nab*Nij,One,Work(I010),1,Work(Ioff),1)
          Endif 
          Call Putlst(Work(Ioff),1,Nab,1,Irrep_ab,List_h)
#ifdef _DEBUG_LVL0
          Call checksum("OOVV    :",Work(Ioff),Nab*Nij)
#endif 
      Enddo

C VVOO block (H(ab,ij)=-delta(ab)(F(ij)+F(ji)))

      Ioff   = Ione
      List_h = 216
      Do Irrep_ij = 1, Nirrep
         Irrep_ab = Dirprd(Irrep_ij,Irrepx)

         Nab = Irpdpd(Irrep_ab,15)
         Nij = Irpdpd(Irrep_ij,14)
         I000 = Ione
         I010 = I000 + Nab*Nij
         Iend = I010 + Nab*Nij
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_1",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nab*Nij)
         Icount = Inull 
         Do Irrep_j = 1, Nirrep
            Irrep_i = Dirprd(Irrep_j,Irrep_ij)
            Do j = 1, Pop(Irrep_j,Ispin)
               Do i = 1, Pop(Irrep_i,Ispin)
                  Do Irrep_b = 1, Nirrep
                     Irrep_a = Dirprd(Irrep_b,Irrep_ab)
                     Do b = 1, Vrt(Irrep_b,Ispin)
                        Do a = 1, Vrt(Irrep_a,Ispin)
                           In = I + Ioffo(Irrep_i)
                           Jn = J + Ioffo(Irrep_j)
                           An = A + Ioffv(Irrep_a)
                           Bn = B + Ioffv(Irrep_b)
                           Work(Ioff+Icount) = 
     +                              -Ivv(a,b)*
     +                               (Htau_pq(in,jn)+Htau_pq(jn,in))
     +                              -Ioo(i,j)*
     +                               (Htau_pq(an,bn)+Htau_pq(bn,an))
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
          If (.Not. ONEP_ONLY) Then
             Call Getlst(Work(I010),1,Nij,1,Irrep_ij,List_h)
             Call Daxpy(Nab*Nij,One,Work(I010),1,Work(Ioff),1)
          Endif 
          Call Putlst(Work(Ioff),1,Nij,1,Irrep_ij,List_h)
#ifdef _DEBUG_LVL0
          Call checksum("VVOO    :",Work(Ioff),Nab*Nij)
#endif 
      Enddo

C No OVOV or VOVO like combinations 

C These blocks are not necessary for orbital optimizations 

#ifdef _NOSKIP
C OOOV block (H(ij,ka)=-delta(ij)(F(k,a)+F(a,k)))

      Ioff   = Ione
      List_h = 210
      Do Irrep_ka = 1, Nirrep
         Irrep_ij = Dirprd(Irrep_ka,Irrepx)

         Nij = Irpdpd(Irrep_ij,14)
         Nka = Irpdpd(Irrep_ka,18)
         I000 = Ione
         I010 = I000 + Nij*Nka
         Iend = I010 + Nij*Nka
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_1",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nij*Nka)
         Icount = Inull
         Do Irrep_a = 1, Nirrep
            Irrep_k = Dirprd(Irrep_a,Irrep_ka)
            Do a = 1, Vrt(Irrep_a,Ispin)
               Do k = 1, Pop(Irrep_k,Ispin)
                  Do Irrep_j = 1, Nirrep
                     Irrep_i = Dirprd(Irrep_j,Irrep_ij)
                     Do j = 1, Pop(Irrep_j,Ispin)
                        Do i = 1, Pop(Irrep_i,Ispin)
                           In = I + Ioffo(Irrep_i)
                           Jn = J + Ioffo(Irrep_j)
                           An = A + Ioffv(Irrep_a)
                           Kn = K + Ioffo(Irrep_k)
                           Work(Ioff+Icount) = 
     +                              -Ioo(i,j)*
     +                               (Htau_pq(kn,an)+Htau_pq(an,kn))
     +                              +Ioo(k,j)*
     +                               (Htau_pq(in,an)+Htau_pq(an,in))
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
          If (.Not. ONEP_ONLY) Then
              Call Getlst(Work(I010),1,Nka,1,Irrep_ka,List_h)
              Call Daxpy(Nka*Nij,One,Work(I010),1,Work(Ioff),1)
          Endif 
          Call Putlst(Work(Ioff),1,Nka,1,Irrep_ka,List_h)
#ifdef _DEBUG_LVL0
          Call checksum("OOOV    :",Work(Ioff),Nij*Nka)
#endif 
      Enddo

C OVOO block (H(ia,kj)=-delta(j,k)(F(i,a)+F(a,i)))

      Ioff   = Ione
      List_h = 211
      Do Irrep_kj = 1, Nirrep
         Irrep_ia = Dirprd(Irrep_kj,Irrepx)

         Nia = Irpdpd(Irrep_ia,18)
         Nkj = Irpdpd(Irrep_kj,14)
         I000 = Ione
         I010 = I000 + Nia*Nkj
         Iend = I010 + Nia*Nkj
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_1",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nia*Nkj)
         Icount = Inull
         Do Irrep_j = 1, Nirrep
            Irrep_k = Dirprd(Irrep_j,Irrep_kj)
            Do j = 1, Pop(Irrep_j,Ispin)
               Do k = 1, Pop(Irrep_k,Ispin)
                  Do Irrep_a = 1, Nirrep
                     Irrep_i = Dirprd(Irrep_a,Irrep_ia)
                     Do a = 1, Vrt(Irrep_a,Ispin)
                        Do i = 1, Pop(Irrep_i,Ispin)
                           In = I + Ioffo(Irrep_i)
                           Jn = J + Ioffo(Irrep_j)
                           An = A + Ioffv(Irrep_a)
                           Kn = K + Ioffo(Irrep_k)
                           Work(Ioff+Icount) = 
     +                              -Ioo(k,j)*
     +                               (Htau_pq(in,an)+Htau_pq(an,in))
     +                              +Ioo(i,j)*
     +                               (Htau_pq(kn,an)+Htau_pq(an,kn))
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
          If (.Not. ONEP_ONLY) Then
             Call Getlst(Work(I010),1,Nkj,1,Irrep_kj,List_h)
             Call Daxpy(Nia*Nkj,One,Work(I010),1,Work(Ioff),1)
          Endif 
          Call Putlst(Work(Ioff),1,Nkj,1,Irrep_kj,List_h)
#ifdef _DEBUG_LVL0
          Call checksum("OVOO    :",Work(Ioff),Nia*Nkj)
#endif 
      Enddo

C VVVO block (H(ab,ci)=Delta(a,b)*(F(c,i)+F(i,c))

      Ioff   = Ione
      List_h = 230
      Do Irrep_ci = 1, Nirrep
         Irrep_ab = Dirprd(Irrep_ci,Irrepx)

         Nci = Irpdpd(Irrep_ci,11)
         Nab = Irpdpd(Irrep_ab,15)
         I000 = Ione
         I010 = I000 + Nab*Nci
         Iend = I010 + Nab*Nci
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_1",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nab*Nci)
         Icount = Inull
         Do Irrep_i = 1, Nirrep
            Irrep_c = Dirprd(Irrep_i,Irrep_ci)
            Do i = 1, Pop(Irrep_i,Ispin)
               Do c = 1, Vrt(Irrep_c,Ispin)
                  Do Irrep_b = 1, Nirrep
                     Irrep_a = Dirprd(Irrep_b,Irrep_ab)
                     Do b = 1, Vrt(Irrep_b,Ispin)
                        Do a = 1, Vrt(Irrep_a,Ispin)
                           In = I + Ioffo(Irrep_i)
                           An = A + Ioffv(Irrep_a)
                           Bn = B + Ioffv(Irrep_b)
                           Cn = C + Ioffv(Irrep_c)
                           Work(Ioff+Icount) = 
     +                              -Ivv(a,b)*
     +                               (Htau_pq(cn,in)+Htau_pq(in,cn))
     +                              +Ivv(b,c)*
     +                               (Htau_pq(an,in)+Htau_pq(in,an))
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
          If (.Not. ONEP_ONLY) Then
             Call Getlst(Work(I010),1,Nci,1,Irrep_ci,List_h)
             Call Daxpy(Nab*Nci,One,Work(I010),1,Work(Ioff),1)
          Endif 
          Call Putlst(Work(Ioff),1,Nci,1,Irrep_ci,List_h)
#ifdef _DEBUG_LVL0
          Call checksum("VVVO    :",Work(Ioff),Nab*Nci)
#endif 
      Enddo

C VOVV block (H(ci,ab)=Delta(a,b)*(F(c,i)+F(i,c))

      Ioff   = Ione
      List_h = 231
      Do Irrep_ab = 1, Nirrep
         Irrep_ci = Dirprd(Irrep_ab,Irrepx)

         Nci = Irpdpd(Irrep_ci,11)
         Nab = Irpdpd(Irrep_ab,15)
         I000 = Ione
         I010 = I000 + Nab*Nci
         Iend = I010 + Nab*Nci
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_1",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nab*Nci)
         Icount = Inull
         Do Irrep_b = 1, Nirrep
            Irrep_a = Dirprd(Irrep_b,Irrep_ab)
            Do b = 1, Vrt(Irrep_b,Ispin)
               Do a = 1, Vrt(Irrep_a,Ispin)
                  Do Irrep_i = 1, Nirrep
                     Irrep_c = Dirprd(Irrep_i,Irrep_ci)
                     Do i = 1, Pop(Irrep_i,Ispin)
                        Do c = 1, Vrt(Irrep_c,Ispin)
                           In = I + Ioffo(Irrep_i)
                           An = A + Ioffv(Irrep_a)
                           Bn = B + Ioffv(Irrep_b)
                           Cn = C + Ioffv(Irrep_c)
                           Work(Ioff+Icount) = 
     +                              -Ivv(a,b)*
     +                               (Htau_pq(cn,in)+Htau_pq(in,cn))
     +                              +Ivv(c,b)*
     +                               (Htau_pq(an,in)+Htau_pq(in,an))
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
          If (.Not. ONEP_ONLY) Then
             Call Getlst(Work(I010),1,Nab,1,Irrep_ab,List_h)
             Call Daxpy(Nab*Nci,One,Work(I010),1,Work(Ioff),1)
          Endif 
          Call Putlst(Work(Ioff),1,Nab,1,Irrep_ab,List_h)
#ifdef _DEBUG_LVL0
          Call checksum("VOVV    :",Work(Ioff),Nab*Nci)
#endif 
      Enddo
#endif 

#ifdef _DEBUG_LVL0
      Call Pccd_check_hess(Work,Maxcor)
#endif 

      Return
      End 


