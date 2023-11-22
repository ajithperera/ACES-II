#include "flags.h"
      Subroutine Pccd_form_htau_1d_pqrs_2(H_pq,D_pq,H_oo,H_vv,H_ov,H_vo,
     +                                    D_oo,D_vv,D_ov,D_vo,Work,
     +                                    Ioffo,Ioffv,Maxcor,Nocc,
     +                                    Nvrt,Nbas)

      Implicit Double Precision(A-H,O-Z)
      Double Precision Ioo,Ivv
      Integer a,b,c,d,An,Bn,Cn,Dn

      Logical Do_oo,Do_vv,Do_vo,Do_ov,Do_hoovv,Do_doovv,
     +        Do_hvvoo,Do_dvvoo,Do_hvoov,Do_dvoov,Do_hovoo,
     +        Do_Dovoo,Do_hvooo,Do_dvooo,Do_hvvov,Do_dvvov,
     +        Do_hvovv,Do_dvovv,Do_hvvvo,Do_dvvvo,Do_hovvv,
     +        Do_dovvv,Do_hovvo,Do_dovvo,Do_hooov,Do_dooov,
     +        Do_hoovo,Do_doovo

      Dimension H_pq(Nbas,Nbas)
      Dimension D_pq(Nbas,Nbas)
      Dimension H_oo(Nocc,Nocc)
      Dimension H_vv(Nvrt,Nvrt)
      Dimension H_ov(Nocc,Nvrt)
      Dimension H_vo(Nvrt,Nocc)
      Dimension D_oo(Nocc,Nocc)
      Dimension D_vv(Nvrt,Nvrt)
      Dimension D_ov(Nocc,Nvrt)
      Dimension D_vo(Nvrt,Nocc)
      Dimension Ioffo(8)
      Dimension Ioffv(8)
      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "flags.com"

      Data One,Ione,Inull /1.0D0,1,0/
 
      Hoo = Ddot(Nocc*Nocc,H_oo,1,H_oo,1)
      Hvv = Ddot(Nvrt*Nvrt,H_vv,1,H_vv,1)
      Hvo = Ddot(Nvrt*Nocc,H_vo,1,H_vo,1)
      Hov = Ddot(Nvrt*Nocc,H_ov,1,H_ov,1)
      Doo = Ddot(Nocc*Nocc,D_oo,1,D_oo,1)
      Dvv = Ddot(Nvrt*Nvrt,D_vv,1,D_vv,1)
      Dvo = Ddot(Nvrt*Nocc,D_vo,1,D_vo,1)
      Dov = Ddot(Nvrt*Nocc,D_ov,1,D_ov,1)

      Cc_conv = Iflags(h_IFLAGS_cc_conv)
      Thres   = Cc_con*1.0D-03

      Do_oo    = .True.
      Do_vv    = .True.
      Do_vo    = .True.
      Do_ov    = .True.  

      Do_hoovv = .True.
      Do_doovv = .True.

      Do_hvvoo = .True.
      Do_dvvoo = .True.

      Do_hvoov = .True.
      Do_dvoov = .True.

      Do_hovvo = .True.
      Do_dovvo = .True.

      Do_hovoo = .True.
      Do_hovoo = .True.

      Do_hvooo = .True.
      Do_dvooo = .True.

      Do_hvvov = .True.
      Do_dvvov = .True.

      Do_hvovv = .True.
      Do_dvovv = .True.

      Do_hvvvo = .True.
      Do_Dvvvo = .True.

      Do_hovvv = .True.
      Do_dovvv = .True.

      Do_hovvo = .True.
      Do_dovvo = .True.

      Do_hooov = .True.
      Do_dooov = .True.

      Do_hoovo = .False.
      Do_doovo = .False. 


      If (Hoo .lt. Thres .OR. Doo .lt. Thres) Do_oo = .False.
      If (Hvv .lt. Thres .OR. Dvv .lt. Thres) Do_vv = .False.
      If (Hvo .lt. Thres .OR. Dvo .lt. Thres) Do_vo = .False.
      If (Hov .lt. Thres .OR. Dov .lt. Thres) Do_ov = .False.

      If (Hoo .lt. Thres .OR. Dvv .lt. Thres) Do_hoovv = .False.
      If (Doo .lt. Thres .OR. Hvv .lt. Thres) Do_doovv = .False.

      If (Hvv .lt. Thres .OR. Doo .lt. Thres) Do_hvvoo = .False.
      If (Dvv .lt. Thres .OR. Doo .lt. Thres) Do_dvvoo = .False.

      If (Hvo .lt. Thres .OR. Dov .lt. Thres) Do_hvoov = .False.
      If (Dvo .lt. Thres .OR. Hov .lt. Thres) Do_dvoov = .False.

      If (Hov .lt. Thres .OR. Doo .lt. Thres) Do_hovoo = .False.
      If (Dov .lt. Thres .OR. Hoo .lt. Thres) Do_Dovoo = .False.

      If (Hvo .lt. Thres .OR. Doo .lt. Thres) Do_hvooo = .False.
      If (Dvo .lt. Thres .OR. Hoo .lt. Thres) Do_dvooo = .False.

      If (Hvv .lt. Thres .OR. Dov .lt. Thres) Do_hvvov = .False.
      If (Dvv .lt. Thres .OR. Hov .lt. Thres) Do_dvvov = .False.

      If (Hvo .lt. Thres .OR. Dvv .lt. Thres) Do_hvovv = .False.
      If (Dvo .lt. Thres .OR. Dvv .lt. Thres) Do_dvovv = .False.

      If (Hvv .lt. Thres .OR. Dvo .lt. Thres) Do_hvvvo = .False.
      If (Dvv .lt. Thres .OR. Hvo .lt. Thres) Do_dvvvo = .False.

      If (Hov .lt. Thres .OR. Dvv .lt. Thres) Do_hovvv = .False.
      If (Dov .lt. Thres .OR. Hvv .lt. Thres) Do_dovvv = .False.

      If (Hov .lt. Thres .OR. Dvo .lt. Thres) Do_hovvo = .False.
      If (Dov .lt. Thres .OR. Hvo .lt. Thres) Do_dovvo = .False.

      If (Hoo .lt. Thres .OR. Dov .lt. Thres) Do_hooov = .False.
      If (Doo .lt. Thres .OR. Hov .lt. Thres) Do_Dooov = .False.

      If (Hoo .lt. Thres .OR. Dvo .lt. Thres) Do_hoovo = .False.
      If (Doo .lt. Thres .OR. Hvo .lt. Thres) Do_doovo = .False.

C OOOO block - H(ij,kl)

      Irrepx = Ione
      Ispin  = Ione
      Write(6,*) 

      If (Do_oo) Then 
      Ioff   = Ione
      List_h = 213
      Do Irrep_kl = 1, Nirrep
         Irrep_ij = Dirprd(Irrep_kl,Irrepx)

         Nij = Irpdpd(Irrep_ij,14)
         Nkl = Irpdpd(Irrep_kl,14)
         I000 = Ione
         I010 = I000 + Nij*Nkl
         Iend = I010 + Nij*Nkl
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_2",
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
     +                                   H_pq(in,jn)*D_pq(ln,kn)+ 
     +                                   H_pq(kn,ln)*D_pq(jn,in)-
     +                                   H_pq(jn,kn)*D_pq(ln,in)-
     +                                   H_pq(in,ln)*D_pq(kn,jn)
                           Icount = Icount + 1
                        Enddo 
                      Enddo
                  Enddo                   
               Enddo
            Enddo 
          Enddo
#ifdef _DEBUG_LVL0
          Call checksum("OOOO    :",Work(Ioff),Nij*Nkl)
#endif
          Call Getlst(Work(I010),1,Nkl,1,Irrep_kl,List_h)
          Call Daxpy(Nij*Nkl,One,Work(I010),1,Work(Ioff),1)
          Call Putlst(Work(Ioff),1,Nkl,1,Irrep_kl,List_h)
      Enddo 
      Endif 

C VVVV block - H(ab,cd)
      If (Do_vv) Then
      Ioff   = Ione
      List_h = 234
      Do Irrep_cd = 1, Nirrep
         Irrep_ab = Dirprd(Irrep_cd,Irrepx)

         Nab = Irpdpd(Irrep_cd,15)
         Ncd = Irpdpd(Irrep_ab,15)
         I000 = Ione
         I010 = I000 + Nab*Ncd
         Iend = I010 + Nab*Ncd
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_2",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nab*Ncd)
         Icount=Inull
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
     +                                   H_pq(an,bn)*D_pq(dn,cn)+
     +                                   H_pq(cn,dn)*D_pq(bn,an)-
     +                                   H_pq(bn,cn)*D_pq(dn,an)-
     +                                   H_pq(an,dn)*D_pq(cn,bn)
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
#ifdef _DEBUG_LVL0
          Call checksum("VVVV    :",Work(Ioff),Nab*Ncd)
#endif 
          Call Getlst(Work(I010),1,Ncd,1,Irrep_cd,List_h)
          Call Daxpy(Nab*Ncd,One,Work(I010),1,Work(Ioff),1)
          Call Putlst(Work(Ioff),1,Ncd,1,Irrep_cd,List_h)
      Enddo
      Endif 

C VOOV block - H(ai,jb) 

      If (Do_vo .or. Do_ov .or. Do_hoovv .or. 
     +    Do_hvvoo .or. Do_doovv .or. Do_dvvoo) Then
      Ioff   = Ione
      List_h = 223
      Do Irrep_jb = 1, Nirrep
         Irrep_ai = Dirprd(Irrep_jb,Irrepx)

         Nai = Irpdpd(Irrep_ai,11)
         Njb = Irpdpd(Irrep_jb,18)
         I000 = Ione
         Iend = I000 + Nai*Njb
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_2",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nai*Njb)
         Icount = Inull
         Do Irrep_b = 1, Nirrep
            Irrep_j = Dirprd(Irrep_b,Irrep_jb)
            Do b = 1, Vrt(Irrep_b,Ispin)
               Do j = 1, Pop(Irrep_j,Ispin)
                  Do Irrep_i = 1, Nirrep
                     Irrep_a = Dirprd(Irrep_i,Irrep_ai)
                     Do i = 1, Pop(Irrep_i,Ispin)
                        Do a = 1, Vrt(Irrep_a,Ispin)
                           In = I + Ioffo(Irrep_i)
                           Jn = J + Ioffo(Irrep_j)
                           An = A + Ioffv(Irrep_a)
                           Bn = B + Ioffv(Irrep_b)
                           Work(Ioff+Icount) = 
     +                                   H_pq(an,in)*D_pq(bn,jn)+
     +                                   H_pq(jn,bn)*D_pq(in,an)-
     +                                   H_pq(in,jn)*D_pq(bn,an)-
     +                                   H_pq(an,bn)*D_pq(jn,in)
#ifdef _NOSKIP
     + 
     +                                   H_pq(an,jn)*D_pq(bn,in)+
     +                                   H_pq(in,bn)*D_pq(jn,an)-
     +                                   H_pq(in,jn)*D_pq(bn,an)-
     +                                   H_pq(an,bn)*D_pq(jn,in)
#endif 
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
          Call Putlst(Work(Ioff),1,Njb,1,Irrep_jb,List_h)
#ifdef _DEBUG_LVL0
          Call checksum("VOOV    :",Work(Ioff),Nai*Njb)
#endif 
      Enddo
      Endif 

C OVVO block - H(ia,bj)

      If (Do_hoovv .or. Do_hvvoo .or. Do_hvoov .or. Do_ov .or.
     +    Do_doovv .or. Do_dvvoo .or. Do_dvoov) Then
      Ioff   = Ione
      List_h = 218
      Do Irrep_bj = 1, Nirrep
         Irrep_ia = Dirprd(Irrep_bj,Irrepx)

         Nia = Irpdpd(Irrep_ia,18)
         Nbj = Irpdpd(Irrep_bj,11)
         I000 = Ione
         Iend = I000 + Nia*Nbj
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_2",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nia*Nbj)
         Icount = Inull
         Do Irrep_b = 1, Nirrep
            Irrep_j = Dirprd(Irrep_b,Irrep_bj)
            Do b = 1, Vrt(Irrep_b,Ispin)
               Do j = 1, Pop(Irrep_j,Ispin)
                  Do Irrep_a = 1, Nirrep
                     Irrep_i = Dirprd(Irrep_a,Irrep_ia)
                     Do i = 1, Pop(Irrep_i,Ispin)
                        Do a = 1, Vrt(Irrep_a,Ispin)
                           In = I + Ioffo(Irrep_i)
                           Jn = J + Ioffo(Irrep_j)
                           An = A + Ioffv(Irrep_a)
                           Bn = B + Ioffv(Irrep_b)
                           Work(Ioff+Icount) =
     +                                   H_pq(in,an)*D_pq(jn,bn)+
     +                                   H_pq(bn,jn)*D_pq(an,in)-
     +                                   H_pq(an,bn)*D_pq(jn,in)-
     +                                   H_pq(in,jn)*D_pq(bn,an)
#ifdef _NOSKIP
     +                                   D_pq(in,bn)*H_pq(jn,an)+
     +                                   D_pq(an,jn)*H_pq(bn,in)-
     +                                   D_pq(an,jn)*H_pq(jn,in)-
     +                                   D_pq(in,jn)*H_pq(bn,an)
#endif 
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
          Call Putlst(Work(Ioff),1,Nbj,1,Irrep_bj,List_h)
#ifdef _DEBUG_LVL0
          Call checksum("OVVO    :",Work(Ioff),Nia*Nbj)
#endif 
      Enddo
      Endif 

C OVOV block - H(ia,jb)

      If (Do_vo .or. Do_ov .or. Do_hoovv .or. Do_hvvoo .or.
     +    Do_doovv .or. Do_dvvoo) Then
      Ioff   = Ione
      List_h = 225
      Do Irrep_jb = 1, Nirrep
         Irrep_ia = Dirprd(Irrep_jb,Irrepx)

         Nia = Irpdpd(Irrep_ia,11)
         Njb = Irpdpd(Irrep_jb,11)
         I000 = Ione
         Iend = I000 + Nia*Njb
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_2",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nia*Njb)
         Icount = Inull
         Do Irrep_b = 1, Nirrep
            Irrep_j = Dirprd(Irrep_b,Irrep_jb)
            Do b = 1, Vrt(Irrep_b,Ispin)
               Do j = 1, Pop(Irrep_j,Ispin)
                  Do Irrep_a = 1, Nirrep
                     Irrep_i = Dirprd(Irrep_a,Irrep_ia)
                     Do i = 1, Pop(Irrep_i,Ispin)
                        Do a = 1, Vrt(Irrep_a,Ispin)
                           In = I + Ioffo(Irrep_i)
                           Jn = J + Ioffo(Irrep_j)
                           An = A + Ioffv(Irrep_a)
                           Bn = B + Ioffv(Irrep_b)
                           Work(Ioff+Icount) =
     +                                   H_pq(in,jn)*D_pq(bn,an)+
     +                                   H_pq(jn,bn)*D_pq(an,in)-
     +                                   H_pq(an,jn)*D_pq(bn,in)-
     +                                   H_pq(in,bn)*D_pq(jn,an)
#ifdef _NOSKIP
     +
     +                                   D_pq(in,jn)*H_pq(bn,an)+
     +                                   D_pq(an,bn)*H_pq(jn,in)-
     +                                   D_pq(an,jn)*H_pq(bn,in)-
     +                                   D_pq(in,bn)*H_pq(jn,an)
#endif 
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
          Call Putlst(Work(Ioff),1,Njb,1,Irrep_jb,List_h)
#ifdef _DEBUG_LVL0
          Call checksum("OVOV    :",Work(Ioff),Nia*Njb)
#endif 
      Enddo
      Endif 
C These contributions do enter into orbital optimizations 
#ifdef _NOSKIP
C OOOV block H(ij,ka)

      If (Do_hooov .or. Do_hovoo .or. Do_hvooo .or.
     +    Do_dooov .or. Do_dovoo .or. Do_dvooo ) Then
      Ioff   = Ione
      List_h = 210
      Do Irrep_ka = 1, Nirrep
         Irrep_ij = Dirprd(Irrep_ka,Irrepx)

         Nka = Irpdpd(Irrep_ka,18)
         Nij = Irpdpd(Irrep_ij,14)
         I000 = Ione
         I010 = I000 + Nij*Nka
         Iend = I010 + Nij*Nka
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_2",
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
     +                                   H_pq(in,jn)*D_pq(an,kn)+
     +                                   H_pq(kn,an)*D_pq(jn,in)-
     +                                   H_pq(jn,kn)*D_pq(an,in)-
     +                                   H_pq(in,an)*D_pq(kn,jn)

#ifdef _NOSKIP
     +                                   D_pq(in,kn)H_pq(an,jn)+
     +                                   D_pq(jn,an)H_pq(kn,in)-
     +                                   D_pq(an,kn)H_pq(in,jn)-
     +                                   D_pq(jn,in)H_pq(kn,an) 
#endif 
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
#ifdef _DEBUG_LVL0
          Call checksum("OOOV    :",Work(Ioff),Nka*Nij)
#endif 
          Call Getlst(Work(I010),1,Nka,1,Irrep_ka,List_h)
          Call Daxpy(Nka*Nij,One,Work(I010),1,Work(Ioff),1)
          Call Putlst(Work(Ioff),1,Nka,1,Irrep_ka,List_h)
      Enddo
      Endif 

C OVOO block H(ia,kj)

      If (Do_hoovo .or. Do_hvooo .or. Do_hovoo .or. 
     +    Do_hooov .or. Do_doovo .or. Do_dvooo .or. 
     +    Do_dovoo .or. Do_dooov) Then
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
     +                                   H_pq(in,an)*D_pq(jn,kn)+
     +                                   H_pq(kn,jn)*D_pq(an,in)-
     +                                   H_pq(an,kn)*D_pq(jn,in)-
     +                                   H_pq(in,jn)*D_pq(kn,an)

#ifdef _NOSKIP     
     +                                   D_pq(in,kn)*H_pq(an,jn)+
     +                                   D_pq(an,jn)*H_pq(kn,in)-
     +                                   D_pq(kn,an)*H_pq(in,jn)-
     +                                   D_pq(in,jn)*H_pq(kn,an)
#endif 
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
#ifdef _DEBUG_LVL0
          Call checksum("OVOO    :",Work(Ioff),Nia*Nkj)
#endif 
          Call Getlst(Work(I010),1,Nkj,1,Irrep_kj,List_h)
          Call Daxpy(Nia*Nkj,One,Work(I010),1,Work(Ioff),1)
          Call Putlst(Work(Ioff),1,Nkj,1,Irrep_kj,List_h)
      Enddo
      Endif

C VVVO block H(ab,ci)

      If (Do_hvvov .or. do_hvovv .or. Do_dvvov .or. do_dvovv ) Then 
      Ioff   = Ione
      List_h = 230
      Do Irrep_ci = 1, Nirrep
         Irrep_ab = Dirprd(Irrep_ci,Irrepx)

         Nci = Irpdpd(Irrep_ci,11)
         Nab = Irpdpd(Irrep_ab,15)
         I000 = Ione
         I010 = I000 + Nci*Nab
         Iend = I010 + Nci*Nab
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_1",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nci*Nab)
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
     +                                   H_pq(an,bn)*D_pq(in,cn)+
     +                                   H_pq(cn,in)*D_pq(bn,an)-
     +                                   H_pq(bn,cn)*D_pq(in,an)-
     +                                   H_pq(an,in)*D_pq(cn,bn)

#ifdef _NOSKIP
     +                         
     +                                   D_pq(an,cn)*H_pq(in,bn)+
     +                                   D_pq(bn,in)*H_pq(cn,an)-
     +                                   D_pq(bn,cn)*H_pq(in,an)-
     +                                   D_pq(an,in)*H_pq(cn,bn)
#endif 
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
#ifdef _DEBUG_LVL0
          Call checksum("VVVO    :",Work(Ioff),Nci*Nab)
#endif 
          Call Getlst(Work(I010),1,Nci,1,Irrep_ci,List_h)
          Call Daxpy(Nci*Nab,One,Work(I010),1,Work(Ioff),1)
          Call Putlst(Work(Ioff),1,Nci,1,Irrep_ci,List_h)
      Enddo
      Endif 

C OVVV block H(ai,bc)

      If (Do_hvvvo .or. Do_hovvv .or. Do_dvvvo .or. Do_dovvv) Then
      Ioff   = Ione
      List_h = 231 
      Do Irrep_bc = 1, Nirrep
         Irrep_ai = Dirprd(Irrep_bc,Irrepx)

         Nbc = Irpdpd(Irrep_bc,15)
         Nai = Irpdpd(Irrep_ai,11)
         I000 = Ione
         I010 = I000 + Nbc*Nai
         Iend = I010 + Nbc*Nai
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_2",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nbc*Nai)
         Icount = Inull
         Do Irrep_c = 1, Nirrep
            Irrep_b = Dirprd(Irrep_c,Irrep_bc)
            Do c = 1, Vrt(Irrep_c,Ispin)
               Do b = 1, Pop(Irrep_b,Ispin)
                  Do Irrep_i = 1, Nirrep
                     Irrep_a = Dirprd(Irrep_i,Irrep_ai)
                     Do i = 1, Pop(Irrep_i,Ispin)
                        Do a = 1, Vrt(Irrep_a,Ispin)
                           In = I + Ioffo(Irrep_i)
                           An = A + Ioffv(Irrep_a)
                           Bn = B + Ioffv(Irrep_b)
                           Cn = C + Ioffv(Irrep_c)
                           Work(Ioff+Icount) =
     +                                   H_pq(an,in)*D_pq(cn,bn)+
     +                                   H_pq(bn,cn)*D_pq(in,an)-
     +                                   H_pq(in,bn)*D_pq(cn,an)-
     +                                   H_pq(an,bn)*D_pq(cn,in)

#ifdef _NOSKIP
     +
     +                                   D_pq(an,bn)*H_pq(cn,in)+
     +                                   D_pq(in,cn)*H_pq(bn,an)-
     +                                   D_pq(in,bn)*H_pq(cn,an)-
     +                                   D_pq(an,bn)*H_pq(cn,in)
#endif 
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
#ifdef _DEBUG_LVL0
          Call checksum("OVVV    :",Work(Ioff),Nia*Nkj)
#endif 
          Call Getlst(Work(I010),1,Nbc,1,Irrep_bc,List_h)
          Call Daxpy(Nbc*Nai,One,Work(I010),1,Work(Ioff),1)
          Call Putlst(Work(Ioff),1,Nbc,1,Irrep_bc,List_h)
      Enddo
      Endif 
#endif 
#ifdef _DEBUG_LVL0
      Call Pccd_check_hess(Work,Maxcor)
#endif
      Return
      End 




