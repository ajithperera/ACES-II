      Subroutine form_dbls_correctns_2rcc(Rcc_eigs,Rcc_vecl,Rcc_vecr,
     &                                    E2_rcc,Work,Length,Iuhf,
     &                                    Len_ov,Irrepx,Imult,
     &                                    Eom_sf_rccd,Eom_sf_drccd)

      Implicit Double Precision(A-H,O-Z)

      Dimension Rcc_vecl(Len_ov,Len_ov),Rcc_vecr(Len_ov,Len_ov)
      Dimension Rcc_eigs(Len_ov)
      Dimension Work(Length),E2_rcc(Len_ov,8)
      Dimension Nroot(8)
      Dimension Idoo(2),Idvv(2)
      Dimension Listijka(2),Listabci(2),Listmbej(3)
      Double Precision Ovlp 
      Logical Eom_sf_rccd,Eom_sf_drccd

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Common /Tdalist/ Listetda,Listvtda 
      Data Listc1,Listt2,Listt2in,Listt2rs /490,444,461,440/
      Data Listijka /107,7/
      Data Listabci /127,27/
      Data Listd0   /447/
      Data Listmbej /54,54,254/
      Data Scale,One,Onem,Dnull/1.0D0,1.0D0,-1.0D0,0.0D0/

      Call Getrec(20,'JOBARC','EESYMINF',NIRREP,NRoot)

#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(a,2(1x,i2),1x,2l)") " IRREPX,IMULT and EOM_SF_D/RCCD:", 
     &                    Irrepx,Imult,Eom_sf_rccd,Eom_sf_drccd
CSSS      call check_ints(Work,Length,Iuhf,0)
      write(6,*)
      lengtha =irpdpd(irrepx,9)
      lengthb =irpdpd(irrepx,10)
      Do Iroot = 1, Nroot(Irrepx)
      joff = 1
      do ispin = 1, 1+iuhf 
      joff = joff+lengtha*(ispin-1)
      if (ispin.eq.1) then
      call output(Rcc_vecr(joff,iroot),1,Lengtha,1,1,Lengtha,1,1)
      call output(Rcc_vecl(joff,iroot),1,Lengtha,1,1,Lengtha,1,1)
      write(6,*)
      call checksum("Rcc_R(P,Q)",Rcc_vecr(joff,Iroot),Lengtha,s)
      call checksum("Rcc_L(P,Q)",Rcc_vecl(joff,Iroot),Lengtha,s)
      ovlp = Ddot(Lengtha,Rcc_vecr(joff,Iroot),1,
     +                    Rcc_vecl(joff,Iroot),1)
      else 
      call output(Rcc_vecr(joff,iroot),1,Lengthb,1,1,Lengthb,1,1)
      call output(Rcc_vecl(joff,iroot),1,Lengthb,1,1,Lengthb,1,1)
      write(6,*)
      ovlp = Ovlp+ Ddot(Lengtha,Rcc_vecr(joff,Iroot),1,
     +                  Rcc_vecl(joff,Iroot),1)
      call checksum("Rcc_R(P,Q)",Rcc_vecr(joff,Iroot),Lengthb,s)
      call checksum("Rcc_L(P,Q)",Rcc_vecl(joff,Iroot),Lengthb,s)
      endif 
      Write(6,"(a,1x,F5.2)") "<L|R> norm: ", Ovlp 
      enddo 
      enddo 
      write(6,*)
#endif 

      Do My_root = 1, Nroot(Irrepx)
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "--------Forming R1*W->R2*L1*W=E(D)----------"
#endif 
         Ioff = 1
         Call Form_d2_4rcc(Work,Rcc_eigs,Len_ov,Length,Iuhf,Irrepx,
     +                     Listd0,1)
         Do Ispin = 3,3-2*Iuhf,-1  
            Call Zerolist(Work,Length,Listt2in-1+Ispin)
         Enddo 
         Call Zero_490_lists(Work,Length,Irrepx,Iuhf)
         E2_rcc_aa = 0.0D0
         E2_rcc_bb = 0.0D0

         Do Ispin = 1, (Iuhf+1)
      
C Lets write the RPA vectors to the same lists that we used to
C store the CIS vectors. That will make the rest of the code's 
C work flow easier.

            Len_ph_aa = Irpdpd(Irrepx,9)
            Len_ph_bb = Irpdpd(Irrepx,10)
            Len_ph_pq = Len_ph_aa + Len_ph_bb * Iuhf

            Ioff = Ioff + (Ispin - 1) * Len_ph_aa

C This was done to compare with the UHF for debugging purposes. 
C            If (iuhf .eq. 0) Call Dscal(Len_ph_aa,1.0D0/Dsqrt(2.0D0),
C     +                                  Rcc_vecr(ioff,My_root),1)

            Call Putlst(Rcc_vecr(ioff,My_root),1,1,1,Ispin,Listc1) 
               
         Enddo 

C The doubles correction to RPA requires doubles target lists
C (lists 461 have already been formed). Therefore, we can proceed
C to built the first part of what is required for doubles correction,

         Listz0 = Listt2in 
         Listw1 = Listijka(2)
         Listw2 = Listabci(2)

         Call C1rpaint2c2a(Work,Length,Iuhf,Irrepx,Listw1,Listz0,
     &                     Listc1,Imult)
         Call C1rpaint2c2b(Work,Length,Iuhf,Irrepx,Listw2,Listz0,
     &                     Listc1,Imult)
#ifdef _DEBUG_LVL0
         Write(6,"(a)") "W(ijka)*R1+W(abci)*R1-> W2 (in list 463)"
         Call Check_vee_lists(Work,Length,Irrepx,Iuhf,1)
         Write(6,*)
#endif 
         Call form_new_c2(Work,Length,Iuhf,Irrepx,Listz0,Listd0)

#ifdef _DEBUG_LVL0
         Write(6,"(a)")"D2(W(ijka)*R1+W(abci)*R1)-> R2 (in list 463)"
         Call Check_vee_lists(Work,Length,Irrepx,Iuhf,1)
         Write(6,*)
#endif 
         Call C2inc1a(Work,Length,Iuhf,Irrepx,Listw2,Listz0,
     &                Listc1,Imult)
         Call C2inc1b(Work,Length,Iuhf,Irrepx,Listw1,Listz0,
     &                Listc1,Imult)
#ifdef _DEBUG_LVL0
         Write(6,"(a)")"R2*W(ab,ci)+R2*W(ika))-> R1 (in list 490)"
         Call Check_vee_lists(Work,Length,Irrepx,Iuhf,1)
#endif 
         Call Dbls_correctns_2rpa(Work,Length,Iuhf,RCC_vecl(1,My_root),
     &                            Len_ph_aa,Len_ph_bb,Len_ph_pq,
     &                            E2_rcc_aa,E2_rcc_bb,Listc1)
#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(a,F15.9)") "E2_rcc_aa: ", E2_rcc_aa*Scale
      if (Iuhf .eq. 1) write(6,"(a,F15.9)") "E2_rcc_bb: ", E2_rcc_bb*
     &                                       Scale
#endif 
         E2_rcc(My_root,Irrepx) = (E2_rcc_aa + E2_rcc_bb)*Scale 
C
#ifdef _CIS
         Call Zero_490_lists(Work,Length,Irrepx,Iuhf)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "--------Forming R1*W(mb,ej)*T*L1*W=E(D)----------"
#endif
C
C This term applies only to CIS. For RPA and DRPA this term is alrready
C in HBAR

         Len_ov = Irpdpd(Irrepx,9) + Irpdpd(Irrepx,10)*Iuhf
         It2c_vec = 1
         Iw2c_vec = It2c_vec + Len_ov
         Iend     = Iw2c_vec + Len_ov
         Length   =  Length - Iend
         If (Iend .Ge. Length) Call Insmem("dbls_correction_2rcc",
     &                                      Iend,Length)

         Call Resort(Work(Iend),Length,Iuhf,1,44,34)
         Call Gt2xf(Work(It2c_vec),1,Irrepx,34,Listc1,0,Work(Iend),
     &              Length,Iuhf,0)
         Call Resort(Work(Iend),Length,Iuhf,1,14,34)
         Call Gt2xf(Work(Iw2c_vec),1,Irrepx,34,Listc1,0,Work(Iend),
     &              Length,Iuhf,0)

         Do Ispin = 1, Iuhf+1
            Length_aa = Irpdpd(Irrepx,9)
            Length_bb = Irpdpd(Irrepx,10)
            Itoff_aa  = It2c_vec
            Iwoff_aa  = Iw2c_vec
            Itoff_bb  = It2c_vec + Length_aa
            Iwoff_bb  = Iw2c_vec + Length_aa

            if (Ispin .Eq. 1) E2_rcc_aa = Ddot(Length_aa,
     &                                         Work(Itoff_aa),
     &                                         1,Work(Iwoff_aa),1)
            if (Ispin .Eq. 2) E2_rcc_bb = Ddot(Length_bb,
     &                                         Work(Itoff_bb),
     &                                         1,Work(Iwoff_bb),1)
         Enddo

#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(a,F15.9)") "E2_rcc_aa: ", E2_rcc_aa
      if (Iuhf .eq. 1) write(6,"(a,F15.9)") "E2_rcc_bb: ", E2_rcc_bb
#endif
         E2_rcc(My_root,Irrepx) = E2_rcc(My_root,Irrepx) + E2_rcc_aa +
     &                            E2_rcc_bb
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "--------Forming R1*W(mn,ef)*T*L1*W=E(D)----------"
#endif
         If (Eom_sf_rccd .OR. Eom_sf_drccd) Then
         E2_rcc_oo = Dnull
         E2_rcc_vv = Dnull
         Itda_vec_off = 1
         Do Ispin = 1, Iuhf+1
            Length_ai = Irpdpd(Irrepx,9)
            Idoo(1) = Iend
            Idoo(2) = Idoo(1) + Nfmi(1)
            Idvv(1) = Idoo(2) + Iuhf*Nfmi(2)
            Idvv(2) = Idvv(1) + Nfea(1)
            Iend    = Idvv(2) + Iuhf*Nfea(1)
            Max_scr = Max(Nfmi(1),Nfmi(2),Nfea(1),Nfea(2))
            Iend    = Iend + Max_scr
            Length  = Length - Iend
            If (Iend .Ge. Length) Call Insmem("dbls_correction",
     &                                         Iend,Length)

            Call Dzero(Work(Iend),Nfmi(Ispin))
            Call Modgij(1,Irrepx,Irrepx,Work(Idoo(Ispin)),Work(Iend),
     &                  Work(Itda_vec_off),Work(Itda_vec_off),Ispin,
     &                  Onem)
            Call Dzero(Work(Iend),Nfea(Ispin))
            Call Modgab(1,Irrepx,Irrepx,Work(Idvv(Ispin)),Work(Iend),
     &                  Work(Itda_vec_off),Work(Itda_vec_off),Ispin,
     &                  One)
            Call Gformg(1,1,44,14,0,Work(Iend),Length,0,One,One,Iuhf)

            Length_oo = Nfmi(Ispin)
            Length_vv = Nfea(Ispin)
            Call Getlst(Work(Iend),1,1,1,Ispin,91)
            E2_rcc_oo = E2_rcc_oo + Ddot(Length_oo,Work(Iend),1,
     &                                   Work(Idoo(Ispin)),1)
            Call Getlst(Work(Iend),1,1,1,Ispin,92)
            E2_rcc_vv = E2_rcc_vv + Ddot(Length_vv,Work(Iend),1,
     &                                   Work(Idvv(Ispin)),1)
            Itda_vec_off = Itda_vec_off +  Length_ai
         Enddo

#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(a,2(1x,F15.9))") "E2_rcc_oo and E2_rcc_vv: ",
     &                            E2_rcc_oo,E2_rcc_vv
#endif
         E2_rcc(My_root,Irrepx) = E2_rcc(My_root,Irrepx) + E2_rcc_oo +
     &                            E2_rcc_vv
         Enddo
#endif 
      Enddo     

#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,"(a)") "EOM(S/F)-RCC(D) correction"
      Do My_root = 1, Nroot(Irrepx)
         write(6,"(a,F15.9)") "E2_rcc   : ", E2_rcc(My_root,Irrepx)
      Enddo 
#endif 
      Return
      End 

      
