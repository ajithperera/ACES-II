      Subroutine Get_energy(T2amps,Work,Maxcor,List1,List2,Listoff,
     +                      T2ln,T2ln_aa,T2ln_bb,Singles,Cct2,Nonhf,
     +                      Iuhf,E,ES)

      Implicit Double Precision(A-H,O-Z)
      Integer Dissiz_trn,Dissiz
      Integer T2ln,T2ln_aa,T2ln_bb,T2off
      Logical Singles,CCt2,Tau,Nonhf

#include "syminf.com"
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "sympop.com"
#include "info.com"

      Dimension Coo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension Cvv(Nfea(1)+Iuhf*Nfea(2))
      Dimension Work(Maxcor)
      Dimension T2amps(T2ln)
      Dimension E(3),ES(3)

      Data Ione,Inul,Two /1,0,2.0D0/

      Irrepx = Ione
      Listw0 = 13
      I000   =  Ione
      Tau    = .False.

      Call Dzero(E,3)
      Call Dzero(ES,3)

      If (Singles .and. Cct2) Then  
         Tau    = .True.
         Ndim_a = Irpdpd(1,9)
         It1a   = Ione 
         Iend   = It1a + Ndim_a
         Call Getlst(Work(It1a),1,1,1,1+Listoff,List1)
#ifdef _DEBUG_LVL0
         call checksum("T1_AA:",Work(It1a),Ndim_a)
#endif 

         If (Iuhf .Ne. 0) Then
           Ndim_b = Irpdpd(1,10)
           It1b   = Iend
           Iend   = It1b + Ndim_b
           Call Getlst(Work(It1b),1,1,1,2+Listoff,List1)
#ifdef _DEBUG_LVL0
         call checksum("T1_BB:",Work(It1b),Ndim_a)
#endif 
         Endif 

         List0 = 93
         If (Nonhf) then
             If1a = Iend
             Iend = If1a + Ndim_a
             Call Getlst(Work(If1a),1,1,1,3,List0)
 
             If (Iuhf .Ne. 0) Then
                If1b   = Iend
                Iend   = If1b + Ndim_b
                Call Getlst(Work(If1b),1,1,1,4,List0)
            Endif 
            E(1)  = Ddot(Ndim_b,Work(If1a),1,Work(It1a),1)
            ES(1) = E(1)
            If (Iuhf .Ne.0) Then
               E(2)  = Ddot(Ndim_b,Work(If1b),1,Work(It1b),1)
               ES(2) = E(2)
            Endif
         Endif 
      Else
         Iend = Ione 
         Tau  = .False.
      Endif 
#ifdef _DEBUG_LVL0
      Write(6,"(a,3F15.10)") "E1AA,E1BB       :", E(1),E(2)
#endif 
      
      Do Ispin =1, 1+Iuhf

         Ioff  = Ione
         Joff  = Ione
         T2off = (Ispin-1)*T2ln_aa + Ione

         Do Irrep = 1, Nirrep

            Irrepr = Irrep
            Irrepl = Dirprd(Irrepr,Irrepx)

            Numdis_trn = Irpdpd(Irrepr,2+Ispin)
            Dissiz_trn = Irpdpd(Irrepl,Ispin)
            
            I000  = Iend
            I010  = I000 + Numdis_trn*Dissiz_trn
            Ilast = I010 + Numdis_trn*Dissiz_trn
            If (Ilast .Gt. Maxcor) Call Insmem("get_energy",Iend,
     +                                          Maxcor)
C T2(A<B,I<J)
            If (Cct2) Then
                Call Getlst(Work(I000),1,Numdis_trn,1,Irrepr,
     +                      List2+Ispin)
            Else
                Call Dcopy(Dissiz_trn*Numdis_trn,T2amps(T2off),Ione,
     +                     Work(I000),Ione)
            Endif

            Listw = Listw0 + Ispin
            If (Ispin .Eq. 1) Then
            Call Energy_aa(Work(It1a),Work(It1a),Work(I000),Work(I010),
     +                     Ndim_a,Ndim_a,Numdis_trn,Dissiz_trn,Listw,
     +                     Irrep,Ispin,Tau,E,ES)
            Elseif (Ispin .Eq. 2) Then
            Call Energy_aa(Work(It1b),Work(It1b),Work(I000),Work(I010),
     +                     Ndim_b,Ndim_b,Numdis_trn,Dissiz_trn,Listw,
     +                     Irrep,Ispin,Tau,E,ES)
            Endif 
            T2off = T2off + Numdis_trn*Dissiz_trn
         Enddo
      Enddo 

      T2off = Ione + T2ln_aa + T2ln_bb

      Do Irrep = 1, Nirrep

         Irrepr = Irrep
         Irrepl = Dirprd(Irrepr,Irrepx)

          Numdis = Irpdpd(Irrepr,14)
          Dissiz = Irpdpd(Irrepl,13)

          I000  = Iend 
          I010  = I000 + Numdis*Dissiz
          Ilast = I010 + Numdis*Dissiz
          If (Ilast .Gt. Maxcor) Call Insmem("get_energy",Iend,
     +                                        Maxcor)
C T2(Ab,Ij)
          If (Cct2) Then
              Call Getlst(Work(I000),1,Numdis,1,Irrepr,List2+3)
          Else
              Call Dcopy(Numdis*Dissiz,T2amps(T2off),Ione,Work(I000),
     +                   Ione)
          Endif

          Listw = 16
          If (Iuhf .Eq. 0) Then
             Call Energy_ab(Work(It1a),Work(It1a),Work(I000),
     +                      Work(I010),Ndim_a,Ndim_a,Dissiz,
     +                      Numdis,Listw,Irrep,Tau,E,ES)
          Else
             Call Energy_ab(Work(It1a),Work(It1b),Work(I000),
     +                      Work(I010),Ndim_a,Ndim_b,Dissiz,
     +                      Numdis,Listw,Irrep,Tau,E,ES)
          Endif 

          T2off = T2off + Numdis*Dissiz
      Enddo

      If (Iuhf .Eq. 0) Then
          E(1)  = Two*E(1)
          ES(1) = Two*(1)
      Endif 

#ifdef _DEBUG_LVL0
      Write(6,"(a,3F15.10)") "E2AA,E2BB,E2AB  :", E(1),E(2),E(3)
#endif 
 
      Return
      End
 
