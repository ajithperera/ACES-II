      Subroutine T2t1_aa2(T2,T1,Work,Maxcor,Pop,Vrt,Ispin,T2ln,T1ln,
     +                    T2ln_aa,T2ln_bb,T2ln_ab,T1ln_aa,T1ln_bb)

      Implicit Double Precision(A-H,O-Z)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"

      Integer Pop(8),Vrt(8)
      Integer T1ln_aa,T1ln_bb
      Integer T2ln_aa,T2ln_bb,T2ln_ab
      Integer T2ln,T1ln
      Integer T2off 

      Dimension Work(Maxcor)
      Dimension T2(T2ln),T1(T1ln)

      Data Izero,Ione,Done,Dzero/0,1,1.0D0,0.0D0/

      Irrepx = Ione
      Listw  = 6 + Ispin
      T2off  = Ione + (Ispin-1)*T2ln_aa
    
      Do Irrep_if = 1, Nirrep
         Irrep_mn = Dirprd(Irrep_if,Irrepx) 
         Irrep_af = Irrep_mn

         Nsqdim = Izero 
         Do Irrep_f = 1, Nirrep
            Irrep_a = Dirprd(Irrep_f,Irrep_af)
            Nsqdim = Nsqdim + Vrt(Irrep_a)*Vrt(Irrep_f)
         Enddo 
       
         Nrow_t = Irpdpd(Irrep_af,Ispin)
         Ncol_t = Irpdpd(Irrep_mn,2+Ispin)
         Nrow_w = Irpdpd(Irrep_mn,2+Ispin)
         Ncol_w = Irpdpd(Irrep_if,8+Ispin)

         I000 = Ione
         I010 = I000 + Nrow_t*Ncol_t
         I020 = I010 + Nrow_t*Nsqdim
         I030 = I020 + Nrow_w*Ncol_w
         I040 = I030 + Max(Nrow_w,Ncol_w,Nrow_t,Ncol_t)
         Iend = I040 + Max(Nrow_w,Ncol_w,Nrow_t,Ncol_t)

         If (Iend .Gt. Maxcor) Call Insmem("t2t1_aa2",Iend,
     +                                       Maxcor)
C T2(a<f,m<n)
         Call Dcopy(Nrow_t*Ncol_t,T2(T2off),1,Work(I000),1)

         T2off = T2off + Nrow_t*Ncol_t

C T2(a<f,m<n) -> T2(m<n,af)
         Call Transp(Work(I000),Work(I010),Ncol_t,Nrow_t)
         Call Symexp(Irrep_af,Vrt,Ncol_t,Work(I010))

C W(m<n,if)
         Call Getlst(Work(I020),1,Ncol_w,1,Irrep_if,Listw)

C W(m<n,if) -> W(m<n,fi)
         Call Symtr1(Irrep_if,Pop,Vrt,Nrow_w,Work(I020),Work(I030),
     +               Work(Iend),Work(I040))

         Ioffs = Ione + (Ispin-1)*T1ln_aa
         Ioffw = Izero 
         Iofft = Izero 
          
C T(m<n,af)*W(m<n,fi)-> T(m<nf,a)*W(m<nf,i)

         Do Irrep_i = 1, Nirrep
            Irrep_f = Dirprd(Irrep_i,Irrep_if)
            Irrep_a = Dirprd(Irrep_f,Irrep_af) 

            Nsum = Nrow_w*Vrt(Irrep_f)
            Nrow = Vrt(Irrep_a)
            Ncol = Pop(Irrep_i)
            Numn = Vrt(Irrep_f)
   
            Lda = Ncol_t*Vrt(Irrep_f)
            Ldb = Nrow_w*Vrt(Irrep_f)
            Ldc = Nrow

            If (Min(Nrow,Ncol,Numn) .Ne. 0) Then
            Call Xgemm("T","N",Nrow,Ncol,Nsum,Done,Work(I010+Iofft),
     +                  Lda,Work(I020+Ioffw),Ldb,Done,T1(Ioffs),Ldc)
            Endif 

            Ioffs = Ioffs + Nrow*Ncol
            Ioffw = Ioffw + Nrow_w*Ncol*Numn
            Iofft = Iofft + Ncol_t*Nrow*Numn 

         Enddo 
      Enddo 
      
      Return
      End 
