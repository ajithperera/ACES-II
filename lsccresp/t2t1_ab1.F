      Subroutine T2t1_ab1(T2,T1,Work,Maxcor,Pop1,Pop2,Vrt1,Vrt2,
     +                    Ispin,Iuhf,T2ln,T1ln,T2ln_aa,T2ln_bb,
     +                    T2lln_ab,T1ln_aa,T1ln_bb)

      Implicit Double Precision(A-H,O-Z)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"

      Integer Pop1(8),Vrt1(8),Pop2(8),Vrt2(8)
      Integer T1ln_aa,T1ln_bb
      Integer T2ln_aa,T2ln_bb,T2ln_ab
      Integer T2ln,T1ln
      Integer T2off 

      Dimension Work(Maxcor)
      Dimension T2(T2ln),T1(T1ln)

      Data Izero,Ione,Done,Onem,Two/0,1,1.0D0,-1.0D0,2.0D0/

      Irrepx = Ione
      Listw  = 31 - Ispin
      T2off  = Ione + T2ln_aa + T2ln_bb
    
      Do Irrep_an = 1, Nirrep
         Irrep_ef = Dirprd(Irrep_an,Irrepx) 
         Irrep_in = Irrep_ef

         Nrow_t = Irpdpd(Irrep_ef,15)
         Ncol_t = Irpdpd(Irrep_in,14)
         Nrow_w = Irpdpd(Irrep_ef,15)
         If (Ispin .Eq. 1) Then
            Ncol_w = Irpdpd(Irrep_an,11)
         Elseif (Ispin .Eq. 2) Then
            Ncol_w = Irpdpd(Irrep_an,18)
         Endif 

         I000 = Ione
         I010 = I000 + Nrow_t*Ncol_t
         I020 = I010 + Nrow_w*Ncol_w
         I030 = I020 + Nrow_t*Ncol_t
         Iend = I030 + Max(Nrow_w,Nrow_t)
         If (Iend .Gt. Maxcor) Call Insmem("t2t1_ab1",Iend,Maxcor)

C T2(Ef,In)
         Call Dcopy(Nrow_t*Ncol_t,T2(T2off),1,Work(I000),1)
         T2off = T2off + Nrow_t*Ncol_t

         If (Iuhf .Eq. 0) Then
            Call Symtra(Irrep_in,Pop1,Pop2,Nrow_t,Work(I000),
     +                  Work(I020))
            Call Dscal(Nrow_t*Ncol_t,Two,Work(I000),1)
            Call Daxpy(Nrow_t*Ncol_t,Onem,Work(I020),1,Work(I000),1)
          Endif 

C T2(Ef,In) -> T2(Ef,nI)

         If (Ispin .Eq. 1) Then
             Call Symtr1(Irrep_in,Pop1,Pop2,Nrow_t,Work(I000),
     +                   Work(I020),Work(Iend),Work(I030))
         Endif 

         Call Getlst(Work(I010),1,Ncol_w,1,Irrep_an,Listw)

C W(Ef,An) -> W(ef,nA) Ispin=1, W(Ef,Na) Ispin=2

         If (Ispin .Eq. 1) Then 
             Call Symtr1(Irrep_an,Vrt1,Pop2,Nrow_w,Work(I010),
     +                   Work(I020),Work(Iend),Work(I030)) 
         Endif 

         Ioffs = Ione + (Ispin-1)*T1ln_aa
         Ioffw = Izero 
         Iofft = Izero 

C T(Ef,nI)*W(Ef,nA)-> W(Efn,A)*T(Efn,I) Isipn=1
C T(Ef,Ni)*W(Ef,Na)-> W(EfN,a)*T(EfN,i) Isipn=2

         Do Irrep_i = 1, Nirrep
            Irrep_n = Dirprd(Irrep_i,Irrep_in)
            Irrep_a = Dirprd(Irrep_n,Irrep_an)

            Nsum = Nrow_w*Pop2(Irrep_n)
            Nrow = Vrt1(Irrep_a)
            Ncol = Pop1(Irrep_i)
            Numn = Pop2(Irrep_n)

            Lda = Nsum
            Ldb = Nrow_t*Pop2(Irrep_n)
            Ldc = Nrow

            Call Xgemm("T","N",Nrow,Ncol,Nsum,Done,Work(I010+Ioffw),
     +                  Lda,Work(I000+iofft),Ldb,Done,T1(Ioffs),Ldc)
            Ioffs = Ioffs + Nrow*Ncol
            Ioffw = Ioffw + Nrow_w*Nrow*Numn
            Iofft = Iofft + Nrow_t*Ncol*Numn

         Enddo
      Enddo

      Return
      End 
