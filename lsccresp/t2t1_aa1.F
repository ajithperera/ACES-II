      Subroutine T2t1_aa1(T2,T1,Work,Maxcor,Pop,Vrt,Ispin,T2ln,T1ln,
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

      Data Izero,Ione,Dzero,Done,Onem /0,1,0.0D0,1.0D0,-1.0D0/

      Irrepx = Ione
      Listw  = 26 + Ispin
      T2off  = Ione + (Ispin-1)*T2ln_aa

      Do Irrep_an = 1, Nirrep

         Irrep_ef = Dirprd(Irrep_an,Irrepx) 

         Nsqdim = IZero 
         Irrep_in = Dirprd(Irrep_ef,Irrepx)
         Do Irrep_n = 1, Nirrep
            Irrep_i = Dirprd(Irrep_n,Irrep_in)
            Nsqdim  = Nsqdim + Pop(Irrep_n)*Pop(Irrep_i)
         Enddo 
       
         Nrow_t = Irpdpd(Irrep_ef,Ispin)
         Ncol_t = Irpdpd(Irrep_in,2+Ispin)
         Nrow_w = Irpdpd(Irrep_ef,Ispin)
         Ncol_w = Irpdpd(Irrep_an,8+Ispin)

         I000 = Ione
         I010 = I000 + Nrow_t*Nsqdim
         I020 = I010 + Nrow_w*Ncol_w
         I030 = I020 + Nrow_w
         Iend = I030 + Nrow_w
         If (Iend .Gt. Maxcor) Call Insmem("t2t1_aa1",Iend,Maxcor)
    
C T2(e<f,i<n)
         Call Dcopy(Nrow_t*Ncol_t,T2(T2off),1,Work(I000),1)

         T2off = T2off + Nrow_t*Ncol_t

C T2(e<f,in)
         Call Symexp(Irrep_in,Pop,Nrow_t,Work(I000))

C W(e<f,an)
         Call Getlst(Work(I010),1,Ncol_w,1,Irrep_an,Listw)

C W(e<f,an) -> W(e<f,na)
         Call Symtr1(Irrep_an,Vrt,Pop,Nrow_w,Work(I010),Work(I020),
     +               Work(Iend),Work(I030)) 

         Ioffs=  Ione + (Ispin-1)*T1ln_aa
         Ioffw = Izero 
         Iofft = Izero 
         
C -T(e<f,in)*W(e<f,na)-> -W(e<fn,a)*T(e<fn,i)

         Do Irrep_i = 1, Nirrep
            Irrep_n = Dirprd(Irrep_i,Irrep_in)
            Irrep_a = Dirprd(Irrep_n,Irrep_an) 
          
            Nsum = Nrow_w*Pop(Irrep_n)
            Nrow = Vrt(Irrep_a)
            Ncol = Pop(Irrep_i)
            Numn = Pop(Irrep_n)
   
            Lda = Nsum
            Ldb = Nrow_t*Pop(Irrep_n)
            Ldc = Nrow
 
            If (Min(Numn,Nrow,Ncol) .Ne. 0) Then
            Call Xgemm("T","N",Nrow,Ncol,Nsum,Onem,Work(I010+Ioffw),
     +                  Lda,Work(I000+iofft),Ldb,Done,T1(Ioffs),Ldc)
            Endif  

            Ioffs = Ioffs + Nrow*Ncol
            Ioffw = Ioffw + Nrow_w*Nrow*Numn
            Iofft = Iofft + Nrow_t*Ncol*Numn 

         Enddo 
      Enddo 
      
      Return
      End 
