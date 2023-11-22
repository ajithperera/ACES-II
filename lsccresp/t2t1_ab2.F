      Subroutine T2t1_ab2(T2,T1,Work,Maxcor,Pop1,Pop2,Vrt1,Vrt2,
     +                    Ispin,Iuhf,T2ln,T1ln,T2ln_aa,T2ln_bb,
     +                    T2ln_ab,T1ln_aa,T1ln_bb)

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

      Data Izero,Ione,Done,Onem,Dzero/0,1,1.0D0,-1.0D0,0.0D0/
      Data Two/2.0D0/

      Irrepx = Ione
      Listw  = 11 - Ispin
      T2off  = Ione + T2ln_aa + T2ln_bb
    
      Do Irrep_if = 1, Nirrep
         Irrep_mn = Dirprd(Irrep_if,Irrepx) 
         Irrep_af = Dirprd(Irrep_mn,Irrepx)

         Nrow_t = Irpdpd(Irrep_af,15)
         Ncol_t = Irpdpd(Irrep_mn,14)
         Nrow_w = Irpdpd(Irrep_mn,14)
         If (Ispin .Eq. 1) Then 
            Ncol_w = Irpdpd(Irrep_if,18)
         Elseif (Ispin .Eq. 2) Then
            Ncol_w = Irpdpd(Irrep_if,11)
         Endif 

         I000 = Ione
         I010 = I000 + Nrow_t*Ncol_t
         I020 = I010 + Nrow_t*Ncol_t
         I030 = I020 + Nrow_w*Ncol_w
         I040 = I030 + Max(Nrow_w,Nrow_t,Ncol_t)
         Iend = I040 + Max(Nrow_w,Nrow_t,Ncol_t)
         If (Iend .Gt. Maxcor) Call Insmem("t2t1_ab2",Iend,Maxcor)

C T2(Af,Mn)
         Call Dcopy(Nrow_t*Ncol_t,T2(T2off),1,Work(I000),1)

         T2off = T2off + Nrow_t*Ncol_t

         If (Iuhf .Eq. 0) Then
            Call Symtra(Irrep_mn,Pop1,Pop2,Nrow_t,Work(I000),
     +                  Work(I010))
            Call Dscal(Nrow_t*Ncol_t,Two,Work(I000),1)
            Call Daxpy(Nrow_t*Ncol_t,Onem,Work(I010),1,Work(I000),1)
         Endif

C T2(Af,Mn) -> T2(Mn,Af)
         Call Transp(Work(I000),Work(I010),Ncol_t,Nrow_t)

C T2(Mn,Af)->T2(Mn,fA)
         If (Ispin .Eq. 1) Then
            Call Symtr1(Irrep_af,Vrt1,Vrt2,Ncol_t,Work(I010),
     +                  Work(I030),Work(Iend),Work(I040)) 
         Endif 
         
C W(Mn,If) Ispin=1, W(Mn,Fi) Ispin=2
         Call Getlst(Work(I020),1,Ncol_w,1,Irrep_if,Listw)

C W(Mn,If) -> W(Mn,fI)
         If (Ispin .Eq. 1) Then
            Call Symtr1(Irrep_If,Pop1,Vrt2,Nrow_w,Work(I020),
     +                  Work(I030),Work(Iend),Work(I040)) 
         Endif
          
         Ioffs = Ione + (Ispin-1)*T1ln_aa
         Ioffw = Izero 
         Iofft = Izero 

C T(Mn,fA)*W(Mn,fI)-> T(Mnf,A)*W(Mnf,I) (Ispin=1)
C T(Mn,Fa)*W(Mn,Fi)-> T(a,FMn)*W(MnF,i) (Ispin=2)

         Do Irrep_i = 1, Nirrep
            Irrep_f = Dirprd(Irrep_i,Irrep_if)
            Irrep_a = Dirprd(Irrep_f,Irrep_af) 

            Nsum = Nrow_w*Vrt2(Irrep_f)
            Nrow = Vrt1(Irrep_a)
            Ncol = Pop1(Irrep_i)
            Numn = Vrt2(Irrep_f)
   
            Lda = Ncol_t*Vrt2(Irrep_f)
            Ldb = Nsum
            Ldc = Nrow

            Call Xgemm("T","N",Nrow,Ncol,Nsum,Onem,Work(I010+Iofft),
     +                  Lda,Work(I020+ioffw),Ldb,Done,T1(Ioffs),Ldc)
         
            Ioffs = Ioffs + Nrow*Ncol
            Ioffw = Ioffw + Nrow_w*Ncol*Numn
            Iofft = Iofft + Ncol_t*Nrow*Numn 

         Enddo 
      Enddo 
      
      Return
      End 
