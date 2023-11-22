#include "flags.h"
      Subroutine Rcc_make_d2(E,D,Iwork,Nbasis,Isize,Laabb,Iuhf)

      Implicit Double Precision (A-H, O-Z)
      Dimension E(2*Nbasis),D(Laabb)
      Dimension Iwork(Isize)
      Integer A,B,A1

      Logical Rle
      COMMON/EXTRAPO/RLE

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "flags.com"
#include "info.com"

      NNP1O2(N)=(N*(N+1))/2
      NNM1O2(N)=(N*(N-1))/2 
      ILRG(IX) = INT(0.5*(1.0D0+SQRT(8.0D0*IX-7))+1.0D-08)+1 

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,"(a)") "--------Entering rcc_make_d2---------"
      Write(6,*)
#endif

      Call Getrec(20,"JOBARC","SCFEVALA",Nbasis*Iintfp,E)
      Call Getrec(20,"JOBARC","SCFEVALB",Nbasis*Iintfp,E(Nbasis+1))

      Irrepx = 1
      Call Newtyp2(Irrepx,48,5,7,.True.)
      Call Newtyp2(Irrepx,64,5,7,.True.)
      If (Iuhf .Ne. 0) Call Newtyp2(Irrepx,49,6,8,.True.)
      If (Iuhf .Ne. 0) Call Newtyp2(Irrepx,65,6,8,.True.)
      
      Do Ispin = 1, Iuhf+1
         Ioff = (Ispin-1) * Nbasis 
         Nrow = Nnp1o2(Nvrto(Ispin))
         Ncol = Nnp1o2(Nocco(Ispin))
         I000 = 1
         I010 = I000 + Nirrep
         I020 = I010 + Ncol 
         I030 = I020 + Nirrep
         I040 = I030 + Nrow 

         If (Ispin .EQ. 1) Then
            Call Getrec(20,"JOBARC",'SOAOA1X ',Nirrep,Iwork(I000))
            Call Getrec(20,"JOBARC",'SOAOA1  ',Ncol,  Iwork(I010))
            Call Getrec(20,"JOBARC",'SVAVA1X ',Nirrep,Iwork(I020))
            Call Getrec(20,"JOBARC",'SVAVA1  ',Nrow,  Iwork(I030))
         Else
            Call Getrec(20,"JOBARC",'SOBOB1X ',Nirrep,Iwork(I000))
            Call Getrec(20,"JOBARC",'SOBOB1  ',Ncol,  Iwork(I010))
            Call Getrec(20,"JOBARC",'SVBVB1X ',Nirrep,Iwork(I020))
            Call Getrec(20,"JOBARC",'SVBVB1  ',Nrow,  Iwork(I030))
         Endif 

#ifdef _DEBUG_LVLM
        Write(6,"(a)") " The symmetry arrays"
        Write(6,"(6(1x,I4))") (Iwork(I000-1+i),i=1,Nirrep)
        Write(6,*) 
        Write(6,"(8(1x,I4))") (Iwork(I010-1+i),i=1,Ncol)
        Write(6,*) 
        Write(6,"(8(1x,I4))") (Iwork(I020-1+i),i=1,Nirrep)
        Write(6,*) 
        Write(6,"(8(1x,I4))") (Iwork(I030-1+i),i=1,Nrow)
#endif 
          Ioffl = 0
          Ioffr = 0
          Do Irrep = 1, Nirrep
             Noo    = Iwork(I000+Irrep-1)
             Nvv    = Iwork(I020+Irrep-1)
             Icount = 0

             Do Ioo = 1, Noo
                Ioo_loc = Ioo + Ioffr
                Ioo_off = Iwork(I010+Ioo_loc-1) 
                I       = Ilrg(Ioo_off) - 1
                I1      = Ilrg(Ioo_off)
                J       = Ioo_off - Nnm1O2(I1-1) 
#ifdef _DEBUG_LVLM
                Write(6,"(a,3(1x,i4))") "I,J: ",i,j
                write(6,*)
#endif 
                Do Ivv = 1, Nvv
                   Icount = Icount + 1
                   Ivv_loc = Ivv + Ioffl
                   Ivv_off = Iwork(I030+Ivv_loc-1) 
                   A       = Nocco(Ispin) + Ilrg(Ivv_off) - 1
                   A1      = Ilrg(Ivv_off)
                   B       = Nocco(Ispin) + Ivv_off - Nnm1O2(A1-1) 
#ifdef _DEBUG_LVLM
                   Write(6,"(a,4(1x,i4))") "I,J,A,B: ",I,J,A,B
#endif 
                   D(Icount) = 1.0D0/(E(Ioff+I)+E(Ioff+J)-
     +                                E(Ioff+A)-E(Ioff+B))
                Enddo
             Enddo 
             Call Putlst(D,1,Noo,1,Irrep,47+ispin)

#ifdef _DEBUG_LVLM
       Write(6,"(a,a,i2,1x,I2)") " The denominator arrays for irrep",
     +                           " and spin: ",irrep,ispin
          Call output(D,1,Nvv,1,Noo,Nvv,Noo,1)
#endif 
             If (Rle) Then
                Do Iinv = 1, Icount
                    D(Iinv) = 1.0D0/D(Iinv)
                Enddo 
                Call Putlst(D,1,Noo,1,Irrep,63+ispin)
             Endif 

             Ioffr = Ioffr + Noo
             Ioffl = Ioffl + Nvv
          Enddo
      Enddo

      Return
      end
