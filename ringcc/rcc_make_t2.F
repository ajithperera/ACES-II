#include "flags.h"
      Subroutine Rcc_make_t2(Work,Maxcor,Iuhf)

      Implicit Double Precision (A-H, O-Z)
      Dimension Work(Maxcor)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "sym.com"
#include "parcc.com"
#include "flags.com"
#include "info.com"

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "--------Entering rcc_make_t2---------"
      Write(6,*)
#endif 

      Irrepx = 1
      Imod   = 0

C Form lists to store <AB|IJ> integrals (ordered as AB;IJ)
C List 1 and 2.

      Call Inipck(Irrepx,19,21,1,Imod,0,1)
      If (Iuhf .Ne.0) Call Inipck(Irrepx,20,22,2,Imod,0,1)

       Eaa = 0.0D0
       Ebb = 0.0D0
       Do Ispin = 1, 1+Iuhf 
       
          Nsize = Idsymsz(Irrepx,Isytyp(1,4+Ispin),Isytyp(2,4+Ispin))
          Nrows = Nvrto(Ispin) * Nvrto(Ispin) +
     +            Nocco(Ispin) * Nocco(Ispin) +
     +            Nvrto(Ispin) * Nocco(Ispin) 

          I000 = 1
          I010 = I000 + Nsize
          I020 = I010 + Nsize 
          I030 = I020 + Nrows 
          Iend = I030
          If (Iend .Ge. Maxcor) Call Insmem("rcc_make_t2",Iend,
     +                                       Maxcor)
C Read in Coulomb integrals as Ispin=1: W(AI,BJ) Ispin=2 (ai,bj)

          Call Getall(Work(I010),Nsize,Irrepx,4+Ispin)

C Form the integrals (AB,IJ)by permuting I and B;Ispin=1: W(AB,IJ)
C Ispin=2 (ab,ij)

          Call Sstgen(Work(I010),Work(I000),Nsize,Vrt(1,Ispin),
     +                Pop(1,Ispin),Vrt(1,Ispin),Pop(1,Ispin),
     +                Work(I020),Irrepx,"1324")

          Call Putall(Work(I000),Nsize,Irrepx,Ispin)

          Do Irrepr = 1, Nirrep
             Irrepl = Dirprd(Irrepr,Irrepx)

             Nrow_tring_ab = Irpdpd(Irrepl,4+Ispin)
             Ncol_tring_ij = Irpdpd(Irrepr,6+Ispin)
             Nrow_squar_ab = Irpdpd(Irrepl,18+Ispin)
             Ncol_squar_ij = Irpdpd(Irrepr,20+Ispin)

             I000 = 1
             I010 = I000 + Nrow_tring_ab * Ncol_tring_ij
             I020 = I010 + Nrow_squar_ab * Ncol_squar_ij
             I030 = I020 + Nrow_squar_ab * Ncol_squar_ij
             I040 = I030 + Nrow_squar_ab * Ncol_squar_ij
             Iend = I040

             If (Iend .Ge. Maxcor) Call Insmem("rcc_make_d2",Iend,
     +                                         Maxcor)
C Retrive D(A<=B,I<=J) and expand,

             Call Getlst(Work(I000),1,Ncol_tring_ij,1,Irrepr,47+Ispin)

C D(A<=B,I<=J) -->D(AB,I<=J)
             
             Call Symexp6(Irrepl,Vrt(1,ispin),Vrt(1,ispin),
     +                    Nrow_squar_ab,Nrow_tring_ab,
     +                    ncol_tring_ij,work(I010),work(I000),
     +                    work(I030))
             Call Transp(work(I010),work(I020),ncol_tring_ij,
     +                        nrow_squar_ab)
C D(IJ,A<=B) -->D(iJ,AB)

             Call Symexp6(Irrepr,Pop(1,ispin),Pop(1,ispin),
     +                    Ncol_squar_ij,Ncol_tring_ij,
     +                    nrow_squar_ab,work(I010),work(I020),
     +                    work(i030))
C D(IJ,AB) -->D(AB,IJ)

             Call Transp(work(I010),work(I020),Nrow_squar_ab,
     +                        Ncol_squar_ij)
#ifdef _DEBUG_LVLM
      Write(6,*) "D<AB|IJ>"
      call output(work(I020),1,Nrow_squar_ab,1,Ncol_squar_ij,
     +            Nrow_squar_ab,Ncol_squar_ij,1)
#endif 
C Retrive integrals I(AB;IJ) array

             Call Getlst(Work(I010),1,Ncol_squar_ij,1,Irrepr,Ispin)

#ifdef _DEBUG_LVLM
      Write(6,*) "I<AB|IJ>"
      call output(work(I010),1,Nrow_squar_ab,1,Ncol_squar_ij,
     +            Nrow_squar_ab,Ncol_squar_ij,1)
#endif
             Do I = 1,Nrow_squar_ab * Ncol_squar_ij
                Work(I030+I-1)  = Work(I020+I-1) * Work(I010+I-1)
             Enddo

             Call Putlst(Work(I030),1,Ncol_squar_ij,1,Irrepr,
     +                   43+Ispin)

             If (Ispin .EQ. 1) Then
             Eaa = Eaa + Ddot(Nrow_squar_ab*Ncol_squar_ij,Work(I030),1,
     +                        Work(I020),1)
             Else
             Ebb = Ebb + Ddot(Nrow_squar_ab*Ncol_squar_ij,Work(I030),1,
     +                        Work(I020),1)
             Endif

#ifdef _DEBUG_LVLM
      write(6,*) " Coulomb only T2 amplitudes and corresponding energy"
      call output(Work(I030),1,Nrow_squar_ab,1,Ncol_squar_ij,
     +            Nrow_squar_ab,Ncol_squar_ij,1)
      call checksum("T2(AA)",Work(I030), Nrow_squar_ab*Ncol_squar_ij)
      If (Ispin .Eq. 1) Write(6,"(a,1x,F15.7)")
     +                                " Coulomb only energy: ", Eaa
      If (Ispin .Eq. 2) Write(6,"(a,1x,F15.7)")
     +                                " Coulomb only energy: ", Ebb
#endif
          Enddo
      Enddo

      Return
      End
c-------------------------------------------------------------------
#ifdef _NOSKIP
C This block is not used. I did not see at the begning (should
C have) symmetric amplitude means T2(ab,ij) means T2(ba,ji). 
C The following assume T2(ab,ij) = T2(ba,ij). Kind if dumb. 

          Do Irrepr = 1, Nirrep
             Irrepl = Dirprd(Irrepr,Irrepx)

             Nrow_tring_ab = Irpdpd(Irrepl,4+Ispin)
             Ncol_tring_ij = Irpdpd(Irrepr,6+Ispin)
             Nrow_squar_ab = Irpdpd(Irrepl,18+Ispin)
             Ncol_squar_ij = Irpdpd(Irrepr,20+Ispin)

             I000 = 1
             I010 = I000 + Nrow_squar_ab * Ncol_squar_ij
             I020 = I010 + Nrow_squar_ab * Ncol_tring_ij
             I030 = I020 + Nrow_tring_ab * Ncol_tring_ij
             I040 = I030 + Nrow_tring_ab * Ncol_tring_ij
             Iend = I040

             If (Iend .Ge. Maxcor) Call Insmem("rcc_make_d2",Iend,
     +                                         Maxcor)
             Call Getlst(Work(I000),1,Ncol_squar_ij,1,Irrepr,Ispin)

#ifdef _DEBUG_LVL0
      write(6,"(a,i2,1x,i2)") "I<AB,IJ> for irrep and spin: ", 
     +                         Irrepr, Ispin
      call output(work(i000),1,Nrow_squar_ab,1,Ncol_squar_ij,
     +            Nrow_squar_ab,Ncol_squar_ij,1)
#endif 
C Convert <AB,IJ> --> <AB;I<=J>

             Call Sq3sym(Irrepr,Pop(1,Ispin),Nrow_squar_ab,
     +                   Work(I010),Work(I000))
#ifdef _DEBUG_LVL0
      Write(6,*) "I<AB,I<=J>"
      call output(work(i010),1,Nrow_squar_ab,1,Ncol_tring_ij,
     +            Nrow_squar_ab,Ncol_tring_ij,1)
#endif 

C Convert <AB,I<=J> --> <A<=B;I<=J>

             Call Sq4sym(Irrepl,Vrt(1,Ispin),Nrow_tring_ab,
     +                   Nrow_squar_ab,Ncol_tring_ij,
     +                   Work(I020),Work(I010))

#ifdef _DEBUG_LVL0
      Write(6,*) "I<A<=B,I<=J>"
      call output(work(i020),1,Nrow_tring_ab,1,Ncol_tring_ij,
     +            Nrow_tring_ab,Ncol_tring_ij,1)
#endif 
C Store  <A<=B;I<=J> integral for future (energy calculations might
C need this)
             Call Putlst(Work(I020),1,Ncol_tring_ij,1,Irrepr,2+Ispin)

C Obtain the D(A<=B;I<=J) denominator array

             Call Getlst(Work(I000),1,Ncol_tring_ij,1,Irrepr,47+Ispin)

#ifdef _DEBUG_LVL0
      Write(6,*) "D<A<=B,I<=J>"
      call output(work(i000),1,Nrow_tring_ab,1,Ncol_tring_ij,
     +            Nrow_tring_ab,Ncol_tring_ij,1)
#endif 
C Store  <A<=B;I<=J> integral for future (energy calculations might
C Construct the Coulomb only T2(A<=B,I<=j) AAAA and BBBB

             Do I = 1,Nrow_tring_ab * Ncol_tring_ij
                Work(I030+I-1)  = Work(I020+I-1) * Work(I000+I-1)
             Enddo 
              
             Call Putlst(Work(I030),1,Ncol_tring_ij,1,Irrepr,
     +                   43+Ispin)

             If (Ispin .EQ. 1) Then
             Eaa = Eaa + Ddot(Nrow_tring_ab*Ncol_tring_ij,Work(I030),1,
     +                        Work(I020),1)
             Else 
             Ebb = Ebb + Ddot(Nrow_tring_ab*Ncol_tring_ij,Work(I030),1,
     +                        Work(I020),1)
             Endif 
#ifdef _DEBUG_LVLM
      write(6,*) " Coulomb only T2 amplitudes and corresponding energy"
      call output(Work(I030),1,Nrow_tring_ab,1,Ncol_tring_ij,
     +            Nrow_tring_ab,Ncol_tring_ij,1)
      If (Ispin .Eq. 1) Write(6,"(a,1x,F15.7)") 
     +                                " Coulomb only energy: ", Eaa
      If (Ispin .Eq. 2) Write(6,"(a,1x,F15.7)") 
     +                                " Coulomb only energy: ", Ebb
#endif 
          Enddo
      Enddo
#ifdef _DEBUG_LVL0
      write(6,"(a)") " Debugging symexp6 "
      do irrepr=1,nirrep
         irrepl=dirprd(irrepr,1)
         nrow_sq=irpdpd(irrepl,19)
         ncol_sq=irpdpd(irrepr,21)
         nrow_tr=irpdpd(irrepl,5)
         ncol_tr=irpdpd(irrepr,7)

         i000=1
         i010=i000+nrow_tr*ncol_tr
         i020=i010+nrow_sq*ncol_sq
         i030=i020+nrow_sq*ncol_sq
         i040=i030+nrow_sq*ncol_sq
         call getlst(work(i000),1,Ncol_tr,1,Irrepr,48)

         call output(work(i000),1,nrow_tr,1,ncol_tr,nrow_tr,
     +               ncol_tr,1)
         ispin=1
         call symexp6(irrepl,vrt(1,ispin),vrt(1,ispin),nrow_sq,
     +                nrow_tr,ncol_tr,work(i020),work(i000),
     +                work(i030))
         call transp(work(i020),work(i030),ncol_tr,nrow_sq)
         call symexp6(irrepr,pop(1,ispin),pop(1,ispin),ncol_sq,
     +                ncol_tr,nrow_sq,work(i010),work(i030),
     +                work(i020))
         call transp(work(i010),work(i020),nrow_sq,ncol_sq)
         write(6,"(a)") " D(AB,IJ)"
         call output(work(i020),1,nrow_sq,1,ncol_sq,nrow_sq,
     +               ncol_sq,1)
      enddo
#endif 
#endif 

