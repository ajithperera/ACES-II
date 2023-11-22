#include "flags.h"
      Subroutine Prep_scfmos(Work,Maxcor,Coo,Cvv,Iuhf)

      Implicit Double Precision(A-H,O-Z)
      Logical Dropmo

#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "flags.com"
#include "flags2.com"
#include "info.com"

      Dimension Work(Maxcor)
      Dimension Coo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension Cvv(Nfea(1)+Iuhf*Nfea(2))
 
      Data Ione,Inull /1,0/

      Call Getrec(20,"JOBARC","NAOBASFN",Ione,Naobfns)
      Call Getrec(20,"JOBARC","NBASTOT ",Ione,Nbfns_tot)
      Call Getrec(0,"JOBARC","NUMDROPA",Length,Ndrop)
  
      Dropmo = (Length .Gt. Inull)
         
      If (Dropmo) Then 
         Call Getrec(20,"JOBARC","NUMDROPA",Ione,Ndrop)
         If (Iuhf .Ne. 0) Call Getrec(20,"JOBARC","NUMDROPB",Ione,
     +                                Ndrop)
      Endif 

      Nbfns = Nbfns_tot - Ndrop 
      If (Nbfns .Ne. Nocco(1)+nvrto(2)) Then
          Write(6,"(2a)") " Internal inconsistency in number of basis",
     +                    " functions."
          Write(6,"(a,I3,a,I3)") " Read: ",Nbfns, " Internally: ", 
     +                            Nocco(1)+nvrto(1)
      Endif 
 
      Nbfns2 = Nbfns_tot*Nbfns
      I000 = Ione
      I010 = 1000 + (Iuhf+1)*Nbfns2
      I020 = 1010 + (Iuhf+1)*Nbfns
      Iend = I020 + Nbfns2
       
      Call Getrec(20,"JOBARC","SCFEVECA",Nbfns2,Work(I000))
      If (Iuhf .Ne.0) Call Getrec(20,"JOBARC","SCFEVECB",NBfns2,
     +                            Work(I000+Nbfns2))

      Call Getrec(20,"JOBARC","SCFEVALA",Nbfns,Work(I010))
      If (Iuhf .Ne.0) Call Getrec(20,"JOBARC","SCFEVALB",NBfns,
     +                            Work(I010+Nbfns))

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " The alpha SCF eigenvalues"
      Write(6,"(5(1x,F15.5))") (Work(I010+i-1),i=1,Nbfns)
      Write(6,"(a)") " The alpha SCF vectors"
      Call output(Work(I000),1,Nbfns_tot,1,Nbfns,Nbfns_tot,Nbfns,1)
      If (Iuhf .Ne. 0) Then
      Write(6,"(a)") " The beta SCF eigenvalues"
        Write(6,"(5(1x,F15.5))") (Work(I010+Nbfns+i-1),i=1,Nbfns)
         Write(6,"(a)") " The beta SCF vectors"
         Call output(Work(I000+Nbfns2),1,Nbfns_tot,1,Nbfns,Nbfns_tot,
     +               Nbfns,1)
      Endif 
#endif

#ifdef _NOSKIP
C Symmetry pack the scf-eigenvectors and create occupied vectors per
C symmetry block and virtual vectors for another block.

C The occupied blocks

      Oocount = Inull
      Do Ispin = 1, Iuhf+1
         Ibegin = (Iuhf-1)*Nbfns2
         Do Irrepr = 1, Nirrep
            Istart = (Irrepr-1)*Nbfns_tot + Ibegin
            Do Iocc = 1, Pop(Irrepr,Ispin)
               Oocount = Oocount +1
               Coo(Oocount) = Work(Istart+Iocc)
            Enddo
         Enddo
      Enddo
      
C The virtual blocks

      Vvcount = Inull
      Do Ispin = 1, Iuhf+1
         Ibegin = Nocco(1)*Nbfns_tot + (Iuhf-1)*Nbfns2
         Do Irrepr = 1, Nirrep
            Istart = (Irrepr-1)*Nbfns_tot + Ibegin
            Do Ivrt = 1, Vrt(Irrepr,1)
               Vvcount = Vvcount +1
               Cvv(Vvcount) = Work(Istart+Ivrt)
            Enddo
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The occ and vrt scf-eigenvectors"
      do irrep =1,Nirrep
      Noo  = Pop(irrep,1)
      Nvv  = vrt(irrep,1)
      Write(6,"(a,I2)") " Irrep: ", Irrep
      call output(Coo(Ioff),1,noo,1,noo,noo,noo,1)
      call output(Cvv(joff),1,nvv,1,nvv,nvv,nvv,1)
      Ioff = Ioff+Noo*Noo
      Joff = Joff+Nvv*Nvv
      enddo
#endif
#endif 

      Return
      End
