        subroutine GetSDInfo(nocca,noccb,nvrta,nvrtb,nbas,vecDim,
     +                       Naobfns)
        integer nocc,nvirt,nbas,vecDim
        integer Naobfns,Nbfirr
        integer Ndrop,Ione
        Dimension Nbfirr(8)
        
#include "syminf.com"
#include "sym.com"
#include "machsp.com"
#include "info.com"

        Ione=1
        Nbas=Nocco(1)+Nvrto(1)
        Call Getrec(20,"JOBARC",'NDROPGEO',1,Idrop)

        If (Idrop .Ne. 0) Then
         Call Getrec(20,"JOBARC","NUMDROPA",1,Ndrop)
#ifdef _DEBUG_LVL0
         Write(6,"(2a,I4)"),'@SDinfo number of dropped basis', 
     +                      ' fxn:', Ndrop
#endif 
CSSS         Nbas = (Nbas - Ndrop)
        Endif

! Build and return all relevant quantities

        vecDim=Nbas*Nbas
        nocca=Nocco(1)
        noccb=Nocco(2)
#ifdef _DEBUG_LVLM
        write(6,"(a,i4)"),'@SDinfo nocca:', nocca
        write(6,"(a,i4)"),'@SDinfo noccb:', noccb
        write(6,"(a,i4)"),'@SDinfo nvrta:', nvrta
        write(6,"(a,i4)"),'@SDinfo nvrtb:', nvrtb
#endif 
        nvrta=Nvrto(1)
        nvrtb=Nvrto(2)
        Naobfns=nocc+nvrt

        end subroutine
