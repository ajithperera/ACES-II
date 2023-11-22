      subroutine addfock(myi,myj,myk,faktor,ispin,kspin,iuhf)
       implicit none
       integer i,j,myi,myj,myk,ispin,kspin,iuhf
       real *8 faktor,faktor2
#include "bwcc.com"
       integer findex
       findex(i,j)=(i*(i-1))/2+j

cjp find out whether the occupation of that spinorbital differs in reference
cjp no.1 and iref, only the contribute and contribute with appropriate sign
c
      if(iocc(myk,iref,kspin).ne.iocc0(myk,kspin)) then
             faktor2=-faktor
             if(iocc(myk,iref,kspin).ne.0) faktor2=faktor
c add the contribution
       fockcontr(findex(myi,myj),ispin)=
     &       fockcontr(findex(myi,myj),ispin)+faktor2
      endif
c
cjp sum to fock2elcontr (diagonal fock, contribution of 2el integrals)
cjp needed for energy evaluation of given fermi vacuum
cjp (we do not have the transformed 1-el integrals here, in this way
cjp dropping will be automatically processed)
c
      if(iocc(myk,iref,kspin).ne.0 .and. myi.eq.myj) then
      fock2elcontr(myi,ispin)=fock2elcontr(myi,ispin)+faktor
      endif
      if(iref.eq.1 .and.
     & iocc0(myk,kspin).ne.0 .and. myi.eq.myj) then
      fock2elcontr0(myi,ispin)=fock2elcontr0(myi,ispin)+faktor
      endif
      return
      end
