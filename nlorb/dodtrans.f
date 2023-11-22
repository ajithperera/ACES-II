      subroutine dodtrans(dtrans,dintrans,mombf,centbf,nbasal,nshell,
     & shell,largel,natoms,angmax,nbastot)
      implicit none

      integer largel, natoms, angmax, nbastot, mombf(nbastot),
     & centbf(nbastot), nbasal(largel+1,natoms), nshell(natoms),
     & shell(largel+1,natoms)

      double precision dtrans(nbastot,nbastot),
     & dintrans(nbastot,nbastot)

      call owntrans(dtrans,mombf,centbf,nbasal,nshell,shell,
     & largel,natoms,angmax,nbastot)
      call dtranspose(dtrans,dintrans,nbastot)
      call test(dtrans,dintrans,nbastot)

      return
      end
