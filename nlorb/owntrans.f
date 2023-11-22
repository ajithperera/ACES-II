      subroutine owntrans(trans,angmom,centbf,nbasal,nshells,
     & shells,largel,natoms,angmax,nbas)
      implicit none
 
      integer nbas, iii, angmax, angmom(nbas), degsph(angmax),
     & atom, centbf(nbas), nbasal(largel+1,natoms), largel,
     & natoms, iiio, jjj, kkk, nshells(natoms),
     & shells(largel+1,natoms), lll, lllsh

      double precision trans(nbas,nbas)

      call sphdeg(degsph,angmax)     

      iii = 1
      iiio = iii

 10   continue

      if (angmom(iii).lt.1) then

         trans(iii,iii) = 1.0D0
	 iii = iii + 1
         iiio = iii

      end if

      if (iii.gt.nbas) then
         go to 20
      end if
            
      if (angmom(iii).ge.1) then

	 atom = centbf(iii)
         do lll = 1, largel+1
            if (shells(lll,atom).eq.angmom(iii)) then
               lllsh = lll
            end if 
         end do

         do jjj = 1, degsph(angmom(iii)+1)
            kkk = iiio+nbasal(lllsh,atom)*(jjj-1)
            trans(iii,kkk) = 1.0D0
            iii = iii + 1
         end do

         iiio = iiio + 1
         if (angmom(iii-1).ne.angmom(iii)) then
            iiio = iii
         end if

      end if

      if (iii.gt.nbas) then
         go to 20
      end if
	
      go to 10

 20   continue

      return
      end
