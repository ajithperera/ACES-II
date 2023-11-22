      subroutine sortmom(centbf,mombf,nshells,shells,nbasal,
     & largel,mxnal,natoms,nbas,angmax,junk,count,count2)
      implicit none

      integer angmax 
 
      integer centbf(nbas), mombf(nbas), nshells(natoms),
     & shells(largel+1,natoms), nbas, iii, jjj, kkk, lll,
     & temp, junk(nbas), count(natoms,angmax), degen(angmax),
     & nbasal(largel+1,natoms), mxnal, largel, natoms,
     & count2(natoms,angmax)

      mxnal = 0
      temp = 0
      kkk = 0
      call izerosec(nshells,natoms)
      call izerosec(junk,nbas)
      call sphdeg(degen,angmax) 

      call countmom(natoms,nbas,centbf,mombf,angmax,count) 
      call countsh(natoms,angmax,count,count2,degen) 

      do iii = 1, natoms
         do jjj = 1, nbas 
            if (centbf(jjj).eq.iii) then
               
               junk(jjj) = mombf(jjj) + 1

               if (junk(jjj).ne.temp) then
                  kkk = kkk + 1
                  nshells(iii) = nshells(iii) + 1
                  shells(kkk,iii) = mombf(jjj)
                  nbasal(kkk,iii) = count2(iii,mombf(jjj) + 1)
                  temp = junk(jjj)
               end if

            end if  
         end do
 
         do kkk = 1, nshells(iii)
            if (nbasal(kkk,iii).gt.mxnal) then
               mxnal = nbasal(kkk,iii)
            end if
         end do
 
         write(*,10) iii, nshells(iii) 
         do kkk = 1, nshells(iii)
            write(*,20) kkk, shells(kkk,iii), nbasal(kkk,iii)
         end do
         write(*,*)
         temp = 0
         kkk = 0

      end do 

 10   format('Atom ',I5,' has ',I5,' shells')
 20   format('Shell ',I5,' has angular momentum ',I5
     & , ' of dimension ',I5)
 
      return
      end
	  
