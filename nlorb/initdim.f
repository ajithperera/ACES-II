      subroutine initdim(nbas,nbastot,nocc,nvir,natoms,
     & nelec)
      implicit none
      
      integer nbas, nbastot, nocc(2), nvir(2), natoms,
     & nelec, uhfrhf, incdummy
      
      call getrec(20,'JOBARC','NAOBASFN',1,nbas)
      call getrec(20,'JOBARC','NBASTOT ',1,nbastot)
      call getrec(20,'JOBARC','NOCCORB ',2,nocc)
      call getrec(20,'JOBARC','NVRTORB ',2,nvir)
      call getrec(20,'JOBARC','NREALATM',1,natoms)
      call getrec(20,'JOBARC','NATOMS  ',1,incdummy)
      call getrec(20,'JOBARC','UHFRHF',1,uhfrhf)

      write(*,*)
      write(*,*)

      if (uhfrhf.eq.0) then
         write(*,21) 
      end if
      if (uhfrhf.eq.1) then
         write(*,22) 
      end if

      write(*,10) nbas
      write(*,17) nbastot
      write(*,11) nocc(1)
      if (uhfrhf.eq.1) then
         write(*,19) nocc(2)
      end if
      write(*,12) nvir(1)
      if (uhfrhf.eq.1) then
         write(*,20) nvir(2)
      end if
      write(*,13) natoms
      if (uhfrhf.eq.0) then
         nelec = nocc(1)*2
         write(*,14) nelec
      end if
      if (uhfrhf.eq.1) then
         nelec = nocc(1)+nocc(2)
         write(*,14) nelec
      end if
      if (incdummy.ne.natoms) then
         write(*,23) incdummy - natoms
      end if

      write(*,*)
      write(*,*)

 10   format('Total Number of Basis Functions         = ', I6)
 11   format('Total Number of Occupied Orbitals Alpha = ', I6)
 12   format('Total Number of Virtual Orbitals Alpha  = ', I6)
 13   format('Total Number of Atoms                   = ', I6)
 14   format('Total Number of Electrons               = ', I6)
 15   format('Total Nuclear Repulsion Energy          = ', F20.12)
 17   format('Total Number of Funcs w/o contaminants  = ', I6)
 19   format('Total Number of Occupied Orbitals Beta  = ', I6)
 20   format('Total Number of Virtual Orbitals Beta   = ', I6)
 21   format('Reference function is RHF')
 22   format('Reference function is UHF')
 23   format('Total Number of Dummy Atoms             = ', I6)

      return
      end
