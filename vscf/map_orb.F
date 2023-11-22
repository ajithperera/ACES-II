      subroutine map_orb(nbas, nirrep, nocc, nbfirr, ivmlsym, imap,
     $     ispin)
c     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                                                               
      integer nocc(nirrep), nbfirr(nirrep), nbas
      integer ivmlsym(nbas), imap(nbas)
      integer drop_irp(8)
c
c  provide a mapping such that imap(iorb_sym) = imap(iorb_energ)
c  iorb_sym is the orbital ordering in correlated calculations 
C  (occ per irrep) then vrt per irrep)
c  iorb_energ is the energy ordering in scf calculation. The original routine
C  is from Marcel Nooijen. 

      if (ispin .eq. 1) then
         call getrec(20,'JOBARC','NUMDROPA', 1, ndrop)
      else
         call getrec(20,'JOBARC','NUMDROPB', 1, ndrop)
      endif
c     
c     determine number of dropped orbitals per irrep
c     
      do irp = 1, nirrep
         ii = 0
         do j = 1, ndrop
            if (ivmlsym(j) .eq. irp) ii = ii +1
c            ivlsym(j) = -1
         enddo
         drop_irp(irp) = ii
      enddo

      ii = 0
      do irp = 1, nirrep
         n1 = 1
         n2 = nocc(irp) - drop_irp(irp)
         i = 0
         do j = ndrop+1, nbas
            if (ivmlsym(j) .eq. irp) then
               i = i + 1        ! ith orbital of proper symmetry
               if (i .ge. n1 .and. i  .le. n2) then
                  ii = ii + 1
                  imap(ii) = j - ndrop
               endif
            endif
         enddo
      enddo
c     
c I have not accounted for drop vrt ...
c
      do irp = 1, nirrep
         n1 = nocc(irp)+1
         n2 = nbfirr(irp)
         i = 0
         do j = 1, nbas
            if (ivmlsym(j) .eq. irp) then
               i = i + 1        ! ith orbital of proper symmetry
               if (i .ge. n1 .and. i  .le. n2) then
                  ii = ii + 1
                  imap(ii) = j - ndrop
               endif
            endif
         enddo
c     
      enddo
c    
#ifdef _DEBUG_LVLM
         write(6,*) ' Sorting of orbitals for vee '
         write(6,*) ' iorb_sym    iorb_energ      symmetry '
         do i = 1, nbas - ndrop
            write(6,880) i, imap(i), ivmlsym(imap(i))
         enddo
 880  format(3x, i6, 5x, i6, 6x, i4)
#endif 
      if (ispin .eq. 1) then
         call putrec(20,'JOBARC', 'ORBMAP_A', nbas-ndrop, imap)
      elseif (ispin .eq. 2) then
         call putrec(20,'JOBARC', 'ORBMAP_B', nbas-ndrop, imap)
      endif
      call putrec(20,'JOBARC', 'ORBSYM  ', nbas, Ivmlsym)
c     
c     imap contains the energy ordered label for the symmetry ordered label
c     
      return
      end
