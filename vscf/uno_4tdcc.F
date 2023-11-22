      Subroutine UNO_4TDCC(Dens,Xform,Eval,Evec,Work,Maxdcor,Ldim1,
     &                     Ldim2,Iuhf,Nbas)

      Implicit Double Precision (A-H,O-Z)

#include "machsp.com"
#include "symm2.com"
#include "maxbasfn.par"

      Dimension Dens((Iuhf+1)*Ldim1), Work(Maxdcor)
      Dimension Xform(Ldim1), Eval((Iuhf+1)*Nbas)
      Dimension Evec((iuhf+1)*ldim2), Itmp(Maxbasfn)

      Logical aotran 
C
      COMMON /POPUL/  NOCC(8,2)

c calculate spatial density matrix.
c
      I000 = 1
      I010 = I000 + 2*Ldim2
      I020 = I010 + 2*Ldim2
      I030 = I020 + 2*Ldim2
      I040 = I030 + 2*Ldim2

      call vadd(Work(I000), dens(itriof(1)), dens(ldim1+itriof(1)),
     $     ldim1, 1.0d0)
c
c get eigenvectors and eigenvalues for average density
c
c The density matrix is given in AO basis. Transform to SO 
c basis before diagonalization
c
      aotran = .true.
      ioff = 1
      do irp = 1, nirrep
         if (nbfirr(irp).eq. 0) go to 100
         nsize=nbfirr(irp)
         call expnd2(Work(I000+itriof(irp)-1), Work(I010), nsize)
         if (aotran) then
c
            call expnd2(xform(itriof(irp)), Work(I030), nsize)
c
c calculate s^(-1/2) by inverting xform. Does this work if 
c linear dependencies? Note that Xform is S^(1/2).
c
            call zero(Work(I020),nsize*nsize)
            do i = 1, nsize
               Work(I020 + (i-1)*nsize + i - 1)= 1.0d0
            enddo
            info = 0

#ifdef _DEBUG_LVLM
      write(6,"(a)") "The total expanded density "
      call output(Work(I000), 1, nsize, 1, nsize, nsize, nsize, 1)
#endif
            call dgesv(nsize, nsize, Work(I030), nsize, Itmp,
     $                 Work(I020), nsize, info)
            if (info .ne. 0) then
               write(6,*) ' @uno_makerhf:something wrong inv. xform'
               call errex
            endif

#ifdef _DEBUG_LVLM
        write(6,*) ' checking inverse'
        I050 = I040 + 2*LDIM2
        call expnd2(xform(itriof(irp)), Work(I040), nsize)
        
        call xgemm('N', 'N', nsize, nsize, nsize, 1.0d0,
     $              Work(I020), nsize, Work(I040), nsize, 0.0d0, 
     $              Work(I050), nsize)
        call output(Work(I050), 1, nsize, 1, nsize, nsize, nsize, 1)
#endif 
            call xgemm('N', 'T', nsize, nsize, nsize, 1.0d0,
     $           Work(I010), nsize, Work(I020), nsize, 0.0d0,
     $           Work(I030), nsize)
            call xgemm('N', 'N', nsize, nsize, nsize, 1.0d0, 
     $           Work(I020), nsize, Work(I030), nsize, 0.0d0,
     $           Work(I010), nsize)
            trace = 0.0
            do i = 1, nsize
               trace = trace + Work(I010 +(i-1)*nsize+i-1)
            enddo

#ifdef _DEBUG_LVLM
      write(6,"(a,F15.8)")' @uno_makerhf: trace of transformed dens ' ,
     $             trace
#endif
         endif

         call eig(Work(i010), Work(i020), nsize, nsize, 0)
c
c copy eigenvalues (obtained on diagonal of densscr) in eval
c and eigenvectors in evec
c
         icount = 1
         do i = 1, nsize
            eval(i+ioff-1) = Work(I010 + (icount-1))
            icount = icount + nsize + 1
         enddo
c
         call dcopy(nsize*nsize, Work(I020), 1, evec(isqrof(irp)), 1)
c
c
#ifdef _DEBUG_LVLM
            write(6,*) ' @uno_makerhf: symmetry block :', irp
            write(6,*) ' natural occupation numbers:'
            call output(eval(ioff), 1, 1, 1, nsize, 1, nsize, 1)
            trace = 0.0d0
            do i = 1, nsize
               trace = trace + eval(ioff+i-1)
            enddo
            write(6,*) ' trace of density matrix block', trace
            write(6,*) ' corresponding symmetry blocked natural',
     $           ' orbitals'
            call output(evec(isqrof(irp)), 1, nsize, 1, nsize,
     $           nsize, nsize, 1)
#endif
c
         ioff = ioff + nsize
c
 100  enddo
c
      factor = 1.0d0
      NoccA  = 0
      call zero(Work(I010),2*ldim2)

      Do irrp = 1, Nirrep
         NoccA = NoccA + NOCC(irrp, 1)
      Enddo
#ifdef _DEBUG_LVLM
      Write(6,"(a,I4)") " The number of occ. orbitals :", NoccA
#endif 
      do iorb = 1, NoccA

         amax = -0.5
         icount = 0
         do irp = 1, nirrep
            do i = 1, nbfirr(irp)
               icount = icount + 1
               if (eval(icount) .gt. amax) then
                  irp0 = irp
                  icount0 = icount
                  i0 = i
                  amax = eval(icount)
               endif
            enddo
         enddo
C
c add factor * c(p, icount) * c_dagger(q,icount) to density matrix
c
         eval(icount0) = -2.0d0
         icount = isqrof(irp0) + (i0-1) * nbfirr(irp0)
         nsize  = nbfirr(irp0)

         call xgemm('N', 'T', nsize, nsize, 1, factor,
     $        evec(icount), nsize, evec(icount), nsize,
     $        1.0d0, Work(I010+isqrof(irp0)-1), nsize)
      enddo
C
c  duplicate alpha to beta density
c
      do i = 1, ldim2
         Work(I010+ldim2+i-1) = Work(I010+i-1)
      enddo

#ifdef _DEBUG_LVLM
      do irp = 1, nirrep
         nsize = nbfirr(irp)
         if (nsize .gt. 0) then
            tracea = 0.0d0
            traceb = 0.0d0
            do i = 1, nsize
               tracea = tracea +
     $                  Work(i010+isqrof(irp)+(i-1)*nsize+i-2)
            enddo
            do i = 1, nsize
                traceb = traceb + Work(i010+ldim2+
     $                   isqrof(irp)+(i-1)*nsize+i-2)
            enddo
            write(6,"(a,2F15.8)")' Partial traces new density matrix ',
     $                  tracea,traceb
         endif
      enddo
#endif
c
c new density matrix is formed in orthogonal basis. transform to ao basis
c
      if (aotran) then
         do irp = 1, nirrep
            if (nbfirr(irp).ne. 0) then
               nsize=nbfirr(irp)
               call expnd2(xform(itriof(irp)), Work(I020), nsize)

c alpha spin
               call xgemm('N', 'T', nsize, nsize, nsize, 1.0d0,
     $              Work(i010+isqrof(irp)-1), nsize, Work(I020), 
     $              nsize, 0.0d0, Work(I030), nsize)
               call xgemm('N', 'N', nsize, nsize, nsize, 1.0d0,
     $              Work(I020), nsize, Work(I030), nsize, 0.0d0, 
     $              Work(I010+isqrof(irp)-1), nsize)
c beta spin
               call xgemm('N', 'T', nsize, nsize, nsize, 1.0d0,
     $              Work(i010+isqrof(irp)-1+ldim2), nsize, Work(I020), 
     $              nsize, 0.0d0, Work(I030), nsize)
               call xgemm('N', 'N', nsize, nsize, nsize, 1.0d0,
     $              Work(I020), nsize, Work(I030), nsize, 0.0d0,
     $              Work(I010+isqrof(irp)-1+ldim2), nsize)
C
            Endif 
         Enddo
      Endif 

      Do Irrp = 1, Nirrep 
         if (nbfirr(irrp).ne. 0) then

C Alpha Density

             Call squez2(Work(I010+isqrof(irrp)-1), 
     &                   Dens(Itriof(irrp)),nbfirr(irrp))

C Beta Density

             Call squez2(Work(I010+isqrof(irrp)-1+ldim2), 
     &                   Dens(Ldim1+Itriof(irrp)),nbfirr(irrp))
             
         Endif 
      Enddo

#ifdef _DEBUG_LVLM
      Write(6,*)
      Do Irrp = 1, Nirrep
         if (nbfirr(irrp).ne. 0) then
         suma = 0.d0
         sumb = 0.d0
         nsize = Itriln(Irrp)
         do i = 1, nsize 
         suma = suma + Dens((Itriof(irrp)+i-1)) * 
     &                 Dens((Itriof(irrp)+i-1))
         sumb = sumb + Dens(Ldim1+(Itriof(irrp)+i-1)) * 
     &                 Dens(Ldim1+(Itriof(irrp)+i-1))
         enddo
         Write(6,"(a,2F15.8)") " Checksum of Alpha/Beta density = ", 
     &                           suma, sumb
         Endif 
      Enddo 
#endif 
C
      Return
      End


