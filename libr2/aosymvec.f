
c This routine creates a "symmetry vector" for the XX AO distributions.

      subroutine aosymvec(iAOMap,nAO)
      IMPLICIT INTEGER (A-Z)
      DIMENSION iAOMap(nAO*nAO)
      COMMON/SYMINF/NSTART,nirrep,IRREPS(255,2),dirprd(8,8)
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &             IRPDPDAO(8),IRPDPDAOS(8),ISTART(8,8),ISTARTMO(8,3)

      do irrep = 1, nirrep
         iCount = 0
         do irp_ket = 1, nirrep
            irp_bra = dirprd(irp_ket,irrep)

            bra0 = IOFFAO(irp_bra)
            bra1 = bra0 - 1 + IAOPOP(irp_bra)
            ket0 = IOFFAO(irp_ket)
            ket1 = ket0 - 1 + IAOPOP(irp_ket)

               ndx = nAO * (ket0-1)
            do ket = ket0-1, ket1-1
            do bra = bra0, bra1
               iCount = iCount + 1
               iAOMap(ndx+bra) = iCount
            end do
               ndx = ndx + nAO
            end do

         end do
      end do

      return
c     end subroutine aosymvec
      end

