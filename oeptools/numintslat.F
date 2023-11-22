      subroutine numintslat (densao, slatao)
c
c This subroutines calculates matrix elements of the Slater           
c orbital-dependent potential by numerical integration over a grid.   
c
c It's restricted to RHF and C1 symmetry
c
c Igor Schweigert, 2003 $Id: numintslat.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c
#include <aces.h>
#include <const.par>
#include <intgrtflags.com>
#include <flags.com>
#include <kscore.com>
#include <mol.com>
#include <grid.com>
#include <int.com>
#include <dftfunc.com>
c#include <energ.com>
c#include <io.par>
c#include <machsp.com>
c
c     .. input parameters ..
c
      double precision
     &     densao (nbastot, nbastot), slatao (nbastot, nbastot)
c
c     .. local variables ..
c
      integer
     &     i, iunq, nunq, iatm, nredatm
c
      double precision
     &     grdwgt, grdpnt (3), denval, sltval, totden
c      
      integer 
     &     zirrtmp, zradpt, zrwt, zrelsiz,
     &     zwtintr, zrtmp, zrint, setptr,
     $     iradpt, iangpt, grid, zcdnt, zrsqrd, zvalprim, 
     &     znull, pnull, zvalao,  
     $     znucatt, zscraxp, zscraxa, zscr1xa
c
      logical evalpt
c
c     .. the offset for the local data ..
c
      znull  =setptr(1,F_REAL,    1)
      pnull  =setptr(1,F_INTEGER, 1)
c
c     .. allocate memory for the local data ..
c
      zvalao    = setptr(1,F_REAL, nbastot)
      zvalprim  = setptr(1,F_REAL, totprim)
      znucatt   = setptr(1,F_REAL, totprim*totprim)
c
c     .. arrays needed for numerical integration grids                ..
c
      zrsqrd    = setptr(1,F_REAL, natoms)
      zrwt      = setptr(1,F_REAL, int_numradpts*maxangpts)
      zirrtmp   = setptr(1,F_REAL, natoms)
      zradpt    = setptr(1,F_REAL, int_numradpts)
      zrelsiz   = setptr(1,F_REAL, natoms*natoms)
      zrtmp     = setptr(1,F_REAL, natoms)
      zrint     = setptr(1,F_REAL, natoms*maxangpts)
      zwtintr   = setptr(1,F_REAL, natoms)
      zcdnt     = setptr(1,F_REAL, natoms*3)
c
c     .. scratch arrays ..
c
      zscraxp   = setptr(1,F_REAL, nbastot*totprim)
      zscraxa   = setptr(1,F_REAL, nbastot*nbastot)
      zscr1xa   = setptr(1,F_REAL, nbastot)
c
c     .. determine relative atomic size for the atomic partitioning   ..
c     .. The relative sizes of atoms are used in INTGRT/oct.F         ..
c     .. to determine the radial integration limits in the part that  ..
c     .. is commented out in OCT.F So unless this "rigid partioning"  ..
c     .. is implemented c the relative sizes are not needed           ..
c
      call dzero (dcore (zrelsiz), natoms*natoms)
c     
c     .. integrate only over symmetry unique atoms                    ..
c
      nunq = compnorb
      if (c1symmet.eq.1) nunq = natoms
      totden = 0.d0
c
      do iunq = 1, nunq
c
c        .. pick up an unique atom ..
c
         iatm = kscore(polist+iunq-1)
         if (c1symmet.eq.1) iatm = iunq
         icntr = iatm
c
c        .. determine the number of symmetry equivalent atoms         ..

         nredatm = kscore(pcomppopv+iunq-1)
         if (c1symmet.eq.1) nredatm = 1
c
c        .. set up the grid for integration over symm. uniq. octants ..
c
         call oct(
     &        kscore(patomchrg), dcore(zatmvc),
     &        dcore(zrij),      kscore(pradgrid),
     &        dcore(zrelsiz),   dcore(zrtmp),
     &        dcore(zrint),     dcore(zradpt),
     &        dcore(zrwt),      dcore(zbslrd),
     &        dcore(zgridxyz),  kscore(pgrdangpts))
c
c        .. loop over radial points ..
c
         do iradpt=1,int_numradpts
c
c           .. loop over angular points ..
c
            grid=kscore(pradgrid+iradpt-1)
            do iangpt=1,kscore(pgrdangpts+grid-1)
c
c              .. determine if the curr. grid point is included,      ..
c              .. its cartesian coordinates, and integration weigth   ..
c
               evalpt=.true.
               call symoct(
     &              kscore(patomchrg), dcore(zatmvc),
     &              dcore(zrij),      dcore(zaij),
     &              dcore(zcdnt),     dcore(zrsqrd),
     &              dcore(zirrtmp),   dcore(zwtintr),
     &              grdwgt,           kscore(pradgrid),
     &              dcore(zradpt),    dcore(zrwt),
     &              kscore(pintegaxis),dcore(zgridxyz),
     &              dcore(zgridwt),   kscore(pgrdangpts),
     &              iradpt,           iangpt,
     &              grid,             evalpt)
c
c              .. if this point is not included, skip the integration ..
c
               if ( .not. evalpt) goto 111
c
c              .. scale this point's weight by the number of symmetry ..
c              .. equivalent atoms                                    ..
c
               grdwgt = dble (nredatm) * grdwgt
c
c               .. calculate values of the primitives at the grid point ..
c
               call gridprmval (
     &              natoms, totprim, maxshell, numangfct, 
     &              kscore (pnshellatom), kscore (pnshellprim),
     &              kscore (pangfct),  dcore (zalpha), 
     &              dcore (zcdnt), dcore (zvalprim) )
c
c              .. calculate the value of the density at the grid point ..
c
               call griddenval (
     &              totprim, nbastot, dcore (zvalprim), 
     &              dcore (zprimcoef), densao, dcore (zvalao), denval)
c
c              .. check the density value ..
c     
               if (denval .lt. thresh) goto 111
               totden = totden + denval*grdwgt
c
c             .. absolute coordinates of the grid point needed for nuc.att...
c     
               do i=1, 3
                  grdpnt (i) =
     $                 dcore (zcdnt + (i-1) * natoms) +
     $                 dcore (zcoord + i-1)
               enddo
c
c             .. calculate the nuclear attraction integrals ..
c
               call nucatt (
     &              natoms, totprim, maxshell, numangfct, 
     &              kscore (pnshellatom), kscore (pnshellprim),
     $              kscore (pangfct), dcore (zalpha), 
     &              grdpnt, dcore (zcoord), dcore (znucatt))
c
c             .. calculate the value of the Slater potential at the point ..
c
               call gridslat (
     $              totprim, nbastot, 
     $              dcore (zprimcoef), dcore (zvalao), densao, denval,  
     $              dcore (znucatt), dcore (zscr1xa), dcore (zscraxp),
     $              dcore (zscraxa), sltval)
c
c              .. add this grid point's contribution to the Slater AO ME ..
c
               call xgemm(
     &              'n', 't', nbastot, nbastot, 1,
     &              sltval*grdwgt, dcore(zvalao), nbastot,
     &              dcore(zvalao), nbastot,
     &              1.d0, slatao, nbastot)
c
 111           continue
            enddo
         enddo
      enddo
c
c      write (6,'(a,f20.12)') 'Density integrates to ', totden
c
c     .. done ..
c
      call relptr(1,F_REAL,   znull)
      call relptr(1,F_INTEGER,pnull)
      return
      end
