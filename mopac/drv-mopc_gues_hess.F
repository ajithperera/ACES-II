      subroutine drv_mopc_gues_hess
C
C The interface to MOPAC program. The purpose of the interface is
C to get an approx. initial Hessian as starting point for geometry
C optimizations. Marshall Cory and Ajith Perera, 12/2004. 
C
      implicit double precision (a-h,o-z)
C
#include "machsp.com"
#include "mxatms.par"
C
      character*2 atype(111)
      character*3 ham(3)
      dimension coords(3*mxatms), IatmNum(mxatms)
C
      data bohrtoa /0.529177249D0/, IOne /1/, Izero /0/

      data atype / ' H', 'He', 'Li', 'Be', ' B', ' C', ' N', ' O', ' F',
     .             'Ne', 'Na', 'Mg', 'Al', 'Si', ' P', ' S', 'Cl', 'Ar',
     .             ' K', 'Ca', 'Sc', 'Ti', ' V', 'Cr', 'Mn', 'Fe', 'Co',
     .             'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
     .             'Rb', 'Sr', ' Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',
     .             'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', ' I', 'Xe',
     .             'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu',
     .             'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
     .             'Ta', ' W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl',
     .             'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
     .             'Pa', ' U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es',
     .             'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs',
     .             'Mt', 'Ds', 'Rg' /
C 
      call getrec(20, 'JOBARC', "NATOMS  ", IOne, Natoms)
      call getrec(20, 'JOBARC', "ATOMCHRG", Natoms, IatmNum) 
      call getrec(20, 'JOBARC', "COORD   ", 3*Natoms*IIntfp, coords)    

      open(unit=10, status='unknown', form='formatted', file='HINP')
      rewind 10
      ham(1) = 'AM1'
      ham(2) = 'UHF'
      ham(3) = 'XYZ'
C
      write(10,'(3(a,1x),a)')ham(1), ham(2), ham(3),
     .          'GEO-OK PRECISE ITRY=10000 FORCE LET'
      write(10,'(a)') ' '
      write(10,'(a)') ' '

      i = IOne
C
      do 100 Iatm=1, Natoms 
         if (IatmNum(Iatm) .Ne. IZero) Then
            write(10,'(a, 5x, 3(f20.12,1x,a))') atype(IatmNum(i)), 
     .      coords((i-1)*3+1)*bohrtoa,' 1 ',coords((i-1)*3+2)
     .      *bohrtoa,' 1 ',coords((i-1)*3+3)*bohrtoa,' 1 '
            I = I + 1 
         Else
            I = I  + 1 
         Endif
 100  continue

      write(10,'(a)') ' '
      write(10,'(a)') ' '
      write(10,'(a)') ' '

      close(unit=10, status='keep')

      return
      end
